import re
import yaml
from pathlib import Path
import sys
import requests
from thefuzz import fuzz
import json
from typing import Optional, Dict, Any, List, Tuple

# --- Configuration ---
YAML_FILE_NAME = "../api_fallback_resolutions.yaml"
CONTEXT_RESOLUTIONS_KEY = "context_specific_resolutions"
GLOBAL_NAME_RESOLUTIONS_KEY = "global_name_resolutions"
RESOLUTIONS_KEY_OLD = "resolved_mapping"
NO_MATCH_MARKER = "__NO_MATCH_CONFIRMED__"

# Global cache for SwissLipids API responses
api_slp_cache: Dict[str, Any] = {}

# --- Regexes (Identical) ---
REGEX_BEST_MATCH = re.compile(
    r"INFO - API Fallback: Best match for '(?P<name>[^']*)' is API_NAME '(?P<api_name>[^']*)' "
    r"\(SLM: (?P<slm_id>[^)]*)\) with score (?P<score>\d+) "
    r"\(SLP (?P<slp_id>[^,]*), Rhea (?P<rhea_id>[^)]*)\)\."
)
REGEX_NO_MATCH = re.compile(
    r"INFO - API Fallback: No matches >= \d+ for '(?P<name>[^']*)' "
    r"in SLP (?P<slp_id>[^,]*), Rhea (?P<rhea_id>[^)]*)\."
)
REGEX_AMBIGUITY = re.compile(
    r"WARNING - API Fallback AMBIGUITY for '(?P<name>[^']*)' "
    r"\(SLP (?P<slp_id>[^,]*), Rhea (?P<rhea_id>[^)]*)\):"
)
REGEX_FAILED = re.compile(
    r"WARNING - API Fallback FAILED for (?:reactant|product) '(?P<name>[^']*)' "
    r"with SLP (?P<slp_id>[^,]*), Rhea (?P<rhea_id>[^)]*)\."
)


def _strip_html(text: str) -> str:
    clean = re.compile("<.*?>")
    return re.sub(clean, "", text)


def resolve_lipid_via_swisslipids_api(
    unresolved_lipid_name: str,
    reaction_context_slp_id: str,
    reaction_context_rhea_id_str: Optional[str],
) -> Optional[str]:
    global api_slp_cache
    api_base_url = "https://www.swisslipids.org/api/index.php/entity/"
    try:
        if reaction_context_slp_id in api_slp_cache:
            slp_data = api_slp_cache[reaction_context_slp_id]
        else:
            print(
                f"INFO: (Live API Call) Fetching {api_base_url + reaction_context_slp_id}"
            )
            response = requests.get(api_base_url + reaction_context_slp_id, timeout=10)
            response.raise_for_status()
            slp_data = response.json()
            api_slp_cache[reaction_context_slp_id] = slp_data
    except requests.exceptions.RequestException as e:
        print(
            f"WARNING: (Live API Call) API request failed for SLP ID {reaction_context_slp_id}: {e}"
        )
        return None
    except json.JSONDecodeError as e:
        print(
            f"WARNING: (Live API Call) Failed to decode JSON from API for SLP ID {reaction_context_slp_id}: {e}"
        )
        return None
    if not slp_data or "reactions" not in slp_data or not slp_data["reactions"]:
        print(
            f"DEBUG: (Live API Call) No reactions found in API response for SLP ID {reaction_context_slp_id}"
        )
        return None
    MIN_ACCEPTABLE_SCORE = 90
    target_rhea_id_int = None
    if reaction_context_rhea_id_str and reaction_context_rhea_id_str.strip().isdigit():
        target_rhea_id_int = int(reaction_context_rhea_id_str.strip())
    candidate_molecules_in_reaction = []
    found_target_reaction = False
    if target_rhea_id_int is not None:
        for _slcr_id, reaction_data_wrapper in slp_data["reactions"].items():
            reaction_details = reaction_data_wrapper.get("reaction", {})
            rhea_info = reaction_details.get("rhea", {})
            api_rhea_id_str = rhea_info.get("rhea_id")
            if (
                api_rhea_id_str
                and api_rhea_id_str.isdigit()
                and int(api_rhea_id_str) == target_rhea_id_int
            ):
                found_target_reaction = True
                for participant_type in ["rhea_reactants", "rhea_products"]:
                    participants = rhea_info.get(participant_type, {})
                    for slm_id, molecule_info in participants.items():
                        api_molecule_name_raw = molecule_info.get("name", "")
                        api_molecule_name_clean = _strip_html(
                            api_molecule_name_raw
                        ).strip()
                        if api_molecule_name_clean:
                            candidate_molecules_in_reaction.append(
                                (slm_id, api_molecule_name_clean)
                            )
                break
        if not found_target_reaction:
            print(
                f"DEBUG: (Live API Call) Target Rhea ID {target_rhea_id_int} not found for SLP {reaction_context_slp_id}."
            )
            return None
    else:
        print(
            f"DEBUG: (Live API Call) No Rhea ID context for SLP {reaction_context_slp_id}. Skipping for '{unresolved_lipid_name}'."
        )
        return None
    if not candidate_molecules_in_reaction:
        print(
            f"DEBUG: (Live API Call) No candidates from API for SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int} for '{unresolved_lipid_name}'."
        )
        return None
    best_match_slm_id = None
    highest_score = -1
    ambiguous_matches = []
    for slm_id, api_name in candidate_molecules_in_reaction:
        score = fuzz.ratio(unresolved_lipid_name.lower(), api_name.lower())
        print(
            f"DEBUG: (Live API Call) Fuzzy: '{unresolved_lipid_name}' vs API_NAME '{api_name}' (SLM: {slm_id}) -> score {score}"
        )
        if score >= MIN_ACCEPTABLE_SCORE:
            if score > highest_score:
                highest_score = score
                best_match_slm_id = slm_id
                ambiguous_matches = [(score, slm_id, api_name)]
            elif score == highest_score:
                ambiguous_matches.append((score, slm_id, api_name))
    if not best_match_slm_id:
        print(
            f"INFO: (Live API Call) No matches >= {MIN_ACCEPTABLE_SCORE} for '{unresolved_lipid_name}' in SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int}."
        )
        return None
    if len(ambiguous_matches) > 1:
        unique_slm_ids_in_ambiguity = set(item[1] for item in ambiguous_matches)
        if len(unique_slm_ids_in_ambiguity) > 1:
            print(
                f"WARNING: (Live API Call) AMBIGUITY for '{unresolved_lipid_name}' (SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int}): "
                f"Multiple distinct SLM IDs score {highest_score}. Matches: {ambiguous_matches}. Skipping."
            )
            return None
        else:
            print(
                f"INFO: (Live API Call) Multiple representations for SLM ID '{best_match_slm_id}' score {highest_score} for '{unresolved_lipid_name}'. "
                f"Selected SLM: {best_match_slm_id}."
            )
            return best_match_slm_id
    print(
        f"INFO: (Live API Call) Best match for '{unresolved_lipid_name}' is API_NAME '{ambiguous_matches[0][2]}' (SLM: {best_match_slm_id}) "
        f"score {highest_score} (SLP {reaction_context_slp_id}, Rhea {target_rhea_id_int})."
    )
    return best_match_slm_id


def load_resolutions(yaml_path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    context_resolutions = {}
    global_name_resolutions = {}
    if yaml_path.exists():
        try:
            with open(
                yaml_path, "r", encoding="utf-8"
            ) as f:  # Ensure reading with UTF-8
                data = yaml.safe_load(f)
                if data:
                    if RESOLUTIONS_KEY_OLD in data and isinstance(
                        data[RESOLUTIONS_KEY_OLD], dict
                    ):
                        context_resolutions = data.get(RESOLUTIONS_KEY_OLD, {})
                        print(
                            f"INFO: Loaded data from old key '{RESOLUTIONS_KEY_OLD}'. Will save with new key '{CONTEXT_RESOLUTIONS_KEY}'."
                        )
                    else:
                        context_resolutions = data.get(CONTEXT_RESOLUTIONS_KEY, {})
                    global_name_resolutions = data.get(GLOBAL_NAME_RESOLUTIONS_KEY, {})
        except yaml.YAMLError:
            print(
                f"Warning: Could not parse {yaml_path}. Starting with empty resolutions."
            )
        except Exception as e:
            print(
                f"Warning: Error loading {yaml_path}: {e}. Starting with empty resolutions."
            )
    return context_resolutions, global_name_resolutions


def save_resolutions(
    yaml_path: Path,
    context_resolutions: Dict[str, str],
    global_name_resolutions: Dict[str, str],
):
    try:
        data_to_save = {
            CONTEXT_RESOLUTIONS_KEY: context_resolutions,
            GLOBAL_NAME_RESOLUTIONS_KEY: global_name_resolutions,
        }
        # Diagnostic print for one of the problematic keys if it exists
        # test_key_part = "β-D-galactosyl-(1→3)-N-acetyl-β-D-galactosaminyl"
        # for key in context_resolutions.keys():
        #     if test_key_part in key:
        #         print(f"DEBUG SAVE: Key being saved: '{key}' (type: {type(key)})")
        # for key in global_name_resolutions.keys():
        #      if test_key_part in key:
        #         print(f"DEBUG SAVE: Global Key being saved: '{key}' (type: {type(key)})")

        with open(yaml_path, "w", encoding="utf-8") as f:  # Ensure writing with UTF-8
            yaml.dump(
                data_to_save,
                f,
                sort_keys=True,
                default_flow_style=False,
                allow_unicode=True,
            )
    except IOError as e:
        print(f"Error: Could not save resolutions to {yaml_path}. {e}")
    except Exception as e:
        print(f"An unexpected error occurred during saving: {e}")


def parse_log_for_fallback_cases(log_content: str) -> List[Dict[str, Any]]:
    cases = {}
    for line_num, line in enumerate(log_content.splitlines()):
        match_best = REGEX_BEST_MATCH.search(line)
        match_no = REGEX_NO_MATCH.search(line)
        match_ambi = REGEX_AMBIGUITY.search(line)
        match_fail = REGEX_FAILED.search(line)
        name, slp, rhea = None, None, None
        log_slm, log_api_name, log_score = None, None, None
        status = None
        if match_best:
            data = match_best.groupdict()
            name, slp, rhea = data["name"], data["slp_id"], data["rhea_id"]
            log_slm, log_api_name, log_score = (
                data["slm_id"],
                data["api_name"],
                int(data["score"]),
            )
            status = "BEST_MATCH"
        elif match_no:
            data = match_no.groupdict()
            name, slp, rhea = data["name"], data["slp_id"], data["rhea_id"]
            status = "NO_MATCH_LOG"
        elif match_fail:
            data = match_fail.groupdict()
            name, slp, rhea = data["name"], data["slp_id"], data["rhea_id"]
            status = "FAILED_LOG"
        elif match_ambi:
            data = match_ambi.groupdict()
            name, slp, rhea = data["name"], data["slp_id"], data["rhea_id"]
            status = "AMBIGUITY_LOG"
        if name and slp and rhea:
            rhea_clean = rhea.strip() if rhea else None
            if rhea_clean == "None":
                rhea_clean = None
            case_key = (name, slp, rhea_clean if rhea_clean else "N/A")
            current_entry = {
                "name": name,
                "slp_id": slp,
                "rhea_id": rhea_clean,
                "log_slm_id": log_slm,
                "log_api_name": log_api_name,
                "log_score": log_score,
                "log_status": status,
                "line_num": line_num + 1,
            }
            if case_key not in cases or (
                status == "BEST_MATCH" and cases[case_key]["log_status"] != "BEST_MATCH"
            ):
                cases[case_key] = current_entry
            elif cases[case_key]["log_status"] != "BEST_MATCH" and status is not None:
                if status == "BEST_MATCH":  # Prioritize best match if it appears later
                    cases[case_key] = current_entry
                # else, keep the first non-best_match status or update if a more "resolved" status comes
                elif cases[case_key]["log_status"] in [
                    "NO_MATCH_LOG",
                    "FAILED_LOG",
                    "AMBIGUITY_LOG",
                ] and status in [
                    "BEST_MATCH"
                ]:  # Should be covered by above, but defensive
                    cases[case_key] = current_entry
    sorted_cases = sorted(list(cases.values()), key=lambda x: x["line_num"])
    return sorted_cases


def main():
    yaml_file = Path(YAML_FILE_NAME)
    context_resolutions, global_name_resolutions = load_resolutions(yaml_file)

    while True:
        log_file_path_str = input(
            "Enter the path to your log file (or type 'exit' to quit): "
        ).strip()
        if log_file_path_str.lower() == "exit":
            print("Exiting program.")
            sys.exit(0)
        log_file_path = Path(log_file_path_str)
        if log_file_path.is_file():
            break
        else:
            print(f"Error: Log file not found at '{log_file_path}'. Please try again.")

    try:
        with open(
            log_file_path, "r", encoding="utf-8"
        ) as f:  # Ensure reading with UTF-8
            log_content = f.read()
    except Exception as e:
        print(f"Error reading log file '{log_file_path}': {e}")
        sys.exit(1)

    all_fallback_cases = parse_log_for_fallback_cases(log_content)
    new_cases_to_review = []
    auto_applied_global_count = 0
    already_resolved_context_count = 0
    updated_by_global_count = 0

    for case in all_fallback_cases:
        molecule_name = case["name"]
        rhea_for_key = case["rhea_id"] if case["rhea_id"] else "N/A"
        case_composite_key = f"{molecule_name}::{case['slp_id']}::{rhea_for_key}"

        # Ensure the molecule_name itself is a clean unicode string from the start
        if not isinstance(molecule_name, str):  # Should already be str from regex
            molecule_name = str(molecule_name)

        if molecule_name in global_name_resolutions:
            global_decision = global_name_resolutions[molecule_name]
            if (
                case_composite_key not in context_resolutions
                or context_resolutions[case_composite_key] != global_decision
            ):
                print(
                    f"INFO: Auto-applying global decision for '{molecule_name}' ({global_decision}) to context {case['slp_id']}/{rhea_for_key}."
                )
                context_resolutions[case_composite_key] = global_decision
                updated_by_global_count += 1
            else:  # Already consistent with global
                auto_applied_global_count += 1
            continue

        if case_composite_key in context_resolutions:
            already_resolved_context_count += 1
            continue

        new_cases_to_review.append(case)

    if updated_by_global_count > 0:
        print(
            f"\n{updated_by_global_count} context-specific case(s) were updated by existing global name resolutions."
        )
        save_resolutions(
            yaml_file, context_resolutions, global_name_resolutions
        )  # Save these auto-updates
    if auto_applied_global_count > 0:
        print(
            f"\n{auto_applied_global_count} case(s) were already consistent with global name resolutions."
        )
    if already_resolved_context_count > 0:
        print(
            f"\n{already_resolved_context_count} case(s) were already resolved for their specific context (and no global rule applied)."
        )

    if not new_cases_to_review:
        print("\nNo new API fallback cases found in the log to review interactively.")
    else:
        print(
            f"\nFound {len(new_cases_to_review)} new API fallback case(s) for interactive review."
        )
        print("-----------------------------------------------------------------")

    try:
        for i, case in enumerate(new_cases_to_review):
            molecule_name = case["name"]
            rhea_for_key = case["rhea_id"] if case["rhea_id"] else "N/A"
            case_composite_key = f"{molecule_name}::{case['slp_id']}::{rhea_for_key}"

            print(
                f"\n--- Reviewing Case {i+1} of {len(new_cases_to_review)} (Log Line: {case['line_num']}) ---"
            )
            print(
                f"  Unresolved Name: '{molecule_name}'"
            )  # Python print handles unicode fine
            print(f"  SLP Context:     {case['slp_id']}")
            print(f"  Rhea Context:    {case['rhea_id'] if case['rhea_id'] else 'N/A'}")
            print(f"  Status in Log:   {case['log_status']}")
            if case["log_slm_id"]:
                print(
                    f"  Logged SLM ID:   {case['log_slm_id']} ('{case['log_api_name']}') with score {case['log_score']}"
                )

            chosen_slm_id_for_yaml = None
            user_made_decision_this_loop = False

            while not user_made_decision_this_loop:
                print("\n  Actions:")
                action_prompt_parts = []
                if case["log_slm_id"] and case["log_status"] == "BEST_MATCH":
                    action_prompt_parts.append("(a)ccept logged")
                action_prompt_parts.extend(
                    [
                        "(q)uery API live",
                        "(m)anual SLM ID",
                        "(n)o match (globally for this name)",
                        "(s)kip session",
                        "(e)xit",
                    ]
                )
                action_prompt = "  " + ", ".join(action_prompt_parts) + ": "

                choice = input(action_prompt).strip().lower()

                if (
                    choice == "a"
                    and case["log_slm_id"]
                    and case["log_status"] == "BEST_MATCH"
                ):
                    chosen_slm_id_for_yaml = case["log_slm_id"]
                    print(f"    -> Accepted logged SLM ID: {chosen_slm_id_for_yaml}")
                    user_made_decision_this_loop = True
                elif choice == "q":
                    print("    -> Running live API query...")
                    live_slm_id = resolve_lipid_via_swisslipids_api(
                        molecule_name, case["slp_id"], case["rhea_id"]
                    )
                    if live_slm_id:
                        print(
                            f"    INFO: Live API query suggested SLM ID: {live_slm_id}"
                        )
                        sub_choice = (
                            input("      Use this SLM ID? (y/n): ").strip().lower()
                        )
                        if sub_choice == "y":
                            chosen_slm_id_for_yaml = live_slm_id
                            user_made_decision_this_loop = True
                    # If 'n' or no live_slm_id, loop continues for new action
                elif choice == "m":
                    manual_id = input(
                        "    Enter SLM ID (e.g., SLM:xxxx or CHEBI:xxxx): "
                    ).strip()
                    if manual_id:
                        chosen_slm_id_for_yaml = manual_id
                        print(f"    -> Manually set SLM ID: {chosen_slm_id_for_yaml}")
                        user_made_decision_this_loop = True
                    else:
                        print("    Invalid SLM ID entered.")
                elif choice == "n":
                    chosen_slm_id_for_yaml = NO_MATCH_MARKER
                    print(
                        f"    -> Marked as '{NO_MATCH_MARKER}' globally for '{molecule_name}'."
                    )
                    user_made_decision_this_loop = True
                elif choice == "s":
                    print(f"    -> Skipped '{molecule_name}' for this session.")
                    user_made_decision_this_loop = True
                elif choice == "e":
                    raise KeyboardInterrupt
                else:
                    print("    Invalid choice. Please try again.")

            if chosen_slm_id_for_yaml is not None:
                # Ensure molecule_name is a str for YAML key
                molecule_name_key = str(molecule_name)

                context_resolutions[case_composite_key] = chosen_slm_id_for_yaml
                global_name_resolutions[molecule_name_key] = chosen_slm_id_for_yaml

                # print(f"DEBUG: Saving context_key='{case_composite_key}', global_key='{molecule_name_key}', value='{chosen_slm_id_for_yaml}'")
                save_resolutions(
                    yaml_file, context_resolutions, global_name_resolutions
                )
                print(
                    f"    Progress saved for '{molecule_name}' (context-specific and global)."
                )

        if new_cases_to_review:
            print("\n-----------------------------------------------------------------")
            print(
                "All new API fallback cases from this log have been processed or skipped for this session."
            )

    except KeyboardInterrupt:
        print("\n\nReview session interrupted by user.")

    finally:
        final_context_res, final_global_res = load_resolutions(yaml_file)
        print(f"\nSummary:")
        print(f"  Resolution decisions saved in: {yaml_file.resolve()}")
        print(f"  Total context-specific decisions: {len(final_context_res)}")
        print(f"  Total global name decisions: {len(final_global_res)}")


if __name__ == "__main__":
    main()
