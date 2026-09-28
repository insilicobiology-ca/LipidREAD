import csv
import itertools
import math
import shutil
import re

import pandas as pd
import psycopg2
import requests
from zipfile import ZipFile
import gzip
from tqdm import tqdm
from _archive.deprecated.compli_utils import get_family
from _archive.deprecated.compli_utils import get_family_generalization
from _archive.deprecated.clean_up import generate_reaction_combinations
from _archive.deprecated.clean_up import split_reaction_text
from _archive.deprecated.clean_up import clean_abbreviation

# import matching_algo
# from matching_algo import get_family

# import matching_algo as m_algo

# conn = psycopg2.connect(
#     "dbname='enznet_db' user='lipidnetwork_user' host='localhost' password='lipidnetwork_pass'")
# cur = conn.cursor()

conn = psycopg2.connect(
    "dbname='postgres' user='postgres' host='104.193.174.179' password='NRC2022!' port=5432 options='-c search_path=lipograph' ")
cur = conn.cursor()



def quick_abv_fix():
    query = "SELECT molecule_id, abbreviation from molecules where sl_lipid_class = 'SLM:000000333'"
    cur.execute(query)
    results = cur.fetchall()
    for num, row in enumerate(results):
        molecule_id = row[0]
        abv = row[1]
        cleaned_abv = clean_abbreviation(abv)
        query2 = "UPDATE molecules SET cleaned_abbreviation = %s where molecule_id = %s"
        cur.execute(query2, (cleaned_abv, molecule_id))
    conn.commit()

# quick_abv_fix()
def add_clean_abbreviation_to_reactions(table_name):
    try:
        query = f"SELECT molecule_id, abbreviation, cleaned_abbreviation FROM {table_name}"
        cur.execute(query)
        all_lipids = cur.fetchall()
        for num, row in enumerate(all_lipids):
            molecule_id = row[0]
            abv = row[1]
            if abv is not None:
                cleaned_abv = clean_abbreviation(abv)
                query = "UPDATE molecules SET cleaned_abbreviation = %s WHERE molecule_id = %s"
                cur.execute(query, (cleaned_abv, molecule_id))
                conn.commit()
                print("Updated this term successfully: " + str(cleaned_abv))
        print(f"Table {table_name} populated successfully.")

    except (Exception, psycopg2.Error) as error:
        print(f"Error while populating table {table_name}: {error}")


# Usage example
# tsv_file_path = r'C:\MastersVault\Project_One\tables_and_raw_data\SL_tables\enzymes_cleaned_fixed.tsv'
# table_name = 'molecules'
# add_clean_abbreviation_to_reactions(table_name)

def upload_organisms_from_tsv(tsv_file_path):
    try:
        # Dictionary to store unique organism IDs and taxon scientific names
        organism_dict = {}

        # Set to store unique organism IDs
        organism_ids = set()

        # Read the TSV file and extract unique organism IDs and taxon_scientific_names
        with open(tsv_file_path, 'r', newline='') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            next(reader)  # Skip the header row
            for row in reader:
                organism_id = int(row[3])  # Assuming organism_id is in the first column
                taxon_scientific_name = row[4]  # Assuming taxon_scientific_name is in the second column

                if organism_id not in organism_dict:
                    organism_dict[organism_id] = taxon_scientific_name

        # Insert the unique organism IDs and taxon scientific names into the organism table
        for organism_id, taxon_scientific_name in organism_dict.items():
            cur.execute("INSERT INTO organism (organism_id, taxon_scientific_name) VALUES (%s, %s)",
                           (organism_id, taxon_scientific_name))

        # Commit the transaction and close the connection
        conn.commit()
        # cur.close()
        conn.close()

        print("Organisms uploaded successfully.")

    except (Exception, psycopg2.Error) as error:
        print(f"Error while uploading organisms: {error}")


# tsv_file_path = r'C:\MastersVault\Project_One\tables_and_raw_data\SL_tables\enzymes_cleaned_fixed.tsv'
# upload_organisms_from_tsv(tsv_file_path)

def convert_header_names(tsv_file, new_tsv_file):
    try:
        # Read the TSV file and extract headers
        with open(tsv_file, 'r', newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            headers = next(reader)

        # Convert header names by replacing spaces with underscores
        converted_headers = [header.replace(' ', '_') for header in headers]

        # Create a new TSV file with fixed headers
        with open(new_tsv_file, 'w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(converted_headers)

        print(f"New TSV file '{new_tsv_file}' created with fixed headers.")

    except FileNotFoundError:
        print("File not found.")
    except csv.Error as e:
        print(f"Error reading or writing TSV file: {e}")


# tsv_file = r'C:\MastersVault\Project_One\tables_and_raw_data\SL_tables\lipids.tsv'
# new_tsv_file = r'C:\MastersVault\Project_One\tables_and_raw_data\SL_tables\lipids_headers.tsv'
# convert_header_names(tsv_file, new_tsv_file)

def copy_temp_sl_enzymes_from_tsv(tsv_file, table_name):
    try:
        cur.execute("SET client_encoding = 'UTF8'")
        create_query = """CREATE TABLE temp_enzymes (
                            swisslipids_id VARCHAR(13),
                            uniprotkb_acs TEXT,
                            gene_name TEXT,
                            protein_taxon TEXT,
                            taxon_scientific_name TEXT,
                            rhea_id INT,
                            reaction_text TEXT,
                            evidence_tag_id TEXT)
                        """
        cur.execute(create_query)
        # Truncate the existing data in the table (optional)
        # ... (code to truncate the table)

        # Copy the TSV file into the table
        with open(tsv_file, 'r', newline='', encoding='utf-8') as file:
            next(file)
            cur.copy_from(file, table_name, sep='\t', null='')
        conn.commit()
        # cur.close()
        print("Data copied successfully!")

    except FileNotFoundError:
        print("File not found.")
    except csv.Error as e:
        print(f"Error reading TSV file: {e}")




def populate_enzymes_from_csv(csv_file_path):
    try:
        # Step 1: Import CSV into temporary table
        temp_table_name = 'temp_enzymes'
        cur.execute(f"CREATE TEMP TABLE {temp_table_name} (enzyme_name TEXT, uniprot_id TEXT, organism_id INT)")

        with open(csv_file_path, 'r', newline='') as csv_file:
            reader = csv.reader(csv_file)
            next(reader)  # Skip the header row
            for row in reader:
                enzyme_name = row[2]  # Assuming enzyme_name is in the first column of the CSV
                uniprot_id = row[1]  # Assuming uniprot_id is in the second column of the CSV
                organism_id = int(row[2])  # Assuming organism_id is in the third column of the CSV

                cur.execute(f"INSERT INTO {temp_table_name} (enzyme_name, uniprot_id, organism_id) "
                               f"VALUES (%s, %s, %s)", (enzyme_name, uniprot_id, organism_id))

        # Step 2: Map organism IDs
        mapped_ids = {}
        cur.execute("SELECT DISTINCT organism_id FROM temp_enzymes")
        organism_ids = cur.fetchall()

        for organism_id in organism_ids:
            cur.execute("SELECT organism_id FROM organism WHERE organism_id = %s", (organism_id,))
            matched_id = cur.fetchone()

            if matched_id:
                mapped_ids[organism_id] = matched_id[0]

        # Step 3: Insert data into enzymes table
        cur.execute("INSERT INTO enzymes (enzyme_name, uniprot_id, organism_id) "
                       "SELECT enzyme_name, uniprot_id, %s FROM temp_enzymes", (mapped_ids,))

        # Commit the transaction and close the connection
        conn.commit()
        # cur.close()
        conn.close()

        print("Enzymes table populated successfully.")

    except (Exception, psycopg2.Error) as error:
        print(f"Error while populating enzymes: {error}")


# Usage example
# csv_file_path = 'path/to/enzymes.csv'
# populate_enzymes_from_csv(csv_file_path)


def copy_temp_sl_lipids_from_tsv(tsv_file, table_name):
    try:
        # Truncate the existing data in the table (optional)
        # ... (code to truncate the table)

        # Copy the TSV file into the table
        with open(tsv_file, 'r', newline='') as file:
            next(file)
            cur.copy_from(file, table_name, sep='\t', null='')
        conn.commit()
        # cur.close()
        print("Data copied successfully!")

    except FileNotFoundError:
        print("File not found.")
    except csv.Error as e:
        print(f"Error reading TSV file: {e}")



def upload_reactions_tsv_to_postgresql(tsv_file_path, table_name):
    try:
        unique_reactions = set()
        # Execute the COPY command to upload the TSV file and insert data into the table
        with open(tsv_file_path, 'r', encoding='utf-8') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            next(reader)  # Skip the header row
            for row in reader:
                reaction_text = row[6].strip()
                # unique_reaction_text
                query = f"INSERT INTO {table_name} (reaction_text, rhea_id) VALUES (%s, %s)"
                cur.execute(query, (row[6].strip(), row[5]))

        # Commit the transaction and close the connection\
        conn.commit()
        # cur.close()
        conn.close()

        print(f"Table {table_name} populated successfully.")

    except (Exception, psycopg2.Error) as error:
        print(f"Error while populating table {table_name}: {error}")




def add_to_reaction_pairs_table(reaction_enzyme_id, combinations):
    """
    Add reactant and product pairs to the reaction_pairs table.
    :param reaction_enzyme_id: The ID of the reaction enzyme associated with the pairs.
    :param reactants: List of reactants.
    :param products: List of products.
    """
    try:
        # Iterate over each reactant and product pair
        for combination in combinations:
            reactant = combination[0]
            product = combination[1]

            # Query the molecules table to get the molecule_id for reactant and product

            reactant_query = "SELECT molecule_id FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
            cur.execute(reactant_query, {'molecule_name': reactant, 'molecule_synonyms': f"%|" + reactant + "|%"})
            # reactant_molecule_id = cur.fetchone()
            reactant_molecule_id_initial = cur.fetchall()
            returned_reactant_count = cur.rowcount
            if returned_reactant_count > 1:
                # TODO this is a problem for Ceramide (d18:0) and Ceramide (d20:0) and Ceramide (t18:0)
                print(f"This molecule as a reactant: {reactant} had multiple matches for molecule_id. The reaction_enzyme_id is {reaction_enzyme_id}")
            elif returned_reactant_count == 1:
                for num, row in enumerate(reactant_molecule_id_initial):
                    reactant_molecule_id = row[0]

            product_query = "SELECT molecule_id FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
            cur.execute(product_query, {'molecule_name': product, 'molecule_synonyms': f"%|" + product + "|%"})
            # cur.execute("SELECT molecule_id FROM molecules WHERE molecule_name = %s", (product,))
            # product_molecule_id = cur.fetchone()
            product_molecule_id_initial = cur.fetchall()
            returned_product_count = cur.rowcount
            # TODO: Need to deal with
            if returned_product_count > 1:
                print(f"This molecule as a product: {product} had multiple matches for molecule_id. The reaction_enzyme_id is {reaction_enzyme_id}")
            elif returned_product_count == 1:
                for num, row in enumerate(product_molecule_id_initial):
                    product_molecule_id = row[0]

            if not reactant_molecule_id_initial:
                print(reactant)
                continue
            elif not product_molecule_id_initial:
                print(product)
                continue

            print(f"Ready to insert this reaction_enzyme_id: {reaction_enzyme_id}, the reactant_id = {reactant_molecule_id} and product_id = {product_molecule_id}")
            # Insert the reactant and product into the reaction_pairs table
            cur.execute(
                "INSERT INTO reaction_pairs (reaction_enzyme_id, reactant_molecule_id, product_molecule_id) VALUES (%s, %s, %s)",
                (reaction_enzyme_id, reactant_molecule_id, product_molecule_id))

        # Commit the changes and close the connection
        conn.commit()
        # conn.close()
    except (Exception, psycopg2.Error) as error:
        print(f"Error while populating table: Reactions: {error}")

def populate_reaction_pairs():
    """
    Populate the reaction_pairs table with data from the reaction_enzyme and reactions tables.
    """
    # Execute the SQL query to retrieve reaction_enzyme_id, reaction_id, and reaction_text
    cur.execute("SELECT re.reaction_enzyme_id, re.reaction_id, r.reaction_text FROM reaction_enzyme AS re JOIN reactions AS r ON re.reaction_id = r.reaction_id")
    # cur.execute("""SELECT re.reaction_enzyme_id, re.reaction_id, r.reaction_text
    # FROM reaction_enzyme AS re JOIN reactions AS r ON re.reaction_id = r.reaction_id
    # WHERE r.reaction_text = 'psychosine + H2O => D-galactose + octadecasphing-4-enine'""")
    # Fetch all the rows from the result
    rows = cur.fetchall()

    # Iterate over each row
    for row in rows:
        reaction_enzyme_id = row[0]
        reaction_id = row[1]
        reaction_text = row[2]

        # Split the reaction text into reactants and products
        reaction_components, reversible = split_reaction_text(reaction_text)
        if (reaction_components, reversible) == (-1, -1):
            continue
        reactants = reaction_components[0]
        products = reaction_components[1]

        # Generate the combinations of reactants and products
        combinations = generate_reaction_combinations(reactants, products, reversible)

        # Add the combinations to the reaction_pairs table
        add_to_reaction_pairs_table(reaction_enzyme_id, combinations)

    # Commit the changes and close the connection
    # conn.commit()
    # conn.close()

# populate_reaction_pairs()

def combinations_check(reactants, products, not_found_list, metabolites):

    for reactant in reactants:
        if reactant.startswith('a '):
            reactant = reactant[2:]
        elif reactant.startswith('an '):
            reactant = reactant[3:]
        reactant_query = "SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
        cur.execute(reactant_query, {'molecule_name': reactant, 'molecule_synonyms': f"%|" + reactant + "|%"})
        reactant_molecule_id = cur.fetchone()
        if reactant_molecule_id is None and matching_algo.metabolite_exclusion_list(reactant) is False:
            if reactant not in not_found_list and reactant not in metabolites:
                metabolites.add(reactant)
                not_found_list.append(reactant)
                print(reactant)

    for product in products:
        if product.startswith('a '):
            product = product[2:]
        elif product.startswith('an '):
            product = product[3:]
        product_query = "SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
        cur.execute(product_query, {'molecule_name': product, 'molecule_synonyms': f"%|" + product + "|%"})
        product_molecule_id = cur.fetchone()
        if product_molecule_id is None and matching_algo.metabolite_exclusion_list(product) is False:
            if product not in not_found_list and product not in metabolites:
                metabolites.add(product)
                not_found_list.append(product)
                print(product)



def find_missing_lipids():
    # excel_file_path = r"C:\MastersVault\Project_One\lipid_not_foundv2_working.xlsx"
    #
    # df = pd.read_excel(excel_file_path)
    metabolites = set()
    # for num, row in df.iterrows():
    #     # if num == 0:
    #     #     continue
    #     molecule_name = row['Molecule_name']
    #     annotation = row['Addressed?']
    #     abbreviation = row['Abbreviation']
    #     if annotation in ['no code', 'not found', 'metabolite']:
    #         metabolites.append(molecule_name)


    not_found_lipids = []
    cur.execute("SELECT reaction_text, rhea_id FROM reactions")

    # Fetch all the rows from the result
    rows = cur.fetchall()

    for num, row in enumerate(rows):
        reaction_text = row[0]
        rhea_id = row[1]
        # Split the reaction text into reactants and products
        reaction_components, reversible = clean_up.split_reaction_text(reaction_text)
        if (reaction_components, reversible) == (-1, -1):
            continue

        reactants = reaction_components[0]
        products = reaction_components[1]

        combinations_check(reactants, products, not_found_lipids, metabolites)

    df = pd.DataFrame(metabolites)
    df.to_csv("lipids_not_found_Oct25.tsv", sep='\t', index=False, encoding='utf-8')



def concatenate_to_synonyms(existing_synonyms, new_synonym):
    if existing_synonyms is None:
        synonym_list = []
        synonym_list.append(new_synonym)
        updated_synonyms = " | ".join(synonym_list)
        return updated_synonyms
    synonym_list = existing_synonyms.split(" | ")
    if new_synonym in synonym_list:
        print("The synonym already exists.")
        return -1
        # You may choose to return or handle the duplicate synonym in some way
    else:
        # Add the new synonym to the list
        synonym_list.append(new_synonym)

        # Join the updated synonym list back into a string
        updated_synonyms = " | ".join(synonym_list)
        return updated_synonyms

def update_molecule_synonyms(molecule_id, new_value):
    try:
        query = f"UPDATE molecules SET molecule_synonyms = %s WHERE molecule_id = %s"
        cur.execute(query, (new_value, molecule_id))
        print(f"Updated this molecule_id: {molecule_id}")
        conn.commit()

    except (Exception, psycopg2.Error) as error:
        print(f"Error while updating table: {table_name}: {error}")

def add_cleaned_name_to_synonyms():
    query = r"select molecule_id, molecule_name, molecule_synonyms from molecules where molecule_name LIKE 'Ganglioside %' and molecule_name LIKE '%(d18:1(4E))'"
    cur.execute(query)

    result = cur.fetchall()

    for num, row in enumerate(result):
        molecule_id = row[0]
        molecule_name = row[1]
        existing_synonyms = row[2]
        clean_name = clean_up.clean_abbreviation(molecule_name)
        print(f"Name: {molecule_name} is now")
        print(f"Name: {clean_name}")
        new_synonyms = concatenate_to_synonyms(existing_synonyms, clean_name)
        print(f"Old synonyms: {existing_synonyms}")
        print(f"New synonyms: {new_synonyms}")
        print("------------")
        if new_synonyms == -1:
            # synonym already exists in list, so go to next lipid
            continue
        else:
            update_molecule_synonyms(molecule_id, new_synonyms)
    print("Done updating!")


# add_cleaned_name_to_synonyms()

def add_abbreviations_from_list():
    abbreviations = [
        ['heptadecasphing-4-enine', 'Sph(d17:1)'],
        ['sphing-4-enine', 'Sph(d18:1)'],
        ['sphinganine', 'Sph(d18:0)'],
        ['sphingosine 1-phosphate', 'S1P(d18:1)'],
        ['sphing-4-enine-phosphocholine', 'SM(d18:1/0:0)'],
        ['1-hexadecanoyl-N-(acetyl)-sphing-4-enine', '1-O-16:0-Cer(d18:1/2:0)'],
        ["ß-D-galactosyl-(11')-sphing-4-enine", "beta-GalSph(d18:1)"]

        # need to add beta- to all Gal and Glc
    ]

    for abbreviation in abbreviations:
        full_name = abbreviation[0]
        abv = abbreviation[1]
        query = "UPDATE molecules SET abbreviation = %s, cleaned_abbreviation = %s WHERE molecule_name = %s"
        cur.execute(query, (abv, abv, full_name,))
    print("Updated abbreviations using list!")
    conn.commit()


# add_abbreviations_from_list()


def show_abbreviation_from_list():
    abbreviations = [
        ['heptadecasphing-4-enine', 'Sph(d17:1)'],
        ['sphing-4-enine', 'Sph(d18:1)'],
        ['sphinganine', 'Sph(d18:0)'],
        ['sphingosine 1-phosphate', 'S1P(d18:1)'],
        ['sphing-4-enine-phosphocholine', 'SM(d18:1/0:0)'],
        ['1-hexadecanoyl-N-(acetyl)-sphing-4-enine', '1-O-16:0-Cer(d18:1/2:0)'],
        ["ß-D-galactosyl-(11')-sphing-4-enine", "beta-GalSph(d18:1)"]
        # need to add beta- to all Gal and Glc
    ]
    for molecule in abbreviations:
        molecule_name = molecule[0]
        abbreviation = molecule[1]
        query = f"SELECT molecule_name, abbreviation FROM molecules where molecule_name = %s"
        cur.execute(query, (molecule_name,))
        result = cur.fetchall()
        for num, row in enumerate(result):
            mol_name = row[0]
            abv = row[1]
            print(f"For this molecule: {mol_name} the abbreviation is: {abv}")


# show_abbreviation_from_list()


def add_synonyms_from_list():
    synonyms = [
        ['sphinganine', 'octadecasphinganine'],
        ['sphingosine 1-phosphate', 'sphing-4-enine 1-phosphate'],
        ['sphing-4-enine', ['octadecasphing-4-enine', 'sphingosine']],
        ['1-hexadecanoyl-N-(acetyl)-sphing-4-enine', '1-hexadecanoyl-N-acetylsphingosine'],
        ['Ganglioside GM1 (d18:1(4E))', 'Ganglioside GM1 (d18:1)'],
        ['Ganglioside GM2 (d18:1(4E))', 'Ganglioside GM2 (d18:1(4E))']
    ]

    for pair in synonyms:
        molecule_name = pair[0]
        synonyms = pair[1]
        query = f"SELECT molecule_name, molecule_synonyms FROM molecules where molecule_name = %s"
        cur.execute(query, (molecule_name,))
        result = cur.fetchall()
        for num, row in enumerate(result):
            molecule_name = row[0]
            existing_synonyms = row[1]
            if isinstance(synonyms, list):
                for synonym in synonyms:
                    new_synonyms = concatenate_to_synonyms(existing_synonyms, synonym)
                    if new_synonyms == -1:
                        continue
                    else:
                        print(f"For this molecule: {molecule_name}, the old synonym is: {existing_synonyms}")
                        print(f"The new synonyms are: {new_synonyms}")
                        query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_name = %s"
                        cur.execute(query, (new_synonyms, molecule_name))

            else:
                new_synonyms = concatenate_to_synonyms(existing_synonyms, synonyms)
                if new_synonyms == -1:
                    continue
                else:
                    print(f"For this molecule: {molecule_name}, the old synonym is: {existing_synonyms}")
                    print(f"The new synonyms are: {new_synonyms}")
                    query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_name = %s"
                    cur.execute(query, (new_synonyms, molecule_name))
    conn.commit()

# add_synonyms_from_list()

def concatenate_synonym(existing_synonyms, new_synonym):
    if existing_synonyms:
        # Remove leading and trailing '|' if present
        existing_synonyms = existing_synonyms.strip('|')
        synonyms_list = existing_synonyms.split('|')
    else:
        # If no existing synonyms, initialize an empty list
        synonyms_list = []

    # Check if the new synonym already exists
    if new_synonym not in synonyms_list:
        # Add the new synonym to the list
        synonyms_list.append(new_synonym)

    # Join the list back into a string with '|' separator
    updated_synonyms = '|' + '|'.join(synonyms_list) + '|'

    return updated_synonyms

def bulk_add_to_synonyms():
    query = "select molecule_id, cleaned_abbreviation, molecule_synonyms from molecules where cleaned_abbreviation LIKE '% | %'"
    cur.execute(query)
    result = cur.fetchall()
    for num, row in enumerate(result):
        molecule_id = row[0]
        cleaned_abv = row[1]
        molecule_synonyms = row[2]
        updated_synonyms = molecule_synonyms

        abbreviations = cleaned_abv.split(" | ")

        if len(abbreviations) > 1:
            new_cleaned_abv = abbreviations[0]
            new_synonyms = abbreviations[1:]

            for synonym in new_synonyms:
                updated_synonyms = concatenate_synonym(updated_synonyms, synonym)

            query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_id = %s"
            cur.execute(query, (updated_synonyms, molecule_id))

            query = f"UPDATE molecules SET cleaned_abbreviation = %s where molecule_id = %s"
            cur.execute(query, (new_cleaned_abv, molecule_id))

            print(f"Committing these synonyms: {updated_synonyms} and this is the new cleaned_abv: {new_cleaned_abv}")
            conn.commit()

# bulk_add_to_synonyms()


def bulk_add_pc_to_synonyms():
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #            WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000000352'"""
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #            WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000000724'"""
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #              WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000501063'"""
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #            WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000758183'"""
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #            WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000049329'"""
    # query = """SELECT molecule_id, cleaned_abbreviation, molecule_synonyms FROM molecules
    #            WHERE cleaned_abbreviation IS NOT NULL AND sl_lipid_class = 'SLM:000399706'"""
    # query = "select molecule_id, cleaned_abbreviation, molecule_synonyms from molecules where cleaned_abbreviation LIKE 'CerP%'"
    # query = "select molecule_id, cleaned_abbreviation, molecule_synonyms from molecules where sl_lipid_class = 'SLM:000043367'"
    # query = "select molecule_id, cleaned_abbreviation, molecule_synonyms from molecules where cleaned_abbreviation LIKE 'PE-Cer%'"
    query = "select molecule_id, cleaned_abbreviation, molecule_synonyms from molecules where cleaned_abbreviation LIKE 'P%(%/0:0)'"
    cur.execute(query)
    result = cur.fetchall()
    for num, row in enumerate(result):
        molecule_id = row[0]
        new_synonym = row[1]
        current_synonyms = row[2]
        updated_synonyms = concatenate_synonym(current_synonyms, new_synonym)

        query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_id = %s"
        cur.execute(query, (updated_synonyms, molecule_id))
        print(f"Updated the synonym for {new_synonym}")

    conn.commit()

# bulk_add_pc_to_synonyms()

def add_to_molecule_synonyms(molecule_id, new_synonym):
    query = r"select molecule_id, molecule_synonyms from molecules where molecule_id = %s"
    cur.execute(query, (molecule_id,))

    result = cur.fetchall()
    for num, row in enumerate(result):
        existing_synonyms = row[1]
        updated_synonyms = concatenate_synonym(existing_synonyms, new_synonym)
        # print(f"The old synonyms were {existing_synonyms} and the new synonyms are {updated_synonyms}")
        query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_id = %s"
        cur.execute(query, (updated_synonyms, molecule_id))

    conn.commit()
    print(f"Updated the synonyms to {updated_synonyms}")

# add_to_molecule_synonyms(28380, "1-O-(13-methyltetradecyl)-2-(13-methyltetradecanoyl)-sn-glycero-3-phosphoethanolamine")
# add_to_molecule_synonyms(49787, "1-O-(1Z-13-methyltetradecenyl)-2-(13-methyltetradecanoyl)-sn-glycero-3-phosphoethanolamine")
# add_to_molecule_synonyms(498417, "beta-D-galactosyl-(1->4)-beta-D-glucosyl-(11)-ceramide")
# add_to_molecule_synonyms(498416, "alpha-D-galactosyl-(1->4)-beta-D-galactosyl-(1->4)-beta-D-glucosyl-(11)-ceramide")
# add_to_molecule_synonyms(394802, "N-hexadecanoyl-sphinganine-1-phosphoethanolamine")
# add_to_molecule_synonyms(387549, "octadecasphing-4-enine-phosphocholine")
# add_to_molecule_synonyms(387669, "beta-D-Gal-(11')-Cer")
# add_to_molecule_synonyms(561, "Ganglioside Gb3 (d18:1)")
# add_to_molecule_synonyms(387478, "Sphingoid base 1-phosphate")
# add_to_molecule_synonyms(385397, "Sphingoid base")
# add_to_molecule_synonyms(779251, "beta-D-glucosyl-(11)-sphing-4-enine")
# add_to_molecule_synonyms(779250, "beta-D-glucosyl-(1<->1')-N-acylsphing-4-enine")
# add_to_molecule_synonyms(387669, "beta-D-galactosyl-(1<->1')-N-acylsphing-4-enine")
# add_to_molecule_synonyms(582, "beta-D-galactosyl-(1<->1')-N-acylsphing-4-enine")
# add_to_molecule_synonyms(605, "beta-D-glucosyl-(1<->1')-N-acylsphing-4-enine")
# add_to_molecule_synonyms(570, "N-acylsphing-4-enine 1-phosphate")
# add_to_molecule_synonyms(385397, "sphingoid base")
# add_to_molecule_synonyms(387669, "beta-D-galactosylceramide")
# add_to_molecule_synonyms(13183, "1-octadecanoyl-2-(9Z,12Z)-octadecadienoyl-sn-glycero-3-phosphocholine")
# add_to_molecule_synonyms(779251, "beta-D-glucosyl-(11)-sphing-4-enine")
# add_to_molecule_synonyms(393408, "N-octanoylsphing-4-enine")
# add_to_molecule_synonyms(392940, "N-hexanoyl-sphing-4-enine-1-phosphocholine")
# add_to_molecule_synonyms(393181, "beta-D-glucosyl-N-(9Z-octadecenoyl)-sphing-4E-enine")
# add_to_molecule_synonyms(385397, "Sphingoid bases")
# add_to_molecule_synonyms(387549, "SM(d18:1/0:0)")
# add_to_molecule_synonyms(55, "fatty acyl-CoA")
# add_to_molecule_synonyms(394802, "N-hexadecanoyl-4-hydroxysphinganine-1-phosphoethanolamin")
# add_to_molecule_synonyms(496020, "cholesteryl (9Z,12Z,15Z-octadecatrienoate)")
# add_to_molecule_synonyms(393197, "beta-D-glucosyl-(11)-N-octadecanoylsphing-4-enine")


def add_to_molecule_synonyms_using_sl_id(swisslipids_id, new_synonym):
    query = r"select molecule_id, molecule_synonyms, molecule_name from molecules where swisslipids_id = %s"
    cur.execute(query, (swisslipids_id,))

    result = cur.fetchall()
    for num, row in enumerate(result):
        molecule_id = row[0]
        molecule_name = row[2]
        existing_synonyms = row[1]
        if molecule_name == new_synonym:
            continue
        if existing_synonyms is not None:
            if new_synonym in existing_synonyms:
                continue
        updated_synonyms = concatenate_synonym(existing_synonyms, new_synonym)
        # print(f"The old synonyms were {existing_synonyms} and the new synonyms are {updated_synonyms}")
        query = f"UPDATE molecules SET molecule_synonyms = %s where molecule_id = %s"
        cur.execute(query, (updated_synonyms, molecule_id))

        conn.commit()
        print(f"Updated the synonyms to {updated_synonyms}")

def query_print_test():
    molecule_name = 'sphingosine'
    molecule_synonym = 'sphingosine'
    query_string = "SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
    print(query_string)
    whatis = cur.mogrify(query_string, {'molecule_name': molecule_name, 'molecule_synonyms': f"%|" + molecule_name + "|%"})
    # print(whatis)




# ---------------------------------------------------------------

def check_molecule_exists(molecule_name):
    new_query = """SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or 
    cleaned_abbreviation LIKE %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"""
    cur.execute(new_query, {'molecule_name': molecule_name, 'molecule_synonyms': f"%|" + molecule_name + "|%"})
    molecule_id = cur.fetchone()
    # TODO: double check that fetchone returns a tuple or integer
    if isinstance(molecule_id, tuple):
        return molecule_id[0]
    return molecule_id


def add_molecule(molecule_name, swisslipids_id, abbreviation, cleaned_abbreviation, molecule_synonyms, sl_lipid_class):
    query = """INSERT INTO molecules (molecule_name, swisslipids_id, abbreviation, 
    cleaned_abbreviation, molecule_synonyms, sl_lipid_class) VALUES (%s, %s, %s, %s, %s, %s) RETURNING molecule_id"""
    cur.execute(query, (molecule_name, swisslipids_id, abbreviation, cleaned_abbreviation, molecule_synonyms, sl_lipid_class))
    new_molecule_id = cur.fetchone()[0]
    conn.commit()
    return new_molecule_id


def main_molecule_addition(molecule_name):
    molecule_info = {
        'molecule_name': "ß-D-galactosyl-(11')-sphingoid base",
        'swisslipids_id': None,
        'abbreviation': 'beta-GalSph',
        'cleaned_abbreviation': 'beta-GalSph',
        'molecule_synonyms': None,
        'sl_lipid_class': 'SLM:000501452'
    }
    molecule_name = molecule_info['molecule_name']
    swisslipids_id = molecule_info['swisslipids_id']
    abbreviation = molecule_info['abbreviation']
    cleaned_abbreviation = molecule_info['cleaned_abbreviation']
    molecule_synonyms = molecule_info['molecule_synonyms']
    sl_lipid_class = molecule_info['sl_lipid_class']

    molecule_id = check_molecule_exists(molecule_name)
    if not molecule_id:
        add_molecule(molecule_name, swisslipids_id, abbreviation, cleaned_abbreviation, molecule_synonyms, sl_lipid_class)
    print(f"Added molecule: {molecule_name} successfully!")

# main_molecule_addition('Beta-D-glucosylceramide')
# main_molecule_addition("alpha-D-Gal-(1->4)-beta-D-Gal-(11')-Cer")
# main_molecule_addition("cholesteryl-beta-D-glucoside")
# main_molecule_addition("cholesteryl 3-beta-D-galactoside")
# main_molecule_addition("N-acyldeoxysphinganine")
# main_molecule_addition("1-deoxysphingosine")
# main_molecule_addition("N-(acyl)-sphingosyl-1,3-cyclic phosphate")
# main_molecule_addition("N-dodecanoyl-sphing-4-enine-1,3-cyclic phosphate")

# main_molecule_addition("ß-D-galactosyl-(11')-sphingoid base")


def check_organism_exists(organism_id):
    query = "SELECT organism_id FROM organism WHERE organism_id = %s"
    cur.execute(query, (organism_id,))
    organism_id = cur.fetchone()
    # TODO: double check that fetchone returns a tuple or integer
    if isinstance(organism_id, tuple):
        return organism_id[0]
    return organism_id


def add_organism(taxon_scientific_name, organism_id):
    # answer = input(f"Do you really want to add this organism: {taxon_scientific_name} with this organism_id: {organism_id}?")
    # if answer.lower() != 'yes':
    #     print(f"Did not add this organism: {taxon_scientific_name} with this organism_id: {organism_id}")
    #     return
    query = "INSERT INTO organism (taxon_scientific_name, organism_id) VALUES (%s, %s) RETURNING organism_id"
    cur.execute(query, (taxon_scientific_name, organism_id,))
    new_organism_id = cur.fetchone()[0]
    conn.commit()
    return new_organism_id




def check_enzyme_exists(uniprot_id):
    query = "SELECT enzyme_id FROM enzymes WHERE uniprot_id = %s"
    cur.execute(query, (uniprot_id,))
    enzyme_id = cur.fetchone()
    if isinstance(enzyme_id, tuple):
        return enzyme_id[0]
    return enzyme_id

def add_enzyme(enzyme_name, uniprot_id, organism_id, swisslipids_p_id):
    # organism_id = check_organism_exists(organism_id)
    # if not organism_id:
        # organism_id = add_organism(taxon_scientific_name, organism_id)
    query = "INSERT INTO enzymes (enzyme_name, uniprot_id, organism_id, swisslipids_p_id) VALUES (%s, %s, %s, %s) RETURNING enzyme_id"
    # TODO: handle the fact that most new enzymes will not have swisslipids_id, so they should just be empty
    cur.execute(query, (enzyme_name, uniprot_id, organism_id, swisslipids_p_id))
    new_enzyme_id = cur.fetchone()[0]
    conn.commit()
    return new_enzyme_id


def check_reaction_exists(reaction_text):
    query = "SELECT reaction_id FROM reactions WHERE reaction_text = %s"
    cur.execute(query, (reaction_text,))
    reaction_id = cur.fetchone()
    if isinstance(reaction_id, tuple):
        return reaction_id[0]
    return reaction_id

def add_reaction(reaction_text, rhea_id, doi):
    query = "INSERT INTO reactions (reaction_text, rhea_id, doi) VALUES (%s, %s, %s) RETURNING reaction_id"
    cur.execute(query, (reaction_text, rhea_id, doi))
    new_reaction_id = cur.fetchone()[0]
    conn.commit()
    return new_reaction_id


def check_reaction_enzyme_exists(reaction_id, enzyme_id):
    query = "SELECT reaction_enzyme_id FROM reaction_enzyme WHERE reaction_id = %s AND enzyme_id = %s"
    cur.execute(query, (reaction_id, enzyme_id,))
    reaction_enzyme_id = cur.fetchone()
    if isinstance(reaction_enzyme_id, tuple):
        return reaction_enzyme_id[0]
    return reaction_enzyme_id

def add_reaction_enzyme(reaction_id, enzyme_id):
    query = "INSERT INTO reaction_enzyme (reaction_id, enzyme_id) VALUES (%s, %s) RETURNING reaction_enzyme_id"
    cur.execute(query, (reaction_id, enzyme_id))
    new_reaction_enzyme_id = cur.fetchone()[0]
    conn.commit()
    return new_reaction_enzyme_id


def add_reaction_enzyme_ortholog(reaction_id, enzyme_id, human_ortholog_gene):
    query = "INSERT INTO reaction_enzyme (reaction_id, enzyme_id, human_ortholog_enzyme_id) VALUES (%s, %s, %s) RETURNING reaction_enzyme_id"
    cur.execute(query, (reaction_id, enzyme_id, human_ortholog_gene))
    new_reaction_enzyme_id = cur.fetchone()[0]

    conn.commit()
    return new_reaction_enzyme_id

def check_pair_exists(reaction_enzyme_id, reactant_molecule_id, product_molecule_id):
    query = "select pair_id from reaction_pairs where reaction_enzyme_id = %s and reactant_molecule_id = %s and product_molecule_id = %s"
    cur.execute(query, (reaction_enzyme_id, reactant_molecule_id, product_molecule_id,))
    pair_id = cur.fetchone()
    if isinstance(pair_id, tuple):
        return pair_id[0]
    return pair_id


# TODO: remove reactants, products, reversible parameters and change into reaction_text
def add_reaction_pairs(reaction_enzyme_id, reaction_text):
    # def add_reaction_pairs(reaction_enzyme_id, reaction_text):
    # def add_reaction_pairs(reaction_enzyme_id, reactants, products, reversible):
    reaction_components, reversible = split_reaction_text(reaction_text)
    if (reaction_components, reversible) == (-1, -1):
        print("Invalid reaction text. It should contain ' + + ', '=>', or '<='")
        return
    reactants = reaction_components[0]
    products = reaction_components[1]

    #uncomment this
    combinations = generate_reaction_combinations(reactants, products, reversible)

    for combination in combinations:
        reactant = combination[0]
        product = combination[1]

        # Query the molecules table to get the molecule_id for reactant and product
        # TODO: figure out what at what level to add the reaction, i.e. what is the reactant and product
        # perhaps, user could specify swisslipids_id of molecules so it can be associated
        # with molecules that are parents
        reactant_query = "SELECT molecule_id FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
        cur.execute(reactant_query, {'molecule_name': reactant, 'molecule_synonyms': f"%|" + reactant + "|%"})
        reactant_molecule_id = cur.fetchone()
        if reactant_molecule_id is not None:
            reactant_molecule_id = reactant_molecule_id[0]
            if reactant_molecule_id in [613, 387864, 55]:
                continue

        product_query = "SELECT molecule_id FROM molecules WHERE molecule_name = %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
        cur.execute(product_query, {'molecule_name': product, 'molecule_synonyms': f"%|" + product + "|%"})
        product_molecule_id = cur.fetchone()
        if product_molecule_id is not None:
            product_molecule_id = product_molecule_id[0]
            if product_molecule_id in [613, 387864, 55]: # triggered if molecule is a fatty acid, fatty alcohol or fatty_acyl-CoA
                continue

        if reactant_molecule_id is None or product_molecule_id is None:
            continue

        # check if pair already exists
        pair_id = check_pair_exists(reaction_enzyme_id, reactant_molecule_id, product_molecule_id)
        if not pair_id:
            print(f"New pair to add: reaction_enzyme_id = {reaction_enzyme_id}, reactant_id = {reactant_molecule_id} and product_id = {product_molecule_id}")
            cur.execute(
                "INSERT INTO reaction_pairs (reaction_enzyme_id, reactant_molecule_id, product_molecule_id) VALUES (%s, %s, %s)",
                (reaction_enzyme_id, reactant_molecule_id, product_molecule_id))
    conn.commit()
        # TODO finish the rest of this by adding product query, look at (add_to_reaction_pairs_table) function


def add_new_reaction_to_database(reaction_info):
    try:
        # Extract reaction information
        reaction_text = reaction_info['reaction_text']
        rhea_id = reaction_info['rhea_id']
        doi = reaction_info['doi']
        enzyme_info = reaction_info['enzyme_info']


        # Check if the organism exists, if not add it and retrieve new organism_id
        organism_id = check_organism_exists(enzyme_info['organism_id'])
        if not organism_id:
            organism_id = add_organism(enzyme_info['taxon_scientific_name'], enzyme_info['organism_id'])

        # Check if the enzyme exists, if not add it and retrieve new enzyme_id
        enzyme_id = check_enzyme_exists(enzyme_info['uniprot_id'])
        if not enzyme_id:
            enzyme_id = add_enzyme(enzyme_info['enzyme_name'], enzyme_info['uniprot_id'],
                               organism_id, enzyme_info['swisslipids_p_id'])

        # Check if the enzyme exists, if not add it and retrieve new reaction_id
        reaction_id = check_reaction_exists(reaction_text)
        if not reaction_id:
            reaction_id = add_reaction(reaction_text, rhea_id, doi)

        # Establish the relationship between reaction and enzyme
        # Check if the reaction_enzyme connection exists, and if not add it and retrieve new_reaction_enzyme_id
        reaction_enzyme_id = check_reaction_enzyme_exists(reaction_id, enzyme_id)
        if not reaction_enzyme_id:
            reaction_enzyme_id = add_reaction_enzyme(reaction_id, enzyme_id)

        # Add the reaction pairs
        add_reaction_pairs(reaction_enzyme_id, reaction_text)

        # Commit the changes to the database
        # conn.commit()

        # Close the cursor and connection
        # cur.close()
        # conn.close()

        print("New reaction successfully added to the database.")
        return True

    except psycopg2.Error as e:
        print("Error adding reaction to the database:", e)
        return False


reaction_info = {
    'reaction_text': "sphing-4-enine + UDP-alpha-D-galactose => UDP + ß-D-galactosyl-(11')-sphing-4-enine + H+",
    'rhea_id': 19488,
    'doi': '10.1007/BF00964821',
    'enzyme_info': {
        'enzyme_name': 'Ugt8',
        'uniprot_id': 'Q64676',
        'taxon_scientific_name': 'Mus musculus',
        'organism_id': 10090,
        'swisslipids_p_id': None
    }
}

# add_new_reaction_to_database(reaction_info)

# ayo, reverse = clean_up.split_reaction_text("N-acyl-sphingoid base + UDP-alpha-D-galactose => UDP + a beta-D-galactosylceramide + H+")







# ---------------------------------------------------------------



def check_if_molecule_present(lipid):
    new_query = "SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or cleaned_abbreviation LIKE %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
    cur.execute(new_query, {'molecule_name': lipid, 'molecule_synonyms': f"%|" + lipid + "|%"})
    returned_row = cur.fetchall()
    for num, row in enumerate(returned_row):
        print(row)

# check_if_molecule_present("1',3'-bis-[1,2-di-(9Z,12Z-octadecadienoyl)-sn-glycero-3-phospho]-glycerol")


def go_through_lipids_not_found():
    # excel_file_path = r"C:\MastersVault\Project_One\lipid_not_foundv2_working.xlsx"

    df = pd.read_excel(excel_file_path)

    for num, row in df.iterrows():
        # if num == 0:
        #     continue
        molecule_name = row['Molecule_name']
        annotation = row['Addressed?']
        abbreviation = row['Abbreviation']
        if annotation in ['no code', 'not found', 'metabolite']:
            continue

        # Check just using current name if it can be found, if yes nothing to be done
        # if no, use the swisslipids_id and add molecule_name to synonyms
        if annotation.startswith('SLM'):
            query_to_find_lipid = "SELECT * FROM molecules WHERE molecule_name = %(molecule_name)s or cleaned_abbreviation LIKE %(molecule_name)s or molecule_synonyms LIKE %(molecule_synonyms)s"
            cur.execute(query_to_find_lipid, {'molecule_name': molecule_name, 'molecule_synonyms': f"%|" + molecule_name + "|%"})
            result = cur.fetchall()
            if result:
                continue
            else:
                split_annotation = annotation.split(' - ')
                swisslipids_id = split_annotation[0]
                add_or_no = split_annotation[1]
                if add_or_no == 'not in lipids':
                    continue
                else:
                    query_to_find_sl_id = "SELECT * from molecules where swisslipids_id = %s"
                    cur.execute(query_to_find_sl_id, (swisslipids_id,))
                    if cur.fetchone() is None:
                        continue
                    else:
                        add_to_molecule_synonyms_using_sl_id(swisslipids_id, molecule_name)
                        print(f"Added {molecule_name} to synonyms of {swisslipids_id}")


# go_through_lipids_not_found()

def fix_missing_synonyms():
    # file_path = 'lipid_not_foundv2_working.xlsx'
    df = pd.read_excel(file_path)

    for index, row in df.iterrows():
        addressed_value = str(row['Addressed?'])
        molecule_name = row['Molecule_name']

        # Skip rows with specified conditions
        if (
                addressed_value in ['metabolite', 'no code', 'not found']
                or 'not in lipids' in addressed_value
        ):
            print(f"Skipping row {index + 2}: {molecule_name} - {addressed_value}")
            continue

        if 'add to synonyms' not in addressed_value:
            continue

        # Extract the SwissLipids ID from the row
        swisslipids_id = addressed_value.split(' - ')[0]

        add_to_molecule_synonyms_using_sl_id(swisslipids_id, molecule_name)

# fix_missing_synonyms()


def repopulate_pairs():
    # query = "SELECT re.reaction_enzyme_id, r.reaction_text FROM reaction_enzyme re JOIN reactions r on re.reaction_id = r.reaction_id"
    query = "SELECT re.reaction_enzyme_id, r.reaction_text FROM reaction_enzyme re JOIN reactions r on re.reaction_id = r.reaction_id where r.reaction_id = 2747"
    # query = "SELECT re.reaction_enzyme_id, r.reaction_text FROM reaction_enzyme re JOIN reactions r on re.reaction_id = r.reaction_id where reaction_enzyme_id = 5160"
    # query = "SELECT re.reaction_enzyme_id, r.reaction_text FROM reaction_enzyme re JOIN reactions r on re.reaction_id = r.reaction_id where r.reaction_id = 2409"
    # query = "SELECT re.reaction_enzyme_id, r.reaction_text FROM reaction_enzyme re JOIN reactions r on re.reaction_id = r.reaction_id where r.reaction_text LIKE '%1-octadecanoyl-2-(9Z,12Z)-octadecadienoyl-sn-glycero-3-phosphocholine%'"
    cur.execute(query)
    all_reaction_enzyme_rows = cur.fetchall()
    for num, row in enumerate(all_reaction_enzyme_rows):
        reaction_enzyme_id = row[0]
        reaction_text = row[1]
        add_reaction_pairs(reaction_enzyme_id, reaction_text)

# repopulate_pairs()

def find_swisslipids_protein_id_if_exists(uniprot_id):
    query = "SELECT swisslipids_p_id FROM enzymes WHERE uniprot_id = %s"
    cur.execute(query, (uniprot_id,))
    swisslipids_p_id = cur.fetchone()
    if isinstance(swisslipids_p_id, tuple):
        return swisslipids_p_id[0]
    return swisslipids_p_id


def uniprot_to_reaction_info():
    # uniprot_file = r"C:\MastersVault\Project_One\tables_and_raw_data\Uniprot_files\s1p_sph_reviewed_rhea27519\uniprotkb_rhea_27519_2023_09_28.tsv"
    df = pd.read_csv(uniprot_file, sep='\t')

    for index, row in df.iterrows():
        uniprot_id = row['Entry']
        gene_name = row['Gene Names (primary)']
        organism_name = row['Organism']
        organism_id = row['Organism (ID)']

        if pd.isnull(gene_name):
            continue

        clean_organism_name = organism_name.split("(", 1)[0].strip()

        swisslipids_p_id = find_swisslipids_protein_id_if_exists(uniprot_id)
        if not swisslipids_p_id:
            swisslipids_p_id = None

        reaction_info = {
            'reaction_text': "1-phosphosphingoid base + H2O => sphingoid base + phosphate",
            'rhea_id': None,
            'enzyme_info': {
                'enzyme_name': gene_name,
                'uniprot_id': uniprot_id,
                'taxon_scientific_name': clean_organism_name,
                'organism_id': organism_id,
                'swisslipids_p_id': swisslipids_p_id
            }
        }

        add_new_reaction_to_database(reaction_info)

# uniprot_to_reaction_info()

def multiple_enzymes_same_reaction_info():
    cers_enzyme_list = [
                            ['ACER1', 'Q8TDN7', 'SLP:000000165'],
                            ['ACER2', 'Q5QJU3', 'SLP:000000164'],
                            ['ACER3', 'Q9NUN7', 'SLP:000000680'],
                            ['ASAH1', 'Q13510', 'SLP:000000162']
                        ]

    for enzyme in cers_enzyme_list:
        gene_name = enzyme[0]
        uniprot_id = enzyme[1]
        swisslipids_id = enzyme[2]


        reaction_info = {
            'reaction_text': "Ceramide + H2O => sphingoid base + a fatty acid",
            'rhea_id': None,
            'enzyme_info': {
                'enzyme_name': gene_name,
                'uniprot_id': uniprot_id,
                'taxon_scientific_name': 'Homo sapiens',
                'organism_id': 9606,
                'swisslipids_p_id': swisslipids_id
            }
        }

        add_new_reaction_to_database(reaction_info)


# multiple_enzymes_same_reaction_info()


def connection_to_postgresql():
    query = "SELECT * from reactions where reaction_text LIKE '%ceramide%'"
    cur.execute(query)
    result = cur.fetchall()
    for num, row in enumerate(result):
        print(row)

# connection_to_postgresql()


def get_download_url_and_date(table_name):
    info_url = 'https://www.swisslipids.org/api/downloadData'
    info_response = requests.get(info_url)


    for file in info_response.json():
        if file['file'] == table_name:
            download_url = file['url']
            date = file['date']
            return download_url, date

# def download_file(url):
#     local_filename = "lipids.tsv"
#     # NOTE the stream=True parameter below
#     output_path = fr'C:\MastersVault\Project_One\AutomatedUpdates\lipids.tsv.gz'
#     with requests.get(url, stream=True) as r:
#         r.raise_for_status()
#         with open(output_path, 'wb') as f:
#             for chunk in r.iter_content(chunk_size=15728640):
#                 # If you have chunk encoded response uncomment if
#                 # and set chunk_size parameter to None.
#                 #if chunk:
#                 f.write(chunk)
#     return local_filename

# download_file("https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv")

def download_table(table_name, output_zip_file_name):
    download_url, date = get_download_url_and_date(table_name)

    # print(download_url)
    # print(date)

    url = download_url
    response = requests.get(url)

    output_path = fr'C:\MastersVault\Project_One\AutomatedUpdates\{output_zip_file_name}'

    if response.status_code == 200:
        file_content = response.content

        with open(output_path, "wb") as file:
            file.write(file_content)
    else:
        print(f"Response status code was {response.status_code}. Could not download file - please check if website is online.")

    return output_path


# def unzip_file(file):
#     with ZipFile(file, 'r') as zObject:
#         zObject.extract('enzymes.tsv', path=r'C:\MastersVault\Project_One\AutomatedUpdates')
#     zObject.close()

def unzip_file(zip_path, extract_path):
    try:
        with ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        print("Extraction successful.")
    except Exception as e:
        print("Extraction failed:", e)


def extract_gz_file(gz_path, extract_path):
    try:
        with gzip.open(gz_path, 'rb') as f_in, open(extract_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print("Extraction successful.")
    except Exception as e:
        print("Extraction failed:", e)



def main_automated_table_updates():
    # Download enzymes and lipids zipped files
    # enzymes_path = download_table('enzymes.tsv', 'enzymes.tsv.gz')
    # lipids_path = download_table('lipids.tsv', 'lipids.tsv.gz')

    # enzymes_output_file_path = r'C:\MastersVault\Project_One\AutomatedUpdates\enzymes.tsv'
    # lipids_output_file_path = r'C:\MastersVault\Project_One\AutomatedUpdates\lipids.tsv'
    # Extract zipped files
    # extract_gz_file(enzymes_path, enzymes_output_file_path)
    # extract_gz_file(lipids_path, lipids_output_file_path)
    # print(enzymes_file)

    # cleaned_enzymes_output_path = r'C:\MastersVault\Project_One\AutomatedUpdates\clean_enzymes.tsv'
    # clean_up.fix_shifted_reaction_text(enzymes_output_file_path, cleaned_enzymes_output_path)
    # clean_up.fix_reaction_text(cleaned_enzymes_output_path, cleaned_enzymes_output_path)
    # columns = clean_up.fix_tsv_file(cleaned_enzymes_output_path, cleaned_enzymes_output_path)
    # print(columns)

    # copy_temp_sl_enzymes_from_tsv(cleaned_enzymes_output_path, 'temp_enzymes')

    temp_table_creation_query = """
    DROP TABLE IF EXISTS transient_table;
    CREATE TEMP TABLE transient_table AS  
    SELECT  
        enzymes.enzyme_id,  
        enzymes.swisslipids_p_id,  
        reactions.rhea_id  
    FROM  
        enzymes  
    JOIN  
        reaction_enzyme ON enzymes.enzyme_id = reaction_enzyme.enzyme_id  
    JOIN  
        reactions ON reaction_enzyme.reaction_id = reactions.reaction_id;
    """

    cur.execute(temp_table_creation_query)
    conn.commit()

    find_new_entries_query = """
    SELECT swisslipids_id, uniprotkb_acs, gene_name, protein_taxon, taxon_scientific_name, te.rhea_id, reaction_text, evidence_tag_id  
    FROM temp_enzymes te  
    LEFT JOIN transient_table ON te.swisslipids_id = transient_table.swisslipids_p_id  
        AND te.rhea_id = transient_table.rhea_id  
    WHERE transient_table.swisslipids_p_id IS NULL AND transient_table.rhea_id IS NULL;
    """

    cur.execute(find_new_entries_query)


    new_reactions = cur.fetchall()
    for num, row in enumerate(new_reactions):
        reaction_info = {
            'reaction_text': row[6],
            'rhea_id': row[5],
            'enzyme_info': {
                'enzyme_name': row[2],
                'uniprot_id': row[1],
                'taxon_scientific_name': row[4],
                'organism_id': row[3],
                'swisslipids_p_id': row[0]
            }
        }

    #     add_new_reaction_to_database(reaction_info)




    # unzipped_file = unzip_file()
    # Step1: Download enzymes.tsv
    # Step2: Download lipid.tsv
    # Step3: Unzip enzymes and lipids
    # Step4: Clean up enzymes and lipids tables
    # Step5: Copy them into temp tables
    # Step6: Left outer join to see new additions and insert that into new table (maybe give user feedback on what is being inserted)



# main_automated_table_updates()

def is_lipid(molecule_name):
    # Define a query to check if a molecule is a lipid
    query = "SELECT molecule_id, molecule_name FROM molecules WHERE molecule_name = %(molecule_name)s OR molecule_synonyms LIKE %(molecule_synonyms)s"

    # Execute the query with the given parameters
    cur.execute(query, {'molecule_name': molecule_name, 'molecule_synonyms': f"%|" + molecule_name + "|%"})

    excluded_metabolite_parts = ["CMP", "-[ACP]", "-CoA", "fatty acids", "long-chain fatty acid", "UDP", "all-trans-retinol"]

    # Fetch the result
    result = cur.fetchone()
    if result:
        molecule = result[1]
        # if molecule.endswith("-CoA") or molecule.lower() == "fatty acids" or molecule.endswith("-[ACP]"):
        #     return False

        # for metabolite in excluded_metabolite_parts:
        #     if metabolite in molecule:
        #         return False

    # Check if the result exists (molecule is a lipid)
    return result is not None

def find_super_reactions():

    query = "SELECT reaction_id, reaction_text from reactions"
    cur.execute(query)
    reactions = cur.fetchall()

    super_reactions = []
    pbar_outer = tqdm(total=len(reactions))
    for num, row in enumerate(reactions):
        reaction_id = row[0]
        reaction_text = row[1]

        reaction_components, reversible = clean_up.split_reaction_text(reaction_text)
        if (reaction_components, reversible) == (-1, -1):
            print("Invalid reaction text. It should contain ' + + ', '=>', or '<='")
            print(f"Reaction doesn't have appropriate signs: {reaction_text}")
            continue

        reactants = reaction_components[0]
        products = reaction_components[1]

        # Count the number of lipids in reactants and products
        lipid_count_reactants = sum(is_lipid(molecule_name) for molecule_name in reactants)
        lipid_count_products = sum(is_lipid(molecule_name) for molecule_name in products)

        # Check if it's a super reaction based on your criteria
        if lipid_count_reactants >= 2 and lipid_count_products >= 1:
            super_reactions.append((reaction_id, reaction_text))

        pbar_outer.update(1)

    # print(super_reactions)
    df = pd.DataFrame(super_reactions).to_csv("super_reactionSept18v2.tsv", sep='\t', header=None, index=None)
    return super_reactions

# find_super_reactions()

def separate_text_into_lists():
    text = """
    Ceramide: SLM:000399814
    Sphingomyelin: SLM:000001000
    Fatty acids: SLM:000000984
    1,2-diacyl-sn-glycero-3-phosphoethanolamine: SLM:000000239 
    Sphingoid base: SLM:000390097
    PE-Cer / Ceramide phosphoethanolamine: SLM:000399706
    1,2-diacyl-sn-glycero-3-phosphocholine: SLM:000000261
    1-acyl-sn-glycero-3-phosphocholine: SLM:000000352
    2-acyl-sn-glycero-3-phosphocholine: SLM:000000724
    1-O-alkyl-sn-glycero-3-phosphocholine: SLM:000000448
    1-O-(1Z-alkenyl)-sn-glycero-3-phosphocholine: SLM:000049329
    sn-glycero-3-phosphoethanolamine: SLM:000000372
    1,2-diacyl-sn-glycero-3-phosphoethanolamine: SLM:000000239
    1-acyl-sn-glycero-3-phosphoethanolamine: SLM:000000354
    2-acyl-sn-glycero-3-phosphoethanolamine: SLM:000043286
    1-O-alkyl-sn-glycero-3-phosphoethanolamine: SLM:000046328
    1-O-(alk-1-enyl)-glycero-3-phosphoethanolamine: SLM:000001155
    1,2-diacyl-sn-glycero-3-phospho-L-serine: SLM:000000336
    1-acyl-sn-glycero-3-phospho-L-serine: SLM:000000335
    2-acylglycero-3-phospho-L-serine: SLM:000001168
    1,2-diacyl-sn-glycero-3-phospho-1D-myo-inositol: SLM:000000324
    1-acyl-sn-glycero-3-phospho-1D-myo-inositol: SLM:000000334
    2-acyl-sn-glycero-3-phospho-1D-myo-inositol: SLM:000000552
    Fatty alcohols: SLM:000390053
    alpha-D-Gal-(1->4)-beta-D-Gal-(11')-Cer
    ß-D-glucosyl-(11')-sphing-4-enine
    Beta-D-glucosylceramide
    """
    lines = text.strip().split('\n')
    result = []
    for line in lines:
        parts = line.strip().split(': ')
        result.append(parts)

    df = pd.DataFrame(result).to_csv("parent_molecules.tsv", sep='\t', index=None, index_label=None, header=False)
    return result

# separate_text_into_lists()
def search_parent_list(parent_list, search_term):
    results = []
    for sublist in parent_list:
        for item in sublist:
            if search_term in item:
                return True
    return False

def get_molecule_info(reactant):
    query = "SELECT molecule_id, molecule_name, sl_lipid_class, cleaned_abbreviation, swisslipids_id FROM molecules WHERE molecule_name = %(molecule_name)s OR molecule_synonyms LIKE %(molecule_synonyms)s"
    cur.execute(query, {'molecule_name': reactant, 'molecule_synonyms': f"%|" + reactant + "|%"})
    molecule_info = cur.fetchone()
    return molecule_info


def get_parent_molecule(swisslipids_id):
    if " | " in swisslipids_id:
        swisslipids_id = swisslipids_id.split(" | ")[0]
    query = "SELECT molecule_id, molecule_name, swisslipids_id, sl_lipid_class FROM molecules WHERE swisslipids_id = %s"
    cur.execute(query, (swisslipids_id,))
    molecule_info = cur.fetchone()
    return molecule_info


def extract_parent_molecules(components, list_of_parent_molecules):
    parent_components = []

    for component in components:
        if is_lipid(component):
            molecule_info = get_molecule_info(component)
            class_sl_id = molecule_info[2]
            parent_info = get_parent_molecule(class_sl_id)
            parent_sl_id = parent_info[2]
            parent_name = parent_info[1]
            parent_sl_lipid_class = parent_info[3]
            grandparent_info = get_parent_molecule(parent_sl_lipid_class)
            grandparent_name = grandparent_info[1]

            if search_parent_list(list_of_parent_molecules, parent_sl_id) or search_parent_list(list_of_parent_molecules, parent_name):
                parent_components.append(parent_name)
                continue
            elif search_parent_list(list_of_parent_molecules, grandparent_name) or search_parent_list(list_of_parent_molecules, parent_sl_lipid_class):
                parent_components.append(grandparent_name)
                continue

            print("=========================")
            print(f"Component: {component}")
            print("Choose the parent molecule:\n")
            print(f"1. Parent: {parent_name}")
            print(f"2. Grandparent: {grandparent_name}")
            print(f"3. Keep component as is: {component}\n")
            choice = input("Enter your choice (1/2/3): ")

            if choice == '1':
                    parent_components.append(parent_name)
                    list_of_parent_molecules.append([parent_name])
                    # decision = input(
                    #     f"Could not find parent of this molecule: {component} - Please check and come back again")
                    # return
            elif choice == '2':
                    parent_components.append(grandparent_name)
                    list_of_parent_molecules.append([grandparent_name])
                    # decision = input(
                    #     f"Could not find grandparent of this molecule: {component} - Please check and come back again")
                    # return
            elif choice == '3':
                parent_components.append(component)
                list_of_parent_molecules.append([component])
            else:
                print("Invalid choice. Please enter 1 or 2.")
                return
        else:
            parent_components.append(component)

    df = pd.DataFrame(list_of_parent_molecules).to_csv("parent_molecules.tsv", sep='\t', index=None, index_label=None, header=False, mode='w+')
    return parent_components

# def write_to_parent_molecule_list(new_parent_list):
#     df = pd.DataFrame(new_parent_list).to_csv("parent_molecules.tsv", sep='\t', index=None, index_label=None, header=False)

def generalize_reactions(reaction_text):
    # TODO: just read available file
    list_of_parent_molecules = separate_text_into_lists()
    components, reversible = clean_up.split_reaction_text(reaction_text)
    reactants, products = components

    parent_reactants = extract_parent_molecules(reactants, list_of_parent_molecules)
    parent_products = extract_parent_molecules(products, list_of_parent_molecules)

    return parent_reactants, parent_products

def combine_components(parent_reactants, parent_products):
    reactants_str = " + ".join(parent_reactants)
    products_str = " + ".join(parent_products)
    reaction_str = f"{reactants_str} => {products_str}"
    return reaction_str

def parse_super_reactions_from_file(file_path):
    super_reactions = []
    with open(file_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                reaction_id = row[0]
                reaction_text = row[1]
                super_reactions.append((reaction_id, reaction_text))
    return super_reactions

def main_generalize():
    generalized_reaction_list = []
    super_reactions = parse_super_reactions_from_file("super_reactionSept18v2.tsv")
    for reaction in super_reactions:
        reaction_text = reaction[1]
        parent_reactants, parent_products = generalize_reactions(reaction_text)
        new_reaction_text = combine_components(parent_reactants, parent_products)
        print(new_reaction_text)
        generalized_reaction_list.append(new_reaction_text)

# main_generalize()

def generalize_v2():
    # go through all reaction_texts, for each molecule, if lipid, get headgroup
    #
    return


def get_molecule_data_using_family(molecule_family):
    query = "SELECT molecule_id, molecule_name, cleaned_abbreviation FROM molecules WHERE cleaned_abbreviation = %s"
    cur.execute(query, (molecule_family,))
    molecule_data = cur.fetchall()
    return molecule_data


def get_sn_chain_location(molecule):
    if molecule.startswith("1-"):
        return "sn1"
    elif molecule.startswith("2-"):
        return "sn2"
    else:
        return "unexpected behaviour"


def generalize_molecules(component, parent_molecules):
    bad_flag = False
    generalized_molescules = []
    number_of_lipid_molecules = 0
    for molecule in component:
        if is_lipid(molecule):
            number_of_lipid_molecules += 1
            molecule_info = get_molecule_info(molecule)
            class_sl_id = molecule_info[2]
            molecule_cleaned_abv = molecule_info[3]

            if class_sl_id is not None:
                parent_info = get_parent_molecule(class_sl_id)
                if parent_info is not None:
                    parent_sl_id = parent_info[2]
                    parent_name = parent_info[1]
                    parent_sl_lipid_class = parent_info[3]
                else:
                    parent_sl_lipid_class = None
                    parent_name = "no parent"
                if parent_sl_lipid_class is not None:
                    grandparent_info = get_parent_molecule(parent_sl_lipid_class)
                else:
                    grandparent_info = None
                if grandparent_info is not None:
                    grandparent_name = grandparent_info[1]
                else:
                    grandparent_name = "no grandparent"

                # if (molecule in parent_molecules) or (parent_name in parent_molecules) or (grandparent_name in parent_molecules):
                if molecule in parent_molecules:
                    generalized_molescules.append(molecule)
                    if molecule in ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs']:
                        number_of_lipid_molecules -= 1
                elif parent_name in parent_molecules:
                    generalized_molescules.append(parent_name)
                    if parent_name in ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs']:
                        number_of_lipid_molecules -= 1
                elif grandparent_name in parent_molecules:
                    generalized_molescules.append(grandparent_name)
                    if grandparent_name in ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs']:
                        number_of_lipid_molecules -= 1

                else:
                    # TODO: if headgroup is None, need to deal with it
                    if molecule_cleaned_abv is None:
                        bad_flag = True
                        generalized_molescules.append(f"{molecule} has no abbreviation")
                        continue
                    else:
                        molecule_family = get_family(molecule_cleaned_abv)
                        if molecule_family is None:
                            bad_flag = True
                            generalized_molescules.append(f"|{molecule} - no family|")
                            continue
                        if molecule_family.startswith("LP") and not (
                                molecule_family.endswith("(O-)") or molecule_family.endswith("(P-)")):
                            if '_sn1' not in molecule_family and '_sn2' not in molecule_family:
                                sn_chain = get_sn_chain_location(molecule)
                                molecule_family = f"{molecule_family}_{sn_chain}"

                        molecule_data = get_molecule_data_using_family(molecule_family)
                        number_of_lipids_found = cur.rowcount
                        if number_of_lipids_found > 1:
                            print(f"Found more than 1 lipid for this headgroup: {molecule_family}")
                            print(f"Please, ensure that only one headgroup is associated with that abbreviation")
                            generalized_molescules.append(f"{molecule_family} has more than one instance in DB")
                            bad_flag = True
                            # return
                        else:
                            for num, row in enumerate(molecule_data):
                                molecule_name = row[1]
                                generalized_molescules.append(molecule_name)
            else:
                generalized_molescules.append(molecule)

        else:
            # molecule is not a lipid (likely metabolite or H2O etc.), add it as is
            generalized_molescules.append(molecule)

    return generalized_molescules, bad_flag, number_of_lipid_molecules



def generalize_reaction(specific_reaction_text):
    parent_molecules = ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs', 'N-acyl-glycinates',
                        'N-acyl-1-acyl-sn-glycero-3-phosphoethanolamine', '1,2-diacyl-sn-glycero-3-phospho-N-acylethanolamine',
                        'N-acylethanolamines', 'triacyl-sn-glycerol', '2,3-diacyl-sn-glycerol', '1,2-diacyl-sn-glycerol',
                        '1-O-alkyl-2-acyl-sn-glycerol']
    components, reversible = split_reaction_text(specific_reaction_text)
    if (components, reversible) == (-1, -1):
        return None
    reactants, products = components
    
    generalized_reactants = generalize_molecules(reactants, parent_molecules)
    generalized_products = generalize_molecules(products, parent_molecules)

    separator = " + "

    reactants = separator.join(generalized_reactants)
    products = separator.join(generalized_products)
    # generalized reactants or products is empty
    if not reactants or not products:
        return None
    reaction_text = f"{reactants} => {products}"
    return reaction_text

# generalized = generalize_reaction("Ceramide (d18:1(4E)) + 1,2-diacyl-sn-glycero-3-phosphocholine => 1,2-diacyl-sn-glycerol + Sphingomyelin (d18:1(4E))")
# print(generalized)
    

def get_combination_molecule_ids(reactants, products):
    combinations = list(itertools.product(reactants, products))

    reaction_ids = []

    for combination in combinations:
        reactant_molecule_id = combination[0]
        product_molecule_id = combination[1]

        query = """select DISTINCT r.reaction_id from reaction_pairs 
                    JOIN reaction_enzyme re on reaction_pairs.reaction_enzyme_id = re.reaction_enzyme_id
                    JOIN reactions r on re.reaction_id = r.reaction_id 
                    where reactant_molecule_id = %s and product_molecule_id = %s"""
        cur.execute(query, (reactant_molecule_id, product_molecule_id,))
        returned_reaction_ids = cur.fetchall()
        if not reaction_ids:
            reaction_ids.extend(row[0] for row in returned_reaction_ids)
        else:
            temp_reaction_ids = []
            temp_reaction_ids.extend(row[0] for row in returned_reaction_ids)
            reaction_ids = list(set(reaction_ids) & set(temp_reaction_ids))
    # print(reaction_ids)
    if reaction_ids:
        return True, reaction_ids[0]
    return False, 0



# TODO: algorithm that checks what components exist in reaction, to make sure no duplicates happen, only for lipid components
def check_reaction_lipid_components_already_in_db(reactants, products):
    # components, reversible = clean_up.split_reaction_text(reaction_text)
    # reactants, products = components

    #reactants
    reactants_molecule_ids = []
    for reactant in reactants:
        if is_lipid(reactant):
            reactant_info = get_molecule_info(reactant)
            reactant_molecule_id = reactant_info[0]
            if reactant_molecule_id in [613, 387864, 55, 210]:
                continue
            reactants_molecule_ids.append(reactant_molecule_id)
    reactants_molecule_ids.sort()

    # if reactants_molecule_ids:
    # print(f"Here's the fucked up reaction text: {reaction_text}")
    # print(f"Here are the reactant_ids: {reactants_molecule_ids}")
    # query = "select DISTINCT reactant_molecule_id from reaction_pairs where reactant_molecule_id in %s"
    # cur.execute(query, (tuple(reactants_molecule_ids),))
    # reactant_result = [row[0] for row in cur.fetchall()]

    #products
    products_molecule_ids = []
    for product in products:
        if is_lipid(product):
            product_info = get_molecule_info(product)
            product_molecule_id = product_info[0]
            if product_molecule_id in [613, 387864, 55, 210]:
                continue
            products_molecule_ids.append(product_molecule_id)
    products_molecule_ids.sort()

    # # if products_molecule_ids:
    # query = "select DISTINCT product_molecule_id from reaction_pairs where product_molecule_id in %s"
    # cur.execute(query, (tuple(products_molecule_ids),))
    # product_result = [row[0] for row in cur.fetchall()]

    # if not reactants_molecule_ids or not products_molecule_ids:

    check, reaction_id = get_combination_molecule_ids(reactants_molecule_ids, products_molecule_ids)

    if check is True:
        return True, reaction_id
    else:
        return False, 0

    # reaction_text_list = set()
    # if set(reactant_result) == set(reactants_molecule_ids):
    #     # print(f"Reactant list is identical")
    #     if set(product_result) == set(products_molecule_ids):
    #         # print(f"Product list is identical")
    #         # print(f"Reaction EXISTS in database")
    #         reaction_id = get_combination_molecule_ids(reactants_molecule_ids, products_molecule_ids)
    #         return True, reaction_id
    #     else:
    #         # print(f"Product list is NOT identical")
    #         # print(f"Reaction does not exist in database")
    #         return False, 0
    # else:
    #     # print(f"Reactant list is NOT identical")
    #     # print(f"Reaction does not exist in database")
    #     return False, 0

    # TODO: check if the components are lipids, then once I have their molecule_id, collect the combination of reactants and products and check if that already exists in reaction pairs


def bequeath_enzyme_ids(specific_reaction_id, general_reaction_id, reactants, products, reversible):
    query = "SELECT * from reaction_enzyme where reaction_id = %s"
    cur.execute(query, (specific_reaction_id,))
    specific_reaction_enzyme_ids = cur.fetchall()
    for num, row in enumerate(specific_reaction_enzyme_ids):
        enzyme_id = row[2]
        reaction_enzyme_id = check_reaction_enzyme_exists(general_reaction_id, enzyme_id)
        if not reaction_enzyme_id:
            reaction_enzyme_id = add_reaction_enzyme(general_reaction_id, enzyme_id)
        add_reaction_pairs(reaction_enzyme_id, reactants, products, reversible)


def read_reactions_from_tsv():
    # general_reactions_file_path = "C:\MastersVault\Project_One\generalized_reactionsOct27v2.tsv"

    df = pd.read_csv(general_reactions_file_path, sep='\t')

    progress_bar = tqdm(total=len(df), desc="Processing Rows")

    for index, row in df.iterrows():
        reaction_id = row['0']
        reaction_text = row['1']

        # if int(index) < 342:
        #     continue

        # reaction_id = '714'
        # reaction_text = '1-O-alkyl-2-acyl-sn-glycero-3-phosphoethanolamine + H2O => 1-O-alkyl-2-acyl-sn-glycero-3-phosphoethanolamine + H+ + Fatty acids'
        reaction_id = '2745'
        reaction_text = 'Beta-D-galactosylceramide + cholesterol => Ceramide + cholesteryl 3-beta-D-galactoside'
        # reaction_id = '368'
        # # reaction_id = '1263'
        # reaction_id = '114'
        # reaction_text = 'H2O + Ganglioside GM2 => Ganglioside GA2 + N-acetylneuraminate'
        # reaction_text = '1,2-diacyl-sn-glycero-3-phosphocholine + 1,2-diacyl-sn-glycero-3-phosphoethanolamine => 2-acyl-sn-glycero-3-phosphocholine + H+ + 1,2-diacyl-sn-glycero-3-phospho-N-acylethanolamine'
        # # reaction_text = '1,2-diacyl-sn-glycero-3-phosphoethanolamine => H+ + 1,2-diacyl-sn-glycero-3-phospho-N-acylethanolamine + 2-acyl-sn-glycero-3-phosphoethanolamine'

        components, reversible = split_reaction_text(reaction_text)
        if (components, reversible) == (-1, -1):
            # could not split reaction, missing arrows, return True since this reaction_text should be ignored
            # return None
            return True
        reactants, products = components

        # check if reaction_text already in DB, if yes, then only to attach all enzyme_ids of the specific reaction to the generalized one
        check, general_reaction_id = check_reaction_lipid_components_already_in_db(reactants, products)
        if check is True:
            bequeath_enzyme_ids(reaction_id, general_reaction_id, reactants, products, reversible)
            # print(f"Bequeathed enzyme_ids to existing reaction")

        # otherwise, add the new general reaction_text, and then attach the reaction_enzyme_ids to it (from the general one)
        else:
            new_general_reaction_id = check_reaction_exists(reaction_text)
            if not new_general_reaction_id:
                new_general_reaction_id = add_reaction(reaction_text, None, None)
                print(f"Added this reaction text: {reaction_text}")
            bequeath_enzyme_ids(reaction_id, new_general_reaction_id, reactants, products, reversible)
            # print(f"Bequeathed enzyme_ids to new reaction")
            # query = "SELECT * from reaction_enzyme where reaction_id = %s"
            # cur.execute(query, (reaction_id,))
            # reaction_enzyme_ids = cur.fetchall()
            # for num, row in enumerate(reaction_enzyme_ids):
            #     enzyme_id = row[2]
            #     reaction_enzyme_id = check_reaction_enzyme_exists(reaction_id, enzyme_id)
            #     if not reaction_enzyme_id:
            #         reaction_enzyme_id = add_reaction_enzyme(reaction_id, enzyme_id)
            #     add_reaction_pairs(reaction_enzyme_id, reactants, products, reversible)
        progress_bar.update(1)

# read_reactions_from_tsv()



# check_reaction_lipid_components_already_in_db("Ceramide + 1,2-diacyl-sn-glycero-3-phosphocholine => Sphingomyelin + 1,2-diacyl-sn-glycerol")


def check_reaction_components_already_in_db(reactants, products):
    pass


def check_all_components_correctly_transferred(reaction_text):
    pass

# go through reactions, for each reaction,generalize it then check if it already exists, if not then add it to list,

def has_lipid_without_headgroup_or_no_lipid_components_found(reaction_text):
    parent_molecules = ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs', 'N-acyl-glycinates']
    components, reversible = split_reaction_text(reaction_text)
    if (components, reversible) == (-1, -1):
        # could not split reaction, missing arrows, return True since this reaction_text should be ignored
        # return None
        return True
    reactants, products = components

    lipid_reactants = []
    for reactant in reactants:
        if is_lipid(reactant):
            molecule_info = get_molecule_info(reactant)
            cleaned_abv = molecule_info[3]
            if cleaned_abv is None:
                return True
            lipid_reactants.append(reactant)

    lipid_products = []
    for product in products:
        if is_lipid(product):
            molecule_info = get_molecule_info(product)
            cleaned_abv = molecule_info[3]
            if cleaned_abv is None:
                return True
            lipid_products.append(product)
    if (len(lipid_reactants) < 1) or (len(lipid_products) < 1):
        return True
    return False



def main_generalizer():
    reaction_text_list = set()
    reaction_text_to_see = set()
    skipped_reactions = set()
    query = "SELECT reaction_id, reaction_text from reactions LIMIT 150"
    cur.execute(query)
    reactions = cur.fetchall()

    # Create a tqdm progress bar
    with tqdm(total=len(reactions), desc="Processing reactions") as pbar:
        for num, row in enumerate(reactions):
            # TODO: need to check if a reaction even has at least 1 lipid component on each side of reaction
            reaction_id = row[0]
            reaction_text = row[1]
            # if reaction_text == '1,2,3-tri-(9Z-octadecenoyl)-glycerol => 1,2,3-tri-(9Z-octadecenoyl)-glycerol':
            #     print(f"hereee!!")
            print(f"Pre-generalized: {reaction_text}")
            if has_lipid_without_headgroup_or_no_lipid_components_found(reaction_text):
                # Skip and collect reactions with lipids without headgroups
                skipped_reactions.add(reaction_text)
            else:
                general_text = generalize_reaction(reaction_text)
                if general_text is None:
                    continue
                reaction_text_to_see.add(general_text)
                if not check_reaction_lipid_components_already_in_db(general_text):
                    if general_text not in reaction_text_list:
                        reaction_text_list.add((reaction_id, general_text))
            # Update the progress bar
            pbar.update(1)

    df = pd.DataFrame(reaction_text_list)
    df.to_csv("Generalized_reactions.tsv", sep='\t')

    df_to_see = pd.DataFrame(reaction_text_to_see)
    df_to_see.to_csv("Generalized_reactions_to_See.tsv", sep='\t')

    skipped_df = pd.DataFrame(skipped_reactions)
    skipped_df.to_csv("skipped_reactions.tsv", sep='\t')

# main_generalizer()


def generalize_reactionv2(reactants, products, reaction_id):
    parent_molecules = ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs', 'N-acyl-glycinates', 'Cholesterol esters'
                        'N-acyl-1-acyl-sn-glycero-3-phosphoethanolamine',
                        '1,2-diacyl-sn-glycero-3-phospho-N-acylethanolamine',
                        'N-acylethanolamines', 'triacyl-sn-glycerol', '2,3-diacyl-sn-glycerol',
                        '1,2-diacyl-sn-glycerol',
                        '1-O-alkyl-2-acyl-sn-glycerol', 'cholesterol']
    bad_reaction_flag = False
    generalized_reactants, bad_reactants_flag, number_of_lipid_reactants = generalize_molecules(reactants, parent_molecules)
    generalized_products, bad_products_flag, number_of_lipid_products = generalize_molecules(products, parent_molecules)

    if bad_reactants_flag is True or bad_products_flag is True:
        bad_reaction_flag = True

    if number_of_lipid_reactants == 0 or number_of_lipid_products == 0:
        bad_reaction_flag = True

    return generalized_reactants, generalized_products, bad_reaction_flag


def generalize_driver():

    generalized_reactions = set()
    reactions_to_troubleshoot = set()

    query = "SELECT reaction_id, reaction_text from reactions where reaction_id = 1578"
    cur.execute(query)
    reactions = cur.fetchall()

    with tqdm(total=len(reactions), desc="Processing reactions") as pbar:
        for num, row in enumerate(reactions):
            reaction_id = row[0]
            reaction_text = row[1]

            components, reversible = split_reaction_text(reaction_text)
            if (components, reversible) == (-1, -1): # triggered if reaction_text is not complete
                continue

            reactants, products = components

            original_reaction_dimension = (len(reactants), len(products))

            generalized_reactants, generalized_products, reaction_flag = generalize_reactionv2(reactants, products, reaction_id)

            generalized_reaction_dimension = (len(generalized_reactants), len(generalized_products))
            if not generalized_reactants or not generalized_products:
                reactions_to_troubleshoot.add((reaction_id, f"{reaction_text} | reactants or products empty!"))
                continue

            if original_reaction_dimension != generalized_reaction_dimension:
                reactions_to_troubleshoot.add((reaction_id, f"{reaction_text} | mismatch between dimensions before and after!"))
                continue

            if reaction_flag is True:
                reactions_to_troubleshoot.add((reaction_id, f"{reaction_text} | either no abv, or >1 family found"))
                continue

            generalized_reactants_text = " + ".join(generalized_reactants)
            generalized_products_text = " + ".join(generalized_products)

            generalized_reaction_text = f"{generalized_reactants_text} => {generalized_products_text}"
            generalized_reactions.add((reaction_id, generalized_reaction_text))

            pbar.update(1)

    troublesome_reactions_df = pd.DataFrame(reactions_to_troubleshoot)
    generalized_reactions_df = pd.DataFrame(generalized_reactions)

    troublesome_reactions_df.to_csv("troublesome_reactionsOct27v2.tsv", sep='\t')
    generalized_reactions_df.to_csv("generalized_reactionsOct27v2.tsv", sep='\t')

# generalize_driver()


def get_family_molecule_id(family_name):
    family_dict = {
        'PC': 60,
        'PC(O-)': 287,
        'PC(P-)': 48017,
        'LPC_sn1': 127,
        'LPC_sn2': 416,
        'LPC(O-)': 193,
        'LPC(P-)': 49538,
        'PE': 52,
        'PE(O-)': 200,
        'PE(P-)': 49571,
        'LPE_sn1': 129,
        'LPE_sn2': 43174,
        'LPE(O-)': 46367,
        'LPE(P-)': 571,
        'PS': 118,
        'PS(O-)': 43316,
        'PS(P-)': 387584,
        'LPS_sn1': 117,
        'LPS_sn2': 43248,
        'LPS(O-)': 870,
        'LPS(P-)': 869,
        'GA1': 498413,
        'GA2': 498421,
        'GD1a': 387337,
        'GD1b': 387359,
        'GD1c': 498507,
        'GD2': 116647,
        'GD3': 116648,
        'GM1': 487255,
        'GM1a': 387338,
        'GM1b': 498424,
        'GM2': 116646,
        'GM3': 116649,
        'GM4': 779043,
        'GQ1b': 498470,
        'GT1a': 498469,
        'GT1b': 387360,
        'GT2': 779080,
        'GT3': 779085,
        '1-O-acyl-Cer': 779108,
        'Cer': 393862,
        'SM': 624,
        'beta-GlcCer': 779250,
        'beta-GalCer': 387669,
        'Sph': 385397,
        'S1P': 387478,
        'C1P': 397565,
        'LacCer': 498417,
        'PE_Cer': 397566,
        'SulfoGalCer': 498624,
        'beta-GalSph': 779262,
        'beta-GlcSph': 779261,
        'PA': 112,
        'PA(O-)': 288,
        'PA(P-)': 46469,
        'LPA_sn1': 115,
        'LPA_sn2': 43118,
        'LPA(O-)': 196,
        'LPA(P-)': 43597,
        'MG': 122,
        'DG': 61,
        'TG': 128,
        'MLCL': 133,
        'CL': 13,
        'PG': 53,
        'PG(O-)': 506011,
        'PG(P-)': 506012,
        'LPG_sn1': 30,
        'LPG_sn2': 46386,
        'LPG(O-)': 506013,
        'LPG(P-)': 506014,
        'FA': 613,
        'NAPE': 83,
        'Chol': 78,
        'beta-GlcChol': 779253,
        'CE': 210,
        'Fuc-GA1': 498659,
        'GalNAc-GD1a': 779097,
        'GB3': 498416,
        'SulfoLacCer': 498488,
        'iGB3': 755350,
        'NAE': 113,
        'GP-NAE': 506063,
        'NALPE': 506039
    }
    # molecule_ids = f""
    # for key in family_dict.keys():
    #     if key.startswith("P") or key.startswith("L"):
    #         continue
    #     molecule_id = family_dict[key]
    #     molecule_ids = molecule_ids + f"{molecule_id}, "
    # print(molecule_ids)
    if family_name in family_dict:
        return family_dict[family_name]
    else:
        return None

# get_family_molecule_id("GA1")

def generalize_pair(reactant_molecule_id, product_molecule_id):
    reactant_family = get_family_generalization(reactant_molecule_id)
    product_family = get_family_generalization(product_molecule_id)

    # if product_family == '24S-hydroxycholesterol':
    #     print("here")

    if reactant_family is None:
        raise ValueError(f"Could not get family for this molecule: {reactant_molecule_id}")
    elif product_family is None:
        raise ValueError(f"Could not get family for this molecule: {product_molecule_id}")

    reactant_family_molecule_id = get_family_molecule_id(reactant_family)
    product_family_molecule_id = get_family_molecule_id(product_family)


    if reactant_family_molecule_id is None:
        raise ValueError(f"Family is not in family_dict: {reactant_family}")
    elif product_family_molecule_id is None:
        raise ValueError(f"Family is not in family_dict: {product_family}")

    return reactant_family_molecule_id, product_family_molecule_id

def make_generalized_reaction_text(molecule_pairs):
    reactants_components = []
    products_components = []

    for pair in molecule_pairs:
        reactant_molecule_id = pair[0]
        product_molecule_id = pair[1]

        reactant_query = "SELECT molecule_name from molecules where molecule_id = %s"
        product_query = "SELECT molecule_name from molecules where molecule_id = %s"

        cur.execute(reactant_query, (reactant_molecule_id,))
        reactant_family_name = cur.fetchone()[0]

        cur.execute(product_query, (product_molecule_id,))
        product_family_name = cur.fetchone()[0]

        if reactant_family_name not in reactants_components:
            reactants_components.append(reactant_family_name)
        if product_family_name not in products_components:
            products_components.append(product_family_name)

    reactants_components.sort()
    products_components.sort()

    generalized_reactants_text = " + ".join(reactants_components)
    generalized_products_text = " + ".join(products_components)

    generalized_reaction_text = f"{generalized_reactants_text} => {generalized_products_text}"

    return generalized_reaction_text




def process_reaction_enzyme(reaction_enzyme_id):
    """
    Function that will process a reaction_enzyme entry, and all its pairs, to make new generalized pairs and
    new generalized reaction text
    :param reaction_enzyme_id:
    :return: pair_list, list of lists, each element is a reactant-product pair, in molecule_id form
    """
    query = ("SELECT mr.cleaned_abbreviation, mp.cleaned_abbreviation FROM "
             "reaction_pairs JOIN molecules mr on reaction_pairs.reactant_molecule_id = mr.molecule_id JOIN molecules mp on "
             "reaction_pairs.product_molecule_id = mp.molecule_id WHERE reaction_enzyme_id = %s "
             "and mr.cleaned_abbreviation is not null and mp.cleaned_abbreviation is not null")
    cur.execute(query, (reaction_enzyme_id,))
    pairs = cur.fetchall()
    pair_list = []
    for num, row in enumerate(pairs):

        reactant = row[0]
        product = row[1]

        reactant_family_molecule_id, product_family_molecule_id = generalize_pair(reactant, product)

        pair_list.append([reactant_family_molecule_id, product_family_molecule_id])

    return pair_list


def add_reaction_pairsv2(new_reaction_enzyme_id, reactant_molecule_id, product_molecule_id):

    if reactant_molecule_id is None:
        raise ValueError(f"reactant_molecule_id is None for this reaction_enzyme_id: {new_reaction_enzyme_id}")
    elif product_molecule_id is None:
        raise ValueError(f"product_molecule_id is None for this reaction_enzyme_id: {new_reaction_enzyme_id}")

    new_pair_id = check_pair_exists(new_reaction_enzyme_id, reactant_molecule_id, product_molecule_id)
    if not new_pair_id:
        print(
            f"New pair to add: reaction_enzyme_id = {new_reaction_enzyme_id}, reactant_id = {reactant_molecule_id} and product_id = {product_molecule_id}")
        new_pair_id = cur.execute(
            "INSERT INTO reaction_pairs (reaction_enzyme_id, reactant_molecule_id, product_molecule_id) VALUES (%s, %s, %s) RETURNING pair_id",
            (new_reaction_enzyme_id, reactant_molecule_id, product_molecule_id))

    conn.commit()
    return new_pair_id


def add_generalized_reaction_pairs(new_reaction_enzyme_id, pair_list):
    """
    This function will add generalized pairs to reaction pairs by going through pairs in pair_list
    :param new_reaction_enzyme_id:
    :param pair_list: list of list, each element is a reactant-product pair
    :return:
    """

    for pair in pair_list:
        reactant_molecule_id = pair[0]
        product_molecule_id = pair[1]

        new_pair_id = add_reaction_pairsv2(new_reaction_enzyme_id, reactant_molecule_id, product_molecule_id)


def bequeath_reaction_enzyme_ids(specific_reaction_id, generalized_reaction_id, pair_list):
    query = "SELECT reaction_enzyme_id, enzyme_id from reaction_enzyme WHERE reaction_id = %s"
    cur.execute(query, (specific_reaction_id,))
    specific_reaction_enzyme_ids = cur.fetchall()

    for num, row in enumerate(specific_reaction_enzyme_ids):
        specific_reaction_enzyme_id = row[0]
        enzyme_id = row[1]

        new_reaction_enzyme_id = check_reaction_enzyme_exists(generalized_reaction_id, enzyme_id)
        if not new_reaction_enzyme_id:
            new_reaction_enzyme_id = add_reaction_enzyme(generalized_reaction_id, enzyme_id)

        add_generalized_reaction_pairs(new_reaction_enzyme_id, pair_list)



def check_pair_exists_no_reaction_enzyme(reactant_molecule_id, product_molecule_id):
    query = "SELECT pair_id, reaction_enzyme_id FROM reaction_pairs WHERE reactant_molecule_id = %s and product_molecule_id = %s"
    cur.execute(query, (reactant_molecule_id, product_molecule_id))
    pairs = cur.fetchall()
    if pairs is not None:
        return


def check_if_family_pairs_exists(pair_list):
    if len(pair_list) == 1:
        reactant_molecule_id = pair_list[0][0]
        product_molecule_id = pair_list[0][1]

        check_pair_exists_no_reaction_enzyme(reactant_molecule_id, product_molecule_id)

    elif len(pair_list) == 2:
        pass


def generalize_reaction_id(specific_reaction_id, gen_reactions, existing_gen_reactions):


    query = ("SELECT reaction_enzyme_id, reactions.reaction_text from reaction_enzyme re "
             "JOIN reactions ON re.reaction_id = reactions.reaction_id WHERE re.reaction_id = %s LIMIT 1")
    cur.execute(query, (specific_reaction_id,))
    result = cur.fetchone()

    reaction_enzyme_id = result[0]
    reaction_text = result[1]

    pair_list = process_reaction_enzyme(reaction_enzyme_id)

    # check_if_family_pairs_exists(pair_list)

    generalized_reaction_text = make_generalized_reaction_text(pair_list)

    # if generalized_reaction_text not in gen_reactions:
    #     gen_reactions[generalized_reaction_text] = specific_reaction_id
    # else:
    #     old_value = gen_reactions[generalized_reaction_text]
    #     new_value = f"{old_value};{specific_reaction_id}"
    #     gen_reactions[generalized_reaction_text] = new_value



    generalized_reaction_id = check_reaction_exists(generalized_reaction_text)
    if not generalized_reaction_id:
        if generalized_reaction_text in existing_gen_reactions:
            generalized_reaction_id = existing_gen_reactions[generalized_reaction_text]
            print(f"reaction already generalized to this id: {generalized_reaction_text}")
        else:
            print(f"Adding new reaction: {generalized_reaction_text}")
            generalized_reaction_id = add_reaction(generalized_reaction_text, None, None)
    bequeath_reaction_enzyme_ids(specific_reaction_id, generalized_reaction_id, pair_list)




    # answer = ""
    # while True:
    #     print("Would you like to generalize the following reaction:")
    #     print(reaction_text)
    #     print(f"into this:")
    #     print(generalized_reaction_text)
    #     answer = input(f"enter yes or no").lower().strip()
    #     if answer not in ["yes", "no"]:
    #         continue
    #     else:
    #         break
    #
    # existing_check = ""
    # if answer == "no":
    #     while True:
    #         print("Would you like to input existing generalized reaction id?")
    #         existing_check = input("enter yes or no").lower().strip()
    #         if answer not in ["yes", "no"]:
    #             continue
    #         else:
    #             break
    #     if existing_check == "yes":
    #         already_existing_general_reaction_id = int(input("enter existing reaction id:").strip())
    #         generalized_reaction_id = already_existing_general_reaction_id
    # elif answer == "yes":
    #     generalized_reaction_id = check_reaction_exists(generalized_reaction_text)
    #     if not generalized_reaction_id:
    #         generalized_reaction_id = add_reaction(generalized_reaction_text, None, None)


    # bequeath_reaction_enzyme_ids(specific_reaction_id, generalized_reaction_id, pair_list)

# generalize_reaction_id(1038)

# generalize_reaction_id(1448)

def generalize_driverv2():
    query = ("""
            SELECT DISTINCT r.reaction_id, reaction_text, mr.cleaned_abbreviation, mr.molecule_name
            FROM reaction_pairs AS rp
                JOIN reaction_enzyme AS re ON rp.reaction_enzyme_id = re.reaction_enzyme_id
                JOIN reactions AS r on re.reaction_id = r.reaction_id
                JOIN molecules AS mr ON rp.reactant_molecule_id = mr.molecule_id
                JOIN molecules AS mp ON rp.product_molecule_id = mp.molecule_id
                JOIN enzymes AS e ON re.enzyme_id = e.enzyme_id
                JOIN lipograph.organism o on e.organism_id = o.organism_id
                    WHERE (rp.reactant_molecule_id IN (SELECT molecule_id FROM molecules WHERE cleaned_abbreviation LIKE 'Cer%' ))
    """)
             # "or rp.product_molecule_id in "
             # "(SELECT molecule_id FROM molecules WHERE cleaned_abbreviation LIKE 'PC(O-%'))")

    cur.execute(query)
    reactions = cur.fetchall()

    existing_gen_reactions = {
        '1,2-diacyl-sn-glycero-3-phosphocholine + H2O => 1-acyl-sn-glycero-3-phosphocholine + fatty acid + H+': 2739,
        # '1,2-diacyl-sn-glycero-3-phosphocholine + H2O => 2-acyl-sn-glycero-3-phosphocholine + fatty acid + H+': 2740
    }

    gen_reactions = {}
    for num, row in enumerate(reactions):
        reaction_id = row[0]
        reaction_text = row[1]
        reaction_id = 1121
        reaction_text = '1,2-diacyl-sn-glycero-3-phosphocholine + H2O => H+ + Fatty acids + acyl-sn-glycero-3-phosphocholine'

        done = []
        skip = [1429, 161, 2088, 525, 1363, 537, 343, 817, 928, 1766, 1325, 1964, 2701, 617, 792, 1769, 2245]
        if reaction_id in skip:
            continue

        if reaction_id not in done:
            generalize_reaction_id(reaction_id, gen_reactions, existing_gen_reactions)
            break

    for reaction_key in gen_reactions:
        print(f"{reaction_key}: {gen_reactions[reaction_key]}")

generalize_driverv2()



def analyze_ortholog_file(human_gene, ortholog_df, human_enzyme_id):
    # file_path = 'human_all_hcop_sixteen_column.txt'
    # data = pd.read_csv(file_path, sep='\t')

    for index, row in ortholog_df.iterrows():
        human_gene_name = row['human_symbol']
        ortholog_taxon_id = row['ortholog_species']
        ortholog_gene_name = row['ortholog_species_symbol']
        human_assert_ids = row['human_assert_ids']
        ortholog_assert_ids = row['ortholog_species_assert_ids']
        human_ensembl_gene_id = row['human_ensembl_gene']
        ortholog_ensembl_gene_id = row['ortholog_species_ensembl_gene']
        ortholog_uniprot_idd = row['uniprotkb_id']

        human_uniprot_id = get_uniprot_id(human_assert_ids)
        ortholog_uniprot_id = get_uniprot_id(ortholog_assert_ids)

        ortholog_enzyme_info = {
            'enzyme_name': ortholog_gene_name,
            'uniprot_id': ortholog_uniprot_id,
            'organism_id': ortholog_taxon_id,
            'swisslipids_p_id': None
        }
        organism_id = ortholog_enzyme_info['organism_id']
        compli_organism_id = check_organism_exists(organism_id)
        if not compli_organism_id:
            continue
        print("Here we would add a new enzyme and them bequeath")
        ortholog_enzyme_id = add_new_enzyme(ortholog_enzyme_info)
        bequeath_reaction_ids(human_enzyme_id, ortholog_enzyme_id)






def get_uniprot_id(assert_ids):
    pattern = r'UniProtKB=([A-Z0-9]+)'

    match = re.search(pattern, assert_ids)

    # Extract the UniProtKB if a match is found
    if match:
        uniprotkb = match.group(1)
        return uniprotkb
    else:
        return None

def get_orthologs(df, gene_name):
    orthologs = df[df['human_symbol'] == gene_name]
    return orthologs


def add_new_enzyme(enzyme_info):
    # Check if the enzyme exists, if not add it and retrieve new enzyme_id
    enzyme_id = check_enzyme_exists(enzyme_info['uniprot_id'])
    if not enzyme_id:
        enzyme_id = add_enzyme(enzyme_info['enzyme_name'], enzyme_info['uniprot_id'],
                               enzyme_info['organism_id'], enzyme_info['swisslipids_p_id'])
    return enzyme_id


def bequeath_reaction_pairs(parent_reaction_enzyme_id, child_reaction_enzyme_id):
    """
    Function will find all pair_ids associated with parent reaction_enzyme_id and make new entries in
    reaction_pairs with child_reaction_enzyme_id and the reactant and product molecule ids of the parent
    :param parent_reaction_enzyme_id:
    :param child_reaction_enzyme_id:
    :return:
    """
    query = "SELECT * from reaction_pairs WHERE reaction_enzyme_id = %s"
    cur.execute(query, (parent_reaction_enzyme_id,))
    reaction_pairs = cur.fetchall()  # all reaction_pairs of the parent
    for num, row in enumerate(reaction_pairs):
        reaction_enzyme_id = row[1]  # will not be used as all entries will have the same one
        reactant_molecule_id = row[2]
        product_molecule_id = row[3]
        # check if pair already exists
        pair_id = check_pair_exists(child_reaction_enzyme_id, reactant_molecule_id, product_molecule_id)
        if not pair_id:
            print(
                f"New pair to add: reaction_enzyme_id = {child_reaction_enzyme_id}, reactant_id = {reactant_molecule_id} and product_id = {product_molecule_id}")
            cur.execute(
                "INSERT INTO reaction_pairs (reaction_enzyme_id, reactant_molecule_id, product_molecule_id) VALUES (%s, %s, %s)",
                (child_reaction_enzyme_id, reactant_molecule_id, product_molecule_id))
    conn.commit()

def check_if_any_reaction_pairs_exist(reaction_enzyme_id):
    query = "SELECT pair_id from reaction_pairs where reaction_enzyme_id = %s"
    cur.execute(query, (reaction_enzyme_id,))
    pairs = cur.fetchall()
    number_of_reaction_pairs_found = cur.rowcount
    if number_of_reaction_pairs_found > 0:
        return True
    return False


def bequeath_reaction_ids(human_enzyme_id, ortholog_enzyme_id):
    """
    Function that will find all reactions associated with human enzyme_id, and then bequeath those reactions
    to the ortholog enzyme_id
    :param human_enzyme_id: human enzyme ID
    :param ortholog_enzyme_id: ortholog enzyme ID
    :return:
    """
    query = "SELECT * from reaction_enzyme where enzyme_id = %s"
    cur.execute(query, (human_enzyme_id,))
    human_reaction_enzyme_ids = cur.fetchall()
    for num, row in enumerate(human_reaction_enzyme_ids):
        parent_reaction_enzyme_id = row[0]
        reaction_id = row[1]
        # check if parent reaction_enzyme has any pairs, if not (False) then skip as there's nothing to bequeath
        if check_if_any_reaction_pairs_exist(parent_reaction_enzyme_id) is False:
            continue
        ortholog_reaction_enzyme_id = check_reaction_enzyme_exists(reaction_id, ortholog_enzyme_id)
        if not ortholog_reaction_enzyme_id:
            ortholog_reaction_enzyme_id = add_reaction_enzyme_ortholog(reaction_id, ortholog_enzyme_id, human_enzyme_id)
        bequeath_reaction_pairs(parent_reaction_enzyme_id, ortholog_reaction_enzyme_id)

def initial_check_for_reaction_pairs_exist(human_enzyme_id):
    query = "SELECT * from reaction_enzyme where enzyme_id = %s"
    cur.execute(query, (human_enzyme_id,))
    human_reaction_enzyme_ids = cur.fetchall()

    if cur.rowcount == 0:
        return False
    for num, row in enumerate(human_reaction_enzyme_ids):
        reaction_enzyme_id = row[0]

        # if at least one reaction_enzyme has pairs, then return True, otherwise False
        if check_if_any_reaction_pairs_exist(reaction_enzyme_id) is True:
            return True
    return False


def human_ortholog_driver(data):

    file_path = 'human_all_hcop_sixteen_column.txt'
    # data = pd.read_csv(file_path, sep='\t', dtype={'human_symbol': str})

    query = "SELECT enzyme_id, enzyme_name from enzyme_ortho where organism_id = 9606"
    cur.execute(query)
    human_enzymes = cur.fetchall()
    with tqdm(total=592, desc="Processing enzymes") as pbar:
        for num, row in enumerate(human_enzymes):
            enzyme_id = row[0]
            gene_name = row[1]
            if initial_check_for_reaction_pairs_exist(enzyme_id) is False:
                continue
            orthologs = get_orthologs(data, gene_name)
            analyze_ortholog_file(gene_name, orthologs, enzyme_id)
            print('done')
            pbar.update(1)

# human_ortholog_driver()



def extract_uniprotkb_id(row):
    pattern = r'UniProtKB=([A-Z0-9]+)'
    matches = re.findall(pattern, row)
    if matches:
        return matches[0]  # Return the first match found
    else:
        return None  # Return None if no match found


def get_new_extracted_uniprotkb_df(df):
    df['uniprotkb_id'] = df['ortholog_species_assert_ids'].apply(extract_uniprotkb_id)


def make_unirpotKB_text_files(df):
    pattern = r'UniProtKB=([A-Z0-9]+)'
    # for num, row in df.iterrows():
    extracted = df['ortholog_species_assert_ids'].str.extractall(pattern)[0].tolist()

    extracted_unique = list(set(extracted))


    # Split into chunks of 100,000 entries
    chunk_size = 100000
    num_chunks = math.ceil(len(extracted_unique) / chunk_size)
    chunks = [extracted_unique[i:i + chunk_size] for i in range(0, len(extracted_unique), chunk_size)]

    # Write each chunk to a separate text file
    for i, chunk in enumerate(chunks):
        with open(f'uniprotKB_orthologs_unique_{i}.txt', 'w') as file:
            file.write('\n'.join(chunk))
    print("done")

def extract_uniprotkb_from_tsv(file_paths):
    uniprotkb_list = []

    df = pd.read_csv(file_paths[0], sep='\t')
    for i in range(1, 3):
        # Read the TSV file
        this_df = pd.read_csv(file_paths[i], sep='\t')
        df = pd.concat([df, this_df])

    # Extract UniProtKB IDs and add them to the list
    uniprotkb_list.extend(df['Entry'].tolist())
    return uniprotkb_list, df



# Function to check if a row is valid
def is_valid_row(row, dict_df):
    return row['uniprotkb_id'] in dict_df and row['ortholog_species'] == dict_df[row['uniprotkb_id']]

def human_ortholog_clean_up():
    file_path = 'human_all_hcop_sixteen_column.txt'
    data = pd.read_csv(file_path, sep='\t', dtype={'human_symbol': str, 'ortholog_species_symbol': str})

    entries_with_gene_name = data[(data['human_symbol'] != '-') & (data['ortholog_species_symbol'] != '-') &
                                  (data['ortholog_species_assert_ids'].str.contains('UniProtKB'))]

    df = entries_with_gene_name.copy()

    get_new_extracted_uniprotkb_df(df)


    # file_paths = ['ortholog_mapping/unique_reviewed/uniprotKB_orthologs_unique_0.tsv', 'ortholog_mapping/unique_reviewed/uniprotKB_orthologs_unique_1.tsv',
    #               'ortholog_mapping/unique_reviewed/uniprotKB_orthologs_unique_2.tsv']
    valid_uniprotkb_ids, swiss_prot_reviewed_df = extract_uniprotkb_from_tsv(file_paths)

    duplicated_uniprotkb_ids = swiss_prot_reviewed_df[swiss_prot_reviewed_df.duplicated(subset='Entry')]

    dict_df = swiss_prot_reviewed_df.set_index('Entry')['Organism (ID)'].to_dict()


    # Filter DataFrame to keep only rows with valid UniProtKB IDs
    df_cleaned = df[df['uniprotkb_id'].isin(valid_uniprotkb_ids)].copy()

    # need to check that every leftover uniprotKB id, has the correct organism id attached to it
    filtered_df = df_cleaned[
        df_cleaned.apply(lambda x: x['uniprotkb_id'] in dict_df and x['ortholog_species'] == dict_df[x['uniprotkb_id']], axis=1)]

    # filtered_df_function = df_cleaned[df_cleaned.apply(is_valid_row, axis=1, dict_df = dict_df)]
    # duplicate_rows = df[df.duplicated(subset='uniprotkb_id')]

    # duplicated_clean = filtered_df[filtered_df.duplicated(subset='uniprotkb_id')]
    human_ortholog_driver(filtered_df)

    # make_unirpotKB_text_files(entries_with_gene_name)
    print(1)

# human_ortholog_clean_up()




def clean_up_hgnc_ortholog_db():
    # file_path = 'human_all_hcop_sixteen_column.txt'
    data = pd.read_csv(file_path, sep='\t', dtype={'human_symbol': str})



# parents = ['Fatty acids', 'Fatty alcohols', 'Fatty acyl-CoAs', 'N-acyl-glycinates']
def populate_cleaned_abbreviation(tsv_file_path, table_name):
    pass

def update_pair_table_abbreviations(old_name, new_name):
    pass


def update_t_lipid_table():
    pass


def add_pairs_to_pair_table():
    pass


def add_reaction_to_reactions_table():
    pass


def add_enzymes_to_enzymes_table():
    pass


def add_reaction_to_db():
    pass


def get_reactions_for_enzyme():
    pass


def get_reactions_for_lipid():
    pass



