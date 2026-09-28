import zipfile
import os
import subprocess
from datetime import datetime


def get_git_files():
    """Get all files tracked by git in the current repository."""
    try:
        result = subprocess.run(
            ["git", "ls-files"], capture_output=True, text=True, check=True
        )
        return result.stdout.strip().split("\n") if result.stdout.strip() else []
    except subprocess.CalledProcessError:
        print("Error: Not in a git repository or git command failed")
        return []
    except FileNotFoundError:
        print("Error: git command not found")
        return []


def zip_git_files(output_zip, exclude_patterns=None):
    """
    Zip all files tracked by git.

    Args:
        output_zip (str): Path to output zip file
        exclude_patterns (list): List of patterns to exclude (optional)
    """
    if exclude_patterns is None:
        exclude_patterns = []

    git_files = get_git_files()

    if not git_files:
        print("No git files found.")
        return

    # Filter out excluded patterns
    filtered_files = []
    for file in git_files:
        if not any(pattern in file for pattern in exclude_patterns):
            filtered_files.append(file)
        else:
            print(f"Excluding: {file}")

    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for file in filtered_files:
            if os.path.isfile(file):
                zipf.write(file)
                print(f"Added: {file}")
            else:
                print(f"Warning: {file} is tracked by git but not found in filesystem")

    print(f"\nSuccessfully zipped {len(filtered_files)} files into {output_zip}")


if __name__ == "__main__":
    # Get current date in YYYY-MM-DD format
    current_date = datetime.now().strftime("%Y-%m-%d")

    # Name of the output zip file
    output_zip_file = f"repo_files_{current_date}.zip"

    # Optional: exclude certain files/patterns (git already handles most via .gitignore)
    exclude_patterns = [
        # Add any additional patterns you want to exclude beyond .gitignore
        # Most filtering is already handled by your comprehensive .gitignore
    ]

    # Zip all git-tracked files
    zip_git_files(output_zip_file, exclude_patterns)
