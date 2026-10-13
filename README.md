Trying to correlate the ATP/ADP ratio and NCT kinetics
For new folder match in folder_pairs.csv

## Analysis scripts

Scripts for the NCT kinetics, cytosol masks and ATP correlation live in [`analysis/`](analysis/).
Only code and READMEs are committed there. This repo is public, so images, TIFFs, CSVs and other
data stay out (`analysis/.gitignore` blocks them).

## Making changes

Each change goes on its own branch and is opened as a pull request against `main`; the repo owner reviews
and merges.

1. `git checkout main && git pull`
2. `git checkout -b <short-description>`
3. Edit the script, then `git add` only the `.py` / `.md` files you changed
4. Commit with a message saying what changed and why, push the branch, and open a PR
