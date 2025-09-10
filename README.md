# Authors: If the script has annotation, please ensure it's IN FULL ENGLISH!!!

## How to ignore the files > 5M

For MG-HFD project authors: Before running `git add .`, please make sure you have ran the following:

### Step 1: Create a Script to Ignore Large Files

Create a script file (e.g., ignore_large_files.sh):

```bash
#!/bin/bash

# Find all files larger than 5MB and add them to .gitignore
find . -type f -size +5M | sed 's|^\./||' >> .gitignore

# Remove duplicate entries from .gitignore
sort -u -o .gitignore .gitignore

echo "All files larger than 5MB have been added to .gitignore."

```

### Step 2: Make the Script Executable

Make sure the script is executable:

```bash
chmod +x ignore_large_files.sh
```

### Step 3: Run the Script

Run the script to update your .`gitignore`:

```bash
./ignore_large_files.sh

```
