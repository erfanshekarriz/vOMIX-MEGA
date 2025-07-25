#!/bin/bash

REPO="holab-hku/vOMIX-MEGA" # Replace with your GitHub repository
PWD=$(pwd)
DIR=$1
DIR_CLEAN="${PWD}/${DIR%/}"

show_help() {
    echo "Usage: vomix_update.sh [vOMIX directory] [options]"
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo ""
    echo "Description:"
    echo "  Script used to update a current vOMIX-MEGA development version."
    echo "  e.g."
    echo "  vomix_update.sh  /home/vOMIX-MEGA"
}

# Check for help option
if [[ "$1" == "-h" || "$1" == "--help" || "$2" == "-h" || "$2" == "--help" ]]; then
    show_help
    exit 0
fi


rm -rf vomix-tmp

echo "Cloning repository $REPO..."
git clone "https://github.com/$REPO.git" vomix-tmp

if [ $? -ne 0 ]; then
    echo "Error cloning repository $REPO."
    exit 1
fi

# 1) Snakefile
echo "[1/5] Updating Snakefile..."
oldfile="${DIR_CLEAN}/workflow/Snakefile"
newfile="vomix-tmp/workflow/Snakefile"
cmp -s "$oldfile" "$newfile"
case $? in
    0) 
        # Do nothing if cmp is successful
        ;;
    1 | 2) 
        echo "Updating Snakefile"
        mv "$newfile" "$oldfile"
        ;;
esac


# 2) config.yml
echo "[2/5] Updating config.yml..."
oldfile="${DIR_CLEAN}/config/config.yml"
newfile="vomix-tmp/config/config.yml"
if [[ -e "$oldfile" ]]; then 
    cp "$oldfile" "${DIR_CLEAN}/config/config_old.yml"
else 
    : 
fi
cmp -s "$oldfile" "$newfile"
case $? in
    0) 
        # Do nothing if cmp is successful
        ;;
    1 | 2) 
        echo "Updating Snakefile"
        mv "$newfile" "$oldfile"
        ;;
esac

# 3) Rules
echo "[3/5] Updating rules..." 
for newfile in $(find vomix-tmp/workflow/rules/ -maxdepth 1 -type f); do
    BASE=$(basename "$newfile")
    oldfile="${DIR_CLEAN}/workflow/rules/${BASE}"
    
    cmp -s "$oldfile" "$newfile"
    case $? in
        0) 
            # Do nothing if cmp is successful
            ;;
        1 | 2) 
            echo "Updating $BASE"
            mv "$newfile" "$oldfile"
            ;;
    esac
done

# 4) Environments
echo "[4/5] Updating envs..." 
for newfile in $(find vomix-tmp/workflow/envs/ -maxdepth 1 -type f); do
    BASE=$(basename "$newfile")
    oldfile="${DIR_CLEAN}/workflow/envs/${BASE}"

    cmp -s "$oldfile" "$newfile"
    case $? in
        0)
            # Do nothing if cmp is successful
            ;;
        1 | 2)
            echo "Updating $BASE"
            mv "$newfile" "$oldfile"
            ;;
    esac
done

# 5) Scripts
echo "[5/5] Updating scripts..." 
for newfile in $(find vomix-tmp/workflow/scripts/ -maxdepth 1 -type f); do
    BASE=$(basename "$newfile")
    oldfile="${DIR_CLEAN}/workflow/scripts/${BASE}"

    cmp -s "$oldfile" "$newfile"
    case $? in
        0)
            # Do nothing if cmp is successful
            ;;
        1 | 2)
            echo "Updating $BASE"
            mv "$newfile" "$oldfile"
            ;;
    esac
done

echo "Done!" 

rm -rf vomix-tmp
