#!/bin/bash

GREEN='\033[0;32m' # Zelená
RED='\033[0;31m' # Červená
NC='\033[0m' # Resetování barev

# Nastavte seznam souborů, ve kterých chcete vyhledávat
files=("StRoot/macros/RunPicoD0AnaMaker.C" "StRoot/macros/RunPicoD0EventMaker.C" "StRoot/macros/RunPicoTowerTest.C")

# Projdeme všechny soubory v seznamu
for file in "${files[@]}"
do
    # Kontrola, zda soubor existuje
    if [ -f "$file" ]; then
        # Vyhledání sekvencí v souboru
        if grep -q 'int pYear = 2014' "$file"; then
            echo -e "The production is set to ${GREEN}2014${NC} in file ${GREEN}$file${NC}."
            read -p "Do you want to change it to 2016? (y/n): " answer
            if [ "$answer" == "y" ]; then
                sed -i 's/int pYear = 2014/int pYear = 2016/g' "$file"
                echo "Succesfully set to ${GREEN}2016${NC}."
            fi
        elif grep -q 'int pYear = 2016' "$file"; then
            echo -e "The production is set to ${GREEN}2016${NC} in file ${GREEN}$file${NC}."
            read -p "Do you want to change it to 2014? (y/n): " answer
            if [ "$answer" == "y" ]; then
                sed -i 's/int pYear = 2016/int pYear = 2014/g' "$file"
                echo -e "Succesfully set to ${GREEN}2014${NC}."
            fi
        else
            echo -e "The production is set to ${RED}unknown year${NC} in $file."
        fi
    else
        echo -e "${RED}File $file does not exist.${NC}"
    fi
done
