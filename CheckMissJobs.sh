#!/bin/bash

# Vložit na začátek skriptu, jednou:
initial_lines_printed=0

shopt -s nullglob
files=(sched*.session.xml)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    echo "❌ Nebyl nalezen žádný soubor typu sched*.session.xml."
    exit 1
fi

# Mapování: zobrazený text → skutečný název souboru
declare -A display_to_file
display_options=()

for file in "${files[@]}"; do
    # Zjisti čas poslední změny ve formátu: YYYY-MM-DD HH:MM
    mod_time=$(stat -c "%y" "$file" | cut -d'.' -f1)

    # Najdi nejvyšší jobIndex v daném souboru
    jobIndex=$(awk '
        $0 ~ /<void property="jobIndex">/ {
            getline
            if ($0 ~ /<int>[0-9]+<\/int>/) {
                match($0, /<int>([0-9]+)<\/int>/, arr)
                if (arr[1] > max) max = arr[1]
            }
        }
        END { print max }
    ' "$file")

    if [ -n "$jobIndex" ]; then
        display="[$mod_time] ($((jobIndex)) jobs)"
    else
        display="[$mod_time] (neznámý počet jobů)"
    fi

    display_to_file["$display"]="$file"
    display_options+=("$display")
done

selected=0

draw_menu() {
    clear
    echo "🕓 Vyber soubor podle času vytvoření (šipky + Enter):"
    for i in "${!display_options[@]}"; do
        if [ "$i" -eq "$selected" ]; then
            echo -e "> \e[1m${display_options[$i]}\e[0m"
        else
            echo "  ${display_options[$i]}"
        fi
    done
}

while true; do
    draw_menu
    IFS= read -rsn1 key
    if [[ $key == $'\x1b' ]]; then
        read -rsn2 -t 0.001 rest
        key+=$rest
    fi

    case "$key" in
        $'\x1b[A')
            ((selected--))
            ((selected < 0)) && selected=$((${#display_options[@]} - 1))
            ;;
        $'\x1b[B')
            ((selected++))
            ((selected >= ${#display_options[@]})) && selected=0
            ;;
        "")
            break
            ;;
    esac
done

# Získání skutečného názvu po výběru
SELECTED_FILE="${display_to_file[${display_options[$selected]}]}"

# Kontrola
if [ -z "$SELECTED_FILE" ]; then
    echo "❌ Není vybrán žádný soubor!"
    exit 1
fi

# Opět spočítej počet jobů pro vybraný soubor (pokud s tím dál pracuješ)
jobIndex=$(awk '
    $0 ~ /<void property="jobIndex">/ {
        getline
        if ($0 ~ /<int>[0-9]+<\/int>/) {
            match($0, /<int>([0-9]+)<\/int>/, arr)
            if (arr[1] > max) max = arr[1]
        }
    }
    END { print max }
' "$SELECTED_FILE")

if [ -z "$jobIndex" ]; then
    echo "⚠️ Nepodařilo se najít žádný jobIndex."
else
    echo "✅ Vybral jsi: $SELECTED_FILE"
    echo "🔢 Počet jobů: $((jobIndex))"
fi
# 🧩 Možnosti k výběru (bez fajfek uvnitř)
options=("Missing root files" "Broken root files" "Missing picoDst" "All")
selected_flags=(false false false false)
checkbox_index=0

draw_checkboxes() {
    # Pokud už něco bylo vypsáno dřív, vrať kurzor na začátek výpisu
    if [ "$initial_lines_printed" -gt 0 ]; then
        tput cuu $initial_lines_printed
    fi

    # Teď vykresli menu
    echo ""
    echo "🧠 Vyber jednu nebo více možností (mezerník = zaškrtnout, Enter = potvrdit):"
    printed=2
    for i in "${!options[@]}"; do
        if [ "${selected_flags[$i]}" = true ]; then
            marker="☑"
        else
            marker="☐"
        fi
        if [ "$i" -eq "$checkbox_index" ]; then
            echo -e "> $marker ${options[$i]}"
        else
            echo "  $marker ${options[$i]}"
        fi
        ((printed++))
    done

    initial_lines_printed=$printed
}


while true; do
    draw_checkboxes
    IFS= read -rsn1 key
    if [[ $key == $'\x1b' ]]; then
        read -rsn2 -t 0.001 rest
        key+=$rest
    fi

    case "$key" in
        $'\x1b[A')
            ((checkbox_index--))
            ((checkbox_index < 0)) && checkbox_index=$((${#options[@]} - 1))
            ;;
        $'\x1b[B')
            ((checkbox_index++))
            ((checkbox_index >= ${#options[@]})) && checkbox_index=0
            ;;
        " ")
            if [ "$checkbox_index" -eq $((${#options[@]} - 1)) ]; then
                # "Vše" byla zaškrtnuta/odškrtnuta
                new_state=$(
                    if [ "${selected_flags[$checkbox_index]}" = true ]; then
                        echo false
                    else
                        echo true
                    fi
                )
                for i in "${!selected_flags[@]}"; do
                    selected_flags[$i]=$new_state
                done
            else
                # Přepni jen danou možnost
                selected_flags[$checkbox_index]=$(
                    if [ "${selected_flags[$checkbox_index]}" = true ]; then
                        echo false
                    else
                        echo true
                    fi
                )


               # Pokud se něco mění mimo "Vše", zruš zaškrtnutí "Vše"
vse_index=$((${#options[@]} - 1))
selected_flags[$vse_index]=false
            fi
            ;;
        "")
            # Zkontroluj, jestli aspoň něco je zaškrtnuté (kromě "Vše")
            any_selected=false
            for i in "${!selected_flags[@]}"; do
                if [ "$i" -lt $((${#options[@]} - 1)) ] && [ "${selected_flags[$i]}" = true ]; then
                    any_selected=true
                    break
                fi
            done

            if [ "$any_selected" = true ]; then
                break
            else
                echo -e "\n⚠️ Musíš vybrat alespoň jednu možnost (kromě 'Vše')."
                read -n1 -s -r -p "Stiskni libovolnou klávesu pro pokračování..."
            fi
            ;;
    esac
done

# 💾 Výsledný výběr
selected_options=()
for i in "${!options[@]}"; do
    # Vynech "Vše" ze seznamu
    if [ "$i" -lt $((${#options[@]} - 1)) ] && [ "${selected_flags[$i]}" = true ]; then
        selected_options+=("${options[$i]}")
    fi
done

# 📤 Výpis
echo
echo "✅ Vybral jsi:"
for opt in "${selected_options[@]}"; do
    echo " - $opt"
done

echo -e "\n🔧 Načítám informace ze souboru..."

# 1️⃣ Najít nejvyšší jobIndex
jobIndex=$(awk '
    $0 ~ /<void property="jobIndex">/ {
        getline
        if ($0 ~ /<int>[0-9]+<\/int>/) {
            match($0, /<int>([0-9]+)<\/int>/, arr)
            if (arr[1] > max) max = arr[1]
        }
    }
    END { print max }
' "$SELECTED_FILE")

# 2️⃣ Najít výstupní souborový vzor s ${JOBINDEX}.root (poslední výskyt)
patterns=$(grep -oE '[^[:space:]]*${JOBINDEX}\.root' "$SELECTED_FILE")

if [ -z "$patterns" ]; then
    patterns=$(grep '&quot;[^&]*${JOBINDEX}\.root' "$SELECTED_FILE" \
        | sed -n 's/.*&quot;\([^&]*${JOBINDEX}\.root\).*/\1/p')
fi

output_pattern=$(echo "$patterns" | tail -n1)

# 💬 Shrnutí
echo "📁 Soubory očekávám ve složce: production/"
echo "📄 Vzor názvu souboru: $output_pattern"
echo "🔢 Job indexy: 1 až $jobIndex"

# 3️⃣ Kontrola chybějících souborů
if printf "%s\n" "${selected_options[@]}" | grep -q "Missing root files"; then
    echo -e "\n🔍 Kontroluji chybějící ROOT soubory..."

    missing_files=()

    for ((i = 1; i <= jobIndex; i++)); do
        output_file="${output_pattern//\$\{JOBINDEX\}/$i}"
        full_path="production/$output_file"

        if [ ! -f "$full_path" ]; then
            missing_files+=("$output_file")
        fi
    done

    if [ ${#missing_files[@]} -eq 0 ]; then
        echo "✅ Všechny ROOT soubory existují."
    else
        echo -e "\n❌ Chybí následující ROOT soubory (${#missing_files[@]}):"
        for f in "${missing_files[@]}"; do
            echo " - production/$f"
        done
    fi
fi

# 4️⃣ Kontrola rozbitých souborů (<1 KB)
if printf "%s\n" "${selected_options[@]}" | grep -q "Broken root files"; then
    echo -e "\n🔍 Kontroluji velikost ROOT souborů..."

    broken_files=()

    for ((i = 1; i <= jobIndex; i++)); do
        output_file="${output_pattern//\$\{JOBINDEX\}/$i}"
        full_path="production/$output_file"

        if [ -f "$full_path" ]; then
            size=$(stat -c%s "$full_path" 2>/dev/null)
            if [ "$size" -lt 1024 ]; then
                broken_files+=("$output_file ($size B)")
            fi
        fi
    done

    if [ ${#broken_files[@]} -eq 0 ]; then
        echo "✅ Žádné podezřele malé ROOT soubory nenalezeny."
    else
        echo -e "\n❌ Podezřelé ROOT soubory (<1KB):"
        for f in "${broken_files[@]}"; do
            echo " - production/$f"
        done
    fi
fi
if printf "%s\n" "${selected_options[@]}" | grep -q "Missing picoDst"; then
    echo -e "\n🔍 Kontroluji logy pro výskyt chyb s přístupem k picoDst souborům..."

    log_pattern=$(grep -oE '[^[:space:]]+\$\{JOBINDEX\}\.log' "$SELECTED_FILE" \
        | grep -v '[&]' \
        | head -n1)

    if [ -z "$log_pattern" ]; then
        echo "⚠️ Nepodařilo se najít žádný logový soubor s \${JOBINDEX}.log."
    else
        echo "📄 Nalezený vzor: $log_pattern"

        missing_logs=()
        missing_picoDst_refs=()

        for ((i = 1; i <= jobIndex; i++)); do
            log_file="${log_pattern//\$\{JOBINDEX\}/$i}"
            full_path="log/$log_file"

            if [ ! -f "$full_path" ]; then
                missing_logs+=("$log_file")
            else
                matches=$(grep 'Error in <TXNetFile::CreateXClient>' "$full_path" | grep -o 'root://[^ ]*picoDst.root')
                if [ -n "$matches" ]; then
                    while IFS= read -r line; do
                        missing_picoDst_refs+=("$line")
                    done <<< "$matches"
                fi
            fi
        done

        if [ ${#missing_logs[@]} -eq 0 ]; then
            echo "✅ Všechny logy existují."
        else
            echo -e "\n❌ Následující logy chybí (${#missing_logs[@]}):"
            for f in "${missing_logs[@]}"; do
                echo " - log/$f"
            done
        fi

        if [ ${#missing_picoDst_refs[@]} -eq 0 ]; then
            echo "✅ V žádném logu nebyla nalezena chyba přístupu k picoDst souborům."
        else
            # 📁 Uložit do souboru pojmenovaného podle sched* souboru
            base_sched_name=$(basename "$SELECTED_FILE" .session.xml)
            output_file="MissingDst_${base_sched_name}.list"

            printf "%s\n" "${missing_picoDst_refs[@]}" > "$output_file"

            echo "❌ Nalezeno ${#missing_picoDst_refs[@]} chybových přístupů k picoDst souborům."
            echo "📄 Uloženo do souboru: $output_file"
        fi
    fi
fi
resubmit_options=("Ano" "Ne")
resubmit_selected=0
resubmit_lines_printed=0

draw_resubmit_menu() {
    # Vrať kurzor nahoru, pokud už bylo menu jednou vykresleno
    if [ "$resubmit_lines_printed" -gt 0 ]; then
        tput cuu "$resubmit_lines_printed"
    fi
    echo ""
    echo "🔁 Chceš resubmitnout nefunkční joby? (šipky + Enter)"
    printed=2
    for i in "${!resubmit_options[@]}"; do
        if [ "$i" -eq "$resubmit_selected" ]; then
            echo -e "> \e[1m${resubmit_options[$i]}\e[0m"
        else
            echo "  ${resubmit_options[$i]}"
        fi
        ((printed++))
    done
    resubmit_lines_printed=$printed
}

while true; do
    draw_resubmit_menu
    IFS= read -rsn1 key
    if [[ $key == $'\x1b' ]]; then
        read -rsn2 -t 0.001 rest
        key+=$rest
    fi

    case "$key" in
        $'\x1b[A'|$'\x1b[D') # šipka nahoru nebo vlevo
            ((resubmit_selected--))
            ((resubmit_selected < 0)) && resubmit_selected=$((${#resubmit_options[@]} - 1))
            ;;
        $'\x1b[B'|$'\x1b[C') # šipka dolů nebo vpravo
            ((resubmit_selected++))
            ((resubmit_selected >= ${#resubmit_options[@]})) && resubmit_selected=0
            ;;
        "")
            break
            ;;
    esac
done

# Vytiskni finální volbu a pokračuj
echo
if [ "${resubmit_options[$resubmit_selected]}" = "Ne" ]; then
    echo "👋 Skript končí. Žádné resubmission nebylo provedeno."
    exit 0
else
    echo "➡️  Pokračujeme s přípravou resubmission fáze..."
fi

if printf "%s\n" "${selected_options[@]}" | grep -q "Missing root files"; then
    echo -e "\n🚀 Spouštím resubmission pro chybějící ROOT soubory..."

    for f in "${missing_files[@]}"; do
        # Najdi číslo mezi posledním podtržítkem a .root (např. _42.root)
        if [[ "$f" =~ _([0-9]+)\.root$ ]]; then
            job_number="${BASH_REMATCH[1]}"
            echo "🔁 Resubmitting job for missing file: $job_number"
            star-submit -r "$job_number" "$SELECTED_FILE" > /dev/null 2>&1
        fi
    done
fi
if printf "%s\n" "${selected_options[@]}" | grep -q "Broken root files"; then
    echo -e "\n🚀 Spouštím resubmission pro ROOT soubory s podezřelou velikostí..."

    for f in "${broken_files[@]}"; do
        # Odstraň velikost v závorce: "D0_Output_..._42.root (271 B)"
        clean_name=$(echo "$f" | cut -d' ' -f1)

        if [[ "$clean_name" =~ _([0-9]+)\.root$ ]]; then
            job_number="${BASH_REMATCH[1]}"
            echo "🔁 Resubmitting job for broken file: $job_number"
            star-submit -r "$job_number" "$SELECTED_FILE" > /dev/null 2>&1
        fi
    done
fi

