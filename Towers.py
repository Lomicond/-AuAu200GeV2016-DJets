import os

folder_path = 'production/'  # Upravte cestu k vaší složce

# Získání seznamu souborů v dané složce
files = os.listdir(folder_path)

# Inicializace proměnné pro nejvyšší číslo
max_number = 0

# Inicializace proměnné pro chybějící čísla
missing_numbers = set()

# Procházení souborů a hledání nejvyššího čísla
for file_name in files:
    if file_name.endswith('.root') and '1_Tower_' in file_name:
        number_str = file_name.split('_')[2].split('.')[0]
        number = int(number_str)
        if number > max_number:
            max_number = number

# Kontrola chybějících čísel
for i in range(1, max_number + 1):
    if i not in [int(file.split('_')[2].split('.')[0]) for file in files if file.endswith('.root') and '1_Tower_' in file]:
        missing_numbers.add(i)

# Počet chybějících souborů
missing_files_count = len(missing_numbers)

# Výpis nejvyššího čísla
print("Nejvyšší číslo:", max_number)

# Výpis chybějících čísel
if missing_numbers:
    print("Chybějící čísla:", sorted(missing_numbers))
else:
    print("Žádná čísla nechybí.")

# Výpis počtu chybějících souborů
print("Počet chybějících souborů:", missing_files_count)

