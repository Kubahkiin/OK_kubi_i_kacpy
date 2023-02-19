import gen

jak_dlugie = 300

instance = gen.dna(jak_dlugie)

gen.save_to(f'dna_{jak_dlugie}.txt', [instance])
