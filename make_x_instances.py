import gen

ile = 10
jak_dlugie = 100
instances = []

for i in range(ile):
    instances.append(gen.dna(jak_dlugie))

gen.save_to("instances.txt", instances)