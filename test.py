import json

def load_species_tax_dict():
    try:
        with open('species_tax_dict.json', 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

species_tax_dict = load_species_tax_dict()
bad = []
no_phylum = []
Pandalanales = []
for key,val in species_tax_dict.items():
    if 'Pandalanales' in key:
        Pandalanales.append((key, val))
    try:
        if 'Streptophyta' not in val['phylum']:
              bad.append((key, val['phylum']))
    except:
        pass
    try:
        if 'phylum' not in val:
            no_phylum.append((key, val))
    except:
        pass

print('No Streptophyta', bad)
for x in bad:
    print(x[0].replace('_', ' '))

# for x in bad:
#     species_tax_dict.pop(x[0], None)
#
# with open('species_tax_dict.json', 'w') as f:
#     json.dump(species_tax_dict, f)

print('No Phylum', no_phylum)
for x in no_phylum:
    print(x[0].replace('_', ' '))
#
# for x in no_phylum:
#     species_tax_dict.pop(x[0], None)
#
# with open('species_tax_dict.json', 'w') as f:
#     json.dump(species_tax_dict, f)

# print('Pandalanales', Pandalanales)
for x in Pandalanales:
    print(x[0].replace('_', ' '))
    print(x[1])

# for x in Pandalanales:
#     species_tax_dict.pop(x[0], None)

# with open('species_tax_dict.json', 'w') as f:
#     json.dump(species_tax_dict, f)