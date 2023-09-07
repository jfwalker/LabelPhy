from species_info import SpeciesInfo
import json

if __name__ == '__main__':

    TREE = open('Rooted.tre', 'r')
    species_tax_dict = {}

    def species_tax(species_list):
        SI = SpeciesInfo(email='ebretz2@uic.edu', key='03085733d8ef177e9eab46c5b9ab1eb84909')
        species_tax = []
        overlap_dict = {}

        # Try to load species_tax_dict from file
        try:
            with open('species_tax_dict.json', 'r') as f:
                species_tax_dict = json.load(f)
        except FileNotFoundError:
            species_tax_dict = {}

        for species in species_list:
            if species not in species_tax_dict:
                species_tax_dict[species] = SI.get_rank_dict(species)
                # Save updated species_tax_dict to file
                with open('species_tax_dict.json', 'w') as f:
                    json.dump(species_tax_dict, f)
            species_tax.append(species_tax_dict[species])
            
        for x in species_tax:
            for dct in species_tax:
                for key, value in dct.items():
                    if key.lower() != 'clade':
                        if key in overlap_dict:
                            overlap_dict[key].append(value)
                        else:
                            overlap_dict[key] = [value]
        
        for key in list(overlap_dict.keys()):
                overlap_dict[key] = list(set.intersection(*map(set, overlap_dict[key])))

        ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus']
        overlap_dict = {key: val for key, val in overlap_dict.items() if val and key.lower() in ranks}
        
            
            
        return overlap_dict

    with TREE as tree:
        start_pos = 0
        count = 0
        pos = 0
        species_start = False
        species_list = []
        prev_char = ''
        for line in tree:
            for char in line:
                pos += 1
                if char == '(' and not start_pos != 0:
                    start_pos = pos
                    count += 1
                elif char == '(' and start_pos:
                    count += 1
                elif char == ')':
                    count -= 1
                    lowest_common = list(species_tax(species_list).items())[-1]
                    print(lowest_common)
                elif char.isalpha() or char == '_' and not prev_char.isnumeric():
                    if not species_start:
                        species_start = True
                        species_name = char
                    else:
                        species_name += char
                elif char == '.' and prev_char.isalpha():
                    if not species_start:
                        species_start = True
                        species_name = char
                    else:
                        species_name += char
                else:
                    if char == '-' and prev_char == 'e':
                        species_name = ''
                        species_start = False

                        continue
                    elif species_start:
                        species_list.append(species_name)
                        species_start = False
                prev_char = char
