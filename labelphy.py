import re
import sys
from species_info import SpeciesInfo
import json
import sys
import argparse
import os


if __name__ == '__main__':

    species_tax_dict = {}
    SI = SpeciesInfo()


    def species_tax(species_list):
        species_tax = []
        overlap_dict = {}

        species_tax_dict = load_species_tax_dict()

        for species in species_list:
            species_tax_dict = update_species_tax_dict(species, species_tax_dict, SI)
            species_tax.append(species_tax_dict[species])
        
        

        overlap_dict = update_overlap_dict(species_tax, overlap_dict)
        
        ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus',]
        overlap_dict = filter_overlap_dict(overlap_dict, ranks)
        return overlap_dict

    def load_species_tax_dict():
        try:
            with open('species_tax_dict.json', 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            return {}

    def update_species_tax_dict(species, species_tax_dict, SI):
        if species not in species_tax_dict:
            sys.stdout.write('\033[F')
            sys.stdout.write('\033[K')
            
            sys.stdout.write(f'\033[1mSearching NCBI for species:\033[0m {species}\n')

            species_tax_dict[species] = SI.get_rank_dict(species)
            with open('species_tax_dict.json', 'w') as f:
                json.dump(species_tax_dict, f)
        return species_tax_dict

    def update_overlap_dict(species_tax, overlap_dict):
        for x in species_tax:
            for dct in species_tax:
                for key, value in dct.items():
                    overlap_dict.setdefault(key, []).append(value)
        for key in list(overlap_dict.keys()):
            overlap_dict[key] = list(set.intersection(*map(set, overlap_dict[key])))
        return overlap_dict

    def filter_overlap_dict(overlap_dict, ranks):
        return {key: val for key, val in overlap_dict.items() if val and key.lower() in ranks}

    def wait_for_keypress():
        input('\033[1mPress any key to continue...\033[0m')
        sys.stdout.write('\033[F')
        sys.stdout.write('\033[K')

    def output(species_list, lowest_common):
        sys.stdout.write(f'\033[0;32m{species_list}\033[34m\n{lowest_common}')
        wait_for_keypress()


    def main():
        stack = []
        species = []
        bracket_positions = []
        close_brackets = []
        lowest_list = []
        with open(f'{args.f}', 'r') as tree:
            for line in tree:
                line = line.rstrip()
                for i, char in enumerate(line):
                    if char == '(':
                        stack.append(i)
                    elif char == ')':
                        if stack:
                            start = stack.pop()
                            bracket_positions.append((start, i))
                            close_brackets.append(i)
                            species.append(re.findall('[a-zA-Z_\.]+(?<=[a-zA-Z])(?!-)', line[start:i]))
      
            species_list = [x for y in species for x in y]
            species_list = list(set(species_list))

            sys.stdout.write(f'\033[1mNumber of Tips:\033[0m {len(species_list)}\n\n')
            for x in species:
                lowest_common = species_tax(x)
                lowest_list.append(list(lowest_common.items())[-1])
            if SI.species_errors:
                sys.stdout.write(f'\033[1mNumber of species errors:\033[0m {len(SI.species_errors)}\n')
            for x,y in reversed(list(enumerate(close_brackets))):
                colon_pos = y
                while line[colon_pos] != ':':
                    if colon_pos == len(line) - 1:
                        break
                    colon_pos += 1
                if args.v:
                    sys.stdout.write(f'\033[1m{species[x]}\033[0m\n')
                    sys.stdout.write(f'\033[32m{lowest_list[x][-1][-1]}\033[m\n')
                    wait_for_keypress()  
                string_to_replace = ')' + lowest_list[x][-1][-1] + ':'
                line = line[:y] + ')' + lowest_list[x][-1][-1] + line[colon_pos:].replace(string_to_replace, '):')

                
            sys.stdout.write(f'\033[1mFinished...\033[0m\nWriting output as {args.o}\n')
            with open(args.o, 'w') as f:
                f.write(line)

    
    parser = argparse.ArgumentParser(description='LabelPhy')
    parser.add_argument('-f', type=str, help='Input .tre File')
    parser.add_argument('-v', action='store_true', help='Verbose / Iterate Every Node')
    parser.add_argument('-e', type=str, help='Email for NCBI API')
    parser.add_argument('-k', type=str, help='API Key for NCBI API')
    parser.add_argument('-s', action='store_true', help='Save NCBI Email and API Key')
    parser.add_argument('-o', type=str, help='Output File Name')
    
    args = parser.parse_args()

    if args.e and args.k:
        SI.email = args.e
        SI.key = args.k
        if args.s:
            with open('api_key.txt', 'w') as f:
                f.write(f'{args.e}\n{args.k}')
        print(SI.email, SI.key)
    else:
        try:
            with open('api_key.txt', 'r') as f:
                SI.email = f.readline().rstrip()
                SI.key = f.readline().rstrip()
        except:
            sys.stdout.write(f'\033[1;31mError:\033[0m No API Key Found\n')
            sys.stdout.write(f'Use -e and -k to specify email and API key\n')
            sys.stdout.write(f'Or use -s to save email and API key\n')
            sys.exit(1)
            
    if args.o:
        dir_name = os.path.dirname(args.o)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)
        
    if args.f:
        main()
    else:
        parser.print_help()
        sys.exit(1)