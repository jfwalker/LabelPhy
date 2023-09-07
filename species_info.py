from Bio import Entrez
import sys

class SpeciesInfo:
    def __init__(self, email: str = '', key: str = '', species_name: str = ''):
        self.key = key
        self.email = email
        self.species_name = species_name
        self.species_id = ''
        self.species_rank = ''
        self.species_rank_dict = {}
        self.species_lineage = ''
        self.species_record = {}
        self.species_record_count = 0
        self.species_taxonomy = {}
        self.species_errors = []

    def count_ranks_in_taxonomy(self):
        rank_counts = {}
        for item in self.species_taxonomy:
            rank = item['Rank']
            scientific_name = item['ScientificName']
            rank_counts.setdefault(rank, []).append(scientific_name)
        return rank_counts

    def search_species_name(self):
        Entrez.email = self.email
        Entrez.api_key = self.key
        try:
            handle = Entrez.esearch(db="taxonomy", term=self.species_name.replace('_', ' '))
            record = Entrez.read(handle)
            return record

        except Exception as e:
            sys.stdout.write(f'\033[1;31mError searching NCBI for species ID:\033[0m {self.species_name}: {e}\n')
            return None

    def search_correction(self):
        species_name_list = self.species_name.split('_')
        sys.stdout.write(f'Retrying {self.species_name} as {" ".join(species_name_list[:-1])}...\n')
        try:
            handle = Entrez.esearch(db="taxonomy", term=' '.join(species_name_list[:-1]))
            record = Entrez.read(handle)
            # sys.stdout.write(f'\033[;32mSuccess.\033[0m\n\n')
            return get_species_id()
        except:
            sys.stdout.write(f'\033[;31mFailed.\033[0m\n')
            sys.stdout.write(f'Retrying {self.species_name} as {" ".join(species_name_list[1:-1])}...\n')
            try:
                handle = Entrez.esearch(db="taxonomy", term=' '.join(species_name_list[1:-1]))
                record = Entrez.read(handle)
                sys.stdout.write(f'\033[;32mSuccess.\033[0m\n')
                return get_species_id()
            except:
                return None

    def get_species_id(self):
        Entrez.email = self.email
        Entrez.api_key = self.key
        try:
            if not self.species_id:
                self.species_id = self.search_species_name()['IdList'][0]
        except Exception as e:
            sys.stdout.write(f'\033[1;31mError getting species ID:\033[0m {self.species_name}: {e}\n')
            sys.stdout.write(f'Retrying {self.species_name}...\n')
            try:
                self.species_id = self.search_species_name()['IdList'][0]
                sys.stdout.write(f'\033[;32mSuccess.\033[0m\n')
            except:
                sys.stdout.write(f'\033[;31mFailed.\033[0m\n')
                self.species_errors.append((self.species_name))
                search_correction = self.search_correction()
                return search_correction

        try:    
            handle = Entrez.efetch(db="taxonomy", id=self.species_id, retmode="xml")
            record = Entrez.read(handle)
            return record
        except Exception as e:
            sys.stdout.write(f'\033[1;31mError searching NCBI for species data:\033[0m {self.species_name}: {e}\n')
            sys.stdout.write(f'Retrying {self.species_name}...\n')
            try:
                handle = Entrez.efetch(db="taxonomy", id=self.species_id, retmode="xml")
                record = Entrez.read(handle)
                sys.stdout.write(f'\033[;32mSuccess.\033[0m\n\n')
                return record
            except:
                sys.stdout.write(f'\033[;31mFailed.\033[0m\n\n')
                return None
            
    def set_species_info(self):
        search_result = self.search_species_name()
        if search_result:
            try:
                self.species_id = search_result['IdList'][0]
            except IndexError:
                self.species_id = None
            species_info = self.get_species_id()
            if species_info:
                self.species_record = species_info
                self.species_record_count = len(species_info)
                self.species_rank = species_info[0]['LineageEx'][0]['Rank']
                self.species_taxonomy = species_info[0]['LineageEx']
                self.species_rank_dict = self.count_ranks_in_taxonomy()

            else:
                sys.stdout.write('No species info found.\n')
        else:
            sys.stdout.write('No search result found.\n')

    def get_rank_dict(self, species_name):
        self.species_name = species_name
        self.set_species_info()
        return self.species_rank_dict

    