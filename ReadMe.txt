LabelPhy

options:
  -h, --help  show this help message and exit
  -f F        Input .tre File
  -v          Verbose / Iterate Every Node
  -e E        Email for NCBI API
  -k K        API Key for NCBI API
  -s          Save NCBI Email and API Key
  -o O        Output File Name

python labelphy.py -f Rooted.tre -o Output_tree.tre -e NCBI_EMAIL -k NCBI_KEY -s

-f,-o 	Pick a tree, Place a tree

-e,-k 	Email and API key for NCBI

-s 	Will save a txt in your directory with your email and api key so you don't have to type it out every time

-v 	Verbose is for if you're insane enough to double check every single node and its children

-c	Same as -v Verbose but with the option to enter a custom node label

The first run of a tree will take awhile because this script will search NCBI for every individual species.
It takes this data and saves it to a .json
So, every following run will be much faster

Script is currently extremely barebones. Think of it more as a proof-of-concept at the moment.
