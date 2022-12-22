#!/urs/bin/python3

import os
from . import tools

this_directory = os.path.abspath(os.path.dirname(__file__))
supplement_path = '{}/supplements/supplementary_data.json'.format(this_directory)

supplementary_data = tools.open_json(supplement_path)
not_r_groups = supplementary_data['not_r_groups']
standard_bond_counts = supplementary_data['standard_bond_counts']
atomic_weights = supplementary_data['atomic_weights']
metal_symbols = supplementary_data['metal_symbols']
index_to_charge = supplementary_data['index_to_charge']
# charge_to_index = supplementary_data['charge_to_index']
