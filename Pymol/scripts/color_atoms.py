import pandas as pd
import numpy as np
from num2words import num2words
from pymol import cmd
from datetime import date, timedelta
import os

path = '/gpfs/group/balch/tasks/CoV2/IR_FR_slider/minVar'
delta = timedelta(days=1)
cov2_structures = pd.read_excel('../CoV2_structures_ITASSER.xlsx',usecols=[2,3])

for val in cov2_structures.index:
	#check path and create if folder doesn't exist
	colored_path = f"../colored_proteins/{cov2_structures.loc[val,'name']}"
	if not os.path.exists(colored_path):
		os.mkdir(colored_path)

	start_date = date(2020, 2, 1)
	end_date = date(2020, 10, 29)
	while start_date <= end_date:
		gpr_vals = pd.read_csv(f"{path}/{cov2_structures.loc[val,'name']}_minVarZ_3M_rgb_pymol_{start_date.strftime('%m_%d_%y')}.csv")
		cmd.load(f"../structures/{cov2_structures.loc[val,'struct']}")

		cmd.select('All_CA','name CA')
		cmd.iterate('All_CA',"cmd.set_color(f'{num2words(resi)}_{name}',list(np.floor(gpr_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
		cmd.alter('All_CA',"color = cmd.get_color_index(f'{num2words(resi)}_{name}')")
		cmd.orient()
		cmd.zoom()
		cmd.deselect()
		cmd.bg_color('white')
		#Fix below for correct path
		cmd.png(f"../colored_proteins/{cov2_structures.loc[val,'name']}/{cov2_structures.loc[val,'name']}_{start_date.strftime('%m_%d_%y')}.png")
		#cmd.png(f"../colored_proteins/test/{cov2_structures.loc[val,'name']}_{start_date.strftime('%m_%d_%y')}.png")
		#cmd.save(f"../sessions/test/{cov2_structures.loc[val,'name']}.pse")
		cmd.delete('all')
		start_date += delta