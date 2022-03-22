import pandas as pd
import numpy as np
from num2words import num2words
from pymol import cmd
from datetime import date, timedelta

import sys,os
sys.path.append(os.getcwd())
import camera_angles

#Path to GRP input values
#path = '/gpfs/group/balch/tasks/CoV2/IR_FR_slider/minVar'
path = '../gpr_values/minvar_7_1_20/minvar_7_1_20'
delta = timedelta(days=1)
cov2_structures = pd.read_excel('../CoV2_structures_ITASSER.xlsx',usecols=[2,3])
structure = '6vxx.pdb'

#Index of single protein to use from excel file
val=15

colored_path = f"../colored_proteins/{cov2_structures.loc[val,'name']}_6vxx_{date.today().strftime('%m_%d')}"
if not os.path.exists(colored_path):
	os.mkdir(colored_path)

start_date = date(2020, 2, 1)
end_date = date(2021, 1, 5)
while start_date <= end_date:
	gpr_vals = pd.read_csv(f"{path}/{cov2_structures.loc[val,'name']}_minVarZ_3M_rgb_pymol_{start_date.strftime('%m_%d_%y')}.csv")
	#There are NA values, replace them with average to smoothe that section
	gpr_vals = gpr_vals.fillna(gpr_vals.mean())
	cmd.load(f"../structures/{structure}")
	cmd.select('All_CA','name CA')
	cmd.iterate('All_CA',"cmd.set_color(f'{num2words(resi)}_{name}',list(np.floor(gpr_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
	cmd.alter('All_CA',"color = cmd.get_color_index(f'{num2words(resi)}_{name}')")

	#Set camera angle from saved matrix
	cmd.set_view (camera_angles.m6vxx)
	cmd.bg_color('white')
	cmd.pseudoatom('forlabel',label=f"{start_date.strftime('%m_%d_%y')}")
	cmd.set('label_relative_mode',1)
	cmd.set('label_screen_point',[-0.75,-0.8,1])
	cmd.set('label_size',25)
	cmd.hide('licorice')
	cmd.deselect()

# 	#Fix below for correct path
	cmd.png(f"{colored_path}/6vxx_{cov2_structures.loc[val,'name']}_{start_date.strftime('%m_%d_%y')}.png",
		width=1200, height=1200, dpi=1000, ray=1)
	cmd.delete('all')
	start_date += delta
