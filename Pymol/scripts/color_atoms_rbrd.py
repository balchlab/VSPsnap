import pandas as pd
import numpy as np
from num2words import num2words
from pymol import cmd
from datetime import date, timedelta

import sys,os
sys.path.append(os.getcwd())
import camera_angles

path = '/gpfs/group/balch/tasks/CoV2/IR_FR_slider/minVar'
delta = timedelta(days=1)
cov2_structures = pd.read_excel('../CoV2_structures_ITASSER.xlsx',usecols=[2,3])


colored_path = f"../colored_proteins/RBRD_{date.today().strftime('%m_%d')}"
if not os.path.exists(colored_path):
	os.mkdir(colored_path)

start_date = date(2020, 2, 1)
end_date = date(2020, 10, 29)
while start_date <= end_date:
	#GPR values
	pol_vals = pd.read_csv(f"{path}/{cov2_structures.loc[10,'name']}_minVarZ_3M_rgb_pymol_{start_date.strftime('%m_%d_%y')}.csv")
	pol_vals = pol_vals.fillna(pol_vals.mean())
	nsp8_vals = pd.read_csv(f"{path}/{cov2_structures.loc[7,'name']}_minVarZ_3M_rgb_pymol_{start_date.strftime('%m_%d_%y')}.csv")
	nsp8_vals = nsp8_vals.fillna(nsp8_vals.mean())
	nsp7_vals = pd.read_csv(f"{path}/{cov2_structures.loc[6,'name']}_minVarZ_3M_rgb_pymol_{start_date.strftime('%m_%d_%y')}.csv")
	nsp7_vals = nsp7_vals.fillna(nsp7_vals.mean())

	#Structure
	cmd.load(f"../structures/6yyt.pdb")

	#Color pol
	cmd.select('All_CA_CHAIN_A','chain A and name CA')
	cmd.iterate('All_CA_CHAIN_A',"cmd.set_color(f'{num2words(resi)}_pol',list(np.floor(pol_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
	cmd.alter('All_CA_CHAIN_A',"color = cmd.get_color_index(f'{num2words(resi)}_pol')")
	#Color nsp8
	cmd.select('All_CA_CHAIN_B_D','(chain B chain D) and name CA')
	cmd.iterate('All_CA_CHAIN_B_D',"cmd.set_color(f'{num2words(resi)}_nsp8',list(np.floor(nsp8_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
	cmd.alter('All_CA_CHAIN_B_D',"color = cmd.get_color_index(f'{num2words(resi)}_nsp8')")
	#Color nsp7
	cmd.select('All_CA_CHAIN_C','chain C and name CA')
	cmd.iterate('All_CA_CHAIN_C',"cmd.set_color(f'{num2words(resi)}_nsp7',list(np.floor(nsp7_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
	cmd.alter('All_CA_CHAIN_C',"color = cmd.get_color_index(f'{num2words(resi)}_nsp7')")

	cmd.set_view (camera_angles.RbRd)
	cmd.center()
	cmd.pseudoatom('forlabel',label=f"{start_date.strftime('%m_%d_%y')}")
	cmd.set('label_relative_mode',1)
	cmd.set('label_screen_point',[-0.75,-0.8,1])
	cmd.set('label_size',25)
	cmd.bg_color('white')
	cmd.deselect()

	cmd.png(f"{colored_path}/rbrd_{start_date.strftime('%m_%d_%y')}.png",
		width=1200, height=1200, dpi=1000, ray=1)
	cmd.delete('all')
	start_date += delta


# #GPR values
# grp_pol = 'Pol_minVarZ_3M_rgb_pymol'
# gpr_nsp8 = 'nsp8_minVarZ_3M_rgb_pymol'
# gpr_nsp7 = 'nsp7_minVarZ_3M_rgb_pymol'
# pol_vals = pd.read_csv(f"../gpr_values/3M/{grp_pol}.csv")
# nsp8_vals = pd.read_csv(f"../gpr_values/3M/{gpr_nsp8}.csv")
# nsp7_vals = pd.read_csv(f"../gpr_values/3M/{gpr_nsp7}.csv")

# #Structure
# structure_file = '6yyt.pdb'
# cmd.load(f"../structures/{structure_file}")

# #resi is for residual number, can't create a color name using an integer hence convert to word.
# #creates rbg color for residual number
# #resi -> is for residual number string
# #resv -> is for residual number as integer
# #Color pol
# cmd.select('All_CA_CHAIN_A','chain A and name CA')
# cmd.iterate('All_CA_CHAIN_A',"cmd.set_color(f'{num2words(resi)}_{grp_pol}',list(np.floor(pol_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
# cmd.alter('All_CA_CHAIN_A',"color = cmd.get_color_index(f'{num2words(resi)}_{grp_pol}')")

# #Color nsp8
# cmd.select('All_CA_CHAIN_B_D','(chain B chain D) and name CA')
# cmd.iterate('All_CA_CHAIN_B_D',"cmd.set_color(f'{num2words(resi)}_{gpr_nsp8}',list(np.floor(nsp8_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
# cmd.alter('All_CA_CHAIN_B_D',"color = cmd.get_color_index(f'{num2words(resi)}_{gpr_nsp8}')")

# #Color nsp7
# cmd.select('All_CA_CHAIN_C','chain C and name CA')
# cmd.iterate('All_CA_CHAIN_C',"cmd.set_color(f'{num2words(resi)}_{gpr_nsp7}',list(np.floor(nsp7_vals.loc[resv-1,['R','G','B']]*255).astype(int)))")
# cmd.alter('All_CA_CHAIN_C',"color = cmd.get_color_index(f'{num2words(resi)}_{gpr_nsp7}')")

# cmd.set_view (camera_angles.RbRd)
# cmd.center()
# cmd.pseudoatom('forlabel',label=f"{start_date.strftime('%m_%d_%y')}")
# cmd.set('label_relative_mode',1)
# cmd.set('label_screen_point',[-0.75,-0.8,1])
# cmd.set('label_size',25)
# cmd.bg_color('white')
# cmd.deselect()


# #Fix below for correct path
# image_out = '6xdc_12_10'
# #cmd.png(f"../colored_proteins/3M/{image_out}.png", dpi=1000, ray=1)
# #cmd.png(f"../colored_proteins/3M/{image_out}.png",width=1200, height=1200, dpi=1000, ray=1)
# #cmd.save(f"../sessions/3M/{image_out}.pse")
# #cmd.delete('all')

