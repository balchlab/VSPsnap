Collection of scripts to run Pymol from the Scripps server. 

### A few pointers.
- Scripts show general outline of selecting atoms within a structure and assigning the color through gpr output files.
- File paths and names (for gpr value files, pdb structures and output diretory) will have to be edited for custom use.
- When running on the server, one first needs to load the Pymol module (module load pymol).
- A few Python libraries may be needed to download. For example num2words (pip install --user num2words, intstalls Python libraries to local user when using server version of Python).
- If wanting to use a custom camera angle/view point. First set the camera viewpoint in Pymol visually through the interactive/desktop version of Pymol. Then export the orientation matrix following https://pymolwiki.org/index.php/Get_View. Copy and paste the camera view matrix into camera_angles.py and assign with cmd.set_view (camera_angles.RbRd) for example.
