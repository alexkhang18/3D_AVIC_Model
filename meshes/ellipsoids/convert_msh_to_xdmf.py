import meshio
import numpy as np 

msh = meshio.read("./process2.msh")
# print(msh.cells)
for cell in msh.cells:
		if cell.type == 'triangle':
			triangle_cells = cell.data
		elif cell.type =='tetra':
			tetra_cells = cell.data


# Get physical labels
for key in msh.cell_data_dict["gmsh:physical"].keys():
	if key == 'triangle':
		triangle_data = msh.cell_data_dict["gmsh:physical"][key]
	elif key == 'tetra':
		tetra_data = msh.cell_data_dict["gmsh:physical"][key]
		tetra_data2 = msh.cell_data_dict["gmsh:geometrical"][key]

gel_num_ele = np.max(np.array(np.where(tetra_data == 100)).shape)
cell_num_ele = np.max(np.array(np.where(tetra_data == 200)).shape)
nuc_num_ele = np.max(np.array(np.where(tetra_data == 300)).shape)
total_num_ele = gel_num_ele + cell_num_ele + nuc_num_ele

tetra_mesh = meshio.Mesh(points=msh.points,
	cells=[("tetra",tetra_cells)],
	cell_data={"gmsh:physical": [tetra_data],"gmsh:geometrical": [tetra_data2]})

triangle_mesh = meshio.Mesh(points=msh.points,
	cells=[("triangle", triangle_cells)],
	cell_data = {"triangle":[triangle_data]})


meshio.write("ellipsoidal_cell.xdmf", tetra_mesh)

with open("mesh_num_ele.txt", "w+") as f:
	f.write("No. Elements Total: {:d}\n".format(total_num_ele))
	f.write("No. Elements Gel: {:d}\n".format(gel_num_ele))
	f.write("No. Elements Cell: {:d}\n".format(cell_num_ele))
	f.write("No. Elements Nucleus: {:d}\n".format(nuc_num_ele))

