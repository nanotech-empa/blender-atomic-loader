def create_sphere(pos, diameter):
    import bpy
    import bmesh
    import mathutils
    
    # Create an empty mesh and the object.
    mesh = bpy.data.meshes.new('Basic_Sphere')
    basic_sphere = bpy.data.objects.new("Basic_Sphere", mesh)
    
    # Add the object into the scene.
    bpy.context.collection.objects.link(basic_sphere)
    
    # Select the newly created object
    bpy.context.view_layer.objects.active = basic_sphere
    basic_sphere.select_set(True)
    
    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()

    # create a location matrix
    mat_loc = mathutils.Matrix.Translation(pos)
    bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, diameter=diameter, matrix=mat_loc)
    bm.to_mesh(mesh)
    bm.free()
    
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.ops.object.shade_smooth()

def cylinder_between(point1, point2, radius):
  
    import math
    import bpy
    import bmesh

    x1=point1[0]
    y1=point1[1]
    z1=point1[2]
    x2=point2[0]
    y2=point2[1]
    z2=point2[2]

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1    
    dist = math.sqrt(dx**2 + dy**2 + dz**2)

    bpy.ops.mesh.primitive_cylinder_add(
        radius = radius, 
        depth = dist,
        location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)   
    ) 

    phi = math.atan2(dy, dx) 
    theta = math.acos(dz/dist) 

    bpy.context.object.rotation_euler[1] = theta 
    bpy.context.object.rotation_euler[2] = phi 

def get_molecule(aseframe):
    import numpy as np
    
    mol=aseframe[np.where(np.asarray(aseframe.get_chemical_symbols())!='Au')[0]]
    # check for possible H passivating gold, far from the molecule
    # first get the z-coordinate of center of mass of the carbon backbone
    cmass_mol=int(mol[np.where(np.asarray(mol.get_chemical_symbols())=='C')[0]].get_center_of_mass()[2])
    # get hydrogens z-coordinate and cast it to an integer 
    # (to round up layers and get rid of fluctuations)
    id_hydrogens=np.where(np.asarray(mol.get_chemical_symbols())=='H')[0]
    hydrogens_z=np.asarray(list(map(int,mol[id_hydrogens].positions[:,2])))
    # discard the hydrogens far away
    id_mol=np.asarray([i for i,at in enumerate(mol) if i not in id_hydrogens[np.where(np.abs(hydrogens_z-cmass_mol)>3)[0]]])
    return mol[id_mol]

def get_substrate(aseframe):
    import numpy as np
    return aseframe[np.where(np.asarray(aseframe.get_chemical_symbols())=='Au')[0]]

def get_neighbours(aseframe):
    # get the neighbours of each atom in the molecule
    from scipy.spatial import cKDTree
    # get the kd-tree for nearest neighbors
    kdtree=cKDTree(aseframe.positions)
    neigh_list={}
    for i,pos in enumerate(aseframe.positions):
        # Query the kd-tree and get the first 7 neighbours
        _,neighs=kdtree.query(pos,k=8)
        # exclude the 1st element (the point itself)
        neigh_list[i]=neighs[1:]
    return neigh_list

def get_bonds(aseframe,cutoff=1.5):
    import numpy as np
    # get the distance matrix
    dist_mat=aseframe.get_all_distances()
    bond_list={}
    for i,dists in enumerate(dist_mat):
        bond_list[i]=[j for j in np.where(dists<1.5)[0] if i!=j]
    return bond_list

def get_bonds(aseframe,cutoff=1.5):
    import numpy as np
    # get the bonds list
    
    # get the distance matrix
    dist_mat=aseframe.get_all_distances()
    bond_list=[]
    for i,dists in enumerate(dist_mat):
        bond_list.append([np.sort([i,j]) for j in np.where(dists<1.5)[0] if i!=j])
    # returne the list of unique bonds
    return np.asarray(list(set([tuple(i) for i in np.concatenate(bond_list)])))

def split_bonds(aseframe,cutoff=1.5):
    import numpy as np
    # split the bonds with hydrogens from those without
    tot_bonds=get_bonds(aseframe,cutoff)
    # find all the bonds with hydrogens
    ll=[]
    for bond_atoms in tot_bonds:
        atom_types=[aseframe[j].symbol for j in bond_atoms]
        ll.append('H' in atom_types)
    ll=np.asarray(ll)
    id_hbonds=np.where(ll==True)[0]
    id_nohbonds=np.where(ll==False)[0]
    return tot_bonds[id_hbonds],tot_bonds[id_nohbonds]

def draw_type(aseframe,atom_type,radious):
    import numpy as np
    # draw ths spheres for a specific atom type
    for pos in aseframe[np.where(np.asarray(aseframe.get_chemical_symbols())==atom_type)[0]].positions:
        create_sphere(pos, radious*2)
