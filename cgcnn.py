import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

import networkx as nx
import seaborn as sns
import os
import ase.io
import pandas as pd
from tqdm import tqdm_notebook
from aflow import *
import itertools
import urllib3
http = urllib3.PoolManager()
#just getting simulations with slight magnetism 
result = search(batch_size=10000).filter(K.spin_atom>0.25).filter(K.natoms>3).select(K.energy_atom, 
                                                                                       K.spin_atom,
                                                                                       K.species)

onic Radii from
# http://abulafia.mt.ic.ac.uk/shannon/ptable.php
# taking approx. average of ionic radii
R_ionic = {"H": -0.16, 
            "He": -1., 
            "Li": 0.9,
            "Be": 0.41,
            "B": 0.19,
            "C": 0.25,
            "N": 1.0,
            "O": 1.38,
            "F": 1.3,
            "Ne": -1,
            "Na": 1.12,
            "Mg": 0.8,
            "Al": 0.62,
            "Si": 0.4,
            "P": 0.43,
            "S": 0.6,
            "Cl": 0.7,
            "Ar": -1.,
            "K": 1.6,
            "Ca": 1.29,
            "Sc": 0.97,
            "Ti": 0.69,
            "V": 0.71,
            "Cr": 0.74,
            "Mn": 0.8,
            "Fe": 0.82,
            "Co": 0.75,
            "Ni": 0.77,
            "Cu": 0.78,
            "Zn": 0.86,
            "Ga": 0.71,
            "Ge": 0.68,
            "As": 0.59,
            "Se": 0.7,
            "Br": 0.8,
            "Kr": -1.,
            "Rb": 1.61,
            "Sr": 1.37,
            "Y": 1.11,
            "Zr": 0.89,
            "Nb": 0.82,
            "Mo": 0.66,
            "Tc": 0.58,
            "Ru": 0.701,
            "Rh": 0.7401,
            "Pd": 0.82,
            "Ag": 1.1,
            "Cd": 1.09,
            "In": 0.94,
            "Sn": 0.83,
            "Sb": 0.87,
            "Te": 1.5,
            "I": 1.05,
            "Xe": 0.49,
            "Cs": 1.811,
            "Ba": 1.51,
            "La":1.35 ,
            "Ce": 1.39,
            "Pr": 1.233,
            "Nd": 1.241,
            "Pm": 1.201,
            "Sm": 1.19,
            "Eu": 1.31,
            "Gd": 1.13,
            "Tb": 1.09,
            "Dy": 1.18,
            "Ho":1.25 ,
            "Er": 1.111,
            "Tm": 1.17,
            "Yb": 1.22,
            "Lu": 1.07,
            "Hf": 0.801,
            "Ta": 0.804,
            "W": 0.74,
            "Re": 0.72,
            "Os": 0.688,
            "Ir": 0.765,
            "Pt": 0.87,
            "Au": 1.1,
            "Hg": 1.12,
            "Tl": 1.4,
            "Pb": 1.3,
            "Bi": 1.105,
            "Po": 1.08,
            "At": 0.76,
            "Rn": -0.1,
            "Fr": 1.94,
            "Ra": 1.7,
            "Ac": 1.26,
            "Th": 1.19,
            "Pa": 1.04,
            "U": 0.99,
            "Np": 1.07,
            "Pu": 1.022,
            "Am": 1.3,
            "Cm": 1.067,
            "Bk": 1.056,
            "Cf": 1.045,
            "Es": -0.1,
            "Fm": -0.1,
            "Md": -0.1,
            "No": 1.24,
            "Lr": -0.1,
            "Rf": -0.1,
            "Db": -0.1,
            "Sg": -0.1,
            "Bh": -0.1,
            "Hs": -0.1,
            "Mt": -0.1,
            "Ds": -0.1,
            "Rg": -0.1,
            "Uub": -0.1,
            "Uut": -0.1,
            "Uuq": -0.1,
            "Uup": -0.1,
            "Uuh": -0.1,
            "Uus": -0.1,
            "Uuo": -0.1}

group_num = {"H": 1, 
            "He": 18., 
            "Li": 1,
            "Be": 2,
            "B": 13,
            "C": 14,
            "N": 15,
            "O": 16,
            "F": 17,
            "Ne": 18,
            "Na": 1,
            "Mg": 2,
            "Al": 13,
            "Si": 14,
            "P": 15,
            "S": 16,
            "Cl": 17,
            "Ar": 18.,
            "K": 1,
            "Ca": 2,
            "Sc": 3,
            "Ti": 4,
            "V": 5,
            "Cr": 6,
            "Mn": 7,
            "Fe": 8,
            "Co": 9,
            "Ni": 10,
            "Cu": 11,
            "Zn": 12,
            "Ga": 13,
            "Ge": 14,
            "As": 15,
            "Se": 16,
            "Br": 17,
            "Kr": 18,
            "Rb": 1,
            "Sr": 2,
            "Y": 3,
            "Zr": 4,
            "Nb": 5,
            "Mo": 6,
            "Tc": 7,
            "Ru": 8,
            "Rh": 9,
            "Pd": 10,
            "Ag": 11,
            "Cd": 12,
            "In": 13,
            "Sn": 14,
            "Sb": 15,
            "Te": 16,
            "I": 17,
            "Xe": 18,
            "Cs": 1,
            "Ba": 2,
            "La": 4,
            "Ce": 5,
            "Pr": 6,
            "Nd": 7,
            "Pm": 8,
            "Sm": 9,
            "Eu": 10,
            "Gd": 11,
            "Tb": 12,
            "Dy": 13,
            "Ho":14 ,
            "Er": 15,
            "Tm": 16,
            "Yb": 17,
            "Lu": 18,
            "Hf": 4,
            "Ta": 5,
            "W": 6,
            "Re": 7,
            "Os": 8,
            "Ir":9,
            "Pt": 10,
            "Au": 11,
            "Hg": 12,
            "Tl": 13,
            "Pb": 14,
            "Bi": 15,
            "Po": 16,
            "At": 17,
            "Rn": 18,
            "Fr": 1,
            "Ra": 2,
            "Ac": 4,
            "Th": 5,
            "Pa": 6,
            "U": 7,
            "Np": 8,
            "Pu": 9,
            "Am": 10,
            "Cm": 11,
            "Bk": 12,
            "Cf": 13,
            "Es": 14,
            "Fm": 15,
            "Md": 16,
            "No": 17,
            "Lr": 18,
            "Rf": 4,
            "Db": 5,
            "Sg": 6,
            "Bh": 7,
            "Hs": 8,
            "Mt": 9,
            "Ds": 10,
            "Rg": 11
            }

period_num = {"H": 1, 
            "He": 1., 
            "Li": 2,
            "Be": 2,
            "B": 2,
            "C": 2,
            "N": 2,
            "O": 2,
            "F": 2,
            "Ne":2,
            "Na": 3,
            "Mg": 3,
            "Al": 3,
            "Si": 3,
            "P": 3,
            "S": 3,
            "Cl": 3,
            "Ar": 3.,
            "K": 4,
            "Ca": 4,
            "Sc": 4,
            "Ti": 4,
            "V": 4,
            "Cr": 4,
            "Mn": 4,
            "Fe": 4,
            "Co": 4,
            "Ni": 4,
            "Cu": 4,
            "Zn": 4,
            "Ga": 4,
            "Ge": 4,
            "As": 4,
            "Se": 4,
            "Br": 4,
            "Kr": 4,
            "Rb": 5,
            "Sr": 5,
            "Y": 5,
            "Zr": 5,
            "Nb": 5,
            "Mo": 5,
            "Tc": 5,
            "Ru": 5,
            "Rh": 5,
            "Pd": 5,
            "Ag": 5,
            "Cd": 5,
            "In": 5,
            "Sn": 5,
            "Sb": 5,
            "Te": 5,
            "I": 5,
            "Xe":5,
            "Cs": 6,
            "Ba": 6,
            "La": 8,
            "Ce": 8,
            "Pr": 8,
            "Nd": 8,
            "Pm": 8,
            "Sm": 8,
            "Eu": 8,
            "Gd": 8,
            "Tb": 8,
            "Dy": 8,
            "Ho": 8,
            "Er": 8,
            "Tm": 8,
            "Yb": 8,
            "Lu": 8,
            "Hf": 6,
            "Ta": 6,
            "W": 6,
            "Re": 6,
            "Os": 6,
            "Ir":6,
            "Pt": 6,
            "Au": 6,
            "Hg": 6,
            "Tl": 6,
            "Pb": 6,
            "Bi": 6,
            "Po": 6,
            "At": 6,
            "Rn": 6,
            "Fr": 7,
            "Ra": 7,
            "Ac": 9,
            "Th": 9,
            "Pa": 9,
            "U": 9,
            "Np": 9,
            "Pu": 9,
            "Am": 9,
            "Cm": 9,
            "Bk": 9,
            "Cf": 9,
            "Es": 9,
            "Fm": 9,
            "Md": 9,
            "No": 9,
            "Lr": 9,
            "Rf": 7,
            "Db": 7,
            "Sg": 7,
            "Bh": 7,
            "Hs": 7,
            "Mt": 7,
            "Ds": 7,
            "Rg": 7
            }

electroneg_dict = {"H": 6, 
            "He": 0, 
            "Li": 1,
            "Be": 3,
            "B": 5,
            "C": 7,
            "N": 8,
            "O": 9,
            "F": 9,
            "Ne":0,
            "Na": 1,
            "Mg": 3,
            "Al": 4,
            "Si": 5,
            "P": 5,
            "S": 7,
            "Cl": 9,
            "Ar": 0,
            "K": 1,
            "Ca": 2,
            "Sc": 3,
            "Ti": 3,
            "V": 4,
            "Cr": 4,
            "Mn": 3,
            "Fe": 4,
            "Co": 4,
            "Ni": 5,
            "Cu": 5,
            "Zn": 4,
            "Ga": 4,
            "Ge": 5,
            "As": 5,
            "Se": 7,
            "Br": 8,
            "Kr": 8,
            "Rb": 1,
            "Sr": 1,
            "Y": 2,
            "Zr": 3,
            "Nb": 4,
            "Mo": 5,
            "Tc": 5,
            "Ru": 6,
            "Rh": 6,
            "Pd": 6,
            "Ag": 5,
            "Cd": 4,
            "In": 4,
            "Sn": 5,
            "Sb": 5,
            "Te": 5,
            "I": 7,
            "Xe":7,
            "Cs": 1,
            "Ba": 1,
            "La": 2,
            "Ce": 2,
            "Pr": 2,
            "Nd": 2,
            "Pm": 2,
            "Sm": 2,
            "Eu": 2,
            "Gd": 2,
            "Tb": 2,
            "Dy": 2,
            "Ho": 2,
            "Er": 2,
            "Tm": 2,
            "Yb": 2,
            "Lu": 2,
            "Hf": 3,
            "Ta": 3,
            "W": 6,
            "Re": 5,
            "Os": 6,
            "Ir":6,
            "Pt": 6,
            "Au": 7,
            "Hg": 5,
            "Tl": 4,
            "Pb": 6,
            "Bi": 5,
            "Po": 5,
            "At": 6,
            "Rn": 0,
            "Fr": 1,
            "Ra": 1,
            "Ac": 2,
            "Th": 3,
            "Pa": 3,
            "U": 3,
            "Np": 3,
            "Pu": 2,
            "Am": 3,
            "Cm": 3,
            "Bk": 3,
            "Cf": 3,
            "Es": 3,
            "Fm": 3,
            "Md": 3,
            "No": 3,
            "Lr": 0,
            "Rf": -1,
            "Db": -1,
            "Sg": -1,
            "Bh": -1,
            "Hs": -1,
            "Mt": -1,
            "Ds": -1,
            "Rg": -1
            }

def length(v):
    return np.linalg.norm(v)

def unit_vector(vector):
    return vector / length(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle_deg_between(v1, v2):
    return np.degrees(angle_between(v1, v2))

def get_lattice_constants(lattice_vectors):
    lat_const_series = pd.Series()
    for i in range(3):
        lat_const_series["lattice_vector_"+str(i+1)+"_ang"] = length(lattice_vectors[i])
    lat_const_series["lattice_angle_alpha_degree"] = angle_deg_between(lattice_vectors[1],lattice_vectors[2])
    lat_const_series["lattice_angle_beta_degree"] = angle_deg_between(lattice_vectors[2],lattice_vectors[0])
    lat_const_series["lattice_angle_gamma_degree"] = angle_deg_between(lattice_vectors[0],lattice_vectors[1])
    return lat_const_series

distance_categories = np.linspace(0.7,5.6,11)
    
electroneg_cats = np.arange(0,10)
group_numcats = np.arange(0,19)
period_num_cats = np.arange(1,10)

def get_crytal_graph(reduced_coords, dists, mat_size=100):
    #matrix = []
    
    natom = len(reduced_coords)
    G = nx.Graph()
    
    #The 'factor' variable determines how likely  
    #it is there will be more vertices in the graph
    if(natom==4):
        factor=2.5
    elif(natom<8):
        factor=2.1
    elif(natom<=32):
        factor=1.1
    elif(natom<80):
        factor=1.05
    
    elif(natom<100):
        factor=1.0
    else:
        factor=0.9
    for i in range(natom):
        symbol_i = reduced_coords[i][1]
        for j in range(i):
            symbol_j = reduced_coords[j][1]
            #if (symbol_i == "O" and symbol_j != "O") or (symbol_i != "O" and symbol_j == "O"):
            node_i = symbol_i + "_" + str(i) 
            node_j = symbol_j + "_" + str(j) 
            R_max = (R_ionic[symbol_i] + R_ionic[symbol_j]) * factor
            if dists[i, j] < R_max:
                G.add_edge(node_i, node_j, weight=dists[i,j])
    
    #print (len(G.nodes()))
    #matrix = -np.ones([len(G.nodes()),900])
    node_list = []
    for i, node in enumerate(G.nodes()):
        
        crdn = G.neighbors(node)
        #print (len(crdn))
        #have to one_hot_encode each node! Boring but need a dataset here.........
        neighbour_list = []
        for j, neighbour in enumerate(crdn):
            distance = G[node][neighbour]['weight']
            
            element = neighbour.split('_')[0]
            one_at = np.searchsorted(distance_categories, distance, side='left')
            if(one_at>=len(distance_categories)):
                one_at-=1
            for idx_1 in range(0,len(distance_categories)):
                if(idx_1 == one_at):
                    neighbour_list.append(1)
                else:
                    neighbour_list.append(0)
                    
            y = electroneg_dict[element]
            one_at = np.searchsorted(electroneg_cats, y, side='left')
            if(one_at>=len(electroneg_cats)):
                one_at-=1
            for idx_2 in range(0,len(electroneg_cats)):
                if(idx_2 == one_at):
                    neighbour_list.append(1)
                else:
                    neighbour_list.append(0)
            
            z = group_num[element]
            one_at = np.searchsorted(group_numcats, z, side='left')
            if(one_at>=len(group_numcats)):
                one_at-=1
            for idx_3 in range(0,len(group_numcats)):
                if(idx_3 == one_at):
                    neighbour_list.append(1)
                else:
                    neighbour_list.append(0)
        
            pp = period_num[element]
            one_at = np.searchsorted(period_num_cats, pp, side='left')
            if(one_at>=len(period_num_cats)):
                one_at-=1
            for idx_4 in range(0,len(period_num_cats)):
                if(idx_4 == one_at):
                    neighbour_list.append(1)
                else:
                    neighbour_list.append(0)
        node_list.append(np.asarray(neighbour_list)) 
        
    return G, np.asarray(node_list)

def get_shortest_distances(reduced_coords, amat):
    natom = len(reduced_coords)
    dists = np.zeros((natom, natom))
    Rij_min = np.zeros((natom, natom, 3))

    for i in range(natom):
        for j in range(i):
            rij = reduced_coords[i][0] - reduced_coords[j][0]
            d_min = np.inf
            R_min = np.zeros(3)
            for l in range(-1, 2):
                for m in range(-1, 2):
                    for n in range(-1, 2):
                        r = rij + np.array([l, m, n])
                        R = np.matmul(amat, r)
                        d = length(R)
                        if d < d_min:
                            d_min = d
                            R_min = R
            dists[i, j] = d_min
            dists[j, i] = dists[i, j]
            Rij_min[i, j] = R_min
            Rij_min[j, i] = -Rij_min[i, j]
    return dists, Rij_min

count=0
latt_vec = []
entries = []
cart_cords = []
n_elms = []
formation_enth = []
mag_mom = []
big_mats = []
graphs = []

for entry in tqdm_notebook(result):
    link = (str(entry)+'/POSCAR.bands')
    r = http.request('GET', link)
    if (r.status==200):
        formation_enth.append(entry.energy_atom)
        mag_mom.append(entry.spin_atom)
        entries.append(entry)
        poscar_count = 0
        poscar = r.data.decode('utf-8')
        #make lists for within loop to be overwritten and appended at end of each iteration
        #over the aflow batch
        _compound = []
        _latt_vec = []
        #_type_o_cords = []
        _cart_cords = []
        _n_elms = []
        compound_yolk = []
        
        for line in poscar.splitlines():
            if(poscar_count==0):
                _compound = entry.species                
            elif(poscar_count==1):
                lattice_constant = float(line)
            elif(poscar_count in [2,3,4]):
                _latt_vec.append([float(i) for i in line.split()])
            elif(poscar_count==5):
                _n_elms.append([int(i) for i in line.split()])
                total_no_elms = np.sum(_n_elms)
                count=0
                for element in _compound:
                    compound_yolk.append([element]*_n_elms[0][count])
                    count+=1
                compound_yolk = list(itertools.chain.from_iterable(compound_yolk))
            elif(poscar_count==6):
                #_type_o_cords.append(line)
                add_element_count=0
                
            elif((poscar_count>=7) and (poscar_count<7+total_no_elms)):
                lil_r = np.array([float(i)*lattice_constant for i in line.split()])
                _cart_cords.append(((np.matmul(np.array(_latt_vec), lil_r)), compound_yolk[add_element_count]))
                add_element_count+=1
            poscar_count +=1
        
        A = np.transpose(_latt_vec)
        B = inv(A)
                    
        crystal_red = [[np.matmul(B, R), symbol] for (R, symbol) in _cart_cords]
        crystal_dist, crystal_Rij = get_shortest_distances(crystal_red, A)

        G, maaa = get_crytal_graph(crystal_red, crystal_dist)
        big_mats.append(maaa)
        graphs.append(G)
        
        latt_vec.append(np.array(_latt_vec))
        n_elms.append(_n_elms[0])
        cart_cords.append(_cart_cords)


#The material and graph should be slightly similar
print(nx.draw_spring(graphs[139]))
print(entries[139])
