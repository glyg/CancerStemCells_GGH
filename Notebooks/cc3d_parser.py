# -*- coding: utf-8 -*-

import os
from xml.etree import ElementTree
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import vtk
import json
from vtk.util.numpy_support import vtk_to_numpy

from skimage.filters import rank
from skimage.morphology import disk

from skimage import filters
from skimage.future import graph

from skimage import measure
from scipy import special

def parse_cc3d(data_dir, simname):

    sim_dict = {}
    sim_str = []

    sim_str.append('## File used: ')
    sim_str.append('Using: `{}`'.format(data_dir))


    xml_file = os.path.join(data_dir+'/Simulation/{}.xml'.format(simname))
    py_file = os.path.join(data_dir+'/Simulation/{}Steppables.py'.format(simname))

    xml_tree = ElementTree.ElementTree(file=xml_file)
    xml_root = xml_tree.getroot()

    for child in xml_root:

        if child.tag == 'Potts':
            steps = int(child.find('Steps').text)
            dim = {k: int(v) for k, v
                   in child.find('Dimensions').items()}
            sim_dict['dim'] = dim
            sim_dict['steps'] = steps
        elif child.get('Name') == 'Contact':
            energies = child.findall('Energy')
            energies = {'{}-{}'.format(energy.attrib['Type1'],
                                       energy.attrib['Type2']):
                        float(energy.text) for energy in energies}

    sim_dict['energies'] = energies

    sim_str.append('<hr/>')
    sim_str.append('## Energies: ')
    for k, v in energies.items():
        sim_str.append('{}: {}'.format(k, v))


    pysettings = []

    with open(py_file) as pf:
        start_parse = False
        for line in pf.readlines():
            if line.startswith('# <parameter settings>'):
                start_parse = True
            elif (start_parse
                  and not line.startswith('#')
                  and len(line) > 3):
                pysettings.append(line[:-1])
            elif line.startswith('# </parameter settings>'):
                break


    sim_dict['pysettings'] = pysettings

    sim_str.append('<hr/>\n')
    sim_str.append('## python set variables:')
    for p in pysettings:
        sim_str.append(p)

    vtk_dir = os.path.join(data_dir,'LatticeData')

    vtk_files = [os.path.join(vtk_dir, f)
                 for f in os.listdir(vtk_dir)
                 if f.endswith('.vtk')]
    vtk_files.sort()
    sim_dict['vtk_files'] = vtk_files

    sim_str.append('<hr/> \n')
    sim_str.append('##  VTK files\n')
    sim_str.append('starts: {}\n stops: {}'.format(vtk_files[0], vtk_files[-1]))
    return sim_dict, sim_str

def parse_vtk(vtk_file, sim_dict, data_fields):

    dim = sim_dict['dim']
    reader = vtk.vtkStructuredPointsReader()
    #reader = vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.Update()
    field_data = reader.GetOutput()
    out_data = {}
    for data_field in data_fields:
        out_data[data_field] = vtk_to_numpy(
            field_data.GetPointData().GetArray(data_field)).reshape((dim['x'],
                                                                     dim['y']))
    reader.CloseVTKFile()
    return out_data


def parse_all_vtks(sim_dict, field_names):
    vtk_files = sim_dict['vtk_files']
    num_steps = len(vtk_files)
    dim = sim_dict['dim']
    data_fields = {field: np.zeros((num_steps, dim['x'], dim['y']))
                   for field in field_names}
    step_values = []
    for i, vtk_file in enumerate(vtk_files):
        step_values.append(int(vtk_file.split('.')[1][-3:]))
        frame_data = parse_vtk(vtk_file, sim_dict, field_names)
        for key, val in frame_data.items():
            data_fields[key][i] = val
    return data_fields, np.array(step_values)


class Tumor:
    '''
    Container class
    '''
    field_names = ['CellType', 'CellId', 'CellAge']

    def __init__(self, sim_dict):
        self.sim_dict = sim_dict
        self.data_fields, self.step_values = parse_all_vtks(sim_dict,
                                                            self.field_names)
        self.get_idxs()
        self.get_cell_df()
        self.get_edge_df()
        self.get_entropy()
        self.cell_averages()

    def get_idxs(self):

        cell_types, cell_ids = (self.data_fields['CellType'],
                                self.data_fields['CellId'])
        e_idx = []
        v_idx = []
        for mcs, cell_type, cell_id in zip(self.step_values,
                                           cell_types, cell_ids):
            rag = graph.rag_mean_color(cell_type, cell_id)
            direct = [(mcs, s, t) for s, t in rag.edges()]
            fliped = [(mcs, t, s) for s, t in rag.edges()]
            e_idx.extend(direct + fliped)
            v_idx.extend([(mcs, cell_id) for cell_id in rag.nodes()])

        self.e_idx = pd.MultiIndex.from_tuples(e_idx,
                                               names=['t', 'srce', 'trgt'])
        self.v_idx = pd.MultiIndex.from_tuples(v_idx,
                                               names=['t', 'cell'])

    def get_cell_df(self):

        cell_data = ['type', 'age', 'area', 'cx', 'cy']
        self.cell_df = pd.DataFrame(index=self.v_idx, columns=cell_data)

        for i, mcs in enumerate(self.step_values):
            cell_type = self.data_fields['CellType'][i]
            cell_id = self.data_fields['CellId'][i]
            cell_age = self.data_fields['CellAge'][i]
            properties = measure.regionprops(cell_id.astype(np.int))
            for p in properties:
                lbl = p['label']
                self.cell_df.loc[(mcs, lbl), 'area'] = p['area']
                cx, cy = p['centroid']
                self.cell_df.loc[(mcs, lbl), 'cx'] = cx
                self.cell_df.loc[(mcs, lbl), 'cy'] = cy

                self.cell_df.loc[(mcs, lbl), 'age'] = cell_age[int(cx),
                                                               int(cy)]
                self.cell_df.loc[(mcs, lbl), 'type'] = cell_type[int(cx),
                                                                 int(cy)]

    def get_edge_df(self):

        edge_columns = ['srce_type', 'trgt_type']
        self.edge_df = pd.DataFrame(index=self.e_idx, columns=edge_columns)

        t_idx = self.e_idx.get_level_values(level='t')
        srce_idx = self.e_idx.get_level_values(level='srce')
        trgt_idx = self.e_idx.get_level_values(level='trgt')

        self.srce_idx = pd.MultiIndex.from_arrays([t_idx, srce_idx],
                                                  names=['t', 'cell'])
        self.trgt_idx = pd.MultiIndex.from_arrays([t_idx, trgt_idx],
                                                  names=['t', 'cell'])

        self.edge_df['srce_type'] = self.cell_df['type'].loc[self.srce_idx].values
        self.edge_df['trgt_type'] = self.cell_df['type'].loc[self.trgt_idx].values

        self.edge_df['type_diff'] = self.edge_df.trgt_type == self.edge_df.srce_type

        self.edge_df = self.edge_df.drop(0, level='srce')
        self.edge_df = self.edge_df.drop(0, level='trgt')


    def get_entropy(self):

        pis = self.edge_df.type_diff.groupby(
            level=['t', 'srce']).apply(
            lambda df: df.sum()/float(df.size))

        pis.index.names = ['t', 'cell']
        self.cell_df['pis'] = pis
        entropy = - pis * np.log2(pis)
        self.cell_df['entropy'] = entropy

    def cell_averages(self):


        cell_df = self.cell_df
        self.is_csc = cell_df.loc[cell_df.type == 1].index
        self.is_npc = cell_df.loc[cell_df.type == 2].index

        self.n_csc = (cell_df.type == 1).groupby(level='t').sum()
        self.n_npc = (cell_df.type == 2).groupby(level='t').sum()

        self.area_csc = cell_df.area.loc[self.is_csc].groupby(
            level='t').apply(np.mean)
        self.area_npc = cell_df.area.loc[self.is_npc].groupby(
            level='t').apply(np.mean)

        self.entropy_csc = cell_df.entropy.loc[self.is_csc].groupby(
            level='t').apply(np.mean)
        self.entropy_npc = cell_df.entropy.loc[self.is_npc].groupby(
            level='t').apply(np.mean)

        self.pis_csc = cell_df.pis.loc[self.is_csc].groupby(
            level='t').apply(np.mean)
        self.pis_npc = cell_df.pis.loc[self.is_npc].groupby(
            level='t').apply(np.mean)


def parse_data_dir(DATA_ROOT):

    data_dirs = os.listdir(DATA_ROOT)
    data_dirs.sort()
    data_dirs = [os.path.join(DATA_ROOT, d) for d in data_dirs]
    data_dirs = [d for d in data_dirs if os.path.isdir(d)]

    sim_strs = {}
    sim_dicts = {}
    for n, ddir in enumerate(data_dirs):
        sim_dicts[n], sim_strs[n] = parse_cc3d(ddir, simname='Sim2')
        exec(' '.join(sim_dicts[n]['pysettings']))
        sim_dicts[n]['py_params'] = params.copy()
        with file(ddir+'.json', 'w') as fp:
            json.dump(sim_dicts[n], fp)

    tumors = {}
    for k, sim_dict in sim_dicts.items():
        tumors[k] = Tumor(sim_dict)

    return tumors
