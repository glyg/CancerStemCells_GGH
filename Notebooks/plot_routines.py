# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import colors



def show_fields(tumors, frame_num, f_size_i):
    fig, axes = plt.subplots(len(tumors), len(tumors[0].data_fields))
    fig.set_size_inches(f_size_i)
    for ax_line, tumor in zip(axes, tumors.values()):
        diff_adh = tumor.sim_dict['energies']['CancerStemCell-NonCancerous']
        for ax, (name, field) in zip(ax_line, tumor.data_fields.items()):
            ax.imshow(field[frame_num], cmap='hot', origin='lower', interpolation='nearest')
            ax.set_title(name, fontdict={'fontsize':12})
    return fig, axes

def show_type(tumors, frame_num, f_size_i):
    fig, axes = plt.subplots(2, len(tumors)//2)
    fig.set_size_inches(f_size_i)
    for ax, tumor in zip(axes.ravel(), tumors.values()):
        diff_adh = tumor.sim_dict['energies']['CancerStemCell-NonCancerous']
        ax.imshow(tumor.data_fields['CellType'][frame_num], cmap='hot',
                  origin='lower', interpolation='nearest')
        ax.set_title(diff_adh, fontdict={'fontsize':12})
        ax.set_xticks([])
        ax.set_yticks([])
    return fig, axes

def show_time_components(tumors, f_size_i):

    columns = [u'area', u'ncells', u'pis']
    ylims = {u'area': (0, 45),
             u'ncells': (0, 1000),
             u'pis': (0, 1.1),}
    ylabels = {u'area': u'Area (px)',
               u'ncells': u'Number of cells',
               u'pis': u'Clustering',}

    csc_cm = plt.get_cmap('afmhot')
    csc_norm  = colors.Normalize(vmin=0, vmax=len(tumors))
    csc_map = plt.cm.ScalarMappable(norm=csc_norm, cmap=csc_cm)

    npc_cm = plt.get_cmap('bone')
    npc_norm  = colors.Normalize(vmin=0, vmax=len(tumors))
    npc_map = plt.cm.ScalarMappable(norm=npc_norm, cmap=npc_cm)

    fig, axes = plt.subplots(1, len(columns))
    fig.set_size_inches(f_size_i)
    font = {'size':12}

    for n, tumor in enumerate(tumors.values()):
        csc_color = csc_map.to_rgba(n)
        npc_color = npc_map.to_rgba(n)
        for ax, col in zip(axes, columns):
            ax.plot(tumor.csc_df.index, tumor.csc_df[col],
                    color=csc_color, marker='o', lw=2, alpha=0.8)
            ax.plot(tumor.npc_df.index, tumor.npc_df[col],
                    color=npc_color, marker='o', lw=2, alpha=0.8)

            ax.set_ylim(ylims[col])
            ax.set_ylabel(ylabels[col], fontdict=font)
            ax.set_xlabel(u'Time (step)', fontdict=font)
    fig.set_tight_layout(True)
    return fig, axes
