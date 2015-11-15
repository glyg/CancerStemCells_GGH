# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import colors

# # matplotlib settings

fontdict = {'fontsize':12}
sim_cm = plt.get_cmap('viridis')
csc_cm = plt.get_cmap('afmhot')
npc_cm = plt.get_cmap('bone')


def show_fields(tumors, frame_num, f_size_i):
    '''Snapshots of the simulation for various data types
    '''
    fig, axes = plt.subplots(len(tumors), len(tumors[0].data_fields))
    fig.set_size_inches(f_size_i)
    for ax_line, tumor in zip(axes, tumors.values()):
        diff_adh = tumor.sim_dict['energies']['CancerStemCell-NonCancerous']
        for ax, (name, field) in zip(ax_line, tumor.data_fields.items()):
            ax.imshow(field[frame_num], cmap='hot', origin='lower', interpolation='nearest')
            ax.set_title(name, fontdict=fontdict)
    return fig, axes

def show_type(tumors, frame_num, f_size_i):
    ''' Image of the simulation with the cells colored by type
    '''
    
    fig, axes = plt.subplots(2, len(tumors)//2)
    fig.set_size_inches(f_size_i)
    for ax, tumor in zip(axes.ravel(), tumors.values()):
        diff_adh = tumor.sim_dict['energies']['CancerStemCell-NonCancerous']
        ax.imshow(tumor.data_fields['CellType'][frame_num], cmap='hot',
                  origin='lower', interpolation='nearest')
        ax.set_title(diff_adh, fontdict=fontdict)
        ax.set_xticks([])
        ax.set_yticks([])
    return fig, axes

def show_time_components(tumors, f_size_i):
    '''
    Plot the variation against time of
    a dicttionnary of tumor data to be ploted is
    specified bellow.
    '''

    columns = [u'area', u'ncells', u'pis']
    ylims = {u'area': (0, 45),
             u'ncells': (0, 1000),
             u'pis': (0, 1.1),}
    ylabels = {u'area': u'Area (px)',
               u'ncells': u'Number of cells',
               u'pis': u'Clustering',}

    fig, axes = plt.subplots(1, len(columns))
    fig.set_size_inches(f_size_i)
    font = fontdict
    csc_norm  = colors.Normalize(vmin=0, vmax=len(tumors))
    csc_map = plt.cm.ScalarMappable(norm=csc_norm, cmap=csc_cm)
    npc_norm  = colors.Normalize(vmin=0, vmax=len(tumors))
    npc_map = plt.cm.ScalarMappable(norm=npc_norm, cmap=npc_cm)

    for n, tumor in enumerate(tumors.values()):
        csc_color = csc_map.to_rgba(n)
        npc_color = npc_map.to_rgba(n)
        for ax, col in zip(axes, columns):
            ax.plot(tumor.csc_df.index, tumor.csc_df[col],
                    color=csc_color, marker='', lw=2, alpha=0.6)
            ax.plot(tumor.npc_df.index, tumor.npc_df[col],
                    color=npc_color, marker='', lw=2, alpha=0.6)

            ax.set_ylim(ylims[col])
            ax.set_ylabel(ylabels[col], fontdict=font)
            ax.set_xlabel(u'Time (step)', fontdict=font)
    fig.set_tight_layout(True)
    return fig, axes

def plot_color_gradients(cmaps, title):
    '''
    Small utility to plot a standalone colorbar
    for colormap cmap
    '''
    nrows = len(cmaps)
    fig, axes = plt.subplots(nrows=nrows)
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    axes[0].set_title(title, fontsize=14)
    for ax, (cmap, ticks, cmap_label) in zip(axes, cmaps.values()):
        gradient = np.vstack((ticks, ticks))
        ax.imshow(gradient, aspect='auto', 
                  cmap=plt.get_cmap(cmap.name),
                  interpolation='nearest')
        ax.set_yticks([])
        ax.set_ylabel(cmap_label)
        ax.set_gid('off')
