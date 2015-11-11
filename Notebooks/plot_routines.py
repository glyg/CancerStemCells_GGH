# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt



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
    fig, axes = plt.subplots(1, len(tumors))
    fig.set_size_inches(f_size_i)
    for ax, tumor in zip(axes, tumors.values()):
        diff_adh = tumor.sim_dict['energies']['CancerStemCell-NonCancerous']
        ax.imshow(tumor.data_fields['CellType'][frame_num], cmap='hot', origin='lower', interpolation='nearest')
        ax.set_title(diff_adh, fontdict={'fontsize':12})
        ax.set_xticks([])
        ax.set_yticks([])
    return fig, axes

def show_time_components(tumors, f_size_i):

    fig, axes = plt.subplots(1, 3)
    fig.set_size_inches(f_size_i)
    cm = plt.get_cmap('afmhot')
    from matplotlib import colors

    cNorm  = colors.Normalize(vmin=0, vmax=len(tumors))
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=cm)

    for n, tumor in enumerate(tumors.values()):

        color = scalarMap.to_rgba(n)
        axes[0].plot(tumor.n_csc.index, tumor.n_csc,
                     color=color, marker='o', lw=2, alpha=0.8)
        axes[0].plot(tumor.n_npc.index, tumor.n_npc,
                     color=color, marker='s', lw=2, alpha=0.8)

        axes[1].plot(tumor.area_csc.index, tumor.area_csc,
                     color=color, marker='o', lw=2, alpha=0.8)
        axes[1].plot(tumor.area_npc.index, tumor.area_npc,
                     color=color, marker='s', lw=2, alpha=0.8)

        axes[2].plot(tumor.pis_csc.index, tumor.pis_csc,
                     color=color, marker='o', lw=2, alpha=0.8)
        axes[2].plot(tumor.pis_npc.index, tumor.pis_npc,
                     color=color, marker='s', lw=2, alpha=0.8)
        axes[2].set_title('Clustering')
    axes[0].set_title('Populations')
    axes[1].set_title('Areas')
    axes[1].set_ylim(0, 45)
    axes[2].set_ylim(0, 1.1)

    return fig, axes
