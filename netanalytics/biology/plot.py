import matplotlib

import numpy as np
import matplotlib.pyplot as plt

def plot_enrichment_results(bars, file, groups, labels, ylabel, title):
    N = len(groups)
    if len(bars) != N:
        raise ValueError("You passed a different number of bars to plot "
                         "then groups")
    if len(bars[0]) != len(labels):
        raise ValueError("You passed %d bars but you passed %d labels"%
                         (len(bars[0]), len(labels)))
        
    fig, ax = plt.subplots(figsize=(15,6))
    
    ind = np.arange(N)    # the x locations for the groups
    width = 0.9/len(bars[0])         # the width of the bars
    
    colors = ['darksalmon', 'darkturquoise', 'palevioletred', 'darkslateblue', 'plum', 'brown']
    bars = np.array(bars)
    ps = []
    for c in range(bars.shape[1]):
        ps.append(ax.bar(ind+c*width, bars[:,c], width, color=colors[c]))

    ax.set_title(title)
    ax.set_xticks(ind + width)
    ax.set_ylabel(ylabel)
    ax.set_xticklabels((groups), rotation=90)

    ax.legend(ps, labels)
    ax.autoscale_view()
    if file[-3:] == 'pdf':
          plt.savefig(file.split('.')[-2]+".png", transparent=True, dpi=200, bbox_inches='tight')
    elif file[-3:] == 'png':
          plt.savefig(file.split('.')[-2]+".pdf", transparent=True, dpi=200, bbox_inches='tight')
    plt.savefig(file, transparent=True, dpi=200, bbox_inches='tight')
    plt.show()
    
