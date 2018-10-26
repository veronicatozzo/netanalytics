import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_degree_distribution_histogram(degrees, ax=None, title="",
                                       xlabel="Degrees",
                                       ylabel="Probability",
                                       to_file=None):
    """
    degrees: dict
        Dictionary (key, value)=(degree, occurence)

    ax: fig.ax, default=None
        The axis to use to plot the figure, if None a new one is generated.

    title:

    xlabel:

    ylabel:

    to_file: string, default=None
        The path to the file in which to save the plot. If None the plot will
        not be saved, just shown.

    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(15,5))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    sns.distplot(degrees, bins=degrees.shape[0],#range(int(np.min(keys)), int(np.max(keys)), 1),
                 ax=ax, kde=True, norm_hist=False,
                 hist=True)
    if to_file is not None:
        plt.savefig(to_file, transparency=True, bbox_inches='tight', dpi=200)

    plt.show()
