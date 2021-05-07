#  Created on: 4 juil. 2018
#      Author: anthony
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import sys
from Config import Config, lcompute, format_float, ticks, proditer
import os

# Experiment parameters #######################################################
# meta_parameters is a dict whose values are list of parameters
meta_parameters = {
    "str_M": ["L95"],
    "algo": ["IEnKS_Q M=I", "IEnKS_bundle"],
    "Na": ["1000"],
    "sigmaQ": ["0.001**0.5", "0.1", "0.1**0.5"],
    "inflation": ["1", "1.02", "1.05", "1.1", "1.15",
                  "1.2", "1.25", "1.3", "1.4", "1.5",
                  "1.75", "2", "2.5", "3", "4"],
    "L": [str(n) for n in range(1, 4)],
    "S": ["L"],
    "G": ["0"],
    "Ndt": ["1"],
    "Nmb": ["19"],
    "Nmq": ["40"],
    "str_S": ["X"],
    "str_lin": ["bundle"],
    "str_stateUp": ["lin"],
    "return": [["errF", "errS", "n_mod"]]
}  # the last key should be "return"
###############################################################################

# Initializations #############################################################

# An experiment directory with the script name is made
directory = sys.argv[0][0:-3]
if not os.path.isdir(directory):
    os.makedirs(directory)
if not os.path.isdir(directory+"/data"):
    os.makedirs(directory+"/data")
if not os.path.isdir(directory+"/config"):
    os.makedirs(directory+"/config")

# Building the config files that will be given to the C++ executable
meta_cles = meta_parameters.keys()
meta_valeurs = meta_parameters.values()
meta_cfg = Config(meta_cles, meta_valeurs)  # homemade Config class
if len(sys.argv) == 1 or sys.argv[1] != "plot":
    meta_lcfg = meta_cfg.croise()  # creates the list of all config files
                                   # having one element of each parameters
                                   # value list
    # Adding default parameters and special parameters
    lcfg, icfg, ncfg = [], 0, len(meta_lcfg)
    for cfg_ in meta_lcfg:
        cfg = Config()
        for cle in meta_cles:
            cfg[cle] = cfg_[cle]
        cfg["id"] = str(icfg)+"/"+str(ncfg)
        icfg += 1
        if len(sys.argv) > 1 and sys.argv[1] == "overwrite":
            cfg["overwrite"] = "y"

        # Parsing special dependent parameters
        Ndt = eval(cfg["Ndt"])
        L = eval(cfg["L"])
        S = eval(cfg["S"])
        sigmaQ = eval(cfg["sigmaQ"])
        cfg["S"] = str(S)
        cfg["sigmaQ"] = str(sigmaQ)
        if cfg["algo"] == "IEnKS_bundle":
            cfg["str_S"] = "0"
        if cfg["algo"] == "IEnKF_Q_bundle":
            cfg["str_S"] = "0"
            cfg["L"] = "1"
            cfg["S"] = "1"
        if cfg["algo"] == "IEnKS_Q M=I":
            cfg["Nmb"] = "40"
            cfg["str_S"] = "svd-I_19"
            cfg["algo"] = "IEnKS_Q"
        if cfg["algo"] == "IEnKS_Q M=0":
            cfg["Nmb"] = "40"
            cfg["str_S"] = "svd-0_19"
            cfg["algo"] = "IEnKS_Q"
        if cfg["algo"] == "IEnKS_Q M=M":
            cfg["Nmb"] = "40"
            cfg["str_S"] = "full"
            cfg["algo"] = "IEnKS_Q"
        # Editing the config name
        label = ""
        for cle in meta_cles[:-1]:
            label += "_"+cle+"="+format_float(cfg[cle], 2)
        label = label[1:]
        cfg["directory"] = directory
        cfg["name"] = label.replace("/", "÷")
        # Appending the current config to the list of configs
        lcfg.append(cfg)

    # Parallel launch of the executable with the configs list
    lcompute(lcfg)
    # Save the results
    meta_cfg.save(lcfg)
###############################################################################

# Post-Treatment ##############################################################
elif len(sys.argv)>1 and sys.argv[1]=="plot":
    params = {'legend.fontsize': '10',
              'axes.labelsize': '10',
              'axes.titlesize': '10',
              'xtick.labelsize': '7',
              'ytick.labelsize': '7',
              'lines.linewidth': '2',
              'figure.subplot.top': '0.99',
              'figure.subplot.bottom': '0.138',
              'figure.subplot.left': '0.058',
              'figure.subplot.right': '0.995',
              'figure.subplot.hspace': '0.035',
              'figure.subplot.wspace': '0.135',
              'figure.figsize': '11, 8',
              'axes.grid': 'True',
              'text.usetex': 'False'
              }
    plt.rcParams.update(params)

    # plot structure
    abscisse = "L"
    ordonnee = "return"
    largeur = "sigmaQ"
    legende = "algo"

    # distant data loading
    if len(sys.argv) > 2 and sys.argv[2] == "kish":
        meta_cfg.load(directory, "filliona@kish:/libre/filliona/IEnKS/Python/" + directory)
    else:
        meta_cfg.load(directory)

    # some transformations of the results
    r = meta_cfg.get_results("errF")
    meta_cfg.set_results("errF", np.ma.log(r))
    r = meta_cfg.get_results("errS")
    meta_cfg.set_results("errS", np.ma.log(r))
    r = meta_cfg.get_results("n_mod")
    r = np.ma.masked_where((r > 2000)+(r <= 0), r)
    meta_cfg.set_results("n_mod", r)

    # optimal inflation
    meta_cfg.opt()

    i_x, n_x = meta_cfg.cles.index(abscisse), len(meta_cfg[abscisse])
    i_y, n_y = meta_cfg.cles.index(ordonnee), len(meta_cfg[ordonnee])
    i_L, n_L = meta_cfg.cles.index(largeur), len(meta_cfg[largeur])
    i_lgd, n_lgd, = meta_cfg.cles.index(legende), len(meta_cfg[legende])

    # plot
    del meta_cfg.valeurs[-1][meta_cfg.valeurs[-1].index("inflation")]
    n_y -= 1
    shape = meta_cfg.shape()
    fig, axs = plt.subplots(n_y, n_L, squeeze=False)
    
    colors = {}
    colors["IEnKS_Q M=M"] = 'C0'
    colors["IEnKS_Q M=I"] = 'C4'
    colors["IEnKS_Q M=0"] = 'C5'
    colors["IEnKF_Q_bundle"] = 'C3'
    colors["IEnKF_Q_transform"] = 'C3'
    colors["IEnKS_bundle"] = 'C2'

    markers = {}
    markers["IEnKS_Q M=M"] = 's'
    markers["IEnKS_Q M=I"] = '>'
    markers["IEnKS_Q M=0"] = '<'
    markers["IEnKF_Q_bundle"] = '^'
    markers["IEnKF_Q_transform"] = '^'
    markers["IEnKS_bundle"] = 'o'

    fancy_algo = {}
    fancy_algo["IEnKS_Q M=M"] = "IEnKS-Q $\widetilde{\mathcal{M}}=\mathcal{M}$"
    fancy_algo["IEnKS_Q M=I"] = "IEnKS-Q $\widetilde{\mathcal{M}}=\mathbf{I}$"
    fancy_algo["IEnKS_Q M=0"] = "IEnKS-Q $\widetilde{\mathcal{M}}=\mathbf{0}$"
    fancy_algo["IEnKF_Q_bundle"] = "IEnKF-Q (Sakov et al., 2018)"
    fancy_algo["IEnKF_Q_transform"] = "IEnKF-Q transform (Sakov et al., 2018)"
    fancy_algo["IEnKS_bundle"] = "IEnKS"

    fancy_return = {}
    fancy_return["errF"] = "filtering RMSE"
    fancy_return["errS"] = "smoothing RMSE"
    fancy_return["n_mod"] = "number of model evaluations"

    for tup in proditer(shape[:i_x]+shape[i_x+1:]):
        tup = tup[:i_x]+(slice(0, None),)+tup[i_x:]
        ax = axs[tup[i_y], tup[i_L]]
        im = ax.plot(
                 meta_cfg.results[tup],
                 label=fancy_algo[meta_cfg[legende][tup[i_lgd]]],
                 marker=markers[meta_cfg[legende][tup[i_lgd]]],
                 color=colors[meta_cfg[legende][tup[i_lgd]]],
                 )
        ax.set_xticks(range(len(meta_cfg.results[tup])))
        ax.set_xticklabels([])
        if meta_cfg[ordonnee][tup[i_y]] == "errS" or meta_cfg[ordonnee][tup[i_y]]=="errF":
            tck = ticks(meta_cfg.results[tup[:i_lgd]+(slice(0, None),)
                                         + tup[i_lgd + 1:]].flatten(), 10)
            ax.set_yticks(tck)
            ax.set_yticklabels(format_float(np.exp(tck), 3))
        if tup[i_y] == n_y-1:
            ax.set_xticklabels(meta_cfg[abscisse], rotation=0)
            ax.set_xlabel(abscisse+"\n"+"q = "
                          + ["0.001", "0.01", "0.1"][tup[i_L]])
        if tup[i_L] == 0:
            ax.set_ylabel(fancy_return[meta_cfg[ordonnee][tup[i_y]]])
    lbl = ""
    for icle in range(len(shape)-1):
        if shape[icle] == 1 and icle != i_x:
            lbl += " "+meta_cfg.cles[icle]+"="+str(meta_cfg.valeurs[icle][0])
    axs[-1, 1].legend(bbox_to_anchor=(0.4, -0.5), loc=8, ncol=5)
    fig.savefig(directory+"/"+directory+".eps")
    plt.show()
