# -*- coding: utf-8 -*-

import itertools
import numpy as np
import multiprocessing as mp
import os
import configparser

default_parameters = {
    "directory": "test",
    "name": "Fake",
    "return": ["errF"],
    "Ns": "40",
    "No": "40",
    "Na": "1",
    "Nmb": "40",
    "Nmq": "40",
    "Nma": "40",
    "Nm": "40",
    "Nq": "1",
    "sigmaB": "1",
    "sigmaR": "1",
    "sigmaQ": "1",
    "algo": "IEnKS",
    "str_Q": "I-1",
    "str_H": "id",
    "str_M": "L95",
    "str_S": "0",
    "str_LR": "none",
    "str_lin": "bundle",
    "str_stateUp": "lin",
    "C": "0",
    "F": "0",
    "S": "1",
    "jmax": "10",
    "Ptol": "-3",
    "Ppas": "-4",
    "inflation": "1",
    "L": "1",
    "K": "1",
    "G": "0",
    "Ndt": "1",
    "Nspin": "5000",
    "burn": "0.05",
    "id": "1",
    "overwrite": "n"
}
cles_def = default_parameters.keys()
valeurs_def = default_parameters.values()


class Config:
    """
    A class of configuration file that is basically a dict with additional
    functions
    """

    def __init__(self, cles=cles_def, valeurs=valeurs_def):
        self.cles = list(cles)  # the list of keys
        self.valeurs = list(valeurs)  # the list of values
        self.results = None  # will contains the exp results

    def add(self, cle, valeur):
        self.cles.append(cle)
        self.valeurs.append(valeur)

    def croise(self):
        """
        Requires: self.valeurs to be a list of lists
        returns a list of Config objects by taking
        the Cartesian product of self.valeurs elements
        """
        lcfg = []
        for valeur in itertools.product(*self.valeurs):
            p = Config(self.cles, list(valeur))
            lcfg.append(p)
        return lcfg

    def shape(self):
        """
        returns the Config shape that is a tuple with the length of all
        elements in self.values
        """
        sh = ()
        for l in self.valeurs:
            sh = sh + (len(l),)
        return sh

    def tostring(self):
        """
        converts a Config object into a string
        """
        s = ""
        for cle in self.cles:
            val = str(self[cle])
            val = val.replace("/", ":")
            s += "_"+cle+"="+val
        return s

    def __getitem__(self, cle):
        """
        to get an element as a dictionary
        """
        return self.valeurs[self.cles.index(cle)]

    def __setitem__(self, cle, valeur):
        """
        to set an element as a dictionary
        """
        if cle in self.cles:
            indice = self.cles.index(cle)
            self.valeurs[indice] = valeur
        else:
            self.cles.append(cle)
            self.valeurs.append(valeur)

    def save(self, lcfg):
        """
        gathers the results of the experiments in lcfg
        and saves them into a single multidimensional numpy array
        """
        rets = self["return"][0]
        self["return"] = rets
        shape = self.shape()
        self.results = np.empty(shape)
        i, n = 0, len(rets)
        for tup in proditer(shape[:-1]):
            cfg = lcfg[i]
            for i_ret in range(n):
                self.results[tup+(i_ret,)] = np.load(cfg["directory"]+"/data/"
                                                     + cfg["name"]+"_return="
                                                     + rets[i_ret]+".npy")
            i += 1
        np.save(cfg["directory"]+"/results.npy", self.results)

    def load(self, loc_dir, dist_dir=""):
        """
        download distant results
        """
        self["return"] = self["return"][0]
        if dist_dir != "":
            os.system("scp "+dist_dir+"/results.npy "+loc_dir)
        self.results = np.ma.masked_array(np.load(loc_dir+"/results.npy"))

    def f(self, f, s):
        """
        mapping a function f on all the results in s
        """
        i = self.valeurs[-1].index(s)
        self.results[Ellipsis, i] = f(self.results[Ellipsis, i])

    def get_results(self, s):
        """
        returns a result by its name
        """
        if type(s) == str:
            s = self.valeurs[-1].index(s)
        return self.results[Ellipsis, s]

    def set_results(self, s, val):
        """
        sets a result by its name
        """
        if type(s) == str:
            s = self.valeurs[-1].index(s)
        self.results[Ellipsis, s] = val

    def opt(self, x="inflation", y='errF'):
        """
        optimizes result y over parameter x
        the parameter key is deleted and the argmin values become results
        """
        i_x = self.cles.index(x)
        i_y = self.valeurs[-1].index(y)
        i_opts = np.argmin(self.results[Ellipsis, i_y], axis=i_x)
        del self.cles[i_x]
        v_x = self.valeurs[i_x]
        del self.valeurs[i_x]
        self.valeurs[-1].append(x)
        results_ = np.ma.empty(self.shape())
        for tup in proditer(results_.shape):
            i_opt = i_opts[tup[:-1]]
            if tup[-1] == results_.shape[-1]-1:
                results_[tup] = v_x[i_opt]
            else:
                results_[tup] = self.results[
                        tup[:i_x]+(i_opt,)+tup[i_x:]
                        ]
        self.results = results_


def format_float(x, n):
    """
    formats a masked array or list of floats
    """
    if isinstance(x, np.ma.core.MaskedArray):
        x = x[~ x.mask]
    if isinstance(x, list) or isinstance(x, np.ndarray):
        if len(x) == 0:
            return []
        else:
            return [format_float(x[0], n)] + format_float(x[1:], n)
    elif isinstance(x, int):
        return str(x)
    elif isinstance(x, str):
        return x
    else:
        return "{:.{n}f}".format(x, n=n)


def format_int(x):
    """
    formats a list or a np array of integers
    """
    if isinstance(x, list) or isinstance(x, np.ndarray):
        if len(x) == 0:
            return []
        else:
            return [format_int(x[0])] + format_int(x[1:])
    elif isinstance(x, int):
        return str(x)
    elif isinstance(x, str):
        return x
    else:
        return "{0:.0f}".format(x)


def proditer(tup):
    """
    returns the product of the elements in tup
    """
    lst = []
    for n in tup:
        lst.append(tuple(range(n)))
    return itertools.product(*lst)


def ticks(t, Nreso):
    """
    sets ticks for nice plots
    """
    if len(t) < Nreso:
        return t
    else:
        t = np.ma.masked_invalid(t)
        t = t[~t.mask]
        t = np.sort(t)
        reso = (t[-1]-t[0])/Nreso
        j, i = 0, 0
        ticks = [t[0]]
        while (i < len(t)):
            while (j < len(t)) and (t[j]-t[i] < reso):
                j += 1
            if j < len(t):
                ticks += [t[j]]
            i = j
            j += 1
        return ticks


def ticks2(t, n):
    """
    another way of setting ticks
    """
    if len(t) <= 2:
        return t
    else:
        mid = 0.5*(t[0]+t[-1])
        i = np.argmin([abs(mid-e) for e in t])
        if n == 0:
            return [t[0], t[i], t[-1]]
        else:
            l1 = ticks2(t[:i], n-1)
            l2 = ticks2(t[i:], n-1)
            return l1[:-1] + l2


def ic(t):
    mid = 0.5*(t[0]+t[-1])
    return np.argmin([abs(mid-e) for e in t])


def compute(cfg):
    """
    Launches the C++ executable on a Config object cfg
    """
    directory = cfg["directory"]
    name = cfg["name"]

    # Building the configfile with ConfigParser libr
    # Warning transforms every var in lowercase
    config = configparser.ConfigParser()
    config.optionxform = str
    config.add_section("config")
    for cle in cfg.cles:
        config.set("config", cle, str(cfg[cle]))
    nameConfig = directory+'/config/'+name+'.cfg'
    # remove an already existing config file
    try:
        os.remove(nameConfig)
    except OSError:
        pass
    with open(nameConfig, 'w') as configfile:
        config.write(configfile)

    # Launch the C++ executable with the config file as parameter
    # if the overwrite parameter is "y" or a result is missing
    exist = True
    for nom_ret in cfg["return"]:
        exist = exist and os.path.isfile(directory+"/data/"+name
                                         + "_return=" + nom_ret + ".npy")
    if (cfg["overwrite"] == "y") or (not exist):
        print('Assimilation n° '+str(cfg["id"])+" : "+name)
        os.system("../Debug/IEnKS " + nameConfig)


def lcompute(lcfg, Nprocess=mp.cpu_count()):
    """
    Launches the C++ executable in parallel over each config object in the list
    lcfg with pool.map
    """
    pool = mp.Pool(processes=Nprocess)
    pool.map(compute, lcfg)
