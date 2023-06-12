from sage.plot.histogram import Histogram
from sage.matrix.constructor import matrix
from sage.plot.matrix_plot import matrix_plot
from sage.plot.plot import list_plot
from sage.rings.rational_field import QQ
from sage.misc.flatten import flatten
from sage.misc.functional import numerical_approx
from sage.functions.other import floor, ceil
from sage.arith.misc import falling_factorial


import numpy as np
from matplotlib.colors import ListedColormap

from rec_sequences.CFiniteSequenceRing import *
from rec_sequences.DFiniteSequenceRing import *
from ore_algebra import *

import pickle
import random
import itertools
import os
from multiprocessing import Pool
import re
import urllib.request
from urllib.error import URLError
import logging
from matplotlib import pyplot as plt


def save_order_degree_plot(recs, folder_name, ens) :
    r"""
    Creates and saves a plot of the orders/degrees of the given
    recurrences in the given folder.
    """
    min_order = 1
    order_bound = 10
    degree_bound = 10
    ord_deg_matrix = matrix(QQ, order_bound-min_order, degree_bound)
    for rec in recs :
        if rec.order() < order_bound and rec.order() >= min_order \
             and rec.degree() < degree_bound :
            ord_deg_matrix[rec.order()-min_order ,rec.degree()] += 1
    total = len(recs)
    ord_deg_matrix = ord_deg_matrix/total
    title = f"Degrees/orders of sequences with ensure={ens}"
    
    N = 256
    vals = np.ones((N, 4))
    vals[:, 0] = np.linspace(plt.cm.tab20(1)[0], 1, N)
    vals[:, 1] = np.linspace(plt.cm.tab20(1)[1], 1, N)
    vals[:, 2] = np.linspace(plt.cm.tab20(1)[2], 1, N)
    newcmp = ListedColormap(vals).reversed()
    
    plot = matrix_plot(ord_deg_matrix, yrange = (0.5,order_bound-1),
                       axes_labels=['degree','order'], colorbar=True, 
                       cmap=newcmp)
    plot.save(folder_name + f"/order_degree_ens{ens}.png")

def create_plot(recs_ensures, folder_name, file_name, func, ensures, total) :
    r"""
    Creates and saves a plot of the comparisons of certain values describing
    recurrences with different ensure values.
    """
    recs_orders = []
    for recs in recs_ensures :
        recs_orders.append([func(rec) for rec in recs])
        
    max_order = max(flatten(recs_orders))
    
    recs_orders_cum = []
    for recs_orders_ensure in recs_orders :
        recs_orders_cum.append([recs_orders_ensure.count(i) 
                                    for i in range(0, max_order)])
    
    max_order_display = 25
    colors = ["#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"]
    percentages = [numerical_approx(100*len(recs)/total, digits=3) 
                    for recs in recs_ensures]
    labels = [f"$e={ens}$, ${perc}\\%$" 
                for ens, perc in zip(ensures, percentages)]
    list_plots = [list_plot(list(zip(range(0,max_order_display),
                                     recs_orders_cum[i][:max_order_display])), 
                            plotjoined=True, color=plt.cm.tab20(i)[:-1], legend_label=labels[i],fontsize=13) 
                    for i in range(len(ensures))]
    plot = sum(list_plots)
    plot.save(folder_name + "/" + file_name + ".png", legend_font_size=13)

def save_orders_plot(recs_ensures, folder_name, ensures, total) :
    r"""
    Creates and saves a plot of the comparisons of the orders
    of the recurrences with different ensure values.
    """
    get_order = lambda rec : rec.order()
    return create_plot(recs_ensures, folder_name, "compare_ensures_orders",
                       get_order, ensures, total)
    
def save_degrees_plot(recs_ensures, folder_name, ensures, total) :
    r"""
    Creates and saves a plot of the comparisons of the orders
    of the recurrences with different ensure values.
    """
    get_degree = lambda rec : rec.degree()
    return create_plot(recs_ensures, folder_name, "compare_ensures_degrees",
                       get_degree, ensures, total)


def print_compare_scraped(seqs_guessed, names_scraped, ensures) :
    r"""
    """
    names_scraped = set(names_scraped)
    names_ens = [set(seqs_guessed_ens) 
                    for seqs_guessed_ens in seqs_guessed]
    print("Comparison between guessed sequences and scraped sequences "
          "from Wiki")
    logging.error("Comparison between guessed sequences and scraped sequences from Wiki")
    print("Given numbers represent: ")
    logging.error("Given numbers represent: ")
    print("1. Sequence found by guessing, but not by scraping")
    logging.error("1. Sequence found by guessing, but not by scraping")
    print("2. Sequence found both by guessing and scraping")
    logging.error("2. Sequence found both by guessing and scraping")
    print("3. Sequence not found by guessing, but by scraping")
    logging.error("3. Sequence not found by guessing, but by scraping")
    for ens, names in zip(ensures, names_ens) :
        not_scraped = len(names.difference(names_scraped))
        not_guessed = len(names_scraped.difference(names))
        both = len(names.intersection(names_scraped))
        text = f"   ensure={ens}: {not_scraped}, {both}, {not_guessed}"
        print(text)
        logging.error(text)

def get_files(directory) :
    r"""
    Given a local directory, return a list of all pickle files.
    
    INPUT:
    
    - ``directory`` -- a path to a local directory 
    
    OUTPUT:
    
    A list of pickle files in the given ``directory``.
    """
    files = []
    for file in os.listdir(directory):
        if file.endswith(".pickle"):
            files.append(os.path.join(directory, file))
    return files

def guess_function(data, order, degree, ensure) :
    D = len(data)-1
    max_order = min(floor((D-ensure)/2), order)
    for r in range(0, max_order+1) :
        max_degree = max(floor(degree/(r+1)-1), 0) # curve
         # linear system overdetermined
        max_degree = min(max(0, floor((D-r+1-ensure)/(r+1)-1)), max_degree)
        for d in range(0, max_degree+1) :
            N = (r+1)*(d+1) + ensure
            if N-1+r > D :
                continue
            M = matrix(QQ, N, (r+1)*(d+1))
            for n in range(0, N) :
                M.set_row(n, flatten(
                [[
                    falling_factorial(n-j+i, i)*data[n-j+i] if j <= n else 0 for i in range(r+1)] 
                    for j in range(d+1)
                ]))
            ker = M.right_kernel()
            if ker.dimension() != 0 :
                R = PolynomialRing(QQ, "n")
                Q = OreAlgebra(R, "Dn")
                return Q("Dn")
    raise ValueError("no deq found")

def guess_rec(start_index, terms, ensure, rec = False, max_order = 10, 
              c_fin = False, max_degree = 0) :
    r"""
    Tries to guess a recurrence for the given terms of a sequence ``seq``.
    
    INPUT:
    
    - ``start_index`` -- the start index of the sequence under consideration,
      i.e., ``terms[0] == seq[start_index]``.
    - ``terms`` -- a list of numbers representing the terms of the sequence.
    - ``ensure`` -- a natural number representing how many more equations than
      variables the corresponding linear system for guessing should have
    - ``rec`` -- if a recurrence operator which annihilates the given terms is 
      given, it is checked whether this annihilator is small enough to make sure
      that the given ``ensure`` level is still met
    - ``max_order`` -- the maximal order of the recurrence that is guessed
    - ``c_fin`` -- if ``True`` only an operator with constant coefficient is 
      guessed
    - ``max_degree`` -- the maximal degree of the recurrence operator that is
      guessed.
    
    OUTPUT:
    
    A pair ``(rec, flag)`` where ``rec`` is the guessed recurrence operator
    that annihilates the given terms (up to a certain index which is determined
    by the ``ensure``-level) and ``flag`` is an integer representing 
    
        - 0 ... no recurrence ``rec`` given
        - 1 ... too few terms given to verify given recurrence ``rec``
        - 2 ... recurrence ``rec`` could be verified
        - 3 ... recurrence ``rec`` is false
    """
    if c_fin :
        R = CFiniteSequenceRing(QQ)
    else :
        PR = PolynomialRing(QQ, "n")
        R = DFiniteSequenceRing(PR)
    
    # 0 ... recurrence not given
    # 1 ... too few terms to verify
    # 2 ... recurrence verified
    # 3 ... recurrence false
    # print(f"Terms: {terms[:5]}, {len(terms)}")
    recurrence_checked = 0
    if c_fin and not isinstance(rec, bool) :
        r = rec.order()
        d = rec.degree()
        # number of equations 
        N = ceil((r+1)*(d+1)*(1+ensure/10))+floor(ensure)
        # recurrence is given, test whether it works 
        # linear system has N rows, last term checked
        # is therefore at N+r-1
        if len(terms) < N+r :
            recurrence_checked = 1
        else :
            data_cut = terms[:(N+r)]
            if any(rec(data_cut)) : # not all zero
                recurrence_checked = 3
            else :
                recurrence_checked = 2
                return rec, recurrence_checked
        
    # if only very few terms are known we assume not C/D-finite
    if len(terms) <= 10 :
        return False
    try :
        if c_fin :
            return R._own_guess(terms, ensure=ensure, max_order=max_order,
                                max_degree = 0, ensure_relative=True,sparse=100), \
                   recurrence_checked
        else :
            #R = PolynomialRing(QQ, "x")
            #Q = OreAlgebra(R, "Dx")
            #return guess(terms[:200], Q, ensure=ensure, order=max_order, 
            #             degree=max_degree), 0
            return R._own_guess(terms, ensure=ensure, max_order=max_order,
                    max_degree = max_degree, sparse=100), 0
    except ValueError:
        return False, recurrence_checked

def guess_recs(seqs, number, ensure, threads, max_order, 
               c_fin = False, max_degree = 0) :
    r"""
    guess a recurrence for every such sequence and check ``ensure``
    many terms on this guessed sequence.
    
    INPUT:
    
    - ``seqs`` -- a dictionary of sequences of the form 
        "A..." : (i, [val[i], val[i+1],...], rec)
      where the recurrence ``rec`` is optional
    - ``number`` -- an integer; if positive, only sequences with identifier
      less or equal ``number`` are considered; otherwise, all sequences are
      considered
    - ``ensure`` -- a natural numbers indicating the certainty-level
      that should be used 
    - ``threads`` -- a natural number indicating how many threads should be used
      for guessing the recurrences
    - ``max_order`` -- the maximal order of the recurrence that is guessed
    - ``c_fin`` -- if ``True`` only an operator with constant coefficient is 
      guessed
    - ``max_degree`` -- the maximal degree of the recurrence operator that is
      guessed.
      
    OUTPUT:
    
    A tuple ``(guessed_seqs, rec_verified, rec_not_verifiable, rec_wrong)``
    where ``guessed`` seqs is a dictionary of the guessed sequences of the form
        "A..." : (i, [val[i], val[i+1],...], rec),
    ``rec_verified`` gives the number of sequences where the given recurrence 
    could be verified, ``rec_not_verifiable`` gives the number of sequences
    where not enough terms are given to verify the sequence and ``rec_wrong``
    gives the number of sequences where the given recurrence could be shown
    to be wrong.
    """
    guessed_seqs = dict()
    
    # prepare multithreading by putting all information in one list
    names = []
    start_indices = []
    terms = []
    recs = []
    for name, seq in seqs.items() :
        if number < 0 or int(name[1:]) <= number :
            start_index, terms_seq = seq[0], seq[1]
            names.append(name)
            start_indices.append(start_index)
            terms.append(terms_seq)
            recs.append(False if len(seq) <= 2 else seq[2])
    
    parameters = list(zip(start_indices,
                          terms,
                          itertools.repeat(ensure),  
                          recs,
                          itertools.repeat(max_order),
                          itertools.repeat(c_fin),
                          itertools.repeat(max_degree)))
    
    rec_verified, rec_not_verifiable, rec_wrong = 0, 0, 0
    with Pool(threads) as p :
        result = p.starmap(guess_rec, parameters)

        for i, rec in enumerate(result) :
            if isinstance(rec,tuple) and len(rec) > 1 :
                # have additional information on recurrence
                if rec[1] == 2 :
                    rec_verified += 1
                elif rec[1] == 1 :
                    rec_not_verifiable += 1
                elif rec[1] == 3 :
                    rec_wrong += 1
                elif rec[1] == 0 :
                    # no recurrence given to check
                    pass 
                else :
                    print("Something went wrong, got wrong flag for recurrence")
                    print(rec)
                    logging.error("Something went wrong, got wrong flag for recurrence")
                    logging.error(f"Recurrence: {rec}")
                rec = rec[0]
            if not isinstance(rec, bool) :
                # recurrence was found
                guessed_seqs[names[i]] = (start_indices[i], terms[i], rec)

    #if len(guessed_seqs) :
    #    print(list(guessed_seqs.keys()))
       
    return guessed_seqs, rec_verified, rec_not_verifiable, rec_wrong

def save_recs(name, number, ensures, threads, path_pickles, name_pickles,
              max_order, c_fin = False, max_degree = 0) :
    r"""
    Guess and save recurrences satisfied by the specified sequences.
    For every given ``ensure`` level a pickle file containing a
    dictionary with the guessed sequences/recurrences in the form
    
        "A......" : (start_index, terms, recurrence) 
        
    is saved.
    
    INPUT:
    
    - ``name`` -- name of the folder containing pickle files with all
      sequences
    - ``number`` -- an integer; if positive, only sequences with identifier
      less or equal ``number`` are considered; otherwise, all sequences are
      considered
    - ``ensures`` -- a list of natural numbers indicating the certainty-levels
      that should be used 
    - ``threads`` -- a natural number indicating how many threads should be used
      for guessing the recurrences
    - ``path_pickles`` -- the path where the recurrences should be saved
    - ``name_pickles`` -- the file names for the saved recurrences
    - ``max_order`` -- the maximal order of the recurrence that is guessed
    - ``c_fin`` -- if ``True`` only an operator with constant coefficient is 
      guessed
    - ``max_degree`` -- the maximal degree of the recurrence operator that is
      guessed.
    """
    # get all pickle files with sequences
    files = get_files(name)
    
    previously_guessed = False
    
    for ens in ensures :
        guessed_seqs = dict()
        rec_verified, rec_not_verifiable, rec_wrong = 0, 0, 0
        if ens == ensures[0] :
            # for first ensure we have have to go through all sequences
            for file in files :
                with open(file, 'rb') as handle:
                    # dictionary "A..." : (i, [val[i], val[i+1],...])
                    seqs = pickle.load(handle)
                    new_guessed = guess_recs(seqs, number, ens, threads, 
                                             max_order, c_fin, max_degree)[0]
                    guessed_seqs.update(new_guessed)
                    
                print(f"Ensure={ens}, file={file} done")
                print(f"{len(guessed_seqs)} sequences found so far")
                
                logging.error(f"Ensure={ens}, file={file} done")
                logging.error(f"{len(guessed_seqs)} sequences found so far")
                
        else :
            guessed, rec_verified, rec_not_verifiable, rec_wrong = \
                guess_recs(previously_guessed, number, ens, threads, max_order,
                           c_fin, max_degree)
            # only have to consider sequences which were already guessed
            guessed_seqs.update(guessed)
        
        previously_guessed = guessed_seqs
        
        filename = path_pickles + "/" + name_pickles + f"ens{ens}"
        
        # temporary difference: do not save all terms in the file
        #seqs_to_save = {key : (seq[0], seq[1][:5], seq[2]) for key, seq in guessed_seqs.items()}
        #save(seqs_to_save, filename)
        
        save(guessed_seqs, filename)
        num_dfin = len(guessed_seqs)
        print("========================================================")  
        print(f"Results for ensure={ens}: ")
        print(f"{rec_verified} recurrences verified to be correct")
        print(f"{rec_not_verifiable} recurrences not verifiable to be correct")
        print(f"{rec_wrong} recurrences found to be wrong")
        print(f"{num_dfin} sequences saved at {filename}")
        print("========================================================")  
        
        
        logging.error(f"======================================================")
        logging.error(f"Results for ensure={ens}: ")
        logging.error(f"{rec_verified} recurrences verified to be correct")
        logging.error(f"{rec_not_verifiable} recurrences not verifiable to be correct")
        logging.error(f"{rec_wrong} recurrences found to be wrong")
        logging.error(f"{num_dfin} sequences saved at {filename}")
        logging.error(f"======================================================")
        
            
def load_pickled_seqs(ensures, path_pickles, name_pickles):
    r"""
    Returns a list of names and recurrences of the saved sequences.
    
    INPUT:

    - ``ensures`` -- a list of natural numbers indicating the certainty-levels
      that were used 
    - ``path_pickles`` -- the path where the recurrences are saved
    - ``name_pickles`` -- the file names for the saved recurrences
    
    OUTPUT:
    
    A pair ``(names_ensures, recs_ensures)`` where the first is a list of
    names ``A......`` of the sequences for which a recurrence could be computed
    and the latter is a list of the computed recurrences.
    """
    names_ensures = []
    recs_ensures = []
    
    for ens in ensures :
        
        filename = path_pickles + "/" + name_pickles + f"ens{ens}"
        with open(filename + ".pickle", 'rb') as handle:
            seqs = pickle.load(handle)
        recs = [seq[2] for seq in seqs.values()]
        recs_ensures.append(recs)
        names_ensures.append(list(seqs.keys()))
        
    return names_ensures, recs_ensures
    
def save(obj, path) :
    r"""
    Pickles and saves an object at a given path.
    
    INPUT:
    
    - ``obj`` -- any Python object that can be pickled
    - ``path`` -- a path where the object should be saved
    """
    with open(path + '.pickle', 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
def make_path(path) :
    r"""
    If ``path`` is not empty, create a folder at the given location.
    
    INPUT:
    
    - ``path`` -- a string representing a path 
    """
    if path != "" :
        try:
            os.makedirs(path)
        except FileExistsError:
            # directory already exists
            pass
    
def parse_file(file) :
    r"""
    Returns a list of indices and a list of the corresponding
    terms of a b-file from the OEIS.
    
    INPUT:
    
    - ``file`` -- the b-file of an OEIS sequence as returned from
      the ``urllib.request.urlopen`` method
      
    OUTPUT:
    
    A tuple ``(indices, terms)`` where ``indices`` and ``terms`` are 
    lists of integers such that ``seq[indices[i]] = terms[i]`` for
    the sequence ``seq`` under consideration.
    """
    indices = []
    terms = []
    for line in file :
        line = line.decode("utf-8").strip()
        if not line or line[0] == "#":
            continue
        line_split = re.search(r"(-?\b\d+)\s+(-?\b\d+)", line)
        if not line_split :
            raise ValueError(f"Something went wrong splitting line {line}")
        else :
            index, term = line_split.groups()
        
        indices.append(int(index))
        terms.append(int(term))
        
    return indices, terms

def check_indices(indices) :
    r"""
    Check if indices are a consecutive subset of the natural numbers.
    
    INPUT:
    
    - ``indices`` -- a list of integers
    
    OUTPUT:
    
    ``True`` if ``indices`` is a consecutive subset of the natural numbers
    and ``False`` otherwise.
    """
    for i in range(1, len(indices)) :
        if indices[i] != indices[i-1]+1 :
            raise ValueError("Indices not consecutive")
    return True

def get_terms(index) :
    r"""
    Given an index of an OEIS sequence (excluding the "A"), returns
    a list of the terms of this sequence as they are given in the
    OEIS (in the corresponding b-file).
    
    INPUT:
    
    - ``index`` -- an integer
    
    OUTPUT:
    
    A triple ``(id, start, terms)`` where ``id`` is the identifier
    of the sequence, i.e., a string of the form ``A......``.
    ``start`` is a natural number which indicates the index of the
    first term of the sequence and ``terms`` is a list of terms
    of the integer sequence.
    
    The method checks whether all terms are consecutive in the
    corresponding b-files.
    """
    url = "http://oeis.org/A{index}/b{index}.txt"
    index_string = str(index).zfill(6)
    current_url = url.format(index=index_string)
    try :
        website_data = urllib.request.urlopen(current_url)
        if website_data == None :
            raise URLError("Return is None")
            
        indices, terms = parse_file(website_data)
        check_indices(indices)
        if indices :
            # otherwise sequence is just the empty sequence
            return "A" + index_string, indices[0], terms
        else :
            raise ValueError("No terms given in the file")
        
    except URLError as e :
        print(f"File {current_url} could not be opened. Error: {e}")
        logging.error(f"File {current_url} could not be opened. Error: {e}")
        
    except ValueError as e :
        print(f"ValueError with file {current_url}. Error: {e}")
        logging.error(f"ValueError with file {current_url}. Error: {e}")
