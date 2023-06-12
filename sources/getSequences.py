# this file is used to use guessing on the sequences from the
# OEIS to determine sequences which are C-finite and D-finite

from sage.all import *
from utility import *

if __name__ == '__main__':

    #---------------------------------------------
    #           change parameters here
    #---------------------------------------------

    # name of the folder containing the OEIS sequences
    name = "all_sequences"
    # the maximal order of the c-finite recurrences which are guessed
    max_order_cfin = 100
    # the maximal order of the d-finite recurrences which are guessed
    max_order_dfin = 200
    # only d-finite recurrences of order r, degree d with 
    # (r+1)(d+1) <= max_degree are guessed
    max_degree = 200
    # confidence levels for the guessing
    checks = [1, 3, 5, 8, 10, 15] 
    # the number of sequences used (-1 for all or e.g. 1000 to only
    # use sequences A000001 up to A001000)
    number = -1
    # path where the pickled D-finite results should be saved
    # leave empty to not execute code for D-finite sequences
    path_pickles_dfin = "dfin"
    # path where the pickled C-finite results should be saved
    # leave empty to not execute code for C-finite sequences
    path_pickles_cfin = "cfin"
    # names of the files the pickled results should have
    name_pickles = ""
    # number of threads used in the computations
    threads = 7
    # should sequences be guessed again; 
    # if False only the evaluation is done
    guess_sequences = True
    # for creating the statistical evaluations: 
    # the total number of OEIS sequences under consideration
    total = 359659
    # filename for the log file
    logfile = "guess.log"
    # file name of scraped sequences from Wiki page
    # leave empty if this comparison should not be done
    file_scraped = "scraped_cfin"
    
    logging.basicConfig(filename=logfile,
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.ERROR)
    
    #---------------------------------------------
    

    make_path(path_pickles_dfin)
    make_path(path_pickles_cfin)
      
    if path_pickles_cfin != "" and guess_sequences :
        print("Start guessing C-finite sequences")
        save_recs(name, number, checks, threads, path_pickles_cfin,
                  name_pickles, max_order_cfin, c_fin = True) 
    
    if path_pickles_cfin != "" :
        seqs_names_ensures, recs_ensures \
            = load_pickled_seqs(checks, path_pickles_cfin, name_pickles)        
        
        save_orders_plot(recs_ensures, path_pickles_cfin, checks, total)
        
        print("Order comparison plot saved")
        
        if file_scraped != "" :
            with open(f"{file_scraped}.pickle", 'rb') as handle:
                scraped = pickle.load(handle)
            
            print_compare_scraped(seqs_names_ensures, scraped, checks)
         
    if path_pickles_dfin != "" and guess_sequences :
        print("Start guessing D-finite sequences")
        save_recs(name, number, checks, threads, path_pickles_dfin,
                  name_pickles, max_order_dfin, c_fin = False,
                  max_degree=max_degree) 
    
    if path_pickles_dfin != "":
        save_name_ord_deg_plot = path_pickles_dfin + "/" + "order_degree_plots"
        make_path(save_name_ord_deg_plot)

        _, recs_ensures = load_pickled_seqs(checks, path_pickles_dfin, 
                                            name_pickles)
        
        for recs, ens in zip(recs_ensures, checks) :
            save_order_degree_plot(recs, save_name_ord_deg_plot, ens)
            
        print("Order/degree matrices plots saved")
        
        save_orders_plot(recs_ensures, path_pickles_dfin, checks, total)
        
        print("Order comparison plot saved")
        
        save_degrees_plot(recs_ensures, path_pickles_dfin, checks, total)
        
        print("Degree comparison plot saved")
    
            
