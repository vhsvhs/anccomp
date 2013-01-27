"""
 USAGE:
 %> python compare_ancs.py --specpath <path1> --modelpath <path2> --window_sizes x y z

 path1 is the filepath to a specification file, formatted as shown below.
 path2 is the filepath to a text file containing a substitution matrix, such as JTT, WAG, LG, etc.
 "x y z" are window sizes for smoothing of the final data

 SPECIFICATION FILE uses the following format:
 msapaths <1> <2> <3> etc.
 seed <taxa name>
 compare <node a> <node b> <msapath>
 compare . . . .

 WINDOW SIZES
 The list of desired sizes can be arbitrarily long.

 OPTIONAL RUNTIME COMMANDS:
 --limstart X, analysis will be limited to sites X and beyond in the alignment

"""

from config import *
from anccomp_tools import *
from pymol_viz import *

show_splash()


"""
Part 1:
Read the input.
"""
specpath = ap.getArg("--specpath")
modelpath = ap.getArg("--modelpath")
winsizes = ap.getList("--window_sizes", type=int)
adir = get_anccomp_dir(ap)
if not os.path.exists(adir):
    os.system("mkdir " + adir)
[msa_nodes, seed, msa_weights] = read_specs( specpath ) #msa_nodes[msapath] = [anc. node 1 path, anc. node 2 path]
m = get_matrix(modelpath) # m is the Markovian substitution matrix.


""" 
Part 2:
Align the alignments.
"""
[longest_msa, msa_refsite_mysite, msa_mysite_refsite] = align_msas(msa_nodes.keys(), seed)


"""
Part 3:
Compare ancestors, store results into msa_hscores
"""
w = 1
msa_hscores = {} # key = msa path, value = data from compare_dat_files
msa_pscores = {}
for msapath in msa_nodes:
    this_ancpath = msa_nodes[msapath][0]
    that_ancpath = msa_nodes[msapath][1]
    print "\n. I'm comparing the ancestor [", this_ancpath, "] to [", that_ancpath, "] with a smoothing window =", w
    hdata = compare_dat_files(this_ancpath, that_ancpath, m, w)    
    msa_hscores[msapath] = hdata
    pdata = compare_dat_files(this_ancpath, that_ancpath, m, w, method="p")    
    msa_pscores[msapath] = pdata
    cdata = compute_cdata(msapath, w)


"""
Part 4:
Integrate the results (i.e. msa_hscores) over multiple alignments.
"""
blended_hdata = blend_msa_data(msa_hscores,msa_refsite_mysite,longest_msa,msa_weights)
blended_pdata = blend_msa_data(msa_pscores,msa_refsite_mysite,longest_msa,msa_weights)

"""
Part 5:
Write output, including text tables and a PDF plot.
"""
write_table(blended_hdata, msa_hscores, msa_refsite_mysite, longest_msa, method="h")
rank_sites(blended_hdata, msa_nodes, msa_refsite_mysite, msa_hscores, method="h")

write_table(blended_pdata, msa_pscores, msa_refsite_mysite, longest_msa, method="p")
rank_sites(blended_pdata, msa_nodes, msa_refsite_mysite, msa_pscores, method="p")

cranpaths = []  # this array will hold paths to R scripts that we'll execute (in R) at the end.    
for w in winsizes:
    blended_hdata = window_analysis(blended_hdata, w, ap)
    blended_pdata = window_analysis(blended_pdata, w, ap)
        
    """ 
    Plot H scores, using R 
    """
    plot_outpath = get_plot_outpath(ap, tag=("Hw=" + w.__str__()) )
    combo_substring = ""
    if False != ap.getOptionalArg("--combo_method"):
        combo_substring = ", " + ap.getOptionalArg("--combo_method")
    weight_substring = ""
    if False != ap.getOptionalArg("--training_weight"):
        weight_substring = ", weight=" + ap.getOptionalArg("--training_weight")
    plot_title = "H score, winsize = " + w.__str__() + combo_substring + weight_substring + ", " + time.asctime().__str__() + ""
    cranpath = plot( blended_hdata, plot_outpath, plot_title, "H score", "blue")
    cranpaths.append(cranpath)    

    """ 
    Plot P scores, using R 
    """
    plot_outpath = get_plot_outpath(ap, tag=("Pw=" + w.__str__()) )
    combo_substring = ""
    if False != ap.getOptionalArg("--combo_method"):
        combo_substring = ", " + ap.getOptionalArg("--combo_method")
    weight_substring = ""
    if False != ap.getOptionalArg("--training_weight"):
        weight_substring = ", weight=" + ap.getOptionalArg("--training_weight")
    plot_title = "P score, winsize = " + w.__str__() + combo_substring + weight_substring + ", " + time.asctime().__str__() + ""
    cranpath = plot( blended_pdata, plot_outpath, plot_title, "P score", "blue")
    cranpaths.append(cranpath)    

"""
Now execute all the R scripts.
"""
for c in cranpaths:
    print "\n. I'm plotting results, using the R script written at [", c, "] . . ."
    os.system("r --no-save < " + c)

"""
Part 6:
If specified, visualize H scores on PyMol structure...
"""
if False != ap.getOptionalArg("--pdb_path"):
    seedseq = get_seq_from_msa(seed, longest_msa)
    do_pymol_viz(ap, blended_hdata, seedseq)

print "\n\n. Finished.  Results were written to the folder", get_plot_outpath(ap)
print "\n. Goodbye."
