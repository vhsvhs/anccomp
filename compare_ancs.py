"""
 USAGE:
 %> python compare_ancs.py --specpath <path1> --modelpath <path2> --window_sizes x y z 
 
 OPTIONAL PARAMETERS:
 --combo_method
 --training_weight
 --limstart
 --metrics c, h, hp, and/or p

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

metrics = ap.getOptionalList("--metrics")
if metrics == None:
    metrics = ["h"]



""" 
Part 2:
Align the alignments.
"""
[longest_msa, msa_refsite_mysite, msa_mysite_refsite] = align_msas(msa_nodes.keys(), seed)
"""
NOTE:
* longest_msa = the path to the longest of the alignments.  We'll use the numbering in the longest MSA
as the reference numbering when we integrate results from multiple alignments.
* msa_refsite_mysite = hashtable. key = msa path, value = hashtable, where key = site number in the longest MSA,
and value = site number in the msa path
* msa_mysite_refsite = hashtable. key = msa path, value = hashtable, where key = site number in msa path, value
= site number in the longest MSA.
"""


"""
Part 3:
Compare the ancestors.
"""
w = 1
metric_data = {} # key = metric ID, value = hashtable of data
if "h" in metrics:
    metric_data["h"] = {}
if "hp" in metrics:
    metric_data["hp"] = {}
if "p" in metrics:
    metric_data["p"] = {}

for msapath in msa_nodes:
    this_ancpath = msa_nodes[msapath][0]
    that_ancpath = msa_nodes[msapath][1]
    print "\n. I'm comparing the ancestor [", this_ancpath, "] to [", that_ancpath, "] with a smoothing window =", w
    for metric in metrics:
        metric_data[metric][msapath] = compare_dat_files(this_ancpath, that_ancpath, m, w, method=metric)    

"""
Part 4:
Integrate the results over multiple alignments.
"""
metric_blendeddata = {} # key = metric ID, value = hashtable of blended (i.e. integrated) data
for metric in metrics:
    metric_blendeddata[metric] = blend_msa_data(metric_data[metric],msa_refsite_mysite,longest_msa,msa_weights)

"""
Part 5:
Write output, including text tables and a PDF plot.
"""
metric_ranked = {}
for metric in metrics:
    write_table(metric_blendeddata[metric], metric_data[metric], msa_refsite_mysite, longest_msa, method=metric)
    metric_ranked[metric] = rank_sites(metric_blendeddata[metric], msa_nodes, msa_refsite_mysite, metric_data[metric], method=metric)

if metrics.__len__() > 1:
    comparisons = []
    for i in range(0, metrics.__len__()):
        for j in range(i, metrics.__len__()):
            comparisons.append( (metrics[i], metrics[j]) )
    for comparison in comparisons:
        this_metric = comparison[0]
        that_metric = comparison[1]
        r = correlate_metrics(metric_ranked[this_metric], metric_ranked[that_metric], metric_blendeddata[this_metric], metric_blendeddata[that_metric], this_metric + "-" + that_metric)
        
cranpaths = []  # this array will hold paths to R scripts that we'll execute (in R) at the end.    
w_metric_blendeddata = {} # key = window size for smoothing, value = hashtable, where key = metric ID, value = blended data
for w in winsizes:    
    w_metric_blendeddata[w] = {}
    for metric in metrics:
        w_metric_blendeddata[w][metric] = {}
        w_metric_blendeddata[w][metric] = window_analysis(metric_blendeddata[metric], w, ap)
    
    if w == 1:
        plot_outpath = get_plot_outpath(ap, tag=("histo" ) )
        cranpath = plot_histogram(w_metric_blendeddata[w], plot_outpath)
        
    for metric in metrics:
        plot_outpath = get_plot_outpath(ap, tag=(metric + ".w=" + w.__str__()) )
        combo_substring = ""
        if False != ap.getOptionalArg("--combo_method"):
            combo_substring = ", " + ap.getOptionalArg("--combo_method")
        weight_substring = ""
        if False != ap.getOptionalArg("--training_weight"):
            weight_substring = ", weight=" + ap.getOptionalArg("--training_weight")
        plot_title = metric + " score, winsize = " + w.__str__() + combo_substring + weight_substring + ", " + time.asctime().__str__() + ""
        cranpath = plot( w_metric_blendeddata[w][metric], plot_outpath, plot_title, metric + " score", "blue")
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
#if False != ap.getOptionalArg("--pdb_path"):
#    seedseq = get_seq_from_msa(seed, longest_msa)
#    do_pymol_viz(ap, blended_hdata, seedseq)

print "\n\n. Finished.  Results were written to the folder", get_plot_outpath(ap)
print "\n. Goodbye."
