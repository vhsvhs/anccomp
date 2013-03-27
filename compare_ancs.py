"""
 USAGE:
 %> python compare_ancs.py --specpath <path1> --modelpath <path2> --window_sizes x y z 
 
 path1 is the filepath to a specification file, with formatting explained below.
 
 path2 is the filepath to a text file containing a substitution matrix, such as JTT, WAG, LG, etc.
 
 "x y z" are window sizes for smoothing of the final data

 The SPECIFICATION FILE uses the following format:
 msapaths <1> <2> <3> etc.
 seed <taxa name>
 compare <node a> <node b> <msapath>
 compare . . . .

 WINDOW SIZES
 The list of desired sizes can be arbitrarily long.

 OPTIONAL PARAMETERS:
 --combo_method // product or sum

 --limstart N // limit the analysis to sites N and beyond

 --limstop M // limit the analysis to all sites up to and including M.

 --restrict_sites <list> // limit the analysis to the sites in the list. This option cannot be used with --limstop or --limstart.

 --metrics <list> // Compare the ancestors using the metrics in the list.  Metric options include c, h, hb, and p.

 --skip_plots True // If True, then R scripts will be written, but not invoked.

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
read_specs( specpath ) #msa_comparisons[msapath] = [anc. node 1 path, anc. node 2 path]
m = get_matrix(modelpath) # m is the Markovian substitution matrix.

metrics = ap.getOptionalList("--metrics")
if metrics == None:
    metrics = ["h"]



""" 
Part 2:
Align the alignments, using the method align_msas().
The mapping of sites from the reference MSA to other MSAs will be written to two hashtables:
    ap.params["msa_refsite2mysite"]
    and ap.params["msa_mysite2refsite"]
The reference alignment (i.e. the longest alignment) will be saved as the string ap.params["longest_msa"].
Invariant sites will be listed in ap.params["invariant_sites"]
"""
align_msas()

"""
NOTE:
* ap.params["msa_refsite2mysite"] = hashtable. key = msa path, value = hashtable, where key = site number in the longest MSA,
and value = site number in the msa path
* ap.params["msa_mysite2refsite"] = hashtable. key = msa path, value = hashtable, where key = site number in msa path, value
= site number in the longest MSA.
"""

"""
Part 2b: Restrict the analysis to a subset of total sites (this is optional, if 
the user specified --limstart, --limstop, or --restruct_sites.).
At the end of this section, ap.params["rsites"][msapath] will contain the appropriate
sites on which to restrict this analysis for the MSA at msapath.

"""
print "\n. I'm building a restriction site library. . ."

limstart = 1
lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
limstop = ap.params["msa_refsite2mysite"][ lnick ].keys().__len__() + 1

ap.params["rsites"] = {}
for msa in ap.params["msa_nick2path"]:
    ap.params["rsites"][msa] = []


if (ap.doesContainArg("--limstart") or ap.doesContainArg("--limstop")) and ap.doesContainArg("--restrict_sites"):
    print "\n. Sorry, I'm confused.  You cannot specific --restrict_sites along with --limstart or --limstop."
    print ". I'm quitting."
    exit()

"""Option 1: the user specified --limstart and/or limstop."""
x = ap.getOptionalArg("--limstart")
if x != False:
    limstart = int(x)
x = ap.getOptionalArg("--limstop")
if x != False:
    limstop = int(x)
for site in range(limstart, limstop):
    if site in ap.params["msa_refsite2mysite"][lnick]:
        #print "adding site:", site
        #print ap.params["rsites"].keys()
        ap.params["rsites"][lnick].append(site) 



""" Option 2: the user specified --restrict_sites."""
x = ap.getOptionalList("--restrict_sites")
if x != None:
    for i in x:
        ap.params["rsites"][lnick].append(i)
else:
    for i in range(limstart, limstop+1):
        ap.params["rsites"][lnick].append(i)

# Map the restriction sites onto each MSA:
for site in ap.params["rsites"][lnick]:
    for msanick in ap.params["msa_nick2path"]:
        if msanick != lnick:
            if site in ap.params["msa_refsite2mysite"][msanick]:
                ap.params["rsites"][msanick].append( ap.params["msa_refsite2mysite"][msanick][site] )

# Cull the invariant sites from our analysis:
for site in ap.params["invariant_sites"]:
    if site in ap.params["rsites"][lnick]:
        ap.params["rsites"][ lnick ].pop(site)
        for msanick in ap.params["msa_nick2path"]:
            if site in ap.params["msa_refsite2mysite"][msanick]:
                ap.params["rsites"][msanick].pop( ap.params["msa_refsite2mysite"][msanick][site] )    

#print "\n. I'm excluding the following invariant sites:", 
#print ap.params["invariant_sites"]

#print "\n. The analysis will focus on the following sites only:"
#print ap.params["rsites"]

"""
Part 2c:
Calculate statistics comparing the ancestral DAT files.
"""
#print "\n. I'm calculating distances between your ancestral .dat files. . ."
#
#old_ancs = []
#new_ancs = []
#ap.params["anc_msasource"] = {}
#for msapath in ap.params["msa_path2nick"]:
#    this_ancpath = ap.params["msa_path2nick"][msapath][0]
#    that_ancpath = ap.params["msa_path2nick"][msapath][1]
#    #old_ancs.append(this_ancpath)
#    #new_ancs.append(that_ancpath)
#    ap.params["anc_msasource"][ this_ancpath ] = msapath
#    ap.params["anc_msasource"][ that_ancpath ] = msapath
#
#fout = open(adir + "/dat_comparison.txt", "w")
#for o in old_ancs:
#    for p in old_ancs:
#        if o != p:
#            #anc_data = {}
#            #anc_data[o] = anc_ppvals[o]
#            #anc_data[p] = anc_ppvals[p]
#            x = compare_dat_files_stats(o,p,ap.params["msa_mysite2refsite"], ap.params["msa_refsite2mysite"])
#            print o, p, x
#            fout.write(o + "\t")
#            fout.write(p + "\t")
#            fout.write( "%.4f"%(float(x[1])/x[0]) + "\t")
#            fout.write( "%.4f"%(float(x[2])/x[0]) + "\t")
#            fout.write( "%.4f"%(float(x[3])/x[0]) + "\t")
#            fout.write( "%.4f"%(float(x[4])/x[0]) + "\t")
#            fout.write( x[0].__str__() + "\n")        
#fout.close()

"""
Part 3:
Compare the ancestors.
"""
print "\n. I'm comparing ancestors using the metrics:", metrics

w = 1
metric_data = {} # key = metric ID, value = hashtable of data
if "h" in metrics:
    metric_data["h"] = {}
if "hb" in metrics:
    metric_data["hb"] = {}
if "p" in metrics:
    metric_data["p"] = {}

for msanick in ap.params["msa_nick2path"]:
    this_ancpath = ap.params["msa_comparisons"][msanick][0]
    that_ancpath = ap.params["msa_comparisons"][msanick][1]
    print "\n. I'm comparing the ancestor [", this_ancpath, "] to [", that_ancpath, "] with a smoothing window =", w
    for metric in metrics:
        metric_data[metric][msanick] = compare_dat_files(this_ancpath, that_ancpath, m, w, method=metric)    

"""
Part 4:
Integrate the results over multiple alignments.
"""
#print metric_data
#exit()
metric_blendeddata = {} # key = metric ID, value = hashtable of blended (i.e. integrated) data
for metric in metrics:
    print "\n. I'm integrating the data for", metric, "from across all the alignments."
    metric_blendeddata[metric] = blend_msa_data(metric_data[metric])

"""
Part 5:
Write output, including text tables and PDF plots.
"""
cranpaths = []  # this array will hold paths to R scripts that we'll execute (in R) at the end.    

metric_ranked = {}
for metric in metrics:
    write_table(metric_blendeddata[metric], metric_data[metric], method=metric)
    metric_ranked[metric] = rank_sites(metric_blendeddata[metric], metric_data[metric], method=metric)

if metrics.__len__() > 1:
    metric_comparisons = []
    for i in range(0, metrics.__len__()):
        for j in range(i, metrics.__len__()):
            metric_comparisons.append( (metrics[i], metrics[j]) )
    for comparison in metric_comparisons:
        this_metric = comparison[0]
        that_metric = comparison[1]
        for cranpath in correlate_metrics(metric_ranked[this_metric], metric_ranked[that_metric], metric_blendeddata[this_metric], metric_blendeddata[that_metric], this_metric + "-" + that_metric):
            cranpaths.append( cranpath )
        
w_metric_blendeddata = {} # key = window size for smoothing, value = hashtable, where key = metric ID, value = blended data
for w in winsizes:    
    w_metric_blendeddata[w] = {}
    for metric in metrics:
        w_metric_blendeddata[w][metric] = {}
        w_metric_blendeddata[w][metric] = window_analysis(metric_blendeddata[metric], w, ap)
    
    if w == 1:
        plot_outpath = get_plot_outpath(ap, tag=("histo" ) )
        for cranpath in plot_histogram(w_metric_blendeddata[w], plot_outpath):
            cranpaths.append( cranpath )
        
        stats_outpath = get_plot_outpath(ap, tag=("stats" ) ) + ".txt"
        fout = open(stats_outpath, "w")
        fout.write("metric\t mean\t s.d.\n")
        for metric in w_metric_blendeddata[w]:
            values = []
            data = w_metric_blendeddata[w][metric]
            for site in data:
                values.append( data[site] )
            mean_val = mean(values)
            sdev = sd(values)
            
            fout.write(metric + "\t %.3f"%mean_val + "\t %.3f"%sdev + "\n")
        fout.close()

        
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
if False == ap.getOptionalArg("--skip_plots"):
    for c in cranpaths:
        print "\n. I'm plotting results, using the R script written at [", c, "] . . ."
        os.system("r --no-save < " + c)

"""
Part 6:
If specified, visualize H scores on PyMol structure...
"""
#if False != ap.getOptionalArg("--pdb_path"):
#    seedseq = get_seq_from_msa(seed, ap.params["longest_msa"])
#    do_pymol_viz(ap, blended_hdata, seedseq)

print "\n\n. Finished.  Results were written to the folder", get_plot_outpath(ap)
print "\n. Goodbye."
