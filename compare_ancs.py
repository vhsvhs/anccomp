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

 --restrict_to_seed

 --metrics <list> // Compare the ancestors using the metrics in the list.  Metric options include Df, p, k

 --skip_plots True // If True, then R scripts will be written, but not invoked.

 --renumber_site True // The output plots and tables will use site numbers relative to the seed sequence, rather than to the absolute alignment.

 --pdb_start_site X // at what site in the seed does the PDB file begin with?

"""

from config import *
from anccomp_tools import *
from pymol_viz import *

show_splash()

"""
Part 1:
Read the input.
"""
read_cli(ap)


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
Part 2b: Restrict the analysis to a subset of total sites.
This part is optional, and will be invoked only if the user specified 
--limstart, --limstop, --restruct_sites, or --restrict_to_seed.
At the end of this section, ap.params["rsites"][msapath] will contain the appropriate
sites on which to restrict this analysis for the multiple sequence alignment at msapath.
"""
build_rsites()   

# Write the meta-alignment to a text file:
#if msapaths.__len__() > 1:
write_meta_alignment(ap)

"""
Part 3:
Compare the ancestors.
"""
if ap.params["metrics"].__len__() > 1: # plural. . .
    print "\n. I'm comparing ancestors using the following metrics:", ap.params["metrics"]
else: # . . . or singular
    print "\n. I'm comparing ancestors using the metric", ap.params["metrics"][0]    

[msa_changes, metric_data] = compare_ancestors() # metric_data[metric][msa nickname][site] = score

"""
Part 4:
Integrate the results over multiple alignments.
"""
metric_blendeddata = {} # key = metric ID, value = hashtable of blended (i.e. integrated) data
for metric in ap.params["metrics"]:
    print "\n. I'm integrating the data for", metric, "from across all the alignments."
    metric_blendeddata[metric] = blend_msa_data(metric_data[metric])


"""
Part 5: Write output tables and plot PDFs.

Part 5a: Write a summary table with all sites and their scores.
"""
metric_ranked = rank_all(metric_data, metric_blendeddata)
write_summary_table(metric_blendeddata, metric_data, metric_ranked)

print "108"

"""
Part 5b: Rank the sites, and correlate metrics (only if multiple metrics were used).
"""
cranpaths = []  # this array will hold paths to R scripts that we'll execute (in R) at the end. 
cranpaths += correlate_all(metric_ranked, metric_blendeddata)

"""
Part 5c: Smooth the data, using a user-specified window.
         Scripts for R plots are also written in this method. . .
"""       
cranpaths += smooth_data(metric_blendeddata)     


"""
Now execute all the R scripts.
"""
if False == ap.getOptionalToggle("--skip_plots"):
    for c in cranpaths:
        print "\n. I'm plotting results, using the R script written at [", c, "] . . ."
        os.system("/usr/bin/Rscript " + c)
else:
    print "\n. I'm skipping the result plots."
    print ". You can build these plots later by invoking the following R scripts:"
    for c in cranpaths:
        print "   ", c


#"""
#Part 6:
#If specified, visualize H scores on PyMol structure...
#"""
if False != ap.getOptionalArg("--pdb_path"):
    #refseq = # the ML sequence of the ancestor
    # ap.params["msa_seedseq"][ap.params["longest_msa"]
    # ap.params["msa_comparisons"][lnick][1]
    
    do_pymol_viz(ap, metric_blendeddata)


#write_html_summary(ap)


print "\n\n. I'm finished.  The results were written to the folder", get_plot_outpath(ap)
print "\n. Goodbye."
