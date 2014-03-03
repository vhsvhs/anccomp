import math, os, re, sys, time
from config import *
from splash import *
import scipy.stats as ss
from argparser import *
ap = ArgParser(sys.argv)

METRIC_COLORS = {}
METRIC_COLORS["Df"] = "black"
METRIC_COLORS["k"] = "red"
METRIC_COLORS["p"] = "blue"

METRIC_WT = {}
METRIC_WT["Df"] = "2.5"
METRIC_WT["k"] = "1.0"
METRIC_WT["p"] = "1.0"

PDIST_WEIGHT = 1.0
KDIST_WEIGHT = 1.0



# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )


def roundup(x, step):
    """Roundup x to the nearest 100"""
    if x > 0:
        return int(math.ceil(x / (1.0*step))) * step
    else:
        return (int(math.ceil(x / (1.0*step))) * step) - step

#####################################################################

def get_output_dir(ap):
    return ap.getArg("--runid")

def get_plot_outpath(ap, tag=""):
    return get_output_dir(ap) + "/" + tag 

def get_table_outpath(ap, tag=""):
    return get_output_dir(ap) + "/" + tag 

def datpath_to_short(msapath):
    tokens = msapath.split("/")
    path = tokens[ tokens.__len__()-1 ]
    path = re.sub(".dat", "", path)
    return path

###################################################
def get_seq_from_msa(seq, msapath):
    fin = open(msapath, "r")
    lines = fin.readlines()
    fin.close()
    for l in lines:
        if l.startswith(seq):
            tokens = l.split()
            return tokens[1]

def fill_missing_states(h):
    for a in AA_ALPHABET:
        if a not in h:
            h[a] = TINYPP
    sum = 0.0
    # normalize to sum of 1.0. . .
    for a in AA_ALPHABET:
        sum += h[a]
    for a in AA_ALPHABET:
        h[a] = h[a] / sum        
    return h

def datline_2_pphash(line):
    tokens = line.split()[1:]
    h = {}
    i = 0
    while i < tokens.__len__():
        if tokens[i].__contains__("-"):
            h["-"] = 1.0
            h = fill_missing_states(h)
            return h
        h[ tokens[i].upper() ] = float(tokens[i+1])
        i += 2
    h = fill_missing_states(h)
    return h


def get_markov_model(ap):
    """Matrix files should contain a 19x19 matrix"""
    filepath = ap.getArg("--modelpath")
    
    fin = open(filepath, "r")
    m = {}
    m_sum = 0.0
    for c in AA_ALPHABET:
        m[c] = {}
    lines = fin.readlines()
    curr_line = 1
    for l in lines:
        if l.__len__() > 3:
            #print "line=",l
            tokens = l.split()
            for i in range(0, tokens.__len__()):
                #print "curr_line=", curr_line
                this = AA_ALPHABET[curr_line]
                that = AA_ALPHABET[i]
                val = float(tokens[i])
                m_sum += val
                m[ this ][ that ] = val
                m[ that ][ this ] = val
                #print this, that, val
            curr_line += 1
        if curr_line > (AA_ALPHABET.__len__() -1):
            break
    fin.close()
    
    """Normalize m."""
    norm_m = {}
    for i in AA_ALPHABET:
        if i not in norm_m:
            norm_m[i] = {}
        for j in AA_ALPHABET:
            if i != j:
                #print i, j
                if j not in norm_m:
                    norm_m[j] = {}
                norm_m[i][j] = m[i][j] / m_sum
                norm_m[j][i] = m[j][i] / m_sum
    #return norm_m
    ap.params["model"] = norm_m

def fill_matrix_rev_values(m):
    for c in AA_ALPHABET:
        for d in AA_ALPHABET:
            if d not in m[c].keys():
                m[c][d] = float(m[d][c])
    return m    

###################################################################
def read_cli(ap):
    """Read parameters specified by the command-line interfance (CLI)"""
    read_specs(ap)
    
    ap.params["winsizes"] = ap.getOptionalList("--window_sizes", type=int)
    if ap.params["winsizes"] == None:
        ap.params["winsizes"] = [1]
    
    adir = get_output_dir(ap)
    if not os.path.exists(adir):
        os.system("mkdir " + adir)

    get_markov_model(ap) # initializes ap.params["model"]    
    
    ap.params["metrics"] = ap.getOptionalList("--metrics")
    if ap.params["metrics"] == None:
        ap.params["metrics"] = ["Df"]

    if ap.doesContainArg("--restrict_sites") or ap.doesContainArg("--restrict_to_seed"):
        ap.params["winsizes"] = [1]
    
def read_specs(ap):
    """Read specifications given in the spec. file (using the --specpath command)."""
    specpath = ap.getArg("--specpath")
    
    ap.params["msa_comparisons"] = {} # key = msa nickname, value = array of *.dat tuples
    ap.params["msa_weights"] = {}
    seed = ""
    ap.params["msa_path2nick"] = {}
    ap.params["msa_nick2path"] = {}
    ap.params["msa_seedseq"] = {}
    fin = open(specpath, "r")
    for l in fin.readlines():
        if l.startswith("#"): # skip lines with comments
            continue
        elif l.startswith("seed"):
            tokens = l.split()
            ap.params["seed"] = tokens[1]
        elif l.startswith("compare"):
            tokens = l.split()
            ap.params["msa_comparisons"][tokens[3]] = (tokens[1],tokens[2])
        elif l.startswith("msaweight"):
            tokens = l.split()
            ap.params["msa_weights"][ tokens[1] ] = float(tokens[2])
        elif l.startswith("msaname"):
            tokens = l.split()
            msapath = tokens[1]
            msanick = tokens[2]
            check_msa(msapath)
            ap.params["msa_path2nick"][ msapath ] = msanick
            ap.params["msa_nick2path"][ msanick ] = msapath
        # "msapaths" is depricated in the spec files, but I'm keeping in the
        # code for legacy-support reasons.
#         if l.startswith("msapaths"):
#             tokens = l.split()
#             for t in tokens[1:]:
#                 ap.params["msa_path2nick"][t] = t
        #if l.startswith("pdb"): # pdb pdb_path homologous_anc_path
        #    tokens = l.split(0)
        #    ap.params["pdbpaint"] = {}
        #    ap.params["pdbpaint"][ tokens[1] ] = tokens[2]
    fin.close()    
    check_specs()
    
    for msapath in ap.params["msa_path2nick"]:
        seedseq = get_seed_seq(msapath, ap.params["seed"])
        ap.params["msa_seedseq"][ msapath ] = seedseq
        
        
def check_specs():
    """This method cross-references the configuration values in the spec. file."""
    for msanick in ap.params["msa_comparisons"]:
        foundit = False
        for m in ap.params["msa_path2nick"]:
            if ap.params["msa_path2nick"][m] == msanick:
                foundit = True
        if foundit == False:
            print "\n. Hmmm, I'm confused."
            print "  In your configuration file you compare nodes from the alignment", msanick
            print "  but I did not find a line with 'msaname' for", msanick
            print "\n. Check your configuration file and try again.\n"
            exit(0)
            #ap.params["msa_path2nick"][ msapath ] = msapath

#    #ap.params["anc_msasource"] = {}
#    for msapath in ap.params["msa_path2nick"]:
#        msanick = ap.params["msa_path2nick"][msapath]
#        this_ancpath = ap.params["msa_comparisons"][msanick][0]
#        that_ancpath = ap.params["msa_comparisons"][msanick][1]
#        #ap.params["anc_msasource"][ this_ancpath ] = msanick
#        #ap.params["anc_msasource"][ that_ancpath ] = msanick


def check_msa(msapath):
    """Verify that the multiple sequence alignment at path 'msapath' is a
    valid Phylip-formatted alignment."""
    if False == os.path.exists(msapath):
        print "\n. Sorry, I can't find your alignment ", msapath
        exit()
    
    
#    fin = open(msapath, "r")
#    for l in fin.xreadlines():
#        tokens = l.split()
#        if tokens.__len__() < 2 or tokens.__len__() > 2:
#            print "\n. Error: I don't think your alignment is in Phylip format."
#            print ". Please check your alignment:", msapath, "\n"
#            exit()
#    fin.close()

def get_msa_len(p):
    #print path
    #os.system("cat " + path)
    #print p
    fin = open(p, "r")
    lines = fin.readlines()
    fin.close()
    return float(lines[0].strip().split()[1])

def get_seed_seq(msapath, seed):
    """Returns the seed sequence, with gaps."""
    seed_seq = ""
    fin = open(msapath, "r")
    for l in fin.readlines():
        if l.__len__() > 1:
            if l.split()[0] == seed:
                tokens = l.split()
                seed_seq = tokens[1]
    fin.close()
    if seed_seq.__len__() == 0:
        print "\n. Hmmm, I couldn't find a seed sequence for ", seed, "in ", msapath
        exit()
    return seed_seq

def align_msas():
    """Aligns each alignment file to the longest alignment.  This allows for
    ancestral inferences to be compared between alignments, despite the fact that
    each alignment may have different lengths, different indel placements, and different
    site numbers.
    At the end of this method, two hashtables will be added to ap.params:
    1. ap.params["msa_refsite2mysite"], key = msa path, value = a nested hashtable, 
        where key = site number in the longest MSA and value = site number in the msa path
    2. ap.params["msa_mysite2refsite"], key = msa path, value = a nested hashtable, 
        where key = site number in msa path and value = site number in the longest MSA.
    """
    
    # Find the longest MSA. . .
    msapaths = ap.params["msa_path2nick"].keys()
    seed = ap.params["seed"]
    longest_len = 0
    longest_msa = None
    for p in msapaths:
        len = get_msa_len(p)
        if len > longest_len:
            longest_len = len
            longest_msa = p
        print "\n. Alignment", ap.params["msa_path2nick"][p], "has", int(len), "sites."
    ap.params["longest_msa"] = longest_msa
    
    if msapaths.__len__() > 1:
        print "\n. I'm aligning the alignments, using taxa", seed, "as the seed."
                    
    # Verify that the seed sequences are actually the SAME sequence in each alignment.
    msa_noindel = {} # key = msa path, value = seed sequence without indels
    for p in msapaths:
        msa_noindel[p] = re.sub("\-", "", ap.params["msa_seedseq"][p])
    
    for p in msapaths:
        if msa_noindel[p] != msa_noindel[ longest_msa ]:
            print "\n. Sorry, I had to stop because the seed sequence is not the same in all your alignments."
            print ". Tip: Check if you used a truncated version of the seed sequence in one of the alignments."
            print re.sub("\-", "", ap.params["msa_seedseq"][p])
            print re.sub("\-", "", ap.params["msa_seedseq"][longest_msa])
            exit()

    # Align the seed sequences. . .
    ap.params["msa_refsite2mysite"] = {} # key = MSA nickname, value = hash, where key = site in longest alg, value = my site (or None)
    ap.params["msa_mysite2refsite"] = {} # the reverse lookup of ap.params["msa_refsite2mysite"]
    for p in msapaths:            
        nick = ap.params["msa_path2nick"][p]
        if p != longest_msa:
            print "\n. . .I'm Aligning ", nick, " to ", ap.params["msa_path2nick"][longest_msa], ""
        ap.params["msa_refsite2mysite"][nick] = {}
        ap.params["msa_mysite2refsite"][nick] = {}
        ref_site = 1
        my_site = 1
        while my_site-1 < ap.params["msa_seedseq"][p].__len__() and ref_site-1 < ap.params["msa_seedseq"][longest_msa].__len__():
            if ap.params["msa_seedseq"][p][my_site-1] == ap.params["msa_seedseq"][longest_msa][ref_site-1]:
                my_state = ap.params["msa_seedseq"][p][my_site-1]
                ref_state = ap.params["msa_seedseq"][longest_msa][ref_site-1]
                if my_state != ref_state:
                    print "\n. Hmmm, something went wrong with the meta-alignment."
                    print " (anccomp_tools.py point 205."
                    print "anccomp 169: matching", my_site, "in", alg, "to", ref_site, "states:", my_state, ref_state    
                    print "\n"
                    exit(1)
                ap.params["msa_refsite2mysite"][nick][ref_site] = my_site
                ap.params["msa_mysite2refsite"][nick][my_site] = ref_site
                ref_site += 1
                my_site += 1
                #continue
            elif ap.params["msa_seedseq"][longest_msa][ref_site-1] == "-":
                while ap.params["msa_seedseq"][longest_msa][ref_site-1] == "-" and ref_site < ap.params["msa_seedseq"][longest_msa].__len__()+1:
                    ref_site += 1
            elif ap.params["msa_seedseq"][p][my_site-1] == "-":
                while ap.params["msa_seedseq"][p][my_site-1] == "-" and my_site < ap.params["msa_seedseq"][p].__len__()+1:
                    my_site += 1
    
    # Identify invariant sites:
    invariant_sites = []
    fin = open(longest_msa, "r")
    lines = fin.readlines()
    taxa_seq = {}
    #last_taxa = ""
    for i in range(1, lines.__len__()):
        if lines[i].__len__() > 1:
            tokens = lines[i]
            taxa = tokens[0]
            seq = tokens[1]
            taxa_seq[taxa] = seq
    fin.close()
    for site in range(0, taxa_seq[ taxa_seq.keys()[0] ].__len__()):
        found_diff = False
        compareto = taxa_seq[ taxa_seq.keys()[0] ][site]
        for taxa in taxa_seq.keys():
            if taxa_seq[taxa][site] != compareto:
                found_diff = True
        if found_diff == False:
            invariant_sites.append( site )
    
    ap.params["invariant_sites"] = invariant_sites
    ap.params["msa_refsite2mysite"] = ap.params["msa_refsite2mysite"]
    ap.params["msa_mysite2refsite"] = ap.params["msa_mysite2refsite"]


def write_meta_alignment(ap):
    """Write a text file with the alignment of alignments."""
    # Write a log file about the meta-alignment,
    # and fill the hash ap.params["msa_refsite2mysite"]
    msa_ids = {} # key = MSA path, value = integer ID
    ids_msa = {} # reverse hash of above
    longest_msa = ap.params["longest_msa"]
    
    msa_ids[longest_msa] = 1
    ids_msa[1] = longest_msa
    
    
    counter = 2
    for path in ap.params["msa_path2nick"]:
        if path != longest_msa:
            msa_ids[ path ] = counter
            ids_msa[ counter ] = path
            counter += 1
    
    pdb_path = ap.getOptionalArg("--pdb_path")
    #if pdb_path != None:
        
    #pdbseqsite2refsite
    
    fout = open(get_output_dir(ap) + "/meta_alignment.txt", "w")
    ids = ids_msa.keys()
    ids.sort()
    fout.write("Alignment Key:\n")
    for i in ids:
        fout.write( "M" + i.__str__() + " : " + ap.params["msa_path2nick"][ ids_msa[i] ]  + "\n")
    fout.write("Seed: ungapped sites in " + ap.params["seed"] + "\n")
    header = ""
    #fout.write("Residues shown in parentheses express the state of taxon " + ap.params["seed"] + " at the corresponding site.\n\n")
    fout.write("The mark 'x' indicates that no corresponding site was found in the alignment.\n\n")
    for i in ids:
        header += "M" + i.__str__() + "\t"
    header += "Seed\t"
    #header += "PDB\t"
    fout.write(header + "\n\n")
    #print ap.params["msa_path2nick"]
    #print ap.params["msa_path2nick"][ longest_msa ]
    
    #print ap.params["msa_mysite2seedsite"]
    #print ap.params["msa_mysite2seedsite"][ ap.params["msa_path2nick"][longest_msa] ]
    
    for site in range(0, ap.params["msa_seedseq"][ longest_msa ].__len__()):
        line = ""
        for i in ids:
            msanick = ap.params["msa_path2nick"][ids_msa[i]]
            if (site+1) in ap.params["msa_refsite2mysite"][ msanick ]:
                mysite = ap.params["msa_refsite2mysite"][ msanick ][ site+1 ]
                line += mysite.__str__() + "\t"
                #state = ap.params["msa_seedseq"][ids_msa[i]][mysite-1]
                #line += mysite.__str__() + " (" + state + ")\t"
            else:
                line += "x\t"
        if (site+1) in ap.params["msa_mysite2seedsite"][ ap.params["msa_path2nick"][longest_msa] ]: 
            seedsite = ap.params["msa_mysite2seedsite"][ ap.params["msa_path2nick"][longest_msa] ][ site+1 ]
            line += seedsite.__str__()
            line += " (" + ap.params["msa_seedseq"][longest_msa][site] + ")"
            line += "\t"
        else:
            line += "x\t"
        fout.write( line + "\n" )
    fout.close()

def build_rsites():
    """Restrict the analysis to a subset of total sites.
    This part is optional, and will be invoked only if the user specified 
    --limstart, --limstop, --restruct_sites, or --restrict_to_seed.
    Each of these parameters specifies the restriction sites in different ways.
    See the user manual for more information.
    At the end of this method, a hashtable will be added to ap.params,
    named ap.params["rsites"], where key = nickname for an MSA, and 
    value is a list of sites that are usable in that MSA.
    """
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    
    ap.params["rsites"] = {}
    for msa in ap.params["msa_nick2path"]:
        ap.params["rsites"][msa] = []
    
    
    if (ap.doesContainArg("--limstart") or ap.doesContainArg("--limstop")) and ap.doesContainArg("--restrict_sites"):
        print "\n. You cannot specific --restrict_sites along with --limstart or --limstop."
        print ". I'm quitting."
        exit()
    
    limstart = 1
    limstop = ap.params["msa_refsite2mysite"][ lnick ].keys().__len__()
    
    # Build the restriction site library, using the site numbers in the longest MSA. . .
    x = ap.getOptionalList("--restrict_sites")
    y = ap.getOptionalToggle("--restrict_to_seed")
    p = ap.getOptionalArg("--limstart")
    q = ap.getOptionalArg("--limstop")
    seed_seq = get_seed_seq(  ap.params["msa_nick2path"][lnick], ap.params["seed"]  )
    
    if x != None:
        for i in x:
            ap.params["rsites"][lnick].append(int(i))
            #print i
    
    elif y == True and p == False and q == False:
        print "\n. I'm restricting the analysis to only those sites found in the seed sequence", ap.params["seed"]
        for site in range(1, limstop+1):
            if seed_seq[site-1] != "-":
                #print "471:", site, seed_seq[site-1]
                ap.params["rsites"][lnick].append(site)
    else:
        if p != False:
            limstart = int(p)
        if q != False:
            limstop = int(q)
        print "\n. I'm restricting the analysis to sites", limstart, "to", limstop, "in", lnick
        if y != False:
            print "  and which are homologous to a non-indel in the seed."
        for site in range(limstart, limstop+1):
            if y != False:
                if seed_seq[site-1] != "-": 
                    ap.params["rsites"][lnick].append(site)
            else:
                ap.params["rsites"][lnick].append(site)                
    
    
    #print "debug 488:", lnick, ap.params["rsites"][lnick]
    #print "debug 489:", ap.params["rsites"][lnick].__len__()
    
    # Map the restriction sites onto the other MSAs. . . 
    for site in ap.params["rsites"][lnick]:
        #print "490:", site
        for msanick in ap.params["msa_nick2path"]:
            if msanick != lnick:
                if site in ap.params["msa_refsite2mysite"][msanick]:
                    #print "debug 492:", site, msanick, lnick, ap.params["msa_refsite2mysite"][msanick][site]
                    ap.params["rsites"][msanick].append( ap.params["msa_refsite2mysite"][msanick][site] )
    # Cull the invariant sites from our analysis:
#    for site in ap.params["invariant_sites"]:
#        if site in ap.params["rsites"][lnick]:
#            ap.params["rsites"][ lnick ].pop(site)
#            for msanick in ap.params["msa_nick2path"]:
#                if site in ap.params["msa_refsite2mysite"][msanick]:
#                    ap.params["rsites"][msanick].pop( ap.params["msa_refsite2mysite"][msanick][site] ) 

    if ap.doesContainArg("--renumber_sites"):
        print ". Renumbering sites..."
        ref2seed = {} # key = reference site, value = corresponding seed sequence site
        count = 0
        for site in ap.params["rsites"][lnick]:
            count += 1
            ref2seed[site] = count
    
        print ". Mapping sites to seed sequence."
        ap.params["msa_mysite2seedsite"] = {}
        for msanick in ap.params["msa_nick2path"]:
            ap.params["msa_mysite2seedsite"][msanick] = {}
            for mysite in ap.params["rsites"][msanick]:
                ap.params["msa_mysite2seedsite"][msanick][mysite] = ref2seed[ ap.params["msa_mysite2refsite"][msanick][mysite] ]
                #print "556:", msanick, mysite, ref2seed[ ap.params["msa_mysite2refsite"][msanick][mysite] ]
    
def compare_dat_files(patha, pathb, msanick, method):
    """Compares two ancestors (at patha and pathb) from the same alignment (with nicknamed msanick).
        Comparison is perfomed using 'method'
        This method returns a hashtable, where key = site, value = comparison score at that site."""
    
    fin = open(patha, "r")
    linesa = fin.readlines()
    fin.close()
    fin = open(pathb, "r")
    linesb = fin.readlines()
    fin.close()
    
    results = {} 
        
    for li in range(0, linesa.__len__()):    
        al = linesa[li]
        bl = linesb[li]        
        asite = int(al.split()[0])
        bsite = int(bl.split()[0])
        if asite != bsite:
            print "\n. Hmm, I've encountered an error. I think your ancestors have a mismatched number of sites."
            print ". I'm quitting."
            exit(1)

        if asite in ap.params["rsites"][msanick]:
            ancx = datline_2_pphash(al)
            ancy = datline_2_pphash(bl) 
            # Compare ancestors X and Y at this site:      
            results[asite] = d(ancx, ancy, method=method)  
    #print "551:", msanick, results.keys().__len__()
    return results

def dat2pp(datpath):
    data = {}
    fin = open(datpath, "r")
    linesa = fin.readlines()
    fin.close()
    
    for li in range(0, linesa.__len__()):    
        l = linesa[li]
        if l.startswith("#") or l.__len__() < 3:
            continue  
        site = int(l.split()[0])
        data[site] = datline_2_pphash(l)
    return data

def dat2seq(datpath, remove_indels=False):
    """Returns the ML sequence from an ancestral *.dat file."""
    fin = open(datpath, "r")
    lines = fin.readlines()
    fin.close()
    mlseq = ""
    for l in lines:
        if l.__len__() > 2:
            tokens = l.split()
            state = tokens[1]
            if (remove_indels == False) or (remove_indels == True and state != "-"):
                mlseq += tokens[1] 
    return mlseq



def compare_ancestors():
    msa_changes = {} # key = msa nickname, value = output from the method 'count_changes_between_ancestors'
    msa_htmlfrags = {}
    metric_data = {}
    for msanick in ap.params["msa_nick2path"]:
        this_ancpath = ap.params["msa_comparisons"][msanick][0]
        that_ancpath = ap.params["msa_comparisons"][msanick][1]
        print "\n. . .I'm comparing the ancestor [", this_ancpath, "] to [", that_ancpath, "]."
    
        # Compare the ancestors in simple terms, by the number of amino acid changes and the number of indel changes, etc.
        [msa_changes[msanick], msa_htmlfrags[msanick]] = count_changes_between_ancestors( this_ancpath, that_ancpath, msanick )
        
        # Compare the ancestors using our statistical metrics
        for metric in ap.params["metrics"]:
            if metric not in metric_data:
                metric_data[metric] = {}
            metric_data[metric][msanick] = compare_dat_files(this_ancpath, that_ancpath, msanick, metric)    
    write_changes_summary(msa_changes)
    return [metric_data, msa_htmlfrags]

def write_changes_summary(msa_changes):

        #print msa, nsites, countindel, countred, countorange, countgreen
    fout = open(get_output_dir(ap) + "/ancestral_changes.txt", "w")
    
    #key = "===============================================================================================\n"
    #key += "A summary of changes between ancestors AncX and AncY. . .\n\n"
    #key += "Key:\n\n"
    #key += "N sites: the total number of sites in AncX and AncY.\n"
    #key += "Indel changes: the number of insertion or deletions between the two ancestors.\n"
    #key += "Type 1: AncY has a different ML state than AncX, and AncX has no support for AncY's state,\n"
    #key += "        and AncY has no support for AncX's state.\n"
    #key += "Type 2: AncY has a different ML state than AncX, but AncX has mild support for AncY's state,\n"
    #key += "        or AncY has mild support for AncX's state.\n"
    #key += "Type 3: AncX and AncY have the same ML state, but either AncX or AncY has strong uncertainty\n"
    #key += "        about this state.\n"
    #key += "==============================================================================================\n\n" 
    #fout.write(key)
    
    header = "Alignment & Model\tN sites\tN indel\tType1\tType2\tType3\n"
    fout.write(header)
    for msa in msa_changes:
        [nsites, countindel, countred, countorange, countgreen] = msa_changes[msa]
        fout.write(msa + "\t" + nsites.__str__() + "\t" + countindel.__str__() + "\t" + countred.__str__() + "\t" + countorange.__str__() + "\t" + countgreen.__str__() + "\n")
    fout.close()
    
def pptransform(ancpp):    
    """Input a PP distribution for one site.
    Output: transformed hashtable with data[pp] = one or more states with this PP."""
    out = {}
    for state in ancpp:
        pp = "%.3f"%ancpp[state]
        if pp not in out:
            out[pp] = [state]
        else:
            out[pp].append(state)
    return out

def get_htmlfrag(ancpp1, ancpp2, rowcolor=None):
    anctrans1 = pptransform(ancpp1)
    anctrans2 = pptransform(ancpp2)
    pps1 = anctrans1.keys()
    pps1.sort(reverse=True)
    pps2 = anctrans2.keys()
    pps2.sort(reverse=True)
    out1 = "<td>1:</td>"
    count1len = 0
    for pp in pps1:
        if float(pp) > 0.0:
            for state in anctrans1[pp]:
                out1 += "<td>" + state + " (" + pp + ")" + "<td>"
                count1len += 1
    count2len = 0
    out2 = "<td>2:</td>"
    for pp in pps2:
        if float(pp) > 0.0:
            for state in anctrans2[pp]:
                out2 += "<td>" + state + " (" + pp + ")" + "<td>"            
                count2len += 1
    
    if count1len < count2len:
        for i in range(0, count2len-count1len):
            out1 += "<td></td>"
    elif count2len < count1len:
        for i in range(0, count1len-count2len):
            out2 += "<td></td>"
    
    out = "<table>"
    style = ""
    #style = "whiterow"
    #if rowcolor == "red":
    #    style = "redrow"
    #elif rowcolor == "orange":
    #    style = "orangerow"
    out += "<tr class=\"" + style + "\">" + out1 + "</tr>"
    out += "<tr class=\"" + style + "\">" + out2 + "</tr>"
    out += "</table>"
    return out

def count_changes_between_ancestors(patha, pathb, msanick):    
    """Counts the number of amino acid differences between two ancestors."""
    anc_data = {}
    anc_data[patha] = dat2pp(patha)
    anc_data[pathb] = dat2pp(pathb)
    
    htmlfrags = {} # output lines of HTML fragments for table1
    
        
    #
    # These variables will be filled with values over the duration
    # of this method. . .
    #
    nsites = 0
    countindel = 0
    countred = 0    # sites with different ML states
    countorange = 0 # sites with different ML states, but same ML+1 states
    countgreen = 0  # sites same ML state, but one ancestor is uncertain
    redsites = []
    orangesites = []
    greensites = []
    
    sites = ap.params["msa_mysite2refsite"][msanick].keys()
    sites.sort()

    for site in sites:
        if site in ap.params["rsites"][msanick]:            
            rowcolor = "None"
            nsites += 1
            indel = False
            red = False
            orange = False
            green = False
            a_state = getmlstate( anc_data[patha][site] )
            b_state = getmlstate( anc_data[pathb][site] )
            if b_state == None and a_state == None:
                continue
            elif b_state != a_state:
                if b_state == None or a_state == None:
                    indel = True
                    #print "indel"
                elif b_state not in anc_data[patha][site]:
                    red = True
                    rowcolor = "red"
                elif a_state not in anc_data[pathb][site]:
                    red = True
                    rowcolor = "red"
                elif False == anc_data[pathb][site].keys().__contains__(a_state):
                    red = True
                    rowcolor = "red"
                    #print "red"
                elif False == anc_data[ patha ][site].keys().__contains__(b_state):
                    red = True
                    rowcolor = "red"
                    #print "red"
                elif anc_data[patha][site][b_state] < 0.05 and anc_data[pathb][site][a_state] < 0.05:
                    red = True
                    rowcolor = "red"
                    #print "red"
                # test for orange:
                else:
                    orange = True
                    rowcolor = "orange"
                    #print "orange"
            # test for green
            elif (anc_data[pathb][site][b_state] < 0.8 and anc_data[patha][site][b_state] > 0.8) or (anc_data[pathb][site][b_state] > 0.8 and anc_data[patha][site][b_state] < 0.8):
                    green = True
                    rowcolor = "green"
            
            
            if indel:
                countindel += 1
            if red:
                countred += 1
                redsites.append(site)
                ##print "site", site, " - case 1"
                ##print "\t",getmlstate( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
                #print "\t",getmlstate( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
            if orange:
                countorange +=1
                orangesites.append(site)
                #print "site", site, " - case 2"
                #print "\t", getmlstate( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
                #print "\t",getmlstate( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
            if green:
                countgreen += 1
                greensites.append(site)
                #print "site", site, " - case 3"
                #print "\t",getmlstate( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
                #print "\t",getmlstate( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
        
            
            #
            # and gather HTML fragments for later, when we write a report.
            #
            htmlfrags[site] = get_htmlfrag( anc_data[patha][site], anc_data[pathb][site], rowcolor )
        
        
    return [ [nsites, countindel, countred, countorange, countgreen], htmlfrags]
    

def compute_cdata(msapath, winsize):
    """Computes basic conservation."""
    c = {} # key = site, value = prportion of sites with same state
    fin = open(msapath, "r")
    lines = fin.readlines()
    fin.close()
    tokens = lines[0].split()
    ntaxa = int(tokens[0])
    fracbit = 1.0 / float(ntaxa)
    nsites = int(tokens[1])
    site_states = {}
    for site in range(1, nsites+1):    
        site_states[site] = {}
    # parse the sequences, put data in site_states
    for l in lines[1:]:
        l = l.strip()
        seq = l.split()[1]
        for site in range(1, nsites+1):
            state = seq[site-1] # because seq is an array with 0-based indexing
            if state in site_states[site]:
                site_states[site][state] += fracbit
            else:
                site_states[site][state] = fracbit
    # calculate C
    for site in range(1, nsites+1):
        
        # percentage-same method
        max_c = 0.0
        for state in site_states[site]:
            if site_states[site][state] > max_c:
                max_c = site_states[site][state]
        c[site] = max_c
        
        # entropy method
        #c[site] = 2.0 - entropy( site_states[site] )
    c = window_analysis(c, winsize, ap)
    return c


def entropy(ppx):
    print "Entropy", ppx
    ret = 0.0
    for i in AA_ALPHABET:
        if ppx[a] > 0:
            ret += (ppx[a] * math.log(ppx[a]))
    return ret
     
# depricated:
#def conditional_entropy(ppx, ppy):
#    """Retuns H(ppy|ppx)"""
#    """Conditional entropy = h(y|x) = h(x,y) - h(x)"""
#    ret = 0.0
#    for a in AA_ALPHABET:
#        for b in AA_ALPHABET:
#            jointp = ppx[a]*ppy[b]
#            if jointp > 0:
#                ret +=  jointp * math.log( jointp )
#    ret = ret - entropy(ppx)
#    return -1 * ret
#
#
#def information_gain(ppx, ppy, m):
#    ret = 0.0
#    ret = entropy(ppy) - conditional_entropy(ppx, ppy)
#    return ret

def geometric_distance(i, j, ppi, ppj, m):
    return abs(ppi - ppj)

def kl_divergence(ppx, ppy):
    """Calculates the one-way Kullback-Leibler divergence between ppx and ppy"""

    # this first part is classic KL...
    ret = TINY    
    
    yvals = ppx
    xvals = ppy
    
    # Set low values to be uncertain (i.e. == 0.05)
#    for a in AA_ALPHABET:
#        if ppy[a] < 0.01:
#            yvals[a] = 0.05
#        else:
#            yvals[a] = ppy[a]
#        if ppx[a] < 0.01:
#            xvals[a] = 0.05
#        else:
#            xvals[a] = ppx[a]
    
    # Normalize. . . .
#    ysum = 0
#    for a in AA_ALPHABET:
#        ysum += yvals[a]
#    xsum = 0
#    for a in AA_ALPHABET:
#        xsum += xvals[a]
#    for a in AA_ALPHABET:
#        yvals[a] = yvals[a]/ysum
#        xvals[a] = xvals[a]/xsum
    
    for a in AA_ALPHABET:
        #if ppy[a] > 0.01 and ppx[a] > 0.01:            
        ret += (xvals[a] * math.log(xvals[a]/yvals[a]))
    return ret

def getmlstate(pp):
    maxa = None
    maxpp = 0.0
    for a in AA_ALPHABET:
        if a in pp:
            if maxpp < pp[a]:
                maxa = a
                maxpp = pp[a]
    return maxa
            
def kdist(ppx, ppy):
    if "-" in ppx.keys() and "-" in ppy.keys(): # indel to indel. . .
        return 0.0

    klxy = kl_divergence(ppx, ppy)
    klyx = kl_divergence(ppy, ppx)
    klsum = klxy + klyx
    return klsum    
            
def pdist(ppx, ppy):
    m = ap.params["model"]
    expected_p = 0.0
    observed_p = 0.0
    if "-" in ppx.keys() and "-" in ppy.keys(): # indel to indel. . .
        return 0.0
    pval = 0.0
    for a in AA_ALPHABET:
        for b in AA_ALPHABET:
            if a != b:
                #observed_p += ppx[a] * ppy[b]
                #expected_p += ppx[a] * math.exp( m[a][b] )# * 0.01 ) # to-do: times branch length?
                ep = ppx[a] * math.exp( m[a][b] )# * 0.01 ) # to-do: times branch length?
                op = ppx[a] * ppy[b]
                pval += abs( ep - op )**2
                #pval += op / ep
                #print "849", a, b, ppx[a], ppy[b], (ppx[a] * ppy[b]), (ppx[a] * math.exp( m[a][b] ))
    #pval = observed_p / expected_p
    return pval

def e_scale(ppx, ppy):
    enty = entropy(ppy)
    entx = entropy(ppx)
    e_dir = 1
    if (entx < enty) and ppy[ getmlstate(ppy) ] < CERTAINTY_CUTOFF and (enty - entx) > ENT_CUTOFF:
        e_dir *= -1
    return e_dir


def entropy(h, ppcull=0.0):
    # first cull low values, if that was requested...
    newh = {}
    for c in h:
        if h[c] >= ppcull:
            newh[c] = h[c]
    # calculate entropy on the filtered array. . .
    epy = 0.0
    for c in newh:    
        epy += newh[c] * math.log( newh[c] )
    return -1.0 * epy

def d(ancx, ancy, method="Df"):
    """Compares a single site in ancestor x to a single site in ancestor Y."""    

    # First, fix the PP vector to contain no zeros.
#    if "-" in ancx.keys() or "X" in ancx.keys():
#        ancx = {}
#        for aa in AA_ALPHABET:
#            ancx[aa] = 0.05
#    if "-" in ancy.keys() or "X" in ancy.keys():
#        ancy = {}
#        for aa in AA_ALPHABET:
#            ancy[aa] = 0.05
            
    pieces = []        
    if method == "Df": # i.e., k times p
        # observed/expectation:
        pieces.append( PDIST_WEIGHT*pdist(ancx, ancy) )
        # KL distance:
        pieces.append( KDIST_WEIGHT*kdist(ancx, ancy) )

    elif method == "k":
        # KL distance only:
        pieces.append( kdist(ancx, ancy) )
    elif method == "p":
        # observed/expectation only:
        pieces.append( pdist(ancx, ancy) )

    score = pieces[0]
    for piece in pieces[1:]:
        score *= piece
        
    score *= e_scale(ancx, ancy)
    return score

def window_analysis(data, winsize, ap):
    winmethod = ap.getOptionalArg("--window_function")
    if winmethod == "geometric":
        return geometric_window_analaysis(data, winsize)
    if winmethod == "dirichlet":
        return dirichlet_window_analysis(data, winsize)
    else:
        return gauss_window_analaysis(data, winsize)


# part of gaussian method
def apply_filter(index,array,window, winsize):
    """apply_filter applies the Gaussian filter."""
    sites = array.keys()
    sites.sort()
    
    #print sites
    
    min = index - winsize
    if min < sites[0]:
        min = sites[0]
        
    max = index + winsize
    if max > sites[ sites.__len__()-1 ]:
        max = sites[ sites.__len__()-1 ]
    sum = 0.0
    #print min, max
    for i in range(min, max):
        if i in array:
            delta = i - min
            sum += array[i] * window[delta]
    #print "site", index, "raw=", array[index], "min=", min, "max=", max, "sum=", sum
    return sum

# part of gaussian methods
def get_window_weights(N):
    """Get Gaussian window scalars."""
    support_points = [(float(3 * i)/float(N))**2.0 for i in range(-N,N + 1)]
    gii_factors = [math.exp(-(i/2.0)) for i in support_points]
    ki = float(sum(gii_factors))
    return [giin/ki for giin in gii_factors]

def gauss_window_analaysis(data, winsize):
    """Smoothing with a Gaussian blur in 1 dimension."""
    window_weights = get_window_weights(winsize)
    ret = {}
    sites = data.keys()
    sites.sort()
    for i in sites:
        ret[i] = apply_filter(i,data,window_weights,winsize)
    return ret

def dirichlet_window_analysis(data, winsize):
    """returns w, w[site] = windows analysis metric for the window starting at site.
    Data to be compared is in h[site]."""
    w = {}
    sites = data.keys()
    sites.sort()
    for s in sites:
        sum = 0.0
        for i in range(s, s+winsize):
            if i <= sites.__len__():
                sum += data[i]
        w[s] = sum / winsize
    return w


def fill_missing_sites(data):
    new_data = {}
    min_s = int( data[data.keys()[0]] ) # lowest site
    max_s = min_s                # highest site
    for s in data.keys():
        if s < min_s:
            min_s = int( s )
        if s > max_s:
            max_s = int( s )
    for i in range(min_s, max_s+1):
        if i in data:
            new_data[i] = data[i]
        else:
            new_data[i] = 0.0
    return new_data

def plot_one_metric(data, outpath, title, ylab, color, image_type="pdf"):
    if False == ap.doesContainArg("--renumber_sites"):
        data = fill_missing_sites(data)
    
    #print "1004:", data
    #exit()
    
    cranpath = outpath + "." + image_type + ".rscript"
    cranout = open(cranpath, "w")
    x = "x <- c("
    y = "y <- c("
    sites = data.keys()
    sites.sort()

    miny = 0.0
    maxy = 0.0
    minx = sites[0]
    allvalues = []
    maxx = sites[ sites.__len__()-1 ]
    if ap.doesContainArg("--renumber_sites"):
        minx = 1
        maxx = sites.__len__()

    
    for s in sites:
        yval = 0
        if data[s] != None:
            if data[s] < 0.00000001:
                yval = 0.0
            else:
                yval = data[s]
        #    else:
        #        yval = math.log( float(data[s]) + 1.0 )
        allvalues.append( yval )
    
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    for s in sites:
        if ap.doesContainArg("--renumber_sites"):
            x += int(ap.params["msa_mysite2seedsite"][lnick][s]).__str__() + ","
        else:
            x += s.__str__() + ","
        if data[s] != None:
            y += data[s].__str__() + ","
            if float(data[s]) > maxy:
                maxy = float(data[s])
            if float(data[s]) < miny:
                miny = float(data[s])
        else:
            y += "0.0,"
    x = re.sub(",$", "", x)
    x += ")"
    cranout.write( x + "\n")
    y = re.sub(",$", "", y)
    y += ")"
    cranout.write( y + "\n")
        
    if image_type == "pdf":
        cranout.write("pdf('" + outpath + ".pdf', width=6.5, height=4);\n")
    elif image_type == "png":
        cranout.write("png('" + outpath + ".png', width=800, height=300);\n")
    cranout.write("plot(c(" + minx.__str__() + "," + roundup(maxx, 100).__str__() + "), c(" + miny.__str__() + "," + maxy.__str__() + "), type='n',xlab='Sites homologous to " + ap.params["seed"] + "', ylab='" + ylab + "', main='" + title + "', lwd=2, las=1, col='" + color + "', bty='n');\n")
    
    # insert any rectangles here, for highlighted sites
    highlight_sites = ap.getOptionalList("--highlight_sites")
    if highlight_sites != None:
        print "\n. For visualization, I will highlight the sites in the range", highlight_sites
        for pair in highlight_sites:
            tokens = pair.split("-")
            cranout.write("rect(" + tokens[0].__str__() + ", " + miny.__str__() + " , " +  tokens[1].__str__()+ ", " + maxy.__str__() + ", col=\"gray92\", lwd=0, lty=\"blank\");\n")
    
    # draw lines for 0, +/- 2 stdev, and +/- 4 stdev
    sdev = sd(allvalues)
    avg = mean(allvalues)
    cranout.write("abline(0,0, lwd=0.75, col='black');\n")
    cranout.write("abline(" + avg.__str__() + ",0, lwd=0.5, lty='dashed', col='gray50');\n")
    cranout.write("abline(" + (avg+2*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='deepskyblue2');\n")
    cranout.write("abline(" + (avg-2*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='deepskyblue2');\n")
    cranout.write("abline(" + (avg+4*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='orange1');\n")
    cranout.write("abline(" + (avg-4*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='orange1');\n")
    cranout.write("abline(" + (avg-6*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='red1');\n")
    cranout.write("abline(" + (avg+6*sdev).__str__() + ",0, lwd=0.5, lty='dashed', col='red1');\n")


    # draw text for stdev lines.
    cranout.write("text(" + maxx.__str__() + "," + (avg).__str__() + ", \"mean\", cex=0.7, col='gray50');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg+2*sdev).__str__() + ", \"+2 sigma\", cex=0.7, col='deepskyblue2');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg-2*sdev).__str__() + ", \"-2 sigma\", cex=0.7, col='deepskyblue2');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg+4*sdev).__str__() + ", \"+4 sigma\", cex=0.7, col='orange1');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg-4*sdev).__str__() + ", \"-4 sigma\", cex=0.7, col='orange1');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg+6*sdev).__str__() + ", \"+6 sigma\", cex=0.7, col='red1');\n")
    cranout.write("text(" + maxx.__str__() + "," + (avg-6*sdev).__str__() + ", \"-6 sigma\", cex=0.7, col='red1');\n")
    
    
    # draw the data
    cranout.write("points(x, y, lwd=2, type='l', col='" + color + "');\n")
    
    cranout.write("dev.off();\n")
    cranout.close()
    return cranpath

def normalize_vector(v):
    """Returns the data normalize to [-1,1]"""
    min = None
    max = None
    norm_v = {}
        
    for site in v.keys():
        i = v[site]
        if i == None:
            i = 0.0
        if min == None:
            min = i
        if max == None:
            max = i
        if i < min:
            min = i
        if i > max:
            max = i
    
    if max > -1*min:
        normalizer = max
    elif min < -1*max:
        normalizer = -1*min
    
    for site in v.keys():
        norm_v[site] = v[site] / normalizer

    return norm_v

def plot_multi_metrics(data, outpath, title, image_type="pdf"):
    """data[metric][site] = value for metric"""
    
    cranpath = outpath + "." + image_type + ".rscript"
    cranout = open(cranpath, "w")

    sites = data[ data.keys()[0] ].keys()
    sites.sort() 
    minx = sites[0]
    maxx = sites[ sites.__len__()-1 ]
    if ap.doesContainArg("--renumber_sites"):
        minx = 1
        maxx = sites.__len__()
        
#
#    depricated code to log-transform the plots
#
#    for s in sites:
#        yval = 0
#        if data[s] != None:
#            if data[s] < 0.00000001:
#                yval = 0.0
#            else:
#                yval = math.log( float(data[s]) + 1.0 )

    if image_type == "pdf":
        cranout.write("pdf('" + outpath + ".pdf', width=8, height=3.5);\n")
    if image_type == "png":
        cranout.write("png('" + outpath + ".png', width=800, height=300);\n")
    cranout.write("plot(c(" + minx.__str__() + "," + maxx.__str__() + "), c(-1,1), type='n',xlab='sites in " + ap.params["seed"] + "', ylab='normalized score', main='" + title + "', lwd=2, col='white');\n")
    metrics = data.keys()
    metrics.sort(reverse=False)
    for metric in metrics:
        if False == ap.doesContainArg("--renumber_sites"):
            this_data = fill_missing_sites(data[metric])
        else:
            this_data = data[metric]
        norm_data = normalize_vector(this_data)
        sites = norm_data.keys()
        sites.sort()

        lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
        x = "x" + metric + " <- c("
        y = "y" + metric + " <- c("
        for s in sites:
            if ap.doesContainArg("--renumber_sites"):
                x += int(ap.params["msa_mysite2seedsite"][lnick][s]).__str__() + ","
            else:
                x += s.__str__() + ","
            if norm_data[s] != None:
                if abs(norm_data[s]) < 0.00000001 and abs(norm_data[s]) > -0.00000001:
                    y += "0.0,"
                else:
                    y += norm_data[s].__str__() + ","
            else:
                y += "0.0,"
        x = re.sub(",$", "", x)
        x += ")"
        cranout.write( x + "\n")
        y = re.sub(",$", "", y)
        y += ")"
        cranout.write( y + "\n")
        cranout.write("points(x" + metric + ", y" + metric + ", lwd=" + METRIC_WT[metric].__str__() + ", type='l', col='" + METRIC_COLORS[metric] + "');\n")
    
    legx = "legx <- c("
    for metric in data:
        legx += "\"" + metric + "\","
    legx = re.sub(",$", "", legx)
    legx += ")"
    cranout.write(legx + "\n")
    legcol = "legcol <- c("
    for metric in data:
        legcol += "\"" + METRIC_COLORS[metric] + "\","
    legcol = re.sub(",$", "", legcol)
    legcol += ")"
    cranout.write(legcol + "\n")
    cranout.write("legend(\"topleft\", legx, pch = 20, col=legcol, title = \"Metrics\");\n")
    
    cranout.write("dev.off();\n")
    cranout.close()
    return cranpath
    
    
def blend_msa_data(msa_data):
    """This method sums the site scores from all the alignments into a single vector of site scores.
    This integration relies on the meta-alignment, previously generated in the method align_msas.
    INPUT: msa_data[msa nickname][site] = score
    OUTPUT: bdata[site] = summed score
    """
    bdata = {} # key = reference site from MSA-MSA, value = blended score
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    
    for ref_site in ap.params["msa_refsite2mysite"][ lnick ].keys():
        # Skip sites not in the restriction sites array.  (If no rsites settings were given by the user,
        # then all sites, by default, should be in the rsites array).
        if ref_site not in ap.params["rsites"][ lnick ]:
            continue
        
        # First check if this site exists in all the alignments.  If not, then skip this site.
        skip_this_site = False        
        for msa in msa_data:
            #print "651: msa = ", msa
            if ref_site not in ap.params["msa_refsite2mysite"][msa]:
                bdata[ref_site] = 0.0
                skip_this_site = True

        # Next, sum the scores from all the alignments.
        if skip_this_site == False:
           for msa in msa_data:
                mysite = ap.params["msa_refsite2mysite"][msa][ref_site]
                #  alg, mysite, ref_site
                myscore = 0.0
                if mysite != None:
                    myscore = msa_data[msa][mysite]
                if mysite != None:
                    if ref_site in bdata:            
                        bdata[ref_site] += (ap.params["msa_weights"][msa] * myscore) / msa_data.keys().__len__()
                    else:
                        bdata[ref_site] = ap.params["msa_weights"][msa] * myscore / msa_data.keys().__len__()
    return bdata

def write_change_summary_table(msa_changes):
    """Writes a table with a summary of how many states changed between ancestors."""
    for msanick in msa_changes:
        this_ancpath = ap.params["msa_comparisons"][msanick][0]
        that_ancpath = ap.params["msa_comparisons"][msanick][1]
        
        # [nsites, countindel, countred, countorange, countgreen]
        

def write_summary_table(data, msa_scores, metric_ranked):
    foutpath = get_table_outpath(ap, tag="summary.txt")
    print "\n. I'm writing a table with all scores to", foutpath
    fout = open(foutpath, "w")
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    refsites = ap.params["msa_refsite2mysite"][ lnick ].keys()
    refsites.sort()
    msapaths = msa_scores[ap.params["metrics"][0]].keys()
    msapaths.sort()
    
    metric_site_rank = {}
    for metric in ap.params["metrics"]:
        metric_site_rank[metric] = {}
        rank = 1
        for tuple in metric_ranked[metric]:
            site = tuple[0]
            metric_site_rank[metric][site] = rank
            rank += 1
    
    header = ""
    header += "site\t"
    for metric in ap.params["metrics"]:
        header += metric + "\t"
    for metric in ap.params["metrics"]:
        header += metric + "_rank\t"
    if msapaths.__len__() > 1:
        for m in msapaths:
            for metric in ap.params["metrics"]:
                header += metric + ":" + m + "\t"
    header += "\n"
    fout.write(header)

    line = ""
    for s in refsites:
        if s in ap.params["rsites"][ lnick ]:
            #print "965:", s 
            if ap.doesContainArg("--renumber_sites"):
                line = ap.params["msa_mysite2seedsite"][lnick][s].__str__() + "\t"
            else:
                line = s.__str__() + "\t"
            for metric in ap.params["metrics"]:
                line += "%.3f"%data[metric][s] + "\t"
            for metric in ap.params["metrics"]:
                line += metric_site_rank[metric][s].__str__() + "\t"
            if msapaths.__len__() > 1:
                for m in msapaths:
                    for metric in ap.params["metrics"]:
                        if s in ap.params["msa_refsite2mysite"][m]:
                            mys = ap.params["msa_refsite2mysite"][m][s]
                            line += "%.3f"%msa_scores[metric][m][mys] + "\t"
                        else:
                            line += "---\t"
            line += "\n"
            fout.write(line)
    fout.close()

def rank_all(metric_data, metric_blendeddata):
    metric_ranked = {}
    for i in ap.params["metrics"]:
        print "\n. I'm ranking the data for", i, "from across all the alignments."
        metric_ranked[i] = rank_sites(metric_blendeddata[i], metric_data[i], method=i)
    return metric_ranked
        
def correlate_all(metric_ranked, metric_blendeddata):
    cranpaths = [] # a list of R scripts that will be executed.
    if False == ap.doesContainArg("--skip_correlation"):
        if ap.params["metrics"].__len__() > 1:
            metric_comparisons = []
            for i in range(0, ap.params["metrics"].__len__()):
                for j in range(i, ap.params["metrics"].__len__()):
                    if i != j:
                        metric_comparisons.append( (ap.params["metrics"][i], ap.params["metrics"][j]) )
            for comparison in metric_comparisons:
                this_metric = comparison[0]
                that_metric = comparison[1]
                for cranpath in correlate_metrics(metric_ranked[this_metric], metric_ranked[that_metric], metric_blendeddata[this_metric], metric_blendeddata[that_metric], this_metric + "-" + that_metric, xunit=this_metric, yunit=that_metric):
                    cranpaths.append( cranpath )
    return cranpaths

def rank_sites(blended_data, msa_scores, method="h", writetable=True):
    """This method takes the site-value data and orders the sites in ascending order.
    Data is written to the *.ranked.txt file."""
    
    score_refsites = {} # key = score, value = array of ref_sites with this score
    sites = blended_data.keys()
    sites.sort()
    for s in sites:
        hval = "%.5f"%blended_data[s]
        hval = blended_data[s]
        if hval in score_refsites:
            score_refsites[hval].append(s)
        else:
            score_refsites[hval] = [ s ]
    
    hscores = score_refsites.keys()
    hscores.sort(reverse=True)

    msa_thatanc_lines = {}
    msa_thisanc_lines = {}
    for msa in ap.params["msa_comparisons"]:
        that_ancpath = ap.params["msa_comparisons"][msa][0]
        fin = open(that_ancpath, "r")
        msa_thatanc_lines[msa] = fin.readlines()
        fin.close()
        this_ancpath = ap.params["msa_comparisons"][msa][1]
        fin = open(this_ancpath, "r")
        msa_thisanc_lines[msa] = fin.readlines() 
        fin.close()
        
    if writetable == True:
        fout = open( get_table_outpath(ap, tag=method+".ranked.txt"), "w")
        outlines = []
        for i in range(0, hscores.__len__()):
            h = hscores[i]
            for refsite in score_refsites[h]:
                lout = "--> "
                lout += " " + method + " = %.3f"%h
                lout += " rank: " + (i+1).__str__()
                #print ap.params["msa_mysite2seedsite"]
                #exit()
                lout += " seed site: " + ap.params["msa_mysite2seedsite"][ ap.params["msa_path2nick"][ ap.params["longest_msa"] ] ][refsite].__str__()              
                lout += "\n" 
                header_lout = ""
                header_lout += "Site: " + ap.params["msa_mysite2seedsite"][ ap.params["msa_path2nick"][ ap.params["longest_msa"] ] ][refsite].__str__()
                header_lout += "\t" + method + "= " + "%.3f"%h
                header_lout += "\trank: " + (i+1).__str__()
                fout.write(header_lout + "\n")
                
                #print "1395:", refsite, ap.params["msa_seedseq"][ ap.params["longest_msa"] ].__len__()
                left = ap.params["msa_seedseq"][ ap.params["longest_msa"] ][refsite-2]
                here = ap.params["msa_seedseq"][ ap.params["longest_msa"] ][refsite-1]
                if refsite < ap.params["msa_seedseq"][ ap.params["longest_msa"] ].__len__():
                    right = ap.params["msa_seedseq"][ ap.params["longest_msa"] ][refsite]
                else:
                    right = "-"
                lout += "\tseed context: " + left + here + right + "\n"
                for msa in ap.params["msa_comparisons"]:
                    if refsite in ap.params["msa_refsite2mysite"][msa]:
                        mys = ap.params["msa_refsite2mysite"][msa][refsite]
                        lout += msa + " site " + mys.__str__() + " (" + method + " = %.3f"%msa_scores[msa][mys] + ")\n"
                        for l in msa_thatanc_lines[msa]:
                            if l.startswith(mys.__str__() + " "):
                                lout += "   " + datpath_to_short(ap.params["msa_comparisons"][msa][0]) + "\tsite " + re.sub("site ", "", l)
                        for l in msa_thisanc_lines[msa]:
                            if l.startswith(mys.__str__() + " "):
                                lout += "   " + datpath_to_short(ap.params["msa_comparisons"][msa][1]) + "\tsite " + re.sub("site ", "", l)
                lout += "\n"
                outlines.append(lout)
        fout.close()
        
        fout = open( get_table_outpath(ap, tag=method+".details.txt"), "w")
        for lout in outlines:
            fout.write(lout)
        fout.close()

    ranked_sites = [] # array of tuples in order, (refsite, value for this site)
    for i in range(0, hscores.__len__()):
        h = hscores[i]
        for refsite in score_refsites[h]:
            ranked_sites.append( (refsite, h) )
    return ranked_sites

def get_bins(min, max, stride):
    # Returns the bin
    bins = []
    if min < 0:
        minbin = 0.0
        while minbin > min:
            minbin -= stride
    elif min >= 0:
        minbin = 0.0
    #print "minbin=", minbin, min
    if max > 0:
        maxbin = 0.0
        while maxbin < max:
            maxbin += stride
    #print "maxbin=", maxbin, max
    i = minbin
    while i < maxbin:
        if i > 0.0 and 0.0 not in bins:
            bins.append(0.0)
        bins.append(i)
        if i >= 0.0 and i < 0.0+stride and TINY not in bins:
            bins.append( TINY )
        i += stride
    return bins

def get_bin(value, bins):
    if value == 0.0:
        return 0.0
    ret = bins[0]
    for i in bins:
        if i <= value:
            ret = i
    return ret

def plot_histogram(metric_data, image_type="pdf"):
    cranpaths = []
    for metric in metric_data:
        maxx = None
        minx = None
        allvalues = []
        data = metric_data[metric]
        for site in data:
            value = data[site]
            allvalues.append(value)
            if maxx == None:
                maxx = value
            if minx == None:
                minx = value
            if maxx < value:
                maxx = value
            if minx > value:
                minx = value
        binwidth = 0.1
        if minx == None or maxx == None:
            print "\n. Something went wrong."
            print metric_data
            exit()
        if maxx - minx > 1000:
            binwidth = 50
        elif maxx - minx > 500:
            binwidth = 25
        elif maxx - minx > 100:
            binwidth = 10
        elif maxx - minx > 50:
            binwidth = 5
        elif maxx - minx > 10:
            binwidth = 0.5
        elif maxx - minx > 5:
            binwidth = 0.1
        elif maxx - minx > 1:
            binwidth = 0.01
        elif maxx - minx > 0.5:
            binwidth = 0.01
        elif maxx - minx > 0.1:
            binwidth = 0.005
        elif maxx - minx > 0.05:
            binwidth = 0.001
        elif maxx - minx > 0.01:
            binwidth = 0.0005
    
        if ap.doesContainArg("--force_bin_width"):
            binwidth = float( ap.getArg("--force_bin_width") )
            #print binwidth, maxx-minx
            #exit()
            
        #print "1111:"
        bins = get_bins(minx, maxx, binwidth)
        #print "1113"
        bin_count= {}
        for b in bins:
            bin_count[b] = 0.0
        
        data = metric_data[metric]
        n = data.__len__()
        for site in data:
            this_bin = get_bin( data[site], bins )
            #print "901:", data[site], this_bin
            if bin_count[this_bin] < 20:
                bin_count[this_bin] += 1.0
        
        """
        Write the R script.
        """        
        if image_type == "pdf":
            pdfpath = get_plot_outpath(ap, tag=metric + "-histogram.pdf")
            cranstr = "pdf(\"" + pdfpath + "\", width=5, height=4);\n"    
        elif image_type == "png":
            pngpath = get_plot_outpath(ap, tag=metric + "-histogram.png")
            cranstr = "png(\"" + pngpath + "\", width=500, height=300);\n" 
    
        #for metric in metric_data:
        cranstr += metric + " <- c("
        for bin in bins:
            if ap.doesContainArg("--plot_histogram_log_axis"):
                if bin_count[bin] == 0.0:
                    bin_count[bin] = 0.5 / metric_data[metric].__len__() 
            cranstr += bin_count[bin].__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ")\n"

        sdev = sd(allvalues)
        avg = mean(allvalues)
    
        cranstr += "bins <- c("
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"

        # New style, drawing each bar with rect()
        xmax = roundup(maxx, 1.0)
        xmin = roundup(minx, 1.0)
        cranstr += "plot(c(" + roundup(xmin,2.0).__str__() + "," + roundup(xmax, 2.0).__str__() + "), c(0,20), type='n',xlab='" + metric + " score', ylab='N sites', main='" + metric + " histogram " + time.asctime().__str__() + "',lwd=2, las=1, col='black', bty='n');\n"
        for bin in bins:
            colstr = "black"
            
            """ depricated.
            if bin < avg:
                if bin <= avg - 6*sdev:
                    colstr = "red1"
                elif bin <= avg - 4*sdev:
                    colstr = "orange1"
                elif bin <= avg - 2*sdev:
                    colstr = "deepskyblue2"
                else:
                    colstr = "black" 
            elif bin > avg:
                if bin >= avg + 6*sdev:
                    colstr = "red1"
                elif bin >= avg + 4*sdev:
                    colstr = "orange1"
                elif bin >= avg + 2*sdev:
                    colstr = "deepskyblue2"
                else:
                    colstr = "black" 
            """
            
            cranstr += "rect(" + bin.__str__() + ", 0.0, " + (bin+binwidth).__str__() + ", " + bin_count[bin].__str__() + ", col=\"" + colstr + "\", lwd=0.1);\n"

        # depricated
        #cranstr += "par(lwd = 0.1);\n"   
        #if ap.doesContainArg("--plot_histogram_log_axis"): # make the Y axis log, or not?
        #    cranstr += "barx = barplot(as.matrix("+ metric + "), xlab=\"" + metric + "\", ylab=\"N sites\", beside=TRUE, log='y', l col=bincolors, las=1, names.arg=bins);\n"
        #else:
        #    cranstr += "barx = barplot(as.matrix("+ metric + "), xlab=\"" + metric + "\", ylab=\"N sites\", beside=TRUE, ylim=range(0,20), col=bincolors, las=1, names.arg=bins);\n"
    
#        # draw lines for 0, +/- 2 stdev, and +/- 4 stdev

        cranstr += "abline(0,0, lwd=0.75, col='black');\n"
        cranstr += "abline(v=" + avg.__str__() + ", lwd=0.5, lty='dashed', col='gray50');\n"
        if roundup(xmax,2.0) > (avg+2*sdev):
            cranstr += "abline(v=" + (avg+2*sdev).__str__() + ", lwd=0.5, lty='dashed', col='deepskyblue2');\n"
            cranstr += "text(" + (avg+2*sdev).__str__() + ", 20, \"+2 sigma\", cex=0.5, col='deepskyblue2');\n"
        if roundup(xmin,2.0) < (avg-2*sdev):
            cranstr += "abline(v=" + (avg-2*sdev).__str__() + ", lwd=0.5, lty='dashed', col='deepskyblue2');\n"
            cranstr += "text(" + (avg-2*sdev).__str__() + ", 20, \"-2 sigma\", cex=0.5, col='deepskyblue2');\n"
        if roundup(xmax,2.0) > (avg+4*sdev):
            cranstr += "abline(v=" + (avg+4*sdev).__str__() + ", lwd=0.5, lty='dashed', col='orange1');\n"
            cranstr += "text(" + (avg+4*sdev).__str__() + ", 20, \"+4 sigma\", cex=0.5, col='orange1');\n"
        if roundup(xmin,2.0) < (avg-4*sdev):
            cranstr += "abline(v=" + (avg-4*sdev).__str__() + ", lwd=0.5, lty='dashed', col='orange1');\n"
            cranstr += "text(" + (avg-4*sdev).__str__() + ", 20, \"-4 sigma\", cex=0.5, col='orange1');\n"
        if roundup(xmax,2.0) > (avg+6*sdev):
            cranstr += "abline(v=" + (avg+6*sdev).__str__() + ", lwd=0.5, lty='dashed', col='red1');\n"
            cranstr += "text(" + (avg+6*sdev).__str__() + ", 20, \"+6 sigma\", cex=0.5, col='red1');\n"
        if roundup(xmin,2.0) < (avg-6*sdev):
            cranstr += "abline(v=" + (avg-6*sdev).__str__() + ", lwd=0.5, lty='dashed', col='red1');\n"
            cranstr += "text(" + (avg-6*sdev).__str__() + ", 20, \"-6 sigma\", cex=0.5, col='red1');\n"
    
        # draw text for stdev lines.
        cranstr += "text(" + (avg).__str__() + ", 20, \"mean\", cex=0.5, col='gray50');\n"

        cranstr += "dev.off();\n"
    
        cranpath = get_plot_outpath(ap, tag=metric + "-" + image_type + "-histogram.rscript")
        fout = open(cranpath, "w")
        fout.write( cranstr )
        fout.close()
        cranpaths.append(cranpath) 
    return cranpaths

#
# ma and mb are sorted lists in score order
#
def correlate_metrics(ma, mb, ma_site_val, mb_site_val, tag, xunit=None, yunit=None):    
    cranpaths = []
    
    print "entereted correlate_metrics"
    
    if  ma.__len__() != mb.__len__():
        print "\n. Hmm, something is wrong.  Anccomp_tools.py point 742."
        exit()
    ma_ranked = [] # A values in rank order
    mb_ranked = [] # V values in rank order
    rank_a_b = [] # array of triples (A val, B val, mark)
    count_skipped = 0
    hrank = 0
    for i in range(0, ma.__len__()):
        hsite = ma[i][0]
        maval = ma[i][1]
        
        psite = mb[i][0]
        mbval = mb[i][1]

        if maval == 0.0 or mbval == 0.0:
            #count_skipped += 1
            #print "skipping", i, maval, mbval
            continue

        hrank += 1

        prank = 0
        j = -1
        while j < mb.__len__():
            j += 1
            prank += 1
            if mb[j][0] == psite:
                j = mb.__len__()
            
        print maval, mbval, hrank, prank
        ma_ranked.append(maval)
        mb_ranked.append(mbval)
        if prank < 0:
            print prank
            print "\n. Something is wrong with prank.  Goodbye."
            exit()
        mark = 1
        if maval == 0.0 and mbval == 0.0:
            mark = 0
        rank_a_b.append( (hrank, prank, mark) )
    spearmans =  ss.spearmanr(ma_ranked, mb_ranked)
    path  = plot_correlation(rank_a_b, get_plot_outpath(ap, tag="corr-rank." + tag), tag + " rank correlation", "blue", xlab=xunit, ylab=yunit, memo=("Spearman: %.3f"%spearmans[0]), image_type="pdf")
    cranpaths.append( path )
    path  = plot_correlation(rank_a_b, get_plot_outpath(ap, tag="corr-rank." + tag), tag + " rank correlation", "blue", xlab=xunit, ylab=yunit, memo=("Spearman: %.3f"%spearmans[0]), image_type="png")
    cranpaths.append( path )

    fout = open( get_plot_outpath(ap, tag="corr-rank." + tag) + ".txt" , "w")
    fout.write( spearmans[0].__str__() + "\n")
    for v in rank_a_b:
        fout.write( v[0].__str__() + "\t" + v[1].__str__() + "\n")
    fout.close()
    
    value_a_b = []
    mavals = []
    mbvals = []
    for site in ma_site_val.keys():
        maval = ma_site_val[site]
        mbval = mb_site_val[site]
        mark = 1
        if maval == 0.0 and mbval == 0.0:
            mark = 0
        value_a_b.append( [ maval,mbval,mark ] )
        mavals.append(maval)
        mbvals.append(mbval)
    pearsons = ss.pearsonr(mavals, mbvals)
    pearsonsT = pearsons[0]
    pearsonsP = pearsons[1]
    path  = plot_correlation(value_a_b, get_plot_outpath(ap, tag="corr-value." + tag), tag + "value correlation", "gray50", xlab=xunit, ylab=yunit, memo=("Pearson: %.3f"%pearsonsT), image_type="pdf")
    cranpaths.append( path )
    path  = plot_correlation(value_a_b, get_plot_outpath(ap, tag="corr-value." + tag), tag + "value correlation", "gray50", xlab=xunit, ylab=yunit, memo=("Pearson: %.3f"%pearsonsT), image_type="png" )
    cranpaths.append( path )
    
    fout = open( get_plot_outpath(ap, tag="corr-value." + tag) + ".txt" , "w")
    fout.write( pearsonsT.__str__() + "\n")
    for v in value_a_b:
        fout.write( v[0].__str__() + "\t" + v[1].__str__() + "\n")
    fout.close()
    
    return cranpaths
    #return [spearmans, pearsons]
    
def plot_correlation(data, outpath, title, color, xlab=None, ylab=None, memo=None, image_type="pdf"):
    cranpath = outpath + "." + image_type + ".rscript"
    cranout = open(cranpath, "w")
    miny = 0.0
    maxy = 0.0
    minx = 0
    maxx = 0
    x = "x <- c("
    y = "y <- c("
    marks = "m <- c("
    colors = "colors <- c("
        
    for tuple in data:
        #print s
        x += tuple[0].__str__() + ","
        y += tuple[1].__str__() + ","
        if tuple[2] == 0:
            marks += "1,"
            colors += "'grey',"
        elif tuple[2] == 1:
            marks += "20,"
            colors += "'" + color + "',"
        
        if float(tuple[0]) > maxx:
            maxx = float(tuple[0])
        if float(tuple[0]) < minx:
            minx = float(tuple[0])

        if float(tuple[1]) > maxy:
            maxy = float(tuple[1])
        if float(tuple[1]) < miny:
            miny = float(tuple[1])

    x = re.sub(",$", "", x)
    x += ")"
    cranout.write( x + "\n")
    y = re.sub(",$", "", y)
    y += ")"
    cranout.write( y + "\n")
    marks = re.sub(",$", "", marks)
    marks += ")"
    cranout.write( marks + "\n")
    colors = re.sub(",$", "", colors)
    colors += ")"
    cranout.write( colors + "\n")
    
    if image_type == "pdf":
        cranout.write("pdf('" + outpath + ".pdf', width=6, height=6);\n")
    elif image_type == "png":
        cranout.write("png('" + outpath + ".png', width=500);\n")
    cranout.write("plot(c(" + minx.__str__() + "," + maxx.__str__() + "), c(" + miny.__str__() + "," + maxy.__str__() + "), type='n',xlab=\"" + xlab + "\", ylab=\"" + ylab + "\", main='" + title + "', lwd=2, col='" + color + "', bty='n');\n")
    if memo != None:
        textx = minx + 0.2*(maxx-minx)
        cranout.write( "text(" + textx.__str__() + "," + maxy.__str__() + ",\"" + memo + "\");\n")
    cranout.write("points(x, y, lwd=2, type='p', col=colors, pch=m);\n")
    cranout.write("dev.off();\n")
    cranout.close()
    return cranpath    

def smooth_data(metric_blendeddata):
    w_metric_blendeddata = {} # key = window size for smoothing, value = hashtable, where key = metric ID, value = blended data
    cranpaths = []
    for w in ap.params["winsizes"]:    
        w_metric_blendeddata[w] = {}
        for metric in ap.params["metrics"]:
            w_metric_blendeddata[w][metric] = {}
            if w == 1:
                w_metric_blendeddata[w][metric] = metric_blendeddata[metric] # pass through
            else:
                w_metric_blendeddata[w][metric] = window_analysis(metric_blendeddata[metric], w, ap)
            
        if w == 1:
            # Plot the histogram of values. . . .
            for cranpath in plot_histogram(w_metric_blendeddata[w], image_type="pdf"):
                cranpaths.append( cranpath )
            for cranpath in plot_histogram(w_metric_blendeddata[w], image_type="png"):
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
    
        #
        # Plot the site-by-metric for each metric
        #
        for metric in ap.params["metrics"]:
            # Plot the data for each site. . .
            plot_outpath = get_plot_outpath(ap, tag=(metric + "-by-site.w=" + w.__str__()) )
            combo_substring = ""
            if False != ap.getOptionalArg("--combo_method"):
                combo_substring = ", " + ap.getOptionalArg("--combo_method")
            plot_title = metric + " score, winsize = " + w.__str__() + combo_substring + ", " + time.asctime().__str__() + ""
            cranpath = plot_one_metric( w_metric_blendeddata[w][metric], plot_outpath, plot_title, metric + " score", "black", image_type="pdf")
            cranpaths.append(cranpath)  
            cranpath = plot_one_metric( w_metric_blendeddata[w][metric], plot_outpath, plot_title, metric + " score", "black", image_type="png")
            cranpaths.append(cranpath) 
        
        #
        # Plot the site-by-metric, normalized, for all metrics on one plot.
        #
        multi_plot_outpath = get_plot_outpath(ap, tag=("all-by-site.w=" + w.__str__()) )
        plot_title = "All metrics, normalized, winsize = " + w.__str__() + combo_substring + ", " + time.asctime().__str__() + ""
        cranpath = plot_multi_metrics(w_metric_blendeddata[w], multi_plot_outpath, plot_title,image_type="pdf")
        cranpaths.append(cranpath)
        cranpath = plot_multi_metrics(w_metric_blendeddata[w], multi_plot_outpath, plot_title, image_type="png")
        cranpaths.append(cranpath)
    return cranpaths




