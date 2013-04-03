import math, os, re, sys, time
from config import *
from splash import *
import scipy.stats as ss
from argparser import *
ap = ArgParser(sys.argv)

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


def get_anccomp_dir(ap):
    return ap.getArg("--runid")

def get_plot_outpath(ap, tag=""):
    return get_anccomp_dir(ap) + "/" + tag 

def get_table_outpath(ap, tag=""):
    return get_anccomp_dir(ap) + "/" + tag 

def datpath_to_short(msapath):
    tokens = msapath.split("/")
    path = tokens[ tokens.__len__()-1 ]
    path = re.sub(".dat", "", path)
    return path

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
            h['-'] = 1.0
            return h
        h[ tokens[i].upper() ] = float(tokens[i+1])
        i += 2
    h = fill_missing_states(h)
    return h

#def get_matrix(filepath):
#    """Matrix files should contain a 19x19 matrix, similar to the exemplar file jtt.txt"""
#    fin = open(filepath, "r")
#    m = {}
#    for c in AA_ALPHABET:
#        m[c] = {}
#    lines = fin.readlines()
#    curr_line = 0
#    for l in lines:
#        if l.__len__() > 2:
#            tokens = l.split()
#            for i in range(0, tokens.__len__()):
#                m[ AA_ALPHABET[curr_line] ][ AA_ALPHABET[i] ] = float(tokens[i])
#        curr_line += 1
#    fin.close()
#    return m

def get_matrix(filepath):
    """Matrix files should contain a 19x19 matrix"""
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
    return norm_m

def fill_matrix_rev_values(m):
    for c in AA_ALPHABET:
        for d in AA_ALPHABET:
            if d not in m[c].keys():
                m[c][d] = float(m[d][c])
    return m    

def read_specs(p):
    ap.params["msa_comparisons"] = {} # key = msa nickname, value = array of *.dat tuples
    ap.params["msa_weights"] = {}
    seed = ""
    ap.params["msa_path2nick"] = {}
    ap.params["msa_nick2path"] = {}
    fin = open(p, "r")
    for l in fin.readlines():
        if l.startswith("#"): # skip lines with comments
            continue
        if l.startswith("msapaths"):
            tokens = l.split()
            for t in tokens[1:]:
                ap.params["msa_path2nick"][t] = t
                #ap.params["msa_comparisons"][t] = None
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
            ap.params["msa_path2nick"][ tokens[1] ] = tokens[2]
            ap.params["msa_nick2path"][ tokens[2] ] = tokens[1]
    fin.close()    
    
    check_specs()
    
def check_specs():
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

    ap.params["anc_msasource"] = {}
    for msapath in ap.params["msa_path2nick"]:
        msanick = ap.params["msa_path2nick"][msapath]
        this_ancpath = ap.params["msa_comparisons"][msanick][0]
        that_ancpath = ap.params["msa_comparisons"][msanick][1]
        ap.params["anc_msasource"][ this_ancpath ] = msanick
        ap.params["anc_msasource"][ that_ancpath ] = msanick


def get_msa_len(p):
    #print path
    #os.system("cat " + path)
    #print p
    fin = open(p, "r")
    lines = fin.readlines()
    fin.close()
    return float(lines[0].strip().split()[1])

def get_seed_seq(msapath, seed):
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
        print "\n. MSA [", p, "] has", int(len), "sites."
    ap.params["longest_msa"] = longest_msa
    
    print "\n. I'm aligning the alignments, using taxa", seed, "as the seed."
    
    # Find the seed sequence in each MSA. . .
    msa_seedseq = {}
    for p in msapaths:
        msa_seedseq[p] = get_seed_seq(p, seed)
        
#        fin = open(p, "r")
#        for l in fin.readlines():
#            if l.__len__() > 1:
#                if l.split()[0] == seed:
#                    tokens = l.split()
#                    msa_seedseq[p] = tokens[1]
#        fin.close()
#        if msa_seedseq[p].__len__() == 0:
#            print "\n. Hmmm, I couldn't find a seed sequence for ", seed, "in ", p
        
    # verify that the focus sequences are actually the SAME sequence 
    for p in msapaths:
        if re.sub("\-", "", msa_seedseq[p]).__len__() != re.sub("\-", "", msa_seedseq[longest_msa]).__len__():
            print re.sub("\-", "", msa_seedseq[p])
            print re.sub("\-", "", msa_seedseq[longest_msa])
            print "\n. Hmmm, I had to stop because the seed sequence is not the same in all your alignments."
            print ". Tip: Check if you used a truncated version of the seed sequence in some of the alignments."
            print ""

    # Align the focus sequences. . .
    ap.params["msa_refsite2mysite"] = {} # key = MSA nickname, value = hash, where key = site in longest alg, value = my site (or None)
    ap.params["msa_mysite2refsite"] = {} # the reverse lookup of ap.params["msa_refsite2mysite"]
    for p in msapaths:
        print "\n. Aligning [", p, "] to [", longest_msa, "]"
        nick = ap.params["msa_path2nick"][p]
        ap.params["msa_refsite2mysite"][nick] = {}
        ap.params["msa_mysite2refsite"][nick] = {}
        ref_site = 1
        my_site = 1
        while my_site-1 < msa_seedseq[p].__len__() and ref_site-1 < msa_seedseq[longest_msa].__len__():
            if msa_seedseq[p][my_site-1] == msa_seedseq[longest_msa][ref_site-1]:
                my_state = msa_seedseq[p][my_site-1]
                ref_state = msa_seedseq[longest_msa][ref_site-1]
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
            elif msa_seedseq[longest_msa][ref_site-1] == "-":
                while msa_seedseq[longest_msa][ref_site-1] == "-" and ref_site < msa_seedseq[longest_msa].__len__()+1:
                    ref_site += 1
            elif msa_seedseq[p][my_site-1] == "-":
                while msa_seedseq[p][my_site-1] == "-" and my_site < msa_seedseq[p].__len__()+1:
                    my_site += 1
    
    # Write a log file about the meta-alignment,
    # and fill the hash ap.params["msa_refsite2mysite"]
    msa_ids = {} # key = MSA path, value = integer ID
    ids_msa = {} # reverse hash of above
    msa_ids[longest_msa] = 1
    ids_msa[1] = longest_msa
    counter = 2
    for p in msapaths:
        if p != longest_msa:
            msa_ids[ p ] = counter
            ids_msa[ counter ] = p
            counter += 1
    fout = open(get_anccomp_dir(ap) + "/meta_alignment.txt", "w")
    ids = ids_msa.keys()
    ids.sort()
    fout.write("Alignment Key:\n")
    for i in ids:
        fout.write( "M" + i.__str__() + " : " + ap.params["msa_path2nick"][ ids_msa[i] ]  + "\n")
    fout.write("\n")
    header = ""
    fout.write("Residues shown in parentheses express the state of taxon " + ap.params["seed"] + " at the corresponding site.\n\n")
    fout.write("The mark 'x' indicates that no corresponding site was found in the alignment.\n")
    for i in ids:
        header += "M" + i.__str__() + "\t"
    fout.write(header + "\n\n")
    for site in range(0, msa_seedseq[longest_msa].__len__()):
        line = ""
        for i in ids:
            msanick = ap.params["msa_path2nick"][ids_msa[i]]
            if (site+1) in ap.params["msa_refsite2mysite"][ msanick ]:
                mysite = ap.params["msa_refsite2mysite"][ msanick ][ site+1 ]
                state = msa_seedseq[ids_msa[i]][mysite-1]
                line += mysite.__str__() + " (" + state + ")\t"
            else:
                line += "x\t"
        fout.write( line + "\n" )
    fout.close()

    
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


def build_rsites():
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    
    ap.params["rsites"] = {}
    for msa in ap.params["msa_nick2path"]:
        ap.params["rsites"][msa] = []
    
    
    if (ap.doesContainArg("--limstart") or ap.doesContainArg("--limstop")) and ap.doesContainArg("--restrict_sites"):
        print "\n. Sorry, I'm confused.  You cannot specific --restrict_sites along with --limstart or --limstop."
        print ". I'm quitting."
        exit()
    
    limstart = 1
    limstop = ap.params["msa_refsite2mysite"][ lnick ].keys().__len__() + 1
    
    # Build the restriction site library, using the site numbers in the longest MSA. . .
    x = ap.getOptionalList("--restrict_sites")
    y = ap.getOptionalArg("--restrict_to_seed")
    if x != None:
        for i in x:
            ap.params["rsites"][lnick].append(int(i))
            #print i
    elif y != None:
        print "\n. I'm restricting the analysis to only those sites found in the seed sequence", ap.params["seed"]
        seed_seq = get_seed_seq(  ap.params["msa_nick2path"][lnick], ap.params["seed"]  )
        for site in range(0, limstop-1):
            if seed_seq[site] != "-":
                ap.params["rsites"][lnick].append(site)
                print site
    else:
        x = ap.getOptionalArg("--limstart")
        if x != False:
            limstart = int(x)
        y = ap.getOptionalArg("--limstop")
        if y != False:
            limstop = int(y)
        for i in range(limstart, limstop+1):
            ap.params["rsites"][lnick].append(i)
    
    # Map the restriction sites onto the other MSAs. . . 
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

    if ap.doesContainArg("--renumber_sites"):
        ref2seed = {} # key = reference site, value = corresponding seed sequence site
        count = 0
        for site in ap.params["rsites"][lnick]:
            count += 1
            ref2seed[site] = count
    
        ap.params["msa_mysite2seedsite"] = {}
        for msanick in ap.params["msa_nick2path"]:
            ap.params["msa_mysite2seedsite"][msanick] = {}
            for mysite in ap.params["rsites"][msanick]:
                ap.params["msa_mysite2seedsite"][msanick][mysite] = ref2seed[ ap.params["msa_mysite2refsite"][msanick][mysite] ]
                #print msanick, mysite, ref2seed[ ap.params["msa_mysite2refsite"][msanick][mysite] ]
    

def compare_dat_files(patha, pathb, m, winsize, method="h"):
    """Compares *.dat files for two ancestors from the same alignment."""
    """Returns windata and consdata"""
    """windata[site] = h score"""
    """consdata[site] = r score"""

    msaone = ap.params["anc_msasource"][patha]
    
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

        if asite in ap.params["rsites"][msaone]:
            ancx = datline_2_pphash(al)
            ancy = datline_2_pphash(bl)       
            results[asite] = d(ancx, ancy, m, method=method)  
    return results

def dat2pp(datpath):
    data = {}
    fin = open(datpath, "r")
    linesa = fin.readlines()
    fin.close()
    
    for li in range(0, linesa.__len__()):    
        l = linesa[li]  
        site = int(l.split()[0])
        data[site] = datline_2_pphash(l)
    return data


def compare_dat_files_stats(datone, dattwo):    
    anc_data = {}
    anc_data[datone] = dat2pp(datone)
    anc_data[dattwo] = dat2pp(dattwo)
    
    msaone = ap.params["anc_msasource"][datone]
    msatwo = ap.params["anc_msasource"][dattwo]
    
    # compare dat files:
    nsites = 0
    countindel = 0
    countred = 0 # different states
    countorange = 0 # diff. states, but pairwise secondary states
    countgreen = 0 # same state, but adds uncertainty
    redsites = []
    orangesites = []
    greensites = []
    
    sites = ap.params["msa_refsite2mysite"][ ap.params["longest_msa"] ].keys()
    
    #print "\n"
    #print ap.params["msa_mysite2refsite"]
    #print sites
    for site in sites:
        #
        # to-do: filter for rsites
        #    rsites = ap.params["rsites"][ap.params["anc_msasource"][datpath]]
        #    longest_msa = ap.params["longest_msa"]
        
        if site not in ap.params["msa_refsite2mysite"][msaone].keys() and site not in ap.params["msa_refsite2mysite"][msatwo].keys():
            continue
        else:
            nsites += 1
        
        indel = False
        red = False
        orange = False
        green = False
        
        if site not in ap.params["msa_refsite2mysite"][msaone] or site not in ap.params["msa_refsite2mysite"][msatwo]:
            indel = True
            countindel += 1
            continue
        
        siteone = ap.params["msa_refsite2mysite"][msaone][site]
        sitetwo = ap.params["msa_refsite2mysite"][msatwo][site]
        a_state = getmlstate( anc_data[datone][siteone] )
        b_state = getmlstate( anc_data[dattwo][sitetwo] )
        if b_state == None:
            indel = True
        # test for indel mismatch:
        elif b_state != a_state:
            if b_state == "-" or a_state == "-":
                indel = True
        # test for red:
        elif b_state != a_state and indel == False:
            if False == anc_data[dattwo][sitetwo].keys().__contains__(a_state):
                red = True
            elif False == anc_data[ datone ][siteone].keys().__contains__(b_state):
                red = True
            elif anc_data[datone][siteone][b_state] < 0.05 and anc_data[dattwo][sitetwo][a_state] < 0.05:
                red = True
        # test for orange:
        elif b_state != a_state and red == False and indel == False:
            orange = True
        # test for green
        elif b_state == a_state:
            #print "372:",  anc, site, b_state
            if (anc_data[dattwo][sitetwo][b_state] < 0.8 and anc_data[datone][siteone][b_state] > 0.8) or (anc_data[dattwo][sitetwo][b_state] > 0.8 and anc_data[datone][siteone][b_state] < 0.8):
                green = True
        if indel:
            countindel += 1
        if red:
            countred += 1
            redsites.append(site)
            #print "site", site, " - case 1"
            #print "\t",getmlstate( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
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
    
    return [nsites, countindel, countred, countorange, countgreen]
    
#    print "============================================================"
#    print "Legend:"
#    print "- case 1: sites with disagreeing ML states, and neither ancestral vector has strong support for the others' state."
#    print "- case 2: sites with disagreeing ML states, but one or both of the ancestral vectors has PP < 0.05 for the others' ML state."
#    print "- case 3: sites with the same ML state, but one of the ancestral vectors strongly supports (PP > 0.8) the state, while the other vector poorly supports (PP < 0.8) the state."
#    print "============================================================"
#    print "\n"
#    
#    print "============================================================"
#    print "Summary:"
#    print "nsites=", nsites
#    print "indel_mismatch=", countindel
#    print "case 1=", countred, "sites: ",  redsites
#    print "case 2=", countorange, "sites: ",  orangesites
#    print "case 3=", countgreen, "sites: ", greensites
#    print "============================================================"
#    print "\n\n"


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
     
        
def conditional_entropy(ppx, ppy):
    """Retuns H(ppy|ppx)"""
    """Conditional entropy = h(y|x) = h(x,y) - h(x)"""
    ret = 0.0
    for a in AA_ALPHABET:
        for b in AA_ALPHABET:
            jointp = ppx[a]*ppy[b]
            if jointp > 0:
                ret +=  jointp * math.log( jointp )
    ret = ret - entropy(ppx)
    return -1 * ret


def information_gain(ppx, ppy, m):
    ret = 0.0
    ret = entropy(ppy) - conditional_entropy(ppx, ppy)
    return ret

def geometric_distance(i, j, ppi, ppj, m):
    return abs(ppi - ppj)

def kl_divergence(ppx, ppy, m):
    # this first part is classic KL...
    ret = 0.0    
    
    yvals = {}
    xvals = {}
    
    # Set low values to be uncertain (i.e. == 0.05)
    for a in AA_ALPHABET:
        if ppy[a] < 0.01:
            yvals[a] = 0.05
        else:
            yvals[a] = ppy[a]
        if ppx[a] < 0.01:
            xvals[a] = 0.05
        else:
            xvals[a] = ppx[a]
    
    # Normalize. . . .
    ysum = 0
    for a in AA_ALPHABET:
        ysum += yvals[a]
    xsum = 0
    for a in AA_ALPHABET:
        xsum += xvals[a]
    for a in AA_ALPHABET:
        yvals[a] = yvals[a]/ysum
        xvals[a] = xvals[a]/xsum
    
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
            
def hdist(ppx, ppy, m, method):
    # 1. Get reciprocal KL distance. . .     
    klxy = kl_divergence(ppx, ppy, m)
    klyx = kl_divergence(ppy, ppx, m)
    klsum = klxy + klyx
            
    if method == "h":
    # 2. Calculate the importance of the a.a. changes. . .
        s = 0.0
        for a in AA_ALPHABET:
            for b in AA_ALPHABET:
                if ppx[a] > 1.0/(AA_ALPHABET.__len__()) and ppy[b] > ppx[b] and a!=b:
                    s_part = (1.0/m[a][b])
                    s_part *= ppx[a] * ppy[b]
                    s += s_part

    elif method == "hb":
        # 3. Calculate observed/expected probability shifts:
        expected_p = 0.0
        observed_p = 0.0
        for a in AA_ALPHABET:
            for b in AA_ALPHABET:
                if a != b:
                    expected_p += ppx[a] * math.exp( m[a][b] )
                    observed_p += ppx[a] * ppy[b]
        s = observed_p / expected_p
    
    # 4. Adjust the sign +/- based on the direction of entropy change
    enty = entropy(ppy)
    entx = entropy(ppx)
    e_dir = 1
    if (entx < enty) and ppy[ getmlstate(ppy) ] < CERTAINTY_CUTOFF:
    #if (entx - enty) <= ENTROPY_CUTOFF:
        e_dir *= -1

    # 5. Combine evidence from 1, 2, and 3
    # Multiple by 100 to make the numbers interpretable
    combostyle = ap.getOptionalArg("--combo_method")
    if combostyle == False:
        combostyle = "product"
    if combostyle == "sum":
        h = klsum + (1-w)*s + e_dir
    elif combostyle == "product":
        h = klsum * s * e_dir
    return h


def pdist(ppx, ppy, m):
    expected_p = 0.0
    observed_p = 0.0
    for a in AA_ALPHABET:
        for b in AA_ALPHABET:
            if a != b:
                expected_p += ppx[a] * math.exp( m[a][b] )
                observed_p += ppx[a] * ppy[b]
    #print "421:", ppx, ppy, expected_p, observed_p
    pval = observed_p / expected_p

    enty = entropy(ppy)
    entx = entropy(ppx)
    e_dir = 1
    if (entx < enty) and ppy[ getmlstate(ppy) ] < CERTAINTY_CUTOFF:
    #if (entx - enty) <= ENTROPY_CUTOFF:
        e_dir *= -1

    return e_dir * pval * 100

def mutual_information(ppx, ppy):
    ret = 0.0
    for a in AA_ALPHABET:
        if ppx[a] > 0.0:
            for b in AA_ALPHABET:
                if ppy[b] > 0.0:
                    ret += ppx[a]*ppy[a]*math.log()

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

def method8(ancx, ancy, m):
    expect = {}
    for a in AA_ALPHABET:
        for b in AA_ALPHABET:
            if b not in expect:
                expect[b] = 0.0
            expect[b] += ancx[a] * m[a][b]
    # normalize expect:
    sum = 0.0
    for a in expect:
        sum += expect[a]
    for a in expect:
        expect[a] = expect[a] / sum

    # compare expect to observation:
    ret = 0.0
    for a in AA_ALPHABET:
        ret += abs( expect[a] - ancy[a] )
    return ret

def d(ancx, ancy, m, method="h"):    
    """First, deal with gaps and X (uncertain) characters...."""
    if "-" in ancx.keys() or "X" in ancx.keys():
        ancx = {}
        for aa in AA_ALPHABET:
            ancx[aa] = 0.05
    if "-" in ancy.keys() or "X" in ancy.keys():
        ancy = {}
        for aa in AA_ALPHABET:
            ancy[aa] = 0.05
    """Second, compare the site..."""
    if method == "h" or method == "hb":
        return hdist(ancx, ancy, m, method=method)        
    elif method == "p":
        return pdist(ancx, ancy, m)
            
def cons_delta(ancx, ancy):
    entx = 0.0
    enty = 0.0
    if "-" not in ancx.keys() and "X" not in ancx.keys():
        entx = entropy(ancx)
    if "-" not in ancy.keys() and "X" not in ancy.keys():
        enty = entropy(ancy)
    #print entx, enty, (entx-enty)
    return (entx - enty)

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
    min_s = data[data.keys()[0]]
    max_s = min_s
    for s in data.keys():
        if s < min_s:
            min_s = s
        if s > max_s:
            max_s = s
    for i in range(1, max_s+1):
        if i in data:
            new_data[i] = data[i]
        else:
            new_data[i] = None
    return new_data

def plot(data, outpath, title, ylab, color):
    #data = fill_missing_sites(data)
    
    cranpath = outpath + ".cran"
    cranout = open(cranpath, "w")
    x = "x <- c("
    y = "y <- c("
    sites = data.keys()
    sites.sort()

    miny = 0.0
    maxy = 0.0
    miny = 0.0
    maxy = 0.0
    
    minx = sites[0]
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
                yval = math.log( float(data[s]) + 1.0 )
    
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
        
    cranout.write("pdf('" + outpath + ".pdf', width=8, height=3.5);\n")
    cranout.write("plot(c(" + minx.__str__() + "," + maxx.__str__() + "), c(" + miny.__str__() + "," + maxy.__str__() + "), type='n',xlab='sites', ylab='" + ylab + "', main='" + title + "', lwd=2, col='" + color + "');\n")
    cranout.write("points(x, y, lwd=2, type='l', col='" + color + "');\n")
    cranout.write("dev.off();\n")
    cranout.close()
    return cranpath

def blend_msa_data(msa_data):
    """This method sums the site scores from all the alignments into a single vector of site scores.
    INPUT: msa_data[msa algorithm][site] = score
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
                        bdata[ref_site] += (ap.params["msa_weights"][msa] * myscore)
                    else:
                        bdata[ref_site] = ap.params["msa_weights"][msa] * myscore
    return bdata

def write_table(hdata, msa_scores, method="h"):
    foutpath = get_table_outpath(ap, tag=method + ".summary.txt")
    fout = open(foutpath, "w")
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    refsites = ap.params["msa_refsite2mysite"][ lnick ].keys()
    refsites.sort()
    msapaths = msa_scores.keys()
    msapaths.sort()
    header = ""
    header += "site\t"
    header += method + "\t"
    if msapaths.__len__() > 1:
        for m in msapaths:    
            header += method + "(" + m + ")\t"
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
            line += "%.3f"%hdata[s] + "\t"
            if msapaths.__len__() > 1:
                for m in msapaths:
                    if s in ap.params["msa_refsite2mysite"][m]:
                        mys = ap.params["msa_refsite2mysite"][m][s]
                        line += "%.3f"%msa_scores[m][mys] + "\t"
                    else:
                        line += "\t"
            line += "\n"
            fout.write(line)
    fout.close()

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
        for i in range(0, hscores.__len__()):
            h = hscores[i]
            for refsite in score_refsites[h]:
                lout = "--> "
                lout += " " + method + " = %.3f"%h
                lout += " rank: " + (i+1).__str__()
                lout += "\n" 
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
                fout.write(lout)
        fout.close()

    ranked_sites = [] # array of tuples, (refsite, value for this site)
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

def plot_histogram(metric_data, tag):
    cranpaths = []
    metric_bin_count = {}
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
        metric_bin_count[metric] = {}
        for b in bins:
            metric_bin_count[metric][b] = 0.0
        
        data = metric_data[metric]
        n = data.__len__()
        for site in data:
            this_bin = get_bin( data[site], bins )
            #print "901:", data[site], this_bin
            metric_bin_count[metric][this_bin] += 1.0/n
        
        """
        Write the R script.
        """        
        pdfpath = tag + "." + metric + ".pdf"
        cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    
        #for metric in metric_data:
        cranstr += metric + " <- c("
        for bin in bins:
            if metric_bin_count[metric][bin] == 0.0:
                metric_bin_count[metric][bin] = 0.5 / metric_data[metric].__len__() 
            cranstr += metric_bin_count[metric][bin].__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ")\n"
    
        cranstr += "bins <- c("
        for bin in bins:
            cranstr += bin.__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"
    
        cranstr += "barx = barplot(as.matrix("+ metric + "), xlab=\"" + metric + "\", beside=TRUE, log='y', col=c(\"blue\"), names.arg=bins);\n"
    
#        avg = mean(allvalues)
#        stdev = sd(allvalues)
#        cranstr += "abline(v=" + avg.__str__() + ", lty=2,lwd=2)\n"
#        cranstr += "abline(v=" + (avg-stdev).__str__() + ", lty=3,lwd=1)\n"
#        cranstr += "abline(v=" + (avg+stdev).__str__() + ", lty=3,lwd=1)\n"
#        cranstr += "abline(v=" + (avg-(2*stdev)).__str__() + ", lty=3,lwd=1)\n"
#        cranstr += "abline(v=" + (avg+(2*stdev)).__str__() + ", lty=3,lwd=1)\n"
    
        cranpath = tag + "." + metric + ".cran"
        fout = open(cranpath, "w")
        fout.write( cranstr )
        fout.close()
        cranpaths.append(cranpath) 
    return cranpaths

def correlate_metrics(ma, mb, ma_site_val, mb_site_val, tag):    
    cranpaths = []
    
    if  ma.__len__() != mb.__len__():
        print "\n. Hmm, something is wrong.  Anccomp_tools.py point 742."
        exit()
    ma_ranked = [] # A values in rank order
    mb_ranked = [] # V values in rank order
    rank_a_b = [] # array of triples (A val, B val, mark)
    count_skipped = 0
    for i in range(0, ma.__len__()):
        maval = ma[i][1]
        mbval = mb[i][1]
        hsite = ma[i][0]
        psite = mb[i][0]
        hrank = i+1-count_skipped
        prank = 0
        for j in range(0, mb.__len__()):
            if mb[j][0] == ma[i][0]:
                prank = j+1-count_skipped
        #print maval, mbval, hrank, prank
        ma_ranked.append(maval)
        mb_ranked.append(mbval)
        if prank < 0:
            print prank
            exit()
        mark = 1
        if maval == 0.0 and mbval == 0.0:
            mark = 0
        rank_a_b.append( (hrank, prank, mark) )
    spearmans =  ss.spearmanr(ma_ranked, mb_ranked)[0]
    path  = plot_correlation(rank_a_b, get_plot_outpath(ap, tag="corr-rank." + tag), tag + " rank correlation", "blue" )
    #os.system("r --no-save < " + path)
    cranpaths.append( path )

    fout = open( get_plot_outpath(ap, tag="corr-rank." + tag) + ".txt" , "w")
    fout.write( spearmans.__str__() + "\n")
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
    path  = plot_correlation(value_a_b, get_plot_outpath(ap, tag="corr-value." + tag), tag + "value correlation", "green" )
    #os.system("r --no-save < " + path)
    cranpaths.append( path )
    
    fout = open( get_plot_outpath(ap, tag="corr-value." + tag) + ".txt" , "w")
    fout.write( pearsons.__str__() + "\n")
    for v in value_a_b:
        fout.write( v[0].__str__() + "\t" + v[1].__str__() + "\n")
    fout.close()
    
    return cranpaths
    #return [spearmans, pearsons]
    
def plot_correlation(data, outpath, title, color):
    cranpath = outpath + ".cran"
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
        
    cranout.write("pdf('" + outpath + ".pdf', width=6, height=6);\n")
    cranout.write("plot(c(" + minx.__str__() + "," + maxx.__str__() + "), c(" + miny.__str__() + "," + maxy.__str__() + "), type='n',xlab='hvals', ylab='', main='" + title + "', lwd=2, col='" + color + "');\n")
    cranout.write("points(x, y, lwd=2, type='p', col=colors, pch=m);\n")
    cranout.write("dev.off();\n")
    cranout.close()
    return cranpath    


