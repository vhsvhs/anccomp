#
# USAGE:
# %> python compare_ancs.py --specpath <path1> --modelpath <path2> --window_sizes x y z
#
# path1 is the filepath to a specification file, formatted as shown below.
# path2 is the filepath to a text file containing a substitution matrix, such as JTT, WAG, LG, etc.
# "x y z" are window sizes for smoothing of the final data
#
# SPECIFICATION FILE uses the following format:
# msapaths <1> <2> <3> etc.
# seed <taxa name>
# compare <node a> <node b> <msapath>
# compare . . . .
#
# WINDOW SIZES
# The list of desired sizes can be arbitrarily long.
#
# OPTIONAL RUNTIME COMMANDS:
# --limstart X, analysis will be limited to sites X and beyond in the alignment
#

import math, os, re, sys, time
from splash import *
from argparser import *
AA_ALPHABET = ["A",  "R",    "N"  , "D" ,  "C"   , "Q"   , "E"  ,  "G"    ,"H"  ,  "I"  , "L"  ,  "K"  ,  "M"   , "F"  ,  "P"  ,  "S" ,   "T" ,  "W" ,  "Y" , "V"]
D_SCALAR = 1000
ap = ArgParser(sys.argv)
COLORS = ["blue", "red", "green", "purple", "orange"]

TINY = 0.001
ENTROPY_CUTOFF = -0.2 


def get_anccomp_dir(ap):
    return ap.getArg("--runid")

def get_plot_outpath(ap, tag=""):
    return get_anccomp_dir(ap) + "/" + tag 

def get_table_outpath(ap, tag=""):
    return get_anccomp_dir(ap) + "/" + tag 

def msapath_to_short(msapath):
    tokens = msapath.split("/")
    return tokens[ tokens.__len__()-1 ]

def fill_missing_states(h):
    for a in AA_ALPHABET:
        if a not in h:
            h[a] = TINY
    sum = 0.0
    # normalize to sum of 1.0. . .
    for a in AA_ALPHABET:
        sum += h[a]
    for a in AA_ALPHABET:
        h[a] = h[a] / sum        
    return h

def line_2_hash(line):
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
    #print h
    return h

def get_matrix(filepath):
    """Matrix files should contain a 20x20 matrix, similar to the exemplar file jtt.txt"""
    fin = open(filepath, "r")
    m = {}
    for c in AA_ALPHABET:
        m[c] = {}
    lines = fin.readlines()
    curr_line = 0
    for l in lines:
        if l.__len__() > 2:
            tokens = l.split()
            for i in range(0, tokens.__len__()):
                m[ AA_ALPHABET[curr_line] ][ AA_ALPHABET[i] ] = float(tokens[i])
        curr_line += 1
    fin.close()
    #m = fill_matrix_rev_values(m)
    return m

def fill_matrix_rev_values(m):
    for c in AA_ALPHABET:
        for d in AA_ALPHABET:
            if d not in m[c].keys():
                m[c][d] = float(m[d][c])
    return m    

def read_specs(p):
    msa_nodes = {} # key = msapath, value = array of *.dat tuples
    msa_weights = {}
    seed = ""
    fin = open(p, "r")
    for l in fin.readlines():
        if l.startswith("#"): # skip lines with comments
            continue
        if l.startswith("msapaths"):
            tokens = l.split()
            for t in tokens[1:]:
                msa_nodes[t] = None
        elif l.startswith("seed"):
            tokens = l.split()
            seed = tokens[1]
        elif l.startswith("compare"):
            tokens = l.split()
            msa_nodes[tokens[3]] = (tokens[1],tokens[2])
        elif l.startswith("msaweight"):
            tokens = l.split()
            msa_weights[ tokens[1] ] = float(tokens[2])
    fin.close()
    return [msa_nodes, seed, msa_weights]

def get_msa_len(p):
    #print path
    #os.system("cat " + path)
    #print p
    fin = open(p, "r")
    lines = fin.readlines()
    fin.close()
    lines[1].strip()
    return lines[1].split()[1].__len__()


def align_msas(msapaths, seed):
    # Find the longest MSA. . .
    longest_len = 0
    longest_msa = None
    for p in msapaths:
        len = get_msa_len(p)
        if len > longest_len:
            longest_len = len
            longest_msa = p
        print "\n. MSA [", p, "] has", len, "sites."
    #exit()
    
    print "\n. I'm aligning the alignments, using taxa", seed, "as the seed."
    
    alg_focus_seq = {}
    for p in msapaths:
        alg_focus_seq[p] = ""
        fin = open(p, "r")
        for l in fin.readlines():
            #print l
            if l.startswith(seed):
                tokens = l.split()
                alg_focus_seq[p] = tokens[1]
                #alg_focus_seq[p] = re.sub("\-", "", alg_focus_seq[p])
        fin.close()
        if alg_focus_seq[p].__len__() == 0:
            print "\n. Hmmm, I couldn't find a seed sequence for ", seed, "in ", p
        
    # verify that the focus sequences are actually the SAME sequence 
    for p in msapaths:
        if re.sub("\-", "", alg_focus_seq[p]).__len__() != re.sub("\-", "", alg_focus_seq[longest_msa]).__len__():
            print "\n. Hmmm, I had to stop because the seed sequence is not the same in all your alignments."
            print ". Tip: Check if you used a truncated version of the seed sequence in some of the alignments."
            print ""

    # Align the focus sequences. . .
    msa_refsite_mysite = {} # key = alg, value = hash, where key = site in longest alg, value = my site (or None)
    msa_mysite_refsite = {} # the reverse lookup of msa_refsite_mysite
    for alg in msapaths:
        print "\n. Aligning [", alg, "] to [", longest_msa, "]"
        msa_refsite_mysite[alg] = {}
        msa_mysite_refsite[alg] = {}
        ref_ptr = 1
        my_ptr = 1
        while my_ptr < alg_focus_seq[alg].__len__()+1 and ref_ptr < alg_focus_seq[longest_msa].__len__()+1:
            if alg_focus_seq[alg][my_ptr-1] == alg_focus_seq[longest_msa][ref_ptr-1]:
                #print "anccomp 169: matching", my_ptr, "in", alg, "to", ref_ptr, "in", longest_msa, "states:", alg_focus_seq[alg][my_ptr-1], alg_focus_seq[longest_msa][ref_ptr-1]   
                msa_refsite_mysite[alg][ref_ptr] = my_ptr
                msa_mysite_refsite[alg][my_ptr] = ref_ptr
                ref_ptr += 1
                my_ptr += 1
                #continue
            elif alg_focus_seq[longest_msa][ref_ptr-1] == "-":
                while alg_focus_seq[longest_msa][ref_ptr-1] == "-" and ref_ptr < alg_focus_seq[longest_msa].__len__()+1:
                    ref_ptr += 1
            elif alg_focus_seq[alg][my_ptr-1] == "-":
                while alg_focus_seq[alg][my_ptr-1] == "-" and my_ptr < alg_focus_seq[alg].__len__()+1:
                    my_ptr += 1
    
    #print msa_refsite_mysite
    #print alg_focus_seq 
    return [longest_msa, msa_refsite_mysite, msa_mysite_refsite]



#################################################################
#
# main part 1. . .
# 
# This block is all about reading the alignments and aligning the alignments
show_splash()
specpath = ap.getArg("--specpath")
modelpath = ap.getArg("--modelpath")
winsizes = ap.getList("--window_sizes", type=int)
adir = get_anccomp_dir(ap)
if not os.path.exists(adir):
    os.system("mkdir " + adir)
[msa_nodes, seed, msa_weights] = read_specs( specpath ) #msa_nodes = hashtable, key = msapath, value = [anc. node 1 path, anc. node 2 path]
m = get_matrix(modelpath)
[longest_msa, msa_refsite_mysite, msa_mysite_refsite] = align_msas(msa_nodes.keys(), seed)
############################################################################

def compare_dat_files(patha, pathb, m, winsize):
    """Compares *.dat files for two ancestors from the same alignment."""
    """Returns windata and consdata"""
    """windata[site] = h score"""
    """consdata[site] = r score"""
    fin = open(patha, "r")
    linesa = fin.readlines()
    fin.close()
    fin = open(pathb, "r")
    linesb = fin.readlines()
    fin.close()
    
    kl = {}
    ent = {}
    cent = {}
    ig = {}
    h = {} # KL distance
    r = {} # key = site, value = conservation difference
    gd = {}
    
    limstart = 0
    x = ap.getOptionalArg("--limstart")
    if x != False:
        limstart = int(x)
    
    curr_site = 1
    aptr = 0
    bptr = 0
    while (aptr < linesa.__len__() and bptr < linesb.__len__()):
        if curr_site >= limstart:
            al = linesa[aptr]
            bl = linesb[bptr]        
            while not al.startswith( curr_site.__str__() ):
                aptr += 1
                al = linesa[aptr]
            while not bl.startswith( curr_site.__str__() ):
                bptr += 1
                bl = linesb[bptr]
            
            ha = line_2_hash(al)
            hb = line_2_hash(bl)
                        
            h[curr_site] = d(ha, hb, m)
        curr_site += 1
        aptr += 1
        bptr += 1    
    windata = window_analysis(h, winsize, ap)
    #print windata
    #consdata = window_analysis(r, winsize, ap)
    return windata

def compute_cdata(msapath, winsize):
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

def normalize_m(m):
    sum = 0.0
    for i in AA_ALPHABET:
        for j in AA_ALPHABET:
            sum += m[i][j]
    new_m = {}
    for i in AA_ALPHABET:
        new_m[i] = {}
        for j in AA_ALPHABET:
            new_m[i][j] = m[i][j] / sum
    return new_m

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
            
def hdist(ppx, ppy, m):
    # 1. Get reciprocal KL distance. . .     
    klxy = kl_divergence(ppx, ppy, m)
    klyx = kl_divergence(ppy, ppx, m)

    klsum = klxy + klyx
        
    # 2. Calculate the importance of the a.a. changes. . .
    #m = normalize_m(m)
    s = 0.0
    for a in AA_ALPHABET:
        for b in AA_ALPHABET:
            if ppx[a] > 1.0/(AA_ALPHABET.__len__()) and ppy[b] > ppx[b]:
                s_part = (1.0/m[a][b])
                s_part *= ppx[a] * ppy[b]
                s += s_part
    #print "\t s = %.3f"%(100*s)

    # 3. Adjust the sign +/- based on the direction of entropy change
    enty = entropy(ppy)
    entx = entropy(ppx)
    e_dir = 1
    if (entx - enty) <= ENTROPY_CUTOFF:
        e_dir *= -1

    # 4. Combine evidence from 1, 2, and 3
    # Multiple by 100 to make the numbers interpretable
    combostyle = ap.getOptionalArg("--combo_method")
    if combostyle == False:
        combostyle = "product"
    if combostyle == "sum":
        w = float(ap.getOptionalArg("--training_weight"))
        if w > 1.0:
            print "\n. Error: the --training_weight must range from 0.0 to 1.0."
        h = w*klsum + (1-w)*s + e_dir * 100
    elif combostyle == "product":
        h = klsum * s * e_dir * 100
    return h

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

def d(ancx, ancy, m):    
    d_sum = 0.0
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
    d_sum = hdist(ancx, ancy, m)       
    return d_sum

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
    for i in range(min, max):
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

def blend_msa_data(msa_data):
    #algweight = 1.0 / msa_data.keys().__len__()

    limstart = 0
    x = ap.getOptionalArg("--limstart")
    if x != False:
        limstart = int(x)

    bdata = {} # key = reference site from MSA-MSA, value = blended score
    for ref_site in msa_refsite_mysite[longest_msa].keys():
        
        
        if ref_site < limstart:
            continue
        
        """Does this ref_site exist in all alignments?"""
        skip_this_site = False        
        for alg in msa_data:
            if False == msa_refsite_mysite[alg].keys().__contains__( ref_site ):
                bdata[ref_site] = 0.0
                skip_this_site = True

        if skip_this_site == False:
           for alg in msa_data:
                mysite = msa_refsite_mysite[alg][ref_site]
                myscore = 0.0
                if mysite != None:
                    myscore = msa_data[alg][mysite]
                if mysite != None:
                    if ref_site in bdata:            
                        bdata[ref_site] += (msa_weights[alg] * myscore)
                    else:
                        bdata[ref_site] = msa_weights[alg] * myscore
    return bdata

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
    data = fill_missing_sites(data)
    
    cranpath = outpath + ".cran"
    cranout = open(cranpath, "w")
    miny = 0.0
    maxy = 0.0
    minx = 0
    maxx = 0
    x = "x <- c("
    y = "y <- c("
    sites = data.keys()
    sites.sort()
    minx = sites[0]
    maxx = sites[ sites.__len__()-1 ]
    
    for s in sites:
        yval = 0
        if data[s] != None:
            if data[s] < 0.00000001:
                yval = 0.0
            else:
                yval = math.log( float(data[s]) + 1.0 )
    
    for s in sites:
        #print s
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

def write_table(hdata, msa_hscores):
    foutpath = get_table_outpath(ap, tag="anccomp.table.txt")
    fout = open(foutpath, "w")
    sites = msa_refsite_mysite[longest_msa].keys()
    sites.sort()
    msapaths = msa_hscores.keys()
    msapaths.sort()
    header = ""
    header += "ref_site\t"
    header += "h\t"
    #header += "r\t"
    #header += "c\t"
    for m in msapaths:    
        header += "h(" + msapath_to_short(m) + ")\t"
        #header += "r(" + msapath_to_short(m) + ")\t"
        #header += "c(" + msapath_to_short(m) + ")\t"
    header += "\n"
    fout.write(header)
    for s in sites:
        line = s.__str__() + "\t"
        line += "%.3f"%hdata[s] + "\t"
        #line += "%.3f"%rdata[s] + "\t"
        #line += "%.3f"%cdata[s] + "\t"
        for m in msapaths:
            if s in msa_refsite_mysite[m]:
                mys = msa_refsite_mysite[m][s]

            #print s, msa_refsite_mysite[m][s]
            #print s, msa_hscores[m][ mys ]

                line += "%.3f"%msa_hscores[m][mys] + "\t"
            #line += "%.3f"%msa_rscores[m][mys] + "\t"
            #line += "%.3f"%msa_cscores[m][mys] + "\t"
            else:
                line += "\t"
        line += "\n"
        #print line
        fout.write(line)
    fout.close()

def rank_sites(blended_hdata):
    hscore_refsite = {}
    #rscore_refsite = {}
    sites = blended_hdata.keys()
    sites.sort()
    for s in sites:
        hval = "%.5f"%blended_hdata[s]
        hval = blended_hdata[s]
        if hval in hscore_refsite:
            hscore_refsite[hval].append(s)
        else:
            hscore_refsite[hval] = [ s ]
    
    hscores = hscore_refsite.keys()
    hscores.sort(reverse=True)

    msa_thatanc_lines = {}
    msa_thisanc_lines = {}
    for msapath in msa_nodes:
        that_ancpath = msa_nodes[msapath][0]
        fin = open(that_ancpath, "r")
        msa_thatanc_lines[msapath] = fin.readlines()
        fin.close()
        this_ancpath = msa_nodes[msapath][1]
        fin = open(this_ancpath, "r")
        msa_thisanc_lines[msapath] = fin.readlines() 
        fin.close()
        
    fout = open( get_table_outpath(ap, tag="h.ranked.txt"), "w")
    for i in range(0, hscores.__len__()):
        h = hscores[i]
        for refsite in hscore_refsite[h]:
            lout = "--> "
            #lout = "\n. Site "
            #lout += refsite.__str__()
            lout += " h= %.3f"%h
            lout += " rank: " + (i+1).__str__()
            #lout += " r= %.3f"%rscores[i]
            lout += "\n" 
            for msapath in msa_nodes:
                if refsite in msa_refsite_mysite[msapath]:
                    mys = msa_refsite_mysite[msapath][refsite]
                    lout += msapath_to_short(msapath) + "\th= %.3f"%msa_hscores[msapath][mys] + "\n"
                    for l in msa_thatanc_lines[msapath]:
                        if l.startswith(mys.__str__() + " "):
                            lout += "   " + msapath_to_short(msa_nodes[msapath][0]) + " site " + l
                    for l in msa_thisanc_lines[msapath]:
                        if l.startswith(mys.__str__() + " "):
                            lout += "   " + msapath_to_short(msa_nodes[msapath][1]) + " site " + l
            lout += "\n"
            fout.write(lout)
    fout.close()



""" Main part 2. . ."""
cranpaths = []  # this array will hold R scripts that we'll execute (in R) at the end.
#for w in winsizes:
w = 1
msa_hscores = {} # key = msa path, value = data from compare_dat_files
msa_cscores = {} # basic conservation scores
for msapath in msa_nodes:
    this_ancpath = msa_nodes[msapath][0]
    that_ancpath = msa_nodes[msapath][1]
    print "\n. I'm comparing the ancestor [", this_ancpath, "] to [", that_ancpath, "] with a smoothing window =", w
    hdata = compare_dat_files(this_ancpath, that_ancpath, m, w)
    
    msa_hscores[msapath] = hdata
    cdata = compute_cdata(msapath, w)
    #msa_cscores[msapath] = cdata

"""Average the data over multiple MSA algorithms...."""
blended_hdata = blend_msa_data(msa_hscores)
#blended_cdata = blend_msa_data(msa_cscores)
    
for w in winsizes:
    #print blended_hdata
    blended_hdata = window_analysis(blended_hdata, w, ap)
    
    #print "Using these site:", blended_hdata.keys()
    
    """Write table data, and rank data, but only for the non-smoothed data (i.e. w=1)."""
    if w == 1:
        write_table(blended_hdata, msa_hscores)
        rank_sites(blended_hdata)
            
    """ Plot H scores """
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
    
""" Finally, execute all the R scripts. . ."""
for c in cranpaths:
    print "\n. I'm plotting results, using the R script written at [", c, "] . . ."
    os.system("r --no-save < " + c)

print "\n\n. Finished.  Results were written to the folder", get_plot_outpath(ap)

print "\n. Goodbye."
