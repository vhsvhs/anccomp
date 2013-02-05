import math, os, re, sys, time
from config import *
from splash import *
import scipy.stats as ss
from argparser import *
ap = ArgParser(sys.argv)

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
    msa_nodes = {} # key = msapath, value = array of *.dat tuples
    msa_weights = {}
    seed = ""
    ap.params["msaname"] = {}
    fin = open(p, "r")
    for l in fin.readlines():
        if l.startswith("#"): # skip lines with comments
            continue
        if l.startswith("msapaths"):
            tokens = l.split()
            for t in tokens[1:]:
                ap.params["msaname"][t] = t
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
        elif l.startswith("msaname"):
            tokens = l.split()
            ap.params["msaname"][ tokens[1] ] = tokens[2]
    fin.close()
    ap.params["msa_nodes"] = msa_nodes
    ap.params["seed"] = seed
    ap.params["msa_weights"] = msa_weights
    
    #check_specs()
    
    return [msa_nodes, seed, msa_weights]

def check_specs():
    for msapath in ap.params["msapaths"]:
        if msapath not in ap.params["msaname"]:
            ap.params["msaname"][ msapath ] = msapath

    for msaname in ap.params["msaname"]:
        if msaname not in ap.params["msapaths"]:
            print "\n. Warning, you defined an 'msaname' for", msaname, ", but you never specified a path to that sequence alignment.\n"
            print "\n. I will continue, but there may be errors. . .\n"

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
    
    msa_seedseq = {}
    for p in msapaths:
        msa_seedseq[p] = ""
        fin = open(p, "r")
        for l in fin.readlines():
            #print l
            if l.startswith(seed):
                tokens = l.split()
                msa_seedseq[p] = tokens[1]
                #msa_seedseq[p] = re.sub("\-", "", msa_seedseq[p])
        fin.close()
        if msa_seedseq[p].__len__() == 0:
            print "\n. Hmmm, I couldn't find a seed sequence for ", seed, "in ", p
        
    # verify that the focus sequences are actually the SAME sequence 
    for p in msapaths:
        if re.sub("\-", "", msa_seedseq[p]).__len__() != re.sub("\-", "", msa_seedseq[longest_msa]).__len__():
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
        ref_site = 1
        my_site = 1
        while my_site-1 < msa_seedseq[alg].__len__() and ref_site-1 < msa_seedseq[longest_msa].__len__():
            if msa_seedseq[alg][my_site-1] == msa_seedseq[longest_msa][ref_site-1]:
                my_state = msa_seedseq[alg][my_site-1]
                ref_state = msa_seedseq[longest_msa][ref_site-1]
                if my_state != ref_state:
                    print "\n. Hmmm, something went wrong with the meta-alignment."
                    print " (anccomp_tools.py point 205."
                    print "anccomp 169: matching", my_site, "in", alg, "to", ref_site, "states:", my_state, ref_state    
                    print "\n"
                    exit(1)
                msa_refsite_mysite[alg][ref_site] = my_site
                msa_mysite_refsite[alg][my_site] = ref_site
                ref_site += 1
                my_site += 1
                #continue
            elif msa_seedseq[longest_msa][ref_site-1] == "-":
                while msa_seedseq[longest_msa][ref_site-1] == "-" and ref_site < msa_seedseq[longest_msa].__len__()+1:
                    ref_site += 1
            elif msa_seedseq[alg][my_site-1] == "-":
                while msa_seedseq[alg][my_site-1] == "-" and my_site < msa_seedseq[alg].__len__()+1:
                    my_site += 1
    
    """Write a log file about the meta-alignment."""
    msa_ids = {}
    msa_ids[ longest_msa ] = 1
    ids_msa = {}
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
    for i in ids:
        fout.write( "M" + i.__str__() + " : " + ap.params["msaname"][ ids_msa[i] ]  + "\n")
    fout.write("\n")
    header = ""
    fout.write("Residues shown in parentheses express the value of taxon " + ap.params["seed"] + " at the corresponding site.\n\n")
    for i in ids:
        header += "M" + i.__str__() + "\t"
    fout.write(header + "\n\n")
    for site in range(0, msa_seedseq[longest_msa].__len__()):
        line = ""
        for i in ids:
            if (site+1) in msa_refsite_mysite[ ids_msa[i] ]:
                mysite = msa_refsite_mysite[ ids_msa[i] ][ site+1 ]
                state = msa_seedseq[ids_msa[i]][mysite-1]
                line += mysite.__str__() + " (" + state + ")\t"
            else:
                line += "-\t"
        fout.write( line + "\n" )
    fout.close()
    return [longest_msa, msa_refsite_mysite, msa_mysite_refsite]

def compare_dat_files(patha, pathb, m, winsize, method="h"):
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
    
    results = {} 
    
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
            
            ancx = line_2_hash(al)
            ancy = line_2_hash(bl)
                        
            results[curr_site] = d(ancx, ancy, m, method=method)
        curr_site += 1
        aptr += 1
        bptr += 1    
    return results


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

    elif method == "hp":
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
        w = float(ap.getOptionalArg("--training_weight"))
        if w > 1.0:
            print "\n. Error: the --training_weight must range from 0.0 to 1.0."
        h = w*klsum + (1-w)*s + e_dir
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
    if method == "h" or method == "hp":
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

def blend_msa_data(msa_data,msa_refsite_mysite,longest_msa,msa_weights):
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
        for msa in msa_data:
            #print "651: msa = ", msa
            if False == msa_refsite_mysite[msa].keys().__contains__( ref_site ):
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

def write_table(hdata, msa_scores,msa_refsite_mysite,longest_msa,method="h"):
    foutpath = get_table_outpath(ap, tag=method + ".summary.txt")
    fout = open(foutpath, "w")
    sites = msa_refsite_mysite[longest_msa].keys()
    sites.sort()
    msapaths = msa_scores.keys()
    msapaths.sort()
    header = ""
    header += "site\t"
    header += method + "\t"
    if msapaths.__len__() > 1:
        for m in msapaths:    
            header += method + "(" + ap.params["msaname"][m] + ")\t"
    header += "\n"
    fout.write(header)
    for s in sites:
        line = s.__str__() + "\t"
        line += "%.3f"%hdata[s] + "\t"
        if msapaths.__len__() > 1:
            for m in msapaths:
                if s in msa_refsite_mysite[m]:
                    mys = msa_refsite_mysite[m][s]
                    line += "%.3f"%msa_scores[m][mys] + "\t"
                else:
                    line += "\t"
        line += "\n"
        fout.write(line)
    fout.close()

def rank_sites(blended_data, msa_nodes, msa_refsite_mysite, msa_scores, method="h", writetable=True):
    """This method takes the site-value data and orders the sites in ascending order.
    Data is written to the *.ranked.txt file."""
    
    #print "696:", msa_nodes
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
    for msapath in msa_nodes:
        that_ancpath = msa_nodes[msapath][0]
        fin = open(that_ancpath, "r")
        msa_thatanc_lines[msapath] = fin.readlines()
        fin.close()
        this_ancpath = msa_nodes[msapath][1]
        fin = open(this_ancpath, "r")
        msa_thisanc_lines[msapath] = fin.readlines() 
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
                for msapath in msa_nodes:
                    if refsite in msa_refsite_mysite[msapath]:
                        mys = msa_refsite_mysite[msapath][refsite]
                        lout += ap.params["msaname"][ msapath ] + " site " + mys.__str__() + " (" + method + " = %.3f"%msa_scores[msapath][mys] + ")\n"
                        for l in msa_thatanc_lines[msapath]:
                            if l.startswith(mys.__str__() + " "):
                                lout += "   " + datpath_to_short(msa_nodes[msapath][0]) + "\tsite " + re.sub("site ", "", l)
                        for l in msa_thisanc_lines[msapath]:
                            if l.startswith(mys.__str__() + " "):
                                lout += "   " + datpath_to_short(msa_nodes[msapath][1]) + "\tsite " + re.sub("site ", "", l)
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
    elif min > 0:
        minbin = 0.0
    if max > 0:
        maxbin = 0.0
        while maxbin < max:
            maxbin += stride
    i = minbin
    while i < maxbin:
        bins.append(i)
        i += stride
    return bins

def get_bin(value, bins):
    ret = bins[0]
    for i in bins:
        if i <= value:
            ret = i
    return ret

def plot_histogram(metric_data, tag):
    maxx = None
    minx = None
    for metric in metric_data:
        data = metric_data[metric]
        for site in data:
            value = data[site]
            if maxx == None:
                maxx = value
            if minx == None:
                minx = value
            if maxx < value:
                maxx = value
            if minx > value:
                minx = value
    binwidth = 0.1
    if maxx - minx > 1000:
        binwidth = 100
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
    elif maxx - minx > 0.1:
        binwidth = 0.005
    
    bins = get_bins(minx, maxx, binwidth)
    metric_bin_count = {}
    
    for metric in metric_data:
        metric_bin_count[metric] = {}
        for b in bins:
            metric_bin_count[metric][b] = 0.0
        
        data = metric_data[metric]
        n = data.__len__()
        for site in data:
            this_bin = get_bin( data[site], bins )
            print data[site], this_bin, bins
            metric_bin_count[metric][this_bin] += 1.0/n
        
    
    """
    Write the R script.
    """        
    pdfpath = tag + ".pdf"
    
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    #cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"
    
    for metric in metric_data:
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
    
    cranpath = tag + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)       
        

def correlate_metrics(ma, mb, ma_site_val, mb_site_val, tag):    
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
    os.system("r --no-save < " + path)

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
    os.system("r --no-save < " + path)
    
    fout = open( get_plot_outpath(ap, tag="corr-value." + tag) + ".txt" , "w")
    fout.write( pearsons.__str__() + "\n")
    for v in value_a_b:
        fout.write( v[0].__str__() + "\t" + v[1].__str__() + "\n")
    fout.close()
    
    return [spearmans, pearsons]
    
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