from config import *
from anccomp_tools import *

from smith_waterman import *

#
# The PyMol visualization requires that you have pdb-tools (i.e. Mike Harm's tools).
#
PDBTOOLSDIR = ap.getOptionalArg("--pdbtoolsdir")
PDBSEQ = "python " + PDBTOOLSDIR + "/pdb_seq.py"

def get_seq_from_pdb(pdbpath):
    outpath = ".pdbseq"
    os.system(PDBSEQ + " " + pdbpath + " > " + outpath)
    fin = open(outpath, "r")
    lines = fin.readlines()
    fin.close()
    pdb_seq = ""
    for l in lines:
        if l.startswith(">") and pdb_seq.__len__() > 1:
            return pdb_seq
        elif l.startswith(">"):
            continue
        pdb_seq += l.strip()
    return pdb_seq


def get_pdb_startsite(pdb_seq, seed_seq):
    for i in range(0, seed_seq.__len__()):
        seed_seq_no_indels = re.sub("-", "", seed_seq[i:])
        if seed_seq_no_indels.startswith(pdb_seq[0:10]):
            n = 11
            while n < pdb_seq.__len__():
                #print n
                if seed_seq_no_indels.startswith(pdb_seq[0:n]):
                    n += 1
                else:
                    return (i,n)
            return (i,n)
    return (-1,-1)


def get_pdb_offset(pdb_path):
    """What is the first site number in the protein sequence in the PDB?"""
    fin = open(pdb_path, "r")
    for l in fin.xreadlines():
        if l.startswith("ATOM"):
            site = int(l.split()[4])
            return site
    fin.close()
    return 1

def get_pdb_sites(pdb_path):
    """Returns a sorted list with the site numbers in the PDB."""
    sites = []
    fin = open(pdb_path, "r")
    for l in fin.xreadlines():
        if l.startswith("ATOM"):
            site = int(l.split()[4])
            if site not in sites:
                sites.append(site)
    fin.close()
    return sites    


def pymol_viz_helper(ap, data, seedseq):
    pass



def getpdbsites(ap):
    # Get the PDB sequence:
    pdb_path = ap.getOptionalArg("--pdb_path")
    print "\n. I'm plotting scores onto the PDB file", pdb_path
    pdb_seq = get_seq_from_pdb(pdb_path)
    #print "\n. The PDB sequence is:\n", pdb_seq
    
    # Does the PDB sequence match the seed sequence?
    #lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    best_match_msanick = None
    best_identity = 0
    align2 = None
    seq = None
    if ap.params["msa_comparisons"].keys().__len__() > 1:
        print "\n. I'm determining which MSA best matches the sequence in your PDB."
    for msanick in ap.params["msa_comparisons"]:
        this_seq = dat2seq( ap.params["msa_comparisons"][msanick][1] )
        (this_identity, this_score, this_align1, this_align2, this_symbol) = needle(this_seq, pdb_seq)
        #print msanick, this_identity
        if this_identity > best_identity: # this sequence is a better match
            seq = this_seq
            best_match_msanick = msanick
            best_identity = this_identity
            align2 = this_align2
            
    pdbsites = get_pdb_sites(pdb_path)
    
    pdbseqsite2structsite = {}
    pdbseqsite2refsite = {} # key = 1-based pdb site number, value = 1-based ref site number
    pdbsite = 0
    for i in range(0, align2.__len__()):
        if align2[i] != "-":
            structsite = pdbsites[ pdbsite ]
            pdbseqsite2structsite[pdbsite] = structsite
            pdbseqsite2refsite[ pdbsite ] = i
            #print "refsite:", i+1, seq[i], "pdbseqsite:", pdbsite+1, pdb_seq[pdbsite], "struct site:", structsite 
            pdbsite += 1
    return pdbseqsite2refsite
    
def do_pymol_viz(ap, mb):    
    if PDBTOOLSDIR == False:
        print "\n. Sorry, but you need pdb-tools installed in order to analyze your PDB file."
        return
    
    # Get the PDB sequence:
    pdb_path = ap.getOptionalArg("--pdb_path")
    print "\n. I'm plotting scores onto the PDB file", pdb_path
    pdb_seq = get_seq_from_pdb(pdb_path)
    #print "\n. The PDB sequence is:\n", pdb_seq
    
    # Does the PDB sequence match the seed sequence?
    #lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    best_match_msanick = None
    best_identity = 0
    align2 = None
    seq = None
    if ap.params["msa_comparisons"].keys().__len__() > 1:
        print "\n. I'm determining which MSA best matches the sequence in your PDB."
    for msanick in ap.params["msa_comparisons"]:
        this_seq = dat2seq( ap.params["msa_comparisons"][msanick][1] )
        (this_identity, this_score, this_align1, this_align2, this_symbol) = needle(this_seq, pdb_seq)
        #print msanick, this_identity
        if this_identity > best_identity: # this sequence is a better match
            seq = this_seq
            best_match_msanick = msanick
            best_identity = this_identity
            align2 = this_align2
            
    pdbsites = get_pdb_sites(pdb_path)
    
    pdbseqsite2structsite = {}
    pdbseqsite2refsite = {} # key = 1-based pdb site number, value = 1-based ref site number
    pdbsite = 0
    for i in range(0, align2.__len__()):
        if align2[i] != "-":
            structsite = pdbsites[ pdbsite ]
            pdbseqsite2structsite[pdbsite] = structsite
            pdbseqsite2refsite[ pdbsite ] = i
            #print "refsite:", i+1, seq[i], "pdbseqsite:", pdbsite+1, pdb_seq[pdbsite], "struct site:", structsite 
            pdbsite += 1
    
    # correct up to here.
    
    for metric in mb:
        data = mb[metric]
        outlines = ""
        outlines += "cmd.load(\"" + os.path.abspath(pdb_path) + "\")\n"
        #pdb_site = pdb_offset
        
        maxh = data[data.keys()[0]]
        minh = data[data.keys()[0]]
        allvals = []
        for site in data:
            allvals.append( data[site] )
            if maxh < data[site]:
                maxh = data[site]
            if minh > data[site]:
                minh = data[site]
        
        red = "red"
        orange = "orange"
        blue = "blue"
        snow = "white"
        
        sdev = sd(allvals)
        avg = mean(allvals)

        tokens = pdb_path.split("/")
        ancname = tokens[ tokens.__len__()-1 ].split(".")[0]
        scriptpath = get_output_dir(ap) + "/pymol_script." + metric + "." + ancname + ".p"        
        fout = open(scriptpath, "w")
        fout.write("cmd.load(\"" + os.path.abspath(pdb_path) + "\")\n")
        
        pdbsites = pdbseqsite2refsite.keys()
        pdbsites.sort()
        for site in pdbsites:
            refsite = pdbseqsite2refsite[site] # get the reference site for this PDB site
            h = -1
            if (refsite+1) in data:
                h = data[refsite+1]
                if h < avg:
                    if h <= avg - 6*sdev:
                        this_color = red
                    elif h <= avg - 4*sdev:
                        this_color =  orange
                    elif h <= avg - 2*sdev:
                        this_color = blue
                    else:
                        this_color = snow 
                elif h > avg:
                    if h >= avg + 6*sdev:
                        this_color = red
                    elif h >= avg + 4*sdev:
                        this_color = orange
                    elif h >= avg + 2*sdev:
                        this_color = blue
                    else:
                        this_color = snow                                
            else:
                this_color = snow
            #print ". Site", refsite+1, seq[refsite], pdbseqsite2structsite[site], pdb_seq[site], h, this_color 

            fout.write("cmd.select(\"x\", \"" + ancname + " and resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n")
            fout.write("cmd.color(\"" + this_color + "\", \"x\")\n")
            #fout.write("cmd.color(\"" + this_color + "\", \"resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n") 
            #fout.write("cmd.color(\"" + this_color + "\", \"resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n")
        fout.close()
        
        print "\n\n. Load PyMol and run the following script: " + os.getcwd() + "/" + scriptpath