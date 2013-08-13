from config import *
from anccomp_tools import *

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
    j = 0
    for i in range(0, seed_seq.__len__()):
        seed_seq_no_indels = re.sub("-", "", seed_seq[i:])
        if seed_seq_no_indels.startswith(pdb_seq[0:20]):
            return i

def get_pdb_offset(pdb_path):
    fin = open(pdb_path, "r")
    for l in fin.xreadlines():
        if l.startswith("ATOM"):
            site = int(l.split()[5])
            return site
    fin.close()
    return 1

def pymol_viz_helper(ap, data, seedseq):
    pass

def do_pymol_viz(ap, mb, seedseq):
    if PDBTOOLSDIR == False:
        print "\n. Sorry, but you need pdb-tools installed in order to analyze your PDB file."
        return
    pdb_path = ap.getOptionalArg("--pdb_path")
    #ref_site = int(ap.getOptionalArg("--pdb_starts_at"))
    #if ref_site == False:
    #    print "\n. You need to specify --pdb_starts_at"
    #    exit()
    print "\n. Plotting data on the PDB file", pdb_path
    
    pdb_seq = get_seq_from_pdb(pdb_path)
    print pdb_seq
    lnick = ap.params["msa_path2nick"][ ap.params["longest_msa"] ]
    seed_seq = get_seed_seq(  ap.params["msa_nick2path"][lnick], ap.params["seed"]  )
    ref_site = get_pdb_startsite(pdb_seq, seed_seq)
    pdb_offset = get_pdb_offset(pdb_path)
    #print "50:"
    #print pdb_seq[0:10]
    #print seed_seq[ref_site:(ref_site+10)]
    
    for metric in mb:
        data = mb[metric]
        #print data
        #print seedseq
        #print seedseq.__len__()
        
        outlines = ""
        outlines += "cmd.load(\"" + pdb_path + "\")\n"
        #outlines += "select all\n"
        #outlines += "color white\n"
        pdb_site = pdb_offset-1
        
        maxh = data[data.keys()[0]]
        minh = data[data.keys()[0]]
        for site in data:
            if maxh < data[site]:
                maxh = data[site]
            if minh > data[site]:
                minh = data[site]
        
        while (ref_site < seedseq.__len__() ):
            if seedseq[ref_site] != "-":
                if ref_site in data:
                    h = data[ref_site]
                    if h > 0:
                        red = 1.0#0.9 + (h/maxh)*0.1
                        green = 1.0 - (h/maxh)#0.9 + (h/maxh)*0.1
                        blue = 1.0 - (h/maxh)
                    if h < 0:
                        blue = 1.0#0.9 + (h/minh)*0.1
                        green = 1.0 - (h/minh)#0.9 + (h/minh)*0.1
                        red = 1.0 - (h/minh)
                    if h == 0:
                        blue = 1.0#0.9
                        green = 1.0#0.9
                        red = 1.0#0.9
                    #print ref_site, seedseq[ref_site], h, blue, green, red
                    this_color = "[" + red.__str__() + "," + green.__str__() + "," + blue.__str__() + "]"
                    outlines += "cmd.set_color('color" + pdb_site.__str__() + "'," + this_color + ")\n"
                    outlines += "cmd.color(\"color" + pdb_site.__str__() + "\", \"resi " + pdb_site.__str__() + "\")\n"
                pdb_site += 1
            ref_site += 1
                
        scriptpath = "pymol_script." + metric + ".p"
        fout = open(scriptpath, "w")
        fout.write(outlines)
        fout.close()
        
        print "\n\n. Run the script at " + os.getcwd() + "/" + scriptpath