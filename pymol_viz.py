from config import *
from anccomp_tools import *
from smith_waterman import *

#
# The PyMol visualization requires that you have pdb-tools (i.e. Mike Harm's tools).
#
PDBTOOLSDIR = ap.getOptionalArg("--pdbtoolsdir")
#PDBTOOLSDIR = "/"
PDBSEQ = "python " + PDBTOOLSDIR + "/pdb_seq.py"


def get_seq_from_pdb(pdbpath):
    outpath = ".pdbseq"
    os.system(PDBSEQ + " " + pdbpath + " > " + outpath)
    fin = open(outpath, "r")
    lines = fin.readlines()
    fin.close()
    os.system("rm " + outpath)
    pdb_seq = ""
    for l in lines:
        if l.startswith(">") and pdb_seq.__len__() > 1:
            return pdb_seq
        elif l.startswith(">"):
            continue
        pdb_seq += l.strip()
    return pdb_seq


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


"""Depricated.
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
            
    pdbsites = get_pdb_sites(pdb_path) # a sorted list of site numbers in the PDB
    
    pdbseqsite2structsite = {} # key = pdb sequence site (1-based), value 
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
"""

def do_pymol_viz(ap): 
    """ This method generates two different PyMOL visualizaitons.
    First, averaged Df, k, and p scores are painted onto the PDB, resulting in three output PyMOL *.pse sessions.
    Second, type 1/2/3 changes for each alignment and model are painted onto the PDB, resulting in NxM output *.pse
    sessions, where N = the number of alignment methods and M equals the number of phylogenetic models. 
    
    This method assumes that ap.params["pdb_path"] contains value.
    mb is the hashtable 'metric_blendeddata' from compare_ancs.py
    
    """ 
    
    # PDB Tools are required, in order to extract the sequence from the PDB file.
    if PDBTOOLSDIR == False:
        print "\n. Sorry, but you need pdb-tools installed in order to analyze your PDB file."
        return
        
    print "\n. I'm plotting scores onto the PDB file", ap.params["pdb"]
    
    #
    # First, build a map between the sites in the PDB file, and the sites in our sequences.
    #
    
    # Get the PDB sequence:
    pdb_seq = get_seq_from_pdb(ap.params["pdb"])
    best_match_msanick = None
    best_identity = 0
    align2 = None # sequence
    seq = None
    if ap.params["msa_comparisons"].keys().__len__() > 1:
        print "\n. I'm determining which MSA best matches the sequence in your PDB."
    # Find the MSA-model combo whose ancestral sequence best matches the PDB sequence. \\
    # We'll use this ancestral sequence to build the map to the PDB
    for msanick in ap.params["msa_comparisons"]:
        this_seq = dat2seq( ap.params["msa_comparisons"][msanick][1] )
        # Run Needleman-Wunch on this sequence versus the PDB sequence
        (this_identity, this_score, this_align1, this_align2, this_symbol) = needle(this_seq, pdb_seq)
        #print msanick, this_identity
        if this_identity > best_identity: # this ancestral sequence is a better match to the PDB.
            seq = this_seq
            best_match_msanick = msanick
            best_identity = this_identity
            align2 = this_align2 # align2 now equals the aligned PDB sequence, with gaps.
    pdbsites = get_pdb_sites(ap.params["pdb"]) # How are site numbers references in the PDB file?    
    pdbseqsite2structsite = {}
    pdbseqsite2refsite = {} # key = 1-based pdb site number, value = 1-based ref site number
    pdbsite = 0 # a counter to sites in the aligned-not-aligned sequence
    for i in range(0, align2.__len__()):
        if align2[i] != "-":
            structsite = pdbsites[ pdbsite ]
            pdbseqsite2structsite[pdbsite] = structsite
            pdbseqsite2refsite[ pdbsite ] = i
            #print "refsite:", i+1, seq[i], "pdbseqsite:", pdbsite+1, pdb_seq[pdbsite], "struct site:", structsite 
            pdbsite += 1

    """At the end of this function, pymol_script_paths will contain multiple PyMOL scripts.
    Running these scripts will result in several output PyMOL *.pse sessions."""
    pymol_script_paths = []

    #
    # VISUALIZATION #1: Paint Df, k, and/or p scores onto the PDB.
    #
    for metric in ap.params["metric_blendeddata"]:
        data = ap.params["metric_blendeddata"][metric]
        
        allvals = []
        mingain = None
        maxloss = None
        for site in data:
            allvals.append( data[site] )
            
        for i in range(0, allvals.__len__()):
            ax = allvals[i]
            if ax > 0:
                if mingain == None:
                    mingain = ax 
                elif mingain > ax:
                    mingain = ax
            elif ax < 0:
                if maxloss == None:
                    maxloss = ax
                elif maxloss < ax:
                    maxloss = ax
        
        red = "red"
        orange = "orange"
        blue = "blue"
        snow = "white"
        
        maxh = max(allvals) 
        minh = min(allvals)

        # Get the name of the ancestor, using some information from the msa_comparison
        that_ancpath = ap.params["msa_comparisons"][msanick][1]
        ancname = that_ancpath.split(".")[ that_ancpath.split(".").__len__()-2 ]
        ancname = ancname.split("/")[1]

        tokens = ap.params["pdb"].split("/")
        pdbmodelname = re.sub( ".pdb", "", tokens[ tokens.__len__()-1 ] )
        scriptpath = get_output_dir(ap) + "/pymol_script.scores." + metric + "." + pdbmodelname + ".p"        
        fout = open(scriptpath, "w")
        fout.write("cmd.load(\"" + os.path.abspath(ap.params["pdb"]) + "\")\n")
        fout.write("cmd.set_name(\"" + pdbmodelname + "\", \"" + ancname + "\")\n")
        
        pdbsites = pdbseqsite2refsite.keys()
        pdbsites.sort()
        for site in pdbsites:
            refsite = pdbseqsite2refsite[site] # get the reference site for this PDB site
            h = -1
            red = 255
            green = 255
            blue = 255
            if (refsite+1) in data:
                h = data[refsite+1]
                if h > 0:
                    #h = math.log(h)
                    red = 255
                    green = (1 - (h-mingain)/(maxh-mingain) ) * 255
                    blue =  (1 - (h-mingain)/(maxh-mingain) ) * 255
                    #green = math.sqrt((1 - (h-mingain)/(maxh-mingain) )) * 255
                    #blue =  math.sqrt((1 - (h-mingain)/(maxh-mingain) )) * 255
                if h < 0:
                    #h = -1 * math.log( abs(h) )
                    blue = 255
                    red =  (1 -(maxloss-h)/(maxloss-minh)) * 255
                    green =  (1 -(maxloss-h)/(maxloss-minh)) * 255
                    #red =  math.sqrt((1 -(maxloss-h)/(maxloss-minh))) * 255
                    #green =  math.sqrt((1 -(maxloss-h)/(maxloss-minh))) * 255
                if h < minh or h > maxh:
                    print "out of range", h
                    exit()
            else:
                this_color = snow
            #print "237:", site, pdbseqsite2structsite[site], red, green, blue
            this_color = "colorsite" + site.__str__()
            fout.write("cmd.set_color(\"" + this_color + "\", (" + int(red).__str__() + "," + int(green).__str__() + "," + int(blue).__str__() + ") )\n")    
            fout.write("cmd.color(\"" + this_color + "\", \"" + ancname + " and resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n")
        
        # Add some style:
        fout.write("cmd.hide('everything')\n")
        fout.write("cmd.show('cartoon', 'all')\n")
        
        # Save a ray-traced PNG
        png_path = get_output_dir(ap) + "/pymol_ray.scores." + metric + "." + ancname + ".png"
        fout.write("cmd.png(\"" + png_path + "\", 500, 500, ray=True)\n")
        
        # Save the PyMOL session as an *.pse file
        pymolsessionpath = get_output_dir(ap) + "/pymol_script.scores." + metric + "." + ancname + ".pse"
        fout.write("cmd.save(\"" + pymolsessionpath + "\")\n")
        fout.close()
        
        pymol_script_paths.append(scriptpath )

    #
    # Visualization #1b - redo #1, but for each msa-model pair
    #
    for metric in ap.params["metric_data"]:
        for msa in ap.params["metric_data"][metric]:
            data = ap.params["metric_data"][metric][msa]
            
            allvals = []
            mingain = None
            maxloss = None
            for site in data:
                allvals.append( data[site] )
                
            for i in range(0, allvals.__len__()):
                ax = allvals[i]
                if ax > 0:
                    if mingain == None:
                        mingain = ax 
                    elif mingain > ax:
                        mingain = ax
                elif ax < 0:
                    if maxloss == None:
                        maxloss = ax
                    elif maxloss < ax:
                        maxloss = ax
            
            red = "red"
            orange = "orange"
            blue = "blue"
            snow = "white"
            
            maxh = max(allvals) 
            minh = min(allvals)
    
            that_ancpath = ap.params["msa_comparisons"][msa][1]
            ancname = that_ancpath.split(".")[ that_ancpath.split(".").__len__()-2 ]
            ancname = ancname.split("/")[1]
    
            tokens = ap.params["pdb"].split("/")
            pdbmodelname = re.sub( ".pdb", "", tokens[ tokens.__len__()-1 ] )
            scriptpath = get_output_dir(ap) + "/pymol_script.scores." + metric + "." + msa + "." + pdbmodelname + ".p"        
            fout = open(scriptpath, "w")
            fout.write("cmd.load(\"" + os.path.abspath(ap.params["pdb"]) + "\")\n")
            fout.write("cmd.set_name(\"" + pdbmodelname + "\", \"" + ancname + "\")\n")
            
            pdbsites = pdbseqsite2refsite.keys()
            pdbsites.sort()
            for site in pdbsites:
                refsite = pdbseqsite2refsite[site] # get the reference site for this PDB site
                h = -1
                red = 255
                green = 255
                blue = 255
                if (refsite+1) in data:
                    h = data[refsite+1]
                    if h > 0:
                        #h = math.log(h)
                        red = 255
                        green = math.sqrt((1 - (h-mingain)/(maxh-mingain) )) * 255
                        blue =  math.sqrt((1 - (h-mingain)/(maxh-mingain) )) * 255
                    if h < 0:
                        #h = -1 * math.log( abs(h) )
                        blue = 255
                        red =  math.sqrt((1 -(maxloss-h)/(maxloss-minh))) * 255
                        green =  math.sqrt((1 -(maxloss-h)/(maxloss-minh))) * 255
                    if h < minh or h > maxh:
                        print "out of range", h
                        exit()
                else:
                    this_color = snow
                #print "237:", site, pdbseqsite2structsite[site], red, green, blue
                this_color = "colorsite" + site.__str__()
                fout.write("cmd.set_color(\"" + this_color + "\", (" + int(red).__str__() + "," + int(green).__str__() + "," + int(blue).__str__() + ") )\n")    
                fout.write("cmd.color(\"" + this_color + "\", \"" + ancname + " and resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n")
            
            # Add some style:
            fout.write("cmd.hide('everything')\n")
            fout.write("cmd.show('cartoon', 'all')\n")
            
            # Save a ray-traced PNG
            png_path = get_output_dir(ap) + "/pymol_ray.scores." + metric + "." + msa + "." + ancname + ".png"
            fout.write("cmd.png(\"" + png_path + "\", 500, 500, ray=True)\n")
            
            # Save the PyMOL session as an *.pse file
            pymolsessionpath = get_output_dir(ap) + "/pymol_script.scores." + metric + "." + msa + "." + ancname + ".pse"
            fout.write("cmd.save(\"" + pymolsessionpath + "\")\n")
            fout.close()
            
            pymol_script_paths.append(scriptpath )
            

    #
    # VISUALIZATION #2: Type 1/2/3/indel onto the PDB
    #
    for msanick in ap.params["msa_comparisons"]:
        # Make a custom alignment of this MSA to the PDB:

        this_seq = dat2seq( ap.params["msa_comparisons"][msanick][1] )
        # Run Needleman-Wunch on this sequence versus the PDB sequence
        (this_identity, this_score, this_align1, this_align2, this_symbol) = needle(this_seq, pdb_seq)   
        pdbsites = get_pdb_sites(ap.params["pdb"])
        pdbseqsite2structsite = {}
        pdbseqsite2mysite = {} # key = 1-based pdb site number, value = 1-based ref site number
        pdbsite = 0 # a counter to sites in the aligned-not-aligned sequence
        for i in range(0, align2.__len__()):
            if align2[i] != "-":
                structsite = pdbsites[ pdbsite ]
                pdbseqsite2structsite[pdbsite] = structsite
                pdbseqsite2mysite[ pdbsite ] = i
                pdbsite += 1

        that_ancpath = ap.params["msa_comparisons"][msanick][1]
        ancname = that_ancpath.split(".")[ that_ancpath.split(".").__len__()-2 ]
        ancname = ancname.split("/")[1]

        tokens = ap.params["pdb"].split("/")
        pdbmodelname = re.sub( ".pdb", "", tokens[ tokens.__len__()-1 ] )
        scriptpath = get_output_dir(ap) + "/pymol_script.types." + msanick +"." + pdbmodelname + ".p"        
        fout = open(scriptpath, "w")
        fout.write("cmd.load(\"" + os.path.abspath(ap.params["pdb"]) + "\")\n")
        fout.write("cmd.set_name(\"" + pdbmodelname + "\", \"" + ancname + "\")\n")

        pdbsites = pdbseqsite2mysite.keys()
        pdbsites.sort()
        for site in pdbsites:
            refsite = pdbseqsite2mysite[site] # get the reference site for this PDB site
            h = -1
            red = 255
            green = 255
            blue = 255
            if refsite in ap.params["redsites"][msanick]:
                red = 255
                green = 0
                blue = 0
            elif refsite in ap.params["orangesites"][msanick]:
                red = 255
                green = 196
                blue = 0
            elif refsite in ap.params["greensites"][msanick]:
                red = 168
                green = 210
                blue = 255
            else:
                red = 255
                green = 255
                blue = 255
            
            #print "237:", site, pdbseqsite2structsite[site], red, green, blue
            this_color = "colorsite" + site.__str__()
            fout.write("cmd.set_color(\"" + this_color + "\", (" + int(red).__str__() + "," + int(green).__str__() + "," + int(blue).__str__() + ") )\n")    
            fout.write("cmd.color(\"" + this_color + "\", \"" + ancname + " and resi " + (pdbseqsite2structsite[site]).__str__() + "\")\n")
            
        # Add some style:
        fout.write("cmd.hide('everything')\n")
        fout.write("cmd.show('cartoon', 'all')\n")
        
        # Save a ray-traced PNG
        png_path = get_output_dir(ap) + "/pymol_ray.types." + msanick + "." +ancname + ".png"
        fout.write("cmd.png(\"" + png_path + "\", 500, 500, ray=True)\n")
        
        # Save the PyMOL session as an *.pse file
        pymolsessionpath = get_output_dir(ap) + "/pymol_script.types." + msanick +"." + ancname + ".pse"
        fout.write("cmd.save(\"" + pymolsessionpath + "\")\n")
        fout.close()
        
        pymol_script_paths.append(scriptpath )
    
    #
    # Finally, invoke all the PyMOL scripts that have been accumulated.
    #
    for path in pymol_script_paths:
        os.system(ap.params["pymol_exe"] + " -c -u " + path)
    
