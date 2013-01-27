from config import *
from anccomp_tools import *

PYMOL = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"


def do_pymol_viz(ap, data, seedseq):
    pdb_path = ap.getOptionalArg("--pdb_path")
    start_site = int(ap.getOptionalArg("--pdb_start_site"))
    if start_site == False:
        print "\n. You need to specify --pdb_start_site"
        exit()
        
    print "\n. Plotting H data on the PDB file", pdb_path
    
    print data
    print seedseq
    
    outlines = ""
    outlines += "cmd.load(\"" + os.getcwd() + "/" + pdb_path + "\")\n"
    pdb_site = 1
    ref_site = start_site
    
    maxh = data[1]
    minh = data[1]
    for site in data:
        if maxh < data[site]:
            maxh = data[site]
        if minh > data[site]:
            minh = data[site]
    
    while (ref_site < seedseq.__len__() ):
        if seedseq[ref_site] != "-":
            h = data[ref_site]
            if h > 0:
                red = 1.0#0.9 + (h/maxh)*0.1
                green = 1.0#0.9 + (h/maxh)*0.1
                blue = 1.0 - (h/maxh)
            if h < 0:
                blue = 1.0#0.9 + (h/minh)*0.1
                green = 1.0#0.9 + (h/minh)*0.1
                red = 1.0 - (h/minh)
            if h == 0:
                blue = 1.0#0.9
                green = 1.0#0.9
                red = 1.0#0.9
            print ref_site, seedseq[ref_site], h, blue, green, red
            this_color = "[" + red.__str__() + "," + green.__str__() + "," + blue.__str__() + "]"
            outlines += "cmd.set_color('color" + pdb_site.__str__() + "'," + this_color + ")\n"
            outlines += "cmd.color(\"color" + pdb_site.__str__() + "\", \"resi " + pdb_site.__str__() + "\")\n"
            pdb_site += 1
        ref_site += 1
            
    scriptpath = "pymol_script.p"
    fout = open(scriptpath, "w")
    fout.write(outlines)
    fout.close()
    
    print "\n\n. Run the script at " + os.getcwd() + "/" + scriptpath
    
    #os.system(PYMOL + " " + scriptpath)