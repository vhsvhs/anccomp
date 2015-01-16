from config import *
from anccomp_tools import *

def write_table1(ap, metric_data):
    for msanick in metric_data[ metric_data.keys()[0] ].keys():
        write_table1_helper(ap, msanick, metric_data)

def write_table1_helper(ap, msanick, metric_data):
    """Returns one big string containing HTML for an interactive table that
    can be embedded in other pages."""
    
    metrics = metric_data.keys()
    sites = metric_data[metrics[0]][msanick].keys()
    sites.sort()
    
    outpath = get_table_outpath(ap, tag=msanick + ".html")
    
    this_ancpath = ap.params["msa_comparisons"][msanick][0]
    that_ancpath = ap.params["msa_comparisons"][msanick][1]
    
    fout = open(outpath, "w")
    fout.write("<table class=\"sortable\">\n")

    allvalues = []
    for site in sites:
        allvalues.append(metric_data["Df"][msanick][site])

    sdev = sd(allvalues)
    avg = mean(allvalues)

    fout.write("<tr class=\"headerrow\"><td>Site</td>")
    for m in metrics:
        fout.write("<td>" + m + "</td>")
    fout.write("<td>Posterior Probabilities</td>")
    fout.write("</tr>")
    for site in sites:        
        value = metric_data["Df"][msanick][site]
        
        # if this site is +2 s.d., then color the row
        style = "whiterow"
        if value < avg:
            if value <= avg - 6*sdev:
                style = "redrow"
            elif value <= avg - 4*sdev:
                style = "orangerow"
            elif value <= avg - 2*sdev:
                style = "bluerow"
            else:
                style = "whiterow" 
        elif value > avg:
            if value >= avg + 6*sdev:
                style = "redrow"
            elif value >= avg + 4*sdev:
                style = "orangerow"
            elif value >= avg + 2*sdev:
                style = "bluerow"
            else:
                style = "whiterow" 
        
        
        fout.write("\n<tr class=\"" + style + "\">\n")
        fout.write("<td>" + site.__str__() + "</td>\n")
        for m in metrics:
            fout.write("<td>" + "%.3f"%metric_data[m][msanick][site] + "</td>\n")
        fout.write("<td>")
        fout.write(htmlfrags[site])
        fout.write("</td>\n")
        

        #fout.write("<td>" + site + "</td>\n")
        #fout.write("<td>" + site + "</td>\n")
        #fout.write("<td>" + site + "</td>\n")
        #fout.write("<td>" + site + "</td>\n")
        #fout.write("<td>" + site + "</td>\n")
        #fout.write("<td>" + site + "</td>\n")
    
        fout.write("</tr>\n")
     
    fout.write("</table>\n")

def get_td_fragments(anc1, anc2):
    fin1 = open(anc1, "r")
    fin1.close()
    