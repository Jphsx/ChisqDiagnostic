import glob
import sys


#specify dir which contains pngs from --plots output
plotdir = '/home/justin/work/research/susy/10-12-21/2Lpngs'

texfile='plots.tex'

pngfiles = glob.glob(plotdir+"/*.png")

#print(pngfiles)
prefitpng = [ x for x in pngfiles if "prefit" in x ]
postfitpng = [ x for x in pngfiles if "fit_b" in x ]
#print("PREFIT")
#print(prefitpng[0])

#print("POSTFIT")
#print(postfitpng)

prefitpng = sorted(prefitpng, key =lambda x: (x.split("_")[2], x.split("_")[1], x.split("_")[3]))

#open tex file to write to
f = open(texfile,"w")

#generate tex
f.write("\\documentclass[10pt,a4paper]{article}\n")
f.write("\\usepackage[utf8]{inputenc}\n")
f.write("\\usepackage{amsmath}\n")
f.write("\\usepackage{amsfonts}\n")
f.write("\\usepackage{amssymb}\n")
f.write("\\usepackage{graphicx}\n")
f.write("\\usepackage[margin=0.13in]{geometry}\n")
f.write("\\usepackage{float}\n")
f.write("\\usepackage{hyperref}\n")
f.write("\\begin{document}\n\n")
##body of tex
#/home/justin/work/research/susy/10-12-21/2Lpngs/Ch2L_mumu_bron_1jS_ge1jISR_PTISR0_CMS_th1x_prefit.png
#f.write("\\includegraphics[scale=0.5]{"+prefitpng[0]+"}\n")
#f.write("\\includegraphics[scale=0.5]{"+postfitpng[0]+"}\n")
for pre in prefitpng:
    #print region name
    RName = pre.split('/')
    RName = RName[-1]
    RName = RName.split("CMS")
    RName = RName[0]
    f.write("\\textbf{\\url{"+RName+"}} (Pre-fit,Post-fit) \\\\")
    f.write("\\includegraphics[scale=0.4]{"+pre+"}\n")
    #identify matching plot in post
    for post in postfitpng:
        if( RName in post ):
            f.write("\\includegraphics[scale=0.4]{"+post+"}\\\\\n")
            break

f.write("\end{document}\n")
f.close()
