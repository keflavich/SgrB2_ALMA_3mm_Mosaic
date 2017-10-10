#!/bin/env python
import subprocess
import shutil
import glob
import argparse
import os
import update_bibentries

assert os.path.exists('bibdesk.bib')

name = 'sgrb2_cores'

parser = argparse.ArgumentParser(description='Make latex files.')
parser.add_argument('--referee', default=False,
                    action='store_true', help='referee style?')
parser.add_argument('--texpath', default='/usr/texbin/',
                    help='path to pdflatex')
parser.add_argument('--infile', default=name+'.tex')
parser.add_argument('--outfile', default='auto')
parser.add_argument('--clean', default=False, action='store_true')
parser.add_argument('--both', default=False, action='store_true')

args = parser.parse_args()


def do_everything():
    if args.clean:
        for globstr in (name+"*.aux", name+"*.bbl", name+"*.blg",
                        name+"*.dvi", name+"*.log", name+"*.lot",
                        name+"*.lof"):
            for fn in glob.glob(globstr):
                os.remove(fn)

    PDFLATEX=os.path.join(args.texpath,'pdflatex')
    pdflatex_args = "-halt-on-error -synctex=1 --interaction=nonstopmode".split()

    BIBTEX = os.path.join(args.texpath, 'bibtex')

    with open('preface.tex','r') as f:
        preface = f.read()

    with open('preface_aa.tex','w') as aa:
        if args.referee:
            aa.write('\documentclass[referee]{aa}\n')
            aa.write(preface)
        else:
            aa.write('\documentclass{aa}\n')
            aa.write(preface)

    with open('solobib.tex','w') as fh:
        fh.write("\\bibliographystyle{aasjournal}\n")
        fh.write("\\bibliography{bibdesk}")

    pdfcmd = [PDFLATEX] + pdflatex_args + [args.infile]
    bibcmd = [BIBTEX, args.infile.replace(".tex","")]

    print("Executing PDF command: {0}".format(pdfcmd))
    subprocess.call(pdfcmd)
    print("Executing bibtex command: {0}".format(bibcmd))
    bibresult = subprocess.call(bibcmd)
    assert bibresult == 0
    print("Executing PDF command a second time: {0}".format(pdfcmd))
    pdfresult = subprocess.call(pdfcmd)
    assert pdfresult == 0

    if not os.system('bibexport')==0:
        print("bibexport is not installed.")
    else:
        assert os.system('bibexport -o extracted.bib sgrb2_cores.aux') == 0
        print("bibexport created extracted.bib")
        
        try:
            update_bibentries.update_bibentries()
            print("Successfully updated bibentries")
        except ImportError:
            print("Could not update bibliography entries because of import error")

    with open('solobib.tex','w') as fh:
        fh.write("\\bibliographystyle{aasjournal}\n")
        fh.write("\\bibliography{extracted}")
    print("Created solobib.tex")


    print("Executing bibtex command a second time: {0}".format(bibcmd))
    bibresult = subprocess.call(bibcmd)
    assert bibresult == 0
    print("Executing PDF command a third time: {0}".format(pdfcmd))
    pdfresult = subprocess.call(pdfcmd)
    assert pdfresult == 0

    if args.outfile == 'auto':
        outprefix = name+'_referee' if args.referee else name
    else:
        outprefix = os.path.splitext(args.outfile)[0]

    # Don't move unnecessarily; messes with Skim.app (maybe?)
    if os.path.split(os.path.basename(outprefix))[0] != name:
        shutil.move(name+".pdf",outprefix+".pdf")


    gscmd = ["gs",
             "-dSAFER",
             "-dBATCH",
             "-dNOPAUSE",
             "-sDEVICE=pdfwrite",
             "-sOutputFile={0}_compressed.pdf".format(outprefix),
             "{0}.pdf".format(outprefix)]

    subprocess.call(gscmd)

if args.both:
    args.referee = True
    do_everything()
    args.referee = False
    do_everything()
else:
    do_everything()
