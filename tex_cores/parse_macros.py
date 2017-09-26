import re

def find_closing_brace(lines):
    nbraces = 0
    first_brace_found = False
    for linenum,line in enumerate(lines):
        for c in line:
            if c == '{':
                nbraces += 1
            elif c == '}':
                nbraces -= 1
            elif c == "%":
                # Skip the rest of the line
                break # break the inner loop
            if nbraces == 0 and first_brace_found:
                return linenum
            elif nbraces == 1:
                first_brace_found=True

def find_matching_braces(lines, nmatch):
    nbraces = 0
    nmatches = 0
    halfmatch = False
    content = ""
    macargs = []
    for linenum,line in enumerate(lines):
        prevchar = ""
        for c in line:

            if c == '{':
                nbraces += 1
                if halfmatch:
                    content = content + "{"
                else:
                    halfmatch = True
                continue
            elif c == '}':
                nbraces -= 1
            elif c == "%" and prevchar != "\\":
                #content = content + "\n"
                # Skip the rest of the line
                break # break the inner loop
            if nbraces == 0 and halfmatch:
                halfmatch = False
                nmatches += 1
                macargs.append(content)
                content = ""
            elif halfmatch:
                content = content + c
            if nmatches == nmatch:
                return linenum, macargs
            prevchar = c

command = re.compile(r"^[^%]*\\def\\([A-Za-z]+)((#[0-9])*){")

#def strip_comments(lines):
#    return [line[:line.find("%")] for line in lines]

def find_macro(lines):
    for linenum,line in enumerate(lines):
        search = command.search(line)
        if search:
            start_num = linenum
            end_num = linenum+find_closing_brace(lines[linenum:])
            return search.group(1),start_num,end_num+1

def find_all_macros(lines):
    prev_macro = find_macro(lines)
    command_name, start_num, end_num = prev_macro
    macros = {command_name: (start_num,end_num)}
    placekeeper = 0
    while prev_macro:
        prev_macro = find_macro(lines[placekeeper+end_num:])
        if prev_macro is None:
            break
        placekeeper += end_num
        command_name, start_num, end_num = prev_macro
        macros[command_name]=(placekeeper+start_num,placekeeper+end_num)

    return macros

def parse_all_macros(lines):
    macros = {mac: parse_macro("".join(lines[start:end]))
              for mac,(start,end) in find_all_macros(lines).items()}
    return macros


def parse_macro(macro):

    cmd_name, input_ids, _ = command.search(macro).groups()
    input_ids = [input_ids[i:i+2] for i in range(0, len(input_ids), 2)]

    new_macro = command.sub("", macro)
    new_macro = re.sub("#([0-9])",
                       lambda m: "{{{0}}}".format(int(m.group(1))-1),
                       new_macro.replace("{","{{").replace("}","}}"))
    new_macro = new_macro[:new_macro.rfind("}")-1]
    return new_macro, len(input_ids)

def test():
    fig, nargs = parse_macro(
        r"""
\def\Figure#1#2#3#4#5{
\begin{figure*}[htp]
\includegraphics[scale=#4,angle=#5]{#1}
\caption{#2}
\label{#3}
\end{figure*}
}
        """)

    assert nargs == 5

    rslt = find_next_macro_in_use(r"""
\Figure{figures/Orion_SgrB2HII_side_by_side.png}
{blah}
{fig:orioncompare}{1}{\textwidth}
""".split("\n"),
                                  fig,
                                  'Figure',
                                  nargs)

    print(rslt)
    assert rslt is not None



    fig2,nargs = parse_macro(
        r"""
\def\FigureTwo#1#2#3#4#5#6{
\begin{figure*}[!htp]
\subfigure[]{ \includegraphics[scale=#5,width=#6]{#1} }
\subfigure[]{ \includegraphics[scale=#5,width=#6]{#2} }
\caption{#3}
\label{#4}
\end{figure*}
}
""")

    assert nargs == 6

    rslt = find_next_macro_in_use(r"""
\FigureTwo
{figures/continuum_peak_DeepSouth_saturated.pdf}
{figures/cores_on_continuum_peak_DeepSouth_saturated.pdf}
{blah }
{fig:continuumDS}{1}{0.5\textwidth}
                                  """.split("\n"),
                                  fig2,
                                  'FigureTwo',
                                  nargs)
    print(rslt)
    assert rslt is not None

    figonecol, nargs = parse_macro(
r"""
\def\FigureOneCol#1#2#3#4#5{
\begin{figure}[!htp]
\includegraphics[scale=#4,width=#5]{#1}
\caption{#2}
\label{#3}
\end{figure}
}
""")

    assert nargs == 5

    rslt = find_next_macro_in_use(
r"""
\FigureOneCol{figures/Orion_SgrB2HII_side_by_side.png}
{blah}
{fig:orioncompare}{1}{0.5\textwidth}
""".split("\n"),
                                  figonecol,
                                  'FigureOneCol',
                                  nargs)
    print(rslt)
    assert rslt is not None



#def find_macro_in_use(lines, macro, macroname, nargs):
#    command = re.compile(r'\\'+macroname)
#    substitutions = {}
#    for linenum,line in enumerate(lines):
#        if command.search(line):
#            end_num, macargs = find_matching_braces(lines[linenum:], nargs)
#            substitutions[(linenum, linenum+end_num-1)] = macargs
#
#    return substitutions
#
#def replace_macro(lines, macro, macroname, nargs):
#    substitutions = find_macro_in_use(lines, macro, macroname, nargs)
#    for (start,end),macargs in substitutions.items():
#        print("replacing {0} with {1}".format(lines[start:end+1],)
#                                              macro.format(*macargs).split("\n"))
#        
#        lines = lines[:start] + macro.format(*macargs).split("\n") + lines[end:]
#    return lines
#
#def replace_all_macros(lines, macros):
#    for macroname, (macro, nargs) in macros.items():
#        lines = replace_macro(lines, macro, macroname, nargs)
#
#    return lines

def find_next_macro_in_use(lines, macro, macroname, nargs):
    defcommand = re.compile(r'\\def')
    usecommand = re.compile(r'\\'+macroname+'($|[^A-Za-z])')
    for linenum,line in enumerate(lines):
        src = usecommand.search(line)
        if src:
            if defcommand.search(line):
                print("Skipped {0} on {2}:{1}".format(macroname, line, linenum))
                continue
            if line.find("%")>=0 and line.find("%") < src.start():
                print("COMMENT Skipped {0} on {2}:{1}".format(macroname, line, linenum))
                continue
            print("Matched {0} on {2}:{1}".format(macroname, line, linenum))
            end_num, macargs = find_matching_braces(lines[linenum:], nargs)
            #substitutions[(linenum, linenum+end_num-1)] = macargs
            return linenum, linenum+end_num, macargs

def replace_all_macros_onepass(lines, macros):
    for macroname, (macro, nargs) in macros.items():
        index = 0

        next_macro = find_next_macro_in_use(lines[index:], macro,
                                            macroname, nargs)
        while next_macro:
            start,end,macargs = next_macro
            print(index, start, end, macroname, macargs)
            newlines = [x+"\n" for x in macro.format(*macargs).split("\n")]

            lines = lines[:start+index] + newlines + lines[end+1+index:]

            index = start+index+len(newlines)

            next_macro = find_next_macro_in_use(lines[index:], macro,
                                                macroname, nargs)

    return lines

if __name__ == "__main__":
    import sys
    import subprocess

    test()
    
    assert sys.argv[0] == 'parse_macros.py'
    name_in = sys.argv[1]
    name_out = sys.argv[2]
    assert name_in[-4:] == '.tex'

    expanded_name = '{0}_expanded.tex'.format(name_in[:-4])

    with open(expanded_name, 'w') as f:
        result = subprocess.call(['latexpand', '--keep-comments', name_in],
                                 stdout=f)

    with open(expanded_name, 'r') as f:
        lines = f.readlines()
    macros = {k:v
              for k,v in parse_all_macros(lines).items()
              if 'Figure' in k or 'Table' in k}
    assert 'Figure' in macros
    newlines = replace_all_macros_onepass(lines, macros)

    for line in newlines:
        if line[-1] != "\n":
            import ipdb;ipdb.set_trace()

    with open(name_out, 'w') as f:
        f.writelines(newlines)
