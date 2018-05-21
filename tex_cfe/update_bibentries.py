def update_bibentries():
    import os
    import ads
    import bibtexparser

    with open('extracted.bib', 'r') as fh:
        bib_database = bibtexparser.load(fh)

    with open('bibdesk.bib', 'r') as fh:
        full_bib_database = bibtexparser.load(fh)

    arxivs = []

    for entry in bib_database.entries:
        if 'doi' in entry:
            match = [x for x in full_bib_database.entries if 'doi' in x and x['doi'] == entry['doi']][0]
            #print(entry['journal'], match['journal'], match['eprint'] if 'eprint' in match else "")
        elif 'journal' in entry and entry['journal'] == 'ArXiv e-prints':
            match = [x for x in full_bib_database.entries if x['ID'] == entry['ID']][0]
            print(entry['journal'], match['journal'], match['eprint'] if 'eprint' in match else "")
            arxivs.append(match['eprint'])
        else:
            pass
            #print(entry['ID'], entry['author'], entry['year'], entry['journal'] if 'journal' in entry else "")
            #print(list(entry.keys()), entry['journal'] if 'journal' in entry else "", entry['volume'] if 'volume' in entry else '')

    # have to execut adsbibdesk + this
    print(" ".join(arxivs))

if __name__ == "__main__":
    assert os.system('bibexport -o extracted.bib sgrb2_cores.aux') == 0
    update_bibentries()
