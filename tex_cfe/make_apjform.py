with open('sgrb2_cores.tex','r') as fhr:
    with open('sgrb2_cores_apjform.tex', 'w') as fhw:
        for line in fhr:
            if 'figures/' in line:
                fhw.write(line.replace("figures/",""))
            else:
                fhw.write(line)
