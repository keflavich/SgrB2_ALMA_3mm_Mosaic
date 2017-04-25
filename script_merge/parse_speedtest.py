from datetime import datetime

def parse_logfile(logfilename):
    results = {}
    with open(logfilename, 'r') as fh:
        for line in fh:
            if 'Begin Task: tclean' in line:
                starttime = datetime.strptime(" ".join(line.split()[:2]),
                                              "%Y-%m-%d %H-%M-%S")
            elif 'imagename=' in line:
                imname = line.split('imagename="')[1].split('"')[0]
            elif 'End Task: tclean' in line:
                endtime = datetime.strptime(" ".join(line.split()[:2]),
                                            "%Y-%m-%d %H-%M-%S")
                results[imname] = (endtime-starttime).seconds
    return results
