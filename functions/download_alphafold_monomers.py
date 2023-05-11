from urllib.error import  HTTPError, URLError
from urllib.request import Request, urlopen
import os.path as op

def download_alphafold_structure(url,outfolder=False,force_rerun = False):
    if not outfolder:
        outfolder = '../GEMPRO/structures/all_alphafold/'
    
    
    file_name = op.basename(url)
    filepath = op.join(outfolder, file_name)
    if op.exists(filepath) and not force_rerun:
        return

    req = Request(url)
    # Open the url
    try:
        f = urlopen(req)
#         print "downloading " + url

        # Open our local file for writing
        with open(filepath, "wb") as local_file:
            local_file.write(f.read())

    #handle errors
    except (HTTPError) as  e:
#         print "HTTP Error:",e.code , url
        pass
    except (URLError) as  e:
#         print "URL Error:",e.reason , url
        pass

def download_alphafold_metrics(url,outfolder=False,force_rerun = False):
    if not outfolder:
        outfolder = '../GEMPRO/structures/all_alphafold-metrics/'
    
    file_name = op.basename(url)
    filepath = op.join(outfolder, file_name)
    if op.exists(filepath) and not force_rerun:
        return
    
    req = Request(url)
    # Open the url
    try:
        f = urlopen(req)
#         print "downloading " + url


        # Open our local file for writing
        with open(filepath, "wb") as local_file:
            local_file.write(f.read())


    #handle errors
    except (HTTPError) as  e:
#         print "HTTP Error:",e.code , url
        pass
    except (URLError) as  e:
#         print "URL Error:",e.reason , url
        pass