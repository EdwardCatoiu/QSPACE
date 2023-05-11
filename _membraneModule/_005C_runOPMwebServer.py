# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

from .OPMserver import *


def run_005C_opmServer(query_structureIds ,
                       chromeDriverLocation,
                       pdbOutputs = {},
                       htmlOutputs = {},
                       not_solved_dict = {},
                       force_rerun = False,
                       to_opm_PDBfolder = qspaceDirs['opmStructuresToSendDir'] ,
                       from_opm_PDBfolder = qspaceDirs['opmOutputStructuresDir'] ,
                       errorsOPM = [],
                       
                       OPMdataFolder = qspaceDirs['opmOutputDataDir'],
                       
                      ):

    #using the already calculated inputs
    if pdbOutputs == {}:
        infile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-pdbOutputs.json')
        if op.exists(infile):
            with open(infile, 'r') as f:
                pdbOutputs = json.load(f)
        else:
            with open(infile, 'w') as f:
                json.dump(pdbOutputs,f)
            
    if htmlOutputs == {}:
        infile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-htmlOutputs.json')
        if op.exists(infile):
            with open(infile, 'r') as f:
                htmlOutputs = json.load(f)
        else:
            with open(infile, 'w') as f:
                json.dump(htmlOutputs,f)
                
    if not_solved_dict == {}:
        NotSolvedInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-not_solved_dict.json')
        if op.exists(NotSolvedInfile):
            with open(NotSolvedInfile, 'r') as f:
                not_solved_dict = json.load(f)
        else:
            with open(NotSolvedInfile, 'w') as f:
                json.dump(not_solved_dict,f)
    ########################################
    
    
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(chromeDriverLocation,chrome_options=chrome_options)

   

    # Figure out which structures Still need to be sent to OPM
    opm_needed = opmNeeded(query_structureIds,
                                     OPMstructureFolder = from_opm_PDBfolder,
                                     force_rerun = force_rerun)
        
    
    for s_id  in tqdm(opm_needed):
    
#         still_waiting = OPMserver.opmNeeded(query_structureIds,
#                                             OPMstructureFolder = from_opm_PDBfolder,
#                                             force_rerun = force_rerun
#                                            )
        
#         if s_id not in still_waiting and not force_rerun:
#             continue

            
        f = "{}.pdb".format(s_id)
        sfile = op.join(to_opm_PDBfolder,f)
        

        NotSolvedInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-not_solved_dict.json')
        if op.exists(NotSolvedInfile):
            with open(NotSolvedInfile, 'r') as f:
                not_solved_dict = json.load(f)

        if s_id in errorsOPM:
            continue

        if s_id in not_solved_dict: # this means that we are waiting for the calculations but we have already sent it to OPM
            continue


        print ('\n',s_id)

        driver.get("https://opm.phar.umich.edu/ppm_server2")
        time.sleep(5)


        #uploadfile
        chooseFile = driver.find_element_by_class_name('file-upload')
        chooseFile.send_keys(os.path.abspath(sfile))
        time.sleep(5)

        instructions = driver.find_elements_by_class_name('instructions')[0]
        submitButton = instructions.find_elements_by_class_name('submit-button')[-1]
        submitButton.click()

        time.sleep(10)
#         driver.save_screenshot("ppm_server2.png")

        try:
            resultsLink = driver.find_elements_by_link_text('open in new tab')[0].get_property('href')
        except IndexError:
            errorsOPM +=[s_id]
            print ('error^')
            continue

        ######### re check the results ##########
        driver.get(resultsLink)
        print (resultsLink )   

        time.sleep(4)
        resultsWaiting = True
        i = 0 
        not_solved = False
#         driver.save_screenshot("ppm_server2.png")

        while resultsWaiting:
            if i > 15:
                not_solved = True
                break

            try:
                i +=1
                print ('.', end=""),
                driver.execute_script("window.scrollTo(0, 0);")
#                 driver.save_screenshot("ppm_server2.png")

                pdbFileLink = driver.find_elements_by_link_text('PDB file')[0].get_property('href')
                resultsWaiting = False
#                 driver.save_screenshot("ppm_server2.png")

            except IndexError:
                time.sleep(6)


        if not_solved:
            #add the URL to the NotSolvedDictionary, these will be checked later....
            with open(NotSolvedInfile, 'r') as f:
                not_solved_dict = json.load(f)
            not_solved_dict.update({s_id:driver.current_url})
            with open(NotSolvedInfile, 'w') as f:
                json.dump(not_solved_dict,f)
            continue

        #scroll up
        driver.execute_script("window.scrollTo(0, 0);")
        time.sleep(1)
        #######################################     PDB

        time.sleep(2)
        r_pdb = requests.get(url = pdbFileLink)
        time.sleep(2)
        outFilePDB = "OPMstructure-{}".format(op.basename(pdbFileLink))
        with open(op.join(from_opm_PDBfolder, outFilePDB), 'w') as f:
            f.write(r_pdb.text)
            f.close()
#         print "Saving ... > {}\n".format(op.join(from_opm_PDBfolder, outFilePDB))
        #######################################     HTML
        time.sleep(2)
        driver.get(resultsLink)
        time.sleep(3)

        r_html = requests.get(url = driver.current_url)

        outFileHtml = "OPMhtml-{}.html".format(s_id)
        with open(op.join(OPMdataFolder, outFileHtml), 'w') as f:
            f.write(r_html.text)
            f.close()


        data = {}
        ##### CLEAN IT UP
        info = r_html.text
        t = info.replace('\n','')
        tt = t.replace('\r','')
        ttt = tt.replace('\t','')

        soup = BeautifulSoup(ttt, "lxml")
        soup = BeautifulSoup(ttt, "lxml")
        # Find all tables in the HTML code
        tables = soup.find_all("table", attrs={"class":"data"})

        table_orientation  = tables[0]
        table_embedded = tables[1]
        table_output = tables[3]
        orientation = scrape_orient(table_orientation)
        Embedded, TransMemSeg = scrape_residues(table_embedded)
        pdbLink = getPDBlink(table_output)

        data = {}
        data.update(orientation)
        data.update({'Embedded' : Embedded})
        data.update({"TransMem" : TransMemSeg})
        data.update({'PDB_link' : pdbLink})
        data.update({'PDB_Filelink' : pdbFileLink})
        data.update({'HTML_Filelink' : resultsLink})
        data = {s_id : data}


        outFileJSON = "OPMjson-{}.json".format(s_id)
        with open(op.join(OPMdataFolder, outFileJSON), 'w') as f:
            json.dump(data, f)

            
        PDBinfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-pdbOutputs.json')
        with open(PDBinfile, 'r') as f:
            pdbOutputs = json.load(f)
        pdbOutputs.update({s_id : pdbFileLink})
        with open(PDBinfile, 'w') as f:
            json.dump(pdbOutputs,f)
            
        htmlInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-htmlOutputs.json')
        with open(htmlInfile, 'r') as f:
            htmlOutputs = json.load(f)
        htmlOutputs.update({s_id : resultsLink})
        with open(htmlInfile, 'w') as f:
            json.dump(htmlOutputs,f)
                    
        time.sleep(3)
    driver.quit()
    return errorsOPM, not_solved_dict


def run_005C_goBackandCheck(not_solved_dict,
                            to_opm_folder = qspaceDirs['opmStructuresToSendDir'] ,
                            from_opm_PDBfolder = qspaceDirs['opmOutputStructuresDir'] ,
                            force_rerun = False,
                            OPMdataFolder = qspaceDirs['opmOutputDataDir'],
                            chromeDriverLocation = op.join(qspaceDirs['root_dir'], 'chromedriver') ,
                           ):
    
    tocheck = copy.deepcopy(not_solved_dict)
    
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(chromeDriverLocation,chrome_options=chrome_options)
    
    
    for s_id in tqdm(tocheck):
        try:
            resultsLink = not_solved_dict[s_id]
        except KeyError:
            continue
        print ('\n',s_id)
        f = s_id + '.pdb'
        sfile = op.join(to_opm_folder, f)

        outfilePDB = "OPMstructure-{}out.pdb".format(s_id.replace('-','_'))
        if op.exists(op.join(from_opm_PDBfolder, outfilePDB)) and not force_rerun:
            continue

        time.sleep(1)
        driver.get(resultsLink)
        print ('\n',resultsLink    )

        time.sleep(3)
        resultsWaiting = True
        i = 0 
        not_solved = False
#         driver.save_screenshot("ppm_server2.png")

        while resultsWaiting:
            if i > 5:
                not_solved = True
                break

            try:
                i +=1
                print ('.', end=""),
                driver.execute_script("window.scrollTo(0, 0);")
                driver.save_screenshot("ppm_server2.png")

                pdbFileLink = driver.find_elements_by_link_text('PDB file')[0].get_property('href')
                resultsWaiting = False
#                 driver.save_screenshot("ppm_server2.png")

            except IndexError:
                time.sleep(6)


        if not_solved:
            #add the URL to the NotSolvedDictionary, these will be checked later....
            NotSolvedInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-not_solved_dict.json')
            with open(NotSolvedInfile, 'r') as f:
                not_solved_dict = json.load(f)
            not_solved_dict.update({s_id:driver.current_url})
            with open(NotSolvedInfile, 'w') as f:
                json.dump(not_solved_dict,f)
            continue

        NotSolvedInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-not_solved_dict.json')
        with open(NotSolvedInfile, 'r') as f:
            not_solved_dict = json.load(f)    
        del not_solved_dict[s_id]
        with open(NotSolvedInfile, 'w') as f:
            json.dump(not_solved_dict,f)
                

        #scroll up
        driver.execute_script("window.scrollTo(0, 0);")
        time.sleep(1)
        #######################################     PDB

        time.sleep(2)
        r_pdb = requests.get(url = pdbFileLink)
        time.sleep(1.5)
        outFilePDB = "OPMstructure-{}".format(op.basename(pdbFileLink))
        with open(op.join(from_opm_PDBfolder, outFilePDB), 'w') as f:
            f.write(r_pdb.text)
            f.close()
#         print "Saving ... > {}\n".format(op.join(from_opm_PDBfolder, outFilePDB))

        #######################################     HTML
    #     time.sleep(2)
        driver.get(resultsLink)
        time.sleep(2)

        r_html = requests.get(url = driver.current_url)

        outFileHtml = "OPMhtml-{}.html".format(s_id)
        with open(op.join(OPMdataFolder, outFileHtml), 'w') as f:
            f.write(r_html.text)
            f.close()


        data = {}
        ##### CLEAN IT UP
        info = r_html.text
        t = info.replace('\n','')
        tt = t.replace('\r','')
        ttt = tt.replace('\t','')

        soup = BeautifulSoup(ttt, "lxml")
        soup = BeautifulSoup(ttt, "lxml")
        # Find all tables in the HTML code
        tables = soup.find_all("table", attrs={"class":"data"})

        table_orientation  = tables[0]
        table_embedded = tables[1]
        table_output = tables[3]
        orientation = scrape_orient(table_orientation)
        Embedded, TransMemSeg = scrape_residues(table_embedded)
        pdbLink = getPDBlink(table_output)

        data = {}
        data.update(orientation)
        data.update({'Embedded' : Embedded})
        data.update({"TransMem" : TransMemSeg})
        data.update({'PDB_link' : pdbLink})
        data.update({'PDB_Filelink' : pdbFileLink})
        data.update({'HTML_Filelink' : resultsLink})
        data = {s_id : data}


        outFileJSON = "OPMjson-{}.json".format(s_id)
        with open(op.join(OPMdataFolder, outFileJSON), 'w') as f:
            json.dump(data, f)


        PDBinfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-pdbOutputs.json')
        with open(PDBinfile, 'r') as f:
            pdbOutputs = json.load(f)
        pdbOutputs.update({s_id : pdbFileLink})
        with open(PDBinfile, 'w') as f:
            json.dump(pdbOutputs,f)
            
        htmlInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-htmlOutputs.json')
        with open(htmlInfile, 'r') as f:
            htmlOutputs = json.load(f)
        htmlOutputs.update({s_id : resultsLink})
        with open(htmlInfile, 'w') as f:
            json.dump(htmlOutputs,f)
            
        time.sleep(3)
    return not_solved_dict
    driver.quit()

