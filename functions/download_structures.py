import os
import os.path as op
import zlib
import urllib
def download_cif(pdb_id,
                 outfolder ,
                 force_download = False,
                 pdb_checked = []):
    
    pdb_id = pdb_id.lower()
    try:
        if len(pdb_id) > 4:
            return False
        #print pdb_id,
        
        
        outfile_name = '%s.cif' %(pdb_id)
        outfile_path = op.join(outfolder, outfile_name)
        if op.exists(outfile_path) and not force_download:
            return outfile_path
        
        
        if pdb_id in pdb_checked and not force_download:
            if op.exists(outfile_path):
                return outfile_path
            else:
                return False
            
        folder = pdb_id[1:3]
#         server = 'https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/%s/' %(folder)
        server = 'https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/%s/' %(folder)
        f = urllib.request.urlopen(op.join(server , '%s.cif.gz' %(pdb_id)))
        decompressed_data=zlib.decompress(f.read(), 16+zlib.MAX_WBITS)
        with open(outfile_path , 'wb') as f:
            f.write(decompressed_data)
            f.close()
                              
    except (urllib.error.URLError) as e:
        return False
    
    return outfile_path


def download_bioassemblies(pdb_id,
                           outfolder ,
                           force_download = False,
                          ):
    pdb_id = pdb_id.lower()
    
    bioassembly_list = []

    folder = pdb_id[1:3]
    server = 'https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/%s/' %(folder)

    for i in range(1,15):
        outfile_name = '%s-assembly%i.cif' %(pdb_id, i)
        outfile_path = op.join(outfolder, outfile_name)
        if op.exists(outfile_path) and not force_download:
            bioassembly_list +=[outfile_path]
            continue

        try:

            f = urllib.request.urlopen(op.join(server , '%s-assembly%i.cif.gz' %(pdb_id,i)))
#                 print op.join(server , '%s-assembly%i.cif.gz' %(pdb_id,i))
            decompressed_data=zlib.decompress(f.read(), 16+zlib.MAX_WBITS)
#                 print 'decomp'
            with open(outfile_path , 'wb') as f:
                f.write(decompressed_data)
                f.close()
            bioassembly_list +=[outfile_path]
#                     print 'write'
#                     log.info(outfile_name)

        except (urllib.error.URLError) as e:
#                 print outfile_name,'<<--N/A'
            break
    
    return bioassembly_list