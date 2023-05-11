#these functions are copied from https://github.com/SBRG/ssbio/ and updated for python 3.7
import ssbio
# print (ssbio.__dict__)
import re
import os
import os.path as op

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest

import logging
log = logging.getLogger(__name__)


def is_valid_uniprot_id(instring):
    """Check if a string is a valid UniProt ID.

    See regex from: http://www.uniprot.org/help/accession_numbers

    Args:
        instring: any string identifier

    Returns: True if the string is a valid UniProt ID

    """
    valid_id = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    if valid_id.match(str(instring)):
        return True
    else:
        return False
    
    
def download_uniprot_file(uniprot_id, filetype, outdir='', force_rerun=False):
    """Download a UniProt file for a UniProt ID/ACC

    Args:
        uniprot_id: Valid UniProt ID
        filetype: txt, fasta, xml, rdf, or gff
        outdir: Directory to download the file

    Returns:
        str: Absolute path to file

    """
    if not is_valid_uniprot_id(uniprot_id):
        raise ValueError('Invalid UniProt ID')

    my_file = '{}.{}'.format(uniprot_id, filetype)
    url = 'http://www.uniprot.org/uniprot/{}'.format(my_file)
    outfile = op.join(outdir, my_file)

    if not op.exists(outfile) or force_rerun:
        urlrequest.urlretrieve(url, outfile)

    return outfile