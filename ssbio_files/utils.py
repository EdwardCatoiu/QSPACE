
# replace the file :
# ssbio/ssbio/protein/sequence/utils/utils.py
# with this notebook


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC


def cast_to_str(obj):
    """Return a string representation of a Seq or SeqRecord.

    Args:
        obj (str, Seq, SeqRecord): Biopython Seq or SeqRecord

    Returns:
        str: String representation of the sequence

    """

    if isinstance(obj, str):
        return obj
    if isinstance(obj, Seq):
        return str(obj)
    if isinstance(obj, SeqRecord):
        return str(obj.seq)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def cast_to_seq(obj):
    """Return a Seq representation of a string or SeqRecord object.

    Args:
        obj (str, Seq, SeqRecord): Sequence string or Biopython SeqRecord object
        alphabet: See Biopython SeqRecord docs -- no longer needed 

    Returns:
        Seq: Seq representation of the sequence

    """

    if isinstance(obj, Seq):
        return obj
    if isinstance(obj, SeqRecord):
        return obj.seq
    if isinstance(obj, str):
        obj = obj.upper()
        return Seq(obj)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def cast_to_seq_record(obj,  id="<unknown id>", name="<unknown name>",
                       description="<unknown description>", dbxrefs=None,
                       features=None, annotations=None,
                       letter_annotations=None):
    """Return a SeqRecord representation of a string or Seq object.

    Args:
        obj (str, Seq, SeqRecord): Sequence string or Biopython Seq object
        alphabet: See Biopython SeqRecord docs  #no longer needed Edited out by Edward Catoiu
        id: See Biopython SeqRecord docs
        name: See Biopython SeqRecord docs
        description: See Biopython SeqRecord docs
        dbxrefs: See Biopython SeqRecord docs
        features: See Biopython SeqRecord docs
        annotations: See Biopython SeqRecord docs
        letter_annotations: See Biopython SeqRecord docs

    Returns:
        SeqRecord: SeqRecord representation of the sequence

    """

    if isinstance(obj, SeqRecord):
        return obj
    if isinstance(obj, Seq):
        return SeqRecord(obj, id, name, description, dbxrefs, features, annotations, letter_annotations)
    if isinstance(obj, str):
        obj = obj.upper()
        return SeqRecord(Seq(obj), id, name, description, dbxrefs, features, annotations, letter_annotations)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')