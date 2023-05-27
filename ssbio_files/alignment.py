#replace the file :
#ssbio/ssbio/protein/sequence/utils/alignment.py  
# with this notebook

"""
Sequence Alignment
==================
"""

import logging
import os.path as op
import subprocess
import tempfile
from collections import defaultdict
from itertools import count, groupby
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import pairwise2
from Bio.Align import MultipleSeqAlignment
# from Bio.SubsMat import MatrixInfo as matlist
import ssbio.utils
import ssbio.protein.sequence.utils

# Quiet the SettingWithCopyWarning when converting dtypes in get_deletions/mutations methods
pd.options.mode.chained_assignment = None

log = logging.getLogger(__name__)


def pairwise_sequence_alignment(a_seq, b_seq, engine, a_seq_id=None, b_seq_id=None,
                                gapopen=10, gapextend=0.5,
                                outfile=None, outdir=None, force_rerun=False):
    """Run a global pairwise sequence alignment between two sequence strings.

    Args:
        a_seq (str, Seq, SeqRecord, SeqProp): Reference sequence
        b_seq (str, Seq, SeqRecord, SeqProp): Sequence to be aligned to reference
        engine (str): `biopython` or `needle` - which pairwise alignment program to use
        a_seq_id (str): Reference sequence ID. If not set, is "a_seq"
        b_seq_id (str): Sequence to be aligned ID. If not set, is "b_seq"
        gapopen (int): Only for `needle` - Gap open penalty is the score taken away when a gap is created
        gapextend (float): Only for `needle` - Gap extension penalty is added to the standard gap penalty for each 
            base or residue in the gap
        outfile (str): Only for `needle` - name of output file. If not set, is {id_a}_{id_b}_align.txt
        outdir (str): Only for `needle` - Path to output directory. Default is the current directory.
        force_rerun (bool): Only for `needle` - Default False, set to True if you want to rerun the alignment 
            if outfile exists.

    Returns:
        MultipleSeqAlignment: Biopython object to represent an alignment

    """
    engine = engine.lower()

    if engine not in ['biopython', 'needle']:
        raise ValueError('{}: invalid engine'.format(engine))

    if not a_seq_id:
        a_seq_id = 'a_seq'
    if not b_seq_id:
        b_seq_id = 'b_seq'

    a_seq = ssbio.protein.sequence.utils.cast_to_str(a_seq)
    b_seq = ssbio.protein.sequence.utils.cast_to_str(b_seq)

    if engine == 'biopython': # this does not work --Eddie
        from Bio.SubsMat import MatrixInfo as matlist #moved import inside function --Eddie
        # TODO: allow different matrices? needle uses blosum62 by default, how to change that?
        # TODO: how to define gap open/extend when using matrix in biopython global alignment?
        log.warning('Gap penalties not implemented in Biopython yet')

        blosum62 = matlist.blosum62
        alignments = pairwise2.align.globaldx(a_seq, b_seq, blosum62)  # TODO: add gap penalties
        best_alignment = alignments[0]

        a = ssbio.protein.sequence.utils.cast_to_seq_record(best_alignment[0], id=a_seq_id)
        b = ssbio.protein.sequence.utils.cast_to_seq_record(best_alignment[1], id=b_seq_id)
        alignment = MultipleSeqAlignment([a, b], annotations={"score": best_alignment[2],
                                                              "start": best_alignment[3],
                                                              "end"  : best_alignment[4]})
        alignment.annotations['percent_identity'] = get_percent_identity(best_alignment[0], best_alignment[1]) * 100

        return alignment

    if engine == 'needle':
        alignment_file = run_needle_alignment(seq_a=a_seq, seq_b=b_seq, gapopen=gapopen, gapextend=gapextend,
                                              write_outfile=True,  # Has to be true, AlignIO parses files on disk
                                              outdir=outdir, outfile=outfile, force_rerun=force_rerun)
        log.debug('Needle alignment at {}'.format(alignment_file))

        if not op.exists(alignment_file):
            raise ValueError('{}: needle alignment file does not exist'.format(alignment_file))

        # Use AlignIO to parse the needle alignment, alignments[0] is the first alignment (the only one in pairwise)
        # alignments = list(AlignIO.parse(alignment_file, "emboss"))
        # alignment = alignments[0]
        alignment = needle_statistics_alignio(alignment_file)

        # Rename the sequence IDs
        alignment[0].id = a_seq_id
        alignment[1].id = b_seq_id

        # # Add needle statistics as annotations in the alignment object
        # stats = needle_statistics(alignment_file)
        # alignment_ids = list(stats.keys())
        # if len(alignment_ids) > 1:
        #     raise ValueError('Needle alignment file contains more than one pairwise alignment')
        # needle_id = alignment_ids[0]
        # alignment.annotations['percent_identity'] = stats[needle_id]['percent_identity']
        # alignment.annotations['percent_similarity'] = stats[needle_id]['percent_similarity']
        # alignment.annotations['percent_gaps'] = stats[needle_id]['percent_gaps']
        # alignment.annotations['score'] = stats[needle_id]['score']

        return alignment


def run_needle_alignment(seq_a, seq_b, gapopen=10, gapextend=0.5, write_outfile=True,
                         outdir=None, outfile=None, force_rerun=False):
    """Run the needle alignment program for two strings and return the raw alignment result.

    More info:
    EMBOSS needle: http://www.bioinformatics.nl/cgi-bin/emboss/help/needle
    Biopython wrapper: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84
    Using strings as input: https://www.biostars.org/p/91124/

    Args:
        id_a: ID of reference sequence
        seq_a (str, Seq, SeqRecord): Reference sequence
        id_b: ID of sequence to be aligned
        seq_b (str, Seq, SeqRecord): String representation of sequence to be aligned
        gapopen: Gap open penalty is the score taken away when a gap is created
        gapextend: Gap extension penalty is added to the standard gap penalty for each base or residue in the gap
        outdir (str, optional): Path to output directory. Default is the current directory.
        outfile (str, optional): Name of output file. If not set, is {id_a}_{id_b}_align.txt
        force_rerun (bool): Default False, set to True if you want to rerun the alignment if outfile exists.

    Returns:
        str: Raw alignment result of the needle alignment in srspair format.

    """
    # TODO: check if needle is installed and raise error if not

    if not outdir:
        outdir = ''

    # TODO: rewrite using utils functions - does not report error if needle is not installed currently
    # TODO: rethink outdir/outfile, also if this should return the tempfile or just a file object or whatever
    if write_outfile:
        seq_a = ssbio.protein.sequence.utils.cast_to_str(seq_a)
        seq_b = ssbio.protein.sequence.utils.cast_to_str(seq_b)

        if not outfile:
            outfile = op.join(tempfile.gettempdir(), 'temp_alignment.needle')
        else:
            outfile = op.join(outdir, outfile)

        if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
            cmd = 'needle -outfile="{}" -asequence=asis::{} -bsequence=asis::{} -gapopen={} -gapextend={}'.format(
                outfile, seq_a, seq_b, gapopen, gapextend)
            command = subprocess.Popen(cmd,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       shell=True)
            out, err = command.communicate()

        return outfile

    else:
        seq_a = ssbio.protein.sequence.utils.cast_to_str(seq_a)
        seq_b = ssbio.protein.sequence.utils.cast_to_str(seq_b)

        cmd = 'needle -auto -stdout -asequence=asis::{} -bsequence=asis::{} -gapopen={} -gapextend={}'.format(seq_a, seq_b, gapopen, gapextend)
        command = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout = command.stdout.read()
        return stdout


def run_needle_alignment_on_files(id_a, faa_a, id_b, faa_b, gapopen=10, gapextend=0.5,
                                  outdir='', outfile='', force_rerun=False):
    """Run the needle alignment program for two fasta files and return the raw alignment result.

    More info:
    EMBOSS needle: http://www.bioinformatics.nl/cgi-bin/emboss/help/needle
    Biopython wrapper: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84

    Args:
        id_a: ID of reference sequence
        faa_a: File path to reference sequence
        id_b: ID of sequence to be aligned
        faa_b: File path to sequence to be aligned
        gapopen: Gap open penalty is the score taken away when a gap is created
        gapextend: Gap extension penalty is added to the standard gap penalty for each base or residue in the gap
        outdir (str, optional): Path to output directory. Default is the current directory.
        outfile (str, optional): Name of output file. If not set, is {id_a}_{id_b}_align.txt
        force_rerun (bool): Default False, set to True if you want to rerun the alignment if outfile exists.

    Returns:
        str: Raw alignment result of the needle alignment in srspair format.

    """

    # TODO: rewrite using utils functions so we can check for needle installation
    # # If you don't want to save the output file, just run the alignment and return the raw results
    # if not outfile and not outdir:
    #     needle_cline = NeedleCommandline(asequence=faa_a, bsequence=faa_b,
    #                                      gapopen=gapopen, gapextend=gapextend,
    #                                      stdout=True, auto=True)
    #     stdout, stderr = needle_cline()
    #     raw_alignment_text = stdout.decode('utf-8')

    # Make a default name if no outfile is set
    if not outfile:
        outfile = op.join(outdir, '{}_{}.needle'.format(id_a, id_b))
    else:
        outfile = op.join(outdir, outfile)

    # Check if the outfile already exists
    if op.exists(outfile) and not force_rerun:
        return outfile
    # If it doesn't exist, or force_rerun=True, run the alignment
    else:
        cmd = 'needle -outfile="{}" -asequence="{}" -bsequence="{}" -gapopen={} -gapextend={}'.format(outfile,
                                                                                                      faa_a,
                                                                                                      faa_b,
                                                                                                      gapopen,
                                                                                                      gapextend)
        command = subprocess.Popen(cmd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True)
        out, err = command.communicate()

    return outfile


def get_percent_identity(a_aln_seq, b_aln_seq):
    """Get the percent identity between two alignment strings"""

    if len(a_aln_seq) != len(b_aln_seq):
        raise ValueError('Sequence lengths not equal - was an alignment run?')

    count = 0
    gaps = 0
    for n in range(0, len(a_aln_seq)):
        if a_aln_seq[n] == b_aln_seq[n]:
            if a_aln_seq[n] != "-":
                count += 1
            else:
                gaps += 1

    return count / float((len(a_aln_seq) - gaps))


def get_alignment_df(a_aln_seq, b_aln_seq, a_seq_id=None, b_seq_id=None):
    """Summarize two alignment strings in a dataframe.

    Args:
        a_aln_seq (str): Aligned sequence string
        b_aln_seq (str): Aligned sequence string
        a_seq_id (str): Optional ID of a_seq
        b_seq_id (str): Optional ID of b_aln_seq

    Returns:
        DataFrame: a per-residue level annotation of the alignment

    """
    if len(a_aln_seq) != len(b_aln_seq):
        raise ValueError('Sequence lengths not equal - was an alignment run?')

    if not a_seq_id:
        a_seq_id = 'a_seq'
    if not b_seq_id:
        b_seq_id = 'b_seq'

    a_aln_seq = ssbio.protein.sequence.utils.cast_to_str(a_aln_seq)
    b_aln_seq = ssbio.protein.sequence.utils.cast_to_str(b_aln_seq)

    a_idx = 1
    b_idx = 1

    appender = []

    for i, (a,b) in enumerate(zip(a_aln_seq, b_aln_seq)):
        to_append = {}

        if a == b and a != '-' and b != '-':
            aa_flag = 'match'
        elif a != b and a == '-' and b != '-':
            aa_flag = 'insertion'
        elif a != b and a != '-' and b == '-':
            aa_flag = 'deletion'
        elif a != b and a != '-' and b == 'X':
            aa_flag = 'unresolved'
        elif a != b and b != '-' and a == 'X':
            aa_flag = 'unresolved'
        elif a != b and a != '-' and b != '-':
            aa_flag = 'mutation'

        to_append['id_a'] = a_seq_id
        to_append['id_b'] = b_seq_id
        to_append['type'] = aa_flag

        if aa_flag == 'match' or aa_flag == 'unresolved' or aa_flag == 'mutation':
            to_append['id_a_aa'] = a
            to_append['id_a_pos'] = int(a_idx)
            to_append['id_b_aa'] = b
            to_append['id_b_pos'] = int(b_idx)
            a_idx += 1
            b_idx += 1

        if aa_flag == 'deletion':
            to_append['id_a_aa'] = a
            to_append['id_a_pos'] = int(a_idx)
            a_idx += 1

        if aa_flag == 'insertion':
            to_append['id_b_aa'] = b
            to_append['id_b_pos'] = int(b_idx)
            b_idx += 1

        appender.append(to_append)

    cols = ['id_a', 'id_b', 'type', 'id_a_aa', 'id_a_pos', 'id_b_aa', 'id_b_pos']
    alignment_df = pd.DataFrame.from_records(appender, columns=cols)
    alignment_df = alignment_df.fillna(value=np.nan)

    return alignment_df


def get_alignment_df_from_file(alignment_file, a_seq_id=None, b_seq_id=None):
    """Get a Pandas DataFrame of the Needle alignment results. Contains all positions of the sequences.

    Args:
        alignment_file:
        a_seq_id: Optional specification of the ID of the reference sequence
        b_seq_id: Optional specification of the ID of the aligned sequence

    Returns:
        Pandas DataFrame: all positions in the alignment

    """
    alignments = list(AlignIO.parse(alignment_file, "emboss"))
    alignment_df = pd.DataFrame(columns=['id_a', 'id_b', 'type', 'id_a_aa', 'id_a_pos', 'id_b_aa', 'id_b_pos'])

    for alignment in alignments:
        if not a_seq_id:
            a_seq_id = list(alignment)[0].id
        a_seq = str(list(alignment)[0].seq)
        if not b_seq_id:
            b_seq_id = list(alignment)[1].id
        b_seq = str(list(alignment)[1].seq)

        df = get_alignment_df(a_seq, b_seq, a_seq_id, b_seq_id)
        alignment_df = alignment_df.append(df).reset_index(drop=True)

    return alignment_df


def get_mutations(aln_df):
    """Get a list of residue numbers (in the original sequence's numbering) that are mutated

    Args:
        aln_df (DataFrame): Alignment DataFrame
        just_resnums: If only the residue numbers should be returned, instead of a list of tuples of
            (original_residue, resnum, mutated_residue)

    Returns:
        list: Residue mutations

    """
    mutation_df = aln_df[aln_df['type'] == 'mutation']
    tuples = []
    if not mutation_df.empty:
        subset = mutation_df[['id_a_aa', 'id_a_pos', 'id_b_aa']]
        subset['id_a_pos'] = subset['id_a_pos'].astype(int)
        tuples = [tuple(x) for x in subset.values]
    return tuples


def get_unresolved(aln_df):
    """Get a list of residue numbers (in the original sequence's numbering) that are unresolved

    Args:
        aln_df (DataFrame): Alignment DataFrame

    Returns:
        list: Residue numbers that are mutated

    """
    unresolved_df = aln_df[aln_df['type'] == 'unresolved']
    unresolved = []
    if not unresolved_df.empty:
        unresolved_df['id_a_pos'] = unresolved_df['id_a_pos'].astype(int)
        unresolved = unresolved_df.id_a_pos.tolist()
    return unresolved


def get_deletions(aln_df):
    """Get a list of tuples indicating the first and last residues of a deletion region, as well as the length of the deletion.

    Examples:
        # Deletion of residues 1 to 4, length 4
        >>> test = {'id_a': {0: 'a', 1: 'a', 2: 'a', 3: 'a'}, 'id_a_aa': {0: 'M', 1: 'G', 2: 'I', 3: 'T'}, 'id_a_pos': {0: 1.0, 1: 2.0, 2: 3.0, 3: 4.0}, 'id_b': {0: 'b', 1: 'b', 2: 'b', 3: 'b'}, 'id_b_aa': {0: np.nan, 1: np.nan, 2: np.nan, 3: np.nan}, 'id_b_pos': {0: np.nan, 1: np.nan, 2: np.nan, 3: np.nan}, 'type': {0: 'deletion', 1: 'deletion', 2: 'deletion', 3: 'deletion'}}
        >>> my_alignment = pd.DataFrame.from_dict(test)
        >>> get_deletions(my_alignment)
        [((1.0, 4.0), 4)]

    Args:
        aln_df (DataFrame): Alignment DataFrame

    Returns:
        list: A list of tuples with the format ((deletion_start_resnum, deletion_end_resnum), deletion_length)

    """

    deletion_df = aln_df[aln_df['type'] == 'deletion']
    if not deletion_df.empty:
        deletion_df['id_a_pos'] = deletion_df['id_a_pos'].astype(int)
    deletions = []

    for k, g in groupby(deletion_df.index, key=lambda n, c=count(): n - next(c)):
        tmp = list(g)
        deletion_indices = (min(tmp), max(tmp))

        deletion_start_ix = deletion_indices[0]
        deletion_end_ix = deletion_indices[1]

        deletion_length = deletion_end_ix - deletion_start_ix + 1

        id_a_pos_deletion_start = aln_df.ix[deletion_start_ix].id_a_pos
        id_a_pos_deletion_end = aln_df.ix[deletion_end_ix].id_a_pos

        deletion_region = (id_a_pos_deletion_start, id_a_pos_deletion_end)

        # Logging where the insertion is
        log.debug('Deletion of length {} at residues {}'.format(deletion_length, deletion_region))

        to_append = (deletion_region, deletion_length)
        deletions.append(to_append)

    return deletions


def get_insertions(aln_df):
    """Get a list of tuples indicating the first and last residues of a insertion region, as well as the length of the insertion.

    If the first tuple is:
        (-1, 1) that means the insertion is at the beginning of the original protein
        (X, Inf) where X is the length of the original protein, that means the insertion is at the end of the protein

    Examples:
        # Insertion at beginning, length 3
        >>> test = {'id_a': {0: 'a', 1: 'a', 2: 'a', 3: 'a'}, 'id_a_aa': {0: np.nan, 1: np.nan, 2: np.nan, 3: 'M'}, 'id_a_pos': {0: np.nan, 1: np.nan, 2: np.nan, 3: 1.0}, 'id_b': {0: 'b', 1: 'b', 2: 'b', 3: 'b'}, 'id_b_aa': {0: 'M', 1: 'M', 2: 'L', 3: 'M'}, 'id_b_pos': {0: 1, 1: 2, 2: 3, 3: 4}, 'type': {0: 'insertion', 1: 'insertion', 2: 'insertion', 3: 'match'}}
        >>> my_alignment = pd.DataFrame.from_dict(test)
        >>> get_insertions(my_alignment)
        [((-1, 1.0), 3)]

    Args:
        aln_df (DataFrame): Alignment DataFrame

    Returns:
        list: A list of tuples with the format ((insertion_start_resnum, insertion_end_resnum), insertion_length)

    """

    insertion_df = aln_df[aln_df['type'] == 'insertion']
    # if not insertion_df.empty: # don't need to do this for insertions
    #     insertion_df['id_a_pos'] = insertion_df['id_a_pos'].astype(int)

    insertions = []

    for k, g in groupby(insertion_df.index, key=lambda n, c=count(): n - next(c)):
        tmp = list(g)
        insertion_indices = (min(tmp), max(tmp))

        insertion_start = insertion_indices[0] - 1
        insertion_end = insertion_indices[1] + 1

        # Checking if insertion is at the beginning or end
        if insertion_start < 0:
            insertion_start = insertion_indices[0]
            insertion_length = insertion_end - insertion_start
        elif insertion_end >= len(aln_df):
            insertion_end = insertion_indices[1]
            insertion_length = insertion_end - insertion_start
        else:
            insertion_length = insertion_end - insertion_start - 1

        id_a_pos_insertion_start = aln_df.ix[insertion_start].id_a_pos
        id_a_pos_insertion_end = aln_df.ix[insertion_end].id_a_pos

        # Checking if insertion is at the beginning or end
        if np.isnan(id_a_pos_insertion_start) and id_a_pos_insertion_end == 1:
            insertion_region = (-1, id_a_pos_insertion_end)
        elif np.isnan(id_a_pos_insertion_end):
            insertion_region = (id_a_pos_insertion_start, float('Inf'))
        else:
            insertion_region = (id_a_pos_insertion_start, id_a_pos_insertion_end)

        # Logging where the insertion is
        if insertion_region[0] == -1:
            log.debug('Insertion of length {} at beginning'.format(insertion_length))
        elif insertion_region[1] == float('Inf'):
            log.debug('Insertion of length {} at end'.format(insertion_length))
        else:
            log.debug('Insertion of length {} at residues {}'.format(insertion_length, insertion_region))

        to_append = (insertion_region, insertion_length)
        insertions.append(to_append)

    return insertions


def map_resnum_a_to_resnum_b(resnums, a_aln, b_aln):
    """Map a residue number in a sequence to the corresponding residue number in an aligned sequence.

    Examples:
    >>> map_resnum_a_to_resnum_b([1,2,3], '--ABCDEF', 'XXABCDEF')
    {1: 3, 2: 4, 3: 5}
    >>> map_resnum_a_to_resnum_b(5, '--ABCDEF', 'XXABCDEF')
    {5: 7}
    >>> map_resnum_a_to_resnum_b(5, 'ABCDEF', 'ABCD--')
    {}
    >>> map_resnum_a_to_resnum_b(5, 'ABCDEF--', 'ABCD--GH')
    {}
    >>> map_resnum_a_to_resnum_b([9,10], '--MKCDLHRLE-E', 'VSNEYSFEGYKLD')
    {9: 11, 10: 13}

    Args:
        resnums (int, list): Residue number or numbers in the first aligned sequence
        a_aln (str, Seq, SeqRecord): Aligned sequence string
        b_aln (str, Seq, SeqRecord): Aligned sequence string

    Returns:
        int: Residue number in the second aligned sequence
    """
    resnums = ssbio.utils.force_list(resnums)

    aln_df = get_alignment_df(a_aln, b_aln)

    maps = aln_df[aln_df.id_a_pos.isin(resnums)]

    able_to_map_to_b = maps[pd.notnull(maps.id_b_pos)]
    successful_map_from_a = able_to_map_to_b.id_a_pos.values.tolist()

    mapping = dict([(int(a), int(b)) for a,b in zip(able_to_map_to_b.id_a_pos, able_to_map_to_b.id_b_pos)])

    cant_map = list(set(resnums).difference(successful_map_from_a))
    if len(cant_map) > 0:
        log.warning('Unable to map residue numbers {} in first sequence to second'.format(cant_map))

    return mapping


def pairwise_alignment_stats(reference_seq_aln, other_seq_aln):
    """Get a report of a pairwise alignment.

    Args:
        reference_seq_aln (str, Seq, SeqRecord): Reference sequence, alignment form
        other_seq_aln (str, Seq, SeqRecord): Other sequence, alignment form

    Returns:
        dict: Dictionary of information on mutations, insertions, sequence identity, etc.

    """
    if len(reference_seq_aln) != len(other_seq_aln):
        raise ValueError('Sequence lengths not equal - was an alignment run?')

    reference_seq_aln = ssbio.protein.sequence.utils.cast_to_str(reference_seq_aln)
    other_seq_aln = ssbio.protein.sequence.utils.cast_to_str(other_seq_aln)

    infodict = {}

    # Percent identity to the reference sequence
    stats_percent_ident = get_percent_identity(a_aln_seq=reference_seq_aln, b_aln_seq=other_seq_aln)
    infodict['percent_identity'] = stats_percent_ident

    # Other alignment results
    aln_df = get_alignment_df(a_aln_seq=reference_seq_aln, b_aln_seq=other_seq_aln)
    infodict['deletions'] = get_deletions(aln_df)
    infodict['insertions'] = get_insertions(aln_df)
    infodict['mutations'] = get_mutations(aln_df)
    infodict['unresolved'] = get_unresolved(aln_df)

    return infodict


def needle_statistics(infile):
    """Reads in a needle alignment file and spits out statistics of the alignment.

    Args:
        infile (str): Alignment file name

    Returns:
        dict: alignment_properties - a dictionary telling you the number of gaps, identity, etc.

    """

    alignments = list(AlignIO.parse(infile, "emboss"))
    alignment_properties = defaultdict(dict)

    with open(infile) as f:
        line = f.readline()

        for i in range(len(alignments)):
            while line.rstrip() != "#=======================================":
                line = f.readline()
                if not line:
                    raise StopIteration

            while line[0] == "#":
                # Read in the rest of this alignment header,
                # try and discover the number of records expected and their length
                parts = line[1:].split(":", 1)
                key = parts[0].lower().strip()
                if key == '1':
                    a_id = parts[1].strip()
                if key == '2':
                    b_id = parts[1].strip()
                if key == 'identity':
                    ident_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    ident_num = int(ident_parse[0].split('/')[0])
                    ident_percent = float(ident_parse[1])
                    alignment_properties[a_id + '_' + b_id]['identity'] = ident_num
                    alignment_properties[a_id + '_' + b_id]['percent_identity'] = ident_percent
                if key == 'similarity':
                    sim_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    sim_num = int(sim_parse[0].split('/')[0])
                    sim_percent = float(sim_parse[1])
                    alignment_properties[a_id + '_' + b_id]['similarity'] = sim_num
                    alignment_properties[a_id + '_' + b_id]['percent_similarity'] = sim_percent
                if key == 'gaps':
                    gap_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    gap_num = int(gap_parse[0].split('/')[0])
                    gap_percent = float(gap_parse[1])
                    alignment_properties[a_id + '_' + b_id]['gaps'] = gap_num
                    alignment_properties[a_id + '_' + b_id]['percent_gaps'] = gap_percent
                if key == 'score':
                    score = float(parts[1].strip())
                    alignment_properties[a_id + '_' + b_id]['score'] = score

                # And read in another line...
                line = f.readline()

    return alignment_properties


def needle_statistics_alignio(infile):
    """Reads in a needle alignment file and returns an AlignIO object with annotations

    Args:
        infile (str): Alignment file name

    Returns:
        AlignIO: annotated AlignIO object

    """

    alignments = list(AlignIO.parse(infile, "emboss"))

    if len(alignments) > 1:
        raise ValueError('Alignment file contains more than one pairwise alignment')

    alignment = alignments[0]

    with open(infile) as f:
        line = f.readline()

        for i in range(len(alignments)):
            while line.rstrip() != "#=======================================":
                line = f.readline()
                if not line:
                    raise StopIteration

            while line[0] == "#":
                # Read in the rest of this alignment header,
                # try and discover the number of records expected and their length
                parts = line[1:].split(":", 1)
                key = parts[0].lower().strip()
                if key == 'identity':
                    ident_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    ident_num = int(ident_parse[0].split('/')[0])
                    ident_percent = float(ident_parse[1])
                    alignment.annotations['identity'] = ident_num
                    alignment.annotations['percent_identity'] = ident_percent
                if key == 'similarity':
                    sim_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    sim_num = int(sim_parse[0].split('/')[0])
                    sim_percent = float(sim_parse[1])
                    alignment.annotations['similarity'] = sim_num
                    alignment.annotations['percent_similarity'] = sim_percent
                if key == 'gaps':
                    gap_parse = parts[1].strip().replace('(','').replace(')','').replace('%','').split()
                    gap_num = int(gap_parse[0].split('/')[0])
                    gap_percent = float(gap_parse[1])
                    alignment.annotations['gaps'] = gap_num
                    alignment.annotations['percent_gaps'] = gap_percent
                if key == 'score':
                    score = float(parts[1].strip())
                    alignment.annotations['score'] = score

                # And read in another line...
                line = f.readline()

    return alignment