import argparse
import time
from functools import wraps
from typing import List, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pro2codon import pn2codon
import sys


def write_result_to_fasta_file(
    file_name: str,
    seq_header_array: List[Tuple[str, str]],
):
    
    list_str = list(sum(seq_header_array, ()))

    for i in range(0, len(list_str), 2):
        list_str[i] = f">{list_str[i]}"

    with open(file_name, "w") as fw:
        fw.writelines("\n".join(list_str))

    print(f"List of reverse-translated NTs written to: {file_name}")

    return None


def read_and_convert_fasta_files(
    aa_file: str,
    nt_file: str,
) -> Tuple[
    List[Tuple[str, str]],
    List[Tuple[str, str]]
]:

    nt_seqs = SimpleFastaParser(open(nt_file))
    aa_seqs = SimpleFastaParser(open(aa_file))

    aas = []
    nts = {}
    nt_headers = {}

    for i, (aa, nt) in enumerate(zip(aa_seqs, nt_seqs)):
        aa_header, aa_seq = aa
        nt_header, nt_seq = nt

        nt_headers[i] = aa_header
        nts[nt_header] = (nt_header, nt_seq)
        aas.append((aa_header, aa_seq))


    nts_in_order = [None for _ in range(len(nts))]

    for i, correct_header in nt_headers.items():
        try:
            nts_in_order[i] = nts[correct_header]
        except:
            print("ERROR CAUGHT: There is a single header in PEP sequence FASTA file that does not exist in NUC sequence FASTA file")
            sys.exit()

    return (aas, nts_in_order)


def convert_to_codon(
    filepath: str,
    table_index: int,
    aa_seqs: List[Tuple[str, str]],
    nt_seqs: List[Tuple[str, str]],
    do_log: bool,
):
    return pn2codon(filepath, table_index, aa_seqs, nt_seqs, do_log)


def init_argparse() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(

        usage="pal2nal.py [OPTIONS] [FILES]...",

        description="Convert Amino Acids to Codons with error checking using the default NCBI table, or your specified table."

    )

    parser.add_argument(
        "-t", "--table", default="genetic_table.json",
        help="""
        Path to the genetics table in JSON format.

        The format of JSON should be:
        [
            {
                "table_id:" 1,
                "M": ["ATG", "ATA"],
                "N": ...
            },
            {
                "table_id:" 2,
                "M": ["ATA"],
                "N": ...
            }
        ]
       
        Consult `genetic_table.json` which ships with the code to get a sense of how the JSON file should be formatted.
        """
    )

    parser.add_argument(
        "-i", "--table_id", default=1,
        help="""
        ID of the genetic codes translation table according to NCBI (if you go by the default table). 
        The IDs are as follows:

        Translation table 1: The standard code
        Translation table 2: The vertebrate mitochondrial code
        Translation table 3: The yeast mitochondrial code
        Translation table 4: The mold, protozoan, and coelenterate mitochondrial code and the mycoplasma/spiroplasma code
        Translation table 5: The invertebrate mitochondrial code
        Translation table 6: The ciliate, dasycladacean and hexamita nuclear code
        Translation table 7: The kinetoplast code; cf. table 4.
        Translation table 8: cf. table 1.
        Translation table 9: The echinoderm and flatworm mitochondrial code
        Translation table 10: The euplotid nuclear code
        Translation table 11: The bacterial, archaeal and plant plastid code
        Translation table 12: The alternative yeast nuclear code
        Translation table 13: The ascidian mitochondrial code
        Translation table 14: The alternative flatworm mitochondrial code
        Translation table 15: The Blepharisma nuclear code
        Translation table 16: The chlorophycean mitochondrial code
        Translation table 21: The trematode mitochondrial code
        Translation table 22: The Scenedesmus obliquus mitochondrial code
        Translation table 23: The Thraustochytrium mitochondrial code
        Translation table 24: The Pterobranchia mitochondrial code
        Translation table 25: The candidate division SR1 and gracilibacteria code
        Translation table 26: The Pachysolen tannophilus nuclear code
        Translation table 27: The karyorelict nuclear code
        Translation table 28: The Condylostoma nuclear code
        Translation table 29: The Mesodinium nuclear code
        Translation table 30: The Peritrich nuclear code
        Translation table 31: The Blastocrithidia nuclear code
        Translation table 32: The Balanophoraceae plastid code
        Translation table 33: The Cephalodiscidae mitochondrial code
        """
    )

    parser.add_argument(
        "-aaf", "--amino_acid_fasta", default="test.aa.fa",
        help="""
        A path to the FASTA file containing Amino Acid Sequences.
        """
    )

    parser.add_argument(
        "-ntf", "--nucleotide_fasta", default="test.nt.fa",
        help="""
        A path to the FASTA file containing source Nucleotide Sequences.
        """
    )

    parser.add_argument(
        "-tr", "--target_files", default="test_result.nt.fa",
        help="""
        Target file to save the final results as a FASTA (.fa) file
        """
    )

    parser.add_argument(
        "-l", "--log", default="0",
        help="""
        Enable logging (separate from errors).
        0 = Disable, 1 = Enable.
        Might be intrusive.
        """
    )

    return parser


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(
            f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


@timeit
def pal_to_nal() -> None:
    print("Starting the operation...")

    parser = init_argparse()
    args = parser.parse_args()

    table_path = args.table
    table_id = args.table_id
    aas_path = args.amino_acid_fasta
    nts_path = args.nucleotide_fasta
    target_path = args.target_files
    log = True if args.log == "1" else False

    aa_seqs, nt_seqs = read_and_convert_fasta_files(
        aas_path,
        nts_path
    )

    write_result_to_fasta_file(target_path, convert_to_codon(
        filepath=table_path,
        table_index=table_id,
        aa_seqs=aa_seqs,
        nt_seqs=nt_seqs,
        do_log=log
    ))

    print("Done! You can see the final time below.")


if __name__ == "__main__":
    pal_to_nal()
