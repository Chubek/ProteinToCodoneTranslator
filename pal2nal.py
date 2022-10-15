import argparse
import os
import time
from functools import wraps
from typing import List, Tuple
from functools import reduce

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pro2codon import pn2codon
from requests import get


def download_and_save_file(url: str) -> str:
    name = url.split("/")[-1]

    if not any([p in url[:8] for p in ["http", "www"]]):
        return name

    if os.path.exists(name):
        print("File already exists, did not download.")
        return name

    resp = get(url)

    with open(name, "w") as fw:
        fw.write(resp.text)

    return name


def write_result_to_fasta_file(
    acc: List,
    tuple_input: Tuple[List[Tuple[str, str]], bool, str],
):
    file_name, append_file, seq_header_array = tuple_input

    list_str = sum(seq_header_array, ())

    if os.path.exists(file_name) and append_file:
        mode = "w+"
    else:
        mode = "w"

    with open(file_name, mode) as fw:
        fw.writelines(list_str)

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
    nts = []

    i = 0

    for aa, nt in zip(aa_seqs, nt_seqs):
        if aa[0] != nt[0]:
            continue

        aas.append((aa[0], aa[1]))
        nts.append((nt[0], nt[1]))
    
    return (aas, nts)


def convert_to_codon(
    filepath: str,
    table_index: int,
    aa_seqs: List[Tuple[str, str]],
    nt_seqs: List[Tuple[str, str]],
):
    return pn2codon(filepath, table_index, aa_seqs, nt_seqs)


def init_argparse() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(

        usage="pal2nal.py [OPTIONS] [FILES]...",

        description="Convert Amino Acids to Codons using the default NCBI table, or your specified table."

    )

    parser.add_argument(
        "-t", "--table", default="genetic_table.json",
        help="""
        Path or a direct download URL to the genetics table.
        Make sure that the URL starts with (http(s)://) or www
        The program checks if the file exists before attempting to download so make sure a file with the same name does not exist in the folder.        
        
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
        "-aaf", "--amino_acid_fasta", default="test.aa.fa", nargs="*",
        help="""
        One or multiple paths or direct download URLs to the FASTA file containing Amino Acid Sequences.
        Number of these arguments must be equal to `-ntf`

        The URL must contain (http(s)://) or www
        The program will check if the file already exists in the library before downloading it, so make sure a file of the same name does not already exist.
        """
    )

    parser.add_argument(
        "-ntf", "--nucleotide_fasta", default="test.nt.fa", nargs="*",
        help="""
        One or multiple paths or direct download URLs to the FASTA file containing source Nucleotide Sequences.
        Number of these arguments must be equal to `-aaf`

        The URL must contain (http(s)://) or www
        The program will check if the file already exists in the library before downloading it, so make sure a file of the same name does not already exist.
        """
    )

    parser.add_argument(
        "-tr", "--target_files", default="test_result.nt.fa", nargs="*",
        help="""
        Target file(s) to save the final results as a FASTA (.fa) file
        
        If one file specifed and multiple NTs and AAs are passed, it will combine them. If an equal number of files are specified, it will map each to each. Otherwise will error out.
        However if no argument is specified, it will save the result as `[source_nt_filename]_result_[num].fa`
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
    print("Starting the operation, if you have not read the help using `pal2nal.py -h` and this is your first time, please do!")

    parser = init_argparse()
    args = parser.parse_args()

    table = args.table
    table_id = args.table_id
    aas_paths_or_urls = sum([[args.amino_acid_fasta]], [])
    nts_paths_or_urls = sum([[args.nucleotide_fasta]], [])
    targets = sum([[args.target_files]], [])

    only_one = False

    if len(targets) != 0:
        if len(targets) == 1:
            if len(aas_paths_or_urls) == len(nts_paths_or_urls) == 1:
                only_one = True
        else:
            if len(aas_paths_or_urls) != len(nts_paths_or_urls):
                raise Exception(
                    "Length of NTs and AAs is not the same --- pass equal number of files, each mapping to the other respectively")
            elif len(aas_paths_or_urls) != len(targets):
                raise Exception(
                    "Length of targets is not equal to the number of inputs. Either pass none to name automatically or pass one to combine. Otherwise pass targets respectively.")
            else:
                print("I/O length check succeeded")

    table_file_name = download_and_save_file(table)

    i = [1]

    def work(tup):
        aa_f, nt_f = tup

        aa_f = download_and_save_file(aa_f)
        nt_f = download_and_save_file(nt_f)

        aa_seqs, nt_seqs = read_and_convert_fasta_files(aa_f, nt_f)

        i[0] += 1

        append_file = False

        if len(targets) == 0:
            fname = nt_f.split(".")[-2].replace(".nt", "")
            fname = f"{fname}_result_{i[0]}.nt.fa"
        elif only_one:
            fname = targets[0]
        elif len(targets) == 1 and not only_one:
            fname = nt_f.split(".")[-2].replace(".nt", "")
            fname = f"{fname}_result.nt.fa"
            append_file = True
        else:
            fname = targets[i[0] - 1]

        return (fname, append_file, convert_to_codon(
            filepath=table_file_name, table_index=table_id, aa_seqs=aa_seqs, nt_seqs=nt_seqs))

    zip_list = list(zip(
        aas_paths_or_urls,
        nts_paths_or_urls
    ))

    reduce(write_result_to_fasta_file, map(work, zip_list), [])

    print("Done! You can see the final time below.")


if __name__ == "__main__":
    pal_to_nal()
