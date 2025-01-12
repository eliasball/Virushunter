#!/usr/bin/env python3
import os
import gzip
import shutil
import logging
from Bio import SeqIO
from hashlib import md5
import ncbi_genome_download as ngd
import argparse 


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)

def extract_descendant_taxids(file_path):
    """
    Extracts the descendant_taxids column from a gimme_taxa output file.

    Parameters:
        file_path (str): Path to the gimme_taxa output file.

    Returns:
        list: A list of descendant taxonomic IDs.
    """
    tax_ids = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                if not line.startswith("parent_taxids"):  # Skip header
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        tax_ids.append(parts[1])  # Extract descendant_taxids
        logging.info(f"Extracted {len(tax_ids)} descendant tax IDs.")
    except Exception as e:
        logging.error(f"Failed to read taxids file: {e}")
    return tax_ids


def download_genomes(
    tax_ids=None,
    groups="viral",
    output="Viral_Genomes",
    file_formats=None,
    assembly_levels=None,
    sections=None,
    dry_run=False,
    flatten_output=False
):
    """
    Downloads genomes from NCBI using ncbi-genome-download and processes downloaded files.

    Parameters:
        tax_ids (list): List of taxonomic IDs to download genomes for.
        groups (str): The taxonomic group (e.g., "viral").
        output (str): Directory to store downloaded files.
        file_formats (list): List of file formats to download (e.g., ["fasta"]).
        assembly_levels (list): List of assembly levels (e.g., ["complete"]).
        sections (list): List of sections to download from (e.g., ["genbank", "refseq"]).
        dry_run (bool): If True, simulate downloads without actual downloading.
        flatten_output (bool): If True, ensures files are moved to the top-level folder.
    """
    file_formats = file_formats or ["fasta"]  # Default to fasta
    assembly_levels = assembly_levels or ["complete"]  # Default to complete
    sections = sections or ["genbank"]  # Default to genbank
    tax_ids = tax_ids or []  # Default to an empty list

    # Resolve output folder relative to script location
    output_path = os.path.abspath(output)  # Ensure absolute path

    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    try:
        for tax_id in tax_ids:
            logging.info(f"Processing tax ID: {tax_id}")
            for section in sections:
                logging.info(f"Starting download for section: {section}, tax ID: {tax_id}, group: {groups}")
                ngd.download(
                    section=section,
                    taxids=tax_id,
                    groups=groups,
                    output=output_path,  # Ensure files are downloaded to the correct folder
                    file_formats=file_formats,
                    assembly_levels=assembly_levels,
                    dry_run=dry_run
                )
                logging.info(f"Download completed for section: {section}, tax ID: {tax_id}")

            # Post-download processing
            process_downloaded_files(output_path, flatten_output, tax_id)
    except Exception as e:
        logging.error(f"An error occurred during genome download: {e}")

def filter_duplicate_genomes(output_path):
    """
    Removes duplicate genome sequences by comparing the sequences in all `.fna` files.

    Parameters:
        output_path (str): Path to the directory containing `.fna` files.
    """
    logging.info("Filtering duplicate genome sequences...")
    seen_sequences = set()  # Tracks unique genome sequences by hash
    files_to_remove = []  # Tracks files to be removed

    for file in os.listdir(output_path):
        if file.endswith(".fna"):
            file_path = os.path.join(output_path, file)
            try:
                # Read the sequences from the FASTA file
                sequences = [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]

                # Compute a combined hash for all sequences in the file
                sequence_hash = md5("".join(sequences).encode('utf-8')).hexdigest()

                # Check for duplicates
                if sequence_hash in seen_sequences:
                    logging.info(f"Duplicate file detected, removing: {file_path}")
                    files_to_remove.append(file_path)
                else:
                    seen_sequences.add(sequence_hash)
            except Exception as e:
                logging.error(f"Error reading file {file_path}: {e}")

    # Remove all files marked for removal
    for file_path in files_to_remove:
        logging.info(f"Removing duplicate file: {file_path}")
        os.remove(file_path)

    logging.info("Duplicate genome filtering completed.")
    
    


def process_downloaded_files(output_path, flatten_output, taxid):
    """
    Processes the downloaded files: unzips `.fna.gz` files, filters duplicates,
    and deletes empty subdirectories.

    Parameters:
        output_path (str): Path to the downloaded files directory.
        flatten_output (bool): If True, ensures the output is flattened to the top-level folder.
        taxid (str): Taxonomic ID to use for renaming files.
    """
    logging.info(f"Processing downloaded files for taxid {taxid}...")

    for root, _, files in os.walk(output_path):
        for file in files:
            file_path = os.path.join(root, file)

            # Unzip `.fna.gz` files
            if file.endswith(".fna.gz"):
                logging.info(f"Unzipping file: {file_path}")
                unzip_file(file_path, output_path, taxid)

            # Remove `MD5SUMS` files
            elif file == "MD5SUMS":
                logging.info(f"Removing file: {file_path}")
                os.remove(file_path)

    # Filter duplicates after unzipping
    filter_duplicate_genomes(output_path)

    # Remove empty subdirectories
    remove_empty_directories(output_path)
    logging.info(f"File processing for taxid {taxid} completed.")




def unzip_file(gzip_path, output_path, taxid):
    """
    Unzips a `.gz` file into the specified directory, processes it to keep the first header 
    and concatenate all sequences, and removes the original zipped file.
    Removes all fasta headers except for the first one and concatenates the sequences.

    Parameters:
        gzip_path (str): Path to the `.gz` file.
        output_path (str): Directory to store the processed file.
        taxid (str): Taxonomic ID to use for naming the file.
    """
    try:
        # Ensure unique naming based on the taxid
        file_name_base = f"{taxid}.fna"
        dest_path = os.path.join(output_path, file_name_base)
        file_counter = 1

        # Handle cases where a file with the same name already exists
        while os.path.exists(dest_path):
            file_name_base = f"{taxid}_{file_counter:02d}.fna"
            dest_path = os.path.join(output_path, file_name_base)
            file_counter += 1

        # Unzip and process the file
        with gzip.open(gzip_path, 'rt') as f_in:  # Read as text
            lines = f_in.readlines()

        first_header = None
        concatenated_sequence = ""

        for line in lines:
            if line.startswith(">"):
                if first_header is None:
                    first_header = line.strip()  # Store the first header
            else:
                concatenated_sequence += line.strip()  # Concatenate sequence without newlines

        # Write the processed content to the destination file
        with open(dest_path, 'w') as f_out:
            f_out.write(f"{first_header}\n{concatenated_sequence}\n")

        # Remove the original .gz file
        os.remove(gzip_path)
        logging.info(f"Unzipped and processed {gzip_path} to {dest_path}")

    except Exception as e:
        logging.error(f"Failed to unzip and process file {gzip_path}: {e}")



def remove_empty_directories(output_path):
    """
    Removes all empty subdirectories within the specified path.

    Parameters:
        output_path (str): Path to the directory to clean up.
    """
    logging.info("Removing empty directories...")
    for root, dirs, _ in os.walk(output_path, topdown=False):
        for directory in dirs:
            dir_path = os.path.join(root, directory)
            if not os.listdir(dir_path):  # Check if the directory is empty
                logging.info(f"Removing empty directory: {dir_path}")
                os.rmdir(dir_path)
    logging.info("Empty directories removed.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Genome extractor")
    parser.add_argument(
        "host_taxid",
        type=str,
        help="Host taxonomic ID to name the output folder"
    )
    args = parser.parse_args()

    # Path to the gimme_taxa output file
    taxid_file = "taxids.txt"

    # Extract descendant tax IDs
    descendant_tax_ids = extract_descendant_taxids(taxid_file)

    # Run genome download
    download_genomes(
        tax_ids=descendant_tax_ids,
        output=args.host_taxid,
        sections=["genbank", "refseq"],
        flatten_output=True
    )
