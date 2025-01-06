import os
import requests
from bs4 import BeautifulSoup

def extract_virus_taxids_from_page(url):
    """
    Extract all virus taxids from multiple tables on the given Virus-Host DB page.

    Args:
        url (str): The URL of the Virus-Host DB page for a specific host taxid.

    Returns:
        list: A list of virus taxids extracted from the page.
    """
    # Fetch the webpage content
    response = requests.get(url)
    response.raise_for_status()  # Raise an error if the request fails

    # Parse the HTML content
    soup = BeautifulSoup(response.text, "html.parser")

    # Find all tables with class "info"
    tables = soup.find_all("table", class_="info")
    if not tables:
        print("No valid tables found on the page.")
        return []

    # Extract taxids from the "Scientific Name" row in each table
    taxids = []
    for table in tables:
        # Find the row containing the Scientific Name
        sci_name_row = table.find("tr")
        if not sci_name_row:
            continue

        # Find the link containing the taxid
        taxid_link = sci_name_row.find("a", href=True)
        if taxid_link:
            href = taxid_link["href"]
            # Extract the taxid (the part after the last "/")
            taxid = href.split("/")[-1]
            taxids.append(taxid)

    return taxids

# List of desired hosts to download

host_taxids = [10090, 9534, 9913, 9823, 9031, 1644094, 9455, 10036, 9544] # 9 major animal hosts
#host_taxids = [9534]
url = "https://www.genome.jp/virushostdb/"
formatted_taxids = []
for host in host_taxids:
    url = f"https://www.genome.jp/virushostdb/{host}"
    taxids = extract_virus_taxids_from_page(url)
    taxids = taxids[1:]
    formatted_taxids.extend(taxids)
    print(f"Extracting all sequences for host {host}")
    formatted_taxids_str = ','.join(map(str, formatted_taxids))
    os.system(f"gimme_taxa.py {formatted_taxids_str} > taxids.txt")
    os.system(f"python ../code/dataPreparation/genomeExtractor.py {host}")
