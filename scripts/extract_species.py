import re
import pandas as pd
from ete3 import NCBITaxa
from nltk.tokenize import word_tokenize


pd.set_option('display.width', 1000)

abstracts_df = pd.read_csv('~/2024-pubmed-mining/results/abstracts.tsv',
                           sep='\t', header=None, names=['Journal_Name', 'Year',
                                                         'Month', 'Title', 'Abstract',
                                                         'PMID', 'xml'])


## Extract species name

ncbi = NCBITaxa()

def add_simple_species(row):
    """
    Extract couple of words corresponding to the following pattern:
     - first word starting with an uppercase letter
     - second starting with a lowercase letter
    
    Args:
        row (dict): Dictionnary containing "Title" as key.

    Returns:
        list : [candidate_species, candidate_genus]
    """
    candidate_species = [] 
    candidate_genus = []
    tokens = word_tokenize(row["Title"])
    for i in range(len(tokens)-1):
        if tokens[i][0].isupper() and tokens[i+1][0].islower():
            candidate_species.append(' '.join([tokens[i], tokens[i+1]]))
            candidate_genus.append(tokens[i])
    return [candidate_species, candidate_genus]


def get_taxids(d):
    """
    Args:
        d (dict): Dictionnary having taxids lists as values

    Returns:
        list: deduplicated taxids list
    """
    tuples = d.values()
    taxids_list = []
    for tuple_el in tuples:
        for list_ids in tuple_el:
            taxids_list.append(list_ids)
    return list(set(taxids_list)) #deduplicate


def extract_tax(taxid_list):
    """
    Extract taxonomy from taxid

    Args:
        taxid_list (list): list of taxid (str)

    Return:
        list : list of dictionnary containing the taxonomy associated to each taxid
    """
    taxo_full = []
    for taxid in taxid_list:
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        lineage_translated = ncbi.get_taxid_translator(lineage)
        dico_taxo = {}
        for r_taxid, rank in ranks.items():
            if rank != 'no rank':
                dico_taxo[rank] = lineage_translated[r_taxid]
                dico_taxo[rank + '_taxid'] = r_taxid
        taxo_full.append(dico_taxo)
    return taxo_full


def iterate_over_df(df):
    """
    Pipeline to obtain species taxid and the complete associated taxonomy
    """
    df[['candidate_species', 'candidate_species_genus']] = df.apply(lambda x: add_simple_species(x), axis=1, result_type='expand')
    print(df['candidate_species'])
    df['tax_ids'] = df['candidate_species'].apply(ncbi.get_name_translator)
    df['tax_ids_genus'] = df['candidate_species_genus'].apply(ncbi.get_name_translator)
    df['tax_ids_simple'] = df['tax_ids'].apply(get_taxids)
    df['tax_ids_simple_genus'] = df['tax_ids_genus'].apply(get_taxids)
    df['full_taxo'] = df['tax_ids_simple'].apply(extract_tax)
    df['full_taxo_genus'] = df['tax_ids_simple_genus'].apply(extract_tax)
    return df


species_title_df = iterate_over_df(abstracts_df)


def count_tax_ids(df):
    df['n_tax_ids']= df['tax_ids_simple'].apply(len)
    return df

species_title_df = count_tax_ids(species_title_df)


## table with one line per taxa 

def return_new_row(PMID, taxo):
    """
    Args:
        PMID (str)
        taxo (dict) : dictionary containing the taxonomy

    Returns:
        row_new_df (dict)
    """
    superkingdom = taxo['superkingdom'] if 'superkingdom' in taxo.keys() else ''
    superkingdom_taxid = taxo['superkingdom_taxid'] if 'superkingdom_taxid' in taxo.keys() else ''
    kingdom = taxo['kingdom'] if 'kingdom' in taxo.keys() else ''
    phylum = taxo['phylum'] if 'phylum' in taxo.keys() else ''
    phylum_taxid = taxo['phylum_taxid'] if 'phylum_taxid' in taxo.keys() else ''
    class_tax = taxo['class'] if 'class' in taxo.keys() else ''
    class_tax_taxid = taxo['class_taxid'] if 'class_taxid' in taxo.keys() else ''
    order = taxo['order'] if 'order' in taxo.keys() else ''
    order_taxid = taxo['order_taxid'] if 'order_taxid' in taxo.keys() else ''
    family = taxo['family'] if 'family' in taxo.keys() else ''
    family_taxid = taxo['family_taxid'] if 'family_taxid' in taxo.keys() else ''
    genus = taxo['genus'] if 'genus' in taxo.keys() else ''
    genus_taxid = taxo['genus_taxid'] if 'genus_taxid' in taxo.keys() else ''
    species = taxo['species'] if 'species' in taxo.keys() else ''
    species_taxid = taxo['species_taxid'] if 'species_taxid' in taxo.keys() else ''
    row_new_df = {
            'PMID': PMID,
            'superkingdom': superkingdom,
            'superkingdom_taxid': superkingdom_taxid,
            'kingdom': kingdom,
            'phylum': phylum,
            'phylum_taxid': phylum_taxid,
            'class': class_tax,
            'class_taxid': class_tax_taxid,
            'order' : order,
            'order_taxid': order_taxid,
            'family': family,
            'family_taxid': family_taxid,
            'genus': genus,
            'genus_taxid': genus_taxid,
            'species': species,
            'species_taxid': species_taxid
        }
    return row_new_df


def get_species_table(df):
    """
    Args:
        df (DataFrame)

    Returns:
        new_df (DataFrame) : species table
    """
    new_df = pd.DataFrame(columns=['PMID', 'superkingdom', 'kingdom', 
                                   'phylum', 'class', 'order', 'family', 
                                   'genus', 'species'])
    for index, row in df.iterrows():
        PMID = row['PMID']
        taxos = row['full_taxo']
        if len(taxos) == 0:
            taxos = row['full_taxo_genus']
        if len(taxos) == 0:
            #PMID without taxonomy
            taxo = {}
            row_new_df = return_new_row(PMID, taxo)
            new_df = new_df.append(row_new_df, ignore_index=True)
        if len(taxos) == 1:
            taxo = taxos[0]
            row_new_df = return_new_row(PMID, taxo)
            new_df = new_df.append(row_new_df, ignore_index=True)
        else:
            for taxo in taxos:
                row_new_df = return_new_row(PMID, taxo)
                new_df = new_df.append(row_new_df, ignore_index=True)
    return new_df



species_df = get_species_table(species_title_df)


# Removal of bacillus-eukaryotes 
species_df = species_df[(species_df['genus'] != 'Bacillus') | (species_df['superkingdom'] != 'Eukaryota')]

# Removal of multiple mentions of the same species in the same title.
species_df = species_df.drop_duplicates().reset_index(drop=True)


species_df_prok = species_df[(species_df['superkingdom']=='Bacteria') | (species_df['superkingdom']== 'Archaea')].reset_index(drop=True)

species_df_prok['pmid_count'] = species_df_prok.groupby('PMID')['PMID'].transform('count')

species_df_prok.to_csv('../results/species_all_prokaryotes.tsv', sep='\t')
