import re
#from collections import Counter
import pandas as pd
#import numpy as np
#from lxml import etree
#from ete3 import NCBITaxa
from nltk.tokenize import word_tokenize
#from sentence_splitter import SentenceSplitter, split_text_into_sentences
#import matplotlib.pyplot as plt
#from matplotlib_venn import venn2
#import matplotlib.patches as mpatches
#import seaborn as sns
#from natsort import natsorted


pd.set_option('display.width', 1000)


abstracts_df = pd.read_csv('~/2024-pubmed-mining/results/abstracts.tsv',
                           sep='\t', header=None, names=['Journal_Name', 'Year',
                                                         'Month', 'Title', 'Abstract',
                                                         'PMID', 'xml'])


# filter on abstracts mentionning a quinone
q_abstracts_df = abstracts_df[abstracts_df['Abstract'].str.contains('quinone|quinol', na=False, regex=True)]


pattern_q_quinone= r"(\b[Dd]?[UuRrMmCc]?[Qq][\s-]?[^aA-zZ]?\d{1,2}[\s-]?(?:(?:H\d{1,2})|[\(\[].*?[\)\]])?)" 
pattern_k_quinone= r"(\b[Dd]?[DdMm]?[Mm][Kk][\s-]?\d{1,2}[\s-]?(?:(?:H\d{1,2})|[\(\[].*?[\)\]])?)"
# no PQ in the dataset


def get_quinone_word(tokens_list):
    """
    Extract string (word) from a list of words that contains "*quinone*/ol*"
    
    Args:
        tokens_list (list): abstract tokenized 
    
    Returns:
        list : [quinone_list]
    """
    pattern_quinone = r"[a-zA-Z]*quinon[a-zA-Z]*"
    quinone_list = []
    print(tokens_list)
    for i in range(len(tokens_list)):
        res = re.match(pattern_quinone, tokens_list[i])
        if res is not None:
            if (i < len(tokens_list)-1 and tokens_list[i+1].isdigit()):
                quinone_list.append(tokens_list[i]+'-'+tokens_list[i+1])
            else:
                quinone_list.append(tokens_list[i])
    return quinone_list


def get_quinone_short(abstract, pattern):
    """
    Extract short mention of quinones (MK, UQ...)

    Args:
        abstract (str)
    
    Returns:
        list : [quinone_list]
    """
    quinone_list = re.findall(pattern, abstract)
    return quinone_list


def make_str(q_list):
    new_q_list = list(set(q_list))
    return ' '.join(new_q_list)#.lower()


def search_main_quinone(abstract):
    """
    Extract the main quinones.

    Args:
        abstract (str)
    
    Returns:
        list : [quinone_list]
    """
    pattern_main_quinone = r"main [a-zA-Z]*quinon[a-zA-Z]*|major [a-zA-Z]*quinon[a-zA-Z]*"
    res = re.findall(pattern_main_quinone, abstract)
    if len(res) != 0:
        match_before = re.search(r'[,.]\s*([^,.]*)\s+'+res[0], abstract)
        if match_before:
            before = match_before.group(1)
        else:
            before = ""
        string_to_extract = " "+before
        match_after = re.search(r''+res[0]+r'\s*([^.]*)\.', abstract)
        if match_after:
            after = match_after.group(1)
        else:
            after = ""
        string_to_extract = string_to_extract+" "+after+" "
        return string_to_extract
    return ''  


def extract_quinone(df):
    """
    Main function to extract quinones in abstract.

    Args :
        df (DataFrame): abstracts dataframe

    Returns:
        df (DataFrame): abstracts dataframe with additionnal columns containing quinones 
                        terms (short : i.e MK-6, UQ-10, MK8(H2); long: i.e ubiquinone, 
                        menaquinone-10)
    """
    # tokenise abstracts
    df['tokens'] = df['Abstract'].apply(word_tokenize)
    # extract word which contains "quinone"; "quinol"
    df['quinone_list'] = df['tokens'].apply(get_quinone_word)
    df['quin'] = df['quinone_list'].apply(make_str)
    df['k_short'] = df.apply(lambda x: get_quinone_short(x["Abstract"], pattern_k_quinone), axis=1)
    df['q_short'] = df.apply(lambda x: get_quinone_short(x["Abstract"], pattern_q_quinone), axis=1)
    df['main_sentence'] = df['Abstract'].apply(search_main_quinone)
    df['main_sentence_tokens'] = df['main_sentence'].apply(word_tokenize)
    df['main_quinone_list'] = df['main_sentence_tokens'].apply(get_quinone_word)
    df['main_quin'] = df['main_quinone_list'].apply(make_str)
    df['main_k_quin'] = df.apply(lambda x: get_quinone_short(x["main_sentence"], pattern_k_quinone), axis=1)
    df['main_q_quin'] = df.apply(lambda x: get_quinone_short(x["main_sentence"], pattern_q_quinone), axis=1)
    return df

quinone_abstract_df = extract_quinone(q_abstracts_df)


def filter_false_hit(q_list):
    """
    Function to filter out false hits.

    Args: 
        q_list (list)

    Returns:
        new_list (list): filtered list 
    """
    # remove if ends with "("
    new_list = [x for x in q_list if not x.endswith('(')]
    # remove if present > 4 times # if number > 30
    for elem in set(q_list):
        if new_list.count(elem) > 3:
            new_list = [x for x in new_list if x != elem]
        if elem.endswith(')'):
            if '(' in elem and elem.endswith(')'):
                continue
            if elem.endswith('))'):
                new_list = [x for x in new_list if x != elem]
                new_list.append(elem[:-2])
            else:
                new_list = [x for x in new_list if x != elem]
                new_list.append(elem[:-1])
        if elem.endswith('-') :
            new_list = [x for x in new_list if x != elem]
    return new_list


def filter_on_df(df):
    df['q_short_filt'] = df["q_short"].apply(filter_false_hit)
    df['k_short_filt'] = df["k_short"].apply(filter_false_hit)
    df['main_q_quin_filt'] = df["main_q_quin"].apply(filter_false_hit)
    df['main_k_quin_filt'] = df["main_k_quin"].apply(filter_false_hit)
    return df

quinone_abstract_filt_df = filter_on_df(quinone_abstract_df)

quinone_abstract_filt_df.to_csv('../results/intermediate_quinones_table.tsv', sep='\t')


def merge_quin_infos(row):
    """
    Merge quinone-related lists into a single filtered list.

    Args:
        row (pandas.Series): A row containing quinone-related columns.

    Returns:
        list: Merged list with filtered quinone entries, excluding certain terms (-quinone)
    """
    final_list = []
    liste_to_rm = ['quinone', 'quinones', 'benzoquinone', 'benzoquinones', 'naphtoquinones',
                    'naphtoquinone', 'lipoquinone', 'lipoquinones']
    row['quinone_list'] = set(row['quinone_list'])
    final_list.extend([x for x in row['quinone_list'] if not "anthraquin" in x and x not in liste_to_rm])
    final_list.extend(set(row['q_short_filt']))
    final_list.extend(set(row['k_short_filt']))
    return final_list


def merge_main_quin_infos(row):
    final_list = []
    liste_to_rm = ['quinone', 'quinones', 'benzoquinone', 'benzoquinones', 'naphtoquinones',
                    'naphtoquinone', 'lipoquinone', 'lipoquinones']
    row['main_quinone_list'] = set(row['main_quinone_list'])
    final_list.extend([x for x in row['main_quinone_list'] if not "anthraquin" in x and x not in liste_to_rm])
    final_list.extend(set(row['main_q_quin_filt']))
    final_list.extend(set(row['main_k_quin_filt']))
    return final_list


def build_dico_quinone(quin_col):
    """
    Build a dictionary of quinones.

    Args:
        quin_col (list of str): List of quinone strings.

    Returns:
        dict: Dictionary with detected quinones as keys and 1 as value.
    """
    liste = [x for x in quin_col if "/" not in x]
    liste.extend([x for x in quin_col if "/" in x])
    sat = r"H[\[\(]?\d{1,2}"
    quinones = {}
    for quin in liste:
        if quin.lower().startswith('ubi'):
            quinones['ubiquinone']=1
        if quin.lower().startswith('rhodo'):
            quinones['rhodoquinone']=1
        if quin.lower().startswith('caldar'):
            quinones['caldariellaquinone']=1
        if quin.lower().startswith('thermoplasma'):
            quinones['thermoplasmaquinone']=1
        if quin.lower().startswith('menathioq'):
            quinones['menathioquinone'] = 1            
        if quin.lower().startswith('methylmen') or quin.lower().startswith('monomet'):
            quinones['methylmenaquinone'] = 1
        if quin.lower().startswith('demethylmen'):
            quinones['demethylmenaquinone'] = 1
        if quin.lower().startswith('dimethylmen'):
            quinones['dimethylmenaquinone'] = 1
        if quin.lower().startswith('menaq') or quin.lower().startswith('mean') or quin.lower().startswith('mk'):
            quinones['menaquinone']=1
        if quin.lower().startswith('q') or  quin.lower().startswith('uq') or  quin.lower().startswith('cq') :
            quinones['ubiquinone']=1
        if quin.upper().startswith('MK') or quin.upper().startswith('MENAQ') or quin.upper().startswith('MEANAQ') :
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 20 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                if bool(re.search(sat, quinone)):
                    sat_match = re.findall(sat, quinone)
                    sat_match = ", ".join(sat_match)
                    sat_list = [int(s) for s in re.findall(r"\d+", sat_match)]
                    for i in range(len(sat_list)):
                        quinones['MK-'+str(length)+'(H'+str(sat_list[i])+')'] = 1
                else:
                    quinones['MK-'+str(length)] = 1
        if quin.upper().startswith('MQ'):
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            length = digits[0]
            quinones['MQ-'+str(length)] = 1
            quinones['methylene-ubiquinone'] = 1
        if quin.upper().startswith('DMMK') or quin.upper().startswith('DIMETHYL') :
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 20 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                if bool(re.search(sat, quinone)):
                    sat_match = re.search(sat, quinone).group(0)
                    sat_list = [int(s) for s in re.findall(r"\d+", sat_match)]
                    for i in range(len(sat_list)):
                        quinones['DMMK-'+str(length)+'(H'+str(sat_list[i])+')'] = 1
                else:
                    quinones['DMMK-'+str(length)] = 1
        if quin.upper().startswith('DMK') or quin.upper().startswith('DEMETHYL') :
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 20 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                if bool(re.search(sat, quinone)):
                    sat_match = re.search(sat, quinone).group(0)
                    sat_list = [int(s) for s in re.findall(r"\d+", sat_match)]
                    for i in range(len(sat_list)):
                        quinones['DMK-'+str(length)+'(H'+str(sat_list[i])+')'] = 1
                else:
                    quinones['DMK-'+str(length)] = 1
        if quin.upper().startswith('MMK') or quin.upper().startswith('MONOMET')  or quin.upper().startswith('METHYL'):
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 20 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                if bool(re.search(sat, quinone)):
                    sat_match = re.search(sat, quinone).group(0)
                    sat_list = [int(s) for s in re.findall(r"\d+", sat_match)]
                    for i in range(len(sat_list)):
                        quinones['MMK-'+str(length)+'(H'+str(sat_list[i])+')'] = 1
                else:
                    quinones['MMK-'+str(length)] = 1
        if quin.upper().startswith('UQ') or quin.upper().startswith('Q')  or quin.upper().startswith('UBIQ') or quin.upper().startswith('CQ') :
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 13 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                if bool(re.search(sat, quinone)):
                    sat_match = re.search(sat, quinone).group(0)
                    sat_list = [int(s) for s in re.findall(r"\d+", sat_match)]
                    for i in range(len(sat_list)):
                        quinones['UQ'+str(length)+'(H'+str(sat_list[i])+')'] = 1
                else:
                    quinones['UQ-'+str(length)] = 1
        if quin.upper().startswith('RHODO') or quin.upper().startswith('RQ'):
            quinone = quin.upper()
            digits = [int(s) for s in re.findall(r"\d+", quinone)]
            if len(digits) == 0:
                pass
            elif digits[0] >= 13 or digits[0] <= 3:
                pass
            else:
                length = digits[0]
                quinones['RQ-'+str(length)] = 1
    return quinones


def get_saturation(complete_quinone):
    """
    Extract saturation information from quinone strings using regex.
    
    Args:
        complete_quinone (list)
    
    return
        (str): Comma-separated string of unique saturation identifiers
    """
    #for i in range(len(complete_quinone)):
    complete_quinone_str = ', '.join(complete_quinone)
    pattern = r"(H)([(]?[(]?)?(\d{1,2})"
    convertion = {"₀":"0", "₁":"1", "₂":"2", "₃":"3", "₄":"4", "₅":"5",
                 "₆":"6", "₇": "7", "₈":"8", "₉": "9"}
    for key, value in convertion.items():
        complete_quinone_str = complete_quinone_str.replace(key, value)
    saturation = re.findall(pattern, complete_quinone_str)
    sat_list = []
    for sat in saturation:
        sat_list.append(sat[0]+sat[2])
    return ', '.join(sorted(set(sat_list)))


def get_chain_length(dico_quin):
    """
    Extract chain length information from quinone dictionary keys using regex.

    Args:
        dico_quin (dict): Dictionary with quinone names as keys.

    Returns:
        (str): Comma-separated string of unique chain lengths (sorted).
    """
    pattern = r"(.K|.Q)-(\d{1,2})"
    chain_length_list = []
    for key in dico_quin.keys():
        chain_lengths = re.findall(pattern, key)
        for cl in chain_lengths:
            chain_length_list.append(cl[1])
    return ', '.join(sorted(set(chain_length_list)))


def get_quinones_list(dico_quin):
    """
    Return a sorted list of quinones.
    Lowercase keys are listed first (alphabetically), followed by uppercase keys (alphabetically).
    
    Args:
        dico_quin (dict): Dictionary with quinone names as keys.
    
    Returns:
        sorted_quin_list (str): Comma-separated string of sorted quinone names.
    """
    upper_case = sorted([x for x in dico_quin.keys() if x[0].isupper()])
    lower_case = sorted([x for x in dico_quin.keys() if x[0].islower()])
    sorted_quin_list = ", ".join(lower_case + upper_case)
    return sorted_quin_list


def final_steps_extraction_quinone(df):
    df['complete_quinone'] = df.apply(merge_quin_infos, axis=1)
    df['main_complete_quinone'] = df.apply(merge_main_quin_infos, axis=1)
    df['dico_quin'] = df['complete_quinone'].apply(build_dico_quinone)
    df['dico_main_quin'] = df['main_complete_quinone'].apply(build_dico_quinone)
    df['saturation'] = df['complete_quinone'].apply(get_saturation)
    df['chain_length'] = df['dico_quin'].apply(get_chain_length)
    df['quinones'] = df['dico_quin'].apply(get_quinones_list)
    df['main_quinones'] = df['dico_main_quin'].apply(get_quinones_list)
    return df

quinone_abstract_filt_df = final_steps_extraction_quinone(quinone_abstract_filt_df)

quinone_abstract_filt_df.to_csv('../results/raw_table_quinones.tsv', sep ='\t')


quinone_abstract_filt_df.PMID.astype(str)
quinones_df = quinone_abstract_filt_df[['PMID', 'Journal_Name', 'Year', 'quinones', 'main_quinones', 'chain_length', 'saturation', 'main_sentence']]


#### merge species tab and quinone tab

## confidence score based on the number of species found in each PMID

species_df = pd.read_csv('../results/species_all_prokaryotes.tsv', sep='\t', index_col=0)

species_df["confidence"] = 1
species_df.loc[species_df["pmid_count"] > 1, "confidence"] = 2

quinones_species_df_prok_merged = quinones_df.merge(species_df, how='left')


cols_order = ['PMID', 'Journal_Name', 'Year', 'species', 'species_taxid',
                'genus', 'genus_taxid', 'family', 'family_taxid', 'order',
                'order_taxid', 'class', 'class_taxid', 'phylum', 'phylum_taxid',
                'kingdom', 'superkingdom', 'metabolism', 'metabolism_in_Title',
                'metabolism_in_Abstract', 'main_sentence', 'quinones',
                'main_quinones', 'saturation', 'chain_length', 'pmid_count', 'confidence']


quinones_species_df_prok_merged = quinones_species_df_prok_merged[cols_order]


quinones_species_df_prok_merged.to_csv("../results/species_quinone_results.tsv", 
                                        sep="\t", index=None)
