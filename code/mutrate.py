import os, sys
import numpy as np
import pandas as pd
import gzip as gz
import glob as glob
import multiprocessing as mp

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_gff(fileIN):
    """
    Parses info from GenCode gff file
    Input: path to gff file
    Output: writte to file
    """
    output = open('../data/processed/gff.txt', 'w')
    output.writelines(["name", '\t', "class", '\t', "type", '\t', "transcript", '\t', "chr", '\t', "start", '\t', "end", '\t', "strand", '\n'])
    i = 0
    with gz.open(fileIN, 'rb') as f: 
        for line in f:
            i += 1
            if i > 7:		# ignore header
                current_line = line.decode('utf-8').rstrip().split('\t')
                #print(current_line)
                if len(current_line) >= 9:
                    current_chr = current_line[0].replace('chr', '')
                    current_type = current_line[2]
                    current_start = current_line[3]
                    current_end = current_line[4]
                    current_strand = current_line[6]

                    current_anno = current_line[8].split(';')
                    anno_dict = {}
                    for j in range(len(current_anno)):
                        current_entry = current_anno[j].split("=")
                        anno_dict[current_entry[0]] = current_entry[1]
                    current_class = anno_dict.get('gene_type', "NA")
                    current_name = anno_dict.get('gene_name', "NA")
                    current_transcript = anno_dict.get('transcript_id', "NA")
                    if current_transcript != "NA":
                        current_transcript = current_transcript.split('.')[0]


                    if current_class == "protein_coding": # and current_type == "gene":
                        #gffDF.loc[len(gffDF)] = (current_name, current_class, current_type, current_chr, current_start, current_end, current_strand)
                        output.writelines([current_name, '\t', current_class, '\t', current_type, '\t', current_transcript, '\t', current_chr, '\t', current_start, '\t', current_end, '\t', current_strand, '\n'])
            if i%1000000 == 0:
                print(i)
   
    output.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_uniprot_ids():
    """
    Match Uniprot ID to Ensembl canonical transcript IDs from Uniprot annotation
    Input: hard-coded table with Ensembl and Uniprot IDs
    Output: written to file
    """

    uniprot = pd.read_csv("../data/uniprotkb_taxonomy_id_9606_AND_reviewed_2024_08_11.tsv", header=0, index_col=False, sep='\t')

    resultDF = pd.DataFrame(columns=['uniprot', 'canonical'])
    for i in range(len(uniprot)):
        current_entry = uniprot.iloc[i]
        current_uni = current_entry['Entry']
        current_ensembl = str(current_entry['Ensembl'])
        current_altern = str(current_entry['Alternative products (isoforms)'])

        if current_ensembl[0:4] == "ENST": 
            if "Sequence=Displayed;" in current_altern:
                current_canonical = current_altern.split("Sequence=Displayed;")[0].split("IsoId=")[-1].replace(';','').split()[0]
                current_ensembl = current_ensembl.split(';')
                for j in range(len(current_ensembl)):
                    j_ensembl = current_ensembl[j]
                    if '[' in j_ensembl:
                        j_ensembl = j_ensembl.split('[')
                        current_key = j_ensembl[1].replace(']', '').split()[0]
                        #print(np.array(current_key.split()[0]), np.array(current_canonical.split())[0])
                        if np.array(current_key.split()[0]) == np.array(current_canonical.split()[0]):        # wrangle weird string formatting
                            canonical_transcript = j_ensembl[0]
                            break
            else:
                list_ensembl = current_ensembl.strip(';').split(';')
                canonical_transcript = list_ensembl[0]

            resultDF.loc[len(resultDF)] = (current_uni, canonical_transcript.split('.')[0])

    resultDF.to_csv("../data/processed/uniprot_canonical.txt", header=True, index=False, sep='\t')

    return resultDF




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cdsseq_dict(fileIN):
    """
    Parse dictionary of transcript coding sequences
    """

    seq = SeqIO.parse(fileIN, "fasta")

    seq_dict = {}
    for s in seq:
        current_name = s.name.split("|")
        current_enst = ""
        current_cds = ""
        for i in range(len(current_name)):
            if current_name[i][0:4]=="ENST":
                current_enst = current_name[i].split('.')[0]
            if current_name[i][0:4]=="CDS:":
                current_cds = np.array(current_name[i].split(':')[1].split('-') , dtype=int)   
        current_seq = s.seq
 
        if len(current_enst) > 0 and len(current_cds) > 0: 
            seq_cds = current_seq[(current_cds[0]-1):(current_cds[1]+2)] 	# force as multiple of 3
            seq_cds = seq_cds[:int(np.floor(len(seq_cds)/3.)*3.)]		 	#
            seq_dict[current_enst] = seq_cds

    return seq_dict




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def prepare_nt_tables(enstIN, tseqs, gff):
    """
    Prepares tables of chromosom and transcript positions for adding mutations, 
    consequence predictions, and counts
    Transcript boundaries from gff file
    Input: list of canonical transcript Ensembl IDs, sequence dictionary, gff file
    Output: one table for each gene, written to files in ../data/nt/
    """

    outputDIR = '../data/nt/'
    if not os.path.exists(outputDIR):
        os.mkdir(outputDIR)

    for i in sorted(list(set(enstIN['canonical']))):
    	if i in list(gff['transcript']) and i in list(tseqs.keys()):
            current_gff = gff[gff['transcript']==i]
            gff_cds = current_gff[current_gff['type']=="exon"]
 
            current_name = list(current_gff['name'])[0]
            current_strand = list(set(current_gff['strand']))[0]
            current_chr = list(set(current_gff['chr']))[0]

            if "start_codon" in list(current_gff['type']):
                start_codon = current_gff[current_gff['type']=="start_codon"]
                start_codon_start = int(list(start_codon['start'])[0])
                start_codon_end = int(list(start_codon['end'])[0])
            else:
                start_codon_start = int(np.min(np.array(current_gff['start']))) 	# if no entries, span whole transcript
                start_codon_end = int(np.max(np.array(current_gff['end'])))			#

            if "stop_codon" in list(current_gff['type']):
                stop_codon = current_gff[current_gff['type']=="stop_codon"]
                stop_codon_start = int(list(stop_codon['start'])[0])
                stop_codon_end = int(list(stop_codon['end'])[0])
            else:
                stop_codon_start = int(np.min(np.array(current_gff['start'])))		# same for negative strand
                stop_codon_end = int(np.max(np.array(current_gff['end'])))			#	

            transcript_position = 1
            current_seq = tseqs.get(i, "NA")

            output = gz.open("../data/nt/"+current_name+".txt.gz", "wt")
            output.writelines("Chr" + '\t' + "Chr_pos" + '\t' + "Pos" + '\t' + "WT" + '\t' + "MUT" + '\n')

            if current_strand == "+":
                for j in np.arange(len(gff_cds)):
                    current_exon = gff_cds.iloc[j]
                    current_start = current_exon['start']
                    current_end = current_exon['end']
                    for chr_position in range(current_start, current_end+1):
                        if (chr_position >= start_codon_start) and (chr_position <= stop_codon_end):
                            for nt in ['A', 'C', 'G', 'T']:
                                output.writelines(str(current_chr) + '\t' + str(chr_position) + '\t' + str(transcript_position) + '\t' + str(current_seq[transcript_position-1]) + '\t' + nt + '\n'  )
                            transcript_position += 1

            elif current_strand == "-":
                for j in np.arange(len(gff_cds)): #[::-1]:
                    current_exon = gff_cds.iloc[j]
                    current_start = current_exon['start']
                    current_end = current_exon['end']
                    for chr_position in range(current_end, current_start-1, -1):
                        if (chr_position <= start_codon_end) and (chr_position >= stop_codon_start):
                            for nt in ['A', 'C', 'G', 'T']:
                                output.writelines(str(current_chr) + '\t' + str(chr_position) + '\t' + str(transcript_position) + '\t' + str(current_seq[transcript_position-1]) + '\t' + nt + '\n'  )
                            transcript_position += 1

            output.close()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def variantpred_by_transcript(fileIN):
    """

    This function loops through the dbNSFP database and copies the scores for 
    each missense variant by the selected variant predictors to a temporary file
    for each gene. The use of temporary files is not efficient, but makes this 
    possible on memory-limited machines.   


    Split select variant predictions from dbNSFP4.7 into per-transcript files
  
    Columns of dbNSFP_variant: 
    1	chr: chromosome number
    2	pos(1-based): physical position on the chromosome 
    3	ref: reference nucleotide allele (as on the + strand)
    4	alt: alternative nucleotide allele (as on the + strand)
    5	aaref: reference amino acid
    6	aaalt: alternative amino acid
    12	aapos: amino acid position as to the protein.
    13	genename: gene name; if the nsSNV can be assigned to multiple genes, gene names are separated by ";"
    14	Ensembl_geneid: Ensembl gene id
    15  Ensembl_transcriptid: Ensembl transcript ids (Multiple entries separated by ";")
    46	Polyphen2_HDIV_pred: Polyphen2 prediction based on HumDiv
    49	Polyphen2_HVAR_pred: Polyphen2 prediction based on HumVar
    122	ESM1b_score: ESM1b scores 
    124	ESM1b_pred: prediction is based on the threshold of -7.5 
    125	EVE_score: The larger the score, the more likely the variant is pathogenic. Detals see https://doi.org/10.1038/s41586-021-04043-8.
    126	EVE_rankscore: The rankscore is the ratio of the rank of the score over the total number of EVE scores in dbNSFP.
    129 EVE_Class25_pred: The EVE classification (B)enign, (U)ncertain, or (P)athogenic when setting 25% as uncertain in a Gaussian mixture model.
    138	AlphaMissense_score: unsupervised model for predicting the pathogenicity of human missense variants
    140	AlphaMissense_pred: likely (B)enign, (A)mbiguous, or likely (P)athogenic with 90% expected precision estimated from ClinVar for likely benign and likely pathogenic classes.
    204	1000Gp3_AC: Alternative allele counts in the whole 1000 genomes phase 3 (1000Gp3) data.
    205	1000Gp3_AF: Alternative allele frequency in the whole 1000Gp3 data.
    352	gnomAD_genomes_AC: Alternative allele count in the whole gnomAD genome samples v4.0.0
    353	gnomAD_genomes_AN: Total allele count in the whole gnomAD genome samples v4.0.0
    429	ALFA_Total_AC: Alternative allele count of the total samples in the Allele Frequency Aggregator
    430	ALFA_Total_AN: Total allele count of the total samples in the Allele Frequency Aggregator
    """


    if not os.path.exists('TMP/'):
        os.mkdir('TMP/')

    
    # columns detailed aboved, adjusted for phython indexing  
    sel_columns = np.array([1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 46, 49, 122, 124, 125, 126, 129, 138, 140, 204, 205, 352, 353, 429, 430]) -1


    i = 0
    with gz.open(fileIN, 'rb') as f: 
        for line in f:
            current_line = line.decode('utf-8').rstrip().split('\t')
            if i == 0:
                #current_header = "\t".join(current_line[:203]) + "\n"
                current_header = "\t".join(list(np.array(current_line)[sel_columns])) + "\n"
            elif i > 0:
                current_ensg = current_line[13].split(';')[0]
                #current_write = "\t".join(current_line[:203]) + "\n"
                current_write = "\t".join(list(np.array(current_line)[sel_columns])) + "\n"
            
                if os.path.exists("TMP/"+current_ensg+".txt.gz"):
                    output = gz.open("TMP/"+current_ensg+".txt.gz", "at")
                else:
                    output = gz.open("TMP/"+current_ensg+".txt.gz", "wt")
                    output.writelines(current_header)

                output.writelines(current_write)
                output.close()

            i += 1

            #if i>100:
            #    break



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_snpdb(fileIN, clean=False):
    """
    THIS IS NOT USED IN THIS WORK AND CAN BE OMITTED

    Split SNPdb file into per-gene temporary files (to work with present hardwar constraints)
    Input file format from NIH SNPdb

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  SAMN10492695    SAMN10492696    SAMN10492697    SAMN10492698    SAMN10492699    SAMN10492700    SAMN10492701    SAMN10492702    SAMN11605645    SAMN10492703    SAMN10492704    SAMN10492705

    #SAMN10492695   "European"
    #SAMN10492696   "African Others"
    #SAMN10492697   "East Asian"
    #SAMN10492698   "African American"
    #SAMN10492699   "Latin American 1"
    #SAMN10492700   "Latin American 2"
    #SAMN10492701   "Other Asian"
    #SAMN10492702   "South Asian"
    #SAMN11605645   "African"
    #SAMN10492703   "Asian"
    #SAMN10492704   "Total"
    #SAMN10492705   "Other"

    #['##fileformat=VCFv4.1']
    #['##build_id=20230706150541']
    #['##Population=https://www.ncbi.nlm.nih.gov/biosample/?term=GRAF-pop']
    #['##FORMAT=<ID=AN,Number=1,Type=Integer,Description="Total allele count for the population, including REF">']
    #['##FORMAT=<ID=AC,Number=A,Type=Integer,Description="Allele count for each ALT allele for the population">']
    #['##FORMAT=<ID=HWEP,Number=1,Type=Integer,Description="int(-log(HWE score test p-value)); -1 indicates that the HWE score test p-value could not be computed">']
    #['##FORMAT=<ID=GR,Number=1,Type=Integer,Description="Genotype homozygous reference allele (AA) count; in rare cases may not be the GRCh reference allele">']
    #['##FORMAT=<ID=GV,Number=1,Type=Integer,Description="Genotype heterozygous ref/alt (A/B) count; reported for the most common two alleles that may or may not include the reference allele">']
    #['##FORMAT=<ID=GA,Number=1,Type=Integer,Description="Genotype homozygous alternate allele (B/B) count; could be any of the non-biallelic variant alleles.">']
    #['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMN10492695', 'SAMN10492696', 'SAMN10492697', 'SAMN10492698', 'SAMN10492699', 'SAMN10492700', 'SAMN10492701', 'SAMN10492702', 'SAMN11605645', 'SAMN10492703', 'SAMN10492704', 'SAMN10492705']
    """

    chrnames = pd.read_csv("../data/human/GCA_000001405.15_GRCh38_assembly_report+ucsc_names.txt", header=59, sep='\t', index_col=False)
    chr_dict = {}
    for i in range(24):
        current_refseq = chrnames.iloc[i]['RefSeq-Accn']
        current_ucsc = chrnames.iloc[i]['UCSC-style-name']
        chr_dict[current_refseq] = current_ucsc


    ## clean previously written files
    if clean:
        genefiles = glob.glob("../data/genes/*.txt")
        for i in genefiles:
            os.remove(i)


    tmplist=[]
    i = 0
    with gz.open(fileIN, 'rb') as f: 
        for line in f:
            current_line = line.decode('utf-8').rstrip() #.split('\t')
            if current_line[0:1] != "#":
                current_line = current_line.split('\t')
                current_refseq = current_line[0]
                current_chr = chr_dict.get(current_refseq, None)

                if current_chr != None:
                    current_chr = current_chr.replace('chr', '')
                    current_pos = int(current_line[1])
              
                    current_wt = current_line[3]
                    current_mut = current_line[4]

                    current_counts_european = current_line[9]
                    current_counts_african = current_line[10]
                    current_counts_eastasian = current_line[11]
                    current_counts_africanamerican = current_line[12]
                    current_counts_latam1 = current_line[13]
                    current_counts_latam2 = current_line[14]
                    current_counts_otherasian = current_line[15]
                    current_counts_southasian = current_line[16]
                    current_counts_african = current_line[17]
                    current_counts_asian = current_line[18]
                    current_counts_total = current_line[19]
                    current_counts_other = current_line[20]

                    current_write = str(current_chr) + '\t' + str(current_pos) + '\t' + str(current_wt) + '\t' + str(current_mut) + '\t' + str(current_counts_european) + '\t' + str(current_counts_african) + '\t' + str(current_counts_eastasian) + '\t' + str(current_counts_africanamerican) + '\t' + str(current_counts_latam1) + '\t' + str(current_counts_latam2) + '\t' + str(current_counts_otherasian) + '\t' + str(current_counts_southasian) + '\t' + str(current_counts_african) + '\t' + str(current_counts_asian) + '\t' + str(current_counts_total) + '\t' + str(current_counts_other) + '\n'
                    current_header = "Chr" + '\t' + "pos" + '\t' + "WT" + '\t' + "MUT" + '\t' + "European" + '\t' + "AfricanOther" + '\t' + "EastAsian" + '\t' + "AfricanAmerican" + '\t' + "Latam1" + '\t' + "Latam2" + '\t' + "OtherAsian" + '\t' + "SouthAsian" + '\t' + "African" + '\t' + "Asian" + '\t' + "Total" + '\t' + "Other" + '\n'  

                    current_gff = gff_dict.get(current_chr, None)
                    current_orfs = orf_dict.get(current_chr, None)
                    current_sel = (current_gff[:,0] <= current_pos) * (current_gff[:,1] >= current_pos)

                    if np.sum(current_sel) > 0:
                        for j in np.where(current_sel == True)[0]:
                            gff_name = current_orfs[j]

                            if os.path.exists("../data/genes/"+gff_name+".txt.gz"):
                                output = gz.open("../data/genes/"+gff_name+".txt.gz", "at")
                            else:
                                output = gz.open("../data/genes/"+gff_name+".txt.gz", "wt")
                                output.writelines(current_header)

                            output.writelines([current_write])
                            output.close()
            i += 1

            if i % 10000000 == 0:
                print(i, "processed")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def add_mut_data_todf(name, ensg, enst, seq):
    """
    Adds predicted variant effet data and allele counts to per-gene nt mutation tables
    Included variant predictors: Polyphen2, EVE, AlphaMissense, ESM1b
    Included population data: SNPdb, ALFA, Gnomad
    For reasons of efficiency/available resources, precomputed per-gene data is merged
    # ex.: add_mut_data_todf("VIT", "ENSG00000205221", 'ENST00000389975')
    Input: name (gene name), ensg (Ensembl gene ID) and enst (Ensemble transcript ID)
    Output is written to ../data/ntfull/
    """
    
    def index_score(scoreIN, idx):
        """
        retrieves transcript specific score at position idx:
        score format: A;B;C;D
        idx: index of the (canonical) transcript of interest
        """
        score = str(scoreIN)
        scoreOUT = "NA"

        if ';' not in score:
            scoreOUT = score
        else:
            score = score.split(';')
            if len(score) > idx:
                scoreOUT = score[idx]
            else:
                scoreOUT = "NA"

        return scoreOUT


    outputDIR = '../data/ntfull/'
    if not os.path.exists(outputDIR):
        os.mkdir(outputDIR)

 
    file_nt = "../data/nt/"+name+".txt.gz"
    file_snpdb = "../data/genes/"+name+".txt.gz"
    file_dbnsfp = "TMP/"+ensg+".txt.gz"

    table = pd.read_csv(file_nt, compression="gzip", header=0, index_col=False, sep='\t', low_memory=False)
    SNPdb = pd.read_csv(file_snpdb, compression="gzip", header=0, index_col=False, sep='\t', low_memory=False)
    dbNSFP = pd.read_csv(file_dbnsfp, compression="gzip", header=0, index_col=False, sep='\t', low_memory=False)
    
    L = len(table)

    current_seq = np.array(seq)
    current_codon_seq = list( np.reshape(current_seq, (int(np.floor(len(current_seq)/3.)), 3 ) ) )

    # initializen additional columns
    table['codon_wt'] = '.'
    table['codon_mut'] = '.'

    table['polyphen2_hdiv'] = '.'
    table['polyphen2_hvar'] = '.'
    table['EVE'] = '.'
    table['AlphaMissense'] = '.'
    table['ESM1b'] = '.'

    table['SNPdb_AN'] = '.'
    table['SNPdb_AC'] = '.'
    table['ALFA_AN'] = '.'
    table['ALFA_AC'] = '.'
    table['gnomad_AN'] = '.'
    table['gnomad_AC'] = '.'


    for i in range(L):           # len(table)
        current_chr = table.iloc[i]['Chr']
        current_chrpos = table.iloc[i]['Chr_pos']
        current_wt = table.iloc[i]['WT']
        current_mut = table.iloc[i]['MUT']

        current_pos = table.iloc[i]['Pos'] - 1			# python indexing
        current_wt = table.iloc[i]['WT']
        current_mut = table.iloc[i]['MUT']
        current_nt = current_seq[current_pos]

        codon_position = int( np.floor( (current_pos) / 3. ) )
        current_codon =  current_codon_seq[codon_position]
        current_phase = int( current_pos % 3. ) 
        current_codon_mut = np.copy(current_codon)
        current_codon_mut[current_phase] = current_mut

        # add codons
        if current_nt == current_wt:
            table.loc[i, 'codon_wt'] = ''.join( list( current_codon ) ) 
            table.loc[i, 'codon_mut'] = ''.join( list( current_codon_mut ) ) 
        #print(i, current_pos, current_nt, current_wt, current_mut, seq_codon_wt, seq_codon_mut)

        # crossref SNPdb data
        SNPdb_chr = np.array(SNPdb['Chr']) == current_chr
        SNPdb_chrpos = np.array(SNPdb['pos']) == current_chrpos
        SNPdb_wt = np.array(SNPdb['WT']) == current_wt
        SNPdb_mut = np.array(SNPdb['MUT']) == current_mut
        SNPdb_sel = SNPdb_chr * SNPdb_chrpos * SNPdb_wt * SNPdb_mut
        if np.sum(SNPdb_sel) > 0:
            SNPdb_total = str(list(SNPdb.iloc[SNPdb_sel]['Total'])[0]).split(':')
            table.loc[i, 'SNPdb_AN'] = int(SNPdb_total[0])
            table.loc[i, 'SNPdb_AC'] = int(SNPdb_total[1])

        # crossref dbNSFP data
        dbNSFP_chr = np.array(dbNSFP['#chr']) == current_chr
        dbNSFP_chrpos = np.array(dbNSFP['pos(1-based)']) == current_chrpos
        dbNSFP_wt = np.array(dbNSFP['ref']) == current_wt
        dbNSFP_mut = np.array(dbNSFP['alt']) == current_mut
        dbNSFP_sel = dbNSFP_chr * dbNSFP_chrpos * dbNSFP_wt * dbNSFP_mut
        current_dbNSFP = dbNSFP.iloc[dbNSFP_sel]
        if len(current_dbNSFP) > 0:
            for j in range(len(current_dbNSFP)):
                dbNSFP_j = current_dbNSFP.iloc[j]
                dbNSFP_enst = dbNSFP_j['Ensembl_transcriptid'].strip().split(';')
        
                if enst in dbNSFP_enst:
                    enst_idx = dbNSFP_enst.index(enst)
        
                    table.loc[i, 'polyphen2_hdiv'] = index_score(list(current_dbNSFP['Polyphen2_HDIV_pred'])[0], enst_idx)
                    table.loc[i, 'polyphen2_hvar'] = index_score(list(current_dbNSFP['Polyphen2_HVAR_pred'])[0], enst_idx)
                    table.loc[i, 'EVE'] = index_score(list(current_dbNSFP['EVE_Class25_pred'])[0], enst_idx)
                    table.loc[i, 'AlphaMissense'] = index_score(list(current_dbNSFP['AlphaMissense_pred'])[0], enst_idx)
                    table.loc[i, 'ESM1b'] = index_score(list(current_dbNSFP['ESM1b_pred'])[0], enst_idx)

                    table.loc[i, 'ALFA_AN'] = list(current_dbNSFP['ALFA_Total_AN'])[0]
                    table.loc[i, 'ALFA_AC'] = list(current_dbNSFP['ALFA_Total_AC'])[0]
                    table.loc[i, 'gnomad_AN'] = list(current_dbNSFP['gnomAD_genomes_AN'])[0]
                    table.loc[i, 'gnomad_AC'] = list(current_dbNSFP['gnomAD_genomes_AC'])[0]

                    break

    table.to_csv("../data/ntfull/"+name+".txt.gz", header=True, index=False, sep='\t', compression='gzip')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mut_data_parallel(packaged_input):
    """
    import/export wrapper for parallelized add_mut function
    """

    input_ensg = packaged_input[0]
    input_enst = packaged_input[1]
    input_name = packaged_input[2]
    input_seq  = packaged_input[3]

    file_nt = "../data/nt/"+input_name+".txt.gz"
    file_snpdb = "../data/genes/"+input_name+".txt.gz"
    file_dbnsfp = "TMP/"+input_ensg+".txt.gz"

    if os.path.exists(file_nt) and os.path.exists(file_snpdb) and os.path.exists(file_dbnsfp):
        add_mut_data_todf(input_name, input_ensg, input_enst, input_seq)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def normalize_mutmat():
    """
    Scale nucleic acid mutation matrix from Lynch (PNAS 2010) to yield
    an overall average mutation rate of 0.46 * 10e-9 mutations per site per generation 
    """

    ccnt = pd.read_csv("../data/processed/composition_nt.txt", header=0, index_col=False, sep='\t')

    #mutmat = np.array([[0, 1.52, 4.29, 1.29],
    #                   [2.58, 0, 2.95, 9.61],
    #                   [9.61, 2.95, 0, 2.58],
    #                   [1.29, 4.29, 1.52, 0]]) #* 10e-9

    mutmat = np.array([[0, 0.078, 0.221, 0.067],
                       [0.110, 0, 0.125, 0.408],
                       [0.408, 0.125 , 0, 0.110],
                       [0.067, 0.221, 0.078, 0]]) # *10e-9 

    n = np.mean(np.sum(mutmat, 1) * ccnt['count'] / np.mean(ccnt['count'])) / 0.46 
    mm2 = mutmat / n

    print(np.mean(np.sum(mm2, 1)) )

    return mm2



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_riskscore(df, predictor, unweighted=False):
    """
    Computes genotypic mutation risk score
    VEPs used: 'polyphen2_hdiv', 'polyphen2_hvar', 'EVE', 'AlphaMissense', 'ESM1b'
    Input:  - per-gene dataframe with all possible mutations and prediction scores
            - choice of prediction algorithm
            - unweighted: if True
    Output: - average score, i.e. normalized by sequence legnth
            - full score
            - neutral score, i.e. all mutations, not just damaging ones
            - coverage, i.e. percentage of variants with prediction scores vs. all possible variants
    """

    # mutation rates from Lynch paper
    # normalized to 0.46 * 10e-9 mutations per site per generation
    mutmat = np.array([[0., 0.07054, 0.19986, 0.06059],
                       [0.09948, 0., 0.11304, 0.36897],
                       [0.36897, 0.11304, 0., 0.09948],
                       [0.06059, 0.19986, 0.07054, 0.]
                        ])

    if unweighted:                                       # not ts/tv weighted, all equal
        mutmat = np.array([[0, 0.15333, 0.15333, 0.15333],
                           [0.15333, 0, 0.15333, 0.15333],
                           [0.15333, 0.15333, 0, 0.15333],
                           [0.15333, 0.15333, 0.15333, 0]])

    mutmat = pd.DataFrame(data=mutmat, columns=['A', 'C', 'G', 'T'], index=['A', 'C', 'G', 'T'])


    if predictor in ['polyphen2_hdiv', 'polyphen2_hvar', 'EVE', 'AlphaMissense', 'ESM1b']:
        if predictor == 'AlphaMissense':
            patho = 'P'							# (B)enign, (A)mbiguous, or likely (P)athogenic
        elif predictor == 'EVE':
            patho = 'P' 						# (B)enign, (U)ncertain, or (P)athogenic
        elif predictor == 'ESM1b':
            patho = 'D' 						# deleterious (D) vs. tolerated (T) variants
        elif predictor == 'polyphen2_hdiv':		
            patho = 'D' 						# "D" ("probably damaging"),"P" ("possibly damaging")  and "B" ("benign")
        elif predictor == 'polyphen2_hvar':			
            patho = 'D' 						# "D" ("probably damaging"), "P" ("possibly damaging") and "B" ("benign")

        A = 0.
        P = 0.

        tot = 0.
        cov = 0.

        N = 0.

        list_positions = sorted(list(set(df['Pos'])))
        for i in list_positions:
            current_df = df[df['Pos']==i]
            current_df = current_df[current_df['WT']!=current_df['MUT']]

            current_wt = list(set(current_df['WT']))[0]
            current_mut = np.array(list(current_df['MUT']))
            current_pred = np.array(list(current_df[predictor]))

            N += np.sum( np.array(mutmat.loc[current_wt]) )

            tot += 3
            cov += np.sum(current_pred != '.')

            sel_mut = current_pred != '.'
            if np.sum(sel_mut) > 0:
                sel_patho = current_pred == patho

                current_rate_mut = np.sum( mutmat.loc[current_wt][current_mut[sel_mut]] )
                current_rate_patho = np.sum( mutmat.loc[current_wt][current_mut[sel_patho]] )

                A += current_rate_mut
                P += current_rate_patho

    #print(P * 10e-9)

    if cov > 0:
        neutral = np.around( (A / cov ) * tot /len(list_positions), 3)
        full = np.around( (P / cov) * tot , 3)
        average = np.around( (P / cov) * tot/len(list_positions) , 3)
        coverage = np.around(cov/tot*100, 3)

    else:
        neutral = np.nan
        full = np.nan 
        average = np.nan 
        coverage = 0

    return average, full, neutral, coverage





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_mutations():
    """
    Compute genotypic mutation risk for all genes with a file in '../data/ntfull'
    Wrapper function to loop through all gene tables and compute the mutation risk
    Output: table with scores, written to file
    """

    genoDF = pd.DataFrame(columns=['Gene', 'PP2hdiv_avrg', 'PP2hdiv_full', 'PP2hdiv_neutral','PP2hdiv_coverage',
    								'EVE_avrg', 'EVE_full', 'EVE_neutral', 'EVE_coverage', 'AlphaMiss_avrg', 'AlphaMiss_full', 'AlphaMiss_neutral', 'AlphaMiss_coverage', 
    								'ESM1b_avrg', 'ESM1b_full', 'ESM1b_neutral', 'ESM1b_coverage'])

    counter = 0
    for i in glob.glob('../data/ntfull/*.txt.gz'):
        current_name = i.split('/')[-1].split('.')[0]

        current_df = pd.read_csv(i, header=0, index_col=None, sep='\t', compression='gzip')
        hdiv_a, hdiv_f, hdiv_n, hdiv_c = geno_riskscore(current_df, 'polyphen2_hdiv')
        #hvar_a, hvar_f, hvar_n, hvar_c = geno_riskscore(current_df, 'polyphen2_hvar')
        eve_a, eve_f, eve_n, eve_c = geno_riskscore(current_df, 'EVE')
        alpha_a, alpha_f, alpha_n, alpha_c = geno_riskscore(current_df, 'AlphaMissense')
        esm1b_a, esm1b_f, esm1b_n, esm1b_c = geno_riskscore(current_df, 'ESM1b')

        genoDF.loc[len(genoDF)] = (current_name, hdiv_a, hdiv_f, hdiv_n, hdiv_c, eve_a, eve_f, eve_n, eve_c, alpha_a, alpha_f, alpha_n, alpha_c, esm1b_a, esm1b_f, esm1b_n, esm1b_c)
        
        counter += 1
        if counter % 1000 == 0:
            print(counter, 'genes processed')
    
    genoDF.to_csv("../data/processed/mutations_genetic2.txt", header=True, index=False, sep='\t', na_rep='NA')

 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_mutations_unweighted():
    """
    Computes the 'unweighted' mutation risk scores, i.e. without nucleotide-specific mutation rates. 
    Output: table with scores, written to file
    """

    genoDF = pd.DataFrame(columns=['Gene', 'PP2hdiv_avrg', 'PP2hdiv_full', 'PP2hdiv_neutral','PP2hdiv_coverage',
                                    'EVE_avrg', 'EVE_full', 'EVE_neutral', 'EVE_coverage', 'AlphaMiss_avrg', 'AlphaMiss_full', 'AlphaMiss_neutral', 'AlphaMiss_coverage', 
                                    'ESM1b_avrg', 'ESM1b_full', 'ESM1b_neutral', 'ESM1b_coverage'])

    counter = 0
    for i in glob.glob('../data/ntfull/*.txt.gz'):
        current_name = i.split('/')[-1].split('.')[0]

        current_df = pd.read_csv(i, header=0, index_col=None, sep='\t', compression='gzip')
        hdiv_a, hdiv_f, hdiv_n, hdiv_c = geno_riskscore(current_df, 'polyphen2_hdiv', unweighted=True)
        eve_a, eve_f, eve_n, eve_c = geno_riskscore(current_df, 'EVE', unweighted=True)
        alpha_a, alpha_f, alpha_n, alpha_c = geno_riskscore(current_df, 'AlphaMissense', unweighted=True)
        esm1b_a, esm1b_f, esm1b_n, esm1b_c = geno_riskscore(current_df, 'ESM1b', unweighted=True)

        genoDF.loc[len(genoDF)] = (current_name, hdiv_a, hdiv_f, hdiv_n, hdiv_c, eve_a, eve_f, eve_n, eve_c, alpha_a, alpha_f, alpha_n, alpha_c, esm1b_a, esm1b_f, esm1b_n, esm1b_c)
        
        counter += 1
        if counter % 1000 == 0:
            print(counter, 'genes processed')
    
    genoDF.to_csv("../data/processed/mutations_genetic_unweighted2.txt", header=True, index=False, sep='\t', na_rep='NA')

 
 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_riskscore_strict(df):
    """
    Compute genetic mutation risk score only for sites that also include EVE predictions in dbSNFP4.7
    'polyphen2_hdiv', 'polyphen2_hvar', 'EVE', 'AlphaMissense', 'ESM1b'
    """

    mutmat = np.array([[0., 0.07054, 0.19986, 0.06059],
                       [0.09948, 0., 0.11304, 0.36897],
                       [0.36897, 0.11304, 0., 0.09948],
                       [0.06059, 0.19986, 0.07054, 0.]]) #* 10e-9
    mutmat = pd.DataFrame(data=mutmat, columns=['A', 'C', 'G', 'T'], index=['A', 'C', 'G', 'T'])

  
    P_alpha = 0.
    P_eve = 0.
    P_esm1b = 0.
    tot = 0.
    cov = 0.

    list_positions = sorted(list(set(df['Pos'])))
    for i in list_positions:
        current_df = df[df['Pos']==i]
        current_df = current_df[current_df['WT']!=current_df['MUT']]    # only mutations

        current_wt = list(set(current_df['WT']))[0]
        current_mut = np.array(list(current_df['MUT']))
        
        current_alpha = np.array(list(current_df['AlphaMissense']))
        current_esm1b = np.array(list(current_df['ESM1b']))
        current_eve = np.array(list(current_df['EVE']))


        sel_mut = (current_alpha != '.') * (current_esm1b != '.') * (current_eve != '.')

        tot += 3
        cov += np.sum(sel_mut)

        if np.sum(sel_mut) > 0:
            
            sel_patho_alpha = (current_alpha == 'P') * sel_mut
            sel_patho_esm1b = (current_esm1b == 'D') * sel_mut
            sel_patho_eve = (current_eve == 'P') * sel_mut

            P_alpha += np.sum( mutmat.loc[current_wt][current_mut[sel_patho_alpha]] )
            P_esm1b += np.sum( mutmat.loc[current_wt][current_mut[sel_patho_esm1b]] )
            P_eve += np.sum( mutmat.loc[current_wt][current_mut[sel_patho_eve]] )


    if cov > 0:
        average_alpha = np.around( (P_alpha / cov) * tot/len(list_positions), 3)
        full_alpha = np.around( (P_alpha / cov) * tot , 3)
        average_eve = np.around( (P_eve / cov) * tot/len(list_positions), 3)
        full_eve = np.around( (P_eve / cov) * tot , 3)
        average_esm1b = np.around( (P_esm1b / cov) * tot/len(list_positions), 3)
        full_esm1b = np.around( (P_esm1b / cov) * tot , 3)
        coverage = np.around(cov/tot*100, 3)

    else:
        average_alpha, average_eve, average_esm1b = np.nan, np.nan, np.nan
        full_alpha, full_eve, full_esm1b = np.nan, np.nan, np.nan 
        coverage = 0


    return average_alpha, full_alpha, average_eve, full_eve, average_esm1b, full_esm1b, coverage


#[35725680. 34142927. 35117474. 29772872.]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_mutations_strict():
    """
    Wrapper function to compute the 'strict' genotypic mutation risk for all files in ntfull
    Output: table of scores, written to file
    """


    genoDF = pd.DataFrame(columns=['Gene', 'coverage', 'EVE_avrg', 'EVE_full', 'AlphaMiss_avrg', 'AlphaMiss_full',
    								'ESM1b_avrg', 'ESM1b_full'])

    counter = 0
    for i in glob.glob('../data/ntfull/*.txt.gz'):
        current_name = i.split('/')[-1].split('.')[0]
        print(current_name)

        current_df = pd.read_csv(i, header=0, index_col=None, sep='\t', compression='gzip')
        alpha_a, alpha_f, eve_a, eve_f, esm1b_a, esm1b_f, coverage = geno_riskscore_strict(current_df)
        
        genoDF.loc[len(genoDF)] = (current_name, coverage, eve_a, eve_f, alpha_a, alpha_f, esm1b_a, esm1b_f)
        
        counter += 1
        if counter % 1000 == 0:
            print(counter, 'genes processed')
    

    print(genoDF)
    genoDF.to_csv("../data/processed/mutations_genetic_strict2.txt", header=True, index=False, sep='\t', na_rep='NA')

 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def count_mutations():
    """
    Get mutation counts by type in data
    """

    mutdict = {}

    counter = 0
    for i in glob.glob('../data/ntfull/*.txt.gz'):
        current_name = i.split('/')[-1].split('.')[0]

        current_df = pd.read_csv(i, header=0, index_col=None, sep='\t', compression='gzip')


        for i in range(len(current_df)):
            current_wt = current_df.iloc[i]['WT']
            current_mut = current_df.iloc[i]['MUT']

            if current_wt != current_mut:
               current_pair = "_".join([current_wt, current_mut])
               if current_pair in list(mutdict.keys()):
                   mutdict[current_pair] += 1
               else:
                   mutdict[current_pair] = 1

        counter += 1
        if counter % 1000 == 0:
            print(counter, 'genes processed')
    
    countDF = pd.DataFrame(columns=['mutation', 'count'])
    countDF.loc[len(countDF)] = ("AT>GC", mutdict['A_G']+mutdict['T_C'])    
    countDF.loc[len(countDF)] = ("GC>AT", mutdict['G_A']+mutdict['C_T'])    
    countDF.loc[len(countDF)] = ("AT>TA", mutdict['A_T']+mutdict['T_A'])    
    countDF.loc[len(countDF)] = ("GC>TA", mutdict['G_C']+mutdict['T_A'])    
    countDF.loc[len(countDF)] = ("AT>CG", mutdict['A_C']+mutdict['T_G'])    
    countDF.loc[len(countDF)] = ("GC>CG", mutdict['G_C']+mutdict['C_G'])    

    print(countDF)
    countDF.to_csv("tmp.count", header=True, index=False, sep='\t')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_alpha(fileIN):
    """
    Parse AlphaMissense AA data into per protein files (slower but runs on available hardware)
    """

    if not os.path.exists('../data/alpha/'):
        os.mkdir('../data/alpha/')

    tmplist=[]
    i = 0
    with gz.open(fileIN, 'rb') as f: 
        for line in f:
            current_line = line.decode('utf-8').rstrip() #.split('\t')
            if current_line[0:1] != "#":
                current_line = current_line.split('\t')
                current_uniprot = current_line[0]
                current_var = current_line[1]
                current_wt = current_var[0]
                current_mut = current_var[-1]
                current_pos = current_var.replace(current_wt,'').replace(current_mut, '')
                current_class = current_line[3]
                #print(current_uniprot, current_wt, current_mut, current_pos, current_class)

                current_write = str(current_uniprot) + '\t' + str(current_pos) + '\t' + str(current_wt) + '\t' + str(current_mut) + '\t' + str(current_class) + '\n'
                current_header = "uniprot" + '\t' + "pos" + '\t' + "wt" + '\t' + "mut" + '\t' + "class" + '\n'  
                
                if os.path.exists("../data/alpha/"+current_uniprot+".txt.gz"):
                    output = gz.open("../data/alpha/"+current_uniprot+".txt.gz", "at")
                else:
                    output = gz.open("../data/alpha/"+current_uniprot+".txt.gz", "wt")
                    output.writelines(current_header)

                output.writelines([current_write])
                output.close()

            i += 1

            if i % 10000000 == 0:
                print(i, "processed")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pheno_riskscore(df, predictor, uniprot):
    """
    Phenotypic mutation risk with 'AlphaMissense' and 'ESM1b'
    Computed are phenotypic mutations with:
    - the mTEL model and yeast ribosome energy parameters (ribo)
    - mutation frequencies derived from experimental peptides (plpl)
    - a random control where all phenotypic mutations are equally likely (rand)

    Input: dataframe with sequence, mutations and predictions; predictor; uniprot ID
    Output: - number of positions (length)
            - phenotypic mutation risk scores for different models as well as for 
            damaging and all (control) mutations
            - additional summary and qc statistics
    """

   
    resm = np.zeros(( 61, 20 ))
    resn = np.zeros(( 61, 20 ))

    if predictor in ['AlphaMissense', 'ESM1b']:
        if predictor == 'AlphaMissense':
            patho = 'P'							# (B)enign, (A)mbiguous, or likely (P)athogenic
            current_uniprot = uniprot
            current_enst = list(cid[cid['UniProtID']==current_uniprot]['TranscriptID'])[0]

            # get matching codons
            if current_enst in list(dict_seq.keys()):
                current_seq = dict_seq[current_enst]
                list_codons = []
                for j in range( int(np.floor(len(current_seq)/3.)) ):
                    list_codons.append( "".join(list(current_seq[(j*3):((j*3)+3)]))  )
                current_seq2 = np.array(list_codons)

            list_pos = []
            if (len(current_seq2) - 1) == np.max(df['pos']):
                aaseq = np.zeros(( len(current_seq2) ), dtype="S1")
                prob_ribo_all = 0.
                prob_ribo_del = 0.
                prob_plpl_all = 0.
                prob_plpl_del = 0.
                prob_rand_all = 0.
                prob_rand_del = 0.

                for i in range(len(df)):
                    current_pos = df.iloc[i]['pos'] - 1 #python indexing!!!
                    current_wt = current_seq2[current_pos]
                    current_AA_wt = df.iloc[i]['wt']
                    aaseq[current_pos] = current_AA_wt
                    if gencode.get(current_wt, '*') == current_AA_wt:
                        current_AA_mut = df.iloc[i]['mut']
                        current_pred = df.iloc[i]['class']
                        list_pos.append(current_pos)
                        #print(current_pos, current_wt, current_AA_wt, current_AA_mut, current_pred)

                        if current_AA_wt != current_AA_mut:
                            prob_ribo_all += mutmat.loc[current_wt][current_AA_mut]
                            prob_plpl_all += pil2.loc[current_wt][current_AA_mut]
                            prob_rand_all += mutmat2.loc[current_wt][current_AA_mut]
                            resn[codons_nonstop.index(current_wt), AA.index(current_AA_mut)] += 1

                            if current_pred == "pathogenic":
                                prob_ribo_del += mutmat.loc[current_wt][current_AA_mut]
                                prob_plpl_del += pil2.loc[current_wt][current_AA_mut]
                                prob_rand_del += mutmat2.loc[current_wt][current_AA_mut]
                                resm[codons_nonstop.index(current_wt), AA.index(current_AA_mut)] += 1

            else:
                prob_ribo_all = np.nan
                prob_ribo_del = np.nan
                prob_plpl_all = np.nan
                prob_plpl_del = np.nan
                prob_rand_all = np.nan
                prob_rand_del = np.nan

                #print(len(set(list_pos)), len(current_seq2), np.around(prob_ribo_all, 5), np.around(prob_ribo_del,5), np.around(prob_plpl_all, 5), np.around(prob_plpl_del,5)  )

        elif predictor == 'ESM1b':
            patho = 'D' 						# deleterious (D) vs. tolerated (T) variants
            current_uniprot = uniprot
            current_enst = list(cid[cid['UniProtID']==current_uniprot]['TranscriptID'])[0]

            # get matching codons
            if current_enst in list(dict_seq.keys()):
                current_seq = dict_seq[current_enst]
                list_codons = []
                for j in range( int(np.floor(len(current_seq)/3.)) ):
                    list_codons.append( "".join(list(current_seq[(j*3):((j*3)+3)]))  )
                current_seq2 = np.array(list_codons)

            list_idx = [x[2:] for x in list(df.columns)]    
            list_aa = [x[0] for x in list(df.columns)]
            df2 = pd.DataFrame(np.transpose(np.array(df)), columns=df.index, index=list_idx)

            list_pos = []
            if (len(current_seq2) - 1) == len(list(df.columns)):
                #print(current_seq2, len(current_seq2))

                aaseq = np.zeros(( len(current_seq2) ), dtype="S1")
                prob_ribo_all = 0.
                prob_ribo_del = 0.
                prob_plpl_all = 0.
                prob_plpl_del = 0.
                prob_rand_all = 0.
                prob_rand_del = 0.
 
                for i in range(len(df2)):
                    current_pos = int(list(df2.index)[i]) - 1 #python indexing!!!
                    current_wt = current_seq2[current_pos]
                    current_AA_wt = list_aa[i]
                    aaseq[current_pos] = current_AA_wt
                    if gencode.get(current_wt, '*') == current_AA_wt:
                        for aa in AA:

                            current_pred = df2.iloc[i][aa]
                            list_pos.append(current_pos)
                            #print(current_pos, current_wt, current_AA_wt, aa, current_pred)

                            if current_AA_wt != aa:

                                prob_ribo_all += mutmat.loc[current_wt][aa]
                                prob_plpl_all += pil2.loc[current_wt][aa]
                                prob_rand_all += mutmat2.loc[current_wt][aa]

                                resn[codons_nonstop.index(current_wt), AA.index(aa)] += 1

                                if current_pred < -7.5:         # threshold from paper
                                    prob_ribo_del += mutmat.loc[current_wt][aa]
                                    prob_plpl_del += pil2.loc[current_wt][aa]
                                    prob_rand_del += mutmat2.loc[current_wt][aa]

                                    resm[codons_nonstop.index(current_wt), AA.index(aa)] += 1

            else:
                prob_ribo_all = np.nan
                prob_ribo_del = np.nan
                prob_plpl_all = np.nan
                prob_plpl_del = np.nan
                prob_rand_all = np.nan
                prob_rand_del = np.nan

    res = np.zeros(( 61, 20 )) * np.nan
    res[resn>0] = resm[resn>0]/resn[resn>0]

    result_len = len(set(list_pos))
    result_probs = [np.around(prob_ribo_all, 3), np.around(prob_ribo_del,3), np.around(prob_plpl_all, 3), np.around(prob_plpl_del,3), np.around(prob_rand_all, 3), np.around(prob_rand_del,3)]
    result_stats = (res, resm, resn)


    return result_len, result_probs, result_stats





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pheno_stats(df, predictor):
    """
    Mean and standard deviation of codon/AA substitution pairs for predicted phenotypioc mutations
    """

    resm = np.zeros(( 61, 20 ))
    resn = np.zeros(( 61, 20 ))

    if predictor in ['polyphen2_hdiv', 'polyphen2_hvar', 'EVE', 'AlphaMissense', 'ESM1b']:
        if predictor == 'AlphaMissense':
            patho = 'P'							# (B)enign, (A)mbiguous, or likely (P)athogenic
        elif predictor == 'EVE':
            patho = 'P' 						# (B)enign, (U)ncertain, or (P)athogenic
        elif predictor == 'ESM1b':
            patho = 'D' 						# deleterious (D) vs. tolerated (T) variants
        elif predictor == 'polyphen2_hdiv':		
            patho = 'D' 						# "D" ("probably damaging"),"P" ("possibly damaging")  and "B" ("benign")
        elif predictor == 'polyphen2_hvar':			
            patho = 'D' 						# "D" ("probably damaging"), "P" ("possibly damaging") and "B" ("benign")

        P = 0.
        tot = 0.
        cov = 0.
        for i in range(len(df)):
            current_wt = df.iloc[i]['codon_wt']
            current_AA_wt = gencode.get(current_wt, "*")
            current_mut = df.iloc[i]['codon_mut']
            current_AA_mut = gencode.get(current_mut, "*")
            current_pred = df.iloc[i][predictor]

            if current_wt != current_mut and current_wt in list(codons_nonstop) and current_AA_mut != "*":
                current_prob = mutmat.loc[current_wt][current_AA_mut]               
                
                if current_pred == patho:
                    resm[codons_nonstop.index(current_wt), AA.index(current_AA_mut)] += 1

                if current_pred != '.':
                    resn[codons_nonstop.index(current_wt), AA.index(current_AA_mut)] += 1


    res = np.zeros(( 61, 20 )) * np.nan
    res[resn>0] = resm[resn>0]/resn[resn>0]

    return res, resn, resm




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pheno_mutations():
    """
    Wrapper function to compute the phenotypic mutation risk computed for all genes in ntfull
    This function loops through the annotated gene table and computes scores for each gene
    Output: a table of scores, written to file
    """

    phenoDF = pd.DataFrame(columns=['alpha_len', 'alpha_ribo_all', 'alpha_ribo_del', 'alpha_plpl_all', 'alpha_plpl_del', 'alpha_rand_all', 'alpha_rand_del',
                                    'esm1b_len', 'esm1b_ribo_all', 'esm1b_ribo_del', 'esm1b_plpl_all', 'esm1b_plpl_del', 'esm1b_rand_all', 'esm1b_rand_del'])

    list_glob = glob.glob('../data/ntfull/*.txt.gz') 	#[0:3]

    if not os.path.exists('../data/esm1b/ALL_hum_isoforms_ESM1b_LLR/'):
        print("WARNING: ESM1b scores have to be downloaded and copied here!")
        break
    

    alpha_all = np.zeros((61, 20, len(list_glob) )) * np.nan
    alpha_all_n = np.zeros(( 61, 20 ))
    alpha_all_m = np.zeros(( 61, 20 ))
   
    esm1b_all = np.zeros((61, 20, len(list_glob) )) * np.nan
    esm1b_all_n = np.zeros(( 61, 20 ))
    esm1b_all_m = np.zeros(( 61, 20 ))

    counter = 0
    for i in list_glob:		
        current_name = i.split('/')[-1].split('.')[0]
        if current_name in list(cid['Name']):
            current_uniprot = list(cid[cid['Name']==current_name]['UniProtID'])[0]

            if os.path.exists('../data/alpha/'+current_uniprot+'.txt.gz'):
                current_alpha = pd.read_csv('../data/alpha/'+current_uniprot+'.txt.gz', header=0, index_col=None, sep='\t', compression='gzip')
                alpha_len, alpha_probs, alpha_stats = pheno_riskscore(current_alpha, 'AlphaMissense', current_uniprot)

                alpha_ribo_all = alpha_probs[0] 
                alpha_ribo_del = alpha_probs[1]
                alpha_plpl_all = alpha_probs[2] 
                alpha_plpl_del = alpha_probs[3] 
                alpha_rand_all = alpha_probs[4] 
                alpha_rand_del = alpha_probs[5]

                alpha_all[:,:,counter] = alpha_stats[0]
                alpha_all_m += alpha_stats[1]
                alpha_all_n += alpha_stats[2]

            else:
                alpha_len = np.nan
                alpha_ribo_all = np.nan
                alpha_ribo_del = np.nan
                alpha_plpl_all = np.nan
                alpha_plpl_del = np.nan
                alpha_rand_all = np.nan 
                alpha_rand_del = np.nan


            if os.path.exists('../data/esm1b/ALL_hum_isoforms_ESM1b_LLR/'+current_uniprot+'_LLR.csv'):
                current_esm1b = pd.read_csv('../data/esm1b/ALL_hum_isoforms_ESM1b_LLR/'+current_uniprot+'_LLR.csv', header=0, index_col=0, sep=',')
                esm1b_len, esm1b_probs, esm1b_stats  = pheno_riskscore(current_esm1b, 'ESM1b', current_uniprot)

                esm1b_ribo_all = esm1b_probs[0]
                esm1b_ribo_del = esm1b_probs[1] 
                esm1b_plpl_all = esm1b_probs[2] 
                esm1b_plpl_del = esm1b_probs[3] 
                esm1b_rand_all = esm1b_probs[4] 
                esm1b_rand_del = esm1b_probs[5]

                esm1b_all[:,:,counter] = alpha_stats[0]
                esm1b_all_m += alpha_stats[1]
                esm1b_all_n += alpha_stats[2]

            else:
                esm1b_len = np.nan
                esm1b_ribo_all = np.nan
                esm1b_ribo_del = np.nan
                esm1b_plpl_all = np.nan
                esm1b_plpl_del = np.nan
                esm1b_rand_all = np.nan 
                esm1b_rand_del = np.nan


        phenoDF.loc[len(phenoDF)] = (alpha_len, alpha_ribo_all, alpha_ribo_del, alpha_plpl_all, alpha_plpl_del, alpha_rand_all, alpha_rand_del,
                                     esm1b_len, esm1b_ribo_all, esm1b_ribo_del, esm1b_plpl_all, esm1b_plpl_del, esm1b_rand_all, esm1b_rand_del)
       

        alpha_stats, alpha_n, alpha_m = pheno_stats(current_df, 'AlphaMissense')
        alpha_all[:,:,counter] = alpha_stats
        alpha_all_n += alpha_n
        alpha_all_m += alpha_m

        esm1b_stats, esm1b_n, esm1b_m = pheno_stats(current_df, 'ESM1b')
        esm1b_all[:,:,counter] = esm1b_stats
        esm1b_all_n += esm1b_n
        esm1b_all_m += esm1b_m

        counter += 1
        if counter % 1000 == 0:
            print(counter, 'genes processed')

    print(phenoDF)
    print(np.shape(alpha_all))

    alpha_mean = pd.DataFrame(data=np.around(np.nanmean(alpha_all, axis=2),5), index=codons_nonstop, columns=AA)
    alpha_mean.to_csv("../data/processed/pheno_alpha_mean.txt", header=True, index=True, sep='\t', na_rep="NA")
    alpha_sd = pd.DataFrame(data=np.around(np.nanstd(alpha_all, axis=2),5), index=codons_nonstop, columns=AA)
    alpha_sd.to_csv("../data/processed/pheno_alpha_sd.txt", header=True, index=True, sep='\t', na_rep="NA")
    np.savetxt("../data/processed/counts_alpha_n.txt", alpha_all_n, fmt='%i')
    np.savetxt("../data/processed/counts_alpha_m.txt", alpha_all_m, fmt='%i')

    esm1b_mean = pd.DataFrame(data=np.around(np.nanmean(esm1b_all, axis=2),5), index=codons_nonstop, columns=AA)
    esm1b_mean.to_csv("../data/processed/pheno_esm1b_mean.txt", header=True, index=True, sep='\t', na_rep="NA")
    esm1b_sd = pd.DataFrame(data=np.around(np.nanstd(esm1b_all, axis=2),5), index=codons_nonstop, columns=AA)
    esm1b_sd.to_csv("../data/processed/pheno_esm1b_sd.txt", header=True, index=True, sep='\t', na_rep="NA")
    np.savetxt("../data/processed/counts_esm1b_n.txt", esm1b_all_n, fmt='%i')
    np.savetxt("../data/processed/counts_esm1b_m.txt", esm1b_all_m, fmt='%i')


    phenoDF.to_csv("../data/processed/mutations_phenotypic.txt", header=True, index=False, sep='\t', na_rep='NA')

 





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def important_features(data, preOUT):
    """
    Fit Random Forest (RF) and linear models (LM), leave-one-out crossvalidated, get feature importances
    first column is the response variable to be fitted/predicted
    """

    data = data.loc[ np.sum(np.isnan(data) == True, axis=1) == 0 ]
    print(np.shape(data))

    y = np.array(data[data.columns[0]], dtype=float)
    X = np.array(data[data.columns[1:]])
    

    # 1. only fullfit/predict
    regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100)
    regr.fit(X, y)
    pred = regr.predict( X )
    resDF = pd.DataFrame({'true': list(y), 'pred': list(pred) })
    resDF.to_csv("tmp_fullfit.txt", header=True, index=False, sep='\t')
    np.savetxt("tmp_fullfit_fa.txt", regr.feature_importances_)

    rr = np.array(resDF)
    sel = (np.isnan(rr[:,0]) == False) * (np.isnan(rr[:,1]) == False)
    cc = np.corrcoef(rr[sel,0], rr[sel,1])[0,1]
    print("R =",np.around(cc, 4), np.shape(rr))
    print(regr.feature_importances_)


    #linear reg
    model = LinearRegression()
    model.fit(X, y)
    lmpred = model.predict( X ) 
    sel = (np.isnan(y) == False) * (np.isnan(lmpred) == False)
    cc = np.corrcoef(y[sel], lmpred[sel])[0,1]
    print("R =",np.around(cc, 4), np.shape(rr))
    #print(model.coef_)


    # 2. leave-one-out cross-validation for predictive power
    y = np.array(data[data.columns[0]], dtype=float)
    X = np.array(data[data.columns[1:]])

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    pred = np.zeros(( len(y) ))
    fa = np.zeros(( len(y), np.shape(X)[1] ))
    pred2 = np.zeros(( len(y) ))
    fa2 = np.zeros(( len(y), np.shape(X)[1] ))
    for i in range( len(y) ):          # len(pdblist)
        sel_train = np.ones(( len(y) ), dtype=bool)
        sel_train[i] = False

        y_train = y[sel_train]
        X_train = X[sel_train,:]
        X_pred = X[i,:].reshape(1, -1)

        regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100)
        regr.fit(X_train, y_train)
        pred[i] = regr.predict( X_pred )

        fa[i,:] = regr.feature_importances_


        model = LinearRegression()
        model.fit(X_train, y_train)
        pred2[i] = model.predict( X_pred )
        fa2[i,:] = model.coef_

    resDF = pd.DataFrame({'True': list(np.around(y,6)), 'RF': list(np.around(pred,6)) , 'LR': list(np.around(pred2,6)) })
    resDF.to_csv("../data/processed/"+preOUT+"_pred.txt", header=True, index=False, sep='\t')

    faRF = pd.DataFrame(data=np.around(fa, 5), columns=list(data.columns[1:]))
    faRF.to_csv("../data/processed/"+preOUT+"_features_RF.txt", header=True, index=False, sep='\t')

    faLR = pd.DataFrame(data=fa2, columns=list(data.columns[1:]))
    faLR.to_csv("../data/processed/"+preOUT+"_features_LR.txt", header=True, index=False, sep='\t')

    #np.savetxt("../data/processed/features_RF_"+preOUT+".txt", fa)
    #np.savetxt("../data/processed/features_LR_"+preOUT+".txt", fa2)

    print(preOUT, "[done]")













if __name__ == '__main__':

    # define codon list as global var
    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   
    gencode = pd.DataFrame(
        {'codon':['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
        'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
        'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
        'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'],
        'aa':['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'] })
    gencode = dict(gencode.values)

    nt = ['A', 'C', 'G', 'T']
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']




    # load gene annotations
    parse_gff('../data/human/gencode.v46.annotation.gff3.gz')
    gff = pd.read_csv("../data/processed/gff.txt", header=0, index_col=False, sep='\t', low_memory=False)

    # load transcript sequences
    dict_seq = cdsseq_dict("../data/human/gencode.v46.pc_transcripts.fa")


    # parse canonical transcript ids
    canonical_enst = extract_uniprot_ids()
    canonical_enst = pd.read_csv("../data/processed/uniprot_canonical.txt", header=0, index_col=False, sep='\t')
    cid = pd.read_csv("../data/human/canonical_id.txt", header=0, index_col=False, sep='\t')


    # parse nt tables
    prepare_nt_tables(canonical_enst, dict_seq, gff)	# DF of pos/seq

    # parse SNPdb data
    parse_snpdb('../variation/freq.vcf.gz')			# SNPdb file from NIH

    # parse dbNSFP variant prediction data (loop through chromosome files)
    for i in glob.glob('/home/sebastian/dbNSFP4.7a/dbNSFP4.7a_variant.chr*.gz'):
        variantpred_by_transcript(i)
    


    # merge: mutation table + dbNSFP + SNPdb
    cid = pd.read_csv("../data/human/canonical_id.txt", header=0, index_col=False, sep='\t')
    print(len(cid))
    list_seq = []
    for ix, i in enumerate(list(cid['TranscriptID'])):
        if i in list(dict_seq.keys()):
            list_seq.append(ix)
    print(len(list_seq))
    cid = cid.iloc[list_seq]

    pool = mp.Pool(processes=20)        # parallelize computation
    for i in range(len(cid)):
        current_ensg = cid.iloc[i]['GeneID']
        current_enst = cid.iloc[i]['TranscriptID']
        current_name = cid.iloc[i]['Name']
        current_seq = dict_seq[current_enst]
        pool.apply_async(mut_data_parallel, [(current_ensg, current_enst, current_name, current_seq)] ) 
    pool.close()
    pool.join()




    # Compute genotypic mutation scores
    normalize_mutmat()
    geno_mutations()
    geno_mutations_unweighted()
    geno_mutations_strict()
    count_mutations()




    # Normalize mutation matrices for phenotypic mutations to 1 in 400 codons
    cccd = pd.read_csv("../data/processed/composition_codons.txt", header=0, index_col=None, sep='\t')
    cccd = np.array(cccd['count'])/np.mean(cccd['count'])
    ccaa = pd.read_csv("../data/processed/composition_AA.txt", header=0, index_col=None, sep='\t')
    ccaa = np.array(ccaa['count'])/np.mean(ccaa['count'])

    # Ribosome mTEL model
    mutmat = pd.read_csv("../data/processed/pheno_human_yeastE.txt", header=0, index_col=0, sep='\t')
    m2 = np.array(mutmat) * (0.0025 / np.mean(cccd * np.sum(mutmat, axis=1)) )                    # 0.0025
    mutmat = pd.DataFrame(data=m2, index=mutmat.index, columns=mutmat.columns)
    mutmat2 = np.copy(np.array(mutmat))     # same
    mutmat2[mutmat2>0] = 1
    mutmat2 = mutmat2 * (0.0025 / np.mean(cccd * np.sum(mutmat2, axis=1)) )                   # 0.0025
    mutmat2 = pd.DataFrame(data=mutmat2, index=mutmat.index, columns=mutmat.columns)

    # Mutation freqs from peptide data
    pilpel = pd.read_csv("../data/pheno/pilpel.txt", header=0, index_col=0, sep='\t')
    pilpel.loc["L"] = np.zeros(( 19 )) * np.nan     # add missing "L"
    pilpel["L"] = np.zeros(( 20 )) * np.nan         # add missing "L"
    pilpel = pilpel.loc[sorted(pilpel.index)]
    pilpel = pilpel[sorted(pilpel.columns)]
    pMAT = np.array(pilpel)
    pMAT[np.isnan(pMAT)==True] = np.percentile(pMAT[np.isnan(pMAT)==False], 25)  #np.nanmedian(pMAT)
    np.fill_diagonal(pMAT, 0)
    pMAT = pMAT * (0.0025 / np.mean(ccaa * np.sum(pMAT, axis=1)) )                        # 0.0025
    pilpel = pd.DataFrame(data=pMAT, index=pilpel.index, columns=pilpel.columns)
    #pilpel = pd.read_csv("normalized_pepAA.txt", header=0, index_col=0, sep='\t')
    pil2 = pd.read_csv("../data/processed/peptides_mutations.txt", header=0, index_col=0, sep='\t')
    pil2 = pil2 * (0.0025 / np.mean(cccd * np.sum(pil2, axis=1)) )                       # 0.0025
    #pil2 = pd.read_csv("normalized_pepCd.txt", header=0, index_col=0, sep='\t')

    print(np.mean(np.sum(mutmat, axis=1)))
    print(np.mean(np.sum(mutmat, axis=1)))
    print(np.mean(np.sum(pilpel, axis=1)))
    print(np.mean(np.sum(pil2, axis=1)))


    # Compute phenotypic mutation risk
    parse_alpha("../data/AlphaMissense_aa_substitutions.tsv.gz")
    pheno_mutations()








    # Compile feature matrix

    # Consensus genotypic mutation risk as average of AlphaMissense and ESM1b scores
    geno = pd.read_csv("../data/processed/mutations_genetic2.txt", header=0, index_col=None, sep='\t')
    geno = geno[['Gene', 'AlphaMiss_avrg', 'AlphaMiss_full', 'AlphaMiss_coverage', 
    				'ESM1b_avrg', 'ESM1b_full', 'ESM1b_coverage']]
    s =  np.absolute(1 - ( geno['ESM1b_avrg']/( (geno['AlphaMiss_avrg']+geno['ESM1b_avrg'])/2 ) ))
    sel = s < 0.5
    geno = geno.loc[sel]
    geno['score_avrg'] = (geno['AlphaMiss_avrg']+geno['ESM1b_avrg'])/2.
    geno['score_full'] = (geno['AlphaMiss_full']+geno['ESM1b_full'])/2.
    geno.to_csv("../data/processed/mutations_genetic_score.txt", header=True, index=False, sep='\t')
    #print([len(s), np.nansum(s > 0), np.nansum(s < 0.5)])  

    # intrinsic disorder
    diso = pd.read_csv("../data/human/mobidb_prediction-disorder-th_50", header=0, index_col=None, sep='\t')
    canonical_id = pd.read_csv("../data/human/canonical_id.txt", header=0, index_col=None, sep='\t')

    disoDF = pd.DataFrame(columns=['Name', 'UniProt', 'diso'])
    for i in range(len(diso)):
        current_uniprot = diso.iloc[i]['UniProt']
        current_diso = np.array([int(x) for x in list(diso.iloc[i]['diso'])])
        if current_uniprot in list(canonical_enst['uniprot']):
            current_enst = list(canonical_enst[canonical_enst['uniprot']==current_uniprot]['canonical'])[0]
            if current_enst in list(canonical_id['TranscriptID']):
                current_name = list(canonical_id[canonical_id['TranscriptID']==current_enst]['Name'])[0]
                current_diso_pct = np.around(np.mean(current_diso), 3)
                disoDF.loc[len(disoDF)] = (current_name, current_uniprot, current_diso_pct)
    disoDF.to_csv("../data/processed/disorder_percent.txt", header=True, index=False, sep='\t')

    geno = pd.read_csv("../data/processed/mutations_genetic.txt", header=0, index_col=None, sep='\t')
    list_diso = []
    for i in range(len(geno)):
        current_name = geno.iloc[i]['Gene']
        if current_name in list(disoDF['Name']):
            current_diso = np.around(list(disoDF[disoDF['Name']==current_name]['diso'])[0], 3)
            list_diso.append(current_diso)
        else:
            list_diso.append('None')

    geno['diso'] = list_diso

    geno2 = geno[['Gene', 'AlphaMiss_avrg', 'ESM1b_avrg', 'diso']]
    print(geno2)
    geno2.to_csv("../data/processed/mutations_genetic_disorder.txt", header=True, index=False, sep='\t', na_rep="NA")


    # average number of hydrophobic contacts in protein structure
    # this requires running CSU on the human AlphaFold database first to generate the input files
    #list_csu = glob.glob('../data/csu/csu*.txt')
    #cmDF = pd.DataFrame(columns=['id', 'avrg_degree'])
    #for i in list_csu:
    #    current_id = i.split('-')[1]
    #    current_data = np.loadtxt(i)
    #    current_avrg_degree = np.around( np.mean(np.sum(current_data, axis=1)), 2)
    #    cmDF.loc[len(cmDF)] = (current_id, current_avrg_degree)   
    #cmDF.to_csv("../data/processed/avrg_contacts.txt", header=True, index=False, sep='\t')

  
    # load mutation risk data
    geno_all = pd.read_csv("../data/processed/mutations_genetic2.txt", header=0, index_col=None, sep='\t')
    geno = pd.read_csv("../data/processed/mutations_genetic_score.txt", header=0, index_col=None, sep='\t')
    pheno = pd.read_csv("../data/processed/mutations_phenotypic.txt", header=0, index_col=None, sep='\t')
    pheno.index = list(geno_all['Gene'])
    cid = pd.read_csv("../data/human/canonical_id.txt", header=0, index_col=False, sep='\t')


    # load protein features
    lenDF = pd.read_csv("../data/processed/enst_length.txt", header=0, index_col=None, sep='\t')
    disoDF = pd.read_csv("../data/processed/disorder_percent.txt", header=0, index_col=None, sep='\t')
    contactDF = pd.read_csv("../data/processed/avrg_contacts.txt", header=0, index_col=None, sep='\t')
    ppiDF = pd.read_csv("../data/processed/ppi_degree.txt", header=0, index_col=None, sep='\t')
    aggDF = pd.read_csv("../data/processed/human_tango.txt", header=0, index_col=None, sep='\t')
    tissueEXP = pd.read_csv("../data/human/rna_tissue_consensus.tsv", header=0, sep='\t', index_col=None)


    # Parse RNA expression data in correct order and format
    #expdat = np.zeros(( len(list(set(tissueEXP['Gene']))), len(list(set(tissueEXP['Tissue'])))   )) * np.nan
    #for ix, i in enumerate(sorted(list(set(tissueEXP['Gene'])))):
    #    d2 = tissueEXP[tissueEXP['Gene']==i]
    #    for jx, j in enumerate(sorted(list(set(tissueEXP['Tissue'])))):
    #        if j in list(d2['Tissue']):
    #            expdat[ix, jx] = list(d2[d2['Tissue']==j]['nTPM'])[0]
    #    if ix%2000 == 0:
    #        print(ix)
    #expDF = pd.DataFrame(data=expdat, index=sorted(list(set(tissueEXP['Gene']))), columns=sorted(list(set(tissueEXP['Tissue']))))
    #s = np.sum(np.isnan(expDF)==True, axis=0) == 0
    #expDF = expDF[expDF.columns[s]]
    #expDF.to_csv("../data/processed/rna_tissue_consensus_mat.txt", header=True, index=True, sep='\t')
    expDF = pd.read_csv("../data/processed/rna_tissue_consensus_mat.txt", header=0, index_col=0, sep='\t')


    # Prase RNA half-life data in correct order and format
    #hlDAT = pd.read_csv("../data/AvgKdegs_genes_v1.1.csv", header=0, index_col=None, sep=',')
    #hldat = np.zeros(( len(list(set(hlDAT['feature_ID']))), len(list(set(hlDAT['cell_line']))))) * np.nan
    #for ix, i in enumerate(sorted(list(set(hlDAT['feature_ID'])))):
    #    d2 = hlDAT[hlDAT['feature_ID']==i]
    #    for jx, j in enumerate(sorted(list(set(hlDAT['cell_line'])))):
    #        if j in list(d2['cell_line']):
    #            hldat[ix, jx] = list(d2[d2['cell_line']==j]['avg_halflife'])[0]
    #hlDF = pd.DataFrame(data=hldat, index=sorted(list(set(hlDAT['feature_ID']))), columns=sorted(list(set(hlDAT['cell_line']))))
    #hlDF.to_csv("../data/processed/rna_halflife_mat.txt", header=True, index=True, sep='\t')
    hlDF = pd.read_csv("../data/processed/rna_halflife_mat.txt", header=0, index_col=0, sep='\t')
    #print(hlDF)


    # Parse feature matrix
    resultDF = pd.DataFrame(columns=['Gene', 'geno_full', 'geno_avrg', 'pheno_full', 'pheno_avrg', 'length', 'diso', 'contacts', 'PPI_HPA', 'PPI_EBI', 'aggregation'])
    for i in list(geno['Gene']):
        if i in list(cid['Name']):
            current_uniprot = list(cid[cid['Name']==i]['UniProtID'])[0]
            current_enst = list(cid[cid['Name']==i]['TranscriptID'])[0]

            # score
            score_full = list(geno[geno['Gene']==i]['score_full'])[0]
            score_avrg = list(geno[geno['Gene']==i]['score_avrg'])[0]
            pheno_full = (pheno.loc[i]['alpha_ribo_del'] + pheno.loc[i]['esm1b_ribo_del'])/2.
            pheno_avrg = ((pheno.loc[i]['alpha_ribo_del']/pheno.loc[i]['alpha_len']) + (pheno.loc[i]['esm1b_ribo_del']/pheno.loc[i]['esm1b_len']) )/2.

            # length
            if current_enst in list(lenDF['enst']):
                current_len = list(lenDF[lenDF['enst']==current_enst]['length'])[0]
            else:
                current_len = np.nan

            # protein disorder (mobiDB)
            if i in list(disoDF['Name']):
                current_diso = list(disoDF[disoDF['Name']==i]['diso'])[0]
            else:
                current_diso = np.nan

            # protein contacts
            if current_uniprot in list(contactDF['id']):
                current_cm = list(contactDF[contactDF['id']==current_uniprot]['avrg_degree'])[0]
            else:
                current_cm = np.nan

            # ppi degree
            if i in list(ppiDF['gene']):
                current_ppi_hpa = list(ppiDF[ppiDF['gene']==i]['HPA'])[0]
                current_ppi_ebi = list(ppiDF[ppiDF['gene']==i]['EBI'])[0]
            else:
                current_ppi_hpa = np.nan
                current_ppi_ebi = np.nan

            # aggregation
            if current_uniprot in list(aggDF['uniprot']):
                current_agg = list(aggDF[aggDF['uniprot']==current_uniprot]['agg'])[0]
            else:
                current_agg = np.nan

            resultDF.loc[len(resultDF)] = (i, score_full, score_avrg, pheno_full, pheno_avrg, current_len,     current_diso, current_cm, current_ppi_hpa, current_ppi_ebi, current_agg)


    result_exp = np.zeros(( len(resultDF), len(list(expDF.columns)) )) * np.nan
    result_hl = np.zeros(( len(resultDF), len(list(hlDF.columns)) )) * np.nan
    # mRNA expression
    # mRNA half-life
    for ix, i in enumerate(list(resultDF['Gene'])):
        current_ensg = list(cid[cid['Name']==i]['GeneID'])[0]
        if current_ensg in list(expDF.index):
            result_exp[ix,:] = np.around(np.array(expDF.loc[current_ensg]), 3)
        if i in list(hlDF.index):
            result_hl[ix,:] = np.around(np.array(hlDF.loc[i]), 3)
        if current_ensg in list(msb_mRNA.index):
            result_msb_mrna[ix,:] = np.around(np.array(msb_mRNA.loc[current_ensg]), 3)
        if current_ensg in list(msb_prot.index):
            result_msb_prot[ix,:] = np.around(np.array(msb_prot.loc[current_ensg]), 3)


    resultDF_exp = pd.DataFrame(data=result_exp, columns=expDF.columns)
    result_exp2 = np.copy(result_exp)
    result_exp2[result_exp2 > 0] = np.log(result_exp2[result_exp2 > 0])     # !! log !!
    resultDF_exp2 = pd.DataFrame(data=result_exp2, columns=expDF.columns)
    resultDF_hl = pd.DataFrame(data=result_hl, columns=hlDF.columns)
    resultDF_msb_mrna = pd.DataFrame(data=result_msb_mrna, columns=msb_mRNA.columns)
    resultDF_msb_mrna.to_csv("../data/processed/all_msb_mrna.txt", header=True, index=False, sep='\t', na_rep="NA")
    resultDF_msb_prot = pd.DataFrame(data=result_msb_prot, columns=msb_prot.columns)
    resultDF_msb_prot.to_csv("../data/processed/all_msb_prot.txt", header=True, index=False, sep='\t', na_rep="NA")

    allDF = pd.concat((resultDF, resultDF_exp, resultDF_hl), axis=1)
    allDF2 = pd.concat((resultDF, resultDF_exp2, resultDF_hl), axis=1)

    print(allDF)
    allDF.to_csv("../data/processed/all_reg_vars.txt", header=True, index=False, sep='\t', na_rep="NA")
    allDF2.to_csv("../data/processed/all_reg_vars_log.txt", header=True, index=False, sep='\t', na_rep="NA")

  

    arv = pd.read_csv("../data/processed/all_reg_vars_log.txt", header=0, index_col=None, sep='\t', na_values='NA')
    #print(arv.columns)
 

    list_el = ['adipose tissue', 'adrenal gland', 'appendix', 'bone marrow', 'breast', 'cerebral cortex', 'cervix', 'choroid plexus', 'colon', 'duodenum', 'endometrium', 'epididymis', 'esophagus', 'fallopian tube', 'gallbladder', 'heart muscle', 'kidney', 'liver', 'lung', 'lymph node', 'ovary', 'pancreas', 'parathyroid gland', 'placenta', 'prostate', 'rectum', 'retina', 'salivary gland', 'seminal vesicle', 'skeletal muscle', 'skin', 'small intestine', 'smooth muscle', 'spleen', 'stomach', 'testis', 'thymus', 'thyroid gland', 'tongue', 'tonsil', 'urinary bladder']
    list_hl = ['CH22', 'Calu3', 'HEK293', 'HeLa', 'K562', 'MCF7', 'MDAMB', 'MOLM', 'MV411', 'Nalm6', 'U2OS']


    features_phenofull_all = ['pheno_full', 'length', 'diso', 'contacts', 'PPI_EBI', 'aggregation'] + list_el + list_hl
    features_genofull_all = ['geno_full', 'length', 'diso', 'contacts', 'PPI_EBI', 'aggregation'] + list_el + list_hl
    features_phenoavrg_all = ['pheno_avrg', 'diso', 'contacts', 'PPI_EBI', 'aggregation'] + list_el + list_hl
    features_genoavrg_all = ['geno_avrg', 'diso', 'contacts', 'PPI_EBI', 'aggregation'] + list_el + list_hl

    features_phenoavrg_prot = ['pheno_avrg', 'diso', 'contacts', 'PPI_EBI', 'aggregation']
    features_genoavrg_prot = ['geno_avrg', 'diso', 'contacts', 'PPI_EBI', 'aggregation']
    features_phenoavrg_el = ['pheno_avrg'] + list_el
    features_genoavrg_el = ['geno_avrg'] + list_el
    features_phenoavrg_hl = ['pheno_avrg'] + list_hl
    features_genoavrg_hl = ['geno_avrg'] + list_hl

    features_phenoavrg_best = ['pheno_avrg',  'contacts', 'cerebral cortex']
    features_genoavrg_best = ['geno_avrg',  'contacts', 'cerebral cortex']


    important_features(arv[features_phenofull_all], "mutrate_phenofull_all")
    important_features(arv[features_genofull_all], "mutrate_genofull_all")
    important_features(arv[features_phenoavrg_all], "mutrate_phenoavrg_all")
    important_features(arv[features_genoavrg_all], "mutrate_genoavrg_all")

    important_features(arv[features_phenoavrg_prot], "mutrate_phenoavrg_prot")
    important_features(arv[features_genoavrg_prot], "mutrate_genoavrg_prot")
    important_features(arv[features_phenoavrg_el], "mutrate_phenoavrg_el")
    important_features(arv[features_genoavrg_el], "mutrate_genoavrg_el")
    important_features(arv[features_phenoavrg_hl], "mutrate_phenoavrg_hl")
    important_features(arv[features_genoavrg_hl], "mutrate_genoavrg_hl")

    important_features(arv[features_phenoavrg_best], "mutrate_phenoavrg_best")
    important_features(arv[features_genoavrg_best], "mutrate_genoavrg_best")



    resultDF_msb_mrna = pd.read_csv("../data/processed/all_msb_mrna.txt", header=0, index_col=None, sep='\t') 
    resultDF_msb_prot = pd.read_csv("../data/processed/all_msb_prot.txt", header=0, index_col=None, sep='\t')
    resultDF_msb_prot /= np.mean(resultDF_msb_prot)

    DF_msb_prot = pd.concat((resultDF, resultDF_msb_prot), axis=1)
    important_features(DF_msb_prot[['pheno_avrg'] + list(resultDF_msb_prot.columns)], "mutrate_phenoavrg_msb_protel")
    important_features(DF_msb_prot[['geno_avrg'] + list(resultDF_msb_prot.columns)], "mutrate_genoavrg_msb_protel")

    DF_msb_mrna = pd.concat((resultDF, resultDF_msb_mrna), axis=1)
    important_features(DF_msb_mrna[['pheno_avrg'] + list(resultDF_msb_mrna.columns)], "mutrate_phenoavrg_msb_el")
    important_features(DF_msb_mrna[['geno_avrg'] + list(resultDF_msb_mrna.columns)], "mutrate_genoavrg_msb_el")

    DF_msb_prot_best = pd.concat((arv, resultDF_msb_prot), axis=1)
    important_features(DF_msb_prot_best[['pheno_avrg', 'contacts', 'Brain'] ], "mutrate_phenoavrg_msb_protel_best")
    important_features(DF_msb_prot_best[['geno_avrg', 'contacts', 'Brain'] ], "mutrate_genoavrg_msb_protel_best")

    DF_msb_mrna_best = pd.concat((arv, resultDF_msb_mrna), axis=1)
    important_features(DF_msb_mrna[['pheno_avrg', 'contacts', 'Brain'] ], "mutrate_phenoavrg_msb_el_best")
    important_features(DF_msb_mrna[['geno_avrg', 'contacts', 'Brain'] ], "mutrate_genoavrg_msb_el_best")
