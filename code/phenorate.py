import os, sys
import numpy as np
import pandas as pd
import pickle
import itertools as it
from Bio import SeqIO



#--------------------------------------------------------------------------------------------
def arrival_probs(ACdf, tRNAdf):
    """
    Compute tRNA arrival probabilities for all codons. 
    Codons without cognate tRNA yield 0. 
    Arrival probability simplifies to relative abundance
    ACdf: codon-anticodon-AA dataframe
    tRNAdf: input tRNA counts and abundances
    """

    lmbda = np.zeros(( len(ACdf) )) * np.nan
    for cx, codon in enumerate(list(ACdf['codon'])):
        current_anticodon = list(ACdf[ACdf['codon']==codon]['anticodon'])[0]
        if current_anticodon in list(tRNAdf['Anticodon']):
            current_tRNA_count = list(tRNAdf[tRNAdf['Anticodon']==current_anticodon]['Count'])[0]
        else:
            current_tRNA_count = 0
        lmbda[cx] = current_tRNA_count 

    pa = lmbda / np.sum(lmbda)
    pa_n = np.zeros(( len(pa) ))
    pa_n[pa>0] = 1 - pa[pa>0]

    resultDF = ACdf.copy()
    resultDF['pa'] = pa
    resultDF['pa_n'] = pa_n
    resultDF['lambda'] = lmbda

    return resultDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_pb(pos_weights, nt_energies):
    """
    Codon - anticodon binding energies computed from 
    nt interaction energy matrix and positional weights
    """

    x = list(it.product((0,1,2,3),(0,1,2,3),(0,1,2,3)))
    y = np.array(list(it.product(x,x)))
    z = ["".join(a) for a in list(np.array(nt)[np.array(x)]) ]

    ENG = np.reshape( np.sum( pos_weights * nt_energies[y[:,0,:], y[:,1,:][:,::-1] ], axis=1), (64, 64) )
  
    delta_A = np.exp(-ENG)  # stability of codon/anticodon interaction
    delta_A = np.transpose(delta_A) # transpose 2x to normalize along axis=1
    delta_A /= np.sum(delta_A, axis=0)
    pbMAT = np.transpose(delta_A)

    E2 = pd.DataFrame(data=pbMAT, columns=z, index=z)

    return E2



#--------------------------------------------------------------------------------------------
def binding_probs(ACdf, E2):
    """
    Rates of phenotypic mutations for each codon and to every amino acid
    Normalized to an overall rate of 1 phenotypic mutation for every 400 codons
    """
    
    er = 0.0025		# 1 in 400 error 

    list_codons = sorted(list(ACdf['codon']))
    list_tRNAs = list(ACdf[ACdf['lambda']>0]['anticodon'])

    list_aa = sorted(list(set(ACdf['AA'])))

    tRNAdf = ACdf[ACdf['lambda']>0]
    tRNAdf.index = tRNAdf['anticodon']

    pa =  tRNAdf['pa']
    pb =  E2.loc[list_codons][list_tRNAs]

    pa_n = 1 - pa

    resultMAT = np.zeros(( len(list_codons), len(list_aa) ))
    for cx, codon in enumerate(list_codons):
        current_AA = list(ACdf[ACdf['codon']==codon]['AA'])[0]

        for ax, aa in enumerate(list_aa):
            current_syn = list(tRNAdf[tRNAdf['AA']==aa]['codon'])

            for sx, c in enumerate(current_syn):
                current_AC = list(ACdf[ACdf['codon']==c]['anticodon'])[0]
                if current_AC in list(tRNAdf['anticodon']) and current_AA != aa: 	# only mutations, not WT
                    codon_pb = pb.loc[codon]

                    current_pb = pb.loc[codon, current_AC]  
                    current_pa = tRNAdf.loc[current_AC, 'pa']
                    current_pa_n = 1-current_pa

                    current_ar = tRNAdf.loc[current_AC, 'lambda']
                    other_ar = np.array(tRNAdf['lambda'])[np.array(tRNAdf['anticodon'])!=current_AC]
                    other_pb = np.array(codon_pb)[np.array(tRNAdf['anticodon'])!=current_AC]
                    current_pb_n = current_pb * np.prod( np.power( (1-other_pb), (other_ar)/current_ar )  )

                    resultMAT[cx,ax] +=  (current_pa * current_pb) + (current_pa_n * current_pb_n)
 
    resultMAT = resultMAT/np.mean(np.mean(resultMAT)) * er 	# normalize to expected error rate
    phenoDF = pd.DataFrame(data=resultMAT, index=list_codons, columns=list_aa)


    return phenoDF




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def canonicalCU(fileIN, canonicalDF):
    """
    Computes codon usage across input fasta
    Input:  - input fasta file
            - dataframe with list of canonical transcripts and IDs
    Output: - dictionary of codon usage
    """

    seq = SeqIO.parse(fileIN, "fasta")

    CU_dict = {}
    c = 0
    for s in seq:
        current_name = s.name.split("|")
        current_enst = ""
        current_ensg = ""
        current_cds = ""
        for i in range(len(current_name)):
            if current_name[i][0:4]=="ENST":
                current_enst = current_name[i].split('.')[0]
            if current_name[i][0:4]=="ENSG":
                current_ensg = current_name[i].split('.')[0]
            if current_name[i][0:4]=="CDS:":
                current_cds = np.array(current_name[i].split(':')[1].split('-') , dtype=int)   

        current_seq = s.seq 
        if len(current_enst) > 0 and len(current_ensg) and len(current_cds) > 0: 
            seq_cds = current_seq[(current_cds[0]-1):(current_cds[1]+2)] 	# force as multiple of 3
            seq_cds = seq_cds[:int(np.floor(len(seq_cds)/3.)*3.)]		 	#

        if current_enst in list(canonicalDF['TranscriptID']):

            current_CU = np.zeros(( len(codons_nonstop) ), dtype=int)
            for pos in range( int(np.floor(len(seq_cds)/3.)) ):
                current_codon = seq_cds[(3*pos):((3*pos) + 3)]
                if current_codon in codons_nonstop:
                    current_CU[codons_nonstop.index(current_codon)] += 1        

            CU_dict[current_ensg] = current_CU

    return CU_dict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mut_fractions(peptideDF):
    """
    Extract mutation frequencies from experimental data on phenotypic mutations
    """

    AAsubmat = np.zeros(( 20, 20 ))
    for i in range(len(peptideDF)):
        current_pep = peptideDF.iloc[i]
        current_pep_ref = current_pep['ref_peptide']
        current_pep_alt = current_pep['alt_peptide']
        current_int_ref = current_pep['ref_int']
        current_int_alt = current_pep['alt_int']
        current_codon = current_pep['codon']

        if len(current_pep_ref) == len(current_pep_alt):
            for j in range(len(current_pep_ref)):
                current_aa_ref = current_pep_ref[j]
                current_aa_alt = current_pep_alt[j]

                if current_aa_ref == current_aa_alt:
                    AAsubmat[AA.index(current_aa_ref), AA.index(current_aa_ref)] += 1 #(current_int_ref + current_int_alt)
                else:
                    AAsubmat[AA.index(current_aa_ref), AA.index(current_aa_alt)] += 1 #current_int_alt
        else:
            print(current_pep_ref)

    for i in range(np.shape(AAsubmat)[0]):
        AAsubmat[i,:] /= AAsubmat[i,i] #np.sum(AAsubmat[i,:])
        AAsubmat[i,i] = np.nan

    return AAsubmat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def codon_mut_fractions(peptideDF, codonDict):
    """
    Extract frequencies of phenotypic mutations from codon to AA
    """

    submat = np.zeros(( 61, 21 ))

    for i in range(len(peptideDF)):		# 
        current_pep = peptideDF.iloc[i]
        current_pep_ref = current_pep['ref_peptide']
        current_pep_alt = current_pep['alt_peptide']
        current_int_ref = current_pep['ref_int']
        current_int_alt = current_pep['alt_int']

        current_codon = current_pep['codon']
        current_refpos = int(current_pep['misProt'])
        current_ensembl = current_pep['ensembl']

        if current_ensembl in list(ensemblID['Protein stable ID']):
            current_enst = list(ensemblID[ensemblID['Protein stable ID']==current_ensembl]['Transcript stable ID'])[0]
        else:
            current_enst = "NA"

        if current_enst != "NA" and current_enst in list(codonDict.keys()):
            current_seq = codonDict[current_enst]
            refpos_codon = current_seq[current_refpos]

            print(current_codon, refpos_codon)

            mutpos = -1
            for j in range(len(current_pep_ref)):
                current_aa_ref = current_pep_ref[j]
                current_aa_alt = current_pep_alt[j]
                
                if current_aa_ref != current_aa_alt:
                    mutpos = j
                    break

            L = len(current_pep_ref)
            pep_start = current_refpos - mutpos

            for j in range(len(current_pep_ref)):
                current_aa_ref = current_pep_ref[j]
                current_aa_alt = current_pep_alt[j]

                if len(current_seq) > pep_start + j:
                    current_pep_codon = current_seq[(pep_start + j)]
                    current_transl = gencode[current_pep_codon]
                else:
                    current_transl = None

                print(j, current_aa_ref, current_aa_alt, current_pep_codon, current_transl)
                    
                if current_aa_ref == current_transl:

                    if current_aa_ref == current_aa_alt:
                        submat[codons_nonstop.index(current_pep_codon), 20] += 1 #last columns is base counts, first 20 is mutations
                    else:
                        submat[codons_nonstop.index(current_pep_codon), AA.index(current_aa_alt)] += 1 #current_int_alt
        else:
            print(current_pep_ref)

    submat2 = submat[:,0:20]
    for i in range(np.shape(submat)[0]):
        submat2[i,:] /= submat[i,20] 
    subDF = pd.DataFrame(data=submat2, index=codons_nonstop, columns=AA)
    subDF.to_csv("tmp.sub2", header=True, index=True, sep='\t')

    return submat2



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
def phenomutations_from_peptides():
    """
    # Project peptides with phenotypic mutations onto coding seqs to derive codon -> AA freqs
    """
    peps = pd.read_csv("../data/human/pilpel_peptides.txt", header=0, index_col=None, sep='\t')
    enstensp = pd.read_csv("../data/human/list.ensp", header=0, index_col=None, sep=' ')

    pepMAT = np.zeros(( len(codons_nonstop), 20 ), dtype=int)
    for i in range(len(peps)):
        current_pep = peps.iloc[i]
        current_ref = np.array(list(current_pep['ref_peptide']))
        current_alt = np.array(list(current_pep['alt_peptide']))
        current_mismatch = np.where(current_ref != current_alt)[0][0]
        current_misProt = current_pep['misProt']
        current_ensp = current_pep['ensembl']
        current_length = len(current_ref)
        current_codon = current_pep['codon']
        current_dest = list(current_pep['destination'])[0]
        current_orig = list(current_pep['origin'])[0]

        if current_ensp in list(enstensp['ENSP']):
            current_enst = list(enstensp[enstensp['ENSP']==current_ensp]['ENST'])[0]
            if current_enst in list(dict_seq.keys()):
                current_seq = dict_seq[current_enst]
                list_codons = []
                for j in range( int(np.floor(len(current_seq)/3.)) ):
                    list_codons.append( "".join(list(current_seq[(j*3):((j*3)+3)]))  )
                current_seq2 = np.array(list_codons)

                current_cpep = current_seq2[ (current_misProt - current_mismatch) : (current_misProt + current_length - current_mismatch)]

                for j in range(len(current_cpep)):
                    if current_cpep[j] in codons_nonstop:

                        if j != current_mismatch:		# no mutation
                            current_ccodon = current_cpep[j]
                            current_caa = gencode[current_ccodon]
                            pepMAT[codons_nonstop.index(current_ccodon), AA.index(current_caa)] += 1
                        else: 							# mutation
                            pepMAT[codons_nonstop.index(current_codon), AA.index(current_dest)] += 1

    pepmutDF = pd.DataFrame(data=pepMAT, index=codons_nonstop, columns=AA)

    for i in range(len(pepmutDF)):
        current_df = pepmutDF.iloc[i]
        current_wt = gencode[pepmutDF.index[i]]
        current_ref = current_df[current_wt]
        pepmutDF.iloc[i] /= current_ref
       
    dd = np.array(pepmutDF)
    m = np.percentile(dd[dd!=1], 25)
    dd[dd==0] = m
    dd[dd==1] = 0
    dd = dd/(np.mean(dd[dd>0])) * 0.0025 	# normalize to exp error rate
    pepmutDF = pd.DataFrame(data=dd, index=pepmutDF.index, columns=pepmutDF.columns)
    pepmutDF.to_csv("../data/processed/peptides_mutations.txt", header=True, index=True, sep='\t')

    return pepmutDF





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


    dict_seq = cdsseq_dict("../data/gencode.v46.pc_transcripts.fa")


    # compute transcript lengths
    ##lenDF = pd.DataFrame(columns=['enst', 'length'])
    ##for i in list(dict_seq.keys()):
    ##    lenDF.loc[len(lenDF)] = (i, len(dict_seq[i]))
    ##lenDF.to_csv("../data/processed/enst_length.txt", header=True, index=False, sep='\t')
    
    # compute background AA, NT, and codon counts in dataset
    codon_count = np.zeros(( len(codons_nonstop) ), dtype=int)
    AA_count = np.zeros(( len(AA) ), dtype=int)
    nt_count = np.zeros(( len(nt) ), dtype=int)
    for i in list(dict_seq.keys()):
        current_seq = dict_seq[i]
        for j in range( int(np.floor(len(current_seq)/3.)) ):
            current_codon = "".join(list(current_seq[(j*3):((j*3)+3)])) 
            if current_codon in codons_nonstop:
                codon_count[codons_nonstop.index(current_codon)] += 1
                current_AA = gencode[current_codon]
                AA_count[AA.index(current_AA)] += 1
                nt_count[nt.index(current_codon[0])] += 1
                nt_count[nt.index(current_codon[1])] += 1
                nt_count[nt.index(current_codon[2])] += 1
    cdDF = pd.DataFrame({'codon': list(codons_nonstop), 'count': list(codon_count)})
    cdDF.to_csv("../data/processed/composition_codons.txt", header=True, index=False, sep='\t')
    aaDF = pd.DataFrame({'AA': list(AA), 'count': list(AA_count)})
    aaDF.to_csv("../data/processed/composition_AA.txt", header=True, index=False, sep='\t')
    ntDF = pd.DataFrame({'nt': list(nt), 'count': list(nt_count)})
    ntDF.to_csv("../data/processed/composition_nt.txt", header=True, index=False, sep='\t')


    # load tRNA and codon-anticodon data
    anticodonDF = pd.read_csv("../data/tRNA/anticodon.txt", header=0, index_col=None, sep='\t')
    trnaDF = pd.read_csv('../data/tRNA/human_tRNA.txt', delimiter=" ", header=0)



    # Comparison of phenotypic mutations computed with yeast ribosome energy parameters to human freqs
    E = pickle.load(open('../data/tRNA/scer2_basic_young_E.pkl', 'rb')) 	# from Pechmann JRSI (2025)
    P = np.loadtxt('../data/tRNA/scer2_basic_young_P.txt')
    CC = np.loadtxt('../data/tRNA/scer2_basic_young_CC.txt')

    CC_order = np.argsort(CC)
    CC_order = CC_order[::-1]

    pheno_top10 = np.zeros(( 61, 20, 10 )) * np.nan
    pheno_top_cor = np.zeros(( 61, 20, 10 )) * np.nan
    for i in range(10):						# representative average of the best 10 fits

        topE = E[:,:,CC_order[i]]
        topP = P[:,CC_order[i]]
        e2 = compute_pb(topP, topE)
        ACdf = arrival_probs(anticodonDF, trnaDF)
        PHENO = binding_probs(ACdf, e2)

        pheno_top10[:,:,i] = np.array(PHENO)
   
    phenoDF = pd.DataFrame(data=np.nanmean(pheno_top10, axis=2), index=PHENO.index, columns=PHENO.columns)
    phenoDF.to_csv("../data/processed/pheno_human_yeastE.txt", header=True, index=True, sep='\t')

    # Mutation freqs and peptides with phenotypic mutations in human from doi:10.1101/.2025.02.17.638604
    pilpel = pd.read_csv("../data/human/pilpel.txt", header=0, index_col=0, sep='\t')
    resDF = pd.DataFrame(columns=['p1', 'p2'])
    for aa in list(set(anticodonDF['AA'])):
        current_codons = list(anticodonDF[anticodonDF['AA']==aa]['codon'])
        current_pheno = phenoDF.loc[current_codons]
        if aa in list(pilpel.columns):
            current_pilpel = pilpel.loc[aa]
            for aa2 in list(set(anticodonDF['AA'])):
                pred = np.mean(np.array(list(current_pheno[aa2])))
                if aa2 in list(pilpel.columns):
                    plpl = pilpel.loc[aa][aa2]
                    r = np.array([pred, plpl])
                    if np.all(np.isnan(r)==False):
                        resDF.loc[len(resDF)] = (pred, plpl)
    resDF.to_csv("../data/processed/pheno_pred_peptides.txt", header=True, index=False, sep='\t')


    # Phenotypic mutation frequencies codon -> AA from experimental data
    pepmutDF = phenomutations_from_peptides() 		# once computed, just read file instead
    #pepmutDF = pd.read_csv("../data/processed/peptides_mutations.txt", header=0, index_col=0, sep='\t')
 