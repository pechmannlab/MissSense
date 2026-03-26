import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import time
import multiprocessing as mp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def components_entropy(Gmatrix):
    """
    Calculation of network entropy
    BFS to compute size of connected components
    """

    N = np.shape(Gmatrix)[0]
    comp_array = np.zeros(( N ), dtype=int)
    visited = np.zeros(( N ), dtype=bool)
    component = 0
    read_pointer = 0
    queue = []

    H = 1

    seeds = np.arange(N)
    while len(seeds) > 0:
        i = int( np.random.choice(seeds, 1))
        component += 1
        read_pointer = 0
        queue = []
        queue.append(i)
        visited[i] = True

        while read_pointer < len(queue):
            current_node = queue[read_pointer]
            comp_array[current_node] = component
            read_pointer += 1

            current_neighbours = np.where( Gmatrix[current_node,:] == 1)[0]
            if len(current_neighbours) > 0:
                current_neighbours = current_neighbours[ visited[current_neighbours] == False ] 
                queue += list(current_neighbours)
                visited[current_neighbours] = True
            else:
                break
                
        seeds = np.where( visited == False )[0]
        
    C = np.array( [np.sum(comp_array == c) for c in list(set(list(comp_array)))] ) 
    p = C / float(N) 
    p = p[p > 0]			# just in case, tho log of zero counts not possible really here
    H = -1 * ( np.sum( p * np.log(p) ) )

    return H



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def interactome_resilience(Gmat):
    """
    Implementation of interactome resilience algorithm
    by Zitnik & Leskovec (PNAS 2019)
    """

    N = np.shape(Gmat)[0]
    G2 = np.copy(Gmat)

    result = np.zeros(( 100 ))
    current_selection = np.ones(( N ), dtype=bool)

    for i in np.arange(100):
        current_fragment = np.min([ int(N / 100. ), np.shape(G2)[0] ])
        current_failure = np.random.choice(np.where(current_selection==True)[0], int(current_fragment) )
        current_selection[current_failure] = False
        G2 = np.copy(Gmat)
        G2[current_selection,:] = G2[:,current_selection] = 0

        current_entro = components_entropy(G2)
        norm_entro = (1/np.log(N)) * current_entro
        result[i] = norm_entro

    resilience = np.around( 1 - ( np.sum( 0.01 * result ) ), 3 ) 

    return resilience

  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_network(fileIN):
    """
    Load graph from edgelist
    Filter nodes to the canonical isoforms and genes in the current dataset
    Output: binary numpy array and nodelist
    """

    cid = pd.read_csv("../data/human/canonical_id.txt", header=0, index_col=False, sep='\t')
    arv = pd.read_csv("../data/processed/all_reg_vars.txt", header=0, index_col=None, sep='\t', na_values='NA')

    G = nx.Graph()
    c = 0
    nodelist = []
    for line in open(fileIN, 'r'):	
        if c > 0:                   		
            current_line = line.split('\t')
            A = current_line[0].split(":")[-1]
            B = current_line[1].split(":")[-1]
            ids = sorted([A, B])
            id1 = ids[0]
            id2 = ids[1]
            if id1 in list(cid['UniProtID']) and id2 in list(cid['UniProtID']):
                n1 = list(cid[cid['UniProtID']==id1]['Name'])[0]
                n2 = list(cid[cid['UniProtID']==id2]['Name'])[0]
                if n1 in list(arv['Gene']) and n2 in list(arv['Gene']):

                    if n1 not in nodelist:
                        nodelist.append(str(n1))
                    if n2 not in nodelist:
                        nodelist.append(str(n2))

                    G.add_edge(n1, n2)
        
        c += 1

        if c % 100000 == 0:
            print(c, "processed")

    nl_sorted = sorted(list(set(nodelist)))
    GG = nx.to_numpy_array( G, nodelist=nl_sorted)
    np.savetxt("../data/ppi/nx_ebippi.txt", GG, fmt='%i')

    nlDF = pd.DataFrame({'nodelist': nl_sorted})
    nlDF.to_csv("../data/ppi/nx_nodelist.txt", header=True, index=False, sep='\t')

    return GG



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pheno_network(Gin, seed, alpha=1):
    """
    Simulation of phenotypic mutations in the protein interaction network:
    - For each node (protein), the fraction of damaged protein interactions
      is computed from the protein mutation risk and a arbitrary scaling factor
      chosen such that the fractions of damaged interactions distributes between 
      0 and 1 in a way that there is some spread, i.e. differences betweenn most proteins
    - a random fraction of protein interactions for that protein is set to 0. E.g. if 
      a protein has 20 interactions, and based on the phenotypic mutation risk it is 
      modelled that 25% of that protein's interactions are lost due to phenotypic mutations, 
      then 5 out of the 20 interactions of this protein are randomly set to 0. 
    Output: new network of same size with removed interactions to simulate phenotypic mutations
    """

    np.random.seed(seed)
    
    G2 = np.zeros( np.shape(Gin), dtype=int)
    arv = pd.read_csv("../data/processed/all_reg_vars.txt", header=0, index_col=None, sep='\t', na_values='NA')
    for ix, i in enumerate(list(nodes['nodelist'])):
        if i in list(arv['Gene']):
            current_p = list(arv[arv['Gene']==i]['pheno_full'])[0]
            current_pct = current_p / 3 #np.nanmax(arv['pheno_full'])
            current_n  = np.where(Gin[ix,:]==1)[0]
            if np.isnan(current_pct) == False:
                current_drop = np.min([int(np.floor( current_pct * alpha * len(current_n) )), len(current_n)])
                current_idx = np.random.choice(current_n, size=current_drop, replace=False)
                G2[ix, current_idx] = G2[current_idx, ix] = 1

    Gout = np.copy(Gin)
    Gout[G2==1] = 0

    return Gout



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def geno_network(Gin, genodf):
    """
    Simulation of genotypic damaging mutations in the protein network
    - The total risk of a damaging mutations joinly for all genes is 
      evaluated against a random number r1
    - If a damaging mutation is to occur, a gene is chosen based on its
      mutation risk by random number r2. All interactions of that node
      are set to 0. 
    Output: a new network of the same size  
    """

    Gout = np.copy(Gin)
  
    z = np.sum(Gout == 0, axis=1) == np.shape(Gout)[1]
    genodf.loc[z,'geno_full'] = 0
    genodf['norm'] = genodf['geno_full']/np.nansum(genodf['geno_full'])
    total = np.sum(genodf['geno_full'])

    r1 = np.random.random_sample()
    if r1 < total:

        r2 = np.random.random_sample()
        cs = np.cumsum(np.array(genodf['norm']))
        sel = np.where( r2 < cs )[0][0]      
        current_gene = genodf.loc[sel]['gene']  #
        idx = list(genodf['gene']).index(current_gene)

        Gout[idx,:] = Gout[:,idx] = 0

    return Gout



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pheno_parallel(packaged_input):
    """
    Wrapper to parallelize computation of phenotypic mutations and network resilience
    """
    G = packaged_input[0]
    seed = packaged_input[1]
    alpha = packaged_input[2]
    G_pheno = pheno_network(G, seed, alpha)
    r  = interactome_resilience(G_pheno)
    return r




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def forward_evolsims():
    """
    Forward simulations of the effect of genotypic and phenotypic mutations on network resilience
    - 3 replicas of genotypic mutations run in parallel
    - cell division times of 2d and 5d are tested
    - for each generation and replica, 1000 runs of phenotypic mutations are run 
    - phenotypic mutations increase with age so that at age 80 there is a 50% increase
    """

    resiliDF = pd.DataFrame(columns=['divs', 'rep', 'n_geno', 'R_geno'])
    resiphDF = pd.DataFrame(columns=['time', 'rep', 'tdiv', 'alpha', 'R_combi'])

    pool = mp.Pool(processes=10)
    output = pool.map(interactome_resilience, [G for iteration in range(10)])
    output = np.array(output)
    pool.close()
    R_geno = np.around(np.mean(output), 3)
    resiliDF.loc[len(resiliDF)] = (0, 0, 0, R_geno)
    resiliDF.loc[len(resiliDF)] = (0, 1, 0, R_geno)
    resiliDF.loc[len(resiliDF)] = (0, 2, 0, R_geno)


    pool = mp.Pool(processes=10)
    output = pool.map(pheno_parallel, [(G, 0, 1) for iteration in range(10)])
    output = np.array(output)
    pool.close()
    R_pheno = np.mean(output)
    resiphDF.loc[len(resiphDF)] = (0, 0, 'd2', 1, np.around(R_pheno, 3))
    resiphDF.loc[len(resiphDF)] = (0, 1, 'd2', 1, np.around(R_pheno, 3))
    resiphDF.loc[len(resiphDF)] = (0, 2, 'd2', 1, np.around(R_pheno, 3))
    resiphDF.loc[len(resiphDF)] = (0, 0, 'd5', 1, np.around(R_pheno, 3))
    resiphDF.loc[len(resiphDF)] = (0, 1, 'd5', 1, np.around(R_pheno, 3))
    resiphDF.loc[len(resiphDF)] = (0, 2, 'd5', 1, np.around(R_pheno, 3))


    R  = interactome_resilience(G)
    print(0, R)

    G1 = np.copy(G)
    G2 = np.copy(G)
    G3 = np.copy(G)
    Gbag = [G1, G2, G3]                 # bag of graphs to iterate through

    tp_d2 = np.arange(11) * 10 * 365/2
    tp_d5 = np.arange(11) * 10 * 365/5
    m = np.zeros(( 3 ), dtype=int) 
    cell_divisions = int(0)
    max_divisions = 18260                   
    while cell_divisions < max_divisions:
        cell_divisions += 1

        for gx, g in enumerate(Gbag):
            gg = geno_network(g, genoDF)

            if np.any(g != gg):
                m[gx] += 1
                Gbag[gx] = np.copy(gg)

                #R_geno  = interactome_resilience(Gbag[gx])
                #print(cell_divisions, gx, m[gx], R_geno)

                pool = mp.Pool(processes=3)
                output = pool.map(interactome_resilience, [Gbag[gx] for iteration in range(3)])
                output = np.array(output)
                pool.close()

                R_geno = np.around(np.mean(output), 3)
                resiliDF.loc[len(resiliDF)] = (cell_divisions, gx, m[gx], R_geno)


            if np.any(tp_d2 == cell_divisions):
                current_t = np.where(tp_d2 == cell_divisions)[0][0] * 10  
                current_alpha = 1 + 0.5/np.exp(np.power(80,0.5)) * np.exp(np.power(current_t,0.5))

                pool = mp.Pool(processes=10)
                output = pool.map(pheno_parallel, [(Gbag[gx], iteration, current_alpha) for iteration in range(10)])
                output = np.array(output)
                pool.close()

                R_pheno = np.mean(output)
                resiphDF.loc[len(resiphDF)] = (current_t, gx, 'd2', np.around(current_alpha, 5), np.around(R_pheno, 3))


            if np.any(tp_d5 == cell_divisions):
                current_t = np.where(tp_d5 == cell_divisions)[0][0] * 10  
                current_alpha = 1 + 0.5/np.exp(np.power(80,0.5)) * np.exp(np.power(current_t,0.5))

                pool = mp.Pool(processes=10)
                output = pool.map(pheno_parallel, [(Gbag[gx], iteration, current_alpha) for iteration in range(10)])
                output = np.array(output)
                pool.close()

                R_pheno = np.mean(output)
                resiphDF.loc[len(resiphDF)] = (current_t, gx, 'd5', np.around(current_alpha, 5), np.around(R_pheno, 3))


    resiliDF.to_csv("../data/processed/resilience_geno.txt", header=True, index=False, sep='\t')
    resiphDF.to_csv("../data/processed/resilience_combi.txt", header=True, index=False, sep='\t')

    resiphenoDF = pd.DataFrame(columns=['time', 'alpha', 'R_pheno'])
    for i in (np.arange(11) * 10):
        current_alpha = 1 + 0.5/np.exp(np.power(80,0.5)) * np.exp(np.power(i,0.5))

        pool = mp.Pool(processes=10)
        output = pool.map(pheno_parallel, [(G, iteration, current_alpha) for iteration in range(10)])
        output = np.array(output)
        pool.close()

        R_pheno = np.mean(output)
        resiphenoDF.loc[len(resiphenoDF)] = (i, np.around(current_alpha, 5), np.around(R_pheno, 3))


    resiphenoDF.to_csv("../data/processed/resilience_pheno.txt", header=True, index=False, sep='\t')
    print(resiphenoDF)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def network_randvsgeno():
    """
    Network resilience upon 
    - random node loss
    - node loss based on mutation risk
    """
   
    resultDF = pd.DataFrame(columns=['drop', 'class', 'R'])
    for drop in [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000 ]:

        # random node loss, average of 10 replicas
        rr = np.zeros(( 10 ))
        for i in range(10):
            Grand = np.copy(G)
            idx = np.random.choice(np.arange(np.shape(Grand)[0]), drop, replace=False) 
            Grand[idx,:] = Grand[:,idx] = 0
            rr[i]  = interactome_resilience(Grand)
        resultDF.loc[len(resultDF)] = (drop, 'rand', np.around(np.mean(rr), 5) )


        # node loss based on mutation risk, average of 10 replicas
        rr = np.zeros(( 10 ))
        for i in range(10):
            Grand = np.copy(G)
            r = np.random.random_sample(drop+50)            # in case of double choice
            cs = np.cumsum(np.array(genoDF['norm']))
            list_idx = []
            for j in range(drop+50):
                sel = np.where( r[j] < cs )[0][0]      
                list_idx.append(sel)

            idx = np.array(list(set(list_idx))[0:drop])     # sel only req
            Grand[idx,:] = Grand[:,idx] = 0
            rr[i]  = interactome_resilience(Grand)
        resultDF.loc[len(resultDF)] = (drop, 'prot', np.around(np.mean(rr), 5) )

    resultDF.to_csv("../data/processed/resilience_rand_geno.txt", header=True, index=False, sep='\t')
    print(resultDF)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hpo_enrichment():
    """
    Evaluate the mutation risk of genes annotated in the Human Phenotype Ontology (HPO)
    - extracts all genes belonging to a phenotype or disease ontology
    - computes the avarage of the gene mutation risks
    - generates a randomized control of average of the same sizes of gene sets (10000 draws each)
    """


    # disease annotations
    disease_genes = []
    diseaseDF = pd.DataFrame(columns=['id', 'disease', 'score', 'n_score', 'genes', 'n_genes'])
    for i in list(set(g2d['disease_id'])):
        current_hpo = g2d[g2d['disease_id']==i]
        current_genes = list(set(current_hpo['gene_symbol']))
        if i in list(set(hpo['database_id'])):
            disease_genes += current_genes
            current_disease = list(hpo[hpo['database_id']==i]['disease_name'])[0]
            current_disease = str(current_disease).replace(" ", "_")

            score = np.zeros(( len(current_genes) )) * np.nan
            if len(current_genes) > 0 and "#" not in current_disease:
                for jx, j in enumerate(current_genes):
                    if j in list(arv['Gene']):
                        score[jx] = float(list(arv[arv['Gene']==j]['geno_full'])[0])

            score = score[np.isnan(score)==False]
            if len(score) > 0:
                ms = np.mean(np.array(score))
                diseaseDF.loc[len(diseaseDF)] = (i, current_disease, ms, len(score), ",".join(current_genes), len(current_genes))

    diseaseDF.to_csv("../data/processed/hpo_disease.txt", header=True, index=False, sep='\t', na_rep="NA")
    diseasegenesDF = pd.DataFrame({'gene': list(set(disease_genes))})
    diseasegenesDF.to_csv("../data/processed/hpo_diseasegenes.txt", header=True, index=False)


    # phenotype annotations
    pheno_genes = []
    p2g = pd.read_csv("../data/tmp/missense/phenotype_to_genes.txt", header=0, index_col=None, sep='\t')
    phenotypeDF = pd.DataFrame(columns=['id', 'phenotype', 'score', 'n_score', 'genes', 'n_genes'])
    for i in list(set(p2g['hpo_id'])):
        current_hpo = p2g[p2g['hpo_id']==i]
        current_genes = list(set(current_hpo['gene_symbol']))
        current_hponame = str(list(current_hpo['hpo_name'])[0]).replace(' ', '_')
        if len(current_genes) > 100:
            pheno_genes += current_genes
            score = np.zeros(( len(current_genes) )) * np.nan
            for jx, j in enumerate(current_genes):
                if j in list(arv['Gene']):
                    score[jx] = float(list(arv[arv['Gene']==j]['geno_full'])[0])

            score = score[np.isnan(score)==False]
            if len(score) > 0:
                ms = np.mean(np.array(score))
            else:
                ms = np.nan

            phenotypeDF.loc[len(phenotypeDF)] = (i, current_hponame, ms, len(score), ",".join(current_genes), len(current_genes))

    phenotypeDF.to_csv("../data/processed/hpo_phenotype.txt", header=True, index=False, sep='\t')
    phenogenesDF = pd.DataFrame({'gene': list(set(pheno_genes))})
    phenogenesDF.to_csv("../data/processed/hpo_phenogenes.txt", header=True, index=False)



    # randomized control
    phenotypeDF = pd.read_csv("../data/processed/hpo_phenotype.txt", header=0, index_col=None, sep='\t')
    randDF = pd.DataFrame(columns=['n', 'median', 'q95', 'q99'])
    gf = np.array(arv['geno_full'])
    gf = gf[np.isnan(gf)==False]    # no np.nan here
    for i in range(1, 300):
        rando = np.zeros(( 10000 )) * np.nan
        for r in range(10000):
            np.random.seed(r)
            current_rand = np.random.choice(gf, i, replace=False)
            rando[r] = np.nanmean(current_rand)

        current_median = np.nanmedian(rando)
        current_q95 = np.percentile(rando[np.isnan(rando)==False], 95)
        current_q99 = np.percentile(rando[np.isnan(rando)==False], 99)

        randDF.loc[len(randDF)] = (i, current_median, current_q95, current_q99)

    randDF.to_csv("../data/processed/hpo_rand.txt", header=True, index=False, sep='\t')





if __name__ == '__main__':


    ## load EBI consensus human PPI network
    #G = load_network("../data/ppi/human.txt")
    G = np.loadtxt("../data/ppi/nx_ebippi.txt")
    nodes = pd.read_csv("../data/ppi/nx_nodelist.txt", header=0, index_col=None)


    # dataframe with mutation risk matching the network nodelist
    arv = pd.read_csv("../data/processed/all_reg_vars.txt", header=0, index_col=None, sep='\t', na_values='NA')
    genoDF = pd.DataFrame(columns=['gene', 'geno_full'])
    for ix, i in enumerate(list(nodes['nodelist'])):
        if i in list(arv['Gene']):
            current_p = float(list(arv[arv['Gene']==i]['geno_full'])[0]) * 10e-9
        else:
            current_p = 0 #np.nan
        genoDF.loc[len(genoDF)] = (i, current_p)

    genoDF['norm'] = genoDF['geno_full']/np.nansum(genoDF['geno_full'])
    genoDF.to_csv("tmp.geno", header=True, index=False, sep='\t')



    # Forward simulation of the effect of genotypic and phenotypic mutations on network resilience
    forward_evolsims()

    # Network resilience for random vs. mutation risk driven node loss
    network_randvsgeno()



    # extract bottleneck nodes through betweenness-centrality
    gg = nx.from_numpy_array(G)
    bc = nx.betweenness_centrality(gg)
    bcDF = pd.DataFrame(columns=['gene', 'betweenness', 'degree'])
    for ix, i in enumerate(list(nodes['nodelist'])):
        print(i, bc[ix])
        bcDF.loc[len(bcDF)] = (i, bc[ix], np.sum(G[ix,:]))
    bcDF.to_csv("../data/processed/betweenness_centrality.txt", header=True, index=False, sep='\t')


    # averages over replicas in human protein half-life data
    hl = pd.read_csv("../data/human/human_proteinHL.txt", header=0, index_col=None, sep='\t')
    list_columns = list(hl.columns)
    hl['Bcells'] = np.nanmean( np.array(hl[list_columns[1:3]]), axis=1)
    hl['NKcells'] = np.nanmean( np.array(hl[list_columns[3:5]]), axis=1)
    hl['Hepatocytes'] = np.nanmean( np.array(hl[list_columns[5:7]]), axis=1)
    hl['Monocytes'] = np.nanmean( np.array(hl[list_columns[7:9]]), axis=1)
    hl2 = hl[['gene_name', 'Bcells', 'NKcells', 'Hepatocytes', 'Monocytes']]
    hl2.to_csv("../data/processed/human_proteinHL.txt", header=True, index=False, sep='\t', na_rep="NA")



    # load HPO phenotype and disease annotations
    arv = pd.read_csv("../data/processed/all_reg_vars.txt", header=0, index_col=None, sep='\t', na_values='NA')
    g2d = pd.read_csv("../data/human/genes_to_disease.txt", header=0, index_col=None, sep='\t')
    hpo = pd.read_csv("../data/human/phenotype.hpoa", header=0, skiprows=4, index_col=None, sep='\t', low_memory=False)

    hpo_enrichment()
