import sys
from Bio import SeqIO
from Bio.PDB import *
from bioservices import UniProt

import io
import glob
import os
import biolib
import pandas as pd
import re
import csv
import mdtraj as md
from contact_map import ContactFrequency
import numpy as np

deeptmhmm = biolib.load('DTU/DeepTMHMM')
working_directory = r'/mnt/c/Users/Jonah/OneDrive/Documents/AlphaFold_Project/contact_map'
df = pd.read_csv('results5.csv',header=0)
print(df)
os.chdir(working_directory)
cutoff_dist = 2.0 #in nm

def find_outside(protein_pdb, subcellular_location):
    """This function detect if there is contact in an inside/outside region.

    Args:
        protein_pdb (trajectorymd.Trajectory): pdb file returned from the function mdtraj.load_pdb
        subcellular_location (pandas.Series): 

    Returns:
        _type_: _description_
    """
    indices_outside = [i for i, elem in enumerate([*subcellular_location]) if 'O' in elem]
    indices_inside = [i for i, elem in enumerate([*subcellular_location]) if 'I' in elem]
    frame_contacts = ContactFrequency(protein_pdb,cutoff=cutoff_dist)
    contact_df = frame_contacts.residue_contacts.df
    len_seq = len(contact_df)
    len_topo = len(subcellular_location)
    print("length df :", len_seq)
    print("length string :", len_topo)

    steps_outside = subcellular_location.count('O')
    steps_inside = subcellular_location.count('I')
    print('length of outside part : ', steps_outside )
    print('length of inside part : ', steps_inside )
    return_df = contact_df[indices_outside]
    return_df = return_df.iloc[indices_inside]
    if 1.0 in return_df.values:
        result_pred = 1
        print("Wrong prediction")
    else:
        result_pred = 0
        print('Right Prediction')
    print('\n')

    contact_df = contact_df.rename_axis(index='idx', columns='cols')
    df_contact_points = contact_df.stack().reset_index(name='value').query('value == 1')
    df_contact_points['contact_tuples'] = list(zip(df_contact_points.idx, df_contact_points.cols))
    list_tuples = df_contact_points['contact_tuples'].tolist()
    return return_df,indices_inside,result_pred,steps_outside,steps_inside, len_seq, len_topo,list_tuples

def get_bfactors(PDBFile):
    p=PDBParser()
    structure=p.get_structure('name', PDBFile)
    bfactors =[]
    for model in structure:   # X-Ray generally only have 1 model, while more in NMR
        for chain in model: 
            for residue in chain:
                for atom in residue:
                    bfactors.append(atom.get_bfactor())
    AF_pred_mean = np.mean(bfactors)
    return 
    
def uniprot_info(protein_code):
    u = UniProt(verbose=False)
    gff = u.search('accession:%s'%protein_code, frmt='gff')
    fasta = u.search('accession:%s'%protein_code, frmt='fasta')
    searchlines = io.StringIO(gff)
    gff_df = pd.read_csv(
        searchlines, 
        sep="\t", 
        header=None,
        comment="#", 
        names=("seqname", "source", "feature", "start",
            "end", "score", "strand", "frame", "attribute", 'other'),
    )

    gff_df['Note'] = gff_df['attribute'].str.extract(r'Note\=(.*?)\;')
    gff_df['Ontology_term'] = gff_df['attribute'].str.extract(r'Ontology_term\=(.*?)\;')
    gff_df['Evidence'] = gff_df['attribute'].str.extract(r'evidence\=(.*?)\;')
    gff_df.drop(columns=['attribute', 'other'],inplace=True)
    return gff_df

def main():   
    for PDBFile in glob.iglob('./afdb_homo_sapiens/**.pdb',recursive=True):
        try:
            regexPDB = re.search('AF-(.*)-F1', str(PDBFile))        #get the code on file name to be able to look in uniprot
            PDBcode = regexPDB.group(1)
            gff_df = uniprot_info(PDBcode)
            if PDBcode in df['SMR'].tolist():
                print('%s is already in the csv'%PDBcode)
                pass
            elif gff_df['feature'].str.contains('Transmembrane').any() == False:
                print('%s has no transmembrane'%PDBcode)
                pass
            elif PDBcode is None:
                pass
            else:
                AF_pred_mean = get_bfactors(PDBFile)
                gff_df = uniprot_info(PDBcode)
                md_pdb = md.load_pdb(PDBFile)
                print("md_traj pdb loaded")
            
                with open(PDBFile, 'r') as pdb_file:
                    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                        to_fasta = '>' + record.id + '\n' + record.seq
                        with open("query.fasta", "w") as text_file:
                            text_file.write(str(to_fasta))
                        deeptmhmm_res = deeptmhmm.cli(args='--fasta query.fasta')
                        deeptmhmm_res.save_files('./outputs_DeepTMHMM/output_DeepTMHMM_%s'%PDBcode)
                        df_topology = pd.read_csv('./outputs_DeepTMHMM/output_DeepTMHMM_%s/predicted_topologies.3line'%PDBcode,delim_whitespace=True,header=None)
                        topology=df_topology.loc[2, 0]
                        print('\n')
                        _,_,result_pred,steps_outside,steps_inside, len_seq, len_topo,list_tuples = find_outside(md_pdb,topology)
                        if 'M' in topology:
                            is_membrane = 1
                        else:
                            is_membrane = 0
                        to_append = [PDBcode,record.seq,topology,result_pred,steps_outside,steps_inside, len_seq, len_topo,is_membrane,AF_pred_mean]
                        print(record.seq)
                        print(str(record.seq))
                        print(topology)
                        print('\n')
                        with open('results5.csv','a') as fd:
                            wr = csv.writer(fd, dialect='excel')
                            wr.writerow(to_append)
        except AttributeError:
            print('attribute error')
            pass        


if __name__ == "__main__":
    main()
    


#restart at O00462