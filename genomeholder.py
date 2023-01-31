import os
#from Bio.Blast.Applications import NcbiblastxCommandline
#from Bio import SeqIO
from os.path import join
#import sqlite3
#from Bio.Blast import NCBIXML
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature, FeatureLocation
#import pickle
from classholder import IsPair
from helperobject import helperobject
#import shutil
approved_gb_file_endings= ["gbk", "gb","gbff"]


class MyContig:
    def __init__(self,home_directory,filepath,approved_gb_file_endings):
        self.filepath=filepath
        self.recordobject=None


class MyGenome:
    def __init__(self,home_directory,approved_gb_file_endings):
        self.dic={}
        self.dic["files"]=[]
        self.directory=""
        self.files=[]
        self.home_directory=home_directory
        self.approved_gb_file_endings=approved_gb_file_endings

    def initialize(self):
        print ("created mygenome")
        if os.path.isdir(self.home_directory):
            self.dic["home_directory"]=self.home_directory
            self.dic["number_of_files"]=0
            self.dic["mycontig_objects"]=None
            for gbfile in os.listdir(self.home_directory):
                if gbfile[0]==".":
                    continue
                if gbfile.split(".")[-1] in self.approved_gb_file_endings and not os.path.isdir(os.path.join(self.home_directory,gbfile)):
                    self.dic["files"].append(gbfile)
                    self.dic["number_of_files"]+=1
                    self.files.append(gbfile)
        return self.dic["number_of_files"]


class Genomelist:
    def __init__(self,approved_gb_file_endings):
        self.approved_gb_file_endings = approved_gb_file_endings
        self.treatfiles=[]
        self.genomes={}

    def treatgbk(self,genbank_file_name):
        gbfilehandle=open(genbank_file_name,"r")
        parsed_genbank=list(SeqIO.parse(fh,"genbank"))
        if len(parsed_genbank)>1:
            return -1
            inputdic["number_of_records_inside"]=len(parsed_genbank)
            inputdic["number_of_gbfiles_in_folder"]=0
            if len(parsed_genbank)>1:
                raw_input("more than one record")
            self.treatfiles.append(genbank_file_name)

        

    def populateList(self,genome_directory):
        if os.path.isdir(genome_directory): 
            for filename in os.listdir(genome_directory):
                if filename[0]==".":
                    continue
                genome_base_dir=os.path.join(genome_directory,filename)
                if os.path.isdir(genome_base_dir):
                    genome=MyGenome(genome_directory,self.approved_gb_file_endings)
                    numrecords=genome.initialize()
                    if not genome_base_dir in self.genomes.keys():
                        self.genomes[genome_base_dir]=genome
                else:
                    print ("was not dir")
        return 1


