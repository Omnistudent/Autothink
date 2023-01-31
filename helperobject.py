import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
#rom Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class helperobject:


    
    def __init__(self):
        self.testvar="here we are, testing"

    @classmethod
    def get_lines_from_file(self,filename,removeheader,separator):
        self.outlines=0
        self.tempfh=open(filename,"r")
        self.cvread=self.tempfh.read()
        print ("read")
        self.cvsplit=self.cvread.split(separator)
        if removeheader==True:
            self.cvsplit= self.cvsplit[1:]
        self.outlines=self.cvsplit
        self.blanksremoved = list(filter(lambda x: x!="", self.cvsplit))

        #irfinderstring="einverted -sequence "+irfile+" -gap "+str(gap_penalty)+" -threshold "+str(min_score_threshold)+" -match "+str(match_score)+" -mismatch "+str(mismatch_score)+ " -outfile "+self.outfile+ " -outseq "+self.outseq+ " -maxrepeat "+str(maxrepeat)
        #print irfinderstring
        #os.system(irfinderstring)
        return self.blanksremoved

    @classmethod
    # Gets the starts and ends of footprints
    def find_starts_ends(self,querylist):
        lookingfor="start"
        results=[]
        fcounter=1
        currentstart=-1
        for i in querylist:
            if lookingfor=="start":
                if i==1:
                    currentstart=fcounter
                    lookingfor="end"
                    fcounter+=1
                    continue
            if lookingfor=="end":
                if i==0:
                    results.append([currentstart,fcounter-1])
                    lookingfor="start"
            fcounter+=1
        if lookingfor=="end":
            results.append([currentstart,len(querylist)])
        return results

    @classmethod
    def maketablemakestring(self,tablename):
        makestring="CREATE TABLE "+tablename+"("
        makestring+="readid"+" "+"TEXT"+","        
        makestring+="seq"+" "+"TEXT"+","        
        makestring+="seqlen"+" "+"INT"+"," 
        makestring=makestring[:-1]
        makestring+=")"
        return makestring

    @classmethod
    def removeFeature_for_curation(self,featureslist,toberemoved):
        newfeatures=[]
        for oldfeat in featureslist:
            if oldfeat.location.start==toberemoved.location.start and  oldfeat.location.end==toberemoved.location.end and oldfeat.qualifiers["is_name"]==toberemoved.qualifiers["is_name"]:
                continue
            else:
                newfeatures.append(oldfeat)
        return newfeatures

    @classmethod
    # Get part of nucleotide sequence
    def getseq(self,start,end,seq):
        return seq[start-1:end]

    @classmethod
    def remove_ending(text):
        rest = text.split(".", 1)[0]
        return res
    
    
