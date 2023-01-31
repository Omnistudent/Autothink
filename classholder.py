import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class InvertFinder:

    
    def __init__(self):
        self.testvar="here we are, testing"
        self.tempfile="/Users/security/invertfindertemp.fa"
        self.outfile="/Users/security/invertfindertemp_outfile.fa"
        self.outseq="/Users/security/invertfindertemp_outseq.fa"

    def getIRs(self,seq,maxrepeat):
        gap_penalty=12
        min_score_threshold=50
        match_score=3
        mismatch_score=-4

        irfile=self.makefasta(seq)
        irfinderstring="einverted -sequence "+irfile+" -gap "+str(gap_penalty)+" -threshold "+str(min_score_threshold)+" -match "+str(match_score)+" -mismatch "+str(mismatch_score)+ " -outfile "+self.outfile+ " -outseq "+self.outseq+ " -maxrepeat "+str(maxrepeat)
        #print irfinderstring
        os.system(irfinderstring)

    def makefasta(self,seq):
        s=Seq(str(seq),IUPAC.IUPACUnambiguousDNA())
        newrec=SeqRecord(s)
        newrec.id="invertfindertemp"
        newrec.description=""
        newrec.name=""
        tempfilefh=open(self.tempfile,"w")
        SeqIO.write([newrec],self.tempfile,"fasta")
        tempfilefh.close()
    
        return self.tempfile
        


class IsPair:

    
    def __init__(self,orig_feature,new_feature,mindistance,shortest_distance):
        self.orig_feature=orig_feature
        self.new_feature=new_feature
        self.mindistance=mindistance
        self.shortest_distance=str(shortest_distance)
        orle=int(float(new_feature.qualifiers["orflength"][0]))
        orle1=int(float(orig_feature.qualifiers["orflength"][0]))
        self.newOrfLenght=int(orle1+orle)/2

        if shortest_distance=="end_newstart":
            if orig_feature.strand==1:
                aareturns=[[self.orig_feature.qualifiers["sbjct_end"][0],"old_end_aa"],[self.new_feature.qualifiers["sbjct_start"][0],"new_start_aa"]]
                self.meeting_end=self.orig_feature.qualifiers["sbjct_end"][0]
                self.meeting_start=self.new_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_start=self.orig_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_end=self.new_feature.qualifiers["sbjct_end"][0]
            
                print ("type"+self.shortest_distance)
                print ("self.meeting_end"+self.meeting_end)
                print ("self.meeting_start"+self.meeting_start)
                print ("self.joined_aa_start"+self.joined_aa_start)
                print ("self.joined_aa_end"+self.joined_aa_end)

            elif orig_feature.strand==-1:
                self.meeting_end=self.new_feature.qualifiers["sbjct_end"][0]
                self.meeting_start=self.orig_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_start=self.new_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_end=self.orig_feature.qualifiers["sbjct_end"][0]
            
                print ("type"+self.shortest_distance)
                print ("self.meeting_end"+self.meeting_end)
                print ("self.meeting_start"+self.meeting_start)
                print ("self.joined_aa_start"+self.joined_aa_start)
                print ("self.joined_aa_end"+self.joined_aa_end)

        elif shortest_distance=="newend_start":
            if orig_feature.strand==1:
                self.meeting_end=self.new_feature.qualifiers["sbjct_end"][0]
                self.meeting_start=self.orig_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_start=self.new_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_end=self.orig_feature.qualifiers["sbjct_end"][0]
            
                print ("type"+self.shortest_distance)
                print ("self.meeting_end"+self.meeting_end)
                print ("self.meeting_start"+self.meeting_start)
                print ("self.joined_aa_start"+self.joined_aa_start)
                print ("self.joined_aa_end"+self.joined_aa_end)
            elif orig_feature.strand==-1:
                self.meeting_end=self.orig_feature.qualifiers["sbjct_end"][0]
                self.meeting_start=self.new_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_start=self.orig_feature.qualifiers["sbjct_start"][0]
                self.joined_aa_end=self.new_feature.qualifiers["sbjct_end"][0]
            
                print ("----------------type"+self.shortest_distance)
                print ("self.meeting_end"+self.meeting_end)
                print ("self.meeting_start"+self.meeting_start)
                print ("self.joined_aa_start"+self.joined_aa_start)
                print ("self.joined_aa_end"+self.joined_aa_end)
                #self.joined_aa_start="NA"
                #self.joined_aa_end="NA"
        else:
            #third typeend_newend
#third typenewstart_start
            raw_input("third type"+shortest_distance)

    def getShortestType(self):
        return self.shortest_distance

    def get_aa_orflength(self):
        return self.newOrfLenght

    def get_joined_aa_start(self):
        return self.joined_aa_start
    def get_joined_aa_end(self):
        return self.joined_aa_end

    def getStrand(self):
        return self.orig_feature.strand

    def getMinDistance(self):
        return self.mindistance

    def getNewCloseFeature(self):
        return self.new_feature

    def getOrigCloseFeature(self):
        return self.orig_feature

    def getOpposingStrands(self):
        if self.orig_feature.strand!=self.new_feature.strand:
            return True
        else:
            return False

    def getClosestDistance(self):
        aareturns=[[],[]]
        aareturns="NOOOOO"
        if self.getShortestType()=="newstart_start":
            returns=[int(self.new_feature.location.start),int(orig_feature.location.start),"newstart_start"]
        elif self.getShortestType()=="start_newstart":
            returns=[int(self.orig_feature.location.start),int(new_feature.location.start),"start_newstart"]
        elif self.getShortestType()=="newstart_end":
            returns=[int(self.new_feature.location.start),int(orig_feature.location.end),"newstart_end"]

        elif self.getShortestType()=="end_newstart":
            if self.orig_feature.strand==1:
                aareturns=[[self.orig_feature.qualifiers["sbjct_end"][0],"old_end_aa"],[self.new_feature.qualifiers["sbjct_start"][0],"new_start_aa"]]
            if self.orig_feature.strand==-1:
                aareturns=[[self.orig_feature.qualifiers["sbjct_start"][0],"old start aa"],[self.new_feature.qualifiers["sbjct_end"][0],"new_end_aa"]]
            returns=[int(self.orig_feature.location.end),int(self.new_feature.location.end),"end_newstart"]

        elif self.getShortestType()=="start_newend":
            returns=[int(self.orig_feature.location.start),int(self.new_feature.location.end),"start_newend"]
        elif self.getShortestType()=="newend_start":
            if self.orig_feature.strand==1:
                aareturns=[[self.new_feature.qualifiers["sbjct_end"][0],"new_end_aa"],[self.orig_feature.qualifiers["sbjct_start"][0],"old_start_aa"]]
            if self.orig_feature.strand==-1:
                aareturns=[[self.new_feature.qualifiers["sbjct_start"][0],"new_end_aa"],[self.orig_feature.qualifiers["sbjct_end"][0],"old_start_aa"]]






            returns=[int(self.new_feature.location.end),int(self.orig_feature.location.start),"newend_start"]
        elif self.getShortestType()=="end_newend":
            returns=[int(self.orig_feature.location.end),int(self.new_feature.location.end),"end_newend"]
        elif self.getShortestType()=="newend_end":
            returns=[int(self.new_feature.location.end),int(self.orig_feature.location.end),"newend_end"]
        else:
            returns=[0,0,"bokren"]
            raw_input("broken")

        return [returns,aareturns]

    def getOrigName(self):
        return str(self.orig_feature.qualifiers["is_name"][0]) 
    def getNewName(self):
        return str(self.new_feature.qualifiers["is_name"][0]) 

    def getName(self):
        nameslist=[str(self.getOrigName()).replace(" ",""),str(self.getNewName()).replace(" ","")] 
        nameslist.sort()
        return nameslist[0]+"|"+nameslist[1]



    def getStartOfBoth(self):

        location_of_start=int(self.orig_feature.location.start)
        location_of_end=int(self.orig_feature.location.end)
        location_of_closest_start=int(self.new_feature.location.start)
        location_of_closest_end=int(self.new_feature.location.end)
        return min(location_of_start,location_of_end,location_of_closest_start,location_of_closest_end)

    def getEndOfBoth(self):
        location_of_start=int(self.orig_feature.location.start)
        location_of_end=int(self.orig_feature.location.end)
        location_of_closest_start=int(self.new_feature.location.start)
        location_of_closest_end=int(self.new_feature.location.end)

        return max(location_of_start,location_of_end,location_of_closest_start,location_of_closest_end)

    def getSubjEnd(self,feature):
        if feature.strand==1:
            return str(feature.qualifiers["sbjct_end"][0])
        if feature.strand==-1:
            return str(feature.qualifiers["sbjct_start"][0])

    def getSubjStart(self,feature):
        if feature.strand==1:
            return str(feature.qualifiers["sbjct_start"][0])
        if feature.strand==-1:
            return str(feature.qualifiers["sbjct_end"][0])



class Drawer():
    def __init__(self):
        self.test="t"

    def getSubjEnd(self,feature):
        if feature.strand==1:
            return str(feature.qualifiers["sbjct_end"][0])
        if feature.strand==-1:
            return str(feature.qualifiers["sbjct_start"][0])

    def getSubjStart(self,feature):
        if feature.strand==1:
            return str(feature.qualifiers["sbjct_start"][0])
        if feature.strand==-1:
            return str(feature.qualifiers["sbjct_end"][0])

    def draw(self,start,end,genome):
        stringlen=170
        outstring=['-']*stringlen
        featurename_outstring=[' ']*stringlen
        featurekind_outstring=[' ']*stringlen
        feature_subj_pos=[' ']*stringlen

        fragment=genome[start:end]
        #print genome.seq[start:end]
        #featurenamelist=[]
        islist=["IS","ISstart","ISend","ISmiddle","wholeIS_0.8","completeIS_0.95","joinedhit"]
        for feature in fragment.features:
            if feature.type in islist:
                print ("start to end of draw:"+str(start)+"-"+str(end))

                st_start=int(float(feature.location.start)*(float(stringlen)/(end-start)))
                #featurenamelist.append(st_start*" "+str(feature.type))
                st_end=int(float(feature.location.end)*(float(stringlen)/(end-start)))
                outstring[st_start:st_end]=(1+st_end-st_start)*["I"]
                featurename_outstring[st_start-1:st_start+len(str(feature.type))-1]=str(feature.type)
                featurekind_outstring[st_start-1:st_start+len(str(feature.qualifiers["is_name"][0]))-1]=str(feature.qualifiers["is_name"][0])
                feature_subj_pos[st_start-1:st_start+len(self.getSubjStart(feature))-1]=self.getSubjStart(feature)
                feature_subj_pos[st_end-len(self.getSubjEnd(feature))-2:st_end-2]=str(self.getSubjEnd(feature))

                if feature.strand==1:
                    outstring[st_end]=">"
                if feature.strand==-1:
                    outstring[st_start-1]="<"
        for feature in fragment.features:
            if feature.type in islist:
                st_start=int(float(feature.location.start)*(float(stringlen)/(end-start)))
                st_end=int(float(feature.location.end)*(float(stringlen)/(end-start)))
                if feature.strand==1:
                    outstring[st_end]=">"
                if feature.strand==-1:
                    outstring[st_start-1]="<"
                #outstring[start:end]=(1+end-start)*["I"]
        print ("".join(featurename_outstring))
        print ("".join(outstring))
        print ("".join(featurekind_outstring))
        print ("".join(feature_subj_pos))
