################################################################
#   
#   Autothink.py
#   Written by theoden.vigil-stenman@su.se
#
#   - Idetifies transposase sequences by repeated blastx searches aginst a database of transposase aa sequences
#   - Joins neighboring hits whose aa sequences complement each other
#   - Output in genbank file
#
#   Needs the file classholder.py,helperobject,genomeholder (in path or same directory as autothink.py)
#
#################################################################

import os
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO
from os.path import join
import sqlite3
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pickle
from classholder import IsPair
from helperobject import helperobject
from genomeholder import Genomelist
import shutil
import argparse

# Removes the old features to replace with "joined" features
def removeFeature_for_curation(featureslist,toberemoved):
    newfeatures=[]
    for oldfeat in featureslist:
        if oldfeat.location.start==toberemoved.location.start and  oldfeat.location.end==toberemoved.location.end and oldfeat.qualifiers["is_name"]==toberemoved.qualifiers["is_name"]:
            continue
        else:
            newfeatures.append(oldfeat)
    return newfeatures

# Checks dictionary to see if position is already treated
def none_in_range_for_curation(filename, featurelocation):
    if filename in none_in_range_dic.keys():
        for entry in none_in_range_dic[filename]:
            if entry[0]==min(int(featurelocation.start),int(featurelocation.end)) and entry[1]==max(int(featurelocation.start),int(featurelocation.end)):
                return "none_in_range"
    return "none_in_dic" 

# Finds the closest feature of the right type
def findClosest_for_curation(myfeature, gbfile,maxdistance):
    pr_report("gbfile to open",gbfile,"findclosest_for_curation",1)

    dist_to_features_fh=open(gbfile,"r")
    features_to_check_distance_to=list(SeqIO.parse(dist_to_features_fh,"genbank"))[0]
    dist_to_features_fh.close()

    closefeaturesObs=[]
    for feature in features_to_check_distance_to.features:
        if feature.type in islist:
            checkstart= int(feature.location.start)
            checkend =int(feature.location.end)
            if checkstart==myfeature.location.start and checkend==myfeature.location.end:
                continue
            checkstart_start=abs(checkstart-int(myfeature.location.start))
            checkstart_end=abs(checkstart-int(myfeature.location.end))
            checkend_start=abs(checkend-int(myfeature.location.start))
            checkend_end=abs(checkend-int(myfeature.location.end))
            mindistance=min(checkstart_start,checkstart_end,checkend_start,checkend_end)
            shortest_distance="wrong_in_choose_shortest_distance"
            if checkstart_start==mindistance:
                if checkstart-int(myfeature.location.start)>0:
                    shortest_distance="start_newstart"
                else:
                    shortest_distance="newstart_start"
            elif checkstart_end==mindistance:
                if checkstart-int(myfeature.location.end)>0:
                    shortest_distance="end_newstart"
                else:
                    shortest_distance="newstart_end"
            elif checkend_start==mindistance:
                if checkend-int(myfeature.location.start)>0:
                    shortest_distance="start_newend"
                else:
                    shortest_distance="newend_start"
            elif checkend_end==mindistance:
                if checkend-int(myfeature.location.end)>0:
                    shortest_distance="newend_end"
                else:
                    shortest_distance="end_newend"
            if mindistance<maxdistance:
                closefeaturesObs.append(IsPair(myfeature,feature,mindistance,shortest_distance))
    sortedclosefeaturesObs=sorted(closefeaturesObs,key=lambda x: x.getMinDistance(), reverse=False)
    return sortedclosefeaturesObs

# Puts transposases with no neighbors in a file, so searches don't have made over and over
def append_to_nir_dic_for_curation(full_file_name,feature_location):
    if full_file_name in none_in_range_dic.keys():
        none_in_range_dic[full_file_name].append([min(int(feature_location.start),int(feature_location.end)),max(int(feature_location.start),int(feature_location.end))])
    else:
        none_in_range_dic[full_file_name]=[[min(int(feature_location.start),int(feature_location.end)),max(int(feature_location.start),int(feature_location.end))]]

def parse_gb_file_for_curation(full_file_name,orgname,read_pickle_dic,aminoacid_distance):
    try:
        gb_fh=open(full_file_name,"r")
        parsed=list(SeqIO.parse(gb_fh,"genbank"))
        gb_fh.close()
    except:
        print (full_file_name)
        input("could not open")
    if len(parsed)!=1:
        print ("num recs in "+full_file_name+":"+str(len(parsed)))
        input("number of records in file is not one")
    else:
        for record in parsed:
            featurelist=[]
            #s=Seq(str(record.seq),IUPAC.IUPACUnambiguousDNA())
            s=Seq(str(record.seq))
            newrec=SeqRecord(s)
            newrec.id=record.id
            newrec.description=record.description
            for feature in record.features:
                if feature.type in islist:
                    closest=findClosest_for_curation(feature,full_file_name,maxdistance)
                    if closest==[]:
                        if read_pickle_dic==False:
                            append_to_nir_dic_for_curation(orgname,feature.location)
                        continue
                    elif closest[0].getOpposingStrands():
                        continue
                    aaboutdistance=closest[0].getClosestDistance()[1]
                    aadistance2=abs(int(aaboutdistance[0][0])-int(aaboutdistance[1][0]))
                    if aadistance2<aminoacid_distance:
                        newfeatures=helperobject.removeFeature_for_curation(record.features,closest[0].getNewCloseFeature())
                        newfeatures2=helperobject.removeFeature_for_curation(newfeatures,closest[0].getOrigCloseFeature())



                        drawstart=closest[0].getStartOfBoth()
                        drawend=closest[0].getEndOfBoth()   
                        newfeature = SeqFeature(FeatureLocation(drawstart,drawend), strand=closest[0].getStrand(),type="joinedhit")

                        new_aa_length = int((abs(int(drawstart)-int(drawend)))/3)
                        print (new_aa_length)

                        newfeature.id="joined"
                        newfeature.qualifiers["is_name"]=[closest[0].getName()]
                        newfeature.qualifiers["sbjct_start"]=[closest[0].get_joined_aa_start()]
                        newfeature.qualifiers["sbjct_end"]=[closest[0].get_joined_aa_end()]
                        newfeature.qualifiers["orflength"]=[closest[0].get_aa_orflength()]
                        newfeature.qualifiers["ntlength"]=[closest[0].getShortestType()]
                        newfeature.qualifiers["perc_of_orf"]=str(round(float(new_aa_length)/float(closest[0].get_aa_orflength()),3))
                        newfeature.qualifiers["family"]="mixed"
                        newfeature.qualifiers["group"]="mixed"
                        newfeature.qualifiers["score"]="mixed"
                        newfeature.qualifiers["expected"]="mixed"
                        newfeature.qualifiers["frame"]="mixed"
                        newfeature.qualifiers["origin"]="mixed"
                        newfeature.qualifiers["numorfs"]="mixed"



                        
                        newrec.features.append(newfeature)

                        newfeatures2.append(newfeature)
                        newrec.features=newfeatures2
                        return newrec
        nodna=Seq("A")
        nonerec=SeqRecord(nodna)
        nonerec.id="none"
        return nonerec

# Get part of nucleotide sequence
def getseq(start,end,seq):
    return seq[start-1:end]

# blasts a footprint (?) against the database aa sequences. Returns the hit with the highest score 
def doblast(isseq,offset):
    # Write query to temporary file
    writequery=open(tempfile,"w")
    writequery.write(">test\n"+str(isseq))
    writequery.close()
    blastline=blast_program_directory+"blastx -db "+args.isblastdbfile+" -out "+tempfileout+" -query "+tempfile+" -query_gencode 11 -num_threads 4 -outfmt 5 -evalue "+ str(args.search1eval)
    os.system(blastline)

    try:
        blastout_handle=open(tempfileout,"r")
        blast_records = list(NCBIXML.parse(blastout_handle))	
    except:
        print ("blast parse trouble in function getbestblast")
        input("trouble")
        return None
    #collect all hsp hits in one list
    hsps=[]
    for record in blast_records:
        for alg in record.alignments:
            for hsp in alg.hsps:
                # Collect hit info in dictionary
                dic={}
                dic["hit_def"]=str(alg.hit_def)
                dic["sbjct_start"]=hsp.sbjct_start
                dic["match"]=hsp.match
                dic["identities"]=hsp.identities
                dic["positives"]=hsp.positives
                dic["sbjct_end"]=hsp.sbjct_end
                dic["expected"]=str(hsp.expect)
                dic["frame"]=hsp.frame
                dic["bits"]=hsp.bits
                dic["query"]=str(hsp.query)
                dic["mod_query_end"]=hsp.query_end+offset-1
                dic["mod_query_start"]=hsp.query_start+offset-1
                dic["query_end"]=hsp.query_end
                dic["query_start"]=hsp.query_start
                dic["sbjct"]=str(hsp.sbjct)
                dic["score"]=hsp.score
                dic["align_length"]=hsp.align_length
                dic["query_length"]=record.query_length
                dic["queried_seq"]=isseq
                dic["hit_seq"]=helperobject.getseq(int(min(int(hsp.query_start),int(hsp.query_end))),int(max(int(hsp.query_start),int(hsp.query_end))),isseq)
                hsps.append(dic)
    # Sort collected list by score
    mysortedhsps=sorted(hsps,key=lambda x: x["score"], reverse=True)
    if len(mysortedhsps)>0:
        # Return the hit with the hightest score, if any
        return mysortedhsps[0]
    else:
        return None


# Gets the starts and ends of footprints
def find_starts_ends(querylist):
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

# Parses the results of the first (genome vs aa database) search.
# Makes "footprints" of hits 
def parse_xml_file(sample,cur,samplename):
    replies=[]
    xmlfh=open(sample,"r")
    try:
        blast_records = list(NCBIXML.parse(xmlfh))	
    except:
        input("no file")
    if len(list(blast_records))<1:
        input("no records")
    parseitrecordcounter=0
    for rec in blast_records:
        parseitrecordcounter+=1
        hit_def="NONE"
        mycsvstring=""
        gbhitslist=[]
        if (len(rec.alignments))<1:
            input ("parseit: no algs for "+samplename)
            #continue
        hit_name=((str(rec.query)).split(" "))[0]
        cur.execute('''SELECT seq FROM genomes WHERE readid=?''', (samplename+hit_name,))
        ans=cur.fetchall()
        queryseq=ans[0][0].encode('ascii','ignore')

        # make genbank record
        #s=Seq(str(queryseq),IUPAC.IUPACUnambiguousDNA())
        s=Seq(str(queryseq))
        newrec=SeqRecord(s)
        newrec.id= str(rec.query).split(" ",1)[0]
        if len(str(newrec.id))>15:
            newrec.id= str(rec.query)[0:16]
        if len(str(rec.query).split(" "))>1:
            newrec.description= str(rec.query).split(" ",1)[1]
        else:
            newrec.description= rec.query
        allfeatures=[]
        querylist=[0]*rec.query_length
        algs=rec.alignments
        # Get coverage of hits (mark them on a string of "0" representing the genome)
        for alg in algs:
            for hsp in alg.hsps:

                if hsp.expect<args.search1eval:
                    querystart=min(hsp.query_start,hsp.query_end)
                    queryend=max(hsp.query_start,hsp.query_end)
                    querylist[querystart-1:queryend]=(1+queryend-querystart)*[1]
                    #feature = SeqFeature(FeatureLocation(hsp.query_start,hsp.query_end), strand=int(hsp.frame[1]),type="firsthit")
                    feature = SeqFeature(FeatureLocation(hsp.query_start,hsp.query_end), type="firsthit")
                    feature.id=alg.hit_id
                    feature.qualifiers["hit"]=alg.hit_def
                    allfeatures.append(feature)
        # Convert the footprints string into start and end numbers
        starts_and_ends=helperobject.find_starts_ends(querylist)
        #starts_and_ends=find_starts_ends(querylist)

        # Blast the footprint sequence, 
        for st_end in starts_and_ends:
            #feature = SeqFeature(FeatureLocation(st_end[0],st_end[1]), strand=1,type="footprint")
            feature = SeqFeature(FeatureLocation(st_end[0],st_end[1]), type="footprint")
            allfeatures.append(feature)

        # Blast footprints and put results in a list
        for res in starts_and_ends:
            footprintseq=queryseq[res[0]-1:res[1]]
            firstsearch=[footprintseq,1,len(footprintseq)]
            remainsearches=[firstsearch]
            while len(remainsearches)>0:
                currentsearch=remainsearches.pop()
                currentseq=currentsearch[0]
                oldstart=currentsearch[1]
                oldend=currentsearch[2]
                if len(currentseq)>4:
                    doblast_results=doblast(currentseq.decode('ASCII'),res[0])
                else:
                    doblast_results=None
                if (doblast_results==None) or (float(doblast_results["expected"])>args.search2eval):
                    pass
                else:
                    blast_hitstart=min(int(doblast_results["query_start"]),int(doblast_results["query_end"]))
                    blast_hitend=max(int(doblast_results["query_start"]),int(doblast_results["query_end"]))
                    recorded_blast_hitstart=blast_hitstart+oldstart-1
                    recorded_blast_hitend=blast_hitend+oldstart-1
                    listToAppend=[recorded_blast_hitstart+res[0],recorded_blast_hitend+res[0],doblast_results["hit_def"],doblast_results,res[0]]
                    gbhitslist.append(listToAppend)
                    if not blast_hitstart==1:
                        leftstart=oldstart
                        leftend=oldstart+blast_hitstart-2
                        leftremains_list=[helperobject.getseq(leftstart,leftend,footprintseq),leftstart,leftend]
                        remainsearches.append(leftremains_list)
                    if not blast_hitend==len(currentseq):
                        rightstart=oldstart+blast_hitend
                        rightend=oldend
                        rightremains_list=[helperobject.getseq(rightstart,rightend,footprintseq),rightstart,rightend]
                        remainsearches.append(rightremains_list)

        # It shouldn't matter if these hits are sorted or not
        sortedhits=sorted(gbhitslist,key=lambda x: x[3]["score"], reverse=True)


        # Make genbank features of the hits
        for hit in sortedhits:		
            if hit[3]["frame"][0]>=1:
                hitstrand=1
            if hit[3]["frame"][0]<0:
                hitstrand=-1
            splitname=hit[2].split("__")
            orflength=splitname[6].replace("orflength:","")
            subjstart=int(hit[3]["sbjct_start"])
            subjend=int(hit[3]["sbjct_end"])
            aahitlen=1+max(subjstart,subjend)-min(subjstart,subjend)
            perc_of_orf=round(aahitlen/float(orflength),3)
            minorf=min(subjstart,subjend)
            maxorf=max(subjstart,subjend)
            orftype="IS"

            # If detailed names of the kind of hit is desired
            if args.detailed_orfnames:
            #if detailed_orfnames:
                if perc_of_orf>=args.complete_cutoff: 
                #if perc_of_orf>=complete_cutoff: 
                    orftype="completeIS"+"_"+str(args.complete_cutoff)
                elif perc_of_orf>=0.8: 
                    orftype="wholeIS"+"_"+str(args.whole_cutoff)
                elif minorf<=float(orflength)*args.isstart_start_cutoff and maxorf<=float(orflength)*args.isstart_end_cutoff:
                    orftype="ISstart"
                elif minorf>=float(orflength)*args.isend_start_cutoff:
                    orftype="ISend"
                elif minorf>=float(orflength)*args.ismiddle_start_cutoff and maxorf<=float(orflength)*args.ismiddle_end_cutoff:
                    orftype="ISmiddle"

            feature = SeqFeature(FeatureLocation(int(hit[0]),int(hit[1])), strand=hitstrand,type=orftype)
            feature.id=hit[3]["hit_def"]
            family=splitname[1].replace("family:","")
            group=splitname[2].replace("group:","")
            origin=splitname[3].replace("origin:","")
            accession=splitname[4].replace("accession:","")
            ntlength=splitname[5].replace("ntlength:","")
            isnumorfs=splitname[7].replace("isnumorfs","")
            feature.qualifiers["is_name"]=splitname[0]

            feature.qualifiers["original_orf"]=str(orflength)
            feature.qualifiers["original_aa_length"]=str(aahitlen)

            feature.qualifiers["family"]=family
            feature.qualifiers["group"]=group
            feature.qualifiers["origin"]=origin
            feature.qualifiers["ntlength"]=ntlength
            feature.qualifiers["orflength"]=orflength
            feature.qualifiers["numorfs"]=isnumorfs
            feature.qualifiers["sbjct_start"]=str(hit[3]["sbjct_start"])
            feature.qualifiers["sbjct_end"]=str(hit[3]["sbjct_end"])
            feature.qualifiers["expected"]=str(hit[3]["expected"])
            feature.qualifiers["score"]=str(hit[3]["score"])
            feature.qualifiers["query_length"]=str(hit[3]["query_length"])
            #feature.qualifiers["match"]=str(hit[3]["match"])
            feature.qualifiers["query_start"]=str(hit[3]["query_start"])
            feature.qualifiers["identities"]=str(hit[3]["identities"])
            feature.qualifiers["align_length"]=str(hit[3]["align_length"])
            feature.qualifiers["positives"]=str(hit[3]["positives"])
            #feature.qualifiers["query"]=str(hit[3]["query"])
            #feature.qualifiers["queried_seq"]=str(hit[3]["queried_seq"])
            feature.qualifiers["query_end"]=str(hit[3]["query_end"])
            feature.qualifiers["frame"]=str(hit[3]["frame"])
            feature.qualifiers["bits"]=str(hit[3]["bits"])
            #feature.qualifiers["sbjct"]=str(hit[3]["sbjct"])
            feature.qualifiers["mod_query_start"]=str(hit[3]["mod_query_start"])
            feature.qualifiers["mod_query_end"]=str(hit[3]["mod_query_end"])
            feature.qualifiers["perc_of_orf"]=str(perc_of_orf)
            allfeatures.append(feature)
        newrec.features=allfeatures
        
        # Return the genome searhed (hit name), a genbank record, and empty string and filename(samplename) 
        replies.append([hit_name,newrec,mycsvstring,samplename])
    return replies
    # End of parsexml function


def treatGenbank(samplename,localgenomedirectory,resultsdir,tempresdir):
    pr_report("resultsdir",resultsdir,"general",1)
    sample=localgenomedirectory+"/"+samplename

    file_end=sample.split(".")[-1]
    if file_end in approved_endings:
    #if ".gbk" in sample or ".gb" in sample or ".gbff" in sample:
        fh=open(sample,"r")
        parsed_genbank=list(SeqIO.parse(fh,"genbank"))

        parse_results_dir=resultsdir

        for record in parsed_genbank:
            record.name=record.name+"_org_gb"
            #record.seq.alphabet=IUPAC.IUPACUnambiguousDNA()

       # save original gb file
        #orig_gb_file_loc=os.path.join(resultsdir,"original_gb_file_"+samplename)
        #orig_gb_file_fh=open(orig_gb_file_loc,"w")

        #SeqIO.write(parsed_genbank,orig_gb_file_fh,"genbank")
        #orig_gb_file_fh.close()
        

        for contigrecord in parsed_genbank:
            #if len(parsed_genbank)>1:
                #input("more than one record in parsed genbank")

            destination2=resultsdir+"/"+"blastres_folder/"+contigrecord.name+".xml"
            blastresfolder=resultsdir+"/"+"blastres_folder"
            
            if not os.path.exists(blastresfolder):
                os.makedirs(blastresfolder)


            SeqIO.write(contigrecord,open(tempfafile,"w"),"fasta")

            pr_report("blasting",sample,"general",1)

            # if xml file doesn't exist, create it by blasting

            if not os.path.isfile(destination2):
                blastx_cline2 = NcbiblastxCommandline(query=tempfafile, db=args.isblastdbfile, evalue=1e-6, outfmt=5, out=destination2,max_target_seqs=10000,num_threads=4,query_gencode=11)
                #blastx_cline2 = NcbiblastxCommandline(query=tempfafile, db=isblastdbfile, evalue=args.search1eval, outfmt=5, out=destination2,max_target_seqs=args.maxtargets,num_threads=4,query_gencode=11)
                print("blastcommand")
                print(blastx_cline2)
                os.system(blast_program_directory+str(blastx_cline2))

            # put genome in database
            genome_dropstring="DROP TABLE IF EXISTS "+sqlite_tablename
            genome_makestring=helperobject.maketablemakestring(sqlite_tablename)
            with makegenome_con:
                    makegenome_cur.execute(genome_dropstring)
                    makegenome_cur.execute(genome_makestring)
            #record=parsed_genbank[0]


            record=contigrecord
            readid=samplename+record.id
            seq= str(record.seq)
            seqlen=len(seq)
            insertline="INSERT INTO "+sqlite_tablename+" VALUES('"
            insertline+=readid+"','"
            insertline+=seq+"',"
            insertline+=str(seqlen)+","
            insertline= insertline[:-1]
            insertline= insertline+")"
            with makegenome_con:
                makegenome_cur.execute(insertline)

            # parse_xml_file returns a list of
            # hit_name,newrec,mycsvstring,samplename
            # it returns ONE genbank record hit
            parseresults2=parse_xml_file(destination2,makegenome_cur,samplename)
           # C:\Users\Eris\Documents\autothinktestfolder\outputdir\bmg51.gb/bmg51.gb/blastres_folder/bmg51_org_gb.xml
#bmg51.gb


            # if no results from parse_xml_file:
            if len(parseresults2)<1:
                input("no replies form parsexml")
                genbank_output=open(resultddir+"/"+str(contigrecord.name)+".gbk","w")
                if len(str(record.id))>15:
                    record.id= str(record.id)[0:16]
                if len(str(record.name))>15:
                    record.name= str(record.name)[0:16]
                #record.seq= Seq(str(record.seq),IUPAC.IUPACUnambiguousDNA())
                record.seq= Seq(str(record.seq))
                #record.name=""
                record.description=""
                record.features=[]
                #test
                record.annotations["molecule_type"] = "DNA"
                #test
                SeqIO.write([record],genbank_output,"genbank")
                genbank_output.close()

            curatedfilename=""
            curatedmarkerword="_Curated"

            for parseresult2 in parseresults2:
                hit_name=parseresult2[0]
                gbrec=parseresult2[1]
                gbrec.annotations["molecule_type"] = "DNA"
                #csvline=parseresult2[2]
                mysamplename=parseresult2[3]

                curatedfilename=resultsdir+"/"+curatedmarkerword+str(samplename)+"/"+str(contigrecord.name)+".gbk"

                altparsegbdir=resultsdir+"/"+curatedmarkerword+str(samplename)
                

                if not os.path.exists(altparsegbdir):
                    os.makedirs(altparsegbdir)
                #C:\Users\Eris\Documents\autothinktestfolder\outputdir\bmg51.gb/bmg51.gb/_Curatedbmg51.gb/bmg51_org_gb.gbk
                
                genbank_output2=open(curatedfilename,"w")
                SeqIO.write([gbrec],genbank_output2,"genbank")
                genbank_output2.close()

            # Make dictionary for curation

            counter=0
            joinparseresults=""
            #while joinparseresults!=None:
            counter+=1
            if counter>1:
                break
            joinparseresults2=parse_gb_file_for_curation(curatedfilename,contigrecord.name,False,aadistance)
            pickle.dump(none_in_range_dic, open(pickledic,"wb"))

            # Perform the actual curation
            counter=0
            infile=os.path.join(resultsdir,str(contigrecord.name)+".gbk")
            joinparseresults=parse_gb_file_for_curation(curatedfilename,contigrecord.name,True,aadistance)

            curated_directory=os.path.join(os.path.join(resultsdir,filename),curatedmarkerword)

            while joinparseresults.id!="none":
                counter+=1
                joinparseresults=parse_gb_file_for_curation(curatedfilename,contigrecord.name,True,aadistance)
                #test
                joinparseresults.annotations["molecule_type"] = "DNA"
                #test is it really dna_
                
                curated_directory=os.path.join(os.path.join(resultsdir,filename),curatedmarkerword)
                if not joinparseresults.id=="none":
                    outfile=curatedfilename

                    #if not os.path.exists(curated_directory+samplename):
                    #    os.makedirs(curated_directory+samplename)
                    #    raw_input("made dirs "+str(curated_directory+samplename))

                    outfh=open(outfile,"w")
                    SeqIO.write([joinparseresults],outfh,"genbank")
                    outfh.close()
                    infile=outfile
                
            if counter==0:
                if not os.path.exists(curated_directory+samplename):
                    os.makedirs(curated_directory+samplename)
                #shutil.copy(parse_results_dir+samplename+"/"+str(contigrecord.name)+".gbk",curated_directory+samplename+"/"+str(contigrecord.name)+".gbk")
                #pr_report("curated samplebname",curated_directory+samplename,"general",2)

        return resultsdir+"/"+curatedmarkerword+samplename+"/"+str(contigrecord.name)+".gbk" 
        #return infile
    else:
        print (sample+" does not have gb, gbk or gbff in file name")
        return resultsdir+"/"+curatedmarkerword+samplename+"/"+str(contigrecord.name)+".gbk" 
        #return infile

def remove_ending(text):
    rest = text.split(".", 1)[0]
    return rest

def remove_footprints_from_gb_file(sentpath,store):
    
    # initialize rangedic
    is_fraction_counts={}
    rangedic={}
    for r in onerange:
        rangedic[r]=0
        is_fraction_counts[r]=0
    
    # parse sent genbank
    fh=open(sentpath,"r")
    internalparsed_genbank=list(SeqIO.parse(fh,"genbank"))
    fh.close()

    if len(internalparsed_genbank)>1:
        pr_report("more than one record in remove_footprint_from_gb_file, length:",str(len(internalparsed_genbank)),"general",1)

    rec=internalparsed_genbank[0]
    allfeats=[]

    # prepare_out_csv
    #csv="name\tstart\tend\tlength\tid\ttype\textract\t"
    #csv+="isfamily"+SEP
    #csv+="isgroup"+SEP
    #csv+="is_score"+SEP+"is_expected"+SEP+"is_frame"+SEP+"is_perc_of_orf"+SEP+"is_origin"+SEP

    csv=""
    for header in ishit_headers: 
        csv+=header+SEP

    for ik in onerange:
        csv+=str(int(ik)-10)+"tolessthan"+str(ik)+SEP

    csv+="\n"

    removedcsv=csv
    hit_csv_string=csv

    remove_firsthit=True
    remove_footprint=True

    cleaned_features=[]


    for feat in rec.features:
        if remove_firsthit==True:
            if str(feat.type).lower()=="firsthit":
                continue

        if remove_footprint==True:
            if str(feat.type).lower()=="footprint":
                continue

        cleaned_features.append(feat)

        perc_of_orf="unknown"

        #insertline="INSERT INTO "+sqlite_tablename+" VALUES('"
        #insertline+=readid+"','"
        #insertline+=seq+"',"
        #insertline+=str(seqlen)+","
        #insertline= insertline[:-1]
        #insertline= insertline+")"
        #with makegenome_con:
        #    makegenome_cur.execute(insertline)

        hit_qualifiers= ["is_name","family","group","score","expected","frame","perc_of_orf","origin","numorfs"]

        hit_csv_string+=str(feat.type)+SEP
        hit_csv_string+=str(feat.location.start)+SEP
        hit_csv_string+=str(feat.location.end)+SEP
        hit_csv_string+=str(1+abs(int(feat.location.start)-int(feat.location.end)))+SEP #length
        hit_csv_string+=str(feat.id)+SEP
        for hit_qual in hit_qualifiers:
            if hit_qual in feat.qualifiers.keys():
                hit_csv_string+=str(feat.qualifiers[hit_qual][0])+SEP
                if hit_qual=="perc_of_orf":
                    perc_of_orf= feat.qualifiers[hit_qual][0]
            else:
                hit_csv_string+="unknown"+SEP

        if perc_of_orf!="unknown":
            for upperlimit in onerange:
                if float(perc_of_orf)<float(upperlimit)*0.01:
                    uppernum=str(upperlimit)
                    lowernum=str(int(upperlimit)-10)
                    is_fraction_counts[upperlimit]+=1
                    break

        for uppervar in onerange:
            hit_csv_string+=str(is_fraction_counts[uppervar])+SEP

        for hit_qual in feat.qualifiers.keys():
            hit_csv_string+=str(hit_qual)+":"+str(feat.qualifiers[hit_qual][0])+SEP

        hit_csv_string+="\n"

    fho3=open(store+"_firsthits_removed.csv","w")
    fho3.write(hit_csv_string)
    fho3.close()

    fho2=open(store+"all.csv","w")
    fho2.write(hit_csv_string)
    fho2.close()
    
    cleanedrec=rec
    cleanedrec.features=cleaned_features
    fho=open(sentpath+"_cleaned.gbk","w")
    SeqIO.write(cleanedrec,fho,"genbank")
    fho.close()


    return is_fraction_counts


def pr_report(var,val,subject,lvl):
    if subject in approvedsubjects:
        if int(lvl)>=int(approvedlevel):
            print (str(var)+ "_ val:"+str(val)+"_")

def get_most_freq_in_dic(dic):
    val=0
    k=""
    for key, value in sorted(dic.items(), key=lambda k,v: (v,k)):
    #for key, value in sorted(dic.iteritems(), key=lambda k,v: (v,k)):
        #for key, value in sorted(mydict.iteritems(), key=lambda (k,v): (v,k)):
        #for k in sorted(d, key=lambda k: d[k], reverse=True)
        val=int(value)
        k=str(key)
    return (val,k)


def getfilestats(mydir,sumcsvfile,rangedic):
    totlen=0

    csvlines=helperobject.get_lines_from_file(mydir,True,"\n")
    namesdic={}
    isnum=0

    numkeylist=[]
    for i in rangedic.keys():
        numkeylist.append(int(i))

    fractionheader=""
    for ik in sorted(numkeylist,reverse=False):
        fractionheader+=str(int(ik)-10)+"tolessthan"+str(ik)+SEP

    for line  in csvlines:
        if line=="":
            continue
        pr_report("line",line,"general",1)
        splitline=line.split("\t")
        totlen+=int(splitline[3])
        name=str(splitline[5])
        print (splitline)
        print (name)
        isnum+=1
        if name in namesdic.keys():
            namesdic[name]+=1
        else:
            namesdic[name]=1
    #num1=get_most_freq_in_dic(namesdic)[1]
    num1=1
    num2=2
    #num2=get_most_freq_in_dic(namesdic)[0]

    outcsv="number_of_is"+SEP+"total_is_nt\tmost_numerous_is"+SEP+"most_numerous_is_count"+"\t"
    #outstuff=write_to_tab_line([totalcount,totalnts,maxkey,dickdic[maxkey]])
    outcsv+=fractionheader
    outcsv+="\n"


    outcsv+=str(isnum)+SEP
    outcsv+=str(totlen)+"\t"+str(num1).replace("[","").replace("[","")+"\t"
    outcsv+=str(num2)+"\t"
    keylist=rangedic.keys()

    for u in sorted(numkeylist):
        outcsv+=str(rangedic[str(u)])+"\t"




    outcsv+="\n"
    sumcsvfh=open(sumcsvfile,"w")
    sumcsvfh.write(outcsv)
    sumcsvfh.close()


def write_to_tab_line(mylist): 
    outstring=""
    for head in mylist:
        outstring+=str(head)+SEP
    return outstring

# Gets total nucleotides for all files in dir, number of nt that is IS, most frequent IS and the number of the most freq is
def getdirstats(mydir,gb,sumcsvfile):
#nerange=['5', '15', '25', '35', '45', '55', '65', '75', '85', '95', '105', '115'
    
    # get list of files in dir, ignore osx system files
    files=list(filter(lambda x: x[0]!=".",os.listdir(mydir)))

    #list(filter(lambda x: x!="", self.cvsplit)

    totalnts=0
    totalcount=0
    totaltot=""
    dickdic={}
    rangedic={}
    
    # initialize dic with upper limit as key, 0 as value
    for uppernum in onerange:
        rangedic[uppernum]=0

    outcline=""
    #outcline="totalis"+SEP
    #outcline="totalisnt"+SEP
    #outcline="most_freq_is"+SEP
    #outcline="most_freq_is_num"+SEP
    headers=write_to_tab_line(["totalis","totalisnt","most_freq_is", "most_freq_is_num"])
    outcline+=headers
    for n in range_list_names_global:
        outcline+=n+SEP
    print ("tt"+outcline+"tt")

    outcline+="\n"

    for f in files:
        if f[0]==".":
            continue
        if os.path.isdir(os.path.join(mydir,f)) and gb==False:

            summarypath=mydir+"/"+f+"/"+f+"_summary.csv"

            cvplines=helperobject.get_lines_from_file(summarypath,True,"\n")

            for u in cvplines:
                tabsp=u.split(SEP)
                totalcount+=int(tabsp[0])
                totalnts+=int(tabsp[1])

                for tabn,ind in zip(range(4,16),onerange):
                    rangedic[ind]+=int(tabsp[tabn])

                mostis= tabsp[2]
                mostisnum= int(tabsp[3])

                if mostis in dickdic.keys():
                    dickdic[mostis]+=mostisnum
                else:
                    dickdic[mostis]=mostisnum

        if gb==True and os.path.exists(mydir+"/"+f+"_summary.csv"):

            gb_summarypath=mydir+"/"+f+"_summary.csv"

            cvplines2=helperobject.get_lines_from_file(gb_summarypath,True,"\n")

            for u in cvplines2:
                tabsp=u.split(SEP)
                totalnts+=int(tabsp[1])
                totalcount+=int(tabsp[0])

                for tabn,ind in zip(range(4,16),onerange):
                    rangedic[ind]+=int(tabsp[tabn])
                    print (ind)

                mostis= tabsp[2]
                mostisnum= int(tabsp[3])

                if mostis in dickdic.keys():
                    dickdic[mostis]+=mostisnum
                else:
                    dickdic[mostis]=mostisnum


    maxkey=max(dickdic, key=dickdic.get)
    sorted(dickdic.items(), key=lambda x: x[1],reverse=True)

    outstuff=write_to_tab_line([totalcount,totalnts,maxkey,dickdic[maxkey]])
    print (outstuff)
    
    #outcline+=str(totalcount)+SEP
    #outcline+=str(totalnts)+"\t"
    #outcline+=str(maxkey)+"\t"
    #outcline+=str(dickdic[maxkey])+"\t"

    outcline+=outstuff

    for ind in onerange:
        outcline+=str(rangedic[ind])+SEP

    print ("tt"+outcline+"tt")

    outcline+="\n"
    outcfh=open(sumcsvfile,"w")
    outcfh.write(outcline)
    outcfh.close()


##############################
#   VARIABLES
##############################
parser = argparse.ArgumentParser()


# Location of blast program

#blast_program_directory="C:/NCBI/blast-BLAST_VERSION+/bin/"

script_folder=os.path.dirname(__file__)+"/"

#parser.add_argument(        "--genomesummaryfile",      help="name of csv file with genome summary, default=prefix+genomesummary.csv", default=script_folder+"genomesummary.csv",   required=False)
#parser.add_argument("-o",   "--outputfile",             help="destionation of summary file, default=prefix+output.csv", default=script_folder+"output.csv",   required=False)
parser.add_argument("-O",   "--outputdir",              help="directory output, default=prefix+outputdir",              default=script_folder+"output",       required=False)
parser.add_argument("-c",   "--files_to_check_for_is",  help="csv for total results, default=basedir+total_output.csv",      default=script_folder+"files_to_check_for_is.csv",       required=False)
parser.add_argument(        "--search1eval",            help="evalue cutoff for first search, default=1e-6",            default=1e-6,                                           type=float)
parser.add_argument(        "--search2eval",            help="evalue cutoff for second search, default=1e-6",           default=1e-6,                                           type=float)
parser.add_argument(        "--maxtargets",             help="max blast records to return, default=10000",              default=10000,                                          type=int)
parser.add_argument(        "--maxdistance",            help="Max distance between hits to be considered for joining, default=1000",              default=1000,                                          type=int)
parser.add_argument(        "--detailed_orfnames",      help="Annotate hits as start, end, middle?, default=True",      default=True,                       required=False,                     action="store_true")
parser.add_argument(        "--complete_cutoff",        help="cutoff for 'complete' annotation, default=0.95",                                                                                          default=0.95,   type=float)
parser.add_argument(        "--whole_cutoff",           help="cutoff for 'whole' annotation, default=0.8",                                                                                              default=0.8,    type=float)
parser.add_argument(        "--isstart_start_cutoff",   help="To be considered to be a 'start', hit start position is below X * protein aa length, default=0.3",                                        default=0.3,    type=float)
parser.add_argument(        "--isstart_end_cutoff",     help="To be considered to be an 'start', hit end position is below X * protein aa length, default=0.5",                                         default=0.5,    type=float)
parser.add_argument(        "--isend_start_cutoff",     help=" To be considered to be an 'end', hit start position is above isend_start_cutoff * protein aa length, default=0.7",                       default=0.7,    type=float)
parser.add_argument(        "--ismiddle_start_cutoff",  help=" To be considered to be a 'middle', hit start and end positions are between X * protein aa length and Y* protein aa length, default=0.3", default=0.3,    type=float)
parser.add_argument(        "--ismiddle_end_cutoff",    help=" To be considered to be a 'middle', hit start and end positions are between X * protein aa length and Y* protein aa length, default=0.7", default=0.7,    type=float)
parser.add_argument(        "--genome_directory",       help="location of folder with genomes to be treated",           default=script_folder+"/"+"testgenomes2/",  required=False)
parser.add_argument(        "--isblastdbfile",          help="formatted database of amino acids",                       default=script_folder+"/"+"is_aa_30_nov2016.fa",  required=False)
parser.add_argument(        "--sqlite_tablename",       help="table name for sql genome databse, default=genomes",      default="genomes")
parser.add_argument(        "--aadistance",             help="for joining, default=50",                                 default=50,    type=int) #old val 0.3
parser.add_argument(        "--blast_location",         help="path to blast executable dir",                                 default="C:/NCBI/blast-BLAST_VERSION+/bin/",required=False) 
args=parser.parse_args()

print (args)


# csv file with genomes to treat

approvedsubjects=["firstparse","blast","curate","posttreat","general"]
approvedlevel=0

results_directory=args.outputdir

output_csv=args.files_to_check_for_is

pickledic=results_directory+"autocurate_pickled_dick.pydic"
genome_db_destination=results_directory+"genomedb.db"
hit_db_destination=results_directory+"hits.db"
tempfile=results_directory+"/temp_blastfile.fa"
tempfileout=results_directory+"temp_idblastresult.xml"
curated_directory=results_directory+"curated/"
tempfafile=script_folder+"/collection.fa"
blast_program_directory=args.blast_location
sqlite_tablename="genomes"
aadistance=args.aadistance
SEP="\t"


# List of feature types to be treated
islist=["IS","ISstart","ISend","ISmiddle","wholeIS_0.8","completeIS_0.95","joinedhit"]

# Max distance between hits to be considered for joining
maxdistance=args.maxdistance

none_in_range_dic={}

onerange=['5', '15', '25', '35', '45', '55', '65', '75', '85', '95', '105', '115']
range_list_names_global=['-5_5', '5_15', '15_25', '25_35', '35_45', '45_55', '55_65', '65_75', '75_85', '85_95', '95_105', '105_115']

ishit_headers=["name","start","end","length","id","type","isfamily","isgroup","is_score","is_expected","is_frame","is_perc_of_orf","is_origin"]
finalcsvheaders=["organism","nt covered by is","number of is"]

#hit_features=["is_name","family","group", "origin", "ntlength", "orflength", "numorfs", "sbjct_start", "sbjct_end", "expected", "score", "query_length", "match", "query_start", "identities", "align_length", "positives", "query", "queried_seq", "query_end", "frame", "bits", "sbjct", "mod_query_start", "mod_query_end","perc_of_orf","original_orf","original_aa_length"]

approved_endings= ["gbk", "gb","gbff"]


##################################
#   START OF SCRIPT EXECUTION
##################################

if not os.path.exists(results_directory):
    os.makedirs(results_directory)

# Make connections for genome database
makegenome_con = sqlite3.connect(genome_db_destination)
makegenome_cur = makegenome_con.cursor()

# Make connections for hit database
makehits_con = sqlite3.connect(hit_db_destination)
makehits_cur = makehits_con.cursor()

# Open the input csv and parse it

outputcsvline=""

first=True
headerholder={}






#################################
# For each file in the genome directory, make results file structure, write info to args.outputcsv
###################################

for filename in os.listdir(args.genome_directory):
    print(filename)
    if filename[0]==".":
        continue
    inputdic={}
    inputdic["filename"]=filename
    inputdic["files_to_treat"]=[]
    inputdic["directory"]=args.genome_directory
    inputdic["resultsdirectory"]=str(os.path.join(results_directory,filename))

    # if directory
    if os.path.isdir(os.path.join(args.genome_directory,filename)):
        inputdic["filetype"]="directory"
        inputdic["number_of_records_inside"]=0
        inputdic["number_of_gbfiles_in_folder"]=0

        genome_dir=os.path.join(args.genome_directory,filename)
        for genomefile in os.listdir(genome_dir):
            if genomefile[0]==".":
                continue
            internalsample=os.path.join(genome_dir,genomefile)
            if os.path.isdir(internalsample):
                input("directory within directory in startgenome folder. Quitting") 
                quit()
            if ".gbk" in genomefile or ".gb" in genomefile or ".gbff" in genomefile:
                
                if not os.path.exists(os.path.join(os.path.join(results_directory,filename),genomefile)):
                    os.makedirs(os.path.join(os.path.join(results_directory,filename),genomefile))
                    pr_report("making dir structure", str(os.path.join(os.path.join(results_directory,filename),genomefile)),"general",1)
                
                inputdic["files_to_treat"].append(str(genomefile))

                fh=open(internalsample,"r")
                internalparsed_genbank=list(SeqIO.parse(fh,"genbank"))
                if len(internalparsed_genbank)!=1:
                    input("multiple records in files in directory in startgenome folder. Quitting") 
                    quit()
                inputdic["number_of_gbfiles_in_folder"]+=1
        filetreattxt=""

        # fill in treat files text
        for f in inputdic["files_to_treat"]:
            filetreattxt+=f+"\n"
        # write treat files text
        filetret_file_dir=os.path.join(results_directory,filename)
        filetret_file2=os.path.join(filetret_file_dir,"treatfiles_"+remove_ending(filename)+".txt")
        filetretfh=open(filetret_file2,"w")
        filetretfh.write(filetreattxt)
        filetretfh.close()
    # if not a directory
    else:
        if ".gbk" in filename or ".gb" in filename or ".gbff" in filename:
            if not os.path.exists(os.path.join(results_directory,filename)):
                os.makedirs(os.path.join(results_directory,filename))
                pr_report("making dir structure", str(os.path.join(results_directory,filename)),"general",1)
            inputdic["filetype"]="gbfile"
            sample=os.path.join(args.genome_directory,filename)
            inputdic["files_to_treat"].append(str(filename))
            fh=open(sample,"r")
            parsed_genbank=list(SeqIO.parse(fh,"genbank"))
            inputdic["number_of_records_inside"]=len(parsed_genbank)
            inputdic["number_of_gbfiles_in_folder"]=0
            inputdic["genome_nucleotide_size"]=len(parsed_genbank[0].seq)
            inputdic["organism"]=len(parsed_genbank[0].seq)
            inputdic["source"]=len(parsed_genbank[0].seq)

            #accesstions=(parsed_genbank[0].annotations["accessions"])

#annotations', 'count', 'dbxrefs', 'description', 'features', 'format', 'id', 'islower', 'isupper', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper


            if len(parsed_genbank)>1:
                input("more than one record")

            filetreattxt=""
            for f in inputdic["files_to_treat"]:
                filetreattxt+=f+"\n"
            filetret_file_dir=os.path.join(results_directory,filename)
            filetret_file2=os.path.join(filetret_file_dir, "treatfiles_"+remove_ending(filename)+".txt")
            filetretfh=open(filetret_file2,"w")
            filetretfh.write(filetreattxt)
            filetretfh.close()
        else:
            pr_report("non-gb file. Quitting.",filename,"general",1)
            input("exit")
            quit()

    outputstats=""

    if first:
        for k in inputdic.keys():
            outputstats+=str(k)+"\t"
        outputstats+="\n"
    first=False

    for k in inputdic.keys():
        outputstats+=str(inputdic[k])+SEP
    outputstats+="\n"
    outputcsvline+=outputstats
    headerholder=inputdic


# Write csv file with file info
outputcsv_fh=open(output_csv,"w")
outputcsv_fh.write(outputcsvline)
outputcsv_fh.close()

print ("wrote to ")
print (output_csv)



#################################
# Reads args.outputcsv
#
#reads genomes from treatfiles_ + the genome name from arts.outputcsv
# For each genome name (allways 1?) it open the genome, adds nt to totalnt in the resultsdic
# for each contigrecord in the genome, it does treatGenbank
# for each contigrecord, treatGenbank extracts the nucleotide sequcence and blasts it ,storing the results in a blastres folder inside the genome folder.
# If there is already a results xml file here, the blast is not performed, the old file is used.
# a folder it makes
# it also stores the sequence in a database
# the results are then parsed by parse_xml_file
###################################


# Read csv file with file info
inputcsv=output_csv
print("reading from inputcsv")
print(inputcsv)

splitinput=helperobject.get_lines_from_file(inputcsv,True,"\n")

allfiles_resultsdics=[]


# for earch line in the input csv
for line in splitinput:
    print("line in splitinput")
    print(line)
    pr_report("line_in_splitinput",line,"splitinput",1)
    one_entry_resultsdic={}
    one_entry_resultsdic["totalnt"]=0
    one_entry_resultsdic["totalfiles"]=0
    one_entry_resultsdic["mineval_first_search"]=args.search1eval
    one_entry_resultsdic["mineval_id_searches"]=args.search2eval
    one_entry_resultsdic["summaryfile"]=""
    
    outline=""
    if line=="":
        continue

    one_entry_resultsdic["inputpart"]=line.replace("\n","")

    tabs=line.split("\t")
    #pr_report("tabs",tabs,"general",1)
    filename=tabs[2]
    genomedir=tabs[2]
    print("filename")
    print(filename)
    #genomedir
    #pr_report("genomedir",genomedir,"general",1)
    
    newfileline=os.path.join(results_directory,tabs[3])
   # C:\Users\Eris\Documents\autothinktestfolder\outputdir\bmg51.gb
    newfileline2=os.path.join(newfileline,"treatfiles_"+remove_ending(tabs[0])+".txt")
    treatfileslist=helperobject.get_lines_from_file(newfileline2,False,"\n")

    genomeresultsdir="" 

    last_finished_file=""
    designated_results_dir=tabs[4]

    for row in treatfileslist:
        if ".gbk" in row or ".gb" in row or ".gbff" in row:
            sample=args.genome_directory+"/"+row

            genomeresultsdir=os.path.join(results_directory,genomedir)
            temp_resultsdir=""
            if tabs[1]=="directory":
                sample=args.genome_directory+genomedir+"/"+row
                sampledir=args.genome_directory+genomedir+"/"

                resultdsdirplusfilenamenoend=os.path.join(results_directory,remove_ending(genomedir))
                genomeresultsdir=os.path.join(resultdsdirplusfilenamenoend,row)
                #genomeresultsdir=os.path.join(os.path.join(results_directory,remove_ending(genomedir)),row)
            else:
                sampledir=args.genome_directory
                genomeresultsdir=os.path.join(results_directory,genomedir)
                genomeresultsdir=os.path.join(results_directory,tabs[0])

                

            fh=open(sample,"r")
            parsed_genbank=list(SeqIO.parse(fh,"genbank"))
            one_entry_resultsdic["totalnt"]+=len(parsed_genbank[0].seq)
            one_entry_resultsdic["totalfiles"]+=1

            if len(parsed_genbank)>1:
                input("more than one record")
            for contigrecord in parsed_genbank:
                last_finished_file=treatGenbank(row,sampledir,genomeresultsdir+"/"+row,temp_resultsdir)
            fractioncount=remove_footprints_from_gb_file(last_finished_file,genomeresultsdir+"/"+row)
            getfilestats(genomeresultsdir+"/"+row+"_firsthits_removed.csv",genomeresultsdir+"/"+row+"_summary.csv",fractioncount)
            print (genomeresultsdir+"/"+row+"_summary.csv")
            one_entry_resultsdic["summaryfile"]=str(genomeresultsdir+"/"+row+"_summary.csv")

    if tabs[1]=="directory":
        #getdirstats(designated_results_dir,False,designated_results_dir+"total_summary.csv")
        print (designated_results_dir)
        #takes dir to treat, if single file, outputlocation
        getdirstats(designated_results_dir,False,args.genomesummaryfile)
    #if tabs[1]=="gbfile":
        #getdirstats(designated_results_dir,True,designated_results_dir+"total_summary.csv")
        #getdirstats(designated_results_dir,True,args.genomesummaryfile)



    allfiles_resultsdics.append(one_entry_resultsdic)


finalcsvline=""
for u in finalcsvheaders:
    finalcsvline+=u +SEP
for u in range_list_names_global:
    finalcsvline+=u +SEP
finalcsvline+="\n"

for line in splitinput: 
    continue
    isnts=0
    isnumber=0
    tempdic={}
    if line=="":
        continue
    tabl=line.split(SEP)
    if len(tabl)>5:
        statscsvf=results_directory+"/"+str(tabl[2])+args.genomesummaryfile
        #statscsvf=results_directory+"/"+str(tabl[2])+"total_summary.csv"
        finalcsvline+=str(tabl[2])+"\t"

        hugsplit=helperobject.get_lines_from_file(statscsvf,False,"\n")
        haders=hugsplit[0]

        hugtabs=hugsplit[1].split(SEP)
        if len(hugtabs)>2:

            isnts=int(hugtabs[1])
            isnumber=int(hugtabs[0])

        else:
            input("no len")

        finalcsvline+=str(isnts)+SEP
        finalcsvline+=str(isnumber)+SEP


        
        for ind in (range(4,16)):
            finalcsvline+=str(hugtabs[ind])+SEP


    else:
        pr_report("short line", tabl,"general",1)
        input("short line")
    finalcsvline+="\n"

#finaloutcsv=results_directory+"/"+"total_summary.csv"
#finaloutcsvfh=open(finaloutcsv,"w")
#finaloutcsvfh.write(finalcsvline)
#finaloutcsvfh.close()

# write resultcsv

first=True
resultstats=""
resultstats2=""
resultcsvline=""


for m in headerholder.keys():
        print (m)
#        raw_input("prant headerholder keys")
        resultstats+=m+SEP
resultstats+='totalis\ttotalisnt\tmost_freq_is\tmost_freq_is_num\t-5_5\t5_15\t15_25\t25_35\t35_45\t45_55\t55_65\t65_75\t75_85\t85_95\t95_105\t105_115\t'
headers_orgkeys=['inputpart', 'totalfiles', 'mineval_id_searches', 'totalnt', 'summaryfile', 'mineval_first_search']

print (resultstats)

#resultstats+="totalis","totalisnt","most_freq_is","most_freq_is_num","-5_5","5_15","15_25","25_35","35_45","45_55","55_65","65_75","75_85","85_95","95_105","105_115"

heds=["myfilestotreat","myfiletype","myfilename","numcrecords","resultsdirI","gendir","numgbfilesinfolder","totalfiles","mineval","totalnt","summaryfile","evalfirstsearch","totalis","totalisnt","most_freq_is","most_freq_is_num","-5_5","5_15","15_25","25_35","35_45","45_55","55_65","65_75","75_85","85_95","95_105","105_115"]

resultstats+="\n"
resultstats2+="\n"
resultstats+=write_to_tab_line(heds)
resultstats2+=write_to_tab_line(heds)
resultstats+="\n"
resultstats2+="\n"
for organism in allfiles_resultsdics:
    resultstats+=""
    myinputpart=organism["inputpart"]
    myinputpart= myinputpart.split(SEP)
    heads2=""
    
    myfilestotreat=myinputpart[0]
    heads2+="myfilestotreat"+SEP
    resultstats2+=myfilestotreat+SEP
    myfiletype=str(myinputpart[1])
    heads2+="myfiletype"+SEP
    resultstats2+=str(myfiletype)+SEP
    myfilename=myinputpart[2]
    heads2+="myfilename"+SEP
    resultstats2+=str(myfilename)+SEP
    mynumberofrecords_inside=myinputpart[3]
    print (mynumberofrecords_inside)
    heads2+="numcrecords"+SEP
    resultstats2+=str(mynumberofrecords_inside)+SEP
    myresultsdirectory=myinputpart[4]
    heads2+="resultsdir"+SEP
    resultstats2+=str(myresultsdirectory)+SEP
    gendirectory=myinputpart[5]
    heads2+="gendir"+SEP
    resultstats2+=str(gendirectory)+SEP
    mynumberofgbfilesinfolder=myinputpart[6]
    heads2+="numgbfilesinfolder"+SEP
    resultstats2+=str(mynumberofgbfilesinfolder)+SEP
    mytotalfiles=organism["totalfiles"]
    heads2+="totalfiles"+SEP
    resultstats2+=str(mytotalfiles)+SEP
    mymineval_id_searches=organism["mineval_id_searches"]
    heads2+="mineval"+SEP
    resultstats2+=str(mymineval_id_searches)+SEP
    mytotalnt=organism["totalnt"]
    heads2+="totalnt"+SEP
    resultstats2+=str(mytotalnt)+SEP
    mytotasummaryfile=organism["summaryfile"]
    heads2+="summaryfile"+SEP
    resultstats2+=str(mytotasummaryfile)+SEP
    mymineval_first_search=organism["mineval_first_search"]
    heads2+="evalfirstsearch"+SEP
    resultstats2+=str(mymineval_first_search)+SEP
    

    print (organism.keys())
    for k in organism.keys():
        #heads2+="key in org"+SEP
        #heads2+="orgnanims k"+SEP

        resultstats+="key:"+str(k)+" "+str(organism[k])+"\t"
        #resultstats2+="key:"+str(k)+" "+str(k)+"\t"
        #resultstats2+="key:"+str(k)+" "+str(organism[k])+"\t"


    #filename=organism["inputpart"].split(SEP)[2]

    #: 'C:\\Users\\Eris\\Documents\\autothinktestfolder\\outputdir\\/C:\\Users\\Eris\\Documents\\autothinktestfolder\\frankiatestgenomes\\C:\\Users\\Eris\\Documents\\autothinktestfolder\\genomesummary.csv'
    #statscsvf=results_directory+"/"+myfilename+args.genomesummaryfile
    statscsvf=results_directory+"/"+myfilestotreat+"/"+myfilestotreat+"_summary.csv"
    hugsplit=helperobject.get_lines_from_file(statscsvf,False,"\n")
    myhugtabs=hugsplit[1]
    print (hugsplit)
    heads =hugsplit[0]
    for t,u in zip(heads,hugsplit[1:]):
        heads2+=t+SEP
        #resultstats+=myhugtabs+SEP
        resultstats+="hugsplit:"+str(t)+" "+myhugtabs+SEP
        resultstats2+=str(t)+" "+myhugtabs+SEP
        #hugtabs=hugsplit[0].split(SEP)

    #if first:
    #    for k in organism.keys():
    #        resultstats+=str(k)+"\t"
    #    resultstats+="\n"
    #first=False
    resultstats+="\n"
    resultstats2+="\n"
resultcsvline+=resultstats2

result2csvpath=os.path.join(results_directory,"result2.csv")
resultcsv_fh=open(result2csvpath,"w")
resultcsv_fh.write(resultcsvline)
resultcsv_fh.close()
