import argparse, regex
parser = argparse.ArgumentParser()
parser.add_argument("-b", required=True, type=str, help="barcode list")
parser.add_argument("-f", required=True, type=str, help="fastq file")
parser.add_argument("-l", required=True, type=str, help="linker sequence")
parser.add_argument("-n", default=2, type=int, help="nwise")
args=parser.parse_args()
barf=args.b
fastq=args.f
linker=args.l.upper()
nwise=args.n
adaptor=fastq.split("_")[1]

#get all the barcodes used in the experiment
def getBarcodes(barf):
 file=open(barf,"r")
 barcodes=[]
 for ln in file:
  barcodes.append(ln.strip("\r\n"))
 return barcodes
barcodes=getBarcodes(barf)

#allowing one nt mismatch in the sample id/barcode in the sequence
def oneMismatch(id,seq):
 restring='('+id+'){e<=1}'
 m=regex.findall(restring,seq,overlapped=True)
 return m

def countCombos(barcodes,combo,dict,nwise,successC):
 key=""
 for i in combo: #for each barcode in the combo
  if i in barcodes: #exact match of the barcode
   key+=str(barcodes.index(i))+"_" #generate the key
  else: #if the barcode is not found in the barcode list
   for b in barcodes: #loop thorough all the barcodes
    m=oneMismatch(i,b) #to see whether there's a barcode with a 1nt mismatch
    if len(m) == 1: #the list can be > 1, which can be ambiguious, so only one 1nt-mistmach map is counted
     key+=str(barcodes.index(b))+"_" #genreate the key
 key=key.strip("_") #generate the key
 if key != "": #if the key is not empty
  skey=len(key.split("_")) #see whether the key contains all the barcodes in all dimensions
  if skey == nwise and key not in dict: # add key to successful count and the dictionary if the key is not in the dictionary
   dict[key]=[1]
   successC+=1
  elif skey == nwise and key in dict: #add the entry to the key and successful count if the key is in the dictionary already
   dict[key].append(dict[key][-1]+1)
   successC+=1
 return (dict,successC)

def chop(nwise,S):
 start=16
 end=start+18
 seq=[S[start:end]]
 for a in range (nwise):
  start=end-2
  end=start+18
  seq.append(S[start:end])
 truncate=[]
 for i in seq:
  m=oneMismatch(linker,i)
  for j in m:
   if len(j) == len(linker):
    truncate.append(j)
 if len(truncate) == 3:
  barcodes=[]
  for link in truncate:
   end=S.find(link)
   barcode=S[end-8:end]
   barcodes.append(barcode)
  return barcodes

def seqtoBarcode(s,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode):
 seq=s[len(adaptor):].split(linker) #looking for exact matches for all linkers
 if len(seq) == nwise+1:  #if all of the linkers were sequenced to be correct
  combo=[seq[0][-8:]] #add the 1st barcode to combo
  for j in range (1,nwise): #for each of the rest of the barcodes
   if len(seq[j]) == lenBarcode: #if after splitting the following barcodes are of the lengths of a barcode =
    combo.append(seq[j]) #add the barcode to the combo
   if len(combo) == nwise: #if there are all n # of barcodes in a combo
    row=adaptor #generate the row added to the sample speparator output
    for b in combo:
     row+=","+b
    sampsepout.write(row+"\r\n")
    seqs.append(combo) #add the combo to list of sequences
    COMB=countCombos(barcodes,combo,dict,nwise,successC)
    dict,successC=COMB[0],COMB[1] #dictionary contains all the keys and # of entries; all successful counts
 else:
  combo=chop(nwise,s[len(adaptor):])
  if combo != None:
   row=adaptor #generate the row added to the sample speparator output
   for b in combo:
    row+=","+b
   sampsepout.write(row+"\r\n")
   seqs.append(combo) #add the combo to list of sequences
   COMB=countCombos(barcodes,combo,dict,nwise,successC)
   dict,successC=COMB[0],COMB[1] #dictionary contains all the keys and # of entries; all successful counts
 return (dict,successC)

def addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC):
 lenBarcode=len(barcodes[0])
 seq=seq.split("\n")[0] #this keeps the 2nd line of the fastq entry
 if "N" not in seq:  #make sure that none of the nt sequenced is uncertain
  if seq.startswith(adaptor) == True:  #if the seqeunce starts with exact match with the sample id (adaptor here)
   runseqtoBarcode=seqtoBarcode(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode)
   dict,successC=runseqtoBarcode[0],runseqtoBarcode[1]
  else: #if the start of the sequence is not an exact match of the sample id (adaptor here)
   m=oneMismatch(adaptor,seq[:len(adaptor)]) #find whether a sequence with 1nt mismatch within the sample id region
   if len(m) > 0: #the max. of this is 1, which indicates that a sequence contains 1nt mismatch from the sample id region
    runseqtoBarcode=seqtoBarcode(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode)
    dict,successC=runseqtoBarcode[0],runseqtoBarcode[1]    
 return (seqs,dict,successC)  

#the running of the program starts here
def readFastq(fastq,linker,nwise,barcodes,adaptor):
 file=open(fastq,"r") #open the fastq input file
 sampsepfn="sample_"+adaptor+".csv" #open the sample separator output file
 sampsepout=open(sampsepfn,"w")
 row="Sample_ID"              #generate the header for the sample separator output
 for a in range (nwise):
  row+=",BC"+str(nwise-a)
 sampsepout.write(row+"\r\n")
 sampseprp="sample_"+adaptor+"_report.csv" #open the report file for the sample separator
 sampseprep=open(sampseprp,"w")
 seq,seqs,totalEntries,dict,successC="",[],0,{},0
 for ln in file:   
  if ln.startswith("@"): #each fastq entry starts with a "@"
   totalEntries+=1 
   if seq != "":   #this adds the seq to the seq list if it is not empty
    R=addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC)
    seqs,dict,successC=R[0],R[1],R[2]
   seq=""
  elif "+" not in ln:  #this add the 2nd (the sequence) and 4th (the quality scores) line of a fastq entry to seq
   seq+=ln
 R=addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC)
 seqs,dict=R[0],R[1]
 sampsepout.close()
 sampseprep.write("Total number of reads: ,"+str(totalEntries)+"\r\n")
 sampseprep.write("Successful reads: ,"+str(len(seqs))+"\r\n")
 sampseprep.write("Wrong reads: ,"+str(totalEntries-len(seqs))+"\r\n")
 sampseprep.write("Discarded reads: ,"+str(totalEntries-len(seqs))+"\r\n")
 sampseprep.close()
 return (seqs,dict,successC)
Fastq=readFastq(fastq,linker,nwise,barcodes,adaptor)
seqs,dict,successC=Fastq[0],Fastq[1],Fastq[2]

BCfn="BCcount_"+adaptor+".csv"
BCout=open(BCfn,"w")
BCout.write("Sample ID: ,"+adaptor+"\r\n")
BCout.write("Input count: ,"+str(len(seqs))+"\r\n")
BCout.write("Combinations: ,"+str(len(dict))+"\r\n")
BCout.write("Success Count: ,"+str(successC)+"\r\n")
BCout.write("Discarded Count: ,"+str(len(seqs)-successC)+"\r\n")
BCout.write(""+"\r\n")
BCout.write("Combinations"+"\r\n")
row=""
for i in range(nwise):
 row+="BC"+str(nwise-i)+","
row=row.strip(",")
BCout.write(row+"\r\n")
def outBCanalyzer(barcodes,dict,BCout):
 for key in dict:
  k=key.split("_")
  row=""
  for b in k:
   row+=barcodes[int(b)]+","
  row+=","+str(dict[key][-1])
  BCout.write(row+"\r\n")
 BCout.close()
outBCanalyzer(barcodes,dict,BCout)
