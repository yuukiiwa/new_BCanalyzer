import sys, regex
barf=sys.argv[1]
fastq=sys.argv[2]
linker=sys.argv[3].upper()
nwise=int(sys.argv[4])
adaptor=fastq.split("_")[1]

def getBarcodes(barf):
 file=open(barf,"r")
 barcodes=[]
 for ln in file:
  barcodes.append(ln.strip("\r\n"))
 return barcodes
barcodes=getBarcodes(barf)

def oneMismatch(id,seq):
 restring='('+id+'){e<=1}'
 m=regex.findall(restring,seq,overlapped=True)
 return m

def countCombos(barcodes,combo,dict,nwise,successC):
 key=""
 for i in combo:
  if i in barcodes:
   key+=str(barcodes.index(i))+"_"
 key=key.strip("_")
 if key != "":
  skey=len(key.split("_"))
  if skey == nwise and key not in dict:
   dict[key]=[1]
   successC+=1
  elif skey == nwise and key in dict:
   dict[key].append(dict[key][-1]+1)
   successC+=1
 return (dict,successC)

def seqtoBarcode(s,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode):
 seq=s[len(adaptor):].split(linker)
 if len(seq) == nwise+1:
  combo=[seq[0][-8:]]
  for j in range (1,nwise):
   if len(seq[j]) == lenBarcode:
    combo.append(seq[j])
   if len(combo) == nwise:
    row=adaptor
    for b in combo:
     row+=","+b
    sampsepout.write(row+"\r\n")
    seqs.append(combo)
    COMB=countCombos(barcodes,combo,dict,nwise,successC)
    dict,successC=COMB[0],COMB[1]
 else:
  print(oneMismatch(linker,s[len(adaptor):]))
 return (dict,successC)

def addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC):
 lenBarcode=len(barcodes[0])
 seq=seq.split("\n")[0]
 if "N" not in seq:
  if seq.startswith(adaptor) == True:
   runseqtoBarcode=seqtoBarcode(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode)
   dict,successC=runseqtoBarcode[0],runseqtoBarcode[1]
  else:
   m=oneMismatch(adaptor,seq[:len(adaptor)])
   if len(m) > 0:
    runseqtoBarcode=seqtoBarcode(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC,lenBarcode)
    dict,successC=runseqtoBarcode[0],runseqtoBarcode[1]    
 return (seqs,dict,successC)  

def readFastq(fastq,linker,nwise,barcodes,adaptor):
 file=open(fastq,"r") 
 sampsepfn="sample_"+adaptor+".csv"
 sampsepout=open(sampsepfn,"w")
 row="Sample_ID"
 for a in range (nwise):
  row+=",BC"+str(nwise-a)
 sampsepout.write(row+"\r\n")
 sampseprp="sample_"+adaptor+"_report.csv"
 sampseprep=open(sampseprp,"w")
 seq,seqs,totalEntries,dict,successC="",[],0,{},0
 for ln in file:
  if ln.startswith("@"):
   totalEntries+=1
   if seq != "":
    R=addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict,sampsepout,successC)
    seqs,dict,successC=R[0],R[1],R[2]
   seq=""
  elif "+" not in ln:
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
