import sys
barf=sys.argv[1]
fastq=sys.argv[2]
linker=sys.argv[3].upper()
nwise=int(sys.argv[4])

def getBarcodes(barf):
 file=open(barf,"r")
 barcodes=[]
 for ln in file:
  barcodes.append(ln.strip("\r\n"))
 return barcodes
barcodes=getBarcodes(barf)

def countCombos(barcodes,combo,dict):
 key=""
 for i in combo:
  if i in barcodes:
   key+=str(barcodes.index(i))+"_"
 key=key.strip("_")
 if key != "" and key not in dict:
  dict[key]=[1]
 elif key != "" and key in dict:
  dict[key].append(dict[key][-1]+1)
 return dict

def addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict):
 lenBarcode=len(barcodes[0])
 seq=seq.split("\n")[0]
 if "N" not in seq:
  seq=seq[len(adaptor):].split(linker)
  if len(seq) == nwise+1:
   combo=[seq[0][-8:]]
   for j in range (1,nwise):
    if len(seq[j]) == lenBarcode:
     combo.append(seq[j])
   if len(combo) == nwise:
    seqs.append(combo)
    dict=countCombos(barcodes,combo,dict)
 return (seqs,dict)  

def readFastq(fastq,linker,nwise,barcodes):
 adaptor=fastq.split("_")[1]
 file=open(fastq,"r") 
 seq,seqs,totalEntries,dict="",[],0,{}
 for ln in file:
  if ln.startswith("@"):
   totalEntries+=1
   if seq != "":
    R=addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict)
    seqs,dict=R[0],R[1]
   seq=""
  elif "+" not in ln:
   seq+=ln
 R=addSeqtoSeqs(seq,seqs,adaptor,linker,nwise,barcodes,dict)
 seqs,dict=R[0],R[1]
 return (seqs,totalEntries,dict)
Fastq=readFastq(fastq,linker,nwise,barcodes)
seqs,totalEntries,dict=Fastq[0],Fastq[1],Fastq[2]
#for i in seqs:
# print(i)
for k in dict:
 print(k,dict[k][-1])
print(len(dict))
