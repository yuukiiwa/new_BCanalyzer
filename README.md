# new BCanalyzer
Before running the program, please check whether regex (not the traditional re module) is installed. To install regex:
```
pip install regex
```
To run the program:
```
python newBCanalyzer.py -b <list of barcodes> -d <directrory of the fastq(s)> [-l linker seq] [-n # of dimensions]
#the argument -l is optional, and the default is CAATTC
#the argument -n is optional, and the default value is 2
```
- 12/30/19: back from Europe! Annotating the code! Adding in the batch processing function!
