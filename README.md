# new BCanalyzer
Before running the program, please check whether regex (not the traditional re module) is installed. To install regex:
```
pip install regex
```
To run the program:
```
python newBCanalyzer.py -b <list of barcodes> -f <fastq NGS file> -l <linker sequence> [-n nwise]
#the argument -n is optional, and the default value is 2
```
- 12/30/19: back from Europe! Annotating the code! Adding in the batch processing function!
