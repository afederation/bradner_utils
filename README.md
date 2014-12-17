BRADNER UTILS
================

This is a set of tools used extensively in the Bradner Lab and Young Lab for the analysis of seq data. Most functions are based around the Locus object that represents a location in the genome

Dependencies
-----------------

Python 2.7.6

samtools (http://samtools.sourceforge.net/)

bamliquidator (https://github.com/BradnerLab/pipeline/wiki/bamliquidator)

Functions
-----------------

###1. Input/Output and file handling functions

def open(file,mode='r')

- replaces open with a version that can handle gzipped files

def parseTable(fn, sep, header = False, excel = False)

- opens standard delimited files

def unParseTable(table, output, sep) 

- writes standard delimited files, opposite of parseTable

def gffToBed(gff,output= '')

- converts standard UCSC gff format files to UCSC bed format files

def formatFolder(folderName,create=False)

- checks for the presence of any folder and makes it if create = True

###2. Gene annotation functions
def makeStartDict(annotFile,geneList = [])

- takes a standard UCSC refseq table and creates a dictionary keyed by refseq ID with info about each transcript

def getTSSs(geneList,refseqTable,refseqDict)

- returns the TSS location of any gene

def importRefseq(refseqFile, returnMultiples = False)

- imports a standard UCSC refseq annotation file into a dictionary

def makeGenes(annotFile,geneList=[],asDict = False) 

- takes a UCSC refseq annotation file and a gene list and makes a list or dictionary of Gene class objects

def makeTranscriptCollection(annotFile,upSearch,downSearch,window = 500,geneList = []) 

- takes a UCSC refseq annotation file and makes a LocusCollection where each locus is a full transcript

def importBoundRegion(boundRegionFile,name) 

- imports a bound region file (a standard bed or macs output bed)

###3. Locus class
class Locus(chr,start,end,sense,ID) 

- standard locus class for tracking genomic loci

class LocusCollection(lociList,windowSize=500) 

- a collection of locus objects used for querying large sets of loci. The window size paramater is outdated and not needed for this class to function

###4. Gene class
class Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds,commonName='') 

- gene class object that contains all annotation information about a given transcript

###5. Locus functions
def locusCollectionToGFF(locusCollection) 

- turns a locus collection into a gff

def gffToLocusCollection(gff,window =500)

- turns a gff into a locus collection (reverse of gff)

def makeTSSLocus(gene,startDict,upstream,downstream)

- from a start dict makes a locus surrounding the tss

def makeSearchLocus(locus,upSearch,downSearch) 

- takes an existing locus and makes a larger flanking locus

def makeSECollection(enhancerFile,name,top=0)

- returns a locus collection from a super enhaner table from ROSE

###6. Bam class

class Bam(bamFile) 

- a class for handling and manipulating bam objects.  requires samtools

###7. Misc. functions
def uniquify(seq, idfun=None)  

- makes a list unique

def order(x, NoneIsLast = True, decreasing = False) 

- returns the ascending or descending order of a list

###8 Sequence functions
def fetchSeq(directory,chrom,start,end,UCSC=False,lineBreaks=True,header = True) 

- grabs sequence from a region

def gffToFasta(species,directory,gff,UCSC = True)

- converts a gff to a fasta

def revComp(seq,rev = True, RNA=False)

###9 sortCoordsBED.py

- a helper function to sort the start and end coordinates in a BED file