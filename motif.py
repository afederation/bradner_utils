#testing the motif class

import sys
sys.path.append('/ark/home/cl512/src/pipeline/')
from utils import *
import math
import random


def parseMotifFile(motifFile):
    '''
    parses a goofy fasta like file for a pwm
    where the header w/ the motif name begins with >
    and each line represents the nucleotide frequencies for a motif
    '''
    motifDict = {}
    motifTable = open(motifFile,'r')

    lines = motifTable.readlines()
    lines = [line.rstrip() for line in lines]
    
    for i in range(0,len(lines),5):
    #for i in range(10,15,5):
        header = lines[i]
        motifName = header[1:]
        print('working on %s' % (motifName))
        motifLines = lines[(i+1):(i+5)]
        #print(motifLines)
        #first character is important so you need to sort
        characterOrder = order([line[0] for line in motifLines])
        
        #now split the motifs into something comprehensible
        #gotta grab shit within the brackets
        motifLines = [re.findall('\[.*\]',line)[0][1:-1] for line in motifLines]
        motifLines = [line.split(' ') for line in motifLines]

        #now get rid of weird annoying spaces
        motifLines = [filter(lambda x: x != '',line) for line in motifLines]

        #now set up the pwm table
        pwm = []
        for j in range(len(motifLines[0])):
            pwm.append([])
        # print('here is the pwm')
        # print(pwm)
        # print('here are the motif lines')
        # print(motifLines)
        # print([len(foo) for foo in motifLines])
        #this should set up a list with L empty lists inside where L is the length of the motif

        #now fill in each list with the counts for each nucleotide
        #problem is we hae to iterate by nucleotide
        for j in characterOrder:
            line = motifLines[j]
            #print(line)
            for x in range(len(line)):
                try:
                    pwm[x].append(float(line[x]))
                except:
                    print('BOOOOO %s ' %(motifName))
                    
                    print(line)
                    print(x)
                    print(len(pwm))
                    print(pwm)


        motif = Motif(pwm,motifName)
        motifDict[motifName] = motif

    return motifDict




#class motif

#takes in a pwm matrix in the form of a table where each row is a position
#with 4 values corresponding to frequencies for ACGT

#e.g. [[.5,.2,.1,.2],[.6,.1,.2,.1],.... L] where L is length of motif
class Motif:

    def __init__(self,pwm,name,background =None,epsilon = 0.01):

        #first make sure the pwm is coerced into a floats
        for i in range(len(pwm)):
            pwm[i] = [float(x) for x in pwm[i]]
            

        #background should be a list or dictionry containing frequencies of ACGT
        if background:
            
            if type(background) == dict:
                background = [background['A'],background['C'],background['G'],background['T']]
                background = [float(x) for x in background]
        else:
            #if no background is given
            #need to figure out background frequencies from the motif itself
            #this is stupid
            background =[0.25,0.25,0.25,0.25]
            # for i in range(4):
            #     background.append(sum([line[i] for line in pwm]))
            
            # background = [x/sum(background) for x in background]

            
            
        pwmDict = {}
        nucleotides = ['A','C','G','T']
        for i in range(len(pwm)):

            pwmDict[i] = defaultdict(float)

            for j in range(4):
                #accounts for frequency
                #and takes the log likelihood
                pwmDict[i][nucleotides[j]] = math.log(pwm[i][j]/sum(pwm[i])/background[j] + epsilon)
        
        #HERE ARE THE VARIABLES THAT ARE SET AFTER INIT
        self._epsilon = epsilon
        self._background = background
        self._pwm = pwm                
        self._name = name
        self._pwmDict = pwmDict
        self._length = len(pwmDict)
        self._corePositions = self.findCorePositions()
        self._maxScore = self.maxScore(False)
        self._maxCoreScore = self.maxScore(True)
        self._bestMotif = self.bestMotif()


        coreMotif = ''
        for i in range(self._length):
            if self._corePositions.count(i) == 1:
                coreMotif+=self._bestMotif[i]
            else:
                coreMotif+='N'
        self._coreMotif = coreMotif

    def changeEpsilon(self,epsilon):

        '''
        allows user to reset epsilon correction
        '''
        self.__init__(self._pwm,self._name,self._background,epsilon)
        
    def resetBackground(self,background):

        '''
        allows user to reset backround and recalculate the log likelihood
        '''
        self.__init__(self._pwm,self._name,background,self._epsilon)


    def findCorePositions(self):

        '''
        finds the top 5 most conserved positions. if motif is less than 5, takes all positions
        '''
        coreLength = 6

        if self._length <= coreLength:
            return range(self._length)
        
        else:
            maxList = []
            nucleotides = ['A','C','G','T']
            for i in range(len(self._pwmDict)):
                maxList.append(max([self._pwmDict[i][x] for x in nucleotides]))

            maxOrder = order(maxList,decreasing=True)
            return maxOrder[0:coreLength]
                    
    def maxScore(self,core=True):

        '''
        figures out the max possible score
        '''
        score = 0
        nucleotides = ['A','C','G','T']
        if core:
            for i in self._corePositions:
                score+= max([self._pwmDict[i][x] for x in nucleotides])
        else:
            for i in range(len(self._pwmDict)):
                score+= max([self._pwmDict[i][x] for x in nucleotides])
            
        return score

    def bestMotif(self):
        
        '''
        returns the best motif
        '''

        bestMotif = ''
        nucleotides = ['A','C','G','T']
        for i in range(len(self._pwmDict)):
            positionScores= [self._pwmDict[i][x] for x in nucleotides]
            positionMax = max(positionScores)
            bestMotif+= nucleotides[positionScores.index(positionMax)]

        return bestMotif

    def scoreSequence(self,sequence,core=True,ratio = True):

        '''
        scores a sequence of same length as motif
        spits a massive error if the sequence ain't the same length
        '''

        sequence = upper(sequence)
        if len(sequence) != len(self._pwmDict):
            print('ERROR:SEQUENCE NOT SAME LENGTH AS MOTIF!!!')
            return
    
        score = 0.0
        if core:
            for i in self._corePositions:
                score+= self._pwmDict[i][sequence[i]]
        else:
            for i in range(len(sequence)):

                score+= self._pwmDict[i][sequence[i]]

        if ratio:
            if core:
                return score/self._maxCoreScore
            else:
                return score/self._maxScore
        else:
            return score
    

    def scoreSequenceTRAP(self, sequence, core=True, Lambda=0.7, ):
        
        sequence = upper(sequence)
        if len(sequence) != len(self._pwmDict):
            print('ERROR:SEQUENCE NOT SAME LENGTH AS MOTIF!!!')
            return

        width = len(self._pwmDict)
        R0 = 0.6*width - 6




    def calculateBackground(self,sequence):
        '''
        given a list of sequences or a single sequence
        recalculate the background
        '''

        total = 0
        backgroundDict = defaultdict(int)
        if type(sequence) == str:
            sequence = upper(sequence)
            for i in sequence:
                backgroundDict[i]+=1
                total+=1
        if type(sequence) == list:
            for line in sequence:
                if type(line) == list:
                    line = line[0]
                if line[0] == '>':
                    continue
                else:
                    #print(line)
                    line = upper(line)
                    for i in line:
                        backgroundDict[i] +=1
                        total+=1

        nucleotides = ['A','C','G','T']
        background = []
        for x in nucleotides:
            background.append(backgroundDict[x]/float(total))
        print('new background model is')
        print(background)
        self.resetBackground(background)
        
    def scoreFastaLine(self,fastaLine,cutoff = 0.8,core=True,localBackground=False):
        '''
        scores a fasta line and returns all instances of the motif above a cutoff
        as a dictionary keyed by position
        with the exact sequence as the entry
        '''
        fastaLine = upper(fastaLine)
        if localBackground:
            self.calculateBackground(fastaLine)
        
        matchDict = {}
        for i in range(0,len(fastaLine)-self._length+1,1):
            sequence = fastaLine[i:(i+self._length)]
            score = self.scoreSequence(sequence)
            if score > cutoff:
                matchDict[i] = sequence

        return matchDict

    def countFasta(self,fasta,cutoff =0.8,core=True,backgroundType = 'global',returnData = 'location'):

        '''
        count occurences of the motif for each fasta line.
        returns a list with the locations
        can use default background (none), a global background (global), or a local background(local)
        '''
        backgroundTypeList = ['global','local','none','default']
        backgroundType = lower(backgroundType)
        if backgroundTypeList.count(backgroundType) == 0:
            print('ERROR, background model not in options.  please use global, local, or none')
            return

        returnDataList = ['location','sequence','score']
        returnData = lower(returnData)
        if returnDataList.count(returnData) == 0:
            print('ERROR: returnData must be either location, sequence, or score')
            return

        #set the background
        if backgroundType == 'global':
            self.calculateBackground(fasta)
            localBackground = False
        elif backgroundType == 'local':
            localBackground = True
        else:
            localBackground = False

        #now go through each line of the fasta
        #output as a new table reusing the header lines in the 
        if type(fasta) == str:
            fasta = parseTable(fasta,'\t')

        motifTable = []
        ticker = 0
        for line in fasta:
            if ticker%1000 == 0:
                print(ticker/2)
            ticker+=1
            if type(line) == list:
                line = line[0]
            if line[0] == '>':
                motifTable.append([line])
                try:
                    start = int(line.split('_')[-3])
                    stop = int(line.split('_')[-2])
                except:
                    print(line)
                    exit()
                
            else:
                senseMatchDict = self.scoreFastaLine(line,cutoff,core,localBackground)
                antiMatchDict = self.scoreFastaLine(revComp(line),cutoff,core,localBackground)
                #matches = [int(x) for x in matchDict.keys()]
                #matches.sort()
                #matches = [str(x) for x in matches]
                if returnData == 'score':
                    matches = [motif.scoreSequence(senseMatchDict[x]) for x in senseMatchDict.keys()]
                    matches += [motif.scoreSequence(antiMatchDict[x]) for x in antiMatchDict.keys()]
                elif returnData == 'sequence':
                    matches = [senseMatchDict[x] for x in senseMatchDict.keys()]
                    matches +=[antiMatchDict[x] for x in antiMatchDict.keys()]
                else:
                    matches = [start + int(x) for x in senseMatchDict.keys()]
                    matches += [stop - int(x) for x in antiMatchDict.keys()]
                    matches.sort()
                matches = [str(x) for x in matches]
                motifTable.append([join(matches,',')])
        return motifTable

    def informationContent(self):

        total = 0
        for line in self._pwm:
            lineTotal = 0
            for base in line:
                if base != 0:
                    IC =  (base)*math.log(base)
                else:
                    IC = 0
                lineTotal += IC

            total += 2 + lineTotal

        return total

    def scoreFASTAline_TRAP(line,self,fastaLine,cutoff = 0.8,core=True,localBackground=False):
        '''
        scores a fasta line and returns all instances of the motif above a cutoff
        as a dictionary keyed by position
        with the exact sequence as the entry
        '''
        fastaLine = upper(fastaLine)
        if localBackground:
            self.calculateBackground(fastaLine)

        matchDict = {}
        for i in range(0,len(fastaLine)-self._length+1,1):
            sequence = fastaLine[i:(i+self._length)]
            score = self.scoreSequenceTRAP(sequence)
            if score > cutoff:
                matchDict[i] = sequence

        return matchDict



def mapMotifGFF(gffFile,genome,motifFile,cutoff=0.75,output='',scramble=False):

    '''
    counts all the motifs in each region
    '''
    if lower(genome) == 'mm9':
        genomeDirectory = '/nfs/genomes/mm9/fasta/'

    if lower(genome) == 'hg18':
        genomeDirectory = '/nfs/genomes/human_gp_mar_06/fasta/'

    gff = parseTable(gffFile,'\t')
    print('getting fasta from %s using genome %s from genome directory %s' % (gffFile,genome,genomeDirectory))
    fasta= gffToFasta(genome,genomeDirectory,gff)


    if scramble:
        scrambledFasta =[]
        ticker =0
        
        for line in fasta:
            
            
            if line[0] == '>':
                scrambledFasta.append(line)
            
            else:
                dinuc = range(0,len(line),2)
                random.shuffle(dinuc)
                newString = ''
                for i in dinuc:
                    newString+=line[i:(i+2)]
                scrambledFasta.append(newString)
            ticker+=1

        fasta = scrambledFasta

        
                
    #fasta = fasta[0:100]
    motifDict = parseMotifFile(motifFile)

    header= ['CHROM','NAME','','START','STOP','','SENSE','','NAME']
    motifList = motifDict.keys()
    motifList.sort()
    header+= [upper(x) for x in motifList]

    newGFF = [header]
    for line in gff:
        newGFF.append(line)

    print('The list of motifs in this analysis are')
    print(motifList)
    for x in motifList:
        print('Finding %s motifs in fasta' % (x))
        motif = motifDict[x]
        motifTable = motif.countFasta(fasta,cutoff,True,'global','location')
        for i in range(1,len(motifTable),2):
            if len(motifTable[i][0]) == 0:
                motifCount = 0
            else:
                motifCount = len(motifTable[i][0].split(','))
            gffIndex = i/2+1
            newGFF[gffIndex].append(motifCount)

    if len(output) == 0:
        return newGFF
    else:
        unParseTable(newGFF,output,'\t')

            
    
#==========================================================================
#================================RANK MOTIF================================
#==========================================================================





def rankMotif(summitBed,genomeDirectory,output,window=100,motifRegex=''):


    '''
    takes a summit bed and ranks eboxes by height from a sequence in a window around the bed

    '''
    
    # if len(motifRegex) == 0:
    #     coreMotif = motif._coreMotif
    #     coreMotif = coreMotif.replace('N','.')
    #     while coreMotif[0] == '.':
    #         coreMotif = coreMotif[1:]
    #     while coreMotif[-1] == '.':
    #         coreMotif = coreMotif[:-1]
    #     motifRegex = coreMotif
    #     print('No motif regular expression provided. Using %s derived from position weight matrix' % (motifRegex))

    print('Using user provided motif regular expression: %s' % (motifRegex))

    if type(summitBed) == str:
        summitBed = parseTable(summitBed,'\t')
    
    motifDict = defaultdict(list)
    ticker= 0
    occurrences = 0
    for line in summitBed:
        if ticker % 1000 == 0:
            print(ticker)
        ticker+=1

        chrom = line[0]
        peakName = line[3]
        sense = '.'


        start = int(line[1])-window
        end = int(line[1])+window
        height = float(line[4])

        sequenceLine = fetchSeq(genomeDirectory,chrom,start,end,True)
        
        motifVector = []
        matches = re.finditer(motifRegex,upper(sequenceLine))
        if matches:
            
            for match in matches:
                motifVector.append(match.group())
                occurrences +=1
        #count only 1 of each motif type per line
        motifVector = uniquify(motifVector)
        for nmer in motifVector:

            motifDict[nmer].append(height)


    motifTable =[]
    motifTableOrdered =[['MOTIF','OCCURENCES','AVG_HEIGHT']]
    for nmer in motifDict.keys():
        newLine = [nmer,len(motifDict[nmer]),mean(motifDict[nmer])]
        motifTable.append(newLine)


    occurenceOrder = order([line[2] for line in motifTable],decreasing=True)
    
    for x in occurenceOrder:
        motifTableOrdered.append(motifTable[x])
    print(motifTableOrdered)
    print('motif occured %s times' % (occurrences))
    unParseTable(motifTableOrdered,output,'\t')

#==========================================================================
#===================RANKING MOTIFS FROM GFF================================
#==========================================================================

def rankMotifGFF(gffFile,species,genomeDirectory,motifRankFile,motifRegex,output):

    '''
    finds all of the motifs in a fasta and writes a table of how often they occur and the avg. score
    for that particular sequence
    '''

    
    fasta = gffToFasta(species,genomeDirectory,gffFile)
    motifDict = defaultdict(float)

    motifTable = parseTable(motifRankFile,'\t')


    for line in motifTable[1:]:
        motifDict[line[0]] = float(line[2])
    print(motifDict)

    occurenceDict = defaultdict(int)
    print(motifRegex)
    print(len(fasta))
    for line in fasta:

        if line[0] == '>':
            continue


        matches = re.finditer(motifRegex,upper(line))
        if matches:
            for match in matches:
                #print(match.group())
                occurenceDict[match.group()]+=1


    strengthOrder = order([line[2] for line in motifTable[1:]],decreasing=True)
    #print(strengthOrder)
    motifList = [motifTable[x+1][0] for x in strengthOrder]

    print(occurenceDict)
    occurenceTable = [['MOTIF','OCCURENCES','MOTIF_SCORE']]

    for nmer in motifList:

        occurenceTable.append([nmer,occurenceDict[nmer],motifDict[nmer]])


    unParseTable(occurenceTable,output,'\t')




# #testing in MES

# #gffFile = '/nfs/young_ata4/pipeline/mes/gff/MM9_ENRICHED_MES_OSN_STRICT_12KB_STITCHED_WITH_TSS_-0_+0.gff'
# #genome = 'mm9'
# motifFile = '/nfs/young_ata4/pipeline/mes/motifs/tfMotifs.txt'


# #foo = mapMotifGFF(gffFile,genome,motifFile,output='/nfs/young_ata4/pipeline/mes/tables/MM9_ENRICHED_MES_OSN_STRICT_12KB_STITCHED_WITH_TSS_-0_+0_MOTIFS.gff')



# #testing motif ranking




# summitBed = '/nfs/young_ata4/DAO/superEnhancers/mES/motif_2E/MES_OCT4_OSNSTRICT_summits.bed'

# genomeDirectory = '/nfs/genomes/mm9/fasta/'

# output = '/nfs/young_ata4/pipeline/mes/motifs/OSN_Rank_test.txt'

# #motifRegex = 'TT.T.ATG.A.A'
# #motifRegex = 'TT.T.ATG.'
# #motifRegex = 'T.ATTT.CAT.'
# motifRegex = 'ATG..AA.'


# #rankMotif(summitBed,genomeDirectory,output,500,motifRegex)    


# gffFile = '/nfs/young_ata4/DAO/superEnhancers/mES/motif_2E/mES_superEnhancers_only.gff'
# species = 'mm9'
# motifRankFile = '/nfs/young_ata4/pipeline/mes/motifs/OSN_Rank_test.txt'

# output = '/nfs/young_ata4/pipeline/mes/motifs/OSN_RANK_SUPER.txt'


# #typicals
# rankMotifGFF(gffFile,species,genomeDirectory,motifRankFile,motifRegex,output)


# gffFile = '/nfs/young_ata4/DAO/superEnhancers/mES/motif_2E/mES_typicalEnhancers_only.gff'
# species = 'mm9'
# motifRankFile = '/nfs/young_ata4/pipeline/mes/motifs/OSN_Rank_test.txt'

# output = '/nfs/young_ata4/pipeline/mes/motifs/OSN_RANK_TYPICAL.txt'

# rankMotifGFF(gffFile,species,genomeDirectory,motifRankFile,motifRegex,output)



#Testing in MM1S

# gffFile = '/nfs/young_ata4/pipeline/mm1s_new/gff/HG18_ENRICHED_MM1S_MED1_DMSO_2_12KB_STITCHED_NO_TSS_-0_+0.gff'
# genome = 'hg18'
# motifFile = '/nfs/young_ata4/pipeline/mm1s_new/motifs/MM1S_MOTIFS.txt'


# foo = mapMotifGFF(gffFile,genome,motifFile,output='/nfs/young_ata4/pipeline/mm1s_new/tables/HG18_ENRICHED_MM1S_MED1_DMSO_2_12KB_STITCHED_NO_TSS_-0_+0_MOTIFS.txt')





#print(foo)
# motifDict = parseMotifFile(motifFile)

# print(motifDict.keys())



# motif = motifDict['Pou5f1']
# #motif.resetBackground([0.25,0.25,0.25,0.25])
# print(motif._pwmDict)
# foo = motif._pwmDict

# nucleotides = ['A','C','G','T']
# maxList = []
# for i in range(len(foo)):
#     print(i)
#     maxList.append(max([foo[i][x] for x in nucleotides]))
# print(maxList)

# print(order(maxList,decreasing=True))



# print(motif._maxScore)
# print(motif._bestMotif)
# print(motif._maxCoreScore)
# print(motif._coreMotif)

# print(motif.scoreSequence('CAGTGTCATGCAAAT'))
# print(motif.scoreSequence('CATTGTCATGCAAAT'))
# print(motif.scoreSequence('CTTTGTCATGCAAAT'))

# print(motif.scoreFastaLine('CTTTGTCATGCAAATAAAAAAAAACTTTGTCATGCAAATAAAAAAAAACTTTGTCATGCAAATAAAAAACATTGTCATGCACAT',0.8,True,False))
# print(motif.scoreFastaLine('CTTTGTCATGCAAATAAAAAAAAACTTTGTCATGCAAATAAAAAAAAACTTTGTCATGCAAATAAAAAACATTGTCATGCACAT',0.8,True,True))

# fasta = parseTable('/nfs/young_ata4/pipeline/mes/fasta/MM9_ENRICHED_MES_OSN_STRICT_12KB_STITCHED_WITH_TSS_-0_+0.fasta','\t')

# fasta = fasta[0:100]



# print(motif.countFasta(fasta,0.75))

# motif = motifDict['Myc']
# print(motif.countFasta(fasta,0.75,True,'global','sequence'))
# #print(motif.maxMotif())
# #print(motif.bestMotif())
