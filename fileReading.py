"""
fileReading.py
---------------------------------------------------------------------------------------------------------------------------
This file controls all the file reading/interface with data.  It handles both reading the information from the primer file, 
 and the alligned data points.
"""
from numpy import asarray
import os

###  TODO: add data saving/reading later capabilities (for faster inference on repeated files)

def primerRead(primerFil):
    with open(primerFil) as f:
        lines = [line.rstrip('\n') for line in f]
        
    names = []; loc = []; seq = []; strand = []
    for lin in lines:
        temp = lin.split('\t')
        names.append(temp[0])
        strand.append(temp[1])
        loc.append(temp[3])
        seq.append(temp[4])    
    loc = asarray(list(map(int, loc)))
    pDict = dict({'names': names, 'loc': loc, 'seq': seq, 'strand': strand})
    return(pDict)

#convert the files of reads into objecs which can be used for the interactive graphing. 
def processFiles(dirName, pDis):
    #each of the different barcodes will store its results in a dict, which will in turn be stored in this array
    barcodeDict = [dict() for x in range(0)]
    for file in os.listdir(dirName):
        if file.endswith('align'):
            #initialize the dict
            nDict = dict()
            for k in range(len(pDis)):
                nDict['seq' + str(k)] = []
            
            #read the file, put it in dict
            with open(dirName + '/' + file) as f:
                lines = [line.rstrip('\n') for line in f]
            for lin in lines:
                temp = lin.split('\t')
                loc = int(temp[3])
                if temp[1] == '-':   #minus strand has to be length adjusted
                    loc += len(temp[4])
                pos = (abs(pDis - loc)).argmin()
                seq = nDict['seq' + str(pos)]
                seq.append(str(temp[4]))
                nDict['seq' + str(pos)] = seq
            barcodeDict.append(nDict)          
    return(barcodeDict)

def barcodeNames(barcName, dirName):
    if barcName == '':
        barcFile = dirName[:dirName.rfind('\\')] + '\\barcodeInfo.txt'
    else:
        barcFile = dirName + '\\' + barcName
    bNames = []
    try:
        with open(barcFile) as f:
            lines = [line.rstrip('\n') for line in f]
        for l in lines:
            bNames.append(l[l.find('\t'):])
    except:
        print('Please either pass in barcodeInfo.txt filename, or place such a file in this directory: '
              + dirName[:dirName.rfind('\\')])
    return(bNames)

def barcodeInit(barcFile):
    """
    barcodeInit will take a barcode file and create a dictionary to store reads for specific barcodes
    :param barcFile: the file containing information about the specific barcodes
    :return bDict: a dictionary from barcode names to an empty array
    """
    bNames = []
    try:
        with open(barcFile) as f:
            lines = [line.rstrip('\n') for line in f]
        for l in lines:
            bNames.append(l[:l.find('\t')])
    except:
        print('Please pass in a valid barcode file to use in alignment')
    bDict = dict({})
    for b in bNames:
        bDict[b] = []
    return(bDict)

def sortStrip(dirName, barcodes, saveFolder, fileEnd = 'sequence.fastq', newName = 'sorted.stripped'):
    """
    sortStrip takes raw files from the sequencer and begins to prepare them for aligning by cutting polyA tails.
        It then saves these sorted and stripped files into a new folder (provided by user)
    :param dirName: the name of the directory with all the files to be processed
    :param barcodes: a dictionary for all desired barcode reads to be used in sorting reads from the file
    :param saveFolder: the folder to save the processed files to
    :param fileEnd: the ending of files to process.  Defaults to 'sequence.fastq'
    :param newName: the new ending to add to files which have been processed. Default is 'sorted.stripped'
    """
    if not os.path.exists(saveFolder): #make sure the folder to save new files into exists
            os.makedirs(saveFolder)

    for file in os.listdir(dirName):
        if file.endswith(fileEnd):
            with open(dirName + '/'+ file) as f:
                lines = [line.rstrip('\n') for line in f]
            for k in range(int(len(lines)/4)):
                start = '#'; end = '/'; l = lines[4*k] #mark where the barcode name starts and ends
                bName = l[l.find(start) + len(start):l.rfind(end)]
                if bName in barcodes.keys():
                    seq = lines[4*k+1]
                    qual = lines[4*k + 3]
                    while seq.endswith('A'): #trim off the polyA tails
                        seq = seq[:-1]
                        qual = qual[:-1]
                    barcodes[bName].append(l[:l.rfind(end)])
                    barcodes[bName].append(seq)
                    barcodes[bName].append('+')
                    barcodes[bName].append(qual)
                    
    for key in barcodes.keys():
        fname = saveFolder + key + newName
        with open(fname, "w+") as f:
            for s in barcodes[key]:
                f.write(s + '\n')
        f.close()

def runBokeh(folder, saveFile, fileEnd = 'sorted.stripped', missFile = ''):
    """
    runBokeh will run the allignment tool on all files within the provided folder.  It will then save the alligned
    sequences to the saveFolder (which is the default).  Also possible to save information regarding mismatch sequences
    
    :param folder: folder containing all reads from the sequencer which have had PolyA tails removed. 
    :param saveFile: where to save alligned reads. 
    :param fileEnd: the expected ending of files to align.  Defaults to the one created in this pipeline (sorted.stripped)
    :param missFile: if provided the program will also save information about reads that couldn't be alligned
    """    
    if missFile == '':
        os.system("python ")
    else:
        os.system("python ")