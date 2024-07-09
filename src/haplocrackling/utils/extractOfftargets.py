#!/usr/bin/python3.6
'''
Purpose:    identify all offtarget sites in the whole genome when considering variants

Input:      FASTA, or multi-FASTA, formatted file

Output:     one file with all the sites

To use:     extractOfftargets.py output-file  (input-files... | input-dir>) --vcf vcf-file

'''

import glob, multiprocessing, os, re, shutil, tempfile, heapq, argparse
from ast import literal_eval
from haplocrackling.Helpers import *
from haplocrackling.Paginator import Paginator

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"

def explodeMultiFastaFile(fpInput, fpOutputTempDir):
    newFilesPaths = []

    with open(fpInput, 'r') as fRead:
        fWrite = None

        for line in fRead:
            line = line.strip()
            
            # just found a new fasta segment. open a new file
            if line[0] == '>':
                fpTemp = tempfile.NamedTemporaryFile(
                    mode = 'w+', 
                    delete = False,
                    dir = fpOutputTempDir
                )
                
                # close the current file if necessary
                if fWrite is not None:
                    fWrite.write('\n')
                    fWrite.close()
                    
                # open a new one
                fWrite = open(fpTemp.name, 'w+')
                newFilesPaths.append(fpTemp.name)
                
                fWrite.write(line)
                fWrite.write('\n')

            if line[0] != '>':
                #fWrite.write(line.upper().strip())
                fWrite.write(line.strip())
        
        # close the last file if necessary
        if fWrite is not None:
            fWrite.close()
            
    return newFilesPaths

def processingNode(fpInputs, fpOutputTempDir, fpVcf):
    # Create a temporary file
    fpTemp = tempfile.NamedTemporaryFile(
        mode = 'w+', 
        delete = False,
        dir = fpOutputTempDir
    )
    
    with open(fpTemp.name, 'w+') as outFile:

        for fpInput in fpInputs:
        
            # key: FASTA header, value: sequence
            seqsByHeader = {}

            with open(fpInput, 'r') as inFile:
                header = fpInput
                
                for line in inFile:
                    if line[0] == '>':
                        header = line[1:].strip()
                        seqsByHeader[header] = []
                    else:
                        if header not in seqsByHeader:
                            # it could be a plain text file, without a header
                            seqsByHeader[header] = []
                        seqsByHeader[header].append(line.strip())#.upper())

            # For each FASTA sequence
            offtargets = [] 
            for header in seqsByHeader:
                seq = ''.join(seqsByHeader[header])
                
                # read VCF file, only keep the records for this sequence
                vcf = []
                with open(fpVcf,'r') as fVcf:
                    for line in fVcf:
                        if line[0] == "#": # comments are ignored
                            continue
                            
                        tempArray = line.rstrip().split("\t")

                        if tempArray[6] != "PASS" or tempArray[0] != header: 
                            # in our case, the filter should always be "PASS"
                            # and the record be relevant to our sequence
                            continue

                        vcf.append(tempArray)
                
                for strand, pattern, seqModifier in [
                    ['positive', pattern_forward_offsite, lambda x : x],
                    ['negative', pattern_reverse_offsite, lambda x : rc(x)]
                ]:
                    match_chr = re.findall(pattern, seq)

                    # Extract off-targets using the exact sequence
                    for i in range(0,len(match_chr)):
                        offtargets.append(
                            seqModifier(match_chr[i])[0:20]
                        )

                    # Extract off-targets considering VCF
                    for record in vcf:
                    
                        pos = literal_eval(record[1])-1 # with -1 because VCF POS is 1-based
                        ref = record[3]
                        alt = record[4]
                        
                        tempArray = alt.split(",")
                        for el in tempArray:
                            alt_seq = seq[pos-22:pos] + el + seq[pos+len(ref):pos+len(ref)+23]

                            match_chr = re.findall(pattern, alt_seq)
                            for i in range(0,len(match_chr)):
                                offtargets.append(
                                    seqModifier(match_chr[i])[0:20]
                                )
                    
            outFile.write(''.join(f'{offTarget}\n' for offTarget in offtargets))
            return len(offtargets)

# Node function that sorts a file for multiprocessing pool
def sortingNode(fileToSort, sortedTempDir):
    # Create a temporary file
    sortedFile = tempfile.NamedTemporaryFile(
        mode = 'w+', 
        delete = False,
        dir = sortedTempDir
    )
    # Sort input file and store in new output dir
    with open(fileToSort, 'r') as input:
        # Read 'page'
        page = input.readlines()
        # Sort Page
        page.sort()
        # Write sorted page to file
        sortedFile.writelines(page)
        # Close sorted file
        sortedFile.close()

def paginatedSort(filesToSort, fpOutput, mpPool, maxNumOpenFiles=400): 
    # Create temp file directory
    sortedTempDir = tempfile.TemporaryDirectory()
    printer(f'Created temp directory {sortedTempDir.name} for sorting')

    # Generate args for sorting 
    args = [
        (
            file,
            sortedTempDir.name 
        ) for file in filesToSort
    ]

    # Submit job to multiprocessing pool
    mpPool.starmap(
        sortingNode,
        args
    )
    
    printer('Sorting of each file completed')
    
    # Collect sorted files to merge
    sortedFiles = glob.glob(
        os.path.join(
            sortedTempDir.name,
            '*'
        )
    )
    
    if len(sortedFiles) == 1:
        shutil.move(sortedFiles[0], fpOutput)
    else:
        # Open all the sorted files to merge
        printer(f'Beginning to merge sorted files, {maxNumOpenFiles:,} at a time')
        while len(sortedFiles) > 1:
            # A file to write the merged sequences to
            mergedFile = tempfile.NamedTemporaryFile(delete = False)
            
            # Select the files to merge
            while True:
                try:
                    sortedFilesPointers = [open(file, 'r') for file in sortedFiles[:maxNumOpenFiles]]
                    break
                except OSError as e:
                    if e.errno == 24:
                        printer(f'Attempted to open too many files at once (OSError errno 24)')
                        maxNumOpenFiles = max(1, int(maxNumOpenFiles / 2))
                        printer(f'Reducing the number of files that can be opened by half to {maxNumOpenFiles}')
                        continue
                    raise e
                        
            printer(f'Merging {len(sortedFilesPointers):,}')
            
            # Merge and write
            with open(mergedFile.name, 'w') as f:
                f.writelines(heapq.merge(*sortedFilesPointers))
            
            # Close all of the open files
            for file in sortedFilesPointers:
                file.close()
              
            # prepare for the next set to be merged
            sortedFiles = sortedFiles[maxNumOpenFiles:] + [mergedFile.name]
        
        shutil.move(mergedFile.name, fpOutput)

def startMultiprocessing(fpInputs, fpOutput, mpPool, numThreads, maxOpenFiles, fpVcf):
    printer('Extracting off-targets using multiprocessing approach')
    
    printer(f'Allowed processes: {numThreads}')
    
    fpTempDir = tempfile.TemporaryDirectory()
    printer(f'Created a temporary directory for intermediate files: {fpTempDir.name}')

    if len(fpInputs) == 1 and os.path.isdir(fpInputs[0]):
        fpInputs = glob.glob(
            os.path.join(
                fpInputs[0], 
                '*'
            )
        )

    if len(fpInputs) == 1:
        printer('Only one input file to process')
        
        fpExplodeTempDir = tempfile.TemporaryDirectory()
        printer('Attempting to explode multi-FASTA file')
        printer(f'Writing each file to a temporary directory: {fpExplodeTempDir.name}')
        
    
        fpInputs = explodeMultiFastaFile(
            fpInputs[0],
            fpExplodeTempDir.name
        )
        
        printer(f'Exploded into {len(fpInputs)} files')

    args = [
        (
            [fpInput],
            fpTempDir.name,
            fpVcf
        ) for fpInput in fpInputs
    ]

    printer(f'Beginning to process {len(args)} files')

    numTargets = mpPool.starmap(
        processingNode,
        args
    )

    printer(f'Processing completed. Found {sum(numTargets):,} targets.')
    
    printer('Preparing for ISSL by sorting all intermediate files')
    
    # sort in batches
    paginatedSort(
        glob.glob(
            os.path.join(
                fpTempDir.name, 
                '*'
            )
        ), 
        fpOutput,
        mpPool,
        maxNumOpenFiles=maxOpenFiles
    )

def main():
    parser = argparse.ArgumentParser(description='Extract CRISPR target sites for Crackling-phased.')
    
    parser.add_argument('output',
        help='A file to write the off-targets to.'
    )
    parser.add_argument('inputs',
        help='A space-separated list of paths or a path containing wildcards to FASTA files.',
        nargs='+'
    )
    parser.add_argument('--vcf', 
        help='Path to VCF file.',
        required=True
    )
    parser.add_argument('--maxOpenFiles', 
        help='The number of files allowed to be opened by a single process. Use `ulimit -n` to find out your system setting.', 
        type=int,
        default=1000, 
        required=False
    )
    parser.add_argument('--threads', 
        help='The number of threads to use. Default is `os.cpu_count()`.', 
        type=int,
        default=os.cpu_count(), 
        required=False
    )

    args = parser.parse_args()

    # Create multiprocessing pool
    # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.starmap
    mpPool = multiprocessing.Pool(args.threads)

    startMultiprocessing(
        args.inputs, 
        args.output, 
        mpPool, 
        args.threads, 
        args.maxOpenFiles,
        args.vcf
    )
    
    # Clean up. Close multiprocessing pool
    mpPool.close()

    printer('Goodbye.')

if __name__ == '__main__':
    main()
