# USTARC (Unitig STitch Advanced constRuction with Colors)
## Overview
USTAR3 is a k-mers set compressor, with colors.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) 
and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given k-mers set.

Additionally, it exploits the possibility of reusing already visited nodes to achieve better compression as shown by [Matchtigs](https://github.com/algbio/matchtigs).

You will find four executables and one bash script:
* `ustar`: the main program
* `ustarx`: a k-mers extractor
* `ustars`: compute k-mers statistics
* `ustar-test`: used for debug
* `validate`: a validation script

## Dependencies
Good news, there are no dependencies. 
However, you'll need [GGCAT](https://github.com/algbio/ggcat) 
in order to compute a colored compacted de Bruijn graph (cdBG) of your multi-fasta file.

## How to download and compile
* `git clone https://github.com/enricorox/USTAR3`.
* `cd USTAR3`
* `cmake . && make -j 4`.

## How to run USTAR3
Run GGCAT first: 
* `./GGCAT -k <kmer-size> -e -c <your-multi-fasta>`

Then USTARC:
* `./ustar -k <kmer-size> -i <GGCAT-output> -D7`

To use the best heuristic, and the optimized RLE for colors, add `-s+u -x-c -e opt_rle`:
* `./ustar -k <kmer-size> -i <GGCAT-output> -D7 -s+c -x-c -e opt_rle`

See the help `./ustar -h` for details and advanced options.

## How to validate the output
You can check that the output file contains the same k-mers of
your bcalm file with your preferred kmer counter.

If you want to check that __kmers and counts__ are correct,
* `./ustarx -k <kmer-size> -i <ustar-fasta> -c <ustar-colors> -s`
* `./validate <kmer-size> <your-multi-fasta> <ustar-kmers-colors>` 

Note that you'll need to install [Jellyfish-2](https://github.com/zippav/Jellyfish-2) in order to use `validate`.

## References

Paper under submission...
