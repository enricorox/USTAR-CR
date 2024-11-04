# USTARC (Unitig STitch Advanced constRuction with Colors)
## Overview
USTAR3 is a k-mers set compressor, with colors.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) 
and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given k-mers set.

Additionally, it exploits the possibility of reusing already visited nodes to achieve better compression as shown by [Matchtigs](https://github.com/algbio/matchtigs).

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


## References

If you are using USTARC in your research, please cite 
`Extremely fast and succint compression of k-mers sets with plain text representation of colored de Bruijn graphs`
(under sumbission).
