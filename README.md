# LZW-Kernel
Implementation of the LZW-Kernel function for biological sequences.

For more information see:
    Reference: 
    Authors: Gleb Filatov, Bruno Bauwens, Attila Kert\'esz-Farkas
    Title: LZW-Kernel: fast kernel utilizing variable length code 
    blocks from LZW compressors for protein sequence classification
    email: akerteszfarkas@hse.ru
    
# Install:
    $ make

# Usage:
    For LZW-Kernel use:
    $ lzw-kernel_all_vs_all input_sequences.fasta > sequnces.kernel_mx
    
    For LZW-NCD use:
    $ lzw-dist_all_vs_all input_sequences.fasta > sequnces.dist_mx
    
# Example:
    $ lzw-kernel_all_vs_all test.fasta > test.kernel_mx
