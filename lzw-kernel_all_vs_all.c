/*********************************************************************
    lzw-kernel_all_vs_all.c (version 1.0)                                          
                                                                    
    Compute the LZW-kernel between all sequence pairs
    from a FASTA format file.
    
    Usage:
    lzw-kernel_all_vs_all sequences.fasta > sequences.kernel_mx

    Reference: 
    Authors: Gleb Filatov, Bruno Bauwens, Attila Kert\'esz-Farkas
    Title: LZW-Kernel: fast kernel utilizing variable length code 
    blocks from LZW compressors for protein sequence classification
    email: akerteszfarkas@hse.ru    
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "seqio.h"
#include <math.h>
#include "malloc.h"
#include "memory.h"

static int* buffer;
static int buflen;
static int bufidx;

#define ALPHABET_NUM 26  //Size of the alphabet

int compress(char *x);
int get_common_code_len(int* tree1, int* tree2, int idx1, int idx2, int depth);

int main(int argc, char *argv[]){
    
    SEQFILE *sfp;
    SEQINFO *sip;
    char *seq1;
    int len;
    double gamma=0;

    if (argc != 2) {
        fprintf(stderr, "lzw-kernel_all_vs_all: compute the lzw-kernel\
        between all pairs of proteins in a database.\n\nUsage: \
        lzw-kernel_all_vs_all sequence_file.fasta [> output_file.kernel_mx]\n");
        exit(1);
    }
    buflen = 1000000;
	buffer = malloc(sizeof(int)*buflen*ALPHABET_NUM);   // Approx. 100 Mb of RAM.

	int i, k, l, seq_num;
    int compressed_len;
    
    // Count the sequences in the input fasta file
    // Open the fasta file and ...
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
     
    // ... and count the sequences
    for (seq_num = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++seq_num) {
      free(seq1);
    }

    // Allocate memory for the kernel matrix, for the protein sequence headers, 
    // and for the code words for every protein sequence, respectively.
    int* C  = malloc(sizeof(int)*seq_num);
    SEQINFO** headers = malloc(sizeof(SEQINFO*)*seq_num);
    float** kernel_mx = malloc(sizeof(float*)*seq_num);
    int** code_tree = malloc(sizeof(int*)*seq_num);

    if (buffer == NULL || C == NULL || headers == NULL || kernel_mx == NULL || code_tree == NULL){
      fprintf(stderr,"Unable to allocate memory \n");
      exit(1);
	}

	for (i = 0; i < ALPHABET_NUM; ++i)
		buffer[i] = (i+1)*ALPHABET_NUM;

	bufidx = buflen;

    // Open the fasta file again and ...
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
    
    // ... and process each sequence to construct the code dictionaries
    for (i = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++i) {
        sip = seqfinfo(sfp,1);  // Get the header
        headers[i] = sip;  // Store the header.
        printf("\t%s", headers[i]->description); // Print header to the output

        compressed_len = compress(seq1);   // Run the LZW compression to get the LZW-code words
        
        // Store the code-words
        code_tree[i] = malloc(sizeof(int)*(bufidx)*ALPHABET_NUM);  // Allocate memory for the code-words  
        memcpy(code_tree[i], buffer, sizeof(int)*(bufidx)*ALPHABET_NUM);   // Store the code-words
        
        // Get the total length of the code words. It will be used for normalization
        C[i] = get_common_code_len(code_tree[i], code_tree[i], 0, 0, 1); 
        gamma += C[i];
    }
    gamma /= i;
    printf("\n"); 
    int common_length;
    // Calculate the LZW kernel for each sequence pair
    for (k = 0; k < seq_num; ++k){
        printf("%s", headers[k]->description);  // Print the sequence header to the output
        
        // The LZW kernel is symmetric, so use k(y,x) := k(x,y) for l < k
        for (l = 0; l < k; ++l){
            printf("\t%g", kernel_mx[l][k]);            
        }
        
        kernel_mx[k] = malloc(sizeof(float)*seq_num);
        // Calculate the LZW-Kernel k(x,y) for two sequences x,y. 
        // The calculation requires only parsing the common branches of the dictionaries simultenously for common code words.
        for (l = k; l < seq_num; ++l){
            common_length = get_common_code_len(code_tree[k], code_tree[l], 0, 0, 1);// Find common codes and return the sum of their lengths
            kernel_mx[k][l] = exp((common_length - (0.5 * ( C[k] + C[l] ))) / gamma);
            printf("\t%g", kernel_mx[k][l]);
        }
        printf("\n");        
    }  
    //We clean up, because we are nice boyz!
    for (i = 0; i < seq_num; ++i){
        free(kernel_mx[i]);
        free(code_tree[i]);
        free(headers[i]);
    }
    free(kernel_mx);
    free(code_tree);
    free(headers);    
    free(buffer);
    free(C);
    
    return 0;
}

// Run the LZW-compression and return the length of the compressed string.
// For more information about LZW-compression use gogle.
// The code words are stored in a prefix tree. 
// The prefix tree is represented by an excellent method published in:
// George Anton Kiraz: "Compressed storage of sparse finite-state transducers", 
//        in International Workshop on Implementing Automata, 1999, p. 109-121

int compress(char* str){

	int i, seq_len;
	int cnt = 0;
	int pointer = 0;
	memset(buffer, 0, sizeof(int)*(bufidx)*ALPHABET_NUM);
	bufidx = 1;
    for (seq_len = 0; str[seq_len] != '\0'; ++seq_len){
        if (str[seq_len] < 95)
            str[seq_len] += 32;
        
        str[seq_len] -= 97;

        if (str[seq_len] < 0 || str[seq_len] > ALPHABET_NUM){
            fprintf(stderr,"Unsupported symbol in sequence, at position: %d, symbol:%c", seq_len, str[seq_len]+97);
            exit(1);
        }
        if (buffer[(int)str[seq_len]] == 0)
            buffer[(int)str[seq_len]] = (bufidx++)*ALPHABET_NUM;
    }
    for (i = 0; i < seq_len; ++i){
        if (buffer[pointer+(int)str[i]] == 0 ){
            buffer[pointer+(int)str[i]] = (bufidx++)*ALPHABET_NUM;
            ++cnt;
            pointer = buffer[(int)str[i]];

            if (bufidx >= buflen){
                int* buf = malloc(sizeof(int)*buflen*ALPHABET_NUM*2);//new int[buflen*ALPHABET_NUM*2];
                memcpy(buf, buffer, buflen*sizeof(int)*ALPHABET_NUM);
                memset(&buf[buflen*ALPHABET_NUM], 0, buflen*ALPHABET_NUM*sizeof(int));
                free(buffer);
                buffer = buf;
                buflen *= 2;
            }
        } else {
            pointer = buffer[pointer+(int)str[i]];
        }
    }
	return ++cnt;
}

// Simultaneously traverses the common branches of two prefix trees recursively, 
// and returns the sum of the length of the common code words in bot prefix-trees.
int get_common_code_len(int* tree1, int* tree2, int idx1, int idx2, int depth){
    int i;
    int sub_common_code_len = 0;
    // For the input tree nodes check all branches 
    for (i = 0; i < ALPHABET_NUM; ++i)
        //Check if both branches exist 
        if (tree1[idx1+i] != 0 && tree2[idx2+i] != 0)
            //Traverse the subtree
            sub_common_code_len += get_common_code_len(tree1, tree2, tree1[idx1+i], tree2[idx2+i], depth+1) + depth;
        
    return sub_common_code_len;
}
