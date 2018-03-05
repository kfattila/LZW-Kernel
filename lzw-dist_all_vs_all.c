/*********************************************************************
    lzw-kernel_all_vs_all.c (version 1.0)                                          
                                                                    
    Compute the LZW-kernel between all sequence pairs
    from a FASTA format file.
    
    Usage:
    lzw-dist_all_vs_all sequences.fasta > sequences.dist_mx

    Reference: 
    Authors: Gleb Filatov, Bruno Bauwens, Attila Kert\'esz-Farkas
    Title: LZW-Kernel: fast kernel utilizing variable length code 
    blocks from LZW compressors for protein sequence classification
    email: akerteszfarkas@hse.ru    
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "seqio.h"
#include "malloc.h"
#include "memory.h"
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))

#define ALPHABET_NUM 26

static int* buffer;
static int buflen;
static int bufidx;

int compress(char *x, char *y);
void normalize_seq(char* seq);

int main(int argc, char *argv[])
{
    SEQFILE *sfp;
    SEQINFO *sip;
    char *seq1;
    int len;

    if (argc != 2) {
        fprintf(stderr, "lzw-dist_all_vs_all: compute the lzw-ncd\
        between all pairs of proteins in a database.\n\nUsage: \
        lzw-dist_all_vs_all sequence_file.fasta [> output_file.dist_mx]\n");
        exit(1);
    }
    buflen = 10000;
	buffer = malloc(sizeof(int)*buflen*ALPHABET_NUM);//new int[buflen*ALPHABET_NUM];
	memset(&buffer[26], 0, sizeof(int)*(buflen-1)*ALPHABET_NUM);

	int i, k, l, seq_num;
    
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

    // Allocate memory for the protein sequences and for their headers, 
    // and for the length of the compressed sequences (used for normalization).
    int *C  = malloc(sizeof(int)*seq_num);
    char** seqs = malloc(sizeof(char*)*seq_num);
    SEQINFO** headers = malloc(sizeof(SEQINFO*)*seq_num);

	if (buffer == 0 || C == 0){
      fprintf(stderr,"Unable to allocate memory \n");
      exit(1);
	}

	for (i = 0; i < ALPHABET_NUM; i++)
		buffer[i] = (i+1)*ALPHABET_NUM;

	bufidx = buflen;

    // Open the fasta file again and ...
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }

    // ... and compress each sequence.  
    for (i = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++i) {
        sip = seqfinfo(sfp,1);
        headers[i] = sip;
        printf("\t%s", headers[i]->description); //Print header
        normalize_seq(seq1);
        seqs[i] = seq1;
        C[i] = compress(seq1, NULL);
    }
    printf("\n");   

    // Calculate the LZW-NCD for each sequence pair
    for (k = 0; k < seq_num; ++k){
        printf("%s", headers[k]->description);
        for (l = 0; l < seq_num; ++l){
            float lzw_kernel_value = (float)(compress(seqs[k], seqs[l]) - min(C[k],C[l])) / max(C[k], C[l]);
            printf("\t%.3f",lzw_kernel_value);
        }
        printf("\n");
    }
    
    //We clean up, because we are nice boyz!    
    for (k = 1; k < seq_num; ++k){
        free(seqs[k]);
        free(headers[k]);
    }
    free(C);
    free(buffer);
    free(seqs);
    free(headers);

    return 0;
}

// Run the LZW-compression and return the length of the compressed string. 
// The two input string is treated as one concatenated string. 
// For more information about LZW-compression use gogle.
// The code words are stored in a prefix tree. 
// The prefix tree is represented by an excellent method published in:
// George Anton Kiraz: "Compressed storage of sparse finite-state transducers", 
//        in International Workshop on Implementing Automata, 1999, p. 109-121
int compress(char* strA, char* strB){

	int i,j;

	int cnt = 1;
	int pointer = 0;
	char* str;
	memset(&buffer[26], 0, sizeof(int)*(bufidx-2)*ALPHABET_NUM);

	bufidx = 28;

	for (j = 0; j < 2; ++j){

		if (j == 0) str = strA;
		if (j == 1) str = strB;
        
        if (str == NULL)
            continue;

		for (i = 0; str[i] != 100; ++i){

			if (buffer[pointer+(int)str[i]] == 0 ){
				buffer[pointer+(int)str[i]] = (bufidx++)*ALPHABET_NUM;
				cnt += 1;
				pointer = buffer[(int)str[i]];

				if (bufidx == buflen){
					int* buf = malloc(sizeof(int)*buflen*ALPHABET_NUM*2);// new int[buflen*ALPHABET_NUM*2];
					memcpy(buf, buffer, buflen*sizeof(int)*ALPHABET_NUM);
					memset(&buf[buflen*ALPHABET_NUM], 0, buflen*ALPHABET_NUM*sizeof(int));
					buffer = buf;
					buflen *= 2;
				}
			} else
				pointer = buffer[pointer+(int)str[i]];
		}
	}
	return cnt;
}

// Preprocess the protein sequences suitable for using prefix-trees efficiently.
// It replaces all characters (a-z) with numbers (1-26).

void normalize_seq(char* seq){
    int j;
    for (j = 0; seq[j]; ++j){
        if (seq[j] < 95)
            seq[j] += 32;
        
        seq[j] -= 97;

        if (seq[j] < 0 || seq[j] > 26){
            fprintf(stderr,"Unsupported symbol in sequence, at position: %d, symbol:%c", j, seq[j]+97);
            exit(1);
        }
    }
    seq[j] = 100;   
}
