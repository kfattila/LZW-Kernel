/*********************************************************************
    lzw-kernel_all_vs_all.c (version 1.0)                                          
                                                                    
    Compute the LZW kernel between all pairs of sequeences
    from a FASTA format file.

    Reference: 
    Authors: Gleb Filatov, Bruno Bauwens, Attila Kert\'esz-Farkas
    Title: LZW-Kernel: fast kernel for protein sequence comparison utilizing 
    variable length code blocks from LZW compressors
    email: akerteszfarkas@hse.ru
    
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "seqio.h"
#include <math.h>
#include <new>
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
    double gamma;

    if (argc != 2) {
        fprintf(stderr, "lzw-kernel_all_vs_all : compute the lzw-kernel between all pairs of proteins in a database.\n\nUsage: lzw-kernel_all_vs_all sequence_file.fasta [> output_file.dmx]\n");
        exit(1);
    }
    buflen = 1000000;
	buffer = new int[buflen*ALPHABET_NUM];   // approx. 100 Mb of memory.

	int i, j, seq_num;
    int compressed_len;
    
    // Count the sequences in the file

    // open the database file and
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
     
    // and count the sequences
    for (seq_num = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++seq_num) {
      free(seq1);
    }

    int* C  = new int[seq_num];
    SEQINFO** headers = new SEQINFO*[seq_num];
    float** kernel_mx = new float*[seq_num];
    int** code_tree = new int*[seq_num];

    if (buffer == NULL || C == NULL || headers == NULL || kernel_mx == NULL || code_tree == NULL){
      fprintf(stderr,"Unable to allocate memory \n");
      exit(1);
	}

	for (i = 0; i < ALPHABET_NUM; ++i)
		buffer[i] = (i+1)*ALPHABET_NUM;

	bufidx = buflen;

    /* open the database file  */
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
    
    // and process each sequences to construct the code dictionaries
    for (i = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++i) {
        sip = seqfinfo(sfp,1);  //get the header
        headers[i] = sip;  // store the header.
        printf("\tHeader"); //Print header

        compressed_len = compress(seq1);   //run the LZW compression
        
        // store the codes
        code_tree[i] = new int[(bufidx)*ALPHABET_NUM];  //Allocate memory for the codes   
        memcpy(code_tree[i], buffer, sizeof(int)*(bufidx)*ALPHABET_NUM);   //Store the codes

        C[i] = get_common_code_len(code_tree[i], code_tree[i], 0, 0, 1); // get the length of the code blocks. It will be used for normalization
        gamma += C[i];
     
    }
    gamma /= i;
    printf("\n"); 
    int common_length;
    //Calculate the LZW kernel for each sequence pair
    for (int k = 0; k < seq_num; ++k){
        printf("Header"); //Print header
        
        //The LZW kernel is symmetric, so use k(y,x) := k(x,y) 
        for (int l = 0; l < k; ++l){
            printf("\t%g", kernel_mx[l][k]);            
        }
        
        kernel_mx[k] = new float[seq_num];
        //Calculate the LZW-Kernel k(x_k,x_l) for two sequences x_k and x_l. The calculation requires only parsing the common branches of the dictionaries simultenously for common code words.
        for (int l = k; l < seq_num; ++l){
            common_length = get_common_code_len(code_tree[k], code_tree[l], 0, 0, 1);//find common codes and return the sum of their lengths
            kernel_mx[k][l] = exp((common_length-(0.5*(C[k]+C[l])))/gamma);
            printf("\t%g", kernel_mx[k][l]);
        }
        printf("\n");        
    }  
    
    //We clean up, because we are nice boyz!
    
    for (i = 0; i < seq_num; ++i){
        delete[] kernel_mx[i];
        delete[] code_tree[i];
        delete[] headers[i];
    }
    
    delete[] kernel_mx;
    delete[] code_tree;
    delete[] headers;
	delete[] C;
	delete[] buffer;
    
    return 0;
}

int compress(char* str){

	int i;

	int cnt = 0;
	int pointer = 0;
	memset(buffer, 0, sizeof(int)*(bufidx)*ALPHABET_NUM);
	bufidx = 1;
    for (i = 0; str[i] != '\0'; ++i){
        if (str[i] < 95)
            str[i] += 32;
        
        str[i] -= 97;

        if (str[i] < 0 || str[i] > ALPHABET_NUM){
            fprintf(stderr,"Unsupported symbol in sequence, at position %ith, symbol:%c", i, str[i]+97);
            exit(1);
        }
        if (buffer[str[i]] == 0)
            buffer[str[i]] = (bufidx++)*ALPHABET_NUM;
    }
    str[i] = 100;
    for (i = 0; str[i] != 100; ++i){
        if (buffer[pointer+str[i]] == 0 ){
            buffer[pointer+str[i]] = (bufidx++)*ALPHABET_NUM;
            ++cnt;
            pointer = buffer[str[i]];

            if (bufidx >= buflen){
                int* buf = new int[buflen*ALPHABET_NUM*2];
                memcpy(buf, buffer, buflen*sizeof(int)*ALPHABET_NUM);
                memset(&buf[buflen*ALPHABET_NUM], 0, buflen*ALPHABET_NUM*sizeof(int));
                delete[] buffer;
                buffer = buf;
                buflen *= 2;
            }
        } else {
            pointer = buffer[pointer+str[i]];
        }
    }
	return ++cnt;
}
int get_common_code_len(int* tree1, int* tree2, int idx1, int idx2, int depth){
    int sub_common_code_len = 0;
    for (int i = 0; i < ALPHABET_NUM; ++i)
        if (tree1[idx1+i] != 0 && tree2[idx2+i] != 0)
            sub_common_code_len += get_common_code_len(tree1, tree2, tree1[idx1+i], tree2[idx2+i], depth+1) + depth;
        
    return sub_common_code_len;
}


