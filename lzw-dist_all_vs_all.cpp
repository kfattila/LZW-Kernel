/*********************************************************************
    lzw-kernel_all_vs_all.c (version 0.1)                                          
                                                                    
    Compute compression based distance between all pairs of sequences
    from a FASTA file using the lzw compressor.

    Reference:
    
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "seqio.h"
#include <new>
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
    fprintf(stderr, "lzw-dist_all_vs_all : compute the compression based distance between all pairs of proteins in a database using the LZW compressor.\n\nUsage: lzw-dist_all_vs_all sequence_file.fasta [> output_file.dmx]\n");
    exit(1);
  }
	//empty string
	char e[2];
	e[0] = 100;
    buflen = 10000;
	buffer = new int[buflen*ALPHABET_NUM];
	memset(&buffer[26], 0, sizeof(int)*(buflen-1)*ALPHABET_NUM);

	int i,j, seq_num;
    
    // Count the sequences in the file

    /* open the database file and search the first query sequenced not processed yet */
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
    
    for (seq_num = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++seq_num) {
      free(seq1);
    }


    int *C  = new int[seq_num];
    char** seqs = new char*[seq_num];
    SEQINFO** headers = new SEQINFO*[seq_num];

	if (buffer == 0 || C == 0){
      fprintf(stderr,"Unable to allocate memory \n");
      exit(1);
	}

	for (i = 0; i < ALPHABET_NUM; i++)
		buffer[i] = (i+1)*ALPHABET_NUM;

	bufidx = buflen;

    /* open the database file and search the first query sequenced not processed yet */
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }
    
    for (i = 0; ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL); ++i) {
        sip = seqfinfo(sfp,1);
        headers[i] = sip;
        printf("\t%s", headers[i]->description); //Print header
        normalize_seq(seq1);
        seqs[i] = seq1;
        C[i] = compress(seq1, NULL);
    }
    printf("\n");   

    for (int k = 0; k < seq_num; ++k){
        printf("%s", headers[k]->description);
        for (int l = 0; l < seq_num; ++l){
            float lzw_kernel_value = (float)(compress(seqs[k], seqs[l]) - min(C[k],C[l])) / max(C[k], C[l]);
            printf("\t%.3f",lzw_kernel_value);
        }
        printf("\n");
    }
    
    //We clean up, because we are nice boyz!    
    for (int k = 1; k < seq_num; ++k){
        free(seqs[k]);
        free(headers[k]);
    }
	delete[] C;
	delete[] buffer;
    delete[] seqs;
    delete[] headers;
    
    return 0;
}

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

			if (buffer[pointer+str[i]] == 0 ){
				buffer[pointer+str[i]] = (bufidx++)*ALPHABET_NUM;
				cnt += 1;
				pointer = buffer[str[i]];

				if (bufidx == buflen){
					int* buf = new int[buflen*ALPHABET_NUM*2];
					memcpy(buf, buffer, buflen*sizeof(int)*ALPHABET_NUM);
					memset(&buf[buflen*ALPHABET_NUM], 0, buflen*ALPHABET_NUM*sizeof(int));
					delete[] buffer;
					buffer = buf;
					buflen *= 2;
				}
			} else
				pointer = buffer[pointer+str[i]];
		}
	}
	return cnt;
}

void normalize_seq(char* seq){
    //_strlwr(seq);
    int j;
    for (j = 0; seq[j]; j++){
        if (seq[j] < 95)
            seq[j] += 32;
        
        seq[j] -= 97;

        if (seq[j] < 0 || seq[j] > 26){
            fprintf(stderr,"Unsupported symbol in sequence, at position %ith, symbol:%c", j, seq[j]+97);
            exit(1);
        }
    }
    seq[j] = 100;   
}
