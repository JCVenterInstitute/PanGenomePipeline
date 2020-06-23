/* 
Program to sort kmers for a set of genomes into file buckets based on the first 9 basepairs of each kmer.
A file with a list of file names one per line is input using the -g option. Each file should contain a multifasta file with the contigs
for an entire genome of the same species. Preferably for haplotype genomes this would be a single haplotype. Results will be improved if
the genomes are as complete and accurate as possible. This is the first program of a set of three to determine unique anchor sequences within
genomes. This program divides the program into smaller problems and preprocesses the genome data. All 23mers from each genome will be output
into a set of 262144 files (4**9). DNA is double stranded and so each 23mer has an equivalent reverse complement 23mer on the other strand.
Since strandedness is not necessarily conserved across genome multifasta files each 23mer will be represented by a canonical 23mer which is the
lesser integer value between the 46 bit encoding of the 23mer and the reverse complement 23mer.
This program creates a two tiered directory structure within the current working directory where at the
first level 64 directories are created for each 3mer (AAA, AAC, AAG, ... TTT). Within each first level directory the same 64 directories are
created. Within each second level directory 64 files are created with the same names to contain the 23mers which begin with the same 9mer as
the concatenation of the first level directory, the second level directory, and the file name. The 23mers are output as binary representations
of a C struct Kmer composed of: a 64 bit unsigned integer to hold the bit encoding of the 23mer which is contained in the lower 46 bits, a
32 bit integer for the position of the last base pair of the 23mer in a contig, a 16 bit unsigned integer for the genome number which corresponds
to the order in the input file names file, and a 16 bit unsigned integer for the contig number within the multifasta file. The next program in
the anchor pipeline will sort each file by 23mer as the primary key and genome number as the secondary key. 23mers which occur more than once in
any genome are not unique and will be ignored for anchors. 23mers which appear in fewer than a user set threshold will be treated as possible
errors and not be used for anchors. The third program in the anchor pipeline will generate anchors based on the reduced set of 23mers produced by
the second program. The four basepairs are encoded as two bits: A 00, C 01, G 10, T 11. Any basepair ambiguity code such as N terminates and restarts
determination of 23mers but is included as part of the position within the contig.
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <stdint.h>

#define MAX_UINT16 65535
#define MAX_UINT32 4294967295
#define MAX_INT32 2147482647
#define KMER_BUFFER_LEN 1000
#define KMER_BUFFER_NUMBER 262144
 
struct Kmer { /* this struct captures the each 23mer and its context and is output to one of the 262,144 bucket files */
  uint64_t kmer; /* the lower 46 bits encodes a 23mer given a two bit encoding of each basepair */
  int32_t  pos; /* the position of the last basepair of the 23mer within the contig, the first basepair of a contig has position 1, a negative value indicates the reverse complement of the 23mer appears in the contig */
  uint16_t  genome; /* the genome number of the 23mer based on the order of the genome in the input file names file beginning with 0 */
  uint16_t   contig; /* the contig number of the 23mer based on the order of the contig in the multifasta genome file beginning with 0 */
};

/*
  This subroutine reads through a genome one basepair at a time to construct 23mers and output them to the buclet files
 */
void kmer_bucket_sort_genome( FILE * fp_fasta, char * genome_file_name, uint16_t genome_number, char ** bucket_file_names, struct Kmer ** kmer_buffers, int * kmer_buffer_indices)
{
  FILE * fp_bucket;
  size_t getline_return;
  char * fasta_line = NULL;
  size_t fasta_line_malloc_len = 0;
  int fgetc_return;
  char cur_char, prev_char;
  uint64_t cur_kmer = 0;
  uint64_t revc_kmer = 0;
  uint64_t cur_bp, revc_bp;
  int num_ambs = 0;
  bool reset_kmer = true;
  bool new_contig = false;
  uint64_t kmer_mask = 01777777777777777; /* This is used to mask only the lower 46 bits used for the 23mer when shifting */
  uint64_t bucket_mask = 01777776000000000; /* This is used to select the upper 18 bits of the lower 46 bits to determine the 9 basepairs used for bucket sorting */
  uint16_t contig_number = 0;
  int kmer_bp_count = 0;
  int contig_pos = 0;
  int kmer_bucket;
  
  if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_fasta)) == -1) {
    fprintf(stderr, "%s appears to be empty.\n%s", genome_file_name);
    exit(EXIT_FAILURE);
  }
  if (fasta_line[0] != '>') {
    fprintf(stderr, "First line of fasta file %s does not begin with a >.\n%s", genome_file_name, fasta_line);
    exit(EXIT_FAILURE);
  }
  prev_char = '\n';
  while ((fgetc_return = fgetc(fp_fasta)) != EOF) {
    cur_char = (char) fgetc_return;
    if (isspace(cur_char)) {
      /* ignore white space except for keeping track of newlines to determine fasta header lines */
      prev_char = cur_char;
    } else {
      switch (cur_char) {
	/* determine the next basepair and bit encoding for it both forward strand and reverse complement (starts at upper end of 46 bits) */
	case 'a':
	case 'A':
	  cur_bp = 0;
	  revc_bp = 06000000000000000;
	  break;
	case 'c':
	case 'C':
	  cur_bp = 01;
	  revc_bp = 04000000000000000;
	  break;
	case 'g':
	case 'G':
	  cur_bp = 02;
	  revc_bp = 02000000000000000;
	  break;
	case 't':
	case 'T':
	  cur_bp = 03;
	  revc_bp = 0;
	  break;
	/* look for DNA ambiguity codes which terminates 23mer determination and restarts it with the next basepair */
	case 'Y':
	case 'R':
	case 'W':
	case 'S':
	case 'K':
	case 'M':
	case 'D':
	case 'V':
	case 'H':
	case 'B':
	case 'N':
	case 'y':
	case 'r':
	case 'w':
	case 's':
	case 'k':
	case 'm':
	case 'd':
	case 'v':
	case 'h':
	case 'b':
	case 'n': /* nucleotide ambiguity codes */
	  num_ambs++;
	  reset_kmer = true;
	  break;
	case '>':
	  if (prev_char == '\n') {
	    /* Start of a new contig. Read in header line but ignore it other than incrementing contig number */
	    if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_fasta)) == -1) {
	      fprintf (stderr, "empty contig fasta header: %d for genome %s!\n", contig_number, genome_file_name);
	      exit(EXIT_FAILURE);
	    }
	    new_contig = true;
	    if (contig_number == 65535) {
	      fprintf (stderr, "maximum number of contigs: %d exceeded for genome %s!\n", contig_number, genome_file_name);
	      exit(EXIT_FAILURE);
	    }
	    contig_number++;
	  } else {
	    fprintf (stderr, "Unexpected > not at beginning of line in fasta file %s.\n", genome_file_name);
	    exit(EXIT_FAILURE);
	  }
	  break;
	default:
	  if (isprint (cur_char)) {
	    fprintf (stderr, "Unexpected character in fasta file %s: `-%c'.\n", genome_file_name, cur_char);
	  } else {
	    fprintf (stderr, "Unexpected unprintable character in fasta file %s: `\\x%x'.\n", genome_file_name, cur_char);
	  }
	  exit(EXIT_FAILURE);
	}
      if (new_contig) {
	prev_char = '\n';
	kmer_bp_count = 0; 
	contig_pos = 0;
	cur_kmer = 0;
	revc_kmer = 0;
      } else {
	prev_char = cur_char;
	contig_pos++; /* increment the contig position before 23mer determination since we start numbering at 1 in the contig */
	if (!reset_kmer) {
	  /* calculate the current 23mer by shifting on the current basepair for both the forward strand and reverse complement */
	  kmer_bp_count++; /* this is used to determine if the 23mer has been filled up */
	  cur_kmer = ((cur_kmer << 2) & kmer_mask) | cur_bp;
	  revc_kmer = (revc_kmer | revc_bp) >> 2;
	  if (kmer_bp_count >= 23) {
	    if (cur_kmer < revc_kmer) { /* use canonical kmer which ever is less */
	      kmer_bucket = (int) ((bucket_mask & cur_kmer) >> 28);
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].kmer = cur_kmer;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].pos = contig_pos;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].genome = genome_number;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].contig = contig_number;
	    } else {
	      kmer_bucket = (int) ((bucket_mask & revc_kmer) >> 28);
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].kmer = revc_kmer;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].pos = -contig_pos;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].genome = genome_number;
	      kmer_buffers[kmer_bucket][kmer_buffer_indices[kmer_bucket]].contig = contig_number;
	    }
	    kmer_buffer_indices[kmer_bucket]++;
	    if (kmer_buffer_indices[kmer_bucket] == 1000) {
	      /* kmer bucket buffer is full so output it and reset index */
	      kmer_buffer_indices[kmer_bucket] = 0;
	      fp_bucket = fopen(bucket_file_names[kmer_bucket], "a");
	      if (fp_bucket == NULL) {
		fprintf (stderr, "Could not open file %s\n", bucket_file_names[kmer_bucket]);
		exit(EXIT_FAILURE);
	      }
	      if (fwrite(kmer_buffers[kmer_bucket], sizeof(struct Kmer), 1000, fp_bucket) != 1000) {
		fprintf (stderr, "Could not complete write to  file %s\n", bucket_file_names[kmer_bucket]);
		exit(EXIT_FAILURE);
	      }
	      fclose(fp_bucket);
	    }
	  }
	} else {
	  kmer_bp_count = 0;
	  cur_kmer = 0;
	  revc_kmer = 0;
	}
      }
    }
  }
  free(fasta_line);
  return;
}

int
main (int argc, char **argv)
{
  char *genomes_file = NULL;
  int index;
  int getopt_return;
  FILE * fp_file_names;
  FILE * fp_bucket;
  char * file_name_line = NULL;
  size_t file_name_line_malloc_len = 0;
  size_t getline_return;
  uint16_t genome_number = 0; /* this gets assigned and used to distinguish different genomes for each k-mer */
  bool max_genomes_exceeded = false; /* set if too many genomes for a uint16_t (65,536) */
  struct Kmer ** kmer_buffers; /* this is the malloced array of arrays for each kmer buffer in the top level buckets used to buffer for output */
  struct Kmer * kmer_buffers_pool; /* this is the pool of malloc storage to use for the k-mer buffer arrays */
  int * kmer_buffer_indices; /* malloced array of indices into the kmer buffers */
  char bps[4] = {'A','C','G','T'};
  char ** bucket_file_names; /* array of malloced fixed length strings for bucket file names */
  char * bucket_file_names_pool; /* this is the pool pf malloc storage for the bucket file names */
  char tmp_file_dir_name[12];
  int i,j,k,l,m,n,o,p,q;

  opterr = 0;

  while ((getopt_return = getopt (argc, argv, "g:")) != -1) {
    switch (getopt_return) {
      case 'g':
        genomes_file = optarg;
        break;
      case '?':
        if (optopt == 'g') {
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        } else if (isprint (optopt)) {
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        } else {
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
	}
        exit(EXIT_FAILURE);
      default:
        exit(EXIT_FAILURE);
    }
  }

  fprintf (stderr, "Genomes File = %s\n", genomes_file);

  for (index = optind; index < argc; index++) {
    fprintf (stderr, "Non-option argument %s\n", argv[index]);
  }

  /* generate and store directories and files names for the 262,144 buckets to sort 23mers into */
  /* this is the array of bucket files names */
  bucket_file_names = (char **) malloc((size_t) (KMER_BUFFER_NUMBER * sizeof(char *)));
  if (bucket_file_names  == NULL) {
      fprintf (stderr, "Could not allocate memory for bucket_file_names\n");
    exit(EXIT_FAILURE);
  }

  /* this is the pool of storage for the bucket file names - relies on the file names being 12 characters including the null '\0' character */
  bucket_file_names_pool = (char *) malloc((size_t) (KMER_BUFFER_NUMBER * 12 * sizeof(char)));
  if (kmer_buffers_pool == NULL) {
    fprintf (stderr, "Could not allocate memory for bucket_file_names_pool\n");
    exit(EXIT_FAILURE);
  }

  /* assign the correct storage pool area for each bucket file name */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
   bucket_file_names[index] = bucket_file_names_pool + (size_t) (index * 12 * sizeof(char));
  }

  /* create two level directory structure based on 3mers and bucket file names including the first and second level directories */
  index = 0;
  for (i = 0; i < 4; i++) {
    tmp_file_dir_name[0] = bps[i];
    for (j = 0; j < 4; j++) {
      tmp_file_dir_name[1] = bps[j];
      for (k = 0; k < 4; k++) {
	tmp_file_dir_name[2] = bps[k];
	tmp_file_dir_name[3] = '\0';
	if (mkdir(tmp_file_dir_name, 0777)) {
	  fprintf (stderr, "Could not make directory %s\n", tmp_file_dir_name);
	  exit(EXIT_FAILURE);
	}
	tmp_file_dir_name[3] = '/';
	for (l = 0; l < 4; l++) {
	  tmp_file_dir_name[4] = bps[l];
	  for (m = 0; m < 4; m++) {
	    tmp_file_dir_name[5] = bps[m];
	    for (n = 0; n < 4; n++) {
	      tmp_file_dir_name[6] = bps[m];
	      tmp_file_dir_name[7] = '\0';
	      if (mkdir(tmp_file_dir_name, 0777)) {
		fprintf (stderr, "Could not make directory %s\n", tmp_file_dir_name);
		exit(EXIT_FAILURE);
	      }
	      tmp_file_dir_name[7] = '/';
	      for (o = 0; o < 4; o++) {
		tmp_file_dir_name[8] = bps[o];
		for (p = 0; p < 4; p++) {
		  tmp_file_dir_name[9] = bps[p];
		  for (q = 0; q < 4; q++) {
		    tmp_file_dir_name[10] = bps[q];
		    tmp_file_dir_name[11] = '\0';
		    /* the first of these bucket file names should be AAA/AAA/AAA\0 */
		    if (strcpy(bucket_file_names[index], tmp_file_dir_name) == NULL) {
		      fprintf (stderr, "strcpy for bucket_file_names index %d failed\n", index);
		      exit(EXIT_FAILURE);
		    }
		    index++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /* allocate the array of indices for each kmer bucket buffer */
  kmer_buffer_indices = (int *) malloc((size_t) (KMER_BUFFER_NUMBER * sizeof(int)));
  if (kmer_buffer_indices == NULL) {
      fprintf (stderr, "Could not allocate memory for kmer_buffer_indices\n");
    exit(EXIT_FAILURE);
  }

  /* initialize the kmer bucket buffer indices to 0 */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
    kmer_buffer_indices[index] = 0;
  }

  /* allocate the array of pointers to the kmer buffers */
  kmer_buffers = (struct Kmer **) malloc((size_t) (KMER_BUFFER_NUMBER * sizeof(struct Kmer *)));
  if (kmer_buffers == NULL) {
      fprintf (stderr, "Could not allocate memory for kmer_buffers\n");
    exit(EXIT_FAILURE);
  }

  /* allocate the pool of storage for the kmer buffers */
  kmer_buffers_pool = (struct Kmer *) malloc((size_t) (KMER_BUFFER_NUMBER * KMER_BUFFER_LEN * sizeof(struct Kmer)));
  if (kmer_buffers_pool == NULL) {
    fprintf (stderr, "Could not allocate memory for kmer_buffers_pool\n");
    exit(EXIT_FAILURE);
  }

  /* set up the kmer buffers to point to the assigned storage buffer within the storage pool */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
    kmer_buffers[index] = kmer_buffers_pool + (size_t) (index * KMER_BUFFER_LEN * sizeof(struct Kmer));
  }
  
  fp_file_names = fopen(genomes_file, "r");
  if (fp_file_names == NULL) {
      fprintf (stderr, "Could not open file %s\n", genomes_file);
    exit(EXIT_FAILURE);
  }

  /* loop through each genome in the genomes file input file */
  while ((getline_return = getline(&file_name_line, &file_name_line_malloc_len, fp_file_names)) != -1) {
    FILE * fp_file_name;
    int file_name_len;
    if (max_genomes_exceeded) {
      fprintf (stderr, "Maximum number of genomes: %d exceeded!\n", genome_number);
      exit(EXIT_FAILURE);
    }
    if (file_name_line == NULL) {
      fprintf (stderr, "file_name_line is NULL!\n");
      exit(EXIT_FAILURE);
    }
    file_name_len = strlen(file_name_line);
    if (file_name_len <= 1) {
      fprintf (stderr, "Unexpected empty line in genomes file: %s\n", genomes_file);
      exit(EXIT_FAILURE);
    }
    file_name_line[file_name_len - 1] = '\0';
    fp_file_name = fopen(file_name_line, "r");
    if (fp_file_name == NULL) {
      fprintf (stderr, "Could not open file %s\n", file_name_line);
      exit(EXIT_FAILURE);
    }
    kmer_bucket_sort_genome(fp_file_name, file_name_line, genome_number, bucket_file_names, kmer_buffers, kmer_buffer_indices);
    if (genome_number == 65535) {
      max_genomes_exceeded = true;
    } else {
      genome_number++;
    }
    fclose(fp_file_name);
  }
  free(file_name_line);

  fclose(fp_file_names);

  /* flush any 23mers remaining in the kmer buffers */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
    if (kmer_buffer_indices[index] != 0) {
      fp_bucket = fopen(bucket_file_names[index], "a");
      if (fp_bucket == NULL) {
	fprintf (stderr, "Could not open file %s\n", bucket_file_names[index]);
	exit(EXIT_FAILURE);
      }
      if (fwrite(kmer_buffers[index], sizeof(struct Kmer), kmer_buffer_indices[index], fp_bucket) != kmer_buffer_indices[index]) {
	fprintf (stderr, "Could not complete write to  file %s\n", bucket_file_names[index]);
	exit(EXIT_FAILURE);
      }
      fclose(fp_bucket);
    }
  }
  exit(EXIT_SUCCESS);
}
