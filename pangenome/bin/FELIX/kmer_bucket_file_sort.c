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
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_UINT16 65535
#define MAX_UINT32 4294967295
#define MAX_INT32 2147482647
#define KMER_SIZE 23
#define KMER_BUFFER_LEN 1000
#define KMER_BUFFER_NUMBER 262144
#define CONTIG_SEQ_BUFFER_LEN 100000
#define MIN_ANCHOR_LEN 100
#define OPT_ANCHOR_LEN 1000
 
struct Kmer { /* this struct captures the each 23mer and its context and is output to one of the 262,144 bucket files */
  uint64_t kmer; /* the lower 46 bits encodes a 23mer given a two bit encoding of each basepair */
  int32_t  pos; /* the position of the last basepair of the 23mer within the contig, the first basepair of a contig has position 1, a negative value indicates the reverse complement of the 23mer appears in the contig */
  uint32_t   contig; /* the contig number of the 23mer based on the order of the contig in the multifasta genome file beginning with 0 */
  uint16_t  genome; /* the genome number of the 23mer based on the order of the genome in the input file names file beginning with 0 */
};

struct red_Kmer { /* this struct is a reduced version of the Kmer struct which does not store the k-mer explicitly */
  int32_t  pos; /* the position of the last basepair of the 23mer within the contig, the first basepair of a contig has position 1, a negative value indicates the reverse complement of the 23mer appears in the contig */
  uint32_t   contig; /* the contig number of the 23mer based on the order of the contig in the multifasta genome file beginning with 0 */
  uint16_t  genome; /* the genome number of the 23mer based on the order of the genome in the input file names file beginning with 0 */
  uint16_t   prevalence; /* the number of genomes containing this k-mer */
};

/*
This subroutine is the comparator function for qsort for an array of k-mers in a bucket file.
*/
static int kmer_sort(const void * kmer_ptr1, const void * kmer_ptr2) {
  if (((struct Kmer *)kmer_ptr1)->kmer < ((struct Kmer *)kmer_ptr2)->kmer) {
    return (-1);
  } else if (((struct Kmer *)kmer_ptr1)->kmer > ((struct Kmer *)kmer_ptr2)->kmer) {
    return (1);
  } else {
    return ((int) (((struct Kmer *)kmer_ptr1)->genome - ((struct Kmer *)kmer_ptr2)->genome));
  }
}

/*
This subroutine is the comparator function for qsort for an array of reduced k-mers in a bucket file.
*/
static int red_kmer_sort(const void * kmer_ptr1, const void * kmer_ptr2) {
  if (((struct red_Kmer *)kmer_ptr1)->contig < ((struct red_Kmer *)kmer_ptr2)->contig) {
    return (-1);
  } else if (((struct red_Kmer *)kmer_ptr1)->contig > ((struct red_Kmer *)kmer_ptr2)->contig) {
    return (1);
  } else {
    return ((int) (((struct red_Kmer *)kmer_ptr1)->pos - ((struct red_Kmer *)kmer_ptr2)->pos));
  }
}

/*
This subroutine reads in a specified number of base pairs into an anchor buffer
*/
int read_fasta_kmer(char * genome_file_name, FILE * fp_fasta, int cur_contig, int cur_pos, int num_basepairs, char * contig_seq, int contig_seq_pos)
{
  static int cur_file_pos = 0;
  static int cur_file_contig = -1;
  static char prev_genome_file_name[1024] = "";
  int fgetc_return, getline_return;
  char cur_char;
  char prev_char = '\n';
  size_t fasta_line_malloc_len = 0;
  char * fasta_line = NULL;
  
  fprintf(stderr, "%s %d:%d:%d\n", genome_file_name, cur_contig, cur_pos, num_basepairs);

  if (strncmp(prev_genome_file_name, genome_file_name, (size_t) 1024) != 0) {
    /*  fprintf(stderr, "Prev Genome: %s contig %d last pos %d\n", prev_genome_file_name, cur_file_contig, (cur_file_pos - 1));
	fprintf(stderr, "Genome: %s Prev Genome: %s\n", genome_file_name, prev_genome_file_name); */
    cur_file_pos = 1;
    cur_file_contig = 0;
    strncpy(prev_genome_file_name, genome_file_name, (size_t) 1024);
  }
  
  while ((num_basepairs > 0) && ((fgetc_return = fgetc(fp_fasta)) != EOF)) {
    cur_char = (char) fgetc_return;
    if (isspace(cur_char)) {
      /* ignore white space except for keeping track of newlines to determine fasta header lines */
      prev_char = cur_char;
    } else if (cur_char == '>') {
      if (prev_char == '\n') {
	/* fprintf(stderr, "Genome: %s contig %d last pos %d\n", genome_file_name, cur_file_contig, (cur_file_pos - 1)); */
	cur_file_contig++;
	cur_file_pos = 1;
	/* Start of a new contig. Read in header line but ignore it other than incrementing contig number */
	if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_fasta)) == -1) {
	  fprintf (stderr, "empty contig fasta header: %d for genome %s!\n", cur_file_contig, genome_file_name);
	  exit(EXIT_FAILURE);
	}
	/* fprintf (stderr, "contig fasta header: %d for genome %s\n%s", cur_file_contig, genome_file_name, fasta_line); */
      } else {
	fprintf (stderr, "Unexpected > not at beginning of line in fasta file %s!\n", genome_file_name);
	exit(EXIT_FAILURE);
      }
    } else { /* just read in a nucleotide either catching up to current anchor or part of the anchor */
      if ((cur_file_contig == cur_contig) && (cur_file_pos >= cur_pos)) {
	switch (cur_char) {
	  /* determine the next basepair and bit encoding for it both forward strand and reverse complement (starts at upper end of 46 bits) */
	case 'a':
	case 'A':
	  cur_char = 'A';
	  break;
	case 'c':
	case 'C':
	  cur_char = 'C';
	  break;
	case 'g':
	case 'G':
	  cur_char = 'G';
	  break;
	case 't':
	case 'T':
	  cur_char = 'T';
	  break;
	default:
	  if (isprint (cur_char)) {
	    fprintf (stderr, "Contig %d position %d\nUnexpected character in fasta file %s: %c.\n", cur_file_contig, cur_file_pos, genome_file_name, cur_char);
	  } else {
	    fprintf (stderr, "Contig %d position %d\nUnexpected unprintable character in fasta file %s: \\x%x.\n", cur_file_contig, cur_file_pos, genome_file_name, cur_char);
	  }
	  exit(EXIT_FAILURE);
	}
	contig_seq[contig_seq_pos++] = cur_char;
	num_basepairs--;
      }
      cur_file_pos++;
    }
  }
  /* if (fgetc_return == EOF) {
    fprintf(stderr, "Genome: %s contig %d last pos %d\n", genome_file_name, cur_file_contig, (cur_file_pos - 1));
  }
  free((void *) fasta_line); */
  return (contig_seq_pos);
}

/*
This subroutine writes anchors/medoid sequences out to a file based on an anchor buffer and returns the number of anchors written
*/
int write_anchors(int anchor_prevalence, FILE * fp_single, FILE * fp_cluster, FILE * fp_match, FILE * fp_pgg, FILE * fp_anchors, int anchor_number, char * contig_seq, int stop_pos, bool anchor_break, bool core_anchor)
{
  int num_anchors_written = 0;
  int i,j;

  if (stop_pos < MIN_ANCHOR_LEN) {
  } else if (stop_pos < ((3 * OPT_ANCHOR_LEN) / 2)) {
    fprintf(fp_anchors, ">medoid_%d\n", anchor_number);
    if (core_anchor) {
      fprintf(fp_single, "%d\n", anchor_number);
    }
    fprintf(fp_cluster, "%d\t%d\n", anchor_number, anchor_prevalence);
    fprintf(fp_match, "%d\tmedoid_%d\n", anchor_number, anchor_number);
    for (j = 0; j < stop_pos; j += 60) {
      int out_len = (stop_pos - j) < 60 ? (stop_pos - j) : 60;
      if (fwrite(&contig_seq[j], sizeof(char), (size_t) out_len, fp_anchors) != out_len) {
	fprintf (stderr, "Could not complete write to anchors file for anchor %d\n", anchor_number);
	exit(EXIT_FAILURE);
      }
      fputc('\n', fp_anchors);
    }
    if (! anchor_break) {
      fprintf(fp_pgg, "(%d_3,%d_5)\t1\n", (anchor_number - 1), anchor_number);
    }
    num_anchors_written++;
  } else {
    int anchors_in_buffer = stop_pos / OPT_ANCHOR_LEN;
    int anchor_len = OPT_ANCHOR_LEN;
    if ((stop_pos % OPT_ANCHOR_LEN) > 0) {
      anchors_in_buffer++;
      anchor_len = stop_pos / anchors_in_buffer;
      if ((stop_pos % anchors_in_buffer) > 0) {
	anchor_len++;
      }
    }
    for (i = 0; i < anchors_in_buffer; i++) {
      int anchor_stop_pos = (((i + 1) * anchor_len) > stop_pos) ? stop_pos : ((i + 1) * anchor_len);
      fprintf(fp_anchors, ">medoid_%d\n", anchor_number);
      if (core_anchor) {
	fprintf(fp_single, "%d\n", anchor_number);
      }
      fprintf(fp_cluster, "%d\t%d\n", anchor_number, anchor_prevalence);
      fprintf(fp_match, "%d\tmedoid_%d\n", anchor_number, anchor_number);
      for (j = i * anchor_len; j < anchor_stop_pos; j += 60) {
	int out_len = (stop_pos - j) < 60 ? (stop_pos - j) : 60;
	if (fwrite(&contig_seq[j], sizeof(char), (size_t) out_len, fp_anchors) != out_len) {
	  fprintf (stderr, "Could not complete write to anchors file for anchor %d\n", anchor_number);
	  exit(EXIT_FAILURE);
	}
	fputc('\n', fp_anchors);
      }
      if (! anchor_break) {
	fprintf(fp_pgg, "(%d_3,%d_5)\t1\n", (anchor_number - 1), anchor_number);
      } else {
	anchor_break = false;
      }
      anchor_number++;
      num_anchors_written++;
    }
  }
  return num_anchors_written;
}

/*
This subroutine reads through a genome one basepair at a time to construct 23mers and output them to the buclet files
*/
int kmer_bucket_sort_genome(FILE * fp_fasta, char * genome_file_name, uint16_t genome_number, char ** bucket_file_names, struct Kmer ** kmer_buffers, int * kmer_buffer_indices, size_t * bucket_sizes)
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
  bool first_amb = true;
  uint64_t kmer_mask = 01777777777777777; /* This is used to mask only the lower 46 bits used for the 23mer when shifting */
  uint64_t bucket_mask = 01777776000000000; /* This is used to select the upper 18 bits of the lower 46 bits to determine the 9 basepairs used for bucket sorting */
  uint32_t contig_number = 0;
  int kmer_bp_count = 0;
  int contig_pos = 0;
  int kmer_bucket;
  int new_red_files = 0;
  
  if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_fasta)) == -1) {
    fprintf(stderr, "%s appears to be empty.\n", genome_file_name);
    exit(EXIT_FAILURE);
  }
  if (fasta_line[0] != '>') {
    fprintf(stderr, "First line of fasta file %s does not begin with a >!\n%s", genome_file_name, fasta_line);
    exit(EXIT_FAILURE);
  }
  /* fprintf(stderr, "Genome: %s contig %d\n%s", genome_file_name, contig_number, fasta_line); */
  new_red_files++; /* this gets incremented for every contig in the first genome, and only once for each other genome */
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
	first_amb = true;
	break;
      case 'c':
      case 'C':
	cur_bp = 01;
	revc_bp = 04000000000000000;
	first_amb = true;
	break;
      case 'g':
      case 'G':
	cur_bp = 02;
	revc_bp = 02000000000000000;
	first_amb = true;
	break;
      case 't':
      case 'T':
	cur_bp = 03;
	revc_bp = 0;
	first_amb = true;
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
	if (first_amb) {
	  first_amb = false;
	  fprintf(stderr, "Amb %d:%d:%d", (int) genome_number, (int) contig_number, contig_pos);
	}
	break;
      case '>':
	if (prev_char == '\n') {
	  /* Start of a new contig. Read in header line but ignore it other than incrementing contig number */
	  /* fprintf(stderr, "Genome: %s contig %d last pos %d\n", genome_file_name, contig_number, contig_pos); */
	  if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_fasta)) == -1) {
	    fprintf (stderr, "empty contig fasta header: %d for genome %s!\n", contig_number, genome_file_name);
	    exit(EXIT_FAILURE);
	  }
	  new_contig = true;
	  if (contig_number == MAX_UINT32) {
	    fprintf (stderr, "maximum number of contigs: %d exceeded for genome %s!\n", contig_number, genome_file_name);
	    exit(EXIT_FAILURE);
	  }
	  if (genome_number == 0) {
	    new_red_files++; /* this gets incremented for every contig in the first genome, and only once for each other genome */
	  }
	  contig_number++;
	  /* fprintf(stderr, "Genome: %s contig %d\n>%s", genome_file_name, contig_number, fasta_line); */
	} else {
	  fprintf (stderr, "Unexpected > not at beginning of line in fasta file %s!\n", genome_file_name);
	  exit(EXIT_FAILURE);
	}
	break;
      default:
	if (isprint (cur_char)) {
	  fprintf (stderr, "Unexpected character in fasta file %s: %c.\n", genome_file_name, cur_char);
	} else {
	  fprintf (stderr, "Unexpected unprintable character in fasta file %s: \\x%x.\n", genome_file_name, cur_char);
	}
	exit(EXIT_FAILURE);
      }
      if (new_contig) {
	new_contig = false;
	prev_char = '\n';
	kmer_bp_count = 0; 
	contig_pos = 0;
	cur_kmer = 0;
	revc_kmer = 0;
	first_amb = true;
      } else {
	prev_char = cur_char;
	contig_pos++; /* increment the contig position before 23mer determination since we start numbering at 1 in the contig */
	if (!reset_kmer) {
	  /* calculate the current 23mer by shifting on the current basepair for both the forward strand and reverse complement */
	  kmer_bp_count++; /* this is used to determine if the 23mer has been filled up */
	  cur_kmer = ((cur_kmer << 2) & kmer_mask) | cur_bp;
	  revc_kmer = (revc_kmer | revc_bp) >> 2;
	  if (kmer_bp_count >= 23) {
	    /* if ((contig_pos % 1000000) == 0) {
	      fprintf(stderr, "%d\n", contig_pos);
	      } */
	    if (cur_kmer < revc_kmer) { /* use canonical kmer which ever is less */
	      kmer_bucket = (int) ((bucket_mask & cur_kmer) >> 28);
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->kmer = cur_kmer;
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->pos = contig_pos - (KMER_SIZE - 1);
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->genome = genome_number;
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->contig = contig_number;
	    } else {
	      kmer_bucket = (int) ((bucket_mask & revc_kmer) >> 28);
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->kmer = revc_kmer;
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->pos = -(contig_pos - (KMER_SIZE - 1));
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->genome = genome_number;
	      ((*(kmer_buffers + kmer_bucket)) + kmer_buffer_indices[kmer_bucket])->contig = contig_number;
	    }
	    kmer_buffer_indices[kmer_bucket]++;
	    if (kmer_buffer_indices[kmer_bucket] == KMER_BUFFER_LEN) {
	      /* kmer bucket buffer is full so output it and reset index */
	      kmer_buffer_indices[kmer_bucket] = 0;
	      fp_bucket = fopen(bucket_file_names[kmer_bucket], "a");
	      if (fp_bucket == NULL) {
		fprintf (stderr, "Could not open file %s\n", bucket_file_names[kmer_bucket]);
		exit(EXIT_FAILURE);
	      }
	      if (fwrite(kmer_buffers[kmer_bucket], sizeof(struct Kmer), KMER_BUFFER_LEN, fp_bucket) != KMER_BUFFER_LEN) {
		fprintf (stderr, "Could not complete write to  file %s\n", bucket_file_names[kmer_bucket]);
		exit(EXIT_FAILURE);
	      }
	      bucket_sizes[kmer_bucket] += KMER_BUFFER_LEN;
	      fclose(fp_bucket);
	    }
	  }
	} else {
	  reset_kmer = false;
	  kmer_bp_count = 0;
	  cur_kmer = 0;
	  revc_kmer = 0;
	}
      }
    }
  }
  /* fprintf(stderr, "Genome: %s contig %d last pos %d\n", genome_file_name, contig_number, contig_pos); */
  free((void *) fasta_line);
  return (new_red_files);
}

int
main (int argc, char **argv)
{
  char * genomes_file = NULL;
  int num_red_files = 0; /* this is the number of files/buckets for the reduced k-mers, one for each of the contigs in the first genome, and one per genome after that */
  int num_first_genome_contigs; /* this is the number of contigs in the first genome */
  int red_file_number; /* index of reduced k-mer file */
  char red_file_name[100]; /* used to convert the red_file_number into a file name */
  int index;
  int getopt_return;
  FILE * fp_file_names;
  FILE * fp_bucket;
  FILE * fp_red_bucket;
  char * file_name_line = NULL;
  size_t file_name_line_malloc_len = 0;
  size_t getline_return;
  uint16_t genome_number = 0; /* this gets assigned and used to distinguish different genomes for each k-mer */
  bool max_genomes_exceeded = false; /* set if too many genomes for a uint16_t (65,536) */
  struct Kmer ** kmer_buffers; /* this is the malloced array of arrays for each kmer buffer in the top level buckets used to buffer for output */
  struct Kmer * kmer_buffers_pool; /* this is the pool of malloc storage to use for the k-mer buffer arrays */
  struct Kmer * kmer_array; /* this is the pool of malloc storage to use for the k-mer array containing one bucket of k-mers */
  int * kmer_buffer_indices; /* malloced array of indices into the kmer buffers */
  struct red_Kmer ** red_kmer_buffers; /* this is the malloced array of arrays for each kmer buffer in the top level buckets used to buffer for output */
  struct red_Kmer * red_kmer_buffers_pool; /* this is the pool of malloc storage to use for the k-mer buffer arrays */
  struct red_Kmer * red_kmer_array; /* this is the pool of malloc storage to use for the reduced k-mer array containing one bucket of reduced k-mers */
  int * red_kmer_buffer_indices; /* malloced array of indices into the kmer buffers */
  size_t * bucket_sizes; /* malloced array of bucket sizes */
  size_t * red_bucket_sizes; /* malloced array of reduced k-mer bucket sizes */
  char bps[4] = {'A','C','G','T'};
  char ** bucket_file_names; /* array of malloced fixed length strings for bucket file names */
  char * bucket_file_names_pool; /* this is the pool pf malloc storage for the bucket file names */
  char tmp_file_dir_name[12];
  int i,j,k,l,m,n,o,p,q;
  FILE * fp_file_name;
  FILE * fp_anchors;
  FILE * fp_pgg;
  FILE * fp_single;
  FILE * fp_cluster;
  FILE * fp_match;
  int file_name_len;
  char * fasta_line = NULL;
  size_t fasta_line_malloc_len = 0;
  char prev_char, cur_car;
  int cur_file_genome, cur_genome;
  int anchor_number = 1; /* the anchor/medoid number initialized here to 1 */
  int kmer_threshold = 1; /* the minimum number of genomes a k-mer must be present in to be used as an anchor */
  int core_threshold = 50; /* percentage of genomes needed to be core */
  
  opterr = 0;

  while ((getopt_return = getopt (argc, argv, "g:t:c:")) != -1) {
    switch (getopt_return) {
    case 'g':
      genomes_file = optarg;
      break;
    case 't':
      kmer_threshold = atoi(optarg);
      if (kmer_threshold <= 0) {
	fprintf(stderr, "k-mer threshold must be greater than 0 not: %s\n", optarg);
	exit(EXIT_FAILURE);
      }
      break;
    case 'c':
      core_threshold = atoi(optarg);
      if ((core_threshold < 0) || (core_threshold > 100)){
	fprintf(stderr, "core threshold must be >= 0 and <= 100, not: %s\n", optarg);
	exit(EXIT_FAILURE);
      }
      break;
    case '?':
      if ((optopt == 'g') || (optopt == 't') || (optopt == 'c')) {
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      } else if (isprint (optopt)) {
	fprintf (stderr, "Unknown option -%c.\n", optopt);
      } else {
	fprintf (stderr, "Unknown option character \\x%x.\n", optopt);
      }
      exit(EXIT_FAILURE);
    default:
      exit(EXIT_FAILURE);
    }
  }

  fprintf (stderr, "Genomes File = %s\nK-mer threshold = %d\n", genomes_file, kmer_threshold);

  for (index = optind; index < argc; index++) {
    fprintf (stderr, "Non-option argument %s\n", argv[index]);
  }


  fprintf (stderr, "Generating directories to store k-mer %d buckets\n", KMER_BUFFER_NUMBER);

  /* generate and store directories and files names for the 262,144 buckets to sort 23mers into */
  /* this is the array of bucket files names */
  bucket_file_names = (char **) malloc((size_t) (KMER_BUFFER_NUMBER * sizeof(char *)));
  if (bucket_file_names  == NULL) {
      fprintf (stderr, "Could not allocate memory for bucket_file_names\n");
    exit(EXIT_FAILURE);
  }

  /* this is the pool of storage for the bucket file names - relies on the file names being 12 characters including the null '\0' character */
  bucket_file_names_pool = (char *) malloc((size_t) (KMER_BUFFER_NUMBER * 12 * sizeof(char)));
  if (bucket_file_names_pool == NULL) {
    fprintf (stderr, "Could not allocate memory for bucket_file_names_pool\n");
    exit(EXIT_FAILURE);
  }

  /* assign the correct storage pool area for each bucket file name */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
   bucket_file_names[index] = bucket_file_names_pool + (index * 12);
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
	      tmp_file_dir_name[6] = bps[n];
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

  fprintf (stderr, "Allocating output buffers for k-mer buckets\n");

  /* allocate the array of indices for each kmer bucket size */
  bucket_sizes = (size_t *) malloc((size_t) (KMER_BUFFER_NUMBER * sizeof(size_t)));
  if (bucket_sizes == NULL) {
    fprintf (stderr, "Could not allocate memory for bucket_sizes\n");
    exit(EXIT_FAILURE);
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
    bucket_sizes[index] = 0;
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
    kmer_buffers[index] = kmer_buffers_pool + (index * KMER_BUFFER_LEN);
  }
  
  fp_file_names = fopen(genomes_file, "r");
  if (fp_file_names == NULL) {
    fprintf (stderr, "Could not open file %s\n", genomes_file);
    exit(EXIT_FAILURE);
  }


  fprintf (stderr, "Reading genome file names from %s\n", genomes_file);

  /* loop through each genome in the genomes file input file */
  while ((getline_return = getline(&file_name_line, &file_name_line_malloc_len, fp_file_names)) != -1) {
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

    fprintf (stderr, "Reading genome from %s\n", file_name_line);

    num_red_files += kmer_bucket_sort_genome(fp_file_name, file_name_line, genome_number, bucket_file_names, kmer_buffers, kmer_buffer_indices, bucket_sizes);
    if (genome_number == 0) {
      num_first_genome_contigs = num_red_files;
    }
    if (genome_number == MAX_UINT16) {
      max_genomes_exceeded = true;
    } else {
      genome_number++;
    }
    fclose(fp_file_name);
  }
  /* free((void *) file_name_line);*/

  fclose(fp_file_names);

  fprintf (stderr, "Flushing k-mer bucket output buffers\n");

  /* flush any 23mers remaining in the kmer buffers */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
    if (kmer_buffer_indices[index] != 0) {
      fp_bucket = fopen(bucket_file_names[index], "a");
      if (fp_bucket == NULL) {
	fprintf (stderr, "Could not open file %s\n", bucket_file_names[index]);
	exit(EXIT_FAILURE);
      }
      if (fwrite((void *) kmer_buffers[index], sizeof(struct Kmer), kmer_buffer_indices[index], fp_bucket) != kmer_buffer_indices[index]) {
	fprintf (stderr, "Could not complete write to file %s\n", bucket_file_names[index]);
	exit(EXIT_FAILURE);
      }
      bucket_sizes[index] += kmer_buffer_indices[index];
      fclose(fp_bucket);
    }
  }

  fprintf (stderr, "Freeing k-mer bucket buffers\n");

  /* free memory used for bucket file write buffers */
  free ((void *) kmer_buffer_indices);
  free ((void *) kmer_buffers);
  free ((void *) kmer_buffers_pool);

  if (kmer_threshold > genome_number) {
    fprintf(stderr, "k-mer threshold (%d) must not be greater than the number of genomes (%d)\n", kmer_threshold, genome_number);
    exit(EXIT_FAILURE);
  }
  
  fprintf (stderr, "Allocating output buffers for reduced k-mer buckets\n");

  /* allocate the array of indices for each reduced kmer bucket buffer */
  red_kmer_buffer_indices = (int *) malloc((size_t) (num_red_files * sizeof(int)));
  if (red_kmer_buffer_indices == NULL) {
    fprintf (stderr, "Could not allocate memory for red_kmer_buffer_indices\n");
    exit(EXIT_FAILURE);
  }

  /* allocate the array of indices for each reducedkmer bucket size */
  red_bucket_sizes = (size_t *) malloc((size_t) (num_red_files * sizeof(size_t)));
  if (red_bucket_sizes == NULL) {
    fprintf (stderr, "Could not allocate memory for reduced k-mer bucket_sizes\n");
    exit(EXIT_FAILURE);
  }

  /* initialize the kmer bucket buffer indices to 0 */
  for (index = 0; index < num_red_files; index++) {
    red_kmer_buffer_indices[index] = 0;
    red_bucket_sizes[index] = 0;
  }

  /* allocate the array of pointers to the kmer buffers */
  red_kmer_buffers = (struct red_Kmer **) malloc((size_t) (num_red_files * sizeof(struct red_Kmer *)));
  if (red_kmer_buffers == NULL) {
    fprintf (stderr, "Could not allocate memory for red_kmer_buffers\n");
    exit(EXIT_FAILURE);
  }

  /* allocate the pool of storage for the kmer buffers */
  red_kmer_buffers_pool = (struct red_Kmer *) malloc((size_t) (num_red_files * KMER_BUFFER_LEN * sizeof(struct red_Kmer)));
  if (red_kmer_buffers_pool == NULL) {
    fprintf (stderr, "Could not allocate memory for red_kmer_buffers_pool\n");
    exit(EXIT_FAILURE);
  }

  /* set up the kmer buffers to point to the assigned storage buffer within the storage pool */
  for (index = 0; index < num_red_files; index++) {
    red_kmer_buffers[index] = red_kmer_buffers_pool + (index * KMER_BUFFER_LEN);
  }
  

  fprintf (stderr, "Sort each k-mer bucket and output to %d reduced buckets files\n", num_red_files);

  /* read in one bucket at a time and sort it */
  for (index = 0; index < KMER_BUFFER_NUMBER; index++) {
    kmer_array = (struct Kmer *) malloc((size_t) (bucket_sizes[index] * sizeof(struct Kmer)));
    if (kmer_array == NULL) {
      fprintf (stderr, "Could not allocate memory for kmer_array\n");
      exit(EXIT_FAILURE);
    }
    fp_bucket = fopen(bucket_file_names[index], "r");
    if (fp_bucket == NULL) {
      /* fprintf (stderr, "WARNING: Could not open file %s for reading\n", bucket_file_names[index]); */
      continue; /* some buckets may not have any data */
      /* exit(EXIT_FAILURE);*/
    }
    if (fread((void *) kmer_array, sizeof(struct Kmer), bucket_sizes[index], fp_bucket) != bucket_sizes[index]) {
      fprintf (stderr, "Could not complete read of file %s\n", bucket_file_names[index]);
      exit(EXIT_FAILURE);
    }
    fclose(fp_bucket);
    unlink(bucket_file_names[index]);
    qsort((void *) kmer_array, bucket_sizes[index], sizeof(struct Kmer), kmer_sort);
    for (i = 0; i <  bucket_sizes[index];) {
      struct Kmer first_array_kmer, prev_array_kmer;
      bool duplicate_kmer = false;
      int kmer_count = 1;
      first_array_kmer = kmer_array[i];
      prev_array_kmer = first_array_kmer;
      for (i++; first_array_kmer.kmer == kmer_array[i].kmer; i++) {
	if (prev_array_kmer.genome == kmer_array[i].genome) {
	  duplicate_kmer = true;
	} else {
	  kmer_count++;
	}
	prev_array_kmer = kmer_array[i];
      }
      if ((! duplicate_kmer) && (kmer_count >= kmer_threshold)) {
	if (first_array_kmer.genome == 0) {
	  red_file_number = first_array_kmer.contig;
	} else {
	  red_file_number = num_first_genome_contigs + (first_array_kmer.genome - 1);
	}
	if (first_array_kmer.pos < 0) {
	  ((*(red_kmer_buffers + red_file_number)) + red_kmer_buffer_indices[red_file_number])->pos = - first_array_kmer.pos;
	} else {
	  ((*(red_kmer_buffers + red_file_number)) + red_kmer_buffer_indices[red_file_number])->pos = first_array_kmer.pos;
	}
	((*(red_kmer_buffers + red_file_number)) + red_kmer_buffer_indices[red_file_number])->genome = first_array_kmer.genome;
	((*(red_kmer_buffers + red_file_number)) + red_kmer_buffer_indices[red_file_number])->contig = first_array_kmer.contig;
	((*(red_kmer_buffers + red_file_number)) + red_kmer_buffer_indices[red_file_number])->prevalence = (uint16_t) kmer_count;
	red_kmer_buffer_indices[red_file_number]++;
	if (red_kmer_buffer_indices[red_file_number] == KMER_BUFFER_LEN) {
	  /* kmer bucket buffer is full so output it and reset index */
	  red_kmer_buffer_indices[red_file_number] = 0;
	  sprintf(red_file_name, "%i", red_file_number);
	  fp_red_bucket = fopen(red_file_name, "a");
	  if (fp_red_bucket == NULL) {
	    fprintf (stderr, "Could not open file %s\n", red_file_name);
	    exit(EXIT_FAILURE);
	  }
	  if (fwrite(red_kmer_buffers[red_file_number], sizeof(struct red_Kmer), KMER_BUFFER_LEN, fp_red_bucket) != KMER_BUFFER_LEN) {
	    fprintf (stderr, "Could not complete write to  file %s\n", red_file_name);
	    exit(EXIT_FAILURE);
	  }
	  red_bucket_sizes[red_file_number] += KMER_BUFFER_LEN;
	  fclose(fp_red_bucket);
	}
      }
    }
    free ((void *) kmer_array);
  }

  free ((void *) bucket_sizes);
  
  fprintf (stderr, "Flushing reduced k-mer bucket output buffers\n");

  /* flush any reduced k-mers remaining in the reduced kmer buffers */
  for (index = 0; index < num_red_files; index++) {
    if (red_kmer_buffer_indices[index] != 0) {
      sprintf(red_file_name, "%i", index);
      fp_red_bucket = fopen(red_file_name, "a");
      if (fp_red_bucket == NULL) {
	fprintf (stderr, "Could not open file %s\n", red_file_name);
	exit(EXIT_FAILURE);
      }
      if (fwrite((void *) red_kmer_buffers[index], sizeof(struct red_Kmer), red_kmer_buffer_indices[index], fp_red_bucket) != red_kmer_buffer_indices[index]) {
	fprintf (stderr, "Could not complete write to file %s\n", red_file_name);
	exit(EXIT_FAILURE);
      }
      red_bucket_sizes[index] += red_kmer_buffer_indices[index];
      fclose(fp_red_bucket);
    }

    fprintf (stderr, "%lu reduced k-mers in %d\n", red_bucket_sizes[index], index);
  }

  fprintf (stderr, "Reading genome file names from %s\n", genomes_file);

  /* loop through each genome in the genomes file input file as needed to determine anchors*/
  cur_file_genome = 0;

  fp_file_names = fopen(genomes_file, "r");
  if (fp_file_names == NULL) {
    fprintf (stderr, "Could not open file %s\n", genomes_file);
    exit(EXIT_FAILURE);
  }

  getline_return = getline(&file_name_line, &file_name_line_malloc_len, fp_file_names);
  if (getline_return == -1) {
    fprintf (stderr, "Unexpected end of genome names file: %s!\n", genomes_file);
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

  fprintf (stderr, "Reading genome %d from %s\n", cur_file_genome, file_name_line);

  fp_file_name = fopen(file_name_line, "r");
  if (fp_file_name == NULL) {
    fprintf (stderr, "Could not open file %s\n", file_name_line);
    exit(EXIT_FAILURE);
  }
  if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_file_name)) == -1) {
    fprintf(stderr, "%s appears to be empty.\n", file_name_line);
    exit(EXIT_FAILURE);
  }
  if (fasta_line[0] != '>') {
    fprintf(stderr, "First line of fasta file %s does not begin with a >!\n%s", file_name_line, fasta_line);
    exit(EXIT_FAILURE);
  }
  /* fprintf (stderr, "contig fasta header: 0 for genome %s\n%s", file_name_line, fasta_line); */

  fprintf (stderr, "Opening medoids.fasta file for anchors and pgg.txt for PGG\n");

  fp_anchors = fopen("medoids.fasta", "w");
  if (fp_anchors == NULL) {
    fprintf (stderr, "Could not open anchors/medoids file medoids.fasta for output\n");
    exit(EXIT_FAILURE);
  }

  fp_pgg = fopen("pgg.txt", "w");
  if (fp_pgg == NULL) {
    fprintf (stderr, "Could not open pgg.txt file for output\n");
    exit(EXIT_FAILURE);
  }

  fp_single = fopen("single_copy_clusters.txt", "w");
  if (fp_single == NULL) {
    fprintf (stderr, "Could not open single_copy_clusters.txt file for output\n");
    exit(EXIT_FAILURE);
  }

  fp_cluster = fopen("cluster_sizes.txt", "w");
  if (fp_cluster == NULL) {
    fprintf (stderr, "Could not open cluster_sizes.txt file for output\n");
    exit(EXIT_FAILURE);
  }

  fp_match = fopen("matchtable.txt", "w");
  if (fp_match == NULL) {
    fprintf (stderr, "Could not open matchtable.txt file for output\n");
    exit(EXIT_FAILURE);
  }

  /* read in reduced k-mers from contig/genome buckets and produce anchors */
  for (index = 0; index < num_red_files; index++) {
    char contig_seq[CONTIG_SEQ_BUFFER_LEN];
    int prev_pos, cur_pos, prev_contig, cur_contig, contig_seq_pos, num_anchors_written;
    bool anchor_break = true; /* this is true if the previous anchor is from a different contig or genome or there is no previous anchor */
    float prev_prevalence;
    int anchor_prevalence;
    bool core_anchor = false;
    int first_red_kmer;
    int last_red_kmer;
    int first_anchor_pos;
    int last_anchor_pos;
    bool reset_prevalence;
    int num_outliers;
    int num_so_far;
    
    sprintf(red_file_name, "%i", index);

    fprintf (stderr, "Reading %lu reduced k-mers from %s\n", red_bucket_sizes[index], red_file_name);

    fp_red_bucket = fopen(red_file_name, "r");
    if (fp_red_bucket == NULL) {
      /* fprintf (stderr, "WARNING: Could not open file %s\n", red_file_name); */
      continue; /* Some reduced buckets may not have any data */
      /* exit(EXIT_FAILURE);*/
    }
    red_kmer_array = (struct red_Kmer *) malloc((size_t) (red_bucket_sizes[index] * sizeof(struct red_Kmer)));
    if (red_kmer_array == NULL) {
      fprintf (stderr, "Could not allocate memory for red_kmer_array %lu\n", (red_bucket_sizes[index] * sizeof(struct red_Kmer)));
      exit(EXIT_FAILURE);
    }
    if (fread((void *) red_kmer_array, sizeof(struct red_Kmer), red_bucket_sizes[index], fp_red_bucket) != red_bucket_sizes[index]) {
      fprintf (stderr, "Could not complete read of file %s\n", red_file_name);
      exit(EXIT_FAILURE);
    }
    qsort((void *) red_kmer_array, red_bucket_sizes[index], sizeof(struct red_Kmer), red_kmer_sort);

    /* first smooth out dips in prevalence due to presumed variability changing k-mers */
    
    prev_pos = -1;
    prev_contig = -1;
    first_red_kmer = 0;
    last_red_kmer = 0;
    reset_prevalence = false;
    num_outliers = 0;
    num_so_far = 0;
    prev_prevalence = 0.0;
    for (i = 0; i < red_bucket_sizes[index]; i++) {
      cur_pos = (int) red_kmer_array[i].pos;
      cur_contig = (int) red_kmer_array[i].contig;
      cur_genome = (int) red_kmer_array[i].genome;
      anchor_prevalence = (int) red_kmer_array[i].prevalence;
      /* fprintf(stderr, "%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%f\n", cur_pos, cur_contig, cur_genome, anchor_prevalence, first_red_kmer, last_red_kmer, num_outliers, num_so_far, prev_pos, prev_contig, prev_prevalence); */
      if (num_so_far == 0) {
	/* start of first contig or prevalence was reset */
	num_so_far = 1;
	num_outliers = 0;
	first_red_kmer = i;
	last_red_kmer = i;
	prev_prevalence = (float) anchor_prevalence;
      } else if (cur_contig != prev_contig) {
	/* start of a new contig so reset prevalence for the old contig */
	reset_prevalence = true;
      } else if ((cur_pos - prev_pos) > KMER_SIZE) {
	/* Break in unique k-mers */
	reset_prevalence = true;
      } else if ((anchor_prevalence > ((110 * prev_prevalence) / 100)) || ((anchor_prevalence < ((90 * prev_prevalence) / 100)))) {
	/* Significant change in anchor prevalence */
	num_outliers++;
	if (num_outliers > (2 * KMER_SIZE)) {
	  reset_prevalence = true;
	}
      } else {
	num_outliers = 0;
	last_red_kmer = i;
	prev_prevalence = ((float)((num_so_far * prev_prevalence) + anchor_prevalence)) / ((float)(num_so_far + 1));
	num_so_far++;
      }
      if (reset_prevalence) {
	int max_prevalence = 0;
	/* fprintf(stderr, "reset:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%f\n", cur_pos, cur_contig, cur_genome, anchor_prevalence, first_red_kmer, last_red_kmer, num_outliers, num_so_far, prev_pos, prev_contig, prev_prevalence); */
	for (j = first_red_kmer; j <= last_red_kmer; j++) {
	  if (red_kmer_array[j].prevalence > max_prevalence) {
	    max_prevalence = (int) red_kmer_array[j].prevalence;
	  }
	}
	for (j = first_red_kmer; j <= last_red_kmer; j++) {
	  red_kmer_array[j].prevalence = (uint16_t) max_prevalence;
	}
	reset_prevalence = false;
	num_so_far = 0;
	num_outliers = 0;
	i = last_red_kmer;
	cur_pos = (int) red_kmer_array[i].pos;
	cur_contig = (int) red_kmer_array[i].contig;
      }
      prev_pos = cur_pos;
      prev_contig = cur_contig;
    }
    { /* reset prevalence at the end */
      int max_prevalence = 0;
      /* fprintf(stderr, "last%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%f\n", cur_pos, cur_contig, cur_genome, anchor_prevalence, first_red_kmer, last_red_kmer, num_outliers, num_so_far, prev_pos, prev_contig, prev_prevalence); */
      for (j = first_red_kmer; j <= last_red_kmer; j++) {
	if (red_kmer_array[j].prevalence > max_prevalence) {
	  max_prevalence = (int) red_kmer_array[j].prevalence;
	}
      }
      for (j = first_red_kmer; j <= last_red_kmer; j++) {
	red_kmer_array[j].prevalence = (uint16_t) max_prevalence;
      }
    }

    /* determine and write out anchors/medoids and associated files */
    prev_pos = -1;
    prev_contig = -1;
    contig_seq_pos = 0;
    prev_prevalence = 0;
    first_anchor_pos = 0;
    last_anchor_pos = -1;
    for (i = 0; i <  red_bucket_sizes[index]; i++) {
      if (((last_anchor_pos - first_anchor_pos) + 1) >= CONTIG_SEQ_BUFFER_LEN) {
	/* Buffer full - output first half of current anchors buffer */
	contig_seq_pos = read_fasta_kmer(file_name_line, fp_file_name, cur_contig, first_anchor_pos, (CONTIG_SEQ_BUFFER_LEN / 2), contig_seq, contig_seq_pos);
	/* having linker issues for libm.a (-lm) with gcc on some platforms so using crude approximation for ceilf function */
	anchor_number += write_anchors((int)(prev_prevalence + 0.999999999999), fp_single, fp_cluster, fp_match, fp_pgg, fp_anchors, anchor_number, contig_seq, contig_seq_pos, anchor_break, core_anchor);
	anchor_break = false;
	/* Reset first_anchor_pos to reflect the portion of the anchor written */
	first_anchor_pos += (CONTIG_SEQ_BUFFER_LEN / 2);
	contig_seq_pos = 0;
      }
      cur_pos = (int) red_kmer_array[i].pos;
      cur_contig = (int) red_kmer_array[i].contig;
      cur_genome = (int) red_kmer_array[i].genome;
      anchor_prevalence = (int) red_kmer_array[i].prevalence;
      /* if (((cur_genome == 0) && (cur_contig == 0) && (cur_pos < 210000)) || ((cur_genome == 0) && (cur_contig == 0) && (cur_pos > 4000000)) || ((cur_genome == 0) && (cur_contig == 1))) {
	fprintf(stderr, "%d:%d:%d:%d:%f:%d:%d\n", cur_pos, cur_contig, cur_genome, anchor_prevalence, prev_prevalence, first_anchor_pos, last_anchor_pos);
	} */
      if (cur_genome != cur_file_genome) {
	fprintf (stderr, "Current genome number %d does not match stored genome number %d\n", cur_genome, cur_file_genome);
	exit(EXIT_FAILURE);
      }
      if (cur_contig != prev_contig) {
	if ((prev_contig != -1) && (((last_anchor_pos - first_anchor_pos) + 1) >= MIN_ANCHOR_LEN)) {
	  /* Switched contigs - output current anchors buffer */
	  contig_seq_pos = read_fasta_kmer(file_name_line, fp_file_name, cur_contig, first_anchor_pos, ((last_anchor_pos - first_anchor_pos) + 1), contig_seq, contig_seq_pos);
	  anchor_number += write_anchors((int)(prev_prevalence + 0.999999999999), fp_single, fp_cluster, fp_match, fp_pgg, fp_anchors, anchor_number, contig_seq, contig_seq_pos, anchor_break, core_anchor);
	}
	anchor_break = true;
	contig_seq_pos = 0;
	first_anchor_pos = cur_pos;
	last_anchor_pos = cur_pos + (KMER_SIZE - 1);
      } else {
	if (((cur_pos - prev_pos) > KMER_SIZE) || ((anchor_prevalence > ((110 * prev_prevalence) / 100)) || ((anchor_prevalence < ((90 * prev_prevalence) / 100))))) {
	  /* Break in unique k-mers or significant change in anchor prevalence - output current anchors buffer */
	  if (((last_anchor_pos - first_anchor_pos) + 1) >= MIN_ANCHOR_LEN) {
	    contig_seq_pos = read_fasta_kmer(file_name_line, fp_file_name, cur_contig, first_anchor_pos, ((last_anchor_pos - first_anchor_pos) + 1), contig_seq, contig_seq_pos);
	    num_anchors_written = write_anchors((int)(prev_prevalence + 0.999999999999), fp_single, fp_cluster, fp_match, fp_pgg, fp_anchors, anchor_number, contig_seq, contig_seq_pos, anchor_break, core_anchor);
	    anchor_number += num_anchors_written;
	    if (num_anchors_written) {
	      anchor_break = false;
	    }
	  }
	  contig_seq_pos = 0;
	  first_anchor_pos = cur_pos;
	  last_anchor_pos = cur_pos + (KMER_SIZE - 1);
	} else {
	  last_anchor_pos = cur_pos + (KMER_SIZE - 1);
	}
      }
      prev_pos = cur_pos;
      prev_contig = cur_contig;
      prev_prevalence = ((float)(((((last_anchor_pos - first_anchor_pos) + 1) - KMER_SIZE) * prev_prevalence) + anchor_prevalence)) / ((float)((((last_anchor_pos - first_anchor_pos) + 1) - KMER_SIZE) + 1));
      if (prev_prevalence >= ((core_threshold * genome_number) / 100)) {
	core_anchor = true;
      } else {
	core_anchor = false;
      }
    }
    if (((last_anchor_pos - first_anchor_pos) + 1) >= MIN_ANCHOR_LEN) {
      /* Flush remaining anchors buffer */
      contig_seq_pos = read_fasta_kmer(file_name_line, fp_file_name, cur_contig, first_anchor_pos, ((last_anchor_pos - first_anchor_pos) + 1), contig_seq, contig_seq_pos);
      anchor_number += write_anchors((int)(prev_prevalence + 0.999999999999), fp_single, fp_cluster, fp_match, fp_pgg, fp_anchors, anchor_number, contig_seq, contig_seq_pos, anchor_break, core_anchor);
    }
    free ((void *) red_kmer_array);
    fclose(fp_red_bucket);
    unlink(red_file_name);

    /* check if we need to switch to a new genome or just a new contig in first genome */
    if ((index < (num_first_genome_contigs - 1)) || (index == (num_red_files - 1))) {
    } else {
      getline_return = getline(&file_name_line, &file_name_line_malloc_len, fp_file_names);
      if (getline_return == -1) {
	fprintf (stderr, "Unexpected end of genome names file: %s!\n", genomes_file);
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
      fclose(fp_file_name);

      cur_file_genome++;
      fprintf (stderr, "Reading genome %d from %s\n", cur_file_genome, file_name_line);

      fp_file_name = fopen(file_name_line, "r");
      if (fp_file_name == NULL) {
	fprintf (stderr, "Could not open file %s\n", file_name_line);
	exit(EXIT_FAILURE);
      }
      if ((getline_return = getline(&fasta_line, &fasta_line_malloc_len, fp_file_name)) == -1) {
	fprintf(stderr, "%s appears to be empty.\n", file_name_line);
	exit(EXIT_FAILURE);
      }
      if (fasta_line[0] != '>') {
	fprintf(stderr, "First line of fasta file %s does not begin with a >!\n%s", file_name_line, fasta_line);
	exit(EXIT_FAILURE);
      }
      /* fprintf (stderr, "contig fasta header: 0 for genome %s\n%s", file_name_line, fasta_line); */
    }
  }

  free((void *) fasta_line);
  free((void *) file_name_line);

  fclose(fp_file_name);
  fclose(fp_file_names);
  fclose(fp_anchors);
  fclose(fp_pgg);
  fclose(fp_single);
  fclose(fp_cluster);
  fclose(fp_match);

  fprintf (stderr, "Freeing reduced k-mer bucket buffers\n");

  /* free memory used for bucket file write buffers */
  free ((void *) red_kmer_buffer_indices);
  free ((void *) red_kmer_buffers);
  free ((void *) red_kmer_buffers_pool);
  free ((void *) red_bucket_sizes);

  /* remove two level directory structure based on 3mers */
  for (i = 0; i < 4; i++) {
    tmp_file_dir_name[0] = bps[i];
    for (j = 0; j < 4; j++) {
      tmp_file_dir_name[1] = bps[j];
      for (k = 0; k < 4; k++) {
	tmp_file_dir_name[2] = bps[k];
	tmp_file_dir_name[3] = '/';
	for (l = 0; l < 4; l++) {
	  tmp_file_dir_name[4] = bps[l];
	  for (m = 0; m < 4; m++) {
	    tmp_file_dir_name[5] = bps[m];
	    for (n = 0; n < 4; n++) {
	      tmp_file_dir_name[6] = bps[n];
	      tmp_file_dir_name[7] = '\0';
	      if (rmdir(tmp_file_dir_name)) {
		fprintf (stderr, "Could not remove directory %s\n", tmp_file_dir_name);
		exit(EXIT_FAILURE);
	      }
	    }
	  }
	}
	tmp_file_dir_name[3] = '\0';
	if (rmdir(tmp_file_dir_name)) {
	  fprintf (stderr, "Could not remove directory %s\n", tmp_file_dir_name);
	  exit(EXIT_FAILURE);
	}
      }
    }
  }

  exit(EXIT_SUCCESS);
}
