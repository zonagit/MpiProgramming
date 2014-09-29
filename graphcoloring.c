#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

char * INPUT_PATH = "/uufs/chpc.utah.edu/common/home/u0082100/cs5965/assignment1/graph_files/";
char * OUTPUT_PATH = "/uufs/chpc.utah.edu/common/home/u0082100/CS6965/u0082100/assignment1/";
char * ws;//weak/strong indicator
int rank,npes, root =0; 
int V,E;//number of vertices and edges

void get_num_v_e(char *filename)
{
  FILE* file = fopen(filename,"rb");
  if (!file) {
    printf("Unable to open file %s\n",filename);
    return;
  }
  char line[1000];
  char *token;
  
  //read file line by line
  while (fgets(line,1000,file) != NULL) {
    //tokenize line
    // printf("Full line %s\n",line);
    strtok(line," ");
    
    //read number of vertices and edges
    if (strcmp(line,"p")==0) {
      token = (char *)strtok(NULL," ");
      token = (char *)strtok(NULL," ");
      V = atoi(token);
      token = (char *)strtok(NULL," ");
      E = atoi(token);
      return;
    }
  }

  fclose(file);
} 

int main(int argc,char** argv)
{
  int num_v_per_p, remainder_v_per_p; 
  int first_v, last_v,range;//ids of the first and last vertices in
  //a given processor
  //Initialize
  MPI_Init(&argc,&argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  if (argc != 4) { 
    if (rank == root)
      printf("Usage: mpirun/mpiexec -np xx ./graphcoloring weak/strong input_filename output_filename\n");

    MPI_Finalize();
    return -1;
  }

  //Read whether we are doing strong/weak scaling and the input and output
  //filenames
  ws = argv[1];

  if (strcmp(ws,"strong")!=0 && strcmp(ws,"weak")!=0) {
    if (rank == root)
      printf("Invalid option. Please enter strong or weak for the second command line option.\n");

    MPI_Finalize();
    return -1;
  }   
  char *input_filename = malloc(900);
  strcpy(input_filename,INPUT_PATH);
  strcat(input_filename,argv[2]);

  char *output_filename = malloc(900);
  strcpy(output_filename,OUTPUT_PATH);
  strcat(output_filename,argv[3]);
  // printf("%s %s\n",input_filename,output_filename);
  
  //all processes read the number of vertices and edges
  //a bit inefficient but...
  get_num_v_e(input_filename);

  //Distribute the number of vertices as uniformly as possible among the
  //processors: that means (V+rank)/npes vertices per process (integer
  //division)
  num_v_per_p = V/npes;
  remainder_v_per_p = (V+rank) % npes;
  first_v =  num_v_per_p * rank + 1 + remainder_v_per_p * (remainder_v_per_p<rank); 
  last_v = (rank+1)* num_v_per_p + (remainder_v_per_p+1) * (remainder_v_per_p<rank);
  range = last_v-first_v+1;
  printf("V=%d rank=%d first=%d last=%d range=%d ideal_range=%d\n",V,rank,first_v,last_v,range,(V+rank)/npes);
  free(input_filename);
  free(output_filename);
  //Finalize
  MPI_Finalize();
  return 0;
}
