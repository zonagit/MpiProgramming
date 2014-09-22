#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#define DEBUG_INITIAL_SORT 0
#define DEBUG_SPLITTERS 0
#define DEBUG_BUCKETS 0
#define DEBUG_FINAL_SORT 1
#define USE_SMALL_ARRAY 1

char * ws;//strong or weak scaling indicator
int rank,npes, root =0; 

//Generate a random number within the range min-max
int rand_interval(int min,int max){
  int r;
  const int range = max - min;
  const int buckets = RAND_MAX/range;
  const int limit = buckets * range;

  do {
    r = rand();
  } while (r >= limit);

  return min + (r/buckets);
}

int compare (const void *a, const void *b)
{
  long int la = *(const long int*) a;
  long int lb = *(const long int*)b;

  return (la>lb)-(la<lb);
}

/****************************************************
 * Sample sort algorithm:
 * 1. Randomly partition input of n elements into n/p groups
 * where p is the number of processors and assign each group
 * to one processor.
 * 2. Sort each group in its own processor using n/p logn/p algo (qsort)
 * 3. Select p-1 splitters (one per process)
 * 4. Gather the splitters in processor p0
 * 5.Sort splitters in p0 and create buckets
 * 6. Exchange data
 * 7. Sort
 */

void samplesort(long int *a,size_t count,size_t size,int (compare)(const void *,const void *),MPI_Comm comm)
{
  long int i,j,k;
  long int *local_splitters,*all_splitters,*global_splitters;
  int *send_counts,*recv_counts,*idx;
  int *send_offsets,*recv_offsets;
  long int * buckets;
#if DEBUG_INITIAL_SORT
  for (i=0;i<count;i++) {
    printf("Initial elements per processor:Rank  %d %lu\n",rank,a[i]);
  }

  MPI_Barrier(comm);
#endif
  //2. Sort locally
  qsort(a,count,size,compare);
  
#if DEBUG_INITIAL_SORT
  for (i=0;i<count;i++) {
    printf("Initial elements per processor after local sort:Rank  %d %lu\n",rank,a[i]);
  }

  MPI_Barrier(comm);
#endif
    //3. Each processor selects npes-1 equally spaced elements from its
  //own list. Equal npes-1 spaces mean psize/npes=n/(npes*npes). We are
  //assuming psize/npes>1 which is true
  local_splitters = (long int *)malloc ((npes-1)*size);
  for (i=0;i<(npes-1);i++) {
    local_splitters[i] = a[(i+1)*count/npes];
  }

#if DEBUG_SPLITTERS
  for (i=0;i<(npes-1);i++) {
    printf("local splitters: Rank %d %lu\n",rank,local_splitters[i]);
  }
  MPI_Barrier(comm);
#endif
  //4a. All the npes*(npes-1) set of splitters are gathered at the root
  all_splitters = (long int *)malloc(npes*(npes-1)*size);
  MPI_Gather(local_splitters,npes-1,MPI_LONG,all_splitters,npes-1,MPI_LONG,root,comm);

  //4b. The root sorts all the npes*(npes-1) splitters and then chooses
  //npes-1 equally spaced elements from among them
  global_splitters = local_splitters;
  if (rank == root) {
    qsort(all_splitters,npes*(npes-1),size,compare);
    for (i=0;i<(npes-1);i++) {
      global_splitters[i]= all_splitters[(i+1)*(npes-1)];
    }    
  }

  //5a. Root sends splitters to each processor (including itself...)
  MPI_Bcast(global_splitters,npes-1,MPI_LONG,root,comm);

#if DEBUG_SPLITTERS
  MPI_Barrier(comm);
  if (rank == root)
    for (i=0;i<(npes-1)*npes;i++) {
      printf("all splitters: Rank %d %lu\n",rank,all_splitters[i]);
    }
  for (i=0;i<(npes-1);i++) {
    printf("global splitters: Rank %d %lu\n",rank,global_splitters[i]);
  }

#endif
  //5b. Each processor forms npes buckets from their sorted blocks and sends
  //appropriate elements to other processors (elements in bucket i go to
  //processor i
  //5b. Each processor buckets its array elements according with the
  //global_splitters
  idx = (int *)malloc((npes+1)*sizeof(int));//idx[j] holds the index+1 of the
  //last element of a (in this processor) that are in the bucket delimited by
  // global_splitters[j-1] and global_splitters[j] (>global_splitters[j-1]
  //but <=global_splitters[j]
  idx[0]=0;
  j=k=0;
  //k goes through all elements in the processor, j goes through the
  //npes-1 global splitters
  while ((k<count) && (j<(npes-1))) {
    if (a[k] > global_splitters[j]) {
      j++;
      idx[j]= k;
    }
    else {
      k++;
    }
  }
  //if j=npes-1, then k=j+1=npes and idx[npes]=count
  for (k=j+1;k<=npes;k++)
      idx[k] = count;
  
  //compute number of elements of a in this processor in each bucket
  send_counts = (int *) malloc(npes*sizeof(int));
  for (i=0;i<npes;i++)
    send_counts[i] = idx[i+1]-idx[i];

#if DEBUG_BUCKETS
  for (i=0;i<count;i++) {
    printf("Element in rank:%d %d %lu\n",rank,i,a[i]);
  }
  MPI_Barrier(comm);
  for (i=0;i<npes;i++){
    printf("Counts in rank:%d %lu %lu\n",rank,global_splitters[i],send_counts[i]);
  }
#endif
  recv_counts = (int *) malloc(npes*sizeof(int));
  //each processor now tells each other how many elements they will recieve
  //i.e. how many elements are in their bucket
  MPI_Alltoall(send_counts,1,MPI_LONG,recv_counts,1,MPI_LONG,comm);
  //imagine npes=3. Then
  //process 0: recv_counts[0]=send_counts[0] from p0,
  //recv_counts[1]=send_counts[0] from p1,recv_counts[2]=send_counts[0]
  //from p2, i.e. recv_counts for p0 has the number of elements in the
  //first bucket in each processor. Similarly for p1 and p2
  
#if DEBUG_BUCKETS
  MPI_Barrier(comm);
  for (i=0;i<npes;i++){
    printf("%d %lu %lu\n",rank,send_counts[i],recv_counts[i]);
}
#endif
  //The processors now have to send the proper elements to each other.
  //For example all elements in bucket 0 (as determined by recv_counts for
  //p0) should be sent to p0 from p1, p2,...
  //Compute offsets on the send and receive sides
  send_offsets = (int *) malloc(npes*sizeof(int));
  recv_offsets = (int *) malloc(npes*sizeof(int));
  k=j=0;
  for (i=0;i< npes;i++) {
    send_offsets[i] = k;//this is really the same as idx so a bit of wasted 
    //memory here
    k += send_counts[i];

    recv_offsets[i] = j;
    j += recv_counts[i];
    printf("recv_counts %d %d\n",rank,recv_counts[i]);
  }
  //buckets have different sizes given by j 
  printf("jj %d %d\n",rank,j);
  for (i =0;i<(npes-1);i++)
    printf("%d %d\n",rank,global_splitters[i]);
  buckets = (long int *)malloc (j*size);
  // MPI_Alltoallv(a,send_counts,send_offsets,MPI_LONG,buckets,recv_counts,recv_offsets,MPI_LONG,comm);
  //alltoallv: Imagine sends for p0, send_counts[0]=5,send_counts[1]=2,
  //send_counts[2]=7. Then a[0-4] from p0 goes to buckets[0-4] of p0
  //a[5-6] from p0 goes to buckets[0-1] of p1 and a[7-13] from p0 goes to
  //buckets[0-6] of p2. Recall that for p1 recv_counts[0]=send_counts[1] of
  //p0 and for p2 recv_counts[0]=send_counts[2] of p0, etc
  //each processor sorts local buckets. 
  qsort(buckets,j,size,compare);

#if DEBUG_FINAL_SORT
  MPI_Barrier(comm);
  for (i=0;i<count;i++) {
    printf("Rank:%d a[i]=%lu\n",rank,a[i]);
  }
  MPI_Barrier(comm);
  for (i=0;i<npes-1;i++) {
    printf("Rank:%d global_splitters[i]=%lu\n",rank,global_splitters[i]);
  }
  MPI_Barrier(comm);
  for (i=0;i<j;i++) {
    printf("Rank:%d buckets[i]=%lu\n",rank,buckets[i]);
  }
#endif
  /*  for (i=0;i<psize;i++) {
    printf("Rank after %d %lu\n",rank,a[i]);
    }*/
  
  free(local_splitters);
  free(all_splitters);
  free(global_splitters);
  free(idx);
  free(send_counts);
  free(recv_counts);
  free(send_offsets);
  free(recv_offsets);
  free(buckets);
}

void sort(long int n)
{
  long int i;
  long int* a;
  long int psize;//size of data in a given processor

  if (strcmp(ws,"strong") == 0) {
    psize = n/npes;
  }   
  if (strcmp(ws,"weak") == 0) {
    psize = n;
  }
 
  //Initialize array with random numbers.
  //Notice that rand generates a max random number of about 2 billion
  //so for n near or above a billion there will be many duplicates
  //To avoid this we would have to use a different generator that
  //would generate longs
  srand(n+rank);
  //assume that n is always divisible by npes
  a = (long int*) malloc (psize * sizeof(long int)); 
  for (i=0;i<psize;i++) {
      a[i] = rand();
  }
#if USE_SMALL_ARRAY
  for (i=0;i<n;i++) {
    a[i] = rand_interval(rank,rank*2+30);
  }   
#endif
  samplesort(a,psize,sizeof(long int),compare,MPI_COMM_WORLD);

  free(a);
}

int main(int argc,char** argv)
{ 
  int return_code;
  long int n;

  //Initialize
  MPI_Init(&argc,&argv); 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  if (argc != 3) { 
    if (rank == root)
      printf("Usage: mpirun/mpiexec -np xx ./samplesort n strong/weak\n");

    MPI_Finalize();
    return -1;
  }

  //Read size of array and whether we are doing strong/weak scaling
  n = atol(argv[1]);
  ws = argv[2];

  if (strcmp(ws,"strong")!=0 && strcmp(ws,"weak")!=0) {
    if (rank == root)
      printf("Invalid option. Please enter strong or weak for the last command line option.\n");

    MPI_Finalize();
    return -1;
  }   

#if USE_SMALL_ARRAY
  n = 32;
  ws = "weak";
#endif  
  sort(n);
  
  //Finalize
  MPI_Finalize();
  return 0;
}

