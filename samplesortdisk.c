#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <mpi.h>

#define GBYTE 1073741824
#define DEBUG_FILE_READ 1
#define DEBUG_INITIAL_SORT 0
#define DEBUG_SPLITTERS 0
#define DEBUG_BUCKETS 0
#define DEBUG_FINAL_SORT 0
#define DEBUG_FILE_WRITE 0

char * INPUT_PATH = "/uufs/chpc.utah.edu/common/home/ci-water4-0/CS6965/data/";
char * OUTPUT_PATH = "/uufs/chpc.utah.edu/common/home/ci-water4-0/CS6965/u0082100/assignment1/";
char * input_filename,*output_filename;//filename from where to read the array 
//elements and where to write final sorted buckets
int rank,npes, root =0;

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
  long int * buckets,*file_offsets;

  MPI_File output_file;
  MPI_Status status;
  MPI_Offset file_offset;
  int file_open_error,file_write_error;
 
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
    printf("rank=%d i=%d a[i]=%lu\n",rank,i,a[i]);
  }
  MPI_Barrier(comm);
  for (i=0;i<npes;i++){
    printf("rank=%d global_splitters[i]=%lu send_counts[i]=%lu\n",rank,global_splitters[i],send_counts[i]);
  }
#endif
  recv_counts = (int *) malloc(npes*sizeof(int));
  //each processor now tells each other how many elements they will recieve
  //i.e. how many elements are in their bucket
  MPI_Alltoall(send_counts,1,MPI_INT,recv_counts,1,MPI_INT,comm);
  //imagine npes=3. Then
  //process 0: recv_counts[0]=send_counts[0] from p0,
  //recv_counts[1]=send_counts[0] from p1,recv_counts[2]=send_counts[0]
  //from p2, i.e. recv_counts for p0 has the number of elements in the
  //first bucket in each processor. Similarly for p1 and p2
  
#if DEBUG_BUCKETS
  MPI_Barrier(comm);
  for (i=0;i<npes;i++){
    printf("rank=%d send_counts=%lu recv_counts=%lu\n",rank,send_counts[i],recv_counts[i]);
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
  }
  //buckets have different sizes given by j 
  buckets = (long int *)malloc (j*size);
  MPI_Alltoallv(a,send_counts,send_offsets,MPI_LONG,buckets,recv_counts,recv_offsets,MPI_LONG,comm);
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
    printf("Rank:%d i=%lu buckets[i]=%lu\n",rank,i,buckets[i]);
  }
  MPI_Barrier(comm);
#endif
  if (rank == root) {
    printf("Finished sorting\n");
  }
  
  //write buckets to file
  int *bucket_sizes = (int*) malloc(npes*sizeof(int));
  //all processes need to know the size of all buckets
  MPI_Allgather(&j, 1, MPI_INT, bucket_sizes, 1, MPI_INT, comm);

  file_offsets = (long int*) malloc(npes *size);
  file_offsets[0]=0;
  for (i=1;i<npes;i++) {
    file_offsets[i] = file_offsets[i-1] + bucket_sizes[i-1];
  }
  //open the file in all processes in create/write only mode 
  file_open_error = MPI_File_open(comm, output_filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
  if (file_open_error != MPI_SUCCESS) {
    char error_string[1000];
    int length_of_error_string, error_class;

    MPI_Error_class(file_open_error,&error_class);
    MPI_Error_string(error_class,error_string,&length_of_error_string);
    printf("%3d: %s\n",rank,error_string);

    MPI_Error_string(file_open_error,error_string,&length_of_error_string);
    printf("%3d: %s\n",rank,error_string);

    printf("Failed to open output file. Error Code=%d\n",file_open_error);
    MPI_Abort(comm, file_open_error);
  }
  file_offset = file_offsets[rank]*size;
  MPI_File_seek(output_file, file_offset, MPI_SEEK_SET);
  double start = MPI_Wtime();
  MPI_File_write(output_file, buckets, j, MPI_LONG, &status);
  double end = MPI_Wtime();
  double write_time = end-start;
  MPI_File_close(&output_file);
  double longest_write_time;
  MPI_Allreduce(&write_time,&longest_write_time,1,MPI_DOUBLE,MPI_MAX,comm);
  if (rank == root)
    printf("Finished writing sorted data to file\n");
  //compute write throughput
  long int total_number_of_bytes;
  for (i=0;i<npes;i++) 
    total_number_of_bytes += bucket_sizes[i];
  total_number_of_bytes = total_number_of_bytes * size;
  printf("Write throughput =%f GB/s\n",total_number_of_bytes/longest_write_time/GBYTE);
#if DEBUG_FILE_WRITE
  printf("My rank is %d. My bucket size was %d and my first element was %ld and the last was %ld\n",rank,j,buckets[0],buckets[j-1]);
  printf("rank %d has bucket size %d\n",rank,bucket_sizes[rank]);
  MPI_Barrier(comm);
   
  //read elements from file anin chunks the size of rank 0 bucket size
  //and check if all elements are in increasing order
  if (rank == root) {
    printf("Reading output file\n");
    FILE *ftest = fopen(output_filename, "rb");
    int upper = 2* bucket_sizes[0];
    long int curr_value, prev_value = 0, file_pos = 0;
    int low_limit=0;
    long int *test_buffer = (long int*)malloc(upper*size);
    int stop = 0;
    do {
      curr_value = fread(test_buffer, size, upper, ftest);
      for(i = 0; i < upper; i++) {
	file_pos++;
	if(test_buffer[i] >= prev_value) {
	  prev_value = test_buffer[i];
	  //print first and last elements of first,second,third,... buckets
	  int sum =0;
	  for (k=0;k<npes;k++) {
	    sum += bucket_sizes[k];
	    if (file_pos== sum) {
	      if (k==0)
		printf("rank %d range %ld %ld\n",k,test_buffer[0],test_buffer[i]);
	      else
		printf("rank %d range %ld %ld. File position %ld\n",k,low_limit,test_buffer[i],file_pos);
	      if (k==(npes-1)) {
		stop = 1;
		break;
	      }
	    }
	    else if (file_pos == sum+1) {
	      low_limit = test_buffer[i];
	    }
	  }
	  continue;
	}
	else if (stop == 0) {
	  printf("value %ld at i=%d is smaller than value %ld at some previous position\n", test_buffer[i],i,prev_value);
	  printf("This happened at location %ld in the file %d\n",file_pos); 
	  break;
	}
      }
    } while (curr_value>0 && stop ==0);
    free(test_buffer);
  }
#endif

  free(local_splitters);
  free(all_splitters);
  free(global_splitters);
  free(idx);
  free(send_counts);
  free(recv_counts);
  free(send_offsets);
  free(recv_offsets);
  free(buckets);
  free(file_offsets);
  free(bucket_sizes);
}

void sort()
{
  long int * a;
  long int psize;//size of array data in a given processor
  MPI_File input_file;
  MPI_Status status;
  //MPI_Offset is long long
  MPI_Offset num_bytes, num_bytes_ll, max_num_bytes_ll, file_offset;
  int file_open_error, file_read_error;
 
  //open the file in all processes in read only mode 
  file_open_error = MPI_File_open(MPI_COMM_WORLD,input_filename,MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file);
  if (file_open_error != MPI_SUCCESS) {
    char error_string[1000];
    int length_of_error_string, error_class;

    MPI_Error_class(file_open_error,&error_class);
    MPI_Error_string(error_class,error_string,&length_of_error_string);
    printf("%3d: %s\n",rank,error_string);

    MPI_Error_string(file_open_error,error_string,&length_of_error_string);
    printf("%3d: %s\n",rank,error_string);

    printf("Failed to open input file. Either the name or the path of the file is wrong. Error Code=%d\n",file_open_error);
    MPI_Abort(MPI_COMM_WORLD, file_open_error);
  } 
  MPI_File_get_size(input_file,&num_bytes);
#if DEBUG_FILE_READ
  printf("%3d: total_number_of_bytes = %lld\n",rank,num_bytes);
#endif
  num_bytes_ll = num_bytes/npes;
  //If npes does not divide num_bytes evenly, the last process will have to 
  //read more data, i.e. to the end of the file
  max_num_bytes_ll = num_bytes_ll + num_bytes % npes;

  if (num_bytes_ll != max_num_bytes_ll) {
    printf("Last process will read more bytes than previous processes\n");
  }
  if (max_num_bytes_ll < INT_MAX) {
    //figure out the number of int elements in the file
    psize = num_bytes_ll/sizeof(long int);
#if DEBUG_FILE_READ      
    printf("Array Size=%d\n",psize);
#endif
    //read psize elements into array a for each process
    file_offset = rank * psize;
    //assume that n is always divisible by npes
    a = (long int*) malloc (psize * sizeof(long int)); 
    MPI_File_seek(input_file,file_offset,MPI_SEEK_SET);
    file_read_error = MPI_File_read(input_file,a,psize,MPI_LONG,&status);
    //inspect content of binary file with 
    //hexdump -v -e '7/4 "%12d "' -e '"\n"' ../../../ci-water4-0/CS6965/data/sort_data_debug.dat > tt
    //and check that the first few digits in tt match a[0],a[1].a[2] returned 
    //by p0. Files have 8 byte (long) numbers

#if DEBUG_FILE_READ
    printf("rank:%d a[0]:%ld psize-1:%ld a[psize-1]:%ld\n", rank,a[0],(psize-1),a[psize-1]);
#endif
    if (rank == root) {
      printf("Finished reading data from file\n");
    }
    if (file_read_error != MPI_SUCCESS) {
      printf("Read fail by rank=%d\n",rank);
      MPI_Abort(MPI_COMM_WORLD,file_read_error);
    }
  }
  else {
    printf("Not enough memory to read the file.\n");
    printf("Consider running on more nodes.\n");
  }
  MPI_File_close(&input_file);
 
  double t1,t2;
  t1 = MPI_Wtime();
  samplesort(a,psize,sizeof(long int),compare,MPI_COMM_WORLD);
  t2 = MPI_Wtime();
  // printf("Time elapsed in rank %d was %f\n",rank, t2-t1);
  free(a);
}

int main(int argc,char** argv)
{ 
  int return_code;

  //Initialize
  MPI_Init(&argc,&argv); 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  if (argc != 3) { 
    if (rank == root)
      printf("Usage: mpirun/mpiexec -np xx ./samplesortdisk input_filename output_filename\n");

    MPI_Finalize();
    return -1;
  }

  //Read whether we are doing strong/weak scaling and the filename
  input_filename = malloc(900);
  strcpy(input_filename,INPUT_PATH);
  strcat(input_filename,argv[1]);

  output_filename = malloc(900);
  strcpy(output_filename,OUTPUT_PATH);
  strcat(output_filename,argv[2]);
  sort();
  printf("%s %s\n",input_filename,output_filename);
  free(input_filename);
  free(output_filename);
  //Finalize
  MPI_Finalize();
  return 0;
}
