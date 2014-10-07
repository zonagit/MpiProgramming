#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

float values[] = {4,1,10,9,2,5,-1,-9,-2,10000,-0.05,-3,-1.1};

int compare (const void *a, const void *b)
{
  float fa = *(const float*) a;
  float fb = *(const float*)b;

  return (fa>fb)-(fa<fb);
}

int main ()
{
  int i;

  qsort(values,13,sizeof(float),compare);

  for (i=0;i<13;i++)
    {
      printf("%f,",values[i]);
    }

  putchar('\n');

  printf("Max int %d\n",INT_MAX);
  printf("Max long int %lu\n",LONG_MAX);
  printf("Max random %lu\n",RAND_MAX);
  return 0;
}
    
