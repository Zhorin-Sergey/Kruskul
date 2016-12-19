#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#define PI 3.14159265

int main(int argc, char* argv[]) {
  MPI_Status status;
  int i = 0;
  double h, l, q = 1;
  int thread_count, rank;
  double *duff;
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  h = (1.0) / n;
  l = (2.0) / m;
  double *a = 0;
  double *a1 = 0;

  MPI_Init(&argc, &argv);

  i = 0;
  for (int uy = 0; uy < 5; uy++) {
    m = m + 100;
    n = n + 100;
    int segmentSize, bufferSize;
    double starttime, endtime, starttime1, endtime1;

    MPI_Comm_size(MPI_COMM_WORLD, &thread_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    segmentSize = n / thread_count;
    bufferSize = segmentSize;


    if (rank == 0) {
      int k = 0;
      a = new double[n*m];
      a1 = new double[n*m];
      while (k < n*m) {
        a[k] = 0;
        k++;
      }
      k = 0;
      for (int i = 1; i < n - 1; i++)
      {

        a[i] = 0.1 * sin(PI * i * h);

      }
      for (int i = 1; i < n - 1; i++)
      {
        a[n + i] = a[i] + 0 * l + (l * l / 2) * ((a[i + 1] - 2 * a[i] + a[i - 1]) / (h * h));
      }
      for (int i = 0; i < n*m; i++)
      {
        a1[i] = a[i];
      }
      i = 2 * (n)+1;
      starttime1 = MPI_Wtime();
      while (i < n*m) {
        for (int j = 0; j < n - 2; j++)
          a1[i + j] = 2 * a1[i + j - n] - a1[i + j - 2 * n] + (q * q * l * l / (h * h)) * (a1[i + j - n + 1] - 2 * a1[i + j - n] + a1[i + j - n - 1]);
        i = i + n;
      }
      endtime1 = MPI_Wtime();
      printf("Sequential algorithm time %1.10f\n", endtime1 - starttime1);
      i = 0;
      k = 0;
      int p = 1;
      while (k < m*n) {
        if (n*m <= 25) {
          printf("%f ", a1[k]);
          if (k == (n*p) - 1) {
            printf(" \n");
            p++;
          }
        }
        k++;
      }
      // printf(" \n");
      i = 0;
      k = 0;
      p = 1;
      while (k < m*n) {
        if (n*m <= 25) {
          printf("%f ", a[k]);
          if (k == (n*p) - 1) {
            printf(" \n");
            p++;
          }
        }
        k++;
      }
      //printf(" \n");




    }

    if (rank == 0) {

      i = 0;
      duff = new double[(segmentSize + 2)*m];
      while (i < (segmentSize + 2)*m) {

        duff[i] = 0;
        i++;
      }
      i = 1;

      while (i < thread_count) {
        MPI_Send(a + (i)*segmentSize - 1, segmentSize + 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        MPI_Send(a + n + (i)*segmentSize - 1, segmentSize + 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        i++;
      }
      i = 0;
    }
    else {
      duff = new double[(segmentSize + 2)*m];
      i = 0;
      while (i < (segmentSize + 2)*m) {
        duff[i] = 0;
        i++;
      }
      i = 0;
      MPI_Recv(duff, segmentSize + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(duff + segmentSize + 2, segmentSize + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }
    i = 0;
    while (i<(segmentSize + 2)*m) {
      i++;
    }
    //printf("rank %i \n", rank);

    if (rank > 0) {
      int dest = rank + 1;
      i = 0;
      double d = 0;
      //printf(" \n");
      if (rank == thread_count - 1)
        dest = 0;
      i = 2 * (segmentSize + 2) + 1;
      while (i < (segmentSize + 2)*m) {
        for (int j = 0; j < segmentSize; j++) {

          duff[i + j] = 2 * duff[i + j - (segmentSize + 2)] - duff[i + j - 2 * (segmentSize + 2)] + (q * q * l * l / (h * h)) * (duff[i + j - (segmentSize + 2) + 1] - 2 * duff[i + j - (segmentSize + 2)] + duff[i + j - (segmentSize + 2) - 1]);
          
        }

        MPI_Send(duff + i, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Send(duff + i + segmentSize - 1, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        MPI_Recv(duff + i - 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(duff + i - 2 + (segmentSize + 2), 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &status);

        i = i + (segmentSize + 2);
      }
      i = 0;
      while (i < (segmentSize + 2)*m) {
        //  printf("rank %i %f \n", rank, duff[i]);
        i++;
      }
      MPI_Send(duff, (segmentSize + 2)*m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      // MPI_Send(&d, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
      //double c = 1, d = 2;
      i = 2 * (n)+1;
      int o = 2 * (n)+thread_count*segmentSize;
      starttime = MPI_Wtime();
      while (i < n*m) {
        for (int j = 0; j < segmentSize - 1; j++) {
          //printf("i+j %f \n",a[i+j-n]);
          a[i + j] = 2 * a[i + j - n] - a[i + j - 2 * n] + (q * q * l * l / (h * h)) * (a[i + j - n + 1] - 2 * a[i + j - n] + a[i + j - n - 1]);
          //printf("rank %i \n", rank);
        }
        // printf("rank1 %i %f \n", rank, a[10 + 2 * n]);
        for (int j = 0; j < n - thread_count*segmentSize - 1; j++){
          a[o + j] = 2 * a[o + j - n] - a[o + j - 2 * n] + (q * q * l * l / (h * h)) * (a[o + j - n + 1] - 2 * a[o + j - n] + a[o + j - n - 1]);
        }
        if (thread_count != 1){

          MPI_Send(a + i + segmentSize - 2, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
          MPI_Send(a + o, 1, MPI_DOUBLE, thread_count - 1, 0, MPI_COMM_WORLD);
          MPI_Recv(a + i - 3 + (segmentSize + 2), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(a + o - 1, 1, MPI_DOUBLE, thread_count - 1, 0, MPI_COMM_WORLD, &status);

        }
        i = i + n;
        o = o + n;
        //  printf("rank1 %i %f \n", rank, a[10 + 2 * n]);
      }
      endtime = MPI_Wtime();
      //printf("rank12 %i %f \n", rank, a[10 + 2 * n]);
      int y;
      for (int u = 1; u < thread_count; u++) {
        y = 1 + 2 * (segmentSize + 2);
        MPI_Recv(duff, (segmentSize + 2)*m, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &status);
        // MPI_Recv(&c, 1, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &status);
        i = 0;
        while (i < (segmentSize + 2)*m) {
          // printf("rank %i %f \n", rank, duff[i]);
          i++;
        }
        //printf(" \n");
        for (int i = segmentSize*u + 2 * n; i < n*m; i = i + n) {
          for (int j = 0; j < segmentSize; j++)
            a[i + j] = duff[y + j];
          y = y + segmentSize + 2;
        }

      }

      printf("Parallel algorithm time on %i X %i %1.10f\n",n , m, endtime - starttime);
      int k = 0;
      int p = 1;
      if (m*n <= 25)
      while (k < m*n) {

        printf("%f ", a[k]);
        if (k == (n*p) - 1) {
          printf(" \n");
          p++;

        }
        k++;
      }

      int flag = 1;
      for (int i = 0; i < m*n; i++){
        if (a[i] != a1[i])
          flag = 0;
      }
      printf("%i \n", flag);

    }
  }
    MPI_Finalize();
  
    return 0;
  
}
