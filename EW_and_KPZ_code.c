#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran2.c"

#define L 1000   // system size
#define TMAX 10000 // number of time steps

long sd=-937176381;

int main()
{
    int i, t,count=0,sum=0,time=0,*avl,*small_avl,total_count=0,count_1=0,b,bundle=1;
    double h[L] = {0},gi[L],v[TMAX] = {0};
    double pinning_force[L],F=0,minimum_load=0,diff=0,force[L]={0},closest_zero;

    FILE*fp;
    fp=fopen("QKPZ.txt","w");

    avl=(int *)malloc(L*sizeof(int));

    for(b=0;b<bundle;b++)
    {
       time=0; total_count=0; count_1=0; count=0; F=0; diff=0;
       for(i=0;i<L;i++){avl[i]=0; h[i] = 0; force[i]=0;}

    for(i=0;i<L;i++)
    {
       h[i]=0;
       pinning_force[i] = 4*ran2(&sd)-2;                       ///uniform distribution between -2 to 2 (edward wilkinson)
    /// pinning_force[i] = 10*ran2(&sd)-5;                        /// uniform distribution between -5 to 5 (kardar parisi)
    }

  F=0.0001;
  while(F < 10)
      {
       total_count=0;
       count_1=1;
       time++;
       count=1;
       while(count>0 && count_1<1000*L)
        {
        count=0;
        count_1++;

        for (i = 0; i < L; i++)
         {
            gi[i] = h[(i+1)%L] + h[(i-1+L)%L] - 2*h[i] + pinning_force[i] + F; ///growth rule for EW equation
            ///gi[i] = h[(i+1)%L] + h[(i-1+L)%L] - 2*h[i] + 0.5*((h[(i-1+L)%L]-h[(i+1)%L])/2)*((h[(i-1+L)%L]-h[(i+1)%L])/2) + pinning_force[i] + F;  /// growth rule for kpz equation
            if (gi[i] >= 0)
             {
                h[i] += 1;
                count++;
                pinning_force[i] = 4*ran2(&sd)-2;  ///edward wilkinson
              ///  pinning_force[i] = 10*ran2(&sd)-5;   ///kardar parisi
             }
         }
          total_count+=count;

        }
///---------- end of second while loop (redistribution loop)-------------///////////////////
        if (count_1>10*L){F=L;} //////escape when depinned
        //printf("\n%d %d\n",count_1, count);

         closest_zero = INFINITY;
         for (int i = 0; i < L; i++)
         {
            diff = gi[i];
            if (diff < 0 && fabs(diff) < fabs(closest_zero))
            {
              closest_zero = diff;
            }
         }

        F += fabs(closest_zero);

        avl[time]=total_count;
        force[time]=F;

    }
 ///----------------end of 1st while (force) loop----------------///
        for(i=2;i<time;i++)
        {
            fprintf(fp,"%d %lf %d %d\n",i,force[i],avl[i],time-i-1);
        }
    }
    free(avl);
}

