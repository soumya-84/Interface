#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ran2.c"
#define L 5001

long sd=-937136753;

int main()
{
  int i,k,n,NRUN=1,bundle=1,t,time,TMAX=10000,count2,breaking,list[L],dis,*stat,count,scan,breaking_small,avc_n=0;
  float *thr,*force,load,velo[TMAX],ava_small[10000],force_store[L],force_before,force_after,sigma_t[TMAX],min_diff,tot_force,diff,gamma,load_share[L];
  float sum,energy,*avalanche_series,*load_series;


  stat=(int *)malloc(L*sizeof(int));
  thr=(float *)malloc(L*sizeof(float));
  force=(float *)malloc(L*sizeof(float));
  avalanche_series=(float *)malloc(TMAX*sizeof(float));
  load_series=(float *)malloc(TMAX*sizeof(float));

  gamma=5;

  for(i=0;i<L;i++) {load_share[i]=0;}
  sum=0;
  for(i=0;i<L/2;i++)
  {
     sum+=1./pow(i+1,gamma);
  }

  for(i=0;i<L/2;i++)
  {
     load_share[i+1]=1./(2*pow(i+1,gamma)*sum);
  }

  for(i=0;i<10000;i++){ava_small[i]=0;}
  for(i=0;i<TMAX;i++) {velo[i]=0; sigma_t[i]=0;}

   load=0.; time=0;
   for(n=0;n<NRUN;n++)
   {
        for(i=0;i<L;i++)
        {
         stat[i]=1;
         thr[i]=ran2(&sd);
         force[i]=0;
        }

        count=0;
        for(i=0;i<L;i++)
        {
          count+=stat[i];
        }

        time=0;
    	for(t=0;t<TMAX;t++)
    	{
            time++;
            min_diff=L;

            for(i=0;i<L;i++)
            {
               diff=thr[i]-force[i];
               if(diff<min_diff) {min_diff=diff;}
            }

            for(i=0;i<L;i++) {force[i]+=min_diff;}

          count2=1;
          breaking=0; scan=0; energy=0;
          while(count2>0)
          {
            force_before=0;
            for(i=0;i<L;i++) {list[i]=0; force_store[i]=0; force_before+=force[i];}
            count=0; scan++; breaking_small=0;
            for(i=0;i<L;i++)
            {
              if(force[i]>=thr[i]) {if(stat[i]==1) {breaking++; energy+=thr[i]*force[i];} stat[i]=0;  list[breaking_small]=i;  breaking_small++;}
              count+=stat[i];
            }


            for(i=0;i<L;i++)
            {
               for(k=0;k<breaking_small;k++)
               {
                  dis=fabs(list[k]-i);
                  if(dis>L/2) {dis=L-dis;}
                  if(dis!=0) {force_store[i]+=force[list[k]]*load_share[dis];}
               }

            }
            force_after=0;
            for(i=0;i<L;i++)
            {
              if(stat[i]==0) {force_after+=force_store[i]; force[i]=force_store[i]; stat[i]=1; thr[i]=ran2(&sd);}
              else {force[i]+=force_store[i]; force_after+=force[i];}
            }
          count2=0;
          for(i=0;i<L;i++)
          {
              if(force[i]>thr[i]) {count2++;}
          }
          ava_small[breaking_small]++; avc_n++;

          if(scan>500)
          {
             count2=0; t=2*TMAX;
          }

         }
//end of while loop (redistribution loop)

          load_series[time]=tot_force/L; avalanche_series[time]= breaking;

          tot_force=0;
          for(i=0;i<L;i++)
          {
            tot_force+=force[i];
          }

          if(t<TMAX)
          {
            sigma_t[t]+=tot_force/L;
            velo[t]+=(float)breaking/L;
          }
        }
//end of time

  for(i=2;i<=time;i++){printf("%d %lf %lf %d\n",i,load_series[i],avalanche_series[i],time-i);}
  }
 //end of realisation NRUN loop

 free(stat); free(avalanche_series); free(force); free(load_series); free(thr);
}
