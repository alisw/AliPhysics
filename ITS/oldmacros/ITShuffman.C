#include "iostream.h"

void ITShuffman_f (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("stream.root");
   if (!file) file = new TFile("stream.root");
   file->ls();


   AliITSInStream  *IStream;

   char *path=gDirectory->GetPath();
   printf("path %s\n",path);
   

   Int_t nkeys=gDirectory->GetNkeys();
   printf("nkeys %d\n",nkeys);
   char namecycle[40];

   AliITSHTable *hufft=new AliITSHTable(256);
   for (int i=1;i<=nkeys;i++) {
     Int_t count=0;
     Int_t count1=0;
     Int_t count2=0;
     Int_t count3=0;
       sprintf(namecycle,"AliITSInStream;%d",i);
       printf("namecycle %s\n",namecycle);
       IStream=(AliITSInStream*)gDirectory->Get(namecycle);
       UChar_t *str=IStream->Stream();
       Int_t len=IStream->StreamLength();
       //printf("str len %p %d\n",str,len);
       for (int k=0;k<len;k++) {
           UChar_t elem=str[k];
	   // test
           if (elem==128) count++;
           if (elem==127) count1++;
           if (elem==129) count2++;
           if (elem==134) count3++;
	   //if (elem != 128) printf("i,k,elem %d %d %d\n",i,k,elem);
       }
       //printf("count 128 127 129 134 %d %d %d %d \n",count,count1,count2,count3);
       hufft->GetFrequencies(len,str);
   }

   hufft->BuildHTable();

   cout<<"END  test for InStream "<<endl;

   file->Close();   
}



