#include "iostream.h"

void ITSstream (Int_t evNumber1=0,Int_t evNumber2=0) 
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

// Get AliRun object from file or create it if not on file

   /*
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   */

   AliITSInStream  *IStream;

   char *path=gDirectory->GetPath();
   printf("path %s\n",path);
   

   Int_t nkeys=gDirectory->GetNkeys();
   printf("nkeys %d\n",nkeys);
   char namecycle[40];

   for (int i=1;i<=nkeys;i++) {
     //sprintf(namecycle,"%s%d",GetName(),i);
       sprintf(namecycle,"AliITSInStream;%d",i);
       printf("namecycle %s\n",namecycle);
       IStream=(AliITSInStream*)gDirectory->Get(namecycle);
       //Char_t *str=IStream->Stream();
       UChar_t *str=IStream->Stream();
       Int_t len=IStream->StreamLength();
       printf("str len %p %d\n",str,len);
       for (int k=0;k<len;k++) {
	 //Char_t elem=str[k];
           UChar_t elem=str[k];
	   printf("i,k,elem %d %d %d\n",i,k,elem);
	   //if (elem > 0) printf("i,k,elem %d %d %d\n",i,k,elem);
       }
   }

     cout<<"END  test for InStream "<<endl;

     file->Close();   
}



