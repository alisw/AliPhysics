// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

/** \class AliHLTTPCBenchmark
<pre>
//_____________________________________________________________
//
// AliHLTTPCBenchmark
//
//   Benchmark class for level3 code
//  
//
</pre>
*/

#ifndef no_root
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TMath.h>
#endif
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCBenchmark.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCBenchmark)

AliHLTTPCBenchmark::AliHLTTPCBenchmark()
  :
  fNbench(0),
  fNmax(20),
  fNames(NULL),
  fTimer(NULL),
  fSum(NULL),
  fMin(NULL),
  fMax(NULL),
  fCount(NULL)
{
  //Constructor
}

AliHLTTPCBenchmark::AliHLTTPCBenchmark(const AliHLTTPCBenchmark&)
  :
  fNbench(0),
  fNmax(20),
  fNames(NULL),
  fTimer(NULL),
  fSum(NULL),
  fMin(NULL),
  fMax(NULL),
  fCount(NULL)
{
  HLTFatal("copy constructor untested");
}

AliHLTTPCBenchmark& AliHLTTPCBenchmark::operator=(const AliHLTTPCBenchmark&)
{ 
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTTPCBenchmark::~AliHLTTPCBenchmark()
{
  //deconstructor
   fNbench   = 0;
   if (fNames)  {delete [] fNames; fNames  = 0;}
   if (fTimer)  {delete [] fTimer; fTimer  = 0;}
   if (fSum)    {delete [] fSum;   fSum   = 0;}
   if (fMin)    {delete [] fMin;   fMin   = 0;}
   if (fMax)    {delete [] fMax;   fMax   = 0;}
   if (fCount)  {delete [] fCount; fCount =0;}
}

Int_t AliHLTTPCBenchmark::GetBench(const Char_t *name)
{
  //get bench with name
   for (Int_t i=0;i<fNbench;i++) {
      if (!strcmp(name,(const Char_t*)fNames[i])) return i;
   }
   return -1;
}


void AliHLTTPCBenchmark::Start(const Char_t *name)
{
  //start the benchmark with name
   if (!fNbench) {
#ifdef no_root
     fNames=new Char_t*[fNmax];
     fTimer = new AliHLTTPCStopwatch[fNmax];
#else
     fNames = new TString[fNmax];
     fTimer = new TStopwatch[fNmax];
#endif

     fSum   = new Float_t[fNmax];
     fMin   = new Float_t[fNmax];
     fMax   = new Float_t[fNmax];
     fCount = new Int_t[fNmax];
     for(Int_t i =0;i<fNmax;i++){
       fSum[i]=0;
       fMin[i]=0;
       fMax[i]=0;
       fCount[i]=0;
     }
   }
   Int_t bench = GetBench(name);
   if (bench < 0 && fNbench < fNmax ) {
      // define a new benchmark to Start
#ifdef no_root
     fNames[fNbench]=new Char_t[strlen(name)+1];
     strcpy(fNames[fNbench],name);
#else
      fNames[fNbench] = name;
#endif
      bench = fNbench;
      fNbench++;
      fTimer[bench].Reset();
      fTimer[bench].Start();
   } else if (bench >=0) {
   // Resume the existent benchmark
      fTimer[bench].Reset();
      fTimer[bench].Start();
   }
   else
     LOG(AliHLTTPCLog::kWarning,"AliHLTTPCBenchmark::Start","Start")
     <<"too many benches"<<ENDLOG;
}

void AliHLTTPCBenchmark::Stop(const char *name)
{
  //stop the benchmark with name
   Int_t bench = GetBench(name);
   if (bench < 0) return;

   fTimer[bench].Stop();
   Float_t val = fTimer[bench].CpuTime();
   
   fSum[bench] += val; 
   fCount[bench]++;
   if(fCount[bench]==1){
     fMin[bench] = val;
     fMax[bench] = val;
   }
   else{
     if(val<fMin[bench])fMin[bench]=val;
     if(val>fMax[bench])fMax[bench]=val;
   }
}

void AliHLTTPCBenchmark::Analyze(const Char_t* name)
{
  //get results of benchmark
  Float_t *x = new Float_t[fNbench]; 
  Float_t *y = new Float_t[fNbench];
  Float_t *eyl = new Float_t[fNbench]; 
  Float_t *eyh = new Float_t[fNbench];
  Char_t filename[256];
  sprintf(filename,"%s.dat",name);
  FILE *f= fopen(filename,"w");
  for (Int_t i=0;i<fNbench;i++) {
    Float_t av =0;
    if(fCount[i]) av = fSum[i]/fCount[i]; 
    x[i]=i+1;
    y[i]=av*1000;
    eyl[i]=(av-fMin[i])*1000;
    eyh[i]=(fMax[i]-av)*1000;
#ifdef no_root
    fprintf(f,"%2d. %s: ",i+1,fNames[i]);
#else
    fprintf(f,"%2d. %s: ",i+1,fNames[i].Data());
#endif
    fprintf(f,"total %4.0f patch %4.0f -%4.0f +%4.0f ms\n",fSum[i],av*1000,eyl[i],eyh[i]);
  }
  fclose(f);
  sprintf(filename,"%s.tmp",name);
/* only a workaround!!
  FILE *f2= fopen(filename,"w");
  for (Int_t i=0;i<fNbench;i++) fprintf(f2,"%f ",x[i]); fprintf(f2,"\n");
  for (Int_t i=0;i<fNbench;i++) fprintf(f2,"%f ",y[i]); fprintf(f2,"\n");
  for (Int_t i=0;i<fNbench;i++) fprintf(f2,"%f ",eyl[i]); fprintf(f2,"\n");
  for (Int_t i=0;i<fNbench;i++) fprintf(f2,"%f ",eyh[i]); fprintf(f2,"\n");
  fclose(f2);
*/
#ifndef no_root
  sprintf(filename,"%s.root",name);
  TFile *file = new TFile(filename,"RECREATE");
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(fNbench,x,y,0,0,eyl,eyh);
  gr->SetTitle("benchmark");
  gr->SetMarkerStyle(8);
  gr->SetMinimum(0);
  //gr->Draw("ALP");
  gr->Write();
  file->Close();
  delete file; 
  file=0;
#endif
  delete[] x;
  delete[] y;
  delete[] eyl;
  delete[] eyh;
}

Double_t AliHLTTPCBenchmark::GetCpuTime()
{
  //get cpu time
  {return (Double_t)(clock()) / CLOCKS_PER_SEC;}
}
