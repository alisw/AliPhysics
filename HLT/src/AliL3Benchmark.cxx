#include <TFile.h>
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TStopwatch.h"
#include "AliL3Benchmark.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliL3Logging.h"

ClassImp(AliL3Benchmark)

AliL3Benchmark::AliL3Benchmark()
{
   fNbench = 0;
   fNmax   = 20;
   fNames  = 0;
   fTimer  = 0;
   fSum    = 0;
   fMin    = 0;
   fMax    = 0;
   fCount  = 0;
//   fStopwatch = 0;
}

AliL3Benchmark::~AliL3Benchmark()
{
   fNbench   = 0;
   if (fNames)  {delete [] fNames; fNames  = 0;}
   if (fTimer)  {delete [] fTimer; fTimer  = 0;}
   if (fSum)    {delete [] fSum;   fSum   = 0;}
   if (fMin)    {delete [] fMin;   fMin   = 0;}
   if (fMax)    {delete [] fMax;   fMax   = 0;}
   if (fCount)  {delete [] fCount; fCount =0;}
//   if(fStopwatch) {delete fStopwatch; fStopwatch =0;}
}

Int_t AliL3Benchmark::GetBench(const char *name)
{
   for (Int_t i=0;i<fNbench;i++) {
      if (!strcmp(name,(const char*)fNames[i])) return i;
   }
   return -1;
}


void AliL3Benchmark::Start(const char *name)
{
   if (!fNbench) {
      fNames = new TString[fNmax];
      fTimer = new TStopwatch[fNmax];
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
      fNames[fNbench] = name;
      bench = fNbench;
      fNbench++;
      fTimer[bench].Reset();
      fTimer[bench].Start();
//      if(fStopwatch) {delete fStopwatch; fStopwatch =0;}
//      fStopwatch = new TStopwatch();
//      fStopwatch->Reset();
//      fStopwatch->Start();
   } else if (bench >=0) {
   // Resume the existen benchmark
      fTimer[bench].Reset();
      fTimer[bench].Start();
//      if(fStopwatch) {delete fStopwatch; fStopwatch =0;}
//      fStopwatch = new TStopwatch();
//      fStopwatch->Reset();
//      fStopwatch->Start();
   }
   else
     LOG(AliL3Log::kWarning,"AliL3Benchmark::Start","Start")
     <<"too many benches"<<ENDLOG;
}

void AliL3Benchmark::Stop(const char *name)
{
   Int_t bench = GetBench(name);
   if (bench < 0) return;

   fTimer[bench].Stop();
   Float_t val = fTimer[bench].CpuTime();
//   fStopwatch->Stop();
//   Float_t val = fStopwatch->CpuTime();
   
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

void AliL3Benchmark::Analyze(const char* name){
  Float_t *x = new Float_t[fNbench]; 
  Float_t *y = new Float_t[fNbench];
  Float_t *eyl = new Float_t[fNbench]; 
  Float_t *eyh = new Float_t[fNbench];
  char filename[256];
  sprintf(filename,"%s.dat",name);
  FILE *f= fopen(filename,"w");
  for (Int_t i=0;i<fNbench;i++) {
    Float_t av =0;
    if(fCount[i]) av = fSum[i]/fCount[i]; 
    x[i]=i+1;
    y[i]=av*1000;
    eyl[i]=(av-fMin[i])*1000;
    eyh[i]=(fMax[i]-av)*1000;
    fprintf(f,"%2d. %s: ",i+1,fNames[i].Data());
    fprintf(f,"%4.0f ms\n",av*1000);
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
  sprintf(filename,"%s.root",name);
  TFile *file = new TFile(filename,"RECREATE");
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(fNbench,x,y,0,0,eyl,eyh);
  gr->SetTitle("benchmark");
  gr->SetMarkerStyle(8);
  gr->SetMinimum(0);
  gr->Draw("ALP");
  gr->Write();
  file->Close();
  delete file; 
  file=0;
  delete[] x;
  delete[] y;
  delete[] eyl;
  delete[] eyh;
}


