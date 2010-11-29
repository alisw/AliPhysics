#include  "AliCaloRawAnalyzerPeakFinder.h"
//#include  "AliCaloRawAnalyzerPeakFinder.h"

#include "AliCaloPeakFinderVectors.h"
#define MAXSTART 3
#define SAMPLERANGE 15
#define SHIF 0.5


void pfocdbmaker()
{
  AliCDBEntry e;
  AliCaloRawAnalyzerPeakFinder *p = new  AliCaloRawAnalyzerPeakFinder();
  e.SetObject(p);
  TFile *f = new TFile("pfvectors.root", "recreate" );
  //Read in the Peak finder vecors from file
  p->Write();
  //  e.Write();
  f->Close();
}






