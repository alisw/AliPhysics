///
/// \file pfocdbmaker.C
/// \ingroup EMCAL_SimRecDB
/// \brief Create peak finder vectors
///
/// Create peak finder vectors
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
///

#if !defined(__CINT__)

#include "TFile.h"

#include "AliCaloRawAnalyzerPeakFinder.h"
#include "AliCaloPeakFinderVectors.h"

#include "AliCDBEntry.h"

#endif

#define MAXSTART 3
#define SAMPLERANGE 15
#define SHIF 0.5

///
/// Main method
///
void pfocdbmaker()
{
  AliCDBEntry e;
  
  AliCaloRawAnalyzerPeakFinder *p = new  AliCaloRawAnalyzerPeakFinder();
  
  e.SetObject(p);
  
  TFile *f = new TFile("pfvectors.root", "recreate" );
  
  // Read in the Peak finder vecors from file
  p->Write();
  
  //  e.Write();
  
  f->Close();
}






