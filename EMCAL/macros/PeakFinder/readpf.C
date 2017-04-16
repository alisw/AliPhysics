///
/// \file readpf.C
/// \ingroup EMCAL_SimRecDB
/// \brief Read peak finder vectors
///
/// Read peak finder vectors
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
///

#if !defined(__CINT__)

#include "TFile.h"

#include "Riostream.h"

#include "AliCaloRawAnalyzerPeakFinder.h"
#include "AliCaloPeakFinderVectors.h"

#endif

void readpf()
{
  TFile *f = new TFile("peakfindervectors2.root", "read");
  
  //   f->ls();
  //   f->Print();
  //   AliCaloRawAnalyzerPeakFinder *p = (AliCaloRawAnalyzerPeakFinder*)f->GetKey("AliCaloRawAnalyzerPeakFinder");
  
  AliCaloPeakFinderVectors *p = (AliCaloPeakFinderVectors *) f->Get("AliCaloPeakFinderVectors");
  
  if ( p == 0 )
  {
    cout << "ERROR, P == 0" << endl;
  }
  else
  {
    cout << "INFO, p = "<< p  << endl;
    p->PrintVectors();
  }
  
  cout << endl <<  "********** !!!!!!!!!!!!!! " << endl;
  
  
  // AliCaloRawAnalyzerPeakFinder *p2 =  new AliCaloRawAnalyzerPeakFinder();
  //  p2->PrintVectors();
}


