#ifndef ALIFMDQATASK_H
#define ALIFMDQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//_________________________________________________________________________
// An analysis task to check the FMD data in simulated data
//
//*-- Hans Hjersing Dalsgaard 
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ;  

class AliFMDQATask : public AliAnalysisTask {

public:
  AliFMDQATask(const char *name) ;
  virtual ~AliFMDQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:

  void   RingSelector(const UShort_t detector, const Char_t ring, const Float_t mult) const ;
  Bool_t TestHisto(TH1D * hTest) const ;
  void   FitAll(TH1D* hTest, Float_t &chiSq, Int_t &ndf, Float_t &mpv, Float_t chiMax, Float_t chiLow ) const ;

  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms

  TH1D * fhFMD1i ;//!
  TH1D * fhFMD2i ;//! 
  TH1D * fhFMD2o ;//! 
  TH1D * fhFMD3i ;//! 
  TH1D * fhFMD3o ;//! 
  
  ClassDef(AliFMDQATask, 0); // a FMD photon analysis task 
};
#endif // ALIFMDQATASK_H
