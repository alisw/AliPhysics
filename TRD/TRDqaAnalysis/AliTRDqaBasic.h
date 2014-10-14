#ifndef ALITRDQABASIC_H
#define ALITRDQABASIC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliTRDqaBasic.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//
// This class is a part of a package of high level QA monitoring for TRD.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class TTree; 
class AliESDEvent; 
class TH1D; 
class TH2D;
class AliExternalTrackParam;
class AliESDtrackCuts;

class AliTRDqaBasic : public AliAnalysisTask {

 public:

  AliTRDqaBasic();
  AliTRDqaBasic(const char *name);
  AliTRDqaBasic(const AliTRDqaBasic & trd);
  AliTRDqaBasic &operator=(const AliTRDqaBasic & /*g*/) { return *this; };
  virtual ~AliTRDqaBasic() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

 private:
 
  TTree        * fChain;             //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD;               //! Declaration of leave types

  TObjArray * fOutputContainer;      //! output data container
  
  // cuts
  AliESDtrackCuts *fTrackCuts;       // Cuts

  // histograms

  TH1D *fStatus;                     // track status

  TH1D *fnTracks;                    // histogram
  TH1D *fPtIn;                       // histogram
  TH1D *fPtOut;                      // histogram
  TH1D *fPtVtx;                      // histogram
  TH1D *fPtVtxSec;                   // histogram

  TH2D *fPtPt;                       // histogram

  ClassDef(AliTRDqaBasic, 0);        // a TRD analysis task 

};
#endif // ALITRDQBASIC_H
