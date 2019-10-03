/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// The task:
// stores TPC TRD matching quantities in a THnSparse
//
//  Author:
//  Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//

#ifndef ALITRDPIDMATCHING_H
#define ALITRDPIDMATCHING_H

#include "AliAnalysisTaskSE.h"
#include "AliTRDPIDTree.h"

class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;
class TFile;
class AliESDEvent;
class AliESDpid;
class AliESD;
class AliESDtrackCuts;
class AliAnalysisTask;
class AliESDInputHandler;
class AliAnalysisManager;
class AliCentrality;
class TTree;
class TSystem;
class TStyle;
class TROOT;
class Riostream;
class TChain;



class AliTRDPIDmatching : public AliTRDPIDTree {
 public:
  AliTRDPIDmatching(const char *name= "trd_pid_matching");
  virtual ~AliTRDPIDmatching();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   ProcessData(AliESDEvent *const esdEvent=0);
  virtual void   Terminate(const Option_t*);
  Int_t          CompareFloat(Float_t f1=1, Float_t f2=0) const;

 private:
  AliESDEvent *fESD;      //! ESD event
  AliMCEvent  *fMC;       //! MC event
  THnSparseF  *fTHntrdmatch;  //! THnSparse containing the data

  TObjArray   *fOutputContainer; //! output data container
  AliESDtrackCuts *fesdTrackCuts; //! ESD track cuts
  Bool_t fHasMC;          //! MC boolean
   
  AliTRDPIDmatching(const AliTRDPIDmatching&); // not implemented
  AliTRDPIDmatching& operator=(const AliTRDPIDmatching&); // not implemented
  
  ClassDef(AliTRDPIDmatching, 1); // example of analysis
};

#endif
