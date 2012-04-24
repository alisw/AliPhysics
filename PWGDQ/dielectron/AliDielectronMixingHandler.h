#ifndef ALIDIELECTRONMIXINGHANDLER_H
#define ALIDIELECTRONMIXINGHANDLER_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronMixingHandler                  #
//#                                                           #
//#  Authors:                                                 #
//#   Jens Wiechula, Uni-Tübingen / Jens.Wiechula@cern.ch     #
//#                                                           #
//#############################################################

#include <TNamed.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliDielectronVarManager.h"

class AliDielectron;
class AliVTrack;
class AliVEvent;

class AliDielectronMixingHandler : public TNamed {
public:
  enum { kMaxCuts=10 };
  enum EMixType {
    kOSonly=0,
    kOSandLS,
    kAll
  };
  AliDielectronMixingHandler();
  AliDielectronMixingHandler(const char*name, const char* title);

  virtual ~AliDielectronMixingHandler();

  void AddVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins,
                   Double_t min, Double_t max, Bool_t log=kFALSE);
  void AddVariable(AliDielectronVarManager::ValueTypes type, const char* binLimitStr);

  void SetDepth(UShort_t depth) { fDepth=depth; }
  UShort_t GetDepth()     const { return fDepth; }

  void SetMixType(EMixType type) { fMixType=type; }
  EMixType GetMixType() const    { return fMixType; }

  void SetMixUncomplete(Bool_t mix) { fMixIncomplete=mix; }
  Bool_t GetMixUncomplete() const { return fMixIncomplete; }

  void SetMoveToSameVertex(Bool_t move) { fMoveToSameVertex=move; }
  Bool_t GetMoveToSameVertex() const { return fMoveToSameVertex; }

  Int_t GetNumberOfBins() const;
  Int_t FindBin(const Double_t values[], TString *dim=0x0);
  void Fill(const AliVEvent *ev, AliDielectron *diele);

  void MixRemaining(AliDielectron *diele);

  void Init(const AliDielectron *diele=0x0);
  static void MoveToSameVertex(AliVTrack * const vtrack, const Double_t vFirst[3], const Double_t vMix[3]);

private:
  UShort_t     fDepth;     //Number of events per bin to start the merging
  TClonesArray fArrPools; //Array of events in bins

  UShort_t  fEventCuts[kMaxCuts]; //cut variables
  TObjArray fAxes;        //Axis descriptions of the event binning

  EMixType fMixType;      // which combinations to include in the mixing

  Bool_t fMixIncomplete;  // whether to mix uncomplete bins at the end of the processing
  Bool_t fMoveToSameVertex; //whether to move the mixed tracks to the same vertex position

  void DoMixing(TClonesArray &pool, AliDielectron *diele);

  AliDielectronMixingHandler(const AliDielectronMixingHandler &c);
  AliDielectronMixingHandler &operator=(const AliDielectronMixingHandler &c);

  
  ClassDef(AliDielectronMixingHandler,1)         // Dielectron MixingHandler
};



#endif
