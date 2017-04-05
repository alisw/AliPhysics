/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
// AliCEPUtils
// for
// AliAnalysisTaskCEP
//
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  rewritten by
//  Paul Buehler <paul.buehler@oeaw.ac.at>

#ifndef ALICEPUTILS_H
#define ALICEPUTILS_H

#include "AliSPDUtils.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMCEvent.h"
#include "AliGeomManager.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDAD.h"
#include "AliVVZERO.h"

#include "AliCEPBase.h"

class AliCEPUtils : public TObject {

  private:

    Int_t fTPCnclsS;        // Maximum number of accepted shared clusters in TPC
    Double_t fTrackDCA;     // Maximum distance of track to vertex
    Double_t fTrackDCAz;    // Maximum distance of track to vertex in z-direction
    Double_t fTrackEtaMin;  // EtaMin, used in second loop
    Double_t fTrackEtaMax;  // EtaMax, used in second loop
    TList *fTrackCutListPrim; // TList with primary track selection cuts

    TArrayI *v0daughters;   //! 

  public:

  	AliCEPUtils();
  	~AliCEPUtils();
    
    static TH1F* GetHistStatsFlow();
    void SetTPCnclsS(Int_t n = 3) { fTPCnclsS = n; } 
    void SetTrackDCA(Double_t d = 500.) {fTrackDCA = d;} 
    void SetTrackDCAz(Double_t d = 6.) {fTrackDCAz = d;} 
    void SetTrackEtaRange(Double_t min=-0.9, Double_t max = 0.9) {
      fTrackEtaMin = min;
      fTrackEtaMax = max;
    }
    
    void InitTrackCuts(Bool_t IsRun1, Int_t clusterCut);
    void AddTrackCut(AliESDtrackCuts* cut) {
      fTrackCutListPrim->AddLast((TObject*)cut);
    }

  	static Int_t GetEventType(const AliVEvent *Event);
  	static void SPDLoadGeom(Int_t run);
    static void DetermineMCprocessType (
      AliMCEvent *MCEvent, TString MCgen, Int_t &MCproc);
    Int_t AnalyzeTracks(AliESDEvent* ESDEvent,
      TObjArray* fTracks,TArrayI* fTrackStatus);

    // QA studies
    static TList* GetQArnumHists(Int_t rnummin, Int_t rnummax);
    
    static TList* GetBBFlagQAHists();
    void BBFlagAnalysis(
      AliVEvent *Event,
      TList *lhh );
    static TList* GetSPDPileupQAHists();
    void SPDVtxAnalysis (
      AliVEvent *Event,
      Int_t minContributors, 
      Double_t minZdist, 
      Double_t nSigmaZdist, 
      Double_t nSigmaDiamXY, 
      Double_t nSigmaDiamZ,
      TList * lhh );
    static TList* GetnClunTraQAHists();
    void SPDClusterVsTrackletBGAnalysis (
      AliVEvent *Event,
      TList *lhh );
    static TList* GetVtxQAHists();
    void VtxAnalysis (
      AliVEvent *Event,
      TList *lhh );

    UInt_t GetVtxPos(AliVEvent *Event, TVector3 *vtxpos);

    Bool_t checkstatus(UInt_t stat, UInt_t mask, UInt_t pattern);
    Int_t countstatus(TArrayI *stats, UInt_t mask, UInt_t pattern);
    Int_t countstatus(TArrayI *stats, UInt_t mask, UInt_t pattern, TArrayI *indices);
    Int_t countstatus(TArrayI *stats,
      TArrayI *masks, TArrayI *patterns);
    Int_t countstatus(TArrayI *stats,
      TArrayI *masks, TArrayI *patterns, TArrayI* indices);

    Int_t GetCEPTracks(AliESDEvent *ESDEvent, TArrayI *stats, TArrayI* indices);
    Int_t GetResiduals(AliESDEvent* fESDEvent);
    Bool_t TestFiredChips(AliESDEvent *esd, TArrayI *indices);
    
  	ClassDef(AliCEPUtils, 1);

};

#endif
