#ifndef ALISINGLETRACKEFFICUTS_H
#define ALISINGLETRACKEFFICUTS_H

/**************************************************************************
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


/*__|______________________________________________________________________________|
 |                              -----Info(i)-----                                  |
 |                                                                                 |
 |  .h of cut class for single track  efficiecy                                    |
 |                                                                                 |
 |   ESDs<-->AODs (ON/OFF)                                                         |
 |                                                      Authors:                   |
 |_____________________________________________________________________________|___*/



#include <TString.h>
#include "TObject.h"

class AliVEvent;
#include "AliSingleTrackEffCuts.h"


class AliSingleTrackEffCuts : public TObject 
{
 public:

  AliSingleTrackEffCuts();
  
  virtual ~AliSingleTrackEffCuts() {;}
  
  AliSingleTrackEffCuts(const AliSingleTrackEffCuts& source);
  AliSingleTrackEffCuts& operator=(const AliSingleTrackEffCuts& source);

  Bool_t IsMCEventSelected(TObject *obj);
  Bool_t IsRecoEventSelected(TObject *obj);

  Bool_t IsMCParticleGenerated(TObject *obj);//truth
  Bool_t IsMCParticleInKineAcceptance(TObject *obj);//truth
  Bool_t IsMCParticleInReconstructable(TObject *obj);//mc truth
  Bool_t IsRecoParticleKineAcceptance(TObject *obj);//reco

    
    
  // Setters
  void SetEtaRange(Float_t etamin, Float_t etamax){ fEtaMin=etamin; fEtaMax=etamax; }
  void SetYRange(Float_t ymin, Float_t ymax){ fYMin=ymin; fYMax=ymax; }
  void SetPtRange(Float_t ptmin, Float_t ptmax){fPtMin=ptmin; fPtMax=ptmax; }

  void SetPdgCode(Int_t pdgCode){ fPdgCode = pdgCode; fIsPdgCode=kTRUE; }
  void SetIsCharged(Bool_t charge){ fIsCharged=charge; }

  void SetMinVtxType(Int_t type=3) {fMinVtxType=type;}  
  void SetUseEventsWithOnlySPDVertex(Bool_t flag=kTRUE){ 
    if(flag) fMinVtxType=1;
    else fMinVtxType=3;
  }

  void SetIsAOD(Bool_t flag){ fisAOD = flag; }
  Bool_t IsAOD(){ return fisAOD; }

  void SetMinVtxContr(Int_t contr=1) {fMinVtxContr=contr;} 
  void SetMaxVtxZ(Float_t z=1e6) {fMaxVtxZ=z;}   

  void SetTriggerMask(ULong64_t mask=0) { fTriggerMask=mask; }
  UInt_t GetTriggerMask(){ return fTriggerMask; }

  void SetRequireVtxCuts(Bool_t vtx=kFALSE) {fRequireVtxCuts=vtx;} // cut values setter
  Bool_t   GetRequireVtxCuts() const {return fRequireVtxCuts;} // cut value getter

  void SetNumberOfClusters(Int_t nITS, Int_t nTPC, Int_t nTOF, Int_t nMUON){
    fnClusITS = nITS; fnClusTPC = nTPC; fnClusTOF = nTOF; fnClusMUON = nMUON;
  }
 
 protected:

  Bool_t IsVertexSelected(AliVEvent *event);

  Bool_t fisAOD;  // flag wether it is AOD:1 or ESD:0 analysis
  
  Bool_t fIsPdgCode;
  Int_t fPdgCode;

  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fYMin;
  Float_t fYMax;
  Float_t fPtMin;
  Float_t fPtMax;
  Bool_t  fIsCharged;
  Bool_t  fRequireVtxCuts ; //The type of trigger to be checked
  AliMCEvent  * fMCinfo;          //! MC event handler

  UInt_t fTriggerMask;    // trigger mask

  Int_t   fMinVtxType;  // 0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks
  Int_t   fMinVtxContr; // minimum vertex contributors
  Float_t fMaxVtxZ;   // maximum |z| of primary vertex

  Int_t fnClusITS;
  Int_t fnClusTPC;
  Int_t fnClusTOF;
  Int_t fnClusMUON;

  Int_t fCutOnZVertexSPD; // 0: no cut, 1: |zvtx-SPD - zvtx-TPC|<0.5cm


  ClassDef(AliSingleTrackEffCuts,1)  // base class for cuts on AOD reconstructed heavy-flavour decays
};

#endif
