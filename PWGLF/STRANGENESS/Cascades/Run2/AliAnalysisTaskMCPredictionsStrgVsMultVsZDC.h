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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskMCPredictions.h from David Dobrigkeit Chinellato
//
// --- Francesca Ercolessi
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskMCPredictionsStrgVsMultVsZDC_H
#define AliAnalysisTaskMCPredictionsStrgVsMultVsZDC_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
class TProfile2D;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC();
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const char *name, Float_t lCenterOfMassEnergy = 13000, Bool_t kDoPythia = kTRUE, Bool_t kDoEPOS = kFALSE);
  virtual ~AliAnalysisTaskMCPredictionsStrgVsMultVsZDC();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Double_t Rapidity(Double_t E, Double_t Pz) const;

  Bool_t CheckIsNotPrimary(Bool_t IsFastGenerator, AliStack* MCstack, Int_t iLabelStack, TParticle* part) const;

  Bool_t IsEPOSLHC() const;

  void SetSelectINELgtZERO ( Bool_t lOpt ) { fkSelectINELgtZERO = lOpt; }

  //---------------------------------------------------------------------------------------
private:
  // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
  // your data member object is created on the worker nodes and streaming is not needed.
  // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14

  TList  *fListHist;  //! List of Cascade histograms

  //Histograms
  TH1D *fHistEventCounter; //! histogram for event counting
  TH1D *fHistV0MMult; //!
  TH1D *fHistV0AMult; //!
  TH1D *fHistV0CMult; //!
  TH1D *fHistMult05; //!
  TH1D *fHistMult08; //!
  TH1D *fHistMult10; //!
  TH1D *fHistMult08to15; //!
  TH1D *fHistSPDClusters; //!
  TH1D *fHistSPDCl0; //!
  TH1D *fHistSPDCl1; //!
  TH1D *fHistNMPI; //!
  TH1D *fHistQ2; //!
  TH1D *fHistb; //!
  TH1D *fHistLeadingE; //!
  TH1D *fHistEffEnergy; //!
  TH1D *fHistRxy; //!
  TH2D *f2DHistINELgt0SPDV0M; //!
  TH2D *f2DHistLeadingESPDV0M; //!
  TH2D *f2DHistLeadingERecoPercSPDV0M; //!
  TH2D *f2DHistEffEnergySPDV0M; //!
  TH2D *f2DHistNchSPDV0M; //!
  TH2D *f2DHistNchRecoPercSPDV0M; //!
  TH2D *f2DHistINELgt0RecoPercSPDV0M; //!
  TH2D *f2DHistNMPISPDV0M; //!
  TH2D *f2DHistINELgt0Nch0815V0M; //!
  TH2D *f2DHistLeadingENch0815V0M; //!
  TH2D *f2DHistEffEnergyNch0815V0M; //!
  TH2D *f2DHistNchNch0815V0M; //!
  TH2D *f2DHistNMPINch0815V0M; //!
  TH2D *f2DHistINELgt0Nch0815ZDC; //!
  TH2D *f2DHistLeadingENch0815ZDC;  //!
  TH2D *f2DHistEffEnergyNch0815ZDC; //!
  TH2D *f2DHistNchNch0815ZDC; //!
  TH2D *f2DHistNMPINch0815ZDC; //!
  TH2D *f2DHistQ2SPDV0M; //!
  TH2D *f2DHistbSPDV0M; //!
  TH1D *fHistPt[22]; //!
  TH1D *fHistDecayVtxPos[22];  //!
  TH2D *f2DHistPartSPDV0M[22]; //!
  TH2D *f2DHistPartRecoPercSPDV0M[22]; //!
  TH2D *f2DHistPartNch0815V0M[22]; //!
  TH2D *f2DHistPartNch0815ZDC[22]; //!
  TH2D *f2DHistAvPtSPDV0M[22]; //!
  TH2D *f2dHistZDCVsLE; //!
  TH2D *f2dHistZDCVsEE; //!
  TH2D *f2dHistZDCVsLEA; //!
  TH2D *f2dHistZDCVsLEC; //!
  TH2D *f2dHistZPVsLP; //!
  TH2D *f2dHistZNVsLN; //!
  TH2D *f2dHistZDCVsLEnoacc; //!
  TH2D *f2dHistZPVsLPnoacc; //!
  TH2D *f2dHistZNVsLNnoacc; //!
  TH2D *f2dHistSPDClRecoVsTrue; //!
  TH2D *f2dHistV0MRecoVsTrue; //!
  TH2D *f2dHistTrueVsRecoSPDCl; //!
  TH2D *f2dHistTrueVsRecoSPDCl0; //!
  TH2D *f2dHistTrueVsRecoSPDCl1; //!
  TH3D *f3dHistPi0SPDMultSPDCl; //!
  TProfile2D *p2dHistPi0SPDMultSPDCl; //!
  TH2D *f2dHistTrueVsRecoNch0815; //!

  //Bool
  Bool_t fkSelectINELgtZERO;
  Bool_t fkDoPythia;
  Bool_t fkDoEPOS;

  //Variables
  Float_t fCenterOfMassEnergy;

  //Reco variables
  Float_t fCentrality_V0M; //!
  Float_t fCentrality_SPDClusters; //!
  Float_t fZNApp; //!
  Float_t fZNCpp; //!
  Float_t fZPApp; //!
  Float_t fZPCpp; //!

  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const AliAnalysisTaskMCPredictionsStrgVsMultVsZDC&);            // not implemented
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC& operator=(const AliAnalysisTaskMCPredictionsStrgVsMultVsZDC&); // not implemented

  ClassDef(AliAnalysisTaskMCPredictionsStrgVsMultVsZDC, 1);
};

#endif

