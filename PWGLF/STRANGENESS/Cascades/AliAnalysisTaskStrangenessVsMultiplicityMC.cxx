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
// This task is meant to explore the possibility of using a VZERO amplitude
// based multiplicity estimator for proton-proton collisions. For this, two 
// main operation methods for this task are foreseen: 
// 
//  1) (under development) it should act as an auxiliary task and provide a 
//     calibrated estimator 
//
//  2) "Debug mode" which will also create a ROOT TTree object with event 
//     by event info potentially used for exploration / calibration. This 
//     includes the following info: 
//    
//      --- All VZERO Amplitudes (saved as Float_t) 
//      --- (optional) time for each channel
//      --- (optional) time width for each channel 
//      --- GetReferenceMultiplicity Estimator, |eta|<0.5 
//      --- GetReferenceMultiplicity Estimator, |eta|<0.8 
//      --- (if MC) True Multiplicity, |eta|<0.5
//      --- (if MC) True Multiplicity,  2.8 < eta < 5.8 (VZEROA region)
//      --- (if MC) True Multiplicity, -3.7 < eta <-1.7 (VZEROC region)
//      --- Run Number
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskStrangenessVsMultiplicityMC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessVsMultiplicityMC)

AliAnalysisTaskStrangenessVsMultiplicityMC::AliAnalysisTaskStrangenessVsMultiplicityMC()
  : AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), 
  fkSaveV0Tree      ( kFALSE ),
  fkSaveCascadeTree ( kTRUE  ),
  fkRunVertexers    ( kTRUE  ), 
  //---> Variables for fTreeEvent
  fAmplitude_V0A   (0),   
  fAmplitude_V0C   (0),   
  fAmplitude_V0AEq (0),   
  fAmplitude_V0CEq (0),  
  fCentrality_V0A(0), 
  fCentrality_V0C(0), 
  fCentrality_V0M(0), 
  fCentrality_V0AEq(0), 
  fCentrality_V0CEq(0), 
  fCentrality_V0MEq(0), 
  fRefMultEta5(0),
  fRefMultEta8(0),
  fTrueMultEta5(0),
  fTrueMultEta8(0),
  fTrueMultVZEROA(0),
  fTrueMultVZEROC(0),
  fRunNumber(0),
  //---> Variables for fTreeCascade
  fTreeCascVarCharge(0), 
	fTreeCascVarMassAsXi(0),
	fTreeCascVarMassAsOmega(0),
	fTreeCascVarPt(0),
	fTreeCascVarPtMC(0),
	fTreeCascVarRapXi(0),
	fTreeCascVarRapOmega(0),
	fTreeCascVarRapMC(0),
	fTreeCascVarNegEta(0),
	fTreeCascVarPosEta(0),
	fTreeCascVarBachEta(0),
	fTreeCascVarDCACascDaughters(0),
	fTreeCascVarDCABachToPrimVtx(0),
	fTreeCascVarDCAV0Daughters(0),
	fTreeCascVarDCAV0ToPrimVtx(0),
	fTreeCascVarDCAPosToPrimVtx(0),
	fTreeCascVarDCANegToPrimVtx(0),
	fTreeCascVarCascCosPointingAngle(0),
	fTreeCascVarCascRadius(0),
	fTreeCascVarV0Mass(0),
	fTreeCascVarV0CosPointingAngle(0),
	fTreeCascVarV0CosPointingAngleSpecial(0),
	fTreeCascVarV0Radius(0),
  fTreeCascVarLeastNbrClusters(0),
	fTreeCascVarDistOverTotMom(0),
	fTreeCascVarNegNSigmaPion(0),
	fTreeCascVarNegNSigmaProton(0),
	fTreeCascVarPosNSigmaPion(0),
	fTreeCascVarPosNSigmaProton(0),
	fTreeCascVarBachNSigmaPion(0),
	fTreeCascVarBachNSigmaKaon(0),
	fTreeCascVarCentV0A(0),
	fTreeCascVarCentV0C(0),
	fTreeCascVarCentV0M(0),
	fTreeCascVarCentV0AEq(0),
	fTreeCascVarCentV0CEq(0),
	fTreeCascVarCentV0MEq(0),
	fTreeCascVarAmpV0A(0),
	fTreeCascVarAmpV0C(0),
	fTreeCascVarAmpV0AEq(0),
	fTreeCascVarAmpV0CEq(0),
  fTreeCascVarRefMultEta8(0),
  fTreeCascVarRefMultEta5(0),
  fTreeCascVarTrueMultEta5(0),
  fTreeCascVarTrueMultEta8(0),
  fTreeCascVarTrueMultVZEROA(0),
  fTreeCascVarTrueMultVZEROC(0),
  fTreeCascVarIsPhysicalPrimary(0), 
  fTreeCascVarPID(0), 
  fTreeCascVarRunNumber(0), 
  //---> Histograms
  fHistEventCounter(0), 
  //---> MC Generated Histo (analysis level) 
	fHistPt_GenXiMinus(0),
	fHistPt_GenXiPlus(0),
	fHistPt_GenOmegaMinus(0),
	fHistPt_GenOmegaPlus(0),

  //VsRefMult
	fHistPtVsRefMultEta5_GenXiMinus(0),
	fHistPtVsRefMultEta5_GenXiPlus(0),
	fHistPtVsRefMultEta5_GenOmegaMinus(0),
	fHistPtVsRefMultEta5_GenOmegaPlus(0),
	fHistPtVsRefMultEta8_GenXiMinus(0),
	fHistPtVsRefMultEta8_GenXiPlus(0),
	fHistPtVsRefMultEta8_GenOmegaMinus(0),
	fHistPtVsRefMultEta8_GenOmegaPlus(0),

  //VsCentralities
	fHistPtVsCentV0A_GenXiMinus(0),
	fHistPtVsCentV0A_GenXiPlus(0),
	fHistPtVsCentV0A_GenOmegaMinus(0),
	fHistPtVsCentV0A_GenOmegaPlus(0),
	fHistPtVsCentV0C_GenXiMinus(0),
	fHistPtVsCentV0C_GenXiPlus(0),
	fHistPtVsCentV0C_GenOmegaMinus(0),
	fHistPtVsCentV0C_GenOmegaPlus(0),
	fHistPtVsCentV0M_GenXiMinus(0),
	fHistPtVsCentV0M_GenXiPlus(0),
	fHistPtVsCentV0M_GenOmegaMinus(0),
	fHistPtVsCentV0M_GenOmegaPlus(0),

  //Equalized
	fHistPtVsCentV0AEq_GenXiMinus(0),
	fHistPtVsCentV0AEq_GenXiPlus(0),
	fHistPtVsCentV0AEq_GenOmegaMinus(0),
	fHistPtVsCentV0AEq_GenOmegaPlus(0),
	fHistPtVsCentV0CEq_GenXiMinus(0),
	fHistPtVsCentV0CEq_GenXiPlus(0),
	fHistPtVsCentV0CEq_GenOmegaMinus(0),
	fHistPtVsCentV0CEq_GenOmegaPlus(0),
	fHistPtVsCentV0MEq_GenXiMinus(0),
	fHistPtVsCentV0MEq_GenXiPlus(0),
	fHistPtVsCentV0MEq_GenOmegaMinus(0),
	fHistPtVsCentV0MEq_GenOmegaPlus(0),

  //VsAmp
	fHistPtVsAmpV0A_GenXiMinus(0),
	fHistPtVsAmpV0A_GenXiPlus(0),
	fHistPtVsAmpV0A_GenOmegaMinus(0),
	fHistPtVsAmpV0A_GenOmegaPlus(0),
	fHistPtVsAmpV0C_GenXiMinus(0),
	fHistPtVsAmpV0C_GenXiPlus(0),
	fHistPtVsAmpV0C_GenOmegaMinus(0),
	fHistPtVsAmpV0C_GenOmegaPlus(0),
	fHistPtVsAmpV0M_GenXiMinus(0),
	fHistPtVsAmpV0M_GenXiPlus(0),
	fHistPtVsAmpV0M_GenOmegaMinus(0),
	fHistPtVsAmpV0M_GenOmegaPlus(0),
  //Equalized Amps
	fHistPtVsAmpV0AEq_GenXiMinus(0),
	fHistPtVsAmpV0AEq_GenXiPlus(0),
	fHistPtVsAmpV0AEq_GenOmegaMinus(0),
	fHistPtVsAmpV0AEq_GenOmegaPlus(0),
	fHistPtVsAmpV0CEq_GenXiMinus(0),
	fHistPtVsAmpV0CEq_GenXiPlus(0),
	fHistPtVsAmpV0CEq_GenOmegaMinus(0),
	fHistPtVsAmpV0CEq_GenOmegaPlus(0),
	fHistPtVsAmpV0MEq_GenXiMinus(0),
	fHistPtVsAmpV0MEq_GenXiPlus(0),
	fHistPtVsAmpV0MEq_GenOmegaMinus(0),
	fHistPtVsAmpV0MEq_GenOmegaPlus(0)  

//------------------------------------------------
// Tree Variables 
{

}

AliAnalysisTaskStrangenessVsMultiplicityMC::AliAnalysisTaskStrangenessVsMultiplicityMC(const char *name) 
  : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), 
  fkSaveV0Tree      ( kFALSE ),
  fkSaveCascadeTree ( kTRUE  ), 
  fkRunVertexers    ( kTRUE  ),
  //---> Variables for fTreeEvent
  fAmplitude_V0A (0),   
  fAmplitude_V0C (0), 
  fAmplitude_V0AEq (0),   
  fAmplitude_V0CEq (0), 
  fCentrality_V0A(0), 
  fCentrality_V0C(0), 
  fCentrality_V0M(0), 
  fCentrality_V0AEq(0), 
  fCentrality_V0CEq(0), 
  fCentrality_V0MEq(0), 
  fRefMultEta5(0),
  fRefMultEta8(0),
  fTrueMultEta5(0),
  fTrueMultEta8(0),
  fTrueMultVZEROA(0),
  fTrueMultVZEROC(0),
  fRunNumber(0),
  //---> Variables for fTreeCascade
  fTreeCascVarCharge(0), 
	fTreeCascVarMassAsXi(0),
	fTreeCascVarMassAsOmega(0),
	fTreeCascVarPt(0),
	fTreeCascVarPtMC(0),
	fTreeCascVarRapXi(0),
	fTreeCascVarRapOmega(0),
	fTreeCascVarRapMC(0),
	fTreeCascVarNegEta(0),
	fTreeCascVarPosEta(0),
	fTreeCascVarBachEta(0),
	fTreeCascVarDCACascDaughters(0),
	fTreeCascVarDCABachToPrimVtx(0),
	fTreeCascVarDCAV0Daughters(0),
	fTreeCascVarDCAV0ToPrimVtx(0),
	fTreeCascVarDCAPosToPrimVtx(0),
	fTreeCascVarDCANegToPrimVtx(0),
	fTreeCascVarCascCosPointingAngle(0),
	fTreeCascVarCascRadius(0),
	fTreeCascVarV0Mass(0),
	fTreeCascVarV0CosPointingAngle(0),
	fTreeCascVarV0CosPointingAngleSpecial(0),
	fTreeCascVarV0Radius(0),
  fTreeCascVarLeastNbrClusters(0),
	fTreeCascVarDistOverTotMom(0),
	fTreeCascVarNegNSigmaPion(0),
	fTreeCascVarNegNSigmaProton(0),
	fTreeCascVarPosNSigmaPion(0),
	fTreeCascVarPosNSigmaProton(0),
	fTreeCascVarBachNSigmaPion(0),
	fTreeCascVarBachNSigmaKaon(0),
	fTreeCascVarCentV0A(0),
	fTreeCascVarCentV0C(0),
	fTreeCascVarCentV0M(0),
	fTreeCascVarCentV0AEq(0),
	fTreeCascVarCentV0CEq(0),
	fTreeCascVarCentV0MEq(0),
	fTreeCascVarAmpV0A(0),
	fTreeCascVarAmpV0C(0),
	fTreeCascVarAmpV0AEq(0),
	fTreeCascVarAmpV0CEq(0),
  fTreeCascVarRefMultEta8(0),
  fTreeCascVarRefMultEta5(0),
  fTreeCascVarTrueMultEta5(0),
  fTreeCascVarTrueMultEta8(0),
  fTreeCascVarTrueMultVZEROA(0),
  fTreeCascVarTrueMultVZEROC(0),
  fTreeCascVarIsPhysicalPrimary(0), 
  fTreeCascVarPID(0), 
  fTreeCascVarRunNumber(0), 
  //---> Histograms
  fHistEventCounter(0), 
  //---> MC Generated Histo (analysis level) 
	fHistPt_GenXiMinus(0),
	fHistPt_GenXiPlus(0),
	fHistPt_GenOmegaMinus(0),
	fHistPt_GenOmegaPlus(0),

  //VsRefMult
	fHistPtVsRefMultEta5_GenXiMinus(0),
	fHistPtVsRefMultEta5_GenXiPlus(0),
	fHistPtVsRefMultEta5_GenOmegaMinus(0),
	fHistPtVsRefMultEta5_GenOmegaPlus(0),
	fHistPtVsRefMultEta8_GenXiMinus(0),
	fHistPtVsRefMultEta8_GenXiPlus(0),
	fHistPtVsRefMultEta8_GenOmegaMinus(0),
	fHistPtVsRefMultEta8_GenOmegaPlus(0),

  //VsCentralities
	fHistPtVsCentV0A_GenXiMinus(0),
	fHistPtVsCentV0A_GenXiPlus(0),
	fHistPtVsCentV0A_GenOmegaMinus(0),
	fHistPtVsCentV0A_GenOmegaPlus(0),
	fHistPtVsCentV0C_GenXiMinus(0),
	fHistPtVsCentV0C_GenXiPlus(0),
	fHistPtVsCentV0C_GenOmegaMinus(0),
	fHistPtVsCentV0C_GenOmegaPlus(0),
	fHistPtVsCentV0M_GenXiMinus(0),
	fHistPtVsCentV0M_GenXiPlus(0),
	fHistPtVsCentV0M_GenOmegaMinus(0),
	fHistPtVsCentV0M_GenOmegaPlus(0),

  //Equalized
	fHistPtVsCentV0AEq_GenXiMinus(0),
	fHistPtVsCentV0AEq_GenXiPlus(0),
	fHistPtVsCentV0AEq_GenOmegaMinus(0),
	fHistPtVsCentV0AEq_GenOmegaPlus(0),
	fHistPtVsCentV0CEq_GenXiMinus(0),
	fHistPtVsCentV0CEq_GenXiPlus(0),
	fHistPtVsCentV0CEq_GenOmegaMinus(0),
	fHistPtVsCentV0CEq_GenOmegaPlus(0),
	fHistPtVsCentV0MEq_GenXiMinus(0),
	fHistPtVsCentV0MEq_GenXiPlus(0),
	fHistPtVsCentV0MEq_GenOmegaMinus(0),
	fHistPtVsCentV0MEq_GenOmegaPlus(0),

  //VsAmp
	fHistPtVsAmpV0A_GenXiMinus(0),
	fHistPtVsAmpV0A_GenXiPlus(0),
	fHistPtVsAmpV0A_GenOmegaMinus(0),
	fHistPtVsAmpV0A_GenOmegaPlus(0),
	fHistPtVsAmpV0C_GenXiMinus(0),
	fHistPtVsAmpV0C_GenXiPlus(0),
	fHistPtVsAmpV0C_GenOmegaMinus(0),
	fHistPtVsAmpV0C_GenOmegaPlus(0),
	fHistPtVsAmpV0M_GenXiMinus(0),
	fHistPtVsAmpV0M_GenXiPlus(0),
	fHistPtVsAmpV0M_GenOmegaMinus(0),
	fHistPtVsAmpV0M_GenOmegaPlus(0),
  //Equalized Amps
	fHistPtVsAmpV0AEq_GenXiMinus(0),
	fHistPtVsAmpV0AEq_GenXiPlus(0),
	fHistPtVsAmpV0AEq_GenOmegaMinus(0),
	fHistPtVsAmpV0AEq_GenOmegaPlus(0),
	fHistPtVsAmpV0CEq_GenXiMinus(0),
	fHistPtVsAmpV0CEq_GenXiPlus(0),
	fHistPtVsAmpV0CEq_GenOmegaMinus(0),
	fHistPtVsAmpV0CEq_GenOmegaPlus(0),
	fHistPtVsAmpV0MEq_GenXiMinus(0),
	fHistPtVsAmpV0MEq_GenXiPlus(0),
	fHistPtVsAmpV0MEq_GenOmegaMinus(0),
	fHistPtVsAmpV0MEq_GenOmegaPlus(0)
{

  //Re-vertex: Will only apply for cascade candidates

  fV0VertexerSels[0] =  33.  ;  // max allowed chi2
  fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
  fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
  fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
  fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
  fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
  fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
  fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
  fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
  fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
  fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
  fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
  fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
  fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        

  DefineOutput(1, TList::Class()); // Event Counter Histo
  DefineOutput(2, TTree::Class()); // Event Tree
  DefineOutput(3, TTree::Class()); // V0 Tree
  DefineOutput(4, TTree::Class()); // Cascade Tree
}


AliAnalysisTaskStrangenessVsMultiplicityMC::~AliAnalysisTaskStrangenessVsMultiplicityMC()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHist){
      delete fListHist;
      fListHist = 0x0;
   }
   if (fTreeEvent){
      delete fTreeEvent;
      fTreeEvent = 0x0;
   }
   if (fTreeV0){
      delete fTreeV0;
      fTreeV0 = 0x0;
   }
   if (fTreeCascade){
      delete fTreeCascade;
      fTreeCascade = 0x0;
   }
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMC::UserCreateOutputObjects()
{

   OpenFile(2);	
   // Called once

//------------------------------------------------

   fTreeEvent = new TTree("fTreeEvent","Event");

//------------------------------------------------
// fTree Branch definitions - Event by Event info
//------------------------------------------------

//-----------BASIC-INFO---------------------------

  //--- VZERO Data (Integrated)
  fTreeEvent->Branch("fAmplitude_V0A",&fAmplitude_V0A,"fAmplitude_V0A/F");
  fTreeEvent->Branch("fAmplitude_V0C",&fAmplitude_V0C,"fAmplitude_V0C/F");
  fTreeEvent->Branch("fAmplitude_V0AEq",&fAmplitude_V0AEq,"fAmplitude_V0AEq/F");
  fTreeEvent->Branch("fAmplitude_V0CEq",&fAmplitude_V0CEq,"fAmplitude_V0CEq/F");

  //Info from AliCentrality (not necessarily 'centrality' per se) 
  fTreeEvent->Branch("fCentrality_V0A",&fCentrality_V0A,"fCentrality_V0A/F");
  fTreeEvent->Branch("fCentrality_V0C",&fCentrality_V0C,"fCentrality_V0C/F");
  fTreeEvent->Branch("fCentrality_V0M",&fCentrality_V0A,"fCentrality_V0M/F");
  fTreeEvent->Branch("fCentrality_V0AEq",&fCentrality_V0AEq,"fCentrality_V0AEq/F");
  fTreeEvent->Branch("fCentrality_V0CEq",&fCentrality_V0CEq,"fCentrality_V0CEq/F");
  fTreeEvent->Branch("fCentrality_V0MEq",&fCentrality_V0AEq,"fCentrality_V0MEq/F");
  
  //Official GetReferenceMultiplicity
  fTreeEvent->Branch("fRefMultEta5",&fRefMultEta5,"fRefMultEta5/I");
  fTreeEvent->Branch("fRefMultEta8",&fRefMultEta8,"fRefMultEta8/I");

  fTreeEvent->Branch("fTrueMultEta5",&fTrueMultEta5,"fTrueMultEta5/I");
  fTreeEvent->Branch("fTrueMultEta8",&fTrueMultEta8,"fTrueMultEta8/I");
  fTreeEvent->Branch("fTrueMultVZEROA",&fTrueMultVZEROA,"fTrueMultVZEROA/I");
  fTreeEvent->Branch("fTrueMultVZEROC",&fTrueMultVZEROC,"fTrueMultVZEROC/I");

  //Run Number
  fTreeEvent->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");

  //Create Basic V0 Output Tree
  fTreeV0 = new TTree( "fTreeV0", "V0 Candidates");

  //Create Cascade output tree
  fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");

//------------------------------------------------
// fTreeCascade Branch definitions - Cascade Tree
//------------------------------------------------

//-----------BASIC-INFO---------------------------
	fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");	
  fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
  fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
  fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
  fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
  fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
  fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
  fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");
  fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
  fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
  fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
//-----------INFO-FOR-CUTS------------------------
  fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
  fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
  fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
  fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
  fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
  fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
  fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
  fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
  fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
//-----------MULTIPLICITY-INFO--------------------
  fTreeCascade->Branch("fTreeCascVarCentV0A",&fTreeCascVarCentV0A,"fTreeCascVarCentV0A/F");
  fTreeCascade->Branch("fTreeCascVarCentV0C",&fTreeCascVarCentV0C,"fTreeCascVarCentV0C/F");
  fTreeCascade->Branch("fTreeCascVarCentV0M",&fTreeCascVarCentV0M,"fTreeCascVarCentV0M/F");
  fTreeCascade->Branch("fTreeCascVarCentV0AEq",&fTreeCascVarCentV0AEq,"fTreeCascVarCentV0AEq/F");
  fTreeCascade->Branch("fTreeCascVarCentV0CEq",&fTreeCascVarCentV0CEq,"fTreeCascVarCentV0CEq/F");
  fTreeCascade->Branch("fTreeCascVarCentV0MEq",&fTreeCascVarCentV0MEq,"fTreeCascVarCentV0MEq/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0A",&fTreeCascVarAmpV0A,"fTreeCascVarAmpV0A/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0C",&fTreeCascVarAmpV0C,"fTreeCascVarAmpV0C/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0AEq",&fTreeCascVarAmpV0AEq,"fTreeCascVarAmpV0AEq/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0CEq",&fTreeCascVarAmpV0CEq,"fTreeCascVarAmpV0CEq/F");
  fTreeCascade->Branch("fTreeCascVarRefMultEta8",&fTreeCascVarRefMultEta8,"fTreeCascVarRefMultEta8/I");
  fTreeCascade->Branch("fTreeCascVarRefMultEta5",&fTreeCascVarRefMultEta5,"fTreeCascVarRefMultEta5/I");
  fTreeCascade->Branch("fTreeCascVarTrueMultEta5",&fTreeCascVarTrueMultEta5,"fTreeCascVarTrueMultEta5/I");
  fTreeCascade->Branch("fTreeCascVarTrueMultEta8",&fTreeCascVarTrueMultEta8,"fTreeCascVarTrueMultEta8/I");
  fTreeCascade->Branch("fTreeCascVarTrueMultVZEROA",&fTreeCascVarTrueMultVZEROA,"fTreeCascVarTrueMultVZEROA/I");
  fTreeCascade->Branch("fTreeCascVarTrueMultVZEROC",&fTreeCascVarTrueMultVZEROC,"fTreeCascVarTrueMultVZEROC/I");
  fTreeCascade->Branch("fTreeCascVarIsPhysicalPrimary",&fTreeCascVarIsPhysicalPrimary,"fTreeCascVarIsPhysicalPrimary/I");
  fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
  fTreeCascade->Branch("fTreeCascVarRunNumber",&fTreeCascVarRunNumber,"fTreeCascVarRunNumber/I");
//-----------DECAY-LENGTH-INFO--------------------
  fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
//------------------------------------------------
  fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
  fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
  fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");

//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------
  
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();

  // Multiplicity
  if(! fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
  }

//------------------------------------------------
// V0 Multiplicity Histograms
//------------------------------------------------

   // Create histograms
   OpenFile(1);
   fListHist = new TList();
   fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",5,0,5); 
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "Phys-Sel");  
    fHistEventCounter->GetXaxis()->SetBinLabel(3, "Has Vtx");  
    fHistEventCounter->GetXaxis()->SetBinLabel(4, "Vtx |z|<10cm");  
    fHistEventCounter->GetXaxis()->SetBinLabel(5, "Isn't Pileup");
    fListHist->Add(fHistEventCounter); 
   }

  //Histograms for Efficiency corrections... a bunch of them 
  //1D Histograms - Fine if efficiency doesn't change vs mult (expected) 
  //---> Always filled for |y|<0.5 
  if(! fHistPt_GenXiMinus ) {
    fHistPt_GenXiMinus    = new TH1D( "fHistPt_GenXiMinus",    "Generated;p_{T} (GeV/c)",200,0,20);   fListHist->Add(fHistPt_GenXiMinus);    }
  if(! fHistPt_GenXiPlus ) {
     fHistPt_GenXiPlus     = new TH1D( "fHistPt_GenXiPlus",    "Generated;p_{T} (GeV/c)",200,0,20);   fListHist->Add(fHistPt_GenXiPlus);     }
  if(! fHistPt_GenOmegaMinus ) {
    fHistPt_GenOmegaMinus = new TH1D( "fHistPt_GenOmegaMinus", "Generated;p_{T} (GeV/c)",200,0,20);   fListHist->Add(fHistPt_GenOmegaMinus); }
  if(! fHistPt_GenOmegaPlus ) {
    fHistPt_GenOmegaPlus  = new TH1D( "fHistPt_GenOmegaPlus",  "Generated;p_{T} (GeV/c)",200,0,20);   fListHist->Add(fHistPt_GenOmegaPlus);  }
  //2D Histos for vs Mult calculation 
  if(! fHistPtVsRefMultEta5_GenXiMinus ) {
    fHistPtVsRefMultEta5_GenXiMinus    = new TH2D( "fHistPtVsRefMultEta5_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta5_GenXiMinus);    }
  if(! fHistPtVsRefMultEta5_GenXiPlus ) {
    fHistPtVsRefMultEta5_GenXiPlus     = new TH2D( "fHistPtVsRefMultEta5_GenXiPlus",        "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta5_GenXiPlus);    }
  if(! fHistPtVsRefMultEta5_GenOmegaMinus ) {
    fHistPtVsRefMultEta5_GenOmegaMinus    = new TH2D( "fHistPtVsRefMultEta5_GenOmegaMinus", "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta5_GenOmegaMinus);    }
  if(! fHistPtVsRefMultEta5_GenOmegaPlus ) {
    fHistPtVsRefMultEta5_GenOmegaPlus     = new TH2D( "fHistPtVsRefMultEta5_GenOmegaPlus",  "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta5_GenOmegaPlus);    }
  if(! fHistPtVsRefMultEta8_GenXiMinus ) {
    fHistPtVsRefMultEta8_GenXiMinus    = new TH2D( "fHistPtVsRefMultEta8_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta8_GenXiMinus);    }
  if(! fHistPtVsRefMultEta8_GenXiPlus ) {
    fHistPtVsRefMultEta8_GenXiPlus     = new TH2D( "fHistPtVsRefMultEta8_GenXiPlus",        "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta8_GenXiPlus);    }
  if(! fHistPtVsRefMultEta8_GenOmegaMinus ) {
    fHistPtVsRefMultEta8_GenOmegaMinus    = new TH2D( "fHistPtVsRefMultEta8_GenOmegaMinus", "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta8_GenOmegaMinus);    }
  if(! fHistPtVsRefMultEta8_GenOmegaPlus ) {
    fHistPtVsRefMultEta8_GenOmegaPlus     = new TH2D( "fHistPtVsRefMultEta8_GenOmegaPlus",  "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   fListHist->Add(fHistPtVsRefMultEta8_GenOmegaPlus);    }

  //Centralities: V0A, V0C, V0M, +Eq
  if(! fHistPtVsCentV0A_GenXiMinus ) {
    fHistPtVsCentV0A_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0A_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0A_GenXiMinus); }
  if(! fHistPtVsCentV0A_GenXiPlus ) {
    fHistPtVsCentV0A_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0A_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0A_GenXiPlus); }
  if(! fHistPtVsCentV0A_GenOmegaMinus ) {
    fHistPtVsCentV0A_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0A_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0A_GenOmegaMinus); }
  if(! fHistPtVsCentV0A_GenOmegaPlus ) {
    fHistPtVsCentV0A_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0A_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0A_GenOmegaPlus); }  

  if(! fHistPtVsCentV0C_GenXiMinus ) {
    fHistPtVsCentV0C_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0C_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0C_GenXiMinus); }
  if(! fHistPtVsCentV0C_GenXiPlus ) {
    fHistPtVsCentV0C_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0C_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0C_GenXiPlus); }
  if(! fHistPtVsCentV0C_GenOmegaMinus ) {
    fHistPtVsCentV0C_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0C_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0C_GenOmegaMinus); }
  if(! fHistPtVsCentV0C_GenOmegaPlus ) {
    fHistPtVsCentV0C_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0C_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0C_GenOmegaPlus); }  

  if(! fHistPtVsCentV0M_GenXiMinus ) {
    fHistPtVsCentV0M_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0M_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0M_GenXiMinus); }
  if(! fHistPtVsCentV0M_GenXiPlus ) {
    fHistPtVsCentV0M_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0M_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0M_GenXiPlus); }
  if(! fHistPtVsCentV0M_GenOmegaMinus ) {
    fHistPtVsCentV0M_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0M_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0M_GenOmegaMinus); }
  if(! fHistPtVsCentV0M_GenOmegaPlus ) {
    fHistPtVsCentV0M_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0M_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0M_GenOmegaPlus); }  

  //Equalized
  if(! fHistPtVsCentV0AEq_GenXiMinus ) {
    fHistPtVsCentV0AEq_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0AEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0AEq_GenXiMinus); }
  if(! fHistPtVsCentV0AEq_GenXiPlus ) {
    fHistPtVsCentV0AEq_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0AEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0AEq_GenXiPlus); }
  if(! fHistPtVsCentV0AEq_GenOmegaMinus ) {
    fHistPtVsCentV0AEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0AEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0AEq_GenOmegaMinus); }
  if(! fHistPtVsCentV0AEq_GenOmegaPlus ) {
    fHistPtVsCentV0AEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0AEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0AEq_GenOmegaPlus); }  

  if(! fHistPtVsCentV0CEq_GenXiMinus ) {
    fHistPtVsCentV0CEq_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0CEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0CEq_GenXiMinus); }
  if(! fHistPtVsCentV0CEq_GenXiPlus ) {
    fHistPtVsCentV0CEq_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0CEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0CEq_GenXiPlus); }
  if(! fHistPtVsCentV0CEq_GenOmegaMinus ) {
    fHistPtVsCentV0CEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0CEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0CEq_GenOmegaMinus); }
  if(! fHistPtVsCentV0CEq_GenOmegaPlus ) {
    fHistPtVsCentV0CEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0CEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0CEq_GenOmegaPlus); }  

  if(! fHistPtVsCentV0MEq_GenXiMinus ) {
    fHistPtVsCentV0MEq_GenXiMinus    = new TH2D( 
    "fHistPtVsCentV0MEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0MEq_GenXiMinus); }
  if(! fHistPtVsCentV0MEq_GenXiPlus ) {
    fHistPtVsCentV0MEq_GenXiPlus    = new TH2D( 
    "fHistPtVsCentV0MEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0MEq_GenXiPlus); }
  if(! fHistPtVsCentV0MEq_GenOmegaMinus ) {
    fHistPtVsCentV0MEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsCentV0MEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0MEq_GenOmegaMinus); }
  if(! fHistPtVsCentV0MEq_GenOmegaPlus ) {
    fHistPtVsCentV0MEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsCentV0MEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,100,0,100);   
    fListHist->Add(fHistPtVsCentV0MEq_GenOmegaPlus); }  

  //AMPLITUDES: V0A, V0C, V0M, +Eq
  Double_t lMaxAmplitude = 2500; 
  Long_t lAmplitudeBins = 10000;
  if(! fHistPtVsAmpV0A_GenXiMinus ) {
    fHistPtVsAmpV0A_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0A_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0A_GenXiMinus); }
  if(! fHistPtVsAmpV0A_GenXiPlus ) {
    fHistPtVsAmpV0A_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0A_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0A_GenXiPlus); }
  if(! fHistPtVsAmpV0A_GenOmegaMinus ) {
    fHistPtVsAmpV0A_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0A_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0A_GenOmegaMinus); }
  if(! fHistPtVsAmpV0A_GenOmegaPlus ) {
    fHistPtVsAmpV0A_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0A_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0A_GenOmegaPlus); }  

  if(! fHistPtVsAmpV0C_GenXiMinus ) {
    fHistPtVsAmpV0C_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0C_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0C_GenXiMinus); }
  if(! fHistPtVsAmpV0C_GenXiPlus ) {
    fHistPtVsAmpV0C_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0C_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0C_GenXiPlus); }
  if(! fHistPtVsAmpV0C_GenOmegaMinus ) {
    fHistPtVsAmpV0C_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0C_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0C_GenOmegaMinus); }
  if(! fHistPtVsAmpV0C_GenOmegaPlus ) {
    fHistPtVsAmpV0C_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0C_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0C_GenOmegaPlus); }  

  if(! fHistPtVsAmpV0M_GenXiMinus ) {
    fHistPtVsAmpV0M_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0M_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0M_GenXiMinus); }
  if(! fHistPtVsAmpV0M_GenXiPlus ) {
    fHistPtVsAmpV0M_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0M_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0M_GenXiPlus); }
  if(! fHistPtVsAmpV0M_GenOmegaMinus ) {
    fHistPtVsAmpV0M_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0M_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0M_GenOmegaMinus); }
  if(! fHistPtVsAmpV0M_GenOmegaPlus ) {
    fHistPtVsAmpV0M_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0M_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0M_GenOmegaPlus); }  

  //Equalized
  if(! fHistPtVsAmpV0AEq_GenXiMinus ) {
    fHistPtVsAmpV0AEq_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0AEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0AEq_GenXiMinus); }
  if(! fHistPtVsAmpV0AEq_GenXiPlus ) {
    fHistPtVsAmpV0AEq_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0AEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0AEq_GenXiPlus); }
  if(! fHistPtVsAmpV0AEq_GenOmegaMinus ) {
    fHistPtVsAmpV0AEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0AEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0AEq_GenOmegaMinus); }
  if(! fHistPtVsAmpV0AEq_GenOmegaPlus ) {
    fHistPtVsAmpV0AEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0AEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0AEq_GenOmegaPlus); }  

  if(! fHistPtVsAmpV0CEq_GenXiMinus ) {
    fHistPtVsAmpV0CEq_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0CEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0CEq_GenXiMinus); }
  if(! fHistPtVsAmpV0CEq_GenXiPlus ) {
    fHistPtVsAmpV0CEq_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0CEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0CEq_GenXiPlus); }
  if(! fHistPtVsAmpV0CEq_GenOmegaMinus ) {
    fHistPtVsAmpV0CEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0CEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0CEq_GenOmegaMinus); }
  if(! fHistPtVsAmpV0CEq_GenOmegaPlus ) {
    fHistPtVsAmpV0CEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0CEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0CEq_GenOmegaPlus); }  

  if(! fHistPtVsAmpV0MEq_GenXiMinus ) {
    fHistPtVsAmpV0MEq_GenXiMinus    = new TH2D( 
    "fHistPtVsAmpV0MEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0MEq_GenXiMinus); }
  if(! fHistPtVsAmpV0MEq_GenXiPlus ) {
    fHistPtVsAmpV0MEq_GenXiPlus    = new TH2D( 
    "fHistPtVsAmpV0MEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0MEq_GenXiPlus); }
  if(! fHistPtVsAmpV0MEq_GenOmegaMinus ) {
    fHistPtVsAmpV0MEq_GenOmegaMinus    = new TH2D( 
    "fHistPtVsAmpV0MEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0MEq_GenOmegaMinus); }
  if(! fHistPtVsAmpV0MEq_GenOmegaPlus ) {
    fHistPtVsAmpV0MEq_GenOmegaPlus    = new TH2D( 
    "fHistPtVsAmpV0MEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);   
    fListHist->Add(fHistPtVsAmpV0MEq_GenOmegaPlus); }  

   //List of Histograms: Normal
   PostData(1, fListHist);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTreeEvent);
   PostData(3, fTreeV0);
   PostData(4, fTreeCascade);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
   AliMCEvent  *lMCevent  = 0x0; 
   AliStack    *lMCstack  = 0x0; 
  
  // Connect to the InputEvent   
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.

   // Appropriate for ESD analysis! 
      
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
   if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
   }

  //Get VZERO Information for multiplicity later
  AliVVZERO* esdV0 = lESDevent->GetVZEROData();
  if (!esdV0) {
    AliError("AliVVZERO not available");
    return;
  }
        
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;   
    return;
  }

  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }

   fRunNumber = lESDevent->GetRunNumber();

  Double_t lMagneticField = -10; 
  lMagneticField = lESDevent->GetMagneticField( );

//------------------------------------------------
// Variable Definition
//------------------------------------------------

//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  fHistEventCounter->Fill(0.5); 

  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
  
  //Standard Min-Bias Selection
  if ( ! isSelected ) {
    PostData(1, fListHist);
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  fHistEventCounter->Fill(1.5);
 
  //------------------------------------------------
  // Primary Vertex Requirements Section:
  //  ---> pp: has vertex, |z|<10cm
  //------------------------------------------------
  
  //classical Proton-proton like selection 
  const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();	
  const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
  const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

  Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
  lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

  //Only accept if Tracking or SPD vertex is fine 
  if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
    AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  //Has SPD or Tracking Vertex
  fHistEventCounter -> Fill(2.5); 

  //Always do Primary Vertex Selection 
  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  //Fill Event selected counter
  fHistEventCounter -> Fill(3.5);

  //------------------------------------------------
  // Check if this isn't pileup
  //------------------------------------------------

  if(lESDevent->IsPileupFromSPDInMultBins() ){
    // minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  
    //-> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
    AliWarning("Pb / Event tagged as pile-up by SPD... return !"); 
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return; 
  }
  //Fill Event isn't pileup counter
  fHistEventCounter -> Fill(4.5);

//------------------------------------------------
// Multiplicity Information Acquistion
//------------------------------------------------

  //Monte Carlo Level information ! 
  //--------- GENERATED NUMBER OF CHARGED PARTICLES
  // ---> Variable Definition

  Long_t lNchEta5   = 0; 
  Long_t lNchEta8   = 0; 
  Long_t lNchVZEROA = 0; 
  Long_t lNchVZEROC = 0; 

  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++) 
  {// This is the begining of the loop on tracks
      TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
      if(!particleOne) continue;
      if(!particleOne->GetPDG()) continue;
      Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
      if(TMath::Abs(lThisCharge)<0.001) continue;
      if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
     
      //Double_t gpt = particleOne -> Pt();
      Double_t geta = particleOne -> Eta(); 

      if( TMath::Abs(geta) < 0.5 ) lNchEta5++; 
      if( TMath::Abs(geta) < 0.8 ) lNchEta8++; 
      if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++; 
      if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++; 
  }//End of loop on tracks

  //Attribution 
  fTrueMultEta5 = lNchEta5; 
  fTrueMultEta8 = lNchEta8; 
  fTrueMultVZEROA = lNchVZEROA; 
  fTrueMultVZEROC = lNchVZEROC; 
  //----- End Loop on Stack ------------------------------------------------------------

  //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
  fRefMultEta5 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);
  fRefMultEta8 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.8);

  // VZERO PART
  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multV0AEq  = 0;          //  multiplicity from V0 reco side A
  Float_t  multV0CEq  = 0;          //  multiplicity from V0 reco side C
  Float_t  multV0ACorr  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0CCorr  = 0;            //  multiplicity from V0 reco side C

  //Non-Equalized Signal: copy of multV0ACorr and multV0CCorr from AliCentralitySelectionTask
  //Getters for uncorrected multiplicity  
  multV0A=esdV0->GetMTotV0A();
  multV0C=esdV0->GetMTotV0C();

  //Get Z vertex position of SPD vertex (why not Tracking if available?) 
  Float_t zvtx = lPrimarySPDVtx->GetZ(); 

  //Acquire Corrected multV0A 
  multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);    
  multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);   
    
  //Copy to Event Tree for extra information 
  fAmplitude_V0A = multV0ACorr; 
  fAmplitude_V0C = multV0CCorr; 

  // Equalized signals // From AliCentralitySelectionTask
  for(Int_t iCh = 4; iCh < 7; ++iCh) {
    Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
    multV0AEq += mult;
  }
  for(Int_t iCh = 0; iCh < 3; ++iCh) {
    Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
    multV0CEq += mult;
  }
  fAmplitude_V0AEq = multV0AEq; 
  fAmplitude_V0CEq = multV0CEq; 

  fCentrality_V0A   = -100; 
  fCentrality_V0C   = -100; 
  fCentrality_V0M   = -100; 
  fCentrality_V0AEq = -100; 
  fCentrality_V0CEq = -100; 
  fCentrality_V0MEq = -100; 

  //AliCentrality... Check if working? 
  AliCentrality* centrality;
  centrality = lESDevent->GetCentrality();
  if ( !(centrality->GetQuality()>1) ){ 
    fCentrality_V0A   = centrality->GetCentralityPercentile( "V0A"   ); 
    fCentrality_V0C   = centrality->GetCentralityPercentile( "V0C"   ); 
    fCentrality_V0M   = centrality->GetCentralityPercentile( "V0M"   ); 
    fCentrality_V0AEq = centrality->GetCentralityPercentile( "V0AEq" ); 
    fCentrality_V0CEq = centrality->GetCentralityPercentile( "V0CEq" ); 
    fCentrality_V0MEq = centrality->GetCentralityPercentile( "V0MEq" ); 
  }
  
  //Event-level fill 
  fTreeEvent->Fill() ;
  
//------------------------------------------------

//------------------------------------------------
// Fill Efficiency Denominators, please 
//------------------------------------------------

   Int_t    lThisPDG  = 0;
   Double_t lThisRap  = 0;
   Double_t lThisPt   = 0;

//----- Loop on Generated CASCADES ---------------
   for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++) 
   {// This is the begining of the loop on tracks
      
      TParticle* lPart = 0x0; 
      lPart = lMCstack->Particle( ilab );
      if(!lPart){
         Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
         continue;
      }

      lThisPDG = lPart->GetPdgCode();	      

      if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) ) 
      {
        lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
        lThisPt    = lPart->Pt();

        //Use Physical Primaries only for filling These Histos
        if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;

        if( lThisPDG ==  3312 && TMath::Abs(lThisRap) < 0.5 ){
          fHistPt_GenXiMinus                -> Fill (lThisPt);          
          fHistPtVsRefMultEta5_GenXiMinus   -> Fill (lThisPt, fRefMultEta5); 
          fHistPtVsRefMultEta8_GenXiMinus   -> Fill (lThisPt, fRefMultEta8);
          //Centralities 
          fHistPtVsCentV0A_GenXiMinus       -> Fill (lThisPt, fCentrality_V0A);  
          fHistPtVsCentV0C_GenXiMinus       -> Fill (lThisPt, fCentrality_V0C);  
          fHistPtVsCentV0M_GenXiMinus       -> Fill (lThisPt, fCentrality_V0M);  
          fHistPtVsCentV0AEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0AEq);  
          fHistPtVsCentV0CEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0CEq);  
          fHistPtVsCentV0MEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0MEq);  
          //Amplitudes 
          fHistPtVsAmpV0A_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0A);  
          fHistPtVsAmpV0C_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0C);  
          fHistPtVsAmpV0M_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);  
          fHistPtVsAmpV0AEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0AEq);  
          fHistPtVsAmpV0CEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0CEq);  
          fHistPtVsAmpV0MEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);  
        }
        if( lThisPDG == -3312 && TMath::Abs(lThisRap) < 0.5 ){
          fHistPt_GenXiPlus                -> Fill (lThisPt);          
          fHistPtVsRefMultEta5_GenXiPlus   -> Fill (lThisPt, fRefMultEta5); 
          fHistPtVsRefMultEta8_GenXiPlus   -> Fill (lThisPt, fRefMultEta8);
          //Centralities 
          fHistPtVsCentV0A_GenXiPlus       -> Fill (lThisPt, fCentrality_V0A);  
          fHistPtVsCentV0C_GenXiPlus       -> Fill (lThisPt, fCentrality_V0C);  
          fHistPtVsCentV0M_GenXiPlus       -> Fill (lThisPt, fCentrality_V0M);  
          fHistPtVsCentV0AEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0AEq);  
          fHistPtVsCentV0CEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0CEq);  
          fHistPtVsCentV0MEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0MEq);  
          //Amplitudes 
          fHistPtVsAmpV0A_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0A);  
          fHistPtVsAmpV0C_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0C);  
          fHistPtVsAmpV0M_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);  
          fHistPtVsAmpV0AEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0AEq);  
          fHistPtVsAmpV0CEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0CEq);  
          fHistPtVsAmpV0MEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);  
        }
        if( lThisPDG ==  3334 && TMath::Abs(lThisRap) < 0.5 ){
          fHistPt_GenOmegaMinus                -> Fill (lThisPt);          
          fHistPtVsRefMultEta5_GenOmegaMinus   -> Fill (lThisPt, fRefMultEta5); 
          fHistPtVsRefMultEta8_GenOmegaMinus   -> Fill (lThisPt, fRefMultEta8);
          //Centralities 
          fHistPtVsCentV0A_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0A);  
          fHistPtVsCentV0C_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0C);  
          fHistPtVsCentV0M_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0M);  
          fHistPtVsCentV0AEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0AEq);  
          fHistPtVsCentV0CEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0CEq);  
          fHistPtVsCentV0MEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0MEq);  
          //Amplitudes 
          fHistPtVsAmpV0A_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0A);  
          fHistPtVsAmpV0C_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0C);  
          fHistPtVsAmpV0M_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);  
          fHistPtVsAmpV0AEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0AEq);  
          fHistPtVsAmpV0CEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0CEq);  
          fHistPtVsAmpV0MEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);  
        }
        if( lThisPDG == -3334 && TMath::Abs(lThisRap) < 0.5 ){
          fHistPt_GenOmegaPlus                -> Fill (lThisPt);          
          fHistPtVsRefMultEta5_GenOmegaPlus   -> Fill (lThisPt, fRefMultEta5); 
          fHistPtVsRefMultEta8_GenOmegaPlus   -> Fill (lThisPt, fRefMultEta8);
          //Centralities 
          fHistPtVsCentV0A_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0A);  
          fHistPtVsCentV0C_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0C);  
          fHistPtVsCentV0M_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0M);  
          fHistPtVsCentV0AEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0AEq);  
          fHistPtVsCentV0CEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0CEq);  
          fHistPtVsCentV0MEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0MEq);  
          //Amplitudes 
          fHistPtVsAmpV0A_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0A);  
          fHistPtVsAmpV0C_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0C);  
          fHistPtVsAmpV0M_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);  
          fHistPtVsAmpV0AEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0AEq);  
          fHistPtVsAmpV0CEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0CEq);  
          fHistPtVsAmpV0MEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);  
        }
      }
   }//End of loop on tracks
//----- End Loop on Cascades ------------------------------------------------------------

//------------------------------------------------
// Rerun cascade vertexer! 
//------------------------------------------------
    
  if( fkRunVertexers ){ 
    lESDevent->ResetCascades();
    lESDevent->ResetV0s();

    AliV0vertexer lV0vtxer;
    AliCascadeVertexer lCascVtxer;
                  
    lV0vtxer.SetDefaultCuts(fV0VertexerSels);
    lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);

    lV0vtxer.Tracks2V0vertices(lESDevent);
    lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
  }

//------------------------------------------------
// MAIN CASCADE LOOP STARTS HERE
//------------------------------------------------
// Code Credit: Antonin Maire (thanks^100)
// ---> This is an adaptation

  Long_t ncascades = 0;
	ncascades = lESDevent->GetNumberOfCascades();
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++){
    //------------------------------------------------
    // Initializations
    //------------------------------------------------	
	  //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
	  //Double_t lBestPrimaryVtxRadius3D = -500.0;

	  // - 1st part of initialisation : variables needed to store AliESDCascade data members
	  Double_t lEffMassXi      = 0. ;
	  //Double_t lChi2Xi         = -1. ;
	  Double_t lDcaXiDaughters = -1. ;
	  Double_t lXiCosineOfPointingAngle = -1. ;
	  Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
	  Double_t lXiRadius = -1000. ;
          
	  // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
	  Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
	  Int_t    lNegTPCClusters    = -1; // For ESD only ...
	  Int_t    lBachTPCClusters   = -1; // For ESD only ...
          		
          // - 3rd part of initialisation : about V0 part in cascades
	  Double_t lInvMassLambdaAsCascDghter = 0.;
	  //Double_t lV0Chi2Xi         = -1. ;
	  Double_t lDcaV0DaughtersXi = -1.;
		
	  Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
	  Double_t lDcaPosToPrimVertexXi  = -1.;
	  Double_t lDcaNegToPrimVertexXi  = -1.;
	  Double_t lV0CosineOfPointingAngleXi = -1. ;
	  Double_t lV0CosineOfPointingAngleXiSpecial = -1. ;
	  Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
	  Double_t lV0RadiusXi = -1000.0;
	  Double_t lV0quality  = 0.;
	
	  // - 4th part of initialisation : Effective masses
	  Double_t lInvMassXiMinus    = 0.;
	  Double_t lInvMassXiPlus     = 0.;
	  Double_t lInvMassOmegaMinus = 0.;
	  Double_t lInvMassOmegaPlus  = 0.;
    
	  // - 6th part of initialisation : extra info for QA
	  Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
	  Double_t lXiTransvMom  = 0. ;
	  //Double_t lXiTransvMomMC= 0. ;
	  Double_t lXiTotMom     = 0. ;
		
	  Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
	  //Double_t lBachTransvMom  = 0.;
	  //Double_t lBachTotMom     = 0.;

    fTreeCascVarNegNSigmaPion   = -100;
    fTreeCascVarNegNSigmaProton = -100;
    fTreeCascVarPosNSigmaPion   = -100;
    fTreeCascVarPosNSigmaProton = -100;
    fTreeCascVarBachNSigmaPion  = -100;
    fTreeCascVarBachNSigmaKaon  = -100;
	
	  Short_t  lChargeXi = -2;
	  //Double_t lV0toXiCosineOfPointingAngle = 0. ;
	
	  Double_t lRapXi   = -20.0, lRapOmega = -20.0, lRapMC = -20;//  lEta = -20.0, lTheta = 360., lPhi = 720. ;
	  //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	    
    // -------------------------------------
    // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
    
	  AliESDcascade *xi = lESDevent->GetCascade(iXi);
	  if (!xi) continue;
	
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)	
		//-------------
	  lV0quality = 0.;
	  xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

	  lEffMassXi  			= xi->GetEffMassXi();
	  //lChi2Xi 			    = xi->GetChi2Xi();
	  lDcaXiDaughters 	= xi->GetDcaXiDaughters();
	  lXiCosineOfPointingAngle 	            = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                                 lBestPrimaryVtxPos[1],
                                                                                 lBestPrimaryVtxPos[2] );
		  // Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
	  xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
	  lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );		

		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
        UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                // Care track label can be negative in MC production (linked with the track quality)
                // However = normally, not the case for track index ...
          
	  // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
	  if(lBachIdx == lIdxNegXi) {
		  AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
	  }
    if(lBachIdx == lIdxPosXi) {
    	AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
	  }
          
	  AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	  AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	  AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );

	  if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
		  AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
		  continue;
	  }

      fTreeCascVarPosEta = pTrackXi->Eta();
      fTreeCascVarNegEta = nTrackXi->Eta();
      fTreeCascVarBachEta = bachTrackXi->Eta();
      
      Double_t lBMom[3], lNMom[3], lPMom[3];
      xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
      xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
      xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
      
      //fTreeCascVarBachTransMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] );
      //fTreeCascVarPosTransMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
      //fTreeCascVarNegTransMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
  
    //------------------------------------------------
    // TPC dEdx information 
    //------------------------------------------------
    fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion   );
    fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
    fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
    fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
    fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
    fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );

    //------------------------------------------------
    // TPC Number of clusters info
    // --- modified to save the smallest number 
    // --- of TPC clusters for the 3 tracks
    //------------------------------------------------
              
	  lPosTPCClusters   = pTrackXi->GetTPCNcls();
	  lNegTPCClusters   = nTrackXi->GetTPCNcls();
	  lBachTPCClusters  = bachTrackXi->GetTPCNcls(); 

    // 1 - Poor quality related to TPCrefit
	  ULong_t pStatus    = pTrackXi->GetStatus();
	  ULong_t nStatus    = nTrackXi->GetStatus();
	  ULong_t bachStatus = bachTrackXi->GetStatus();

    //fTreeCascVarkITSRefitBachelor = kTRUE; 
    //fTreeCascVarkITSRefitNegative = kTRUE; 
    //fTreeCascVarkITSRefitPositive = kTRUE; 

    if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
    if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
    if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }

    //Extra Debug Information: booleans for ITS refit
    //if ((pStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitPositive = kFALSE; }
    //if ((nStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitNegative = kFALSE; }
    //if ((bachStatus&AliESDtrack::kITSrefit) == 0) { fTreeCascVarkITSRefitBachelor = kFALSE; }

	  // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
    if(lPosTPCClusters  < 70) { AliWarning("Pb / V0 Pos. track has less than 70 TPC clusters ... continue!"); continue; }
	  if(lNegTPCClusters  < 70) { AliWarning("Pb / V0 Neg. track has less than 70 TPC clusters ... continue!"); continue; }
	  if(lBachTPCClusters < 70) { AliWarning("Pb / Bach.   track has less than 70 TPC clusters ... continue!"); continue; }
	  Int_t leastnumberofclusters = 1000; 
	  if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
	  if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
	  if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;

	  lInvMassLambdaAsCascDghter	= xi->GetEffMass();
	  // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
	  lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters(); 
	  //lV0Chi2Xi 			= xi->GetChi2V0();
	
	  lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
									    lBestPrimaryVtxPos[1],
									    lBestPrimaryVtxPos[2] );
    //Modification: V0 CosPA wrt to Cascade decay vertex
	  lV0CosineOfPointingAngleXiSpecial 	= xi->GetV0CosineOfPointingAngle( lPosXi[0],
									    lPosXi[1],
									    lPosXi[2] );

	  lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0], 
						      lBestPrimaryVtxPos[1], 
						      lBestPrimaryVtxPos[2] );
		
	  lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0], 
						       		lBestPrimaryVtxPos[1], 
						       		lMagneticField  ) ); 
					  // Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
	  xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
	  lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
	
	  lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
						   		lBestPrimaryVtxPos[1], 
						   		lMagneticField  )     ); 
	
	  lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
					      			lBestPrimaryVtxPos[1], 
					      			lMagneticField  )     ); 
		
	  // - II.Step 4 : around effective masses (ESD)
	  // ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		
	  if( bachTrackXi->Charge() < 0 )	{
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , 3312); 	
			  // Calculate the effective mass of the Xi- candidate. 
			  // pdg code 3312 = Xi-
		  lInvMassXiMinus = xi->GetEffMassXi();
		
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , 3334); 	
			  // Calculate the effective mass of the Xi- candidate. 
			  // pdg code 3334 = Omega-
		  lInvMassOmegaMinus = xi->GetEffMassXi();
					
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
	  }// end if negative bachelor
	
	
	  if( bachTrackXi->Charge() >  0 ){
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , -3312); 	
			  // Calculate the effective mass of the Xi+ candidate. 
			  // pdg code -3312 = Xi+
		  lInvMassXiPlus = xi->GetEffMassXi();
		
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , -3334); 	
			  // Calculate the effective mass of the Xi+ candidate. 
			  // pdg code -3334  = Omega+
		  lInvMassOmegaPlus = xi->GetEffMassXi();
		
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
	  }// end if positive bachelor
		  // - II.Step 6 : extra info for QA (ESD)
		  // miscellaneous pieces of info that may help regarding data quality assessment.
		  //-------------

	  xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
		  lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
		  lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
	  xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
		  //lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
		  //lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	  lChargeXi = xi->Charge();

	  //lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
	
	  lRapXi    = xi->RapXi();
	  lRapOmega = xi->RapOmega();
	  //lEta      = xi->Eta();
	  //lTheta    = xi->Theta() *180.0/TMath::Pi();
	  //lPhi      = xi->Phi()   *180.0/TMath::Pi();
	  //lAlphaXi  = xi->AlphaXi();
	  //lPtArmXi  = xi->PtArmXi();

//------------------------------------------------
// Associate Cascade Candidates to Monte Carlo!
//------------------------------------------------	

//Warning: Not using Continues... Need to fill tree later!

  Double_t lXiTransvMomMC= 0. ;	
	Int_t lPDGCodeCascade = 0;	
	Int_t lPID_BachMother = 0;
	Int_t lPID_NegMother = 0;
	Int_t lPID_PosMother = 0;
  fTreeCascVarIsPhysicalPrimary = 0; // 0: not defined, any candidate may have this

	if(fDebug > 5)
		cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent() 
			<< " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
	
//----------------------------------------
// Regular MC ASSOCIATION STARTS HERE
//----------------------------------------

	  Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		  // Abs value = needed ! question of quality track association ...
	  Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
	  Int_t lblBach        = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );

	  TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	  TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	  TParticle* mcBach        = lMCstack->Particle( lblBach );
      
    //fTreeCascVarPosTransMomMC = mcPosV0Dghter->Pt();
    //fTreeCascVarNegTransMomMC = mcNegV0Dghter->Pt();

	  //fTreeCascVarPIDPositive = mcPosV0Dghter -> GetPdgCode();
	  //fTreeCascVarPIDNegative = mcNegV0Dghter -> GetPdgCode();
	  //fTreeCascVarPIDBachelor = mcBach->GetPdgCode();

	  // - Step 4.2 : level of the Xi daughters
		
	  Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
	  Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
	
	  //Rather uncivilized: Open brackets for each 'continue'
	  if(! (lblMotherPosV0Dghter != lblMotherNegV0Dghter) ) { // same mother
	  if(! (lblMotherPosV0Dghter < 0) ) { // mother != primary (!= -1)
	  if(! (lblMotherNegV0Dghter < 0) ) {
					
		// mothers = Lambda candidate ... a priori
	
	  TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
	  TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );
			
	  // - Step 4.3 : level of Xi candidate
	
	  Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetFirstMother() ;
	  Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetFirstMother() ;
				
		if(! (lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter) ) {
		if(! (lblGdMotherPosV0Dghter < 0) ) { // primary lambda ...
		if(! (lblGdMotherNegV0Dghter < 0) ) { // primary lambda ...

		  // Gd mothers = Xi candidate ... a priori
	
	  TParticle* mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
	  TParticle* mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );
					
	  Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother()  );
	
  //		if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
		  if(!(lblMotherBach != lblGdMotherPosV0Dghter)) { //same mother for bach and V0 daughters
	
	  TParticle* mcMotherBach = lMCstack->Particle( lblMotherBach );
	
    lPID_BachMother = mcMotherBach->GetPdgCode();
	  lPID_NegMother = mcGdMotherPosV0Dghter->GetPdgCode();
	  lPID_PosMother = mcGdMotherNegV0Dghter->GetPdgCode();
   
	  if(lPID_BachMother==lPID_NegMother && lPID_BachMother==lPID_PosMother){ 
		  lPDGCodeCascade = lPID_BachMother; 
      lXiTransvMomMC = mcMotherBach->Pt();
      if( lMCstack->IsPhysicalPrimary       (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 1; //Is Primary!
      if( lMCstack->IsSecondaryFromWeakDecay(lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 2; //Weak Decay!
      if( lMCstack->IsSecondaryFromMaterial (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 3; //From Material!
      if ( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) !=0 ){
        lRapMC = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
      }
	  }

  }}}}}}} //Ends all conditionals above...

//----------------------------------------
// Regular MC ASSOCIATION ENDS HERE
//----------------------------------------


  //------------------------------------------------
  // Set Variables for adding to tree
  //------------------------------------------------		
	
          fTreeCascVarCharge	= lChargeXi;
          fTreeCascVarPID = lPDGCodeCascade; 
          if(lInvMassXiMinus!=0)    fTreeCascVarMassAsXi = lInvMassXiMinus;
          if(lInvMassXiPlus!=0)     fTreeCascVarMassAsXi = lInvMassXiPlus;
          if(lInvMassOmegaMinus!=0) fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
          if(lInvMassOmegaPlus!=0)  fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
          fTreeCascVarPt = lXiTransvMom;
          fTreeCascVarPtMC = lXiTransvMomMC;
          fTreeCascVarRapXi = lRapXi ;
          fTreeCascVarRapMC = lRapMC ;
          fTreeCascVarRapOmega = lRapOmega ;
          fTreeCascVarDCACascDaughters = lDcaXiDaughters;
          fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
          fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
          fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
          fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
          fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
          fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle;
          fTreeCascVarCascRadius = lXiRadius;
          fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
          fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi;
          fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
          fTreeCascVarV0Radius = lV0RadiusXi;
          fTreeCascVarLeastNbrClusters = leastnumberofclusters;

          //Copy Multiplicity information 
          fTreeCascVarCentV0A = fCentrality_V0A; 
          fTreeCascVarCentV0C = fCentrality_V0C; 
          fTreeCascVarCentV0M = fCentrality_V0M; 
          fTreeCascVarCentV0AEq = fCentrality_V0AEq; 
          fTreeCascVarCentV0CEq = fCentrality_V0CEq; 
          fTreeCascVarCentV0MEq = fCentrality_V0MEq; 
          fTreeCascVarAmpV0A = fAmplitude_V0A; 
          fTreeCascVarAmpV0C = fAmplitude_V0C; 
          fTreeCascVarAmpV0AEq = fAmplitude_V0AEq; 
          fTreeCascVarAmpV0CEq = fAmplitude_V0CEq; 
          fTreeCascVarRefMultEta8 = fRefMultEta8;
          fTreeCascVarRefMultEta5 = fRefMultEta5;
          fTreeCascVarRunNumber = fRunNumber; 

          fTreeCascVarDistOverTotMom = TMath::Sqrt(
						TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
					);
          fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);

//All vars not specified here: specified elsewhere!

//------------------------------------------------
// Fill Tree! 
//------------------------------------------------

// The conditional is meant to decrease excessive
// memory usage! Be careful when loosening the 
// cut!

  //Xi    Mass window: 150MeV wide
  //Omega mass window: 150MeV wide

  if( fkSaveCascadeTree && ( (fTreeCascVarMassAsXi<1.32+0.075&&fTreeCascVarMassAsXi>1.32-0.075) ||
      (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075) ) ){
      fTreeCascade->Fill();
  }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

	}// end of the Cascade loop (ESD or AOD)

  // Post output data.
  PostData(1, fListHist); 
  PostData(2, fTreeEvent);
  PostData(3, fTreeV0);
  PostData(4, fTreeCascade);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMC::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMC : ouput data container list not available\n");
      return;
   }	
	
   fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
   if (!fHistEventCounter) {
      Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMC : fHistEventCounter not available");
      return;
   }
  
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityMC","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();

   fHistEventCounter->SetMarkerStyle(22);
   fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskStrangenessVsMultiplicityMC::MyRapidity(Double_t rE, Double_t rPz) const
{
   // Local calculation for rapidity
   Double_t ReturnValue = -100;
   if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
      ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
   }
   return ReturnValue;
} 
