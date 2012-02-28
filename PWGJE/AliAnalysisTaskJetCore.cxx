
// ******************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// the distance from track to jet axis and a 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
// *******************************************


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


#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliAnalysisTaskFastEmbedding.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetCore.h"

ClassImp(AliAnalysisTaskJetCore)

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore() :
AliAnalysisTaskSE(),
fESD(0x0),
fAOD(0x0),
fAODExtension(0x0),
fBackgroundBranch(""),
fNonStdFile(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.),
fVtxZMax(10.),
fEvtClassMin(0),
fEvtClassMax(4),
fFilterMask(0),
fRadioFrac(0.2),
fMinDist(0.1),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fAngStructCloseTracks(0),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod3C10(0x0),
fh2JetCoreMethod1C20(0x0),
fh2JetCoreMethod2C20(0x0),
fh2JetCoreMethod3C20(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod3C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2JetCoreMethod3C60(0x0),
fh2JetCoreMethod3C10lead(0x0),
fh2JetCoreMethod3C20lead(0x0),
fh2JetCoreMethod3C30lead(0x0),
fh2JetCoreMethod3C60lead(0x0),
fh2JetCoreMethod3C10sublead(0x0),
fh2JetCoreMethod3C20sublead(0x0),
fh2JetCoreMethod3C30sublead(0x0),
fh2JetCoreMethod3C60sublead(0x0),
fh2SumPtInC10(0x0),
fh2SumPtInC20(0x0),
fh2SumPtInC30(0x0),
fh2SumPtInC60(0x0),
fh2SumPtInC10lead(0x0),
fh2SumPtInC20lead(0x0),
fh2SumPtInC30lead(0x0),
fh2SumPtInC60lead(0x0),
fh2SumPtInC10sublead(0x0),
fh2SumPtInC20sublead(0x0),
fh2SumPtInC30sublead(0x0),
fh2SumPtInC60sublead(0x0),
fh2SumPtOutC10(0x0),
fh2SumPtOutC20(0x0),
fh2SumPtOutC30(0x0),
fh2SumPtOutC60(0x0),
fh2SumPtOutC10lead(0x0),
fh2SumPtOutC20lead(0x0),
fh2SumPtOutC30lead(0x0),
fh2SumPtOutC60lead(0x0),
fh2SumPtOutC10sublead(0x0),
fh2SumPtOutC20sublead(0x0),
fh2SumPtOutC30sublead(0x0),
fh2SumPtOutC60sublead(0x0),
fh2SumPtInC10bkg(0x0),
fh2SumPtInC20bkg(0x0),
fh2SumPtInC30bkg(0x0),
fh2SumPtInC60bkg(0x0),
fh2SumPtInC10bkglead(0x0),
fh2SumPtInC20bkglead(0x0),
fh2SumPtInC30bkglead(0x0),
fh2SumPtInC60bkglead(0x0),
fh2SumPtInC10bkgsublead(0x0),
fh2SumPtInC20bkgsublead(0x0),
fh2SumPtInC30bkgsublead(0x0),
fh2SumPtInC60bkgsublead(0x0),
fh2SumPtOutC10bkg(0x0),
fh2SumPtOutC20bkg(0x0),
fh2SumPtOutC30bkg(0x0),
fh2SumPtOutC60bkg(0x0),
fh2SumPtOutC10bkglead(0x0),
fh2SumPtOutC20bkglead(0x0),
fh2SumPtOutC30bkglead(0x0),
fh2SumPtOutC60bkglead(0x0),
fh2SumPtOutC10bkgsublead(0x0),
fh2SumPtOutC20bkgsublead(0x0),
fh2SumPtOutC30bkgsublead(0x0),
fh2SumPtOutC60bkgsublead(0x0),
fh2DeltaRC10pt1(0x0),
fh2DeltaRC20pt1(0x0),
fh2DeltaRC30pt1(0x0),
fh2DeltaRC60pt1(0x0),  
fh2DeltaRC10pt2(0x0),
fh2DeltaRC20pt2(0x0),
fh2DeltaRC30pt2(0x0),
fh2DeltaRC60pt2(0x0),  
fh2DeltaRC10pt3(0x0),
fh2DeltaRC20pt3(0x0),
fh2DeltaRC30pt3(0x0),
fh2DeltaRC60pt3(0x0),  
fh2DeltaRC10pt4(0x0),
fh2DeltaRC20pt4(0x0),
fh2DeltaRC30pt4(0x0),
fh2DeltaRC60pt4(0x0),
fh2DeltaEtaC10pt1(0x0),
fh2DeltaEtaC20pt1(0x0),
fh2DeltaEtaC30pt1(0x0),
fh2DeltaEtaC60pt1(0x0),  
fh2DeltaEtaC10pt2(0x0),
fh2DeltaEtaC20pt2(0x0),
fh2DeltaEtaC30pt2(0x0),
fh2DeltaEtaC60pt2(0x0),  
fh2DeltaEtaC10pt3(0x0),
fh2DeltaEtaC20pt3(0x0),
fh2DeltaEtaC30pt3(0x0),
fh2DeltaEtaC60pt3(0x0),  
fh2DeltaEtaC10pt4(0x0),
fh2DeltaEtaC20pt4(0x0),
fh2DeltaEtaC30pt4(0x0),
fh2DeltaEtaC60pt4(0x0),
fh2DeltaPhiC10pt1(0x0),
fh2DeltaPhiC20pt1(0x0),
fh2DeltaPhiC30pt1(0x0),
fh2DeltaPhiC60pt1(0x0),  
fh2DeltaPhiC10pt2(0x0),
fh2DeltaPhiC20pt2(0x0),
fh2DeltaPhiC30pt2(0x0),
fh2DeltaPhiC60pt2(0x0),  
fh2DeltaPhiC10pt3(0x0),
fh2DeltaPhiC20pt3(0x0),
fh2DeltaPhiC30pt3(0x0),
fh2DeltaPhiC60pt3(0x0),  
fh2DeltaPhiC10pt4(0x0),
fh2DeltaPhiC20pt4(0x0),
fh2DeltaPhiC30pt4(0x0),
fh2DeltaPhiC60pt4(0x0),
fh2DeltaRC10pt1lead(0x0),
fh2DeltaRC20pt1lead(0x0),
fh2DeltaRC30pt1lead(0x0),
fh2DeltaRC60pt1lead(0x0),  
fh2DeltaRC10pt2lead(0x0),
fh2DeltaRC20pt2lead(0x0),
fh2DeltaRC30pt2lead(0x0),
fh2DeltaRC60pt2lead(0x0),  
fh2DeltaRC10pt3lead(0x0),
fh2DeltaRC20pt3lead(0x0),
fh2DeltaRC30pt3lead(0x0),
fh2DeltaRC60pt3lead(0x0),  
fh2DeltaRC10pt4lead(0x0),
fh2DeltaRC20pt4lead(0x0),
fh2DeltaRC30pt4lead(0x0),
fh2DeltaRC60pt4lead(0x0),
fh2DeltaEtaC10pt1lead(0x0),
fh2DeltaEtaC20pt1lead(0x0),
fh2DeltaEtaC30pt1lead(0x0),
fh2DeltaEtaC60pt1lead(0x0),  
fh2DeltaEtaC10pt2lead(0x0),
fh2DeltaEtaC20pt2lead(0x0),
fh2DeltaEtaC30pt2lead(0x0),
fh2DeltaEtaC60pt2lead(0x0),  
fh2DeltaEtaC10pt3lead(0x0),
fh2DeltaEtaC20pt3lead(0x0),
fh2DeltaEtaC30pt3lead(0x0),
fh2DeltaEtaC60pt3lead(0x0),  
fh2DeltaEtaC10pt4lead(0x0),
fh2DeltaEtaC20pt4lead(0x0),
fh2DeltaEtaC30pt4lead(0x0),
fh2DeltaEtaC60pt4lead(0x0),
fh2DeltaPhiC10pt1lead(0x0),
fh2DeltaPhiC20pt1lead(0x0),
fh2DeltaPhiC30pt1lead(0x0),
fh2DeltaPhiC60pt1lead(0x0),  
fh2DeltaPhiC10pt2lead(0x0),
fh2DeltaPhiC20pt2lead(0x0),
fh2DeltaPhiC30pt2lead(0x0),
fh2DeltaPhiC60pt2lead(0x0),  
fh2DeltaPhiC10pt3lead(0x0),
fh2DeltaPhiC20pt3lead(0x0),
fh2DeltaPhiC30pt3lead(0x0),
fh2DeltaPhiC60pt3lead(0x0),  
fh2DeltaPhiC10pt4lead(0x0),
fh2DeltaPhiC20pt4lead(0x0),
fh2DeltaPhiC30pt4lead(0x0),
fh2DeltaPhiC60pt4lead(0x0),
fh2DeltaRC10pt1sublead(0x0),
fh2DeltaRC20pt1sublead(0x0),
fh2DeltaRC30pt1sublead(0x0),
fh2DeltaRC60pt1sublead(0x0),  
fh2DeltaRC10pt2sublead(0x0),
fh2DeltaRC20pt2sublead(0x0),
fh2DeltaRC30pt2sublead(0x0),
fh2DeltaRC60pt2sublead(0x0),  
fh2DeltaRC10pt3sublead(0x0),
fh2DeltaRC20pt3sublead(0x0),
fh2DeltaRC30pt3sublead(0x0),
fh2DeltaRC60pt3sublead(0x0),  
fh2DeltaRC10pt4sublead(0x0),
fh2DeltaRC20pt4sublead(0x0),
fh2DeltaRC30pt4sublead(0x0),
fh2DeltaRC60pt4sublead(0x0),
fh2DeltaEtaC10pt1sublead(0x0),
fh2DeltaEtaC20pt1sublead(0x0),
fh2DeltaEtaC30pt1sublead(0x0),
fh2DeltaEtaC60pt1sublead(0x0),  
fh2DeltaEtaC10pt2sublead(0x0),
fh2DeltaEtaC20pt2sublead(0x0),
fh2DeltaEtaC30pt2sublead(0x0),
fh2DeltaEtaC60pt2sublead(0x0),  
fh2DeltaEtaC10pt3sublead(0x0),
fh2DeltaEtaC20pt3sublead(0x0),
fh2DeltaEtaC30pt3sublead(0x0),
fh2DeltaEtaC60pt3sublead(0x0),  
fh2DeltaEtaC10pt4sublead(0x0),
fh2DeltaEtaC20pt4sublead(0x0),
fh2DeltaEtaC30pt4sublead(0x0),
fh2DeltaEtaC60pt4sublead(0x0),
fh2DeltaPhiC10pt1sublead(0x0),
fh2DeltaPhiC20pt1sublead(0x0),
fh2DeltaPhiC30pt1sublead(0x0),
fh2DeltaPhiC60pt1sublead(0x0),  
fh2DeltaPhiC10pt2sublead(0x0),
fh2DeltaPhiC20pt2sublead(0x0),
fh2DeltaPhiC30pt2sublead(0x0),
fh2DeltaPhiC60pt2sublead(0x0),  
fh2DeltaPhiC10pt3sublead(0x0),
fh2DeltaPhiC20pt3sublead(0x0),
fh2DeltaPhiC30pt3sublead(0x0),
fh2DeltaPhiC60pt3sublead(0x0),  
fh2DeltaPhiC10pt4sublead(0x0),
fh2DeltaPhiC20pt4sublead(0x0),
fh2DeltaPhiC30pt4sublead(0x0),
fh2DeltaPhiC60pt4sublead(0x0),
fh2AngStructpt1C10(0x0),
fh2AngStructpt2C10(0x0),
fh2AngStructpt3C10(0x0),
fh2AngStructpt4C10(0x0),
fh2AngStructpt1C20(0x0),
fh2AngStructpt2C20(0x0),
fh2AngStructpt3C20(0x0),
fh2AngStructpt4C20(0x0),    
fh2AngStructpt1C30(0x0),
fh2AngStructpt2C30(0x0),
fh2AngStructpt3C30(0x0),
fh2AngStructpt4C30(0x0),   
fh2AngStructpt1C60(0x0),
fh2AngStructpt2C60(0x0),
fh2AngStructpt3C60(0x0),
fh2AngStructpt4C60(0x0)
{
   // default Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAOD(0x0),
fAODExtension(0x0),
fBackgroundBranch(""),
fNonStdFile(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.),
fVtxZMax(10.),
fEvtClassMin(0),
fEvtClassMax(4),
fFilterMask(0),
fRadioFrac(0.2),
fMinDist(0.1),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fAngStructCloseTracks(0),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod3C10(0x0),
fh2JetCoreMethod1C20(0x0),
fh2JetCoreMethod2C20(0x0),
fh2JetCoreMethod3C20(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod3C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2JetCoreMethod3C60(0x0),
fh2JetCoreMethod3C10lead(0x0),
fh2JetCoreMethod3C20lead(0x0),
fh2JetCoreMethod3C30lead(0x0),
fh2JetCoreMethod3C60lead(0x0),
fh2JetCoreMethod3C10sublead(0x0),
fh2JetCoreMethod3C20sublead(0x0),
fh2JetCoreMethod3C30sublead(0x0),
fh2JetCoreMethod3C60sublead(0x0),
fh2SumPtInC10(0x0),
fh2SumPtInC20(0x0),
fh2SumPtInC30(0x0),
fh2SumPtInC60(0x0),
fh2SumPtInC10lead(0x0),
fh2SumPtInC20lead(0x0),
fh2SumPtInC30lead(0x0),
fh2SumPtInC60lead(0x0),
fh2SumPtInC10sublead(0x0),
fh2SumPtInC20sublead(0x0),
fh2SumPtInC30sublead(0x0),
fh2SumPtInC60sublead(0x0),
fh2SumPtOutC10(0x0),
fh2SumPtOutC20(0x0),
fh2SumPtOutC30(0x0),
fh2SumPtOutC60(0x0),
fh2SumPtOutC10lead(0x0),
fh2SumPtOutC20lead(0x0),
fh2SumPtOutC30lead(0x0),
fh2SumPtOutC60lead(0x0),
fh2SumPtOutC10sublead(0x0),
fh2SumPtOutC20sublead(0x0),
fh2SumPtOutC30sublead(0x0),
fh2SumPtOutC60sublead(0x0),
fh2SumPtInC10bkg(0x0),
fh2SumPtInC20bkg(0x0),
fh2SumPtInC30bkg(0x0),
fh2SumPtInC60bkg(0x0),
fh2SumPtInC10bkglead(0x0),
fh2SumPtInC20bkglead(0x0),
fh2SumPtInC30bkglead(0x0),
fh2SumPtInC60bkglead(0x0),
fh2SumPtInC10bkgsublead(0x0),
fh2SumPtInC20bkgsublead(0x0),
fh2SumPtInC30bkgsublead(0x0),
fh2SumPtInC60bkgsublead(0x0),
fh2SumPtOutC10bkg(0x0),
fh2SumPtOutC20bkg(0x0),
fh2SumPtOutC30bkg(0x0),
fh2SumPtOutC60bkg(0x0),
fh2SumPtOutC10bkglead(0x0),
fh2SumPtOutC20bkglead(0x0),
fh2SumPtOutC30bkglead(0x0),
fh2SumPtOutC60bkglead(0x0),
fh2SumPtOutC10bkgsublead(0x0),
fh2SumPtOutC20bkgsublead(0x0),
fh2SumPtOutC30bkgsublead(0x0),
fh2SumPtOutC60bkgsublead(0x0),
fh2DeltaRC10pt1(0x0),
fh2DeltaRC20pt1(0x0),
fh2DeltaRC30pt1(0x0),
fh2DeltaRC60pt1(0x0),  
fh2DeltaRC10pt2(0x0),
fh2DeltaRC20pt2(0x0),
fh2DeltaRC30pt2(0x0),
fh2DeltaRC60pt2(0x0),  
fh2DeltaRC10pt3(0x0),
fh2DeltaRC20pt3(0x0),
fh2DeltaRC30pt3(0x0),
fh2DeltaRC60pt3(0x0),  
fh2DeltaRC10pt4(0x0),
fh2DeltaRC20pt4(0x0),
fh2DeltaRC30pt4(0x0),
fh2DeltaRC60pt4(0x0),
fh2DeltaEtaC10pt1(0x0),
fh2DeltaEtaC20pt1(0x0),
fh2DeltaEtaC30pt1(0x0),
fh2DeltaEtaC60pt1(0x0),  
fh2DeltaEtaC10pt2(0x0),
fh2DeltaEtaC20pt2(0x0),
fh2DeltaEtaC30pt2(0x0),
fh2DeltaEtaC60pt2(0x0),  
fh2DeltaEtaC10pt3(0x0),
fh2DeltaEtaC20pt3(0x0),
fh2DeltaEtaC30pt3(0x0),
fh2DeltaEtaC60pt3(0x0),  
fh2DeltaEtaC10pt4(0x0),
fh2DeltaEtaC20pt4(0x0),
fh2DeltaEtaC30pt4(0x0),
fh2DeltaEtaC60pt4(0x0),
fh2DeltaPhiC10pt1(0x0),
fh2DeltaPhiC20pt1(0x0),
fh2DeltaPhiC30pt1(0x0),
fh2DeltaPhiC60pt1(0x0),  
fh2DeltaPhiC10pt2(0x0),
fh2DeltaPhiC20pt2(0x0),
fh2DeltaPhiC30pt2(0x0),
fh2DeltaPhiC60pt2(0x0),  
fh2DeltaPhiC10pt3(0x0),
fh2DeltaPhiC20pt3(0x0),
fh2DeltaPhiC30pt3(0x0),
fh2DeltaPhiC60pt3(0x0),  
fh2DeltaPhiC10pt4(0x0),
fh2DeltaPhiC20pt4(0x0),
fh2DeltaPhiC30pt4(0x0),
fh2DeltaPhiC60pt4(0x0),
fh2DeltaRC10pt1lead(0x0),
fh2DeltaRC20pt1lead(0x0),
fh2DeltaRC30pt1lead(0x0),
fh2DeltaRC60pt1lead(0x0),  
fh2DeltaRC10pt2lead(0x0),
fh2DeltaRC20pt2lead(0x0),
fh2DeltaRC30pt2lead(0x0),
fh2DeltaRC60pt2lead(0x0),  
fh2DeltaRC10pt3lead(0x0),
fh2DeltaRC20pt3lead(0x0),
fh2DeltaRC30pt3lead(0x0),
fh2DeltaRC60pt3lead(0x0),  
fh2DeltaRC10pt4lead(0x0),
fh2DeltaRC20pt4lead(0x0),
fh2DeltaRC30pt4lead(0x0),
fh2DeltaRC60pt4lead(0x0),
fh2DeltaEtaC10pt1lead(0x0),
fh2DeltaEtaC20pt1lead(0x0),
fh2DeltaEtaC30pt1lead(0x0),
fh2DeltaEtaC60pt1lead(0x0),  
fh2DeltaEtaC10pt2lead(0x0),
fh2DeltaEtaC20pt2lead(0x0),
fh2DeltaEtaC30pt2lead(0x0),
fh2DeltaEtaC60pt2lead(0x0),  
fh2DeltaEtaC10pt3lead(0x0),
fh2DeltaEtaC20pt3lead(0x0),
fh2DeltaEtaC30pt3lead(0x0),
fh2DeltaEtaC60pt3lead(0x0),  
fh2DeltaEtaC10pt4lead(0x0),
fh2DeltaEtaC20pt4lead(0x0),
fh2DeltaEtaC30pt4lead(0x0),
fh2DeltaEtaC60pt4lead(0x0),
fh2DeltaPhiC10pt1lead(0x0),
fh2DeltaPhiC20pt1lead(0x0),
fh2DeltaPhiC30pt1lead(0x0),
fh2DeltaPhiC60pt1lead(0x0),  
fh2DeltaPhiC10pt2lead(0x0),
fh2DeltaPhiC20pt2lead(0x0),
fh2DeltaPhiC30pt2lead(0x0),
fh2DeltaPhiC60pt2lead(0x0),  
fh2DeltaPhiC10pt3lead(0x0),
fh2DeltaPhiC20pt3lead(0x0),
fh2DeltaPhiC30pt3lead(0x0),
fh2DeltaPhiC60pt3lead(0x0),  
fh2DeltaPhiC10pt4lead(0x0),
fh2DeltaPhiC20pt4lead(0x0),
fh2DeltaPhiC30pt4lead(0x0),
fh2DeltaPhiC60pt4lead(0x0),
fh2DeltaRC10pt1sublead(0x0),
fh2DeltaRC20pt1sublead(0x0),
fh2DeltaRC30pt1sublead(0x0),
fh2DeltaRC60pt1sublead(0x0),  
fh2DeltaRC10pt2sublead(0x0),
fh2DeltaRC20pt2sublead(0x0),
fh2DeltaRC30pt2sublead(0x0),
fh2DeltaRC60pt2sublead(0x0),  
fh2DeltaRC10pt3sublead(0x0),
fh2DeltaRC20pt3sublead(0x0),
fh2DeltaRC30pt3sublead(0x0),
fh2DeltaRC60pt3sublead(0x0),  
fh2DeltaRC10pt4sublead(0x0),
fh2DeltaRC20pt4sublead(0x0),
fh2DeltaRC30pt4sublead(0x0),
fh2DeltaRC60pt4sublead(0x0),
fh2DeltaEtaC10pt1sublead(0x0),
fh2DeltaEtaC20pt1sublead(0x0),
fh2DeltaEtaC30pt1sublead(0x0),
fh2DeltaEtaC60pt1sublead(0x0),  
fh2DeltaEtaC10pt2sublead(0x0),
fh2DeltaEtaC20pt2sublead(0x0),
fh2DeltaEtaC30pt2sublead(0x0),
fh2DeltaEtaC60pt2sublead(0x0),  
fh2DeltaEtaC10pt3sublead(0x0),
fh2DeltaEtaC20pt3sublead(0x0),
fh2DeltaEtaC30pt3sublead(0x0),
fh2DeltaEtaC60pt3sublead(0x0),  
fh2DeltaEtaC10pt4sublead(0x0),
fh2DeltaEtaC20pt4sublead(0x0),
fh2DeltaEtaC30pt4sublead(0x0),
fh2DeltaEtaC60pt4sublead(0x0),
fh2DeltaPhiC10pt1sublead(0x0),
fh2DeltaPhiC20pt1sublead(0x0),
fh2DeltaPhiC30pt1sublead(0x0),
fh2DeltaPhiC60pt1sublead(0x0),  
fh2DeltaPhiC10pt2sublead(0x0),
fh2DeltaPhiC20pt2sublead(0x0),
fh2DeltaPhiC30pt2sublead(0x0),
fh2DeltaPhiC60pt2sublead(0x0),  
fh2DeltaPhiC10pt3sublead(0x0),
fh2DeltaPhiC20pt3sublead(0x0),
fh2DeltaPhiC30pt3sublead(0x0),
fh2DeltaPhiC60pt3sublead(0x0),  
fh2DeltaPhiC10pt4sublead(0x0),
fh2DeltaPhiC20pt4sublead(0x0),
fh2DeltaPhiC30pt4sublead(0x0),
fh2DeltaPhiC60pt4sublead(0x0),
fh2AngStructpt1C10(0x0),
fh2AngStructpt2C10(0x0),
fh2AngStructpt3C10(0x0),
fh2AngStructpt4C10(0x0),
fh2AngStructpt1C20(0x0),
fh2AngStructpt2C20(0x0),
fh2AngStructpt3C20(0x0),
fh2AngStructpt4C20(0x0),    
fh2AngStructpt1C30(0x0),
fh2AngStructpt2C30(0x0),
fh2AngStructpt3C30(0x0),
fh2AngStructpt4C30(0x0),   
fh2AngStructpt1C60(0x0),
fh2AngStructpt2C60(0x0),
fh2AngStructpt3C60(0x0),
fh2AngStructpt4C60(0x0)    

 {
   // Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;

   DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetCore::~AliAnalysisTaskJetCore()
{
   delete fListJets[0];
   delete fListJets[1];
}

void AliAnalysisTaskJetCore::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetCore::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetCore::UserCreateOutputObjects()
{
   // Create histograms
   // Called once
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList;
   fOutputList->SetOwner(kTRUE);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);


   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");

   fHistJetSelection = new TH1I("fHistJetSelection", "jet selection", 8, -0.5, 7.5);
   fHistJetSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistJetSelection->GetXaxis()->SetBinLabel(2,"probes IN");
   fHistJetSelection->GetXaxis()->SetBinLabel(3,"no matching jet");
   fHistJetSelection->GetXaxis()->SetBinLabel(4,"not in list");
   fHistJetSelection->GetXaxis()->SetBinLabel(5,"fraction cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(6,"acceptance cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(7,"p_{T} cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(8,"trigger exclude mask");

   fh2JetSelection = new TH2F("fh2JetSelection", "jet selection", 8, -0.5, 7.5,100,0.,200.);
   fh2JetSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fh2JetSelection->GetXaxis()->SetBinLabel(2,"probes IN");
   fh2JetSelection->GetXaxis()->SetBinLabel(3,"no matching jet");
   fh2JetSelection->GetXaxis()->SetBinLabel(4,"not in list");
   fh2JetSelection->GetXaxis()->SetBinLabel(5,"fraction cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(6,"acceptance cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(7,"p_{T} cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(8,"trigger exclude mask");


   //UInt_t entries = 0; // bit coded, see GetDimParams() below
   //UInt_t opt = 0;  // bit coded, default (0) or high resolution (1)

   //  Int_t bins[5]={10,10,150,150,15};
   //Double_t xmin[5]={0.,0.,0.,0.,0.};
   //Double_t xmax[5]={100.,6.5,150.,1.5,1.5};     
   //fhnDeltaRjets = new THnSparseF("fhnDeltaRjets", "fhnDeltaRjets",5,bins,xmin,xmax);



    fh2JetCoreMethod1C10 = new TH2F("JetCoreMethod1C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C10 = new TH2F("JetCoreMethod2C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C10 = new TH2F("JetCoreMethod3C10","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod1C20 = new TH2F("JetCoreMethod1C20","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C20 = new TH2F("JetCoreMethod2C20","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C20 = new TH2F("JetCoreMethod3C20","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod1C30 = new TH2F("JetCoreMethod1C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C30 = new TH2F("JetCoreMethod2C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C30 = new TH2F("JetCoreMethod3C30","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod1C60 = new TH2F("JetCoreMethod1C60","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C60 = new TH2F("JetCoreMethod2C60","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C60 = new TH2F("JetCoreMethod3C60","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod3C10lead = new TH2F("JetCoreMethod3C10lead","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod3C20lead = new TH2F("JetCoreMethod3C20lead","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C30lead = new TH2F("JetCoreMethod3C30lead","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C60lead = new TH2F("JetCoreMethod3C60lead","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod3C10sublead = new TH2F("JetCoreMethod3C10sublead","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod3C20sublead = new TH2F("JetCoreMethod3C20sublead","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C30sublead = new TH2F("JetCoreMethod3C30sublead","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C60sublead = new TH2F("JetCoreMethod3C60sublead","",150, 0., 150.,100, 0., 1.5); 
   
    fh2SumPtInC10 = new TH2F("PtSumInC10","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20 = new TH2F("PtSumInC20","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30 = new TH2F("PtSumInC30","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60 = new TH2F("PtSumInC60","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC10lead = new TH2F("PtSumInC10lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20lead = new TH2F("PtSumInC20lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30lead = new TH2F("PtSumInC30lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60lead = new TH2F("PtSumInC60lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC10sublead = new TH2F("PtSumInC10sublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20sublead = new TH2F("PtSumInC20sublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30sublead = new TH2F("PtSumInC30sublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60sublead = new TH2F("PtSumInC60sublead","",150, 0., 150.,5000, 0., 25.); 

    fh2SumPtOutC10 = new TH2F("PtSumOutC10","",150, 0., 150.,5000, 0., 25.);
    fh2SumPtOutC20 = new TH2F("PtSumOutC20","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30 = new TH2F("PtSumOutC30","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60 = new TH2F("PtSumOutC60","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC10lead = new TH2F("PtSumOutC10lead","",150, 0., 150.,5000, 0., 25.);
    fh2SumPtOutC20lead = new TH2F("PtSumOutC20lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30lead = new TH2F("PtSumOutC30lead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60lead = new TH2F("PtSumOutC60lead","",150, 0., 150.,5000, 0., 25.);      
    fh2SumPtOutC10sublead = new TH2F("PtSumOutC10sublead","",150, 0., 150.,5000, 0., 25.);
    fh2SumPtOutC20sublead = new TH2F("PtSumOutC20sublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30sublead = new TH2F("PtSumOutC30sublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60sublead = new TH2F("PtSumOutC60sublead","",150, 0., 150.,5000, 0., 25.);      


    fh2SumPtInC10bkg = new TH2F("PtSumInC10bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20bkg = new TH2F("PtSumInC20bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30bkg = new TH2F("PtSumInC30bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60bkg = new TH2F("PtSumInC60bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC10bkglead = new TH2F("PtSumInC10bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20bkglead = new TH2F("PtSumInC20bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30bkglead = new TH2F("PtSumInC30bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60bkglead = new TH2F("PtSumInC60bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC10bkgsublead = new TH2F("PtSumInC10bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC20bkgsublead = new TH2F("PtSumInC20bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC30bkgsublead = new TH2F("PtSumInC30bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtInC60bkgsublead = new TH2F("PtSumInC60bkgsublead","",150, 0., 150.,5000, 0., 25.); 


    fh2SumPtOutC10bkg = new TH2F("PtSumOutC10bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC20bkg = new TH2F("PtSumOutC20bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30bkg = new TH2F("PtSumOutC30bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60bkg = new TH2F("PtSumOutC60bkg","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC10bkglead = new TH2F("PtSumOutC10bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC20bkglead = new TH2F("PtSumOutC20bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30bkglead = new TH2F("PtSumOutC30bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60bkglead = new TH2F("PtSumOutC60bkglead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC10bkgsublead = new TH2F("PtSumOutC10bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC20bkgsublead = new TH2F("PtSumOutC20bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC30bkgsublead = new TH2F("PtSumOutC30bkgsublead","",150, 0., 150.,5000, 0., 25.); 
    fh2SumPtOutC60bkgsublead = new TH2F("PtSumOutC60bkgsublead","",150, 0., 150.,5000, 0., 25.); 



    fh2DeltaRC10pt1 = new TH2F("DeltaRC10pt1","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt1 = new TH2F("DeltaRC20pt1","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt1 = new TH2F("DeltaRC30pt1","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt1 = new TH2F("DeltaRC60pt1","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt2 = new TH2F("DeltaRC10pt2","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt2 = new TH2F("DeltaRC20pt2","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt2 = new TH2F("DeltaRC30pt2","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt2 = new TH2F("DeltaRC60pt2","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt3 = new TH2F("DeltaRC10pt3","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt3 = new TH2F("DeltaRC20pt3","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt3 = new TH2F("DeltaRC30pt3","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt3 = new TH2F("DeltaRC60pt3","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt4 = new TH2F("DeltaRC10pt4","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt4 = new TH2F("DeltaRC20pt4","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt4 = new TH2F("DeltaRC30pt4","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt4 = new TH2F("DeltaRC60pt4","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaEtaC10pt1 = new TH2F("DeltaEtaC10pt1","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt1 = new TH2F("DeltaEtaC20pt1","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt1 = new TH2F("DeltaEtaC30pt1","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt1 = new TH2F("DeltaEtaC60pt1","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt2 = new TH2F("DeltaEtaC10pt2","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt2 = new TH2F("DeltaEtaC20pt2","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt2 = new TH2F("DeltaEtaC30pt2","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt2 = new TH2F("DeltaEtaC60pt2","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt3 = new TH2F("DeltaEtaC10pt3","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt3 = new TH2F("DeltaEtaC20pt3","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt3 = new TH2F("DeltaEtaC30pt3","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt3 = new TH2F("DeltaEtaC60pt3","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt4 = new TH2F("DeltaEtaC10pt4","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt4 = new TH2F("DeltaEtaC20pt4","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt4 = new TH2F("DeltaEtaC30pt4","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt4 = new TH2F("DeltaEtaC60pt4","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaPhiC10pt1 = new TH2F("DeltaPhiC10pt1","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt1 = new TH2F("DeltaPhiC20pt1","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt1 = new TH2F("DeltaPhiC30pt1","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt1 = new TH2F("DeltaPhiC60pt1","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt2 = new TH2F("DeltaPhiC10pt2","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt2 = new TH2F("DeltaPhiC20pt2","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt2 = new TH2F("DeltaPhiC30pt2","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt2 = new TH2F("DeltaPhiC60pt2","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt3 = new TH2F("DeltaPhiC10pt3","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt3 = new TH2F("DeltaPhiC20pt3","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt3 = new TH2F("DeltaPhiC30pt3","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt3 = new TH2F("DeltaPhiC60pt3","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt4 = new TH2F("DeltaPhiC10pt4","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt4 = new TH2F("DeltaPhiC20pt4","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt4 = new TH2F("DeltaPhiC30pt4","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt4 = new TH2F("DeltaPhiC60pt4","",150, 0., 15.,100,-6.5,6.5); 

    fh2DeltaRC10pt1lead = new TH2F("DeltaRC10pt1lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt1lead = new TH2F("DeltaRC20pt1lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt1lead = new TH2F("DeltaRC30pt1lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt1lead = new TH2F("DeltaRC60pt1lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt2lead = new TH2F("DeltaRC10pt2lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt2lead = new TH2F("DeltaRC20pt2lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt2lead = new TH2F("DeltaRC30pt2lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt2lead = new TH2F("DeltaRC60pt2lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt3lead = new TH2F("DeltaRC10pt3lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt3lead = new TH2F("DeltaRC20pt3lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt3lead = new TH2F("DeltaRC30pt3lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt3lead = new TH2F("DeltaRC60pt3lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt4lead = new TH2F("DeltaRC10pt4lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt4lead = new TH2F("DeltaRC20pt4lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt4lead = new TH2F("DeltaRC30pt4lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt4lead = new TH2F("DeltaRC60pt4lead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaEtaC10pt1lead = new TH2F("DeltaEtaC10pt1lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt1lead = new TH2F("DeltaEtaC20pt1lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt1lead = new TH2F("DeltaEtaC30pt1lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt1lead = new TH2F("DeltaEtaC60pt1lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt2lead = new TH2F("DeltaEtaC10pt2lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt2lead = new TH2F("DeltaEtaC20pt2lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt2lead = new TH2F("DeltaEtaC30pt2lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt2lead = new TH2F("DeltaEtaC60pt2lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt3lead = new TH2F("DeltaEtaC10pt3lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt3lead = new TH2F("DeltaEtaC20pt3lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt3lead = new TH2F("DeltaEtaC30pt3lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt3lead = new TH2F("DeltaEtaC60pt3lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt4lead = new TH2F("DeltaEtaC10pt4lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt4lead = new TH2F("DeltaEtaC20pt4lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt4lead = new TH2F("DeltaEtaC30pt4lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt4lead = new TH2F("DeltaEtaC60pt4lead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaPhiC10pt1lead = new TH2F("DeltaPhiC10pt1lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt1lead = new TH2F("DeltaPhiC20pt1lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt1lead = new TH2F("DeltaPhiC30pt1lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt1lead = new TH2F("DeltaPhiC60pt1lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt2lead = new TH2F("DeltaPhiC10pt2lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt2lead = new TH2F("DeltaPhiC20pt2lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt2lead = new TH2F("DeltaPhiC30pt2lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt2lead = new TH2F("DeltaPhiC60pt2lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt3lead = new TH2F("DeltaPhiC10pt3lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt3lead = new TH2F("DeltaPhiC20pt3lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt3lead = new TH2F("DeltaPhiC30pt3lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt3lead = new TH2F("DeltaPhiC60pt3lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt4lead = new TH2F("DeltaPhiC10pt4lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt4lead = new TH2F("DeltaPhiC20pt4lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt4lead = new TH2F("DeltaPhiC30pt4lead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt4lead = new TH2F("DeltaPhiC60pt4lead","",150, 0., 15.,100,-6.5,6.5); 


    fh2DeltaRC10pt1sublead = new TH2F("DeltaRC10pt1sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt1sublead = new TH2F("DeltaRC20pt1sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt1sublead = new TH2F("DeltaRC30pt1sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt1sublead = new TH2F("DeltaRC60pt1sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt2sublead = new TH2F("DeltaRC10pt2sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt2sublead = new TH2F("DeltaRC20pt2sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt2sublead = new TH2F("DeltaRC30pt2sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt2sublead = new TH2F("DeltaRC60pt2sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt3sublead = new TH2F("DeltaRC10pt3sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt3sublead = new TH2F("DeltaRC20pt3sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt3sublead = new TH2F("DeltaRC30pt3sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt3sublead = new TH2F("DeltaRC60pt3sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC10pt4sublead = new TH2F("DeltaRC10pt4sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC20pt4sublead = new TH2F("DeltaRC20pt4sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC30pt4sublead = new TH2F("DeltaRC30pt4sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaRC60pt4sublead = new TH2F("DeltaRC60pt4sublead","",150, 0., 15.,100,0.,1.5); 
    fh2DeltaEtaC10pt1sublead = new TH2F("DeltaEtaC10pt1sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt1sublead = new TH2F("DeltaEtaC20pt1sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt1sublead = new TH2F("DeltaEtaC30pt1sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt1sublead = new TH2F("DeltaEtaC60pt1sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt2sublead = new TH2F("DeltaEtaC10pt2sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt2sublead = new TH2F("DeltaEtaC20pt2sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt2sublead = new TH2F("DeltaEtaC30pt2sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt2sublead = new TH2F("DeltaEtaC60pt2sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt3sublead = new TH2F("DeltaEtaC10pt3sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt3sublead = new TH2F("DeltaEtaC20pt3sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt3sublead = new TH2F("DeltaEtaC30pt3sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt3sublead = new TH2F("DeltaEtaC60pt3sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC10pt4sublead = new TH2F("DeltaEtaC10pt4sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC20pt4sublead = new TH2F("DeltaEtaC20pt4sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC30pt4sublead = new TH2F("DeltaEtaC30pt4sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaEtaC60pt4sublead = new TH2F("DeltaEtaC60pt4sublead","",150, 0., 15.,100,-1.5,1.5); 
    fh2DeltaPhiC10pt1sublead = new TH2F("DeltaPhiC10pt1sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt1sublead = new TH2F("DeltaPhiC20pt1sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt1sublead = new TH2F("DeltaPhiC30pt1sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt1sublead = new TH2F("DeltaPhiC60pt1sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt2sublead = new TH2F("DeltaPhiC10pt2sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt2sublead = new TH2F("DeltaPhiC20pt2sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt2sublead = new TH2F("DeltaPhiC30pt2sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt2sublead = new TH2F("DeltaPhiC60pt2sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt3sublead = new TH2F("DeltaPhiC10pt3sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt3sublead = new TH2F("DeltaPhiC20pt3sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt3sublead = new TH2F("DeltaPhiC30pt3sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt3sublead = new TH2F("DeltaPhiC60pt3sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC10pt4sublead = new TH2F("DeltaPhiC10pt4sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC20pt4sublead = new TH2F("DeltaPhiC20pt4sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC30pt4sublead = new TH2F("DeltaPhiC30pt4sublead","",150, 0., 15.,100,-6.5,6.5); 
    fh2DeltaPhiC60pt4sublead = new TH2F("DeltaPhiC60pt4sublead","",150, 0., 15.,100,-6.5,6.5); 





      if(fAngStructCloseTracks>0){
    fh2AngStructpt1C10 = new TH2F("Ang struct pt1 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C10 = new TH2F("Ang struct pt2 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C10 = new TH2F("Ang struct pt3 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C10 = new TH2F("Ang struct pt4 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C20 = new TH2F("Ang struct pt1 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C20 = new TH2F("Ang struct pt2 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C20 = new TH2F("Ang struct pt3 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C20 = new TH2F("Ang struct pt4 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C30 = new TH2F("Ang struct pt1 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C30 = new TH2F("Ang struct pt2 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C30 = new TH2F("Ang struct pt3 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C30 = new TH2F("Ang struct pt4 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C60 = new TH2F("Ang struct pt1 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C60 = new TH2F("Ang struct pt2 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C60 = new TH2F("Ang struct pt3 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C60 = new TH2F("Ang struct pt4 C60","",15,0.,1.5,150,0.,10.); }

   fOutputList->Add(fHistEvtSelection);
   fOutputList->Add(fHistJetSelection);
   fOutputList->Add(fh2JetSelection);

      fOutputList->Add(fh2JetCoreMethod1C10);
      fOutputList->Add(fh2JetCoreMethod2C10);
      fOutputList->Add(fh2JetCoreMethod3C10);
      fOutputList->Add(fh2JetCoreMethod1C20);
      fOutputList->Add(fh2JetCoreMethod2C20);
      fOutputList->Add(fh2JetCoreMethod3C20);
      fOutputList->Add(fh2JetCoreMethod1C30);
      fOutputList->Add(fh2JetCoreMethod2C30);
      fOutputList->Add(fh2JetCoreMethod3C30);
      fOutputList->Add(fh2JetCoreMethod1C60);
      fOutputList->Add(fh2JetCoreMethod2C60);
      fOutputList->Add(fh2JetCoreMethod3C60);
      
      fOutputList->Add(fh2JetCoreMethod3C10lead);
      fOutputList->Add(fh2JetCoreMethod3C20lead);
      fOutputList->Add(fh2JetCoreMethod3C30lead);
      fOutputList->Add(fh2JetCoreMethod3C60lead);
      fOutputList->Add(fh2JetCoreMethod3C10sublead);
      fOutputList->Add(fh2JetCoreMethod3C20sublead);
      fOutputList->Add(fh2JetCoreMethod3C30sublead);
      fOutputList->Add(fh2JetCoreMethod3C60sublead);

      fOutputList->Add(fh2SumPtInC10);
      fOutputList->Add(fh2SumPtInC20);
      fOutputList->Add(fh2SumPtInC30);
      fOutputList->Add(fh2SumPtInC60);
      fOutputList->Add(fh2SumPtInC10lead);
      fOutputList->Add(fh2SumPtInC20lead);
      fOutputList->Add(fh2SumPtInC30lead);
      fOutputList->Add(fh2SumPtInC60lead);
      fOutputList->Add(fh2SumPtInC10sublead);
      fOutputList->Add(fh2SumPtInC20sublead);
      fOutputList->Add(fh2SumPtInC30sublead);
      fOutputList->Add(fh2SumPtInC60sublead); 

      fOutputList->Add(fh2SumPtOutC10);
      fOutputList->Add(fh2SumPtOutC20);
      fOutputList->Add(fh2SumPtOutC30);
      fOutputList->Add(fh2SumPtOutC60);
      fOutputList->Add(fh2SumPtOutC10lead);
      fOutputList->Add(fh2SumPtOutC20lead);
      fOutputList->Add(fh2SumPtOutC30lead);
      fOutputList->Add(fh2SumPtOutC60lead);
      fOutputList->Add(fh2SumPtOutC10sublead);
      fOutputList->Add(fh2SumPtOutC20sublead);
      fOutputList->Add(fh2SumPtOutC30sublead);
      fOutputList->Add(fh2SumPtOutC60sublead);
      
      fOutputList->Add(fh2SumPtInC10bkg);
      fOutputList->Add(fh2SumPtInC20bkg);
      fOutputList->Add(fh2SumPtInC30bkg);
      fOutputList->Add(fh2SumPtInC60bkg);
      fOutputList->Add(fh2SumPtInC10bkglead);
      fOutputList->Add(fh2SumPtInC20bkglead);
      fOutputList->Add(fh2SumPtInC30bkglead);
      fOutputList->Add(fh2SumPtInC60bkglead);
      fOutputList->Add(fh2SumPtInC10bkgsublead);
      fOutputList->Add(fh2SumPtInC20bkgsublead);
      fOutputList->Add(fh2SumPtInC30bkgsublead);
      fOutputList->Add(fh2SumPtInC60bkgsublead);

      fOutputList->Add(fh2SumPtOutC10bkg);
      fOutputList->Add(fh2SumPtOutC20bkg);
      fOutputList->Add(fh2SumPtOutC30bkg);
      fOutputList->Add(fh2SumPtOutC60bkg);     
      fOutputList->Add(fh2SumPtOutC10bkglead);
      fOutputList->Add(fh2SumPtOutC20bkglead);
      fOutputList->Add(fh2SumPtOutC30bkglead);
      fOutputList->Add(fh2SumPtOutC60bkglead);     
      fOutputList->Add(fh2SumPtOutC10bkgsublead);
      fOutputList->Add(fh2SumPtOutC20bkgsublead);
      fOutputList->Add(fh2SumPtOutC30bkgsublead);
      fOutputList->Add(fh2SumPtOutC60bkgsublead);     
      
      fOutputList->Add(fh2DeltaRC10pt1);
      fOutputList->Add(fh2DeltaRC20pt1);
      fOutputList->Add(fh2DeltaRC30pt1);
      fOutputList->Add(fh2DeltaRC60pt1); 
      fOutputList->Add(fh2DeltaRC10pt2);
      fOutputList->Add(fh2DeltaRC20pt2);
      fOutputList->Add(fh2DeltaRC30pt2);
      fOutputList->Add(fh2DeltaRC60pt2); 
      fOutputList->Add(fh2DeltaRC10pt3);
      fOutputList->Add(fh2DeltaRC20pt3);
      fOutputList->Add(fh2DeltaRC30pt3);
      fOutputList->Add(fh2DeltaRC60pt3); 
      fOutputList->Add(fh2DeltaRC10pt4);
      fOutputList->Add(fh2DeltaRC20pt4);
      fOutputList->Add(fh2DeltaRC30pt4);
      fOutputList->Add(fh2DeltaRC60pt4);
      fOutputList->Add(fh2DeltaEtaC10pt1);
      fOutputList->Add(fh2DeltaEtaC20pt1);
      fOutputList->Add(fh2DeltaEtaC30pt1);
      fOutputList->Add(fh2DeltaEtaC60pt1); 
      fOutputList->Add(fh2DeltaEtaC10pt2);
      fOutputList->Add(fh2DeltaEtaC20pt2);
      fOutputList->Add(fh2DeltaEtaC30pt2);
      fOutputList->Add(fh2DeltaEtaC60pt2); 
      fOutputList->Add(fh2DeltaEtaC10pt3);
      fOutputList->Add(fh2DeltaEtaC20pt3);
      fOutputList->Add(fh2DeltaEtaC30pt3);
      fOutputList->Add(fh2DeltaEtaC60pt3); 
      fOutputList->Add(fh2DeltaEtaC10pt4);
      fOutputList->Add(fh2DeltaEtaC20pt4);
      fOutputList->Add(fh2DeltaEtaC30pt4);
      fOutputList->Add(fh2DeltaEtaC60pt4); 
      fOutputList->Add(fh2DeltaPhiC10pt1);
      fOutputList->Add(fh2DeltaPhiC20pt1);
      fOutputList->Add(fh2DeltaPhiC30pt1);
      fOutputList->Add(fh2DeltaPhiC60pt1); 
      fOutputList->Add(fh2DeltaPhiC10pt2);
      fOutputList->Add(fh2DeltaPhiC20pt2);
      fOutputList->Add(fh2DeltaPhiC30pt2);
      fOutputList->Add(fh2DeltaPhiC60pt2); 
      fOutputList->Add(fh2DeltaPhiC10pt3);
      fOutputList->Add(fh2DeltaPhiC20pt3);
      fOutputList->Add(fh2DeltaPhiC30pt3);
      fOutputList->Add(fh2DeltaPhiC60pt3); 
      fOutputList->Add(fh2DeltaPhiC10pt4);
      fOutputList->Add(fh2DeltaPhiC20pt4);
      fOutputList->Add(fh2DeltaPhiC30pt4);
      fOutputList->Add(fh2DeltaPhiC60pt4);      

      fOutputList->Add(fh2DeltaRC10pt1lead);
      fOutputList->Add(fh2DeltaRC20pt1lead);
      fOutputList->Add(fh2DeltaRC30pt1lead);
      fOutputList->Add(fh2DeltaRC60pt1lead); 
      fOutputList->Add(fh2DeltaRC10pt2lead);
      fOutputList->Add(fh2DeltaRC20pt2lead);
      fOutputList->Add(fh2DeltaRC30pt2lead);
      fOutputList->Add(fh2DeltaRC60pt2lead); 
      fOutputList->Add(fh2DeltaRC10pt3lead);
      fOutputList->Add(fh2DeltaRC20pt3lead);
      fOutputList->Add(fh2DeltaRC30pt3lead);
      fOutputList->Add(fh2DeltaRC60pt3lead); 
      fOutputList->Add(fh2DeltaRC10pt4lead);
      fOutputList->Add(fh2DeltaRC20pt4lead);
      fOutputList->Add(fh2DeltaRC30pt4lead);
      fOutputList->Add(fh2DeltaRC60pt4lead);
      fOutputList->Add(fh2DeltaEtaC10pt1lead);
      fOutputList->Add(fh2DeltaEtaC20pt1lead);
      fOutputList->Add(fh2DeltaEtaC30pt1lead);
      fOutputList->Add(fh2DeltaEtaC60pt1lead); 
      fOutputList->Add(fh2DeltaEtaC10pt2lead);
      fOutputList->Add(fh2DeltaEtaC20pt2lead);
      fOutputList->Add(fh2DeltaEtaC30pt2lead);
      fOutputList->Add(fh2DeltaEtaC60pt2lead); 
      fOutputList->Add(fh2DeltaEtaC10pt3lead);
      fOutputList->Add(fh2DeltaEtaC20pt3lead);
      fOutputList->Add(fh2DeltaEtaC30pt3lead);
      fOutputList->Add(fh2DeltaEtaC60pt3lead); 
      fOutputList->Add(fh2DeltaEtaC10pt4lead);
      fOutputList->Add(fh2DeltaEtaC20pt4lead);
      fOutputList->Add(fh2DeltaEtaC30pt4lead);
      fOutputList->Add(fh2DeltaEtaC60pt4lead); 
      fOutputList->Add(fh2DeltaPhiC10pt1lead);
      fOutputList->Add(fh2DeltaPhiC20pt1lead);
      fOutputList->Add(fh2DeltaPhiC30pt1lead);
      fOutputList->Add(fh2DeltaPhiC60pt1lead); 
      fOutputList->Add(fh2DeltaPhiC10pt2lead);
      fOutputList->Add(fh2DeltaPhiC20pt2lead);
      fOutputList->Add(fh2DeltaPhiC30pt2lead);
      fOutputList->Add(fh2DeltaPhiC60pt2lead); 
      fOutputList->Add(fh2DeltaPhiC10pt3lead);
      fOutputList->Add(fh2DeltaPhiC20pt3lead);
      fOutputList->Add(fh2DeltaPhiC30pt3lead);
      fOutputList->Add(fh2DeltaPhiC60pt3lead); 
      fOutputList->Add(fh2DeltaPhiC10pt4lead);
      fOutputList->Add(fh2DeltaPhiC20pt4lead);
      fOutputList->Add(fh2DeltaPhiC30pt4lead);
      fOutputList->Add(fh2DeltaPhiC60pt4lead);      



      fOutputList->Add(fh2DeltaRC10pt1sublead);
      fOutputList->Add(fh2DeltaRC20pt1sublead);
      fOutputList->Add(fh2DeltaRC30pt1sublead);
      fOutputList->Add(fh2DeltaRC60pt1sublead); 
      fOutputList->Add(fh2DeltaRC10pt2sublead);
      fOutputList->Add(fh2DeltaRC20pt2sublead);
      fOutputList->Add(fh2DeltaRC30pt2sublead);
      fOutputList->Add(fh2DeltaRC60pt2sublead); 
      fOutputList->Add(fh2DeltaRC10pt3sublead);
      fOutputList->Add(fh2DeltaRC20pt3sublead);
      fOutputList->Add(fh2DeltaRC30pt3sublead);
      fOutputList->Add(fh2DeltaRC60pt3sublead); 
      fOutputList->Add(fh2DeltaRC10pt4sublead);
      fOutputList->Add(fh2DeltaRC20pt4sublead);
      fOutputList->Add(fh2DeltaRC30pt4sublead);
      fOutputList->Add(fh2DeltaRC60pt4sublead);
      fOutputList->Add(fh2DeltaEtaC10pt1sublead);
      fOutputList->Add(fh2DeltaEtaC20pt1sublead);
      fOutputList->Add(fh2DeltaEtaC30pt1sublead);
      fOutputList->Add(fh2DeltaEtaC60pt1sublead); 
      fOutputList->Add(fh2DeltaEtaC10pt2sublead);
      fOutputList->Add(fh2DeltaEtaC20pt2sublead);
      fOutputList->Add(fh2DeltaEtaC30pt2sublead);
      fOutputList->Add(fh2DeltaEtaC60pt2sublead); 
      fOutputList->Add(fh2DeltaEtaC10pt3sublead);
      fOutputList->Add(fh2DeltaEtaC20pt3sublead);
      fOutputList->Add(fh2DeltaEtaC30pt3sublead);
      fOutputList->Add(fh2DeltaEtaC60pt3sublead); 
      fOutputList->Add(fh2DeltaEtaC10pt4sublead);
      fOutputList->Add(fh2DeltaEtaC20pt4sublead);
      fOutputList->Add(fh2DeltaEtaC30pt4sublead);
      fOutputList->Add(fh2DeltaEtaC60pt4sublead); 
      fOutputList->Add(fh2DeltaPhiC10pt1sublead);
      fOutputList->Add(fh2DeltaPhiC20pt1sublead);
      fOutputList->Add(fh2DeltaPhiC30pt1sublead);
      fOutputList->Add(fh2DeltaPhiC60pt1sublead); 
      fOutputList->Add(fh2DeltaPhiC10pt2sublead);
      fOutputList->Add(fh2DeltaPhiC20pt2sublead);
      fOutputList->Add(fh2DeltaPhiC30pt2sublead);
      fOutputList->Add(fh2DeltaPhiC60pt2sublead); 
      fOutputList->Add(fh2DeltaPhiC10pt3sublead);
      fOutputList->Add(fh2DeltaPhiC20pt3sublead);
      fOutputList->Add(fh2DeltaPhiC30pt3sublead);
      fOutputList->Add(fh2DeltaPhiC60pt3sublead); 
      fOutputList->Add(fh2DeltaPhiC10pt4sublead);
      fOutputList->Add(fh2DeltaPhiC20pt4sublead);
      fOutputList->Add(fh2DeltaPhiC30pt4sublead);
      fOutputList->Add(fh2DeltaPhiC60pt4sublead);      
     
       



        if(fAngStructCloseTracks>0){
       fOutputList->Add(fh2AngStructpt1C10);
       fOutputList->Add(fh2AngStructpt2C10);
       fOutputList->Add(fh2AngStructpt3C10);
       fOutputList->Add(fh2AngStructpt4C10); 
       fOutputList->Add(fh2AngStructpt1C20);
       fOutputList->Add(fh2AngStructpt2C20);
       fOutputList->Add(fh2AngStructpt3C20);
       fOutputList->Add(fh2AngStructpt4C20); 
       fOutputList->Add(fh2AngStructpt1C30);
       fOutputList->Add(fh2AngStructpt2C30);
       fOutputList->Add(fh2AngStructpt3C30);
       fOutputList->Add(fh2AngStructpt4C30);
       fOutputList->Add(fh2AngStructpt1C60);
       fOutputList->Add(fh2AngStructpt2C60);
       fOutputList->Add(fh2AngStructpt3C60);
       fOutputList->Add(fh2AngStructpt4C60);}  




   // =========== Switch on Sumw2 for all histos ===========
   for (Int_t i=0; i<fOutputList->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if (h1){
         h1->Sumw2();
         continue;
      }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if (hn){
         hn->Sumw2();
      }	  
   }
   TH1::AddDirectory(oldStatus);

   PostData(1, fOutputList);
}

void AliAnalysisTaskJetCore::UserExec(Option_t *)
{
   

   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
      return;
   }

   fESD=dynamic_cast<AliESDEvent*>(InputEvent());
   if (!fESD) {
      AliError("ESD not available");
      fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
   } else {
      fAOD = dynamic_cast<AliAODEvent*>(AODEvent());
   }
 
    if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }}
    




   // -- event selection --
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fOutputList);
      return;
   }

   // vertex selection
   AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
   Int_t nTracksPrim = primVtx->GetNContributors();
   if ((nTracksPrim < fMinContribVtx) ||
         (primVtx->GetZ() < fVtxZMin) ||
         (primVtx->GetZ() > fVtxZMax) ){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }

   // event class selection (from jet helper task)
   Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   if(fDebug) Printf("Event class %d", eventClass);
   if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
   }

   // centrality selection
   AliCentrality *cent = 0x0;
   Float_t centValue = 0.; 
   if(fESD) cent = fESD->GetCentrality();
   if(cent) centValue = cent->GetCentralityPercentile("V0M");
   if(fDebug) printf("centrality: %f\n", centValue);
   if (centValue < fCentMin || centValue > fCentMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
   }


   // multiplicity due to input tracks
   //Int_t nInputTracks = GetNInputTracks();
   //if (nInputTracks < fNInputTracksMin || (fNInputTracksMax > -1 && nInputTracks > fNInputTracksMax)){
   //   fHistEvtSelection->Fill(5);
   //   PostData(1, fOutputList);
   //   return;
   // }

   
   fHistEvtSelection->Fill(0); 
   // accepted events  
   // -- end event selection --
  
   // get background
   AliAODJetEventBackground* externalBackground = 0;
   if(fAOD&&!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAOD->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   if(fAODExtension&&!externalBackground&&fBackgroundBranch.Length()){
     externalBackground =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   
   Float_t rho = 0;
   if(externalBackground)rho = externalBackground->GetBackground(0);


   // fetch jets
   TClonesArray *aodJets[2];
   aodJets[0]=0;
   if(fAOD&&!aodJets[0]){
   aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data()));  }
   if(fAODExtension && !aodJets[0]){ 
     aodJets[0] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[0].Data())); 
     aodJets[1] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[1].Data()));  }

   Double_t ptsub[aodJets[0]->GetEntriesFast()];
   Int_t inord[aodJets[0]->GetEntriesFast()];
   for(Int_t n=0;n<aodJets[0]->GetEntriesFast();n++){
     ptsub[n]=0;
     inord[n]=0;}   

   TList ParticleList;
   Int_t nT = GetListOfTracks(&ParticleList);
     for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
      fListJets[iJetType]->Clear();
      if (!aodJets[iJetType]) continue;

      if(fDebug) Printf("%s: %d jets",fJetBranchName[iJetType].Data(),aodJets[iJetType]->GetEntriesFast());
      
   
      for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
         if (jet) fListJets[iJetType]->Add(jet);
         if(iJetType==0){
	   ptsub[iJet]=jet->Pt()-rho*jet->EffectiveAreaCharged();
	 }}}
   
   Double_t etabig=0;
   Double_t ptbig=0;
   Double_t areabig=0;
   Double_t phibig=0.;
   Double_t etasmall=0;
   Double_t ptsmall=0;
   Double_t areasmall=0;
   Double_t distr=0.;
   Double_t phismall=0.;
   Int_t indexlead=-1;
   Int_t indexsublead=-1;
   Int_t indexstop=-1;
  
   if(fListJets[0]->GetEntries()>0) TMath::Sort(fListJets[0]->GetEntries(),ptsub,inord);

   for(Int_t jj=0;jj<fListJets[0]->GetEntries();jj++){
   AliAODJet* jetlead = (AliAODJet*)(fListJets[0]->At(inord[jj]));
   if(jetlead->Pt()-rho*jetlead->EffectiveAreaCharged()<=0) continue;
   if((jetlead->Eta()<fJetEtaMin)||(jetlead->Eta()>fJetEtaMax)) continue;
   indexlead=inord[jj];
   indexstop=jj;
   break;}
   if((indexstop>-1)&&(indexstop+1<fListJets[0]->GetEntries()-1)){
   for(Int_t k=indexstop+1;k<fListJets[0]->GetEntries();k++){
   AliAODJet* jetsublead = (AliAODJet*)(fListJets[0]->At(inord[k]));
   if(jetsublead->Pt()-rho*jetsublead->EffectiveAreaCharged()<=0) continue;
   if((jetsublead->Eta()<fJetEtaMin)||(jetsublead->Eta()>fJetEtaMax)) continue;
   indexsublead=inord[k];
   break;}}
         

   Double_t up1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};



   for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){
           AliAODJet* jetbig = (AliAODJet*)(fListJets[0]->At(i));
           etabig  = jetbig->Eta();
           phibig  = jetbig->Phi();
           ptbig   = jetbig->Pt();
           if(ptbig==0) continue; 
           areabig = jetbig->EffectiveAreaCharged();
           Double_t ptcorr=ptbig-rho*areabig;
           if(ptcorr<=0) continue;
      	   if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue;
                   Double_t dismin=100.;
                   Double_t ptmax=-10.; 
                   Int_t index1=-1;
                   Int_t index2=-1;
                   Double_t fracin=0.;
                   Double_t sumPtIn=0.;
                   Double_t sumPtOut=0.; 
	   //compute sum of the pt of the tracks in a concentric cone
           TRefArray *genTrackList = jetbig->GetRefTracks();
           Int_t nTracksGenJet = genTrackList->GetEntriesFast();
           AliAODTrack* genTrack;
             for(Int_t ir=0; ir<nTracksGenJet; ++ir){
             genTrack = (AliAODTrack*)(genTrackList->At(ir));
	     Float_t etr=genTrack->Eta();
             Float_t phir=genTrack->Phi();
             distr=(etr-etabig)*(etr-etabig)+(phir-phibig)*(phir-phibig);
             distr=TMath::Sqrt(distr);
	     if(distr<=fRadioFrac){ fracin=fracin+genTrack->Pt();}}
    
	     if(centValue<10) {fh2JetCoreMethod3C10->Fill(ptcorr,fracin/ptbig);
	       if(i==indexlead) fh2JetCoreMethod3C10lead->Fill(ptcorr,fracin/ptbig);
               if(i==indexsublead) fh2JetCoreMethod3C10sublead->Fill(ptcorr,fracin/ptbig);}

	      
             if((centValue>20)&&(centValue<40)) {fh2JetCoreMethod3C20->Fill(ptcorr,fracin/ptbig);
                                                 if(i==indexlead) fh2JetCoreMethod3C20lead->Fill(ptcorr,fracin/ptbig);
						 if(i==indexsublead) fh2JetCoreMethod3C20sublead->Fill(ptcorr,fracin/ptbig);}
             if((centValue>30)&&(centValue<60)){ fh2JetCoreMethod3C30->Fill(ptcorr,fracin/ptbig);
                                                if(i==indexlead) fh2JetCoreMethod3C30lead->Fill(ptcorr,fracin/ptbig);
                                                if(i==indexsublead) fh2JetCoreMethod3C30sublead->Fill(ptcorr,fracin/ptbig);}
             if(centValue>60){fh2JetCoreMethod3C60->Fill(ptcorr,fracin/ptbig);
                                                if(i==indexlead) fh2JetCoreMethod3C60lead->Fill(ptcorr,fracin/ptbig);
		                        	if(i==indexsublead) fh2JetCoreMethod3C60sublead->Fill(ptcorr,fracin/ptbig);}
   
	     

                for(Int_t j=0; j<fListJets[1]->GetEntries(); ++j){
                  AliAODJet* jetsmall = (AliAODJet*)(fListJets[1]->At(j));
                  etasmall  = jetsmall->Eta();
                  phismall = jetsmall->Phi();
                  ptsmall   = jetsmall->Pt();
                  areasmall = jetsmall->EffectiveAreaCharged();
                  Double_t tmpDeltaR=(phismall-phibig)*(phismall-phibig)+(etasmall-etabig)*(etasmall-etabig);
		  tmpDeltaR=TMath::Sqrt(tmpDeltaR);
		     //Fraction in the jet core  
                    if((ptsmall>ptmax)&&(tmpDeltaR<=fRadioFrac)){ptmax=ptsmall;  
		    index2=j;}  
                    if(tmpDeltaR<=dismin){ dismin=tmpDeltaR;
		      index1=j;}} //en of loop over R=0.2 jets
                //method1:most concentric jet=core 
		if(dismin<fMinDist){ AliAODJet* jetmethod1 = (AliAODJet*)(fListJets[1]->At(index1));       
      		  if(centValue<10) fh2JetCoreMethod1C10->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if((centValue>20)&&(centValue<40)) fh2JetCoreMethod1C20->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod1C30->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod1C60->Fill(ptcorr,jetmethod1->Pt()/ptbig); }
                //method2:hardest contained jet=core   
		if(index2!=-1){ 
                  AliAODJet* jetmethod2 = (AliAODJet*)(fListJets[1]->At(index2));
                  if(centValue<10) fh2JetCoreMethod2C10->Fill(ptcorr,jetmethod2->Pt()/ptbig);
                  if((centValue>20)&&(centValue<40)) fh2JetCoreMethod2C20->Fill(ptcorr,jetmethod2->Pt()/ptbig); 
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod2C30->Fill(ptcorr,jetmethod2->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod2C60->Fill(ptcorr,jetmethod2->Pt()/ptbig); }  

     Double_t R=fRadioFrac*2.;      
  
     for(int it = 0;it<nT;++it){

	  AliVParticle *part = (AliVParticle*)ParticleList.At(it);
	  Float_t deltaR = jetbig->DeltaR(part);
          Float_t deltaEta = part->Eta()-etabig;
          Float_t deltaPhi = part->Phi()-phibig;  
      	  if((deltaR>=0.4)&&(deltaR<=0.6))sumPtIn=sumPtIn+part->Pt();                     
          if((deltaR>=0.8)&&(deltaR<=1.))sumPtOut=sumPtOut+part->Pt();     
	  if(centValue<10.){
	    
	       if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC10pt1->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC10pt1->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt1->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC10pt2->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt2->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt2->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC10pt3->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt3->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt3->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC10pt4->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt4->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt4->Fill(part->Pt(),deltaPhi); }

	       if(i==indexlead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC10pt1lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC10pt1lead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt1lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC10pt2lead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt2lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt2lead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC10pt3lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt3lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt3lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC10pt4lead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt4lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt4lead->Fill(part->Pt(),deltaPhi); }}


               if(i==indexsublead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC10pt1sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC10pt1sublead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt1sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC10pt2sublead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt2sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt2sublead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC10pt3sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt3sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC10pt3sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC10pt4sublead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC10pt4sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC10pt4sublead->Fill(part->Pt(),deltaPhi); }} 

 


	  }

          if((centValue>20.)&&(centValue<40.)){
              if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC20pt1->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC20pt1->Fill(part->Pt(),deltaEta);
	      if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt1->Fill(part->Pt(),deltaPhi);}
       	      if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC20pt2->Fill(part->Pt(),deltaR);
              if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt2->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt2->Fill(part->Pt(),deltaPhi); }
	      if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC20pt3->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt3->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt3->Fill(part->Pt(),deltaPhi);}
	      if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC20pt4->Fill(part->Pt(),deltaR); 
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt4->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt4->Fill(part->Pt(),deltaPhi); }

               if(i==indexlead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC20pt1lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC20pt1lead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt1lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC20pt2lead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt2lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt2lead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC20pt3lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt3lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt3lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC20pt4lead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt4lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt4lead->Fill(part->Pt(),deltaPhi); }}
                             
               if(i==indexsublead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC20pt1sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC20pt1sublead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt1sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC20pt2sublead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt2sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt2sublead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC20pt3sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt3sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC20pt3sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC20pt4sublead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC20pt4sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC20pt4sublead->Fill(part->Pt(),deltaPhi); }} 

	  }

	  if((centValue>30.)&&(centValue<60.)){

              if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC30pt1->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC30pt1->Fill(part->Pt(),deltaEta);
	      if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt1->Fill(part->Pt(),deltaPhi);}
	      if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC30pt2->Fill(part->Pt(),deltaR);
              if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt2->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt2->Fill(part->Pt(),deltaPhi); }
	      if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC30pt3->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt3->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt3->Fill(part->Pt(),deltaPhi);}
	      if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC30pt4->Fill(part->Pt(),deltaR); 
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt4->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt4->Fill(part->Pt(),deltaPhi); }

               if(i==indexlead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC30pt1lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC30pt1lead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt1lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC30pt2lead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt2lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt2lead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC30pt3lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt3lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt3lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC30pt4lead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt4lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt4lead->Fill(part->Pt(),deltaPhi); }}


                if(i==indexsublead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC30pt1sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC30pt1sublead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt1sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC30pt2sublead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt2sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt2sublead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC30pt3sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt3sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC30pt3sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC30pt4sublead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC30pt4sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC30pt4sublead->Fill(part->Pt(),deltaPhi); }} 





	  }
                 

	  if(centValue>60.){
              if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC60pt1->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC60pt1->Fill(part->Pt(),deltaEta);
	      if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt1->Fill(part->Pt(),deltaPhi);}
	      if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC60pt2->Fill(part->Pt(),deltaR);
              if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt2->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt2->Fill(part->Pt(),deltaPhi); }
	      if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC60pt3->Fill(part->Pt(),deltaR);
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt3->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt3->Fill(part->Pt(),deltaPhi);}
	      if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC60pt4->Fill(part->Pt(),deltaR); 
	      if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt4->Fill(part->Pt(),deltaEta);
              if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt4->Fill(part->Pt(),deltaPhi);}

               if(i==indexlead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC60pt1lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC60pt1lead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt1lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC60pt2lead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt2lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt2lead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC60pt3lead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt3lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt3lead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC60pt4lead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt4lead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt4lead->Fill(part->Pt(),deltaPhi); }}

                if(i==indexsublead){
               if((ptcorr>=70.)&&(ptcorr<85.)) {fh2DeltaRC60pt1sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R))fh2DeltaEtaC60pt1sublead->Fill(part->Pt(),deltaEta);
	       if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt1sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=85.)&&(ptcorr<100.)) {fh2DeltaRC60pt2sublead->Fill(part->Pt(),deltaR);
               if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt2sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt2sublead->Fill(part->Pt(),deltaPhi); }
	       if((ptcorr>=100.)&&(ptcorr<120.)) {fh2DeltaRC60pt3sublead->Fill(part->Pt(),deltaR);
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt3sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R)) fh2DeltaPhiC60pt3sublead->Fill(part->Pt(),deltaPhi);}
	       if((ptcorr>=120.)&&(ptcorr<140.)) {fh2DeltaRC60pt4sublead->Fill(part->Pt(),deltaR); 
	       if((part->Phi()<=phibig+R)&&(part->Phi()>=phibig-R)) fh2DeltaEtaC60pt4sublead->Fill(part->Pt(),deltaEta);
               if((part->Eta()<=etabig+R)&&(part->Eta()>=etabig-R))fh2DeltaPhiC60pt4sublead->Fill(part->Pt(),deltaPhi); }} 

	  }
                 
     } 
     //end of track loop
     Double_t coronain=rho*TMath::Pi()*(1.-0.8*0.8);
     Double_t coronaout=rho*TMath::Pi()*(0.6*0.6-0.4*0.4);
     if(centValue<10.){  
     fh2SumPtInC10bkg->Fill(ptcorr,coronain/ptbig);       
     fh2SumPtOutC10bkg->Fill(ptcorr,coronaout/ptbig);           
     fh2SumPtInC10->Fill(ptcorr,sumPtIn/ptbig);
     fh2SumPtOutC10->Fill(ptcorr,sumPtOut/ptbig); 
     if(i==indexlead){  fh2SumPtInC10bkglead->Fill(ptcorr,coronain/ptbig);       
                        fh2SumPtOutC10bkglead->Fill(ptcorr,coronaout/ptbig);           
                        fh2SumPtInC10lead->Fill(ptcorr,sumPtIn/ptbig);
                        fh2SumPtOutC10lead->Fill(ptcorr,sumPtOut/ptbig);}    
     if(i==indexsublead){ fh2SumPtInC10bkgsublead->Fill(ptcorr,coronain/ptbig);       
                         fh2SumPtOutC10bkgsublead->Fill(ptcorr,coronaout/ptbig);           
                         fh2SumPtInC10sublead->Fill(ptcorr,sumPtIn/ptbig);
                         fh2SumPtOutC10sublead->Fill(ptcorr,sumPtOut/ptbig);}
     }
     if((centValue>20.)&&(centValue<40.)){             
     fh2SumPtInC20bkg->Fill(ptcorr,coronain/ptbig);
     fh2SumPtOutC20bkg->Fill(ptcorr,coronaout/ptbig);    
     fh2SumPtInC20->Fill(ptcorr,sumPtIn/ptbig);
     fh2SumPtOutC20->Fill(ptcorr,sumPtOut/ptbig); 
      if(i==indexlead){ fh2SumPtInC20bkglead->Fill(ptcorr,coronain/ptbig);       
                        fh2SumPtOutC20bkglead->Fill(ptcorr,coronaout/ptbig);           
                        fh2SumPtInC20lead->Fill(ptcorr,sumPtIn/ptbig);
                        fh2SumPtOutC20lead->Fill(ptcorr,sumPtOut/ptbig);}    
     if(i==indexsublead){ fh2SumPtInC20bkgsublead->Fill(ptcorr,coronain/ptbig);       
                         fh2SumPtOutC20bkgsublead->Fill(ptcorr,coronaout/ptbig);           
                         fh2SumPtInC20sublead->Fill(ptcorr,sumPtIn/ptbig);
                         fh2SumPtOutC20sublead->Fill(ptcorr,sumPtOut/ptbig);}



     }
     if((centValue>30.)&&(centValue<60.)){
     fh2SumPtInC30bkg->Fill(ptcorr,coronain/ptbig);
     fh2SumPtOutC30bkg->Fill(ptcorr,coronaout/ptbig);                    
     fh2SumPtInC30->Fill(ptcorr,sumPtIn/ptbig);
     fh2SumPtOutC30->Fill(ptcorr,sumPtOut/ptbig); 

     if(i==indexlead){  fh2SumPtInC30bkglead->Fill(ptcorr,coronain/ptbig);       
                        fh2SumPtOutC30bkglead->Fill(ptcorr,coronaout/ptbig);           
                        fh2SumPtInC30lead->Fill(ptcorr,sumPtIn/ptbig);
                        fh2SumPtOutC30lead->Fill(ptcorr,sumPtOut/ptbig);}    
     if(i==indexsublead){ fh2SumPtInC30bkgsublead->Fill(ptcorr,coronain/ptbig);       
                         fh2SumPtOutC30bkgsublead->Fill(ptcorr,coronaout/ptbig);           
                         fh2SumPtInC30sublead->Fill(ptcorr,sumPtIn/ptbig);
                         fh2SumPtOutC30sublead->Fill(ptcorr,sumPtOut/ptbig);}

     }
     if(centValue>60.){         
     fh2SumPtInC60bkg->Fill(ptcorr,coronain/ptbig);
     fh2SumPtOutC60bkg->Fill(ptcorr,coronaout/ptbig);      
     fh2SumPtInC60->Fill(ptcorr,sumPtIn/ptbig);
     fh2SumPtOutC60->Fill(ptcorr,sumPtOut/ptbig);

     if(i==indexlead){  fh2SumPtInC60bkglead->Fill(ptcorr,coronain/ptbig);       
                        fh2SumPtOutC60bkglead->Fill(ptcorr,coronaout/ptbig);           
                        fh2SumPtInC60lead->Fill(ptcorr,sumPtIn/ptbig);
                        fh2SumPtOutC60lead->Fill(ptcorr,sumPtOut/ptbig);}    
     if(i==indexsublead){ fh2SumPtInC60bkgsublead->Fill(ptcorr,coronain/ptbig);       
                         fh2SumPtOutC60bkgsublead->Fill(ptcorr,coronaout/ptbig);           
                         fh2SumPtInC60sublead->Fill(ptcorr,sumPtIn/ptbig);
                         fh2SumPtOutC60sublead->Fill(ptcorr,sumPtOut/ptbig);}

     }      
  
     //////////////////ANGULAR STRUCTURE//////////////////////////////////////

     //tracks up to R=0.8 distant from the jet axis
     if(fAngStructCloseTracks==1){
      TList CloseTrackList;
      Int_t nn=GetListOfTracksCloseToJet(&CloseTrackList,jetbig);
      Double_t difR=0.04;
      for(Int_t l=0;l<15;l++){
      Double_t rr=l*0.1+0.1;
       for(int it = 0;it<nn;++it){
           AliVParticle *part1 = (AliVParticle*)CloseTrackList.At(it);
           for(int itu=it+1;itu<CloseTrackList.GetEntries();itu++){      
           AliVParticle *part2 = (AliVParticle*)CloseTrackList.At(itu);  
           Double_t ptm=part1->Pt();
           Double_t ptn=part2->Pt();	
           Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
                      Rnm=TMath::Sqrt(Rnm);
                      Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));      
                      Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR)));
                      if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
                      if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}
     }
    
     //only jet constituents
      if(fAngStructCloseTracks==2){

      Double_t difR=0.04;
      for(Int_t l=0;l<15;l++){
      Double_t rr=l*0.1+0.1;

    
      AliAODTrack* part1;
      AliAODTrack* part2;
          for(Int_t it=0; it<nTracksGenJet; ++it){
             part1 = (AliAODTrack*)(genTrackList->At(it));
           for(Int_t itu=0; itu<nTracksGenJet; ++itu){
             part2 = (AliAODTrack*)(genTrackList->At(itu));
           Double_t ptm=part1->Pt();
           Double_t ptn=part2->Pt();	
           Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
                      Rnm=TMath::Sqrt(Rnm);
                      Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));
                      Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR))); 
                      if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
                      if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}}


   












   }


     //end loop over R=0.4 jets	
      if(fAngStructCloseTracks>0){
     for(Int_t l=0;l<15;l++){
     Double_t rr=l*0.1+0.1;
        if(down1[l]!=0){  
	if(centValue<10.)fh2AngStructpt1C10->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt1C20->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt1C30->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>60.) fh2AngStructpt1C60->Fill(rr,rr*up1[l]/down1[l]);}
        if(down2[l]!=0){  
	if(centValue<10.) fh2AngStructpt2C10->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt2C20->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt2C30->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>60.) fh2AngStructpt2C60->Fill(rr,rr*up2[l]/down2[l]);}
        if(down3[l]!=0){  
	if(centValue<10.) fh2AngStructpt3C10->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt3C20->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt3C30->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>60.) fh2AngStructpt3C60->Fill(rr,rr*up3[l]/down3[l]);}
        if(down4[l]!=0){  
	if(centValue<10.) fh2AngStructpt4C10->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt4C20->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt4C30->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>60.) fh2AngStructpt4C60->Fill(rr,rr*up4[l]/down4[l]);}}}

    





   PostData(1, fOutputList);
     }

void AliAnalysisTaskJetCore::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
}











Int_t  AliAnalysisTaskJetCore::GetListOfTracks(TList *list){

    Int_t iCount = 0;
 
    
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   list->Sort();
  return iCount;
 
}

 Int_t  AliAnalysisTaskJetCore::GetListOfTracksCloseToJet(TList *list,AliAODJet *jetbig){

    Int_t iCount = 0;
 
  
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      Double_t disR=jetbig->DeltaR(tr);
      if(disR>0.8)  continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   list->Sort();
   return iCount;

}











Int_t AliAnalysisTaskJetCore::GetNInputTracks()
{

   Int_t nInputTracks = 0;

   TString jbname(fJetBranchName[1]);
   //needs complete event, use jets without background subtraction
   for(Int_t i=1; i<=3; ++i){
      if(jbname.Contains(Form("B%d",i))) jbname.ReplaceAll(Form("B%d",i),"B0");
   }
   // use only HI event
   if(jbname.Contains("AODextraonly")) jbname.ReplaceAll("AODextraonly","AOD");
   if(jbname.Contains("AODextra")) jbname.ReplaceAll("AODextra","AOD");

   if(fDebug) Printf("Multiplicity from jet branch %s", jbname.Data());
   TClonesArray *tmpAODjets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(jbname.Data()));
   if(!tmpAODjets){
      Printf("Jet branch %s not found", jbname.Data());
      Printf("AliAnalysisTaskJetCore::GetNInputTracks FAILED");
      return -1;
   }
   
   for (Int_t iJet=0; iJet<tmpAODjets->GetEntriesFast(); iJet++){
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*tmpAODjets)[iJet]);
      if(!jet) continue;
      TRefArray *trackList = jet->GetRefTracks();
      Int_t nTracks = trackList->GetEntriesFast();
      nInputTracks += nTracks;
      if(fDebug) Printf("#jet%d: %d tracks", iJet, nTracks);
   }
   if(fDebug) Printf("---> input tracks: %d", nInputTracks);

   return nInputTracks;  
}



Double_t AliAnalysisTaskJetCore::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}



