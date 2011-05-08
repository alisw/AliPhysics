/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
// Class AliAnalysisTaskSPDdNdEta                                       //
// Analysis task for dN/dEta reconstruction with the SPD                //
//                                                                      //
// Author:  M. Nicassio (INFN Bari)                                     //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it         //
//////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"

#include "TH1F.h"
#include "TH2F.h" 
#include "TH3F.h" 
#include "TList.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"

#include "AliESDInputHandlerRP.h"
#include "AliESDVZERO.h"
#include "../ITS/AliITSRecPoint.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackReference.h"

#include "AliGenHijingEventHeader.h" 
#include "AliGenDPMjetEventHeader.h"

#include "AliLog.h"

#include "AliTriggerAnalysis.h" 
#include "AliPhysicsSelection.h"
#include "AliTrackletAlg.h" 
#include "AliAnalysisTaskSPDdNdEta.h"


ClassImp(AliAnalysisTaskSPDdNdEta)
//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::AliAnalysisTaskSPDdNdEta(const char *name) 
  : AliAnalysisTaskSE(name), 

  fmyESD(0), 
  fOutput(0), 
  
  fUseMC(kFALSE), 
  fTrigger(AliTriggerAnalysis::kAcceptAll),
  fTR(kFALSE),
  fRecoTracklets(kFALSE),

  fMCCentralityBin(AliAnalysisTaskSPDdNdEta::kall),
  fCentrLowLim(0),
  fCentrUpLim(0),
  fCentrEst(kFALSE),
  fMinClMultLay2(0),
  fMaxClMultLay2(0),
  fMinMultV0(0),
  fVtxLim(0),
  fPartSpecies(0),

  fPhiWindow(0),
  fThetaWindow(0),
  fPhiShift(0),
  fRemoveClustersFromOverlaps(0),
  fPhiOverlapCut(0),
  fZetaOverlapCut(0),
  fPhiRotationAngle(0),
  fPhiWindowAna(0),

  fMultReco(0), 

  fV0Ampl(0),
  fHistSPDRAWMultvsZ(0),
  fHistSPDRAWMultvsZTriggCentrEvts(0),
  fHistSPDRAWMultvsZCentrEvts(0),
  fHistSPDRAWEtavsZ(0),
 
  fHistSPDmultEtacut(0),
  fHistSPDmult(0),
  fHistSPDmultcl1(0),
  fHistSPDmultcl2(0),
  fHistSPDmultcl1vscl2(0),
  fHistSPDmultvscl1(0),
  fHistSPDmultvscl2(0),

  fHistSPDeta(0),
  fHistSPDphi(0),
  fHistSPDtheta(0),
  fHistSPDdePhi(0),
  fHistSPDphivsSPDeta(0),  
  fHistSPDdeTheta(0),
  fHistSPDvtx(0),
  fHistSPDvtxAnalysis(0),
  fHistSPDdePhideTheta(0),

  fHistSPDphicl1(0),
  fHistSPDphicl2(0),

  fHistBkgCorrDen(0),
  fHistReconstructedProtons(0),
  fHistReconstructedKaons(0),
  fHistReconstructedPions(0),
  fHistReconstructedOthers(0), 
  fHistReconstructedSec(0),
  fHistBkgCorrDenPrimGen(0),
  fHistBkgCombLabels(0),
  fHistBkgCorrNum(0),
  fHistAlgEffNum(0),
  fHistNonDetectableCorrDen(0),

  fHistNonDetectableCorrNum(0),
  fHistProtons(0),
  fHistKaons(0),
  fHistPions(0),
  fHistOthers(0),
  fHistAllPrimaries(0),
  fHistTrackCentrEvts(0),
  fHistTrackTrigCentrEvts(0),

  fHistAllEvts(0),
  fHistCentrEvts(0),
  fHistTrigCentrEvts(0),
  fHistSelEvts(0),

  fHistMCmultEtacut(0),
  fHistMCmultEtacutvsSPDmultEtacut(0),
  fHistMCmultEtacutvsSPDmultcl1(0),
  fHistMCmultEtacutvsSPDmultcl2(0),
 
  fHistMCvtxx(0),
  fHistMCvtxy(0),
  fHistMCvtxz(0),

  fHistRecvsGenImpactPar(0),
  fHistMCNpart(0), 

  fHistdPhidThetaPP(0),
  fHistdPhidThetaSS(0),
  fHistdPhidThetaComb(0),

  fHistDeVtx(0)

{

  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
//  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  DefineOutput(1, TList::Class());  

}
//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::~AliAnalysisTaskSPDdNdEta()
{

  // Destructor

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  delete fMultReco;

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }

}
//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::UserCreateOutputObjects() 
{
  if (fRecoTracklets) {
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
    if (fUseMC) man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
    else man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
    man->SetRun(137161); 
    AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
    if(!obj) AliFatal("Unable to load geometry from CDB!");
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    AliGeomManager::GetNalignable("ITS");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  }

  fMultReco = new AliTrackletAlg();

  // Create histograms
  fOutput = new TList();
  fOutput->SetOwner(); 

  Int_t nBinVtx = 40;
  Double_t MaxVtx = 20.;

  Int_t nBinMultCorr = 200;
  Float_t nMaxMult = 20000.;
  Double_t binLimMultCorr[nBinMultCorr+1];
  binLimMultCorr[0] = 0.;
  binLimMultCorr[1] = 1.;
  for (Int_t i = 2; i<=nBinMultCorr;++i) {
    binLimMultCorr[i] = (i-1)*nMaxMult/nBinMultCorr;
  }

  Int_t nBinEtaCorr = 60; 
  Float_t etaMax = 3.;
  Float_t etaMin = -3.;
  Double_t *binLimEtaCorr = new Double_t[nBinEtaCorr+1];
  for (Int_t i = 0; i<nBinEtaCorr+1; ++i) {
    binLimEtaCorr[i] = (Double_t) etaMin+i*(etaMax*2.)/nBinEtaCorr;
  } 

  // Data to be corrected
  // ...event level    
  
  fV0Ampl = new TH1F("fV0Ampl","",500,0.,30000);
  fOutput->Add(fV0Ampl);

  fHistSPDRAWMultvsZ = new TH2F("fHistSPDRAWMultvsZ", "",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDRAWMultvsZ->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fOutput->Add(fHistSPDRAWMultvsZ); 

  fHistSPDRAWMultvsZTriggCentrEvts = new TH2F("fHistSPDRAWMultvsZTriggCentrEvts", "",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDRAWMultvsZTriggCentrEvts->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZTriggCentrEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fOutput->Add(fHistSPDRAWMultvsZTriggCentrEvts);

  fHistSPDRAWMultvsZCentrEvts = new TH2F("fHistSPDRAWMultvsZCentrEvts", "",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDRAWMultvsZCentrEvts->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZCentrEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fOutput->Add(fHistSPDRAWMultvsZCentrEvts);

  // ...track level
  fHistSPDRAWEtavsZ = new TH2F("fHistSPDRAWEtavsZ", "Tracklet pseudorapidity distribution", nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDRAWEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDRAWEtavsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fOutput->Add(fHistSPDRAWEtavsZ);
  
  fHistSPDmultEtacut = new TH1F("fHistSPDmultEtacut", "Tracklet multiplicity distribution",nBinMultCorr,binLimMultCorr);
  fHistSPDmultEtacut->GetXaxis()->SetTitle("Tracklet multiplicity (|#eta|<1.4)");
  fHistSPDmultEtacut->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultEtacut);

  fHistSPDmult = new TH1F("fHistSPDmult", "Tracklet multiplicity distribution", nBinMultCorr,binLimMultCorr);
  fHistSPDmult->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDmult->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmult);

  fHistSPDmultcl1 = new TH1F("fHistSPDmultcl1", "Inner layer cluster multiplicity", nBinMultCorr,binLimMultCorr);
  fHistSPDmultcl1->GetXaxis()->SetTitle("Inner layer cluster multiplicity");
  fHistSPDmultcl1->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultcl1);

  fHistSPDmultcl2 = new TH1F("fHistSPDmultcl2", "Outer layer cluster multiplicity", nBinMultCorr,binLimMultCorr);
  fHistSPDmultcl2->GetXaxis()->SetTitle("Outer layer cluster multiplicity");
  fHistSPDmultcl2->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultcl2);

  fHistSPDmultcl1vscl2 = new TH2F("fHistSPDmultcl1vscl2", "Inner layer cluster vs outer multiplicity", nBinMultCorr,binLimMultCorr, nBinMultCorr,binLimMultCorr);
  fHistSPDmultcl1vscl2->GetXaxis()->SetTitle("Inner layer cluster multiplicity");
  fHistSPDmultcl1vscl2->GetYaxis()->SetTitle("Outer layer cluster multiplicity");
  fOutput->Add(fHistSPDmultcl1vscl2);

  fHistSPDmultvscl1 = new TH2F("fHistSPDmultvscl1", "Tracklet vs inner layer cluster multiplicity", nBinMultCorr,binLimMultCorr, nBinMultCorr,binLimMultCorr);
  fHistSPDmultvscl1->GetXaxis()->SetTitle("Inner layer cluster multiplicity");
  fHistSPDmultvscl1->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvscl1);

  fHistSPDmultvscl2 = new TH2F("fHistSPDmultvscl2", "Tracklet vs outer layer cluster multiplicity", nBinMultCorr,binLimMultCorr, nBinMultCorr,binLimMultCorr);
  fHistSPDmultvscl2->GetXaxis()->SetTitle("Outer layer cluster multiplicity");
  fHistSPDmultvscl2->GetYaxis()->SetTitle("Tracklet multiplicity");
  fOutput->Add(fHistSPDmultvscl2);

  fHistSPDphi = new TH1F("fHistSPDphi", "Tracklet #phi  distribution", 360, 0.,2*TMath::Pi());
  fHistSPDphi->GetXaxis()->SetTitle("#varphi [rad]");
  fHistSPDphi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDphi);

  fHistSPDtheta = new TH1F("fHistSPDtheta", "Tracklet #theta  distribution", 180, 0.,TMath::Pi());
  fHistSPDtheta->GetXaxis()->SetTitle("#theta [rad]");
  fHistSPDtheta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDtheta);

  fHistSPDdePhi= new TH1F("fHistSPDdePhi", "Tracklet #Delta#varphi  distribution",400,-0.1,.1);
  fHistSPDdePhi->GetXaxis()->SetTitle("#Delta#varphi [rad]");
  fHistSPDdePhi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDdePhi);

  fHistSPDphivsSPDeta= new TH2F("fHistSPDphivsSPDeta", "Tracklets - #varphi vs #eta",120,-3.,3,360,0.,2*TMath::Pi());
  fHistSPDphivsSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDphivsSPDeta->GetYaxis()->SetTitle("#varphi [rad]");
  fOutput->Add(fHistSPDphivsSPDeta);

  fHistSPDdeTheta= new TH1F("fHistSPDdeTheta", "Tracklet #Delta#theta distribution",100,-0.05,.05);
  fHistSPDdeTheta->GetXaxis()->SetTitle("#Delta#theta [rad]");
  fHistSPDdeTheta->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDdeTheta);

  fHistSPDvtx = new TH1F("fHistSPDvtx", "SPD vertex distribution",nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDvtx->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtx->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtx);

  fHistSPDvtxAnalysis = new TH1F("fHistSPDvtxAnalysis", "SPD vertex distribution: events selected for the analysis",nBinVtx,-MaxVtx,MaxVtx);
  fHistSPDvtxAnalysis->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtxAnalysis->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtxAnalysis);

  fHistSPDdePhideTheta= new TH2F("fHistSPDdePhideTheta", "Tracklet #Delta#varphi  distribution",2000,-1.,1.,1000,-0.25,.25);
  fHistSPDdePhideTheta->GetXaxis()->SetTitle("#Delta#varphi [rad]");
  fHistSPDdePhideTheta->GetYaxis()->SetTitle("#Delta#theta [rad]");
  fOutput->Add(fHistSPDdePhideTheta);

  fHistSPDphicl1 = new TH1F("fHistSPDphicl1", "Tracklet #phi  distribution", 360, 0.,2*TMath::Pi());
  fHistSPDphicl1->GetXaxis()->SetTitle("#varphi [rad]");
  fHistSPDphicl1->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDphicl1);

  fHistSPDphicl2 = new TH1F("fHistSPDphicl2", "Tracklet #phi  distribution", 360, 0.,2*TMath::Pi());
  fHistSPDphicl2->GetXaxis()->SetTitle("#varphi [rad]");
  fHistSPDphicl2->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDphicl2);


  if (fUseMC) {

    // Track level correction histograms 
    fHistBkgCorrDen = new TH2F("fHistBkgCorrDen","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistBkgCorrDen);
   
    fHistReconstructedProtons = new TH2F("fHistReconstructedProtons","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistReconstructedProtons);
    fHistReconstructedKaons = new TH2F("fHistReconstructedKaons","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistReconstructedKaons);
    fHistReconstructedPions = new TH2F("fHistReconstructedPions","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistReconstructedPions);
    fHistReconstructedOthers = new TH2F("fHistReconstructedOthers","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistReconstructedOthers);
 
    fHistBkgCorrDenPrimGen = new TH2F("fHistBkgCorrDenPrimGen","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistBkgCorrDenPrimGen);

    fHistBkgCombLabels = new TH2F("fHistBkgCombLabels","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistBkgCombLabels);

    if (fTR) {
      fHistBkgCorrNum = new TH2F("fHistBkgCorrNum","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
      fOutput->Add(fHistBkgCorrNum);
 
      fHistAlgEffNum = new TH2F("fHistAlgEffNum","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
      fOutput->Add(fHistAlgEffNum);  

      fHistNonDetectableCorrDen = new TH2F("fHistNonDetectableCorrDen","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
      fOutput->Add(fHistNonDetectableCorrDen);

    }

    fHistNonDetectableCorrNum = new TH2F("fHistNonDetectableCorrNum","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistNonDetectableCorrNum);

    fHistProtons = new TH2F("fHistProtons","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistProtons);
    fHistKaons = new TH2F("fHistKaons","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistKaons);
    fHistPions = new TH2F("fHistPions","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistPions);
    fHistOthers = new TH2F("fHistOthers","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistOthers);
    fHistReconstructedSec = new TH2F("fHistReconstructedSec","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistReconstructedSec); 

 
    fHistAllPrimaries = new TH2F("fHistAllPrimaries","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistAllPrimaries);

    fHistTrackCentrEvts = new TH2F("fHistTrackCentrEvts","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistTrackCentrEvts);

    fHistTrackTrigCentrEvts = new TH2F("fHistTrackTrigCentrEvts","",nBinEtaCorr,binLimEtaCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistTrackTrigCentrEvts);


    // Event level correction histograms  
    fHistAllEvts = new TH2F("fHistAllEvts","",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistAllEvts);

    fHistCentrEvts = new TH2F("fHistCentrEvts","",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistCentrEvts);

    fHistTrigCentrEvts = new TH2F("fHistTrigCentrEvts","",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistTrigCentrEvts);

    fHistSelEvts = new TH2F("fHistSelEvts","",nBinMultCorr,binLimMultCorr,nBinVtx,-MaxVtx,MaxVtx);
    fOutput->Add(fHistSelEvts);


    // MC distributions
    fHistMCmultEtacut = new TH1F("fHistMCmultEtacut","Generated multiplicity",nBinMultCorr,binLimMultCorr);
    fHistMCmultEtacut->GetXaxis()->SetTitle("Generated multiplicity |#eta|<1.4");
    fHistMCmultEtacut->GetYaxis()->SetTitle("Entries");
    fOutput->Add(fHistMCmultEtacut);

    fHistMCmultEtacutvsSPDmultEtacut = new TH2F("fHistMCmultEtacutvsSPDmultEtacut","",nBinMultCorr,binLimMultCorr,nBinMultCorr,binLimMultCorr);
    fHistMCmultEtacutvsSPDmultEtacut->GetXaxis()->SetTitle("Generated multiplicity |#eta|<1.4");
    fHistMCmultEtacutvsSPDmultEtacut->GetYaxis()->SetTitle("Tracklet multiplicity |#eta|<1.4");
    fOutput->Add(fHistMCmultEtacutvsSPDmultEtacut);

    fHistMCmultEtacutvsSPDmultcl1 = new TH2F("fHistMCmultEtacutvsSPDmultcl1","",nBinMultCorr,binLimMultCorr,nBinMultCorr,binLimMultCorr);
    fHistMCmultEtacutvsSPDmultcl1->GetXaxis()->SetTitle("Generated multiplicity |#eta|<1.4");
    fHistMCmultEtacutvsSPDmultcl1->GetYaxis()->SetTitle("Cluster inner layer multiplicity");
    fOutput->Add(fHistMCmultEtacutvsSPDmultcl1);

    fHistMCmultEtacutvsSPDmultcl2 = new TH2F("fHistMCmultEtacutvsSPDmultcl2","",nBinMultCorr,binLimMultCorr,nBinMultCorr,binLimMultCorr);
    fHistMCmultEtacutvsSPDmultcl2->GetXaxis()->SetTitle("Generated multiplicity |#eta|<1.4");
    fHistMCmultEtacutvsSPDmultcl2->GetYaxis()->SetTitle("Cluster outer layer multiplicity");
    fOutput->Add(fHistMCmultEtacutvsSPDmultcl2);


    fHistMCvtxx = new TH1F("fHistMCvtxx", "MC vertex distribution - x",100,-.5,.5);
    fOutput->Add(fHistMCvtxx);
    fHistMCvtxy = new TH1F("fHistMCvtxy", "MC vertex distribution - y",100,-.5,.5);
    fOutput->Add(fHistMCvtxy);
    fHistMCvtxz = new TH1F("fHistMCvtxz", "MC vertex distribution - z",500,-50.,50.);
    fOutput->Add(fHistMCvtxz);

    fHistRecvsGenImpactPar= new TH2F("fHistRecvsGenImpactPar","",20,0.,20.,20,0.,20.);
    fHistRecvsGenImpactPar->GetXaxis()->SetTitle("b_{gen} [fm]");
    fHistRecvsGenImpactPar->GetYaxis()->SetTitle("b_{rec} [fm]");
    fOutput->Add(fHistRecvsGenImpactPar);

    fHistMCNpart = new TH1F("fHistMCNpart","Number of participants",450,0,450);
    fHistMCNpart->GetXaxis()->SetTitle("N_{part}");
    fHistMCNpart->GetYaxis()->SetTitle("Entries"); 
    fOutput->Add(fHistMCNpart);
 
    fHistdPhidThetaPP = new TH2F("fHistdPhidThetaPP","",2000,-1.,1.,1000,-0.25,.25);
    fOutput->Add(fHistdPhidThetaPP);
    fHistdPhidThetaSS = new TH2F("fHistdPhidThetaSS","",2000,-1.,1.,1000,-0.25,.25);
    fOutput->Add(fHistdPhidThetaSS);
    fHistdPhidThetaComb = new TH2F("fHistdPhidThetaComb","",2000,-1.,1.,1000,-0.25,.25);
    fOutput->Add(fHistdPhidThetaComb);

    fHistDeVtx = new TH2F("fHistDeVtx","",80,-20.,20.,5000,-0.5,0.5);
    fHistDeVtx->GetXaxis()->SetTitle("z_{MCvtx} [cm]");
    fHistDeVtx->GetYaxis()->SetTitle("z_{MCvtx}-z_{SPDvtx} [cm]");
    fOutput->Add(fHistDeVtx);

  }
  PostData(1, fOutput);
//  Printf("Output objects created...");
}

//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::UserExec(Option_t *) 
{
  // Main loop

  // Called for each event
//  Printf("User exec..........");

  AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!hand) { printf("No RP handler\n"); return; }

  fmyESD = dynamic_cast<AliESDEvent*>(InputEvent()); 
  if (!fmyESD) {
    Printf("ERROR: fESD not available");
    return;
  }
    
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field && !fmyESD->InitMagneticField()) {Printf("Failed to initialize the B field\n");return;}


  // Trigger selection
  Bool_t eventTriggered = kTRUE;
  static AliTriggerAnalysis* triggerAnalysis = 0; 
  if (eventTriggered) // applying an offline trigger
    eventTriggered = triggerAnalysis->IsTriggerFired(fmyESD, fTrigger);
  if (!eventTriggered)
    Printf("No trigger");


  // Centrality selection 
  Bool_t eventInCentralityBin = kFALSE;
/*  AliESDCentrality *centrality = fmyESD->GetCentrality();
  if (fCentrEst=="") eventInCentralityBin = kTRUE;
  else {
    if(!centrality) {
      AliError("Centrality object not available"); 
    }  else {
      if (centrality->IsEventInCentralityClass(fCentrLowLim,fCentrUpLim,fCentrEst.Data())) eventInCentralityBin = kTRUE;
    }
  }
*/
  if (fCentrEst) {
    AliESDVZERO* esdV0 = fmyESD->GetVZEROData();
    Float_t multV0A=esdV0->GetMTotV0A();
    Float_t multV0C=esdV0->GetMTotV0C();
    fV0Ampl->Fill(multV0A+multV0C);
    if (multV0A+multV0C>=fMinMultV0) eventInCentralityBin = kTRUE; 
  } else if (!fCentrEst) {
    eventInCentralityBin = kTRUE;
  }

  const AliMultiplicity* multESD = fmyESD->GetMultiplicity();

  // ESD vertex
  Bool_t eventWithVertex = kFALSE;
  const AliESDVertex* vtxESD = fmyESD->GetVertex();
  const AliESDVertex* vtxTPC = fmyESD->GetPrimaryVertexTPC();
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  Int_t nContrib = vtxESD->GetNContributors();
  Int_t nContribTPC = vtxTPC->GetNContributors(); 
  if (nContrib>0&&nContribTPC>0) {
    if (vtxESD->GetDispersion()<0.04) 
      if (vtxESD->GetZRes()<0.25) 
        if (nContribTPC>(-10.+0.25*multESD->GetNumberOfITSClusters(0))) {
          fHistSPDvtx->Fill(esdvtx[2]);
          if (TMath::Abs(esdvtx[2])<fVtxLim)  
            eventWithVertex = kTRUE;   
        }
  } 

  // Reconstructing or loading tracklets...
  Int_t multSPD = multESD->GetNumberOfTracklets();
  Int_t nSingleCl1 = multESD->GetNumberOfSingleClusters();
  Int_t multSPDcl1 = nSingleCl1 + multSPD;  
  Int_t multSPDEtacut = 0;
  Int_t multSPDAna = 0;

  Int_t multSPDcl2 = multESD->GetNumberOfITSClusters(1);

  Float_t trackletCoord[multSPDcl1][5];
  Float_t trackletLab[multSPDcl1][2];

  Bool_t eventSelected = kFALSE;
  // Event selection: in centrality bin, triggered with vertex
  if (eventTriggered&&eventWithVertex&&eventInCentralityBin&&multSPDcl2>fMinClMultLay2&&multSPDcl2<fMaxClMultLay2) {
    eventSelected = kTRUE; 
    fHistSPDmultcl1->Fill(multSPDcl1);
    fHistSPDmultcl2->Fill(multSPDcl2);
    fHistSPDmultcl1vscl2->Fill(multSPDcl1,multSPDcl2);
    fHistSPDmultvscl1->Fill(multSPDcl1,multSPD);
    fHistSPDmultvscl2->Fill(multSPDcl2,multSPD);
    if (fRecoTracklets) {
      // Load ITS rec points tree
      TTree *itsClusterTree = hand->GetTreeR("ITS");
      if (!itsClusterTree) {
        AliError(" Invalid ITS cluster tree !\n");
        return;
      }
      float vtxf[3] = {vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ()};

      fMultReco->SetPhiWindow(fPhiWindow);
      fMultReco->SetThetaWindow(fThetaWindow); 
      fMultReco->SetPhiShift(fPhiShift); 
      fMultReco->SetRemoveClustersFromOverlaps(fRemoveClustersFromOverlaps);
      fMultReco->SetPhiOverlapCut(fPhiOverlapCut);
      fMultReco->SetZetaOverlapCut(fZetaOverlapCut); 
      fMultReco->SetPhiRotationAngle(fPhiRotationAngle);
      fMultReco->SetHistOn(kFALSE); 

      fMultReco->Reconstruct(itsClusterTree,vtxf,vtxf);
//      Printf("cl 1 from alg %d",fMultReco->GetNClustersLayer1());
//      Printf("cl 2 from alg %d",fMultReco->GetNClustersLayer2());
      multSPD = fMultReco->GetNTracklets(); 
//      Printf("tracklets found...%d",multSPD);
      nSingleCl1 = fMultReco->GetNSingleClusters();
//      Printf("singl found...%d",nSingleCl1);
      for (Int_t itr = 0; itr<multSPD;++itr) {
        trackletCoord[itr][3] = fMultReco->GetTracklet(itr)[0]; //theta
        trackletCoord[itr][2] = fMultReco->GetTracklet(itr)[1]; //phi
        trackletCoord[itr][1] = fMultReco->GetTracklet(itr)[3]; //deTheta
        trackletCoord[itr][0] = fMultReco->GetTracklet(itr)[2]; //dephi
        trackletCoord[itr][4] = -TMath::Log(TMath::Tan(trackletCoord[itr][3]/2.)); //eta
        trackletLab[itr][0] = static_cast<Int_t>(fMultReco->GetTracklet(itr)[4]);  //label lay1
        trackletLab[itr][1] = static_cast<Int_t>(fMultReco->GetTracklet(itr)[5]);  //label lay2
      }
      for (Int_t icl1 = 0; icl1<multSPDcl1;++icl1) {
        fHistSPDphicl1->Fill(fMultReco->GetClusterLayer1(icl1)[1]); 
      }
      for (Int_t icl2 = 0; icl2<multSPDcl2; ++icl2) {
        fHistSPDphicl2->Fill(fMultReco->GetClusterLayer2(icl2)[1]);          
      }

    } else {
      // set variables from ESD
      for (Int_t itr = 0; itr<multSPD;++itr) {
        trackletCoord[itr][3] = multESD->GetTheta(itr);       //theta
        trackletCoord[itr][2] = multESD->GetPhi(itr);         //phi
        trackletCoord[itr][1] = multESD->GetDeltaTheta(itr);  //deTheta
        trackletCoord[itr][0] = multESD->GetDeltaPhi(itr);    //dePhi
        trackletCoord[itr][4] = multESD->GetEta(itr);         //eta 
        trackletLab[itr][0] = multESD->GetLabel(itr,0);       //label lay1
        trackletLab[itr][1] = multESD->GetLabel(itr,1);       //label lay2
      }
    }
//    Printf("tracklets in ESD...%d",multESD->GetNumberOfTracklets());

    for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
      if (TMath::Abs(trackletCoord[itracklet][0])<fPhiWindowAna) { // select tracklets with tighter cuts
        fHistSPDRAWEtavsZ->Fill(trackletCoord[itracklet][4],esdvtx[2]);
        fHistSPDphi->Fill(trackletCoord[itracklet][2]);
        fHistSPDtheta->Fill(trackletCoord[itracklet][3]);
        fHistSPDdePhi->Fill(trackletCoord[itracklet][0]);
        fHistSPDdeTheta->Fill(trackletCoord[itracklet][1]);
 
        fHistSPDphivsSPDeta->Fill(trackletCoord[itracklet][4],trackletCoord[itracklet][2]);
        multSPDAna++;
        // Calculate multiplicity in etacut
        if (TMath::Abs(trackletCoord[itracklet][4])<1.4) multSPDEtacut++;
      }
      // for combinatorial background normalization
      fHistSPDdePhideTheta->Fill(trackletCoord[itracklet][0], trackletCoord[itracklet][1]);

    }
    if (multSPDAna!=0) fHistSPDvtxAnalysis->Fill(esdvtx[2]);
    fHistSPDmultEtacut->Fill(multSPDEtacut);
    fHistSPDmult->Fill(multSPDAna); 

    fHistSPDRAWMultvsZ->Fill(multSPDAna,esdvtx[2]);

  } // End selected events

  if (eventInCentralityBin) {
    fHistSPDRAWMultvsZCentrEvts->Fill(multSPDAna,esdvtx[2]);  
    if (eventTriggered) fHistSPDRAWMultvsZTriggCentrEvts->Fill(multSPDAna,esdvtx[2]);
  }

  if (fUseMC) {

    if (!fMCEvent) {
      AliError("No MC info found");
      return;
    } 
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliStack* stack = fMCEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    Float_t impactParameter = 0.;
    Int_t npart = 0;

    AliGenHijingEventHeader* hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(fMCEvent->Header()->GenEventHeader());
    AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(fMCEvent->Header()->GenEventHeader());

    if (hijingHeader) impactParameter = hijingHeader->ImpactParameter();
    else if (dpmHeader) impactParameter = dpmHeader->ImpactParameter();

    Bool_t IsEventInMCCentralityBin = kFALSE;
    switch (fMCCentralityBin) {

      case 1: if (impactParameter<3) IsEventInMCCentralityBin = kTRUE;
      break;
      case 2: if (impactParameter>3&&impactParameter<6) IsEventInMCCentralityBin = kTRUE;
      break;
      case 3: if (impactParameter>6&&impactParameter<9) IsEventInMCCentralityBin = kTRUE;
      break;
      case 4: if (impactParameter>9&&impactParameter<12) IsEventInMCCentralityBin = kTRUE;
      break;
      case 5: if (impactParameter>12&&impactParameter<15) IsEventInMCCentralityBin = kTRUE;
      break;
      case 6: if (impactParameter>15) IsEventInMCCentralityBin = kTRUE;
      break;
      case 7:  IsEventInMCCentralityBin = kTRUE;
      break;

    }

    if (IsEventInMCCentralityBin) {
      // MC vertex
      TArrayF vtxMC(3);
      fMCEvent->GenEventHeader()->PrimaryVertex(vtxMC);
      fHistMCvtxx->Fill(vtxMC[0]);   
      fHistMCvtxy->Fill(vtxMC[1]);
      fHistMCvtxz->Fill(vtxMC[2]);

//      Printf("Impact parameter gen: %f", impactParameter);
       if (hijingHeader) npart = hijingHeader->TargetParticipants()+hijingHeader->ProjectileParticipants();
       else if (dpmHeader)npart = dpmHeader->TargetParticipants()+dpmHeader->ProjectileParticipants();

      //Rec centrality vs gen centrality
      AliESDZDC *zdcRec = fmyESD->GetESDZDC();
      Double_t impactParameterZDC = zdcRec->GetImpactParameter();
//      Printf("Impact parameter rec: %f", impactParameterZDC);

      fHistRecvsGenImpactPar->Fill(impactParameter,impactParameterZDC);
      fHistMCNpart->Fill(npart);
  
      // Tracks from MC
      Int_t  multMCCharged = 0;
      Int_t  multMCChargedEtacut = 0;
//      Int_t  nMCPart = stack->GetNprimary();
      Int_t  nMCPart = stack->GetNtrack();  // decay products of D and B mesons are also primaries and produced in HIJING during transport
      Float_t* etagen = new Float_t[nMCPart];  
      Int_t* stackIndexOfPrimaryParts = new Int_t[nMCPart];
      Bool_t* reconstructedPrimaryPart = new Bool_t[nMCPart];
      Bool_t* detectedPrimaryPartLay1 = new Bool_t[nMCPart];
      Bool_t* detectedPrimaryPartLay2 = new Bool_t[nMCPart];
      Bool_t* detectablePrimaryPart = new Bool_t[nMCPart];

      TTree* tRefTree;  
      if (fTR) {
        tRefTree = eventHandler->TreeTR(); 
        fMCEvent->ConnectTreeTR(tRefTree); 
      }

      // Loop over MC particles
      for (Int_t imc=0; imc<nMCPart; imc++) {
        AliMCParticle *mcpart  = (AliMCParticle*)fMCEvent->GetTrack(imc);

        Bool_t isPrimary = stack->IsPhysicalPrimary(imc);
        if (!isPrimary)            continue;
        if (mcpart->Charge() == 0) continue;
        Float_t theta = mcpart->Theta();
        if (theta==0 || theta==TMath::Pi())    continue;
        Float_t eta = mcpart->Eta();
//        Float_t pt = mcpart->Pt();
        etagen[multMCCharged] = eta; 
        stackIndexOfPrimaryParts[multMCCharged] = imc;

        reconstructedPrimaryPart[multMCCharged]=kFALSE;
        detectedPrimaryPartLay1[multMCCharged]=kFALSE;
        detectedPrimaryPartLay2[multMCCharged]=kFALSE;
        detectablePrimaryPart[multMCCharged]=kFALSE;

        if (fTR) {
          Int_t nref = mcpart->GetNumberOfTrackReferences();
          if (nref==0) {
            detectablePrimaryPart[multMCCharged]= kTRUE;
          } else if (nref>0) {
            // Check if the primary is detectable 
            detectablePrimaryPart[multMCCharged] = IsDetectablePrimary(nref,mcpart);
            // Check if the primary is detected (if SPD 100% efficient) 
            if (detectablePrimaryPart[multMCCharged]) {
              detectedPrimaryPartLay1[multMCCharged] = IsDetectedPrimary(nref,mcpart,0);
              detectedPrimaryPartLay2[multMCCharged] = IsDetectedPrimary(nref,mcpart,1);
            }
          }  
        }
        if (eventSelected&&fPartSpecies) {
          if (TMath::Abs(mcpart->PdgCode())==2212) fHistProtons->Fill(etagen[multMCCharged],vtxMC[2]);
          else if (TMath::Abs(mcpart->PdgCode())==321) fHistKaons->Fill(etagen[multMCCharged],vtxMC[2]);
          else if (TMath::Abs(mcpart->PdgCode())==211) fHistPions->Fill(etagen[multMCCharged],vtxMC[2]);
          else fHistOthers->Fill(etagen[multMCCharged],vtxMC[2]);  //includes leptons pdg->11,13
        }
        multMCCharged++;
        if (TMath::Abs(eta)<1.4) multMCChargedEtacut++;
      } // end of MC particle loop

      fHistMCmultEtacut->Fill(multMCChargedEtacut);

      // Event selection: in centrality bin, triggered with vertex 
      if (eventTriggered&&eventWithVertex&&eventInCentralityBin&&multSPDcl2>fMinClMultLay2&&multSPDcl2<fMaxClMultLay2) {
    
        fHistDeVtx->Fill(vtxMC[2],vtxMC[2]-esdvtx[2]);      

        for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
          if (TMath::Abs(trackletCoord[itracklet][0])<fPhiWindowAna) 
            fHistBkgCorrDen->Fill(trackletCoord[itracklet][4],esdvtx[2]);

          if (trackletLab[itracklet][0]==trackletLab[itracklet][1]) {
            Bool_t trakletByPrim = kFALSE; 
            for (Int_t imc=0; imc<multMCCharged; imc++) {
              if (trackletLab[itracklet][0]==stackIndexOfPrimaryParts[imc]) {
                fHistdPhidThetaPP->Fill(trackletCoord[itracklet][0],trackletCoord[itracklet][1]);
                if (TMath::Abs(trackletCoord[itracklet][0])<fPhiWindowAna) {
                  fHistBkgCorrDenPrimGen->Fill(etagen[imc],vtxMC[2]);
                  if (fPartSpecies) {
                    if (TMath::Abs(((AliMCParticle*)fMCEvent->GetTrack(stackIndexOfPrimaryParts[imc]))->PdgCode())==2212) 
                      fHistReconstructedProtons->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else if (TMath::Abs(((AliMCParticle*)fMCEvent->GetTrack(stackIndexOfPrimaryParts[imc]))->PdgCode())==321) 
                      fHistReconstructedKaons->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else if (TMath::Abs(((AliMCParticle*)fMCEvent->GetTrack(stackIndexOfPrimaryParts[imc]))->PdgCode())==211) 
                      fHistReconstructedPions->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else fHistReconstructedOthers->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                  }     
                }
                trakletByPrim = kTRUE; 
                if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) {
                  reconstructedPrimaryPart[imc]=kTRUE; // tracklet by prim and tr (=rp) on both layers 
                }
                break;
              }  
            }
            if (!trakletByPrim) {
              fHistdPhidThetaSS->Fill(trackletCoord[itracklet][0],trackletCoord[itracklet][1]);
              if (TMath::Abs(trackletCoord[itracklet][0])<fPhiWindowAna) {

                fHistBkgCorrDenPrimGen->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                if (fPartSpecies) {
                  Int_t motherlab = ((AliMCParticle*)fMCEvent->GetTrack((Int_t)trackletLab[itracklet][0]))->GetMother(); 
                  if (motherlab>-1) { 
                    if (TMath::Abs(fMCEvent->GetTrack(motherlab)->PdgCode())==2212) fHistReconstructedProtons->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else if (TMath::Abs(fMCEvent->GetTrack(motherlab)->PdgCode())==321) fHistReconstructedKaons->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else if (TMath::Abs(fMCEvent->GetTrack(motherlab)->PdgCode())==211) fHistReconstructedPions->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                    else fHistReconstructedOthers->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                  } else fHistReconstructedSec->Fill(trackletCoord[itracklet][4],esdvtx[2]);
                }
              }
            }
          } else {
            fHistdPhidThetaComb->Fill(trackletCoord[itracklet][0],trackletCoord[itracklet][1]);
            if (TMath::Abs(trackletCoord[itracklet][0])<fPhiWindowAna) {
              fHistBkgCorrDenPrimGen->Fill(trackletCoord[itracklet][4],esdvtx[2]);
              fHistBkgCombLabels->Fill(trackletCoord[itracklet][4],esdvtx[2]); 
            }
          }
        } // end loop tracklets

        for (Int_t imc=0; imc<multMCCharged; imc++) {
          if (fTR) {
            if (reconstructedPrimaryPart[imc]) fHistBkgCorrNum->Fill(etagen[imc],vtxMC[2]); // only for separate corrections
            if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) fHistAlgEffNum->Fill(etagen[imc],vtxMC[2]);
            if (detectablePrimaryPart[imc]) fHistNonDetectableCorrDen->Fill(etagen[imc],vtxMC[2]); 
          }
          fHistNonDetectableCorrNum->Fill(etagen[imc],vtxMC[2]);
        }

        fHistSelEvts->Fill(multSPDAna,vtxMC[2]);
      } // end of selected events

      fHistMCmultEtacutvsSPDmultEtacut->Fill(multMCChargedEtacut,multSPDEtacut);
      fHistMCmultEtacutvsSPDmultcl1->Fill(multMCChargedEtacut,multSPDcl1);
      fHistMCmultEtacutvsSPDmultcl2->Fill(multMCChargedEtacut,multSPDcl2);
      fHistAllEvts->Fill(multSPDAna,vtxMC[2]);

      if (eventInCentralityBin) {
        fHistCentrEvts->Fill(multSPDAna,vtxMC[2]);
        if (eventTriggered) {
          fHistTrigCentrEvts->Fill(multSPDAna,vtxMC[2]);
        }
      }

      for (Int_t imc=0; imc<multMCCharged; imc++) { 
        fHistAllPrimaries->Fill(etagen[imc],vtxMC[2]);
        if (eventInCentralityBin) {
          fHistTrackCentrEvts->Fill(etagen[imc],vtxMC[2]);
          if (eventTriggered) fHistTrackTrigCentrEvts->Fill(etagen[imc],vtxMC[2]); 
        }
      } 
    
      delete[] etagen;
      delete[] stackIndexOfPrimaryParts;
      delete[] reconstructedPrimaryPart;
      delete[] detectedPrimaryPartLay1;
      delete[] detectedPrimaryPartLay2;
      delete[] detectablePrimaryPart;
    }  
  }
  PostData(1, fOutput);
}      
//________________________________________________________________________
Bool_t AliAnalysisTaskSPDdNdEta::IsDetectedPrimary(Int_t nref, AliMCParticle* mcpart, Int_t Layer) {
  Double_t rMinLay1 = 3.4; //min rad for track ref first SPD layer 
  Double_t rMaxLay1 = 4.4; //max rad for track ref first SPD layer 
  Double_t rMinLay2 = 6.9; //min rad for track ref second SPD layer  
  Double_t rMaxLay2 = 7.9; //max rad for track ref second SPD layer 
  Double_t rmin = (Layer > 0 ? rMinLay2 : rMinLay1);
  Double_t rmax = (Layer > 0 ? rMaxLay2 : rMaxLay1);
  AliTrackReference *tref=0x0;

  for (Int_t iref=0;iref<nref;iref++) {
     tref = mcpart->GetTrackReference(iref);
   if (tref) {
       if (tref->R()>rmin&&tref->R()<rmax) {
         if (tref->DetectorId()==0) {
           return kTRUE;
         }
       }
     }
  }
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSPDdNdEta::IsDetectablePrimary(Int_t nref, AliMCParticle* mcpart) {

  Double_t rMinLay2 = 6.9; //min rad for track ref second SPD layer  
  Double_t rMaxLay2 = 7.9; //max rad for track ref second SPD layer 

  AliTrackReference *tref= mcpart->GetTrackReference(nref-1);
  if (tref->DetectorId()!=-1) {
    return kTRUE;
  } else { //last is -1 -> particle disappeared. Where?
    if (tref->R()>rMaxLay2) {
      return kTRUE;  
    } else if (tref->R()>=rMinLay2&&tref->R()<=rMaxLay2) { // this last tr is in lay 2
      for (Int_t iref=0;iref<nref;iref++) {  // look for other tr in lay 2
        tref = mcpart->GetTrackReference(iref);
        if (tref) 
          if (tref->R()>=rMinLay2&&tref->R()<=rMaxLay2) 
            if (tref->DetectorId()==0) return kTRUE; 
      } 
    } // else particle disappeared before lay2
  }
  return kFALSE; 
}
//________________________________________________________________________
void AliAnalysisTaskSPDdNdEta::Terminate(Option_t *) 
{
  Printf("Terminating...");
//  AliAnalysisTaskSE::Terminate();

}
