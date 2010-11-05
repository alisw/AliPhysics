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

#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"

#include "AliESDInputHandlerRP.h"
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

#include "AliLog.h"

#include "AliTriggerAnalysis.h" 
#include "AliPhysicsSelection.h"
#include "AliESDCentrality.h" 
#include "AliTrackletAlg.h" 
#include "AliAnalysisTaskSPDdNdEta.h"


ClassImp(AliAnalysisTaskSPDdNdEta)
//________________________________________________________________________
AliAnalysisTaskSPDdNdEta::AliAnalysisTaskSPDdNdEta(const char *name) 
  : AliAnalysisTaskSE(name), 

  fmyESD(0), 
  fOutput(0), 

  fMCCentralityBin(AliAnalysisTaskSPDdNdEta::kall),
  fCentrLowLim(0),
  fCentrUpLim(0),
  fCentrEst(""),
  
  fUseMC(kFALSE), 
  fPbPb(kFALSE),
  fTrigger(AliTriggerAnalysis::kAcceptAll),
  fTR(kFALSE),
  fRecoTracklets(kFALSE),

  fHistSPDRAWMultvsZ(0),
  fHistSPDRAWMultvsZTriggCentrEvts(0),
  fHistSPDRAWMultvsZCentrEvts(0),
  fHistSPDRAWEtavsZ(0),
 
  fHistSPDmultEtacut(0),
  fHistSPDmult(0),
  fHistSPDmultcl1(0),
  fHistSPDeta(0),
  fHistSPDphi(0),
  fHistSPDtheta(0),
  fHistSPDdePhi(0),
  fHistSPDphivsSPDeta(0),  
  fHistSPDdeTheta(0),
  fHistSPDvtxAnalysis(0),
  fHistSPDdePhideTheta(0),

  fHistBkgCorrDen(0),
  fHistBkgCorrDenPrimGen(0),
  fHistBkgCorrNum(0),
  fHistAlgEffNum(0),
  fHistNonDetectableCorrDen(0),

  fHistNonDetectableCorrNum(0),
  fHistAllPrimaries(0),
  fHistTrackCentrEvts(0),
  fHistTrackTrigCentrEvts(0),

  fHistAllEvts(0),
  fHistCentrEvts(0),
  fHistTrigCentrEvts(0),
  fHistSelEvts(0),

  fHistMCmultEtacut(0),
  fHistMCmultEtacutvsSPDmultEtacut(0),
 
  fHistMCvtxx(0),
  fHistMCvtxy(0),
  fHistMCvtxz(0),

  fHistRecvsGenImpactPar(0),
  fHistMCNpart(0), 

  fHistdPhiPP(0),
  fHistdPhiSS(0),
  fHistdPhiComb(0),

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
//  man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
    man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
//  man->SetSpecificStorage("GRP/Geometry/Data","local://");
//  man->SetSpecificStorage("ITS/Align/Data","local://");
    man->SetRun(130844); 
    AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
////  AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130845,7);
////  AliCDBEntry* obj = man->Get(AliCDBPath("ITS", "Align", "Data"));
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    AliGeomManager::GetNalignable("ITS");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
////  AliGeomManager::ApplyAlignObjsToGeom("ITS",130845,5,-1);
  }
  // Create histograms
  fOutput = new TList();
  fOutput->SetOwner(); 

  Int_t nBinMultCorr = 200;
  Float_t nMaxMult = 20000.;
  Double_t binLimMultCorr[nBinMultCorr+1];
  binLimMultCorr[0] = 0.;
  binLimMultCorr[1] = 1.;
  for (Int_t i = 2; i<=nBinMultCorr;++i) {
    binLimMultCorr[i] = (i-1)*nMaxMult/nBinMultCorr;
  }

  if(!fPbPb) { 
    for (Int_t i = 0; i<nBinMultCorr+1; ++i) {
      if (i<21) binLimMultCorr[i] = (Double_t) i;
      else if (i>=21&&i<27) binLimMultCorr[i] = (Double_t) 20 +5*(i-20);
      else if (i==27) binLimMultCorr[i] = 100;
      else if (i==28) binLimMultCorr[i] = 200; 
      else if (i==29) binLimMultCorr[i] = 400;
      else if (i==30) binLimMultCorr[i] = 600;
    }
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
  fHistSPDRAWMultvsZ = new TH2F("fHistSPDRAWMultvsZ", "",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  fHistSPDRAWMultvsZ->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWMultvsZ->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZ); 

  fHistSPDRAWMultvsZTriggCentrEvts = new TH2F("fHistSPDRAWMultvsZTriggCentrEvts", "",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  fHistSPDRAWMultvsZTriggCentrEvts->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZTriggCentrEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWMultvsZTriggCentrEvts->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZTriggCentrEvts);

  fHistSPDRAWMultvsZCentrEvts = new TH2F("fHistSPDRAWMultvsZCentrEvts", "",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  fHistSPDRAWMultvsZCentrEvts->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDRAWMultvsZCentrEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWMultvsZCentrEvts->Sumw2();
  fOutput->Add(fHistSPDRAWMultvsZCentrEvts);

  // ...track level
  fHistSPDRAWEtavsZ = new TH2F("fHistSPDRAWEtavsZ", "Tracklet pseudorapidity distribution", nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  fHistSPDRAWEtavsZ->GetXaxis()->SetTitle("Pseudorapidity #eta");
  fHistSPDRAWEtavsZ->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDRAWEtavsZ->Sumw2();
  fOutput->Add(fHistSPDRAWEtavsZ);
  

  fHistSPDmultEtacut = new TH1F("fHistSPDmultEtacut", "Tracklet multiplicity distribution",nBinMultCorr,binLimMultCorr);
  fHistSPDmultEtacut->GetXaxis()->SetTitle("Tracklet multiplicity (|#eta|<1.4)");
  fHistSPDmultEtacut->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultEtacut);

  fHistSPDmult = new TH1F("fHistSPDmult", "Tracklet multiplicity distribution", nBinMultCorr,binLimMultCorr);
  fHistSPDmult->GetXaxis()->SetTitle("Tracklet multiplicity");
  fHistSPDmult->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmult);

  fHistSPDmultcl1 = new TH1F("fHistSPDmultcl1", "Inner layer cluster multiplicity distribution", nBinMultCorr,binLimMultCorr);
  fHistSPDmultcl1->GetXaxis()->SetTitle("Inner layer cluster multiplicity");
  fHistSPDmultcl1->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDmultcl1);

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

  fHistSPDvtxAnalysis = new TH1F("fHistSPDvtxAnalysis", "SPD vertex distribution: events selected for the analysis",20,-20.,20.);
  fHistSPDvtxAnalysis->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistSPDvtxAnalysis->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fHistSPDvtxAnalysis);

  fHistSPDdePhideTheta= new TH2F("fHistSPDdePhideTheta", "Tracklet #Delta#varphi  distribution",2000,-1.,1.,1000,-0.25,.25);
  fHistSPDdePhideTheta->GetXaxis()->SetTitle("#Delta#varphi [rad]");
  fHistSPDdePhideTheta->GetYaxis()->SetTitle("#Delta#theta [rad]");
  fOutput->Add(fHistSPDdePhideTheta);


  if (fUseMC) {

    // Track level correction histograms 
    fHistBkgCorrDen = new TH2F("fHistBkgCorrDen","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
    fOutput->Add(fHistBkgCorrDen);
    
    fHistBkgCorrDenPrimGen = new TH2F("fHistBkgCorrDenPrimGen","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
    fOutput->Add(fHistBkgCorrDenPrimGen);

    if (fTR) {
      fHistBkgCorrNum = new TH2F("fHistBkgCorrNum","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
      fOutput->Add(fHistBkgCorrNum);
 
      fHistAlgEffNum = new TH2F("fHistAlgEffNum","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
      fOutput->Add(fHistAlgEffNum);  

      fHistNonDetectableCorrDen = new TH2F("fHistNonDetectableCorrDen","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
      fOutput->Add(fHistNonDetectableCorrDen);

    }

    fHistNonDetectableCorrNum = new TH2F("fHistNonDetectableCorrNum","",nBinEtaCorr,binLimEtaCorr,20,-20.,
20.);
    fOutput->Add(fHistNonDetectableCorrNum);
    
    fHistAllPrimaries = new TH2F("fHistAllPrimaries","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
    fOutput->Add(fHistAllPrimaries);

    fHistTrackCentrEvts = new TH2F("fHistTrackCentrEvts","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
    fOutput->Add(fHistTrackCentrEvts);

    fHistTrackTrigCentrEvts = new TH2F("fHistTrackTrigCentrEvts","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
    fOutput->Add(fHistTrackTrigCentrEvts);


    // Event level correction histograms  
    fHistAllEvts = new TH2F("fHistAllEvts","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistAllEvts);

    fHistCentrEvts = new TH2F("fHistCentrEvts","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistCentrEvts);

    fHistTrigCentrEvts = new TH2F("fHistTrigCentrEvts","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
    fOutput->Add(fHistTrigCentrEvts);

    fHistSelEvts = new TH2F("fHistSelEvts","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
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
 
    fHistdPhiPP = new TH1F("fHistdPhiPP","",400,-0.1,.1);
    fOutput->Add(fHistdPhiPP);
    fHistdPhiSS = new TH1F("fHistdPhiSS","",400,-0.1,.1);
    fOutput->Add(fHistdPhiSS);
    fHistdPhiComb = new TH1F("fHistdPhiComb","",400,-0.1,.1);
    fOutput->Add(fHistdPhiComb);

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

  fmyESD = dynamic_cast<AliESDEvent*>(InputEvent()); 
  if (!fmyESD) {
    Printf("ERROR: fESD not available");
    return;
  }


  // Trigger selection
  Bool_t eventTriggered = kTRUE;
  static AliTriggerAnalysis* triggerAnalysis = 0; 
  if (eventTriggered) // applying an offline trigger
    eventTriggered = triggerAnalysis->IsTriggerFired(fmyESD, fTrigger);
  if (!eventTriggered)
    Printf("No trigger");


  // Centrality selection 
  Bool_t eventInCentralityBin = kFALSE;
  // Centrality selection
  AliESDCentrality *centrality = fmyESD->GetCentrality();
  if (fCentrEst=="") eventInCentralityBin = kTRUE;
  else {
    if(!centrality) {
      AliError("Centrality object not available"); 
    }  else {
      if (centrality->IsEventInCentralityClass(fCentrLowLim,fCentrUpLim,fCentrEst.Data())) eventInCentralityBin = kTRUE;
    }
  }


  // ESD vertex
  Bool_t eventWithVertex = kFALSE;
  const AliESDVertex* vtxESD = fmyESD->GetVertex();
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  //... check resolution
//  Double_t zRes = vtxESD->GetZRes();
//  Double_t zDisp = vtxESD->GetDispersion();
  if (esdvtx[2]!=0) eventWithVertex = kTRUE; //vertex selection in PbPb  


  // Reconstructing or loading tracklets...
  const AliMultiplicity* multESD = fmyESD->GetMultiplicity();
  Int_t multSPD = multESD->GetNumberOfTracklets();
  Int_t nSingleCl1 = multESD->GetNumberOfSingleClusters();
  Int_t multSPDcl1 = nSingleCl1 + multSPD;  
  Int_t multSPDEtacut = 0;
  Int_t multSPDAna = 0;

  fHistSPDmultcl1->Fill(multSPDcl1);

  Float_t trackletCoord[multSPDcl1][5];
  Float_t trackletLab[multSPDcl1][2];

  // Selected events: in centrality bin, triggered with vertex
  if (eventTriggered&&eventWithVertex&&eventInCentralityBin) {

    if (fRecoTracklets) {
      // Load ITS rec points tree
      TTree *itsClusterTree = hand->GetTreeR("ITS");
      if (!itsClusterTree) {
        AliError(" Invalid ITS cluster tree !\n");
        return;
      }
      float vtxf[3] = {vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ()};

      AliTrackletAlg *fMultReco = new AliTrackletAlg();

      fMultReco->SetPhiWindow(0.8);
//      fMultReco->SetPhiWindow(.08); //def
      fMultReco->SetThetaWindow(0.025); //def
//      fMultReco->SetThetaWindow(0.25);
      fMultReco->SetPhiShift(0.0045); //def
      fMultReco->SetRemoveClustersFromOverlaps(kFALSE); // def
      fMultReco->SetPhiOverlapCut(0.005); //def
      fMultReco->SetZetaOverlapCut(0.05); //def
      fMultReco->SetPhiRotationAngle(0.0); // def
//      fMultReco->SetPhiRotationAngle(TMath::Pi());
      fMultReco->SetHistOn(kFALSE); //def
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
        trackletCoord[itr][4] = -TMath::Log(TMath::Tan(trackletCoord[itr][3]/2.));
        trackletLab[itr][0] = static_cast<Int_t>(fMultReco->GetTracklet(itr)[4]);
        trackletLab[itr][1] = static_cast<Int_t>(fMultReco->GetTracklet(itr)[5]);
      }
      delete  fMultReco;
    } else {
      // set variables from ESD
      for (Int_t itr = 0; itr<multSPD;++itr) {
        trackletCoord[itr][3] = multESD->GetTheta(itr);       //theta
        trackletCoord[itr][2] = multESD->GetPhi(itr);         //phi
        trackletCoord[itr][1] = multESD->GetDeltaTheta(itr);  //deTheta
        trackletCoord[itr][0] = multESD->GetDeltaPhi(itr);    //dePhi
        trackletCoord[itr][4] = multESD->GetEta(itr);         //Eta 
        trackletLab[itr][0] = multESD->GetLabel(itr,0);       //label lay1
        trackletLab[itr][1] = multESD->GetLabel(itr,1);       //label lay2
      }
    }
//    Printf("tracklets in ESD...%d",multESD->GetNumberOfTracklets());

    for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
      if (TMath::Abs(trackletCoord[itracklet][0])<0.08) { // select tracklets with tighter cuts
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

    AliGenHijingEventHeader* hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(fMCEvent->Header()->GenEventHeader());
    if (!hijingHeader) {
      Printf("Unknown header type. \n");
      return ;
    }
    Float_t impactParameter = hijingHeader->ImpactParameter();
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
//    Printf("Impact parameter gen: %f", impactParameter);
    Int_t npart = hijingHeader->TargetParticipants()+hijingHeader->ProjectileParticipants();
    //Rec centrality vs gen centrality
    // Centrality (reconstructed)
    AliESDZDC *zdcRec = fmyESD->GetESDZDC();
    Double_t impactParameterZDC = zdcRec->GetImpactParameter();
//    Printf("Impact parameter rec: %f", impactParameterZDC);

    fHistRecvsGenImpactPar->Fill(impactParameter,impactParameterZDC);
    fHistMCNpart->Fill(npart);
    // Tracks from MC
    Int_t  multMCCharged = 0;
    Int_t  multMCChargedEtacut = 0;
    Int_t  nMCPart = stack->GetNprimary();
    Float_t* etagen = new Float_t[nMCPart];  
    Int_t* stackIndexOfPrimaryParts = new Int_t[nMCPart];
    Bool_t* reconstructedPrimaryPart = new Bool_t[nMCPart];
    Bool_t* detectedPrimaryPartLay1 = new Bool_t[nMCPart];
    Bool_t* detectedPrimaryPartLay2 = new Bool_t[nMCPart];
    Bool_t* detectablePrimaryPart = new Bool_t[nMCPart];
    Bool_t* primCounted = new Bool_t[nMCPart];

    TTree* tRefTree;  
    if (fTR) {
      tRefTree = eventHandler->TreeTR(); 
      fMCEvent->ConnectTreeTR(tRefTree); 
    }

    // Loop over MC particles
    for (Int_t imc=0; imc<nMCPart; imc++) {
      AliMCParticle *mcpart  = (AliMCParticle*)fMCEvent->GetTrack(imc);
      Bool_t isPrimary = stack->IsPhysicalPrimary(imc);
      if (!isPrimary)                        continue;
      if (mcpart->Charge() == 0) continue;
      Float_t theta = mcpart->Theta();
      if (theta==0 || theta==TMath::Pi())    continue;
      Float_t eta = mcpart->Eta();
//      Float_t pt = mcpart->Pt();
      etagen[multMCCharged] = eta; 
      stackIndexOfPrimaryParts[multMCCharged] = imc;

      reconstructedPrimaryPart[multMCCharged]=kFALSE;
      detectedPrimaryPartLay1[multMCCharged]=kFALSE;
      detectedPrimaryPartLay2[multMCCharged]=kFALSE;
      detectablePrimaryPart[multMCCharged]=kFALSE;
      primCounted[multMCCharged]=kFALSE;

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
      multMCCharged++;
      if (TMath::Abs(eta)<1.4) multMCChargedEtacut++;
    } // end of MC particle loop
    fHistMCmultEtacut->Fill(multMCChargedEtacut);

    // Event selection
    if (eventTriggered&&eventWithVertex&&eventInCentralityBin) {
    
      fHistDeVtx->Fill(vtxMC[2],vtxMC[2]-esdvtx[2]);      

      for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
        if (TMath::Abs(trackletCoord[itracklet][0])<0.08) {
          fHistBkgCorrDen->Fill(trackletCoord[itracklet][4],esdvtx[2]);

          if (trackletLab[itracklet][0]==trackletLab[itracklet][1]) {
            Bool_t trakletByPrim = kFALSE; 
            for (Int_t imc=0; imc<multMCCharged; imc++) {
              if (trackletLab[itracklet][0]==stackIndexOfPrimaryParts[imc]) {
                if (!primCounted[imc]) {
                  primCounted[imc] = kTRUE;
                }
                fHistdPhiPP->Fill(trackletCoord[itracklet][0]);
                fHistBkgCorrDenPrimGen->Fill(etagen[imc],vtxMC[2]);
                trakletByPrim = kTRUE; 
                if (detectedPrimaryPartLay1[imc]&&detectedPrimaryPartLay2[imc]) {
                  reconstructedPrimaryPart[imc]=kTRUE; // tracklet by prim and tr (=rp) on both layers 
                }
                break;
              }  
            }
            if (!trakletByPrim) {
              fHistdPhiSS->Fill(trackletCoord[itracklet][0]);
              fHistBkgCorrDenPrimGen->Fill(trackletCoord[itracklet][4],esdvtx[2]);
            }
          } else {
            fHistdPhiComb->Fill(trackletCoord[itracklet][0]);
            fHistBkgCorrDenPrimGen->Fill(trackletCoord[itracklet][4],esdvtx[2]);
          }
        }  
      }

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
Bool_t AliAnalysisTaskSPDdNdEta::IsDetectedPrimary(Int_t nref, AliMCParticle* mcpart, Int_t Layer)
{
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
Bool_t AliAnalysisTaskSPDdNdEta::IsDetectablePrimary(Int_t nref, AliMCParticle* mcpart)
{
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
