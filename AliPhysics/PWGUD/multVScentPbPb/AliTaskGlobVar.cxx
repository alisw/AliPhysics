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

///////////////////////////////////////////////////////////////////////////
// Analysis task to extract global variables to the tree                 //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TList.h>
#include <TObjArray.h>
#include "AliAnalysisManager.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliPhysicsSelection.h"
#include "AliCentrality.h" 
#include "AliESDtrackCuts.h"
#include "AliTaskGlobVar.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliTriggerAnalysis.h" 


ClassImp(AliTaskGlobVar)


//________________________________________________________________________
AliTaskGlobVar::AliTaskGlobVar(const char *name) 
  : AliAnalysisTaskSE(name), 
    //
    fUseMC(kFALSE),
    fOutput(0), 
    fOutTree(0),
    fTrackCuts(0),
  fTrackCuts1(0),
  fGlobVars()
{
  // Constructor
  DefineOutput(1, TList::Class());
  //
}

//________________________________________________________________________
AliTaskGlobVar::~AliTaskGlobVar()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
  }
  //
}

//________________________________________________________________________
void AliTaskGlobVar::UserCreateOutputObjects() 
{
  // create ouptut
  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  fOutTree = new TTree("globTree","globTree");
  fOutTree->Branch("g",&fGlobVars,"runID/I:timeStamp/i:zdcNA/F:zdcPA/F:zdcNC/F:zdcPC/F:zdcNAC/F:zdcNCC/F:zem1/F:zem2/F:zvSPD/F:zvTPC/F:chunk/S:flags/S"
		   ":spd1/S:spd2/S:ncontSPDV/S:ncontTPCV/S:nTrTPC/S:nTrTPCITS/S:nTracklets/S:v0A/S:v0C/S:v0Corr/S", 16777216);

  if (fUseMC) {
    fOutTree->Branch("gmc",&fGlobVars.mcZV,"mcZV/F:mcdNdEta/S:mcNPart/S:mcNBColl/S", 16777216);
  }
  fOutput->Add(fOutTree);
  //
  //  fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTrackCuts  = CreatedNdPtTrackCuts(23);
  fTrackCuts->SetEtaRange(-0.8,0.8);
  fTrackCuts1 = CreatedNdPtTrackCuts(72);
  fTrackCuts1->SetEtaRange(-0.8,0.8);
  //
  PostData(1, fOutput);
  printf("CreatingXXX\n");
  //
}

//________________________________________________________________________
void AliTaskGlobVar::UserExec(Option_t *) 
{
  // Main loop
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  AliESDInputHandler *hand = (AliESDInputHandler*)anMan->GetInputEventHandler();
  if (!hand) { printf("No ESD handler\n"); return; }
  AliESDEvent *esd  = hand->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  //
  if (!fUseMC && !ZDCTimeTrigger(esd)) return;
  //
  const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
  Bool_t vtxOK = kTRUE;
  if (vtxESD->GetNContributors()<1) vtxOK = kFALSE; //return;
  if (vtxESD->GetDispersion()>0.04) vtxOK = kFALSE; //return;
  if (vtxESD->GetZRes()>0.25) vtxOK = kFALSE; //return;
  const AliMultiplicity* multESD  = esd->GetMultiplicity();
  const AliESDVertex*   vtxESDTPC = esd->GetPrimaryVertexTPC();
  //  if (vtxESDTPC->GetNContributors()<1) return;
  //  if (vtxESDTPC->GetNContributors()<(-10.+0.25*multESD->GetNumberOfITSClusters(0))) return;
  //
  fGlobVars.runID = esd->GetRunNumber();
  TString rid = "";
  rid +=  fGlobVars.runID;
  TString flname        = hand->GetTree()->GetCurrentFile()->GetName();
  fGlobVars.chunk        = 0;
  int id = 0;
  while ( (id=flname.Index(rid))>=0 ) flname = flname.Data()+id+rid.Length();
  id = flname.First('.');
  if (id>=0) {
    flname = flname.Data() + id+1;
    fGlobVars.chunk = (Short_t)flname.Atoi();
  }
  //  printf("%d %s\n",fGlobVars.chunk,hand->GetTree()->GetCurrentFile()->GetName());
  //
  fGlobVars.timeStamp    = esd->GetTimeStamp();
  fGlobVars.timeStamp    = esd->GetTimeStamp();
  fGlobVars.zvSPD        = vtxESD->GetZ(); 
  fGlobVars.zvTPC        = vtxESDTPC->GetZ(); 
  fGlobVars.ncontSPDV = (Short_t)vtxESD->GetNContributors();
  fGlobVars.ncontTPCV = (Short_t)vtxESDTPC->GetNContributors();
  //  fGlobVars.nTrTPC     = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(esd,kTRUE):-1;
  fGlobVars.nTrTPC       = fTrackCuts  ? (Short_t)fTrackCuts->CountAcceptedTracks(esd):-1;
  fGlobVars.nTrTPCITS    = fTrackCuts1 ? (Short_t)fTrackCuts1->CountAcceptedTracks(esd):-1;
  fGlobVars.nTracklets = multESD->GetNumberOfTracklets();
  //
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  if (esdV0) {
    fGlobVars.v0A = (Short_t)esdV0->GetMTotV0A();
    fGlobVars.v0C = (Short_t)esdV0->GetMTotV0C();
  }
  else fGlobVars.v0A = fGlobVars.v0C = 0;
  //
  fGlobVars.spd1 = (Short_t)multESD->GetNumberOfITSClusters(0);
  fGlobVars.spd2 = (Short_t)multESD->GetNumberOfITSClusters(1);
  //
  float v0Corr,v0CorrR;
  //  v0Corr = GetCorrV0(esd,v0CorrR); MF: Deprecated: V0corr isnot needed any more
  v0Corr = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
  fGlobVars.v0Corr = (Short_t)v0Corr;
  //  fGlobVars.v0CorrResc = (Short_t)v0CorrR;

  //---------------------------------------------
  AliESDZDC *esdZDC = esd->GetESDZDC();
  // --- ZDC offline trigger ---
  // Returns if ZDC triggered, based on TDC information
  Bool_t tdc[32] = {kFALSE};
  for(Int_t itdc=0; itdc<32; itdc++){
    for(Int_t i=0; i<4; i++){
      if (0.025*esdZDC->GetZDCTDCData(itdc, i) != 0){
	tdc[itdc] = kTRUE;
      }
    }
  }
  fGlobVars.flags = 0;
  if ( tdc[12] ) fGlobVars.flags |= GloVars_t::kTDCNA; //  Bool_t zdcNA = tdc[12];
  if ( tdc[10] ) fGlobVars.flags |= GloVars_t::kTDCNC; //  Bool_t zdcNC = tdc[10];
  if ( tdc[13] ) fGlobVars.flags |= GloVars_t::kTDCPA; //  Bool_t zdcPA = tdc[13];
  if ( tdc[11] ) fGlobVars.flags |= GloVars_t::kTDCPC; //  Bool_t zdcPC = tdc[11];
  if ( vtxOK   ) fGlobVars.flags |= GloVars_t::kSPDVTXOK;
  //
  
  const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
  const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
  const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
  const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();

  fGlobVars.zdcNC = (Float_t) (esdZDC->GetZDCN1Energy());
  fGlobVars.zdcPC = (Float_t) (esdZDC->GetZDCP1Energy());
  fGlobVars.zdcNA = (Float_t) (esdZDC->GetZDCN2Energy());
  fGlobVars.zdcPA = (Float_t) (esdZDC->GetZDCP2Energy());

  fGlobVars.zdcNCC = (Float_t) (towZNC[0]);
  fGlobVars.zdcNAC = (Float_t) (towZNA[0]);


  fGlobVars.zem1  = (Float_t) (esdZDC->GetZDCEMEnergy(0)) /8.;
  fGlobVars.zem2  = (Float_t) (esdZDC->GetZDCEMEnergy(1)) /8.;



  //-----------------------------------------------
  //
  // ---------------------- MC ONLY -------------------------------
  AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
  AliStack*    stack=0;
  AliMCEvent*  mcEvent=0;
  fGlobVars.mcdNdEta = 0;
  fGlobVars.mcNPart  = 0;
  fGlobVars.mcNBColl = 0;
  if (fUseMC && eventHandler && (mcEvent=eventHandler->MCEvent()) && (stack=mcEvent->Stack())) {
    int ntr = stack->GetNtrack();
    for (int itr=ntr;itr--;) {
      if (!stack->IsPhysicalPrimary(itr)) continue;
      AliMCParticle *part  = (AliMCParticle*)mcEvent->GetTrack(itr);
      if (!part || !part->Charge()) continue;
      Float_t eta = part->Eta();
      if (TMath::Abs(eta)>0.5) continue;
      fGlobVars.mcdNdEta++;
    }
    //
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())) {
      AliGenHijingEventHeader* hHijing = (AliGenHijingEventHeader*)mcGenH;
      fGlobVars.mcNPart  = (hHijing->ProjectileParticipants()+hHijing->TargetParticipants())/2.;
      fGlobVars.mcNBColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
    }
    else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
      AliGenDPMjetEventHeader* hDpmJet = (AliGenDPMjetEventHeader*)mcGenH;
      fGlobVars.mcNPart  = (hDpmJet->ProjectileParticipants()+hDpmJet->TargetParticipants())/2.;
      fGlobVars.mcNBColl = hDpmJet->NN()+hDpmJet->NNw()+hDpmJet->NwN()+hDpmJet->NwNw();
    }
    else {} // unknown generator
    //
    TArrayF vtxMC;
    mcGenH->PrimaryVertex(vtxMC);
    fGlobVars.mcZV = vtxMC[2];
  }
  //
  fOutTree->Fill();
  //
}      
//________________________________________________________________________
void AliTaskGlobVar::Terminate(Option_t *) 
{
  // print itself
  Printf("Terminating...");
  //  AliAnalysisTaskSE::Terminate();
}


//________________________________________________________________________
Float_t AliTaskGlobVar::GetCorrV0(const AliESDEvent* esd, float &v0CorrResc) const
{
  // correct V0 non-linearity, prepare a version rescaled to SPD2 corr
  const Double_t par0[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			      5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			      5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			      6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			      5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			      6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			      2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			      2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			      4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			      4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			      4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  const Double_t par1[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			      -2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			      -2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			      -2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			      -1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			      -2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			      1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			      3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			      -1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			      -3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			      -1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  const Double_t par2[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			      5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			      4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			      5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			      2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			      5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			      3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			      1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			      1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			      1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			      1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };
  //
  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  for(Int_t i = 0; i < 64; ++i) {
    Double_t b = (esdV0->GetMultiplicity(i)*par1[i]-par0[i]);
    Double_t s = (b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i));
    Double_t n;
    if (s<0) {
      printf("FPE %d %.2f %.2f %.2e\n",i,esdV0->GetMultiplicity(i),b,(b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i)));
      n = -b;
    }
    else {
      n = (-b + TMath::Sqrt(s));
    }
    multChCorr[i] = 2.*esdV0->GetMultiplicity(i)/n*par0[i];
    multCorr += multChCorr[i];
    multCorr2 += (multChCorr[i]/par0[i]/64.);
  }
  v0CorrResc =  multCorr2;
  return multCorr;
}

//_________________________________________________________________
Bool_t AliTaskGlobVar::ZDCTimeTrigger(const AliESDEvent *aEsd) const
{
  // This method implements a selection
  // based on the timing in both sides of zdcN
  // It can be used in order to eliminate
  // parasitic collisions
  Bool_t zdcAccept = kFALSE;

  AliESDZDC *esdZDC = aEsd->GetESDZDC();

  const Float_t refSum = -568.5;
  const Float_t refDelta = -2.1;
  const Float_t sigmaSum = 3.25;
  const Float_t sigmaDelta = 2.25;
  for(Int_t i = 0; i < 4; ++i) {
    if (esdZDC->GetZDCTDCData(10,i) != 0) {
      Float_t tdcC = 0.025*(esdZDC->GetZDCTDCData(10,i)-esdZDC->GetZDCTDCData(14,i)); 
      for(Int_t j = 0; j < 4; ++j) {
	if (esdZDC->GetZDCTDCData(12,j) != 0) {
	  Float_t tdcA = 0.025*(esdZDC->GetZDCTDCData(12,j)-esdZDC->GetZDCTDCData(14,j));
	  if (((tdcC-tdcA-refDelta)*(tdcC-tdcA-refDelta)/(sigmaDelta*sigmaDelta) +
	       (tdcC+tdcA-refSum)*(tdcC+tdcA-refSum)/(sigmaSum*sigmaSum))< 1.0)
	    zdcAccept = kTRUE;
	}
      }
    }
  }
  return zdcAccept;
}


AliESDtrackCuts* AliTaskGlobVar::CreatedNdPtTrackCuts(Int_t cutMode, Bool_t fieldOn)
{
  // copy of PWG0/dNdPt/macros/CreatedNdPtTrackCuts.C
  //
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");


  Double_t cov1, cov2, cov3, cov4, cov5;
  Double_t nSigma;
  Double_t maxDCAtoVertex, maxDCAtoVertexXY, maxDCAtoVertexZ;
  Int_t minNClustersTPC;
  Double_t maxChi2PerClusterTPC;
  Double_t minPt, maxPt;
  TString tag;

  // default cuts for ITS+TPC
  if (cutMode == 0) 
  {
    cov1 = 2;
    cov2 = 2;
    cov3 = 0.5;
    cov4 = 0.5;
    cov5 = 2;
    nSigma = 3;
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
    //    esdTrackCuts->SetMinNsigmaToVertex(nSigma);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "Global tracking";
  }

  // TPC-only cuts (vertex n sigma cut)
  if (cutMode == 1) 
  {
    // beta cuts (still under investigation)
    //cov1 = 4;
    //cov2 = 4;
    cov1 = 2;
    cov2 = 2;
    cov3 = 0.5;
    cov4 = 0.5;
    cov5 = 2;
    nSigma = 4;
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
    //    esdTrackCuts->SetMinNsigmaToVertex(nSigma);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only cuts (vertex maxDCAtoVertex cut)
  if (cutMode == 2) 
  {
    // beta cuts (still under investigation)
    maxDCAtoVertex = 3.0; // cm
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);    
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);    
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only no vertex cuts
  if (cutMode == 3) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only no cuts at all 
  if (cutMode == 4) 
  {

    // beta cuts (still under investigation)
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);

    tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no chi2 
  if (cutMode == 5) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only no kink removal 
  if (cutMode == 6) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no minNClustersTPC 
  if (cutMode == 7) 
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }
  // TPC-only no kink removal no minNClustersTPC 
  if (cutMode == 8) 
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertex = 3.0; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);    
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);    
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only no kink removal no minNClustersTPC no maxChi2PerClusterTPC
  if (cutMode == 9) 
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertex = 3.0; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);    
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);    
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    tag = "TPC-only tracking";
  }

  // TPC-only (loose cuts, absolute DCA cut) 
  if (cutMode == 10) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertex = 2.8; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertex);    
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertex);    
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }


  // TPC-only (loose cuts, no DCA cut) 
  if (cutMode == 11) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts, no DCA cut) 
  if (cutMode == 12) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.2;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (tight cuts, no DCA cut) 
  if (cutMode == 13) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.3;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (loose cuts, no pt cut) 
  if (cutMode == 14) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts, no pt cut) 
  if (cutMode == 15) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (tight cuts, no pt cuts) 
  if (cutMode == 16) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.5;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }
  // TPC-only (loose cuts)
  if (cutMode == 17) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.6; // cm
    maxDCAtoVertexZ  = 2.1; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (standard cuts)
  if (cutMode == 18) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 96;
    maxChi2PerClusterTPC = 3.5;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 1.8; // cm
    minPt=0.2;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

 // TPC-only (tight cuts)
  if (cutMode == 19) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 120;
    maxChi2PerClusterTPC = 3.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 1.8; // cm
    minPt=0.3;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts included)
  if (cutMode == 20) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts excluded)
  if (cutMode == 21) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only (arb. cuts, kink cuts excluded, no chi2, no DCA)
  if (cutMode == 22) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 1.e10;
    maxDCAtoVertexXY = 1.e10; // cm
    maxDCAtoVertexZ  = 1.e10; // cm
    minPt=0.15;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  // TPC-only 
  if (cutMode == 23) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    tag = "TPC-only tracking";
  }

  // TPC-only (no pt cut, no eta cut)
  if (cutMode == 24) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  //
  // systematic errors DCA cut studies
  //
  // TPC-only
  if (cutMode == 25) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 2.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 26) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.6; // cm
    maxDCAtoVertexZ  = 2.4; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  //
  // systematic errors cut studies
  //
  // TPC-only
  if (cutMode == 27) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.8; // cm
    maxDCAtoVertexZ  = 2.6; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 28) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.0; // cm
    maxDCAtoVertexZ  = 2.8; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 29) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.2; // cm
    maxDCAtoVertexZ  = 3.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 30) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 31) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.6; // cm
    maxDCAtoVertexZ  = 3.4; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }


  if (cutMode == 32) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.8; // cm
    maxDCAtoVertexZ  = 3.6; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 33) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.0; // cm
    maxDCAtoVertexZ  = 3.8; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 34) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.2; // cm
    maxDCAtoVertexZ  = 4.0; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 35) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.4; // cm
    maxDCAtoVertexZ  = 4.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

//
// cut stability systematics
//

  if (cutMode == 36) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

 if (cutMode == 37) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 90;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 38) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 39) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 5.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 40) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.4; // cm
    maxDCAtoVertexZ  = 2.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 41) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 3.4; // cm
    maxDCAtoVertexZ  = 4.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }

  if (cutMode == 42) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    minPt=0.0;
    maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    //esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    esdTrackCuts->SetPtRange(minPt,maxPt);

    tag = "TPC-only tracking";
  }
  // test
  if (cutMode == 43) 
  {
    // beta cuts (still under investigation)
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 4.0;
    //maxDCAtoVertexXY = 2.4; // cm
    //maxDCAtoVertexZ  = 3.2; // cm
    //minPt=0.15;
    //maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    //esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    //esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    //esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    tag = "TPC-only tracking";
  }

  // TPC-only + pt cut + eta cut 
  if (cutMode == 45) 
  {
    // beta cuts (still under investigation)
    //minNClustersTPC = 50;
    //maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 2.4; // cm
    maxDCAtoVertexZ  = 3.2; // cm
    //minPt=0.15;
    //maxPt=1.e10;

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    //esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    //esdTrackCuts->SetPtRange(minPt,maxPt);
    //esdTrackCuts->SetEtaRange(minEta,maxEta);

    tag = "TPC-only tracking";
  }

  // TPC-tracks + SPD point + ITS refit
  if (cutMode == 50) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    //Double_t maxEtaInAcc=0.8;
    Double_t maxdcaxyITSTPC=0.2;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //esdTrackCuts->SetEtaRange(-maxEtaInAcc,maxEtaInAcc);

    tag = "TPC-tracks + ITS refit + >1 SPD cluster";
  }

  // TPC-tracks + SPD point + ITS refit
  if (cutMode == 60) 
  {
    Int_t    minclsITS=4;
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcaxyITSTPC=0.2;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMinNClustersITS(minclsITS);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);

    tag = "Global tracking: TPC refit + ITS refit + >3 ITS clusters + >=1 SPD cluster";
  }

  /*
  // TPC-tracks + SPD point + ITS refit + DCAr(pt)
  if (cutMode == 70) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcaxyITSTPC=1.e9;
    Double_t maxdcazITSTPC=1.e9;

    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt)";
  }
  */

  // TPC-tracks + SPD point + ITS refit + DCAr(pt)
  if (cutMode == 70) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt)";
  }

  // TPC+ITS combine tracking + DCAr(pt) + DCAz(pt)
  if (cutMode == 71) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.01+0.011/pt^0.72)
    esdTrackCuts->SetMaxDCAToVertexZPtDep("0.07+0.077/pt^0.72");

    tag = "TPC+ITS combine tracking + DCAr(pt) + DCAz(pt)";
  }

  // TPC+ITS combine tracking + DCAr(pt) (2010)
  if (cutMode == 72) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=2.0;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // DCArphi parametrization (LHC10c pass2)
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

    tag = "TPC+ITS combine tracking + DCAr(pt) (2010)";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) 4-sigma
  if (cutMode == 75) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 4*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.02+0.024/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) 4-sigma";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) 10-sigma
  if (cutMode == 80) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 10*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.05+0.06/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) 10 sigma";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 60 TPCclust
  if (cutMode == 85) 
  {
    Int_t    minclsTPC=60;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 80 clusters
  if (cutMode == 90) 
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 80 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + TPCchi2=3.5
  if (cutMode == 95) 
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=3.5;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + TPCchi2 3.5";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + TPCchi2=4.5
  if (cutMode == 100) 
  {
    Int_t    minclsTPC=80;
    Double_t maxchi2perTPCcl=4.5;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + TPCchi2 4.5";
  }

  // TPC-tracks
  if (cutMode == 110) 
  {

    minNClustersTPC = 70;
    maxChi2PerClusterTPC = 4.0;
    maxDCAtoVertexXY = 1.e9; // cm
    maxDCAtoVertexZ  = 1.e9; // cm

    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxDCAtoVertexXY);
    esdTrackCuts->SetMaxDCAToVertexZ(maxDCAtoVertexZ);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);

    tag = "TPC-tracks loose criteria";
  }


  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 50 TPCclust
  if (cutMode == 120) 
  {
    Int_t    minclsTPC=50;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 70 TPCclust + accept kink daughters
  if (cutMode == 130) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // TPC-tracks + SPD point + ITS refit + DCAr(pt) + 30 TPCclust + accept kink daughters
  if (cutMode == 140) 
  {
    Int_t    minclsTPC=30;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=1.e9;

    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "TPC-tracks + ITS refit + >1 SPD cluster + DCAr(Pt) + 60 TPCclust";
  }

  // Adam Kisiel track selectiion
  if (cutMode == 150) 
  {
    Int_t    minclsTPC=70;
    Double_t maxchi2perTPCcl=4.;
    Double_t maxdcazITSTPC=0.25;
    Double_t maxdcaxyITSTPC=0.2;

    //
    // TPC
    //
    //esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minclsTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    // primary selection
    //
    //esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
    esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyITSTPC);

    // 7*(0.0050+0.0060/pt^0.9)
    //esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");

    tag = "Adam Kisiel track selection";
  }

  // TPC+ITS refit
  // for cut studies
  if (cutMode == 151) 
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //

    tag = "TPC+ITS refit required - for cut studies";
  }

  // TPC+ITS
  // for cut studies
  if (cutMode == 152) 
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    //esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    
    tag = "TPC refit required - for cut studies";
  }

  // TPC
  // for cut studies
  if (cutMode == 153) 
  {
    //
    // TPC
    //
    esdTrackCuts->SetRequireTPCRefit(kFALSE);
    esdTrackCuts->SetRequireITSRefit(kFALSE);
    esdTrackCuts->SetRequireTPCStandAlone(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //
    // ITS
    //
    //esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //
    
    tag = "TPC stand alone - for cut studies";
  }






  // cuts for data without field
  if (!fieldOn)
  {
    cov5 = 1e10;
    tag += " without field";
  }

  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}
