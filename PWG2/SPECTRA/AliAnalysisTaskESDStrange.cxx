
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

// macro to study V0s
// loops over ESD files, and creates AliAODv0
// Author: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"

#include "AliAnalysisTaskESDStrange.h"


ClassImp(AliAnalysisTaskESDStrange)

//________________________________________________________________________
AliAnalysisTaskESDStrange::AliAnalysisTaskESDStrange(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fListHist(), 
    fHistTrackPerEvent(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),
    fHistDcaPosToPrimVertex(0),
    fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0),
    fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0),
    fHistDecayLengthV0(0),
    fHistDcaV0Daughters(0),
    fHistChi2(0),
    fHistCosPointAngle(0),
    fHistCosPointAngleZoom(0),
    fHistPtVsYK0s(0),
    fHistPtVsYK0sMI(0),
    fHistPtVsYLambda(0),
    fHistPtVsYLambdaMI(0),
    fHistPtVsYAntiLambda(0),
    fHistPtVsYAntiLambdaMI(0),
    fHistMassK0(0),
    fHistMassK0MI(0),
    fHistMassLambda(0),
    fHistMassLambdaMI(0),
    fHistMassAntiLambda(0),
    fHistMassAntiLambdaMI(0),
    fHistArmenterosPodolanski(0),
    fHistArmenterosPodolanskiMI(0)
    
    
    
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskESDStrange::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    tree->SetBranchStatus("fV0s.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskESDStrange::CreateOutputObjects() 
{
  // Create histograms
  // Called once

  fListHist = new TList();

  // multiplicity
  fHistTrackPerEvent        = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",50,0,50);
  fListHist->Add(fHistTrackPerEvent);

  // Primary Vertex:
  fHistPrimaryVertexX       = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",40,-1,1);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY       = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",40,-1,1);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ       = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",60,-5,5);
  fListHist->Add(fHistPrimaryVertexZ);

  // Cut checks:
  fHistDcaPosToPrimVertex   = new TH2F("h2DcaPosToPrimVertex", "Positive V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertex);

  fHistDcaNegToPrimVertex   = new TH2F("h2DcaNegToPrimVertex", "Negative V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertex);

  fHistDcaPosToPrimVertexZoom  = new TH2F("h2DcaPosToPrimVertexZoom", "Positive V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexZoom);

  fHistDcaNegToPrimVertexZoom  = new TH2F("h2DcaNegToPrimVertexZoom", "Negative V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexZoom);

  fHistRadiusV0             = new TH2F("h2RadiusV0", "Radius;Radius(cm);Status",1000,0,100,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0);

  fHistDecayLengthV0        = new TH2F("h2DecayLengthV0", "V0s decay Length;decay length(cm);Status", 200, 0, 100,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0);

  fHistDcaV0Daughters       = new TH2F("h2DcaV0Daughters", "DCA between daughters;dca(cm);Status", 160, 0, 4,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0Daughters);

  fHistChi2                 = new TH2F("h2Chi2", "V0s chi2;chi2;Status", 33, 0, 33,2,-0.5,1.5);
  fListHist->Add(fHistChi2);

  fHistCosPointAngle           = new TH2F("h2CosPointAngle", "Cosine of V0's pointing angle", 100,0,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngle);

  fHistCosPointAngleZoom       = new TH2F("h2CosPointAngleZoom", "Cosine of V0's pointing angle", 100,0.9,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleZoom);


  // Pt and rapidity distribution:
  fHistPtVsYK0s             = new TH2F("h2PtVsYK0s", "K^{0} candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0s);
  fHistPtVsYK0sMI           = new TH2F("h2PtVsYK0sMI", "K^{0} MI candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0sMI);

  fHistPtVsYLambda          = new TH2F("h2PtVsYLambda", "#Lambda^{0} candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambda);
  fHistPtVsYLambdaMI        = new TH2F("h2PtVsYLambdaMI", "#Lambda^{0} MI candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambdaMI);

  fHistPtVsYAntiLambda      = new TH2F("h2PtVsYAntiLambda", "#bar{#Lambda}^{0} candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambda);
  fHistPtVsYAntiLambdaMI    = new TH2F("h2PtVsYAntiLambdaMI", "#bar{#Lambda}^{0} MI candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambdaMI);


  // Mass:
  fHistMassK0               = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0);
  fHistMassK0MI             = new TH1F("h1MassK0MI", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0MI);

  fHistMassLambda           = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);
  fHistMassLambdaMI         = new TH1F("h1MassLambdaMI", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambdaMI);

  fHistMassAntiLambda       = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);
  fHistMassAntiLambdaMI     = new TH1F("h1MassAntiLambdaMI", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambdaMI);


  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fHistArmenterosPodolanskiMI   = new TH2F("h2ArmenterosPodolanskiMI","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);

  
}

//________________________________________________________________________
void AliAnalysisTaskESDStrange::Exec(Option_t *) 
{
  //*********************************************
  // Called for each event
  //*********************************************

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  //Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  if (!(fESD->GetNumberOfTracks())) {
    Printf("No ESD track in the event");
    return;
  }
  fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());

  // Primary Vertex
  Double_t  PrimaryVtxPosition[3];
  Double_t  PrimaryVtxCov[6];

  const AliESDVertex *primaryVtx = fESD->GetPrimaryVertex();

  primaryVtx->GetXYZ(PrimaryVtxPosition);
  primaryVtx->GetCovMatrix(PrimaryVtxCov); 

  AliAODVertex *primary = new AliAODVertex(PrimaryVtxPosition, PrimaryVtxCov, primaryVtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);

  fHistPrimaryVertexX->Fill(primary->GetX());
  fHistPrimaryVertexY->Fill(primary->GetY());
  fHistPrimaryVertexZ->Fill(primary->GetZ());

  // V0 variables:
  // to get info from ESD files and fill AliAODVertex:
  Float_t   tdcaPosToPrimVertexXYZ[2], tdcaNegToPrimVertexXYZ[2]; // ..[0] = Impact parameter in XY plane and ..[1] = Impact parameter in Z            
  Double_t  tdcaDaughterToPrimVertex[2];                          // ..[0] = Pos and ..[1] = Neg
  Double_t  tdcaV0Daughters     = 0, tdcaV0ToPrimVertex   = 0;
  Double_t  tMomPos[3];
  Double_t  tMomNeg[3];
  Double_t  V0Position[3];
  Double_t  V0Cov[6];

  // to fill AliAODtrack:
  Double_t  TrackP[3];
  Double_t  TrackPosition[3];
  Double_t  TrackcovTr[21];
  Double_t  Trackpid[10];


  Double_t rcPosX        = 0,  rcPosY  = 0, rcPosZ  = 0;
  Double_t rcPosR        = 0;
  Double_t cosPointAngle = 0;
  Double_t deltaPos2     = 0;
  Double_t deltaPos[3];

  Int_t    myStatus = 0;
  Double_t pt       = 0;

  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;
  UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg       = 0;

  AliAODTrack  *myPosAodTrack  = new AliAODTrack();
  AliAODTrack  *myNegAodTrack  = new AliAODTrack();
  AliAODVertex *myAODVertex    = new AliAODVertex();
  AliAODv0     *myAODv0        = new AliAODv0();

  //cout<<"number V0s:"<<fESD->GetNumberOfV0s()<<endl;


  //*************************
  // V0 loop
  //*************************
  for (Int_t iV0 = 0; iV0 < fESD->GetNumberOfV0s(); iV0++) {
  
    AliESDv0 *v0 = fESD->GetV0(iV0);
    if (!v0) continue;
    myAODv0->ResetV0();

    // AliAODVertex
    v0->GetXYZ(V0Position[0], V0Position[1], V0Position[2]);
    v0->GetPosCov(V0Cov);
    myAODVertex->SetPosition(V0Position[0],V0Position[1],V0Position[2]);
    myAODVertex->SetCovMatrix(V0Cov);
    myAODVertex->SetChi2perNDF(v0->GetChi2V0());
    myAODVertex->SetID((Short_t)iV0);
    myAODVertex->SetParent(primary);
    myAODVertex->SetType(AliAODVertex::kV0);

    // AliAODtrack (V0 Daughters)
    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *TrackPos = fESD->GetTrack(lIndexTrackPos);
    AliESDtrack *TrackNeg = fESD->GetTrack(lIndexTrackNeg);
    if (!TrackPos || !TrackNeg) {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }
    lLabelTrackPos = (UInt_t)TMath::Abs(TrackPos->GetLabel());
    lLabelTrackNeg = (UInt_t)TMath::Abs(TrackNeg->GetLabel());

    myPosAodTrack->SetID((Short_t)(TrackPos->GetID()));  
    myPosAodTrack->SetLabel(lLabelTrackPos);
    TrackPos->GetPxPyPz(TrackP);
    myPosAodTrack->SetP(TrackP);
    TrackPos->GetXYZ(TrackPosition);
    myPosAodTrack->SetPosition(TrackPosition,kFALSE);
    TrackPos->GetCovarianceXYZPxPyPz(TrackcovTr);
    myPosAodTrack->SetCovMatrix(TrackcovTr);
    TrackPos->GetESDpid(Trackpid);
    myPosAodTrack->SetPID(Trackpid);
    myPosAodTrack->SetCharge((Short_t)(TrackPos->Charge()));
    myPosAodTrack->SetITSClusterMap(TrackPos->GetITSClusterMap());
    myPosAodTrack->SetProdVertex(myAODVertex);
    myPosAodTrack->SetUsedForVtxFit(kTRUE);
    myPosAodTrack->SetUsedForPrimVtxFit(kFALSE);
    myPosAodTrack->SetType(AliAODTrack::kSecondary);
    myPosAodTrack->ConvertAliPIDtoAODPID();

    myNegAodTrack->SetID((Short_t)(TrackNeg->GetID()));
    myNegAodTrack->SetLabel(lLabelTrackNeg);
    TrackNeg->GetPxPyPz(TrackP);
    myNegAodTrack->SetP(TrackP);
    TrackNeg->GetXYZ(TrackPosition);
    myNegAodTrack->SetPosition(TrackPosition,kFALSE);
    TrackNeg->GetCovarianceXYZPxPyPz(TrackcovTr);
    myNegAodTrack->SetCovMatrix(TrackcovTr);
    TrackNeg->GetESDpid(Trackpid);
    myNegAodTrack->SetPID(Trackpid);
    myNegAodTrack->SetCharge((Short_t)(TrackNeg->Charge()));
    myNegAodTrack->SetITSClusterMap(TrackPos->GetITSClusterMap());
    myNegAodTrack->SetProdVertex(myAODVertex);
    myNegAodTrack->SetUsedForVtxFit(kTRUE);
    myNegAodTrack->SetUsedForPrimVtxFit(kFALSE);
    myNegAodTrack->SetType(AliAODTrack::kSecondary);
    myNegAodTrack->ConvertAliPIDtoAODPID();
   
    myAODVertex->AddDaughter(myPosAodTrack);
    myAODVertex->AddDaughter(myNegAodTrack);

    // filling myAODv0
    tdcaV0Daughters    = v0->GetDcaV0Daughters();
    tdcaV0ToPrimVertex = v0->GetD();

    if (TrackPos) TrackPos->GetImpactParameters(tdcaPosToPrimVertexXYZ[0],tdcaPosToPrimVertexXYZ[1]);
    if (TrackNeg) TrackNeg->GetImpactParameters(tdcaNegToPrimVertexXYZ[0],tdcaNegToPrimVertexXYZ[1]);
    tdcaDaughterToPrimVertex[0] = TMath::Sqrt(tdcaPosToPrimVertexXYZ[0]*tdcaPosToPrimVertexXYZ[0]+tdcaPosToPrimVertexXYZ[1]*tdcaPosToPrimVertexXYZ[1]);
    tdcaDaughterToPrimVertex[1] = TMath::Sqrt(tdcaNegToPrimVertexXYZ[0]*tdcaNegToPrimVertexXYZ[0]+tdcaNegToPrimVertexXYZ[1]*tdcaNegToPrimVertexXYZ[1]);

    v0->GetPPxPyPz(tMomPos[0],tMomPos[1],tMomPos[2]); 
    v0->GetNPxPyPz(tMomNeg[0],tMomNeg[1],tMomNeg[2]); 

    myAODv0->Fill(myAODVertex, tdcaV0Daughters, tdcaV0ToPrimVertex, tMomPos, tMomNeg, tdcaDaughterToPrimVertex);
    

    // Reconstructed V0 position and Cos pointing angle:
    rcPosX   = myAODv0->DecayVertexV0X();
    rcPosY   = myAODv0->DecayVertexV0Y();
    rcPosZ   = myAODv0->DecayVertexV0Z();
    rcPosR   = TMath::Sqrt(rcPosX*rcPosX+rcPosY*rcPosY);

    deltaPos[0]   = rcPosX - PrimaryVtxPosition[0];
    deltaPos[1]   = rcPosY - PrimaryVtxPosition[1];
    deltaPos[2]   = rcPosZ - PrimaryVtxPosition[2];
    deltaPos2     = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
    cosPointAngle = (deltaPos[0]*(myAODv0->MomV0X()) + deltaPos[1]*(myAODv0->MomV0Y()) + deltaPos[2]*(myAODv0->MomV0Z()))/TMath::Sqrt((myAODv0->Ptot2V0())*deltaPos2);


    myStatus = v0->GetOnFlyStatus();
    pt       = TMath::Sqrt(myAODv0->Ptot2V0());

    // filling histograms
    fHistDcaPosToPrimVertex->Fill(myAODv0->DcaPosToPrimVertex(),myStatus);
    fHistDcaNegToPrimVertex->Fill(myAODv0->DcaNegToPrimVertex(),myStatus);
    fHistDcaPosToPrimVertexZoom->Fill(myAODv0->DcaPosToPrimVertex(),myStatus);
    fHistDcaNegToPrimVertexZoom->Fill(myAODv0->DcaNegToPrimVertex(),myStatus);
    fHistRadiusV0->Fill(myAODv0->RadiusV0(),myStatus);
    fHistDecayLengthV0->Fill(myAODv0->DecayLengthV0(PrimaryVtxPosition),myStatus);
    fHistDcaV0Daughters->Fill(myAODv0->DcaV0Daughters(),myStatus);
    fHistChi2->Fill(myAODv0->Chi2V0(),myStatus);
    fHistCosPointAngle->Fill(cosPointAngle,myStatus);
    if (cosPointAngle >= 0.9) fHistCosPointAngleZoom->Fill(cosPointAngle,myStatus);
    if (!myStatus) {
      fHistPtVsYK0s->Fill(pt,myAODv0->RapK0Short());
      fHistPtVsYLambda->Fill(pt,myAODv0->RapLambda());
      fHistPtVsYAntiLambda->Fill(pt,myAODv0->RapLambda());
      fHistArmenterosPodolanski->Fill(myAODv0->AlphaV0(),myAODv0->PtArmV0());
    }
    else {
      fHistPtVsYK0sMI->Fill(pt,myAODv0->RapK0Short());
      fHistPtVsYLambdaMI->Fill(pt,myAODv0->RapLambda());
      fHistPtVsYAntiLambdaMI->Fill(pt,myAODv0->RapLambda());
      fHistArmenterosPodolanskiMI->Fill(myAODv0->AlphaV0(),myAODv0->PtArmV0());
    }
    // K0s histograms:
    if (TMath::Abs(myAODv0->RapK0Short()) < 1) {  
	if (!myStatus) fHistMassK0->Fill(myAODv0->MassK0Short());	
	else fHistMassK0MI->Fill(myAODv0->MassK0Short());	
      }

    // Lambda and AntiLambda histograms
    if (TMath::Abs(myAODv0->RapLambda()) < 1) {
      if (!myStatus) {
	fHistMassLambda->Fill(myAODv0->MassLambda());
	fHistMassAntiLambda->Fill(myAODv0->MassAntiLambda());
      }
      else {
	fHistMassLambdaMI->Fill(myAODv0->MassLambda());
	fHistMassAntiLambdaMI->Fill(myAODv0->MassAntiLambda());
      }
    }
 
  } // end V0 loop

  if (primary) delete primary;
  if (myPosAodTrack) delete myPosAodTrack;
  if (myNegAodTrack) delete myNegAodTrack;
  if (myAODv0) delete myAODv0;
  
  
  // Post output data.
  PostData(0, fListHist);
}      

//________________________________________________________________________
void AliAnalysisTaskESDStrange::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  /*
  fHistTrackPerEvent = dynamic_cast<TH1F*> (((TList*)GetOutputData(0))->FindObject("fHistTrackPerEvent"));
  if (!fHistTrackPerEvent) {
    Printf("ERROR: fHistTrackPerEvent not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskESDStrange","TrackPerEvent",10,10,510,510);
  c1->cd(1);
  fHistTrackPerEvent->DrawCopy("E");
*/
}

//----------------------------------------------------------------------------

Double_t myRap(Double_t rE, Double_t rPz)
{
  Double_t lRapidity = 1.e30;
  if(rPz && rE && (rPz != rE) && (1.+(2./((rE/rPz)-1.))>0))
    lRapidity = 0.5*log(1.+(2./((rE/rPz)-1.)));
  return lRapidity;
  } 

//----------------------------------------------------------------------------

