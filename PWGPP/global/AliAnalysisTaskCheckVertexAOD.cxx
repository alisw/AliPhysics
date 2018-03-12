#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODTracklets.h"
#include <AliAODMCParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliAnalysisTaskCheckVertexAOD.h"


/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskCheckVertexAOD
// AliAnalysisTaskSE to extract QA plots for vertices in AODs
// 
//
// Authors: 
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

ClassImp(AliAnalysisTaskCheckVertexAOD)
//______________________________________________________________________________
AliAnalysisTaskCheckVertexAOD::AliAnalysisTaskCheckVertexAOD() : 
  AliAnalysisTaskSE("ITSsa resolution"), 
  fOutput{nullptr},
  fHistNEvents{nullptr},
  fHistAllVtxType{nullptr},
  fHistPrimVtxType{nullptr},
  fHistXspdVsContrib{nullptr},
  fHistYspdVsContrib{nullptr},
  fHistZspdVsContrib{nullptr},
  fHistXtrkVsContrib{nullptr},
  fHistYtrkVsContrib{nullptr},
  fHistZtrkVsContrib{nullptr},
  fHistXtpcVsContrib{nullptr},
  fHistYtpcVsContrib{nullptr},
  fHistZtpcVsContrib{nullptr},  
  fHistoNOfPileupVertSPD{nullptr},
  fHistoNOfSelPileupVertSPD{nullptr},
  fHistoNOfPileupVertMV{nullptr},
  fHistoNOfSelPileupVertMV{nullptr},
  fUsePhysSel(kTRUE),
  fTriggerMask(AliVEvent::kAnyINT),
  fMaxMult(500.),
  fSPDContributorsCut(5),
  fSPDZDiffCut(0.8),
  fMVContributorsCut(5),
  fMVCChi2Cut(5.),
  fMVWeiZDiffCut(15.),
  fMVCheckPlpFromDifferentBC(kFALSE),
  fReadMC{kFALSE}
{
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckVertexAOD::~AliAnalysisTaskCheckVertexAOD(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistAllVtxType;
    delete fHistPrimVtxType;
    delete fHistXspdVsContrib;
    delete fHistYspdVsContrib;
    delete fHistZspdVsContrib;
    delete fHistXtrkVsContrib;
    delete fHistYtrkVsContrib;
    delete fHistZtrkVsContrib;
    delete fHistXtpcVsContrib;
    delete fHistYtpcVsContrib;
    delete fHistZtpcVsContrib;  
    delete fHistoNOfPileupVertSPD;
    delete fHistoNOfSelPileupVertSPD;
    delete fHistoNOfPileupVertMV;
    delete fHistoNOfSelPileupVertMV;
  }
  delete fOutput;
}
 
//___________________________________________________________________________
void AliAnalysisTaskCheckVertexAOD::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel"); 
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Good vertex"); 
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Pass zSPD-zTrk vert sel"); 
  fHistNEvents->GetXaxis()->SetBinLabel(5,"|zvert|<10"); 
  fOutput->Add(fHistNEvents);

  fHistAllVtxType = new TH1F("hAllVtxType"," ; Vertex Type ; Entries",12,0.5,12.5);
  fHistAllVtxType->GetXaxis()->SetBinLabel(1,"kPrimaryInvalid");
  fHistAllVtxType->GetXaxis()->SetBinLabel(2,"kUndef");
  fHistAllVtxType->GetXaxis()->SetBinLabel(3,"kPrimary");
  fHistAllVtxType->GetXaxis()->SetBinLabel(4,"kKink");
  fHistAllVtxType->GetXaxis()->SetBinLabel(5,"kV0");
  fHistAllVtxType->GetXaxis()->SetBinLabel(6,"kCascade");
  fHistAllVtxType->GetXaxis()->SetBinLabel(7,"kMulti");
  fHistAllVtxType->GetXaxis()->SetBinLabel(8,"kMainSPD");
  fHistAllVtxType->GetXaxis()->SetBinLabel(9,"kPileupSPD");
  fHistAllVtxType->GetXaxis()->SetBinLabel(10,"kPileupTracks");
  fHistAllVtxType->GetXaxis()->SetBinLabel(11,"kMainTPC");
  fHistAllVtxType->GetXaxis()->SetBinLabel(12,"kPrimaryTPC");
  fOutput->Add(fHistAllVtxType);

  fHistPrimVtxType = new TH1F("hPrimVtxType"," ; Vertex Type ; Entries",8,0.5,8.5);
  fHistPrimVtxType->GetXaxis()->SetBinLabel(1,"kPrimaryInvalid");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(2,"kUndef");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(3,"TrackVertex");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(4,"SPD3DVertex");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(5,"SPDZVertex");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(6,"kPrimaryTPC");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(7,"TPCVertex (old AOD)");
  fHistPrimVtxType->GetXaxis()->SetBinLabel(8,"Other");
  fOutput->Add(fHistPrimVtxType);

  fHistXspdVsContrib=new TH2F("hXspdVsContrib"," ; n_{Contributors} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYspdVsContrib=new TH2F("hYspdVsContrib"," ; n_{Contributors} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZspdVsContrib=new TH2F("hZspdVsContrib"," ; n_{Contributors} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fHistXtrkVsContrib=new TH2F("hXtrkVsContrib"," ; n_{Contributors} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYtrkVsContrib=new TH2F("hYtrkVsContrib"," ; n_{Contributors} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZtrkVsContrib=new TH2F("hZtrkVsContrib"," ; n_{Contributors} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fHistXtpcVsContrib=new TH2F("hXtpcVsContrib"," ; n_{Contributors} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYtpcVsContrib=new TH2F("hYtpcVsContrib"," ; n_{Contributors} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZtpcVsContrib=new TH2F("hZtpcVsContrib"," ; n_{Contributors} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fOutput->Add(fHistXspdVsContrib);
  fOutput->Add(fHistYspdVsContrib);
  fOutput->Add(fHistZspdVsContrib);
  fOutput->Add(fHistXtrkVsContrib);
  fOutput->Add(fHistYtrkVsContrib);
  fOutput->Add(fHistZtrkVsContrib);
  fOutput->Add(fHistXtpcVsContrib);
  fOutput->Add(fHistYtpcVsContrib);
  fOutput->Add(fHistZtpcVsContrib);

  fHistXspdVsMult=new TH2F("hXspdVsMult"," ; n_{Tracklets} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYspdVsMult=new TH2F("hYspdVsMult"," ; n_{Tracklets} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZspdVsMult=new TH2F("hZspdVsMult"," ; n_{Tracklets} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fHistXtrkVsMult=new TH2F("hXtrkVsMult"," ; n_{Tracklets} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYtrkVsMult=new TH2F("hYtrkVsMult"," ; n_{Tracklets} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZtrkVsMult=new TH2F("hZtrkVsMult"," ; n_{Tracklets} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fHistXtpcVsMult=new TH2F("hXtpcVsMult"," ; n_{Tracklets} ; x_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistYtpcVsMult=new TH2F("hYtpcVsMult"," ; n_{Tracklets} ; y_{Vertex} (cm)",100,0.,fMaxMult,300,-1.,1.);
  fHistZtpcVsMult=new TH2F("hZtpcVsMult"," ; n_{Tracklets} ; z_{Vertex} (cm)",100,0.,fMaxMult,300,-20.,20.);
  fOutput->Add(fHistXspdVsMult);
  fOutput->Add(fHistYspdVsMult);
  fOutput->Add(fHistZspdVsMult);
  fOutput->Add(fHistXtrkVsMult);
  fOutput->Add(fHistYtrkVsMult);
  fOutput->Add(fHistZtrkVsMult);
  fOutput->Add(fHistXtpcVsMult);
  fOutput->Add(fHistYtpcVsMult);
  fOutput->Add(fHistZtpcVsMult);

  fHistoNOfPileupVertSPD = new TH1F("hNOfPileupVertSPD","",11,-0.5,10.5);
  fHistoNOfSelPileupVertSPD = new TH1F("hNOfSelPileupVertSPD","",11,-0.5,10.5);
  fHistoNOfPileupVertMV = new TH1F("hNOfPileupVertMV","",11,-0.5,10.5);
  fHistoNOfSelPileupVertMV = new TH1F("hNOfSelPileupVertMV","",11,-0.5,10.5);
  fOutput->Add(fHistoNOfPileupVertSPD);
  fOutput->Add(fHistoNOfSelPileupVertSPD);
  fOutput->Add(fHistoNOfPileupVertMV);
  fOutput->Add(fHistoNOfSelPileupVertMV);

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckVertexAOD::UserExec(Option_t *)
{
  //

  AliAODEvent *aod = (AliAODEvent*) (InputEvent());
  if(!aod) {
    printf("AliAnalysisTaskCheckVertexAOD::UserExec(): bad AOD\n");
    return;
  } 

  TClonesArray *arrayMC=0;

  if(fReadMC){
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      Printf("ERROR: MC particles branch not found!\n");
      return;
    }
  }


  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  for(Int_t jv=0; jv<aod->GetNumberOfVertices(); jv++){
    AliAODVertex *v=(AliAODVertex*)aod->GetVertex(jv);
    Int_t typ=v->GetType();
    Int_t val=typ+3;
    if(typ==AliAODVertex::kPrimaryInvalid) val=1;
    else if(typ==AliAODVertex::kUndef) val=2;
    fHistAllVtxType->Fill(val);
  }

  const AliAODVertex* vtPrim = (AliAODVertex*)aod->GetPrimaryVertex();
  const AliAODVertex *vtTPC  = aod->GetPrimaryVertexTPC();

  Int_t typp=vtPrim->GetType();
  Int_t val=8;
  if(typp==AliAODVertex::kPrimaryInvalid) val=1;
  else if(typp==AliAODVertex::kUndef) val=2;
  else if(typp==AliAODVertex::kPrimaryTPC) val=6;
  else{
    TString title=vtPrim->GetTitle();
    if(title.Contains("VertexerTracks")){
      // track vertex (TPC or global)
      if (TMath::Abs(vtPrim->GetZ()-vtTPC->GetZ())<1e-6 &&
	  TMath::Abs(vtPrim->GetChi2perNDF()-vtTPC->GetChi2perNDF())<1e-6) {
	// TPC vertex
	val=7;
      }else{
	// track vertex
	val=3;
      }
    }else{
      // SPD vertex
      if(title.Contains("ertexer: Z")) val=5;
      if(title.Contains("ertexer: 3D")) val=4;
    }
  }
  fHistPrimVtxType->Fill(val);

  Int_t ntracklets = 0;
  AliAODTracklets *mult=aod->GetTracklets();
  if(mult) ntracklets=mult->GetNumberOfTracklets();

  const AliVVertex* vtSPD = aod->GetPrimaryVertexSPD();
  Int_t ct=0;
  Float_t zt=-999.;
  if(val==3){
    Float_t xt=vtPrim->GetX();
    Float_t yt=vtPrim->GetY();
    zt=vtPrim->GetZ();
    ct=vtPrim->GetNContributors();
    fHistXtrkVsContrib->Fill(ct,xt);
    fHistYtrkVsContrib->Fill(ct,yt);
    fHistZtrkVsContrib->Fill(ct,zt);
    fHistXtrkVsMult->Fill(ntracklets,xt);
    fHistYtrkVsMult->Fill(ntracklets,yt);
    fHistZtrkVsMult->Fill(ntracklets,zt);
  }
  Int_t cs=0;
  Float_t zs=-999.;
  if(vtSPD){
    Float_t xs=vtSPD->GetX();
    Float_t ys=vtSPD->GetY();
    zs=vtSPD->GetZ();
    cs=vtSPD->GetNContributors();
    fHistXspdVsContrib->Fill(cs,xs);
    fHistYspdVsContrib->Fill(cs,ys);
    fHistZspdVsContrib->Fill(cs,zs);
    fHistXspdVsMult->Fill(ntracklets,xs);
    fHistYspdVsMult->Fill(ntracklets,ys);
    fHistZspdVsMult->Fill(ntracklets,zs);
  }
  if(vtTPC){
    Float_t xtpc=vtTPC->GetX();
    Float_t ytpc=vtTPC->GetY();
    Float_t ztpc=vtTPC->GetZ();
    Int_t ctpc=vtTPC->GetNContributors();
    fHistXtpcVsContrib->Fill(ctpc,xtpc);
    fHistYtpcVsContrib->Fill(ctpc,ytpc);
    fHistZtpcVsContrib->Fill(ctpc,ztpc);
    fHistXtpcVsMult->Fill(ntracklets,xtpc);
    fHistYtpcVsMult->Fill(ntracklets,ytpc);
    fHistZtpcVsMult->Fill(ntracklets,ztpc);
  }

  Int_t nPileupVertSPD=aod->GetNumberOfPileupVerticesSPD();
  fHistoNOfPileupVertSPD->Fill(nPileupVertSPD);
  Int_t nSelPileupVertSPD=0;
  for(Int_t iv=0; iv<nPileupVertSPD; iv++){
    const AliAODVertex *spdvp=aod->GetPileupVertexSPD(iv);
    if(spdvp->GetNContributors()>=fSPDContributorsCut){
      Double_t zpile=spdvp->GetZ();
      Double_t zdiff=zpile-zs;
      if(TMath::Abs(zdiff)>fSPDZDiffCut) ++nSelPileupVertSPD;
    }
  }
  fHistoNOfSelPileupVertSPD->Fill(nSelPileupVertSPD);

  if(val==3){
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(fMVContributorsCut);
    utils.SetMaxPlpChi2MV(fMVCChi2Cut);
    utils.SetMinWDistMV(fMVWeiZDiffCut);
    utils.SetCheckPlpFromDifferentBCMV(fMVCheckPlpFromDifferentBC);
    Int_t nPileupVertMV=aod->GetNumberOfPileupVerticesTracks();
    fHistoNOfPileupVertMV->Fill(nPileupVertMV);
    Int_t nSelPileupVertMV=0;
    Int_t bcPrim=vtPrim->GetBC();
    for(Int_t iv=0; iv<nPileupVertMV; iv++){
      const AliAODVertex *trvp=aod->GetPileupVertexTracks(iv);
      Bool_t accept=kFALSE;
      if(trvp->GetNContributors()>=fMVContributorsCut && trvp->GetChi2perNDF() < fMVCChi2Cut){
	if(fMVCheckPlpFromDifferentBC){
	  Int_t bcPile = trvp->GetBC();
	  if (bcPile!=AliVTrack::kTOFBCNA && TMath::Abs(bcPile-bcPrim)>2){
	    accept=kTRUE;
	  }
	}
	Double_t wDst = utils.GetWDist(vtPrim,trvp);
	if (wDst>fMVWeiZDiffCut) accept=kTRUE;
	if(accept) ++nSelPileupVertMV;
      }
    }
    fHistoNOfSelPileupVertMV->Fill(nSelPileupVertMV);
  }

  if (ct<2 || cs<1) return; // one of vertices is missing
  fHistNEvents->Fill(2);
  
  double covPrim[6],covSPD[6];
  vtPrim->GetCovarianceMatrix(covPrim);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = vtPrim->GetZ()-vtSPD->GetZ();
  double errTot = TMath::Sqrt(covPrim[5]+covSPD[5]);
  double errPrim = TMath::Sqrt(covPrim[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigPrim = TMath::Abs(dz)/errPrim;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigPrim>20) return; // bad vertexing
  fHistNEvents->Fill(3);

  if(TMath::Abs(vtPrim->GetZ())>10) return;
  fHistNEvents->Fill(4);

  PostData(1,fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskCheckVertexAOD::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskCheckVertexAOD::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(5));
  return;
}





