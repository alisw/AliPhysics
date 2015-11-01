// ******************************************
// This task computes several jet observables like
// the fraction of energy in inner and outer coronnas,
// jet-track correlations,triggered jet shapes and
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
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TCanvas.h"
#include "TRandom3.h"
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
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskFastEmbedding.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetAntenna.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetAntenna)

AliAnalysisTaskJetAntenna::AliAnalysisTaskJetAntenna() :
AliAnalysisTaskSE(),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
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
fFilterMaskBestPt(0),
fFilterType(0),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fRequireITSRefit(0),
fApplySharedClusterCut(0),
fTrackTypeRec(kTrackUndef),
fRPAngle(0),
fNRPBins(50),
fSemigoodCorrect(0),
fDoMatching(kFALSE),
fHolePos(4.71),
fHoleWidth(0.2),
fCutTM(0.15),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fNevents(0),
fTindex(0),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2JetEntries(0x0),
fh2Circularity(0x0),
fh2JetAxisPhi(0x0),
fhnJetTM(0x0)
{
   // default Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskJetAntenna::AliAnalysisTaskJetAntenna(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
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
fFilterMaskBestPt(0),
fFilterType(0),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fRequireITSRefit(0),
fApplySharedClusterCut(0),
fTrackTypeRec(kTrackUndef),
fRPAngle(0),
fNRPBins(50),
fSemigoodCorrect(0),
fDoMatching(kFALSE),
fHolePos(4.71),
fHoleWidth(0.2),
fCutTM(0.15),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fNevents(0),
fTindex(0),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2JetEntries(0x0),
fh2Circularity(0x0),
fh2JetAxisPhi(0x0),
fhnJetTM(0x0)
 {
   // Constructor


   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;

   DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetAntenna::~AliAnalysisTaskJetAntenna()
{
   delete fListJets[0];
   delete fListJets[1];
}

void AliAnalysisTaskJetAntenna::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetAntenna::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetAntenna::UserCreateOutputObjects()
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
  fOutputList->Add(fHistEvtSelection);
  fh2JetEntries=new TH2F("JetEntries","",150,0,150,10,-0.5,9.5);
  fOutputList->Add(fh2JetEntries);
  fh2Circularity=new TH2F("Circcularity","",10,0,1,150,0,150);
  fOutputList->Add(fh2Circularity);
  fh2JetAxisPhi=new TH2F("JetAxisSmearPhi","",9,0,TMath::Pi(),10,-0.5,9.5);
  fOutputList->Add(fh2JetAxisPhi);
  

  Int_t nbinsJet[10]={3,9,75,9,36,5,7,10,50,2};
  Double_t binlowJet[10]= {0,0, 0, 0,-0.5*TMath::Pi(),0,0,-0.5,0,0};
  Double_t binupJet[10]= {100,0.9, 150,150,1.5*TMath::Pi(),1,150,9.5,200,2};
  fhnJetTM = new THnSparseF("fhnJetTM", "fhnJetTM; cent;dr;pt_jet;pt_track;phi;circ;ptrue;pthard;ptlead;isemebed",10,nbinsJet,binlowJet,binupJet);
 
 Double_t *xPt3=new Double_t[10];
  xPt3[0] = 0.;
  for(Int_t i = 1;i<=9;i++){
    if(xPt3[i-1]<2)xPt3[i] = xPt3[i-1] + 0.4; // 1 - 5
    else if(xPt3[i-1]<11)xPt3[i] = xPt3[i-1] + 3; // 5 - 12
    else xPt3[i] = xPt3[i-1] + 150.; // 18
  }
  fhnJetTM->SetBinEdges(3,xPt3);

  // Double_t *xPt2=new Double_t[10];
  // xPt2[0] = 0.;
  // xPt2[1]=20;
  // xPt2[2]=40;
  // xPt2[3]=60;
  // xPt2[4]=80;
  // xPt2[5]=100;
  // xPt2[6]=120;
  // xPt2[7]=150; 
  
  //fhnJetTM->SetBinEdges(2,xPt2);
 

  Double_t *xPt4=new Double_t[4];
  xPt4[0] = 0.;
  xPt4[1]=10;
  xPt4[2]=30;
  xPt4[3]=50;

  
  fhnJetTM->SetBinEdges(0,xPt4);

  Double_t *xPt5=new Double_t[10];
  xPt5[0] = 0.;
  xPt5[1]=20;
  xPt5[2]=40;
  xPt5[3]=60;
  xPt5[4]=80;
  xPt5[5]=100;
  xPt5[6]=120;
  xPt5[7]=150; 
  
  fhnJetTM->SetBinEdges(6,xPt5);


   fOutputList->Add(fhnJetTM);
   delete [] xPt3;
   // delete [] xPt2;
   delete [] xPt4;
   delete [] xPt5;
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

void AliAnalysisTaskJetAntenna::UserExec(Option_t *)
{


  if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
    AliError("Jet branch name not set.");
    return;
  }

  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    AliError("ESD not available");
    fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  }
  fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());

  static AliAODEvent* aod = 0;
  // take all other information from the aod we take the tracks from
  if(!aod){
    if(!fESD)aod = fAODIn;
    else aod = fAODOut;}



  if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }
  }


  // -- event selection --
  fHistEvtSelection->Fill(1); // number of events before event selection


  Bool_t selected=kTRUE;
  selected = AliAnalysisHelperJetTasks::Selected();
  if(!selected){
    // no selection by the service task, we continue
    PostData(1,fOutputList);
    return;}



  // physics selection: this is now redundant, all should appear as accepted after service task selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  std::cout<<inputHandler->IsEventSelected()<<" "<<fOfflineTrgMask<<std::endl;
  if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
    if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
    fHistEvtSelection->Fill(2);
    PostData(1, fOutputList);
    return;
  }



  // vertex selection
  if(!aod){
    if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
    fHistEvtSelection->Fill(3);
    PostData(1, fOutputList);
  }
  AliAODVertex* primVtx = aod->GetPrimaryVertex();

  if(!primVtx){
    if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
    fHistEvtSelection->Fill(3);
    PostData(1, fOutputList);
    return;
  }

  Int_t nTracksPrim = primVtx->GetNContributors();
  if ((nTracksPrim < fMinContribVtx) ||
      (primVtx->GetZ() < fVtxZMin) ||
      (primVtx->GetZ() > fVtxZMax) ){
    if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
    fHistEvtSelection->Fill(3);
    PostData(1, fOutputList);
    return;
  }



  // centrality selection
  AliCentrality *cent = 0x0;
  Double_t centValue = 0.;
  if(fIsPbPb){
    if(fESD) {cent = fESD->GetCentrality();
      if(cent) centValue = cent->GetCentralityPercentile("V0M");}
    else     centValue=((AliVAODHeader*)aod->GetHeader())->GetCentrality();

    if(fDebug) printf("centrality: %f\n", centValue);
    if (centValue < fCentMin || centValue > fCentMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
    }}


  fHistEvtSelection->Fill(0);
  // accepted events
  // -- end event selection --
  // pt hard
  Double_t pthard=0;
  Double_t pthardbin=0;
  if(fDoMatching){
  pthard = AliAnalysisTaskFastEmbedding::GetPtHard();
  pthardbin = GetPtHardBin(pthard);}
  // get background
  AliAODJetEventBackground* externalBackground = 0;
  if(fAODOut&&!externalBackground&&fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }
  if(fAODExtension&&!externalBackground&&fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }

  if(fAODIn&&!externalBackground&&fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }

  Float_t rho = 0;


  if(externalBackground)rho = externalBackground->GetBackground(0);

  // fetch jets
  TClonesArray *aodJets[2];
  aodJets[0]=0;
  if(fAODOut&&!aodJets[0]){
    aodJets[0] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[0].Data()));
    aodJets[1] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[1].Data()));
  }
  if(fAODExtension && !aodJets[0]){
    aodJets[0] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[0].Data()));
    aodJets[1] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[1].Data()));
  }
  if(fAODIn&&!aodJets[0]){
    aodJets[0] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[0].Data()));
    aodJets[1] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[1].Data()));
  }



  Int_t nT=0;
  TList ParticleList;
  if(!fDoMatching) nT = GetListOfTracks(&ParticleList);
  if(fDoMatching)nT=GetListOfTracksExtra(&ParticleList);
  if(nT<0){
    PostData(1, fOutputList);
    return;
  }

  for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
    fListJets[iJetType]->Clear();
    if (!aodJets[iJetType]) continue;
    if(fDebug) Printf("%s: %d jets",fJetBranchName[iJetType].Data(),aodJets[iJetType]->GetEntriesFast());
    for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
      if (jet) fListJets[iJetType]->Add(jet);
    }
  }

// jet matching 
  static TArrayI aMatchIndex(fListJets[0]->GetEntries());
  static TArrayF aPtFraction(fListJets[0]->GetEntries());
  if(aMatchIndex.GetSize()<fListJets[0]->GetEntries()) aMatchIndex.Set(fListJets[0]->GetEntries());
  if(aPtFraction.GetSize()<fListJets[0]->GetEntries()) aPtFraction.Set(fListJets[0]->GetEntries());

  if(fDoMatching){
 

  // stores matched jets in 'aMatchIndex' and fraction of pT in 'aPtFraction'
  AliAnalysisHelperJetTasks::GetJetMatching(fListJets[0], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[0]->GetEntries()),
					    fListJets[1], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[1]->GetEntries()),
					    aMatchIndex, aPtFraction, fDebug, fMatchMaxDist, fIsPbPb?1:2);}

  



 
 // loop over matched jets
  Int_t ir = -1;  // index of matched reconstruced jet
  Float_t fraction = -1.;
   AliAODJet *jetmatched =0x0;
   
  
   
   
  
    for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){

    Double_t etabig=0;
    Double_t ptbig=0;
    Double_t areabig=0;
    Double_t phibig=0.;
    Double_t pxbig,pybig,pzbig;
    Double_t phitrue=0;
    Double_t etatrue=0;
    Double_t smearphi=0;
    Double_t ptrue=0;
    AliAODJet* jetbig = (AliAODJet*)(fListJets[0]->At(i));
    etabig  = jetbig->Eta();
    phibig  = jetbig->Phi();
    phitrue=phibig;
    etatrue=etabig;
    ptbig   = jetbig->Pt();
    if(ptbig==0) continue;
    areabig = jetbig->EffectiveAreaCharged();
    if(!fDoMatching) ptbig=ptbig-rho*areabig;
    if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue;

    if(fSemigoodCorrect){
      if((phibig>fHolePos-fHoleWidth) && (phibig<fHolePos+fHoleWidth)) continue;
    }
    pxbig=jetbig->Px();
    pybig=jetbig->Py();
    pzbig=jetbig->Pz();

   if(fDoMatching){
     Bool_t jetAccepted=kTRUE;
    ir = aMatchIndex[i];
    if(ir>=0) jetmatched = (AliAODJet*)(fListJets[1]->At(ir));
    else    continue;
    fraction  = aPtFraction[i];
    // minimum fraction required
    
    if(fraction<fJetPtFractionMin) jetAccepted = kFALSE;
    if(jetAccepted){
    // jet acceptance + minimum pT check
    if(jetmatched->Eta()>fJetEtaMax || jetmatched->Eta()<fJetEtaMin)  jetAccepted = kFALSE;
      
    }
    ptrue=ptbig;
    if(!jetAccepted) continue;
    
    etabig=jetmatched->Eta();
    phibig=jetmatched->Phi();
    pxbig=jetmatched->Px();
    pybig=jetmatched->Py();
    pzbig=jetmatched->Pz();
    ptbig=jetmatched->Pt()-rho*jetmatched->EffectiveAreaCharged();


   
      

    smearphi=RelativePhi(phitrue,jetmatched->Phi());
    smearphi=TMath::Abs(smearphi);
   

}




    TVector3  ppJ1(pxbig, pybig, pzbig);
    TVector3  ppJ3(- pxbig * pzbig, - pybig * pzbig, pxbig * pxbig + pybig * pybig);
    ppJ3.SetMag(1.);
    TVector3  ppJ2(-pybig, pxbig, 0);
    ppJ2.SetMag(1.);




    Float_t mxx    = 0.;
    Float_t myy    = 0.;
    Float_t mxy    = 0.;
    Int_t   nc     = 0;
    Float_t sump2  = 0.;
    Float_t ptmax=-10;
    for(int it = 0;it<ParticleList.GetEntries();++it){
      AliVParticle *track = (AliVParticle*)ParticleList.At(it);
      TVector3 pp(track->Px(), track->Py(), track->Pz());
      Float_t phi = track->Phi();
      Float_t eta = track->Eta();
      Float_t pt  = track->Pt();
      Float_t deta = eta - etabig;
      Float_t dphi = RelativePhi(phi,phibig);
      if(TMath::Abs(dphi)>=0.5*TMath::Pi()) continue;
      Float_t r = TMath::Sqrt(dphi * dphi + deta * deta);
     
       if (r < 0.4 && pt>fCutTM) {
	//longitudinal and perpendicular component of the track pT in the
	//local frame
        if(pt>ptmax) ptmax=pt;
	TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
	TVector3 pPerp = pp - pLong;
	//projection onto the two perpendicular vectors defined above
	Float_t ppjX = pPerp.Dot(ppJ2);
	Float_t ppjY = pPerp.Dot(ppJ3);
	Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
	//components of the 2D symmetrical sphericity matrix
	mxx += (ppjX * ppjX / ppjT);
	myy += (ppjY * ppjY / ppjT);
	mxy += (ppjX * ppjY / ppjT);
	nc++;
	sump2 += ppjT;}
      
      if(nc<2) continue;

    } // 1st Track Loop


    // Sphericity Matrix
    const Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
    TMatrixDSym m0(2,ele);
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    Int_t jev = 0;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    TVector2 evec(evec0[0], evec0[1]);
    Float_t circ=0;
    if(jev==1) circ=2*eval[0];
    if(jev==0) circ=2*eval[1];
    fh2Circularity->Fill(circ,ptbig);
    fh2JetEntries->Fill(ptbig,pthardbin);


    if(fDoMatching) fh2JetAxisPhi->Fill(smearphi,pthardbin);
    

    
    for (Int_t ip = 0; ip < ParticleList.GetEntries(); ip++) {
      AliVParticle *track = (AliVParticle*)ParticleList.At(ip);
      Float_t isembed=0.5;
      if(fDoMatching) if(ip<nT) isembed=1.5;
      TVector3 pp(track->Px(), track->Py(), track->Pz());
      Float_t phi = track->Phi();
      Float_t eta = track->Eta();
      Float_t pt  = track->Pt();
    
      Float_t deta = eta - etabig;
      Float_t dphi = RelativePhi(phi,phibig);
      if(TMath::Abs(dphi)>=0.5*TMath::Pi()) continue;

      Float_t dRR = TMath::Sqrt(dphi * dphi + deta * deta);
      TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
      TVector3 pPerp = pp - pLong;
      Float_t ppjX = pPerp.Dot(ppJ2);
      Float_t ppjY = pPerp.Dot(ppJ3);
      TVector2 vr(ppjX, ppjY) ;
      //and this is the angle between the particle and the TM axis.
      float phistr=evec.Phi()-vr.Phi();
      if(phistr>2*TMath::Pi()) phistr -= 2*TMath::Pi();
      if(phistr<-2*TMath::Pi()) phistr += 2*TMath::Pi();
      if(phistr<-0.5*TMath::Pi()) phistr += 2*TMath::Pi();
      if(phistr>1.5*TMath::Pi()) phistr -= 2*TMath::Pi();

      double jetEntries[10] = {centValue,dRR,ptbig,pt,phistr,circ,ptrue,pthardbin,ptmax,isembed};
      fhnJetTM->Fill(jetEntries);

    } // 2nd Track loop
  }//jet loop

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetAntenna::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
}


Int_t  AliAnalysisTaskJetAntenna::GetListOfTracks(TList *list){

  Int_t iCount = 0;
  AliAODEvent *aod = 0;

  if(!fESD)aod = fAODIn;
  else aod = fAODOut;
  if(!aod)return 0;


    for(int it = 0;it < aod->GetNumberOfTracks();++it){
    AliAODTrack *tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) AliFatal("Not a standard AOD");
    Bool_t bGood = false;
    if(fFilterType == 0)bGood = true;
    else if(fFilterType == 1)bGood = tr->IsHybridTPCConstrainedGlobal();
    else if(fFilterType == 2)bGood = tr->IsHybridGlobalConstrainedGlobal();
    if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
    if(fRequireITSRefit==1){if((tr->GetStatus()&AliESDtrack::kITSrefit)==0)continue;}
    if(bGood==false) continue;
    if (fApplySharedClusterCut) {
      Double_t frac = Double_t(tr->GetTPCnclsS()) /Double_t(tr->GetTPCncls());
      if (frac > 0.4) continue;
    }
    if(TMath::Abs(tr->Eta())>0.9)continue;
    if(tr->Pt()<0.15)continue;
    list->Add(tr);
    iCount++;
    }
  return iCount;
}




Int_t  AliAnalysisTaskJetAntenna::GetListOfTracksExtra(TList *list){

  Int_t iCount = 0;
  Int_t nEmbed=0;
  AliAODEvent *aod = 0;

  if(!fESD)aod = fAODIn;
  else aod = fAODOut;
  if(!aod)return 0;
 

      TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(aod->FindListObject("aodExtraTracks"));
      if(!aodExtraTracks)return iCount;
      for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
	AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
	if (!track) continue;
	AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*> (track);
	if(!trackAOD)continue;
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = trackAOD->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = trackAOD->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!trackAOD->TestFilterBit(fFilterMask)||(!bGood))))continue;
        if(fRequireITSRefit){if((trackAOD->GetStatus()&AliESDtrack::kITSrefit)==0)continue;}
         if (fApplySharedClusterCut) {
           Double_t frac = Double_t(trackAOD->GetTPCnclsS()) /Double_t(trackAOD->GetTPCncls());
           if (frac > 0.4) continue;
	 }

	if(TMath::Abs(trackAOD->Eta())>0.9) continue;
	if(trackAOD->Pt()<0.15) continue;
	if(fDebug) printf("pt extra track %.2f \n", trackAOD->Pt());
	list->Add(trackAOD);
	iCount++;
        }

      nEmbed=iCount-1;




    for(int it = 0;it < aod->GetNumberOfTracks();++it){
    AliAODTrack *tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) AliFatal("Not a standard AOD");
    Bool_t bGood = false;
    if(fFilterType == 0)bGood = true;
    else if(fFilterType == 1)bGood = tr->IsHybridTPCConstrainedGlobal();
    else if(fFilterType == 2)bGood = tr->IsHybridGlobalConstrainedGlobal();
    if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
    if(fRequireITSRefit==1){if((tr->GetStatus()&AliESDtrack::kITSrefit)==0)continue;}
    if(bGood==false) continue;
    if (fApplySharedClusterCut) {
      Double_t frac = Double_t(tr->GetTPCnclsS()) /Double_t(tr->GetTPCncls());
      if (frac > 0.4) continue;
    }
    if(TMath::Abs(tr->Eta())>0.9)continue;
    if(tr->Pt()<0.15)continue;
    list->Add(tr);
    iCount++;
    } 
      
     return nEmbed;
}

Double_t AliAnalysisTaskJetAntenna::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}

Int_t AliAnalysisTaskJetAntenna::GetPhiBin(Double_t phi)
{
    Int_t phibin=-1;
    if(!(TMath::Abs(phi)<=2*TMath::Pi())){AliError("phi w.r.t. RP out of defined range");return -1;}
    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));
    phibin=Int_t(fNRPBins*phiwrtrp/(0.5*TMath::Pi()));
    if(phibin<0||phibin>=fNRPBins){AliError("Phi Bin not defined");}
    return phibin;
}

Int_t AliAnalysisTaskJetAntenna::GetPtHardBin(Double_t ptHard){

   const Int_t nBins = 10;
  Double_t binLimits[nBins] = { 5., 11., 21., 36., 57., 84., 117., 156., 200., 249. }; // lower limits
   
  Int_t bin = -1;
  while(bin<nBins-1 && binLimits[bin+1]<ptHard){
    bin++;
  }
   
  return bin;
}
