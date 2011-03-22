
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
//
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "AliRunLoader.h"
#include "AliStack.h"



#include <TPDGCode.h>
#include <TStyle.h>
#include "TLinearFitter.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"



#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliMathBase.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliTrackPointArray.h"
#include "TCint.h"
#include "AliTPCcalibV0.h"
#include "AliV0.h"
#include "TRandom.h"
#include "TTreeStream.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCorrection.h"
#include "AliGRPObject.h"
#include "AliTPCTransform.h"





ClassImp(AliTPCcalibV0)


AliTPCcalibV0::AliTPCcalibV0() : 
   AliTPCcalibBase(),
   fV0Tree(0),
   fHPTTree(0),
   fStack(0),
   fESD(0),
   fPdg(0),
   fParticles(0),
   fV0s(0),
   fGammas(0)
{
  
}   
AliTPCcalibV0::AliTPCcalibV0(const Text_t *name, const Text_t *title):
   AliTPCcalibBase(),
   fV0Tree(0),
   fHPTTree(0),
   fStack(0),
   fESD(0),
   fPdg(0),
   fParticles(0),
   fV0s(0),
   fGammas(0)
{
  fPdg = new TDatabasePDG;       
  // create output histograms
  SetName(name);
  SetTitle(title);
}   

AliTPCcalibV0::~AliTPCcalibV0(){
  //
  //
  //
  delete fV0Tree;
  delete fHPTTree;
}





void  AliTPCcalibV0::ProcessESD(AliESDEvent *esd, AliStack *stack){
  //
  //
  //
  fESD = esd;
  AliKFParticle::SetField(esd->GetMagneticField());
  if (TMath::Abs(AliTracker::GetBz())<1) return;  
  DumpToTree(esd);
  DumpToTreeHPT(esd);
  //
  if (stack) {
    MakeV0s();
    fStack = stack;
    MakeMC();
  }else{
    fStack =0;
  }
}

void  AliTPCcalibV0::DumpToTreeHPT(AliESDEvent *esd){
  //
  // Dump V0s fith full firend information to the 
  // 
  if (TMath::Abs(AliTracker::GetBz())<1) return;
  const Int_t kMinCluster=110;
  const Float_t kMinPt   =3.;
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(esd->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    Printf("ERROR: esdFriend not available");
    return;
  }
  //
  Int_t ntracks=esd->GetNumberOfTracks();
  for (Int_t i=0;i<ntracks;++i) {
    Bool_t isOK=kFALSE;
    AliESDtrack *track = esd->GetTrack(i);
    if (track->GetTPCncls()<kMinCluster) continue;
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(i);
    if (!friendTrack) continue;
    if (TMath::Abs(AliTracker::GetBz())>1){ // cut on momenta if measured
      if (track->Pt()>kMinPt) isOK=kTRUE;
    }
    if (TMath::Abs(AliTracker::GetBz())<1){  // require primary track for the B field OFF data
      Bool_t isAccepted=kTRUE;
      if (!track->IsOn(AliESDtrack::kITSrefit)) isAccepted=kFALSE;
      if (!track->IsOn(AliESDtrack::kTPCrefit)) isAccepted=kFALSE;
      if (!track->IsOn(AliESDtrack::kTOFout))   isAccepted=kFALSE;
      Float_t dvertex[2],cvertex[3]; 
      track->GetImpactParametersTPC(dvertex,cvertex);
      if (TMath::Abs(dvertex[0]/TMath::Sqrt(cvertex[0]+0.01))>20) isAccepted=kFALSE;
      if (TMath::Abs(dvertex[1]/TMath::Sqrt(TMath::Abs(cvertex[2]+0.01)))>20) isAccepted=kFALSE;
      track->GetImpactParameters(dvertex,cvertex);
      if (TMath::Abs(dvertex[0]/TMath::Sqrt(cvertex[0]+0.01))>10) isAccepted=kFALSE;
      if (TMath::Abs(dvertex[1]/TMath::Sqrt(TMath::Abs(cvertex[2]+0.01)))>10) isAccepted=kFALSE;
      if (!isAccepted) isOK=kFALSE;
    }
    
    if (!isOK) continue;
    //

    TObject *calibObject;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (!seed) continue;
      if (!fHPTTree) {
      fHPTTree = new TTree("HPT","HPT");
      fHPTTree->SetDirectory(0);
    }
    if (fHPTTree->GetEntries()==0){
      //
      fHPTTree->SetDirectory(0);
      fHPTTree->Branch("t.",&track);
      fHPTTree->Branch("ft.",&friendTrack);
      fHPTTree->Branch("s.",&seed);
    }else{
      fHPTTree->SetBranchAddress("t.",&track);
      fHPTTree->SetBranchAddress("ft.",&friendTrack);
      fHPTTree->SetBranchAddress("s.",&seed);
    }
    fHPTTree->Fill();
    //
  }
}



void  AliTPCcalibV0::DumpToTree(AliESDEvent *esd){
  //
  // Dump V0s fith full firend information to the 
  // 
  Int_t nV0s  = fESD->GetNumberOfV0s();
  const Int_t kMinCluster=110;
  const Double_t kDownscale=0.01;
  const Float_t kMinPt   =1.0;
  const Float_t kMinMinPt   =0.7;
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(esd->FindListObject("AliESDfriend"));
  if (!esdFriend) {
    Printf("ERROR: esdFriend not available");
    return;
  }
  //
  
  for (Int_t ivertex=0; ivertex<nV0s; ivertex++){
    Bool_t isOK=kFALSE;
    AliESDv0 * v0 = (AliESDv0*) esd->GetV0(ivertex);
    AliESDtrack * track0 = fESD->GetTrack(v0->GetIndex(0)); // negative track
    AliESDtrack * track1 = fESD->GetTrack(v0->GetIndex(1)); // positive track 
    if (track0->GetTPCNcls()<kMinCluster) continue;
    if (track0->GetKinkIndex(0)>0) continue;    
    if (track1->GetTPCNcls()<kMinCluster) continue;
    if (track1->GetKinkIndex(0)>0) continue;
    if (v0->GetOnFlyStatus()==kFALSE) continue;
    //
    if (TMath::Min(track0->Pt(),track1->Pt())<kMinMinPt) continue;
    //
    //
    if (TMath::Max(track0->Pt(),track1->Pt())>kMinPt) isOK=kTRUE;
    if (gRandom->Rndm()<kDownscale) isOK=kTRUE;  
    if (!isOK) continue;
    //
    AliESDfriendTrack *ftrack0 = esdFriend->GetTrack(v0->GetIndex(0));
    if (!ftrack0) continue;
    AliESDfriendTrack *ftrack1 = esdFriend->GetTrack(v0->GetIndex(1));
    if (!ftrack1) continue;
    //
    TObject *calibObject;
    AliTPCseed *seed0 = 0;
    AliTPCseed *seed1 = 0;
    for (Int_t l=0;(calibObject=ftrack0->GetCalibObject(l));++l) {
      if ((seed0=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    for (Int_t l=0;(calibObject=ftrack1->GetCalibObject(l));++l) {
      if ((seed1=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
    if (!seed0) continue;
    if (!seed1) continue;
    AliExternalTrackParam * paramIn0 = (AliExternalTrackParam *)track0->GetInnerParam();
    AliExternalTrackParam * paramIn1 = (AliExternalTrackParam *)track1->GetInnerParam();
    if (!paramIn0) continue;
    if (!paramIn1) continue;
    //
    //
    if (!fV0Tree) {
      fV0Tree = new TTree("V0s","V0s");
      fV0Tree->SetDirectory(0);
    }
    if (fV0Tree->GetEntries()==0){
      //
      fV0Tree->SetDirectory(0);
      fV0Tree->Branch("v0.",&v0);
      fV0Tree->Branch("t0.",&track0);
      fV0Tree->Branch("t1.",&track1);
      fV0Tree->Branch("ft0.",&ftrack0);
      fV0Tree->Branch("ft1.",&ftrack1);
      fV0Tree->Branch("s0.",&seed0);
      fV0Tree->Branch("s1.",&seed1);
    }else{
      fV0Tree->SetBranchAddress("v0.",&v0);
      fV0Tree->SetBranchAddress("t0.",&track0);
      fV0Tree->SetBranchAddress("t1.",&track1);
      fV0Tree->SetBranchAddress("ft0.",&ftrack0);
      fV0Tree->SetBranchAddress("ft1.",&ftrack1);
      fV0Tree->SetBranchAddress("s0.",&seed0);
      fV0Tree->SetBranchAddress("s1.",&seed1);
    }
    fV0Tree->Fill();
  }
}


Long64_t AliTPCcalibV0::Merge(TCollection *const li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibV0* cal = 0;

  while ((cal = (AliTPCcalibV0*)iter->Next())) {
    if (cal->fV0Tree){
      if (!fV0Tree) {
	fV0Tree = new TTree("V0s","V0s");
	fV0Tree->SetDirectory(0);
      }
      if (cal->fV0Tree->GetEntries()>0) AliTPCcalibV0::AddTree(cal->fV0Tree);
      if (cal->fHPTTree->GetEntries()>0) AliTPCcalibV0::AddTreeHPT(cal->fHPTTree);
    }    
  }
  return 0;
}


void AliTPCcalibV0::AddTree(TTree * treeInput){
  //
  // Add the content of tree: 
  // Notice automatic copy of tree in ROOT does not work for such complicated tree
  //  
  AliESDv0 * v0 = new AliESDv0;
  Double_t kMinPt=0.8;
  AliESDtrack * track0 = 0; // negative track
  AliESDtrack * track1 = 0; // positive track 
  AliESDfriendTrack *ftrack0 = 0;
  AliESDfriendTrack *ftrack1 = 0;
  AliTPCseed *seed0 = 0;
  AliTPCseed *seed1 = 0;
  treeInput->SetBranchStatus("ft0.",kFALSE);
  treeInput->SetBranchStatus("ft1.",kFALSE);
  TDatabasePDG pdg;
  Double_t massK0= pdg.GetParticle("K0")->Mass();
  Double_t massLambda= pdg.GetParticle("Lambda0")->Mass();

  Int_t entries= treeInput->GetEntries();
  for (Int_t i=0; i<entries; i++){
    treeInput->SetBranchAddress("v0.",&v0);
    treeInput->SetBranchAddress("t0.",&track0);
    treeInput->SetBranchAddress("t1.",&track1);
    treeInput->SetBranchAddress("ft0.",&ftrack0);
    treeInput->SetBranchAddress("ft1.",&ftrack1);
    treeInput->SetBranchAddress("s0.",&seed0);
    treeInput->SetBranchAddress("s1.",&seed1);
    if (fV0Tree->GetEntries()==0){
      fV0Tree->SetDirectory(0);
      fV0Tree->Branch("v0.",&v0);
      fV0Tree->Branch("t0.",&track0);
      fV0Tree->Branch("t1.",&track1);
      fV0Tree->Branch("ft0.",&ftrack0);
      fV0Tree->Branch("ft1.",&ftrack1);
      fV0Tree->Branch("s0.",&seed0);
      fV0Tree->Branch("s1.",&seed1);
    }else{
      fV0Tree->SetBranchAddress("v0.",&v0);
      fV0Tree->SetBranchAddress("t0.",&track0);
      fV0Tree->SetBranchAddress("t1.",&track1);
      fV0Tree->SetBranchAddress("ft0.",&ftrack0);
      fV0Tree->SetBranchAddress("ft1.",&ftrack1);
      fV0Tree->SetBranchAddress("s0.",&seed0);
      fV0Tree->SetBranchAddress("s1.",&seed1);
    }
    //
    treeInput->GetEntry(i);
    //ftrack0->GetCalibContainer()->SetOwner(kTRUE);
    //ftrack1->GetCalibContainer()->SetOwner(kTRUE);
    Bool_t isOK=kTRUE;
    if (v0->GetOnFlyStatus()==kFALSE) isOK=kFALSE;
    if (track0->GetTPCncls()<100) isOK=kFALSE;
    if (track1->GetTPCncls()<100) isOK=kFALSE;    
    if (TMath::Min(seed0->Pt(),seed1->Pt())<kMinPt) isOK=kFALSE;
    if (TMath::Min(track0->Pt(),track1->Pt())<kMinPt) isOK=kFALSE;
    Bool_t isV0=kFALSE;    
    if (TMath::Abs(v0->GetEffMass(2,2)-massK0)<0.05)     isV0=kTRUE;
    if (TMath::Abs(v0->GetEffMass(4,2)-massLambda)<0.05) isV0=kTRUE; 
    if (TMath::Abs(v0->GetEffMass(2,4)-massLambda)<0.05) isV0=kTRUE;
    if (TMath::Abs(v0->GetEffMass(0,0))<0.02) isV0=kFALSE; //reject electrons
    if (!isV0) isOK=kFALSE;
    if (isOK) fV0Tree->Fill();
    delete v0;
    delete track0;
    delete track1;
    delete ftrack0;
    delete ftrack1;
    delete seed0;
    delete seed1;
    v0=0;
    track0=0;
    track1=0;
    ftrack0=0;
    ftrack1=0;
    seed0=0;
    seed1=0;
  }
}

void AliTPCcalibV0::AddTreeHPT(TTree * treeInput){
  //
  // Add the content of tree: 
  // Notice automatic copy of tree in ROOT does not work for such complicated tree
  //  
  AliESDtrack *track = 0;
  AliESDfriendTrack *friendTrack = 0;
  AliTPCseed *seed = 0;
  if (!treeInput) return;
  if (treeInput->GetEntries()==0) return;
  //
  Int_t entries= treeInput->GetEntries();  
  //
  for (Int_t i=0; i<entries; i++){
    track=0;
    friendTrack=0;
    seed=0;
    //
    treeInput->SetBranchAddress("t.",&track);
    treeInput->SetBranchAddress("ft.",&friendTrack);
    treeInput->SetBranchAddress("s.",&seed);
    treeInput->GetEntry(i);
    //
    if (fHPTTree->GetEntries()==0){
      fHPTTree->SetDirectory(0);
      fHPTTree->Branch("t.",&track);
      fHPTTree->Branch("ft.",&friendTrack);
      fHPTTree->Branch("s.",&seed);
    }else{
      fHPTTree->SetBranchAddress("t.",&track);
      fHPTTree->SetBranchAddress("ft.",&friendTrack);
      fHPTTree->SetBranchAddress("s.",&seed);
    }    
    Bool_t isOK=kTRUE;
    if (!track->IsOn(AliESDtrack::kITSrefit)) isOK=kFALSE;
    if (!track->IsOn(AliESDtrack::kTOFout)) isOK=kFALSE;
    if (isOK) fHPTTree->Fill();
    //
    delete track;
    delete friendTrack;
    delete seed;
  }
}


void AliTPCcalibV0::MakeMC(){
  //
  // MC comparison 
  //    1. Select interesting particles
  //    2. Assign the recosntructed particles 
  //
  //1. Select interesting particles
  const Float_t kMinP   = 0.2;
  const Float_t kMinPt  = 0.1;
  const Float_t kMaxR   = 0.5;
  const Float_t kMaxTan = 1.2;
  const Float_t kMaxRad = 150;
  //
  if (!fParticles) fParticles = new TObjArray;
  TParticle *part=0;
  //  
  Int_t entries = fStack->GetNtrack();
  for (Int_t ipart=0; ipart<entries; ipart++){
    part = fStack->Particle(ipart);
    if (!part) continue;
    if (part->P()<kMinP) continue;
    if (part->R()>kMaxR) continue;
    if (TMath::Abs(TMath::Tan(part->Theta()-TMath::Pi()*0.5))>kMaxTan) continue;
    Bool_t isInteresting = kFALSE;
    if (part->GetPdgCode()==22) isInteresting =kTRUE;
    if (part->GetPdgCode()==310) isInteresting =kTRUE;
    if (part->GetPdgCode()==111) isInteresting =kTRUE;
    if (TMath::Abs(part->GetPdgCode()==3122)) isInteresting =kTRUE;

    //
    if (!isInteresting) continue;    
    fParticles->AddLast(new TParticle(*part));
  }
  if (fParticles->GetEntries()<1) {
    return;
  }
  //
  //
  //
  Int_t sentries=fParticles->GetEntries();;
  for (Int_t ipart=0; ipart<sentries; ipart++){
    part = (TParticle*)fParticles->At(ipart);
    TParticle *p0 = 0;
    TParticle *p1 = 0;

    Int_t nold =0;
    Int_t nnew =0;
    Int_t id0  = part->GetDaughter(0);
    Int_t id1  = part->GetDaughter(1);    
    if (id0>=fStack->GetNtrack() ) id0*=-1;
    if (id1>=fStack->GetNtrack() ) id1*=-1;
    Bool_t findable = kTRUE;
    if (id0<0 || id1<0) findable = kFALSE;
    Int_t charge =0; 
    if (findable){
      p0 = fStack->Particle(id0);
      if (p0->R()>kMaxRad) findable = kFALSE;
      if (p0->Pt()<kMinPt) findable = kFALSE;
      if (p0->Vz()>250) findable= kFALSE;
      if (TMath::Abs(TMath::Tan(p0->Theta()-TMath::Pi()*0.5))>2) findable=kFALSE;
      if (fPdg->GetParticle(p0->GetPdgCode())==0) findable =kFALSE;
      else
	if (fPdg->GetParticle(p0->GetPdgCode())->Charge()==0) charge++;
	  
      p1 = fStack->Particle(id1);
      if (p1->R()>kMaxRad) findable = kFALSE;
      if (p1->Pt()<kMinPt) findable = kFALSE;
      if (TMath::Abs(p1->Vz())>250) findable= kFALSE;
      if (TMath::Abs(TMath::Tan(p1->Theta()-TMath::Pi()*0.5))>2) findable=kFALSE;
      if (fPdg->GetParticle(p1->GetPdgCode())==0) findable = kFALSE;
      else
	if (fPdg->GetParticle(p1->GetPdgCode())->Charge()==0) charge++;
			  
    }
    //   (*fDebugStream)<<"MC0"<<
    //       "P.="<<part<<
    //       "findable="<<findable<<
    //       "id0="<<id0<<
    //       "id1="<<id1<<
    //       "\n";
    if (!findable) continue;
    Float_t minpt = TMath::Min(p0->Pt(), p1->Pt());
    Int_t type=-1;
    
    //
    // 
    AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
    for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
      AliESDv0 * v0 = fESD->GetV0(ivertex);
      // select coresponding track
      AliESDtrack * trackN = fESD->GetTrack(v0->GetIndex(0));
      if (TMath::Abs(trackN->GetLabel())!=id0 && TMath::Abs(trackN->GetLabel())!=id1) continue;
      AliESDtrack * trackP = fESD->GetTrack(v0->GetIndex(1));
      if (TMath::Abs(trackP->GetLabel())!=id0 && TMath::Abs(trackP->GetLabel())!=id1) continue;
      TParticle *pn = fStack->Particle(TMath::Abs(trackN->GetLabel()));
      TParticle *pp = fStack->Particle(TMath::Abs(trackP->GetLabel()));
      //
      //
      if ( v0->GetOnFlyStatus()) nnew++;
      if (!v0->GetOnFlyStatus()) nold++;
      if (part->GetPdgCode()==22 && TMath::Abs(pn->GetPdgCode())==11 &&  TMath::Abs(pp->GetPdgCode())==11) 
	type =1;
      if (part->GetPdgCode()==310 && TMath::Abs(pn->GetPdgCode())==211 &&  TMath::Abs(pp->GetPdgCode())==211) 
	type =0;
      if (part->GetPdgCode()==3122){
	if (TMath::Abs(pn->GetPdgCode())==210 ) type=2;
	else type=3;
      }
      AliKFParticle *v0kf       = Fit(primVtx,v0,pn->GetPdgCode(),pp->GetPdgCode());
      v0kf->SetProductionVertex( primVtx );
      AliKFParticle *v0kfc       = Fit(primVtx,v0,pn->GetPdgCode(),pp->GetPdgCode());
      v0kfc->SetProductionVertex( primVtx );
      v0kfc->SetMassConstraint(fPdg->GetParticle(part->GetPdgCode())->Mass());
      Float_t chi2 = v0kf->GetChi2();
      Float_t chi2C = v0kf->GetChi2();
      //
      //
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"MCRC"<<
	  "P.="<<part<<
	  "type="<<type<<
	  "chi2="<<chi2<<
	  "chi2C="<<chi2C<<
	  "minpt="<<minpt<<
	  "id0="<<id0<<
	  "id1="<<id1<<
	  "Pn.="<<pn<<
	  "Pp.="<<pp<<
	  "tn.="<<trackN<<
	  "tp.="<<trackP<<
	  "nold.="<<nold<<
	  "nnew.="<<nnew<<
	  "v0.="<<v0<<
	  "v0kf.="<<v0kf<<
	  "v0kfc.="<<v0kfc<<
	  "\n";     
	delete v0kf;
	delete v0kfc; 
	//
      }
      
      if (cstream){
	(*cstream)<<"MC"<<
	  "P.="<<part<<
	  "charge="<<charge<<
	  "type="<<type<<
	  "minpt="<<minpt<<
	  "id0="<<id0<<
	  "id1="<<id1<<
	  "P0.="<<p0<<
	  "P1.="<<p1<<
	  "nold="<<nold<<
	  "nnew="<<nnew<<
	  "\n";
      }
    }
    fParticles->Delete(); 
  }
}

void AliTPCcalibV0::MakeFitTreeTrack(const TObjArray * corrArray, Double_t ptCut, Int_t /*run*/){
  //
  // Make a fit tree
  //
  // 0. Loop over selected tracks
  // 1. Loop over all transformation - refit the track with and without the
  //    transformtation
  // 2. Dump the matching paramaeters to the debugStremer
  //
  
  //Connect input
  const Int_t kMinNcl=120;
  TFile f("TPCV0Objects.root");
  AliTPCcalibV0 *v0TPC = (AliTPCcalibV0*) f.Get("v0TPC");
  TTree * treeInput = v0TPC->GetHPTTree();
  TTreeSRedirector *pcstream = new TTreeSRedirector("fitHPT.root");
  AliESDtrack *track = 0;
  AliESDfriendTrack *friendTrack = 0;
  AliTPCseed *seed = 0;
  if (!treeInput) return;
  if (treeInput->GetEntries()==0) return;
  //
  treeInput->SetBranchAddress("t.",&track);
  treeInput->SetBranchAddress("ft.",&friendTrack);
  treeInput->SetBranchAddress("s.",&seed);
  //
  Int_t ncorr=0;
  if (corrArray) ncorr = corrArray->GetEntries();
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
 //  AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();
//   AliGRPObject*  grp = AliTPCcalibDB::Instance()->GetGRP(run);
//   Double_t time=0.5*(grp->GetTimeStart() +grp->GetTimeEnd());
  //
  //
  //  
  Int_t ntracks= treeInput->GetEntries();
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    treeInput->GetEntry(itrack);
    if (!track) continue;
    if (seed->Pt()<ptCut) continue;
    if (track->Pt()<ptCut) continue;
    if (track->GetTPCncls()<kMinNcl) continue;
    //
    // Reapply transformation
    //
    for (Int_t irow=0; irow<159; irow++){
      AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
      if (cluster &&cluster->GetX()>10){
        Double_t x0[3]={cluster->GetRow(),cluster->GetPad(),cluster->GetTimeBin()};
        Int_t index0[1]={cluster->GetDetector()};
        transform->Transform(x0,index0,0,1);
        cluster->SetX(x0[0]);
        cluster->SetY(x0[1]);
        cluster->SetZ(x0[2]);
        //
      }
    }    
    //
    AliExternalTrackParam* paramInner=0;
    AliExternalTrackParam* paramOuter=0;
    AliExternalTrackParam* paramIO=0;
    Bool_t isOK=kTRUE;
    for (Int_t icorr=-1; icorr<ncorr; icorr++){
      //
      AliTPCCorrection *corr = 0;
      if (icorr>=0) corr = (AliTPCCorrection*)corrArray->At(icorr);
      AliExternalTrackParam * trackInner = RefitTrack(seed, corr,85,134,0.1);      
      AliExternalTrackParam * trackIO = RefitTrack(seed, corr,245,85,0.1);      
      AliExternalTrackParam * trackOuter = RefitTrack(seed, corr,245,134,0.1 ); 
      if (trackInner&&trackOuter&&trackIO){
	trackOuter->Rotate(trackInner->GetAlpha());
	trackOuter->PropagateTo(trackInner->GetX(),AliTracker::GetBz());
	if (icorr<0) {
	  paramInner=trackInner;
	  paramOuter=trackOuter;
	  paramIO=trackIO;
	  paramIO->Rotate(seed->GetAlpha());
	  paramIO->PropagateTo(seed->GetX(),AliTracker::GetBz());
	}
      }else{
	isOK=kFALSE;
      }
      
    }
    if (paramOuter&& paramInner) {
      //      Bool_t isOK=kTRUE;
      if (paramInner->GetSigmaY2()>0.01) isOK&=kFALSE;
      if (paramOuter->GetSigmaY2()>0.01) isOK&=kFALSE;
      if (paramInner->GetSigmaZ2()>0.01) isOK&=kFALSE;
      if (paramOuter->GetSigmaZ2()>0.01) isOK&=kFALSE;      
      (*pcstream)<<"fit"<<
	"s.="<<seed<<
	"io.="<<paramIO<<
	"pIn.="<<paramInner<<
	"pOut.="<<paramOuter;      
    }
    //
  }
  delete pcstream;
  /*
    .x ~/rootlogon.C
    Int_t run=117112;
    .x ../ConfigCalibTrain.C(run)
    .L ../AddTaskTPCCalib.C
    ConfigOCDB(run)
    TFile fexb("../../RegisterCorrectionExB.root");
    AliTPCComposedCorrection *cc=  (AliTPCComposedCorrection*) fexb.Get("ComposedExB");
    cc->Init();
    cc->Print("DA"); // Print used correction classes
    TObjArray *array = cc->GetCorrections()
    AliTPCcalibV0::MakeFitTreeTrack(array,2,run);
   
   */
}

void AliTPCcalibV0::MakeFitTreeV0(const TObjArray * corrArray, Double_t ptCut, Int_t run){
  //
  // Make a fit tree
  //
  // 0. Loop over selected tracks
  // 1. Loop over all transformation - refit the track with and without the
  //    transformtation
  // 2. Dump the matching paramaeters to the debugStremer
  //
  
  //Connect input
  TFile f("TPCV0Objects.root");
  AliTPCcalibV0 *v0TPC = (AliTPCcalibV0*) f.Get("v0TPC");
  TTree * treeInput = v0TPC->GetV0Tree();
  TTreeSRedirector *pcstream = new TTreeSRedirector("fitV0.root");
  AliESDv0 *v0 = 0;
  AliESDtrack *track0 = 0;
  AliESDfriendTrack *friendTrack0 = 0;
  AliTPCseed *seed0 = 0;
  AliTPCseed *s0 = 0;
  AliESDtrack *track1 = 0;
  AliESDfriendTrack *friendTrack1 = 0;
  AliTPCseed *s1 = 0;
  AliTPCseed *seed1 = 0;
  if (!treeInput) return;
  if (treeInput->GetEntries()==0) return;
  //
  treeInput->SetBranchAddress("v0.",&v0);
  treeInput->SetBranchAddress("t0.",&track0);
  treeInput->SetBranchAddress("ft0.",&friendTrack0);
  treeInput->SetBranchAddress("s0.",&s0);
  treeInput->SetBranchAddress("t1.",&track1);
  treeInput->SetBranchAddress("ft1.",&friendTrack1);
  treeInput->SetBranchAddress("s1.",&s1);
  //
  TDatabasePDG pdg;
  Int_t ncorr=0;
  if (corrArray) ncorr = corrArray->GetEntries();
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  Double_t massK0= pdg.GetParticle("K0")->Mass();
  Double_t massLambda= pdg.GetParticle("Lambda0")->Mass();
  Double_t massPion=pdg.GetParticle("pi+")->Mass();
  Double_t massProton=pdg.GetParticle("proton")->Mass();
  Int_t pdgPion=pdg.GetParticle("pi+")->PdgCode();
  Int_t pdgProton=pdg.GetParticle("proton")->PdgCode();
  Double_t rmass0=0;
  Double_t rmass1=0;
  Double_t massV0=0;
  Int_t    pdg0=0;
  Int_t    pdg1=0;
  //
  //
  //  
  Int_t nv0s= treeInput->GetEntries();
  for (Int_t iv0=0; iv0<nv0s; iv0++){
    Int_t  v0Type=0;
    Int_t isK0=0;
    Int_t isLambda=0;
    Int_t isAntiLambda=0;
    treeInput->GetEntry(iv0);
    if (TMath::Abs(v0->GetEffMass(2,2)-massK0)<0.03) {isK0=1; v0Type=1;} //select K0s    
    if (TMath::Abs(v0->GetEffMass(4,2)-massLambda)<0.01) {isLambda=1; v0Type=2;} //select Lambda   
    if (TMath::Abs(v0->GetEffMass(2,4)-massLambda)<0.01) {isAntiLambda=1;v0Type=3;} //select Anti Lambda
    if (isK0+isLambda+isAntiLambda!=1) continue;
    rmass0=massPion;
    rmass1=massPion;
    pdg0=pdgPion;
    pdg1=pdgPion;
    if (isLambda) {rmass0=massProton; pdg0=pdgProton;}
    if (isAntiLambda) {rmass1=massProton; pdg1=pdgProton;}
    massV0=massK0;
    if (isK0==0) massV0=massLambda;
    //
    if (!s0) continue;
    seed0=(s0->GetSign()>0)?s0:s1;
    seed1=(s0->GetSign()>0)?s1:s0;
    if (seed0->GetZ()*seed1->GetZ()<0) continue; //remove membrane crossed tracks
    if (seed0->Pt()<ptCut) continue;
    if (seed1->Pt()<ptCut) continue;
    //
    // Reapply transformation
    //
    for  (Int_t itype=0; itype<2; itype++){
      AliTPCseed * seed = (itype==0) ? seed0: seed1;      
      for (Int_t irow=0; irow<159; irow++){
	AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
	if (cluster &&cluster->GetX()>10){
	  Double_t x0[3]={cluster->GetRow(),cluster->GetPad(),cluster->GetTimeBin()};
	  Int_t index0[1]={cluster->GetDetector()};
	  transform->Transform(x0,index0,0,1);
	  cluster->SetX(x0[0]);
	  cluster->SetY(x0[1]);
	  cluster->SetZ(x0[2]);
	  //
	}
      }
    }   
    Bool_t isOK=kTRUE;
    Double_t radius = v0->GetRr();
    Double_t xyz[3];
    v0->GetXYZ(xyz[0],xyz[1],xyz[2]);
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    TObjArray arrayV0in(ncorr+1);
    TObjArray arrayV0io(ncorr+1);
    TObjArray arrayT0(ncorr+1);
    TObjArray arrayT1(ncorr+1);
    arrayV0in.SetOwner(kTRUE);
    arrayV0io.SetOwner(kTRUE);
    //
    for (Int_t icorr=-1; icorr<ncorr; icorr++){
      AliTPCCorrection *corr =0;
      if (icorr>=0) corr = (AliTPCCorrection*)corrArray->At(icorr);
      //
      AliExternalTrackParam * trackInner0 = RefitTrack(seed0, corr,160,85,rmass0);      
      AliExternalTrackParam * trackIO0    = RefitTrack(seed0, corr,245,85,rmass0);      
      AliExternalTrackParam * trackInner1 = RefitTrack(seed1, corr,160,85,rmass1);      
      AliExternalTrackParam * trackIO1    = RefitTrack(seed1, corr,245,85,rmass1);      
      if (!trackInner0) isOK=kFALSE;
      if (!trackInner1) isOK=kFALSE;
      if (!trackIO0)    isOK=kFALSE;
      if (!trackIO1)    isOK=kFALSE;
      if (isOK){
	if (!trackInner0->Rotate(alpha)) isOK=kFALSE;
	if (!trackInner1->Rotate(alpha)) isOK=kFALSE;
	if (!trackIO0->Rotate(alpha)) isOK=kFALSE;
	if (!trackIO1->Rotate(alpha)) isOK=kFALSE;
	//
	if (!AliTracker::PropagateTrackToBxByBz(trackInner0, radius, rmass0, 1, kFALSE)) isOK=kFALSE; 
	if (!AliTracker::PropagateTrackToBxByBz(trackInner1, radius, rmass1, 1, kFALSE)) isOK=kFALSE; 
	if (!AliTracker::PropagateTrackToBxByBz(trackIO0, radius, rmass0, 1, kFALSE)) isOK=kFALSE; 
	if (!AliTracker::PropagateTrackToBxByBz(trackIO1, radius, rmass1, 1, kFALSE)) isOK=kFALSE; 
	if (!isOK) continue;
	arrayT0.AddAt(trackIO0->Clone(),icorr+1);
	arrayT1.AddAt(trackIO1->Clone(),icorr+1);
	Int_t charge=TMath::Nint(trackIO0->GetSign());
	AliKFParticle pin0( *trackInner0,  pdg0*charge);
	AliKFParticle pin1( *trackInner1, -pdg1*charge);
	AliKFParticle pio0( *trackIO0,  pdg0*charge);
	AliKFParticle pio1( *trackIO1, -pdg1*charge);
	AliKFParticle v0in;
	AliKFParticle v0io;
	v0in+=pin0;
	v0in+=pin1;
	v0io+=pio0;
	v0io+=pio1;
	arrayV0in.AddAt(v0in.Clone(),icorr+1);
	arrayV0io.AddAt(v0io.Clone(),icorr+1);
      }
    }
    if (!isOK) continue;
    //
    //AliKFParticle* pin0= (AliKFParticle*)arrayV0in.At(0);
    AliKFParticle* pio0= (AliKFParticle*)arrayV0io.At(0);
    AliExternalTrackParam *param0=(AliExternalTrackParam *)arrayT0.At(0);
    AliExternalTrackParam *param1=(AliExternalTrackParam *)arrayT1.At(0);
    Double_t mass0=0, mass0E=0; 
    pio0->GetMass( mass0,mass0E);
    //
    Double_t mean=mass0-massV0;
    if (TMath::Abs(mean)>0.05) continue;
    Double_t mass[10000];
    //
    Int_t dtype=30;  // id for V0
    Int_t ptype=5;   // id for invariant mass
    //    Int_t id=TMath::Nint(100.*(param0->Pt()-param1->Pt())/(param0->Pt()+param1->Pt()));      // K0s V0 asymetry
    Int_t id=Int_t(1000.*(param0->Pt()-param1->Pt()));      // K0s V0 asymetry
    Double_t gx,gy,gz, px,py,pz;
    Double_t pt = v0->Pt();
    v0->GetXYZ(gx,gy,gz);
    v0->GetPxPyPz(px,py,pz);
    Double_t theta=pz/TMath::Sqrt(px*px+py*py);
    Double_t phi=TMath::ATan2(py,px);
    Double_t snp=0.5*(seed0->GetSnp()+seed1->GetSnp());
    Double_t sector=9*phi/TMath::Pi();
    Double_t dsec=sector-TMath::Nint(sector);
    Double_t refX=TMath::Sqrt(gx*gx+gy*gy);
    //Int_t nentries=v0Type;
    Double_t bz=AliTracker::GetBz();
    Double_t dRrec=0;
    (*pcstream)<<"fitDebug"<< 
      "id="<<id<<
      "v0.="<<v0<<
      "mean="<<mean<<
      "rms="<<mass0E<<
      "pio0.="<<pio0<<
      "p0.="<<param0<<
      "p1.="<<param1;
    (*pcstream)<<"fit"<<  // dump valus for fit
      "run="<<run<<       //run number
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type 30
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<mass0E<<       // rms
      "sector="<<sector<<
      "dsec="<<dsec<<
      //
      "refX="<<refX<<      // reference radius
      "gx="<<gx<<         // global position
      "gy="<<gy<<         // global position
      "gz="<<gz<<         // global position
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "pt="<<pt<<         // pt of the particle
      "id="<<id<<     //delta of the momenta      
      "entries="<<v0Type;//  type of the V0
    for (Int_t icorr=0; icorr<ncorr; icorr++){
      AliTPCCorrection *corr =0;
      if (icorr>=0) corr = (AliTPCCorrection*)corrArray->At(icorr);
      //      AliKFParticle* pin= (AliKFParticle*)arrayV0in.At(icorr+1);
      AliKFParticle* pio= (AliKFParticle*)arrayV0io.At(icorr+1);
      AliExternalTrackParam *par0=(AliExternalTrackParam *)arrayT0.At(icorr+1);
      AliExternalTrackParam *par1=(AliExternalTrackParam *)arrayT1.At(icorr+1);
      Double_t massE=0; 
      pio->GetMass( mass[icorr],massE);
      mass[icorr]-=mass0;
      (*pcstream)<<"fit"<< 
	Form("%s=",corr->GetName())<<mass[icorr];
      (*pcstream)<<"fitDebug"<< 
	Form("%s=",corr->GetName())<<mass[icorr]<<
	Form("%sp0.=",corr->GetName())<<par0<<
	Form("%sp1=",corr->GetName())<<par1;
    }
    (*pcstream)<<"fit"<< "isOK="<<isOK<<"\n";        
    (*pcstream)<<"fitDebug"<< "isOK="<<isOK<<"\n";        
  }
  delete pcstream;
}

AliExternalTrackParam * AliTPCcalibV0::RefitTrack(AliTPCseed *seed, AliTPCCorrection * corr, Double_t xstart, Double_t xstop, Double_t mass){
  //
  // Refit the track:
  //    seed   - tpc track with cluster
  //    corr   - distrotion/correction class  - apllied to the points
  //    xstart - radius to start propagate/update
  //    xstop  - radius to stop propagate/update
  // 
  const Double_t kResetCov=20.;
  const Double_t kSigma=5.;
  Double_t covar[15];
  for (Int_t i=0;i<15;i++) covar[i]=0;
  covar[0]=kSigma*kSigma;
  covar[2]=kSigma*kSigma;
  covar[5]=kSigma*kSigma/Float_t(150*150);
  covar[9]=kSigma*kSigma/Float_t(150*150);
  covar[14]=1*1;
  // 
  AliExternalTrackParam *refit  = new AliExternalTrackParam(*seed);
  refit->PropagateTo(xstart, AliTracker::GetBz());
  refit->AddCovariance(covar);
  refit->ResetCovariance(kResetCov);
  Double_t xmin = TMath::Min(xstart,xstop);
  Double_t xmax = TMath::Max(xstart,xstop);
  Int_t ncl=0;
  //
  Bool_t isOK=kTRUE;
  for (Int_t index0=0; index0<159; index0++){
    Int_t irow= (xstart<xstop)? index0:159-index0;
    AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);  //cluster in local system
    if (!cluster) continue;
    if (cluster->GetX()<xmin) continue;
    if (cluster->GetX()>xmax) continue;
    Double_t alpha = TMath::Pi()*(cluster->GetDetector()+0.5)/9.;
    if (!refit->Rotate(alpha)) isOK=kFALSE;
    Double_t x     = cluster->GetX();
    Double_t y     = cluster->GetY();
    Double_t z     = cluster->GetZ();
    if (corr){      
      Float_t xyz[3]={cluster->GetX(),cluster->GetY(),cluster->GetZ()};  // original position
      corr->DistortPointLocal(xyz,cluster->GetDetector());
      x=xyz[0];
      y=xyz[1];
      z=xyz[2];
    }
    if (!AliTracker::PropagateTrackToBxByBz(refit, x,mass,1.,kFALSE)) isOK=kFALSE;
    if (!isOK) continue;
    Double_t cov[3]={0.01,0.,0.01};
    Double_t yz[2]={y,z};
    if (!refit->Update(yz,cov)) isOK=kFALSE;    
    ncl++;
  }
  if (!AliTracker::PropagateTrackToBxByBz(refit, xstop, mass,1.,kTRUE)) isOK=kFALSE;
  //  
  if (!isOK) {
    delete refit;
    return 0;
  }  
  return refit;
}





//
// Obsolete part
//


void AliTPCcalibV0::MakeV0s(){
  //
  // 
  //
  const Int_t kMinCluster=40;
  const Float_t kMinR    =0;
  if (! fV0s) fV0s = new TObjArray(10);
  fV0s->Clear();
  //
  // Old V0 finder
  //
  for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
    AliESDv0 * v0 = fESD->GetV0(ivertex);
    if (v0->GetOnFlyStatus()) continue;
    fV0s->AddLast(v0);
  }
  ProcessV0(0);
  fV0s->Clear(0);
  //
  // MI V0 finder
  //
  for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
    AliESDv0 * v0 = fESD->GetV0(ivertex);
    if (!v0->GetOnFlyStatus()) continue;
    fV0s->AddLast(v0);
  }
  ProcessV0(1);
  fV0s->Clear();
  //
  // combinatorial
  //
  Int_t ntracks = fESD->GetNumberOfTracks();
  for (Int_t itrack0=0; itrack0<ntracks; itrack0++){
    AliESDtrack * track0 = fESD->GetTrack(itrack0);
    if (track0->GetSign()>0) continue;
    if ( track0->GetTPCNcls()<kMinCluster) continue;
    if (track0->GetKinkIndex(0)>0) continue;    
    //
    for (Int_t itrack1=0; itrack1<ntracks; itrack1++){
      AliESDtrack * track1 = fESD->GetTrack(itrack1);
      if (track1->GetSign()<0) continue;
      if ( track1->GetTPCNcls()<kMinCluster) continue;
      if (track1->GetKinkIndex(0)>0) continue;
      //
      //      AliExternalTrackParam param0(*track0);
      // AliExternalTrackParam param1(*track1);
      AliV0 vertex;
      vertex.SetParamN(*track0);
      vertex.SetParamP(*track1);
      Float_t xyz[3];
      xyz[0] = fESD->GetPrimaryVertex()->GetXv();
      xyz[1] = fESD->GetPrimaryVertex()->GetYv();
      xyz[2] = fESD->GetPrimaryVertex()->GetZv();
      vertex.Update(xyz);
      if (vertex.GetRr()<kMinR) continue;
      if (vertex.GetDcaV0Daughters()>1.) continue;
      if (vertex.GetDcaV0Daughters()>0.3*vertex.GetRr()) continue;
      // if (vertex.GetPointAngle()<0.9) continue;
      vertex.SetIndex(0,itrack0);
      vertex.SetIndex(1,itrack1);      
      fV0s->AddLast(new AliV0(vertex));
    }
  }
  ProcessV0(2);
  for (Int_t i=0;i<fV0s->GetEntries(); i++) delete fV0s->At(i);
  fV0s->Clear();
}








void AliTPCcalibV0::ProcessV0(Int_t ftype){
  //
  // Obsolete
  //  
  if (! fGammas) fGammas = new TObjArray(10);
  fGammas->Clear();
  Int_t nV0s  = fV0s->GetEntries();
  if (nV0s==0) return;
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  //
  for (Int_t ivertex=0; ivertex<nV0s; ivertex++){
    AliESDv0 * v0 = (AliESDv0*)fV0s->At(ivertex);
    AliESDtrack * trackN = fESD->GetTrack(v0->GetIndex(0)); // negative track
    AliESDtrack * trackP = fESD->GetTrack(v0->GetIndex(1)); // positive track
    
    const AliExternalTrackParam * paramInNeg = trackN->GetInnerParam();
    const AliExternalTrackParam * paramInPos = trackP->GetInnerParam();
  
    if (!paramInPos) continue; // in case the inner paramters do not exist
    if (!paramInNeg) continue;
    // 
    // 
    //
    AliKFParticle *v0K0       = Fit(primVtx,v0,-211,211);
    AliKFParticle *v0Gamma    = Fit(primVtx,v0,11,-11);
    AliKFParticle *v0Lambda42 = Fit(primVtx,v0,-2212,211);
    AliKFParticle *v0Lambda24 = Fit(primVtx,v0,-211,2212);
    //Set production vertex
    v0K0->SetProductionVertex( primVtx );
    v0Gamma->SetProductionVertex( primVtx );
    v0Lambda42->SetProductionVertex( primVtx );
    v0Lambda24->SetProductionVertex( primVtx );
    Double_t massK0, massGamma, massLambda42,massLambda24, massSigma;
    v0K0->GetMass( massK0,massSigma);
    v0Gamma->GetMass( massGamma,massSigma);
    v0Lambda42->GetMass( massLambda42,massSigma);
    v0Lambda24->GetMass( massLambda24,massSigma);
    Float_t chi2K0       = v0K0->GetChi2()/v0K0->GetNDF();
    Float_t chi2Gamma    = v0Gamma->GetChi2()/v0Gamma->GetNDF();
    Float_t chi2Lambda42 = v0Lambda42->GetChi2()/v0Lambda42->GetNDF();
    Float_t chi2Lambda24 = v0Lambda24->GetChi2()/v0Lambda24->GetNDF();
    //
    // Mass Contrained params
    //
    AliKFParticle *v0K0C       = Fit(primVtx,v0,-211,211);
    AliKFParticle *v0GammaC    = Fit(primVtx,v0,11,-11);
    AliKFParticle *v0Lambda42C = Fit(primVtx,v0,-2212,211); //lambdaBar
    AliKFParticle *v0Lambda24C = Fit(primVtx,v0,-211,2212); //lambda
    //   
    v0K0C->SetProductionVertex( primVtx );
    v0GammaC->SetProductionVertex( primVtx );
    v0Lambda42C->SetProductionVertex( primVtx );
    v0Lambda24C->SetProductionVertex( primVtx );

    v0K0C->SetMassConstraint(fPdg->GetParticle(310)->Mass());
    v0GammaC->SetMassConstraint(0);
    v0Lambda42C->SetMassConstraint(fPdg->GetParticle(-3122)->Mass());
    v0Lambda24C->SetMassConstraint(fPdg->GetParticle(3122)->Mass());
    //    
    Double_t timeK0, sigmaTimeK0;  
    Double_t timeLambda42, sigmaTimeLambda42;  
    Double_t timeLambda24, sigmaTimeLambda24;  
    v0K0C->GetLifeTime(timeK0, sigmaTimeK0);
    //v0K0Gamma->GetLifeTime(timeK0, sigmaTimeK0);
    v0Lambda42C->GetLifeTime(timeLambda42, sigmaTimeLambda42);
    v0Lambda24C->GetLifeTime(timeLambda24, sigmaTimeLambda24);
    

    //
    Float_t chi2K0C       = v0K0C->GetChi2()/v0K0C->GetNDF();
    if (chi2K0C<0) chi2K0C=100;
    Float_t chi2GammaC    = v0GammaC->GetChi2()/v0GammaC->GetNDF();
    if (chi2GammaC<0) chi2GammaC=100;
    Float_t chi2Lambda42C = v0Lambda42C->GetChi2()/v0Lambda42C->GetNDF();
    if (chi2Lambda42C<0) chi2Lambda42C=100;
    Float_t chi2Lambda24C = v0Lambda24C->GetChi2()/v0Lambda24C->GetNDF();
    if (chi2Lambda24C<0) chi2Lambda24C=100;
    //
    Float_t  minChi2C=99;
    Int_t   type   =-1;
    if (chi2K0C<minChi2C) { minChi2C= chi2K0C; type=0;}
    if (chi2GammaC<minChi2C) { minChi2C= chi2GammaC; type=1;}
    if (chi2Lambda42C<minChi2C) { minChi2C= chi2Lambda42C; type=2;}
    if (chi2Lambda24C<minChi2C) { minChi2C= chi2Lambda24C; type=3;}
    Float_t  minChi2=99;
    Int_t   type0   =-1;
    if (chi2K0<minChi2) { minChi2= chi2K0; type0=0;}
    if (chi2Gamma<minChi2) { minChi2= chi2Gamma; type0=1;}
    if (chi2Lambda42<minChi2) { minChi2= chi2Lambda42; type0=2;}
    if (chi2Lambda24<minChi2) { minChi2= chi2Lambda24; type0=3;}
    
     // 0 is  negative particle; 1 is positive particle
    Float_t betaGamma0 = 0;
    Float_t betaGamma1 = 0;
    
    switch (type) {
     case 0:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-211)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(211)->Mass();
      break;
     case 1:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(11)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(-11)->Mass();
      break;
     case 2:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-2212)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(211)->Mass();
      break;
     case 3:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-211)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(2212)->Mass();
      break;
    }
 
    // cuts and histogram filling
    Int_t numCand = 0; // number of particle types which have a chi2 < 10*minChi2
        
    if (minChi2C < 2 && ftype == 1) {
     //
     if (chi2K0C < 10*minChi2C) numCand++;
     if (chi2GammaC < 10*minChi2C) numCand++;
     if (chi2Lambda42C < 10*minChi2C) numCand++;
     if (chi2Lambda24C < 10*minChi2C) numCand++;
     //
    }
    
    //
    //
    // write output tree
    if (minChi2>50) continue;
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      (*cstream)<<"V0"<<
	"ftype="<<ftype<<
	"v0.="<<v0<<
	"trackN.="<<trackN<<
	"trackP.="<<trackP<<
	//
	"betaGamma0="<<betaGamma0<<
	"betaGamma1="<<betaGamma1<<
	//
	"type="<<type<<
	"chi2C="<<minChi2C<<
	"v0K0.="<<v0K0<<
	"v0Gamma.="<<v0Gamma<<
	"v0Lambda42.="<<v0Lambda42<<
	"v0Lambda24.="<<v0Lambda24<<
	//
	"chi20K0.="<<chi2K0<<
	"chi2Gamma.="<<chi2Gamma<<
	"chi2Lambda42.="<<chi2Lambda42<<
	"chi2Lambda24.="<<chi2Lambda24<<
	//
	"chi20K0c.="<<chi2K0C<<
	"chi2Gammac.="<<chi2GammaC<<
	"chi2Lambda42c.="<<chi2Lambda42C<<
	"chi2Lambda24c.="<<chi2Lambda24C<<
	//
	"v0K0C.="<<v0K0C<<
	"v0GammaC.="<<v0GammaC<<
	"v0Lambda42C.="<<v0Lambda42C<<
	"v0Lambda24C.="<<v0Lambda24C<<
	//
	"massK0="<<massK0<<
	"massGamma="<<massGamma<<
	"massLambda42="<<massLambda42<<
	"massLambda24="<<massLambda24<<
	//
	"timeK0="<<timeK0<<
	"timeLambda42="<<timeLambda42<<
	"timeLambda24="<<timeLambda24<<
	"\n";
    }
    if (type==1) fGammas->AddLast(v0); 
    //
    //
    //
    delete v0K0;
    delete v0Gamma;
    delete v0Lambda42;
    delete v0Lambda24;    
    delete v0K0C;
    delete v0GammaC;
    delete v0Lambda42C;
    delete v0Lambda24C; 
  }
  ProcessPI0(); 
}



void AliTPCcalibV0::ProcessPI0(){
  //
  //
  //
  Int_t nentries = fGammas->GetEntries();
  if (nentries<2) return;
  // 
  Double_t m0[3], m1[3];
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  for (Int_t i0=0; i0<nentries; i0++){
    AliESDv0 *v00 = (AliESDv0*)fGammas->At(i0); 
    v00->GetPxPyPz (m0[0], m0[1], m0[2]);
    AliKFParticle *p00 = Fit(primVtx, v00, 11,-11);
    p00->SetProductionVertex( primVtx );
    p00->SetMassConstraint(0);
    //
    for (Int_t i1=i0; i1<nentries; i1++){
      AliESDv0 *v01 = (AliESDv0*)fGammas->At(i1);
      v01->GetPxPyPz (m1[0], m1[1], m1[2]);
      AliKFParticle *p01 = Fit(primVtx, v01, 11,-11);
      p01->SetProductionVertex( primVtx );
      p01->SetMassConstraint(0);
      if (v00->GetIndex(0) != v01->GetIndex(0) && 
	  v00->GetIndex(1) != v01->GetIndex(1)){
	AliKFParticle pi0( *p00,*p01); 
	pi0.SetProductionVertex(primVtx);
	Double_t n1 = TMath::Sqrt (m0[0]*m0[0] + m0[1]*m0[1] + m0[2]*m0[2]);
        Double_t n2 = TMath::Sqrt (m1[0]*m1[0] + m1[1]*m1[1] + m1[2]*m1[2]);
        Double_t mass = TMath::Sqrt(2.*(n1*n2 - (m0[0]*m1[0] + m0[1]*m1[1] + m0[2]*m1[2])));
	TTreeSRedirector *cstream = GetDebugStreamer();
	if (cstream){
	  (*cstream)<<"PI0"<<
	    "v00.="<<v00<<
	    "v01.="<<v01<<
	    "mass="<<mass<<
	    "p00.="<<p00<<
	    "p01.="<<p01<<
	    "pi0.="<<&pi0<<
	    "\n";	
	}
      }
      delete p01;
    }
    delete p00;
  }
}





AliKFParticle * AliTPCcalibV0::Fit(AliKFVertex & /*primVtx*/, AliESDv0 *v0, Int_t PDG1, Int_t PDG2){
  //
  // Make KF Particle
  //
  AliKFParticle p1( *(v0->GetParamN()), PDG1 );
  AliKFParticle p2( *(v0->GetParamP()), PDG2 );
  AliKFParticle *V0 = new AliKFParticle;
  Double_t x, y, z;
  v0->GetXYZ(x,y,z );
  V0->SetVtxGuess(x,y,z);
  *(V0)+=p1;
  *(V0)+=p2;
  return V0;  
}




void AliTPCcalibV0::BinLogX(TH2F *h) {
  //
  //
  //
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Double_t from = axis->GetXmin();
   Double_t to = axis->GetXmax();
   Double_t *new_bins = new Double_t[bins + 1];
   
   new_bins[0] = from;
   Double_t factor = pow(to/from, 1./bins);
  
   for (int i = 1; i <= bins; i++) {
     new_bins[i] = factor * new_bins[i-1];
   }
   axis->Set(bins, new_bins);
   delete [] new_bins;
   
}




