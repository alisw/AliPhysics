
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
//#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDVertex.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVfriendTrack.h"
#include "AliVfriendEvent.h"

#include "AliMathBase.h" 
#include "AliTPCseed.h"
#include "AliTPCreco.h"
#include "AliTPCclusterMI.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliTrackPointArray.h"
#include "AliTPCcalibV0.h"
#include "AliV0.h"
#include "TRandom.h"
#include "TTreeStream.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCorrection.h"
#include "AliGRPObject.h"
#include "AliTPCTransform.h"
#include "AliAnalysisManager.h"



ClassImp(AliTPCcalibV0)


AliTPCcalibV0::AliTPCcalibV0() : 
   AliTPCcalibBase(),
   fV0Tree(0),
   fHPTTree(0),
   fStack(0),
   fEvent(0),
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
   fEvent(0),
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





void  AliTPCcalibV0::ProcessESD(AliVEvent *event){
  //
  //
  //
  fEvent = event;
  AliKFParticle::SetField(event->GetMagneticField());
  if (TMath::Abs(AliTracker::GetBz())<1) return;  
  DumpToTree(event);
  DumpToTreeHPT(event);
}

void  AliTPCcalibV0::DumpToTreeHPT(AliVEvent *event){
  //
  // Dump V0s fith full firend information to the 
  // 
  if (TMath::Abs(AliTracker::GetBz())<1) return;
  const Int_t kMinCluster=110;
  const Float_t kMinPt   =4.;
  AliVfriendEvent *friendEvent=event->FindFriend();
//   if (!esdFriend) {
//     Printf("ERROR: esdFriend not available");
//     return;
//   }
  //
  Int_t ntracks=event->GetNumberOfTracks();
  for (Int_t i=0;i<ntracks;++i) {
    Bool_t isOK=kFALSE;
    AliVTrack *track = event->GetVTrack(i);
    if(!track) continue;
    if (track->GetTPCncls()<kMinCluster) continue;
    if (TMath::Abs(AliTracker::GetBz())>1){ // cut on momenta if measured
      if (track->Pt()>kMinPt) isOK=kTRUE;
    }
    if (TMath::Abs(AliTracker::GetBz())<1){  // require primary track for the B field OFF data
      Bool_t isAccepted=kTRUE;
      if (!track->IsOn(AliVTrack::kITSrefit)) isAccepted=kFALSE;
      if (!track->IsOn(AliVTrack::kTPCrefit)) isAccepted=kFALSE;
      if (!track->IsOn(AliVTrack::kTOFout))   isAccepted=kFALSE;
      Float_t dvertex[2],cvertex[3]; 
      track->GetImpactParametersTPC(dvertex,cvertex);
      if (TMath::Abs(dvertex[0]/TMath::Sqrt(cvertex[0]+0.01))>20) isAccepted=kFALSE;
      if (TMath::Abs(dvertex[1]/TMath::Sqrt(TMath::Abs(cvertex[2]+0.01)))>20) isAccepted=kFALSE;
      track->GetImpactParameters(dvertex,cvertex);
      if (TMath::Abs(dvertex[0]/TMath::Sqrt(cvertex[0]+0.01))>10) isAccepted=kFALSE;
      if (TMath::Abs(dvertex[1]/TMath::Sqrt(TMath::Abs(cvertex[2]+0.01)))>10) isAccepted=kFALSE;
      if (!isAccepted) isOK=kFALSE;
    } 
    if ( track->GetTPCsignal()>100 && track->GetInnerParam()->Pt()>1 ){
      if (track->IsOn(AliVTrack::kITSin)||track->IsOn(AliVTrack::kTRDout)||track->IsOn(AliVTrack::kTOFin))
	isOK=kTRUE;
      if (isOK){
	TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
	Int_t eventNumber = event->GetEventNumberInFile();
	Bool_t hasFriend=(friendEvent) ? (friendEvent->GetTrack(i)!=0):0;
	Bool_t hasITS=(track->GetNcls(0)>2);
	printf("DUMPIONTrack:%s|%f|%d|%d|%d\n",filename.Data(),track->GetInnerParam()->Pt()*track->GetTPCsignal()/50., eventNumber,hasFriend,hasITS);
      }
    }
    if (!isOK) continue;
    TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
    Int_t eventNumber = event->GetEventNumberInFile();
    Bool_t hasFriend=(friendEvent) ? (friendEvent->GetTrack(i)!=0):0;
    Bool_t hasITS=(track->GetNcls(0)>2);    
    printf("DUMPHPTTrack:%s|%f|%d|%d|%d\n",filename.Data(),track->Pt(), eventNumber,hasFriend,hasITS);
    //
    if (!friendEvent) continue;
    const AliVfriendTrack *friendTrack = friendEvent->GetTrack(i);
    if (!friendTrack) continue;

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

      //**********************TEMPORARY!!*******************************************
      // more investigation is needed with Tree ///!!!
      //all dummy stuff here is just for code to compile and work with ESD

      AliESDfriendTrack *dummyfriendTrack=(AliESDfriendTrack*)friendTrack;
      AliESDtrack *dummytrack=(AliESDtrack*)track;


    if (fHPTTree->GetEntries()==0){
      //
      fHPTTree->SetDirectory(0);
      fHPTTree->Branch("t.",&dummytrack);
      fHPTTree->Branch("ft.",&dummyfriendTrack);
      fHPTTree->Branch("s.",&seed);
    }else{
      fHPTTree->SetBranchAddress("t.",&dummytrack);
      fHPTTree->SetBranchAddress("ft.",&dummyfriendTrack);
      fHPTTree->SetBranchAddress("s.",&seed);
    }
    fHPTTree->Fill();
    //
  }
}



void  AliTPCcalibV0::DumpToTree(AliVEvent *event){
  //
  // Dump V0s fith full firend information to the 
  // 
  Int_t nV0s  = fEvent->GetNumberOfV0s();
  const Int_t kMinCluster=110;
  const Double_t kDownscale=0.01;
  const Float_t kMinPt   =1.0;
  const Float_t kMinMinPt   =0.7;
  AliVfriendEvent *friendEvent=event->FindFriend();
  //
  
  for (Int_t ivertex=0; ivertex<nV0s; ivertex++){
    Bool_t isOK=kFALSE;
    AliESDv0 dummyv0;
    event->GetV0(dummyv0,ivertex);
    AliESDv0 *v0=&dummyv0;

    AliVTrack * track0 = fEvent->GetVTrack(v0->GetIndex(0)); // negative track
    AliVTrack * track1 = fEvent->GetVTrack(v0->GetIndex(1)); // positive track
    if(!track0) continue;
    if (track0->GetTPCNcls()<kMinCluster) continue;
    if (track0->GetKinkIndex(0)>0) continue;
    if(!track1) continue;
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
    TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
    Int_t eventNumber = event->GetEventNumberInFile();
    Bool_t hasITS=(track0->GetNcls(0)+ track1->GetNcls(0)>4);
    printf("DUMPHPTV0:%s|%f|%d|%d|%d\n",filename.Data(), (TMath::Min(track0->Pt(),track1->Pt())), eventNumber,(friendEvent!=0), hasITS);
    //
    if (!friendEvent) continue;
    //
    
    //
    const AliVfriendTrack *ftrack0 = friendEvent->GetTrack(v0->GetIndex(0));
    if (!ftrack0) continue;
    const AliVfriendTrack *ftrack1 = friendEvent->GetTrack(v0->GetIndex(1));
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

    //**********************TEMPORARY!!*******************************************
    // more investigation is needed with Tree ///!!!
    //all dummy stuff here is just for code to compile and work with ESD

    AliESDfriendTrack *dummyftrack0=(AliESDfriendTrack*)ftrack0;
    AliESDfriendTrack *dummyftrack1=(AliESDfriendTrack*)ftrack1;
    AliESDtrack *dummytrack0=(AliESDtrack*)track0;
    AliESDtrack *dummytrack1=(AliESDtrack*)track1;

    if (fV0Tree->GetEntries()==0){
      //
      fV0Tree->SetDirectory(0);
      fV0Tree->Branch("v0.",&v0);
      fV0Tree->Branch("t0.",&dummytrack0);
      fV0Tree->Branch("t1.",&dummytrack1);
      fV0Tree->Branch("ft0.",&dummyftrack0);
      fV0Tree->Branch("ft1.",&dummyftrack1);
      fV0Tree->Branch("s0.",&seed0);
      fV0Tree->Branch("s1.",&seed1);
    }else{
      fV0Tree->SetBranchAddress("v0.",&v0);
      fV0Tree->SetBranchAddress("t0.",&dummytrack0);
      fV0Tree->SetBranchAddress("t1.",&dummytrack1);
      fV0Tree->SetBranchAddress("ft0.",&dummyftrack0);
      fV0Tree->SetBranchAddress("ft1.",&dummyftrack1);
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
  return ;
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
  return ;
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
    for (Int_t irow=0; irow<kMaxRow; irow++){
      AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
      if (cluster &&cluster->GetX()>10){
        Double_t x0[3]={ static_cast<Double_t>(cluster->GetRow()),cluster->GetPad(),cluster->GetTimeBin()};
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
      for (Int_t irow=0; irow<kMaxRow; irow++){
	AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
	if (cluster &&cluster->GetX()>10){
	  Double_t x0[3]={ static_cast<Double_t>(cluster->GetRow()),cluster->GetPad(),cluster->GetTimeBin()};
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
  for (Int_t index0=0; index0<kMaxRow; index0++){
    Int_t irow= (xstart<xstop)? index0:kMaxRow-1-index0;
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



void AliTPCcalibV0::FilterV0s(AliVEvent *event){
  //
  // 
  TDatabasePDG pdg;  
  const Double_t kChi2Cut=20;
  const Double_t kMinR=2;
  const Double_t ptCut=0.2;
  const Int_t kMinNcl=110;
  //
  Int_t nv0 = event->GetNumberOfV0s();
  //AliESDVertex *vertex= (AliESDVertex *)event->GetPrimaryVertex();
  //AliKFVertex kfvertex=*vertex;

  //AliESDVertex vtx;
  //event->GetPrimaryVertex(vtx);
  //AliESDVertex *vertex=&vtx;
  //AliKFVertex *kfvertex=(AliKFVertex*)vertex;

  AliESDVertex vtx;
  event->GetPrimaryVertex(vtx);
  AliESDVertex *vertex=&vtx;
  AliKFVertex kfvertex=*vertex;

  //
  for (Int_t iv0=0;iv0<nv0;iv0++){
    AliESDv0 dummyv0;
    event->GetV0(dummyv0,iv0);
    AliESDv0 *v0=&dummyv0;

    if (!v0) continue;
    if (v0->GetPindex()<0) continue;
    if (v0->GetNindex()<0) continue;
    if (TMath::Max(v0->GetPindex(), v0->GetNindex())>event->GetNumberOfTracks()) continue;
    //
    //   
    AliExternalTrackParam pp=(v0->GetParamP()->GetSign()>0) ? (*(v0->GetParamP())):(*(v0->GetParamN()));
    AliExternalTrackParam pn=(v0->GetParamP()->GetSign()>0) ? (*(v0->GetParamN())):(*(v0->GetParamP()));
    AliKFParticle kfp1( pp, 211 );
    AliKFParticle kfp2( pn, -211 );
    AliKFParticle *v0KFK0 = new AliKFParticle(kfp1,kfp2);
    AliKFParticle *v0KFK0CV = new AliKFParticle(*v0KFK0);
    v0KFK0CV->SetProductionVertex(kfvertex);
    v0KFK0CV->TransportToProductionVertex();
    AliKFParticle *v0KFK0CVM = new AliKFParticle(*v0KFK0CV);
    v0KFK0CVM->SetMassConstraint(pdg.GetParticle("K_S0")->Mass());
    Double_t chi2K0 = v0KFK0CV->GetChi2();
    //    Double_t chi2K0M= v0KFK0CVM->GetChi2();    
    //Double_t massK0=v0KFK0CV->GetMass();
    if (chi2K0>kChi2Cut) continue;
    if (v0->GetRr()>kMinR) continue;
    //
    Double_t maxPt = TMath::Max(v0->GetParamP()->Pt(), v0->GetParamN()->Pt());
    Double_t effMass22=v0->GetEffMass(2,2);
    Double_t effMass42=v0->GetEffMass(4,2);
    Double_t effMass24=v0->GetEffMass(2,4);
    Bool_t isV0= kFALSE;
    isV0|=TMath::Abs(effMass22-pdg.GetParticle("K_S0")->Mass())<0.1;
    isV0|=TMath::Abs(effMass42-pdg.GetParticle("Lambda0")->Mass())<0.1;
    isV0|=TMath::Abs(effMass24-pdg.GetParticle("Lambda0")->Mass())<0.1;

    Double_t sign= v0->GetParamP()->GetSign()* v0->GetParamN()->GetSign();
    if (sign<0&&v0->GetOnFlyStatus()>0.5&&maxPt>ptCut&&isV0){
      AliVTrack * trackP = event->GetVTrack(v0->GetPindex());
      AliVTrack * trackN = event->GetVTrack(v0->GetNindex());
      if (!trackN) continue;
      if (!trackP) continue;
      Int_t nclP= (Int_t)trackP->GetTPCClusterInfo(2,1);
      Int_t nclN= (Int_t)trackN->GetTPCClusterInfo(2,1);
      if (TMath::Min(nclP,nclN)<kMinNcl) continue;

      AliExternalTrackParam trkprmP;
      trackP->GetTrackParam(trkprmP);
      AliExternalTrackParam trkprmN;
      trackN->GetTrackParam(trkprmN);

      Double_t eta = TMath::Max(TMath::Abs(trkprmP.Eta()), TMath::Abs(trkprmN.Eta()));
      Double_t ncls = TMath::Min(TMath::Abs(trackP->GetNcls(0)), TMath::Abs(trackN->GetNcls(0)));
      if (eta<0.8&&ncls>2){
	//	printf("%d\t%f\t%f\t%d\t%d\t%f\t%f\n",i, v0->Pt(), maxPt, v0->GetNindex(),v0->GetPindex(),v0->GetRr(),effMass22);	
	(*fDebugStreamer)<<"v0tree"<<
	  "v0.="<<v0<<
	  "tp.="<<trackP<<
	  "tm.="<<trackN<<
	  //
      //"v.="<<vertex<<
	  "ncls="<<ncls<<
	  "maxPt="<<maxPt<<
	  "\n";        
      }
    }
  }
}
