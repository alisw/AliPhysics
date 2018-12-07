/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
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
//
// The task:
// stores TRD PID quantities in a Tree
// output can then be used for e.g. parameter creation
//
// Author:
// Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//


#include "AliTRDPIDTree.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
//#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliTRDgeometry.h"
#include "TTree.h"

class TCanvas;
class TAxis;
class TFile;
class TStyle;
class TString;
class TH1F;
class TH2D;
class THnSparse;
class TLegend;
class TVirtualFitter;
class AliESDtrackCuts;
class AliStack;
class AliMCParticle;


using namespace std;

ClassImp(AliTRDPIDTree)

Double_t AliTRDPIDTree::fgMinLayer = 1.;

//________________________________________________________________________
AliTRDPIDTree::AliTRDPIDTree(const char *name)
    : AliAnalysisTaskSE(name), fV0tags(0x0), fV0cuts(0x0), fV0electrons(0x0), fV0pions(0x0), fV0protons(0x0),
    fESDEvent(0), fMCEvent(0), fMCStack(0), fTreeTRDPID(0), fPIDResponse(0), fOutputContainer(0), fESDtrackCuts(0),
    fESDtrackCutsV0(0), fListQATRD(0x0), fListQATRDV0(0x0),
    fNumTagsStored(0), fCollisionSystem(3),
    fpdg(0), frun(0), frunnumber(0), fcentrality(0), fTRDNcls(0), fTRDntracklets(0), fTRDntrackletsPID(0),
    fTRDtheta(0), fTRDTPCtgl(0), fTRDsignal(0), fTRDnclsdEdx(0), fTRDnch(0), fPDG(0), fTrackCharge(0), fPDGTRUE(0), fChi2(0),
    fhtrackCuts(0), fhEventCount(0), fhArmenteros(0), fUseExtraPileupCut(0), fHistV0MvsTPCoutBeforePileUpCuts(0), fHistV0MvsTPCoutAfterPileUpCuts(0)
{

  //
  // Constructor
  //
     
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());

}


//_________________________________________________
AliTRDPIDTree::~AliTRDPIDTree()
{

  //
  // Destructor
  //

    delete fV0cuts;
    delete fV0electrons;
    delete fV0pions;
    delete fV0protons;
    delete fV0tags;
    fV0tags = 0;
    fNumTagsStored = 0;
}


//________________________________________________________________________
void AliTRDPIDTree::UserCreateOutputObjects()
{
    //
    // Definition of user output ntuple and histogram file
    //


    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandler)
    printf("Inputhandler not available \n");
  else
    fPIDResponse = inputHandler->GetPIDResponse();

   // V0 Kine cuts 
  fV0cuts = new AliESDv0KineCuts();

  // V0 PID Obj arrays
  fV0electrons = new TObjArray;
  fV0pions     = new TObjArray;
  fV0protons   = new TObjArray;

 

  OpenFile(1);
    fTreeTRDPID = new TTree("TreeTRDPID", "TRD PID");
    fTreeTRDPID->Branch("run", &frunnumber);
    fTreeTRDPID->Branch("centrality", &fcentrality);
    fTreeTRDPID->Branch("TRDslices[48]",fTRDslices);
    fTreeTRDPID->Branch("TRDMomentum[6]",fTRDMomentum);
    fTreeTRDPID->Branch("TRDNcls",&fTRDNcls);
    fTreeTRDPID->Branch("TRDntracklets",&fTRDntracklets);
    fTreeTRDPID->Branch("TRDntrackletsPID",&fTRDntrackletsPID); 
    fTreeTRDPID->Branch("TRDphi[6]",fTRDphi); // values in different layers slightly different, clear correlation -> keep all layers for the moment
    fTreeTRDPID->Branch("TRDtheta",&fTRDtheta); // value in different layers identical -> keep layer 0 only
    fTreeTRDPID->Branch("TRDY[6]",fTRDY); // raus
    fTreeTRDPID->Branch("TRDTPCtgl",&fTRDTPCtgl);
    fTreeTRDPID->Branch("TRDsignal",&fTRDsignal); 
    fTreeTRDPID->Branch("TRDnclsdEdx",&fTRDnclsdEdx); 
    fTreeTRDPID->Branch("TRDnch",&fTRDnch);  
    fTreeTRDPID->Branch("NSigmaTPC[3]",fNSigmaTPC);
    fTreeTRDPID->Branch("NSigmaTOF[3]",fNSigmaTOF);
    fTreeTRDPID->Branch("PDG",&fPDG);
    fTreeTRDPID->Branch("TrackCharge",&fTrackCharge);
    fTreeTRDPID->Branch("PDGTRUE",&fPDGTRUE);
    fTreeTRDPID->Branch("DCA[2]",fDCA);
    fTreeTRDPID->Branch("Chi2",&fChi2);
    fTreeTRDPID->Branch("TRDsigma[5]",fsigmaTRD);
    fTreeTRDPID->Branch("TRDdelta[5]",fdeltaTRD);
    fTreeTRDPID->Branch("TRDratio[5]",fratioTRD);

    PostData(1,fTreeTRDPID);


  OpenFile(2);
  fListQATRD=new TList;
  fListQATRD->SetOwner();
  fListQATRDV0=new TList;
  fListQATRDV0->SetOwner();
  fListQATRDV0->SetName("V0decay");
  fListQATRD->Add(fListQATRDV0);

  SetupV0qa();

  switch(fUseExtraPileupCut){
    case kLHC15o :
      printf("Using extra pile-up cut on TPCout <-> VZEROTotalMult correlation with parameters for LHC15o.\n");
      break;

    case kLHC18q :
      printf("Using extra pile-up cut on TPCout <-> VZEROTotalMult correlation with parameters for LHC18q.\n");
      break;
  }

  fHistV0MvsTPCoutBeforePileUpCuts = new TH2F("fHistV0MvsTPCoutBeforePileUpCuts","V0M amplitude vs TPCout tracks; TPCout tracks; V0M amplitude;",1000,0,20000,1000,0,40000);
  fListQATRD->Add(fHistV0MvsTPCoutBeforePileUpCuts);
  fHistV0MvsTPCoutAfterPileUpCuts = new TH2F("fHistV0MvsTPCoutAfterPileUpCuts","V0M amplitude vs TPCout tracks; TPCout tracks; V0M amplitude;",1000,0,20000,1000,0,40000);
  fListQATRD->Add(fHistV0MvsTPCoutAfterPileUpCuts);  
  fhtrackCuts  = new TH1F("fhtrackCuts","TrackEventCuts QA",10,-0.5,9.5);
  fListQATRD->Add(fhtrackCuts);
  fhEventCount  = new TH1F("fhEventCount","Event Count",5,-0.5,4.5);
  fListQATRD->Add(fhEventCount);

  fHistV0MvsTPCoutBeforePileUpCuts->SetMarkerStyle(20);
  fHistV0MvsTPCoutBeforePileUpCuts->SetMarkerSize(1);
  fHistV0MvsTPCoutAfterPileUpCuts->SetMarkerStyle(20);
  fHistV0MvsTPCoutAfterPileUpCuts->SetMarkerSize(1);

  
  PostData(2,fListQATRD);


}


//_____________________________________________________________________________
void AliTRDPIDTree::UserExec(Option_t *)
{
    //
    //calls the Process function
    //

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      printf("ERROR: Could not get ESDInputHandler \n");
    }
    else fESDEvent = (AliESDEvent*)esdH->GetEvent();
    
    // If MC, forward MCevent
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    fMCStack=NULL;
    if(fMCEvent){
	fMCStack = fMCEvent->Stack();
    }

    FillV0PIDlist();
    //
    Process(fESDEvent, fMCEvent);
    //

    PostData(1,fTreeTRDPID);
    PostData(2,fListQATRD);

    // Clear the V0 PID arrays
    ClearV0PIDlist();

}



//________________________________________________________________________
void AliTRDPIDTree::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
{
    //
    //called for each event
    //

  if (!esdEvent) {
    Printf("ERROR: esdEvent not available"); 
    return;
  }

  if (!mcEvent) {
//    Printf("ERROR: mcEvent not available");
//    return;
  }

  if (!fPIDResponse ) {
    Printf("ERROR: No PIDresponse and/or v0KineCuts!");
    return;
  }

 // Float_t centralityFper=99;

 // AliCentrality *esdCentrality = esdEvent->GetCentrality();
 // centralityFper = esdCentrality->GetCentralityPercentile("V0M");

  Float_t centralityFper = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) esdEvent->FindListObject("MultSelection");
  if( !MultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
   AliWarning("AliMultSelection object not found!");
  } else{
      centralityFper = MultSelection->GetMultiplicityPercentile("V0M");
  }



  const AliESDVertex* fESDEventvertex = esdEvent->GetPrimaryVertexTracks(); 
  if (!fESDEventvertex)
    return;
  
  Int_t ncontr = fESDEventvertex->GetNContributors();




  if (ncontr <= 0) return;
  frunnumber = fESDEvent->GetRunNumber();

  if(fESDEvent) fhEventCount->Fill(1,1);

  //Extra Pile-up Cut
  if(fUseExtraPileupCut!=0){
    Int_t ntrkTPCout = 0;
    for (int it = 0; it < esdEvent->GetNumberOfTracks(); it++) {
      AliESDtrack* ESDTrk = (AliESDtrack*)esdEvent->GetTrack(it);
      if ((ESDTrk->GetStatus() & AliESDtrack::kTPCout) && ESDTrk->GetID() > 0)
        ntrkTPCout++;
     }

    Double_t multVZERO =0;
    AliVVZERO *vzero = (AliVVZERO*)esdEvent->GetVZEROData();
    if(vzero) {
      for(int ich=0; ich < 64; ich++)
      multVZERO += vzero->GetMultiplicity(ich);
    }


    fHistV0MvsTPCoutBeforePileUpCuts->Fill(ntrkTPCout, multVZERO);

    switch(fUseExtraPileupCut){
      case kLHC15o :
        if (multVZERO < (-2200 + 2.5*ntrkTPCout + 1.2e-5*ntrkTPCout*ntrkTPCout))  {
          return;
        }
        break;

      case kLHC18q :
        if (multVZERO < (-2800 + 3.165*ntrkTPCout + 2.5e-5*ntrkTPCout*ntrkTPCout))  {
          return;
        }
        break;
    }

    fHistV0MvsTPCoutAfterPileUpCuts->Fill(ntrkTPCout, multVZERO);
  }
  
  // - Begin: track loop for electrons from V0 -
    for(Int_t itrack = 0; itrack < fV0electrons->GetEntries(); itrack++){
    //  AliVTrack *track=(AliVTrack*)fV0electrons->At(itrack);
      AliESDtrack *track=(AliESDtrack*)fV0electrons->At(itrack);

      fpdg=11;
      FillTree(track, fpdg, frun, centralityFper);

    } // - End: track loop for electrons from V0 -

  // - Begin: track loop for pions from V0 -
    for(Int_t itrack = 0; itrack < fV0pions->GetEntries(); itrack++){
	//       AliVTrack *track=(AliVTrack*)fV0pions->At(itrack);
	AliESDtrack *track=(AliESDtrack*)fV0pions->At(itrack);
	fpdg=211;
	FillTree(track, fpdg, frun, centralityFper);
    }

    // - Begin: track loop for protons from V0 -
    for(Int_t itrack = 0; itrack < fV0protons->GetEntries(); itrack++){
	//       AliVTrack *track=(AliVTrack*)fV0protons->At(itrack);
	AliESDtrack *track=(AliESDtrack*)fV0protons->At(itrack);
	fpdg=2212;
	FillTree(track, fpdg, frun, centralityFper);
    }

   

    PostData(1,fTreeTRDPID);
    PostData(2,fListQATRD);}
  

//________________________________________________________________________
void AliTRDPIDTree::FillTree(AliESDtrack *track, Int_t pdgfromv0, Int_t runnumber, Int_t centralityvalue)
{
    //
    // Fill tree with track and event variables for e, pi, proton from V0 candidates
    //

    // global event and track cuts
    if(!PassTrackCuts(track)) return;

    // track cuts:
    Int_t ntrackletstracking=track->GetTRDntracklets();
    if(ntrackletstracking<fgMinLayer) return;
    Int_t trdnch  = track->GetTRDNchamberdEdx(); 
//    if(trdnch<1) return;



    // check if TOF can provide PID (if not true return)
    UInt_t status = track->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout;
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    if(!(hasTOFout && hasTOFtime)) return;

    
    Double_t nSigmaTPC = -999;
    Double_t nSigmaTOF = -999;

    if(TMath::Abs(pdgfromv0)==11){
	nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
	nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
    }
    if(TMath::Abs(pdgfromv0)==211){
	nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
	nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
    }
    if(TMath::Abs(pdgfromv0)==2212){
	nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
	nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
    }

    if(nSigmaTPC>3||nSigmaTPC<-3){
	return;
    }

    if(nSigmaTOF>3||nSigmaTOF<-3){
	return;
    }

    frun=runnumber;

    Int_t charge = track->Charge();
    fPDG=pdgfromv0;
    fTrackCharge=charge;
    fcentrality = centralityvalue;
    fTRDntracklets = ntrackletstracking;
    fTRDntrackletsPID = track->GetTRDntrackletsPID();

    fTRDNcls=track->GetTRDncls();
    fTRDTPCtgl=track->GetTPCTgl();
    fChi2=track->GetTRDchi2();

    fTRDsignal=track->GetTRDsignal(); // truncated mean signal
    fTRDnclsdEdx=track->GetTRDNclusterdEdx(); // number of clusters dedx
    fTRDnch=trdnch;  // number of chambers dedx

    track->GetImpactParameters(fDCA[0],fDCA[1]);

    fNSigmaTPC[0]=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
    fNSigmaTPC[1]=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
    fNSigmaTPC[2]=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
    fNSigmaTOF[0]=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
    fNSigmaTOF[1]=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
    fNSigmaTOF[2]=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);

    Int_t ntracklets=0;
    for(Int_t iPl=0;iPl<AliVTrack::kTRDnPlanes;iPl++){

	fTRDphi[iPl]=GetPhi(track,iPl,fTRDY[iPl],fTRDthetalayer[iPl]);
	//  printf("testing %i %f %f %f \n",iPl,fTRDtheta,(-TMath::Log(TMath::Tan(0.5 * fTRDtheta))),fTRDeta[iPl]);
	Float_t dEdx=0;
	fTRDMomentum[iPl]= track->GetTRDmomentum(iPl);
	for(int isl= 0; isl<= 7;isl++){
	    fTRDslices[iPl*8+isl]=track->GetTRDslice(iPl,isl);
	}
    }
    fTRDtheta=fTRDthetalayer[0];
    

    fPDGTRUE=0;
    if(fMCStack){
	Int_t label=track->GetLabel();
	if(label>=0){
	    TParticle *part=fMCStack->Particle(label);
	    if(part){
		fPDGTRUE=part->GetPdgCode();
		if(TMath::Abs(fPDGTRUE)==11)fPDGTRUE=-fPDGTRUE; // switch signs (particle code says 11 for e- and -11 for e+)
	    }
	}
    }


    AliPID::EParticleType types[]={AliPID::kElectron, AliPID::kMuon, AliPID::kPion, AliPID::kKaon, AliPID::kProton};
    for(Int_t itrdpid=0; itrdpid<5; itrdpid++){
	fsigmaTRD[itrdpid]= fPIDResponse->NumberOfSigmas(AliPIDResponse::kTRD, track, types[itrdpid]); //  (signal-expsig)/(sigma), where sigma = mean*resolution
	fdeltaTRD[itrdpid]= fPIDResponse->GetSignalDelta(AliPIDResponse::kTRD, track, types[itrdpid]);   // signal/(expsig)
	fratioTRD[itrdpid]= fPIDResponse->GetSignalDelta(AliPIDResponse::kTRD, track, types[itrdpid], kTRUE);  // difference signal to theor. value divided by resolution
    }


    fTreeTRDPID->Fill();

}

//________________________________________________________________________
Int_t  AliTRDPIDTree::CompareFloat(Float_t f1, Float_t f2) const
{
    //
    //compares if the Float_t f1 is equal with f2 and returns 1 if true and 0 if false
    //
    Float_t precision = 0.00001;
    if (((f1 - precision) < f2) &&
	((f1 + precision) > f2))
    {
	return 1;
    }
    else
    {
	return 0;
    }
}


//________________________________________________________________________
void AliTRDPIDTree::Terminate(const Option_t *)
{
    //
    // Terminate function
    //


}


//______________________________________________________________________________
void AliTRDPIDTree::FillV0PIDlist(){

  //
  // Fill the PID object arrays holding the pointers to identified particle tracks
  //

  // Dynamic cast to ESD events (DO NOTHING for AOD events)
  AliESDEvent *event = dynamic_cast<AliESDEvent *>(InputEvent());
  if ( !event )  return;

  if(IsPbPb()) {
      fV0cuts->SetMode(AliESDv0KineCuts::kPurity,AliESDv0KineCuts::kPbPb);
  }
  else {
      fV0cuts->SetMode(AliESDv0KineCuts::kPurity,AliESDv0KineCuts::kPP);
  }


  // V0 selection
  // set event
  fV0cuts->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++)
    fV0tags[i] = 0;

  fNumTagsStored = numTracks;

  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){


    AliESDv0 *v0 = (AliESDv0 *) event->GetV0(iv0);
 
    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue; 
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    // v0 Armenteros plot (QA)
    Float_t armVar[2] = {0.0,0.0};
    fV0cuts->Armenteros(v0, armVar);

    if(fListQATRDV0&&fhArmenteros)fhArmenteros->Fill(armVar[0],armVar[1]);

    // fill the Object arrays
    // positive particles
    if( pdgP == -11){
	fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackP));
        fV0tags[iTrackP] = 11;
    }
    else if( pdgP == 211){
	fV0pions->Add((AliVTrack*)event->GetTrack(iTrackP));
        fV0tags[iTrackP] = 211;
    }
    else if( pdgP == 2212){
	fV0protons->Add((AliVTrack*)event->GetTrack(iTrackP));
	fV0tags[iTrackP] = 2212;
    }

    // negative particles
    if( pdgN == 11){
	fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackN));
        fV0tags[iTrackN] = -11;
    }
    else if( pdgN == -211){
	fV0pions->Add((AliVTrack*)event->GetTrack(iTrackN));
        fV0tags[iTrackN] = -211;
    }
    else if( pdgN == -2212){
	fV0protons->Add((AliVTrack*)event->GetTrack(iTrackN));
        fV0tags[iTrackN] = -2212;
    }


  }
}

//______________________________________________________________________________
Int_t AliTRDPIDTree::GetV0tag(Int_t trackIndex) const
{
  //
  // Get the tag for the corresponding trackIndex. Returns -99 in case of invalid index/tag list.
  //



  if (trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0tags) return -99;
  else
  {
      return fV0tags[trackIndex];
  }
}

//______________________________________________________________________________
void AliTRDPIDTree::ClearV0PIDlist(){

  //
  // Clear the PID object arrays
  //

  fV0electrons->Clear();
  fV0pions->Clear();
  fV0protons->Clear();

  delete[] fV0tags;
  fV0tags = 0;

  fNumTagsStored = 0;
}


//______________________________________________________________________________
void AliTRDPIDTree::SetupV0qa()
{
  //
  // Create the qa objects for V0 Kine cuts
  //
  
  fhArmenteros  = new TH2F("fhArmenteros","Armenteros plot",200,-1.,1.,200,0.,0.4);
  fListQATRDV0->Add(fhArmenteros);
 
}

//________________________________________________________________________
Double_t AliTRDPIDTree::GetPhi(AliESDtrack *const fTrack,Int_t iPl, Double_t &yposlayer, Double_t &thetalayer)
{
    //
    // extrapolate track to TRD radii and convert global phi angle to local coordinate system
    //

    Double_t phi=-999;
    yposlayer=-999;
    thetalayer=-999;
    // Phi at entrance of TRD
    Double_t xtrdbeg=AliTRDgeometry::GetXtrdBeg();
    Double_t xtrdend=AliTRDgeometry::GetXtrdEnd();
    Double_t x=xtrdbeg+iPl*(xtrdend-xtrdbeg)/6;

    const AliExternalTrackParam *tempparam = NULL;
    if(fTrack->GetOuterParam()) tempparam = fTrack->GetOuterParam();
    else if(fTrack->GetInnerParam()) tempparam = fTrack->GetInnerParam();



    if(tempparam){
	AliExternalTrackParam param(*tempparam);
	Bool_t isOk=param.PropagateTo(x,fESDEvent->GetMagneticField());
	if(isOk){
	  phi=param.Phi()-param.GetAlpha();
	  yposlayer= param.GetY();
	  thetalayer= param.Eta();
	}	
    }
    if(phi!=-999){
      if(phi<-TMath::Pi())phi+=2*TMath::Pi();
      if(phi>TMath::Pi())phi-=2*TMath::Pi();
    }

    return phi;
}


//________________________________________________________________________
Bool_t AliTRDPIDTree::PassTrackCuts(AliESDtrack *fESDTrack)
{
    //
    // check if tracks pass minimum quality critieria
    //

    if(!fESDTrack) return kFALSE;

    if(fhtrackCuts) fhtrackCuts->Fill(0);


    // DCA to PV
    Float_t dca[2];
    fESDTrack->GetImpactParameters(dca[0],dca[1]);
    if(dca[0]>5||dca[1]>10){
	if(fhtrackCuts) fhtrackCuts->Fill(1);
	return kFALSE;
    }

    // remove kinks
    if(fESDTrack->GetKinkIndex(0)>0){
      if(fhtrackCuts) fhtrackCuts->Fill(2);
      return kFALSE;
    }

    // eta cut
    if((TMath::Abs(fESDTrack->Eta()))>1.0){
        if(fhtrackCuts) fhtrackCuts->Fill(3);
	return kFALSE;
    }

    // Chi2
    Float_t tpcchi2=99;
    Int_t tpcnclusF=fESDTrack->GetTPCNclsF();
    if(tpcnclusF!=0) tpcchi2=(Float_t)fESDTrack->GetTPCchi2()/tpcnclusF;
    else tpcchi2=1000;
    if(tpcchi2 > 4){
	if(fhtrackCuts) fhtrackCuts->Fill(4);
	return kFALSE;
    }

    //TRD out
    if((fESDTrack->GetStatus()&AliVTrack::kTRDout)==0){
	if(fhtrackCuts) fhtrackCuts->Fill(5);
	return kFALSE;
    }

    // QA #TRD PID tracklets (no cut, just for QA)
    if(fESDTrack->GetTRDntrackletsPID()<4){
	if(fhtrackCuts) fhtrackCuts->Fill(6);
    }

    // QA Missing layers (no cut, just for QA)
    if(HasMissingLayer(fESDTrack)){
	if(fhtrackCuts) fhtrackCuts->Fill(7);
    }

    return kTRUE;
}

//_______________________________________________________________________
Bool_t AliTRDPIDTree::HasMissingLayer(const AliVTrack *fESDTrack){

    //
    // returns if slices missing
    //

    Int_t ntrl=0;
    for(Int_t jPl=0;jPl<AliVTrack::kTRDnPlanes;jPl++){
	Double_t signal=0;
	for(int isl= 0; isl<= 8;isl++){
	    Double_t sigsl=fESDTrack->GetTRDslice(jPl,isl);
	    if(sigsl>0)signal+=sigsl;
	}
        // if signal is missing, stop counting
	if(signal<=0||fESDTrack->GetTRDmomentum(jPl)<=0)break;
	ntrl++;
    }
    //compare with number of layers
    if(ntrl!=GetNTrackletsPID(fESDTrack))return kTRUE;

    return kFALSE;
}

//_______________________________________________________________________
Int_t AliTRDPIDTree::GetNTrackletsPID(const AliVTrack *fESDTrack) const
{
    //
    // returns the number of PID tracklets, re-calculated from dE/dx stored in slices
    //

    Int_t ntracklets=0;

    for(Int_t iPl=0;iPl<AliVTrack::kTRDnPlanes;iPl++){
	Double_t signal=0;
	for(int isl= 0; isl<= 8;isl++){
	    Double_t sigsl=fESDTrack->GetTRDslice(iPl,isl);
	    if(sigsl>0)signal+=sigsl;
	}
	if(signal>0){
	    ntracklets++;
	}
    }
    return ntracklets;
}
