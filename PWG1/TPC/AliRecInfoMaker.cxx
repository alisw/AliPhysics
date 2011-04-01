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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber                                                  //
//  Comparison macro for ESD                                                 //
//  responsible: 
//  marian.ivanov@cern.ch                                                    //
//
//

/* 
marian.ivanov@cern.ch 
Usage:
 


gSystem->Load("libPWG1.so");
//
AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
t2->Exec();


TFile f("cmpESDTracks.root");
TTree * tree = (TTree*)f.Get("ESDcmpTracks");

AliTreeDraw comp;
comp.SetTree(tree)



//
//some cuts definition
TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.01&&abs(MC.fVDist[2])<0.01")
//TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.5&&abs(MC.fVDist[2])<0.5")
//TCut citsin("citsin","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<3.9");
TCut citsin("citsin","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<5");
TCut csec("csec","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)>0.5");


TCut crec("crec","fReconstructed==1");
TCut cteta1("cteta1","abs(MC.fParticle.Theta()/3.1415-0.5)<0.25");
TCut cteta05("cteta05","abs(MC.fParticle.Theta()/3.1415-0.5)<0.1");

TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
TCut csens("csens","abs(sqrt(fVDist[0]**2+fVDist[1]**2)-170)<50");
TCut cmuon("cmuon","abs(MC.fParticle.fPdgCode==-13)");
TCut cchi2("cchi2","fESDtrack.fITSchi2MIP[0]<7.&&fESDtrack.fITSchi2MIP[1]<5.&&fESDtrack.fITSchi2MIP[2]<7.&&fESDtrack.fITSchi2MIP[3]<7.5&&fESDtrack.fITSchi2MIP[4]<6.")


//
//example
comp.T()->SetAlias("radius","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)");
comp.T()->SetAlias("direction","MC.fParticle.fVx*MC.fParticle.fPx+MC.fParticle.fVy*MC.fParticle.fPy");
comp.T()->SetAlias("decaydir","MC.fTRdecay.fX*MC.fTRdecay.fPx+MC.fTRdecay.fY*MC.fTRdecay.fPy");
comp.T()->SetAlias("theta","MC.fTrackRef.Theta()");
comp.T()->SetAlias("primdca","sqrt(RC.fITStrack.fD[0]**2+RC.fITStrack.fD[1]**2)");
comp.T()->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
comp.T()->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");


TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);
TH1F hfake("hfake","hfake",1000,0,150);
TProfile profp0("profp0","profp0",20,-0.4,0.9)

comp.DrawXY("fTPCinP0[3]","fTPCDelta[4]/fTPCinP1[3]","fReconstructed==1"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();
comp.fMean->Draw();  

comp.DrawXY("fITSinP0[3]","fITSDelta[4]/fITSinP1[3]","fReconstructed==1&&fITSOn"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();

comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn",20,0.2,1.5)
comp.fRes->Draw();

comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn&&fITSOn&&fESDtrack.fITSFakeRatio<0.1",10,0.2,1.5)
comp.fRes->Draw();
comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn&&fITSOn&&fESDtrack.fITSFakeRatio>0.1",10,0.2,1.5)
comp.fRes->Draw();

comp.T()->Draw("fESDtrack.fITSsignal/fESDtrack.fTPCsignal","fITSOn&&fTPCOn&&fESDtrack.fITSFakeRatio==0") 

TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);

TH2F * hdedx0 = new TH2F("dEdx0","dEdx0",100, 0,2,200,0,550); hdedx0->SetMarkerColor(1);
TH2F * hdedx1 = new TH2F("dEdx1","dEdx1",100, 0,2,200,0,550); hdedx1->SetMarkerColor(4);
TH2F * hdedx2 = new TH2F("dEdx2","dEdx2",100, 0,2,200,0,550); hdedx2->SetMarkerColor(3);
TH2F * hdedx3 = new TH2F("dEdx3","dEdx3",100, 0,2,200,0,550); hdedx3->SetMarkerColor(2);

comp.T()->Draw("fESDtrack.fITSsignal:MC.fParticle.P()>>dEdx0","fITSOn&&abs(fPdg)==211&&fITStrack.fN==6"+cprim) 
comp.T()->Draw("fESDtrack.fITSsignal:MC.fParticle.P()>>dEdx1","fITSOn&&abs(fPdg)==2212&&fITStrack.fN==6"+cprim) 
comp.T()->Draw("fESDtrack.fITSsignal:MC.fParticle.P()>>dEdx2","fITSOn&&abs(fPdg)==321&&fITStrack.fN==6"+cprim) 
comp.T()->Draw("fESDtrack.fITSsignal:MC.fParticle.P()>>dEdx3","fITSOn&&abs(fPdg)==11&&fITStrack.fN==6"+cprim) 


comp.T()->Draw("fESDtrack.fTRDsignal:MC.fParticle.P()>>dEdx0","fTRDOn&&abs(fPdg)==211&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.T()->Draw("fESDtrack.fTRDsignal:MC.fParticle.P()>>dEdx1","fTRDOn&&abs(fPdg)==2212&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.T()->Draw("fESDtrack.fTRDsignal:MC.fParticle.P()>>dEdx2","fTRDOn&&abs(fPdg)==321&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.T()->Draw("fESDtrack.fTRDsignal:MC.fParticle.P()>>dEdx3","fTRDOn&&abs(fPdg)==11&&fTRDtrack.fN>40&&fStatus[2]>1") 

comp.T()->Draw("fESDtrack.fTPCsignal:fTPCinP0[4]>>dEdx0","fTPCOn&&abs(fPdg)==211&&fESDtrack.fTPCncls>180&&fESDtrack.fTPCsignal>10"+cteta1); 
comp.T()->Draw("fESDtrack.fTPCsignal:fTPCinP0[4]>>dEdx1","fTPCOn&&abs(fPdg)==2212&&fESDtrack.fTPCncls>180&&fESDtrack.fTPCsignal>10"+cteta1); 
comp.T()->Draw("fESDtrack.fTPCsignal:fTPCinP0[4]>>dEdx2","fTPCOn&&abs(fPdg)==321&&fESDtrack.fTPCncls>180&&fESDtrack.fTPCsignal>10"+cteta1); 
comp.T()->Draw("fESDtrack.fTPCsignal:fTPCinP0[4]>>dEdx3","fTPCOn&&abs(fPdg)==11&&fESDtrack.fTPCncls>180&&fESDtrack.fTPCsignal>10"+cteta1); 

hdedx3->SetXTitle("P(GeV/c)");
hdedx3->SetYTitle("dEdx(unit)");
hdedx3->Draw(); hdedx1->Draw("same"); hdedx2->Draw("same"); hdedx0->Draw("same");

comp.DrawXY("fITSinP0[3]","fITSPools[4]","fReconstructed==1&&fPdg==-211&&fITSOn"+cprim,"1",4,0.2,1.0,-8,8)

TProfile prof("prof","prof",10,0.5,5);




*/


#include <stdio.h>
#include <string.h>
//ROOT includes
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TGeoManager.h"
//#include "Getline.h"
//
//ALIROOT includes
//
#include "AliRun.h"
#include "AliESDtrack.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTPCseed.h"
#include "AliITStrackMI.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDkink.h"
#include "AliESDv0.h"
#include "AliV0.h"
//
#include "AliTreeDraw.h"
#include "AliMCInfo.h"
#include "AliGenKinkInfo.h"
#include "AliGenV0Info.h"


#include "AliESDRecInfo.h"
#include "AliESDRecV0Info.h"
#include "AliESDRecKinkInfo.h"
#include "AliRecInfoMaker.h"



ClassImp(AliRecInfoMaker)




AliTPCParam * AliRecInfoMaker::GetTPCParam(){
  //
  // create TPC param
  //
  AliTPCParamSR * par = new AliTPCParamSR;
  par->Update();
  return par;
}



void  AliRecInfoMaker::MakeAliases(TTree * tree)
{
  //
  // aliases definition
  //
  tree->SetAlias("radius","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)");
  tree->SetAlias("direction","MC.fParticle.fVx*MC.fParticle.fPx+MC.fParticle.fVy*MC.fParticle.fPy");
  tree->SetAlias("decaydir","MC.fTRdecay.fX*MC.fTRdecay.fPx+MC.fTRdecay.fY*MC.fTRdecay.fPy");
  tree->SetAlias("theta","MC.fTrackRef.Theta()");
  tree->SetAlias("primdca","sqrt(RC.fITStrack.fD[0]**2+RC.fITStrack.fD[1]**2)");
  tree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
  tree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
  
  tree->SetAlias("trddedx","(RC.fESDtrack.fTRDsignals[0]+RC.fESDtrack.fTRDsignals[1]+RC.fESDtrack.fTRDsignals[2]+RC.fESDtrack.fTRDsignals[3]+RC.fESDtrack.fTRDsignals[4]+RC.fESDtrack.fTRDsignals[5])/6.");
  
  tree->SetAlias("dtofmc2","fESDtrack.fTrackTime[2]-(10^12*MC.fTOFReferences[0].fTime)");
  tree->SetAlias("dtofrc2","(fESDtrack.fTrackTime[2]-fESDtrack.fTOFsignal)");

  tree->SetAlias("psum","fESDtrack.fTOFr[4]+fESDtrack.fTOFr[3]+fESDtrack.fTOFr[2]+fESDtrack.fTOFr[1]+fESDtrack.fTOFr[0]");
  tree->SetAlias("P0","fESDtrack.fTOFr[0]/psum");
  tree->SetAlias("P1","fESDtrack.fTOFr[1]/psum");
  tree->SetAlias("P2","fESDtrack.fTOFr[2]/psum");
  tree->SetAlias("P3","fESDtrack.fTOFr[3]/psum");
  tree->SetAlias("P4","fESDtrack.fTOFr[4]/psum");
  tree->SetAlias("MaxP","max(max(max(P0,P1),max(P2,P3)),P4)");
}


////////////////////////////////////////////////////////////////////////
AliRecInfoMaker::AliRecInfoMaker(const char* fnGenTracks,
				 const char* fnCmp,
				 const char* fnGalice,
				 Int_t nEvents, Int_t firstEvent):
  
  fEventNr(0),                 //! current event number
  fNEvents(0),                 //! number of events to process
  fFirstEventNr(0),            //! first event to process
  fFileCmp(0),                //! output file with cmp tracks
  fTreeCmp(0),                //! output tree with cmp tracks
  fTreeCmpKinks(0),                //! output tree with cmp Kinks
  fTreeCmpV0(0),                //! output tree with cmp V0
  //
  fFileGenTracks(0),                //! input files with generated tracks   
  fTreeGenTracks(0),           //! tree with generated tracks
  fTreeGenKinks(0),            // tree with gen kinks
  fTreeGenV0(0),            // tree with gen V0
  //
  fLoader(0),         //! pointer to the run loader
  //
  fIndexRecTracks(0),         //! index of particle label in the TreeT_ESD
  fFakeRecTracks(0),          //! number of fake tracks
  fMultiRecTracks(0),         //! number of multiple reconstructions
  //
  fIndexRecKinks(0),         //! index of particle label in treeesd
  fMultiRecKinks(0),         //! number of multiple reconstructions
  fSignedKinks(0),           //! indicator that kink was not fake
  //
  fIndexRecV0(0),         //! index of particle label in treeesd
  fMultiRecV0(0),         //! number of multiple reconstructions
  fSignedV0(0),                //! indicator that kink was not fake
  //
  fRecArray(0),           // container with rec infos
  fEvent(0),             //!event
  fESDfriend(0),              //!event friend
  //
  fParamTPC(0),         //! AliTPCParam
  fNParticles(0),              //! number of particles in the input tree genTracks
  fDebug(0),                   //! debug flag  
  fNextTreeGenEntryToRead(0),    //! last entry already read from genTracks tree
  fNextKinkToRead(0),            //! last entry already read from genKinks tree
  fNextV0ToRead(0),            //! last entry already read from genV0 tree
  //
  fMCInfo(0),           //! MC information writen per particle
  fGenKinkInfo(0),      //! MC information writen per Kink
  fGenV0Info(0),      //! MC information writen per Kink
  fRecInfo(0),          //! Rec. information writen per particle
  fFriend(0),          //! friend track
  fRecKinkInfo(0),    //! reconstructed kink info
  fRecV0Info(0)    //! reconstructed kink info
{
  // AliRecInfoMaker - connencts the MC information with reconstructed information
  // fnGenTracks  - file with MC to be created before using AliGenInfoMaker
  // fnCmp        - file name  to be created  
  // fnGalice     - file with Loaders - usualy galice.root 
  //  
  // nEvent       - number of event s to be analyzed
  // AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
  //


  Reset();
  //  fFnGenTracks = fnGenTracks;
  //  fFnCmp = fnCmp;

  memset(fFnGenTracks,0,sizeof(fFnGenTracks));
  memset(fFnCmp,0,sizeof(fFnCmp));

  snprintf(fFnGenTracks,1000,"%s",fnGenTracks);
  snprintf(fFnCmp,1000,"%s",fnCmp);

  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    //delete AliRunLoader::Instance();
    delete gAlice;
    gAlice = 0x0;
  }
  if (fLoader->LoadgAlice()){
    cerr<<"Error occured while l"<<endl;
  }
  Int_t nall = fLoader->GetNumberOfEvents();
  if (nEvents==0) {
    nEvents =nall;
    fNEvents=nall;
    fFirstEventNr=0;
  }    

  if (nall<=0){
    cerr<<"no events available"<<endl;
    fEventNr = 0;
    return;
  }
  if (firstEvent+nEvents>nall) {
    fEventNr = nall-firstEvent;
    cerr<<"restricted number of events availaible"<<endl;
  }
  TGeoManager::Import("geometry.root");


}
////////////////////////////////////////////////////////////////////////
AliRecInfoMaker::AliRecInfoMaker(const AliRecInfoMaker& /*info*/):
  
  fEventNr(0),                 //! current event number
  fNEvents(0),                 //! number of events to process
  fFirstEventNr(0),            //! first event to process
  fFileCmp(0),                //! output file with cmp tracks
  fTreeCmp(0),                //! output tree with cmp tracks
  fTreeCmpKinks(0),                //! output tree with cmp Kinks
  fTreeCmpV0(0),                //! output tree with cmp V0
  //
  fFileGenTracks(0),                //! input files with generated tracks   
  fTreeGenTracks(0),           //! tree with generated tracks
  fTreeGenKinks(0),            // tree with gen kinks
  fTreeGenV0(0),            // tree with gen V0
  //
  fLoader(0),         //! pointer to the run loader
  //
  fIndexRecTracks(0),         //! index of particle label in the TreeT_ESD
  fFakeRecTracks(0),          //! number of fake tracks
  fMultiRecTracks(0),         //! number of multiple reconstructions
  //
  fIndexRecKinks(0),         //! index of particle label in treeesd
  fMultiRecKinks(0),         //! number of multiple reconstructions
  fSignedKinks(0),           //! indicator that kink was not fake
  //
  fIndexRecV0(0),         //! index of particle label in treeesd
  fMultiRecV0(0),         //! number of multiple reconstructions
  fSignedV0(0),                //! indicator that kink was not fake
  //
  fRecArray(0),           // container with rec infos
  fEvent(0),             //!event
  fESDfriend(0),              //!event friend
  //
  fParamTPC(0),         //! AliTPCParam
  fNParticles(0),              //! number of particles in the input tree genTracks
  fDebug(0),                   //! debug flag  
  fNextTreeGenEntryToRead(0),    //! last entry already read from genTracks tree
  fNextKinkToRead(0),            //! last entry already read from genKinks tree
  fNextV0ToRead(0),            //! last entry already read from genV0 tree
  //
  fMCInfo(0),           //! MC information writen per particle
  fGenKinkInfo(0),      //! MC information writen per Kink
  fGenV0Info(0),      //! MC information writen per Kink
  fRecInfo(0),          //! Rec. information writen per particle
  fFriend(0),          //! friend track
  fRecKinkInfo(0),    //! reconstructed kink info
  fRecV0Info(0)    //! reconstructed kink info
{
  //
  // Dummy copu constructor
  //
  memset(fFnGenTracks,0,sizeof(fFnGenTracks));
  memset(fFnCmp,0,sizeof(fFnCmp));
}




////////////////////////////////////////////////////////////////////////
AliRecInfoMaker::~AliRecInfoMaker()
{
  //
  // Destructor
  //
  if (fLoader) {
    delete fLoader;
  }
}

//////////////////////////////////////////////////////////////
Int_t AliRecInfoMaker::SetIO()
{
  //
  // SetIO  - Create the input trees
  //
  CreateTreeCmp();
  if (!fTreeCmp) return 1;
  fParamTPC = GetTPCParam();
  //
  if (!ConnectGenTree()) {
    cerr<<"Cannot connect tree with generated tracks"<<endl;
    return 1;
  }
  return 0;
}

//////////////////////////////////////////////////////////////

Int_t AliRecInfoMaker::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  //
  TFile f("AliESDs.root");
  //
 
  TTree* tree = (TTree*) f.Get("esdTree");
  tree->SetBranchStatus("*",1);
  fEvent = new AliESDEvent;
  
  if (tree->GetBranch("ESD")){
    //    tree->SetBranchAddress("ESD", &fEvent);
    // tree->SetBranchAddress("ESDfriend.",&fESDfriend);
    // tree->GetEntry(eventNr);
    // fEvent->SetESDfriend(fESDfriend);    
  }else{
    fEvent->ReadFromTree(tree);
    fESDfriend = (AliESDfriend*)fEvent->FindListObject("AliESDfriend"); 
    tree->GetEntry(eventNr);
    fEvent->SetESDfriend(fESDfriend); 
  }    
  


  if (!fEvent) return 1;

  return 0;
}



////////////////////////////////////////////////////////////////////////
void AliRecInfoMaker::Reset()
{
  //
  // Reset the class
  //
  fEventNr = 0;
  fNEvents = 0;
  fTreeCmp = 0;
  fTreeCmpKinks =0;
  fTreeCmpV0 =0;
  //  fFnCmp = "cmpTracks.root";
  fFileGenTracks = 0;
  fDebug = 0;
  //
  fParamTPC = 0;
  fEvent =0;
}

////////////////////////////////////////////////////////////////////////
Int_t AliRecInfoMaker::Exec(Int_t nEvents, Int_t firstEventNr)
{
  //
  // Exec comparison for subrange of events
  //
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t AliRecInfoMaker::Exec()
{
  //
  // Exec comparison
  //
  TStopwatch timer;
  timer.Start();

  if (SetIO()==1) 
    return 1;
   
  fNextTreeGenEntryToRead = 0;
  fNextKinkToRead = 0;
  fNextV0ToRead   =0;
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    fEventNr = eventNr;
    SetIO(fEventNr);
    fNParticles = gAlice->GetEvent(fEventNr);    

    fIndexRecTracks = new Short_t[fNParticles*20];  //write at maximum 4 tracks corresponding to particle
    fIndexRecKinks  = new Short_t[fNParticles*20];  //write at maximum 20 tracks corresponding to particle
    fIndexRecV0  = new Short_t[fNParticles*20];  //write at maximum 20 tracks corresponding to particle

    fFakeRecTracks = new Short_t[fNParticles];
    fMultiRecTracks = new Short_t[fNParticles];
    fMultiRecKinks = new Short_t[fNParticles];
    fMultiRecV0 = new Short_t[fNParticles];

    for (Int_t i = 0; i<fNParticles; i++) {
      for (Int_t j=0;j<20;j++){
	fIndexRecTracks[20*i+j] = -1;
	fIndexRecKinks[20*i+j]  = -1;
	fIndexRecV0[20*i+j]  = -1;
      }
      fFakeRecTracks[i] = 0;
      fMultiRecTracks[i] = 0;
      fMultiRecKinks[i] = 0;
      fMultiRecV0[i] = 0;      
    }
  
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\tStart loop over TreeT"<<endl;
    if (TreeTLoop()>0) return 1;

    if (fDebug>2) cout<<"\tStart loop over tree genTracks"<<endl;
    if (TreeGenLoop(eventNr)>0) return 1;
    //BuildKinkInfo0(eventNr);
    BuildV0Info(eventNr); // no V0 info for a moment
    fRecArray->Delete();

    if (fDebug>2) cout<<"\tEnd loop over tree genTracks"<<endl;

    delete [] fIndexRecTracks;
    delete [] fIndexRecKinks;
    delete [] fIndexRecV0;
    delete [] fFakeRecTracks;
    delete [] fMultiRecTracks;
    delete [] fMultiRecKinks;
    delete [] fMultiRecV0;
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  timer.Stop();
  timer.Print();
  return 0;

}
////////////////////////////////////////////////////////////////////////
Bool_t AliRecInfoMaker::ConnectGenTree()
{
//
// connect all branches from the genTracksTree
// use the same variables as for the new cmp tree, it may work
//
  fFileGenTracks = TFile::Open(fFnGenTracks,"READ");
  if (!fFileGenTracks) {
    cerr<<"Error in ConnectGenTree: cannot open file "<<fFnGenTracks<<endl;
    return kFALSE;
  }
  fTreeGenTracks = (TTree*)fFileGenTracks->Get("genTracksTree");
  if (!fTreeGenTracks) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    return kFALSE;
  }
  //
  fMCInfo = new AliMCInfo;
  fTreeGenTracks->SetBranchAddress("MC",&fMCInfo);
  //
  //
  fTreeGenKinks = (TTree*)fFileGenTracks->Get("genKinksTree");
  if (!fTreeGenKinks) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    //return kFALSE;
  }
  else{
    fGenKinkInfo = new AliGenKinkInfo;
    fTreeGenKinks->SetBranchAddress("MC",&fGenKinkInfo);
  }

  fTreeGenV0 = (TTree*)fFileGenTracks->Get("genV0Tree");
  if (!fTreeGenV0) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    //return kFALSE;
  }
  else{
    fGenV0Info = new AliGenV0Info;
    fTreeGenV0->SetBranchAddress("MC",&fGenV0Info);
  }
  //
  if (fDebug > 1) {
    cout<<"Number of gen. tracks with TR: "<<fTreeGenTracks->GetEntries()<<endl;
  }
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
void AliRecInfoMaker::CreateTreeCmp() 
{
  //
  // Create file and tree with comparison information 
  //
  fFileCmp = TFile::Open(fFnCmp,"RECREATE");
  if (!fFileCmp) {
    cerr<<"Error in CreateTreeCmp: cannot open file "<<fFnCmp<<endl;
    return;
  }
  //
  //
  fTreeCmp    = new TTree("ESDcmpTracks","ESDcmpTracks");
  fMCInfo = new AliMCInfo;
  fRecInfo = new AliESDRecInfo;
  AliESDtrack * esdTrack = new AliESDtrack; 
  //  AliITStrackMI * itsTrack = new AliITStrackMI;  
  fTreeCmp->Branch("MC","AliMCInfo",&fMCInfo,256000);
  fTreeCmp->Branch("RC","AliESDRecInfo",&fRecInfo,256000);
  //  fTreeCmp->Branch("ITS","AliITStrackMI",&itsTrack);
  delete esdTrack;
  //
  //
  fTreeCmpKinks    = new TTree("ESDcmpKinks","ESDcmpKinks"); 
  fGenKinkInfo     = new AliGenKinkInfo;
  fRecKinkInfo     = new AliESDRecKinkInfo;
  fTreeCmpKinks->Branch("MC.","AliGenKinkInfo",&fGenKinkInfo,256000);
  fTreeCmpKinks->Branch("RC.","AliESDRecKinkInfo",&fRecKinkInfo,256000);
  //
  //
  fTreeCmpV0       = new TTree("ESDcmpV0","ESDcmpV0"); 
  fGenV0Info     = new AliGenV0Info;
  fRecV0Info     = new AliESDRecV0Info;
  fTreeCmpV0->Branch("MC.","AliGenV0Info",   &fGenV0Info,256000);
  fTreeCmpV0->Branch("RC.","AliESDRecV0Info",&fRecV0Info,256000);
  //
  fTreeCmp->AutoSave(); 
  fTreeCmpKinks->AutoSave(); 
  fTreeCmpV0->AutoSave(); 
}

////////////////////////////////////////////////////////////////////////
void AliRecInfoMaker::CloseOutputFile()  
{
  //
  // Close output file
  //

  if (!fFileCmp) {
    cerr<<"File "<<fFnCmp<<" not found as an open file."<<endl;
    return;
  }
  fFileCmp->cd();
  fTreeCmp->Write();    
  delete fTreeCmp;
  
  fFileCmp->Close();
  delete fFileCmp;
  return;
}
////////////////////////////////////////////////////////////////////////

TVector3 AliRecInfoMaker::TR2Local(AliTrackReference *trackRef,
			    AliTPCParam *paramTPC) {

  //
  // Transform position to the local coord frame
  //
  
  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2Ideal(x,index);
  return TVector3(x);
}
////////////////////////////////////////////////////////////////////////

Int_t AliRecInfoMaker::TreeTLoop()
{
  //
  // loop over all ESD reconstructed tracks and store info in memory
  //
  // + loop over all reconstructed kinks
  TStopwatch  timer;
  timer.Start();
  //  
  Int_t nEntries = (Int_t)fEvent->GetNumberOfTracks();  
  Int_t nKinks = (Int_t) fEvent->GetNumberOfKinks();
  Int_t nV0MIs = (Int_t) fEvent->GetNumberOfV0s();
  fSignedKinks = new Short_t[nKinks];
  fSignedV0    = new Short_t[nV0MIs];
  //
  // load kinks to the memory
  for (Int_t i=0; i<nKinks;i++){
    //    AliESDkink * kink =
    fEvent->GetKink(i);
    fSignedKinks[i]=0;
  }
  //
  for (Int_t i=0; i<nV0MIs;i++){
    //AliV0 * v0MI = 
    (AliV0*)fEvent->GetV0(i);
    fSignedV0[i]=0;
  }
  
  //
  //
  AliESDtrack * track=0;
  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    track = (AliESDtrack*)fEvent->GetTrack(iEntry);
    //
    Int_t label = track->GetLabel();
    Int_t absLabel = abs(label);
    if (absLabel < fNParticles) {
      //      fIndexRecTracks[absLabel] =  iEntry;
      if (label < 0) fFakeRecTracks[absLabel]++;      
      if (fMultiRecTracks[absLabel]>0){
	if (fMultiRecTracks[absLabel]<20)
	  fIndexRecTracks[absLabel*20+fMultiRecTracks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecTracks[absLabel*20] =  iEntry;
      fMultiRecTracks[absLabel]++;
    }
  }
  // sort reconstructed kinks  
  //
  AliESDkink * kink=0;
  for (Int_t iEntry=0; iEntry<nKinks;iEntry++){
    kink = (AliESDkink*)fEvent->GetKink(iEntry);
    if (!kink) continue;
    //
    Int_t label0 = TMath::Abs(kink->GetLabel(0));
    Int_t label1 = TMath::Abs(kink->GetLabel(1));
    Int_t absLabel = TMath::Min(label0,label1);
    if (absLabel < fNParticles) {
      if (fMultiRecKinks[absLabel]>0){
	if (fMultiRecKinks[absLabel]<20)
	  fIndexRecKinks[absLabel*20+fMultiRecKinks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecKinks[absLabel*20] =  iEntry;
      fMultiRecKinks[absLabel]++;
    }
  }  
  // --sort reconstructed V0
  //
  AliV0 * v0MI=0;
  for (Int_t iEntry=0; iEntry<nV0MIs;iEntry++){
    v0MI = (AliV0*)fEvent->GetV0(iEntry);
    if (!v0MI) continue;
    //
    //
    //
    //Int_t label0 = TMath::Abs(v0MI->GetLabel(0));
    //Int_t label1 = TMath::Abs(v0MI->GetLabel(1));
    AliESDtrack * trackn = fEvent->GetTrack((v0MI->GetNindex()));
    AliESDtrack * trackp = fEvent->GetTrack((v0MI->GetPindex()));
    Int_t labels[2]={-1,-1};
    labels[0] = (trackn==0) ? -1 : TMath::Abs(trackn->GetLabel()); 
    labels[1] = (trackp==0) ? -1 : TMath::Abs(trackp->GetLabel()); 
    //
    for (Int_t i=0;i<2;i++){
      Int_t absLabel =  TMath::Abs(labels[i]);
      if (absLabel < fNParticles) {
	if (fMultiRecV0[absLabel]>0){
	  if (fMultiRecV0[absLabel]<20)
	    fIndexRecV0[absLabel*20+fMultiRecV0[absLabel]] =  iEntry; 	
	}
	else      
	  fIndexRecV0[absLabel*20] =  iEntry;
	fMultiRecV0[absLabel]++;
      }
    }
  }  


  printf("Time spended in TreeTLoop\n");
  timer.Print();
  
  if (fDebug > 2) cerr<<"end of TreeTLoop"<<endl;  
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t AliRecInfoMaker::TreeGenLoop(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextTreeGenEntryToRead;
  Double_t nParticlesTR = fTreeGenTracks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextTreeGenEntryToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextTreeGenEntryToRead<<endl;
  TBranch * branch = fTreeCmp->GetBranch("RC");
  //  TBranch * branchF = fTreeCmp->GetBranch("F");
  
  branch->SetAddress(&fRecInfo); // set all pointers
  fRecArray = new TObjArray(fNParticles);
  AliESDtrack dummytrack;  //
  AliESDfriendTrack dummytrackF;  //

  while (entry < nParticlesTR) {
    fTreeGenTracks->GetEntry(entry);
    entry++;
    if (eventNr < fMCInfo->fEventNr) continue;
    if (eventNr > fMCInfo->fEventNr) continue;
    if (fMCInfo->GetCharge()==0) continue;
    //
    fNextTreeGenEntryToRead = entry-1;
    if (fDebug > 2 && fMCInfo->fLabel < 10) {
      cerr<<"Fill track with a label "<<fMCInfo->fLabel<<endl;
    }
    //    if (fMCInfo->fNTPCRef<1) continue; // not TPCref
    //
    fRecInfo->Reset();
    AliESDtrack * track=0;
    fRecInfo->fReconstructed =0;
    TVector3 local = TR2Local(&(fMCInfo->fTrackRef),fParamTPC);
    local.GetXYZ(fRecInfo->fTRLocalCoord);	
    //
    if (fIndexRecTracks[fMCInfo->fLabel*20] >= 0) {
      track= (AliESDtrack*)fEvent->GetTrack(fIndexRecTracks[fMCInfo->fLabel*20]);
      if(!track) continue;
      //
      //
      // find nearest track if multifound
      //Int_t sign = Int_t(track->GetSign()*fMCInfo->fCharge);
      //
      Int_t status = 0;
      if  ((track->GetStatus()&AliESDtrack::kITSrefit)>0) status++;
      if  ((track->GetStatus()&AliESDtrack::kTPCrefit)>0) status++;
      if  ((track->GetStatus()&AliESDtrack::kTRDrefit)>0) status++;

      //
      if (fIndexRecTracks[fMCInfo->fLabel*20+1]>0){
	//
	Double_t p[3];
	track->GetInnerPxPyPz(p);
	Float_t maxp = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
	//
	for (Int_t i=1;i<20;i++){
	  if (fIndexRecTracks[fMCInfo->fLabel*20+i]>=0){
	    AliESDtrack * track2 = (AliESDtrack*)fEvent->GetTrack(fIndexRecTracks[fMCInfo->fLabel*20+i]);
	    if (!track2) continue;
	    //Int_t sign2 = track->GetSign()*fMCInfo->fCharge; //	    
	    //if (sign2<0) continue;
	    track2->GetInnerPxPyPz(p);
	    Float_t mom = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
	    /*
	    if (sign<0){
	      sign = sign2;
	      track = track2;
	      maxp = mom;
	      continue;
	    }
	    */
	    //
	    Int_t status2 = 0;
	    if  ((track2->GetStatus()&AliESDtrack::kITSrefit)>0) status2++;
	    if  ((track2->GetStatus()&AliESDtrack::kTPCrefit)>0) status2++;
	    if  ((track2->GetStatus()&AliESDtrack::kTRDrefit)>0) status2++;
	    if (status2<status) continue;
	    //
	    if (mom<maxp) continue;
	    maxp = mom;
	    track = track2;
	    //
	  }
	}
      }	
      //
      if (track) {
	fRecInfo->SetESDtrack(track);
      }else{
	fRecInfo->SetESDtrack(&dummytrack);
      }
      //

      fRecInfo->fReconstructed = 1;
      fRecInfo->fFake     = fFakeRecTracks[fMCInfo->fLabel];
      fRecInfo->fMultiple = fMultiRecTracks[fMCInfo->fLabel];
      //
      fRecInfo->Update(fMCInfo,fParamTPC,kTRUE);          
    }
    else{
      fRecInfo->SetESDtrack(&dummytrack);
      fRecInfo->Update(fMCInfo,fParamTPC,kFALSE);
    }
    fRecArray->AddAt(new AliESDRecInfo(*fRecInfo),fMCInfo->fLabel);
    fTreeCmp->Fill();
  }
  fTreeCmp->AutoSave();
  printf("Time spended in TreeGenLoop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of TreeGenLoop"<<endl;

  return 0;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Int_t AliRecInfoMaker::BuildKinkInfo0(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextKinkToRead;
  Double_t nParticlesTR = fTreeGenKinks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextKinkToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextKinkToRead<<endl;
  //
  TBranch * branch = fTreeCmpKinks->GetBranch("RC.");
  branch->SetAddress(&fRecKinkInfo); // set all pointers
  
  //
  while (entry < nParticlesTR) {
    fTreeGenKinks->GetEntry(entry);
    entry++;
    if (eventNr < fGenKinkInfo->GetMinus().fEventNr) continue;
    if (eventNr > fGenKinkInfo->GetMinus().fEventNr) continue;;
    //
    fNextKinkToRead = entry-1;
    //
    //
    AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(fGenKinkInfo->GetMinus().fLabel);
    if(!fRecInfo1) continue;
    AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(fGenKinkInfo->GetPlus().fLabel);
    if(!fRecInfo2) continue;
    fRecKinkInfo->fT1 = (*fRecInfo1);
    fRecKinkInfo->fT2 = (*fRecInfo2);
    fRecKinkInfo->fStatus =0;
    if (fRecInfo1 && fRecInfo1->fTPCOn) fRecKinkInfo->fStatus+=1;
    if (fRecInfo2 && fRecInfo2->fTPCOn) fRecKinkInfo->fStatus+=2;
    if (fRecKinkInfo->fStatus==3&&fRecInfo1->fSign!=fRecInfo2->fSign) fRecKinkInfo->fStatus*=-1;
    
    if (fRecKinkInfo->fStatus==3){
      fRecKinkInfo->Update();    
    }
    Int_t label =  TMath::Min(fGenKinkInfo->GetMinus().fLabel,fGenKinkInfo->GetPlus().fLabel);
    Int_t label2 = TMath::Max(fGenKinkInfo->GetMinus().fLabel,fGenKinkInfo->GetPlus().fLabel);
    
    AliESDkink *kink=0;
    fRecKinkInfo->fRecStatus   =0;
    fRecKinkInfo->fMultiple    = fMultiRecKinks[label];
    fRecKinkInfo->fKinkMultiple=0;
    //
    if (fMultiRecKinks[label]>0){

      //      for (Int_t j=0;j<TMath::Min(fMultiRecKinks[label],100);j++){
      for (Int_t j=TMath::Min(fMultiRecKinks[label],Short_t(20))-1;j>=0;j--){
	Int_t index = fIndexRecKinks[label*20+j];
	//AliESDkink *kink2  = (AliESDkink*)fKinks->At(index);
	AliESDkink *kink2  = (AliESDkink*)fEvent->GetKink(index);
	if (TMath::Abs(kink2->GetLabel(0))==label &&TMath::Abs(kink2->GetLabel(1))==label2) {
	  fRecKinkInfo->fKinkMultiple++;
	  fSignedKinks[index]=1;
	  Int_t c0=0;
	  if (kink){
	    //	    if (kink->fTRDOn) c0++;
	    //if (kink->fITSOn) c0++;
	    if (kink->GetStatus(2)>0) c0++;
	    if (kink->GetStatus(0)>0) c0++;
	  }
	  Int_t c2=0;
	  //	  if (kink2->fTRDOn) c2++;
	  //if (kink2->fITSOn) c2++;
	  if (kink2->GetStatus(2)>0) c2++;
	  if (kink2->GetStatus(0)>0) c2++;

	  if (c2<c0) continue;
	  kink =kink2;
	}
	if (TMath::Abs(kink2->GetLabel(1))==label &&TMath::Abs(kink2->GetLabel(0))==label2) {
	  fRecKinkInfo->fKinkMultiple++;
	  fSignedKinks[index]=1;
	  Int_t c0=0;
	  if (kink){
	    //if (kink->fTRDOn) c0++;
	    //if (kink->fITSOn) c0++;
	    if (kink->GetStatus(2)>0) c0++;
	    if (kink->GetStatus(0)>0) c0++;

	  }
	  Int_t c2=0;
	  //	  if (kink2->fTRDOn) c2++;
	  //if (kink2->fITSOn) c2++;
	  if (kink2->GetStatus(2)>0) c2++;
	  if (kink2->GetStatus(0)>0) c2++;

	  if (c2<c0) continue;
	  kink =kink2;
	}
      }
    }
    if (kink){
      fRecKinkInfo->fKink = *kink;
      fRecKinkInfo->fRecStatus=1;
    }
    fTreeCmpKinks->Fill();
  }
  //  Int_t nkinks = fKinks->GetEntriesFast();
  Int_t nkinks = fEvent->GetNumberOfKinks();
  for (Int_t i=0;i<nkinks;i++){
    if (fSignedKinks[i]==0){
      //      AliESDkink *kink  = (AliESDkink*)fKinks->At(i);
      AliESDkink *kink  = (AliESDkink*)fEvent->GetKink(i);
      if (!kink) continue;
      //
      fRecKinkInfo->fKink = *kink;
      fRecKinkInfo->fRecStatus =-2;
      //
      AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(kink->GetLabel(0)));
      AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(kink->GetLabel(1)));
      if (fRecInfo1 && fRecInfo2){
	fRecKinkInfo->fT1 = (*fRecInfo1);
	fRecKinkInfo->fT2 = (*fRecInfo2);
	fRecKinkInfo->fRecStatus =-1;
      }
      fTreeCmpKinks->Fill();
    }
  }


  fTreeCmpKinks->AutoSave();
  printf("Time spended in BuilKinkInfo Loop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of BuildKinkInfo Loop"<<endl;
  return 0;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



Int_t AliRecInfoMaker::BuildV0Info(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  static TDatabasePDG pdgtable;

  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextV0ToRead;
  Double_t nParticlesTR = fTreeGenV0->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextV0ToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextV0ToRead<<endl;
  //
  TBranch * branch = fTreeCmpV0->GetBranch("RC.");
  branch->SetAddress(&fRecV0Info); // set all pointers
  const AliESDVertex * esdvertex = fEvent->GetVertex();
  Float_t vertex[3]= {esdvertex->GetXv(), esdvertex->GetYv(),esdvertex->GetZv()};
  
  //
  while (entry < nParticlesTR) {
    fTreeGenV0->GetEntry(entry);
    entry++;
    fRecV0Info->Reset();  //reset all variables
    if (eventNr < fGenV0Info->GetMinus().fEventNr) continue;
    if (eventNr > fGenV0Info->GetMinus().fEventNr) continue;;
    //
    fNextV0ToRead = entry-1;
    //
    //
    AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(fGenV0Info->GetMinus().fLabel);
    AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(fGenV0Info->GetPlus().fLabel);
    if (fGenV0Info->GetMinus().fCharge*fGenV0Info->GetPlus().fCharge>0) continue;  // interactions
    if (!fRecInfo1 || !fRecInfo2) continue;
    fRecV0Info->fT1 = (*fRecInfo1);
    fRecV0Info->fT2 = (*fRecInfo2);
    fRecV0Info->fV0Status =0;
    if (fRecInfo1 && fRecInfo1->fStatus[1]>0) fRecV0Info->fV0Status+=1;
    if (fRecInfo2 && fRecInfo2->fStatus[1]>0) fRecV0Info->fV0Status+=2;

    if (fRecV0Info->fV0Status==3&&fRecInfo1->fSign==fRecInfo2->fSign) fRecV0Info->fV0Status*=-1;


    if (abs(fRecV0Info->fV0Status)==3){
      fRecV0Info->Update(vertex);
      {
	//
	// TPC V0 Info
	Double_t x,alpha, param[5],cov[15];
	if ( fRecV0Info->fT1.GetESDtrack()->GetInnerParam() && fRecV0Info->fT2.GetESDtrack()->GetInnerParam()){
	  fRecV0Info->fT1.GetESDtrack()->GetInnerExternalParameters(alpha,x,param);
	  fRecV0Info->fT1.GetESDtrack()->GetInnerExternalCovariance(cov);
	  AliExternalTrackParam paramP(x,alpha,param,cov);
	  //
	  fRecV0Info->fT2.GetESDtrack()->GetInnerExternalParameters(alpha,x,param);
	  fRecV0Info->fT2.GetESDtrack()->GetInnerExternalCovariance(cov);
	  AliExternalTrackParam paramM(x,alpha,param,cov);
	  //
	  fRecV0Info->fV0tpc->SetParamN(paramM);
	  fRecV0Info->fV0tpc->SetParamP(paramP);
	  Double_t pid1[5],pid2[5];
	  fRecV0Info->fT1.GetESDtrack()->GetESDpid(pid1);
	  fRecV0Info->fT1.GetESDtrack()->GetESDpid(pid2);
	  //
	  //fRecV0Info->fV0tpc.UpdatePID(pid1,pid2);
	  fRecV0Info->fV0tpc->Update(vertex);
	
	  //
	  //
	  fRecV0Info->fT1.GetESDtrack()->GetExternalParameters(x,param);
	  fRecV0Info->fT1.GetESDtrack()->GetExternalCovariance(cov);
	  alpha = fRecV0Info->fT1.GetESDtrack()->GetAlpha();
	  new (&paramP) AliExternalTrackParam(x,alpha,param,cov);
	  //
	  fRecV0Info->fT2.GetESDtrack()->GetExternalParameters(x,param);
	  fRecV0Info->fT2.GetESDtrack()->GetExternalCovariance(cov);
	  alpha = fRecV0Info->fT2.GetESDtrack()->GetAlpha();
	  new (&paramM) AliExternalTrackParam(x,alpha,param,cov);
	  //
	  fRecV0Info->fV0its->SetParamN(paramM);
	  fRecV0Info->fV0its->SetParamP(paramP);
	  //	  fRecV0Info->fV0its.UpdatePID(pid1,pid2);
	  fRecV0Info->fV0its->Update(vertex);
	}
      }
      //
      // ????
      // 
      if (TMath::Abs(fGenV0Info->GetMinus().fPdg)==11 &&TMath::Abs(fGenV0Info->GetPlus().fPdg)==11){
	if (fRecV0Info->fDist2>10){
	  fRecV0Info->Update(vertex);
	}
	if (fRecV0Info->fDist2>10){
	  fRecV0Info->Update(vertex);
	}
      }
    }   
    //
    // take the V0 from reconstruction
 
    Int_t label =  TMath::Min(fGenV0Info->GetMinus().fLabel,fGenV0Info->GetPlus().fLabel);
    Int_t label2 = TMath::Max(fGenV0Info->GetMinus().fLabel,fGenV0Info->GetPlus().fLabel);    
    AliV0 *v0MI=0;
    AliV0 *v0MIOff=0;
    fRecV0Info->fRecStatus   =0;
    fRecV0Info->fMultiple    = fMultiRecV0[label];
    fRecV0Info->fV0MultipleOn=0;
    fRecV0Info->fV0MultipleOff=0;
    //
    if (fMultiRecV0[label]>0 || fMultiRecV0[label2]>0){

      //      for (Int_t j=0;j<TMath::Min(fMultiRecV0s[label],100);j++){
      for (Int_t j=TMath::Min(fMultiRecV0[label],Short_t(20))-1;j>=0;j--){
	Int_t index = fIndexRecV0[label*20+j];
	if (index<0) continue;
	AliV0 *v0MI2  = (AliV0*)fEvent->GetV0(index);
	// get track labels
	AliESDtrack * trackn = fEvent->GetTrack((v0MI2->GetNindex()));
	AliESDtrack * trackp = fEvent->GetTrack((v0MI2->GetPindex()));
	Int_t vlabeln = (trackn==0) ? -1 : trackn->GetLabel(); 
	Int_t vlabelp = (trackp==0) ? -1 : trackp->GetLabel(); 
	fRecV0Info->fLab[0]=TMath::Abs(vlabelp);
	fRecV0Info->fLab[1]=TMath::Abs(vlabeln); 
	//
	if (TMath::Abs(vlabeln)==label &&TMath::Abs(vlabelp)==label2) {
	  if (v0MI2->GetOnFlyStatus()) {
	    v0MI =v0MI2;
	    fRecV0Info->fV0MultipleOn++;
	  }else  {
	    v0MIOff = v0MI2;
	    fRecV0Info->fV0MultipleOff++;
	  }
	  fSignedV0[index]=1;
	}
	if (TMath::Abs(vlabelp)==label &&TMath::Abs(vlabeln)==label2) {
	  if (v0MI2->GetOnFlyStatus()){
	    v0MI =v0MI2;
	    fRecV0Info->fV0MultipleOn++;
	  }else  {
	    v0MIOff = v0MI2;
	    fRecV0Info->fV0MultipleOff++;
	  }
	  fSignedV0[index]=1;
	}
      }
    }
    if (v0MI){
      new (fRecV0Info->fV0rec) AliV0(*v0MI);
      fRecV0Info->fRecStatus=1;
    }
    if (v0MIOff){
      new (fRecV0Info->fV0recOff) AliV0(*v0MIOff);
      fRecV0Info->fRecStatus=1;
    }
    Int_t mpdg = fGenV0Info->GetMother().GetPdgCode();
    Float_t mass = ( pdgtable.GetParticle(mpdg)==0) ? 0 :pdgtable.GetParticle(mpdg)->Mass();
    fRecV0Info->UpdateKF(*esdvertex,
			 fGenV0Info->GetPlus().GetPdg(),
			 fGenV0Info->GetMinus().GetPdg(),
			 mass);
    fTreeCmpV0->Fill();
  }
  //
  // write fake v0s
  //
  Int_t nV0MIs = fEvent->GetNumberOfV0s();
  for (Int_t i=0;i<nV0MIs;i++){
    if (fSignedV0[i]==0){
      AliV0 *v0MI  = (AliV0*)fEvent->GetV0(i);
      if (!v0MI) continue;
      fRecV0Info->Reset();  //reset all variables
      //
      new (fRecV0Info->fV0rec) AliV0(*v0MI);
      fRecV0Info->fV0Status  =-10;
      fRecV0Info->fRecStatus =-2;
      //
      AliESDtrack * trackn = fEvent->GetTrack((v0MI->GetNindex()));
      AliESDtrack * trackp = fEvent->GetTrack((v0MI->GetPindex()));
      Int_t vlabeln = (trackn==0) ? -1 : trackn->GetLabel(); 
      Int_t vlabelp = (trackp==0) ? -1 : trackp->GetLabel(); 
      fRecV0Info->fLab[0]=TMath::Abs(vlabelp);
      fRecV0Info->fLab[1]=TMath::Abs(vlabeln);      
      if (TMath::Abs(fRecV0Info->fLab[0] - fRecV0Info->fLab[1])<2) continue;
      AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(vlabeln));
      AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(vlabelp));
      if (fRecInfo1 && fRecInfo2){
	fRecV0Info->fT1 = (*fRecInfo1);
	fRecV0Info->fT2 = (*fRecInfo2);
	fRecV0Info->fRecStatus =-1;
      }
      fRecV0Info->Update(vertex);
      fRecV0Info->UpdateKF(*esdvertex,211,211,0.49767);
      fTreeCmpV0->Fill();
    }
  }



  fTreeCmpV0->AutoSave();
  printf("Time spended in BuilV0Info Loop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of BuildV0Info Loop"<<endl;
  return 0;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


