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


///////////////////////////////////////////////////////////////////////////
/*

Origin: marian.ivanov@cern.ch

Generate complex MC information - used for Comparison later on
How to use it?


---Usage outside of the analysis framework

gSystem->Load("libANALYSIS.so")
gSystem->Load("libPWG1.so")
AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0)
t->Exec();

*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "TROOT.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

//ALIROOT includes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include "AliRun.h"
#include "AliStack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliTrackPointArray.h"

#endif
#include "AliMCInfo.h" 
#include "AliGenV0Info.h" 
#include "AliGenKinkInfo.h" 
#include "AliGenInfoMaker.h" 
//
// 

ClassImp(AliGenInfoMaker)






  
////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker(): 
  fGenTracksArray(0),          //clones array with filtered particles
  fGenKinkArray(0),            //clones array with filtered Kinks
  fGenV0Array(0),              //clones array with filtered V0s
  fDebug(0),                   //! debug flag  
  fEventNr(0),                 //! current event number
  fLabel(0),                   //! track label
  fNEvents(0),                 //! number of events to process
  fFirstEventNr(0),            //! first event to process
  fNParticles(0),              //! number of particles in TreeK
  fTreeGenTracks(0),           //! output tree with generated tracks
  fTreeKinks(0),               //!  output tree with Kinks
  fTreeV0(0),                  //!  output tree with V0
  fFileGenTracks(0),           //! output file with stored fTreeGenTracks
  fLoader(0),                  //! pointer to the run loader
  fTreeD(0),                   //! current tree with digits
  fTreeTR(0),                  //! current tree with TR
  fStack(0),                   //! current stack
  fGenInfo(0),                 //! array with pointers to gen info
  fNInfos(0),                  //! number of tracks with infos
  fParamTPC(0),                //! AliTPCParam
  fTPCPtCut(0.1),              //  TPC pt cut            
  fITSPtCut(0.1),              //  ITS pt cut
  fTRDPtCut(0.1),              //  TRD pt cut
  fTOFPtCut(0.1)               //  TOF pt cut
{    
  snprintf(fFnRes,1000,"%s","genTracks.root");
  for(Int_t i=0; i<3; i++) { fVPrim[i] = 0; };  
}




////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker(const char * fnGalice, const char* fnRes,
				 Int_t nEvents, Int_t firstEvent):
  fGenTracksArray(0),          //clones array with filtered particles
  fGenKinkArray(0),            //clones array with filtered Kinks
  fGenV0Array(0),              //clones array with filtered V0s
  fDebug(0),                   //! debug flag  
  fEventNr(0),                 //! current event number
  fLabel(0),                   //! track label
  fNEvents(0),                 //! number of events to process
  fFirstEventNr(0),            //! first event to process
  fNParticles(0),              //! number of particles in TreeK
  fTreeGenTracks(0),           //! output tree with generated tracks
  fTreeKinks(0),               //!  output tree with Kinks
  fTreeV0(0),                  //!  output tree with V0
  fFileGenTracks(0),           //! output file with stored fTreeGenTracks
  fLoader(0),                  //! pointer to the run loader
  fTreeD(0),                   //! current tree with digits
  fTreeTR(0),                  //! current tree with TR
  fStack(0),                   //! current stack
  fGenInfo(0),                 //! array with pointers to gen info
  fNInfos(0),                  //! number of tracks with infos
  fParamTPC(0),                //! AliTPCParam 
  fTPCPtCut(0.1),  
  fITSPtCut(0.1),  
  fTRDPtCut(0.1), 
  fTOFPtCut(0.1)
{
  //
  // 
  //
  for(Int_t i=0; i<3; i++) { fVPrim[i] = 0; };  
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  snprintf(fFnRes,1000,"%s",fnRes);
  //
  //
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    delete AliRunLoader::Instance();
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
}




////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::~AliGenInfoMaker()
{
  //
  // Destructor
  //
  
  if (fLoader){
    fLoader->UnloadgAlice();
    gAlice = 0;
    delete fLoader;
  }
}


AliMCInfo * AliGenInfoMaker::MakeInfo(UInt_t i)
{
  // 
  // Make info structure for given particle index
  //
  if (i<fNParticles) {
    if (fGenInfo[i]) return  fGenInfo[i];
    fGenInfo[i] = new AliMCInfo;  
    fNInfos++;
    return fGenInfo[i];
  }
  else 
    return 0;  
}


Int_t  AliGenInfoMaker::ProcessEvent(AliMCEventHandler* mcinfo){
  //
  // Process MC info from the task
  //
  if (!fParamTPC) {
    SetIO();
    fParamTPC = GetTPCParam();
  }
  fStack  = mcinfo->MCEvent()->Stack();
  fTreeTR = mcinfo->TreeTR();
  fTreeD  = 0;
  // array with preprocessedprocessed information
  fGenTracksArray = new TObjArray;
  fGenKinkArray   = new TObjArray;
  fGenV0Array     = new TObjArray;
  //
  ProcessEvent();
  fEventNr++;
  return 0;
}



Int_t AliGenInfoMaker::ProcessEvent(){
  //
  // Process Event 
  //
  fNParticles = fStack->GetNtrack();
  //
  fGenInfo = new AliMCInfo*[fNParticles];
  for (UInt_t i = 0; i<fNParticles; i++) {
    fGenInfo[i]=0; 
  }
  //
  cout<<"Start to process event "<<fEventNr<<endl;
  cout<<"\tfNParticles = "<<fNParticles<<endl;
  if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeTR"<<endl;
  if (TreeTRLoop()>0) return 1;
  //
  if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeD"<<endl;
  if (TreeDLoop()>0) return 1;
  //
  if (fDebug>2) cout<<"\n\n\n\tStart loop over TreeK"<<endl;
  if (TreeKLoop()>0) return 1;
  if (fDebug>2) cout<<"\tEnd loop over TreeK"<<endl;
  //
  //  if (BuildKinkInfo()>0) return 1;
  if (BuildV0Info()>0) return 1;
  if (fDebug>2) cout<<"\tEnd loop over TreeK"<<endl;
  //
  for (UInt_t i = 0; i<fNParticles; i++) {
    if (fGenInfo[i]) delete fGenInfo[i]; 
  }
  delete []fGenInfo;
  return 0;
}







Int_t  AliGenInfoMaker::SetIO()
{
  //
  // Set IO for given event 
  //
  CreateTreeGenTracks();
  if (!fTreeGenTracks) return 1;
 
  fParamTPC = GetTPCParam();
  //
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  fLoader->SetEventNumber(eventNr);
  //
  fLoader->LoadHeader();
  fLoader->LoadKinematics();  
  fStack = fLoader->Stack();
  //
  fLoader->LoadTrackRefs();
  fLoader->LoadHits();
  fTreeTR = fLoader->TreeTR();
  //
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->LoadDigits();
  fTreeD = tpcl->TreeD();  
  return 0;
}

Int_t AliGenInfoMaker::CloseIOEvent()
{
  //
  // Close IO for current event
  //
  fLoader->UnloadHeader();
  fLoader->UnloadKinematics();
  fLoader->UnloadTrackRefs();
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->UnloadDigits();
  return 0;
}

Int_t AliGenInfoMaker::CloseIO()
{
  fLoader->UnloadgAlice();
  return 0;
}



////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::Exec()  
{
  //
  // Make a comparision MC tree
  // Connect MC information -TPArticle - AliTrackRefernces ...
  //
  TStopwatch timer;
  timer.Start();
  Int_t status =SetIO();
  if (status>0) return status;
  //
  for (fEventNr = fFirstEventNr; fEventNr < fFirstEventNr+fNEvents;
       fEventNr++) {
    SetIO(fEventNr);
    ProcessEvent();
    CloseIOEvent();
  }
  //
  CloseIO();
  CloseOutputFile();
  cerr<<"Exec finished"<<endl;
  timer.Stop();
  timer.Print();
  return 0;
}






////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::CreateTreeGenTracks() 
{
  //
  //
  //
  fFileGenTracks = TFile::Open(fFnRes,"RECREATE");
  if (!fFileGenTracks) {
    cerr<<"Error in CreateTreeGenTracks: cannot open file "<<fFnRes<<endl;
    return;
  }
  fTreeGenTracks = new TTree("genTracksTree","genTracksTree");  
  AliMCInfo * info = new AliMCInfo;
  fTreeGenTracks->Branch("MC","AliMCInfo",&info);
  delete info; 
  //
  AliGenKinkInfo *kinkinfo = new AliGenKinkInfo;
  fTreeKinks = new TTree("genKinksTree","genKinksTree");  
  fTreeKinks->Branch("MC","AliGenKinkInfo",&kinkinfo);
  delete kinkinfo;
  //
  AliGenV0Info *v0info = new AliGenV0Info;
  fTreeV0 = new TTree("genV0Tree","genV0Tree");  
  fTreeV0->Branch("MC","AliGenV0Info",&v0info);
  delete v0info;
  //
  //
  //
  fTreeGenTracks->AutoSave();
  fTreeKinks->AutoSave();
  fTreeV0->AutoSave();
}



////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::CloseOutputFile() 
{
  //
  // Close Output files
  //
  if (!fFileGenTracks) {
    cerr<<"File "<<fFnRes<<" not found as an open file."<<endl;
    return;
  }
  fFileGenTracks->cd();
  fTreeGenTracks->Write();  
  delete fTreeGenTracks;
  fTreeKinks->Write();
  delete fTreeKinks;
  fTreeV0->Write();
  delete fTreeV0;

  fFileGenTracks->Close();
  delete fFileGenTracks;
  return;
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeKLoop()
{
//
// open the file with treeK
// loop over all entries there and save information about some tracks
//

  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  Int_t  ipdg = 0;
  //TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  

  TParticle *particle = stack->ParticleFromTreeK(0);
  fVPrim[0] = particle->Vx();
  fVPrim[1] = particle->Vy();
  fVPrim[2] = particle->Vz();
  Int_t accepted =0;
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    //////////////////////////////////////////////////////////////////////
    info->fLabel = iParticle;
    //
    info->fParticle = *(stack->Particle(iParticle));
    info->fVDist[0] = info->fParticle.Vx()-fVPrim[0];
    info->fVDist[1] = info->fParticle.Vy()-fVPrim[1];
    info->fVDist[2] = info->fParticle.Vz()-fVPrim[2]; 
    info->fVDist[3] = TMath::Sqrt(info->fVDist[0]*info->fVDist[0]+
				  info->fVDist[1]*info->fVDist[1]+info->fVDist[2]*info->fVDist[2]);
    //
    //
    ipdg = info->fParticle.GetPdgCode();
    info->fPdg = ipdg;
    //ppdg = pdg.GetParticle(ipdg);   	   
    info->fEventNr = fEventNr;
    info->Update();
    if (fGenTracksArray){
      fGenTracksArray->AddLast(info->Clone());
    }
    accepted++;
  }
  //
  // write the results to the tree - if specified
  //
  TBranch * br = fTreeGenTracks->GetBranch("MC");
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    //////////////////////////////////////////////////////////////////////    
    br->SetAddress(&info);    
    fTreeGenTracks->Fill();    
  }
  fTreeGenTracks->AutoSave();
  //
  //
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////
Int_t  AliGenInfoMaker::BuildKinkInfo()
{
  //
  // Build tree with Kink Information
  //
  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  //  Int_t  ipdg = 0;
  //TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeKinks->GetBranch("MC");
  //
  AliGenKinkInfo * kinkinfo = new AliGenKinkInfo;
  //
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {
    // load only particles with TR
    AliMCInfo * info = GetInfo(iParticle);
    if (!info) continue;
    if (info->fCharge==0) continue;  
    if (info->fTRdecay.R()>500)  continue;  //R cut - decay outside barrel
    TParticle & particle = info->fParticle;
    Int_t first = particle.GetDaughter(0);
    if (first ==0) continue;

    Int_t last = particle.GetDaughter(1);
    if (last ==0) last=first;
    AliMCInfo * info2 =  0;
    AliMCInfo * dinfo =  0;
    
    for (Int_t id2=first;id2<=last;id2++){
      info2 = GetInfo(id2);
      if (!info2) continue;
      TParticle &p2 = info2->fParticle;
      Double_t r = TMath::Sqrt(p2.Vx()*p2.Vx()+p2.Vy()*p2.Vy());
      if ( TMath::Abs(info->fTRdecay.R()-r)>0.01) continue;
      if (!(p2.GetPDG())) continue;
      dinfo =info2;
     
      kinkinfo->SetInfoMother(*info);
      kinkinfo->SetInfoDaughter(*dinfo);
      if (kinkinfo->GetMinus().GetParticle().GetPDG()==0 ||  kinkinfo->GetPlus().GetParticle().GetPDG()==0) continue;
      kinkinfo->Update();
      br->SetAddress(&kinkinfo);    
      fTreeKinks->Fill();
    }
  }
  fTreeKinks->AutoSave();
  if (fDebug > 2) cerr<<"end of Kink Loop"<<endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////
Int_t  AliGenInfoMaker::BuildV0Info()
{
  //
  // Build tree with V0 Information
  //
  AliStack * stack = fStack;
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}
  
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }  
  //
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeV0->GetBranch("MC");
  //
  AliGenV0Info * v0info = new AliGenV0Info;
  //
  //
  for (UInt_t iParticle = 0; iParticle < fNParticles; iParticle++) {    
    TParticle * mParticle = fStack->Particle(iParticle);
    if (!mParticle) continue;
    if (mParticle->GetPDG()==0) continue;
    if (mParticle->GetPDG()->Charge()!=0) continue;  //only neutral particle
    //
    Int_t first = mParticle->GetDaughter(0);
    Int_t last  = mParticle->GetDaughter(1);
    if (first-last==0) continue;
    //
    AliMCInfo * info0 = GetInfo(first);
    AliMCInfo * info1 = GetInfo(first+1);
    if (info0 && info1){
      if (info0->GetPdg()==0) continue;
      if (info1->GetPdg()==0) continue;
      v0info->SetInfoP(*info0);
      v0info->SetInfoM(*info1);
      v0info->SetMother(*mParticle);
      v0info->Update(fVPrim);
      if (fGenV0Array) fGenV0Array->AddLast(v0info);
      br->SetAddress(&v0info);    
      fTreeV0->Fill();
    }
  }
  fTreeV0->AutoSave();
  if (fDebug > 2) cerr<<"end of V0 Loop"<<endl;
  return 0;
}





////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeDLoop()
{
  //
  // open the file with treeD
  // loop over all entries there and save information about some tracks
  //
  if (!fTreeD) return 0;
  Int_t nInnerSector = fParamTPC->GetNInnerSector();
  Int_t rowShift = 0;
  Int_t zero=fParamTPC->GetZeroSup()+6;  
  //
  //
  AliSimDigits digitsAddress, *digits=&digitsAddress;
  fTreeD->GetBranch("Segment")->SetAddress(&digits);
  
  Int_t sectorsByRows=(Int_t)fTreeD->GetEntries();
  if (fDebug > 1) cout<<"\tsectorsByRows = "<<sectorsByRows<<endl;
  for (Int_t i=0; i<sectorsByRows; i++) {
    if (!fTreeD->GetEvent(i)) continue;
    Int_t sec,row;
    fParamTPC->AdjustSectorRow(digits->GetID(),sec,row);
    if (fDebug > 1) cout<<sec<<' '<<row<<"                          \r";
    // here I expect that upper sectors follow lower sectors
    if (sec > nInnerSector) rowShift = fParamTPC->GetNRowLow();
    //
    digits->ExpandTrackBuffer();
    digits->First();        
    do {
      Int_t iRow=digits->CurrentRow();
      Int_t iColumn=digits->CurrentColumn();
      Short_t digitValue = digits->CurrentDigit();
      if (digitValue >= zero) {
	Int_t label;
	for (Int_t j = 0; j<3; j++) {
	  //	  label = digits->GetTrackID(iRow,iColumn,j);
	  label = digits->GetTrackIDFast(iRow,iColumn,j)-2; 	  
	  if (label >= (Int_t)fNParticles) { //don't label from bakground event
	    continue;
	  }
	  if (label >= 0 && label <= (Int_t)fNParticles) {
	    if (fDebug > 6 ) {
	      cout<<"Inner loop: sector, iRow, iColumn, label, value, row "
		  <<sec<<" "
		  <<iRow<<" "<<iColumn<<" "<<label<<" "<<digitValue
		  <<" "<<row<<endl;
	    }	
	    AliMCInfo * info = GetInfo(label);
	    if (info){
	      info->fTPCRow.SetRow(row+rowShift);
	    }
	  }
	}
      }
    } while (digits->Next());
  }
  
  if (fDebug > 2) cerr<<"end of TreeDLoop"<<endl;  
  return 0;
}




////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeTRLoop()
{
  //
  // loop over TrackReferences and store the first one for each track
  //  
  TTree * treeTR = fTreeTR;
  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  if (fDebug > 1) cout<<"There are "<<nPrimaries<<" entries in TreeTR"<<endl;
  //
  //
  //
  TClonesArray* runArrayTR = new TClonesArray("AliTrackReference");
  //
  if (treeTR->GetBranch("TrackReferences")) treeTR->GetBranch("TrackReferences")->SetAddress(&runArrayTR);
  //
  //
  //
  for (Int_t ipart = 0; ipart<fStack->GetNtrack(); ipart++) {
    TParticle * part = fStack->Particle(ipart);
    if (!part) continue;
    if (part->GetPDG()==0) continue;
    if (part->GetPDG()->Charge()==0) continue;
    
    if (part->Pt()<fITSPtCut) continue;
    if (part->R()>250.) continue;
    if (TMath::Abs(part->Vz())>250.) continue;
    if (TMath::Abs(part->Pz()/part->Pt())>2.5) continue;     
    //AliMCInfo * info = ;
    MakeInfo(ipart); 
  }
  //
  //

  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    treeTR->GetEntry(iPrimPart);
    //
    // gettrack references
    //
    for (Int_t iTrackRef = 0; iTrackRef < runArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)runArrayTR->At(iTrackRef);      
      //
      Int_t label = trackRef->GetTrack();
      AliMCInfo * cinfo = GetInfo(label);
      if (cinfo) cinfo->CalcTPCrows(runArrayTR);
      //
      // TPC
      //
      if (trackRef->DetectorId()== AliTrackReference::kTPC){
	//
	Int_t label1 = trackRef->GetTrack();      
	AliMCInfo * info = GetInfo(label1);
	if (!info && trackRef->Pt()<fTPCPtCut) continue;
	if (!info) {
	  info = MakeInfo(label1);
	}
	if (!info) continue;
	info->fPrimPart =  iPrimPart;
	TClonesArray & arr = *(info->fTPCReferences);
	new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
      }
      //
      // ITS
      //
      if (trackRef->DetectorId()== AliTrackReference::kITS){
	//
	Int_t label1 = trackRef->GetTrack();      
	AliMCInfo * info = GetInfo(label1);
	if (!info && trackRef->Pt()<fITSPtCut) continue;
	if (!info) info = MakeInfo(label1);
	if (!info) continue;
	info->fPrimPart =  iPrimPart;
	TClonesArray & arr = *(info->fITSReferences);
	new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
      }
      //
      // TRD
      //
      if (trackRef->DetectorId()== AliTrackReference::kTRD){
	//
	Int_t label1 = trackRef->GetTrack();      
	AliMCInfo * info = GetInfo(label1);
	if (!info&&trackRef->P()<fTRDPtCut) continue;
	if (!info) info = MakeInfo(label1);
	if (!info) continue;
	info->fPrimPart =  iPrimPart;
	TClonesArray & arr = *(info->fTRDReferences);
	new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
      }
      //
      // TOF
      //
      if (trackRef->DetectorId()== AliTrackReference::kTOF){
	//
	Int_t label1 = trackRef->GetTrack();      
	AliMCInfo * info = GetInfo(label1);
	if (!info&&trackRef->P()<fTOFPtCut) continue;
	if (!info) info = MakeInfo(label1);
	if (!info) continue;
	info->fPrimPart =  iPrimPart;
	TClonesArray & arr = *(info->fTOFReferences);
	new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
      }
      //
      // decay case
      //
      if (trackRef->DetectorId()== AliTrackReference::kDisappeared){
	//
	AliMCInfo * info = GetInfo(label);
	if (!info) continue;
	info->fTRdecay = *trackRef;      	
      }
    }
  }
  //
  //
  return 0;
}

////////////////////////////////////////////////////////////////////////
Float_t AliGenInfoMaker::TR2LocalX(AliTrackReference *trackRef,
				    AliTPCParam *paramTPC) const {

  //
  // Transformation from Gloabal to local tacking system
  //
  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return x[0];
}
////////////////////////////////////////////////////////////////////////



AliTPCParam * AliGenInfoMaker::GetTPCParam(){
  //
  // create default TPC parameters
  //
  AliTPCParamSR * par = new AliTPCParamSR;
  par->Update();
  return par;
}


