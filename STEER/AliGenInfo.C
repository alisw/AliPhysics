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

Macro to generate comples MC information - used for Comparison later on
How to use it?

.L $ALICE_ROOT/STEER/AliGenInfo.C+
AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",1,0)
t->Exec();

*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

//ALIROOT includes
#include "AliRun.h"
#include "AliStack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliMagF.h"
#endif
#include "AliGenInfo.h" 
//
//

AliTPCParam * GetTPCParam(){
  AliTPCParamSR * par = new AliTPCParamSR;
  par->Update();
  return par;
}


//_____________________________________________________________________________
Float_t TPCBetheBloch(Float_t bg)
{
  //
  // Bethe-Bloch energy loss formula
  //
  const Double_t kp1=0.76176e-1;
  const Double_t kp2=10.632;
  const Double_t kp3=0.13279e-4;
  const Double_t kp4=1.8631;
  const Double_t kp5=1.9479;

  Double_t dbg = (Double_t) bg;

  Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);

  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1./dbg,kp5);

  bb=TMath::Log(kp3+bb);
  
  return ((Float_t)((kp2-aa-bb)*kp1/aa));
}





////////////////////////////////////////////////////////////////////////
AliMCInfo::AliMCInfo()
{
  fTPCReferences  = new TClonesArray("AliTrackReference",10);
  fITSReferences  = new TClonesArray("AliTrackReference",10);
  fTRDReferences  = new TClonesArray("AliTrackReference",10);
  fTRdecay.SetTrack(-1);
}

AliMCInfo::~AliMCInfo()
{
  if (fTPCReferences) {
    delete fTPCReferences;
  }
  if (fITSReferences){
    delete fITSReferences;
  }
  if (fTRDReferences){    
    delete fTRDReferences;  
  }
}



void AliMCInfo::Update()
{
  //
  //
  fMCtracks =1;
  if (!fTPCReferences) {
    fNTPCRef =0;
    return;
  }
  Float_t direction=1;
  //Float_t rlast=0;
  fNTPCRef = fTPCReferences->GetEntriesFast();
  fNITSRef = fITSReferences->GetEntriesFast();

  for (Int_t iref =0;iref<fTPCReferences->GetEntriesFast();iref++){
    AliTrackReference * ref = (AliTrackReference *) fTPCReferences->At(iref);
    //Float_t r = (ref->X()*ref->X()+ref->Y()*ref->Y());
    Float_t newdirection = ref->X()*ref->Px()+ref->Y()*ref->Py(); //inside or outside
    if (iref==0) direction = newdirection;
    if ( newdirection*direction<0){
      //changed direction
      direction = newdirection;
      fMCtracks+=1;
    }
    //rlast=r;			    
  }
  //
  // decay info
  fTPCdecay=kFALSE;
  if (fTRdecay.GetTrack()>0){
    fDecayCoord[0] = fTRdecay.X();
    fDecayCoord[1] = fTRdecay.Y();
    fDecayCoord[2] = fTRdecay.Z();
    if ( (fTRdecay.R()<250)&&(fTRdecay.R()>85) && (TMath::Abs(fTRdecay.Z())<250) ){
      fTPCdecay=kTRUE;     
    }
    else{
      fDecayCoord[0] = 0;
      fDecayCoord[1] = 0;
      fDecayCoord[2] = 0;
    }
  }
  //
  //
  //digits information update
  fRowsWithDigits    = fTPCRow.RowsOn();    
  fRowsWithDigitsInn = fTPCRow.RowsOn(63); // 63 = number of inner rows
  fRowsTrackLength   = fTPCRow.Last() - fTPCRow.First();
  //
  //
  // calculate primary ionization per cm  
  if (fParticle.GetPDG()){
    fMass = fParticle.GetMass();  
    fCharge = fParticle.GetPDG()->Charge();
    if (fTPCReferences->GetEntriesFast()>0){
      fTrackRef = *((AliTrackReference*)fTPCReferences->At(0));
    }
    if (fMass>0){
      Float_t p = TMath::Sqrt(fTrackRef.Px()*fTrackRef.Px()+
			      fTrackRef.Py()*fTrackRef.Py()+
			      fTrackRef.Pz()*fTrackRef.Pz());    
      if (p>0.001){
	Float_t betagama = p /fMass;
	fPrim = TPCBetheBloch(betagama);
      }else fPrim=0;
    }
  }else{
    fMass =0;
    fPrim =0;
  }  
}

/////////////////////////////////////////////////////////////////////////////////
void AliGenV0Info::Update()
{
  fMCPd[0] = fMCd.fParticle.Px();
  fMCPd[1] = fMCd.fParticle.Py();
  fMCPd[2] = fMCd.fParticle.Pz();
  fMCPd[3] = fMCd.fParticle.P();
  fMCX[0]  = fMCd.fParticle.Vx();
  fMCX[1]  = fMCd.fParticle.Vy();
  fMCX[2]  = fMCd.fParticle.Vz();
  fMCR       = TMath::Sqrt( fMCX[0]*fMCX[0]+fMCX[1]*fMCX[1]);
  fPdg[0]    = fMCd.fParticle.GetPdgCode();
  fPdg[1]    = fMCm.fParticle.GetPdgCode();
  //
  fLab[0]    = fMCd.fParticle.GetUniqueID();
  fLab[1]    = fMCm.fParticle.GetUniqueID();

}




  
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// End of implementation of the class AliMCInfo
//
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
digitRow::digitRow()
{
  Reset();
}
////////////////////////////////////////////////////////////////////////
digitRow & digitRow::operator=(const digitRow &digOld)
{
  for (Int_t i = 0; i<kgRowBytes; i++) fDig[i] = digOld.fDig[i];
  return (*this);
}
////////////////////////////////////////////////////////////////////////
void digitRow::SetRow(Int_t row) 
{
  if (row >= 8*kgRowBytes) {
    cerr<<"digitRow::SetRow: index "<<row<<" out of bounds."<<endl;
    return;
  }
  Int_t iC = row/8;
  Int_t iB = row%8;
  SETBIT(fDig[iC],iB);
}

////////////////////////////////////////////////////////////////////////
Bool_t digitRow::TestRow(Int_t row)
{
//
// return kTRUE if row is on
//
  Int_t iC = row/8;
  Int_t iB = row%8;
  return TESTBIT(fDig[iC],iB);
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::RowsOn(Int_t upto)
{
//
// returns number of rows with a digit  
// count only rows less equal row number upto
//
  Int_t total = 0;
  for (Int_t i = 0; i<kgRowBytes; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (i*8+j > upto) return total;
      if (TESTBIT(fDig[i],j))  total++;
    }
  }
  return total;
}
////////////////////////////////////////////////////////////////////////
void digitRow::Reset()
{
//
// resets all rows to zero
//
  for (Int_t i = 0; i<kgRowBytes; i++) {
    fDig[i] <<= 8;
  }
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::Last()
{
//
// returns the last row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = kgRowBytes-1; i>=0; i--) {
    for (Int_t j = 7; j >= 0; j--) {
      if TESTBIT(fDig[i],j) return i*8+j;
    }
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::First()
{
//
// returns the first row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = 0; i<kgRowBytes; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (TESTBIT(fDig[i],j)) return i*8+j;
    }
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////
//
// end of implementation of a class digitRow
//
////////////////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::AliGenInfoMaker(const char * fnGalice, const char* fnRes,
				   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  //   fFnRes = fnRes;
  sprintf(fFnRes,"%s",fnRes);
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    delete gAlice->GetRunLoader();
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
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf);
}


AliMCInfo * AliGenInfoMaker::MakeInfo(UInt_t i)
{
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

////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeGenTracks = 0;
  fFileGenTracks = 0;
  fGenInfo = 0;
  fNInfos  = 0;
  //
  //
  fDebug = 0;
  fVPrim[0] = -1000.;
  fVPrim[1] = -1000.;
  fVPrim[2] = -1000.;
  fParamTPC = 0;
}
////////////////////////////////////////////////////////////////////////
AliGenInfoMaker::~AliGenInfoMaker()
{
  
  if (fLoader){
    fLoader->UnloadgAlice();
    gAlice = 0;
    delete fLoader;
  }
}

Int_t  AliGenInfoMaker::SetIO()
{
  //
  // 
  CreateTreeGenTracks();
  if (!fTreeGenTracks) return 1;
  //  AliTracker::SetFieldFactor(); 
 
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
  fTreeTR = fLoader->TreeTR();
  //
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->LoadDigits();
  fTreeD = tpcl->TreeD();  
  return 0;
}

Int_t AliGenInfoMaker::CloseIOEvent()
{
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
Int_t AliGenInfoMaker::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::Exec()  
{
  TStopwatch timer;
  timer.Start();
  Int_t status =SetIO();
  if (status>0) return status;
  //

  for (fEventNr = fFirstEventNr; fEventNr < fFirstEventNr+fNEvents;
       fEventNr++) {
    SetIO(fEventNr);
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
    for (UInt_t i = 0; i<fNParticles; i++) {
      if (fGenInfo[i]) delete fGenInfo[i]; 
    }
    delete []fGenInfo;
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
  fFileGenTracks = TFile::Open(fFnRes,"RECREATE");
  if (!fFileGenTracks) {
    cerr<<"Error in CreateTreeGenTracks: cannot open file "<<fFnRes<<endl;
    return;
  }
  fTreeGenTracks = new TTree("genTracksTree","genTracksTree");  
  AliMCInfo * info = new AliMCInfo;
  //
  fTreeGenTracks->Branch("MC","AliMCInfo",&info);
  delete info; 
  fTreeGenTracks->AutoSave();
}
////////////////////////////////////////////////////////////////////////
void AliGenInfoMaker::CloseOutputFile() 
{
  if (!fFileGenTracks) {
    cerr<<"File "<<fFnRes<<" not found as an open file."<<endl;
    return;
  }
  fFileGenTracks->cd();
  fTreeGenTracks->Write();  
  delete fTreeGenTracks;
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
  TParticlePDG *ppdg = 0;
  // not all generators give primary vertex position. Take the vertex
  // of the particle 0 as primary vertex.
  TDatabasePDG  pdg; //get pdg table  
  //thank you very much root for this
  TBranch * br = fTreeGenTracks->GetBranch("MC");
  TParticle *particle = stack->ParticleFromTreeK(0);
  fVPrim[0] = particle->Vx();
  fVPrim[1] = particle->Vy();
  fVPrim[2] = particle->Vz();
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
    ppdg = pdg.GetParticle(ipdg);   	   
    info->fEventNr = fEventNr;
    info->Update();
    //////////////////////////////////////////////////////////////////////    
    br->SetAddress(&info);    
    fTreeGenTracks->Fill();
  }
  fTreeGenTracks->AutoSave();
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;
  return 0;
}




////////////////////////////////////////////////////////////////////////
Int_t AliGenInfoMaker::TreeDLoop()
{
  //
  // open the file with treeD
  // loop over all entries there and save information about some tracks
  //
  
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
  //track references for TPC
  TClonesArray* TPCArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* ITSArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* TRDArrayTR = new TClonesArray("AliTrackReference");
  TClonesArray* RunArrayTR = new TClonesArray("AliTrackReference");
  //
  if (treeTR->GetBranch("TPC"))    treeTR->GetBranch("TPC")->SetAddress(&TPCArrayTR);
  if (treeTR->GetBranch("ITS"))    treeTR->GetBranch("ITS")->SetAddress(&ITSArrayTR);
  if (treeTR->GetBranch("TRD"))    treeTR->GetBranch("TRD")->SetAddress(&TRDArrayTR);
  if (treeTR->GetBranch("AliRun")) treeTR->GetBranch("AliRun")->SetAddress(&RunArrayTR);
  //
  //
  //
  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    treeTR->GetEntry(iPrimPart);
    //
    // Loop over TPC references
    //
    for (Int_t iTrackRef = 0; iTrackRef < TPCArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TPCArrayTR->At(iTrackRef);            
      //
      if (trackRef->Pt()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      TClonesArray & arr = *(info->fTPCReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // Loop over ITS references
    //
    for (Int_t iTrackRef = 0; iTrackRef < ITSArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)ITSArrayTR->At(iTrackRef);            
      //
      if (trackRef->Pt()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if ( (!info) && trackRef->Pt()<fgITSPtCut) continue;
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      TClonesArray & arr = *(info->fITSReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // Loop over TRD references
    //
    for (Int_t iTrackRef = 0; iTrackRef < TRDArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TRDArrayTR->At(iTrackRef);            
      //
      if (trackRef->Pt()<fgTPCPtCut) continue;
      Int_t label = trackRef->GetTrack();      
      AliMCInfo * info = GetInfo(label);
      if ( (!info) && trackRef->Pt()<fgTRDPtCut) continue;
      if (!info) info = MakeInfo(label);
      if (!info) continue;
      TClonesArray & arr = *(info->fTRDReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // get dacay position
    for (Int_t iTrackRef = 0; iTrackRef < RunArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)RunArrayTR->At(iTrackRef);      
      //
      Int_t label = trackRef->GetTrack();
      AliMCInfo * info = GetInfo(label);
      if (!info) continue;
      if (!trackRef->TestBit(BIT(2))) continue;  //if not decay
      //      if (TMath::Abs(trackRef.X());
      info->fTRdecay = *trackRef;      
    }
  }
  //
  TPCArrayTR->Delete();
  delete  TPCArrayTR;
  TRDArrayTR->Delete();
  delete  TRDArrayTR;
  ITSArrayTR->Delete();
  delete  ITSArrayTR;
  RunArrayTR->Delete();
  delete  RunArrayTR;
  //
  return 0;
}

////////////////////////////////////////////////////////////////////////
Float_t AliGenInfoMaker::TR2LocalX(AliTrackReference *trackRef,
				    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return x[0];
}
////////////////////////////////////////////////////////////////////////



TH1F * AliComparisonDraw::DrawXY(const char * chx, const char *chy, const char* selection, 
		const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy)
{
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  Int_t nBinsRes = 30;
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, minx, maxx, nBinsRes, miny, maxy);
  char cut[1000];
  sprintf(cut,"%s&&%s",selection,quality);
  char expression[1000];
  sprintf(expression,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}

TH1F * AliComparisonDraw::DrawLogXY(const char * chx, const char *chy, const char* selection, 
		const char * quality, Int_t nbins, Float_t minx, Float_t maxx, Float_t miny, Float_t maxy)
{
  //
  Double_t* bins = CreateLogBins(nbins, minx, maxx);
  Int_t nBinsRes = 30;
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nbins, bins, nBinsRes, miny, maxy);
  char cut[1000];
  sprintf(cut,"%s&&%s",selection,quality);
  char expression[1000];
  sprintf(expression,"%s:%s>>hRes2",chy,chx);
  fTree->Draw(expression, cut, "groff");
  TH1F* hMean=0;  
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, chx, chy);
  //
  delete hRes2;
  delete[] bins;
  fRes  = hRes;
  fMean = hMean;
  return hRes;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliComparisonDraw::Eff(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, min, max);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, min, max);
  char inputGen[1000];  
  sprintf(inputGen,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && %s", selection, quality);
  char inputRec[1000];  
  sprintf(inputRec,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  return hEff;
}



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
TH1F * AliComparisonDraw::EffLog(const char *variable, const char* selection, const char * quality, 
			      Int_t nbins, Float_t min, Float_t max)
{
  //
  //
  Double_t* bins = CreateLogBins(nbins, min, max);
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nbins, bins);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nbins, bins);
  char inputGen[1000];  
  sprintf(inputGen,"%s>>hGen", variable);
  fTree->Draw(inputGen, selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && %s", selection, quality);
  char inputRec[1000];  
  sprintf(inputRec,"%s>>hRec", variable);
  fTree->Draw(inputRec, selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, variable, "#epsilon [%]");
  fRes = hEff;
  delete hRec;
  delete hGen;
  delete[] bins;
  return hEff;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

Double_t* AliComparisonDraw::CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax)
{
  Double_t* bins = new Double_t[nBins+1];
  bins[0] = xMin;
  Double_t factor = pow(xMax/xMin, 1./nBins);
  for (Int_t i = 1; i <= nBins; i++)
    bins[i] = factor * bins[i-1];
  return bins;
}




void AliComparisonDraw::AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle)
{
  //
  histo->GetXaxis()->SetTitle(xAxisTitle);
  histo->GetYaxis()->SetTitle(yAxisTitle);
}


TH1F* AliComparisonDraw::CreateEffHisto(TH1F* hGen, TH1F* hRec)
{
  //
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  //
  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGen = hGen->GetBinContent(iBin);
    Double_t nRec = hRec->GetBinContent(iBin);
    if (nGen > 0) {
      Double_t eff = nRec/nGen;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt((eff*(1.-eff)+0.01) / nGen);      
      //      if (error == 0) error = sqrt(0.1/nGen);
      //
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);
    } else {
      hEff->SetBinContent(iBin, 100. * 0.5);
      hEff->SetBinError(iBin, 100. * 0.5);
    }
  }
  return hEff;
}



TH1F* AliComparisonDraw::CreateResHisto(TH2F* hRes2, TH1F **phMean,  Bool_t drawBinFits, 
		     Bool_t overflowBinFits)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  //TKFitGaus* fitFunc = new TKFitGaus("resFunc");
  //   TF1 * fitFunc = new TF1("G","[3]+[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  Int_t dBin = ((overflowBinFits) ? 1 : 0);
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin, bin);
    //    
    if (hBin->GetEntries() > 5) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	sprintf(name, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	sprintf(name, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	sprintf(name, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}



