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
//  Comparison macro for TPC                                                 //
//  responsible: 
//  marian.ivanov@cern.ch                                                    //
//
//
// usage:

//
//  TO CREATE OF TREE WITH MC INFO

//  .L AliTPCComparisonMI.C+
//  TPCFindGenTracks *t = new TPCFindGenTracks("galice.root","genTracks.root",1,0)
//  t->Exec()
//  
//  TO CREATE COMPARISON TREE
//  TPCCmpTr *t2 = new TPCCmpTr("genTracks.root","cmpTracks.root","galice.root",1,0);
//  t2->Exec()
//
//  EXAMPLE OF COMPARISON VISUALIZATION SESSION 
//
// .L AliTPCComparisonMI.C+
// TCut cprim("cprim","MC.fVDist[3]<1");
// TCut cprim("cprim","MC.fVDist[3]<1");
// TCut crec("crec","fReconstructed==1");
// TCut cteta1("cteta1","abs(MC.fTrackRef.Theta()/3.1415-0.5)<0.25");
// TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
// TCut csens("csens","abs(sqrt(fVDist[0]**2+fVDist[1]**2)-170)<50");
// TCut cmuon("cmuon","abs(MC.fParticle.fPdgCode==-13)");
//
//
// AliTPCComparisonDraw comp;
// comp.SetIO();
// (comp.EffVsPt("MC.fRowsWithDigits>120"+cteta1+cprim,"1"))->Draw()
// (comp.EffVsRM("MC.fRowsWithDigits>20"+cteta1+cprim,"1"))->Draw()
// (comp.EffVsRS("MC.fRowsWithDigits>20"+cteta1+cpos1,"1"))->Draw()
//
// (comp.ResPtvsPt("MC.fRowsWithDigits>20"+cteta1+cpos1,"1",0.15,2.))->Draw()
// (comp.MeanPtvsPt("MC.fRowsWithDigits>20"+cteta1+cpos1,"1",0.15,2.))->Draw()
// (comp.ResdEdxvsN("RC.fReconstructed==1&&MC.fPrim<1.5&&abs(MC.fParticle.fPdgCode)==211&&MC.fParticle.P()>0.35"+cteta1+cpos1+cprim,"1",100,160,4))->Draw() 


///////////////////////////////////////////////////////////////////////////////




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
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TText.h"
#include "Getline.h"
#include "TStyle.h"

//ALIROOT includes
#include "AliRun.h"
#include "AliStack.h"
#include "AliTPCtrack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "AliRun.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliComplexCluster.h"
#include "AliTPCComparisonMI.h"
#include "AliMagF.h"
#endif







 
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
AliTPCGenInfo::AliTPCGenInfo()
{
  fReferences  = 0;
  fTRdecay.SetTrack(-1);
}

AliTPCGenInfo::~AliTPCGenInfo()
{
  if (fReferences) delete fReferences;
  fReferences =0;
  
}

void AliTPCGenInfo::Update()
{
  //
  //
  fMCtracks =1;
  if (!fReferences) return;
  Float_t direction=1;
  //Float_t rlast=0;
  for (Int_t iref =0;iref<fReferences->GetEntriesFast();iref++){
    AliTrackReference * ref = (AliTrackReference *) fReferences->At(iref);
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
}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// End of implementation of the class AliTPCGenInfo
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
TPCFindGenTracks::TPCFindGenTracks()
{
  fMCInfo = new AliTPCGenInfo;
  fMCInfo->fReferences = new TClonesArray("AliTrackReference");

  Reset();
}

////////////////////////////////////////////////////////////////////////
TPCFindGenTracks::TPCFindGenTracks(const char * fnGalice, const char* fnRes,
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

////////////////////////////////////////////////////////////////////////
void TPCFindGenTracks::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeGenTracks = 0;
  fFileGenTracks = 0;
  fContainerDigitRow = 0;
  //
  fReferences = 0;
  fReferenceIndex0 = 0;
  fReferenceIndex1 = 0;
  fDecayRef   = 0;
  //
  fDebug = 0;
  fVPrim[0] = -1000.;
  fVPrim[1] = -1000.;
  fVPrim[2] = -1000.;
  fParamTPC = 0;
}
////////////////////////////////////////////////////////////////////////
TPCFindGenTracks::~TPCFindGenTracks()
{
  
  if (fLoader){
    fLoader->UnloadgAlice();
    gAlice = 0;
    delete fLoader;
  }
}


Int_t  TPCFindGenTracks::SetIO()
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
Int_t TPCFindGenTracks::SetIO(Int_t eventNr)
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

Int_t TPCFindGenTracks::CloseIOEvent()
{
  fLoader->UnloadHeader();
  fLoader->UnloadKinematics();
  fLoader->UnloadTrackRefs();
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->UnloadDigits();
  return 0;
}

Int_t TPCFindGenTracks::CloseIO()
{
  fLoader->UnloadgAlice();
  return 0;
}



////////////////////////////////////////////////////////////////////////
Int_t TPCFindGenTracks::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t TPCFindGenTracks::Exec()  
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
    fContainerDigitRow = new digitRow[fNParticles];
    //
    fReferences      = new AliTrackReference[fgMaxTR];
    fReferenceIndex0 = new Int_t[fNParticles];
    fReferenceIndex1 = new Int_t[fNParticles];
    fDecayRef        = new AliTrackReference[fNParticles];

    for (Int_t i = 0; i<fNParticles; i++) {
      fReferenceIndex0[i] = -1;
      fReferenceIndex1[i] = -1;
    }
    
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\tStart loop over TreeD"<<endl;
    if (TreeDLoop()>0) return 1;
    if (fDebug>2) cout<<"\tStart loop over TreeTR"<<endl;
    if (TreeTRLoop()>0) return 1;
    if (fDebug>2) cout<<"\tStart loop over TreeK"<<endl;
    if (TreeKLoop()>0) return 1;
    if (fDebug>2) cout<<"\tEnd loop over TreeK"<<endl;

    delete [] fContainerDigitRow;
    delete [] fReferences;
    delete [] fReferenceIndex0;
    delete [] fReferenceIndex1;    
    delete [] fDecayRef;
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
void TPCFindGenTracks::CreateTreeGenTracks() 
{
  fFileGenTracks = TFile::Open(fFnRes,"RECREATE");
  if (!fFileGenTracks) {
    cerr<<"Error in CreateTreeGenTracks: cannot open file "<<fFnRes<<endl;
    return;
  }
  fTreeGenTracks = new TTree("genTracksTree","genTracksTree");
  


  fMCInfo = new AliTPCGenInfo;
  fMCInfo->fReferences = new TClonesArray("AliTrackReference");

  fTreeGenTracks->Branch("MC","AliTPCGenInfo",&fMCInfo);

  fTreeGenTracks->Branch("fEventNr",&fEventNr,"fEventNr/I");

  fTreeGenTracks->AutoSave();
}
////////////////////////////////////////////////////////////////////////
void TPCFindGenTracks::CloseOutputFile() 
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
Int_t TPCFindGenTracks::TreeKLoop()
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
  for (Int_t iParticle = 0; iParticle < fNParticles; iParticle++) {

// load only particles with TR
    if (fReferenceIndex0[iParticle]<0) continue;
    //////////////////////////////////////////////////////////////////////
    fMCInfo->fLabel = iParticle;

    fMCInfo->fRow = (fContainerDigitRow[iParticle]);
    fMCInfo->fRowsWithDigits = fMCInfo->fRow.RowsOn();    
    if (fMCInfo->fRowsWithDigits<10) continue;
    fMCInfo->fRowsWithDigitsInn = fMCInfo->fRow.RowsOn(63); // 63 = number of inner rows
    fMCInfo->fRowsTrackLength = fMCInfo->fRow.Last() - fMCInfo->fRow.First();
    fMCInfo->fDigitsInSeed = 0;
    if (fMCInfo->fRow.TestRow(seedRow11) && fMCInfo->fRow.TestRow(seedRow12)) 
      fMCInfo->fDigitsInSeed = 1;
    if (fMCInfo->fRow.TestRow(seedRow21) && fMCInfo->fRow.TestRow(seedRow22)) 
      fMCInfo->fDigitsInSeed += 10;
    //
    //
    //
    fMCInfo->fParticle = *(stack->Particle(iParticle));
    //
    //
    fMCInfo->fLabel = iParticle;
    fMCInfo->fVDist[0] = fMCInfo->fParticle.Vx()-fVPrim[0];
    fMCInfo->fVDist[1] = fMCInfo->fParticle.Vy()-fVPrim[1];
    fMCInfo->fVDist[2] = fMCInfo->fParticle.Vz()-fVPrim[2]; 
    fMCInfo->fVDist[3] = TMath::Sqrt(fMCInfo->fVDist[0]*fMCInfo->fVDist[0]+
				      fMCInfo->fVDist[1]*fMCInfo->fVDist[1]+fMCInfo->fVDist[2]*fMCInfo->fVDist[2]);

    //
    Int_t index = fReferenceIndex0[iParticle];
    AliTrackReference  ref  = fReferences[index];
    // if (ref.GetTrack()!=iParticle)
    //  printf("problem2\n");
    fMCInfo->fTrackRef = ref;
  
    Int_t rfindex =0;
    if (fMCInfo->fReferences !=0) delete fMCInfo->fReferences;
    fMCInfo->fReferences = 
      new TClonesArray("AliTrackReference");
    fMCInfo->fReferences->ExpandCreateFast(fReferenceIndex1[iParticle]-fReferenceIndex0[iParticle]+1);
    //
    for (Int_t i = fReferenceIndex0[iParticle];i<=fReferenceIndex1[iParticle];i++){
      AliTrackReference  ref  = fReferences[i];
      AliTrackReference *ref2 = (AliTrackReference*) fMCInfo->fReferences->At(rfindex);
      if (ref.GetTrack()!=iParticle){
	//printf("problem5\n");	
	continue;
      }
      *ref2 = ref;
      rfindex++;
    }   
    //
    //
    ipdg = fMCInfo->fParticle.GetPdgCode();
    fMCInfo->fPdg = ipdg;
    ppdg = pdg.GetParticle(ipdg);
    // calculate primary ionization per cm
    if (ppdg){
      Float_t mass = ppdg->Mass();
      Float_t p = TMath::Sqrt(fMCInfo->fTrackRef.Px()*fMCInfo->fTrackRef.Px()+
			      fMCInfo->fTrackRef.Py()*fMCInfo->fTrackRef.Py()+
			      fMCInfo->fTrackRef.Pz()*fMCInfo->fTrackRef.Pz());
      
      //      Float_t betagama = fMCInfo->fParticle.P()/mass;
      Float_t betagama = p /mass;
      fMCInfo->fPrim = TPCBetheBloch(betagama);
    }	
    fMCInfo->fTPCdecay=kFALSE;
    if (fDecayRef[iParticle].GetTrack()>0){
      fMCInfo->fTRdecay  = fDecayRef[iParticle];
      fMCInfo->fDecayCoord[0] = fMCInfo->fTRdecay.X();
      fMCInfo->fDecayCoord[1] = fMCInfo->fTRdecay.Y();
      fMCInfo->fDecayCoord[2] = fMCInfo->fTRdecay.Z();
      if ( (fMCInfo->fTRdecay.R()<250)&&(fMCInfo->fTRdecay.R()>85) && (TMath::Abs(fMCInfo->fTRdecay.Z())<250) ){
	fMCInfo->fTPCdecay=kTRUE;
      }
    }
    else{
      fMCInfo->fTRdecay.SetTrack(-1);
      fMCInfo->fDecayCoord[0] = 0;
      fMCInfo->fDecayCoord[1] = 0;
      fMCInfo->fDecayCoord[2] = 0;
    }
    fMCInfo->Update();
    //////////////////////////////////////////////////////////////////////
    
    br->SetAddress(&fMCInfo);
    fTreeGenTracks->Fill();

  }
  fTreeGenTracks->AutoSave();

  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;

  return 0;
}




////////////////////////////////////////////////////////////////////////
Int_t TPCFindGenTracks::TreeDLoop()
{
//
// open the file with treeD
// loop over all entries there and save information about some tracks
//

  Int_t nInnerSector = fParamTPC->GetNInnerSector();
  Int_t rowShift = 0;
  //  Int_t zero=fParamTPC->GetZeroSup();
  Int_t zero=fParamTPC->GetZeroSup()+6;

  
  //char treeDName[100]; 
  //sprintf(treeDName,"TreeD_75x40_100x60_150x60_%d",fEventNr);
  //TTree *treeD=(TTree*)fFileTreeD->Get(treeDName);
  //
  //
  AliSimDigits digitsAddress, *digits=&digitsAddress;
  fTreeD->GetBranch("Segment")->SetAddress(&digits);

  Int_t sectorsByRows=(Int_t)fTreeD->GetEntries();
  if (fDebug > 1) cout<<"\tsectorsByRows = "<<sectorsByRows<<endl;
  for (Int_t i=0; i<sectorsByRows; i++) {
//  for (Int_t i=5720; i<sectorsByRows; i++) {
    if (!fTreeD->GetEvent(i)) continue;
    Int_t sec,row;
    fParamTPC->AdjustSectorRow(digits->GetID(),sec,row);
    if (fDebug > 1) cout<<sec<<' '<<row<<"                          \r";
//    cerr<<sec<<' '<<row<<endl;

// here I expect that upper sectors follow lower sectors
    if (sec > nInnerSector) rowShift = fParamTPC->GetNRowLow();
    //
    //digits->ExpandTrackBuffer();     
    //Int_t *tracks = digits->GetTracks();
    digits->First();    
    do {
      Int_t iRow=digits->CurrentRow();
      Int_t iColumn=digits->CurrentColumn();
      Short_t digitValue = digits->CurrentDigit();
//      cout<<"Inner loop: sector, iRow, iColumn "
//	  <<sec<<" "<<iRow<<" "<<iColumn<<endl;
      if (digitValue >= zero) {
	Int_t label;
	for (Int_t j = 0; j<3; j++) {
	  label = digits->GetTrackID(iRow,iColumn,j); 
	  //label = digits->GetTrackIDFast(iRow,iColumn,j); 
	  
	  if (label >= fNParticles) { //don't label from bakground event
	    continue;
	  }
	  if (label >= 0 && label <= fNParticles) {
//	  if (label >= 0 && label <= fDebug) {
	    if (fDebug > 6 ) {
	      cout<<"Inner loop: sector, iRow, iColumn, label, value, row "
		  <<sec<<" "
		  <<iRow<<" "<<iColumn<<" "<<label<<" "<<digitValue
		  <<" "<<row<<endl;
	    }	
	    fContainerDigitRow[label].SetRow(row+rowShift);
	  }
	}
      }
    } while (digits->Next());
  }

  if (fDebug > 2) cerr<<"end of TreeDLoop"<<endl;

  return 0;
}


////////////////////////////////////////////////////////////////////////
Int_t TPCFindGenTracks::TreeTRLoop()
{
  //
  // loop over TrackReferences and store the first one for each track
  //
  
  TTree * treeTR = fTreeTR;
  if (!treeTR) {
    cerr<<"TreeTR not found"<<endl;
    return 1;
  }
  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  if (fDebug > 1) cout<<"There are "<<nPrimaries<<" entries in TreeTR"<<endl;
  //
  //
  //track references for TPC
  TBranch *TPCBranchTR  = treeTR->GetBranch("TPC");
  if (!TPCBranchTR) {
    cerr<<"TPC branch in TR not found"<<endl;
    return 1;
  }
  TClonesArray* TPCArrayTR = new TClonesArray("AliTrackReference");
  TPCBranchTR->SetAddress(&TPCArrayTR);
  //get decay point if exist
  TBranch *runbranch  = treeTR->GetBranch("AliRun");
  if (!runbranch) {
    cerr<<"Run branch in TR not found"<<endl;
    return 1;
  }
  TClonesArray* RunArrayTR = new TClonesArray("AliTrackReference");
  runbranch->SetAddress(&RunArrayTR);
  //
  //
  //
  Int_t index     =  0;
  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    TPCBranchTR->GetEntry(iPrimPart);
    Float_t ptstart = 0;
    for (Int_t iTrackRef = 0; iTrackRef < TPCArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TPCArrayTR->At(iTrackRef);            
      //
      Int_t label = trackRef->GetTrack();
      if (label<0  || label > fNParticles) {
	continue;
      }
      if (fReferenceIndex0[label]<0) ptstart = trackRef->Pt();  //store pt at the TPC entrance
      if (ptstart<fgPtCut) continue;

      if (index>=fgMaxTR) continue;     //restricted size of buffer for TR
      fReferences[index] = *trackRef;
      fReferenceIndex1[label] = index;  // the last ref with given label
      if (fReferenceIndex0[label]==-1) fReferenceIndex0[label] = index;   //the first index with label
      index++;           
    }
    // get dacay position
    runbranch->GetEntry(iPrimPart);    
    for (Int_t iTrackRef = 0; iTrackRef < RunArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)RunArrayTR->At(iTrackRef);      
      //
      if (trackRef->Pt() < fgPtCut) continue;      
      Int_t label = trackRef->GetTrack();
      if (label<0  || label > fNParticles) {
	continue;
      }
      if (trackRef->R()>450) continue;   //not decay  in TPC
      if (trackRef->Z()>450) continue;   //not decay  in TPC
      if (!trackRef->TestBit(BIT(2))) continue;  //if not decay
      
      if (label>=fgMaxTR) continue;     //restricted size of buffer for TR
      fDecayRef[label] = *trackRef;      
    }

  }
  TPCArrayTR->Delete();
  delete  TPCArrayTR;
  RunArrayTR->Delete();
  delete  RunArrayTR;

  return 0;
}

////////////////////////////////////////////////////////////////////////
Float_t TPCFindGenTracks::TR2LocalX(AliTrackReference *trackRef,
				    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return x[0];
}
////////////////////////////////////////////////////////////////////////



  
////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr(const char* fnGenTracks,
		   const char* fnCmp,
		   const char* fnGalice,
		   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  //  fFnGenTracks = fnGenTracks;
  //  fFnCmp = fnCmp;
  sprintf(fFnGenTracks,"%s",fnGenTracks);
  sprintf(fFnCmp,"%s",fnCmp);

  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    //delete gAlice->GetRunLoader();
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


////////////////////////////////////////////////////////////////////////
TPCCmpTr::~TPCCmpTr()
{
  if (fLoader) {
    delete fLoader;
  }
}

//////////////////////////////////////////////////////////////
Int_t TPCCmpTr::SetIO()
{
  //
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

Int_t TPCCmpTr::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  //gAlice->GetEvent(eventNr);
  fLoader->SetEventNumber(eventNr);  
  //
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->LoadTracks();
  fTreeRecTracks = tpcl->TreeT();
  //
  return 0;
}



////////////////////////////////////////////////////////////////////////
void TPCCmpTr::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeCmp = 0;
  //  fFnCmp = "cmpTracks.root";
  fFileGenTracks = 0;
  fDebug = 0;
  //
  fParamTPC = 0;
  fTreeRecTracks = 0;
  fTreePoints =0;

  fTPCTrack = 0; 
  fTracks   = 0;
  fTrackPoints =0;
}

////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::Exec()
{
  TStopwatch timer;
  timer.Start();

  if (SetIO()==1) 
    return 1;
   
  fNextTreeGenEntryToRead = 0;
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    SetIO(fEventNr);
    fNParticles = gAlice->GetEvent(fEventNr);    

    fIndexRecTracks = new Int_t[fNParticles*4];  //write at maximum 4 tracks corresponding to particle
    fFakeRecTracks = new Int_t[fNParticles];
    fMultiRecTracks = new Int_t[fNParticles];
    for (Int_t i = 0; i<fNParticles; i++) {
      fIndexRecTracks[4*i] = -1;
      fIndexRecTracks[4*i+1] = -1;
      fIndexRecTracks[4*i+2] = -1;
      fIndexRecTracks[4*i+3] = -1;

      fFakeRecTracks[i] = 0;
      fMultiRecTracks[i] = 0;
    }
  
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\tStart loop over TreeT"<<endl;
    if (TreeTLoop(eventNr)>0) return 1;

    if (fDebug>2) cout<<"\tStart loop over tree genTracks"<<endl;
    if (TreeGenLoop(eventNr)>0) return 1;
    if (fDebug>2) cout<<"\tEnd loop over tree genTracks"<<endl;

    delete fTreeRecTracks;
    delete [] fIndexRecTracks;
    delete [] fFakeRecTracks;
    delete [] fMultiRecTracks;
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  timer.Stop();
  timer.Print();
  return 0;

}
////////////////////////////////////////////////////////////////////////
Bool_t TPCCmpTr::ConnectGenTree()
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
  fMCInfo = new AliTPCGenInfo;
  fMCInfo->fReferences = new TClonesArray("AliTrackReference");  
  fTreeGenTracks->SetBranchAddress("MC",&fMCInfo);

  //
  //fTreeGenTracks->SetBranchAddress("fEventNr",&fEventNr);


  if (fDebug > 1) {
    cout<<"Number of gen. tracks with TR: "<<fTreeGenTracks->GetEntries()<<endl;
  }
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
void TPCCmpTr::CreateTreeCmp() 
{
  fFileCmp = TFile::Open(fFnCmp,"RECREATE");
  if (!fFileCmp) {
    cerr<<"Error in CreateTreeCmp: cannot open file "<<fFnCmp<<endl;
    return;
  }


  fTreeCmp    = new TTree("TPCcmpTracks","TPCcmpTracks");
  //
  fMCInfo = new AliTPCGenInfo;
  fMCInfo->fReferences = new TClonesArray("AliTrackReference");
  fRecInfo = new AliTPCRecInfo;
  //
  fTPCTrack = new AliTPCtrack;
   //
  fTreeCmp->Branch("MC","AliTPCGenInfo",&fMCInfo);
  fTreeCmp->Branch("RC","AliTPCRecInfo",&fRecInfo);
  fTreeCmp->Branch("fEventNr",&fEventNr,"fEventNr/I");
  fTreeCmp->Branch("fTPCTrack","AliTPCtrack",&fTPCTrack);
  //
  fTreeCmp->AutoSave(); 

}
////////////////////////////////////////////////////////////////////////
void TPCCmpTr::CloseOutputFile()  
{
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

TVector3 TPCCmpTr::TR2Local(AliTrackReference *trackRef,
			    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return TVector3(x);
}
////////////////////////////////////////////////////////////////////////

Int_t TPCCmpTr::TreeTLoop(Int_t eventNr)
{
  //
  // loop over all TPC reconstructed tracks and store info in memory
  //
  TStopwatch  timer;
  timer.Start();
  
  if (!fTreeRecTracks) {
    cerr<<"Can't get a tree with TPC rec. tracks  "<<endl;
    return 1;
  }
  //fTreePoints=(TTree*)fFileRecTracks->Get("trackDebug");
  
  Int_t nEntries = (Int_t) fTreeRecTracks->GetEntries();
  if (fDebug > 2) cout<<"Event, rec. tracks: "<<eventNr<<" "
		      <<nEntries<<endl;
  TBranch * br= fTreeRecTracks->GetBranch("tracks");
  br->SetAddress(&fTPCTrack);
  TBranch *brp = 0;
  if (fTreePoints) brp = fTreePoints->GetBranch("debug");

  if (fTracks){
    fTracks->Delete();    
    delete fTracks;
  }
  if (fTrackPoints){
    fTrackPoints->Delete();
    delete fTrackPoints;
    fTrackPoints =0;
  }
  fTracks      = new TObjArray(nEntries);
  if (brp){
    fTrackPoints = new TObjArray(nEntries);
  }
  else fTrackPoints = 0;

  //
  //load tracks to the memory
  for (Int_t i=0; i<nEntries;i++){
    AliTPCtrack * track = new AliTPCtrack;
    br->SetAddress(&track);
    br->GetEntry(i);
    fTracks->AddAt(track,i);
  }
  //
  //load track points to the memory
  if (brp) for (Int_t i=0; i<nEntries;i++){
    TClonesArray * arr = new TClonesArray("AliTPCTrackPoint2");
    brp->SetAddress(&arr);
    brp->GetEntry(i);
    if (arr!=0)
      for (Int_t j=0;j<arr->GetEntriesFast();j++){
	AliTPCTrackPoint2 * point = (AliTPCTrackPoint2*)arr->UncheckedAt(j);
	if (point && point->fID>=0){
	  fTrackPoints->AddAt(arr,point->fID);
	  break;
	}
      }    
  }
  //

  //
  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    //br->GetEntry(iEntry);
    fTPCTrack = (AliTPCtrack*)fTracks->UncheckedAt(iEntry);
    //
    Int_t label = fTPCTrack->GetLabel();
    Int_t absLabel = abs(label);
    if (absLabel < fNParticles) {
      //      fIndexRecTracks[absLabel] =  iEntry;
      if (label < 0) fFakeRecTracks[absLabel]++;
      
      if (fMultiRecTracks[absLabel]>0){
	if (fMultiRecTracks[absLabel]<4)
	  fIndexRecTracks[absLabel*4+fMultiRecTracks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecTracks[absLabel*4] =  iEntry;
      fMultiRecTracks[absLabel]++;
    }
  }  
  printf("Time spended in TreeTLoop\n");
  timer.Print();
  
  if (fDebug > 2) cerr<<"end of TreeTLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::TreeGenLoop(Int_t eventNr)
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
  branch->SetAddress(&fRecInfo); // set all pointers

  while (entry < nParticlesTR) {
    fTreeGenTracks->GetEntry(entry);
    entry++;
    if (fEventNr < eventNr) continue;
    if (fEventNr > eventNr) break;
    fNextTreeGenEntryToRead = entry-1;
    if (fDebug > 2 && fMCInfo->fLabel < 10) {
      cerr<<"Fill track with a label "<<fMCInfo->fLabel<<endl;
    }

    fRecInfo->Reset();


    if (fIndexRecTracks[fMCInfo->fLabel*4] >= 0) {
      fTPCTrack= (AliTPCtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4]);
      
      //      if (nBytes > 0) {
      if (fTPCTrack) {
	//
	TVector3 local = TR2Local(&(fMCInfo->fTrackRef),fParamTPC);
	local.GetXYZ(fRecInfo->fTRLocalCoord);
	//
	// find nearest track if multifound
	if (fIndexRecTracks[fMCInfo->fLabel*4]+1){
	  Float_t dz = TMath::Abs(local.Z()-fTPCTrack->GetZ());
	  for (Int_t i=1;i<4;i++){
	    if (fIndexRecTracks[fMCInfo->fLabel*4+i]>=0){
	      AliTPCtrack * track = (AliTPCtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4+i]);
	      if  (TMath::Abs(local.Z()-track->GetZ())<dz)
		fTPCTrack = track;		   
	    }
	  }
	}
	fRecInfo->fTP=0;
	if (fTrackPoints){
	  Int_t id = fTPCTrack->GetUniqueID();
	  if (fTrackPoints->UncheckedAt(id)){
	    fRecInfo->fTP = (TClonesArray*)fTrackPoints->UncheckedAt(id);
	    //	    fTrackPoints->AddAt(0,id);   //not owner anymore
	  }
	}
	fRecInfo->fTPCTrack =*fTPCTrack; 
	fRecInfo->fReconstructed = 1;
	fRecInfo->fFake = fFakeRecTracks[fMCInfo->fLabel];
	fRecInfo->fMultiple = fMultiRecTracks[fMCInfo->fLabel];
	fRecInfo->fdEdx = fTPCTrack->GetdEdx();
	//
	fRecInfo->fTPCTrack.PropagateTo(local.X());
	Double_t par[5];
	//
	Double_t localX = local.X();
	fTPCTrack->GetExternalParameters(localX,par);
	fRecInfo->fRecPhi=TMath::ASin(par[2]) + fTPCTrack->GetAlpha();
	if (fRecInfo->fRecPhi<0) fRecInfo->fRecPhi+=2*TMath::Pi();
	if (fRecInfo->fRecPhi>=2*TMath::Pi()) fRecInfo->fRecPhi-=2*TMath::Pi();
//	  fRecInfo->fRecPhi = (fRecInfo->fRecPhi)*kRaddeg;
	fRecInfo->fLambda = TMath::ATan(par[3]);
	fRecInfo->fRecPt_1 = TMath::Abs(par[4]);
      }

    } 

    fTreeCmp->Fill();
  }
  fTreeCmp->AutoSave();
  fTracks->Delete();
  fTPCTrack =0;
  if (fTrackPoints){
    fTrackPoints->Delete();
    delete fTrackPoints;
    fTrackPoints =0;
  } 
  printf("Time spended in TreeKLoop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void AliTPCComparisonDraw::SetIO(const char *fname)
{
  //
   TFile* file = TFile::Open(fname);
   if (!file) {
     printf("Could not open file  - generated new one\n"); 
     TFile* filegen = TFile::Open("genTracks.root");
     if (!filegen){
       printf("FILE with MC information is generated\n"); 
       TPCFindGenTracks *t = new TPCFindGenTracks("galice.root","genTracks.root",0);
       t->Exec();
       delete t;
     }
     filegen = TFile::Open("genTracks.root");
     if (!filegen){
       printf("COMPARISON  FILE COULDNT BE GENERATED \n");
       return;
     }
     printf("COMPARISON  FILE IS GENERATED \n");
     TPCCmpTr *t2 = new TPCCmpTr("genTracks.root","cmpTracks.root","galice.root",0);
     t2->Exec();
   }
   file = TFile::Open(fname);
   if (!file){
     printf("Comparison file couldn't be generated\n"); 
     return;
   }
   //
   fTree = (TTree*) file->Get("TPCcmpTracks");
   if (!fTree) {
    printf("no track comparison tree found\n");
    file->Close();
    delete file;
  }
}

void AliTPCComparisonDraw::ResPt()
{
 //
  //
  gStyle->SetOptFit();
  TCanvas *c1=new TCanvas("TPC pt resolution","TPC pt resolution",0,0,700,850);
  c1->Draw(); c1->cd();
  TPad *p1=new TPad("p1","p1",0.01,0.51,.99,.99); 
  p1->Draw();p1->cd(); 
  p1->SetGridx(); p1->SetGridy();
  //
  c1->cd();
  TPad *p2=new TPad("p2","p2",0.01,0.01,.99,.49); 
  p2->Draw();p2->cd();
  p2->SetGridx(); p2->SetGridy();
  // 
  //Default cuts
  TCut cprim("cprim","MC.fVDist[3]<1");
  TCut cnprim("cnprim","MC.fVDist[3]>1");
  TCut cteta1("cteta1","abs(MC.fTrackRef.Theta()/3.1415-0.5)<0.25");
  TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
  //
  c1->cd();  p1->cd();
  TH1F* hisp = ResPtvsPt("MC.fRowsWithDigits>100"+cteta1+cpos1,"1",0.15,2.,6);
  c1->cd(); 
  c1->Draw();
  p1->cd();
  p1->Draw();
  hisp->Draw(); 
  //
  c1->cd();  
  c1->Draw();
  p2->cd();
  p2->Draw();
  TH1F* his2 =  new TH1F("Ptresolution","Ptresolution",40,-5.,5.);
  fTree->Draw("100.*(abs(1./fTPCTrack.Get1Pt())-MC.fTrackRef.Pt())/MC.fTrackRef.Pt()>>Ptresolution","MC.fRowsWithDigits>100&&RC.fTPCTrack.fN>50&&RC.fMultiple==1"+cteta1+cpos1+cprim);
  AliLabelAxes(his2, "#Delta p_{t} / p_{t} [%]", "entries");
  his2->Fit("gaus");
  his2->Draw();
 
}

void AliTPCComparisonDraw::Eff()
{
  //
  //
  TCanvas *c1=new TCanvas("TPC efficiency","TPC efficiency",0,0,700,850);
  c1->Draw(); c1->cd();
  TPad *p1=new TPad("p1","p1",0.01,0.51,.99,.99); 
  p1->Draw();p1->cd(); 
  p1->SetGridx(); p1->SetGridy();
  //
  c1->cd();
  TPad *p2=new TPad("p2","p2",0.01,0.01,.99,.49); 
  p2->Draw();p2->cd();
  p2->SetGridx(); p2->SetGridy();
  // 
  //Default cuts
  TCut cprim("cprim","MC.fVDist[3]<1");
  TCut cnprim("cnprim","MC.fVDist[3]>1");
  TCut cteta1("cteta1","abs(MC.fTrackRef.Theta()/3.1415-0.5)<0.25");
  TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
  //
  c1->cd();  p1->cd();
  TH1F* hisp = EffVsPt("MC.fRowsWithDigits>100"+cteta1+cpos1+cprim,"RC.fTPCTrack.fN>50");
  hisp->Draw();
  //hisp->DrawClone(); 
  TText * text = new TText(0.25,102.,"Primary particles");
  text->SetTextSize(0.05);
  text->Draw();
  //
  c1->cd();  p2->cd();
  TH1F* hiss = EffVsPt("MC.fRowsWithDigits>100"+cteta1+cpos1+cnprim,"RC.fTPCTrack.fN>50");
  hiss->Draw();
  //hiss->DrawClone();
  text = new TText(0.25,102.,"Secondary particles");
  text->SetTextSize(0.05);
  text->Draw();

 
}


TH1F * AliTPCComparisonDraw::EffVsPt(const char* selection, const char * quality, Float_t min, Float_t max)
{
  //
  //
  Int_t nBins = 10;
  Double_t* bins = CreateLogBins(nBins, min, max);
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nBins, bins);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nBins, bins);
  
  fTree->Draw("MC.fParticle.Pt()>>hGen", selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && RC.fReconstructed && %s", selection, quality);
  fTree->Draw("MC.fParticle.Pt()>>hRec", selectionRec, "groff");

  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, "p_{t} [GeV/c]", "#epsilon [%]");

  delete hRec;
  delete hGen;
  delete[] bins;
  return hEff;
}


TH1F * AliTPCComparisonDraw::EffVsRM(const char* selection, const char * quality, Float_t min, Float_t max, Int_t nBins)
{
  //
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nBins, min, max);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nBins, min, max);
  //  
  fTree->Draw("sqrt(MC.fDecayCoord[0]*MC.fDecayCoord[0] + MC.fDecayCoord[1]*MC.fDecayCoord[1])>>hGen", selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && RC.fReconstructed && %s", selection, quality);
  fTree->Draw("sqrt(MC.fDecayCoord[0]*MC.fDecayCoord[0] + MC.fDecayCoord[1]*MC.fDecayCoord[1])>>hRec", selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, "r_{vertex} [cm]", "#epsilon [%]");
  //
  delete hRec;
  delete hGen;
  return hEff;
}

TH1F * AliTPCComparisonDraw::EffVsRS(const char* selection, const char * quality, Float_t min, Float_t max, Int_t nBins)
{
  //
  TH1F* hGen = new TH1F("hGen", "gen. tracks", nBins, min, max);
  TH1F* hRec = new TH1F("hRec", "rec. tracks", nBins, min, max);
  //  
  fTree->Draw("sqrt(MC.fVDist[0]*MC.fVDist[0] + MC.fVDist[1]*MC.fVDist[1])>>hGen", selection, "groff");
  char selectionRec[256];
  sprintf(selectionRec, "%s && RC.fReconstructed && %s", selection, quality);
  fTree->Draw("sqrt(MC.fVDist[0]*MC.fVDist[0] + MC.fVDist[1]*MC.fVDist[1])>>hRec", selectionRec, "groff");
  //
  TH1F* hEff = CreateEffHisto(hGen, hRec);
  AliLabelAxes(hEff, "r_{vertex} [cm]", "#epsilon [%]");
  //
  delete hRec;
  delete hGen;
  return hEff;

}

TH1F * AliTPCComparisonDraw::ResPtvsPt(const char* selection, const char * quality, Float_t min, Float_t max, Int_t nBins)
{
  //
  Double_t* bins = CreateLogBins(nBins, min, max);
  Int_t nBinsRes = 30;
  Double_t maxRes = 10.;
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nBins, bins, nBinsRes, -maxRes, maxRes);
  
  fTree->Draw("100.*(abs(1./fTPCTrack.Get1Pt())-MC.fTrackRef.Pt())/MC.fTrackRef.Pt():MC.fTrackRef.Pt()>>hRes2", selection, "groff");

  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, "p_{t} [GeV/c]", "#Delta p_{t} / p_{t} [%]");
  //
  delete hRes2;
  delete[] bins;
  return hRes;

}

TH1F * AliTPCComparisonDraw::MeanPtvsPt(const char* selection, const char * quality, Float_t min, Float_t max, Int_t nBins)
{
  //
  Double_t* bins = CreateLogBins(nBins, min, max);
  Int_t nBinsRes = 30;
  Double_t maxRes = 10.;
  TH2F* hRes2 = new TH2F("hRes2", "residuals", nBins, bins, nBinsRes, -maxRes, maxRes);
  
  fTree->Draw("100.*(1./fTPCTrack.Get1Pt()-MC.fTrackRef.Pt())/MC.fTrackRef.Pt():MC.fTrackRef.Pt()>>hRes2", selection, "groff");

  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, "p_{t} [GeV/c]", "mean p_{t} / p_{t} [%]");
  //
  delete hRes2;
  delete[] bins;
  if (!hMean) return 0;
  return hMean;

}


TH1F * AliTPCComparisonDraw::ResdEdxvsN(const char* selection, const char * quality, Float_t min, Float_t max, Int_t nBins)
{
  //
  Int_t nBinsRes = 15;

  TH2F* hRes2 = new TH2F("hRes2", "residuals", nBins, min, max, nBinsRes, 34, 60);
  
  fTree->Draw("RC.fTPCTrack.fdEdx/MC.fPrim:RC.fTPCTrack.fN>>hRes2", selection, "groff");

  TH1F* hMean=0;
  TH1F* hRes = CreateResHisto(hRes2, &hMean);
  AliLabelAxes(hRes, "N points", "sigma dEdx/Nprim [%]");
  //
  delete hRes2;
  return hRes;

}



Double_t* AliTPCComparisonDraw::CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax)
{
  Double_t* bins = new Double_t[nBins+1];
  bins[0] = xMin;
  Double_t factor = pow(xMax/xMin, 1./nBins);
  for (Int_t i = 1; i <= nBins; i++)
    bins[i] = factor * bins[i-1];
  return bins;
}




void AliTPCComparisonDraw::AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle)
{
  //
  histo->GetXaxis()->SetTitle(xAxisTitle);
  histo->GetYaxis()->SetTitle(yAxisTitle);
}


TH1F* AliTPCComparisonDraw::CreateEffHisto(TH1F* hGen, TH1F* hRec)
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



TH1F* AliTPCComparisonDraw::CreateResHisto(TH2F* hRes2, TH1F **phMean,  Bool_t drawBinFits, 
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
    if (hBin->GetEntries() > 10) {
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



