#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include <stdio.h>
#include <string.h>
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliTPCtrack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "TParticle.h"
#include "AliTPC.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
//
#include "iostream.h"  
#include "Rtypes.h"
#include "TSystem.h"
#include "TTimer.h"
#include "Getline.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "AliRun.h"
#include "AliTPCComparisonMI.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliComplexCluster.h"


#endif
 
//
//
Bool_t ImportgAlice(TFile *file);
TFile* OpenAliceFile(char *fn, Bool_t importgAlice = kFALSE, char *mode = "read");

AliTPCParam * GetTPCParam(){
  AliTPCParamSR * par = new AliTPCParamSR;
  par->Update();
  return par;
}



////////////////////////////////////////////////////////////////////////
Bool_t ImportgAlice(TFile *file) {
// read in gAlice object from the file
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice)  return kFALSE;
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
TFile* OpenAliceFile(char *fn, Bool_t importgAlice, char *mode) {
  TFile *file = TFile::Open(fn,mode);
  if (!file->IsOpen()) {
    cerr<<"OpenAliceFile: cannot open file "<<fn<<" in mode "
	<<mode<<endl;
    return 0;
  }
  if (!importgAlice) return file;
  if (ImportgAlice(file)) return file;
  return 0;
}
////////////////////////////////////////////////////////////////////////

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
}

AliTPCGenInfo::~AliTPCGenInfo()
{
  if (fReferences) delete fReferences;
  fReferences =0;
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
TPCFindGenTracks::TPCFindGenTracks(char * fnGalice, char* fnRes,
				   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  fFnRes = fnRes;
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
void TPCFindGenTracks::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeGenTracks = 0;
  fFnRes = "genTracks.root";
  fFileGenTracks = 0;
  fContainerDigitRow = 0;
  //
  fReferences = 0;
  fReferenceIndex0 = 0;
  fReferenceIndex1 = 0;
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
  ;
}


Int_t  TPCFindGenTracks::SetIO()
{
  //
  // 
  CreateTreeGenTracks();
  if (!fTreeGenTracks) return 1;
  AliTracker::SetFieldFactor(); 
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
  gAlice->GetEvent(eventNr);
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
    fReferences = new AliTrackReference[fgMaxTR];
    fReferenceIndex0 = new Int_t[fNParticles];
    fReferenceIndex1 = new Int_t[fNParticles];

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
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  delete gAlice;

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
    if (fMCInfo->fParticle.GetLastDaughter()>0&&fMCInfo->fParticle.GetLastDaughter()<fNParticles){
      TParticle * dparticle0 = 0;
      if (fMCInfo->fParticle.GetDaughter(0)>0) dparticle0 = stack->Particle(fMCInfo->fParticle.GetDaughter(0));
      TParticle * dparticle1 = 0;
      if (fMCInfo->fParticle.GetDaughter(1)>0) dparticle1 = stack->Particle(fMCInfo->fParticle.GetDaughter(1));
      if ((dparticle1)&&(dparticle1->P()>0.15)) {
	fMCInfo->fDecayCoord[0] = dparticle1->Vx();
	fMCInfo->fDecayCoord[1] = dparticle1->Vy();
	fMCInfo->fDecayCoord[2] = dparticle1->Vz();
      }
      else 
	fMCInfo->fDecayCoord[0]=1000000;
    }else
      fMCInfo->fDecayCoord[0]=1000000;
    fMCInfo->fParticle = *(stack->ParticleFromTreeK(iParticle));
    //
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
      if (ref.GetTrack()!=iParticle)
	printf("problem5\n");	
      *ref2 = ref;
      rfindex++;
    }   
    //
    //
    ipdg = fMCInfo->fParticle.GetPdgCode();
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
  TBranch *TPCBranchTR  = treeTR->GetBranch("TPC");
  if (!TPCBranchTR) {
    cerr<<"TPC branch in TR not found"<<endl;
    return 1;
  }
  TClonesArray* TPCArrayTR = new TClonesArray("AliTrackReference");
  TPCBranchTR->SetAddress(&TPCArrayTR);
  //
  //
  //
  Int_t index     =  0;
  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    TPCBranchTR->GetEntry(iPrimPart);
    for (Int_t iTrackRef = 0; iTrackRef < TPCArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TPCArrayTR->At(iTrackRef);      
      //
      if (trackRef->Pt() < fgPtCut) continue;      
      Int_t label = trackRef->GetTrack();
      if (label<0  || label > fNParticles) {
	continue;
      }
      if (index>=fgMaxTR) continue;     //restricted size of buffer for TR
      fReferences[index] = *trackRef;
      fReferenceIndex1[label] = index;  // the last ref with given label
      if (fReferenceIndex0[label]==-1) fReferenceIndex0[label] = index;   //the first index with label
      index++;
           
    }
  }
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
TPCCmpTr::TPCCmpTr(char* fnGenTracks,
		   char* fnCmp,
		   char* fnGalice,
		   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFnGenTracks = fnGenTracks;
  fFnCmp = fnCmp;
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  
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

//////////////////////////////////////////////////////////////
Int_t TPCCmpTr::SetIO()
{
  //
  // 
  CreateTreeCmp();
  if (!fTreeCmp) return 1;
  AliTracker::SetFieldFactor(); 
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
  gAlice->GetEvent(eventNr);
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
  fFnCmp = "cmpTracks.root";
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
TPCCmpTr::~TPCCmpTr()
{
  ;
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
  delete gAlice;

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
  
  Int_t entry = fNextTreeGenEntryToRead;
  Double_t nParticlesTR = fTreeGenTracks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextTreeGenEntryToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextTreeGenEntryToRead<<endl;
  TBranch * branch = fTreeCmp->GetBranch("RC");
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
	if (fRecInfo->fRecPhi<0) fRecInfo->fRecPhi+=2*kPI;
	if (fRecInfo->fRecPhi>=2*kPI) fRecInfo->fRecPhi-=2*kPI;
//	  fRecInfo->fRecPhi = (fRecInfo->fRecPhi)*kRaddeg;
	fRecInfo->fLambda = TMath::ATan(par[3]);
	fRecInfo->fRecPt_1 = TMath::Abs(par[4]);
      }

    } 

    branch->SetAddress(&fRecInfo); // set all pointers
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
  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
