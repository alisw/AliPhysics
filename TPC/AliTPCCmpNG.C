/// \file AliTPCCmpNG.C
///
/// version: 1.0
/// description:
///        define a class TPCGenTrack
///        save TPC related properties of tracks into a single tree
///
/// input:
///        Int_t nEvents      ... nr of events to process
///        Int_t firstEventNr ... first event number (starts from 0)
///        char* fnRecTracks .. name of file with reconstructed tracks
///        char* fnHits ... name of file with hits and Kine Tree
///        char* fnDigits  ... name of file with digits
///        char* fnTracks .. output file name, default genTracks.root
///
/// How to use:
///  Typical usage:
///    .L AliTPCCmpNG.C+
///    TPCFindGenTracks *t = new TPCFindGenTracks("galice.root","tpc.digits.root")
///    t->Exec();
///    .q
///    aliroot
///    .L  AliTPCCmpNG.C+
///    TPCCmpTr *t2 = new TPCCmpTr("tpc.tracks.root","genTracks.root","cmpTracks.root");
///    t2->Exec();
///
///  Details:
///
///  Step 1 - summurize information from simulation
///
///  Compile macro with ACLIC:
///     .L AliTPCCmpNG.C+
///  create an object TPCFindGenTracks, which processes information
///  from simulations. As input it needs:
///     object gAlice: to get magnetic field
///     TreeK: to get parameters of generated particles
///     TreeTR: to get track parameters at the TPC entry
///     TreeD: to get number of digits and digits pattern
///                for a given track in TPC
///  These input objects can be in different files, gAlice, TreeK and
///  TreeTR are in the file fnHits, TreeD in the file fnDigits (can be
///  the same as fnHits. Output is written to the file fnRes
///  ("genTracks.root" by default). Use can specify number of
///  events to process and the first event number:
///    TPCFindGenTracks *t = new TPCFindGenTracks("galice.root","tpc.digits.root","genTracks.root",1,0)
///  The filenames in the example on previous line are defaults, user can
///  specify just the file name with gAlice object (and TreeTR and TreeK),
///  so equivalent shorter initialization is:
///    TPCFindGenTracks *t = new TPCFindGenTracks("galice.root")
///  The task is done by calling Exec() method:
///    t->Exec();
///  User can set different debug levels by invoking:
///    t->SetDebug(debugLevel)
///  Number of events to process and the first event number can be
///  specified as parameters to Exec:
///    t->Exec(nEvents, firstEvent)
///  Then you have to quit root to get rid of problems with deleting gAlice
///  object (it is not deleted, but read again in the following step):
///
///  Step 2 - compare reconstructed tracks with simulated
///
///  Load (and compile) the macro:
///   .L AliTPCCmpNG.C+
///  Create object TPCCmpTr, which does the comparison. As input it requires
///  name of the file with reconstructed TPC tracks. You can specify
///  name of the file with genTrack tree (summarized info about simulation),
///  file with gAlice object, output file name, number of events to process
///  and first event number:
///  TPCCmpTr *t2 = new TPCCmpTr("tpc.tracks.root","genTracks.root","cmpTracks.root","galice.root",1,0);
///  The interface is quite similar to the TPCFindGenTracks class.
///  Then just invoke Exec() method:
///  t2->Exec();
///
///  Step 3 - study the results
///
///  Load the outoput TTree and you can do Draw(), Scan() or other
///  usual things to do with TTree:
///  TFile *f = new TFile("cmpTracks.root")
///  TTree *t = (TTree*)f->Get("TPCcmpTracks")
///  t->Draw("fVDist[3]","fReconstructed")
///
/// History:
///
/// 24.09.02 - first version
/// 24.01.03 - v7, before change from TPC Special Hits to TrackReferences
/// 26.01.03 - change from TPC Special Hits to TrackReferences
///            (loop over TreeTR instead of TreeH)
/// 28.01.03 - v8 last version before removing TPC special point
/// 28.01.03 - remove TPC special point, loop over TreeH
///            store TParticle and AliTrack
/// 29.01.03 - v9 last version before moving the loop over rec. tracks
///            into separate step
/// 03.02.03 - rename to AliTPCCmpNG.C, remove the part with rec. tracks
///            (will be addded in a macro AliTPCCmpTr.C
///
/// \author Jiri Chudoba
/// \date 24.09.2002

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
// #include "AliConst.h"
#endif

#include "AliTPCTracking.C"

// include ML.C:
////////////////////////////////////////////////////////////////////////
//
// name: ML.C  (Macro Library - collection of functions used in diff macros
// date: 08.05.2002
// last update: 08.05.2002
// author: Jiri Chudoba
// version: 1.0
// description: 
//
// History:
//
// 08.05.02 - first version
//
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "Rtypes.h"
#include "TSystem.h"
#include "TTimer.h"
#include "Getline.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "AliRun.h"

void WaitForReturn();
Bool_t ImportgAlice(TFile *file);
TFile* OpenAliceFile(char *fn, Bool_t importgAlice = kFALSE, char *mode = "read");
Int_t Chain(TString baseDir, TString subDirNameMask, TString fn, TChain& chain);

#endif

////////////////////////////////////////////////////////////////////////
void WaitForReturn() 
{
/// wait until user press return;

  char    *input;
  Bool_t done = kFALSE;
  TTimer  *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);

  do {
    timer->TurnOn();
    timer->Reset();
    // Now let's read the input, we can use here any
    // stdio or iostream reading methods. like std::cin >> myinputl;
    input = Getline("Type <return> to continue: "); 
    timer->TurnOff();
    if (input) done = kTRUE;
  } while (!done);
}

////////////////////////////////////////////////////////////////////////

Int_t Chain(TString baseDir, TString subDirNameMask, TString fn, TChain& chain) {
/// chain all files fn in the subdirectories which match subDirName
/// return number of chained files

// open baseDir, loop over subdirs

  void *dirp = gSystem->OpenDirectory(baseDir); 
  if (!dirp) {
    cerr<<"Could not open base directory "<<baseDir.Data()<<endl;
    return 0;
  }
  const char *afile;
  Long_t id, size, flag, modTime;
  Int_t rc = 0;
  char relName[100];
  while((afile = gSystem->GetDirEntry(dirp))) {
//    printf("file: %s\n",afile);
    if (strstr(afile,subDirNameMask.Data())) {
      sprintf(relName,"%s/%s/%s",baseDir.Data(),afile,fn.Data());
//      cerr<<"relName = "<<relName<<endl;
      if(!gSystem->GetPathInfo(relName, &id, &size, &flag, &modTime)) { 
	rc += chain.Add(relName);      
      }
    }
  }
  gSystem->FreeDirectory(dirp);
//  cerr<<"rc = "<<rc<<endl;
  return rc;
}
////////////////////////////////////////////////////////////////////////
Bool_t ImportgAlice(TFile *file) {
/// read in gAlice object from the file

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


////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class TPCGenInfo
//
////////////////////////////////////////////////////////////////////////
class TPCGenInfo: public TObject {

public:
  TPCGenInfo();
  AliTrackReference *fTrackRef;   ///< track reference saved in the output tree
  TParticle *fParticle;           ///< generated particle
  Int_t fLabel;                   ///< track label

  Int_t fRowsWithDigitsInn;    ///< number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;       ///< number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;      ///< last - first row with digit
  Int_t fDigitsInSeed;         ///< digits in the default seed rows

  ClassDef(TPCGenInfo,1)  // container for 
};
ClassImp(TPCGenInfo)
////////////////////////////////////////////////////////////////////////
TPCGenInfo::TPCGenInfo()
{
  fTrackRef = 0;
  fParticle = 0;
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// End of implementation of the class TPCGenInfo
//
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class digitRow
//
////////////////////////////////////////////////////////////////////////
const Int_t kgRowBytes = 32;

class digitRow: public TObject {

public:
  digitRow();
//  digitRow(){;}
  virtual ~digitRow(){;}
  void SetRow(Int_t row);
  Bool_t TestRow(Int_t row);
  digitRow & operator=(const digitRow &digOld);
  Int_t RowsOn(Int_t upto=8*kgRowBytes);
  Int_t Last();
  Int_t First();
  void Reset();

//private:
  UChar_t fDig[kgRowBytes];

  ClassDef(digitRow,1)  // container for digit pattern
};
ClassImp(digitRow)
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
/// return kTRUE if row is on

  Int_t iC = row/8;
  Int_t iB = row%8;
  return TESTBIT(fDig[iC],iB);
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::RowsOn(Int_t upto)
{
/// returns number of rows with a digit
/// count only rows less equal row number upto

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
/// resets all rows to zero

  for (Int_t i = 0; i<kgRowBytes; i++) {
    fDig[i] <<= 8;
  }
}
////////////////////////////////////////////////////////////////////////
Int_t digitRow::Last()
{
/// returns the last row number with a digit
/// returns -1 if now digits

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
/// returns the first row number with a digit
/// returns -1 if now digits

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
// 
// Start of implementation of the class TPCFindGenTracks
//
////////////////////////////////////////////////////////////////////////

class TPCFindGenTracks {

public:
  TPCFindGenTracks();
  TPCFindGenTracks(char* fnHits,
		   char* fnDigits ="tpc.digits.root",
		   char* fnRes    ="genTracks.root",
		   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~TPCFindGenTracks();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeGenTracks();
  void CloseOutputFile();
  Int_t TreeKLoop();
  Int_t TreeTRLoop();
  Int_t TreeDLoop();
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}

  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);

public:
  Int_t fDebug;                   //!< debug flag
  Int_t fEventNr;                 //!< current event number
  Int_t fLabel;                   //!< track label
  Int_t fNEvents;                 //!< number of events to process
  Int_t fFirstEventNr;            //!< first event to process
  Int_t fNParticles;              //!< number of particles in TreeK
  TTree *fTreeGenTracks;          //!< output tree with generated tracks
  char *fFnRes;                   //!< output file name with stored tracks
  char *fFnHits;                  //!< input file name with hits
  char *fFnDigits;                //!< input file name with digits
  TFile *fFileGenTracks;             //!< output file with stored fTreeGenTracks
  TFile *fFileHits;               //!< input file with hits
  TFile *fFileTreeD;              //!< input file with digits
  digitRow *fDigitRow;            //!< pointer to the object saved in Branch
  digitRow *fContainerDigitRow;   //!< big container for partial information
  AliTrackReference *fTrackRef;   //!< track reference saved in the output tree
  AliTrackReference *fContainerTR;//!< big container for partial information
  Int_t *fIndexTR;                //!< index of particle label in the fContainerTR
  Int_t fLastIndexTR;             //!< last used index in fIndexTR

  AliTPCParam* fParamTPC;         //!< AliTPCParam

  Double_t fVPrim[3];             //!< primary vertex position
  Double_t fVDist[4];             //!< distance of the particle vertex from primary vertex
                                  // the fVDist[3] contains size of the 3-vector
  TParticle *fParticle;           //!< generated particle

  Int_t fRowsWithDigitsInn;       //!< number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;          //!< number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;         //!< last - first row with digit
  Int_t fDigitsInSeed;            //!< digits in the default seed rows

private:

// some constants for the original non-pareller tracking (by Y.Belikov)
  static const Int_t seedRow11 = 158;  ///< nRowUp - 1
  static const Int_t seedRow12 = 139;  ///< nRowUp - 1 - (Int_t) 0.125*nRowUp
  static const Int_t seedRow21 = 149;  ///< seedRow11 - shift
  static const Int_t seedRow22 = 130;  ///< seedRow12 - shift
  static const Double_t kRaddeg = 180./TMath::Pi();

  static const Int_t fgMaxIndexTR = 50000; ///< maximum number of tracks with a track ref
  static const Int_t fgMaxParticles = 2000000; ///< maximum number of generated particles
  static const Double_t fgPtCut = .001; ///< do not store particles with generated pT less than this
  static const Float_t fgTrackRefLocalXMax = 82.95;
  static const Float_t fgTrackRefLocalXMaxDelta = 500.; 

  ClassDef(TPCFindGenTracks,1)    // class which creates and fills tree with TPCGenTrack objects
};
ClassImp(TPCFindGenTracks)
  
////////////////////////////////////////////////////////////////////////
TPCFindGenTracks::TPCFindGenTracks()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
TPCFindGenTracks::TPCFindGenTracks(char* fnHits, char* fnDigits, 
				   char* fnRes,
				   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  fFnRes = fnRes;
  fFnHits = fnHits;
  fFnDigits = fnDigits;  
}
////////////////////////////////////////////////////////////////////////
void TPCFindGenTracks::Reset()
{
  fDigitRow = 0;
  fEventNr = 0;
  fNEvents = 0;
  fTreeGenTracks = 0;
  fFnRes = "genTracks.root";
  fFnHits = "rfio:galice.root";
  fFnDigits = "rfio:its.tpc.trd.digits.root";
  fFileGenTracks = 0;
  fFileHits =0;
  fFileTreeD =0;
  fContainerDigitRow = 0;
  fContainerTR = 0;
  fIndexTR = 0;
  fLastIndexTR = -1;
  fParticle = 0;
  fTrackRef = 0;
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

  fDigitRow = new digitRow();
  CreateTreeGenTracks();
  if (!fTreeGenTracks) return 1;
  fFileHits = OpenAliceFile(fFnHits,kTRUE,"read"); //gAlice is read here
  if (!fFileHits) return 1;  
  fFileHits->cd();
  SetFieldFactor(); 

  fFileTreeD = TFile::Open(fFnDigits,"read");
  if (!fFileTreeD->IsOpen()) {
    cerr<<"Cannot open file "<<fFnDigits<<endl;
    return 1;
  }

  fParamTPC = LoadTPCParam(fFileTreeD);
  if (!fParamTPC) {
    cerr<<"TPC parameters not found and could not be created"<<endl;
    return 1;
  }

  for (fEventNr = fFirstEventNr; fEventNr < fFirstEventNr+fNEvents;
       fEventNr++) {
    fNParticles = gAlice->GetEvent(fEventNr);
    fContainerDigitRow = new digitRow[fNParticles];
    fContainerTR = new AliTrackReference[fgMaxIndexTR];
    fIndexTR = new Int_t[fNParticles];
    for (Int_t i = 0; i<fNParticles; i++) {
      fIndexTR[i] = -1;
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
    delete [] fContainerTR;
    delete [] fIndexTR;
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  fFileHits->cd();
  delete gAlice;
  fFileHits->Close();
  delete fFileHits;

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
  TBranch *branchBits = fTreeGenTracks->Branch("bitsBranch", "digitRow", &fDigitRow, 4000, 0);
  if (!branchBits) {
    cerr<<"Error in CreateTreeGenTracks: cannot create branch."<<endl;
    return;
  }
  fTreeGenTracks->Branch("fEventNr",&fEventNr,"fEventNr/I");
  fTreeGenTracks->Branch("fLabel",&fLabel,"fLabel/I");
  fTreeGenTracks->Branch("fRowsWithDigitsInn",&fRowsWithDigitsInn,"fRowsWithDigitsInn/I");
  fTreeGenTracks->Branch("fRowsWithDigits",&fRowsWithDigits,"fRowsWithDigits/I");
  fTreeGenTracks->Branch("fRowsTrackLength",&fRowsTrackLength,"fRowsTrackLength/I");
  fTreeGenTracks->Branch("fDigitsInSeed",&fDigitsInSeed,"fDigitsInSeed/I");

  fTreeGenTracks->Branch("Particle","TParticle",&fParticle);
  fTreeGenTracks->Branch("fVDist",&fVDist,"fVDist[4]/D");
  fTreeGenTracks->Branch("TR","AliTrackReference",&fTrackRef);

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
/// open the file with treeK
/// loop over all entries there and save information about some tracks

  fFileHits->cd();
  if (fDebug > 0) {
    cout<<"There are "<<fNParticles<<" primary and secondary particles in event "
	<<fEventNr<<endl;
  }
  AliStack * stack = gAlice->Stack();
  if (!stack) {cerr<<"Stack was not found!\n"; return 1;}

// not all generators give primary vertex position. Take the vertex
// of the particle 0 as primary vertex.
  fParticle = stack->ParticleFromTreeK(0);
  fVPrim[0] = fParticle->Vx();
  fVPrim[1] = fParticle->Vy();
  fVPrim[2] = fParticle->Vz();

  for (Int_t iParticle = 0; iParticle < fNParticles; iParticle++) {
//  for (Int_t iParticle = 0; iParticle < fDebug; iParticle++) {

// load only particles with TR
    if (fIndexTR[iParticle] < 0) continue;
    fParticle = stack->ParticleFromTreeK(iParticle);
    if (fDebug > 3 && iParticle < 10) {
      cout<<"processing particle "<<iParticle<<" ";
      fParticle->Print();
    }

// fill the tree

    fLabel = iParticle;
    fVDist[0] = fParticle->Vx()-fVPrim[0];
    fVDist[1] = fParticle->Vy()-fVPrim[1];
    fVDist[2] = fParticle->Vz()-fVPrim[2];
    fVDist[3] = TMath::Sqrt(fVDist[0]*fVDist[0]+fVDist[1]*fVDist[1]+fVDist[2]*fVDist[2]);
    fDigitRow = &(fContainerDigitRow[iParticle]);
    fRowsWithDigitsInn = fDigitRow->RowsOn(63); // 63 = number of inner rows
    fRowsWithDigits = fDigitRow->RowsOn();    
    fRowsTrackLength = fDigitRow->Last() - fDigitRow->First();
    if (fDebug > 2 && iParticle < 10) {
      cerr<<"Fill track with a label "<<iParticle<<endl;
    }
    fDigitsInSeed = 0;
    if (fDigitRow->TestRow(seedRow11) && fDigitRow->TestRow(seedRow12)) 
      fDigitsInSeed = 1;
    if (fDigitRow->TestRow(seedRow21) && fDigitRow->TestRow(seedRow22)) 
      fDigitsInSeed += 10;

    if (fIndexTR[iParticle] >= 0) {
      fTrackRef = &(fContainerTR[fIndexTR[iParticle]]);
//      cerr<<"Debug: fTrackRef->X() = "<<fTrackRef->X()<<endl;
    } else {
      fTrackRef->SetTrack(-1);
    }

    fTreeGenTracks->Fill();

  }
  fTreeGenTracks->AutoSave();
//  delete gAlice; gAlice = 0;
//  fFileHits->Close();

  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindGenTracks::TreeDLoop()
{
/// open the file with treeD
/// loop over all entries there and save information about some tracks


//  Int_t nrow_up=fParamTPC->GetNRowUp();
//  Int_t nrows=fParamTPC->GetNRowLow()+nrow_up;
  Int_t nInnerSector = fParamTPC->GetNInnerSector();
  Int_t rowShift = 0;
  Int_t zero=fParamTPC->GetZeroSup();
//  Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
  
  char treeDName[100]; 
  sprintf(treeDName,"TreeD_75x40_100x60_150x60_%d",fEventNr);
  TTree *treeD=(TTree*)fFileTreeD->Get(treeDName);
  AliSimDigits digitsAddress, *digits=&digitsAddress;
  treeD->GetBranch("Segment")->SetAddress(&digits);

  Int_t sectorsByRows=(Int_t)treeD->GetEntries();
  if (fDebug > 1) cout<<"\tsectorsByRows = "<<sectorsByRows<<endl;
  for (Int_t i=0; i<sectorsByRows; i++) {
//  for (Int_t i=5720; i<sectorsByRows; i++) {
    if (!treeD->GetEvent(i)) continue;
    Int_t sec,row;
    fParamTPC->AdjustSectorRow(digits->GetID(),sec,row);
    if (fDebug > 1) cout<<sec<<' '<<row<<"                          \r";
//    cerr<<sec<<' '<<row<<endl;

// here I expect that upper sectors follow lower sectors
    if (sec > nInnerSector) rowShift = fParamTPC->GetNRowLow();
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
	  if (label >= fNParticles) {
	    cerr<<"particle label too big: fNParticles, label "
		<<fNParticles<<" "<<label<<endl;
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
/// loop over TrackReferences and store the first one for each track

  TTree *treeTR=gAlice->TreeTR();
  if (!treeTR) {
    cerr<<"TreeTR not found"<<endl;
    return 1;
  }
  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  if (fDebug > 1) cout<<"There are "<<nPrimaries<<" entries in TreeTR"<<endl;
  TBranch *TPCBranchTR  = treeTR->GetBranch("TPC");
  if (!TPCBranchTR) {
    cerr<<"TPC branch in TR not found"<<endl;
    return 1;
  }
  TClonesArray* TPCArrayTR = new TClonesArray("AliTrackReference");
  TPCBranchTR->SetAddress(&TPCArrayTR);
  for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
    TPCBranchTR->GetEntry(iPrimPart);
    for (Int_t iTrackRef = 0; iTrackRef < TPCArrayTR->GetEntriesFast(); iTrackRef++) {
      AliTrackReference *trackRef = (AliTrackReference*)TPCArrayTR->At(iTrackRef);
      Int_t label = trackRef->GetTrack();
      
// save info in the fContainerTR
      if (label<0  || label > fNParticles) {
	cerr<<"Wrong label: "<<label<<endl;
	continue;
//	return 1;
      }

// store only if we do not have any track reference yet for this label
      if (fIndexTR[label] >= 0) continue;

// store only references with localX < fgTrackRefLocalXMax +- fgTrackRefLocalXMaxDelta
// and the pT > fgPtCut
      Float_t localX = TR2LocalX(trackRef,fParamTPC);
      if (TMath::Abs(localX-fgTrackRefLocalXMax)>fgTrackRefLocalXMaxDelta) continue;
      if (trackRef->Pt() < fgPtCut) continue;


//      cout<<"label, xg "<<label<<" "<<xg<<endl;
      if (fLastIndexTR >= fgMaxIndexTR-1) {
	cerr<<"Too many tracks with track reference. Increase the constant"
	    <<" fgMaxIndexTR"<<endl;
	return 1;
      }
      fIndexTR[label] =  ++fLastIndexTR;

// someone was lazy to implement copy ctor in AliTrackReference, so we have to do
//  it here "manually"
      fContainerTR[fIndexTR[label]].SetPosition(trackRef->X(),trackRef->Y(),trackRef->Z());

      fContainerTR[fIndexTR[label]].SetMomentum(trackRef->Px(),trackRef->Py(),trackRef->Pz());
      fContainerTR[fIndexTR[label]].SetTrack(trackRef->GetTrack());
      fContainerTR[fIndexTR[label]].SetLength(trackRef->GetLength());
//      cerr<<"Debug: trackRef->X(), stored: "<<trackRef->X()<<" "
//	  << fContainerTR[fIndexTR[label]].X()<<endl;

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
// 
// Start of implementation of the class TPCCmpTr
//
////////////////////////////////////////////////////////////////////////

class TPCCmpTr {

public:
  TPCCmpTr();
  TPCCmpTr(char* fnRecTracks,
	   char* fnGenTracks   ="genTracks.root",
	   char* fnCmpRes      ="cmpTracks.root", 
	   char* fnGalice      ="galice.root",
	   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~TPCCmpTr();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeCmp();
  void CloseOutputFile();
  Bool_t ConnectGenTree();
  Int_t TreeGenLoop(Int_t eventNr);
  Int_t TreeTLoop(Int_t eventNr);
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}

// tmp method, should go to TrackReferenceTPC
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);

private:
  digitRow *fDigitRow;            //!< pointer to the object saved in Branch
  Int_t fEventNr;                 //!< current event number
  Int_t fLabel;                   //!< track label
  Int_t fNEvents;                 //!< number of events to process
  Int_t fFirstEventNr;            //!< first event to process

  char *fFnCmp;                   //!< output file name with cmp tracks
  TFile *fFileCmp;                //!< output file with cmp tracks
  TTree *fTreeCmp;                //!< output tree with cmp tracks

  char *fFnGenTracks;             //!< input file name with gen tracks
  TFile *fFileGenTracks;
  TTree *fTreeGenTracks;

  char *fFnHits;                  //!< input file name with gAlice object (needed for B)
  TFile *fFileHits;               //!< input file with gAlice

  char *fFnRecTracks;             //!< input file name with tpc rec. tracks
  TFile *fFileRecTracks;          //!< input file with reconstructed tracks
  TTree *fTreeRecTracks;          //!< tree with reconstructed tracks

  AliTPCtrack *fTPCTrack;         //!< pointer to TPC track to connect branch
  Int_t *fIndexRecTracks;         //!< index of particle label in the TreeT_TPC

  Int_t fRowsWithDigitsInn;       //!< number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;          //!< number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;         //!< last - first row with digit
  Int_t fDigitsInSeed;            //!< digits in the default seed rows
  TParticle *fParticle;           //!< generated particle
  Double_t fVDist[4];             //!< distance of the particle vertex from primary vertex
                                  // the fVDist[3] contains size of the 3-vector
  AliTrackReference *fTrackRef;   //!< track reference saved in the output tree
  Int_t   fReconstructed;         //!< flag if track was reconstructed
  AliTPCParam* fParamTPC;         //!< AliTPCParam

  Int_t fNParticles;              //!< number of particles in the input tree genTracks
  Int_t fDebug;                   //!< debug flag

  Int_t fNextTreeGenEntryToRead;    //!< last entry already read from genTracks tree
  TPCGenInfo *fGenInfo;           //!< container for all the details

  Double_t fRecPhi;         ///< reconstructed phi angle (0;2*kPI)
  Double_t fLambda;         ///< reconstructed
  Double_t fRecPt_1;        ///< reconstructed
  Float_t fdEdx;           ///< reconstructed  dEdx


  ClassDef(TPCCmpTr,1)    // class which creates and fills tree with TPCGenTrack objects
};
/// \cond CLASSIMP
ClassImp(TPCCmpTr)
/// \endcond
  
////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr(char* fnRecTracks,
		   char* fnGenTracks,
		   char* fnCmp,
		   char* fnGalice,
		   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  fFnRecTracks = fnRecTracks;
  fFnGenTracks = fnGenTracks;
  fFnCmp = fnCmp;
  fFnHits = fnGalice;
  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
}
////////////////////////////////////////////////////////////////////////
void TPCCmpTr::Reset()
{
  fDigitRow = 0;
  fEventNr = 0;
  fNEvents = 0;
  fTreeCmp = 0;
  fFnCmp = "cmpTracks.root";
  fFnHits = "galice.root";
  fFileGenTracks = 0;
  fFileHits =0;
  fParticle = 0;
  fTrackRef = 0;

  fRowsWithDigitsInn = 0;
  fRowsWithDigits = 0;
  fRowsTrackLength = 0;
  fDigitsInSeed = 0;

  fDebug = 0;

  fParamTPC = 0;
  fFnRecTracks = "tpc.tracks.root";
  fTreeRecTracks = 0;
  fFileRecTracks = 0;
  fTPCTrack = 0; 
  fGenInfo = 0; 
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

  fDigitRow = new digitRow();
  CreateTreeCmp();
  if (!fTreeCmp) {
    cerr<<"output tree not created"<<endl;
    return 1;
  }
  fFileHits = OpenAliceFile(fFnHits,kTRUE,"read"); //gAlice is read here
  if (!fFileHits) {
    cerr<<"Cannot open file with gAlice object "<<fFnHits<<endl;
    return 1;  
  }
  fFileHits->cd();
  SetFieldFactor(); 

  fParamTPC = LoadTPCParam(fFileHits);
  if (!fParamTPC) {
    cerr<<"TPC parameters not found and could not be created"<<endl;
    return 1;
  }

  if (!ConnectGenTree()) {
    cerr<<"Cannot connect tree with generated tracks"<<endl;
    return 1;
  }
  fFileHits->cd();
  fNextTreeGenEntryToRead = 0;
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    fNParticles = gAlice->GetEvent(fEventNr);    
    fIndexRecTracks = new Int_t[fNParticles];
    for (Int_t i = 0; i<fNParticles; i++) {
      fIndexRecTracks[i] = -1;
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
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  fFileHits->cd();
  delete gAlice;
  fFileHits->Close();
  delete fFileHits;

  timer.Stop();
  timer.Print();
  return 0;
}
////////////////////////////////////////////////////////////////////////
Bool_t TPCCmpTr::ConnectGenTree()
{
/// connect all branches from the genTracksTree
/// use the same variables as for the new cmp tree, it may work

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
  fTreeGenTracks->SetBranchAddress("fEventNr",&fEventNr);
  fTreeGenTracks->SetBranchAddress("fLabel",&fLabel);
  fTreeGenTracks->SetBranchAddress("fRowsWithDigitsInn",&fRowsWithDigitsInn);
  fTreeGenTracks->SetBranchAddress("fRowsWithDigits",&fRowsWithDigits);
  fTreeGenTracks->SetBranchAddress("fRowsTrackLength",&fRowsTrackLength);
  fTreeGenTracks->SetBranchAddress("fDigitsInSeed",&fDigitsInSeed);
  fTreeGenTracks->SetBranchAddress("Particle",&fParticle);
  fTreeGenTracks->SetBranchAddress("fVDist",fVDist);
  fTreeGenTracks->SetBranchAddress("TR",&fTrackRef);

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
  fTreeCmp = new TTree("TPCcmpTracks","TPCcmpTracks");
  TBranch *branchBits = fTreeCmp->Branch("bitsBranch", "digitRow", &fDigitRow, 4000, 0);
  if (!branchBits) {
    cerr<<"Error in CreateTreeCmp: cannot create branch."<<endl;
    return;
  }
  fTreeCmp->Branch("fEventNr",&fEventNr,"fEventNr/I");
  fTreeCmp->Branch("fLabel",&fLabel,"fLabel/I");
  fTreeCmp->Branch("fRowsWithDigitsInn",&fRowsWithDigitsInn,"fRowsWithDigitsInn/I");
  fTreeCmp->Branch("fRowsWithDigits",&fRowsWithDigits,"fRowsWithDigits/I");
  fTreeCmp->Branch("fRowsTrackLength",&fRowsTrackLength,"fRowsTrackLength/I");
  fTreeCmp->Branch("fDigitsInSeed",&fDigitsInSeed,"fDigitsInSeed/I");

  fTreeCmp->Branch("fReconstructed",&fReconstructed,"fReconstructed/I");
  fTreeCmp->Branch("fTPCTrack","AliTPCtrack",&fTPCTrack);

  fTreeCmp->Branch("Particle","TParticle",&fParticle);
  fTreeCmp->Branch("fVDist",&fVDist,"fVDist[4]/D");
  fTreeCmp->Branch("TR","AliTrackReference",&fTrackRef);

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

Float_t TPCCmpTr::TR2LocalX(AliTrackReference *trackRef,
			    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return x[0];
}
////////////////////////////////////////////////////////////////////////

Int_t TPCCmpTr::TreeTLoop(Int_t eventNr)
{
/// loop over all TPC reconstructed tracks and store info in memory

  
  if (!fFileRecTracks) fFileRecTracks = TFile::Open(fFnRecTracks,"read");
  if (!fFileRecTracks->IsOpen()) {
    cerr<<"Cannot open file "<<fFnRecTracks<<endl;
    return 1;
  }

  char treeNameBase[11] = "TreeT_TPC_";
  char treeName[20];
  sprintf(treeName,"%s%d",treeNameBase,eventNr);

  fTreeRecTracks=(TTree*)fFileRecTracks->Get(treeName);
  if (!fTreeRecTracks) {
    cerr<<"Can't get a tree with TPC rec. tracks named "<<treeName<<endl;
    return 1;
  }
  
  Int_t nEntries = (Int_t) fTreeRecTracks->GetEntries();
  if (fDebug > 2) cout<<"Event, rec. tracks: "<<eventNr<<" "
		      <<nEntries<<endl;
  TBranch * br= fTreeRecTracks->GetBranch("tracks");
  br->SetAddress(&fTPCTrack);

  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    br->GetEntry(iEntry);
    Int_t label = fTPCTrack->GetLabel();
    fIndexRecTracks[label] =  iEntry; 
  }  

  if (fDebug > 2) cerr<<"end of TreeTLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::TreeGenLoop(Int_t eventNr)
{
/// loop over all entries for a given event, find corresponding
/// rec. track and store in the fTreeCmp

  fFileGenTracks->cd();
  Int_t entry = fNextTreeGenEntryToRead;
  Double_t nParticlesTR = fTreeGenTracks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextTreeGenEntryToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextTreeGenEntryToRead<<endl;
  while (entry < nParticlesTR) {
    fTreeGenTracks->GetEntry(entry);
    entry++;
    if (fEventNr < eventNr) continue;
    if (fEventNr > eventNr) break;
    fNextTreeGenEntryToRead = entry-1;
    if (fDebug > 2 && fLabel < 10) {
      cerr<<"Fill track with a label "<<fLabel<<endl;
    }

    fReconstructed = 0;
    fdEdx = 0.;
    fRecPhi = fLambda = fRecPt_1 = 0.;

    if (fDebug > 2) {
      cerr<<"fLabel, fIndexRecTracks[fLabel] "<<fLabel<<" "<<fIndexRecTracks[fLabel]<<endl;
    }
    if (fIndexRecTracks[fLabel] >= 0) {
      Int_t nBytes = fTreeRecTracks->GetEvent(fIndexRecTracks[fLabel]);
      if (nBytes > 0) {
	fReconstructed = 1;
	fdEdx = fTPCTrack->GetdEdx();
	Double_t Localx = TR2LocalX(fTrackRef,fParamTPC);
	if (fDebug > 3) {
	  cerr<<"Track local X before prolongation: "<<fTPCTrack->GetX()<<endl;
	}
	fTPCTrack->PropagateTo(Localx);
	Double_t par[5];
	if (fDebug > 3) {
	  cerr<<"Track local X after prolongation: "<<fTPCTrack->GetX()<<endl;
	}
	fTPCTrack->GetExternalParameters(Localx,par);
	fRecPhi=TMath::ASin(par[2]) + fTPCTrack->GetAlpha();
	if (fRecPhi<0) fRecPhi+=2*TMath::Pi();
	if (fRecPhi>=2*TMath::Pi()) fRecPhi-=2*TMath::Pi();
//	  fRecPhi = (fRecPhi)*kRaddeg;
	fLambda = TMath::ATan(par[3]);
	fRecPt_1 = TMath::Abs(par[4]);
      }
    }


    fTreeCmp->Fill();

  }
  fTreeCmp->AutoSave();
//  delete gAlice; gAlice = 0;
//  fFileHits->Close();

  if (fDebug > 2) cerr<<"end of TreeKLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
