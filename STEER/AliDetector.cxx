///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base class for ALICE modules. Both sensitive modules (detectors) and      //
// non-sensitive ones are described by this base class. This class           //
// supports the hit and digit trees produced by the simulation and also      //
// the objects produced by the reconstruction.                               //
//                                                                           //
// This class is also responsible for building the geometry of the           //
// detectors.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliDetectorClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliRun.h"
#include "AliHit.h"
#include "AliPoints.h"
#include <TClass.h>
#include <TNode.h>
#include <TRandom.h>

// Static variables for the hit iterator routines
static Int_t sMaxIterHit=0;
static Int_t sCurIterHit=0;

ClassImp(AliDetector)
 
//_____________________________________________________________________________
AliDetector::AliDetector()
{
  //
  // Default constructor for the AliDetector class
  //
  fNhits      = 0;
  fNdigits    = 0;
  fPoints     = 0;
  fHits       = 0;
  fDigits     = 0;
  fTimeGate   = 200.e-9;
  fBufferSize = 16000;
}
 
//_____________________________________________________________________________
AliDetector::AliDetector(const char* name,const char *title):AliModule(name,title)
{
  //
  // Normal constructor invoked by all Detectors.
  // Create the list for detector specific histograms
  // Add this Detector to the global list of Detectors in Run.
  //

  fTimeGate   = 200.e-9;
  fActive     = kTRUE;
  fNhits      = 0;
  fHits       = 0;
  fDigits     = 0;
  fNdigits    = 0;
  fPoints     = 0;
  fBufferSize = 16000;
}
 
//_____________________________________________________________________________
AliDetector::~AliDetector()
{
  //
  // Destructor
  //
  fNhits      = 0;
  fNdigits    = 0;
  //
  // Delete space point structure
  if (fPoints) fPoints->Delete();
  delete fPoints;
  fPoints     = 0;
}
 
//_____________________________________________________________________________
void AliDetector::Browse(TBrowser *b)
{
  //
  // Insert Detector objects in the list of objects to be browsed
  //
  char name[64];
  if( fHits == 0) return;
  TObject *obj;
  Int_t i, nobjects;
  //
  nobjects = fHits->GetEntries();
  for (i=0;i<nobjects;i++) {
    obj = fHits->At(i);
    sprintf(name,"%s_%d",obj->GetName(),i);
    b->Add(obj, &name[0]);
  }
}

//_____________________________________________________________________________
void AliDetector::FinishRun()
{
  //
  // Procedure called at the end of a run.
  //
}

//_____________________________________________________________________________
AliHit* AliDetector::FirstHit(Int_t track)
{
  //
  // Initialise the hit iterator
  // Return the address of the first hit for track
  // If track>=0 the track is read from disk
  // while if track<0 the first hit of the current
  // track is returned
  // 
  if(track>=0) {
    gAlice->ResetHits();
    gAlice->TreeH()->GetEvent(track);
  }
  //
  sMaxIterHit=fHits->GetEntriesFast();
  sCurIterHit=0;
  if(sMaxIterHit) return (AliHit*) fHits->UncheckedAt(0);
  else            return 0;
}

//_____________________________________________________________________________
AliHit* AliDetector::NextHit()
{
  //
  // Return the next hit for the current track
  //
  if(sMaxIterHit) {
    if(++sCurIterHit<sMaxIterHit) 
      return (AliHit*) fHits->UncheckedAt(sCurIterHit);
    else        
      return 0;
  } else {
    printf("* AliDetector::NextHit * Hit Iterator called without calling FistHit before\n");
    return 0;
  }
}

//_____________________________________________________________________________
void AliDetector::LoadPoints(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  //
  if (fHits == 0) return;
  //
  if (fPoints == 0) fPoints = new TObjArray(gAlice->GetNtrack());
  Int_t nhits = fHits->GetEntriesFast();
  if (nhits == 0) return;
  AliHit *ahit;
  //
  AliPoints *points = 0;
  Int_t trko=-99, trk;
  //
  // Loop over all the hits and store their position
  for (Int_t hit=0;hit<nhits;hit++) {
    ahit = (AliHit*)fHits->UncheckedAt(hit);
    if(trko!=(trk=ahit->GetTrack())) {
      //
      // Initialise a new track
      trko=trk;
      points = new AliPoints(nhits);
      fPoints->AddAt(points,trk);
      points->SetMarkerColor(GetMarkerColor());
      points->SetMarkerStyle(GetMarkerStyle());
      points->SetMarkerSize(GetMarkerSize());
      points->SetDetector(this);
      points->SetParticle(trk);
    }
    points->SetPoint(hit,ahit->fX,ahit->fY,ahit->fZ);
  }
}

//_____________________________________________________________________________
void AliDetector::MakeBranch(Option_t *option)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  //
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  //
  // Get the pointer to the header
  char *H = strstr(option,"H");
  //
  if (fHits   && gAlice->TreeH() && H) {
    gAlice->TreeH()->Branch(branchname,&fHits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for hits\n",branchname);
  }	
}

//_____________________________________________________________________________
void AliDetector::ResetDigits()
{
  //
  // Reset number of digits and the digits array
  //
  fNdigits   = 0;
  if (fDigits)   fDigits->Clear();
}

//_____________________________________________________________________________
void AliDetector::ResetHits()
{
  //
  // Reset number of hits and the hits array
  //
  fNhits   = 0;
  if (fHits)   fHits->Clear();
}

//_____________________________________________________________________________
void AliDetector::ResetPoints()
{
  //
  // Reset array of points
  //
  if (fPoints) {
    fPoints->Delete();
    delete fPoints;
    fPoints = 0;
  }
}

//_____________________________________________________________________________
void AliDetector::SetTreeAddress()
{
  //
  // Set branch address for the Hits and Digits Trees
  //
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  //
  // Branch address for hit tree
  TTree *treeH = gAlice->TreeH();
  if (treeH && fHits) {
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fHits);
  }
  //
  // Branch address for digit tree
  TTree *treeD = gAlice->TreeD();
  if (treeD && fDigits) {
    branch = treeD->GetBranch(branchname);
    if (branch) branch->SetAddress(&fDigits);
  }
}

//_____________________________________________________________________________
void AliDetector::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class Detector.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    AliModule::Streamer(R__b);
    R__b >> fTimeGate;
    R__b >> fIshunt;
    //R__b >> fNhits;
    //
    // Stream the pointers but not the TClonesArrays
    R__b >> fHits; // diff
    R__b >> fDigits; // diff
    
  } else {
    R__b.WriteVersion(AliDetector::IsA());
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    AliModule::Streamer(R__b);
    R__b << fTimeGate;
    R__b << fIshunt;
    //R__b << fNhits;
    //
    // Stream the pointers but not the TClonesArrays
    R__b << fHits; // diff
    R__b << fDigits; // diff
  }
}
 
