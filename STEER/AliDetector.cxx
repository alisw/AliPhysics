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
<img src="picts/AliDetectorClass.gif">
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
  Int_t nhits = fHits->GetEntriesFast();
  if (nhits == 0) return;
  Int_t tracks = gAlice->GetNtrack();
  if (fPoints == 0) fPoints = new TObjArray(tracks);
  AliHit *ahit;
  //
  Int_t *ntrk=new Int_t[tracks];
  Int_t *limi=new Int_t[tracks];
  Float_t **coor=new Float_t*[tracks];
  for(Int_t i=0;i<tracks;i++) {
    ntrk[i]=0;
    coor[i]=0;
    limi[i]=0;
  }
  //
  AliPoints *points = 0;
  Float_t *fp=0;
  Int_t trk;
  Int_t chunk=nhits/4+1;
  //
  // Loop over all the hits and store their position
  for (Int_t hit=0;hit<nhits;hit++) {
    ahit = (AliHit*)fHits->UncheckedAt(hit);
    trk=ahit->GetTrack();
    if(ntrk[trk]==limi[trk]) {
      //
      // Initialise a new track
      fp=new Float_t[3*(limi[trk]+chunk)];
      if(coor[trk]) {
	memcpy(fp,coor[trk],sizeof(Float_t)*3*limi[trk]);
	delete [] coor[trk];
      }
      limi[trk]+=chunk;
      coor[trk] = fp;
    } else {
      fp = coor[trk];
    }
    fp[3*ntrk[trk]  ] = ahit->fX;
    fp[3*ntrk[trk]+1] = ahit->fY;
    fp[3*ntrk[trk]+2] = ahit->fZ;
    ntrk[trk]++;
  }
  //
  for(trk=0; trk<tracks; ++trk) {
    if(ntrk[trk]) {
      points = new AliPoints();
      points->SetMarkerColor(GetMarkerColor());
      points->SetMarkerSize(GetMarkerSize());
      points->SetDetector(this);
      points->SetParticle(trk);
      points->SetPolyMarker(ntrk[trk],coor[trk],GetMarkerStyle());
      fPoints->AddAt(points,trk);
      delete [] coor[trk];
      coor[trk]=0;
    }
  }
  delete [] coor;
  delete [] ntrk;
  delete [] limi;
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
 
