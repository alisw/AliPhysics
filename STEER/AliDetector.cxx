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

/* $Id$ */

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

#include <assert.h>

#include <Riostream.h>
#include <TBrowser.h>
#include <TFile.h>
#include <TFolder.h>
#include <TROOT.h>
#include <TTree.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliHit.h"
#include "AliPoints.h"
#include "AliRun.h"
#include "AliTrackReference.h"


// Static variables for the hit iterator routines
static Int_t sMaxIterHit=0;
static Int_t sCurIterHit=0;


ClassImp(AliDetector)
 
//_______________________________________________________________________
AliDetector::AliDetector():
  fTimeGate(200.e-9),
  fIshunt(0),
  fNhits(0),
  fNdigits(0),
  fBufferSize(1600),
  fHits(0),
  fDigits(0),
  fDigitsFile(0),
  fPoints(0),
  fTrackReferences(0),
  fMaxIterTrackRef(0),
  fCurrentIterTrackRef(0)
{
  //
  // Default constructor for the AliDetector class
  //
}
 
//_______________________________________________________________________
AliDetector::AliDetector(const AliDetector &det):
  AliModule(det),
  fTimeGate(200.e-9),
  fIshunt(0),
  fNhits(0),
  fNdigits(0),
  fBufferSize(1600),
  fHits(0),
  fDigits(0),
  fDigitsFile(0),
  fPoints(0),
  fTrackReferences(0),
  fMaxIterTrackRef(0),
  fCurrentIterTrackRef(0)
{
  det.Copy(*this);
}

//_____________________________________________________________________________
AliDetector::AliDetector(const char* name,const char *title):
  AliModule(name,title),
  fTimeGate(200.e-9),
  fIshunt(0),
  fNhits(0),
  fNdigits(0),
  fBufferSize(1600),
  fHits(0),
  fDigits(0),
  fDigitsFile(0),
  fPoints(0),
  fTrackReferences(new TClonesArray("AliTrackReference", 100)),
  fMaxIterTrackRef(0),
  fCurrentIterTrackRef(0)
{
  //
  // Normal constructor invoked by all Detectors.
  // Create the list for detector specific histograms
  // Add this Detector to the global list of Detectors in Run.
  //

  fActive     = kTRUE;
  AliConfig::Instance()->Add(this);

}
 
//_______________________________________________________________________
AliDetector::~AliDetector()
{
  //
  // Destructor
  //

  // Delete space point structure
  if (fPoints) {
    fPoints->Delete();
    delete fPoints;
    fPoints     = 0;
  }
  // Delete digits structure
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits     = 0;
  }
  if (fDigitsFile) delete [] fDigitsFile;
}

//_______________________________________________________________________
void AliDetector::Publish(const char *dir, void *address, const char *name)
{
  //
  // Register pointer to detector objects. 
  // 
  TFolder *topFolder = dynamic_cast<TFolder *>(gROOT->FindObjectAny("/Folders"));
  if  (topFolder) { 
    TFolder *folder = dynamic_cast<TFolder *>(topFolder->FindObjectAny(dir));
    // TFolder *folder = dynamic_cast<TFolder *>(gROOT->FindObjectAny(dir));
    if (!folder)  {
      cerr << "Cannot register: Missing folder: " << dir << endl;
    } else {
      TFolder *subfolder = dynamic_cast<TFolder *>(folder->FindObjectAny(this->GetName())); 

      if(!subfolder)
         subfolder = folder->AddFolder(this->GetName(),this->GetTitle());
      if (address) {
        TObject **obj = static_cast<TObject **>(address);
        if ((*obj)->InheritsFrom(TCollection::Class())) {
           TCollection *collection = dynamic_cast<TCollection *>(*obj); 
           if (name)
             collection->SetName(name);
        } 
        subfolder->Add(*obj);
      } 
    }  
  }
}

//_______________________________________________________________________
TBranch* AliDetector::MakeBranchInTree(TTree *tree, const char* name, 
                                       void* address, Int_t size,
                                       const char *file)
{ 
    return(MakeBranchInTree(tree,name,0,address,size,99,file));
}

//_______________________________________________________________________
TBranch* AliDetector::MakeBranchInTree(TTree *tree, const char* name, 
                                       const char *classname, 
                                       void* address,Int_t size, 
                                       Int_t splitlevel, const char *file)
{ 
    //
    // Makes branch in given tree and diverts them to a separate file
    //  
    if (GetDebug()>1)
      printf("* MakeBranch * Making Branch %s \n",name);
      
    TDirectory *cwd = gDirectory;
    TBranch *branch = 0;
    
    if (classname) {
      branch = tree->Branch(name,classname,address,size,splitlevel);
    } else {
      branch = tree->Branch(name,address,size);
    }
       
    if (file) {
        char * outFile = new char[strlen(gAlice->GetBaseFile())+strlen(file)+2];
        sprintf(outFile,"%s/%s",gAlice->GetBaseFile(),file);
        branch->SetFile(outFile);
        TIter next( branch->GetListOfBranches());
        while ((branch=dynamic_cast<TBranch*>(next()))) {
           branch->SetFile(outFile);
        } 
       delete outFile;
        
       cwd->cd();
        
       if (GetDebug()>1)
           printf("* MakeBranch * Diverting Branch %s to file %s\n",name,file);
    }
    const char *folder = 0;
    TString folderName(name);  
    
    if (!strncmp(tree->GetName(),"TreeE",5)) folder = "RunMC/Event/Data";
    if (!strncmp(tree->GetName(),"TreeK",5)) folder = "RunMC/Event/Data";
    if (!strncmp(tree->GetName(),"TreeH",5)) {
      folder     = "RunMC/Event/Data/Hits";
      folderName = "Hits" ; 
    }
    if (!strncmp(tree->GetName(),"TreeTrackReferences",5)) {
      folder     = "RunMC/Event/Data/TrackReferences";
      folderName = "TrackReferences" ; 
    }

    if (!strncmp(tree->GetName(),"TreeD",5)) {
      folder     = "Run/Event/Data";
      folderName = "Digits" ; 
    }
    if (!strncmp(tree->GetName(),"TreeS",5)) {
      folder     = "RunMC/Event/Data/SDigits";
      folderName = "SDigits" ; 
    }
    if (!strncmp(tree->GetName(),"TreeR",5)) folder = "Run/Event/RecData";

    if (folder) {
      if (GetDebug())
          printf("%15s: Publishing %s to %s\n",ClassName(),name,folder);
      Publish(folder,address, folderName.Data());
    }  
    return branch;
}

//_______________________________________________________________________
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

//_______________________________________________________________________
void AliDetector::Copy(AliDetector &) const
{
  //
  // Copy *this onto det -- not implemented
  //
  Fatal("Copy","Not implemented\n");
}

//_______________________________________________________________________
void AliDetector::FinishRun()
{
  //
  // Procedure called at the end of a run.
  //
}

//_______________________________________________________________________
void AliDetector::RemapTrackReferencesIDs(Int_t *map)
{
  // 
  // Remapping track reference
  // Called at finish primary
  //
  if (!fTrackReferences) return;
  for (Int_t i=0;i<fTrackReferences->GetEntries();i++){
    AliTrackReference * ref = dynamic_cast<AliTrackReference*>(fTrackReferences->UncheckedAt(i));
    if (ref) {
      Int_t newID = map[ref->GetTrack()];
      if (newID>=0) ref->SetTrack(newID);
      else ref->SetTrack(-1);
      
    }
  }
}

//_______________________________________________________________________
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
  if(sMaxIterHit) return dynamic_cast<AliHit*>(fHits->UncheckedAt(0));
  else            return 0;
}


//_______________________________________________________________________
AliTrackReference* AliDetector::FirstTrackReference(Int_t track)
{
  //
  // Initialise the hit iterator
  // Return the address of the first hit for track
  // If track>=0 the track is read from disk
  // while if track<0 the first hit of the current
  // track is returned
  // 
  if(track>=0) {
    gAlice->ResetTrackReferences();
    gAlice->TreeTR()->GetEvent(track);
  }
  //
  fMaxIterTrackRef     = fTrackReferences->GetEntriesFast();
  fCurrentIterTrackRef = 0;
  if(fMaxIterTrackRef) return dynamic_cast<AliTrackReference*>(fTrackReferences->UncheckedAt(0));
  else            return 0;
}

//_______________________________________________________________________
AliHit* AliDetector::NextHit()
{
  //
  // Return the next hit for the current track
  //
  if(sMaxIterHit) {
    if(++sCurIterHit<sMaxIterHit) 
      return dynamic_cast<AliHit*>(fHits->UncheckedAt(sCurIterHit));
    else        
      return 0;
  } else {
    printf("* AliDetector::NextHit * Hit Iterator called without calling FistHit before\n");
    return 0;
  }
}

//_______________________________________________________________________
AliTrackReference* AliDetector::NextTrackReference()
{
  //
  // Return the next hit for the current track
  //
  if(fMaxIterTrackRef) {
    if(++fCurrentIterTrackRef<fMaxIterTrackRef) 
      return dynamic_cast<AliTrackReference*>(fTrackReferences->UncheckedAt(fCurrentIterTrackRef));
    else        
      return 0;
  } else {
    printf("* AliDetector::NextTrackReference * TrackReference  Iterator called without calling FistTrackReference before\n");
    return 0;
  }
}

//_______________________________________________________________________
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
    ahit = dynamic_cast<AliHit*>(fHits->UncheckedAt(hit));
    trk=ahit->GetTrack();
    assert(trk<=tracks);
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
    fp[3*ntrk[trk]  ] = ahit->X();
    fp[3*ntrk[trk]+1] = ahit->Y();
    fp[3*ntrk[trk]+2] = ahit->Z();
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

//_______________________________________________________________________
void AliDetector::MakeBranch(Option_t *option, const char *file)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  //
 
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  //
  // Get the pointer to the header
  const char *cH = strstr(option,"H");
  //
  if (fHits && gAlice->TreeH() && cH) {
    MakeBranchInTree(gAlice->TreeH(), 
                     branchname, &fHits, fBufferSize, file) ;              
  }	
  
  const char *cD = strstr(option,"D");

  if (cD) {
    if (file) {
       fDigitsFile = new char[strlen (file)];
       strcpy(fDigitsFile,file);
    }
  }
}
//_______________________________________________________________________
void AliDetector::MakeBranchTR(Option_t *option, const char *file)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  //
 
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  //
  // Get the pointer to the header
  const char *cTR = strstr(option,"T");
  //
  if (fTrackReferences && gAlice->TreeTR() && cTR) {
    MakeBranchInTree(gAlice->TreeTR(), 
                     branchname, &fTrackReferences, fBufferSize, file) ;              
  }	  
}

//_______________________________________________________________________
void AliDetector::ResetDigits()
{
  //
  // Reset number of digits and the digits array
  //
  fNdigits   = 0;
  if (fDigits)   fDigits->Clear();
}

//_______________________________________________________________________
void AliDetector::ResetHits()
{
  //
  // Reset number of hits and the hits array
  //
  fNhits   = 0;
  if (fHits)   fHits->Clear();
}

//_______________________________________________________________________
void AliDetector::ResetTrackReferences()
{
  //
  // Reset number of hits and the hits array
  //
  fMaxIterTrackRef   = 0;
  if (fTrackReferences)   fTrackReferences->Clear();
}

//_______________________________________________________________________
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

//_______________________________________________________________________
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

  // Branch address for tr  tree
  TTree *treeTR = gAlice->TreeTR();
  if (treeTR && fTrackReferences) {
    branch = treeTR->GetBranch(branchname);
    if (branch) branch->SetAddress(&fTrackReferences);
  }
}

 
