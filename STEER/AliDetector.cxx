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

#include <TBrowser.h>
#include <TTree.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliHit.h"
#include "AliPoints.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliMC.h"


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
  fPoints(0),
  fLoader(0x0)
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
  fPoints(0),
  fLoader(0x0)
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
  fPoints(0),
  fLoader(0x0)
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

  if (fLoader)
   {
    fLoader->GetModulesFolder()->Remove(this);
   }

}

//_______________________________________________________________________
void AliDetector::Publish(const char */*dir*/, void */*address*/, const char */*name*/) const
{
//
// Register pointer to detector objects. 
// 
  MayNotUse("Publish");
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
                                       Int_t splitlevel, const char */*file*/)
{ 
//
// Makes branch in given tree and diverts them to a separate file
// 
//
//
// if (GetDebug()>1)
 if(GetDebug()) Info("MakeBranch","Making Branch %s",name);
 if (tree == 0x0) 
  {
   Error("MakeBranch","Making Branch %s Tree is NULL",name);
   return 0x0;
  }
 TBranch *branch = tree->GetBranch(name);
 if (branch) 
  {  
    if(GetDebug()) Info("MakeBranch","Branch %s is already in tree.",name);
    return branch;
  }
    
 if (classname) 
  {
    branch = tree->Branch(name,classname,address,size,splitlevel);
  } 
 else 
  {
    branch = tree->Branch(name,address,size);
  }
 if(GetDebug()) Info("MakeBranch","Branch %s returning branch %#x",name,branch);
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
void AliDetector::Copy(TObject &) const
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
    gAlice->ResetHits(); //stupid = if N detector this method is called N times
    TreeH()->GetEvent(track); //skowron
  }
  //
  fMaxIterHit=fHits->GetEntriesFast();
  fCurIterHit=0;
  if(fMaxIterHit) return dynamic_cast<AliHit*>(fHits->UncheckedAt(0));
  else            return 0;
}

//_______________________________________________________________________
AliHit* AliDetector::NextHit()
{
  //
  // Return the next hit for the current track
  //
  if(fMaxIterHit) {
    if(++fCurIterHit<fMaxIterHit) 
      return dynamic_cast<AliHit*>(fHits->UncheckedAt(fCurIterHit));
    else        
      return 0;
  } else {
    printf("* AliDetector::NextHit * Hit Iterator called without calling FistHit before\n");
    return 0;
  }
}

//_______________________________________________________________________
void AliDetector::LoadPoints(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  //
  if (fHits == 0) 
   {
    Error("LoadPoints","fHits == 0. Name is %s",GetName());
    return;
   }
  //
  Int_t nhits = fHits->GetEntriesFast();
  if (nhits == 0) 
   {
//    Error("LoadPoints","nhits == 0. Name is %s",GetName());
    return;
   }
  Int_t tracks = gAlice->GetMCApp()->GetNtrack();
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
    if(trk>tracks) Fatal("LoadPoints","Found track number %d, max track %d\n",trk, tracks);
    if(ntrk[trk]==limi[trk])
     {
      //
      // Initialise a new track
      fp=new Float_t[3*(limi[trk]+chunk)];
      if(coor[trk]) 
       {
          memcpy(fp,coor[trk],sizeof(Float_t)*3*limi[trk]);
          delete [] coor[trk];
       }
      limi[trk]+=chunk;
      coor[trk] = fp;
     } 
    else 
     {
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
void AliDetector::MakeBranch(Option_t *option)
{
  //
  // Create a new branch for this detector in its treeH
  //

  if(GetDebug()) Info("MakeBranch"," for %s",GetName());
  const char *cH = strstr(option,"H");

  if (fHits && TreeH() && cH) 
   {
     MakeBranchInTree(TreeH(), GetName(), &fHits, fBufferSize, 0);
   }	
}

//_______________________________________________________________________
void AliDetector::ResetDigits()
{
  //
  // Reset number of digits and the digits array
  //
  fNdigits   = 0;
  if (fDigits) fDigits->Clear();
}

//_______________________________________________________________________
void AliDetector::ResetHits()
{
  //
  // Reset number of hits and the hits array
  //
  fNhits   = 0;
  if (fHits) fHits->Clear();
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
  //
  // Branch address for hit tree
  
  TTree *tree = TreeH();
  if (tree && fHits) {
    branch = tree->GetBranch(GetName());
    if (branch) 
     {
       if(GetDebug()) Info("SetTreeAddress","(%s) Setting for Hits",GetName());
       branch->SetAddress(&fHits);
     }
    else
     { //can be invoked before branch creation
       if(GetDebug()) Warning("SetTreeAddress","(%s) Failed for Hits. Can not find branch in tree.",GetName());
     }
  }
  
  //
  // Branch address for digit tree
  TTree *treeD = fLoader->TreeD();
  if (treeD && fDigits) {
    branch = treeD->GetBranch(GetName());
    if (branch) branch->SetAddress(&fDigits);
  }
  
  AliModule::SetTreeAddress();
}

//_______________________________________________________________________
void AliDetector::MakeTree(Option_t *option)
 {
 //makes a tree (container) for the data defined in option
 //"H" - hits
 //"D" - digits
 //"S" - summable digits
 //"R" - recontructed points and tracks
 
    AliLoader* loader = GetLoader();
    if (loader == 0x0)
     {
       Error("MakeTree","Can not get loader for %s",GetName());
       return;
     }
    loader->MakeTree(option); //delegate this job to getter
 }

//_______________________________________________________________________
AliLoader* AliDetector::MakeLoader(const char* topfoldername)
{ 
//builds standard getter (AliLoader type)
//if detector wants to use castomized getter, it must overload this method

 if (GetDebug())
   Info("MakeLoader",
        "Creating standard getter for detector %s. Top folder is %s.",
         GetName(),topfoldername);
     
 fLoader = new AliLoader(GetName(),topfoldername);
 return fLoader;
}

//_______________________________________________________________________
TTree* AliDetector::TreeH() const
{
//Get the hits container from the folder
  if (GetLoader() == 0x0) 
    {
    //sunstitude this with make getter when we can obtain the event folder name 
     Error("TreeH","Can not get the getter");
     return 0x0;
    }
 
  TTree* tree = (TTree*)GetLoader()->TreeH();
  return tree;
}

