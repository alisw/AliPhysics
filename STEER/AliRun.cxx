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

/*
$Log$
Revision 1.46  2000/12/14 19:29:27  fca
galice.cuts was not read any more

Revision 1.45  2000/11/30 07:12:49  alibrary
Introducing new Rndm and QA classes

Revision 1.44  2000/10/26 13:58:59  morsch
Add possibility to choose the lego generator (of type AliGeneratorLego or derived) when running
RunLego(). Default is the base class AliGeneratorLego.

Revision 1.43  2000/10/09 09:43:17  fca
Special remapping of hits for TPC and TRD. End-of-primary action introduced

Revision 1.42  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.41  2000/07/13 16:19:09  fca
Mainly coding conventions + some small bug fixes

Revision 1.40  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

Revision 1.39  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

Revision 1.38  2000/06/20 13:05:45  fca
Writing down the TREE headers before job starts

Revision 1.37  2000/06/09 20:05:11  morsch
Introduce possibility to chose magnetic field version 3: AliMagFDM + field02.dat

Revision 1.36  2000/06/08 14:03:58  hristov
Only one initializer for a default argument

Revision 1.35  2000/06/07 10:13:14  hristov
Delete only existent objects.

Revision 1.34  2000/05/18 10:45:38  fca
Delete Particle Factory properly

Revision 1.33  2000/05/16 13:10:40  fca
New method IsNewTrack and fix for a problem in Father-Daughter relations

Revision 1.32  2000/04/27 10:38:21  fca
Correct termination of Lego Run and introduce Lego getter in AliRun

Revision 1.31  2000/04/26 10:17:32  fca
Changes in Lego for G4 compatibility

Revision 1.30  2000/04/18 19:11:40  fca
Introduce variable Config.C function signature

Revision 1.29  2000/04/07 11:12:34  fca
G4 compatibility changes

Revision 1.28  2000/04/05 06:51:06  fca
Workaround for an HP compiler problem

Revision 1.27  2000/03/22 18:08:07  fca
Rationalisation of the virtual MC interfaces

Revision 1.26  2000/03/22 13:42:26  fca
SetGenerator does not replace an existing generator, ResetGenerator does

Revision 1.25  2000/02/23 16:25:22  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.24  2000/01/19 17:17:20  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.23  1999/12/03 11:14:31  fca
Fixing previous wrong checking

Revision 1.21  1999/11/25 10:40:08  fca
Fixing daughters information also in primary tracks

Revision 1.20  1999/10/04 18:08:49  fca
Adding protection against inconsistent Euclid files

Revision 1.19  1999/09/29 07:50:40  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Control class for Alice C++                                              //
//  Only one single instance of this class exists.                           //
//  The object is created in main program aliroot                            //
//  and is pointed by the global gAlice.                                     //
//                                                                           //
//   -Supports the list of all Alice Detectors (fModules).                 //
//   -Supports the list of particles (fParticles).                           //
//   -Supports the Trees.                                                    //
//   -Supports the geometry.                                                 //
//   -Supports the event display.                                            //
//Begin_Html
/*
<img src="picts/AliRunClass.gif">
*/
//End_Html
//Begin_Html
/*
<img src="picts/alirun.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#include <TFile.h>
#include <TRandom.h>
#include <TBRIK.h> 
#include <TNode.h> 
#include <TCint.h> 
#include <TSystem.h>
#include <TObjectTable.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TROOT.h>
#include "TBrowser.h"

#include "TParticle.h"
#include "AliRun.h"
#include "AliDisplay.h"
#include "AliMC.h"
#include "AliLego.h"
#include "AliMagFC.h"
#include "AliMagFCM.h"
#include "AliMagFDM.h"
#include "AliHit.h"
#include "TRandom3.h"
#include "AliMCQA.h"
#include "AliGenerator.h"
#include "AliLegoGenerator.h"

#include "AliDetector.h"

AliRun *gAlice;

static AliHeader *gAliHeader;

ClassImp(AliRun)

//_____________________________________________________________________________
AliRun::AliRun()
{
  //
  // Default constructor for AliRun
  //
  gAliHeader=&fHeader;
  fRun       = 0;
  fEvent     = 0;
  fCurrent   = -1;
  fModules = 0;
  fGenerator = 0;
  fTreeD     = 0;
  fTreeK     = 0;
  fTreeH     = 0;
  fTreeE     = 0;
  fTreeR     = 0;
  fParticles = 0;
  fGeometry  = 0;
  fDisplay   = 0;
  fField     = 0;
  fMC       = 0;
  fNdets     = 0;
  fImedia    = 0;
  fTrRmax    = 1.e10;
  fTrZmax    = 1.e10;
  fInitDone  = kFALSE;
  fLego      = 0;
  fPDGDB     = 0;        //Particle factory object!
  fHitLists  = 0;
  fConfigFunction    = "\0";
  fRandom = 0;
  fMCQA = 0;
  fTransParName = "\0";
}

//_____________________________________________________________________________
AliRun::AliRun(const char *name, const char *title)
  : TNamed(name,title)
{
  //
  //  Constructor for the main processor.
  //  Creates the geometry
  //  Creates the list of Detectors.
  //  Creates the list of particles.
  //
  Int_t i;
  
  gAlice     = this;
  fTreeD     = 0;
  fTreeK     = 0;
  fTreeH     = 0;
  fTreeE     = 0;
  fTreeR     = 0;
  fTrRmax    = 1.e10;
  fTrZmax    = 1.e10;
  fGenerator = 0;
  fInitDone  = kFALSE;
  fLego      = 0;
  fField     = 0;
  fConfigFunction    = "Config();";

  // Set random number generator
  gRandom = fRandom = new TRandom3();
  
  gROOT->GetListOfBrowsables()->Add(this,name);
  //
  // create the support list for the various Detectors
  fModules = new TObjArray(77);
  //
  // Create the TNode geometry for the event display
  
  BuildSimpleGeometry();
  

  fNtrack=0;
  fHgwmk=0;
  fCurrent=-1;
  gAliHeader=&fHeader;
  fRun       = 0;
  fEvent     = 0;
  //
  // Create the particle stack
  fParticles = new TClonesArray("TParticle",100);
  
  fDisplay = 0;
  //
  // Create default mag field
  SetField();
  //
  fMC      = gMC;
  //
  // Prepare the tracking medium lists
  fImedia = new TArrayI(1000);
  for(i=0;i<1000;i++) (*fImedia)[i]=-99;
  //
  // Make particles
  fPDGDB     = TDatabasePDG::Instance();        //Particle factory object!
  //
  // Create HitLists list
  fHitLists  = new TList();
  //
  SetTransPar();
}


//_____________________________________________________________________________
AliRun::~AliRun()
{
  //
  // Defaullt AliRun destructor
  //
  delete fImedia;
  delete fField;
  delete fMC;
  delete fGeometry;
  delete fDisplay;
  delete fGenerator;
  delete fLego;
  delete fTreeD;
  delete fTreeK;
  delete fTreeH;
  delete fTreeE;
  delete fTreeR;
  if (fModules) {
    fModules->Delete();
    delete fModules;
  }
  if (fParticles) {
    fParticles->Delete();
    delete fParticles;
  }
  delete fHitLists;
  delete fPDGDB;
}

//_____________________________________________________________________________
void AliRun::AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const
{
  //
  //  Add a hit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) ((AliModule*) dets[id])->AddHit(track,vol,hits);
}

//_____________________________________________________________________________
void AliRun::AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const
{
  //
  // Add digit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) ((AliModule*) dets[id])->AddDigit(tracks,digits);
}

//_____________________________________________________________________________
void AliRun::Browse(TBrowser *b)
{
  //
  // Called when the item "Run" is clicked on the left pane
  // of the Root browser.
  // It displays the Root Trees and all detectors.
  //
  if (fTreeK) b->Add(fTreeK,fTreeK->GetName());
  if (fTreeH) b->Add(fTreeH,fTreeH->GetName());
  if (fTreeD) b->Add(fTreeD,fTreeD->GetName());
  if (fTreeE) b->Add(fTreeE,fTreeE->GetName());
  if (fTreeR) b->Add(fTreeR,fTreeR->GetName());
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    b->Add(detector,detector->GetName());
  }
  b->Add(fMCQA,"AliMCQA");
}

//_____________________________________________________________________________
void AliRun::Build()
{
  //
  // Initialize Alice geometry
  // Dummy routine
  //
}
 
//_____________________________________________________________________________
void AliRun::BuildSimpleGeometry()
{
  //
  // Create a simple TNode geometry used by Root display engine
  //
  // Initialise geometry
  //
  fGeometry = new TGeometry("AliceGeom","Galice Geometry for Hits");
  new TMaterial("void","Vacuum",0,0,0);  //Everything is void
  TBRIK *brik = new TBRIK("S_alice","alice volume","void",2000,2000,3000);
  brik->SetVisibility(0);
  new TNode("alice","alice","S_alice");
}

//_____________________________________________________________________________
void AliRun::CleanDetectors()
{
  //
  // Clean Detectors at the end of event
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->FinishEvent();
  }
}

//_____________________________________________________________________________
void AliRun::CleanParents()
{
  //
  // Clean Particles stack.
  // Set parent/daughter relations
  //
  TClonesArray &particles = *(gAlice->Particles());
  TParticle *part;
  int i;
  for(i=0; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    if(!part->TestBit(kDaughtersBit)) {
      part->SetFirstDaughter(-1);
      part->SetLastDaughter(-1);
    }
  }
}

//_____________________________________________________________________________
Int_t AliRun::DistancetoPrimitive(Int_t, Int_t)
{
  //
  // Return the distance from the mouse to the AliRun object
  // Dummy routine
  //
  return 9999;
}

//_____________________________________________________________________________
void AliRun::DumpPart (Int_t i) const
{
  //
  // Dumps particle i in the stack
  //
  TClonesArray &particles = *fParticles;
  ((TParticle*) particles[i])->Print();
}

//_____________________________________________________________________________
void AliRun::DumpPStack () const
{
  //
  // Dumps the particle stack
  //
  TClonesArray &particles = *fParticles;
  printf(
	 "\n\n=======================================================================\n");
  for (Int_t i=0;i<fNtrack;i++) 
    {
      printf("-> %d ",i); ((TParticle*) particles[i])->Print();
      printf("--------------------------------------------------------------\n");
    }
  printf(
	 "\n=======================================================================\n\n");
}

void  AliRun::SetField(AliMagF* magField)
{
    // Set Magnetic Field Map
    fField = magField;
    fField->ReadField();
}

//_____________________________________________________________________________
void AliRun::SetField(Int_t type, Int_t version, Float_t scale,
		      Float_t maxField, char* filename)
{
  //
  //  Set magnetic field parameters
  //  type      Magnetic field transport flag 0=no field, 2=helix, 3=Runge Kutta
  //  version   Magnetic field map version (only 1 active now)
  //  scale     Scale factor for the magnetic field
  //  maxField  Maximum value for the magnetic field

  //
  // --- Sanity check on mag field flags
  if(fField) delete fField;
  if(version==1) {
    fField = new AliMagFC("Map1"," ",type,scale,maxField);
  } else if(version<=2) {
    fField = new AliMagFCM("Map2-3",filename,type,scale,maxField);
    fField->ReadField();
  } else if(version==3) {
    fField = new AliMagFDM("Map4",filename,type,scale,maxField);
    fField->ReadField();
  } else {
    Warning("SetField","Invalid map %d\n",version);
  }
}

//_____________________________________________________________________________
void AliRun::FillTree()
{
  //
  // Fills all AliRun TTrees
  //
  if (fTreeK) fTreeK->Fill();
  if (fTreeH) fTreeH->Fill();
  if (fTreeD) fTreeD->Fill();
  if (fTreeR) fTreeR->Fill();
}
 
//_____________________________________________________________________________
void AliRun::PreTrack()
{
     TObjArray &dets = *fModules;
     AliModule *module;

     for(Int_t i=0; i<=fNdets; i++)
       if((module = (AliModule*)dets[i]))
	 module->PreTrack();

     fMCQA->PreTrack();
}

//_____________________________________________________________________________
void AliRun::PostTrack()
{
     TObjArray &dets = *fModules;
     AliModule *module;

     for(Int_t i=0; i<=fNdets; i++)
       if((module = (AliModule*)dets[i]))
	 module->PostTrack();
}

//_____________________________________________________________________________
void AliRun::FinishPrimary()
{
  //
  // Called  at the end of each primary track
  //
  
  //  static Int_t count=0;
  //  const Int_t times=10;
  // This primary is finished, purify stack
  PurifyKine();

  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->FinishPrimary();
  }

  // Write out hits if any
  if (gAlice->TreeH()) {
    gAlice->TreeH()->Fill();
  }
  
  // Reset Hits info
  gAlice->ResetHits();

  //
  //  if(++count%times==1) gObjectTable->Print();
}

//_____________________________________________________________________________
void AliRun::FinishEvent()
{
  //
  // Called at the end of the event.
  //
  
  //
  if(fLego) fLego->FinishEvent();

  //Update the energy deposit tables
  Int_t i;
  for(i=0;i<fEventEnergy.GetSize();i++) {
    fSummEnergy[i]+=fEventEnergy[i];
    fSum2Energy[i]+=fEventEnergy[i]*fEventEnergy[i];
  }
  fEventEnergy.Reset();
  
  // Clean detector information
  CleanDetectors();
  
  // Write out the kinematics
  if (fTreeK) {
    CleanParents();
    fTreeK->Fill();
  }
  
  // Write out the digits
  if (fTreeD) {
    fTreeD->Fill();
    ResetDigits();
  }
  
  // Write out reconstructed clusters  
  if (fTreeR) {
    fTreeR->Fill();
  }

  // Write out the event Header information
  if (fTreeE) fTreeE->Fill();
  
  // Reset stack info
  ResetStack();
  
  // Write Tree headers
  //  Int_t ievent = fHeader.GetEvent();
  //  char hname[30];
  //  sprintf(hname,"TreeK%d",ievent);
  if (fTreeK) fTreeK->Write(0,TObject::kOverwrite);
  //  sprintf(hname,"TreeH%d",ievent);
  if (fTreeH) fTreeH->Write(0,TObject::kOverwrite);
  //  sprintf(hname,"TreeD%d",ievent);
  if (fTreeD) fTreeD->Write(0,TObject::kOverwrite);
  //  sprintf(hname,"TreeR%d",ievent);
  if (fTreeR) fTreeR->Write(0,TObject::kOverwrite);

  ++fEvent;
}

//_____________________________________________________________________________
void AliRun::FinishRun()
{
  //
  // Called at the end of the run.
  //

  //
  if(fLego) fLego->FinishRun();

  // Clean detector information
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->FinishRun();
  }
  
  //Output energy summary tables
  EnergySummary();
  
  // file is retrieved from whatever tree
  TFile *file = 0;
  if (fTreeK) file = fTreeK->GetCurrentFile();
  if ((!file) && (fTreeH)) file = fTreeH->GetCurrentFile();
  if ((!file) && (fTreeD)) file = fTreeD->GetCurrentFile();
  if ((!file) && (fTreeE)) file = fTreeE->GetCurrentFile();
  if( NULL==file ) {
    Error("FinishRun","There isn't root file!");
    exit(1);
  }
  file->cd();
  fTreeE->Write(0,TObject::kOverwrite);
  
  // Clean tree information
  if (fTreeK) {
    delete fTreeK; fTreeK = 0;
  }
  if (fTreeH) {
    delete fTreeH; fTreeH = 0;
  }
  if (fTreeD) {
    delete fTreeD; fTreeD = 0;
  }
  if (fTreeR) {
    delete fTreeR; fTreeR = 0;
  }
  if (fTreeE) {
    delete fTreeE; fTreeE = 0;
  }
  
  // Write AliRun info and all detectors parameters
  Write();
  
  // Close output file
  file->Write();
}

//_____________________________________________________________________________
void AliRun::FlagTrack(Int_t track)
{
  //
  // Flags a track and all its family tree to be kept
  //
  int curr;
  TParticle *particle;

  curr=track;
  while(1) {
    particle=(TParticle*)fParticles->UncheckedAt(curr);
    
    // If the particle is flagged the three from here upward is saved already
    if(particle->TestBit(kKeepBit)) return;
    
    // Save this particle
    particle->SetBit(kKeepBit);
    
    // Move to father if any
    if((curr=particle->GetFirstMother())==-1) return;
  }
}
 
//_____________________________________________________________________________
void AliRun::EnergySummary()
{
  //
  // Print summary of deposited energy
  //

  Int_t ndep=0;
  Float_t edtot=0;
  Float_t ed, ed2;
  Int_t kn, i, left, j, id;
  const Float_t kzero=0;
  Int_t ievent=fHeader.GetEvent()+1;
  //
  // Energy loss information
  if(ievent) {
    printf("***************** Energy Loss Information per event (GEV) *****************\n");
    for(kn=1;kn<fEventEnergy.GetSize();kn++) {
      ed=fSummEnergy[kn];
      if(ed>0) {
	fEventEnergy[ndep]=kn;
	if(ievent>1) {
	  ed=ed/ievent;
	  ed2=fSum2Energy[kn];
	  ed2=ed2/ievent;
	  ed2=100*TMath::Sqrt(TMath::Max(ed2-ed*ed,kzero))/ed;
	} else 
	  ed2=99;
	fSummEnergy[ndep]=ed;
	fSum2Energy[ndep]=TMath::Min((Float_t) 99.,TMath::Max(ed2,kzero));
	edtot+=ed;
	ndep++;
      }
    }
    for(kn=0;kn<(ndep-1)/3+1;kn++) {
      left=ndep-kn*3;
      for(i=0;i<(3<left?3:left);i++) {
	j=kn*3+i;
        id=Int_t (fEventEnergy[j]+0.1);
	printf(" %s %10.3f +- %10.3f%%;",gMC->VolName(id),fSummEnergy[j],fSum2Energy[j]);
      }
      printf("\n");
    }
    //
    // Relative energy loss in different detectors
    printf("******************** Relative Energy Loss per event ********************\n");
    printf("Total energy loss per event %10.3f GeV\n",edtot);
    for(kn=0;kn<(ndep-1)/5+1;kn++) {
      left=ndep-kn*5;
      for(i=0;i<(5<left?5:left);i++) {
	j=kn*5+i;
        id=Int_t (fEventEnergy[j]+0.1);
	printf(" %s %10.3f%%;",gMC->VolName(id),100*fSummEnergy[j]/edtot);
      }
      printf("\n");
    }
    for(kn=0;kn<75;kn++) printf("*"); 
    printf("\n");
  }
  //
  // Reset the TArray's
  //  fEventEnergy.Set(0);
  //  fSummEnergy.Set(0);
  //  fSum2Energy.Set(0);
}

//_____________________________________________________________________________
AliModule *AliRun::GetModule(const char *name) const
{
  //
  // Return pointer to detector from name
  //
  return (AliModule*)fModules->FindObject(name);
}
 
//_____________________________________________________________________________
AliDetector *AliRun::GetDetector(const char *name) const
{
  //
  // Return pointer to detector from name
  //
  return (AliDetector*)fModules->FindObject(name);
}
 
//_____________________________________________________________________________
Int_t AliRun::GetModuleID(const char *name) const
{
  //
  // Return galice internal detector identifier from name
  //
  Int_t i=-1;
  TObject *mod=fModules->FindObject(name);
  if(mod) i=fModules->IndexOf(mod);
  return i;
}
 
//_____________________________________________________________________________
Int_t AliRun::GetEvent(Int_t event)
{
  //
  // Connect the Trees Kinematics and Hits for event # event
  // Set branch addresses
  //

  // Reset existing structures
  ResetStack();
  ResetHits();
  ResetDigits();
  
  // Delete Trees already connected
  if (fTreeK) delete fTreeK;
  if (fTreeH) delete fTreeH;
  if (fTreeD) delete fTreeD;
  if (fTreeR) delete fTreeR;

  // Get header from file
  if(fTreeE) fTreeE->GetEntry(event);
  else Error("GetEvent","Cannot file Header Tree\n");
  
  // Get Kine Tree from file
  char treeName[20];
  sprintf(treeName,"TreeK%d",event);
  fTreeK = (TTree*)gDirectory->Get(treeName);
  if (fTreeK) fTreeK->SetBranchAddress("Particles", &fParticles);
  else    Error("GetEvent","cannot find Kine Tree for event:%d\n",event);
  
  // Get Hits Tree header from file
  sprintf(treeName,"TreeH%d",event);
  fTreeH = (TTree*)gDirectory->Get(treeName);
  if (!fTreeH) {
    Error("GetEvent","cannot find Hits Tree for event:%d\n",event);
  }
  
  // Get Digits Tree header from file
  sprintf(treeName,"TreeD%d",event);
  fTreeD = (TTree*)gDirectory->Get(treeName);
  if (!fTreeD) {
    Warning("GetEvent","cannot find Digits Tree for event:%d\n",event);
  }
  
  
  // Get Reconstruct Tree header from file
  sprintf(treeName,"TreeR%d",event);
  fTreeR = (TTree*)gDirectory->Get(treeName);
  if (!fTreeR) {
    //    printf("WARNING: cannot find Reconstructed Tree for event:%d\n",event);
  }
 
  // Set Trees branch addresses
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->SetTreeAddress();
  }
  
  if (fTreeK) fTreeK->GetEvent(0);
  fNtrack = Int_t (fParticles->GetEntries());
  return fNtrack;
}

//_____________________________________________________________________________
TGeometry *AliRun::GetGeometry()
{
  //
  // Import Alice geometry from current file
  // Return pointer to geometry object
  //
  if (!fGeometry) fGeometry = (TGeometry*)gDirectory->Get("AliceGeom");
  //
  // Unlink and relink nodes in detectors
  // This is bad and there must be a better way...
  //
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->SetTreeAddress();
    TList *dnodes=detector->Nodes();
    Int_t j;
    TNode *node, *node1;
    for ( j=0; j<dnodes->GetSize(); j++) {
      node = (TNode*) dnodes->At(j);
      node1 = fGeometry->GetNode(node->GetName());
      dnodes->Remove(node);
      dnodes->AddAt(node1,j);
    }
  }
  return fGeometry;
}

//_____________________________________________________________________________
void AliRun::GetNextTrack(Int_t &mtrack, Int_t &ipart, Float_t *pmom,
			  Float_t &e, Float_t *vpos, Float_t *polar,
			  Float_t &tof)
{
  //
  // Return next track from stack of particles
  //
  TVector3 pol;
  fCurrent=-1;
  TParticle *track;
  for(Int_t i=fNtrack-1; i>=0; i--) {
    track=(TParticle*) fParticles->UncheckedAt(i);
    if(!track->TestBit(kDoneBit)) {
      //
      // The track has not yet been processed
      fCurrent=i;
      ipart=track->GetPdgCode();
      pmom[0]=track->Px();
      pmom[1]=track->Py(); 
      pmom[2]=track->Pz();
      e     =track->Energy();
      vpos[0]=track->Vx();
      vpos[1]=track->Vy();
      vpos[2]=track->Vz();
      track->GetPolarisation(pol);
      polar[0]=pol.X();
      polar[1]=pol.Y();
      polar[2]=pol.Z();
      tof=track->T();
      track->SetBit(kDoneBit);
      break;
    }
  }
  mtrack=fCurrent;
  //
  // stop and start timer when we start a primary track
  Int_t nprimaries = fHeader.GetNprimary();
  if (fCurrent >= nprimaries) return;
  if (fCurrent < nprimaries-1) {
    fTimer.Stop();
    track=(TParticle*) fParticles->UncheckedAt(fCurrent+1);
    //    track->SetProcessTime(fTimer.CpuTime());
  }
  fTimer.Start();
}

//_____________________________________________________________________________
Int_t AliRun::GetPrimary(Int_t track) const
{
  //
  // return number of primary that has generated track
  //
  int current, parent;
  TParticle *part;
  //
  parent=track;
  while (1) {
    current=parent;
    part = (TParticle *)fParticles->UncheckedAt(current);
    parent=part->GetFirstMother();
    if(parent<0) return current;
  }
}
 
//_____________________________________________________________________________
void AliRun::InitMC(const char *setup)
{
  //
  // Initialize the Alice setup
  //

  if(fInitDone) {
    Warning("Init","Cannot initialise AliRun twice!\n");
    return;
  }

  gROOT->LoadMacro(setup);
  gInterpreter->ProcessLine(fConfigFunction.Data());

  gMC->DefineParticles();  //Create standard MC particles

  TObject *objfirst, *objlast;

  fNdets = fModules->GetLast()+1;

  //
  //=================Create Materials and geometry
  gMC->Init();

  // Added also after in case of interactive initialisation of modules
  fNdets = fModules->GetLast()+1;

   TIter next(fModules);
   AliModule *detector;
   while((detector = (AliModule*)next())) {
      detector->SetTreeAddress();
      objlast = gDirectory->GetList()->Last();
      
      // Add Detector histograms in Detector list of histograms
      if (objlast) objfirst = gDirectory->GetList()->After(objlast);
      else         objfirst = gDirectory->GetList()->First();
      while (objfirst) {
	detector->Histograms()->Add(objfirst);
	objfirst = gDirectory->GetList()->After(objfirst);
      }
   }
   ReadTransPar(); //Read the cuts for all materials
   
   MediaTable(); //Build the special IMEDIA table
   
   //Initialise geometry deposition table
   fEventEnergy.Set(gMC->NofVolumes()+1);
   fSummEnergy.Set(gMC->NofVolumes()+1);
   fSum2Energy.Set(gMC->NofVolumes()+1);
   
   //Compute cross-sections
   gMC->BuildPhysics();
   
   //Write Geometry object to current file.
   fGeometry->Write();
   
   fInitDone = kTRUE;

   fMCQA = new AliMCQA(fNdets);

   //
   // Save stuff at the beginning of the file to avoid file corruption
   Write();
}

//_____________________________________________________________________________
void AliRun::MediaTable()
{
  //
  // Built media table to get from the media number to
  // the detector id
  //
  Int_t kz, nz, idt, lz, i, k, ind;
  //  Int_t ibeg;
  TObjArray &dets = *gAlice->Detectors();
  AliModule *det;
  //
  // For all detectors
  for (kz=0;kz<fNdets;kz++) {
    // If detector is defined
    if((det=(AliModule*) dets[kz])) {
        TArrayI &idtmed = *(det->GetIdtmed()); 
        for(nz=0;nz<100;nz++) {
	// Find max and min material number
	if((idt=idtmed[nz])) {
	  det->LoMedium() = det->LoMedium() < idt ? det->LoMedium() : idt;
	  det->HiMedium() = det->HiMedium() > idt ? det->HiMedium() : idt;
	}
      }
      if(det->LoMedium() > det->HiMedium()) {
	det->LoMedium() = 0;
	det->HiMedium() = 0;
      } else {
	if(det->HiMedium() > fImedia->GetSize()) {
	  Error("MediaTable","Increase fImedia from %d to %d",
		fImedia->GetSize(),det->HiMedium());
	  return;
	}
	// Tag all materials in rage as belonging to detector kz
	for(lz=det->LoMedium(); lz<= det->HiMedium(); lz++) {
	  (*fImedia)[lz]=kz;
	}
      }
    }
  }
  //
  // Print summary table
  printf(" Traking media ranges:\n");
  for(i=0;i<(fNdets-1)/6+1;i++) {
    for(k=0;k< (6<fNdets-i*6?6:fNdets-i*6);k++) {
      ind=i*6+k;
      det=(AliModule*)dets[ind];
      if(det)
	printf(" %6s: %3d -> %3d;",det->GetName(),det->LoMedium(),
	       det->HiMedium());
      else
	printf(" %6s: %3d -> %3d;","NULL",0,0);
    }
    printf("\n");
  }
}

//____________________________________________________________________________
void AliRun::SetGenerator(AliGenerator *generator)
{
  //
  // Load the event generator
  //
  if(!fGenerator) fGenerator = generator;
}

//____________________________________________________________________________
void AliRun::ResetGenerator(AliGenerator *generator)
{
  //
  // Load the event generator
  //
  if(fGenerator)
    if(generator)
      Warning("ResetGenerator","Replacing generator %s with %s\n",
	      fGenerator->GetName(),generator->GetName());
    else
      Warning("ResetGenerator","Replacing generator %s with NULL\n",
	      fGenerator->GetName());
  fGenerator = generator;
}

//____________________________________________________________________________
void AliRun::SetTransPar(char *filename)
{
  fTransParName = filename;
}

//____________________________________________________________________________
void AliRun::ReadTransPar()
{
  //
  // Read filename to set the transport parameters
  //


  const Int_t kncuts=10;
  const Int_t knflags=11;
  const Int_t knpars=kncuts+knflags;
  const char kpars[knpars][7] = {"CUTGAM" ,"CUTELE","CUTNEU","CUTHAD","CUTMUO",
			       "BCUTE","BCUTM","DCUTE","DCUTM","PPCUTM","ANNI",
			       "BREM","COMP","DCAY","DRAY","HADR","LOSS",
			       "MULS","PAIR","PHOT","RAYL"};
  char line[256];
  char detName[7];
  char* filtmp;
  Float_t cut[kncuts];
  Int_t flag[knflags];
  Int_t i, itmed, iret, ktmed, kz;
  FILE *lun;
  //
  // See whether the file is there
  filtmp=gSystem->ExpandPathName(fTransParName.Data());
  lun=fopen(filtmp,"r");
  delete [] filtmp;
  if(!lun) {
    Warning("ReadTransPar","File %s does not exist!\n",fTransParName.Data());
    return;
  }
  //
  printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
  printf(" *%59s\n","*");
  printf(" *       Please check carefully what you are doing!%10s\n","*");
  printf(" *%59s\n","*");
  //
  while(1) {
    // Initialise cuts and flags
    for(i=0;i<kncuts;i++) cut[i]=-99;
    for(i=0;i<knflags;i++) flag[i]=-99;
    itmed=0;
    for(i=0;i<256;i++) line[i]='\0';
    // Read up to the end of line excluded
    iret=fscanf(lun,"%[^\n]",line);
    if(iret<0) {
      //End of file
      fclose(lun);
      printf(" *%59s\n","*");
      printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
      return;
    }
    // Read the end of line
    fscanf(lun,"%*c");
    if(!iret) continue;
    if(line[0]=='*') continue;
    // Read the numbers
    iret=sscanf(line,"%s %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d",
		detName,&itmed,&cut[0],&cut[1],&cut[2],&cut[3],&cut[4],&cut[5],&cut[6],&cut[7],&cut[8],
		&cut[9],&flag[0],&flag[1],&flag[2],&flag[3],&flag[4],&flag[5],&flag[6],&flag[7],
		&flag[8],&flag[9],&flag[10]);
    if(!iret) continue;
    if(iret<0) {
      //reading error
      Warning("ReadTransPar","Error reading file %s\n",fTransParName.Data());
      continue;
    }
    // Check that the module exist
    AliModule *mod = GetModule(detName);
    if(mod) {
      // Get the array of media numbers
      TArrayI &idtmed = *mod->GetIdtmed();
      // Check that the tracking medium code is valid
      if(0<=itmed && itmed < 100) {
	ktmed=idtmed[itmed];
	if(!ktmed) {
	  Warning("ReadTransPar","Invalid tracking medium code %d for %s\n",itmed,mod->GetName());
	  continue;
	}
	// Set energy thresholds
	for(kz=0;kz<kncuts;kz++) {
	  if(cut[kz]>=0) {
	    printf(" *  %-6s set to %10.3E for tracking medium code %4d for %s\n",
		   kpars[kz],cut[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,kpars[kz],cut[kz]);
	  }
	}
	// Set transport mechanisms
	for(kz=0;kz<knflags;kz++) {
	  if(flag[kz]>=0) {
	    printf(" *  %-6s set to %10d for tracking medium code %4d for %s\n",
		   kpars[kncuts+kz],flag[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,kpars[kncuts+kz],Float_t(flag[kz]));
	  }
	}
      } else {
	Warning("ReadTransPar","Invalid medium code %d *\n",itmed);
	continue;
      }
    } else {
      Warning("ReadTransPar","Module %s not present\n",detName);
      continue;
    }
  }
}

//_____________________________________________________________________________
void AliRun::MakeTree(Option_t *option)
{
  //
  //  Create the ROOT trees
  //  Loop on all detectors to create the Root branch (if any)
  //

  char hname[30];
  //
  // Analyse options
  char *oK = strstr(option,"K");
  char *oH = strstr(option,"H");
  char *oE = strstr(option,"E");
  char *oD = strstr(option,"D");
  char *oR = strstr(option,"R");
  //
  if (oK && !fTreeK) {
    sprintf(hname,"TreeK%d",fEvent);
    fTreeK = new TTree(hname,"Kinematics");
    //  Create a branch for particles
    fTreeK->Branch("Particles",&fParticles,4000);
    fTreeK->Write();
  }
  if (oH && !fTreeH) {
    sprintf(hname,"TreeH%d",fEvent);
    fTreeH = new TTree(hname,"Hits");
    fTreeH->SetAutoSave(1000000000); //no autosave
    fTreeH->Write();
  }
  if (oD && !fTreeD) {
    sprintf(hname,"TreeD%d",fEvent);
    fTreeD = new TTree(hname,"Digits");
    fTreeD->Write();
  }
  if (oR && !fTreeR) {
    sprintf(hname,"TreeR%d",fEvent);
    fTreeR = new TTree(hname,"Reconstruction");
    fTreeR->Write();
  }
  if (oE && !fTreeE) {
    fTreeE = new TTree("TE","Header");
    //  Create a branch for Header
    fTreeE->Branch("Header","AliHeader",&gAliHeader,4000);
    fTreeE->Write();
  }
  //
  // Create a branch for hits/digits for each detector
  // Each branch is a TClonesArray. Each data member of the Hits classes
  // will be in turn a subbranch of the detector master branch
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     if (oH || oD || oR) detector->MakeBranch(option);
  }
}

//_____________________________________________________________________________
Int_t AliRun::PurifyKine(Int_t lastSavedTrack, Int_t nofTracks)
{
  //
  // PurifyKine with external parameters
  //
  fHgwmk = lastSavedTrack;
  fNtrack = nofTracks;
  PurifyKine();
  return fHgwmk;
}

//_____________________________________________________________________________
void AliRun::PurifyKine()
{
  //
  // Compress kinematic tree keeping only flagged particles
  // and renaming the particle id's in all the hits
  //
  TClonesArray &particles = *fParticles;
  int nkeep=fHgwmk+1, parent, i;
  TParticle *part, *partnew, *father;
  int *map = new int[particles.GetEntries()];

  // Save in Header total number of tracks before compression
  fHeader.SetNtrack(fHeader.GetNtrack()+fNtrack-fHgwmk);

  // First pass, invalid Daughter information
  for(i=0; i<fNtrack; i++) {
    // Preset map, to be removed later
    if(i<=fHgwmk) map[i]=i ; else map[i] = -99;
    ((TParticle *)particles.UncheckedAt(i))->ResetBit(kDaughtersBit);
  }
  // Second pass, build map between old and new numbering
  for(i=fHgwmk+1; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    if(part->TestBit(kKeepBit)) {
      
      // This particle has to be kept
      map[i]=nkeep;
      if(i!=nkeep) {
	
	// Old and new are different, have to copy
	partnew = (TParticle *)particles.UncheckedAt(nkeep);
	// Change due to a bug in the HP compiler
	//	*partnew = *part;
	memcpy(partnew,part,sizeof(TParticle));
      } else partnew = part;
      
      // as the parent is always *before*, it must be already
      // in place. This is what we are checking anyway!
      if((parent=partnew->GetFirstMother())>fHgwmk) {
	if(map[parent]==-99) printf("map[%d] = -99!\n",parent);
	partnew->SetFirstMother(map[parent]);
      }
      nkeep++;
    }
  }
  fNtrack=nkeep;
  
  // Fix daughters information
  for (i=0; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    parent = part->GetFirstMother();
    if(parent>=0) {
      father = (TParticle *)particles.UncheckedAt(parent);
      if(father->TestBit(kDaughtersBit)) {
      
	if(i<father->GetFirstDaughter()) father->SetFirstDaughter(i);
	if(i>father->GetLastDaughter())  father->SetLastDaughter(i);
      } else {
	// Iitialise daughters info for first pass
	father->SetFirstDaughter(i);
	father->SetLastDaughter(i);
	father->SetBit(kDaughtersBit);
      }
    }
  }
  
#ifdef old
  // Now loop on all detectors and reset the hits
  AliHit *OneHit;
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    if (!detector->Hits()) continue;
    TClonesArray &vHits=*(detector->Hits());
    if(vHits.GetEntries() != detector->GetNhits())
      printf("vHits.GetEntries()!=detector->GetNhits(): %d != %d\n",
	     vHits.GetEntries(),detector->GetNhits());
    for (i=0; i<detector->GetNhits(); i++) {
      OneHit = (AliHit *)vHits.UncheckedAt(i);
      OneHit->SetTrack(map[OneHit->GetTrack()]);
    }
  }
#else

  // Now loop on all registered hit lists
  TIter next(fHitLists);
  TCollection *hitList;
  while((hitList = (TCollection*)next())) {
    TIter nexthit(hitList);
    AliHit *hit;
    while((hit = (AliHit*)nexthit())) {
      hit->SetTrack(map[hit->GetTrack()]);
    }
  }
#endif

  // 
  // This for detectors which have a special mapping mechanism
  // for hits, such as TPC and TRD
  //

   TIter nextmod(fModules);
   AliModule *detector;
   while((detector = (AliModule*)nextmod())) {
     detector->RemapTrackHitIDs(map);
   }
  

  fHgwmk=nkeep-1;
  particles.SetLast(fHgwmk);
  delete [] map;
}

//_____________________________________________________________________________
void AliRun::BeginEvent()
{
  //
  //  Reset all Detectors & kinematics & trees
  //
  char hname[30];
  //

  //
  if(fLego) {
    fLego->BeginEvent();
    return;
  }

  //
  ResetStack();
  ResetHits();
  ResetDigits();

  // Initialise event header
  fHeader.Reset(fRun,fEvent);

  if(fTreeK) {
    fTreeK->Reset();
    sprintf(hname,"TreeK%d",fEvent);
    fTreeK->SetName(hname);
  }
  if(fTreeH) {
    fTreeH->Reset();
    sprintf(hname,"TreeH%d",fEvent);
    fTreeH->SetName(hname);
  }
  if(fTreeD) {
    fTreeD->Reset();
    sprintf(hname,"TreeD%d",fEvent);
    fTreeD->SetName(hname);
  }
  if(fTreeR) {
    fTreeR->Reset();
    sprintf(hname,"TreeR%d",fEvent);
    fTreeR->SetName(hname);
  }
}

//_____________________________________________________________________________
void AliRun::ResetDigits()
{
  //
  //  Reset all Detectors digits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetDigits();
  }
}

//_____________________________________________________________________________
void AliRun::ResetHits()
{
  //
  //  Reset all Detectors hits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetHits();
  }
}

//_____________________________________________________________________________
void AliRun::ResetPoints()
{
  //
  // Reset all Detectors points
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetPoints();
  }
}

//_____________________________________________________________________________
void AliRun::RunMC(Int_t nevent, const char *setup)
{
  //
  // Main function to be called to process a galice run
  // example
  //    Root > gAlice.Run(); 
  // a positive number of events will cause the finish routine
  // to be called
  //

  // check if initialisation has been done
  if (!fInitDone) InitMC(setup);
  
  // Create the Root Tree with one branch per detector
  MakeTree("KHDER");

  gMC->ProcessRun(nevent);

  // End of this run, close files
  if(nevent>0) FinishRun();
}

//_____________________________________________________________________________
void AliRun::RunLego(const char *setup, Int_t nc1, Float_t c1min,
		     Float_t c1max,Int_t nc2,Float_t c2min,Float_t c2max,
		     Float_t rmin,Float_t rmax,Float_t zmax, AliLegoGenerator* gener)
{
  //
  // Generates lego plots of:
  //    - radiation length map phi vs theta
  //    - radiation length map phi vs eta
  //    - interaction length map
  //    - g/cm2 length map
  //
  //  ntheta    bins in theta, eta
  //  themin    minimum angle in theta (degrees)
  //  themax    maximum angle in theta (degrees)
  //  nphi      bins in phi
  //  phimin    minimum angle in phi (degrees)
  //  phimax    maximum angle in phi (degrees)
  //  rmin      minimum radius
  //  rmax      maximum radius
  //  
  //
  //  The number of events generated = ntheta*nphi
  //  run input parameters in macro setup (default="Config.C")
  //
  //  Use macro "lego.C" to visualize the 3 lego plots in spherical coordinates
  //Begin_Html
  /*
    <img src="picts/AliRunLego1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego3.gif">
  */
  //End_Html
  //

  // check if initialisation has been done
  if (!fInitDone) InitMC(setup);
  //Save current generator
  AliGenerator *gen=Generator();

  // Set new generator
  if (!gener) gener  = new AliLegoGenerator();
  ResetGenerator(gener);
  //
  // Configure Generator
  gener->SetRadiusRange(rmin, rmax);
  gener->SetZMax(zmax);
  gener->SetCoor1Range(nc1, c1min, c1max);
  gener->SetCoor2Range(nc2, c2min, c2max);
  
  
  //Create Lego object  
  fLego = new AliLego("lego",gener);

  //Prepare MC for Lego Run
  gMC->InitLego();
  
  //Run Lego Object

  gMC->ProcessRun(nc1*nc2+1);
  
  // Create only the Root event Tree
  MakeTree("E");
  
  // End of this run, close files
  FinishRun();
  // Restore current generator
  ResetGenerator(gen);
  // Delete Lego Object
  delete fLego; fLego=0;
}

//_____________________________________________________________________________
void AliRun::SetConfigFunction(const char * config) 
{
  //
  // Set the signature of the function contained in Config.C to configure
  // the run
  //
  fConfigFunction=config;
}

//_____________________________________________________________________________
void AliRun::SetCurrentTrack(Int_t track)
{ 
  //
  // Set current track number
  //
  fCurrent = track; 
}
 
//_____________________________________________________________________________
void AliRun::SetTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
		      Float_t *vpos, Float_t *polar, Float_t tof,
		      AliMCProcess mech, Int_t &ntr, Float_t weight)
{ 
  //
  // Load a track on the stack
  //
  // done     0 if the track has to be transported
  //          1 if not
  // parent   identifier of the parent track. -1 for a primary
  // pdg    particle code
  // pmom     momentum GeV/c
  // vpos     position 
  // polar    polarisation 
  // tof      time of flight in seconds
  // mecha    production mechanism
  // ntr      on output the number of the track stored
  //
  TClonesArray &particles = *fParticles;
  TParticle *particle;
  Float_t mass;
  const Int_t kfirstdaughter=-1;
  const Int_t klastdaughter=-1;
  const Int_t kS=0;
  //  const Float_t tlife=0;
  
  //
  // Here we get the static mass
  // For MC is ok, but a more sophisticated method could be necessary
  // if the calculated mass is required
  // also, this method is potentially dangerous if the mass
  // used in the MC is not the same of the PDG database
  //
  mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  Float_t e=TMath::Sqrt(mass*mass+pmom[0]*pmom[0]+
			pmom[1]*pmom[1]+pmom[2]*pmom[2]);
  
  //printf("Loading particle %s mass %f ene %f No %d ip %d pos %f %f %f mom %f %f %f kS %d m %s\n",
  //pname,mass,e,fNtrack,pdg,vpos[0],vpos[1],vpos[2],pmom[0],pmom[1],pmom[2],kS,mecha);
  
  particle=new(particles[fNtrack]) TParticle(pdg,kS,parent,-1,kfirstdaughter,
					     klastdaughter,pmom[0],pmom[1],pmom[2],
					     e,vpos[0],vpos[1],vpos[2],tof);
  ((TParticle*)particles[fNtrack])->SetPolarisation(TVector3(polar[0],polar[1],polar[2]));
  ((TParticle*)particles[fNtrack])->SetWeight(weight);
  particle->SetUniqueID(mech);
  if(!done) particle->SetBit(kDoneBit);
  //Declare that the daughter information is valid
  ((TParticle*)particles[fNtrack])->SetBit(kDaughtersBit);
  
  if(parent>=0) {
    particle=(TParticle*) fParticles->UncheckedAt(parent);
    particle->SetLastDaughter(fNtrack);
    if(particle->GetFirstDaughter()<0) particle->SetFirstDaughter(fNtrack);
  } else { 
    //
    // This is a primary track. Set high water mark for this event
    fHgwmk=fNtrack;
    //
    // Set also number if primary tracks
    fHeader.SetNprimary(fHgwmk+1);
    fHeader.SetNtrack(fHgwmk+1);
  }
  ntr = fNtrack++;
}

//_____________________________________________________________________________
void AliRun::KeepTrack(const Int_t track)
{ 
  //
  // flags a track to be kept
  //
  TClonesArray &particles = *fParticles;
  ((TParticle*)particles[track])->SetBit(kKeepBit);
}
 
//_____________________________________________________________________________
void AliRun::StepManager(Int_t id) 
{
  //
  // Called at every step during transport
  //

  //
  // --- If lego option, do it and leave 
  if (fLego)
    fLego->StepManager();
  else {
    Int_t copy;
    //Update energy deposition tables
    AddEnergyDeposit(gMC->CurrentVolID(copy),gMC->Edep());
  
    //Call the appropriate stepping routine;
    AliModule *det = (AliModule*)fModules->At(id);
    if(det) {
      fMCQA->StepManager(id);
      det->StepManager();
    }
  }
}

//_____________________________________________________________________________
void AliRun::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliRun.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    if (!gAlice) gAlice = this;
    gROOT->GetListOfBrowsables()->Add(this,"Run");
    fTreeE = (TTree*)gDirectory->Get("TE");
    if (fTreeE) fTreeE->SetBranchAddress("Header", &gAliHeader);
    else    Error("Streamer","cannot find Header Tree\n");
    R__b >> fNtrack;
    R__b >> fHgwmk;
    R__b >> fDebug;
    fHeader.Streamer(R__b);
    R__b >> fModules;
    R__b >> fParticles;
    R__b >> fField; 
    //    R__b >> fMC;
    R__b >> fNdets;
    R__b >> fTrRmax;
    R__b >> fTrZmax;
    R__b >> fGenerator;
    if(R__v>1) {
      R__b >> fPDGDB;        //Particle factory object!
      fTreeE->GetEntry(0);
    } else {
      fHeader.SetEvent(0);
      fPDGDB     = TDatabasePDG::Instance();        //Particle factory object!
    }
    if(R__v>2) {
      fConfigFunction.Streamer(R__b);
    } else {
      fConfigFunction="Config();";
    }
    if(R__v>3) {
      R__b >> fRandom;
      gRandom = fRandom;
      R__b >> fMCQA;
    } else {
      fRandom = gRandom = new TRandom3();
      fMCQA = new AliMCQA();
    }
  } else {
    R__b.WriteVersion(AliRun::IsA());
    TNamed::Streamer(R__b);
    R__b << fNtrack;
    R__b << fHgwmk;
    R__b << fDebug;
    fHeader.Streamer(R__b);
    R__b << fModules;
    R__b << fParticles;
    R__b << fField;
    //    R__b << fMC;
    R__b << fNdets;
    R__b << fTrRmax;
    R__b << fTrZmax;
    R__b << fGenerator;
    R__b << fPDGDB;        //Particle factory object!
    fConfigFunction.Streamer(R__b);
    R__b << fRandom;
    R__b << fMCQA;
  }
} 
