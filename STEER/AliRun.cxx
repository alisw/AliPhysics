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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "Riostream.h"
#include "TBRIK.h" 
#include "TBrowser.h"
#include "TCint.h" 
#include "TFile.h"
#include "TFolder.h"
#include "TGeometry.h"
#include "TNode.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliDisplay.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliMCQA.h"
#include "AliMagFC.h"
#include "AliMagFCM.h"
#include "AliMagFDM.h"
#include "AliPDG.h"
#include "AliRun.h"
#include "AliStack.h"

AliRun *gAlice;

ClassImp(AliRun)

//_______________________________________________________________________
AliRun::AliRun():
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fDebug(0),
  fHeader(0),
  fTreeD(0),
  fTreeS(0),
  fTreeH(0),
  fTreeTR(0),
  fTreeE(0),
  fTreeR(0),
  fModules(0),
  fGeometry(0),
  fDisplay(0),
  fTimer(),
  fField(0),
  fMC(0),
  fImedia(0),
  fNdets(0),
  fTrRmax(1.e10),
  fTrZmax(1.e10),
  fGenerator(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(0),  //Particle factory object
  fHitLists(0),
  fEventEnergy(0),
  fSummEnergy(0),
  fSum2Energy(0),
  fConfigFunction("\0"),
  fRandom(0),
  fMCQA(0),
  fTransParName("\0"),
  fBaseFileName("\0"),
  fStack(0),
  fTreeDFileName(""),
  fTreeDFile(0),
  fTreeSFileName(""),
  fTreeSFile(0),
  fTreeRFileName(""),
  fTreeRFile(0)
{
  //
  // Default constructor for AliRun
  //
}

//_______________________________________________________________________
AliRun::AliRun(const AliRun& arun):
  TVirtualMCApplication(arun),
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fDebug(0),
  fHeader(0),
  fTreeD(0),
  fTreeS(0),
  fTreeH(0),
  fTreeTR(0),
  fTreeE(0),
  fTreeR(0),
  fModules(0),
  fGeometry(0),
  fDisplay(0),
  fTimer(),
  fField(0),
  fMC(0),
  fImedia(0),
  fNdets(0),
  fTrRmax(1.e10),
  fTrZmax(1.e10),
  fGenerator(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(0),  //Particle factory object
  fHitLists(0),
  fEventEnergy(0),
  fSummEnergy(0),
  fSum2Energy(0),
  fConfigFunction("\0"),
  fRandom(0),
  fMCQA(0),
  fTransParName("\0"),
  fBaseFileName("\0"),
  fStack(0),
  fTreeDFileName(""),
  fTreeDFile(0),
  fTreeSFileName(""),
  fTreeSFile(0),
  fTreeRFileName(""),
  fTreeRFile(0)
{
  //
  // Copy constructor for AliRun
  //
  arun.Copy(*this);
}

//_____________________________________________________________________________
AliRun::AliRun(const char *name, const char *title):
  TVirtualMCApplication(name,title),
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fDebug(0),
  fHeader(new AliHeader()),
  fTreeD(0),
  fTreeS(0),
  fTreeH(0),
  fTreeTR(0),
  fTreeE(0),
  fTreeR(0),
  fModules(new TObjArray(77)), // Support list for the Detectors
  fGeometry(0),
  fDisplay(0),
  fTimer(),
  fField(0),
  fMC(gMC),
  fImedia(new TArrayI(1000)),
  fNdets(0),
  fTrRmax(1.e10),
  fTrZmax(1.e10),
  fGenerator(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(TDatabasePDG::Instance()),        //Particle factory object!
  fHitLists(new TList()),                  // Create HitLists list
  fEventEnergy(0),
  fSummEnergy(0),
  fSum2Energy(0),
  fConfigFunction("Config();"),
  fRandom(new TRandom3()),
  fMCQA(0),
  fTransParName("\0"),
  fBaseFileName("\0"),
  fStack(new AliStack(10000)),        //Particle stack
  fTreeDFileName(""),
  fTreeDFile(0),
  fTreeSFileName(""),
  fTreeSFile(0),
  fTreeRFileName(""),
  fTreeRFile(0)
{
  //
  //  Constructor for the main processor.
  //  Creates the geometry
  //  Creates the list of Detectors.
  //  Creates the list of particles.
  //

  gAlice     = this;

  // Set random number generator
  gRandom = fRandom;

  if (gSystem->Getenv("CONFIG_SEED")) {
     gRandom->SetSeed(static_cast<UInt_t>(atoi(gSystem->Getenv("CONFIG_SEED"))));
  }

  // Add to list of browsable  
  gROOT->GetListOfBrowsables()->Add(this,name);

  // Create the TNode geometry for the event display
  BuildSimpleGeometry();
  
  // Create default mag field
  SetField();

  // Prepare the tracking medium lists
  for(Int_t i=0;i<1000;i++) (*fImedia)[i]=-99;

  // Add particle list to configuration
  AliConfig::Instance()->Add(fPDGDB); 

  // Set transport parameters
  SetTransPar();
}


//_______________________________________________________________________
AliRun::~AliRun()
{
  //
  // Default AliRun destructor
  //
  TFile *curfil =0;
  if(fTreeE)curfil=fTreeE->GetCurrentFile();
  delete fImedia;
  delete fField;
  // delete fMC;
  delete gMC; gMC=0;
  delete fGeometry;
  delete fDisplay;
  delete fGenerator;
  delete fLego;
  delete fTreeD;
  delete fTreeH;
  delete fTreeTR;
  delete fTreeE;
  delete fTreeR;
  delete fTreeS;
  if (fModules) {
    fModules->Delete();
    delete fModules;
  }
  delete fStack;
  delete fHitLists;
  delete fPDGDB;
  delete fMCQA;
  delete fHeader;
  // avoid to delete TFile objects not owned by this object
  // avoid multiple deletions
  if(curfil == fTreeDFile) fTreeDFile=0;
  if(curfil == fTreeSFile) fTreeSFile=0;
  if(curfil == fTreeRFile) fTreeRFile=0;
  if(fTreeSFile == fTreeDFile) fTreeSFile=0;
  if(fTreeRFile == fTreeDFile) fTreeRFile=0;
  if(fTreeRFile == fTreeSFile) fTreeRFile=0;
  if(fTreeDFile){
    if(fTreeDFile->IsOpen())fTreeDFile->Close();
    delete fTreeDFile;
  }
  if(fTreeSFile){
    if(fTreeSFile->IsOpen())fTreeSFile->Close();
    delete fTreeSFile;
  }
  if(fTreeRFile){
    if(fTreeRFile->IsOpen())fTreeRFile->Close();
    delete fTreeRFile;
  }
  if (gROOT->GetListOfBrowsables())
    gROOT->GetListOfBrowsables()->Remove(this);

  gAlice=0;
}

//_______________________________________________________________________
void AliRun::Copy(AliRun &) const
{
  //
  //  Copy method ... not implemented
  //
  Fatal("Copy","Not implemented!\n");
}

//_______________________________________________________________________
void AliRun::AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const
{
  //
  //  Add a hit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) dynamic_cast<AliModule*>(dets[id])->AddHit(track,vol,hits);
}

//_______________________________________________________________________
void AliRun::AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const
{
  //
  // Add digit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) dynamic_cast<AliModule*>(dets[id])->AddDigit(tracks,digits);
}

//_______________________________________________________________________
void AliRun::Browse(TBrowser *b)
{
  //
  // Called when the item "Run" is clicked on the left pane
  // of the Root browser.
  // It displays the Root Trees and all detectors.
  //
  if(!fStack) fStack=fHeader->Stack();
  TTree*  pTreeK = fStack->TreeK();
    
  if (pTreeK) b->Add(pTreeK,pTreeK->GetName());
  if (fTreeH) b->Add(fTreeH,fTreeH->GetName());
  if (fTreeTR) b->Add(fTreeTR,fTreeH->GetName());
  if (fTreeD) b->Add(fTreeD,fTreeD->GetName());
  if (fTreeE) b->Add(fTreeE,fTreeE->GetName());
  if (fTreeR) b->Add(fTreeR,fTreeR->GetName());
  if (fTreeS) b->Add(fTreeS,fTreeS->GetName());
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    b->Add(detector,detector->GetName());
  }
  b->Add(fMCQA,"AliMCQA");
}

//_______________________________________________________________________
void AliRun::Build()
{
  //
  // Initialize Alice geometry
  // Dummy routine
  //
}
 
//_______________________________________________________________________
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

//_______________________________________________________________________
void AliRun::CleanDetectors()
{
  //
  // Clean Detectors at the end of event
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->FinishEvent();
  }
}

//_______________________________________________________________________
Int_t AliRun::DistancetoPrimitive(Int_t, Int_t) const
{
  //
  // Return the distance from the mouse to the AliRun object
  // Dummy routine
  //
  return 9999;
}

//_______________________________________________________________________
void AliRun::DumpPart (Int_t i) const
{
  //
  // Dumps particle i in the stack
  //
    fStack->DumpPart(i);
}

//_______________________________________________________________________
void AliRun::DumpPStack () const
{
  //
  // Dumps the particle stack
  //
    fStack->DumpPStack();
}

//_______________________________________________________________________
void  AliRun::SetField(AliMagF* magField)
{
    // Set Magnetic Field Map
    fField = magField;
    fField->ReadField();
}

//_______________________________________________________________________
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
 
//_______________________________________________________________________
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
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->FinishRun();
  }
  
  //Output energy summary tables
  EnergySummary();

  TFile *file = fTreeE->GetCurrentFile();

  file->cd();
  
  fTreeE->Write(0,TObject::kOverwrite);
  
  // Write AliRun info and all detectors parameters
  Write(0,TObject::kOverwrite);

  // Clean tree information

  fStack->FinishRun();
  
  if (fTreeH) {
    delete fTreeH; fTreeH = 0;
  }
  if (fTreeTR) {
    delete fTreeTR; fTreeTR = 0;
  }
  if (fTreeD) {
    delete fTreeD; fTreeD = 0;
  }
  if (fTreeR) {
    delete fTreeR; fTreeR = 0;
  }
//   if (fTreeE) {
//     delete fTreeE; fTreeE = 0;
//   }
  if (fTreeS) {
    delete fTreeS; fTreeS = 0;
  }
  fGenerator->FinishRun();
  
  // Close output file
  file->Write();
}

//_______________________________________________________________________
void AliRun::FlagTrack(Int_t track)
{
  // Delegate to stack
  //
    fStack->FlagTrack(track);
}
 
//_______________________________________________________________________
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
  Int_t ievent=fHeader->GetEvent()+1;
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
	fSum2Energy[ndep]=TMath::Min(static_cast<Float_t>(99.),TMath::Max(ed2,kzero));
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

//_______________________________________________________________________
void AliRun::Announce() const
{
  //
  // Announce the current version of AliRoot
  //
  printf("%70s",
	 "****************************************************************\n");
  printf("%6s","*");printf("%64s","*\n");

  printf("%6s","*");
  printf("    You are running AliRoot version v3-09-07\n");

  printf("%6s","*");
  printf("    The cvs tag for the current program is $Name$\n");

  printf("%6s","*");printf("%64s","*\n");
  printf("%70s",
	 "****************************************************************\n");
}

//_______________________________________________________________________
AliModule *AliRun::GetModule(const char *name) const
{
  //
  // Return pointer to detector from name
  //
  return dynamic_cast<AliModule*>(fModules->FindObject(name));
}
 
//_______________________________________________________________________
AliDetector *AliRun::GetDetector(const char *name) const
{
  //
  // Return pointer to detector from name
  //
  return dynamic_cast<AliDetector*>(fModules->FindObject(name));
}
 
//_______________________________________________________________________
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
 
//_______________________________________________________________________
Int_t AliRun::GetEvent(Int_t event)
{
  //
  // Connect the Trees Kinematics and Hits for event # event
  // Set branch addresses
  //

  // Reset existing structures
  ResetHits();
  ResetTrackReferences();
  ResetDigits();
  ResetSDigits();
  
  // Delete Trees already connected
  if (fTreeH) { delete fTreeH; fTreeH = 0;}
  if (fTreeTR) { delete fTreeTR; fTreeTR = 0;}
  if (fTreeD) { delete fTreeD; fTreeD = 0;}
  if (fTreeR) { delete fTreeR; fTreeR = 0;}
  if (fTreeS) { delete fTreeS; fTreeS = 0;}

 // Create the particle stack
  if (fHeader) delete fHeader; 
  fHeader = 0;
  
  // Get header from file
  if(fTreeE) {
      fTreeE->SetBranchAddress("Header", &fHeader);

      if (!fTreeE->GetEntry(event)) {
	Error("GetEvent","Cannot find event:%d\n",event);
	return -1;
      }
  }  
  else {
      Error("GetEvent","Cannot find Header Tree (TE)\n");
      return -1;
  }

  // Get the stack from the header, set fStack to 0 if it 
  // fails to get event
  //
  TFile *file = fTreeE->GetCurrentFile();
  char treeName[20];
  
  file->cd();

  if (fStack) delete fStack;
  fStack = fHeader->Stack();
  if (fStack) {
    if (!fStack->GetEvent(event)) fStack = 0;
  }

  // Get Hits Tree header from file
  sprintf(treeName,"TreeH%d",event);
  fTreeH = dynamic_cast<TTree*>(gDirectory->Get(treeName));
  if (!fTreeH) {
      Warning("GetEvent","cannot find Hits Tree for event:%d\n",event);
  }

  // Get TracReferences Tree header from file
  sprintf(treeName,"TreeTR%d",event);
  fTreeTR = dynamic_cast<TTree*>(gDirectory->Get(treeName));
  if (!fTreeTR) {
    Warning("GetEvent","cannot find TrackReferences Tree for event:%d\n",event);
  }

  // get current file name and compare with names containing trees S,D,R
  TString curfilname=static_cast<TString>(fTreeE->GetCurrentFile()->GetName());
  if(fTreeDFileName==curfilname)fTreeDFileName="";
  if(fTreeSFileName==curfilname)fTreeSFileName="";
  if(fTreeRFileName==curfilname)fTreeRFileName="";

  // Get Digits Tree header from file
  sprintf(treeName,"TreeD%d",event);

  if (!fTreeDFile && fTreeDFileName != "") {
    InitTreeFile("D",fTreeDFileName);
  }    
  if (fTreeDFile) {    
    fTreeD = dynamic_cast<TTree*>(fTreeDFile->Get(treeName));
  } else {
    fTreeD = dynamic_cast<TTree*>(file->Get(treeName));
  }
  if (!fTreeD) {
    // Warning("GetEvent","cannot find Digits Tree for event:%d\n",event);
  }
  if(fTreeDFileName != ""){
    if(fTreeDFileName==fTreeSFileName) {
      fTreeSFileName = "";
      fTreeSFile = fTreeDFile;
    }
    if(fTreeDFileName==fTreeRFileName) {
      fTreeRFileName = "";
      fTreeRFile = fTreeDFile;
    }
  }

  file->cd();

  // Get SDigits Tree header from file
  sprintf(treeName,"TreeS%d",event);
  if (!fTreeSFile && fTreeSFileName != "") {
    InitTreeFile("S",fTreeSFileName);
  } 
  if (fTreeSFile) {
    fTreeS = dynamic_cast<TTree*>(fTreeSFile->Get(treeName));
  } else {
    fTreeS = dynamic_cast<TTree*>(gDirectory->Get(treeName));
  }
  if (!fTreeS) {
    // Warning("GetEvent","cannot find SDigits Tree for event:%d\n",event);
  }

  if(fTreeSFileName != ""){
    if(fTreeSFileName==fTreeRFileName){
      fTreeRFileName = "";
      fTreeRFile = fTreeSFile;
    }
  }

  file->cd();
  
  // Get Reconstruct Tree header from file
  sprintf(treeName,"TreeR%d",event);
  if (!fTreeRFile && fTreeRFileName != "") {
    InitTreeFile("R",fTreeRFileName);
  } 
  if(fTreeRFile) {
    fTreeR = dynamic_cast<TTree*>(fTreeRFile->Get(treeName));
  } else {
    fTreeR = dynamic_cast<TTree*>(gDirectory->Get(treeName));
  }
  if (!fTreeR) {
    //    printf("WARNING: cannot find Reconstructed Tree for event:%d\n",event);
  }

  file->cd();
 
  // Set Trees branch addresses
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->SetTreeAddress();
  }

  fEvent=event; //MI change

  return fHeader->GetNtrack();
}

//_______________________________________________________________________
TGeometry *AliRun::GetGeometry()
{
  //
  // Import Alice geometry from current file
  // Return pointer to geometry object
  //
  if (!fGeometry) fGeometry = dynamic_cast<TGeometry*>(gDirectory->Get("AliceGeom"));
  //
  // Unlink and relink nodes in detectors
  // This is bad and there must be a better way...
  //
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    TList *dnodes=detector->Nodes();
    Int_t j;
    TNode *node, *node1;
    for ( j=0; j<dnodes->GetSize(); j++) {
      node = dynamic_cast<TNode*>(dnodes->At(j));
      node1 = fGeometry->GetNode(node->GetName());
      dnodes->Remove(node);
      dnodes->AddAt(node1,j);
    }
  }
  return fGeometry;
}

//_______________________________________________________________________
Int_t AliRun::GetPrimary(Int_t track) const
{
  //
  // return number of primary that has generated track
  //
    return fStack->GetPrimary(track);
}
 
//_______________________________________________________________________
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
    if((det=dynamic_cast<AliModule*>(dets[kz]))) {
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
      det=dynamic_cast<AliModule*>(dets[ind]);
      if(det)
	printf(" %6s: %3d -> %3d;",det->GetName(),det->LoMedium(),
	       det->HiMedium());
      else
	printf(" %6s: %3d -> %3d;","NULL",0,0);
    }
    printf("\n");
  }
}

//_______________________________________________________________________
void AliRun::SetGenerator(AliGenerator *generator)
{
  //
  // Load the event generator
  //
  if(!fGenerator) fGenerator = generator;
}

//_______________________________________________________________________
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

//_______________________________________________________________________
void AliRun::SetTransPar(const char *filename)
{
  //
  // Sets the file name for transport parameters
  //
  fTransParName = filename;
}

//_______________________________________________________________________
void AliRun::SetBaseFile(const char *filename)
{
  fBaseFileName = filename;
}

//_______________________________________________________________________
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
  if(fDebug) {
    printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
    printf(" *%59s\n","*");
    printf(" *       Please check carefully what you are doing!%10s\n","*");
    printf(" *%59s\n","*");
  }
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
      if(fDebug){
	printf(" *%59s\n","*");
	printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
      }
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
	    if(fDebug) printf(" *  %-6s set to %10.3E for tracking medium code %4d for %s\n",
		   kpars[kz],cut[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,kpars[kz],cut[kz]);
	  }
	}
	// Set transport mechanisms
	for(kz=0;kz<knflags;kz++) {
	  if(flag[kz]>=0) {
	    if(fDebug) printf(" *  %-6s set to %10d for tracking medium code %4d for %s\n",
		   kpars[kncuts+kz],flag[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,kpars[kncuts+kz],Float_t(flag[kz]));
	  }
	}
      } else {
	Warning("ReadTransPar","Invalid medium code %d *\n",itmed);
	continue;
      }
    } else {
      if(fDebug) printf("%s::ReadTransParModule: %s not present\n",ClassName(),detName);
      continue;
    }
  }
}


//_______________________________________________________________________
void AliRun::MakeTree(Option_t *option, const char *file)
{
  //
  //  Create the ROOT trees
  //  Loop on all detectors to create the Root branch (if any)
  //

  char hname[30];
  //
  // Analyse options
  const char *oK = strstr(option,"K");
  const char *oH = strstr(option,"H");
  const char *oTR = strstr(option,"T");
  const char *oE = strstr(option,"E");
  const char *oD = strstr(option,"D");
  const char *oR = strstr(option,"R");
  const char *oS = strstr(option,"S");
  //

  TDirectory *cwd = gDirectory;

  TBranch *branch = 0;
  
  if (oK) fStack->MakeTree(fEvent, file);

  if (oE && !fTreeE) {
    fTreeE = new TTree("TE","Header");
    //    branch = fTreeE->Branch("Header", "AliHeader", &fHeader, 4000, 0);         
    branch = fTreeE->Branch("Header", "AliHeader", &fHeader, 4000, 0);         
    branch->SetAutoDelete(kFALSE);	     
    TFolder *folder = dynamic_cast<TFolder *>(gROOT->FindObjectAny("/Folders/RunMC/Event/Header"));
    if (folder) folder->Add(fHeader);
//    branch = fTreeE->Branch("Stack","AliStack", &fStack, 4000, 0);
//    branch->SetAutoDelete(kFALSE);	     
//    if (folder) folder->Add(fStack);
    fTreeE->Write(0,TObject::kOverwrite);
  }
  
  if (file && branch) {
    char * outFile = new char[strlen(gAlice->GetBaseFile())+strlen(file)+2];
    sprintf(outFile,"%s/%s",GetBaseFile(),file);
    branch->SetFile(outFile);
    TIter next( branch->GetListOfBranches());
    while ((branch=dynamic_cast<TBranch*>(next()))) {
       branch->SetFile(outFile);
    } 
    if (GetDebug()>1)
        printf("* MakeBranch * Diverting Branch %s to file %s\n", branch->GetName(),file);
    cwd->cd();
    delete outFile;
  }
  
  if (oH && !fTreeH) {
    sprintf(hname,"TreeH%d",fEvent);
    fTreeH = new TTree(hname,"Hits");
    fTreeH->SetAutoSave(1000000000); //no autosave
    fTreeH->Write(0,TObject::kOverwrite);
  }

  if (oTR && !fTreeTR) {
    sprintf(hname,"TreeTR%d",fEvent);
    fTreeTR = new TTree(hname,"TrackReferences");
    fTreeTR->SetAutoSave(1000000000); //no autosave
    fTreeTR->Write(0,TObject::kOverwrite);
  }

  if (oD && !fTreeD) {
    sprintf(hname,"TreeD%d",fEvent);
    fTreeD = new TTree(hname,"Digits");
    fTreeD->Write(0,TObject::kOverwrite);
  }
  if (oS && !fTreeS) {
    sprintf(hname,"TreeS%d",fEvent);
    fTreeS = new TTree(hname,"SDigits");
    fTreeS->Write(0,TObject::kOverwrite);
  }
  if (oR && !fTreeR) {
    sprintf(hname,"TreeR%d",fEvent);
    fTreeR = new TTree(hname,"Reconstruction");
    fTreeR->Write(0,TObject::kOverwrite);
  }

  //
  // Create a branch for hits/digits for each detector
  // Each branch is a TClonesArray. Each data member of the Hits classes
  // will be in turn a subbranch of the detector master branch
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     if (oH) detector->MakeBranch(option,file);
     if (oTR) detector->MakeBranchTR(option,file);
  }
}

//_______________________________________________________________________
TParticle* AliRun::Particle(Int_t i) const
{
  //
  // Returns particle i on the simulation stack
  //
  return fStack->Particle(i);
}

//_______________________________________________________________________
void AliRun::ResetDigits()
{
  //
  //  Reset all Detectors digits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     detector->ResetDigits();
  }
}

//_______________________________________________________________________
void AliRun::ResetSDigits()
{
  //
  //  Reset all Detectors digits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     detector->ResetSDigits();
  }
}

//_______________________________________________________________________
void AliRun::ResetHits()
{
  //
  //  Reset all Detectors hits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     detector->ResetHits();
  }
}

//_______________________________________________________________________
void AliRun::ResetTrackReferences()
{
  //
  //  Reset all Detectors hits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     detector->ResetTrackReferences();
  }
}

//_______________________________________________________________________
void AliRun::ResetPoints()
{
  //
  // Reset all Detectors points
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
     detector->ResetPoints();
  }
}

//_______________________________________________________________________
void AliRun::InitMC(const char *setup)
{
  //
  // Initialize the Alice setup
  //

  Announce();

  if(fInitDone) {
    Warning("Init","Cannot initialise AliRun twice!\n");
    return;
  }
    
  gROOT->LoadMacro(setup);
  gInterpreter->ProcessLine(fConfigFunction.Data());

  // Register MC in configuration 
  AliConfig::Instance()->Add(gMC);
  gMC->SetStack(fStack);

  gMC->DefineParticles();  //Create standard MC particles
  AliPDG::AddParticlesToPdgDataBase();  

  TObject *objfirst, *objlast;

  fNdets = fModules->GetLast()+1;

  //
  //=================Create Materials and geometry
  gMC->Init();

  // Added also after in case of interactive initialisation of modules
  fNdets = fModules->GetLast()+1;

   TIter next(fModules);
   AliModule *detector;
   while((detector = dynamic_cast<AliModule*>(next()))) {
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

   AliConfig::Instance();
   //
   // Save stuff at the beginning of the file to avoid file corruption
   Write();
}

//_______________________________________________________________________
void AliRun::RunMC(Int_t nevent, const char *setup)
{
  //
  // Main function to be called to process a galice run
  // example
  //    Root > gAlice.Run(); 
  // a positive number of events will cause the finish routine
  // to be called
  //
    fEventsPerRun = nevent;
  // check if initialisation has been done
  if (!fInitDone) InitMC(setup);
  
  // Create the Root Tree with one branch per detector

   MakeTree("ESDRT");

  if (gSystem->Getenv("CONFIG_SPLIT_FILE")) {
     MakeTree("K","Kine.root");
     MakeTree("H","Hits.root");
  } else {
     MakeTree("KH");
  }

  gMC->ProcessRun(nevent);

  // End of this run, close files
  if(nevent>0) FinishRun();
}

//_______________________________________________________________________
void AliRun::RunReco(const char *selected, Int_t first, Int_t last)
{
  //
  // Main function to be called to reconstruct Alice event
  // 
   cout << "Found "<< gAlice->TreeE()->GetEntries() << "events" << endl;
   Int_t nFirst = first;
   Int_t nLast  = (last < 0)? static_cast<Int_t>(gAlice->TreeE()->GetEntries()) : last;
   
   for (Int_t nevent = nFirst; nevent <= nLast; nevent++) {
     cout << "Processing event "<< nevent << endl;
     GetEvent(nevent);
     // MakeTree("R");
     Digits2Reco(selected);
   }
}

//_______________________________________________________________________

void AliRun::Hits2Digits(const char *selected)
{

   // Convert Hits to sumable digits
   // 
   for (Int_t nevent=0; nevent<gAlice->TreeE()->GetEntries(); nevent++) {
     GetEvent(nevent);
    // MakeTree("D");
     Hits2SDigits(selected);
     SDigits2Digits(selected);
   }  
}


//_______________________________________________________________________

void AliRun::Tree2Tree(Option_t *option, const char *selected)
{
  //
  // Function to transform the content of
  //  
  // - TreeH to TreeS (option "S")
  // - TreeS to TreeD (option "D")
  // - TreeD to TreeR (option "R")
  // 
  // If multiple options are specified ("SDR"), transformation will be done in sequence for
  // selected detector and for all detectors if none is selected (detector string 
  // can contain blank separated list of detector names). 


   const char *oS = strstr(option,"S");
   const char *oD = strstr(option,"D");
   const char *oR = strstr(option,"R");
   
   TObjArray *detectors = Detectors();

   TIter next(detectors);

   AliDetector *detector = 0;

   TDirectory *cwd = gDirectory;

   TObject *obj;

   char outFile[32];
   
   while((obj = next())) {
     if (!dynamic_cast<AliModule*>(obj)) 
       Fatal("Tree2Tree","Wrong type in fModules array\n");
     if (!(detector = dynamic_cast<AliDetector*>(obj))) continue;
     if (selected) 
       if (strcmp(detector->GetName(),selected)) continue;
     if (!detector->IsActive()) continue; 
     if (gSystem->Getenv("CONFIG_SPLIT_FILE")) {          
       if (oS) {
	 sprintf(outFile,"SDigits.%s.root",detector->GetName());
	 detector->MakeBranch("S",outFile);
       }    
       if (oD) {
	 sprintf(outFile,"Digits.%s.root",detector->GetName());
	 detector->MakeBranch("D",outFile);
       }    
       if (oR) {
	 sprintf(outFile,"Reco.%s.root",detector->GetName());
	 detector->MakeBranch("R",outFile);
       }    
     } else {
       detector->MakeBranch(option);
     }
     
     cwd->cd(); 
     
     if (oS) {
       cout << "Hits2SDigits: Processing " << detector->GetName() << "..." << endl;
       detector->Hits2SDigits(); 
     }  
     if (oD) {
       cout << "SDigits2Digits: Processing " << detector->GetName() << "..." << endl;
       detector->SDigits2Digits();
     } 
     if (oR) {
       cout << "Digits2Reco: Processing " << detector->GetName() << "..." << endl;
       detector->Digits2Reco(); 
     }
     
     cwd->cd();        
   }   
}

//_______________________________________________________________________
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

  //gMC->ProcessRun(nc1*nc2+1);
  gMC->ProcessRun(nc1*nc2);
  
  // Create only the Root event Tree
  MakeTree("E");
  
  // End of this run, close files
  FinishRun();
  // Restore current generator
  ResetGenerator(gen);
  // Delete Lego Object
  delete fLego; fLego=0;
}

//_______________________________________________________________________
void AliRun::SetConfigFunction(const char * config) 
{
  //
  // Set the signature of the function contained in Config.C to configure
  // the run
  //
  fConfigFunction=config;
}

//_______________________________________________________________________
void AliRun::SetCurrentTrack(Int_t track)
{ 
  //
  // Set current track number
  //
    fStack->SetCurrentTrack(track); 
}
 
//_______________________________________________________________________
void AliRun::SetTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
		      Float_t *vpos, Float_t *polar, Float_t tof,
		      TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{ 
// Delegate to stack
//

     fStack->SetTrack(done, parent, pdg, pmom, vpos, polar, tof,
 		     mech, ntr, weight, is);
}

//_______________________________________________________________________
void AliRun::SetTrack(Int_t done, Int_t parent, Int_t pdg,
  	              Double_t px, Double_t py, Double_t pz, Double_t e,
  		      Double_t vx, Double_t vy, Double_t vz, Double_t tof,
		      Double_t polx, Double_t poly, Double_t polz,
		      TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{ 
  // Delegate to stack
  //
    fStack->SetTrack(done, parent, pdg, px, py, pz, e, vx, vy, vz, tof,
		   polx, poly, polz, mech, ntr, weight, is);
    
}

//_______________________________________________________________________
void AliRun::SetHighWaterMark(const Int_t nt)
{
    //
    // Set high water mark for last track in event
    fStack->SetHighWaterMark(nt);
}

//_______________________________________________________________________
void AliRun::KeepTrack(const Int_t track)
{ 
  //
  // Delegate to stack
  //
    fStack->KeepTrack(track);
}
 
// 
// MC Application
// 

//_______________________________________________________________________
void  AliRun::ConstructGeometry() 
{
  //
  // Create modules, materials, geometry
  //

    TStopwatch stw;
    TIter next(fModules);
    AliModule *detector;
    printf("Geometry creation:\n");
    while((detector = dynamic_cast<AliModule*>(next()))) {
      stw.Start();
      // Initialise detector materials and geometry
      detector->CreateMaterials();
      detector->CreateGeometry();
      printf("%10s R:%.2fs C:%.2fs\n",
	     detector->GetName(),stw.RealTime(),stw.CpuTime());
    }
}

//_______________________________________________________________________
void  AliRun::InitGeometry()
{ 
  //
  // Initialize detectors and display geometry
  //

   printf("Initialisation:\n");
    TStopwatch stw;
    TIter next(fModules);
    AliModule *detector;
    while((detector = dynamic_cast<AliModule*>(next()))) {
      stw.Start();
      // Initialise detector and display geometry
      detector->Init();
      detector->BuildGeometry();
      printf("%10s R:%.2fs C:%.2fs\n",
	     detector->GetName(),stw.RealTime(),stw.CpuTime());
    }
 
}

//_______________________________________________________________________
void  AliRun::GeneratePrimaries()
{ 
  //
  // Generate primary particles and fill them in the stack.
  //

  Generator()->Generate();
}

//_______________________________________________________________________
void AliRun::BeginEvent()
{
    // Clean-up previous event
    // Energy scores  
    fEventEnergy.Reset();  
    // Clean detector information
    CleanDetectors();
    // Reset stack info
    fStack->Reset();

 
  //
  //  Reset all Detectors & kinematics & trees
  //
  char hname[30];
  //
  // Initialise event header
  fHeader->Reset(fRun,fEvent,fEventNrInRun);
  //
  fStack->BeginEvent(fEvent);

  //
  if(fLego) {
    fLego->BeginEvent();
    return;
  }

  //

  ResetHits();
  ResetTrackReferences();
  ResetDigits();
  ResetSDigits();


  if(fTreeH) {
    fTreeH->Reset();
    sprintf(hname,"TreeH%d",fEvent);
    fTreeH->SetName(hname);
  }

  if(fTreeTR) {
    fTreeTR->Reset();
    sprintf(hname,"TreeTR%d",fEvent);
    fTreeTR->SetName(hname);
  }

  if(fTreeD) {
    fTreeD->Reset();
    sprintf(hname,"TreeD%d",fEvent);
    fTreeD->SetName(hname);
    fTreeD->Write(0,TObject::kOverwrite);
  }
  if(fTreeS) {
    fTreeS->Reset();
    sprintf(hname,"TreeS%d",fEvent);
    fTreeS->SetName(hname);
    fTreeS->Write(0,TObject::kOverwrite);
  }
  if(fTreeR) {
    fTreeR->Reset();
    sprintf(hname,"TreeR%d",fEvent);
    fTreeR->SetName(hname);
    fTreeR->Write(0,TObject::kOverwrite);
  }
}

//_______________________________________________________________________
void AliRun::BeginPrimary()
{
  //
  // Called  at the beginning of each primary track
  //
  
  // Reset Hits info
  gAlice->ResetHits();
  gAlice->ResetTrackReferences();

}

//_______________________________________________________________________
void AliRun::PreTrack()
{
  //
  // Method called before each track
  //
  TObjArray &dets = *fModules;
  AliModule *module;
  
  for(Int_t i=0; i<=fNdets; i++)
    if((module = dynamic_cast<AliModule*>(dets[i])))
      module->PreTrack();
  
  fMCQA->PreTrack();
}

//_______________________________________________________________________
void AliRun::Stepping() 
{
  //
  // Called at every step during transport
  //

  Int_t id = DetFromMate(gMC->GetMedium());
  if (id < 0) return;

  //
  // --- If lego option, do it and leave 
  if (fLego)
    fLego->StepManager();
  else {
    Int_t copy;
    //Update energy deposition tables
    AddEnergyDeposit(gMC->CurrentVolID(copy),gMC->Edep());
  
    //Call the appropriate stepping routine;
    AliModule *det = dynamic_cast<AliModule*>(fModules->At(id));
    if(det && det->StepManagerIsEnabled()) {
      fMCQA->StepManager(id);
      det->StepManager();
    }
  }
}

//_______________________________________________________________________
void AliRun::PostTrack()
{
  //
  // Called after a track has been trasported
  //
  TObjArray &dets = *fModules;
  AliModule *module;
  
  for(Int_t i=0; i<=fNdets; i++)
    if((module = dynamic_cast<AliModule*>(dets[i])))
      module->PostTrack();
}

//_______________________________________________________________________
void AliRun::FinishPrimary()
{
  //
  // Called  at the end of each primary track
  //
  
  //  static Int_t count=0;
  //  const Int_t times=10;
  // This primary is finished, purify stack
  fStack->PurifyKine();

  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->FinishPrimary();
  }

  // Write out hits if any
  if (gAlice->TreeH()) {
    gAlice->TreeH()->Fill();
  }

  // Write out hits if any
  if (gAlice->TreeTR()) {
    gAlice->TreeTR()->Fill();
  }
  
  //
  //  if(++count%times==1) gObjectTable->Print();
}

//_______________________________________________________________________
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


  
  // Update Header information 

  fHeader->SetNprimary(fStack->GetNprimary());
  fHeader->SetNtrack(fStack->GetNtrack());  

  
  // Write out the kinematics
  fStack->FinishEvent();
   
  // Write out the event Header information
  if (fTreeE) {
      fHeader->SetStack(fStack);
      fTreeE->Fill();
  }
  
  
  // Write Tree headers
  TTree* pTreeK = fStack->TreeK();
  if (pTreeK) pTreeK->Write(0,TObject::kOverwrite);
  if (fTreeH) fTreeH->Write(0,TObject::kOverwrite);
  if (fTreeTR) fTreeTR->Write(0,TObject::kOverwrite);
  
  ++fEvent;
  ++fEventNrInRun;
}

//_______________________________________________________________________
void AliRun::Field(const Double_t* x, Double_t *b) const
{
  //
  // Returns the magnetic field at point x[3]
  // Units are kGauss
  //
  Float_t xfloat[3];
  for (Int_t i=0; i<3; i++) xfloat[i] = x[i]; 
  
  if (Field()) {
    Float_t bfloat[3];
    Field()->Field(xfloat,bfloat);
    for (Int_t j=0; j<3; j++) b[j] = bfloat[j]; 
  } 
  else {
    printf("No mag field defined!\n");
    b[0]=b[1]=b[2]=0.;
  }

}      

// 
// End of MC Application
// 

//_______________________________________________________________________
void AliRun::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliRun.

  if (R__b.IsReading()) {
    if (!gAlice) gAlice = this;

    AliRun::Class()->ReadBuffer(R__b, this);
    //
    gROOT->GetListOfBrowsables()->Add(this,"Run");

    fTreeE = dynamic_cast<TTree*>(gDirectory->Get("TE"));
    if (fTreeE) {
	  fTreeE->SetBranchAddress("Header", &fHeader);
    }      
    else    Error("Streamer","cannot find Header Tree\n");
    
    fTreeE->GetEntry(0);
    gRandom = fRandom;
  } else {
    AliRun::Class()->WriteBuffer(R__b, this);
  }
}


//_______________________________________________________________________
Int_t AliRun::CurrentTrack() const {
  //
  // Returns current track
  //
  return fStack->CurrentTrack();
}

//_______________________________________________________________________
Int_t AliRun::GetNtrack() const {
  //
  // Returns number of tracks in stack
  //
  return fStack->GetNtrack();
}

//_______________________________________________________________________
TObjArray* AliRun::Particles() const {
  //
  // Returns pointer to Particles array
  //
  return fStack->Particles();
}

//_______________________________________________________________________
TTree* AliRun::TreeK() const {
  //
  // Returns pointer to the TreeK array
  //
  return fStack->TreeK();
}


//_______________________________________________________________________
void AliRun::SetGenEventHeader(AliGenEventHeader* header)
{
    fHeader->SetGenEventHeader(header);
}

//_______________________________________________________________________
TFile* AliRun::InitFile(TString fileName)
{
// 
// create the file where the whole tree will be saved
//
  TDirectory *wd = gDirectory;
  TFile* file = TFile::Open(fileName,"update");
  gDirectory = wd;
  if (!file->IsOpen()) {
    Error("Cannot open file, %s\n",fileName);
    return 0;
  }
  return file;
}

//_______________________________________________________________________
TFile* AliRun::InitTreeFile(Option_t *option, TString fileName)
{
  //
  // create the file where one of the following trees will be saved
  // trees:   S,D,R
  // WARNING: by default these trees are saved on the file on which
  // hits are stored. If you divert one of these trees, you cannot restore
  // it to the original file (usually galice.root) in the same aliroot session
  Bool_t oS = (strstr(option,"S")!=0);
  Bool_t oR = (strstr(option,"R")!=0);
  Bool_t oD = (strstr(option,"D")!=0);
  Int_t choice[3]; 
  for (Int_t i=0; i<3; i++) choice[i] = 0;
  if(oS)choice[0] = 1; 
  if(oD)choice[1] = 1;
  if(oR)choice[2] = 1;

  TFile *ptr=0;

  if(!(oS || oR || oD))return ptr;

  Int_t active[3];
  for (Int_t i=0; i<3; i++) active[i] = 0;
  if(fTreeSFileName != "") active[0] = 1;
  if(fTreeDFileName != "") active[1] = 1;
  if(fTreeDFileName != "") active[2] = 1;

  Bool_t alreadyopen1 = kFALSE;
  Bool_t alreadyopen2 = kFALSE;

  if(oS){
    // if already active and same name with non-null ptr
    if(active[0]==1 && fileName == fTreeSFileName && fTreeSFile){
      Warning("InitTreeFile","File %s already opened",fTreeSFileName.Data());
      ptr = fTreeSFile;
    }
    else {
      // if already active with different name with non-null ptr
      if(active[0]==1 && fileName != fTreeSFileName && fTreeSFile){
        // close the active files and also the other possible files in option
        CloseTreeFile(option);
      }
      fTreeSFileName = fileName;
      alreadyopen1 = 
        (active[1] == 1 && fTreeDFileName == fTreeSFileName && fTreeDFile);
      alreadyopen2 =
        (active[2] == 1 && fTreeRFileName == fTreeSFileName && fTreeRFile);
      if(!(alreadyopen1 || alreadyopen2)){
        ptr = InitFile(fileName);
        fTreeSFile = ptr;
      }
      else {
        if(alreadyopen1){fTreeSFile = fTreeDFile; ptr = fTreeSFile;}
        if(alreadyopen2){fTreeSFile = fTreeRFile; ptr = fTreeSFile;}
      }
      if(choice[1] == 1) { fTreeDFileName = fileName; fTreeDFile = ptr;}
      if(choice[2] == 1) { fTreeRFileName = fileName; fTreeRFile = ptr;}
    }
    return ptr;
  }

  if(oD){
    // if already active and same name with non-null ptr
    if(active[1]==1 && fileName == fTreeDFileName && fTreeDFile){
      Warning("InitTreeFile","File %s already opened",fTreeDFileName.Data());
      ptr = fTreeDFile;
    }
    else {
      // if already active with different name with non-null ptr
      if(active[1]==1 && fileName != fTreeDFileName && fTreeDFile){
        // close the active files and also the other possible files in option
        CloseTreeFile(option);
      }
      fTreeDFileName = fileName;
      alreadyopen1 = 
        (active[0] == 1 && fTreeSFileName == fTreeDFileName && fTreeSFile);
      alreadyopen2 =
        (active[2] == 1 && fTreeRFileName == fTreeDFileName && fTreeRFile);
      if(!(alreadyopen1 || alreadyopen2)){
        ptr = InitFile(fileName);
        fTreeDFile = ptr;
      }
      else {
        if(alreadyopen1){fTreeDFile = fTreeSFile; ptr = fTreeDFile;}
        if(alreadyopen2){fTreeDFile = fTreeRFile; ptr = fTreeDFile;}
      }
      if(choice[2] == 1) { fTreeRFileName = fileName; fTreeRFile = ptr;}
    }
    return ptr;
  }

  if(oR){
    // if already active and same name with non-null ptr
    if(active[2]==1 && fileName == fTreeRFileName && fTreeRFile){
      Warning("InitTreeFile","File %s already opened",fTreeRFileName.Data());
      ptr = fTreeRFile;
    }
    else {
      // if already active with different name with non-null ptr
      if(active[2]==1 && fileName != fTreeRFileName && fTreeRFile){
        // close the active files and also the other possible files in option
        CloseTreeFile(option);
      }
      fTreeRFileName = fileName;
      alreadyopen1 = 
        (active[1] == 1 && fTreeDFileName == fTreeRFileName && fTreeDFile);
      alreadyopen2 =
        (active[0]== 1 && fTreeSFileName == fTreeRFileName && fTreeSFile);
      if(!(alreadyopen1 || alreadyopen2)){
        ptr = InitFile(fileName);
        fTreeRFile = ptr;
      }
      else {
        if(alreadyopen1){fTreeRFile = fTreeDFile; ptr = fTreeRFile;}
        if(alreadyopen2){fTreeRFile = fTreeSFile; ptr = fTreeRFile;}
      }
    }
    return ptr;
  }
  return 0;
}

//_______________________________________________________________________
void AliRun::PrintTreeFile()
{
  //
  // prints the file names and pointer associated to S,D,R trees
  //
  cout<<"===================================================\n";
  TFile *file = fTreeE->GetCurrentFile();
  TString curfilname="";
  if(file)curfilname=static_cast<TString>(file->GetName());
  cout<<" Current tree file name: "<<curfilname<<endl;
  cout<<"Pointer: "<<file<<endl;
  cout<<" Tree S File name: "<<fTreeSFileName<<endl;
  cout<<"Pointer: "<<fTreeSFile<<endl<<endl;
  cout<<" Tree D File name: "<<fTreeDFileName<<endl;
  cout<<"Pointer: "<<fTreeDFile<<endl<<endl;
  cout<<" Tree R File name: "<<fTreeRFileName<<endl;
  cout<<"Pointer: "<<fTreeRFile<<endl<<endl;
  cout<<"===================================================\n";
}
//_______________________________________________________________________
void AliRun::CloseTreeFile(Option_t *option)
{
  // 
  // closes the file containing the tree specified in option
  // (S,D,R)
  //
  Bool_t oS = (strstr(option,"S")!=0);
  Bool_t oR = (strstr(option,"R")!=0);
  Bool_t oD = (strstr(option,"D")!=0);
  Bool_t none = !(oS || oR || oD);
  if(none)return;
  if(oS){
    fTreeSFileName = "";
    if(fTreeSFile){
      if(!((fTreeSFile == fTreeDFile) || (fTreeSFile == fTreeRFile)) &&
           fTreeSFile->IsOpen()){
        fTreeSFile->Close();
        delete fTreeSFile;
      }
    }
    fTreeSFile = 0;
  }
  if(oD){
    fTreeDFileName = "";
    if(fTreeDFile){
      if(!((fTreeDFile == fTreeRFile) || (fTreeDFile == fTreeSFile)) &&
         fTreeDFile->IsOpen()){
        fTreeDFile->Close();
        delete fTreeDFile;
      }
    }
    fTreeDFile = 0;
  }
  if(oR){
    fTreeRFileName = "";
    if(fTreeRFile){
      if(!((fTreeRFile == fTreeSFile) || (fTreeRFile == fTreeDFile)) &&
         fTreeRFile->IsOpen()){
        fTreeRFile->Close();
        delete fTreeRFile;
      }
    }
    fTreeRFile = 0;
  }
}

//_______________________________________________________________________
void AliRun::MakeTree(Option_t *option, TFile *file)
{
  //
  //  Create some trees in the separate file
  //
  const char *oD = strstr(option,"D");
  const char *oR = strstr(option,"R");
  const char *oS = strstr(option,"S");

  TDirectory *cwd = gDirectory;
  char hname[30];

  if (oD) {
    delete fTreeD;
    sprintf(hname,"TreeD%d",fEvent);
    file->cd();
    fTreeD = static_cast<TTree*>(file->Get("hname"));
    if (!fTreeD) {
      fTreeD = new TTree(hname,"Digits");
      fTreeD->Write(0,TObject::kOverwrite);
    }
    cwd->cd();
  }
  if (oS) {
    delete fTreeS;
    sprintf(hname,"TreeS%d",fEvent);
    file->cd();
    fTreeS = static_cast<TTree*>(file->Get("hname"));
    if (!fTreeS) {
      fTreeS = new TTree(hname,"SDigits");
      fTreeS->Write(0,TObject::kOverwrite);
    }
    cwd->cd();
  }

  if (oR) {
    delete fTreeR;
    sprintf(hname,"TreeR%d",fEvent);
    file->cd();
    fTreeR = static_cast<TTree*>(file->Get("hname"));
    if (!fTreeR) {
      fTreeR = new TTree(hname,"RecPoint");
      fTreeR->Write(0,TObject::kOverwrite);
    }
    cwd->cd();
  }
}
