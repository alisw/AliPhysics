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

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include <TBRIK.h> 
#include <TBrowser.h>
#include <TCint.h> 
#include <TFile.h>
#include <TFolder.h>
#include <TGeometry.h>
#include <TKey.h>
#include <TNode.h>
#include <TObjectTable.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliConfig.h"
#include "AliDetector.h"
#include "AliDisplay.h"
#include "AliGenEventHeader.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliHit.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliLoader.h"
#include "AliMCQA.h"
#include "AliMagFC.h"
#include "AliMagFCM.h"
#include "AliMagFDM.h"
#include "AliPDG.h"
#include "AliRun.h"
#include "AliRunLoader.h"
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
  fRunLoader(0x0)
{
  //
  // Default constructor for AliRun
  //
  AliConfig::Instance();//skowron 29 Feb 2002
                        //ensures that the folder structure is build
}

//_______________________________________________________________________
AliRun::AliRun(const AliRun& arun):
  TVirtualMCApplication(arun),
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fDebug(0),
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
  fRunLoader(0x0)
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
  fRunLoader(0x0)
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
  gROOT->GetListOfBrowsables()->Remove(this);

  if (fRunLoader)
   {
    TFolder* evfold = fRunLoader->GetEventFolder();
    TFolder* modfold = dynamic_cast<TFolder*>(evfold->FindObjectAny(AliConfig::GetModulesFolderName()));
    TIter next(fModules);
    AliModule *mod;
    while((mod = (AliModule*)next()))
     { 
       modfold->Remove(mod);
     }
   }
   
  delete fImedia;
  delete fField;
  // delete fMC;
  delete gMC; gMC=0;
  delete fGeometry;
  delete fDisplay;
  delete fGenerator;
  delete fLego;
  if (fModules) {
    fModules->Delete();
    delete fModules;
  }
  
  delete fHitLists;
  delete fPDGDB;
  delete fMCQA;
}

//_______________________________________________________________________
void AliRun::Copy(AliRun &) const
{
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
  //detectors are in folders anyway
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
   fRunLoader->CleanDetectors();
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
   if (fRunLoader->Stack())
    fRunLoader->Stack()->DumpPart(i);
}

//_______________________________________________________________________
void AliRun::DumpPStack () const
{
  //
  // Dumps the particle stack
  //
   if (fRunLoader->Stack())
    fRunLoader->Stack()->DumpPStack();
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

//_____________________________________________________________________________

void AliRun::InitLoaders()
{
  //creates list of getters
  if (GetDebug()) Info("InitLoaders","");
  TIter next(fModules);
  AliModule *mod;
  while((mod = (AliModule*)next()))
   { 
     AliDetector *det = dynamic_cast<AliDetector*>(mod);
     if (det) 
      {
        if (GetDebug()) Info("InitLoaders"," Adding %s ",det->GetName());
        fRunLoader->AddLoader(det);
      }
   }
  if (GetDebug()) Info("InitLoaders","Done");
}
//_____________________________________________________________________________

void AliRun::FinishRun()
{
  //
  // Called at the end of the run.
  //
  
  if(fLego) 
   {
    if (GetDebug()) Info("FinishRun"," Finish Lego");
    fRunLoader->CdGAFile();
    fLego->FinishRun();
   }
  
  // Clean detector information
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    if (GetDebug()) Info("FinishRun"," %s->FinishRun()",detector->GetName());
    detector->FinishRun();
  }
  
  //Output energy summary tables
  if (GetDebug()) Info("FinishRun"," EnergySummary()");
  EnergySummary();
  
  if (GetDebug()) Info("FinishRun"," fRunLoader->WriteHeader(OVERWRITE)");
  fRunLoader->WriteHeader("OVERWRITE");

  // Write AliRun info and all detectors parameters
  fRunLoader->CdGAFile();
  Write(0,TObject::kOverwrite);//write AliRun
  fRunLoader->Write(0,TObject::kOverwrite);//write RunLoader itself
  
  // Clean tree information
  if (GetDebug()) Info("FinishRun"," fRunLoader->Stack()->FinishRun()");
  fRunLoader->Stack()->FinishRun();

  // Clean detector information
  if (GetDebug()) Info("FinishRun"," fGenerator->FinishRun()");
  fGenerator->FinishRun();
  
  fRunLoader->Synchronize();
}

//_______________________________________________________________________
void AliRun::FlagTrack(Int_t track)
{
  // Delegate to stack
  //
    fRunLoader->Stack()->FlagTrack(track);
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
  Int_t ievent=fRunLoader->GetHeader()->GetEvent()+1;
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
  printf("    You are running AliRoot version NewIO\n");

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
// Reloads data containers in folders # event
// Set branch addresses
//
  if (fRunLoader == 0x0)
   {
     Error("GetEvent","RunLoader is not set. Can not load data.");
     return -1;
   }
/*****************************************/ 
/****   P R E    R E L O A D I N G    ****/
/*****************************************/ 
// Reset existing structures
  ResetHits();
  ResetTrackReferences();
  ResetDigits();
  ResetSDigits();

/*****************************************/ 
/****       R  E  L  O  A  D          ****/
/*****************************************/

  fRunLoader->GetEvent(event);

/*****************************************/ 
/****  P O S T    R E L O A D I N G   ****/
/*****************************************/ 

  // Set Trees branch addresses
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) 
   {
     detector->SetTreeAddress();
   }
 
  return fRunLoader->GetHeader()->GetNtrack();
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
    return fRunLoader->Stack()->GetPrimary(track);
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
//_____________________________________________________________________________

void AliRun::BeginEvent()
{
  //
  // Clean-up previous event
  // Energy scores
  if (GetDebug()) 
   {
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
     Info("BeginEvent","          BEGINNING EVENT               ");
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
     Info("BeginEvent",">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
   }
    
  /*******************************/    
  /*   Clean after eventual      */
  /*   previous event            */
  /*******************************/    

  
  //Set the next event in Run Loader -> Cleans trees (TreeK and all trees in detectors),
  fRunLoader->SetEventNumber(++fEventNrInRun);// sets new files, cleans the previous event stuff, if necessary, etc.,  
  if (GetDebug()) Info("BeginEvent","EventNr is %d",fEventNrInRun);
     
  fEventEnergy.Reset();  
    // Clean detector information
  
  if (fRunLoader->Stack())
    fRunLoader->Stack()->Reset();//clean stack -> tree is unloaded
  else
    fRunLoader->MakeStack();//or make a new one
    
  if (GetDebug()) Info("BeginEvent","  fRunLoader->MakeTree(K)");
  fRunLoader->MakeTree("K");
  if (GetDebug()) Info("BeginEvent","  gMC->SetStack(fRunLoader->Stack())");
  gMC->SetStack(fRunLoader->Stack());//Was in InitMC - but was moved here 
                                     //because we don't have guarantee that 
                                     //stack pointer is not going to change from event to event
	                 //since it bellobgs to header and is obtained via RunLoader
  //
  //  Reset all Detectors & kinematics & make/reset trees
  //
    
  fRunLoader->GetHeader()->Reset(fRun,fEvent,fEventNrInRun);
//  fRunLoader->WriteKinematics("OVERWRITE");  is there any reason to rewrite here since MakeTree does so

  if (GetDebug()) Info("BeginEvent","  fRunLoader->MakeTrackRefsContainer()");
  fRunLoader->MakeTrackRefsContainer();//for insurance

  if (GetDebug()) Info("BeginEvent","  ResetHits()");
  ResetHits();
  if (GetDebug()) Info("BeginEvent","  fRunLoader->MakeTree(H)");
  fRunLoader->MakeTree("H");

  //
  if(fLego) 
   {
    fLego->BeginEvent();
    return;
   }

  //create new branches and SetAdresses
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next()))
   {
    if (GetDebug()) Info("BeginEvent","  %s->MakeBranch(H)",detector->GetName());
    detector->MakeBranch("H"); 
    if (GetDebug()) Info("BeginEvent","  %s->MakeBranchTR()",detector->GetName());
    detector->MakeBranchTR();
    if (GetDebug()) Info("BeginEvent","  %s->SetTreeAddress()",detector->GetName());
    detector->SetTreeAddress();
   }
}

//_______________________________________________________________________
TParticle* AliRun::Particle(Int_t i) const
{
  if (fRunLoader)
   if (fRunLoader->Stack())
    return fRunLoader->Stack()->Particle(i);
  return 0x0;   
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

  InitLoaders();

  fRunLoader->MakeTree("E");
  fRunLoader->LoadKinematics("RECREATE");
  fRunLoader->LoadTrackRefs("RECREATE");
  fRunLoader->LoadHits("all","RECREATE");
  
  
  fRunLoader->CdGAFile();

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
  while((detector = dynamic_cast<AliModule*>(next()))) 
   {
     objlast = gDirectory->GetList()->Last();
      
     // Add Detector histograms in Detector list of histograms
     if (objlast) objfirst = gDirectory->GetList()->After(objlast);
     else         objfirst = gDirectory->GetList()->First();
     while (objfirst) 
      {
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
   fRunLoader->WriteGeometry();
   
   fInitDone = kTRUE;

   fMCQA = new AliMCQA(fNdets);

   //
   // Save stuff at the beginning of the file to avoid file corruption
   Write();
   fEventNrInRun = -1; //important - we start Begin event from increasing current number in run
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
  //Hits moved to begin event -> now we are crating separate tree for each event

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
   Int_t nev = fRunLoader->GetNumberOfEvents();
   if (GetDebug()) Info("RunReco","Found %d events",nev);
   Int_t nFirst = first;
   Int_t nLast  = (last < 0)? nev : last;
   
   for (Int_t nevent = nFirst; nevent <= nLast; nevent++) {
     if (GetDebug()) Info("RunReco","Processing event %d",nevent);
     GetEvent(nevent);
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

   while((detector = dynamic_cast<AliDetector*>(next()))) {
     if (selected) 
       if (strcmp(detector->GetName(),selected)) continue;
     if (detector->IsActive())
      { 
       
       AliLoader* loader = detector->GetLoader();
       if (loader == 0x0) continue;
       
       if (oS) 
        {
          if (GetDebug()) Info("Tree2Tree","Processing Hits2SDigits for %s ...",detector->GetName());
          loader->LoadHits("read");
          if (loader->TreeS() == 0x0) loader->MakeTree("S");
          detector->MakeBranch(option);
          detector->SetTreeAddress();
          detector->Hits2SDigits();
          loader->UnloadHits();
          loader->UnloadSDigits();
        }  
       if (oD) 
        {
          if (GetDebug()) Info("Tree2Tree","Processing SDigits2Digits for %s ...",detector->GetName());
          loader->LoadSDigits("read");
          if (loader->TreeD() == 0x0) loader->MakeTree("D");
          detector->MakeBranch(option);
          detector->SetTreeAddress();
          detector->SDigits2Digits();
          loader->UnloadSDigits();
          loader->UnloadDigits();
        } 
       if (oR) 
        {
          if (GetDebug()) Info("Tree2Tree","Processing Digits2Reco for %s ...",detector->GetName());
          loader->LoadDigits("read");
          if (loader->TreeR() == 0x0) loader->MakeTree("R");
          detector->MakeBranch(option);
          detector->SetTreeAddress();
          detector->Digits2Reco(); 
          loader->UnloadDigits();
          loader->UnloadRecPoints();

        }
     }   
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
  fRunLoader->MakeTree("E");
  
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
    fRunLoader->Stack()->SetCurrentTrack(track); 
}
 
//_______________________________________________________________________
void AliRun::PushTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
                      Float_t *vpos, Float_t *polar, Float_t tof,
                      TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{ 
// Delegate to stack
//
    fRunLoader->Stack()->PushTrack(done, parent, pdg, pmom, vpos, polar, tof,
		     mech, ntr, weight, is);
}

//_______________________________________________________________________
void AliRun::PushTrack(Int_t done, Int_t parent, Int_t pdg,
  	              Double_t px, Double_t py, Double_t pz, Double_t e,
  		      Double_t vx, Double_t vy, Double_t vz, Double_t tof,
		      Double_t polx, Double_t poly, Double_t polz,
		      TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{ 
  // Delegate to stack
  //
  fRunLoader->Stack()->PushTrack(done, parent, pdg, px, py, pz, e, vx, vy, vz, tof,
                                polx, poly, polz, mech, ntr, weight, is);
}

//_______________________________________________________________________
void AliRun::SetHighWaterMark(const Int_t nt)
{
    //
    // Set high water mark for last track in event
    fRunLoader->Stack()->SetHighWaterMark(nt);
}

//_______________________________________________________________________
void AliRun::KeepTrack(const Int_t track)
{ 
  //
  // Delegate to stack
  //
    fRunLoader->Stack()->KeepTrack(track);
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
    if (GetDebug()) Info("ConstructGeometry","Geometry creation:");
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
  fRunLoader->Stack()->PurifyKine();

  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->FinishPrimary();
    if(detector->GetLoader())
     {
       detector->GetLoader()->TreeH()->Fill();
     }
  }

  // Write out track references if any
  if (fRunLoader->TreeTR()) 
   {
    fRunLoader->TreeTR()->Fill();
   }
}

//_______________________________________________________________________
void AliRun::FinishEvent()
{
  //
  // Called at the end of the event.
  //
  
  //
  if(fLego) fLego->FinishEvent();

  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    detector->FinishEvent();
  }

  //Update the energy deposit tables
  Int_t i;
  for(i=0;i<fEventEnergy.GetSize();i++) 
   {
    fSummEnergy[i]+=fEventEnergy[i];
    fSum2Energy[i]+=fEventEnergy[i]*fEventEnergy[i];
   }

  AliHeader* header = fRunLoader->GetHeader();
  AliStack* stack = fRunLoader->Stack();
  if ( (header == 0x0) || (stack == 0x0) )
   {//check if we got header and stack. If not cry and exit aliroot
    Fatal("AliRun","Can not get the stack or header from LOADER");
    return;//never reached
   }  
  // Update Header information 
  header->SetNprimary(stack->GetNprimary());
  header->SetNtrack(stack->GetNtrack());  

  
  // Write out the kinematics
  stack->FinishEvent();
   
  // Write out the event Header information
  TTree* treeE = fRunLoader->TreeE();
  if (treeE) 
   {
      header->SetStack(stack);
      treeE->Fill();
   }
  else
   {
    Error("FinishEvent","Can not get TreeE from RL");
   }
  
  fRunLoader->WriteKinematics("OVERWRITE");
  fRunLoader->WriteTrackRefs("OVERWRITE");
  fRunLoader->WriteHits("OVERWRITE");

  if (GetDebug()) 
   { 
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
     Info("FinishEvent","          FINISHING EVENT               ");
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
     Info("FinishEvent","<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
   }
}

//_______________________________________________________________________
void AliRun::Field(const Double_t* x, Double_t *b) const
{
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
    gROOT->GetListOfBrowsables()->Add(this,"Run");

    gRandom = fRandom;
  } else {
    AliRun::Class()->WriteBuffer(R__b, this);
  }
}


//_______________________________________________________________________
Int_t AliRun::GetCurrentTrackNumber() const {
  //
  // Returns current track
  //
  return fRunLoader->Stack()->GetCurrentTrackNumber();
}

//_______________________________________________________________________
Int_t AliRun::GetNtrack() const {
  //
  // Returns number of tracks in stack
  //
  return fRunLoader->Stack()->GetNtrack();
}
//_______________________________________________________________________

//_______________________________________________________________________
TObjArray* AliRun::Particles() const {
  //
  // Returns pointer to Particles array
  //
  if (fRunLoader)
   if (fRunLoader->Stack())
    return fRunLoader->Stack()->Particles();
  return 0x0;
}

//___________________________________________________________________________

//_______________________________________________________________________
void AliRun::SetGenEventHeader(AliGenEventHeader* header)
{
  fRunLoader->GetHeader()->SetGenEventHeader(header);
}

//___________________________________________________________________________

Int_t AliRun::GetEvNumber() const
{ 
//Returns number of current event  
  if (fRunLoader == 0x0)
   {
     Error("GetEvent","RunLoader is not set. Can not load data.");
     return -1;
   }

  return fRunLoader->GetEventNumber();
}

void AliRun::SetRunLoader(AliRunLoader* rloader)
{
  fRunLoader = rloader;
  if (fRunLoader == 0x0) return;
  
  TString evfoldname;
  TFolder* evfold = fRunLoader->GetEventFolder();
  if (evfold) evfoldname = evfold->GetName();
  else Warning("SetRunLoader","Did not get Event Folder from Run Loader");
  
  if ( fRunLoader->GetAliRun() )
   {//if alrun already exists in folder
    if (fRunLoader->GetAliRun() != this )
     {//and is different than this - crash
       Fatal("AliRun","AliRun is already in Folder and it is not this object");
       return;//pro forma
     }//else do nothing
   }
  else
   {
     evfold->Add(this);//Post this AliRun to Folder
   }
  
  TIter next(fModules);
  AliModule *module;
  while((module = (AliModule*)next())) 
   {
     if (evfold) AliConfig::Instance()->Add(module,evfoldname);
     AliDetector* detector = dynamic_cast<AliDetector*>(module);
     if (detector)
      {
        AliLoader* loader = fRunLoader->GetLoader(detector);
        if (loader == 0x0)
         {
           Error("SetRunLoader","Can not get loader for detector %s",detector->GetName());
         }
        else
         {
           if (GetDebug()) Info("SetRunLoader","Setting loader for detector %s",detector->GetName());
           detector->SetLoader(loader);
         }
      }
   }
}

void AliRun::AddModule(AliModule* mod)
{
  if (mod == 0x0) return;
  if (strlen(mod->GetName()) == 0) return;
  if (GetModuleID(mod->GetName()) >= 0) return;
  
  if (GetDebug()) Info("AddModule","%s",mod->GetName());
  if (fRunLoader == 0x0) AliConfig::Instance()->Add(mod);
  else AliConfig::Instance()->Add(mod,fRunLoader->GetEventFolder()->GetName());

  Modules()->Add(mod);
}
