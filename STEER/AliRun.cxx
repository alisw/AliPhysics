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

#include <TBRIK.h> 
#include <TCint.h> 
#include <TDatabasePDG.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
// 
#include "AliLog.h"
#include "AliDetector.h"
#include "AliDisplay.h"
#include "AliHeader.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliMC.h"
#include "AliMagFC.h"
#include "AliMagFCM.h"
#include "AliMagFDM.h"
#include "AliPDG.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliAlignObj.h"

AliRun *gAlice;

ClassImp(AliRun)

//_______________________________________________________________________
AliRun::AliRun():
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fModules(0),
  fGeometry(0),
  fMCApp(0),
  fDisplay(0),
  fField(0),
  fNdets(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(0),  //Particle factory object
  fConfigFunction(""),
  fRandom(0),
  fBaseFileName(""),
  fIsRootGeometry(kFALSE),
  fGeometryFileName(""),
  fTriggerDescriptor(""),
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
  TNamed(arun),
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fModules(0),
  fGeometry(0),
  fMCApp(0),
  fDisplay(0),
  fField(0),
  fNdets(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(0),  //Particle factory object
  fConfigFunction("\0"),
  fRandom(0),
  fBaseFileName(""),
  fIsRootGeometry(kFALSE),
  fGeometryFileName(""),
  fTriggerDescriptor(""),
  fRunLoader(0x0)
{
  //
  // Copy constructor for AliRun
  //
  arun.Copy(*this);
}

//_____________________________________________________________________________
AliRun::AliRun(const char *name, const char *title):
  TNamed(name,title),
  fRun(0),
  fEvent(0),
  fEventNrInRun(0),
  fEventsPerRun(0),
  fModules(new TObjArray(77)), // Support list for the Detectors
  fGeometry(0),
  fMCApp(0),
  fDisplay(0),
  fField(0),
  fNdets(0),
  fInitDone(kFALSE),
  fLego(0),
  fPDGDB(TDatabasePDG::Instance()),        //Particle factory object!
  fConfigFunction("Config();"),
  fRandom(new TRandom3()),
  fBaseFileName(""),
  fIsRootGeometry(kFALSE),
  fGeometryFileName(""),
  fTriggerDescriptor(""),
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
  
  // Create default mag field
  SetField();

  // Add particle list to configuration
  AliConfig::Instance()->Add(fPDGDB); 

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
  
  
  delete fField;
  delete fMCApp;
  delete gMC; gMC=0;
  delete fGeometry;
  delete fDisplay;
  delete fLego;
  if (fModules) {
    fModules->Delete();
    delete fModules;
  }
  
  delete fPDGDB;
}

//_______________________________________________________________________
void AliRun::Copy(TObject &) const
{
  AliFatal("Not implemented!");
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
void AliRun::ResetHits() 
{
  fMCApp->ResetHits();
}

//_______________________________________________________________________
AliGenerator* AliRun::Generator() const 
{
  return fMCApp->Generator();
}

//_______________________________________________________________________
void  AliRun::SetField(AliMagF* magField)
{
  //
  // Set Magnetic Field Map
  //
  fField = magField;
  fField->ReadField();
}

//_______________________________________________________________________
void AliRun::SetRootGeometry(Bool_t flag)
{
// Instruct application that the geometry is to be retreived from a root file.
   fIsRootGeometry = flag;
   if (flag) gMC->SetRootGeometry();
}
//_______________________________________________________________________
void AliRun::SetField(Int_t type, Int_t version, Float_t scale,
		      Float_t maxField, const char* filename)
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
    AliWarning(Form("Invalid map %d",version));
  }
}

//_____________________________________________________________________________

void AliRun::InitLoaders()
{
  //creates list of getters
  AliDebug(1, "");
  TIter next(fModules);
  AliModule *mod;
  while((mod = (AliModule*)next()))
   { 
     mod->SetRunLoader(fRunLoader);
     AliDetector *det = dynamic_cast<AliDetector*>(mod);
     if (det) 
      {
        AliDebug(2, Form("Adding %s", det->GetName()));
        fRunLoader->AddLoader(det);
      }
   }
  AliDebug(1, "Done");
}
//_____________________________________________________________________________

void AliRun::FinishRun()
{
  //
  // Called at the end of the run.
  //
  
  if(fLego) 
   {
    AliDebug(1, "Finish Lego");
    fRunLoader->CdGAFile();
    fLego->FinishRun();
   }
  
  // Clean detector information
  TIter next(fModules);
  AliModule *detector;
  while((detector = dynamic_cast<AliModule*>(next()))) {
    AliDebug(2, Form("%s->FinishRun()", detector->GetName()));
    detector->FinishRun();
  }
  
  AliDebug(1, "fRunLoader->WriteHeader(OVERWRITE)");
  fRunLoader->WriteHeader("OVERWRITE");

  // Write AliRun info and all detectors parameters
  fRunLoader->CdGAFile();
  Write(0,TObject::kOverwrite);//write AliRun
  fRunLoader->Write(0,TObject::kOverwrite);//write RunLoader itself
  
  // Clean tree information
  AliDebug(1, "fRunLoader->Stack()->FinishRun()");
  fRunLoader->Stack()->FinishRun();

  if(fMCApp) fMCApp->FinishRun();

  fRunLoader->Synchronize();
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
     AliError("RunLoader is not set. Can not load data.");
     return -1;
   }
/*****************************************/ 
/****   P R E    R E L O A D I N G    ****/
/*****************************************/ 
// Reset existing structures
  fMCApp->ResetHits();
  fMCApp->ResetTrackReferences();
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
  // Create the TNode geometry for the event display
  if (!fGeometry) { 
    BuildSimpleGeometry();
    //
    // Unlink and relink nodes in detectors
    // This is bad and there must be a better way...
    //
  
    TIter next(fModules);
    AliModule *detector;
    while((detector = dynamic_cast<AliModule*>(next()))) {
      detector->BuildGeometry();
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
  }
  return fGeometry;
}

//_______________________________________________________________________
void AliRun::SetBaseFile(const char *filename)
{
  fBaseFileName = filename;
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
  // Initialize ALICE Simulation run
  //
  Announce();

  if(fInitDone) {
    AliWarning("Cannot initialise AliRun twice!");
    return;
  }
    
  if (!fMCApp)  
    fMCApp=new AliMC(GetName(),GetTitle());
    
  gROOT->LoadMacro(setup);
  gInterpreter->ProcessLine(fConfigFunction.Data());

  fRunLoader->CdGAFile();

  AliPDG::AddParticlesToPdgDataBase();  

  fNdets = fModules->GetLast()+1;

  TIter next(fModules);
  for(Int_t i=0; i<fNdets; ++i)
   {
     TObject *objfirst, *objlast;
     AliModule *detector=dynamic_cast<AliModule*>(fModules->At(i));
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
   
   fMCApp->Init();
   
   //Must be here because some MCs (G4) adds detectors here and not in Config.C
   InitLoaders();
   fRunLoader->MakeTree("E");
   if (fLego == 0x0)
    {
      fRunLoader->LoadKinematics("RECREATE");
      fRunLoader->LoadTrackRefs("RECREATE");
      fRunLoader->LoadHits("all","RECREATE");
    }
   fInitDone = kTRUE;
   //
   // Save stuff at the beginning of the file to avoid file corruption
   fRunLoader->CdGAFile();
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
   AliDebug(1, Form("Found %d events", nev));
   Int_t nFirst = first;
   Int_t nLast  = (last < 0)? nev : last;
   
   for (Int_t nevent = nFirst; nevent <= nLast; nevent++) {
     AliDebug(1, Form("Processing event %d", nevent));
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
          AliDebug(1, Form("Processing Hits2SDigits for %s ...", detector->GetName()));
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
          AliDebug(1, Form("Processing SDigits2Digits for %s ...", detector->GetName()));
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
          AliDebug(1, Form("Processing Digits2Reco for %s ...", detector->GetName()));
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
  // If runloader has been initialized, set the number of events per file to nc1 * nc2

  // Set new generator
  if (!gener) gener  = new AliLegoGenerator();
  //
  // Configure Generator
  gener->SetRadiusRange(rmin, rmax);
  gener->SetZMax(zmax);
  gener->SetCoor1Range(nc1, c1min, c1max);
  gener->SetCoor2Range(nc2, c2min, c2max);
  
  
  //Create Lego object  
  fLego = new AliLego("lego",gener);

  if (!fInitDone) InitMC(setup);
  //Save current generator
  
  AliGenerator *gen=fMCApp->Generator();
  fMCApp->ResetGenerator(gener);
  //Prepare MC for Lego Run
  gMC->InitLego();
  
  //Run Lego Object

  if (fRunLoader) fRunLoader->SetNumberOfEventsPerFile(nc1 * nc2);
  //gMC->ProcessRun(nc1*nc2+1);
  gMC->ProcessRun(nc1*nc2);
  
  // End of this run, close files
  FinishRun();
  // Restore current generator
  fMCApp->ResetGenerator(gen);
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

// 
// MC Application
// 

//_______________________________________________________________________
void AliRun::Field(const Double_t* x, Double_t *b) const
{
  //
  // Return the value of the magnetic field
  //
  Float_t xfloat[3];
  for (Int_t i=0; i<3; i++) xfloat[i] = x[i]; 
  
  if (Field()) {
    Float_t bfloat[3];
    Field()->Field(xfloat,bfloat);
    for (Int_t j=0; j<3; j++) b[j] = bfloat[j]; 
  } 
  else {
    AliError("No mag field defined!");
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

void AliRun::SetGenEventHeader(AliGenEventHeader* header)
{
  fRunLoader->GetHeader()->SetGenEventHeader(header);
}
//_______________________________________________________________________

Int_t AliRun::GetEvNumber() const
{ 
//Returns number of current event  
  if (fRunLoader == 0x0)
   {
     AliError("RunLoader is not set. Can not load data.");
     return -1;
   }

  return fRunLoader->GetEventNumber();
}
//_______________________________________________________________________

void AliRun::SetRunLoader(AliRunLoader* rloader)
{
  //
  // Set the loader of the run
  //
  fRunLoader = rloader;
  if (fRunLoader == 0x0) return;
  
  TString evfoldname;
  TFolder* evfold = fRunLoader->GetEventFolder();
  if (evfold) evfoldname = evfold->GetName();
  else AliWarning("Did not get Event Folder from Run Loader");
  
  if ( fRunLoader->GetAliRun() )
   {//if alrun already exists in folder
    if (fRunLoader->GetAliRun() != this )
     {//and is different than this - crash
       AliFatal("AliRun is already in Folder and it is not this object");
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
     module->SetRunLoader(fRunLoader);
     AliDetector* detector = dynamic_cast<AliDetector*>(module);
     if (detector)
      {
        AliLoader* loader = fRunLoader->GetLoader(detector);
        if (loader == 0x0)
         {
           AliError(Form("Can not get loader for detector %s", detector->GetName()));
         }
        else
         {
           AliDebug(1, Form("Setting loader for detector %s", detector->GetName()));
           detector->SetLoader(loader);
         }
      }
   }
}

void AliRun::AddModule(AliModule* mod)
{
  //
  // Add a module to the module list
  //
  if (mod == 0x0) return;
  if (strlen(mod->GetName()) == 0) return;
  if (GetModuleID(mod->GetName()) >= 0) return;
  
  AliDebug(1, mod->GetName());
  if (fRunLoader == 0x0) AliConfig::Instance()->Add(mod);
  else AliConfig::Instance()->Add(mod,fRunLoader->GetEventFolder()->GetName());

  Modules()->Add(mod);
  
  fNdets++;
}

// added by Alberto Colla
//_____________________________________________________________________________
/*inline*/ Bool_t AliRun::IsFileAccessible(const char* fnam, EAccessMode mode)
{ return !gSystem->AccessPathName(fnam,mode);}

//______________________________________________________
/*inline*/ Bool_t AliRun::IsFileAccessible(Char_t* name,EAccessMode mode)
{
  TString str = name; gSystem->ExpandPathName(str);
  return !gSystem->AccessPathName(str.Data(),mode);
}
