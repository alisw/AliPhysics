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

#include <TCint.h> 
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TTree.h>
// 
#include "AliLog.h"
#include "AliDetector.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliPDG.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliCDBManager.h"
#include "AliAlignObj.h"
#include "AliSimulation.h"
#include "AliLego.h"

AliRun *gAlice;

ClassImp(AliRun)

//_______________________________________________________________________
AliRun::AliRun():
//  fRun(-1),
  fEventNrInRun(-1),
  fModules(0),
  fMCApp(0),
  fNdets(0),
  fConfigFunction(""),
  fBaseFileName(""),
  fRunLoader(0x0)
{
  //
  // Default constructor for AliRun
  //
  AliConfig::Instance();//skowron 29 Feb 2002
                        //ensures that the folder structure is build

}

//_____________________________________________________________________________
AliRun::AliRun(const char *name, const char *title):
  TNamed(name,title),
  fEventNrInRun(-1),
  fModules(new TObjArray(77)), // Support list for the Detectors
  fMCApp(new AliMC(GetName(),GetTitle())),
  fNdets(0),
  fConfigFunction("Config();"),
  fBaseFileName(""),
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
  gRandom = new TRandom3();

  if (gSystem->Getenv("CONFIG_SEED")) {
     gRandom->SetSeed(static_cast<UInt_t>(atoi(gSystem->Getenv("CONFIG_SEED"))));
  }

  // Add to list of browsable  
  gROOT->GetListOfBrowsables()->Add(this,name);
  
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
    if(!modfold) AliFatal(Form("Folder %s not found\n",AliConfig::GetModulesFolderName().Data()));
    TIter next(fModules);
    AliModule *mod;
    while((mod = (AliModule*)next()))
     { 
       modfold->Remove(mod);
     }
   }
    
  delete fMCApp;
  delete gMC; gMC=0;
  if (fModules) {
    fModules->Delete();
    delete fModules;
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
  printf("    The SVN version for the current program is $Id$\n");

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
  fMCApp->ResetDigits();
  fMCApp->ResetSDigits();

/*****************************************/ 
/****       R  E  L  O  A  D          ****/
/*****************************************/

  AliRunLoader::Instance()->GetEvent(event);

/*****************************************/ 
/****  P O S T    R E L O A D I N G   ****/
/*****************************************/ 

  // Set Trees branch addresses
  TIter next(fModules);
  AliDetector *detector;
  while((detector = dynamic_cast<AliDetector*>(next()))) 
   {
     detector->SetTreeAddress();
   }
 
  return AliRunLoader::Instance()->GetHeader()->GetNtrack();
}

//_______________________________________________________________________
void AliRun::SetBaseFile(const char *filename)
{
  fBaseFileName = filename;
}


//_______________________________________________________________________
void AliRun::Hits2Digits(const char *selected)
{

   // Convert Hits to sumable digits
   // 
   for (Int_t nevent=0; nevent<AliRunLoader::Instance()->TreeE()->GetEntries(); nevent++) {
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
void AliRun::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliRun.

  if (R__b.IsReading()) {
    if (!gAlice) gAlice = this;
    AliRun::Class()->ReadBuffer(R__b, this);
    gROOT->GetListOfBrowsables()->Add(this,"Run");
    gRandom = new TRandom3();
  } else {
    AliRun::Class()->WriteBuffer(R__b, this);
  }
}

//_______________________________________________________________________
void AliRun::SetGenEventHeader(AliGenEventHeader* header)
{
  AliRunLoader::Instance()->GetHeader()->SetGenEventHeader(header);
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
  else AliFatal("Did not get Event Folder from Run Loader");
  
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

//_______________________________________________________________________
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

