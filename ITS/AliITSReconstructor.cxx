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
// class for ITS reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "AliITSReconstructor.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliITSDetTypeRec.h"
#include "AliITSLoader.h"
#include "AliITStrackerMI.h"
#include "AliITStrackerSA.h"
#include "AliITSVertexerIons.h"
#include "AliITSVertexerFast.h"
#include "AliITSVertexerPPZ.h"
#include "AliITSVertexerZ.h"
#include "AliESD.h"
#include "AliITSpidESD.h"
#include "AliITSpidESD1.h"
#include "AliITSpidESD2.h"
#include "AliITSInitGeometry.h"

ClassImp(AliITSReconstructor)

//___________________________________________________________________________
AliITSReconstructor::AliITSReconstructor() : AliReconstructor(),
fItsPID(0)
{
  // Default constructor
}
 //___________________________________________________________________________
AliITSReconstructor::~AliITSReconstructor(){
// destructor
  delete fItsPID;
} 
//______________________________________________________________________
AliITSReconstructor::AliITSReconstructor(const AliITSReconstructor &ob) :AliReconstructor(ob),
fItsPID(ob.fItsPID) 
{
  // Copy constructor
}

//______________________________________________________________________
AliITSReconstructor& AliITSReconstructor::operator=(const AliITSReconstructor&  ob ){
  // Assignment operator
  this->~AliITSReconstructor();
  new(this) AliITSReconstructor(ob);
  return *this;
}
//______________________________________________________________________
void AliITSReconstructor::Init(AliRunLoader *runLoader) const{
    // Initalize this constructor bet getting/creating the objects
    // nesseary for a proper ITS reconstruction.
    // Inputs:
    //    AliRunLoader *runLoader   Pointer to the run loader to allow
    //                              the getting of files/folders data
    //                              needed to do reconstruction
    // Output:
    //   none.
    // Return:
    //   none.

    AliITSInitGeometry *initgeom = new AliITSInitGeometry("AliITSvPPRasymmFMD",
							  2);
    AliITSgeom *geom = initgeom->CreateAliITSgeom();
    delete initgeom; // once created, do not need initgeom any more.
    AliITSLoader* loader = static_cast<AliITSLoader*>
	(runLoader->GetLoader("ITSLoader"));
    if (!loader) {
	Error("Init", "ITS loader not found");
	return;
    }
    loader->SetITSgeom(geom);
    return;
}
//_____________________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters


  AliITSLoader* loader = static_cast<AliITSLoader*>(runLoader->GetLoader("ITSLoader"));
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  AliITSDetTypeRec* rec = new AliITSDetTypeRec(loader);
  rec->SetDefaults();

  loader->LoadRecPoints("recreate");
  loader->LoadDigits("read");
  runLoader->LoadKinematics();
  TString option = GetOption();
  Bool_t clusfinder=kTRUE;   // Default: V2 cluster finder
  if(option.Contains("OrigCF"))clusfinder=kFALSE;

  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    runLoader->GetEvent(iEvent);
    if(loader->TreeR()==0x0) loader->MakeTree("R");
    rec->MakeBranch("R");
    rec->SetTreeAddress();
    rec->DigitsToRecPoints(iEvent,0,"All",clusfinder);    
  }

  loader->UnloadRecPoints();
  loader->UnloadDigits();
  runLoader->UnloadKinematics();
}

//_________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRunLoader* runLoader, 
                                      AliRawReader* rawReader) const
{
// reconstruct clusters

 
  AliITSLoader* loader = static_cast<AliITSLoader*>(runLoader->GetLoader("ITSLoader"));
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }

  AliITSDetTypeRec* rec = new AliITSDetTypeRec(loader);
  rec->SetDefaults();
  rec->SetDefaultClusterFindersV2(kTRUE);

  loader->LoadRecPoints("recreate");

  Int_t iEvent = 0;

  while(rawReader->NextEvent()) {
    runLoader->GetEvent(iEvent++);
    if(loader->TreeR()==0x0) loader->MakeTree("R");
    rec->DigitsToRecPoints(rawReader);
  }

  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
AliTracker* AliITSReconstructor::CreateTracker(AliRunLoader* runLoader)const
{
// create a ITS tracker

  
  AliITSgeom* geom = GetITSgeom(runLoader);
  TString selectedTracker = GetOption();
  AliTracker* tracker;    
  if (selectedTracker.Contains("MI")) {
    tracker = new AliITStrackerMI(geom);
  }
  else {
    tracker =  new AliITStrackerSA(geom);  // inherits from AliITStrackerMI
  }

  TString selectedPIDmethod = GetOption();
  AliITSLoader *loader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("CreateTracker", "ITS loader not found");
  }
  if(selectedPIDmethod.Contains("LandauFitPID")){
    loader->AdoptITSpid(new AliITSpidESD2((AliITStrackerMI*)tracker,loader));
  }
  else{
    Double_t parITS[] = {34., 0.15, 10.};
    loader->AdoptITSpid(new AliITSpidESD1(parITS));
  }
  return tracker;
  
}

//_____________________________________________________________________________
AliVertexer* AliITSReconstructor::CreateVertexer(AliRunLoader* /*runLoader*/) const
{
// create a ITS vertexer

  TString selectedVertexer = GetOption();
  if(selectedVertexer.Contains("ions") || selectedVertexer.Contains("IONS")){
    Info("CreateVertexer","a AliITSVertexerIons object has been selected\n");
    return new AliITSVertexerIons("null");
  }
  if(selectedVertexer.Contains("smear") || selectedVertexer.Contains("SMEAR")){
    Double_t smear[3]={0.005,0.005,0.01};
    Info("CreateVertexer","a AliITSVertexerFast object has been selected\n"); 
    return new AliITSVertexerFast(smear);
  }
  if(selectedVertexer.Contains("ppz") || selectedVertexer.Contains("PPZ")){
    Info("CreateVertexer","a AliITSVertexerPPZ object has been selected\n");
    return new AliITSVertexerPPZ("null");
  }
  // by default an AliITSVertexerZ object is instatiated
  Info("CreateVertexer","a AliITSVertexerZ object has been selected\n");
  return new AliITSVertexerZ("null");
}

//_____________________________________________________________________________
void AliITSReconstructor::FillESD(AliRunLoader* runLoader, 
				  AliESD* esd) const
{
// make PID, find V0s and cascade
  AliITSLoader *loader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  AliITSpidESD *pidESD = 0;
  TString selectedPIDmethod = GetOption();
  if(selectedPIDmethod.Contains("LandauFitPID")){
    Info("FillESD","ITS LandauFitPID option has been selected\n");
    pidESD=loader->GetITSpid();
  }
  else{
    Info("FillESD","ITS default PID\n");
    pidESD=loader->GetITSpid();
  }
  if(pidESD!=0){
    pidESD->MakePID(esd);
  }
  else {
    Error("FillESD","!! cannot do the PID !!\n");
  }
}


//_____________________________________________________________________________
AliITSgeom* AliITSReconstructor::GetITSgeom(AliRunLoader* runLoader) const
{
// get the ITS geometry

  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetITSgeom", "couldn't get AliRun object");
    return NULL;
  }
  AliITSLoader *loader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  AliITSgeom* geom = (AliITSgeom*)loader->GetITSgeom();
  if(!geom){
    Error("GetITSgeom","no ITS geometry available");
    return NULL;
  }
  
  return geom;
}
