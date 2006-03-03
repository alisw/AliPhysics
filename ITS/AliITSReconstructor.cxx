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
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"

ClassImp(AliITSReconstructor)

//___________________________________________________________________________
AliITSReconstructor::AliITSReconstructor() : AliReconstructor(){
  // Default constructor
  fItsPID=0;
}
 //___________________________________________________________________________
AliITSReconstructor::~AliITSReconstructor(){
// destructor
  delete fItsPID;
} 
//______________________________________________________________________
AliITSReconstructor::AliITSReconstructor(const AliITSReconstructor &ob) :AliReconstructor(ob) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSpidESD2","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSReconstructor& AliITSReconstructor::operator=(const AliITSReconstructor& /* ob */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//_____________________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters


  AliLoader* loader = runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  gAlice=runLoader->GetAliRun();
  TDirectory* olddir = gDirectory;
  runLoader->CdGAFile();
  AliITSgeom* geom = (AliITSgeom*)gDirectory->Get("AliITSgeom");  
  olddir->cd();
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetLoader((AliITSLoader*)loader);
  rec->SetITSgeom(geom);
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

 
  AliLoader* loader = runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  gAlice=runLoader->GetAliRun();
  TDirectory* olddir = gDirectory;
  runLoader->CdGAFile();
  AliITSgeom* geom = (AliITSgeom*)gDirectory->Get("AliITSgeom");  
  olddir->cd();

  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetLoader((AliITSLoader*)loader);
  rec->SetITSgeom(geom);
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
  // V0 finding
  Double_t cuts[]={33,  // max. allowed chi2
		   0.16,// min. allowed negative daughter's impact parameter 
		   0.05,// min. allowed positive daughter's impact parameter 
		   0.08,// max. allowed DCA between the daughter tracks
		   0.99,// max. allowed cosine of V0's pointing angle
		   0.9,  // min. radius of the fiducial volume
		   2.9   // max. radius of the fiducial volume
  };
  AliV0vertexer vtxer(cuts);
  Double_t vtx[3], cvtx[6];
  esd->GetVertex()->GetXYZ(vtx);
  esd->GetVertex()->GetSigmaXYZ(cvtx);
  vtxer.SetVertex(vtx);
  vtxer.Tracks2V0vertices(esd);

  // cascade finding
  Double_t cts[]={33.,    // max. allowed chi2
		  0.05,   // min. allowed V0 impact parameter 
		  0.008,  // window around the Lambda mass 
		  0.035,  // min. allowed bachelor's impact parameter 
		  0.10,   // max. allowed DCA between a V0 and a track
		  0.9985, //max. allowed cosine of the cascade pointing angle
		  0.9,    // min. radius of the fiducial volume
		  2.9     // max. radius of the fiducial volume
  };
  AliCascadeVertexer cvtxer=AliCascadeVertexer(cts);
  cvtxer.SetVertex(vtx);
  cvtxer.V0sTracks2CascadeVertices(esd);
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
  
  TDirectory * olddir = gDirectory;
  runLoader->CdGAFile();
  AliITSgeom* geom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
  olddir->cd();
  if(!geom){
    Error("GetITSgeom","no ITS geometry available");
    return NULL;
  }
  
  return geom;
}

