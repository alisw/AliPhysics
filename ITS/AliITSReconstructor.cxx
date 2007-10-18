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
#include "AliRawReader.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITStrackerMI.h"
#include "AliITStrackerV2.h"
#include "AliITStrackerSA.h"
#include "AliITSVertexerIons.h"
#include "AliITSVertexerFast.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
#include "AliITSVertexerCosmics.h"
#include "AliESDEvent.h"
#include "AliITSpidESD.h"
#include "AliITSpidESD1.h"
#include "AliITSpidESD2.h"
#include "AliITSInitGeometry.h"


ClassImp(AliITSReconstructor)

AliITSRecoParam *AliITSReconstructor::fgkRecoParam =0;  // reconstruction parameters

//___________________________________________________________________________
AliITSReconstructor::AliITSReconstructor() : AliReconstructor(),
fItsPID(0),
fDetTypeRec(0)
{
  // Default constructor
  if (!fgkRecoParam) {
    AliError("The Reconstruction parameters nonitialized - Used default one");
    fgkRecoParam = AliITSRecoParam::GetHighFluxParam();
  }
}
 //___________________________________________________________________________
AliITSReconstructor::~AliITSReconstructor(){
// destructor
  delete fItsPID;
  if(fgkRecoParam) delete fgkRecoParam;
  if(fDetTypeRec) delete fDetTypeRec;
} 
//______________________________________________________________________
AliITSReconstructor::AliITSReconstructor(const AliITSReconstructor &ob) :AliReconstructor(ob),
									 fItsPID(ob.fItsPID),
									 fDetTypeRec(ob.fDetTypeRec)

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
void AliITSReconstructor::Init() {
    // Initalize this constructor bet getting/creating the objects
    // nesseary for a proper ITS reconstruction.
    // Inputs:
    //   none.
    // Output:
    //   none.
    // Return:
    //   none.

    AliITSInitGeometry initgeom;
    AliITSgeom *geom = initgeom.CreateAliITSgeom();
    AliInfo(Form("Geometry name: %s",(initgeom.GetGeometryName()).Data()));

    fDetTypeRec = new AliITSDetTypeRec();
    fDetTypeRec->SetITSgeom(geom);
    fDetTypeRec->SetDefaults();
    
    return;
}

//_____________________________________________________________________________
void AliITSReconstructor::Reconstruct(TTree *digitsTree, TTree *clustersTree) const
{
// reconstruct clusters

  TString option = GetOption();
  Bool_t clusfinder=kTRUE;   // Default: V2 cluster finder
  if(option.Contains("OrigCF"))clusfinder=kFALSE;

  fDetTypeRec->SetTreeAddressD(digitsTree);
  fDetTypeRec->MakeBranch(clustersTree,"R");
  fDetTypeRec->SetTreeAddressR(clustersTree);
  fDetTypeRec->DigitsToRecPoints(digitsTree,clustersTree,0,"All",clusfinder);    
}

//_________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRawReader* rawReader, TTree *clustersTree) const
{
  // reconstruct clusters from raw data
 
  fDetTypeRec->SetDefaultClusterFindersV2(kTRUE);
  fDetTypeRec->DigitsToRecPoints(rawReader,clustersTree);
}

//_____________________________________________________________________________
AliTracker* AliITSReconstructor::CreateTracker() const
{
// create a ITS tracker

  
  TString selectedTracker = GetOption();
  AliTracker* tracker;    
  if (selectedTracker.Contains("MI")) {
    tracker = new AliITStrackerMI(0);
  }  
  else if (selectedTracker.Contains("V2")) {
    tracker = new AliITStrackerV2(0);
  }
  else {
    tracker =  new AliITStrackerSA(0);  // inherits from AliITStrackerMI
    AliITStrackerSA *sat=(AliITStrackerSA*)tracker;
    if(selectedTracker.Contains("onlyITS"))sat->SetSAFlag(kTRUE);
    if(sat->GetSAFlag())AliDebug(1,"Tracking Performed in ITS only\n");
    if(selectedTracker.Contains("cosmics")||selectedTracker.Contains("COSMICS"))
      sat->SetOuterStartLayer(AliITSgeomTGeo::GetNLayers()-2);
  }

  TString selectedPIDmethod = GetOption();
  AliITSReconstructor* nc = const_cast<AliITSReconstructor*>(this);
  if(selectedPIDmethod.Contains("LandauFitPID")){
    Info("FillESD","ITS LandauFitPID option has been selected\n");
    nc->fItsPID = new AliITSpidESD2((AliITStrackerMI*)tracker);
  }
  else{
    Info("FillESD","ITS default PID\n");
    Double_t parITS[] = {79.,0.13, 10.}; //IB: this is  "pp tuning"
    nc->fItsPID = new AliITSpidESD1(parITS);
  }
 
  return tracker;
  
}

//_____________________________________________________________________________
AliVertexer* AliITSReconstructor::CreateVertexer() const
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
  if(selectedVertexer.Contains("3d") || selectedVertexer.Contains("3D")){
    Info("CreateVertexer","a AliITSVertexer3D object has been selected\n");
    return new AliITSVertexer3D("null");
  }
  if(selectedVertexer.Contains("cosmics") || selectedVertexer.Contains("COSMICS")){
    Info("CreateVertexer","a AliITSVertexerCosmics object has been selected\n");
    return new AliITSVertexerCosmics();
  }
  // by default an AliITSVertexerZ object is instatiated
  Info("CreateVertexer","a AliITSVertexerZ object has been selected\n");
  return new AliITSVertexerZ("null");
}

//_____________________________________________________________________________
void AliITSReconstructor::FillESD(TTree * /*digitsTree*/, TTree *clustersTree, 
				  AliESDEvent* esd) const
{
// make PID, find V0s and cascade
  if(fItsPID!=0) {
    TString selectedPIDmethod = GetOption();
    if(selectedPIDmethod.Contains("LandauFitPID")){
      fItsPID->MakePID(clustersTree,esd);
    }
    else{
      fItsPID->MakePID(esd);
    }
  }
  else {
    Error("FillESD","!! cannot do the PID !!\n");
  }
}
