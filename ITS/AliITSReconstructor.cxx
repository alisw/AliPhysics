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
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITStrackerMI.h"
#include "AliITStrackerV2.h"
#include "AliITStrackerSA.h"
#include "AliITSVertexerFast.h"
#include "AliITSVertexerFixed.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
#include "AliITSVertexerCosmics.h"
#include "AliITSInitGeometry.h"
#include "AliITSTrackleterSPDEff.h"
#include "AliITSMultReconstructor.h"

ClassImp(AliITSReconstructor)

//___________________________________________________________________________
AliITSReconstructor::AliITSReconstructor() : AliReconstructor(),
fDetTypeRec(0)
{
  // Default constructor
}
 //___________________________________________________________________________
AliITSReconstructor::~AliITSReconstructor(){
// destructor
  if(fDetTypeRec) delete fDetTypeRec;
} 
//____________________________________________________________________________
void AliITSReconstructor::GetPidSettings(AliESDpid *ESDpid) {
  //
  // pass PID settings from AliITSRecoParam to AliESDpid
  //
  Int_t pidOpt = GetRecoParam()->GetPID();

  if(pidOpt==1){
    AliDebug(1,"ITS LandauFitPID option has been selected\n");
    ESDpid->SetITSPIDmethod(AliESDpid::kITSLikelihood);
  }
  else{
    AliDebug(1,"ITS default PID\n");
    ESDpid->SetITSPIDmethod(AliESDpid::kITSTruncMean);
  }
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

  Int_t cluFindOpt = GetRecoParam()->GetClusterFinder();

  fDetTypeRec->SetTreeAddressD(digitsTree);
  fDetTypeRec->MakeBranch(clustersTree,"R");
  fDetTypeRec->SetTreeAddressR(clustersTree);
  fDetTypeRec->DigitsToRecPoints(digitsTree,clustersTree,0,"All",cluFindOpt);    
}

//_________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRawReader* rawReader, TTree *clustersTree) const
{
  // reconstruct clusters from raw data
 
  fDetTypeRec->SetDefaultClusterFindersV2(kTRUE);
  fDetTypeRec->DigitsToRecPoints(rawReader,clustersTree);
}

//_____________________________________________________________________________
AliTrackleter* AliITSReconstructor::CreateMultFinder() const
{
  // create the SPD trackeleter for mult. reconstruction
  AliITSMultReconstructor* multReco = new AliITSMultReconstructor();
  multReco->SetHistOn(kFALSE);
  multReco->SetDetTypeRec(fDetTypeRec);
  return multReco;
}

//_____________________________________________________________________________
AliTracker* AliITSReconstructor::CreateTrackleter() const
{
// create the SPD trackeleter (for SPD PlaneEfficiency evaluation)
  if(!GetRecoParam()->GetComputePlaneEff() || !GetRecoParam()->GetUseTrackletsPlaneEff()) return NULL;
  //Int_t trackerOpt = GetRecoParam()->GetTracker();
  AliTracker* trackleter;
  trackleter = new AliITSTrackleterSPDEff();
  AliITSTrackleterSPDEff *spdtrackleter=(AliITSTrackleterSPDEff*)trackleter;
  // here set cuts (from RecoParam)
  if(GetRecoParam()->GetBkgTrackletsPlaneEff()) spdtrackleter->SetReflectClusterAroundZAxisForLayer(1,kTRUE);
  if(GetRecoParam()->GetMCTrackletsPlaneEff()) spdtrackleter->SetMC();
  spdtrackleter->SetHistOn();
  spdtrackleter->SetPhiWindowL2(GetRecoParam()->GetTrackleterPhiWindowL2());
  spdtrackleter->SetZetaWindowL2(GetRecoParam()->GetTrackleterZetaWindowL2());
  spdtrackleter->SetPhiWindowL1(GetRecoParam()->GetTrackleterPhiWindowL1());
  spdtrackleter->SetZetaWindowL1(GetRecoParam()->GetTrackleterZetaWindowL1());
  if(GetRecoParam()->GetUpdateOncePerEventPlaneEff()) spdtrackleter->SetUpdateOncePerEventPlaneEff();
  spdtrackleter->SetMinContVtx(GetRecoParam()->GetMinContVtxPlaneEff());
  return trackleter;
}

//_____________________________________________________________________________
AliTracker* AliITSReconstructor::CreateTracker() const
{
// create a ITS tracker

  Int_t trackerOpt = GetRecoParam()->GetTracker();
  AliTracker* tracker;    
  if (trackerOpt==1) {
    tracker = new AliITStrackerMI(0);
    AliITStrackerMI *mit=(AliITStrackerMI*)tracker;
    mit->SetDetTypeRec(fDetTypeRec);
  }  
  else if (trackerOpt==2) {
    tracker = new AliITStrackerV2(0);
  }
  else {
    tracker =  new AliITStrackerSA(0);  // inherits from AliITStrackerMI
    AliITStrackerSA *sat=(AliITStrackerSA*)tracker;
    sat->SetDetTypeRec(fDetTypeRec);
    if(GetRecoParam()->GetTrackerSAOnly()) sat->SetSAFlag(kTRUE);
    if(sat->GetSAFlag())AliDebug(1,"Tracking Performed in ITS only\n");
	if(GetRecoParam()->GetInwardFindingSA()){
       sat->SetInwardFinding();
       sat->SetInnerStartLayer(GetRecoParam()->GetInnerStartLayerSA());
	}else{
      sat->SetOutwardFinding();
      sat->SetOuterStartLayer(GetRecoParam()->GetOuterStartLayerSA());
	}
    sat->SetMinNPoints(GetRecoParam()->GetMinNPointsSA());
  }

  return tracker;
  
}

//_____________________________________________________________________________
AliVertexer* AliITSReconstructor::CreateVertexer() const
{
// create a ITS vertexer

  AliITSVertexer *vptr = NULL;
  Int_t vtxOpt = GetRecoParam()->GetVertexer();
  if(vtxOpt==3){
    AliFatal("Option AliITSVertexerIons is no longer supported");
    return vptr;
  }
  else if(vtxOpt==4){
    Double_t smear[3]={GetRecoParam()->GetVertexerFastSmearX(),
		       GetRecoParam()->GetVertexerFastSmearY(),
		       GetRecoParam()->GetVertexerFastSmearZ()};
    AliDebug(1,"AliITSVertexerFast has been selected"); 
    vptr = new AliITSVertexerFast(smear);
  }
  else if(vtxOpt==1){
    AliDebug(1,"AliITSVertexerZ has been selected");
    AliITSVertexerZ* vtxrz =  new AliITSVertexerZ();
    Int_t pileupAlgo=GetRecoParam()->GetSPDVertexerPileupAlgo();
    if(pileupAlgo==3) vtxrz->SetSearchForPileup(kFALSE);
    vptr = vtxrz;
  }
  else if(vtxOpt==2){
    AliDebug(1,"AliITSVertexerCosmics has been selected");
    vptr = new AliITSVertexerCosmics();
  }
  else if(vtxOpt==5){ 
    AliDebug(1,"Vertex is fixed in the position of the TDI\n");
    vptr = new AliITSVertexerFixed("TDI");
  }
  else if(vtxOpt==6){ 
    AliDebug(1,"Vertex is fixed in the position of the TED\n");
    vptr = new AliITSVertexerFixed("TED");
  }
  else {
  // by default an AliITSVertexer3D object is instatiated
    AliITSVertexer3D*  vtxr = new AliITSVertexer3D();
    Float_t dzw=GetRecoParam()->GetVertexer3DWideFiducialRegionZ();
    Float_t drw=GetRecoParam()->GetVertexer3DWideFiducialRegionR();
    vtxr->SetWideFiducialRegion(dzw,drw);
    Float_t dzn=GetRecoParam()->GetVertexer3DNarrowFiducialRegionZ();
    Float_t drn=GetRecoParam()->GetVertexer3DNarrowFiducialRegionR();
    vtxr->SetNarrowFiducialRegion(dzn,drn);
    Float_t dphil=GetRecoParam()->GetVertexer3DLooseDeltaPhiCut();
    Float_t dphit=GetRecoParam()->GetVertexer3DTightDeltaPhiCut();
    vtxr->SetDeltaPhiCuts(dphil,dphit);
    Float_t dcacut=GetRecoParam()->GetVertexer3DDCACut();
    vtxr->SetDCACut(dcacut);
    Int_t pileupAlgo=GetRecoParam()->GetSPDVertexerPileupAlgo();
    vtxr->SetPileupAlgo(pileupAlgo);
    UChar_t highmultAlgo=GetRecoParam()->GetSPDVertexerHighMultAlgo();
    vtxr->SetHighMultAlgo(highmultAlgo);
    AliDebug(1,Form("AliITSVertexer3D with pileup algo %d has been selected",pileupAlgo));
    vptr = vtxr;
  }
  vptr->SetDetTypeRec(fDetTypeRec);
  return vptr;
}

//_____________________________________________________________________________
void AliITSReconstructor::FillESD(TTree * /*digitsTree*/, TTree * /*clustersTree*/, 
				  AliESDEvent* /* esd */) const
{
// make PID, find V0s and cascade
/* Now done in AliESDpid
  if(fItsPID!=0) fItsPID->MakePID(esd);
  else AliError("!! cannot do the PID !!");
*/
}
