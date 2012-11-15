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

/* $Id: AliITSUReconstructor.cxx 58442 2012-09-04 17:17:06Z masera $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ITS reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "AliITSUReconstructor.h"
#include "AliRun.h"
#include "AliRawReader.h"
#include "AliESDEvent.h"

#include "AliTracker.h"
#include "AliITStrackerMI.h"

#include "AliITSUGeomTGeo.h"
#include "AliITSUSegmentationPix.h"
#include "AliITSUDigitPix.h"
#include "AliITSUClusterizer.h"
#include "AliITSUClusterPix.h"
#include "AliMagF.h"

ClassImp(AliITSUReconstructor)

//___________________________________________________________________________
AliITSUReconstructor::AliITSUReconstructor() 
:  AliReconstructor()
  ,fGeom(0)
  ,fClusterFinders(0)
  ,fClusters(0)
{
  // Default constructor

}

//___________________________________________________________________________
AliITSUReconstructor::~AliITSUReconstructor()
{
  // destructor
  //
  if (!fGeom) return; // was not initialized
  //
  // same cluster finders and recpoint arrays might be attached to different layers
  for (int i=fGeom->GetNLayers();i--;) {
    TObject* clFinder = fClusterFinders.At(i);
    if (clFinder) {
      while (fClusterFinders.Remove(clFinder)) {}
      delete clFinder;
    }
    //
    delete[] fClusters;
  }
  //
  delete fGeom;
} 

//______________________________________________________________________
void AliITSUReconstructor::Init() 
{
  // Initalize this constructor 
  AliInfo("Initializing");
  if (fGeom) AliFatal("was already done, something is wrong...");
  //
  fGeom = new AliITSUGeomTGeo(kTRUE,kTRUE);
  AliITSUClusterPix::SetGeom(fGeom);
  //  
  AliITSUClusterizer* clusPIX = 0;
  fClusters = new TClonesArray*[fGeom->GetNLayers()];
  //
  for (int ilr=fGeom->GetNLayers();ilr--;) {
    fClusters[ilr] = 0;
    int tpDet = fGeom->GetLayerDetTypeID(ilr)/AliITSUGeomTGeo::kMaxSegmPerDetType;
    if (tpDet == AliITSUGeomTGeo::kDetTypePix) {
      if (!clusPIX)    clusPIX    = new AliITSUClusterizer();
      fClusterFinders.AddAtAndExpand(clusPIX, ilr);
      fClusters[ilr] = new TClonesArray(AliITSUClusterPix::Class());
      //
      // to expand the buffers to max.size
      clusPIX->SetSegmentation((AliITSUSegmentationPix*)fGeom->GetSegmentation(ilr)); 
      continue;
    }
    else {
      AliFatal(Form("ClusterFinder for detector type %d is not defined",tpDet));
    }
  }
  //
  return;

}

//_____________________________________________________________________________
void AliITSUReconstructor::Reconstruct(TTree *digitsTree, TTree *clustersTree) const
{
  // reconstruct clusters. If clustersTree is provided, write the tree
  if (!digitsTree) return;
  AliDebug(1,"ITSU Cluster finder (from digits tree) is initiated here \n");
  //
  // At the moment only pixel digits
  TClonesArray *digArrPix = 0;
  digitsTree->SetBranchAddress("ITSDigitsPix",&digArrPix);
  //
  // a new tree is created for each event: add each layer as separate branch
  TBranch *lrBranch[fGeom->GetNLayers()];
  //
  for (int ilr=0;ilr<fGeom->GetNLayers();ilr++) {
    if (clustersTree) { // do we write clusters tree?
      int tp = fGeom->GetLayerDetTypeID(ilr)/AliITSUGeomTGeo::kMaxSegmPerDetType;
      if (tp==AliITSUGeomTGeo::kDetTypePix) {
	lrBranch[ilr] = clustersTree->Bronch(Form("ITSRecPoints%d",ilr),"TClonesArray",&fClusters[ilr]);
      }
      else {
	AliFatal(Form("Detector type %d is not defined",tp));
      }
    }
  }
  //
  AliITSUClusterizer* clFinder = 0;
  AliMagF* field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
  double bz = 0;
  if (field == 0) AliError("Cannot get magnetic field from TGeoGlobalMagField");
  else bz = field->SolenoidField();
  //
  for (int ilr=0;ilr<fGeom->GetNLayers();ilr++) {
    //
    fClusters[ilr]->Clear();
    clFinder = (AliITSUClusterizer*)fClusterFinders[ilr];
    clFinder->SetSegmentation((AliITSUSegmentationPix*)fGeom->GetSegmentation(ilr));
    clFinder->SetLayerID(ilr);
    clFinder->SetClusters(fClusters[ilr]);
    clFinder->SetRecoParam(GetRecoParam()); // RS: Do we need to set it for every event?
    clFinder->PrepareLorentzAngleCorrection(bz);
    //
    int modF=fGeom->GetFirstModIndex(ilr);
    int modL=fGeom->GetLastModIndex(ilr)+1;
    for (int imod=modF;imod<modL;imod++) {
      digitsTree->GetEntry(imod);   
      int ndig  = digArrPix->GetEntries();
      if (!ndig) continue;
      clFinder->SetVolID(imod);
      clFinder->SetDigits(digArrPix);
      clFinder->Clusterize();
    }
    //
    AliITSUClusterPix::SetSortMode( AliITSUClusterPix::SortModeIdTrkYZ());
    fClusters[ilr]->Sort();
    AliDebug(1,Form(" -> Lr%d : %d Cluster",ilr,fClusters[ilr]->GetEntries()));
    if (clustersTree) lrBranch[ilr]->Fill();
  }
  if (clustersTree) clustersTree->SetEntries();
  //
}

//_____________________________________________________________________________
AliTracker* AliITSUReconstructor::CreateTracker() const
{
  // create a ITS tracker

  AliDebug(1,"ITSU tracking initialization will be done here\n");
 
  return 0;

  /* // from Current ITS
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
  */
  
}

//_____________________________________________________________________________
AliVertexer* AliITSUReconstructor::CreateVertexer() const
{
  // create a ITS vertexer
  // to be implemented for the upgrade

  AliDebug(1,"ITSU vertexer should be initiated here\n");
  return 0;

}

//_____________________________________________________________________________
AliTrackleter* AliITSUReconstructor::CreateMultFinder() const
{
  // create the SPD trackeleter for mult. reconstruction
  // to be implemented for the upgrade

  AliDebug(1,"ITSU MultFinder  should be initiated here\n");
  return 0;

}

//_____________________________________________________________________________
AliTracker* AliITSUReconstructor::CreateTrackleter() const
{
  // create the SPD trackeleter (for SPD PlaneEfficiency evaluation)
  // to be implemented for the upgrade

  AliDebug(1,"ITSU Trackleter  should be initiated here\n");
  return 0;

}

//_____________________________________________________________________________
Int_t AliITSUReconstructor::LoadClusters(TTree* treeRP) 
{
  // read clusters from the tree, if it is provided
  if (!treeRP) return 0;
  for (int ilr=fGeom->GetNLayers();ilr--;) {
    if (!fClusters[ilr]) AliFatal(Form("Clusters array for layer %d is not defined",ilr)); 
    TBranch* br = treeRP->GetBranch(Form("ITSRecPoints%d",ilr));
    if (!br) AliFatal(Form("Provided cluster tree does not contain branch for layer %d",ilr));
    br->SetAddress(&fClusters[ilr]);
  }
  treeRP->GetEntry(0); // we are still in 1 ev/tree mode...
  return 1;
}
