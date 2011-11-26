/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
// 
//        Implementation of the ITS-SPD trackleter class
//
// It retrieves clusters in the pixels (theta and phi) and finds tracklets.
// These can be used to extract charged particle multiplicity from the ITS.
//
// A tracklet consists of two ITS clusters, one in the first pixel layer and 
// one in the second. The clusters are associated if the differences in 
// Phi (azimuth) and Theta (polar angle) are within fiducial windows.
// In case of multiple candidates the candidate with minimum
// distance is selected. 
//
// Two methods return the number of tracklets and the number of unassociated 
// clusters (i.e. not used in any tracklet) in the first SPD layer
// (GetNTracklets and GetNSingleClusters)
//
// The cuts on phi and theta depend on the interacting system (p-p or Pb-Pb)
// and can be set via AliITSRecoParam class
// (SetPhiWindow and SetThetaWindow)  
// 
// Origin: Tiziano Virgili 
//
// Current support and development: 
//         Domenico Elia, Maria Nicassio (INFN Bari) 
//         Domenico.Elia@ba.infn.it, Maria.Nicassio@ba.infn.it
//
// Most recent updates:
//     - multiple association forbidden (fOnlyOneTrackletPerC2 = kTRUE)    
//     - phi definition changed to ALICE convention (0,2*TMath::pi()) 
//     - cluster coordinates taken with GetGlobalXYZ()
//     - fGeometry removed
//     - number of fired chips on the two layers
//     - option to cut duplicates in the overlaps
//     - options and fiducial cuts via AliITSRecoParam
//     - move from DeltaZeta to DeltaTheta cut
//     - update to the new algorithm by Mariella and Jan Fiete
//     - store also DeltaTheta in the ESD 
//     - less new and delete calls when creating the needed arrays
//
//     - RS: to decrease the number of new/deletes the clusters data are stored 
//           not in float[6] attached to float**, but in 1-D array.
//     - RS: Clusters are sorted in Z in roder to have the same numbering as in the ITS reco
//     - RS: Clusters used by ESDtrack are flagged, this information is passed to AliMulitiplicity object 
//           when storing the tracklets and single cluster info
//     - MN: first MC label of single clusters stored  
//_________________________________________________________________________

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBits.h>
#include <TArrayI.h>
#include <string.h>

#include "AliITSMultReconstructor.h"
#include "AliITSReconstructor.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeom.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDetTypeRec.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliLog.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliESDv0.h"
#include "AliV0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliRefArray.h"

//____________________________________________________________________
ClassImp(AliITSMultReconstructor)


//____________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor():
fDetTypeRec(0),fESDEvent(0),fTreeRP(0),fTreeRPMix(0),
fTracklets(0),
fSClusters(0),
fNTracklets(0),
fNSingleCluster(0),
fDPhiWindow(0),
fDThetaWindow(0),
fPhiShift(0),
fRemoveClustersFromOverlaps(0),
fPhiOverlapCut(0),
fZetaOverlapCut(0),
fPhiRotationAngle(0),
fScaleDTBySin2T(0),
fNStdDev(1.0),
fNStdDevSq(1.0),
//
fCutPxDrSPDin(0.1),
fCutPxDrSPDout(0.15),
fCutPxDz(0.2),
fCutDCArz(0.5),
fCutMinElectronProbTPC(0.5),
fCutMinElectronProbESD(0.1),
fCutMinP(0.05),
fCutMinRGamma(2.),
fCutMinRK0(1.),
fCutMinPointAngle(0.98),
fCutMaxDCADauther(0.5),
fCutMassGamma(0.03),
fCutMassGammaNSigma(5.),
fCutMassK0(0.03),
fCutMassK0NSigma(5.),
fCutChi2cGamma(2.),
fCutChi2cK0(2.),
fCutGammaSFromDecay(-10.),
fCutK0SFromDecay(-10.),
fCutMaxDCA(1.),
//
fHistOn(0),
fhClustersDPhiAcc(0),
fhClustersDThetaAcc(0),
fhClustersDPhiAll(0),
fhClustersDThetaAll(0),
fhDPhiVsDThetaAll(0),
fhDPhiVsDThetaAcc(0),
fhetaTracklets(0),
fhphiTracklets(0),
fhetaClustersLay1(0),
fhphiClustersLay1(0),
//
  fDPhiShift(0),
  fDPhiWindow2(0),
  fDThetaWindow2(0),
  fPartners(0),
  fAssociatedLay1(0),
  fMinDists(0),
  fBlackList(0),
//
  fCreateClustersCopy(0),
  fClustersLoaded(0),
  fRecoDone(0),
  fBuildRefs(kTRUE),
  fSPDSeg()
{
  // default c-tor
  for (int i=0;i<2;i++) {
    fNFiredChips[i] = 0;
    fClArr[i] = 0;
    for (int j=0;j<2;j++) fUsedClusLay[i][j] = 0;
    fDetectorIndexClustersLay[i] = 0;
    fOverlapFlagClustersLay[i] = 0;
    fNClustersLay[i] = 0;
    fClustersLay[i] = 0;
  }
  // Method to reconstruct the charged particles multiplicity with the 
  // SPD (tracklets).
  
  SetHistOn();

  if (AliITSReconstructor::GetRecoParam()) { 
    SetPhiWindow(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiWindow());
    SetThetaWindow(AliITSReconstructor::GetRecoParam()->GetTrackleterThetaWindow());
    SetPhiShift(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiShift());
    SetRemoveClustersFromOverlaps(AliITSReconstructor::GetRecoParam()->GetTrackleterRemoveClustersFromOverlaps());
    SetPhiOverlapCut(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiOverlapCut());
    SetZetaOverlapCut(AliITSReconstructor::GetRecoParam()->GetTrackleterZetaOverlapCut());
    SetPhiRotationAngle(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiRotationAngle());
    SetNStdDev(AliITSReconstructor::GetRecoParam()->GetTrackleterNStdDevCut());
    SetScaleDThetaBySin2T(AliITSReconstructor::GetRecoParam()->GetTrackleterScaleDThetaBySin2T());
    SetBuildRefs(AliITSReconstructor::GetRecoParam()->GetTrackleterBuildCl2TrkRefs());
    //
    SetCutPxDrSPDin(AliITSReconstructor::GetRecoParam()->GetMultCutPxDrSPDin());
    SetCutPxDrSPDout(AliITSReconstructor::GetRecoParam()->GetMultCutPxDrSPDout());
    SetCutPxDz(AliITSReconstructor::GetRecoParam()->GetMultCutPxDz());
    SetCutDCArz(AliITSReconstructor::GetRecoParam()->GetMultCutDCArz());
    SetCutMinElectronProbTPC(AliITSReconstructor::GetRecoParam()->GetMultCutMinElectronProbTPC());
    SetCutMinElectronProbESD(AliITSReconstructor::GetRecoParam()->GetMultCutMinElectronProbESD());
    SetCutMinP(AliITSReconstructor::GetRecoParam()->GetMultCutMinP());
    SetCutMinRGamma(AliITSReconstructor::GetRecoParam()->GetMultCutMinRGamma());
    SetCutMinRK0(AliITSReconstructor::GetRecoParam()->GetMultCutMinRK0());
    SetCutMinPointAngle(AliITSReconstructor::GetRecoParam()->GetMultCutMinPointAngle());
    SetCutMaxDCADauther(AliITSReconstructor::GetRecoParam()->GetMultCutMaxDCADauther());
    SetCutMassGamma(AliITSReconstructor::GetRecoParam()->GetMultCutMassGamma());
    SetCutMassGammaNSigma(AliITSReconstructor::GetRecoParam()->GetMultCutMassGammaNSigma());
    SetCutMassK0(AliITSReconstructor::GetRecoParam()->GetMultCutMassK0());
    SetCutMassK0NSigma(AliITSReconstructor::GetRecoParam()->GetMultCutMassK0NSigma());
    SetCutChi2cGamma(AliITSReconstructor::GetRecoParam()->GetMultCutChi2cGamma());
    SetCutChi2cK0(AliITSReconstructor::GetRecoParam()->GetMultCutChi2cK0());
    SetCutGammaSFromDecay(AliITSReconstructor::GetRecoParam()->GetMultCutGammaSFromDecay());
    SetCutK0SFromDecay(AliITSReconstructor::GetRecoParam()->GetMultCutK0SFromDecay());
    SetCutMaxDCA(AliITSReconstructor::GetRecoParam()->GetMultCutMaxDCA());
    //
  } else {
    SetPhiWindow();
    SetThetaWindow();
    SetPhiShift();
    SetRemoveClustersFromOverlaps();
    SetPhiOverlapCut();
    SetZetaOverlapCut();
    SetPhiRotationAngle();

    //
    SetCutPxDrSPDin();
    SetCutPxDrSPDout();
    SetCutPxDz();
    SetCutDCArz();
    SetCutMinElectronProbTPC();
    SetCutMinElectronProbESD();
    SetCutMinP();
    SetCutMinRGamma();
    SetCutMinRK0();
    SetCutMinPointAngle();
    SetCutMaxDCADauther();
    SetCutMassGamma();
    SetCutMassGammaNSigma();
    SetCutMassK0();
    SetCutMassK0NSigma();
    SetCutChi2cGamma();
    SetCutChi2cK0();
    SetCutGammaSFromDecay();
    SetCutK0SFromDecay();
    SetCutMaxDCA();
  } 
  //
  fTracklets                 = 0;
  fSClusters                 = 0;
  //
  // definition of histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fhClustersDPhiAcc   = new TH1F("dphiacc",  "dphi",  100,-0.1,0.1);
  fhClustersDThetaAcc = new TH1F("dthetaacc","dtheta",100,-0.1,0.1);

  fhDPhiVsDThetaAcc = new TH2F("dphiVsDthetaAcc","",100,-0.1,0.1,100,-0.1,0.1);

  fhClustersDPhiAll   = new TH1F("dphiall",  "dphi",  100,0.0,0.5);
  fhClustersDThetaAll = new TH1F("dthetaall","dtheta",100,0.0,0.5);

  fhDPhiVsDThetaAll = new TH2F("dphiVsDthetaAll","",100,0.,0.5,100,0.,0.5);

  fhetaTracklets  = new TH1F("etaTracklets",  "eta",  100,-2.,2.);
  fhphiTracklets  = new TH1F("phiTracklets",  "phi",  100, 0., 2*TMath::Pi());
  fhetaClustersLay1  = new TH1F("etaClustersLay1",  "etaCl1",  100,-2.,2.);
  fhphiClustersLay1  = new TH1F("phiClustersLay1", "phiCl1", 100, 0., 2*TMath::Pi());
  for (int i=2;i--;) fStoreRefs[i][0] =  fStoreRefs[i][1] = kFALSE;
  TH1::AddDirectory(oldStatus);
}

//______________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor(const AliITSMultReconstructor &mr) : 
AliTrackleter(mr),
fDetTypeRec(0),fESDEvent(0),fTreeRP(0),fTreeRPMix(0),
fTracklets(0),
fSClusters(0),
fNTracklets(0),
fNSingleCluster(0),
fDPhiWindow(0),
fDThetaWindow(0),
fPhiShift(0),
fRemoveClustersFromOverlaps(0),
fPhiOverlapCut(0),
fZetaOverlapCut(0),
fPhiRotationAngle(0),
fScaleDTBySin2T(0),
fNStdDev(1.0),
fNStdDevSq(1.0),
//
fCutPxDrSPDin(0.1),
fCutPxDrSPDout(0.15),
fCutPxDz(0.2),
fCutDCArz(0.5),
fCutMinElectronProbTPC(0.5),
fCutMinElectronProbESD(0.1),
fCutMinP(0.05),
fCutMinRGamma(2.),
fCutMinRK0(1.),
fCutMinPointAngle(0.98),
fCutMaxDCADauther(0.5),
fCutMassGamma(0.03),
fCutMassGammaNSigma(5.),
fCutMassK0(0.03),
fCutMassK0NSigma(5.),
fCutChi2cGamma(2.),
fCutChi2cK0(2.),
fCutGammaSFromDecay(-10.),
fCutK0SFromDecay(-10.),
fCutMaxDCA(1.),
//
fHistOn(0),
fhClustersDPhiAcc(0),
fhClustersDThetaAcc(0),
fhClustersDPhiAll(0),
fhClustersDThetaAll(0),
fhDPhiVsDThetaAll(0),
fhDPhiVsDThetaAcc(0),
fhetaTracklets(0),
fhphiTracklets(0),
fhetaClustersLay1(0),
fhphiClustersLay1(0),
fDPhiShift(0),
fDPhiWindow2(0),
fDThetaWindow2(0),
fPartners(0),
fAssociatedLay1(0),
fMinDists(0),
fBlackList(0),
//
fCreateClustersCopy(0),
fClustersLoaded(0),
fRecoDone(0),
fBuildRefs(kTRUE),
fSPDSeg()
 {
  // Copy constructor :!!! RS ATTENTION: old c-tor reassigned the pointers instead of creating a new copy -> would crash on delete
   AliError("May not use");
}

//______________________________________________________________________
AliITSMultReconstructor& AliITSMultReconstructor::operator=(const AliITSMultReconstructor& mr){
  // Assignment operator
  if (this != &mr) {
    this->~AliITSMultReconstructor();
    new(this) AliITSMultReconstructor(mr);
  }
  return *this;
}

//______________________________________________________________________
AliITSMultReconstructor::~AliITSMultReconstructor(){
  // Destructor

  // delete histograms
  delete fhClustersDPhiAcc;
  delete fhClustersDThetaAcc;
  delete fhClustersDPhiAll;
  delete fhClustersDThetaAll;
  delete fhDPhiVsDThetaAll;
  delete fhDPhiVsDThetaAcc;
  delete fhetaTracklets;
  delete fhphiTracklets;
  delete fhetaClustersLay1;
  delete fhphiClustersLay1;
  //
  // delete arrays    
  for(Int_t i=0; i<fNTracklets; i++) delete [] fTracklets[i];
    
  for(Int_t i=0; i<fNSingleCluster; i++) delete [] fSClusters[i];
  
  // 
  for (int i=0;i<2;i++) {
    delete[] fClustersLay[i];
    delete[] fDetectorIndexClustersLay[i];
    delete[] fOverlapFlagClustersLay[i];
    delete   fClArr[i];
    for (int j=0;j<2;j++) if (fUsedClusLay[i][j]) delete fUsedClusLay[i][j];
  }
  delete [] fTracklets;
  delete [] fSClusters;
  //
  delete[] fPartners;      fPartners = 0;
  delete[] fMinDists;      fMinDists = 0;
  delete   fBlackList;     fBlackList = 0;
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::Reconstruct(AliESDEvent* esd, TTree* treeRP) 
{  
  if (!treeRP) { AliError(" Invalid ITS cluster tree !\n"); return; }
  if (!esd) {AliError("ESDEvent is not available, use old reconstructor"); return;}
  // reset counters
  if (fMult) delete fMult; fMult = 0;
  fNClustersLay[0] = 0;
  fNClustersLay[1] = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;
  //
  fESDEvent = esd;
  fTreeRP = treeRP;
  //
  // >>>> RS: this part is equivalent to former AliITSVertexer::FindMultiplicity
  //
  // see if there is a SPD vertex 
  Bool_t isVtxOK=kTRUE, isCosmics=kFALSE;
  AliESDVertex* vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexSPD();
  if (!vtx || vtx->GetNContributors()<1) isVtxOK = kFALSE;
  if (vtx && strstr(vtx->GetTitle(),"cosmics")) {
    isVtxOK = kFALSE;
    isCosmics = kTRUE;
  }
  //
  if (!isVtxOK) {
    if (!isCosmics) {
      AliDebug(1,"Tracklets multiplicity not determined because the primary vertex was not found");
      AliDebug(1,"Just counting the number of cluster-fired chips on the SPD layers");
    }
    vtx = 0;
  }
  if(vtx){
    float vtxf[3] = {vtx->GetX(),vtx->GetY(),vtx->GetZ()};
    FindTracklets(vtxf);
  }
  else {
    FindTracklets(0);
  }
  //
  CreateMultiplicityObject();
}

//____________________________________________________________________
void AliITSMultReconstructor::Reconstruct(TTree* clusterTree, Float_t* vtx, Float_t* /* vtxRes*/) {
  //
  // RS NOTE - this is old reconstructor invocation, to be used from VertexFinder and in analysis mode

  if (fMult) delete fMult; fMult = 0;
  fNClustersLay[0] = 0;
  fNClustersLay[1] = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;
  //
  if (!clusterTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  //
  fESDEvent = 0;
  SetTreeRP(clusterTree);
  //
  FindTracklets(vtx);
  //
}


//____________________________________________________________________
void AliITSMultReconstructor::ReconstructMix(TTree* clusterTree, TTree* clusterTreeMix, const Float_t* vtx, Float_t*) 
{
  //
  // RS NOTE - this is old reconstructor invocation, to be used from VertexFinder and in analysis mode

  if (fMult) delete fMult; fMult = 0;
  fNClustersLay[0] = 0;
  fNClustersLay[1] = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;
  //
  if (!clusterTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  if (!clusterTreeMix) { AliError(" Invalid ITS cluster tree 2nd event !\n"); return; }
  //
  fESDEvent = 0;
  SetTreeRP(clusterTree);
  SetTreeRPMix(clusterTreeMix);
  //
  FindTracklets(vtx);
  //
}


//____________________________________________________________________
void AliITSMultReconstructor::FindTracklets(const Float_t *vtx) 
{
  // - calls LoadClusterArrays that finds the position of the clusters
  //   (in global coord) 

  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). Clusters in the inner layer can be now
  //   rotated for combinatorial studies 
  // - makes an array of tracklets 
  //   
  // After this method has been called, the clusters of the two layers
  // and the tracklets can be retrieved by calling the Get'er methods.


  // Find tracklets converging to vertex
  //
  LoadClusterArrays(fTreeRP,fTreeRPMix);
  // flag clusters used by ESD tracks
  if (fESDEvent) ProcessESDTracks();
  fRecoDone = kTRUE;

  if (!vtx) return;

  InitAux();
  
  // find the tracklets
  AliDebug(1,"Looking for tracklets... ");  
  
  ClusterPos2Angles(vtx); // convert cluster position to angles wrt vtx
  //
  // Step1: find all tracklets allowing double assocation: 
  int found = 1;
  while (found > 0) {
    found = 0;
    for (Int_t iC1=0; iC1<fNClustersLay[0]; iC1++) found += AssociateClusterOfL1(iC1);
  }
  //
  // Step2: store tracklets; remove used clusters 
  for (Int_t iC2=0; iC2<fNClustersLay[1]; iC2++) StoreTrackletForL2Cluster(iC2);
  //
  // store unused single clusters of L1
  StoreL1Singles();
  //
  AliDebug(1,Form("%d tracklets found", fNTracklets));
}

//____________________________________________________________________
void AliITSMultReconstructor::CreateMultiplicityObject()
{
  // create AliMultiplicity object and store it in the ESD event
  //
  TBits fastOrFiredMap,firedChipMap;
  if (fDetTypeRec) {
   fastOrFiredMap  = fDetTypeRec->GetFastOrFiredMap();
   firedChipMap    = fDetTypeRec->GetFiredChipMap(fTreeRP);
  }
  //
  fMult = new AliMultiplicity(fNTracklets,fNSingleCluster,fNFiredChips[0],fNFiredChips[1],fastOrFiredMap);
  fMult->SetMultTrackRefs( fBuildRefs );
  // store some details of reco:
  fMult->SetScaleDThetaBySin2T(fScaleDTBySin2T);
  fMult->SetDPhiWindow2(fDPhiWindow2);
  fMult->SetDThetaWindow2(fDThetaWindow2);
  fMult->SetDPhiShift(fDPhiShift);
  fMult->SetNStdDev(fNStdDev);
  //
  fMult->SetFiredChipMap(firedChipMap);
  AliITSRecPointContainer* rcont = AliITSRecPointContainer::Instance();
  fMult->SetITSClusters(0,rcont->GetNClustersInLayer(1,fTreeRP));
  for(Int_t kk=2;kk<=6;kk++) fMult->SetITSClusters(kk-1,rcont->GetNClustersInLayerFast(kk));
  //
  UInt_t shared[100]; 
  AliRefArray *refs[2][2] = {{0,0},{0,0}};
  if (fBuildRefs) {
    for (int il=2;il--;) 
      for (int it=2;it--;)  // tracklet_clusters->track references to stor
	if (fStoreRefs[il][it]) refs[il][it] = new AliRefArray(fNTracklets,0);
  }
  //
  for (int i=fNTracklets;i--;)  {
    float* tlInfo = fTracklets[i];
    fMult->SetTrackletData(i,tlInfo);
    //
    if (!fBuildRefs) continue; // do we need references?
    for (int itp=0;itp<2;itp++) {	
      for (int ilr=0;ilr<2;ilr++) {
	if (!fStoreRefs[ilr][itp]) continue; // nothing to store
	int clID = int(tlInfo[ilr ? kClID2:kClID1]);
	int nref = fUsedClusLay[ilr][itp]->GetReferences(clID,shared,100);
	if (!nref) continue;
	else if (nref==1) refs[ilr][itp]->AddReference(i,shared[0]);
	else refs[ilr][itp]->AddReferences(i,shared,nref);
      }
    }
  }
  if (fBuildRefs) fMult->AttachTracklet2TrackRefs(refs[0][0],refs[0][1],refs[1][0],refs[1][1]); 
  //
  AliRefArray *refsc[2] = {0,0};
  if (fBuildRefs) for (int it=2;it--;) if (fStoreRefs[0][it]) refsc[it] = new AliRefArray(fNClustersLay[0]);
  for (int i=fNSingleCluster;i--;) {
    float* clInfo = fSClusters[i];
    fMult->SetSingleClusterData(i,clInfo); 
    //
    if (!fBuildRefs) continue; // do we need references?
    int clID = int(clInfo[kSCID]);
    for (int itp=0;itp<2;itp++) {
      if (!fStoreRefs[0][itp]) continue;
      int nref = fUsedClusLay[0][itp]->GetReferences(clID,shared,100);
      if (!nref) continue;
      else if (nref==1) refsc[itp]->AddReference(i,shared[0]);
      else refsc[itp]->AddReferences(i,shared,nref);
    }
  }
  if (fBuildRefs) fMult->AttachCluster2TrackRefs(refsc[0],refsc[1]); 
  fMult->CompactBits();
  //
}


//____________________________________________________________________
void AliITSMultReconstructor::LoadClusterArrays(TTree* tree, TTree* treeMix)
{
  // load cluster info and prepare tracklets arrays
  //
  if (AreClustersLoaded()) {AliInfo("Clusters are already loaded"); return;}
  LoadClusterArrays(tree,0);
  LoadClusterArrays(treeMix ? treeMix:tree,1);
  int nmaxT = TMath::Min(fNClustersLay[0], fNClustersLay[1]);
  if (fTracklets) delete[] fTracklets;
  fTracklets = new Float_t*[nmaxT];
  memset(fTracklets,0,nmaxT*sizeof(Float_t*));
  //
  if (fSClusters) delete[] fSClusters;
  fSClusters = new Float_t*[fNClustersLay[0]]; 
  memset(fSClusters,0,fNClustersLay[0]*sizeof(Float_t*));
  //
  AliDebug(1,Form("(clusters in layer 1 : %d,  layer 2: %d)",fNClustersLay[0],fNClustersLay[1]));
  AliDebug(1,Form("(cluster-fired chips in layer 1 : %d,  layer 2: %d)",fNFiredChips[0],fNFiredChips[1]));
  SetClustersLoaded();
}

//____________________________________________________________________
void AliITSMultReconstructor::LoadClusterArrays(TTree* itsClusterTree, int il) 
{
  // This method
  // - gets the clusters from the cluster tree for layer il
  // - convert them into global coordinates 
  // - store them in the internal arrays
  // - count the number of cluster-fired chips
  //
  // RS: This method was strongly modified wrt original. In order to have the same numbering 
  // of clusters as in the ITS reco I had to introduce sorting in Z
  // Also note that now the clusters data are stored not in float[6] attached to float**, but in 1-D array
  AliDebug(1,Form("Loading clusters and cluster-fired chips for layer %d",il));
  //
  fNClustersLay[il] = 0;
  fNFiredChips[il]  = 0;
  for (int i=2;i--;) fStoreRefs[il][i] = kFALSE;
  //
  AliITSRecPointContainer* rpcont = 0;
  static TClonesArray statITSrec("AliITSRecPoint");
  static TObjArray clArr(100);
  TBranch* branch = 0;
  TClonesArray* itsClusters = 0;
  //
  if (!fCreateClustersCopy) {
    rpcont=AliITSRecPointContainer::Instance();
    itsClusters = rpcont->FetchClusters(0,itsClusterTree);
    if(!rpcont->IsSPDActive()){
      AliWarning("No SPD rec points found, multiplicity not calculated");
      return;
    } 
  }
  else {
    itsClusters = &statITSrec;
    branch = itsClusterTree->GetBranch("ITSRecPoints");
    branch->SetAddress(&itsClusters);
    if (!fClArr[il]) fClArr[il] = new TClonesArray("AliITSRecPoint",100);
  }    
  //
  // count clusters
  // loop over the SPD subdetectors
  int nclLayer = 0;
  int detMin = TMath::Max(0,AliITSgeomTGeo::GetModuleIndex(il+1,1,1));
  int detMax = AliITSgeomTGeo::GetModuleIndex(il+2,1,1);
  for (int idt=detMin;idt<detMax;idt++) {
    if (!fCreateClustersCopy) itsClusters = rpcont->UncheckedGetClusters(idt);
    else                      branch->GetEvent(idt); 
    int nClusters = itsClusters->GetEntriesFast();
    if (!nClusters) continue;
    Int_t nClustersInChip[5] = {0,0,0,0,0};
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
      if (!cluster) continue;
      if (fCreateClustersCopy) 	cluster = new ((*fClArr[il])[nclLayer]) AliITSRecPoint(*cluster);
      clArr.AddAtAndExpand(cluster,nclLayer++);
      Int_t chipNo = fSPDSeg.GetChipFromLocal(0,cluster->GetDetLocalZ());
      if(chipNo>=0)nClustersInChip[ chipNo ]++; 
    }
    for(Int_t ifChip=5;ifChip--;) if (nClustersInChip[ifChip]) fNFiredChips[il]++;
  }
  // sort the clusters in Z (to have the same numbering as in ITS reco
  Float_t *z     = new Float_t[nclLayer];
  Int_t   *index = new Int_t[nclLayer];
  for (int ic=0;ic<nclLayer;ic++) z[ic] = ((AliITSRecPoint*)clArr[ic])->GetZ();
  TMath::Sort(nclLayer,z,index,kFALSE);
  Float_t*   clustersLay              = new Float_t[nclLayer*kClNPar];
  Int_t*     detectorIndexClustersLay = new Int_t[nclLayer];
  Bool_t*    overlapFlagClustersLay   = new Bool_t[nclLayer];
  //
  for (int ic=0;ic<nclLayer;ic++) {
    AliITSRecPoint* cluster = (AliITSRecPoint*)clArr[index[ic]];
    float* clPar = &clustersLay[ic*kClNPar];
    //      
    cluster->GetGlobalXYZ( clPar );
    detectorIndexClustersLay[ic] = cluster->GetDetectorIndex(); 
    overlapFlagClustersLay[ic]   = kFALSE;
    for (Int_t i=3;i--;) clPar[kClMC0+i] = cluster->GetLabel(i);
  }
  clArr.Clear();
  delete[] z;
  delete[] index;
  //
  if (fOverlapFlagClustersLay[il]) delete[] fOverlapFlagClustersLay[il];
  fOverlapFlagClustersLay[il]   = overlapFlagClustersLay;
  //
  if (fDetectorIndexClustersLay[il]) delete[] fDetectorIndexClustersLay[il]; 
  fDetectorIndexClustersLay[il] = detectorIndexClustersLay;
  //
  if (fBuildRefs) {
    for (int it=0;it<2;it++) {
      if (fUsedClusLay[il][it]) delete fUsedClusLay[il][it];
      fUsedClusLay[il][it] = new AliRefArray(nclLayer);
    }
  }
  //
  if (fClustersLay[il]) delete[] fClustersLay[il]; 
  fClustersLay[il] = clustersLay;
  fNClustersLay[il] = nclLayer;
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::LoadClusterFiredChips(TTree* itsClusterTree) {
  // This method    
  // - gets the clusters from the cluster tree 
  // - counts the number of (cluster)fired chips
  
  AliDebug(1,"Loading cluster-fired chips ...");
  
  fNFiredChips[0] = 0;
  fNFiredChips[1] = 0;
  
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray* itsClusters=NULL;
  rpcont->FetchClusters(0,itsClusterTree);
  if(!rpcont->IsSPDActive()){
    AliWarning("No SPD rec points found, multiplicity not calculated");
    return;
  } 

  // loop over the its subdetectors
  Int_t nSPDmodules=AliITSgeomTGeo::GetModuleIndex(3,1,1);
  for (Int_t iIts=0; iIts < nSPDmodules; iIts++) {
    itsClusters=rpcont->UncheckedGetClusters(iIts);
    Int_t nClusters = itsClusters->GetEntriesFast();

    // number of clusters in each chip of the current module
    Int_t nClustersInChip[5] = {0,0,0,0,0};
    Int_t layer = 0;
    Int_t ladder=0;
    Int_t det=0;
    AliITSgeomTGeo::GetModuleId(iIts,layer,ladder,det);
    --layer;  // layer is from 1 to 6 in AliITSgeomTGeo, but from 0 to 5 here
    if(layer<0 || layer >1)continue;
    
    // loop over clusters
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
          
      // find the chip for the current cluster
      Float_t locz = cluster->GetDetLocalZ();
      Int_t iChip = fSPDSeg.GetChipFromLocal(0,locz);
      if (iChip>=0) nClustersInChip[iChip]++; 
      
    }// end of cluster loop

    // get number of fired chips in the current module
    for(Int_t ifChip=0; ifChip<5; ifChip++) {
      if(nClustersInChip[ifChip] >= 1)  fNFiredChips[layer]++;
    }

  } // end of its "subdetector" loop  
  

  AliDebug(1,Form("(cluster-fired chips in layer 1 : %d,  layer 2: %d)",fNFiredChips[0],fNFiredChips[1]));
}
//____________________________________________________________________
void
AliITSMultReconstructor::SaveHists() {
  // This method save the histograms on the output file
  // (only if fHistOn is TRUE). 
  
  if (!fHistOn)
    return;

  fhClustersDPhiAll->Write();
  fhClustersDThetaAll->Write();
  fhDPhiVsDThetaAll->Write();

  fhClustersDPhiAcc->Write();
  fhClustersDThetaAcc->Write();
  fhDPhiVsDThetaAcc->Write();

  fhetaTracklets->Write();
  fhphiTracklets->Write();
  fhetaClustersLay1->Write();
  fhphiClustersLay1->Write();
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagClustersInOverlapRegions (Int_t iC1, Int_t iC2WithBestDist) 
{
  // Flags clusters in the overlapping regions
  Float_t distClSameMod=0.;
  Float_t distClSameModMin=0.;
  Int_t   iClOverlap =0;
  Float_t meanRadiusLay1 = 3.99335; // average radius inner layer
  Float_t meanRadiusLay2 = 7.37935; // average radius outer layer;

  Float_t zproj1=0.;
  Float_t zproj2=0.;
  Float_t deZproj=0.;
  Float_t* clPar1  = GetClusterLayer1(iC1);
  Float_t* clPar2B = GetClusterLayer2(iC2WithBestDist);
  // Loop on inner layer clusters
  for (Int_t iiC1=0; iiC1<fNClustersLay[0]; iiC1++) {
    if (!fOverlapFlagClustersLay[0][iiC1]) {
      // only for adjacent modules
      if ((TMath::Abs(fDetectorIndexClustersLay[0][iC1]-fDetectorIndexClustersLay[0][iiC1])==4)||
         (TMath::Abs(fDetectorIndexClustersLay[0][iC1]-fDetectorIndexClustersLay[0][iiC1])==76)) {
	Float_t *clPar11 = GetClusterLayer1(iiC1);
        Float_t dePhi=TMath::Abs(clPar11[kClPh]-clPar1[kClPh]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay1/TMath::Tan(clPar1[kClTh]);
        zproj2=meanRadiusLay1/TMath::Tan(clPar11[kClTh]);

        deZproj=TMath::Abs(zproj1-zproj2);

        distClSameMod = TMath::Sqrt(TMath::Power(deZproj/fZetaOverlapCut,2)+TMath::Power(dePhi/fPhiOverlapCut,2));
        if (distClSameMod<=1.) fOverlapFlagClustersLay[0][iiC1]=kTRUE;

//        if (distClSameMod<=1.) {
//          if (distClSameModMin==0. || distClSameMod<distClSameModMin) {
//            distClSameModMin=distClSameMod;
//            iClOverlap=iiC1;
//          } 
//        }


      } // end adjacent modules
    } 
  } // end Loop on inner layer clusters

//  if (distClSameModMin!=0.) fOverlapFlagClustersLay[0][iClOverlap]=kTRUE;

  distClSameMod=0.;
  distClSameModMin=0.;
  iClOverlap =0;
  // Loop on outer layer clusters
  for (Int_t iiC2=0; iiC2<fNClustersLay[1]; iiC2++) {
    if (!fOverlapFlagClustersLay[1][iiC2]) {
      // only for adjacent modules
      Float_t *clPar2 = GetClusterLayer2(iiC2);
      if ((TMath::Abs(fDetectorIndexClustersLay[1][iC2WithBestDist]-fDetectorIndexClustersLay[1][iiC2])==4) ||
         (TMath::Abs(fDetectorIndexClustersLay[1][iC2WithBestDist]-fDetectorIndexClustersLay[1][iiC2])==156)) {
        Float_t dePhi=TMath::Abs(clPar2[kClPh]-clPar2B[kClPh]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay2/TMath::Tan(clPar2B[kClTh]);
        zproj2=meanRadiusLay2/TMath::Tan(clPar2[kClTh]);

        deZproj=TMath::Abs(zproj1-zproj2);
        distClSameMod = TMath::Sqrt(TMath::Power(deZproj/fZetaOverlapCut,2)+TMath::Power(dePhi/fPhiOverlapCut,2));
        if (distClSameMod<=1.) fOverlapFlagClustersLay[1][iiC2]=kTRUE;

//        if (distClSameMod<=1.) {
//          if (distClSameModMin==0. || distClSameMod<distClSameModMin) {
//            distClSameModMin=distClSameMod;
//            iClOverlap=iiC2;
//          }
//        }

      } // end adjacent modules
    }
  } // end Loop on outer layer clusters

//  if (distClSameModMin!=0.) fOverlapFlagClustersLay[1][iClOverlap]=kTRUE;

}

//____________________________________________________________________
void AliITSMultReconstructor::InitAux()
{
  // init arrays/parameters for tracklet reconstruction
  
  // dPhi shift is field dependent, get average magnetic field
  Float_t bz = 0;
  AliMagF* field = 0;
  if (TGeoGlobalMagField::Instance()) field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
  if (!field) {
    AliError("Could not retrieve magnetic field. Assuming no field. Delta Phi shift will be deactivated in AliITSMultReconstructor.");
  }
  else bz = TMath::Abs(field->SolenoidField());
  fDPhiShift = fPhiShift / 5 * bz; 
  AliDebug(1, Form("Using phi shift of %f", fDPhiShift));
  //
  if (fPartners) delete[] fPartners; fPartners = new Int_t[fNClustersLay[1]];
  if (fMinDists) delete[] fMinDists; fMinDists = new Float_t[fNClustersLay[1]];
  if (fAssociatedLay1) delete[] fAssociatedLay1; fAssociatedLay1 = new Int_t[fNClustersLay[0]];
  //
  if (fBlackList) delete fBlackList; fBlackList = new AliRefArray(fNClustersLay[0]);
  //
  //  Printf("Vertex in find tracklets...%f %f %f",vtx[0],vtx[1],vtx[2]);
  for (Int_t i=0; i<fNClustersLay[1]; i++) {
    fPartners[i] = -1;
    fMinDists[i] = 2*fNStdDev;
  }
  memset(fAssociatedLay1,0,fNClustersLay[0]*sizeof(Int_t));
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::ClusterPos2Angles(const Float_t *vtx)
{
  // convert cluster coordinates to angles wrt vertex
  for (int ilr=0;ilr<2;ilr++) {
    for (Int_t iC=0; iC<fNClustersLay[ilr]; iC++) {    
      float* clPar = GetClusterOfLayer(ilr,iC);
      CalcThetaPhi(clPar[kClTh]-vtx[0],clPar[kClPh]-vtx[1],clPar[kClZ]-vtx[2],clPar[kClTh],clPar[kClPh]);
      if (ilr==0) {
	clPar[kClPh] = clPar[kClPh] + fPhiRotationAngle;   // rotation of inner layer for comb studies  
	if (fHistOn) {
	  Float_t eta = clPar[kClTh];
	  eta= TMath::Tan(eta/2.);
	  eta=-TMath::Log(eta);
	  fhetaClustersLay1->Fill(eta);    
	  fhphiClustersLay1->Fill(clPar[kClPh]);
	}
      }      
    }
  }
  //
}

//____________________________________________________________________
Int_t AliITSMultReconstructor::AssociateClusterOfL1(Int_t iC1)
{
  // search association of cluster iC1 of L1 with all clusters of L2
  if (fAssociatedLay1[iC1] != 0) return 0;
  Int_t  iC2WithBestDist = -1;   // reset
  Double_t minDist       =  2*fNStdDev;   // reset
  float* clPar1 = GetClusterLayer1(iC1);
  for (Int_t iC2=0; iC2<fNClustersLay[1]; iC2++) {
    //
    if (fBlackList->IsReferred(iC1,iC2)) continue;
    float* clPar2 = GetClusterLayer2(iC2);
    //
    // find the difference in angles
    Double_t dTheta = TMath::Abs(clPar2[kClTh] - clPar1[kClTh]); 
    Double_t dPhi   = TMath::Abs(clPar2[kClPh] - clPar1[kClPh]);
    //        Printf("detheta %f  dephi %f", dTheta,dPhi);
    //
    if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;     // take into account boundary condition
    //
    if (fHistOn) {
      fhClustersDPhiAll->Fill(dPhi);
      fhClustersDThetaAll->Fill(dTheta);    
      fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
    }
    Float_t d = CalcDist(dPhi,dTheta,clPar1[kClTh]);     // make "elliptical" cut in Phi and Theta! 
    // look for the minimum distance: the minimum is in iC2WithBestDist
    if (d<fNStdDev && d<minDist) { minDist=d; iC2WithBestDist = iC2; }
  }
  //
  if (minDist<fNStdDev) { // This means that a cluster in layer 2 was found that matches with iC1
    //
    if (fMinDists[iC2WithBestDist] > minDist) {
      Int_t oldPartner = fPartners[iC2WithBestDist];
      fPartners[iC2WithBestDist] = iC1;
      fMinDists[iC2WithBestDist] = minDist;
      //
      fAssociatedLay1[iC1] = 1;      // mark as assigned
      //
      if (oldPartner != -1) {
	// redo partner search for cluster in L0 (oldPartner), putting this one (iC2WithBestDist) on its fBlackList
	fBlackList->AddReference(oldPartner,iC2WithBestDist);
	fAssociatedLay1[oldPartner] = 0;       // mark as free   
      }
    } else {
      // try again to find a cluster without considering iC2WithBestDist 
      fBlackList->AddReference(iC1,iC2WithBestDist);
    } 
    //
  }
  else fAssociatedLay1[iC1] = 2;// cluster has no partner; remove
  //
  return 1;
}

//____________________________________________________________________
Int_t AliITSMultReconstructor::StoreTrackletForL2Cluster(Int_t iC2)
{
  // build tracklet for cluster iC2 of layer 2
  if (fPartners[iC2] == -1) return 0;
  if (fRemoveClustersFromOverlaps) FlagClustersInOverlapRegions (fPartners[iC2],iC2);
  // Printf("saving tracklets");
  if (fOverlapFlagClustersLay[0][fPartners[iC2]] || fOverlapFlagClustersLay[1][iC2]) return 0;
  float* clPar2 = GetClusterLayer2(iC2);
  float* clPar1 = GetClusterLayer1(fPartners[iC2]);
  //
  Float_t* tracklet = fTracklets[fNTracklets] = new Float_t[kTrNPar]; // RS Add also the cluster id's
  //
  tracklet[kTrTheta] = clPar1[kClTh];    // use the theta from the clusters in the first layer
  tracklet[kTrPhi]   = clPar1[kClPh];    // use the phi from the clusters in the first layer
  tracklet[kTrDPhi] = clPar1[kClPh] - clPar2[kClPh];  // store the difference between phi1 and phi2
  //
  // define dphi in the range [0,pi] with proper sign (track charge correlated)
  if (tracklet[kTrDPhi] > TMath::Pi())   tracklet[kTrDPhi] = tracklet[kTrDPhi]-2.*TMath::Pi();
  if (tracklet[kTrDPhi] < -TMath::Pi())  tracklet[kTrDPhi] = tracklet[kTrDPhi]+2.*TMath::Pi();
  //
  tracklet[kTrDTheta] = clPar1[kClTh] - clPar2[kClTh]; // store the theta1-theta2
  //
  if (fHistOn) {
    fhClustersDPhiAcc->Fill(tracklet[kTrDPhi]); 
    fhClustersDThetaAcc->Fill(tracklet[kTrDTheta]);    
    fhDPhiVsDThetaAcc->Fill(tracklet[kTrDTheta],tracklet[kTrDPhi]);
  }
  //
  // find label
  // if equal label in both clusters found this label is assigned
  // if no equal label can be found the first labels of the L1 AND L2 cluster are assigned
  Int_t label1=0,label2=0;
  while (label2 < 3) {
    if ( int(clPar1[kClMC0+label1])!=-2 && int(clPar1[kClMC0+label1])==int(clPar2[kClMC0+label2])) break;
    if (++label1 == 3) { label1 = 0; label2++; }
  }
  if (label2 < 3) {
    AliDebug(AliLog::kDebug, Form("Found label %d == %d for tracklet candidate %d\n", 
				  (Int_t) clPar1[kClMC0+label1], (Int_t) clPar1[kClMC0+label2], fNTracklets));
    tracklet[kTrLab1] = tracklet[kTrLab2] = clPar1[kClMC0+label1];
  } else {
    AliDebug(AliLog::kDebug, Form("Did not find label %d %d %d %d %d %d for tracklet candidate %d\n", 
				  (Int_t) clPar1[kClMC0], (Int_t) clPar1[kClMC1], (Int_t) clPar1[kClMC2], 
				  (Int_t) clPar2[kClMC0], (Int_t) clPar2[kClMC1], (Int_t) clPar2[kClMC2], fNTracklets));
    tracklet[kTrLab1] = clPar1[kClMC0];
    tracklet[kTrLab2] = clPar2[kClMC0];
  }
  //
  if (fHistOn) {
    Float_t eta = tracklet[kTrTheta];
    eta= TMath::Tan(eta/2.);
    eta=-TMath::Log(eta);
    fhetaTracklets->Fill(eta);
    fhphiTracklets->Fill(tracklet[kTrPhi]);
  }
  //
  tracklet[kClID1] = fPartners[iC2];
  tracklet[kClID2] = iC2;
  //
  // Printf("Adding tracklet candidate");
  AliDebug(1,Form(" Adding tracklet candidate %d ", fNTracklets));
  AliDebug(1,Form(" Cl. %d of Layer 1 and %d of Layer 2", fPartners[iC2], iC2));
  fNTracklets++;
  fAssociatedLay1[fPartners[iC2]] = 1;
  // 
  return 1;
}

//____________________________________________________________________
void AliITSMultReconstructor::StoreL1Singles()
{
  // Printf("saving single clusters...");
  for (Int_t iC1=0; iC1<fNClustersLay[0]; iC1++) {
    float* clPar1 = GetClusterLayer1(iC1);
    if (fAssociatedLay1[iC1]==2||fAssociatedLay1[iC1]==0) { 
      fSClusters[fNSingleCluster] = new Float_t[kClNPar];
      fSClusters[fNSingleCluster][kSCTh] = clPar1[kClTh];
      fSClusters[fNSingleCluster][kSCPh] = clPar1[kClPh];
      fSClusters[fNSingleCluster][kSCLab] = clPar1[kClMC0]; 
      fSClusters[fNSingleCluster][kSCID] = iC1;
      AliDebug(1,Form(" Adding a single cluster %d (cluster %d  of layer 1)",
		      fNSingleCluster, iC1));
      fNSingleCluster++;
    }
  }
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::ProcessESDTracks()
{
  // Flag the clusters used by ESD tracks
  // Flag primary tracks to be used for multiplicity counting 
  //
  if (!fESDEvent || !fBuildRefs) return;
  AliESDVertex* vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
  if (!vtx || vtx->GetNContributors()<1) vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexSPD();
  if (!vtx || vtx->GetNContributors()<1) {
    AliDebug(1,"No primary vertex: cannot flag primary tracks");
    return;
  }
  Int_t ntracks = fESDEvent->GetNumberOfTracks();
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack* track = fESDEvent->GetTrack(itr);
    if (!track->IsOn(AliESDtrack::kITSin)) continue; // use only tracks propagated in ITS to vtx
    FlagTrackClusters(itr);
    FlagIfSecondary(track,vtx);
  }
  FlagV0s(vtx);
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagTrackClusters(Int_t id)
{
  // RS: flag the SPD clusters of the track if it is useful for the multiplicity estimation
  //
  const AliESDtrack* track = fESDEvent->GetTrack(id);
  Int_t idx[12];
  if ( track->GetITSclusters(idx)<3 ) return; // at least 3 clusters must be used in the fit
  Int_t itsType = track->IsOn(AliESDtrack::kITSpureSA) ? 1:0;
  
  for (int i=6/*AliESDfriendTrack::kMaxITScluster*/;i--;) { // ignore extras: note: i>=6 is for extra clusters
    if (idx[i]<0) continue;
    int layID= (idx[i] & 0xf0000000) >> 28; 
    if (layID>1) continue; // SPD only
    int clID = (idx[i] & 0x0fffffff);
    fUsedClusLay[layID][itsType]->AddReference(clID,id);
    fStoreRefs[layID][itsType] = kTRUE;
  }
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagIfSecondary(AliESDtrack* track, const AliVertex* vtx)
{
  // RS: check if the track is primary and set the flag
  double cut = (track->HasPointOnITSLayer(0)||track->HasPointOnITSLayer(1)) ? fCutPxDrSPDin:fCutPxDrSPDout;
  float xz[2];
  track->GetDZ(vtx->GetX(),vtx->GetY(),vtx->GetZ(), fESDEvent->GetMagneticField(), xz);
  if (TMath::Abs(xz[0]*track->P())>cut || TMath::Abs(xz[1]*track->P())>fCutPxDz ||
      TMath::Abs(xz[0])>fCutDCArz   || TMath::Abs(xz[1])>fCutDCArz) 
    track->SetStatus(AliESDtrack::kMultSec);
  else track->ResetStatus(AliESDtrack::kMultSec);
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagV0s(const AliESDVertex *vtx)
{
  // flag tracks belonging to v0s
  //
  const double kK0Mass = 0.4976;
  //
  AliV0 pvertex;
  AliKFVertex vertexKF;
  AliKFParticle epKF0,epKF1,pipmKF0,piKF0,piKF1,gammaKF,k0KF;
  Double_t mass,massErr,chi2c;
  enum {kKFIni=BIT(14)};
  //
  double recVtx[3];
  float recVtxF[3];
  vtx->GetXYZ(recVtx);
  for (int i=3;i--;) recVtxF[i] = recVtx[i];
  //
  int ntracks = fESDEvent->GetNumberOfTracks();
  if (ntracks<2) return;
  //
  vertexKF.X() = recVtx[0];
  vertexKF.Y() = recVtx[1];
  vertexKF.Z() = recVtx[2];
  vertexKF.Covariance(0,0) = vtx->GetXRes()*vtx->GetXRes();
  vertexKF.Covariance(1,2) = vtx->GetYRes()*vtx->GetYRes();
  vertexKF.Covariance(2,2) = vtx->GetZRes()*vtx->GetZRes();
  //
  AliESDtrack *trc0,*trc1;
  for (int it0=0;it0<ntracks;it0++) {
    trc0 = fESDEvent->GetTrack(it0);
    if (trc0->IsOn(AliESDtrack::kMultInV0)) continue;
    if (!trc0->IsOn(AliESDtrack::kITSin)) continue;
    Bool_t isSAP = trc0->IsPureITSStandalone();
    Int_t  q0 = trc0->Charge();
    Bool_t testGamma = CanBeElectron(trc0);
    epKF0.ResetBit(kKFIni);
    piKF0.ResetBit(kKFIni);
    double bestChi2=1e16;
    int bestID = -1;
    //    
    for (int it1=it0+1;it1<ntracks;it1++) {
      trc1 = fESDEvent->GetTrack(it1);
      if (trc1->IsOn(AliESDtrack::kMultInV0)) continue;
      if (!trc1->IsOn(AliESDtrack::kITSin)) continue;
      if (trc1->IsPureITSStandalone() != isSAP) continue; // pair separately ITS_SA_Pure tracks and TPC/ITS+ITS_SA
      if ( (q0+trc1->Charge())!=0 ) continue;             // don't pair like signs
      //
      pvertex.SetParamN(q0<0 ? *trc0:*trc1);
      pvertex.SetParamP(q0>0 ? *trc0:*trc1);
      pvertex.Update(recVtxF);
      if (pvertex.P()<fCutMinP) continue;
      if (pvertex.GetV0CosineOfPointingAngle()<fCutMinPointAngle) continue;
      if (pvertex.GetDcaV0Daughters()>fCutMaxDCADauther) continue;
      double d = pvertex.GetD(recVtx[0],recVtx[1],recVtx[2]);
      if (d>fCutMaxDCA) continue;
      double dx=recVtx[0]-pvertex.Xv(), dy=recVtx[1]-pvertex.Yv();
      double rv = TMath::Sqrt(dx*dx+dy*dy);
      //
      // check gamma conversion hypothesis ----------------------------------------------------------->>>
      Bool_t gammaOK = kFALSE;
      while (testGamma && CanBeElectron(trc1)) {
	if (rv<fCutMinRGamma) break;
	if (!epKF0.TestBit(kKFIni)) {
	  new(&epKF0) AliKFParticle(*trc0,q0>0 ? kPositron:kElectron);
	  epKF0.SetBit(kKFIni);
	}
	new(&epKF1) AliKFParticle(*trc1,q0<0 ? kPositron:kElectron);
	gammaKF.Initialize();
	gammaKF += epKF0;
	gammaKF += epKF1;      
	gammaKF.SetProductionVertex(vertexKF);
	gammaKF.GetMass(mass,massErr);
	if (mass>fCutMassGamma || (massErr>0&&(mass>massErr*fCutMassGammaNSigma))) break;
	if (gammaKF.GetS()<fCutGammaSFromDecay) break;
	gammaKF.SetMassConstraint(0.,0.001);
	chi2c = (gammaKF.GetNDF()!=0) ? gammaKF.GetChi2()/gammaKF.GetNDF() : 1000;
	if (chi2c>fCutChi2cGamma) break;
	gammaOK = kTRUE;
	if (chi2c>bestChi2) break;
	bestChi2 = chi2c;
	bestID = it1;
	break;
      }
      if (gammaOK) continue;
      // check gamma conversion hypothesis -----------------------------------------------------------<<<
      // check K0 conversion hypothesis    ----------------------------------------------------------->>>
      while (1) {
	if (rv<fCutMinRK0) break;
	if (!piKF0.TestBit(kKFIni)) {
	  new(&piKF0) AliKFParticle(*trc0,q0>0 ? kPiPlus:kPiMinus);
	  piKF0.SetBit(kKFIni);
	}
	new(&piKF1) AliKFParticle(*trc1,q0<0 ? kPiPlus:kPiMinus);
	k0KF.Initialize();
	k0KF += piKF0;
	k0KF += piKF1;      
	k0KF.SetProductionVertex(vertexKF);
	k0KF.GetMass(mass,massErr);
	mass -= kK0Mass;
	if (TMath::Abs(mass)>fCutMassK0 || (massErr>0&&(abs(mass)>massErr*fCutMassK0NSigma))) break;
	if (k0KF.GetS()<fCutK0SFromDecay) break;
	k0KF.SetMassConstraint(kK0Mass,0.001);
	chi2c = (k0KF.GetNDF()!=0) ? k0KF.GetChi2()/k0KF.GetNDF() : 1000;
	if (chi2c>fCutChi2cK0) break;
	if (chi2c>bestChi2) break;
	bestChi2 = chi2c;
	bestID = it1;
	break;
      }
      // check K0 conversion hypothesis    -----------------------------------------------------------<<<
    }
    //
    if (bestID>=0) {
      trc0->SetStatus(AliESDtrack::kMultInV0);
      fESDEvent->GetTrack(bestID)->SetStatus(AliESDtrack::kMultInV0);
    }
  }
  //
}

//____________________________________________________________________
Bool_t AliITSMultReconstructor::CanBeElectron(const AliESDtrack* trc) const
{
  // check if the track can be electron
  Double_t pid[AliPID::kSPECIES];
  if (!trc->IsOn(AliESDtrack::kESDpid)) return kTRUE;
  trc->GetESDpid(pid);
  return (trc->IsOn(AliESDtrack::kTPCpid)) ? 
    pid[AliPID::kElectron]>fCutMinElectronProbTPC : 
    pid[AliPID::kElectron]>fCutMinElectronProbESD;
  //
}
