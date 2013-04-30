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

//-----------------------------------------------------------------
//   Implementation of the alignment steering class
//   It provides an access to the track space points
//   written along the esd tracks. The class enables
//   the user to plug any track fitter (deriving from
//   AliTrackFitter class) and minimization fo the
//   track residual sums (deriving from the AliTrackResiduals).
//-----------------------------------------------------------------

#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TSystem.h>

#include "AliAlignmentTracks.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackResidualsChi2.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDfriend.h"

ClassImp(AliAlignmentTracks)

//______________________________________________________________________________
AliAlignmentTracks::AliAlignmentTracks():
  fESDChain(0),
  fPointsFilename("AliTrackPoints.root"),
  fPointsFile(0),
  fPointsTree(0),
  fLastIndex(0),
  fArrayIndex(0),
  fIsIndexBuilt(kFALSE),
  fAlignObjs(0),
  fMisalignObjs(0),
  fTrackFitter(0),
  fMinimizer(0),
  fDoUpdate(kTRUE),
  fCovIsUsed(kFALSE)
{
  // Default constructor
  InitIndex();
  InitAlignObjs();
}

//______________________________________________________________________________
AliAlignmentTracks::AliAlignmentTracks(TChain *esdchain):
  fESDChain(esdchain),
  fPointsFilename("AliTrackPoints.root"),
  fPointsFile(0),
  fPointsTree(0),
  fLastIndex(0),
  fArrayIndex(0),
  fIsIndexBuilt(kFALSE),
  fAlignObjs(0),
  fMisalignObjs(0),
  fTrackFitter(0),
  fMinimizer(0),
  fDoUpdate(kTRUE),
  fCovIsUsed(kFALSE)
{
  // Constructor in the case
  // the user provides an already
  // built TChain with ESD trees
  InitIndex();
  InitAlignObjs();
}


//______________________________________________________________________________
AliAlignmentTracks::AliAlignmentTracks(const char *esdfilename, const char *esdtreename):
  fESDChain(new TChain(esdtreename)),
  fPointsFilename("AliTrackPoints.root"),
  fPointsFile(0),
  fPointsTree(0),
  fLastIndex(0),
  fArrayIndex(0),
  fIsIndexBuilt(kFALSE),
  fAlignObjs(0),
  fMisalignObjs(0),
  fTrackFitter(0),
  fMinimizer(0),
  fDoUpdate(kTRUE),
  fCovIsUsed(kFALSE)
{
  // Constructor in the case
  // the user provides a single ESD file
  // or a directory containing ESD files
  fESDChain->Add(esdfilename);

  InitIndex();
  InitAlignObjs();
}


//______________________________________________________________________________
AliAlignmentTracks::~AliAlignmentTracks()
{
  // Destructor
  if (fESDChain) delete fESDChain;

  DeleteIndex();
  DeleteAlignObjs();

  delete fTrackFitter;
  delete fMinimizer;

  if (fPointsFile) fPointsFile->Close();
}

//______________________________________________________________________________
void AliAlignmentTracks::AddESD(TChain *esdchain)
{
  // Add a chain with ESD files
  if (fESDChain)
    fESDChain->Add(esdchain);
  else
    fESDChain = esdchain;
}

//______________________________________________________________________________
void AliAlignmentTracks::AddESD(const char *esdfilename, const char *esdtreename)
{
  // Add a single file or
  // a directory to the chain
  // with the ESD files
  if (fESDChain)
    fESDChain->AddFile(esdfilename,TChain::kBigNumber,esdtreename);
  else {
    fESDChain = new TChain(esdtreename);
    fESDChain->Add(esdfilename);
  }
}


//________________________________________________________________________
void AliAlignmentTracks::ProcessESD(Bool_t onlyITS,
				    Int_t minITSpts,
				    Bool_t cuts,
				    Float_t minAngleWrtITSModulePlanes,
				    Float_t minMom,Float_t maxMom,
				    Float_t minAbsSinPhi,Float_t maxAbsSinPhi,
				    Float_t minSinTheta,Float_t maxSinTheta)
{
  // Analyzes and filters ESD tracks
  // Stores the selected track space points
  // into the output file

  if (!fESDChain) return;

  AliESDEvent *esd = new AliESDEvent();
  esd->ReadFromTree(fESDChain);
  AliESDfriend *esdf = 0; 
  fESDChain->SetBranchStatus("ESDfriend*",1);
  fESDChain->SetBranchAddress("ESDfriend.",&esdf);

  // Open the output file
  if (fPointsFilename.IsNull()) {
    AliWarning("Incorrect output filename!");
    return;
  }

  TFile *pointsFile = TFile::Open(fPointsFilename,"RECREATE");
  if (!pointsFile || !pointsFile->IsOpen()) {
    AliWarning(Form("Can't open %s !",fPointsFilename.Data()));
    return;
  }

  TTree *pointsTree = new TTree("spTree", "Tree with track space point arrays");
  const AliTrackPointArray *array = 0;
  AliTrackPointArray *array2 = 0;
  if(onlyITS) {   // only ITS AliTrackPoints 
    pointsTree->Branch("SP","AliTrackPointArray", &array2);
  } else {
    pointsTree->Branch("SP","AliTrackPointArray", &array);
  } 

  Int_t ievent = 0;
  while (fESDChain->GetEntry(ievent++)) {
    if (!esd) break;

    esd->SetESDfriend(esdf); //Attach the friend to the ESD

    Int_t ntracks = esd->GetNumberOfTracks();
    for (Int_t itrack=0; itrack < ntracks; itrack++) {
      AliESDtrack * track = esd->GetTrack(itrack);
      if (!track) continue;

      if(track->GetNcls(0) < minITSpts) continue;
      if(cuts) {
	if(track->GetP()<minMom || track->GetP()>maxMom) continue;
	Float_t abssinphi = TMath::Abs(TMath::Sin(track->GetAlpha()+TMath::ASin(track->GetSnp())));
	if(abssinphi<minAbsSinPhi || abssinphi>maxAbsSinPhi) continue;
	Float_t sintheta = TMath::Sin(0.5*TMath::Pi()-TMath::ATan(track->GetTgl()));
	if(sintheta<minSinTheta || sintheta>maxSinTheta) continue;
      } 

      AliTrackPoint point;
      array = track->GetTrackPointArray();

      if(onlyITS) {
	Bool_t layerOK[6]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
	Int_t ipt,volId,modId,layerId;
	Int_t jpt=0;
	for(ipt=0; ipt<array->GetNPoints(); ipt++) {
	  array->GetPoint(point,ipt);
	  volId = point.GetVolumeID();
	  layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	  if(layerId>6) continue;
	  // check minAngleWrtITSModulePlanes
	  if(cuts) {
	    Double_t p[3]; track->GetDirection(p);
	    TVector3 pvec(p[0],p[1],p[2]);
	    Double_t rot[9]; AliGeomManager::GetOrigRotation(volId,rot);
	    TVector3 normvec(rot[1],rot[4],rot[7]);
	    Double_t angle = pvec.Angle(normvec);
	    if(angle>0.5*TMath::Pi()) angle = TMath::Pi()-angle;
	    angle = 0.5*TMath::Pi()-angle;
	    if(angle<minAngleWrtITSModulePlanes) {
	      layerOK[layerId-1]=kFALSE;
	      continue;
	    }
	  }
	  jpt++;
	}
	if(jpt < minITSpts) continue;      
	array2 = new AliTrackPointArray(jpt);
	jpt=0;
	for(ipt=0; ipt<array->GetNPoints(); ipt++) {
	  array->GetPoint(point,ipt);
	  volId = point.GetVolumeID();
	  layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	  if(layerId>6 || !layerOK[layerId-1]) continue;
	  array2->AddPoint(jpt,&point);
	  jpt++;
	}
      } // end if(onlyITS)
 
      pointsTree->Fill();
    }
  }

  if (!pointsTree->Write()) {
    AliWarning("Can't write the tree with track point arrays!");
    return;
  }

  pointsFile->Close();

  return;
}

//_____________________________________________________________________________
void AliAlignmentTracks::ProcessESDCosmics(Bool_t onlyITS,
				   Int_t minITSpts,Float_t maxMatchingAngle,
				   Bool_t cuts,
				   Float_t minAngleWrtITSModulePlanes,
				   Float_t minMom,Float_t maxMom,
				   Float_t minAbsSinPhi,Float_t maxAbsSinPhi,
			     	   Float_t minSinTheta,Float_t maxSinTheta)
{
  // Analyzes and filters ESD tracks
  // Merges inward and outward tracks in one single track
  // Stores the selected track space points
  // into the output file

  if (!fESDChain) return;

  AliESDEvent *esd = new AliESDEvent();
  esd->ReadFromTree(fESDChain);
  AliESDfriend *esdf = 0; 
  fESDChain->SetBranchStatus("ESDfriend*",1);
  fESDChain->SetBranchAddress("ESDfriend.",&esdf);

  // Open the output file
  if (fPointsFilename.IsNull()) {
    AliWarning("Incorrect output filename!");
    return;
  }

  TFile *pointsFile = TFile::Open(fPointsFilename,"RECREATE");
  if (!pointsFile || !pointsFile->IsOpen()) {
    AliWarning(Form("Can't open %s !",fPointsFilename.Data()));
    return;
  }

  TTree *pointsTree = new TTree("spTree", "Tree with track space point arrays");
  const AliTrackPointArray *array = 0;
  AliTrackPointArray *array2 = 0;
  pointsTree->Branch("SP","AliTrackPointArray", &array2);

  Int_t ievent = 0;
  while (fESDChain->GetEntry(ievent++)) {
    if (!esd) break;

    esd->SetESDfriend(esdf); //Attach the friend to the ESD

    Int_t ntracks = esd->GetNumberOfTracks();
    if(ntracks<2) continue;
    Int_t *goodtracksArray = new Int_t[ntracks];
    Float_t *phiArray = new Float_t[ntracks];
    Float_t *thetaArray = new Float_t[ntracks];
    Int_t ngt=0;
    for (Int_t itrack=0; itrack < ntracks; itrack++) {
      AliESDtrack * track = esd->GetTrack(itrack);
      if (!track) continue;

      if(track->GetNcls(0) < minITSpts) continue;
      Float_t phi = track->GetAlpha()+TMath::ASin(track->GetSnp());
      Float_t theta = 0.5*TMath::Pi()-TMath::ATan(track->GetTgl());
      if(cuts) {
	if(track->GetP()<minMom || track->GetP()>maxMom) continue;
	Float_t abssinphi = TMath::Abs(TMath::Sin(phi));
	if(abssinphi<minAbsSinPhi || abssinphi>maxAbsSinPhi) continue;
	Float_t sintheta = TMath::Sin(theta);
	if(sintheta<minSinTheta || sintheta>maxSinTheta) continue;
      } 
      goodtracksArray[ngt]=itrack;
      phiArray[ngt]=phi;
      thetaArray[ngt]=theta;
      ngt++;
    }

    if(ngt<2) {
      delete [] goodtracksArray; goodtracksArray=0;
      delete [] phiArray; phiArray=0;
      delete [] thetaArray; thetaArray=0;
      continue;
    }

    // check matching of the two tracks from the muon
    Float_t min = 10000000.;
    Int_t good1 = -1, good2 = -1;
    for(Int_t itr1=0; itr1<ngt-1; itr1++) {
      for(Int_t itr2=itr1+1; itr2<ngt; itr2++) {
	Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArray[itr1]-thetaArray[itr2]);
	if(deltatheta>maxMatchingAngle) continue;
	Float_t deltaphi = TMath::Abs(TMath::Abs(phiArray[itr1]-phiArray[itr2])-TMath::Pi());
	if(deltaphi>maxMatchingAngle) continue;
	//printf("%f  %f     %f  %f\n",deltaphi,deltatheta,thetaArray[itr1],thetaArray[itr2]);
	if(deltatheta+deltaphi<min) {
	  min=deltatheta+deltaphi;
	  good1 = goodtracksArray[itr1];
	  good2 = goodtracksArray[itr2];
	}
      }
    }

    delete [] goodtracksArray; goodtracksArray=0;
    delete [] phiArray; phiArray=0;
    delete [] thetaArray; thetaArray=0;

    if(good1<0) continue;

    AliESDtrack * track1 = esd->GetTrack(good1);
    AliESDtrack * track2 = esd->GetTrack(good2);

    AliTrackPoint point;
    Int_t ipt,volId,modId,layerId;
    Int_t jpt=0;
    Bool_t layerOK[6][2]; 
    for(Int_t l1=0;l1<6;l1++) for(Int_t l2=0;l2<2;l2++) layerOK[l1][l2]=kTRUE; 
    array = track1->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      if(onlyITS) {
	volId = point.GetVolumeID();
	layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	if(layerId>6) continue;
	// check minAngleWrtITSModulePlanes
	if(cuts) {
	  Double_t p[3]; track1->GetDirection(p);
	  TVector3 pvec(p[0],p[1],p[2]);
	  Double_t rot[9]; AliGeomManager::GetOrigRotation(volId,rot);
	  TVector3 normvec(rot[1],rot[4],rot[7]);
	  Double_t angle = pvec.Angle(normvec);
	  if(angle>0.5*TMath::Pi()) angle = TMath::Pi()-angle;
	  angle = 0.5*TMath::Pi()-angle;
	  if(angle<minAngleWrtITSModulePlanes) {
	    layerOK[layerId-1][0]=kFALSE;
	    continue;
	  }
	}
      }
      jpt++;
    }
    array = track2->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      if(onlyITS) {
	volId = point.GetVolumeID();
	layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	if(layerId>6) continue;
	// check minAngleWrtITSModulePlanes
	if(cuts) {
	  Double_t p[3]; track2->GetDirection(p);
	  TVector3 pvec(p[0],p[1],p[2]);
	  Double_t rot[9]; AliGeomManager::GetOrigRotation(volId,rot);
	  TVector3 normvec(rot[1],rot[4],rot[7]);
	  Double_t angle = pvec.Angle(normvec);
	  if(angle>0.5*TMath::Pi()) angle = TMath::Pi()-angle;
	  angle = 0.5*TMath::Pi()-angle;
	  if(angle<minAngleWrtITSModulePlanes) {
	    layerOK[layerId-1][0]=kFALSE;
	    continue;
	  }
	}
      }
      jpt++;
    }

    if(jpt < 2*minITSpts) continue;
    array2 = new AliTrackPointArray(jpt);
    jpt=0;
    array = track1->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      if(onlyITS) {
	volId = point.GetVolumeID();
	layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	if(layerId>6 || !layerOK[layerId-1][0]) continue;
      }
      array2->AddPoint(jpt,&point);
      jpt++;
    }
    array = track2->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      if(onlyITS) {
	volId = point.GetVolumeID();
	layerId = AliGeomManager::VolUIDToLayer(volId,modId);
	if(layerId>6 || !layerOK[layerId-1][1]) continue;
      }
      array2->AddPoint(jpt,&point);
      jpt++;
    }

    pointsTree->Fill();
  }

  if (!pointsTree->Write()) {
    AliWarning("Can't write the tree with track point arrays!");
    return;
  }

  pointsFile->Close();
  return;
}

//______________________________________________________________________________
void AliAlignmentTracks::ProcessESD(TSelector *selector)
{
  AliWarning(Form("ESD processing based on selector is not yet implemented (%p) !",selector));
}

//______________________________________________________________________________
void AliAlignmentTracks::BuildIndex()
{
  // Build index of points tree entries
  // Used for access based on the volume IDs
  if (fIsIndexBuilt) return;

  fIsIndexBuilt = kTRUE;

  // Dummy object is created in order
  // to initialize the volume paths
  AliAlignObjParams alobj;

  fPointsFile = TFile::Open(fPointsFilename);
  if (!fPointsFile || !fPointsFile->IsOpen()) {
    AliWarning(Form("Can't open %s !",fPointsFilename.Data()));
    return;
  }
  
  //  AliTrackPointArray* array = new AliTrackPointArray;
  AliTrackPointArray* array = 0;
  fPointsTree = (TTree*) fPointsFile->Get("spTree");
  if (!fPointsTree) {
    AliWarning("No pointsTree found!");
    return;
  }
  fPointsTree->SetBranchAddress("SP", &array);

  Int_t nArrays = (Int_t)fPointsTree->GetEntries();
  for (Int_t iArray = 0; iArray < nArrays; iArray++)
    {
      fPointsTree->GetEvent(iArray);
      if (!array) continue;
      for (Int_t ipoint = 0; ipoint < array->GetNPoints(); ipoint++) {
	UShort_t volId = array->GetVolumeID()[ipoint];
	// check if the volId is valid
	if (!AliGeomManager::SymName(volId)) {
	  AliError(Form("The volume id %d has no default volume name !",
			volId));
	  continue;
	}
	Int_t modId;
	Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId)
	              - AliGeomManager::kFirstLayer;
	if (!fArrayIndex[layerId][modId]) {
	  //first entry for this volume
	  fArrayIndex[layerId][modId] = new TArrayI(1000);
	}
	else {
	  Int_t size = fArrayIndex[layerId][modId]->GetSize();
	  // If needed allocate new size
	  if (fLastIndex[layerId][modId] >= size)
	    fArrayIndex[layerId][modId]->Set(size + 1000);
	}

	// Check if the index is already filled
	Bool_t fillIndex = kTRUE;
	if (fLastIndex[layerId][modId] != 0) {
	  if ((*fArrayIndex[layerId][modId])[fLastIndex[layerId][modId]-1] == iArray)
	    fillIndex = kFALSE;
	}
	// Fill the index array and store last filled index
	if (fillIndex) {
	  (*fArrayIndex[layerId][modId])[fLastIndex[layerId][modId]] = iArray;
	  fLastIndex[layerId][modId]++;
	}
      }
    }

}

//______________________________________________________________________________
void AliAlignmentTracks::InitIndex()
{
  // Initialize the index arrays
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fLastIndex = new Int_t*[nLayers];
  fArrayIndex = new TArrayI**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fLastIndex[iLayer] = new Int_t[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    fArrayIndex[iLayer] = new TArrayI*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      fLastIndex[iLayer][iModule] = 0;
      fArrayIndex[iLayer][iModule] = 0;
    }
  }
}

//______________________________________________________________________________
void AliAlignmentTracks::ResetIndex()
{
  // Reset the value of the last filled index
  // Do not realocate memory

  fIsIndexBuilt = kFALSE;
  
  for (Int_t iLayer = 0; iLayer < AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      fLastIndex[iLayer][iModule] = 0;
    }
  }
}

//______________________________________________________________________________
void AliAlignmentTracks::DeleteIndex()
{
  // Delete the index arrays
  // Called by the destructor
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      if (fArrayIndex[iLayer][iModule]) {
	delete fArrayIndex[iLayer][iModule];
	fArrayIndex[iLayer][iModule] = 0;
      }
    }
    delete [] fLastIndex[iLayer];
    delete [] fArrayIndex[iLayer];
  }
  delete [] fLastIndex;
  delete [] fArrayIndex;
}

//______________________________________________________________________________
Bool_t AliAlignmentTracks::ReadAlignObjs(const char *alignObjFileName, const char* arrayName)
{
  // Read alignment object from a file: update the alignobj already present with the one in the file
  // To be replaced by a call to CDB
  
  if(gSystem->AccessPathName(alignObjFileName,kFileExists)){
    printf("Wrong AlignObjs File Name \n");
    return kFALSE;
  } 

  TFile *fRealign=TFile::Open(alignObjFileName);
  if (!fRealign || !fRealign->IsOpen()) {
    AliError(Form("Could not open Align Obj File file %s !",alignObjFileName));
    return kFALSE;
  }  
  printf("Getting TClonesArray \n");
  TClonesArray *clnarray=(TClonesArray*)fRealign->Get(arrayName);
  Int_t size=clnarray->GetSize();
  UShort_t volid;

  for(Int_t ivol=0;ivol<size;ivol++){
    AliAlignObjParams *a=(AliAlignObjParams*)clnarray->At(ivol);
    volid=a->GetVolUID();
    Int_t iModule;
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
    if(iLayer<AliGeomManager::kFirstLayer||iLayer>AliGeomManager::kSSD2)continue;
    printf("Updating volume: %d ,layer: %d module: %d \n",volid,iLayer,iModule);
    *fAlignObjs[iLayer-AliGeomManager::kFirstLayer][iModule] *= *a;
  }
 
  delete clnarray;
  fRealign->Close();
  return kTRUE;
}

//______________________________________________________________________________
void AliAlignmentTracks::InitAlignObjs()
{
  // Initialize the alignment objects array
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fAlignObjs = new AliAlignObj**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fAlignObjs[iLayer] = new AliAlignObj*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer+ AliGeomManager::kFirstLayer,iModule);
      fAlignObjs[iLayer][iModule] = new AliAlignObjParams(AliGeomManager::SymName(volid),volid,0,0,0,0,0,0,kTRUE);
    }
  }
}

//______________________________________________________________________________
void AliAlignmentTracks::ResetAlignObjs()
{
  // Reset the alignment objects array
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++)
      fAlignObjs[iLayer][iModule]->SetPars(0,0,0,0,0,0);
  }
}

//______________________________________________________________________________
void AliAlignmentTracks::DeleteAlignObjs()
{
  // Delete the alignment objects array
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++)
      if (fAlignObjs[iLayer][iModule])
	delete fAlignObjs[iLayer][iModule];
    delete [] fAlignObjs[iLayer];
  }
  delete [] fAlignObjs;
  fAlignObjs = 0;
}

Bool_t AliAlignmentTracks::AlignDetector(AliGeomManager::ELayerID firstLayer,
					 AliGeomManager::ELayerID lastLayer,
					 AliGeomManager::ELayerID layerRangeMin,
					 AliGeomManager::ELayerID layerRangeMax,
					 Int_t iterations)
{
  // Align detector volumes within
  // a given layer range
  // (could be whole detector).
  // Tracks are fitted only within
  // the range defined by the user.
  Int_t nModules = 0;
  for (Int_t iLayer = firstLayer; iLayer <= lastLayer; iLayer++)
    nModules += AliGeomManager::LayerSize(iLayer);
  TArrayI volIds(nModules);

  Int_t modnum = 0;
  for (Int_t iLayer = firstLayer; iLayer <= lastLayer; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {
      UShort_t volId = AliGeomManager::LayerToVolUID(iLayer,iModule);
      volIds.AddAt(volId,modnum);
      modnum++;
    }
  }

  Bool_t result = kFALSE;
  while (iterations > 0) {
    if (!(result = AlignVolumes(&volIds,0x0,layerRangeMin,layerRangeMax))) break;
    iterations--;
  }
  return result;
}

//______________________________________________________________________________
Bool_t AliAlignmentTracks::AlignLayer(AliGeomManager::ELayerID layer,
				      AliGeomManager::ELayerID layerRangeMin,
				      AliGeomManager::ELayerID layerRangeMax,
				      Int_t iterations)
{
  // Align detector volumes within
  // a given layer.
  // Tracks are fitted only within
  // the range defined by the user.
  Int_t nModules = AliGeomManager::LayerSize(layer);
  TArrayI volIds(nModules);
  for (Int_t iModule = 0; iModule < nModules; iModule++) {
    UShort_t volId = AliGeomManager::LayerToVolUID(layer,iModule);
    volIds.AddAt(volId,iModule);
  }

  Bool_t result = kFALSE;
  while (iterations > 0) {
    if (!(result = AlignVolumes(&volIds,0x0,layerRangeMin,layerRangeMax))) break;
    iterations--;
  }
  return result;
}

//______________________________________________________________________________
Bool_t AliAlignmentTracks::AlignVolume(UShort_t volId, UShort_t volIdFit,
				     Int_t iterations)
{
  // Align single detector volume to
  // another volume.
  // Tracks are fitted only within
  // the second volume.
  TArrayI volIds(1);
  volIds.AddAt(volId,0);
  TArrayI volIdsFit(1);
  volIdsFit.AddAt(volIdFit,0);

  Bool_t result = kFALSE;
  while (iterations > 0) {
    if (!(result = AlignVolumes(&volIds,&volIdsFit))) break;
    iterations--;
  }
  return result;
}

//______________________________________________________________________________
Bool_t AliAlignmentTracks::AlignVolumes(const TArrayI *volids, const TArrayI *volidsfit,
				     AliGeomManager::ELayerID layerRangeMin,
				     AliGeomManager::ELayerID layerRangeMax,
				     Int_t iterations)
{
  // Align a set of detector volumes.
  // Tracks are fitted only within
  // the range defined by the user
  // (by layerRangeMin and layerRangeMax)
  // or within the set of volidsfit
  // Repeat the procedure 'iterations' times

  Int_t nVolIds = volids->GetSize();
  if (nVolIds == 0) {
    AliError("Volume IDs array is empty!");
    return kFALSE;
  }
 
  // Load only the tracks with at least one
  // space point in the set of volume (volids)
  BuildIndex();
  AliTrackPointArray **points;
  Int_t pointsdim;
  // Start the iterations
  Bool_t result = kFALSE;
  while (iterations > 0) {
    Int_t nArrays = LoadPoints(volids, points,pointsdim);
    if (nArrays == 0) {
      UnloadPoints(pointsdim, points);
      return kFALSE;
    }

    AliTrackResiduals *minimizer = CreateMinimizer();
    minimizer->SetNTracks(nArrays);
    minimizer->InitAlignObj();
    AliTrackFitter *fitter = CreateFitter();
    for (Int_t iArray = 0; iArray < nArrays; iArray++) {
      if (!points[iArray]) continue;
      fitter->SetTrackPointArray(points[iArray], kFALSE);
      if (fitter->Fit(volids,volidsfit,layerRangeMin,layerRangeMax) == kFALSE) continue;
      AliTrackPointArray *pVolId,*pTrack;
      fitter->GetTrackResiduals(pVolId,pTrack);
      minimizer->AddTrackPointArrays(pVolId,pTrack);
    }
    if (!(result = minimizer->Minimize())) {
      UnloadPoints(pointsdim, points);
      break;
    }

    // Update the alignment object(s)
    if (fDoUpdate) for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
      UShort_t volid = (*volids)[iVolId];
      Int_t iModule;
      AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
      AliAlignObj *alignObj = fAlignObjs[iLayer-AliGeomManager::kFirstLayer][iModule];      
      *alignObj *= *minimizer->GetAlignObj();
      if(iterations==1)alignObj->Print("");
    }

    UnloadPoints(pointsdim, points);
    
    iterations--;
  }
  return result;
}
  
//______________________________________________________________________________
Int_t AliAlignmentTracks::LoadPoints(const TArrayI *volids, AliTrackPointArray** &points,Int_t &pointsdim)
{
  // Load track point arrays with at least
  // one space point in a given set of detector
  // volumes (array volids).
  // Use the already created tree index for
  // fast access.

  if (!fPointsTree) {
    AliError("Tree with the space point arrays not initialized!");
    points = 0;
    return 0;
  }

  Int_t nVolIds = volids->GetSize();
  if (nVolIds == 0) {
    AliError("Volume IDs array is empty!");
    points = 0;
    return 0;
  }

  Int_t nArrays = 0;
  for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
    UShort_t volid = (*volids)[iVolId];
    Int_t iModule;
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);

    // In case of empty index
    if (fLastIndex[iLayer-AliGeomManager::kFirstLayer][iModule] == 0) {
      AliWarning(Form("There are no space-points belonging to the volume which is to be aligned (Volume ID =%d)!",volid));
      continue;
    }
    nArrays += fLastIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
  }

  if (nArrays == 0) {
    AliError("There are no space-points belonging to all of the volumes which are to be aligned!");
    points = 0x0;
    return 0;
  }

  AliTrackPointArray* array = 0;
  fPointsTree->SetBranchAddress("SP", &array);

  // Allocate the pointer to the space-point arrays
  pointsdim=nArrays;
  points = new AliTrackPointArray*[nArrays];
  for (Int_t i = 0; i < nArrays; i++) points[i] = 0x0;

  // Init the array used to flag already loaded tree entries
  Bool_t *indexUsed = new Bool_t[(UInt_t)fPointsTree->GetEntries()];
  for (Int_t i = 0; i < fPointsTree->GetEntries(); i++)
    indexUsed[i] = kFALSE;

  // Start the loop over the volume ids
  Int_t iArray = 0;
  for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
    UShort_t volid = (*volids)[iVolId];
    Int_t iModule;
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);

    Int_t nArraysId = fLastIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
    TArrayI *index = fArrayIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
    AliTrackPoint p;

    for (Int_t iArrayId = 0; iArrayId < nArraysId; iArrayId++) {

      // Get tree entry
      Int_t entry = (*index)[iArrayId];
      if (indexUsed[entry] == kTRUE) {
	nArrays--;
	continue;
      }
      fPointsTree->GetEvent(entry);
      if (!array) {
	AliWarning("Wrong space point array index!");
	continue;
      }
      indexUsed[entry] = kTRUE;

      // Get the space-point array
      Int_t nPoints = array->GetNPoints();
      points[iArray] = new AliTrackPointArray(nPoints);
      for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {
	array->GetPoint(p,iPoint);
	Int_t modnum;
	AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(p.GetVolumeID(),modnum);
	// check if the layer id is valid
	if ((layer < AliGeomManager::kFirstLayer) ||
	    (layer >= AliGeomManager::kLastLayer)) {
	  AliError(Form("Layer index is invalid: %d (%d -> %d) !",
			layer,AliGeomManager::kFirstLayer,AliGeomManager::kLastLayer-1));
	  continue;
	}
	if ((modnum >= AliGeomManager::LayerSize(layer)) ||
	    (modnum < 0)) {
	  AliError(Form("Module number inside layer %d is invalid: %d (0 -> %d)",
			layer,modnum,AliGeomManager::LayerSize(layer)));
	  continue;
	}

	// Misalignment is introduced here
	// Switch it off in case of real
	// alignment job!
	if (fMisalignObjs) {
	  AliAlignObj *misalignObj = fMisalignObjs[layer-AliGeomManager::kFirstLayer][modnum];
	  if (misalignObj)
	    misalignObj->Transform(p);
	}
	// End of misalignment


	AliAlignObj *alignObj = fAlignObjs[layer-AliGeomManager::kFirstLayer][modnum];
	UShort_t volp=p.GetVolumeID();
	Bool_t found=kFALSE;
	if(fCovIsUsed){
	  for (Int_t iVol = 0; iVol < nVolIds; iVol++) {
	    UShort_t vol = (*volids)[iVol];
	    if(volp==vol){
	      alignObj->Transform(p,kFALSE);
	      found=kTRUE;
	      break;
	    }
	  }
	}
	if(!found)alignObj->Transform(p,fCovIsUsed);
	points[iArray]->AddPoint(iPoint,&p);
      }
      iArray++;
    }
  }


  delete [] indexUsed;

  return nArrays;
}

//______________________________________________________________________________
void AliAlignmentTracks::UnloadPoints(Int_t n, AliTrackPointArray **points)
{
  // Unload track point arrays for a given
  // detector volume
  for (Int_t iArray = 0; iArray < n; iArray++)
    delete points[iArray];
  delete [] points;
}

//______________________________________________________________________________
AliTrackFitter *AliAlignmentTracks::CreateFitter()
{
  // Check if the user has already supplied
  // a track fitter object.
  // If not, create a default one.
  if (!fTrackFitter)
    fTrackFitter = new AliTrackFitterRieman;

  return fTrackFitter;
}

//______________________________________________________________________________
AliTrackResiduals *AliAlignmentTracks::CreateMinimizer()
{
  // Check if the user has already supplied
  // a track residuals minimizer object.
  // If not, create a default one.
  if (!fMinimizer)
    fMinimizer = new AliTrackResidualsChi2;

  return fMinimizer;
}

//______________________________________________________________________________
Bool_t AliAlignmentTracks::Misalign(const char *misalignObjFileName, const char* arrayName)
{
  // The method reads from a file a set of AliAlignObj which are
  // then used to apply misalignments directly on the track
  // space-points. The method is supposed to be used only for
  // fast development and debugging of the alignment algorithms.
  // Be careful not to use it in the case of 'real' alignment
  // scenario since it will bias the results.

  // Initialize the misalignment objects array
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fMisalignObjs = new AliAlignObj**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fMisalignObjs[iLayer] = new AliAlignObj*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++)
      fMisalignObjs[iLayer][iModule] = 0x0;
  }

  // Open the misliagnment file and load the array with
  // misalignment objects
  TFile* inFile = TFile::Open(misalignObjFileName,"READ");
  if (!inFile || !inFile->IsOpen()) {
    AliError(Form("Could not open misalignment file %s !",misalignObjFileName));
    return kFALSE;
  }

  TClonesArray* array = ((TClonesArray*) inFile->Get(arrayName));
  if (!array) {
    AliError(Form("Could not find misalignment array %s in the file %s !",arrayName,misalignObjFileName));
    inFile->Close();
    return kFALSE;
  }
  inFile->Close();

  // Store the misalignment objects for further usage  
  Int_t nObjs = array->GetEntriesFast();
  AliGeomManager::ELayerID layerId; // volume layer
  Int_t modId; // volume ID inside the layer
  for(Int_t i=0; i<nObjs; i++)
    {
      AliAlignObj* alObj = (AliAlignObj*)array->UncheckedAt(i);
      alObj->GetVolUID(layerId,modId);
      if(layerId<AliGeomManager::kFirstLayer) {
	AliWarning(Form("Alignment object is ignored: %s",alObj->GetSymName()));
	continue;
      }
      fMisalignObjs[layerId-AliGeomManager::kFirstLayer][modId] = alObj;
    }
  return kTRUE;
}


//________________________________________________
void AliAlignmentTracks::WriteRealignObjArray(TString outfilename,AliGeomManager::ELayerID layerRangeMin,AliGeomManager::ELayerID layerRangeMax){
  
  Int_t last=0;
  TClonesArray *clonesArray=new TClonesArray("AliAlignObjParams",2200);
  TClonesArray &alo=*clonesArray;
  for (Int_t iLayer = layerRangeMin-AliGeomManager::kFirstLayer; iLayer <= (layerRangeMax - AliGeomManager::kFirstLayer);iLayer++) {
  
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
     
      AliAlignObj *alignObj = fAlignObjs[iLayer][iModule]; 
      new(alo[last])AliAlignObjParams(*alignObj);
      last++;
    }
  }
  TFile *file=new TFile(outfilename.Data(),"RECREATE");
  file->cd();
 
  alo.Write("ITSAlignObjs",TObject::kSingleKey);
  file->Close();    
 
     
  delete clonesArray;
  return;
}
