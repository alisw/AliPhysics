/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTask to extract from ESD tracks the AliTrackPointArrays
// with ITS points for selected tracks. This are the input data for alignment
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMap.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TRandom.h>

#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"
#include "AliITSReconstructor.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSgeomTGeo.h"
#include "AliTrackPointArray.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAlignmentDataFilterITS.h"


ClassImp(AliAlignmentDataFilterITS)


//________________________________________________________________________
AliAlignmentDataFilterITS::AliAlignmentDataFilterITS():
AliAnalysisTaskSE(),
fOnlySPDFO(kFALSE),
fDownsamplelowpt(kFALSE),
fGeometryFileName("geometry.root"),
fITSRecoParam(0),
fESD(0),
fListOfHistos(0),
fspTree(0),
fHistNevents(0),
fHistNpoints(0),
fHistPt(0),
fHistLayer0(0),
fHistLayer1(0),
fHistLayer2(0),
fHistLayer3(0),
fHistLayer4(0),
fHistLayer5(0),
fntExtra(0),
fntCosmicMatching(0)
{
  // Constructor
}

//________________________________________________________________________
AliAlignmentDataFilterITS::AliAlignmentDataFilterITS(const char *name):
AliAnalysisTaskSE(name),
fOnlySPDFO(kFALSE),
fDownsamplelowpt(kFALSE),
fGeometryFileName("geometry.root"),
fITSRecoParam(0),
fESD(0),
fListOfHistos(0),
fspTree(0),
fHistNevents(0),
fHistNpoints(0),
fHistPt(0),
fHistLayer0(0),
fHistLayer1(0),
fHistLayer2(0),
fHistLayer3(0),
fHistLayer4(0),
fHistLayer5(0),
fntExtra(0),
fntCosmicMatching(0)
{
  // Constructor

  // Define output slots here

  // Output slot #1 writes into a TTree
  DefineOutput(1,TTree::Class());  //My private output
  // Output slot #2 writes into a TList
  DefineOutput(2,TList::Class());  //My private output
}

//________________________________________________________________________
AliAlignmentDataFilterITS::~AliAlignmentDataFilterITS()
{
  // Destructor
  if (fListOfHistos) {
    delete fListOfHistos;
    fListOfHistos = 0;
  }
  if (fspTree) {
    delete fspTree;
    fspTree = 0;
  }
  if (fHistNevents) {
    delete fHistNevents;
    fHistNevents = 0;
  }
  if (fHistNpoints) {
    delete fHistNpoints;
    fHistNpoints = 0;
  }
  if (fHistPt) {
    delete fHistPt;
    fHistPt = 0;
  }
  if (fHistLayer0) {
    delete fHistLayer0;
    fHistLayer0 = 0;
  }
  if (fHistLayer1) {
    delete fHistLayer1;
    fHistLayer1 = 0;
  }
  if (fHistLayer2) {
    delete fHistLayer2;
    fHistLayer2 = 0;
  }
  if (fHistLayer3) {
    delete fHistLayer3;
    fHistLayer3 = 0;
  }
  if (fHistLayer4) {
    delete fHistLayer4;
    fHistLayer4 = 0;
  }
  if (fHistLayer5) {
    delete fHistLayer5;
    fHistLayer5 = 0;
  }
  if (fntExtra) {
    delete fntExtra;
    fntExtra = 0;
  }
  if (fntCosmicMatching) {
    delete fntCosmicMatching;
    fntCosmicMatching = 0;
  }
}  


//________________________________________________________________________
void AliAlignmentDataFilterITS::UserCreateOutputObjects()
{
  // Create the output container
  //
  
  // load the geometry  
  if(!gGeoManager) {    
    printf("AliAlignmentDataFilterITS::CreateOutputObjects(): loading geometry from %s\n",fGeometryFileName.Data());
    AliGeomManager::LoadGeometry(fGeometryFileName.Data());
    if(!gGeoManager) { 
      printf("AliAlignmentDataFilterITS::CreateOutputObjects(): no geometry loaded \n");
      return;
    }
  }

  // Several histograms are more conveniently managed in a TList
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();

  fHistNevents = new TH1F("fHistNevents", "Number of processed events; N events; bin",5,-0.5,4.5);
  fHistNevents->Sumw2();
  fHistNevents->SetMinimum(0);
  fListOfHistos->Add(fHistNevents);

  fHistNpoints = new TH1F("fHistNpoints", "Number of AliTrackPoints per track; N points; tracks",25,-0.5,24.5);
  fHistNpoints->Sumw2();
  fHistNpoints->SetMinimum(0);
  fListOfHistos->Add(fHistNpoints);

  fHistPt = new TH1F("fHistPt", "p_{t} of tracks; p_{t} [GeV/c]; tracks",100,0,50);
  fHistPt->Sumw2();
  fHistPt->SetMinimum(0);
  fListOfHistos->Add(fHistPt);


  Float_t zmax=14.;
  Int_t nbinsphi=20,nbinsz=4;
  fHistLayer0 = new TH2F("fHistLayer0","Points in layer inner SPD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer0);
  zmax=14.;
  nbinsphi=40;nbinsz=4;
  fHistLayer1 = new TH2F("fHistLayer1","Points in layer outer SPD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer1);
  zmax=22.;
  nbinsphi=14;nbinsz=6;
  fHistLayer2 = new TH2F("fHistLayer2","Points in layer inner SDD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer2);
  zmax=29.5;
  nbinsphi=22;nbinsz=8;
  fHistLayer3 = new TH2F("fHistLayer3","Points in layer outer SDD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer3);
  zmax=45.;
  nbinsphi=34;nbinsz=23;
  fHistLayer4 = new TH2F("fHistLayer4","Points in layer inner SSD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer4);
  zmax=51.;
  nbinsphi=38;nbinsz=26;
  fHistLayer5 = new TH2F("fHistLayer5","Points in layer outer SSD; global   #phi; global z [cm]",nbinsphi,-3.14,3.14,nbinsz,-zmax,zmax);
  fListOfHistos->Add(fHistLayer5);


  fntExtra = new TNtuple("fntExtra","extra clusters in ITS","ncls:layer:ladder:volid:phi:x:y:z:xloc:zloc:dxy:dz:d0mu:z0mu:pt");
  fListOfHistos->Add(fntExtra);

  fntCosmicMatching = new TNtuple("fntCosmicMatching","cosmic tracks matching in ITS","ncls1:ncls2:pt1:pt2:sigmad01:sigmad02:sigmaz01:sigmaz02:dxy:dz:phimu:thetamu:d0mu:z0mu");
  fListOfHistos->Add(fntCosmicMatching);

  fspTree = new TTree("spTree","Tree with ITS track points");
  AliTrackPointArray *array = 0;
  AliESDVertex *vertex = 0;
  Float_t curv=0,curverr=0,runNumber=0;
  TObjString *itsaligndata = 0;
  TObjString *itscalibrespsdd = 0;
  fspTree->Branch("SP","AliTrackPointArray",&array);
  fspTree->Branch("vertex","AliESDVertex",&vertex);
  fspTree->Branch("curv",&curv);
  fspTree->Branch("curverr",&curverr);
  fspTree->Branch("run",&runNumber);
  fspTree->Branch("ITSAlignData",&itsaligndata);
  fspTree->Branch("ITSCalibRespSDD",&itscalibrespsdd);

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // write ITS AliTrackPoints for selected tracks to fspTree
  
  // check the geometry  
  if(!gGeoManager) { 
    printf("AliAlignmentDataFilterITS::Exec(): no geometry loaded \n");
    return;
  }

  // check if we have AliITSRecoParam
  if(!GetRecoParam()) {
    if(!fITSRecoParam) {
      printf("AliAlignmentDataFilterITS::Exec(): no AliITSRecoParam\n");
      return;
    }
  }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD) {
    printf("AliAlignmentDataFilterITS::Exec(): no ESD \n");
    return;
  } 

  //AliESDfriend *esdfriend = (AliESDfriend*)(fESD->FindListObject("AliESDfriend"));

  //if(!esdfriend) printf("AliAlignmentDataFilterITS::Exec(): no ESDfriend \n");
  //fESD->SetESDfriend(esdfriend);

  // Post the data for slot 0
  fHistNevents->Fill(0);


  // write field value to spTree UserInfo
  if(!((fspTree->GetUserInfo())->FindObject("BzkGauss"))) {
    Double_t bz=fESD->GetMagneticField();
    TString bzString; bzString+=bz;
    TObjString *bzObjString = new TObjString(bzString);
    TList *bzList = new TList();	 
    bzList->SetOwner(1);	 
    bzList->SetName("BzkGauss");	 
    bzList->Add(bzObjString);
    fspTree->GetUserInfo()->Add(bzList);
  }

  // write OCDB info to spTree UserInfo
  if(!((fspTree->GetUserInfo())->FindObject("cdbList"))) {
    TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
    if(!tree) {
      printf("ERROR: Could not read chain from input slot 0\n");
    } else {
      // Get the OCDB path and the list of OCDB objects used for reco 
      TMap *cdbMap = (TMap*)(tree->GetTree()->GetUserInfo())->FindObject("cdbMap");
      TList *cdbList = (TList*)(tree->GetTree()->GetUserInfo())->FindObject("cdbList");
      
      //cdbList->Print();
      // write the list to the user info of the output tree
      if(!fspTree) {
	printf("ERROR: fspTree does not exist\n");
      } else {
	TMap *cdbMapCopy = new TMap(cdbMap->GetEntries());	 
	cdbMapCopy->SetOwner(1);	 
	cdbMapCopy->SetName("cdbMap");	 
	TIter iter1(cdbMap->GetTable());	 
	
	TPair* pair = 0;	 
	while((pair = dynamic_cast<TPair*> (iter1.Next()))){	 
	  TObjString* keyStr = dynamic_cast<TObjString*> (pair->Key());	 
	  TObjString* valStr = dynamic_cast<TObjString*> (pair->Value());	 
	  cdbMapCopy->Add(new TObjString(keyStr->GetName()), new TObjString(valStr->GetName()));	 
	}	 
	
	TList *cdbListCopy = new TList();	 
	cdbListCopy->SetOwner(1);	 
	cdbListCopy->SetName("cdbList");	 
	
	TIter iter2(cdbList);	 
	
	TObjString* cdbEntry=0;
	while((cdbEntry =(TObjString*)(iter2.Next()))) {
	  cdbListCopy->Add(new TObjString(*cdbEntry));
	}	 
	cdbListCopy->Print();


	fspTree->GetUserInfo()->Add(cdbMapCopy);	 
	fspTree->GetUserInfo()->Add(cdbListCopy);
      }
    }
  }



  // Process event as Cosmic or Collision
  if(fESD->GetEventSpecie()<=1) {
    printf("AliAlignmentDataFilterITS::Exec(): event specie not set !\n");
    if(GetRecoParam()->GetAlignFilterCosmics()) {
      FilterCosmic(fESD);
    } else {
      FilterCollision(fESD);
    }
  } else if(fESD->GetEventSpecie()==8) {
    FilterCosmic(fESD);
  } else {
    FilterCollision(fESD);
  }

  PostData(2,fListOfHistos);

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::FilterCosmic(const AliESDEvent *esd)
{
  // Extract ITS AliTrackPoints for Cosmics (check angular matching
  // of top and bottom track, merge the two tracks, if requested)
  //

  // Set branch addresses for space points tree
  AliTrackPointArray *arrayForTree=0;
  AliESDVertex *vertexForTree=0;
  Float_t curv,curverr,runNumber;
  TObjString *itsaligndata=0;
  TObjString *itscalibrespsdd = 0;
  fspTree->SetBranchAddress("SP",&arrayForTree);
  fspTree->SetBranchAddress("vertex",&vertexForTree);
  fspTree->SetBranchAddress("curv",&curv);
  fspTree->SetBranchAddress("curverr",&curverr);
  fspTree->SetBranchAddress("run",&runNumber);
  fspTree->SetBranchAddress("ITSAlignData",&itsaligndata);
  fspTree->SetBranchAddress("ITSCalibRespSDD",&itscalibrespsdd);


  runNumber = (Float_t)esd->GetRunNumber();
  Int_t uid=10000+esd->GetEventNumberInFile();

  TTree* esdTree = dynamic_cast<TTree*> (GetInputData(0));
  // Get the list of OCDB objects used for reco 
  TList *cdbList = (TList*)(esdTree->GetTree()->GetUserInfo())->FindObject("cdbList");
  TIter iter2(cdbList);	     
  TObjString* cdbEntry=0;
  TString cdbEntryString;
  while((cdbEntry =(TObjString*)(iter2.Next()))) {
  cdbEntryString = cdbEntry->GetString();
  if(cdbEntryString.Contains("ITS/Align/Data")) 
    itsaligndata = new TObjString(*cdbEntry);
  if(cdbEntryString.Contains("ITS/Calib/RespSDD")) 
    itscalibrespsdd = new TObjString(*cdbEntry);
  }	 


  TString triggeredClass = esd->GetFiredTriggerClasses(); 
  if(fOnlySPDFO && !triggeredClass.Contains("C0SCO-ABCE-NOPF-CENT")) return;


  Int_t ntracks = esd->GetNumberOfTracks();
  if(ntracks<2) return;

  if(esd->GetPrimaryVertexSPD()->GetNContributors()<0) return;

  Double_t vtxpos[3]; esd->GetPrimaryVertexSPD()->GetXYZ(vtxpos);

  Int_t *goodtracksArray = new Int_t[ntracks];
  Float_t *phiArray = new Float_t[ntracks];
  Float_t *thetaArray = new Float_t[ntracks];
  Int_t *nclsArray = new Int_t[ntracks];
  Int_t ngt=0;
  Int_t itrack=0;
  for (itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack *track = esd->GetTrack(itrack);
    if (!track) continue;


    if(track->GetNcls(0)<GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    if((GetRecoParam()->GetAlignFilterOnlyITSSATracks()) && track->GetNcls(1)>0) continue;
    if((GetRecoParam()->GetAlignFilterOnlyITSTPCTracks()) && track->GetNcls(1)==0) continue;

    Float_t phi = track->GetAlpha()+TMath::ASin(track->GetSnp());
    Float_t theta = 0.5*TMath::Pi()-TMath::ATan(track->GetTgl());

    if(track->Pt()<GetRecoParam()->GetAlignFilterMinPt() || 
       track->Pt()>GetRecoParam()->GetAlignFilterMaxPt()) continue;

    goodtracksArray[ngt] = itrack;
    phiArray[ngt]        = phi;
    thetaArray[ngt]      = theta;
    nclsArray[ngt]       = track->GetNcls(0);
    ngt++;
  }

  if(ngt<2) {
    delete [] goodtracksArray; goodtracksArray=0;
    delete [] phiArray; phiArray=0;
    delete [] thetaArray; thetaArray=0;
    delete [] nclsArray; nclsArray=0;
    return;
  }

  // check matching of the two tracks from the muon
  Float_t min = 10000000.;
  Int_t maxCls = 0;
  Int_t good1 = -1, good2 = -1;
  for(Int_t itr1=0; itr1<ngt-1; itr1++) {
    for(Int_t itr2=itr1+1; itr2<ngt; itr2++) {
      Float_t deltatheta = TMath::Abs(TMath::Pi()-thetaArray[itr1]-thetaArray[itr2]);
      if(deltatheta>GetRecoParam()->GetAlignFilterMaxMatchingAngle()) continue;
      Float_t deltaphi = TMath::Abs(TMath::Abs(phiArray[itr1]-phiArray[itr2])-TMath::Pi());
      if(deltaphi>GetRecoParam()->GetAlignFilterMaxMatchingAngle()) continue;
      if(nclsArray[itr1]+nclsArray[itr2] > maxCls) {
	maxCls = nclsArray[itr1]+nclsArray[itr2];
	min = deltatheta+deltaphi;
	good1 = goodtracksArray[itr1];
	good2 = goodtracksArray[itr2];
      } else if(nclsArray[itr1]+nclsArray[itr2] == maxCls) {
	if(deltatheta+deltaphi < min) {
	  min = deltatheta+deltaphi;
	  good1 = goodtracksArray[itr1];
	  good2 = goodtracksArray[itr2];
	}
      }
    }
  }
  
  delete [] goodtracksArray; goodtracksArray=0;
  delete [] phiArray; phiArray=0;
  delete [] thetaArray; thetaArray=0;
  delete [] nclsArray; nclsArray=0;

  if(good1<0) return;
  AliDebug(2,"ok track matching");
  
  // track1 will be the inward track (top)
  // track2 the outward (bottom)
  AliESDtrack *track1=0; 
  AliESDtrack *track2=0;
  AliESDtrack *track = esd->GetTrack(good1);
  if(track->Py()>0) { 
    track1 = esd->GetTrack(good1);
    track2 = esd->GetTrack(good2);
  } else {
    track1 = esd->GetTrack(good2);
    track2 = esd->GetTrack(good1);
  }

  AliTrackPoint point;
  const AliTrackPointArray *array=0;
  Int_t ipt,volId,modId,layerId,lay,lad,det;
  Int_t jpt=0;
  Bool_t layerOK[6][2]; 
  Int_t nclsTrk[2]={0,0};

  for(Int_t l1=0;l1<6;l1++) for(Int_t l2=0;l2<2;l2++) layerOK[l1][l2]=kFALSE;
    
  for(itrack=0; itrack<2; itrack++) {
    if(itrack==0) {
      track = track1;
    } else {
      track = track2;
    }
    array = track->GetTrackPointArray();
    if(!array) {
      AliWarning("No tracks points avaialble");
      continue;
    }
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      if(volId<=0) continue;
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      AliDebug(2,Form("%d %d %d  %f\n",ipt,layerId-1,volId,TMath::Sqrt(point.GetX()*point.GetX()+point.GetY()*point.GetY())));
      if(point.IsExtra() && 
	 (GetRecoParam()->GetAlignFilterSkipExtra())) continue;
      if(layerId<1 || layerId>6) continue;
      if(!GetRecoParam()->GetAlignFilterUseLayer(layerId-1)) continue;
      // check minAngleWrtITSModulePlanes
      Double_t p[3]; track->GetDirection(p);
      TVector3 pvec(p[0],p[1],p[2]);
      Double_t rot[9]; AliGeomManager::GetOrigRotation(volId,rot);
      TVector3 normvec(rot[1],rot[4],rot[7]);
      Double_t angle = pvec.Angle(normvec);
      if(angle>0.5*TMath::Pi()) angle = TMath::Pi()-angle;
      angle = 0.5*TMath::Pi()-angle;
      if(angle<GetRecoParam()->GetAlignFilterMinAngleWrtModulePlanes()) continue;
      layerOK[layerId-1][itrack]=kTRUE;
      jpt++;
      nclsTrk[itrack]++;
    }
  }
  AliDebug(2,Form("nClsTrk1 %d nClsTrk2 %d\n",nclsTrk[0],nclsTrk[1]));
    
  // read ITS cluster maps
  Int_t map1[6],map2[6];
  for(Int_t ilay=0;ilay<6;ilay++) {
    map1[ilay]=0; map2[ilay]=0;
    if(track1->HasPointOnITSLayer(ilay)) map1[ilay]=1;
    if(track2->HasPointOnITSLayer(ilay)) map2[ilay]=1;
  }
  AliDebug(2,Form("ITS map 1: %d %d %d %d %d %d pt %f\n",map1[0],map1[1],map1[2],map1[3],map1[4],map1[5],track1->Pt()));
  AliDebug(2,Form("ITS map 2: %d %d %d %d %d %d pt %f\n",map2[0],map2[1],map2[2],map2[3],map2[4],map2[5],track2->Pt()));
  Int_t idx1[12],idx2[12];
  track1->GetITSclusters(idx1);
  track2->GetITSclusters(idx2);
  AliDebug(2,Form("cls idx 1 %d %d %d %d %d %d %d %d %d %d %d %d\n",idx1[0],idx1[1],idx1[2],idx1[3],idx1[4],idx1[5],idx1[6],idx1[7],idx1[8],idx1[9],idx1[10],idx1[11]));
  AliDebug(2,Form("cls idx 2 %d %d %d %d %d %d %d %d %d %d %d %d\n",idx2[0],idx2[1],idx2[2],idx2[3],idx2[4],idx2[5],idx2[6],idx2[7],idx2[8],idx2[9],idx2[10],idx2[11]));
  

  if(jpt<GetRecoParam()->GetAlignFilterMinITSPointsMerged()) return;
  AliDebug(2,Form(" Total points %d, accepted\n",jpt));  
  fHistNpoints->Fill(jpt);
  fHistPt->Fill(0.5*(track1->Pt()+track2->Pt()));
  
  Float_t d0z0mu[2];
  track1->GetDZ(0,0,0,esd->GetMagneticField(),d0z0mu);
  //printf("d0mu %f  z0mu %f\n",d0z0mu[0],d0z0mu[1]);

  vertexForTree = new AliESDVertex(*(esd->GetPrimaryVertexSPD()));

  Float_t dzOverlap[2];
  Float_t curvArray[2],curverrArray[2];
  Double_t globExtra[3],locExtra[3];
  if(GetRecoParam()->GetAlignFilterCosmicMergeTracks()) {
    arrayForTree = new AliTrackPointArray(jpt);
    arrayForTree->SetUniqueID(uid);
  }
  jpt=0;
  for(itrack=0; itrack<2; itrack++) {
    if(itrack==0) {
      track = track1;
    } else {
      track = track2;
    }
    curvArray[itrack] = track->GetC(esd->GetMagneticField());
    curverrArray[itrack] = TMath::Sqrt(track->GetSigma1Pt2())*track->GetC(esd->GetMagneticField())/track->OneOverPt();

    if(!(GetRecoParam()->GetAlignFilterCosmicMergeTracks())) {
      jpt=0;
      arrayForTree = new AliTrackPointArray(nclsTrk[itrack]);
      arrayForTree->SetUniqueID(uid);
    }
    array = track->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      if(volId<=0) continue;
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<1 || layerId>6 || !layerOK[layerId-1][itrack]) continue;
      arrayForTree->AddPoint(jpt,&point);
      jpt++;
      switch(layerId) {
      case 1:
	fHistLayer0->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 2:
	fHistLayer1->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 3:
	fHistLayer2->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 4:
	fHistLayer3->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 5:
	fHistLayer4->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 6:
	fHistLayer5->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      }
      // Post the data for slot 2
      if(jpt==1) PostData(2,fListOfHistos); // only if this is the first points
      if(!point.IsExtra() || 
	 !(GetRecoParam()->GetAlignFilterFillQANtuples())) continue;
      nclsTrk[itrack]--;
      for(Int_t ll=1;ll<layerId;ll++) modId+=AliITSgeomTGeo::GetNLadders(ll)*AliITSgeomTGeo::GetNDetectors(ll);
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      globExtra[0]=point.GetX();
      globExtra[1]=point.GetY();
      globExtra[2]=point.GetZ();
      AliITSgeomTGeo::GlobalToLocal(lay,lad,det,globExtra,locExtra);
      //printf("%d %d %d %d %d  %f %f %f\n",volId,modId,lay,lad,det,locExtra[0],locExtra[1],locExtra[2]);
      track->GetDZ(point.GetX(),point.GetY(),point.GetZ(),esd->GetMagneticField(),dzOverlap);
      AliTrackPoint pointT;
      Float_t radius,radiusT,phiv,phivT,thetav,thetavT;
      for(Int_t lll=0;lll<ipt;lll++) {
	array->GetPoint(pointT,lll);
	if(pointT.GetVolumeID()<=0) continue;
	Int_t layerIdT = AliGeomManager::VolUIDToLayer(pointT.GetVolumeID(),modId);
	if(layerIdT!=layerId) continue;
	radius=TMath::Sqrt((point.GetX()-vtxpos[0])*(point.GetX()-vtxpos[0])+(point.GetY()-vtxpos[1])*(point.GetY()-vtxpos[1]));
	radiusT=TMath::Sqrt((pointT.GetX()-vtxpos[0])*(pointT.GetX()-vtxpos[0])+(pointT.GetY()-vtxpos[1])*(pointT.GetY()-vtxpos[1]));
	phiv=TMath::ATan2(point.GetY()-vtxpos[1],point.GetX()-vtxpos[0]);
	phivT=TMath::ATan2(pointT.GetY()-vtxpos[1],pointT.GetX()-vtxpos[0]);
	if(TMath::Abs(point.GetZ()-vtxpos[2])<0.00001 || TMath::Abs(pointT.GetZ()-vtxpos[2])<0.00001) continue;
	thetav=TMath::ATan(radius/(point.GetZ()-vtxpos[2]));
	thetavT=TMath::ATan(radiusT/(pointT.GetZ()-vtxpos[2]));
	dzOverlap[0]=(Float_t)((phivT-phiv)*0.5*(radiusT+radius));
	if(TMath::Abs(TMath::Tan(0.5*(thetav+thetavT)))<0.00001) continue;
	dzOverlap[1]=(Float_t)((pointT.GetZ()-point.GetZ())-(radiusT-radius)/TMath::Tan(0.5*(thetav+thetavT)));
	fntExtra->Fill((Float_t)nclsTrk[itrack],(Float_t)(layerId-1),lad,volId,TMath::ATan2(point.GetY(),point.GetX()),point.GetX(),point.GetY(),point.GetZ(),locExtra[0],locExtra[2],dzOverlap[0],dzOverlap[1],d0z0mu[0],d0z0mu[1],track->Pt());
      }
    }

    if(!(GetRecoParam()->GetAlignFilterCosmicMergeTracks())) {
      curv = curvArray[itrack];
      curverr = curverrArray[itrack];
      fspTree->Fill();
    }
  }

  if(GetRecoParam()->GetAlignFilterCosmicMergeTracks()) {
    curv = 0.5*(curvArray[0]-curvArray[1]);  // the "-" is because the two tracks have opposite curvature!
    curverr = 0.5*TMath::Sqrt(curverrArray[0]*curverrArray[0]+curverrArray[1]*curverrArray[1]);
    fspTree->Fill();
  }
  PostData(1,fspTree);

  if(!(GetRecoParam()->GetAlignFilterFillQANtuples())) return; 
  // fill ntuple with track-to-track matching
  Float_t phimu,thetamu,phiout,thetaout,dphi,dtheta,rotymu,rotyout,droty;    
  Float_t d0[2],z0[2];
  Float_t sigmad0[2],sigmaz0[2];
  phimu = track1->GetAlpha()+TMath::ASin(track1->GetSnp());
  thetamu = 0.5*TMath::Pi()-TMath::ATan(track1->GetTgl());
  phiout = track2->GetAlpha()+TMath::ASin(track2->GetSnp());
  thetaout = 0.5*TMath::Pi()-TMath::ATan(track2->GetTgl());
  rotymu = TMath::ATan2(track1->Px(),track1->Pz());
  rotyout = TMath::ATan2(track2->Px(),track2->Pz());

  dphi = phimu - (phiout+TMath::Pi());
  dtheta = thetamu - (TMath::Pi()-thetaout);
  if(rotymu>0) {
    droty = rotymu - (rotyout+TMath::Pi());
  } else {
    droty = rotymu - (rotyout-TMath::Pi());
  }

  Double_t alpha = TMath::ATan2(track1->Py(),track1->Px());

  track1->Propagate(alpha,0.,esd->GetMagneticField());
  track2->Propagate(alpha,0.,esd->GetMagneticField());
  d0[0] = track1->GetY();
  z0[0] = track1->GetZ();
  d0[1] = track2->GetY();
  z0[1] = track2->GetZ();
  Float_t dxy = -(d0[0]-d0[1]);
  Float_t dz  = z0[0]-z0[1];
  sigmad0[0] = TMath::Sqrt(track1->GetSigmaY2());
  sigmaz0[0] = TMath::Sqrt(track1->GetSigmaZ2());
  sigmad0[1] = TMath::Sqrt(track2->GetSigmaY2());
  sigmaz0[1] = TMath::Sqrt(track2->GetSigmaZ2());
  /*  
  Double_t xyz1atxl0[3],xyz1atxl1[3],xyz2atxl0[3],xyz2atxl1[3];
  track1->GetXYZAt(0.,esd->GetMagneticField(),xyz1atxl0);
  track1->GetXYZAt(1.,esd->GetMagneticField(),xyz1atxl1);
  track2->GetXYZAt(0.,esd->GetMagneticField(),xyz2atxl0);
  track2->GetXYZAt(1.,esd->GetMagneticField(),xyz2atxl1);
  Float_t x1aty0 = (xyz1atxl0[0]*xyz1atxl1[1]-xyz1atxl0[1]*xyz1atxl1[0])/(xyz1atxl1[1]-xyz1atxl0[1]);
  Float_t x2aty0 = (xyz2atxl0[0]*xyz2atxl1[1]-xyz2atxl0[1]*xyz2atxl1[0])/(xyz2atxl1[1]-xyz2atxl0[1]);
  Float_t dxaty0 = x1aty0-x2aty0;
  */
  fntCosmicMatching->Fill((Float_t)nclsTrk[0],(Float_t)nclsTrk[1],track1->Pt(),track2->Pt(),sigmad0[0],sigmad0[1],sigmaz0[0],sigmaz0[1],dxy,dz,phimu,thetamu,TMath::Abs(d0z0mu[0]),d0z0mu[1]);
  
  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::FilterCollision(const AliESDEvent *esd)
{
  // Extract ITS AliTrackPoints for Cosmics (check angular matching
  // of top and bottom track, merge the two tracks, if requested)
  //

  // Set branch addresses for space points tree
  AliTrackPointArray *arrayForTree=0;
  AliESDVertex *vertexForTree=0;
  Float_t curv,curverr,runNumber;
  TObjString *itsaligndata=0;
  TObjString *itscalibrespsdd = 0;
  fspTree->SetBranchAddress("SP",&arrayForTree);
  fspTree->SetBranchAddress("vertex",&vertexForTree);
  fspTree->SetBranchAddress("curv",&curv);
  fspTree->SetBranchAddress("curverr",&curverr);
  fspTree->SetBranchAddress("run",&runNumber);
  fspTree->SetBranchAddress("ITSAlignData",&itsaligndata);
  fspTree->SetBranchAddress("ITSCalibRespSDD",&itscalibrespsdd);


  runNumber = (Float_t)esd->GetRunNumber();
  Int_t uid=20000+esd->GetEventNumberInFile();

  TTree* esdTree = dynamic_cast<TTree*> (GetInputData(0));
  // Get the list of OCDB objects used for reco 
  TList *cdbList = (TList*)(esdTree->GetTree()->GetUserInfo())->FindObject("cdbList");
  TIter iter2(cdbList);	     
  TObjString* cdbEntry=0;
  TString cdbEntryString;
  while((cdbEntry =(TObjString*)(iter2.Next()))) {
  cdbEntryString = cdbEntry->GetString();
  if(cdbEntryString.Contains("ITS/Align/Data")) 
    itsaligndata = new TObjString(*cdbEntry);
  if(cdbEntryString.Contains("ITS/Calib/RespSDD")) 
    itscalibrespsdd = new TObjString(*cdbEntry);
  }	 

  Int_t ntracks = esd->GetNumberOfTracks();

  if(ntracks==0) return;

  const AliESDVertex *vertexTracks = esd->GetPrimaryVertexTracks();
  if(!vertexTracks);
  if(vertexTracks->GetNContributors()<=0) return;

  Double_t vtxpos[3]; vertexTracks->GetXYZ(vtxpos);

  Int_t ncls=0;
  Double_t pt=-10000.;
  Double_t d0z0[2],covd0z0[3];
  const AliTrackPointArray *array = 0;

  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack);
    if (!track) continue;

    if(fDownsamplelowpt && TMath::Abs(esd->GetMagneticField())>0.01 &&
       track->Pt()<gRandom->Rndm()) continue;

    if(track->GetNcls(0)<GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    if((GetRecoParam()->GetAlignFilterOnlyITSSATracks()) && track->GetNcls(1)>0) continue;
    if((GetRecoParam()->GetAlignFilterOnlyITSTPCTracks()) && track->GetNcls(1)==0) continue;

    if(track->Pt()<GetRecoParam()->GetAlignFilterMinPt() || 
       track->Pt()>GetRecoParam()->GetAlignFilterMaxPt()) continue;

    pt = track->Pt();
    ncls = track->GetNcls(0);
    Double_t maxd=10000.;
    track->PropagateToDCA(vertexTracks,esd->GetMagneticField(),maxd,d0z0,covd0z0);

    // read ITS cluster map
    Int_t map[6];
    for(Int_t ilay=0;ilay<6;ilay++) {
      map[ilay]=0;
      if(track->HasPointOnITSLayer(ilay)) map[ilay]=1;
    }
    AliDebug(2,Form("ITS map : %d %d %d %d %d %d pt %f\n",map[0],map[1],map[2],map[3],map[4],map[5],track->Pt()));
    Int_t idx[12];
    track->GetITSclusters(idx);
    AliDebug(2,Form("cls idx %d %d %d %d %d %d %d %d %d %d %d %d\n",idx[0],idx[1],idx[2],idx[3],idx[4],idx[5],idx[6],idx[7],idx[8],idx[9],idx[10],idx[11]));
    

    AliTrackPoint point;
    Int_t ipt,volId,modId,layerId,lay,lad,det;
    Int_t jpt=0;
    Bool_t layerOK[6]; for(Int_t l1=0;l1<6;l1++) layerOK[l1]=kFALSE;
    
    array = track->GetTrackPointArray();
    if(!array) {printf("no track points\n"); continue;}
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      if(volId<=0) continue;
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<1 || layerId>6) continue;
      if(point.IsExtra() && 
	 (GetRecoParam()->GetAlignFilterSkipExtra())) continue;
      layerOK[layerId-1]=kTRUE;
      jpt++;
    }

    if(jpt < GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    fHistNpoints->Fill(jpt);
    fHistPt->Fill(pt);
    PostData(2,fListOfHistos);

    Float_t dzOverlap[2];
    Double_t globExtra[3],locExtra[3];
    arrayForTree = new AliTrackPointArray(jpt);
    arrayForTree->SetUniqueID(uid);
    jpt=0;
    array = track->GetTrackPointArray();
    if(!array) continue;
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<1 || layerId>6 || !layerOK[layerId-1]) continue;
      arrayForTree->AddPoint(jpt,&point);
      switch(layerId) {
      case 1:
	fHistLayer0->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 2:
	fHistLayer1->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 3:
	fHistLayer2->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 4:
	fHistLayer3->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 5:
	fHistLayer4->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      case 6:
	fHistLayer5->Fill(TMath::ATan2(point.GetY(),point.GetX()),point.GetZ());
	break;
      }
      jpt++;
      if(!point.IsExtra() || 
	 !(GetRecoParam()->GetAlignFilterFillQANtuples())) continue;
      ncls--;
      for(Int_t ll=1;ll<layerId;ll++) modId+=AliITSgeomTGeo::GetNLadders(ll)*AliITSgeomTGeo::GetNDetectors(ll);
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      globExtra[0]=point.GetX();
      globExtra[1]=point.GetY();
      globExtra[2]=point.GetZ();
      AliITSgeomTGeo::GlobalToLocal(lay,lad,det,globExtra,locExtra);
      track->GetDZ(point.GetX(),point.GetY(),point.GetZ(),esd->GetMagneticField(),dzOverlap);
      AliTrackPoint pointT;
      Float_t radius,radiusT,phiv,phivT,thetav,thetavT;
      for(Int_t lll=0;lll<ipt;lll++) {
	array->GetPoint(pointT,lll);
	Int_t layerIdT = AliGeomManager::VolUIDToLayer(pointT.GetVolumeID(),modId);
	if(layerIdT!=layerId) continue;
	radius=TMath::Sqrt((point.GetX()-vtxpos[0])*(point.GetX()-vtxpos[0])+(point.GetY()-vtxpos[1])*(point.GetY()-vtxpos[1]));
	radiusT=TMath::Sqrt((pointT.GetX()-vtxpos[0])*(pointT.GetX()-vtxpos[0])+(pointT.GetY()-vtxpos[1])*(pointT.GetY()-vtxpos[1]));
	phiv=TMath::ATan2(point.GetY()-vtxpos[1],point.GetX()-vtxpos[0]);
	phivT=TMath::ATan2(pointT.GetY()-vtxpos[1],pointT.GetX()-vtxpos[0]);
	if(TMath::Abs(point.GetZ()-vtxpos[2])<0.00001 || TMath::Abs(pointT.GetZ()-vtxpos[2])<0.00001) continue;
	thetav=TMath::ATan(radius/(point.GetZ()-vtxpos[2]));
	thetavT=TMath::ATan(radiusT/(pointT.GetZ()-vtxpos[2]));
	dzOverlap[0]=(Float_t)((phivT-phiv)*0.5*(radiusT+radius));
	if(TMath::Abs(TMath::Tan(0.5*(thetav+thetavT)))<0.00001) continue;
	dzOverlap[1]=(Float_t)((pointT.GetZ()-point.GetZ())-(radiusT-radius)/TMath::Tan(0.5*(thetav+thetavT)));
	fntExtra->Fill((Float_t)ncls,(Float_t)(layerId-1),lad,volId,TMath::ATan2(point.GetY(),point.GetX()),point.GetX(),point.GetY(),point.GetZ(),locExtra[0],locExtra[2],dzOverlap[0],dzOverlap[1],d0z0[0],d0z0[1],track->Pt());
      }
    }

    curv = track->GetC(esd->GetMagneticField());
    curverr = TMath::Sqrt(track->GetSigma1Pt2())*track->GetC(esd->GetMagneticField())/track->OneOverPt();

    vertexForTree = new AliESDVertex(*vertexTracks);
    if(vertexTracks->UsesTrack(track->GetID())) {
      vertexForTree->SetID(1);
    } else {
      vertexForTree->SetID(0);
    }

    fspTree->Fill();
 
  } // end of tracks loop

  PostData(1,fspTree);

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(2,"AliITSAlignmentDataFiler: Terminate() \n");

  fspTree = dynamic_cast<TTree*> (GetOutputData(1));
  if (!fspTree) {     
    printf("ERROR: fspTree not available\n");
    return;
  }

  fListOfHistos = dynamic_cast<TList*> (GetOutputData(2));
  if (!fListOfHistos) {     
    printf("ERROR: fListOfHistos not available\n");
    return;
  }

  fHistNevents = dynamic_cast<TH1F*>(fListOfHistos->FindObject("fHistNevents"));
  fHistNpoints = dynamic_cast<TH1F*>(fListOfHistos->FindObject("fHistNpoints"));
  fHistPt = dynamic_cast<TH1F*>(fListOfHistos->FindObject("fHistPt"));
  fHistLayer0 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer0"));
  fHistLayer1 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer1"));
  fHistLayer2 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer2"));
  fHistLayer3 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer3"));
  fHistLayer4 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer4"));
  fHistLayer5 = dynamic_cast<TH2F*>(fListOfHistos->FindObject("fHistLayer5"));
  fntExtra = dynamic_cast<TNtuple*>(fListOfHistos->FindObject("fntExtra"));
  fntCosmicMatching = dynamic_cast<TNtuple*>(fListOfHistos->FindObject("fntCosmicMatching"));



  return;
}
//-------------------------------------------------------------------------------
const AliITSRecoParam *AliAlignmentDataFilterITS::GetRecoParam() const 
{
  //
  // Return the ITSRecoParam object
  //
  if(AliITSReconstructor::GetRecoParam()) {
    return AliITSReconstructor::GetRecoParam();
  } else if(fITSRecoParam) {
    return fITSRecoParam;
  } else return NULL;
}
//--------------------------------------------------------------------------------
Int_t AliAlignmentDataFilterITS::WriteTrackPointsInIdealGeom(Char_t *fin, 
							     Char_t *fout,
							     Char_t *fmis,
							     Char_t *fgeo,
							     Bool_t prn)
{
  //
  // Convert AliTrackPoints in fin, reconstructed with fmis, back
  // to ideal geometry
  //
  // M. Lunardon
  //


  TGeoHMatrix deltahm;

  // Load geometry
  if (gSystem->AccessPathName(fgeo)) {
    printf("couldn't find geometry file %s - skipping...\n",fmis);
    return -1;
  }
  
  TFile *geofile=TFile::Open(fgeo);
  TGeoManager *fgGeometry=NULL;

  fgGeometry=(TGeoManager*)geofile->Get("ALICE");

  if (!fgGeometry)
    fgGeometry=(TGeoManager*)geofile->Get("Geometry");

  if (!fgGeometry) {
    AliCDBEntry *entry = (AliCDBEntry*)geofile->Get("AliCDBEntry");
    if (entry)
      fgGeometry = (TGeoManager*)entry->GetObject();
  }

  if (!fgGeometry) return -1;
  AliGeomManager::SetGeometry(fgGeometry);
  if(!AliGeomManager::GetGeometry()) return -1;
  
  
  // open alignment file
  if (gSystem->AccessPathName(fmis)) {
    printf("couldn't open alignment file %s - skipping...\n",fmis);
    return -2;
  }
  TFile *pref = TFile::Open(fmis);
  if (!pref->IsOpen()) return -2;
  
  
  /// apply alignment to ideal geometry
  TClonesArray *prea=(TClonesArray*)pref->Get("ITSAlignObjs");
  if (!prea) {
    if (pref->Get("AliCDBEntry"))
      prea = (TClonesArray*) ((AliCDBEntry*)pref->Get("AliCDBEntry"))->GetObject();
  }
  if (!prea) return -3;  
  Int_t nprea=prea->GetEntriesFast();
  printf("Array of input misalignments with %d entries\n",nprea);
  AliGeomManager::ApplyAlignObjsToGeom(*prea); // apply all levels of objs
  
  AliTrackPointArray *tpain=NULL;
  TFile *tpainfile=NULL;
  TTree *treein=NULL;
  AliTrackPoint point; 
  AliITSAlignMille2Module *m2[2200];
  for (Int_t i=0; i<2198; i++)
    m2[i]=new AliITSAlignMille2Module(AliITSAlignMille2Module::GetVolumeIDFromIndex(i));  
  
  // open input file
  if (gSystem->AccessPathName(fin)) {
    printf("couldn't open file %s - skipping...\n",fin);
    return -4;
  }
  tpainfile = TFile::Open(fin);
  if (!tpainfile->IsOpen()) return -4;
  
  treein=(TTree*)tpainfile->Get("spTree");
  if (!treein) return -5;
  Float_t curv,curverr,runNumber;
  TObjString *itsaligndata=0;
  TObjString *itscalibrespsdd = 0;
  treein->SetBranchAddress("SP", &tpain);
  treein->SetBranchAddress("curv", &curv);
  treein->SetBranchAddress("curverr", &curverr);
  treein->SetBranchAddress("run",&runNumber);
  treein->SetBranchAddress("ITSAlignData",&itsaligndata);
  treein->SetBranchAddress("ITSCalibRespSDD",&itscalibrespsdd);

  int ntrks=treein->GetEntries();
  printf("Reading %d tracks from %s\n",ntrks,fin);
  
  
  // open output file
  TFile *pointsFile = TFile::Open(fout,"RECREATE");
  if (!pointsFile || !pointsFile->IsOpen()) {
    printf("Can't open output file %s !",fout);
    return -6;
  }
  AliTrackPointArray *array = new AliTrackPointArray();
  
  // new!
  TTree *treeout=(TTree*)treein->Clone("spTree");
  treeout->Reset();
  treeout->SetBranchAddress("SP", &array);
  treeout->SetBranchAddress("curv", &curv);
  treeout->SetBranchAddress("curverr", &curverr);
  treeout->SetBranchAddress("run",&runNumber);
  treeout->SetBranchAddress("ITSAlignData",&itsaligndata);
  treeout->SetBranchAddress("ITSCalibRespSDD",&itscalibrespsdd);

  // tracks main loop
  for (Int_t it=0; it<ntrks; it++) {    
    if (!(it%5000) ) printf("...processing track n. %d\n",it);
    
    treein->GetEvent(it);
    
    //////////////////////////////
    
    AliTrackPointArray *atp=tpain;
    AliTrackPointArray *atps=NULL;
    Int_t npts=atp->GetNPoints();
    
    AliTrackPoint p;
    // check points in specific places
    
    // build a new track
    atps=new AliTrackPointArray(npts);
        
    Int_t npto=0;
    for (int i=0; i<npts; i++) {
      atp->GetPoint(p,i);
      
      UShort_t volid=atp->GetVolumeID()[i];
      Int_t index=AliITSAlignMille2Module::GetIndexFromVolumeID(volid);
      
      
      // dealign point
      // get MODIFIED matrix
      TGeoHMatrix *svMatrix = m2[index]->GetSensitiveVolumeMatrix(p.GetVolumeID());
      //TGeoHMatrix *svOrigMatrix = mm->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
      
      Double_t pg[3],pl[3];
      pg[0]=p.GetX();
      pg[1]=p.GetY();
      pg[2]=p.GetZ();
      if (prn) printf("Global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]);
      svMatrix->MasterToLocal(pg,pl);
      
      // check that things went OK: local y should be 0.
      if(TMath::Abs(pl[1])>1.e-6) {
	printf("AliAlignmentDataFilterITS::WriteTrackPointsInIdealGeom: ERROR, local y = %f (should be zero)\n",pl[1]);
	return -7;
      }

      if (prn) printf("Local coordinates of measured point : X=%f  Y=%f  Z=%f \n",pl[0],pl[1],pl[2]);
      
      // update covariance matrix
      TGeoHMatrix hcov;
      Double_t hcovel[9];
      hcovel[0]=(Double_t)(p.GetCov()[0]);
      hcovel[1]=(Double_t)(p.GetCov()[1]);
      hcovel[2]=(Double_t)(p.GetCov()[2]);
      hcovel[3]=(Double_t)(p.GetCov()[1]);
      hcovel[4]=(Double_t)(p.GetCov()[3]);
      hcovel[5]=(Double_t)(p.GetCov()[4]);
      hcovel[6]=(Double_t)(p.GetCov()[2]);
      hcovel[7]=(Double_t)(p.GetCov()[4]);
      hcovel[8]=(Double_t)(p.GetCov()[5]);
      hcov.SetRotation(hcovel);
      // now rotate in local system
      hcov.Multiply(svMatrix);
      hcov.MultiplyLeft(&svMatrix->Inverse());
      // now hcov is LOCAL COVARIANCE MATRIX
      
      /// get original matrix of sens. vol.
      TGeoHMatrix *svOrigMatrix = m2[index]->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
      // modify global coordinates according with pre-aligment
      svOrigMatrix->LocalToMaster(pl,pg);
      // now rotate in local system
      hcov.Multiply(&svOrigMatrix->Inverse());
      hcov.MultiplyLeft(svOrigMatrix);
      // hcov is back in GLOBAL RF
      Float_t pcov[6];
      pcov[0]=hcov.GetRotationMatrix()[0];
      pcov[1]=hcov.GetRotationMatrix()[1];
      pcov[2]=hcov.GetRotationMatrix()[2];
      pcov[3]=hcov.GetRotationMatrix()[4];
      pcov[4]=hcov.GetRotationMatrix()[5];
      pcov[5]=hcov.GetRotationMatrix()[8];
      
      p.SetXYZ(pg[0],pg[1],pg[2],pcov);
      if (prn) printf("New global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]);
      atps->AddPoint(npto,&p);
      if (prn) printf("Adding point[%d] = ( %f , %f , %f )     volid = %d\n",npto,atps->GetX()[npto],atps->GetY()[npto],atps->GetZ()[npto],atps->GetVolumeID()[npto] );
      if (prn) p.Print("");
      
      npto++;
    }
    
    
    ////////////////////////////////////////////////////////////
    array = atps;
    treeout->Fill();
    
    delete atps;
    atps=NULL;
    
  } // end loop on tracks
  
  pointsFile->Write();
  pointsFile->Close();
  tpainfile->Close();
  
  return 0;
}
