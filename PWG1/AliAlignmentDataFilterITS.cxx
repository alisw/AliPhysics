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
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>
#include <TGeoManager.h>

#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliITSReconstructor.h"
#include "AliITSgeomTGeo.h"
#include "AliTrackPointArray.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAlignmentDataFilterITS.h"


ClassImp(AliAlignmentDataFilterITS)


//________________________________________________________________________
AliAlignmentDataFilterITS::AliAlignmentDataFilterITS(const char *name):
AliAnalysisTask(name,"task"),
fESD(0),
fESDfriend(0),
fListOfHistos(0),
fspTree(0),
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

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TTree
  DefineOutput(0,TTree::Class());  //My private output
  // Output slot #1 writes into a TList
  DefineOutput(1,TList::Class());  //My private output
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
void AliAlignmentDataFilterITS::ConnectInputData(Option_t *) 
{
  // Connect ESD
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if(!tree) {
    printf("ERROR: Could not read chain from input slot 0\n");
  } else {
    // Disable all branches and enable only the needed ones

    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);

    tree->SetBranchStatus("ESDfriend*", 1);
    tree->SetBranchAddress("ESDfriend.",&fESDfriend);

    tree->SetBranchStatus("fSPDMult*", 1);
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(!esdH) {
      printf("ERROR: Could not get ESDInputHandler\n");
    } else {
      fESD = esdH->GetEvent();
    }
  }
  
  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::Init()
{
  // Initialization

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::CreateOutputObjects()
{
  // Create the output container
  //

  // Several histograms are more conveniently managed in a TList
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();

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


  fntExtra = new TNtuple("fntExtra","extra clusters in ITS","ncls:layer:ladder:volid:phi:x:y:z:xloc:zloc:dxy:dz:d0mu:z0mu");
  fListOfHistos->Add(fntExtra);

  fntCosmicMatching = new TNtuple("fntCosmicMatching","cosmic tracks matching in ITS","ncls1:ncls2:pt1:pt2:sigmad01:sigmad02:sigmaz01:sigmaz02:dxy:dz:phimu:thetamu:d0mu:z0mu");
  fListOfHistos->Add(fntCosmicMatching);

  fspTree = new TTree("spTree","Tree with ITS track points");
  const AliTrackPointArray *array = 0;
  Float_t curv,curverr;
  fspTree->Branch("SP","AliTrackPointArray",&array);
  fspTree->Branch("curv",&curv);
  fspTree->Branch("curverr",&curverr);

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::Exec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // write ITS AliTrackPoints for selected tracks to fspTree
  
  if(!gGeoManager) {
    printf("AliAlignmentDataFilterITS::Exec(): no geometry loaded \n");
    return;
  }

  if(!fESD) {
    printf("AliAlignmentDataFilterITS::Exec(): no ESD \n");
    return;
  } 
  if(!fESDfriend) {
    printf("AliAlignmentDataFilterITS::Exec(): no ESDfriend \n");
    return;
  } 
  // attach ESDfriend
  fESD->SetESDfriend(fESDfriend);


  // Process event as Cosmic or Collision
  //if(esd->GetEventType()== ???? ) {
  printf("AliAlignmentDataFilterITS::Exec(): MOVE ASAP TO esd->GetEventType() !\n");
  if(AliITSReconstructor::GetRecoParam()->GetAlignFilterCosmics()) {
    FilterCosmic(fESD);
  } else {
    FilterCollision(fESD);
  }

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
  Float_t curv,curverr;
  fspTree->SetBranchAddress("SP",&arrayForTree);
  fspTree->SetBranchAddress("curv",&curv);
  fspTree->SetBranchAddress("curverr",&curverr);


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


    if(track->GetNcls(0)<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    if(AliITSReconstructor::GetRecoParam()->GetAlignFilterOnlyITSSATracks() && track->GetNcls(1)>0) continue;
    if(AliITSReconstructor::GetRecoParam()->GetAlignFilterOnlyITSTPCTracks() && track->GetNcls(1)==0) continue;

    Float_t phi = track->GetAlpha()+TMath::ASin(track->GetSnp());
    Float_t theta = 0.5*TMath::Pi()-TMath::ATan(track->GetTgl());

    if(track->Pt()<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinPt() || 
       track->Pt()>AliITSReconstructor::GetRecoParam()->GetAlignFilterMaxPt()) continue;

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
      if(deltatheta>AliITSReconstructor::GetRecoParam()->GetAlignFilterMaxMatchingAngle()) continue;
      Float_t deltaphi = TMath::Abs(TMath::Abs(phiArray[itr1]-phiArray[itr2])-TMath::Pi());
      if(deltaphi>AliITSReconstructor::GetRecoParam()->GetAlignFilterMaxMatchingAngle()) continue;
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
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      AliDebug(2,Form("%d %d\n",ipt,layerId-1));
      if(point.IsExtra() && 
	 AliITSReconstructor::GetRecoParam()->GetAlignFilterSkipExtra()) continue;
      if(layerId>6) continue;
      if(!AliITSReconstructor::GetRecoParam()->GetAlignFilterUseLayer(layerId-1)) continue;
      // check minAngleWrtITSModulePlanes
      Double_t p[3]; track->GetDirection(p);
      TVector3 pvec(p[0],p[1],p[2]);
      Double_t rot[9]; AliGeomManager::GetOrigRotation(volId,rot);
      TVector3 normvec(rot[1],rot[4],rot[7]);
      Double_t angle = pvec.Angle(normvec);
      if(angle>0.5*TMath::Pi()) angle = TMath::Pi()-angle;
      angle = 0.5*TMath::Pi()-angle;
      if(angle<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinAngleWrtModulePlanes()) continue;
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
  

  if(jpt<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinITSPointsMerged()) return;
  AliDebug(2,Form(" Total points %d, accepted\n",jpt));  
  fHistNpoints->Fill(jpt);
  fHistPt->Fill(0.5*(track1->Pt()+track2->Pt()));
  
  Float_t d0z0mu[2];
  track1->GetDZ(0,0,0,esd->GetMagneticField(),d0z0mu);
  //printf("d0mu %f  z0mu %f\n",d0z0mu[0],d0z0mu[1]);

  Float_t dzOverlap[2];
  Float_t curvArray[2],curverrArray[2];
  Double_t globExtra[3],locExtra[3];
  if(AliITSReconstructor::GetRecoParam()->GetAlignFilterCosmicMergeTracks()) 
    arrayForTree = new AliTrackPointArray(jpt);
  
  jpt=0;
  for(itrack=0; itrack<2; itrack++) {
    if(itrack==0) {
      track = track1;
    } else {
      track = track2;
    }
    curvArray[itrack] = track->GetC(esd->GetMagneticField());
    curverrArray[itrack] = TMath::Sqrt(track->GetSigma1Pt2())*track->GetC(esd->GetMagneticField())/track->OneOverPt();

    if(!AliITSReconstructor::GetRecoParam()->GetAlignFilterCosmicMergeTracks()) {
      jpt=0;
      arrayForTree = new AliTrackPointArray(nclsTrk[itrack]);
    }
    array = track->GetTrackPointArray();
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId>6 || !layerOK[layerId-1][itrack]) continue;
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
      // Post the data for slot 0
      if(jpt==1) PostData(1,fListOfHistos); // only if this is the first points
      if(!point.IsExtra() || 
	 !AliITSReconstructor::GetRecoParam()->GetAlignFilterFillQANtuples()) continue;
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
	fntExtra->Fill((Float_t)nclsTrk[itrack],(Float_t)(layerId-1),lad,volId,TMath::ATan2(point.GetY(),point.GetX()),point.GetX(),point.GetY(),point.GetZ(),locExtra[0],locExtra[2],dzOverlap[0],dzOverlap[1],d0z0mu[0],d0z0mu[1]);
      }
    }

    if(!AliITSReconstructor::GetRecoParam()->GetAlignFilterCosmicMergeTracks()) {
      curv = curvArray[itrack];
      curverr = curverrArray[itrack];
      fspTree->Fill();
    }
  }

  if(AliITSReconstructor::GetRecoParam()->GetAlignFilterCosmicMergeTracks()) {
    curv = 0.5*(curvArray[0]+curvArray[1]);
    curverr = 0.5*TMath::Sqrt(curverrArray[0]*curverrArray[0]+curverrArray[1]*curverrArray[1]);
    fspTree->Fill();
  }
  PostData(0,fspTree);

  if(!AliITSReconstructor::GetRecoParam()->GetAlignFilterFillQANtuples()) return; 
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
  Float_t curv,curverr;
  fspTree->SetBranchAddress("SP",&arrayForTree);
  fspTree->SetBranchAddress("curv",&curv);
  fspTree->SetBranchAddress("curverr",&curverr);

  Int_t ntracks = esd->GetNumberOfTracks();

  if(ntracks==0) return;

  if(esd->GetPrimaryVertexTracks()->GetNContributors()<=0) return;

  Double_t vtxpos[3]; esd->GetPrimaryVertexTracks()->GetXYZ(vtxpos);

  Int_t ncls=0;
  Double_t pt=-10000.;
  Double_t d0z0[2],covd0z0[3];
  const AliTrackPointArray *array = 0;

  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack);
    if (!track) continue;

    if(track->GetNcls(0)<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    if(AliITSReconstructor::GetRecoParam()->GetAlignFilterOnlyITSSATracks() && track->GetNcls(1)>0) continue;
    if(AliITSReconstructor::GetRecoParam()->GetAlignFilterOnlyITSTPCTracks() && track->GetNcls(1)==0) continue;

    if(track->Pt()<AliITSReconstructor::GetRecoParam()->GetAlignFilterMinPt() || 
       track->Pt()>AliITSReconstructor::GetRecoParam()->GetAlignFilterMaxPt()) continue;

    pt = track->Pt();
    ncls = track->GetNcls(0);
    Double_t maxd=10000.;
    track->PropagateToDCA(esd->GetPrimaryVertex(),esd->GetMagneticField(),maxd,d0z0,covd0z0);

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
    if(!array) continue;
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<1 || layerId>6) continue;
      if(point.IsExtra() && 
	 AliITSReconstructor::GetRecoParam()->GetAlignFilterSkipExtra()) continue;
      layerOK[layerId-1]=kTRUE;
      jpt++;
    }

    if(jpt < AliITSReconstructor::GetRecoParam()->GetAlignFilterMinITSPoints()) continue;

    fHistNpoints->Fill(jpt);
    fHistPt->Fill(pt);
    PostData(1,fListOfHistos);

    Float_t dzOverlap[2];
    Double_t globExtra[3],locExtra[3];
    arrayForTree = new AliTrackPointArray(jpt);
    jpt=0;
    array = track->GetTrackPointArray();
    if(!array) continue;
    for(ipt=0; ipt<array->GetNPoints(); ipt++) {
      array->GetPoint(point,ipt);
      volId = point.GetVolumeID();
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<1 || layerId>6 || !layerOK[layerId-1]) continue;
      if(!point.IsExtra() || 
	 !AliITSReconstructor::GetRecoParam()->GetAlignFilterFillQANtuples()) continue;
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
	fntExtra->Fill((Float_t)ncls,(Float_t)(layerId-1),lad,volId,TMath::ATan2(point.GetY(),point.GetX()),point.GetX(),point.GetY(),point.GetZ(),locExtra[0],locExtra[2],dzOverlap[0],dzOverlap[1],d0z0[0],d0z0[1]);
      }
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
    }

    curv = track->GetC(esd->GetMagneticField());
    curverr = TMath::Sqrt(track->GetSigma1Pt2())*track->GetC(esd->GetMagneticField())/track->OneOverPt();

    fspTree->Fill();
 
  } // end of tracks loop

  PostData(0,fspTree);

  return;
}

//________________________________________________________________________
void AliAlignmentDataFilterITS::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(2,"AliITSAlignmentDataFiler: Terminate() \n");

  fspTree = dynamic_cast<TTree*> (GetOutputData(0));
  if (!fspTree) {     
    printf("ERROR: fspTree not available\n");
    return;
  }

  fListOfHistos = dynamic_cast<TList*> (GetOutputData(1));
  if (!fListOfHistos) {     
    printf("ERROR: fListOfHistos not available\n");
    return;
  }

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

