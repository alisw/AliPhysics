/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////
//  Stand alone ITS tracker class                        //
//  Origin:  Elisabetta Crescio - crescio@to.infn.it     //
//  Updated: Francesco Prino    - prino@to.infn.it       //
///////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TArrayI.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TTree.h>

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliITSVertexer.h"
#include "AliITSclusterTable.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITStrackSA.h"
#include "AliITStrackerSA.h"
#include "AliITSReconstructor.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliITStrackerSA)

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA():AliITStrackerMI(),
fPhiEstimate(0),
fITSStandAlone(0),
fLambdac(0),
fPhic(0),
fCoef1(0),
fCoef2(0),
fCoef3(0),
fNloop(0),
fPhiWin(0),
fLambdaWin(0),
fVert(0),
fVertexer(0),
fListOfTracks(0),
fListOfSATracks(0),
fITSclusters(0),
fInwardFlag(0),
fOuterStartLayer(0),
fInnerStartLayer(5),
fMinNPoints(0),
fMinQ(0.),
fCluLayer(0),
fCluCoord(0){
  // Default constructor
  Init();
 
}
//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(const Char_t *geom):AliITStrackerMI(0),
fPhiEstimate(0),
fITSStandAlone(0),
fLambdac(0),
fPhic(0),
fCoef1(0),
fCoef2(0),
fCoef3(0),
fNloop(0),
fPhiWin(0),
fLambdaWin(0),
fVert(0),
fVertexer(0),
fListOfTracks(0),
fListOfSATracks(0),
fITSclusters(0),
fInwardFlag(0),
fOuterStartLayer(0),
fInnerStartLayer(5),
fMinNPoints(0),
fMinQ(0.),
fCluLayer(0),
fCluCoord(0) 
{
  // Standard constructor (Vertex is known and passed to this obj.)
  if (geom) {
    AliWarning("\"geom\" is actually a dummy argument !");
  }

  Init();
  fVert = 0;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(const Char_t *geom, AliESDVertex *vert):AliITStrackerMI(0),
fPhiEstimate(0),
fITSStandAlone(0),
fLambdac(0),
fPhic(0),
fCoef1(0),
fCoef2(0),
fCoef3(0),
fNloop(0),
fPhiWin(0),
fLambdaWin(0),
fVert(vert),
fVertexer(0),
fListOfTracks(0),
fListOfSATracks(0),
fITSclusters(0),
fInwardFlag(0),
fOuterStartLayer(0),
fInnerStartLayer(5),
fMinNPoints(0),
fMinQ(0.),
fCluLayer(0),
fCluCoord(0)
{
  // Standard constructor (Vertex is known and passed to this obj.)
  if (geom) {
    AliWarning("\"geom\" is actually a dummy argument !");
  }
  Init();
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(const Char_t *geom, AliITSVertexer *vertexer):AliITStrackerMI(0),
fPhiEstimate(0),
fITSStandAlone(0),
fLambdac(0),
fPhic(0),
fCoef1(0),
fCoef2(0),
fCoef3(0),
fNloop(0),
fPhiWin(0),
fLambdaWin(0),
fVert(),
fVertexer(vertexer),
fListOfTracks(0),
fListOfSATracks(0),
fITSclusters(0),
fInwardFlag(0),
fOuterStartLayer(0),
fInnerStartLayer(5),
fMinNPoints(0),
fMinQ(0.),
fCluLayer(0),
fCluCoord(0)
{
  // Standard constructor (Vertex is unknown - vertexer is passed to this obj)
  if (geom) {
    AliWarning("\"geom\" is actually a dummy argument !");
  }
  Init();
  fVertexer = vertexer;
 
}
/*
//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(const AliITStrackerSA& tracker):AliITStrackerMI(),
fPhiEstimate(tracker.fPhiEstimate),
fITSStandAlone(tracker.fITSStandAlone),
fLambdac(tracker.fLambdac),
fPhic(tracker.fPhic),
fCoef1(tracker.fCoef1),
fCoef2(tracker.fCoef2),
fCoef3(tracker.fCoef3),
fNloop(tracker.fNloop),
fPhiWin(tracker.fPhiWin),
fLambdaWin(tracker.fLambdaWin),
fVert(tracker.fVert),
fVertexer(tracker.fVertexer),
fListOfTracks(tracker.fListOfTracks),
fListOfSATracks(tracker.fListOfSATracks),
fITSclusters(tracker.fITSclusters),
fInwardFlag(tracker.fInwardFlag),
fOuterStartLayer(tracker.fOuterStartLayer),
fInnerStartLayer(tracker.fInnerStartLayer),
fMinNPoints(tracker.fMinNPoints),
fMinQ(tracker.fMinQ),
fCluLayer(tracker.fCluLayer),
fCluCoord(tracker.fCluCoord) {
  // Copy constructor
  for(Int_t i=0;i<2;i++){
    fPoint1[i]=tracker.fPoint1[i];
    fPoint2[i]=tracker.fPoint2[i];
    fPoint3[i]=tracker.fPoint3[i];
    fPointc[i]=tracker.fPointc[i];
  }
  if(tracker.fVertexer && tracker.fVert){
    fVert = new AliESDVertex(*tracker.fVert);
  }
  else {
    fVert = tracker.fVert;
  }
  for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++){
    fCluLayer[i] = tracker.fCluLayer[i];
    fCluCoord[i] = tracker.fCluCoord[i];
  } 
}
//______________________________________________________________________
AliITStrackerSA& AliITStrackerSA::operator=(const AliITStrackerSA& source){
    // Assignment operator. 
  this->~AliITStrackerSA();
  new(this) AliITStrackerSA(source);
  return *this;
 
}
*/
//____________________________________________________________________________
AliITStrackerSA::~AliITStrackerSA(){
  // destructor
  // if fVertexer is not null, the AliESDVertex obj. is owned by this class
  // and is deleted here
  if(fVertexer){
    if(fVert)delete fVert;
  }
  fVert = 0;
  fVertexer = 0;
 
  if(fPhiWin)delete []fPhiWin;
  if(fLambdaWin)delete []fLambdaWin;
  fListOfTracks->Delete();
  delete fListOfTracks;
  fListOfSATracks->Delete();
  delete fListOfSATracks;
  if(fCluLayer){
    for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++){
      if(fCluLayer[i]){
	fCluLayer[i]->Delete();
	delete fCluLayer[i];
      }
    }
    delete [] fCluLayer;
  }
  if(fCluCoord){
    for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++){
      if(fCluCoord[i]){
	fCluCoord[i]->Delete();
	delete fCluCoord[i];
      }
    }
    delete [] fCluCoord;
  }
  
}

//____________________________________________________________________________
Int_t AliITStrackerSA::Clusters2Tracks(AliESDEvent *event){
// This method is used to find and fit the tracks. By default the corresponding
// method in the parent class is invoked. In this way a combined tracking
// TPC+ITS is performed. If the flag fITSStandAlone is true, the tracking
// is done in the ITS only. In the standard reconstruction chain this option
// can be set via AliReconstruction::SetOption("ITS","onlyITS")
  Int_t rc=0;

  if(!fITSStandAlone){
    rc=AliITStrackerMI::Clusters2Tracks(event);
  }
  else {
    AliDebug(1,"Stand Alone flag set: doing tracking in ITS alone\n");
  }
  if(!rc){ 
    rc=FindTracks(event,kFALSE);
    Int_t nSPDcontr=0;
    const AliESDVertex *spdv = event->GetPrimaryVertexSPD();
    if(spdv) nSPDcontr = spdv->GetNContributors();
    if(AliITSReconstructor::GetRecoParam()->GetSAUseAllClusters()==kTRUE && 
       nSPDcontr<=AliITSReconstructor::GetRecoParam()->GetMaxSPDcontrForSAToUseAllClusters()) {
      rc=FindTracks(event,kTRUE);
    }
  }
  return rc;
}

//____________________________________________________________________________
void AliITStrackerSA::Init(){
  //  Reset all data members
    fPhiEstimate=0;
    for(Int_t i=0;i<3;i++){fPoint1[i]=0;fPoint2[i]=0;fPoint3[i]=0;}
    fLambdac=0;
    fPhic=0;
    fCoef1=0;
    fCoef2=0;
    fCoef3=0;
    fPointc[0]=0;
    fPointc[1]=0;
    fVert = 0;
    fVertexer = 0;
    Int_t nLoops=AliITSReconstructor::GetRecoParam()->GetNLoopsSA();
    if(nLoops==33){
      SetFixedWindowSizes();
    }else{
      Double_t phimin=AliITSReconstructor::GetRecoParam()->GetMinPhiSA();
      Double_t phimax=AliITSReconstructor::GetRecoParam()->GetMaxPhiSA();
      Double_t lambmin=AliITSReconstructor::GetRecoParam()->GetMinLambdaSA();
      Double_t lambmax=AliITSReconstructor::GetRecoParam()->GetMaxLambdaSA();
      SetCalculatedWindowSizes(nLoops,phimin,phimax,lambmin,lambmax);
    }
    fMinQ=AliITSReconstructor::GetRecoParam()->GetSAMinClusterCharge();
    fITSclusters = 0;
    SetOuterStartLayer(1);
    SetSAFlag(kFALSE);
    fListOfTracks=new TClonesArray("AliITStrackMI",100);
    fListOfSATracks=new TClonesArray("AliITStrackSA",100);
    fCluLayer = 0;
    fCluCoord = 0;
    fMinNPoints = 3;
 }
//_______________________________________________________________________
void AliITStrackerSA::ResetForFinding(){
  //  Reset data members used in all loops during track finding
    fPhiEstimate=0;
    for(Int_t i=0;i<3;i++){fPoint1[i]=0;fPoint2[i]=0;fPoint3[i]=0;}
    fLambdac=0;
    fPhic=0;
    fCoef1=0;
    fCoef2=0;
    fCoef3=0;
    fPointc[0]=0;
    fPointc[1]=0;
    fListOfTracks->Clear();
    fListOfSATracks->Clear();
}

 

//______________________________________________________________________
Int_t AliITStrackerSA::FindTracks(AliESDEvent* event, Bool_t useAllClusters){

// Track finder using the ESD object

  AliDebug(2,Form(" field is %f",event->GetMagneticField()));
  AliDebug(2,Form("SKIPPING %d %d %d %d %d %d",ForceSkippingOfLayer(0),ForceSkippingOfLayer(1),ForceSkippingOfLayer(2),ForceSkippingOfLayer(3),ForceSkippingOfLayer(4),ForceSkippingOfLayer(5)));

  if(!fITSclusters){
    Fatal("FindTracks","ITS cluster tree is not accessed - Abort!!!\n Please use method SetClusterTree to pass the pointer to the tree\n");
    return -1;
  }
  //Reads event and mark clusters of traks already found, with flag kITSin
  Int_t nentr=event->GetNumberOfTracks();
  if(!useAllClusters) {
    while (nentr--) {
      AliESDtrack *track=event->GetTrack(nentr);
      if ((track->GetStatus()&AliESDtrack::kITSin) == AliESDtrack::kITSin){
	Int_t idx[12];
	Int_t ncl = track->GetITSclusters(idx);
	for(Int_t k=0;k<ncl;k++){
	  AliITSRecPoint* cll = (AliITSRecPoint*)GetCluster(idx[k]);
	  cll->SetBit(kSAflag);
	}
      }
    }
  }else{
    while (nentr--) {
      AliESDtrack *track=event->GetTrack(nentr);
      if ((track->GetStatus()&AliESDtrack::kITSin) == AliESDtrack::kITSin){
	Int_t idx[12];
	Int_t ncl = track->GetITSclusters(idx);
	for(Int_t k=0;k<ncl;k++){
	  AliITSRecPoint* cll = (AliITSRecPoint*)GetCluster(idx[k]);
	  cll->ResetBit(kSAflag);
	}
      }
    }
  }
  //Get primary vertex
  Double_t primaryVertex[3];
  event->GetVertex()->GetXYZ(primaryVertex);
  //Creates TClonesArray with clusters for each layer. The clusters already used
  //by AliITStrackerMI are not considered
  Int_t nclusters[AliITSgeomTGeo::kNLayers]={0,0,0,0,0,0};
  Int_t dmar[AliITSgeomTGeo::kNLayers]={0,0,0,0,0,0};
  if (fCluLayer == 0) {
    fCluLayer = new TClonesArray*[AliITSgeomTGeo::kNLayers];
    fCluCoord = new TClonesArray*[AliITSgeomTGeo::kNLayers];
    for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++) {
      fCluLayer[i]=0;
      fCluCoord[i]=0;
    }
  }
  for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++){
    AliITSlayer &layer=fgLayers[i];
    if (!ForceSkippingOfLayer(i)) {
      for(Int_t cli=0;cli<layer.GetNumberOfClusters();cli++){
	AliITSRecPoint* cls = (AliITSRecPoint*)layer.GetCluster(cli);
	if(cls->TestBit(kSAflag)==kTRUE) continue; //clusters used by TPC prol.
	if(cls->GetQ()==0) continue; //fake clusters dead zones
	if(i>1 && cls->GetQ()<=fMinQ) continue; // cut on SDD and SSD cluster charge
	nclusters[i]++;
      }
    }
    dmar[i]=0;
    if(!fCluLayer[i]){
      fCluLayer[i] = new TClonesArray("AliITSRecPoint",nclusters[i]);
    }else{
      fCluLayer[i]->Delete();
      fCluLayer[i]->Expand(nclusters[i]);
    }
    if(!fCluCoord[i]){
      fCluCoord[i] = new TClonesArray("AliITSclusterTable",nclusters[i]);
    }else{
      fCluCoord[i]->Delete();
      fCluCoord[i]->Expand(nclusters[i]);
    }
  }

  for(Int_t ilay=0;ilay<AliITSgeomTGeo::GetNLayers();ilay++){
    TClonesArray &clulay = *fCluLayer[ilay];
    TClonesArray &clucoo = *fCluCoord[ilay];
    AliITSlayer &layer=fgLayers[ilay];
    if (!ForceSkippingOfLayer(ilay)) {
      for(Int_t cli=0;cli<layer.GetNumberOfClusters();cli++){
	AliITSRecPoint* cls = (AliITSRecPoint*)layer.GetCluster(cli);
	if(cls->TestBit(kSAflag)==kTRUE) continue;
	if(cls->GetQ()==0) continue;
	if(ilay>1 && cls->GetQ()<=fMinQ) continue; 
	Double_t phi=0;Double_t lambda=0;
	Float_t x=0;Float_t y=0;Float_t z=0;
	Float_t sx=0;Float_t sy=0;Float_t sz=0;
	GetCoorAngles(cls,phi,lambda,x,y,z,primaryVertex);
	GetCoorErrors(cls,sx,sy,sz);
	new (clulay[dmar[ilay]]) AliITSRecPoint(*cls);
	new (clucoo[dmar[ilay]]) AliITSclusterTable(x,y,z,sx,sy,sz,phi,lambda,cli);
	dmar[ilay]++;
      }
    }
  }
   
  // track counter
  Int_t ntrack=0;

  static Int_t nClusLay[AliITSgeomTGeo::kNLayers];//counter for clusters on each layer
  Int_t startLayForSeed=0;
  Int_t lastLayForSeed=fOuterStartLayer;
  Int_t nSeedSteps=lastLayForSeed-startLayForSeed;
  Int_t seedStep=1;
  if(fInwardFlag){
    startLayForSeed=AliITSgeomTGeo::GetNLayers()-1;
    lastLayForSeed=fInnerStartLayer;
    nSeedSteps=startLayForSeed-lastLayForSeed;
    seedStep=-1;
  }

  // loop on minimum number of points
  for(Int_t iMinNPoints=AliITSgeomTGeo::GetNLayers(); iMinNPoints>=fMinNPoints; iMinNPoints--) {

    // loop on starting layer for track finding 
    for(Int_t iSeedLay=0; iSeedLay<=nSeedSteps; iSeedLay++) {
      Int_t theLay=startLayForSeed+iSeedLay*seedStep;
      if(ForceSkippingOfLayer(theLay)) continue;
      Int_t minNPoints=iMinNPoints-theLay;
      if(fInwardFlag) minNPoints=iMinNPoints-(AliITSgeomTGeo::GetNLayers()-1-theLay);
      for(Int_t i=theLay+1;i<AliITSgeomTGeo::GetNLayers();i++)
	if(ForceSkippingOfLayer(i)) 
	  minNPoints--;
      if(minNPoints<fMinNPoints) continue;

      // loop on phi and lambda window size
      for(Int_t nloop=0;nloop<fNloop;nloop++){
	Int_t nclTheLay=fCluLayer[theLay]->GetEntries();
	while(nclTheLay--){ 
	  ResetForFinding();
	  Bool_t useRP=SetFirstPoint(theLay,nclTheLay,primaryVertex);
	  if(!useRP) continue;	    
	  AliITStrackSA trs;
	    
	  Int_t pflag=0;	    
	  Int_t kk;
	  for(kk=0;kk<AliITSgeomTGeo::GetNLayers();kk++) nClusLay[kk] = 0;
	    
	  kk=0;
	  nClusLay[kk] = SearchClusters(theLay,fPhiWin[nloop],fLambdaWin[nloop],
					&trs,primaryVertex[2],pflag);
	  Int_t nextLay=theLay+seedStep;
	  Bool_t goon=kTRUE;
	  if(nextLay<0 || nextLay == 6) goon = kFALSE;
	  while(goon){
	    kk++;
	    nClusLay[kk] = SearchClusters(nextLay,fPhiWin[nloop],fLambdaWin[nloop],
					    &trs,primaryVertex[2],pflag);
	    if(nClusLay[kk]!=0){
	      pflag=1;
	      if(kk==1) {
		fPoint3[0]=fPointc[0];
		fPoint3[1]=fPointc[1];
	      } else {
		UpdatePoints();
	      }
	    }
	    nextLay+=seedStep;
	    if(nextLay<0 || nextLay==6) goon=kFALSE;
	  }

	    
	  Int_t layOK=0;
	  if(!fInwardFlag){
	    for(Int_t nnp=0;nnp<AliITSgeomTGeo::GetNLayers()-theLay;nnp++){
	      if(nClusLay[nnp]!=0) layOK+=1;
	    }
	  }else{
	    for(Int_t nnp=theLay; nnp>=0; nnp--){
	      if(nClusLay[nnp]!=0) layOK+=1;
	    }
	  }
	  if(layOK>=minNPoints){ 
	    AliDebug(2,Form("---NPOINTS: %d; MAP: %d %d %d %d %d %d\n",layOK,nClusLay[0],nClusLay[1],nClusLay[2],nClusLay[3],nClusLay[4],nClusLay[5]));
	    AliITStrackV2* tr2 = 0;
	    tr2 = FitTrack(&trs,primaryVertex);
	    if(!tr2){ 
	      continue;
	    }
	    AliDebug(2,Form("---NPOINTS fit: %d\n",tr2->GetNumberOfClusters()));
	      
	    StoreTrack(tr2,event,useAllClusters);
	    ntrack++;
	      
	  }   
	  
	}//end loop on clusters of theLay
      } //end loop on window sizes
    } //end loop on theLay
  }//end loop on min points

  // search for 1-point tracks in SPD, only for cosmics
  // (A.Dainese 21.03.08)
  if(AliITSReconstructor::GetRecoParam()->GetSAOnePointTracks() && 
     TMath::Abs(event->GetMagneticField())<0.01) {
    Int_t outerLayer=1; // only SPD
    for(Int_t innLay=0; innLay<=TMath::Min(1,fOuterStartLayer); innLay++) {
      //   counter for clusters on each layer  

      for(Int_t nloop=0;nloop<fNloop;nloop++){
	Int_t nclInnLay=fCluLayer[innLay]->GetEntries();
	while(nclInnLay--){ //loop starting from layer innLay
	  ResetForFinding();
	  Bool_t useRP=SetFirstPoint(innLay,nclInnLay,primaryVertex);
	  if(!useRP) continue;
	  AliITStrackSA trs;
	    
	  Int_t pflag=0;	    
	  Int_t kk;
	  for(kk=0;kk<AliITSgeomTGeo::GetNLayers();kk++) nClusLay[kk] = 0;
	  
	  kk=0;
	  nClusLay[kk] = SearchClusters(innLay,fPhiWin[nloop],fLambdaWin[nloop],
				  &trs,primaryVertex[2],pflag);
	  for(Int_t nextLay=innLay+1; nextLay<=outerLayer; nextLay++) {
	    kk++;
	    nClusLay[kk] = SearchClusters(nextLay,fPhiWin[nloop],fLambdaWin[nloop],
				    &trs,primaryVertex[2],pflag);
	    if(nClusLay[kk]!=0){
	      pflag=1;
	      if(kk==1) {
		fPoint3[0]=fPointc[0];
		fPoint3[1]=fPointc[1];
	      } else {
		UpdatePoints();
	      }
	    }
	  }
	  
	  Int_t layOK=0;
	  for(Int_t nnp=0;nnp<AliITSgeomTGeo::GetNLayers()-innLay;nnp++){
	    if(nClusLay[nnp]!=0) layOK+=1;
	  }
	  if(layOK==1) {
	    AliDebug(2,Form("----NPOINTS: %d; MAP: %d %d %d %d %d %d\n",layOK,nClusLay[0],nClusLay[1],nClusLay[2],nClusLay[3],nClusLay[4],nClusLay[5]));
	    AliITStrackV2* tr2 = 0;
	    Bool_t onePoint = kTRUE;
	    tr2 = FitTrack(&trs,primaryVertex,onePoint);
	    if(!tr2){
	      continue;
	    }
	    AliDebug(2,Form("----NPOINTS fit: %d\n",tr2->GetNumberOfClusters()));
	    
	    StoreTrack(tr2,event,useAllClusters);
	    ntrack++;
	    
	  }   
	  
	}//end loop on clusters of innLay
      } //end loop on window sizes
      
    } //end loop on innLay
  } // end search 1-point tracks
  
  if(!useAllClusters) AliInfo(Form("Number of found tracks: %d",event->GetNumberOfTracks()));
  ResetForFinding();
  return 0;

}
 
//________________________________________________________________________

AliITStrackV2* AliITStrackerSA::FitTrack(AliITStrackSA* tr,Double_t *primaryVertex,Bool_t onePoint) {
  //fit of the found track (most general case, also <6 points, layers missing)
  // A.Dainese 16.11.07 

  
  const Int_t kMaxClu=AliITStrackSA::kMaxNumberOfClusters;

  static Int_t firstmod[AliITSgeomTGeo::kNLayers];  
  static Int_t clind[AliITSgeomTGeo::kNLayers][kMaxClu];
  static Int_t clmark[AliITSgeomTGeo::kNLayers][kMaxClu];
  static Int_t end[AliITSgeomTGeo::kNLayers];
  static Int_t indices[AliITSgeomTGeo::kNLayers];

  static AliITSRecPoint *listlayer[AliITSgeomTGeo::kNLayers][kMaxClu];

  for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++) {
    firstmod[i]=AliITSgeomTGeo::GetModuleIndex(i+1,1,1);
    end[i]=0;
    for(Int_t j=0;j<kMaxClu; j++){
      clind[i][j]=0;
      clmark[i][j]=0;
      listlayer[i][j]=0;
   }
  }
  

  Int_t nclusters = tr->GetNumberOfClustersSA();
  for(Int_t ncl=0;ncl<nclusters;ncl++){
    Int_t index = tr->GetClusterIndexSA(ncl); 
    AliITSRecPoint* cl = (AliITSRecPoint*)GetCluster(index);
    Int_t lay = (index & 0xf0000000) >> 28;
    Int_t nInLay=end[lay];
    listlayer[lay][nInLay]=cl;
    clind[lay][nInLay]=index;
    end[lay]++;
  }

  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){
    for(Int_t ncl=0;ncl<tr->GetNumberOfMarked(nlay);ncl++){
      Int_t mark = tr->GetClusterMark(nlay,ncl);
      clmark[nlay][ncl]=mark;
    }
  }


  Int_t firstLay=-1,secondLay=-1;
  for(Int_t i=0;i<AliITSgeomTGeo::GetNLayers();i++) {
    if(end[i]==0) {
      end[i]=1;
    }else{
      if(firstLay==-1) {
	firstLay=i;
      } else if(secondLay==-1) {
	secondLay=i;
      }
    }
  }

  if(firstLay==-1 || (secondLay==-1 && !onePoint)) return 0;
  TClonesArray &arrMI= *fListOfTracks;
  TClonesArray &arrSA= *fListOfSATracks;
  Int_t nFoundTracks=0;


  for(Int_t l0=0;l0<end[0];l0++){ //loop on layer 1
    indices[0]=l0;
    for(Int_t l1=0;l1<end[1];l1++){ //loop on layer 2
      indices[1]=l1;
      for(Int_t l2=0;l2<end[2];l2++){  //loop on layer 3
	indices[2]=l2;
        for(Int_t l3=0;l3<end[3];l3++){ //loop on layer 4   
	  indices[3]=l3;
          for(Int_t l4=0;l4<end[4];l4++){ //loop on layer 5
	    indices[4]=l4;
            for(Int_t l5=0;l5<end[5];l5++){ //loop on layer 6  
	      indices[5]=l5;

	      // estimate curvature from 2 innermost points (or innermost point + vertex)

	      Int_t iFirstLay=indices[firstLay];
	      Int_t mrk1=clmark[firstLay][iFirstLay];

	      AliITSRecPoint* p1=(AliITSRecPoint*)listlayer[firstLay][iFirstLay];
	      Int_t module1 = p1->GetDetectorIndex()+firstmod[firstLay]; 
	      Int_t layer,ladder,detector;
	      AliITSgeomTGeo::GetModuleId(module1,layer,ladder,detector);
	      Float_t yclu1 = p1->GetY();
	      Float_t zclu1 = p1->GetZ();

	      Double_t x1,y1,z1;
	      Double_t x2,y2,z2;
	      Double_t cv=0,tgl2=0,phi2=0;
	      AliITSclusterTable* arr1 = (AliITSclusterTable*)GetClusterCoord(firstLay,mrk1);
	      x1 = arr1->GetX();
	      y1 = arr1->GetY();
	      z1 = arr1->GetZ();

	      if(secondLay>0) {
		Int_t iSecondLay=indices[secondLay];	      
		Int_t mrk2=clmark[secondLay][iSecondLay];
		AliITSclusterTable* arr2 = (AliITSclusterTable*)GetClusterCoord(secondLay,mrk2);
		x2 = arr2->GetX();
		y2 = arr2->GetY();
		z2 = arr2->GetZ();
		cv = Curvature(primaryVertex[0],primaryVertex[1],x1,y1,x2,y2);
		tgl2 = (z2-z1)/TMath::Sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		phi2 = TMath::ATan2((y2-y1),(x2-x1));
	      } else { // special case of 1-point tracks, only for cosmics (B=0)
		x2 = primaryVertex[0];
		y2 = primaryVertex[1];
		z2 = primaryVertex[2];
		cv = 0;
		tgl2 = (z1-z2)/TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		phi2 = TMath::ATan2((y1-y2),(x1-x2));
	      }

	      // create track and attach it the RecPoints
              AliITStrackSA trac(layer,ladder,detector,yclu1,zclu1,phi2,tgl2,cv,1);
	      for(Int_t iLay=5; iLay>=0; iLay--){
		Int_t iInLay=indices[iLay];
		AliITSRecPoint* cl=(AliITSRecPoint*)listlayer[iLay][iInLay];
		if(cl!=0){
		  trac.AddClusterV2(iLay,(clind[iLay][iInLay] & 0x0fffffff)>>0);
		  trac.AddClusterMark(iLay,clmark[iLay][iInLay]);
		}
	      }

              //fit with Kalman filter using AliITStrackerMI::RefitAt()
	      AliITStrackSA ot(trac);

              ot.ResetCovariance(10.);
              ot.ResetClusters();
              
	      // Propagate inside the innermost layer with a cluster 
	      if(ot.Propagate(ot.GetX()-0.1*ot.GetX())) {

		if(RefitAt(AliITSRecoParam::GetrInsideITSscreen(),&ot,&trac)){ //fit from layer 1 to layer 6
		  AliITStrackMI otrack2(ot);
		  otrack2.ResetCovariance(10.); 
		  otrack2.ResetClusters();
		  //fit from layer 6 to layer 1
		  if(RefitAt(AliITSRecoParam::GetrInsideSPD1(),&otrack2,&ot)) {
		    new(arrMI[nFoundTracks]) AliITStrackMI(otrack2);
		    new(arrSA[nFoundTracks]) AliITStrackSA(trac);
		    ++nFoundTracks;
		  }
                              
		}       
	      }
            }//end loop layer 6
          }//end loop layer 5
        }//end loop layer 4        
      }//end loop layer 3
    }//end loop layer 2 
  }//end loop layer 1




  if(fListOfTracks->GetEntries()==0) return 0;

  Int_t lowchi2 = FindTrackLowChiSquare();
  AliITStrackV2* otrack =(AliITStrackV2*)fListOfTracks->At(lowchi2);
  AliITStrackSA* trsa = (AliITStrackSA*)fListOfSATracks->At(lowchi2);
 
  if(otrack==0) return 0;

  CookLabel(otrack,0.); //MI change - to see fake ratio
  Int_t label=FindLabel(otrack);
  otrack->SetLabel(label);  
  Double_t low=0.;
  Double_t up=0.51;    
  otrack->CookdEdx(low,up);

  //remove clusters of found track
  for(Int_t nlay=0;nlay<AliITSgeomTGeo::GetNLayers();nlay++){
    for(Int_t cln=0;cln<trsa->GetNumberOfMarked(nlay);cln++){
      Int_t index = trsa->GetClusterMark(nlay,cln);
      fCluLayer[nlay]->RemoveAt(index);
      RemoveClusterCoord(nlay,index);
      fCluLayer[nlay]->Compress();
    }    
  }

  return otrack;

}

//_______________________________________________________
void AliITStrackerSA::StoreTrack(AliITStrackV2 *t,AliESDEvent *event, Bool_t pureSA) const 
{
  //
  // Add new track to the ESD
  //
  AliESDtrack outtrack;
  outtrack.UpdateTrackParams(t,AliESDtrack::kITSin);
  if(pureSA) outtrack.SetStatus(AliESDtrack::kITSpureSA);
  for(Int_t i=0;i<12;i++) {
    outtrack.SetITSModuleIndex(i,t->GetModuleIndex(i));
  }
  Double_t sdedx[4]={0.,0.,0.,0.};
  for(Int_t i=0; i<4; i++) sdedx[i]=t->GetSampledEdx(i);
  outtrack.SetITSdEdxSamples(sdedx);


  if(AliITSReconstructor::GetRecoParam()->GetSAUsedEdxInfo()){
    Double_t mom=t->P();
    Double_t ppid[AliPID::kSPECIES];
    for(Int_t isp=0;isp<AliPID::kSPECIES;isp++) ppid[isp]=0.;
    ppid[AliPID::kPion]=1.;
    if(mom<0.7){
      Double_t truncmean=t->GetdEdx();
      Int_t ide=fITSPid->GetParticleIdFromdEdxVsP(mom,truncmean,kTRUE);
      if(ide==AliPID::kProton){
	ppid[AliPID::kProton]=1.;
	ppid[AliPID::kPion]=0.;
      }
      else if(ide==AliPID::kKaon){ 
	ppid[AliPID::kKaon]=1.; 
	ppid[AliPID::kPion]=0.;
      }
    }
    outtrack.SetITSpid(ppid);
    outtrack.SetESDpid(ppid);    
  }
  event->AddTrack(&outtrack);

  return;
}


//_______________________________________________________
Int_t AliITStrackerSA::SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, AliITStrackSA* trs,Double_t /*zvertex*/,Int_t pflag){
  //function used to to find the clusters associated to the track

  if(ForceSkippingOfLayer(layer)) return 0;

  Int_t nc=0;
  AliITSlayer &lay = fgLayers[layer];
  Double_t r=lay.GetR();
  if(pflag==1){      
    Float_t cx1,cx2,cy1,cy2;
    FindEquation(fPoint1[0],fPoint1[1],fPoint2[0],fPoint2[1],fPoint3[0],fPoint3[1],fCoef1,fCoef2,fCoef3);
    if (FindIntersection(fCoef1,fCoef2,fCoef3,-r*r,cx1,cy1,cx2,cy2)==0)
       return 0;
    Double_t fi1=TMath::ATan2(cy1-fPoint1[1],cx1-fPoint1[0]);
    Double_t fi2=TMath::ATan2(cy2-fPoint1[1],cx2-fPoint1[0]);
    fPhiEstimate=ChoosePoint(fi1,fi2,fPhic);
  }

 
  Double_t phiExpect=fPhiEstimate;
  Double_t lamExpect=fLambdac;

  Int_t ncl = fCluLayer[layer]->GetEntriesFast();
  for (Int_t index=0; index<ncl; index++) {
    AliITSRecPoint *c = (AliITSRecPoint*)fCluLayer[layer]->UncheckedAt(index);
    if (!c) continue;
    
    AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(layer,index);

    Double_t lambda = arr->GetLambda();
    if (TMath::Abs(lambda-lamExpect)>lambdawindow) continue;

    Double_t phi = arr->GetPhi();
    Double_t deltaPhi = phi-phiExpect;
    if(deltaPhi>TMath::Pi()) deltaPhi-=2*TMath::Pi();
    else if(deltaPhi<-TMath::Pi()) deltaPhi+=2*TMath::Pi();
    if (TMath::Abs(deltaPhi)>phiwindow) continue;
    
    if(trs->GetNumberOfClustersSA()==trs->GetMaxNumberOfClusters()) return 0;
    if(trs->GetNumberOfMarked(layer)==trs->GetMaxNMarkedPerLayer()) return 0;
    Int_t orind = arr->GetOrInd();
    trs->AddClusterSA(layer,orind);
    trs->AddClusterMark(layer,index);
    nc++;
    fLambdac=lambda;
    fPhiEstimate=phi;
    
    fPointc[0]=arr->GetX();
    fPointc[1]=arr->GetY();
    
  }
  return nc;
}

//________________________________________________________________
Bool_t AliITStrackerSA::SetFirstPoint(Int_t lay, Int_t clu, Double_t* primaryVertex){
  // Sets the first point (seed) for tracking

  AliITSRecPoint* cl = (AliITSRecPoint*)fCluLayer[lay]->UncheckedAt(clu);
  if(!cl) return kFALSE;
  if (cl->GetQ()<=0) return kFALSE;
  if(lay>1 && cl->GetQ()<=fMinQ) return kFALSE;

  AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(lay,clu);
  fPhic = arr->GetPhi();
  fLambdac = arr->GetLambda();
  fPhiEstimate = fPhic;
  fPoint1[0]=primaryVertex[0];
  fPoint1[1]=primaryVertex[1];
  fPoint2[0]=arr->GetX();
  fPoint2[1]=arr->GetY();
  return kTRUE; 
}

//________________________________________________________________
void AliITStrackerSA::UpdatePoints(){
  //update of points for the estimation of the curvature  

  fPoint2[0]=fPoint3[0];
  fPoint2[1]=fPoint3[1];
  fPoint3[0]=fPointc[0];
  fPoint3[1]=fPointc[1];

  
}

//___________________________________________________________________
Int_t AliITStrackerSA::FindEquation(Float_t x1, Float_t y1, Float_t x2, Float_t y2, Float_t x3, Float_t y3,Float_t& a, Float_t& b, Float_t& c){

   //given (x,y) of three recpoints (in global coordinates) 
   //returns the parameters a,b,c of circonference x*x + y*y +a*x + b*y +c

   Float_t den = (x3-x1)*(y2-y1)-(x2-x1)*(y3-y1);
   if(den==0) return 0;
   a = ((y3-y1)*(x2*x2+y2*y2-x1*x1-y1*y1)-(y2-y1)*(x3*x3+y3*y3-x1*x1-y1*y1))/den;
   b = -(x2*x2-x1*x1+y2*y2-y1*y1+a*(x2-x1))/(y2-y1);
   c = -x1*x1-y1*y1-a*x1-b*y1;
   return 1;
 }
//__________________________________________________________________________
 Int_t AliITStrackerSA::FindIntersection(Float_t a1, Float_t b1, Float_t c1, Float_t c2,Float_t& x1,Float_t& y1, Float_t& x2, Float_t& y2){
 
 //Finds the intersection between the circonference of the track and the circonference centered in (0,0) represented by one layer
 //c2 is -rlayer*rlayer

  if(a1==0) return 0;
 Double_t m = c2-c1; 
 Double_t aA = (b1*b1)/(a1*a1)+1;
 Double_t bB = (-2*m*b1/(a1*a1));
 Double_t cC = c2+(m*m)/(a1*a1);
 Double_t dD = bB*bB-4*aA*cC;
 if(dD<0) return 0;
 
 y1 = (-bB+TMath::Sqrt(dD))/(2*aA); 
 y2 = (-bB-TMath::Sqrt(dD))/(2*aA); 
 x1 = (c2-c1-b1*y1)/a1;
 x2 = (c2-c1-b1*y2)/a1;

 return 1; 
}
//____________________________________________________________________
Double_t AliITStrackerSA::Curvature(Double_t x1,Double_t y1,Double_t 
x2,Double_t y2,Double_t x3,Double_t y3){

  //calculates the curvature of track  
  Double_t den = (x3-x1)*(y2-y1)-(x2-x1)*(y3-y1);
  if(den==0) return 0;
  Double_t a = ((y3-y1)*(x2*x2+y2*y2-x1*x1-y1*y1)-(y2-y1)*(x3*x3+y3*y3-x1*x1-y1*y1))/den;
  Double_t b = -(x2*x2-x1*x1+y2*y2-y1*y1+a*(x2-x1))/(y2-y1);
  Double_t c = -x1*x1-y1*y1-a*x1-b*y1;
  Double_t xc=-a/2.;

  if((a*a+b*b-4*c)<0) return 0;
  Double_t rad = TMath::Sqrt(a*a+b*b-4*c)/2.;
  if(rad==0) return 0;
  
  if((x1>0 && y1>0 && x1<xc)) rad*=-1;
  if((x1<0 && y1>0 && x1<xc)) rad*=-1;
  //  if((x1<0 && y1<0 && x1<xc)) rad*=-1;
  // if((x1>0 && y1<0 && x1<xc)) rad*=-1;
  
  return 1/rad;
 
}


//____________________________________________________________________
Double_t AliITStrackerSA::ChoosePoint(Double_t p1, Double_t p2, Double_t pp){

  //Returns the point closest to pp

  Double_t diff1 = p1-pp;
  Double_t diff2 = p2-pp;
  
  if(TMath::Abs(diff1)<TMath::Abs(diff2)) fPhiEstimate=p1;
  else fPhiEstimate=p2;  
  return fPhiEstimate;
  
}


//_________________________________________________________________
Int_t AliITStrackerSA::FindTrackLowChiSquare() const {
  // returns track with lowest chi square  
  Int_t dim=fListOfTracks->GetEntries();
  if(dim<=1) return 0;
  AliITStrackV2* trk = (AliITStrackV2*)fListOfTracks->At(0);
  Double_t minChi2=trk->GetChi2();
  Int_t index=0;
  for(Int_t i=1;i<dim;i++){
    trk = (AliITStrackV2*)fListOfTracks->At(i);
    Double_t chi2=trk->GetChi2();
    if(chi2<minChi2){
      minChi2=chi2;
      index=i;
    }
  }
  return index;
}

//__________________________________________________________
Int_t AliITStrackerSA::FindLabel(AliITStrackV2* track){
  //
  
  Int_t labl[AliITSgeomTGeo::kNLayers][3];
  Int_t cnts[AliITSgeomTGeo::kNLayers][3];
  for(Int_t j=0;j<AliITSgeomTGeo::GetNLayers();j++){
    for(Int_t k=0;k<3;k++){
      labl[j][k]=-2;
      cnts[j][k]=1;
    }
  }
  Int_t iNotLabel=0;
  for(Int_t i=0;i<track->GetNumberOfClusters(); i++) {
    Int_t indexc = track->GetClusterIndex(i);
    AliITSRecPoint* cl = (AliITSRecPoint*)GetCluster(indexc);
    Int_t iLayer=cl->GetLayer();
    for(Int_t k=0;k<3;k++){
      labl[iLayer][k]=cl->GetLabel(k);
      if(labl[iLayer][k]<0) iNotLabel++;
    }
  }
  if(iNotLabel==3*track->GetNumberOfClusters()) return -2;

  for(Int_t j1=0;j1<AliITSgeomTGeo::kNLayers; j1++) {
    for(Int_t j2=0; j2<j1;  j2++){
      for(Int_t k1=0; k1<3; k1++){
	for(Int_t k2=0; k2<3; k2++){
	  if(labl[j1][k1]>=0 && labl[j1][k1]==labl[j2][k2] && cnts[j2][k2]>0){
	    cnts[j2][k2]++;
	    cnts[j1][k1]=0;
	  }
	}
      }
    }
  }


  Int_t cntMax=0;
  Int_t label=-1;
  for(Int_t j=0;j<AliITSgeomTGeo::kNLayers;j++){
    for(Int_t k=0;k<3;k++){
      if(cnts[j][k]>cntMax && labl[j][k]>=0){
	cntMax=cnts[j][k];
	label=labl[j][k];
      }
    }
  }

  Int_t lflag=0;
  for(Int_t i=0;i<AliITSgeomTGeo::kNLayers;i++)
    if(labl[i][0]==label || labl[i][1]==label || labl[i][2]==label) lflag++;
  
  if(lflag<track->GetNumberOfClusters()) label = -label;
  return label;
}
//_____________________________________________________________________________
void AliITStrackerSA::SetCalculatedWindowSizes(Int_t n, Float_t phimin, Float_t phimax, Float_t lambdamin, Float_t lambdamax){
  // Set sizes of the phi and lambda windows used for track finding
  fNloop = n;
  if(fPhiWin) delete [] fPhiWin;
  if(fLambdaWin) delete [] fLambdaWin;
  fPhiWin = new Double_t[fNloop];
  fLambdaWin = new Double_t[fNloop];
  Float_t stepPhi=(phimax-phimin)/(Float_t)(fNloop-1);
  Float_t stepLambda=(lambdamax-lambdamin)/(Float_t)(fNloop-1);
  for(Int_t k=0;k<fNloop;k++){
    Float_t phi=phimin+k*stepPhi;
    Float_t lam=lambdamin+k*stepLambda;
    fPhiWin[k]=phi;
    fLambdaWin[k]=lam;
  }
}
//_____________________________________________________________________________
void AliITStrackerSA::SetFixedWindowSizes(Int_t n, Double_t *phi, Double_t *lam){
  // Set sizes of the phi and lambda windows used for track finding
  fNloop = n;
  if(phi){ // user defined values
    fPhiWin = new Double_t[fNloop];
    fLambdaWin = new Double_t[fNloop];
    for(Int_t k=0;k<fNloop;k++){
      fPhiWin[k]=phi[k];
      fLambdaWin[k]=lam[k];
    }
  }
  else {  // default values
            
    Double_t phid[33]   = {0.002,0.003,0.004,0.0045,0.0047,
			   0.005,0.0053,0.0055,
			   0.006,0.0063,0.0065,0.007,0.0073,0.0075,0.0077,
			   0.008,0.0083,0.0085,0.0087,0.009,0.0095,0.0097,
			   0.01,0.0105,0.011,0.0115,0.012,0.0125,0.013,0.0135,0.0140,0.0145};
    Double_t lambdad[33] = {0.003,0.004,0.005,0.005,0.005,
			    0.005,0.005,0.006,
			    0.006,0.006,0.006,0.007,0.007,0.007,0.007,
			    0.007,0.007,0.007,0.007,0.007,0.007,0.007,
			    0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008};
    
    if(fNloop!=33){
      fNloop = 33;
    }
    
    
    fPhiWin = new Double_t[fNloop];
    fLambdaWin = new Double_t[fNloop];

    Double_t factor=AliITSReconstructor::GetRecoParam()->GetFactorSAWindowSizes(); // possibility to enlarge windows for cosmics reco with large misalignments (A.Dainese)
  
    for(Int_t k=0;k<fNloop;k++){
      fPhiWin[k]=phid[k]*factor;
      fLambdaWin[k]=lambdad[k]*factor;
    }
  
  }

}
//_______________________________________________________________________
void AliITStrackerSA::GetCoorAngles(AliITSRecPoint* cl,Double_t &phi,Double_t &lambda, Float_t &x, Float_t &y,Float_t &z,Double_t* vertex){
  //Returns values of phi (azimuthal) and lambda angles for a given cluster
/*  
  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; fGeom->GetModuleId(module,lay,lad,det);
  Float_t tx,ty,tz;  fGeom->GetTrans(lay,lad,det,tx,ty,tz);     

  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi1=TMath::Pi()/2+alpha;
  if (lay==1) phi1+=TMath::Pi();

  Float_t cp=TMath::Cos(phi1), sp=TMath::Sin(phi1);
  Float_t r=tx*cp+ty*sp;

  xyz= r*cp - cl->GetY()*sp;
  y= r*sp + cl->GetY()*cp;
  z=cl->GetZ();
*/
  Float_t xyz[3];
  cl->GetGlobalXYZ(xyz);
  x=xyz[0];
  y=xyz[1];
  z=xyz[2];
 
  phi=TMath::ATan2(y-vertex[1],x-vertex[0]);
  lambda=TMath::ATan2(z-vertex[2],TMath::Sqrt((x-vertex[0])*(x-vertex[0])+(y-vertex[1])*(y-vertex[1])));
}

//________________________________________________________________________
void AliITStrackerSA::GetCoorErrors(AliITSRecPoint* cl,Float_t &sx,Float_t &sy, Float_t &sz){

  //returns sigmax, y, z of cluster in global coordinates
/*
  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; 
  AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
 
  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi=TMath::Pi()/2+alpha;
  if (lay==1) phi+=TMath::Pi();

  Float_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);
*/
  Float_t covm[6];
  cl->GetGlobalCov(covm);
  sx=TMath::Sqrt(covm[0]);
  sy=TMath::Sqrt(covm[3]);
  sz=TMath::Sqrt(covm[5]);
/*
  sx = TMath::Sqrt(sp*sp*cl->GetSigmaY2());
  sy = TMath::Sqrt(cp*cp*cl->GetSigmaY2());
  sz = TMath::Sqrt(cl->GetSigmaZ2());
*/
}

