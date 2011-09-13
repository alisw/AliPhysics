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

////////////////////////////////////////////////////
//  ITS Upgrade Stand alone tracker class         //
//  Authors: A.Mastroserio C.Terrevoli            //
//  e-mail:annalisa.mastroserio@cern.ch           //
//         cristina.terrevoli@ba.infn.it          //
////////////////////////////////////////////////////

#include <stdlib.h>

#include <TArrayI.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TTree.h>
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliITSVertexer.h"
#include "AliITSclusterTable.h"
#include "AliITSRecPointU.h"
#include "AliITStrackU.h"
#include "AliITStrackerMI.h"
#include "AliITSlayerUpgrade.h"
#include "AliITSsegmentationUpgrade.h"
#include "AliITStrackerUpgrade.h"
#include "AliITSReconstructor.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliITStrackerUpgrade)

//____________________________________________________________________________
  AliITStrackerUpgrade::AliITStrackerUpgrade():
    AliITStrackerMI(),
    fNLayers(), fPhiEstimate(0), fITSStandAlone(0), fLambdac(0),
    fPhic(0), fCoef1(0), fCoef2(0), fCoef3(0), fNloop(0), 
    fPhiWin(0),  fLambdaWin(0), fVert(0), fVertexer(0),
    fListOfTracks(0), fListOfUTracks(0), fITSclusters(0),
    fInwardFlag(0), fOuterStartLayer(0),  fInnerStartLayer(0),
    fMinNPoints(0),  fMinQ(0), fTrackToFollow(), fLayers(0), fSegmentation(0x0),
    fCluLayer(0), fCluCoord(0)
{
  // Default constructor
  Init();
 
}
//____________________________________________________________________________
AliITStrackerUpgrade::AliITStrackerUpgrade(Int_t nLay):
  AliITStrackerMI(),
  fNLayers(nLay), fPhiEstimate(0), fITSStandAlone(0), fLambdac(0),
  fPhic(0), fCoef1(0), fCoef2(0), fCoef3(0), fNloop(0), 
  fPhiWin(0),  fLambdaWin(0), fVert(0), fVertexer(0),
  fListOfTracks(0), fListOfUTracks(0), fITSclusters(0),
  fInwardFlag(0), fOuterStartLayer(0),  fInnerStartLayer(nLay),
  fMinNPoints(0),  fMinQ(0),  fTrackToFollow(), fLayers(0), fSegmentation(0x0),
  fCluLayer(0), fCluCoord(0)
{
  //
  // constructor
  //
  Init();
 
}
  
//________________________________________________________________________
AliITStrackerUpgrade::AliITStrackerUpgrade(const AliITStrackerUpgrade& tracker):
  AliITStrackerMI(),
  fNLayers(tracker.fNLayers), fPhiEstimate(tracker.fPhiEstimate),
  fITSStandAlone(tracker.fITSStandAlone), fLambdac(tracker.fLambdac),
  fPhic(tracker.fPhic), fCoef1(tracker.fCoef1), fCoef2(tracker.fCoef2),
  fCoef3(tracker.fCoef3), fNloop(tracker.fNloop), fPhiWin(tracker.fPhiWin),
  fLambdaWin(tracker.fLambdaWin), fVert(tracker.fVert), fVertexer(tracker.fVertexer),
  fListOfTracks(tracker.fListOfTracks), fListOfUTracks(tracker.fListOfUTracks),
  fITSclusters(tracker.fITSclusters), fInwardFlag(tracker.fInwardFlag),
  fOuterStartLayer(tracker.fOuterStartLayer), fInnerStartLayer(tracker.fInnerStartLayer),
  fMinNPoints(tracker.fMinNPoints), fMinQ(tracker.fMinQ), fTrackToFollow(tracker.fTrackToFollow), fLayers(tracker.fLayers),
  fSegmentation(tracker.fSegmentation), fCluLayer(tracker.fCluLayer), fCluCoord(tracker.fCluCoord) 
{
  //									
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
  for(Int_t i=0;i<fNLayers;i++){
    fCluLayer[i] = tracker.fCluLayer[i];
    fCluCoord[i] = tracker.fCluCoord[i];
  }
}
//______________________________________________________________________
AliITStrackerUpgrade& AliITStrackerUpgrade::operator=(const AliITStrackerUpgrade& source){
  // Assignment operator. 
  this->~AliITStrackerUpgrade();
  new(this) AliITStrackerUpgrade(source);
  return *this;
 
}

//____________________________________________________________________________
AliITStrackerUpgrade::~AliITStrackerUpgrade(){
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
  fListOfUTracks->Delete();
  delete fListOfUTracks;
  if(fCluLayer){
    for(Int_t i=0;i<fNLayers;i++){
      if(fCluLayer[i]){
	fCluLayer[i]->Delete();
	delete fCluLayer[i];
      }
    }
    delete [] fCluLayer;
  }
  if(fCluCoord){
    for(Int_t i=0;i<fNLayers;i++){
      if(fCluCoord[i]){
	fCluCoord[i]->Delete();
	delete fCluCoord[i];
      }
    }
    delete [] fCluCoord;
  }
  
  if(fSegmentation) delete fSegmentation;
}

//____________________________________________________________________________
Int_t AliITStrackerUpgrade::Clusters2Tracks(AliESDEvent *event){
  // This method is used to find and fit the tracks. By default the corresponding
  // method in the parent class is invoked. In this way a combined tracking
  // TPC+ITS is performed. If the flag fITSStandAlone is true, the tracking
  // is done in the ITS only. In the standard reconstruction chain this option
  // can be set via AliReconstruction::SetOption("ITS","onlyITS")
  Int_t rc=0;
  if(!fITSStandAlone){
    rc=AliITStrackerMI::Clusters2Tracks(event); // should not be used !!
    // Note standalone tracking is implememted in "AliITStrackerU"
  }
  else {
    AliDebug(1,"Stand Alone flag set: doing tracking in ITS alone\n");
  }
  if(!rc){
    rc=FindTracks(event,kFALSE);
    if(AliITSReconstructor::GetRecoParam()->GetSAUseAllClusters()==kTRUE) {
      rc=FindTracks(event,kTRUE); // eventually a 2nd time ??????
    }
  }
  return rc;
}
//_________________________________________________-
void AliITStrackerUpgrade::Init(){
  //  Reset all data members
  const Int_t nL = fNLayers;
  fPhiEstimate=0;
  for(Int_t i=0;i<2;i++){fPoint1[i]=0;fPoint2[i]=0;fPoint3[i]=0;}
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
  fListOfTracks=new TClonesArray("AliITStrackU",100);
  fListOfUTracks=new TClonesArray("AliITStrackU",100);
  fLayers=new AliITSlayerUpgrade*[nL];//to be fixed
  fCluLayer = 0;
  fCluCoord = 0;
  fMinNPoints = 3;
  for(Int_t layer=0; layer<fNLayers; layer++){
    Double_t p=0.;
    Double_t zC= 0.;
    fLayers[layer] = new AliITSlayerUpgrade(p,zC);
  }

  fSegmentation = new AliITSsegmentationUpgrade();
	
}
//_______________________________________________________________________
Int_t AliITStrackerUpgrade::LoadClusters(TTree *clusTree){ 
  //
  // Load clusters for tracking
  // 

  TClonesArray statITSCluster("AliITSRecPointU");
  TClonesArray *ITSCluster = &statITSCluster;

  TBranch* itsClusterBranch=clusTree->GetBranch("ITSRecPoints");
  if (!itsClusterBranch){
    AliError("can't get the branch with the ITS clusters ! \n");
    return 1;
  }
  itsClusterBranch->SetAddress(&ITSCluster);
  clusTree->GetEvent(0);
  Int_t nCluster = ITSCluster->GetEntriesFast();
  for(Int_t i=0; i<nCluster; i++){
    AliITSRecPointU *recp = (AliITSRecPointU*)ITSCluster->UncheckedAt(i);
    fLayers[recp->GetLayer()]->InsertCluster(new AliITSRecPointU(*recp));
  }//loop clusters

  SetClusterTree(clusTree);
  return 0;
}

//_______________________________________________________
void AliITStrackerUpgrade::ResetForFinding(){
  //  Reset data members used in all loops during track finding
  fPhiEstimate=0;
  for(Int_t i=0;i<2;i++){fPoint1[i]=0;fPoint2[i]=0;fPoint3[i]=0;}
  fLambdac=0;
  fPhic=0;
  fCoef1=0;
  fCoef2=0;
  fCoef3=0;
  fPointc[0]=0;
  fPointc[1]=0;
  fListOfTracks->Clear();
  fListOfUTracks->Clear();
}

//______________________________________________________________________
Int_t AliITStrackerUpgrade::FindTracks(AliESDEvent* event,Bool_t useAllClusters){

  // Track finder using the ESD object
  AliDebug(2,Form(" field is %f",event->GetMagneticField()));

  if(!fITSclusters){
    Fatal("FindTracks","ITS cluster tree is not accessed - Abort!!!\n Please use method SetClusterTree to pass the pointer to the tree\n");
    return -1;
  }

  // Reads event and mark clusters of tracks already found, with flag kITSin
  Int_t nentr=event->GetNumberOfTracks();
  if(!useAllClusters) {
    while (nentr--) {
      AliESDtrack *track=event->GetTrack(nentr);
      if ((track->GetStatus()&AliESDtrack::kITSin) == AliESDtrack::kITSin){
        Int_t idx[12];
        Int_t ncl = track->GetITSclusters(idx);
        for(Int_t k=0;k<ncl;k++){
          AliITSRecPointU* cll = (AliITSRecPointU*)GetCluster(idx[k]);
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
          AliITSRecPointU* cll = (AliITSRecPointU*)GetCluster(idx[k]);
          cll->ResetBit(kSAflag);
        }
      }
    }
  }
  //Get primary vertex
  Double_t primaryVertex[3];
  event->GetVertex()->GetXYZ(primaryVertex);
  Int_t nclusters[fNLayers];
  for(Int_t i=0; i<fNLayers; i++){
    nclusters[i]=0;
  }
  
  Int_t dmar[fNLayers];
  for(Int_t i=0; i<fNLayers; i++){ 
    dmar[i]=0;
  }

  if (fCluLayer == 0) {
    fCluLayer = new TClonesArray*[fNLayers];
    fCluCoord = new TClonesArray*[fNLayers];
    for(Int_t i=0;i<fNLayers;i++) {
      fCluLayer[i]=0;
      fCluCoord[i]=0;
    }
  }
  for(Int_t i=0;i<fNLayers;i++){
    if (!ForceSkippingOfLayer(i)) {
      for(Int_t cli=0;cli<fLayers[i]->GetNumberOfClusters();cli++){
        AliITSRecPointU* cls = (AliITSRecPointU*)fLayers[i]->GetCluster(cli);
        if(cls->TestBit(kSAflag)==kTRUE) continue; //clusters used by TPC prol.
        if(cls->GetQ()==0) continue; //fake clusters dead zones
        if(i>1 && cls->GetQ()<=fMinQ) continue; // cut on SDD and SSD cluster charge
        nclusters[i]++;
      }
    } 
    dmar[i]=0;
    if(!fCluLayer[i]){
      fCluLayer[i] = new TClonesArray("AliITSRecPointU",nclusters[i]);
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
  for(Int_t ilay=0;ilay<fNLayers;ilay++){
    TClonesArray &clulay = *fCluLayer[ilay];
    TClonesArray &clucoo = *fCluCoord[ilay];
    if (!ForceSkippingOfLayer(ilay)){
      AliDebug(2,Form("number of clusters in layer %i : %i",ilay,fLayers[ilay]->GetNumberOfClusters()));
      for(Int_t cli=0;cli<fLayers[ilay]->GetNumberOfClusters();cli++){
        AliITSRecPointU* cls = (AliITSRecPointU*)fLayers[ilay]->GetCluster(cli);
        if(cls->TestBit(kSAflag)==kTRUE) continue;
        if(cls->GetQ()==0) continue;
        Double_t phi=0;Double_t lambda=0;
        Double_t x=0;Double_t y=0;Double_t z=0;
        Float_t sx=0;Float_t sy=0;Float_t sz=0;
        GetCoorAngles(cls,phi,lambda,x,y,z,primaryVertex);
        GetCoorErrors(cls,sx,sy,sz);
        new (clulay[dmar[ilay]]) AliITSRecPoint(*cls);
        new (clucoo[dmar[ilay]]) AliITSclusterTable(x,y,z,sx,sy,sz,phi,lambda,cli);
        dmar[ilay]++;
      }
    } else AliDebug(2,Form("Force skipping layer %i",ilay));
  }


  // track counter
  Int_t ntrack=0;

  static Int_t nClusLay[fgMaxNLayer];//counter for clusters on each layer
  Int_t startLayForSeed=0;
  Int_t lastLayForSeed=fOuterStartLayer;
  Int_t nSeedSteps=lastLayForSeed-startLayForSeed;
  Int_t seedStep=1;
  if(fInwardFlag){
    startLayForSeed=fNLayers-1;
    lastLayForSeed=fInnerStartLayer;
    nSeedSteps=startLayForSeed-lastLayForSeed;
    seedStep=-1;
  }
  // loop on minimum number of points
  for(Int_t iMinNPoints=fNLayers; iMinNPoints>=fMinNPoints; iMinNPoints--) {
    // loop on starting layer for track finding
    for(Int_t iSeedLay=0; iSeedLay<=nSeedSteps; iSeedLay++) {
      Int_t theLay=startLayForSeed+iSeedLay*seedStep;
      if(ForceSkippingOfLayer(theLay)) continue;
      Int_t minNPoints=iMinNPoints-theLay;
      if(fInwardFlag) minNPoints=iMinNPoints-(fNLayers-1-theLay);
      for(Int_t i=theLay+1;i<fNLayers;i++)
        if(ForceSkippingOfLayer(i))
          minNPoints--;
      if(minNPoints<fMinNPoints) continue;

      // loop on phi and lambda window size
      for(Int_t nloop=0;nloop<fNloop;nloop++){
	Int_t nclTheLay=fCluLayer[theLay]->GetEntriesFast();
        while(nclTheLay--){
          ResetForFinding();
          Bool_t useRP=SetFirstPoint(theLay,nclTheLay,primaryVertex);
          if(!useRP) continue;
          AliITStrackU trs(fNLayers);

          Int_t pflag=0;
          Int_t kk;
          for(kk=0;kk<fNLayers;kk++) nClusLay[kk] = 0;
          kk=0;
          nClusLay[kk] = SearchClusters(theLay,fPhiWin[nloop],fLambdaWin[nloop],
                                        &trs,primaryVertex[2],pflag);
          Int_t nextLay=theLay+seedStep;
          if(nextLay<0){
          AliWarning("Negative layer id. Please check!");
          continue;
          }
          Bool_t goon=kTRUE;
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
            if(nextLay<0 || nextLay==fNLayers) goon=kFALSE;
          }


          Int_t layOK=0;
          if(!fInwardFlag){
            for(Int_t nnp=0;nnp<fNLayers-theLay;nnp++){
              if(nClusLay[nnp]!=0) layOK+=1;
            }
          }else{
            for(Int_t nnp=theLay; nnp>=0; nnp--){
              if(nClusLay[nnp]!=0) layOK+=1;
            }
          }
	  if(layOK>=minNPoints){
	    //    AliDebug(2,Form("---NPOINTS: %d; MAP: %d %d %d %d %d %d %d\n",layOK,nClusLay[0],nClusLay[1],nClusLay[2],nClusLay[3],nClusLay[4],nClusLay[5],nClusLay[6]));
            AliITStrackV2* tr2 = 0;
            Bool_t onePoint = kFALSE;
	    tr2 = FitTrack(&trs,primaryVertex,onePoint);
	    if(!tr2){
	      continue;
	    }
	    //	    AliDebug(2,Form("---NPOINTS fit: %d\n",tr2->GetNumberOfClusters()));

	    StoreTrack(tr2,event,useAllClusters);
	    ntrack++;

	  }

        }//end loop on clusters of theLay
      } //end loop on window sizes
    } //end loop on theLay
  }//end loop on min points
  if(!useAllClusters) AliInfo(Form("Number of found tracks: %d",event->GetNumberOfTracks()));

  ResetForFinding();
  return 0;

}
 
//________________________________________________________________________

AliITStrackV2* AliITStrackerUpgrade::FitTrack(AliITStrackU* tr,Double_t *primaryVertex,Bool_t onePoint) {

  const Int_t kMaxClu=AliITStrackU::kMaxNumberOfClusters;
  
  static Int_t clind[fgMaxNLayer][kMaxClu];
  static Int_t clmark[fgMaxNLayer][kMaxClu];
  static Int_t end[fgMaxNLayer];
  static AliITSRecPointU *listlayer[fgMaxNLayer][kMaxClu];
  static Int_t indices[fgMaxNLayer];

  for(Int_t k=0;k<fgMaxNLayer; k++)end[k]=0;

  for(Int_t i=0;i<fNLayers;i++) {
    for(Int_t j=0;j<kMaxClu; j++){
      clind[i][j]=0;
      clmark[i][j]=0;
      listlayer[i][j]=0;
    }
  }
  Int_t inx[fNLayers]; 
  for (Int_t k=0; k<fNLayers; k++) inx[k]=-1;
  Int_t nclusters = tr->GetNumberOfClustersU();
  for(Int_t ncl=0;ncl<nclusters;ncl++){
    Int_t index = tr->GetClusterIndexU(ncl); 
    Int_t lay = (index & 0xf0000000) >> 28;
    //Int_t cli = index&0x0fffffff;
    //AliITSRecPoint* cl = (AliITSRecPoint*)fLayers[lay]->GetCluster(cli);
    AliITSRecPointU* cl = (AliITSRecPointU*)GetCluster(index);
    if(cl->TestBit(kSAflag)==kTRUE) cl->ResetBit(kSAflag);  
    inx[lay]=index;
 
    Int_t nInLay=end[lay];
    listlayer[lay][nInLay]=cl;
    clind[lay][nInLay]=index;
    end[lay]++;
  }

  for(Int_t nlay=0;nlay<fNLayers;nlay++){
    for(Int_t ncl=0;ncl<tr->GetNumberOfMarked(nlay);ncl++){
      Int_t mark = tr->GetClusterMark(nlay,ncl);
      clmark[nlay][ncl]=mark;
    }
  }


  Int_t firstLay=-1,secondLay=-1;
  for(Int_t i=0;i<fNLayers;i++) {
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

  //---adding initialization for non existent layers
  for(Int_t i=fNLayers;i<fgMaxNLayer;i++) {
    if(end[i]==0) {
      end[i]=1;
    }
  }



  if(firstLay==-1 || (secondLay==-1 && !onePoint)) return 0;

  TClonesArray &arrSA= *fListOfUTracks;
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
              for(Int_t l6=0;l6<end[6];l6++){//loop on layer 7
		indices[6]=l6;
	        for(Int_t l7=0;l7<end[7];l7++){ //loop on layer 8
                  indices[7]=l7;

		  Int_t iFirstLay=indices[firstLay];
		  Int_t mrk1=clmark[firstLay][iFirstLay];

		  Double_t x1,y1,z1,sx1,sy1,sz1;
		  Double_t x2,y2,z2,sx2,sy2,sz2;
		  AliITSRecPointU* p1=(AliITSRecPointU*)listlayer[firstLay][iFirstLay];
		  Int_t layer = 0;
		  Double_t radius = 0.;
		  Double_t yclu1 = p1->GetY();
		  Double_t zclu1 = p1->GetZ();
		  layer=p1->GetLayer();
		  radius = fSegmentation->GetRadius(layer);
		  Double_t cv=0.,tgl2=0.,phi2=0.;
		  AliITSclusterTable* arr1 = (AliITSclusterTable*)GetClusterCoord(firstLay,mrk1);
		  x1 = arr1->GetX();
		  y1 = arr1->GetY();
		  z1 = arr1->GetZ();
		  sx1 = arr1->GetSx();
		  sy1 = arr1->GetSy();
		  sz1 = arr1->GetSz();

		  if(secondLay>0) {
		    Int_t iSecondLay=indices[secondLay];
		    Int_t mrk2=clmark[secondLay][iSecondLay];
		    AliITSclusterTable* arr2 = (AliITSclusterTable*)GetClusterCoord(secondLay,mrk2);
		    x2 = arr2->GetX();
		    y2 = arr2->GetY();
		    z2 = arr2->GetZ();
		    sx2 = arr2->GetSx();
		    sy2 = arr2->GetSy();
		    sz2 = arr2->GetSz();
		    cv = Curvature(primaryVertex[0],primaryVertex[1],x1,y1,x2,y2);
		    tgl2 = (z2-z1)/TMath::Sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		    phi2 = TMath::ATan2((y2-y1),(x2-x1));
		  }

		  Double_t x,y,z = 0.;
		  Double_t xz[2];
		  xz[0]= p1->GetDetLocalX();
		  xz[1]= p1->GetDetLocalZ();
		  fSegmentation->DetToGlobal(layer,p1->GetModule(),xz[0], xz[1],x,y,z);


		  Double_t phiclrad;
		  phiclrad=TMath::ATan2(y,x); 
		  if (phiclrad<0) phiclrad+=TMath::TwoPi();
		  else if (phiclrad>=TMath::TwoPi()) phiclrad-=TMath::TwoPi();

		  Double_t alpha = fSegmentation->GetAlpha(fSegmentation->GetModule(phiclrad));
   
                  Int_t clIndice[fNLayers]; 
		  for (Int_t k=0; k<fNLayers; k++) clIndice[k]=-1;

		  AliITStrackU trac(alpha,radius,yclu1,zclu1,phi2,tgl2,cv,1,fNLayers);
		  for(Int_t iLay=fNLayers-1; iLay>=0; iLay--){          
		    Int_t iInLay=indices[iLay];
		    AliITSRecPointU* cl=(AliITSRecPointU*)listlayer[iLay][iInLay];
		    if(cl!=0){
		      trac.AddClusterV2(iLay,(clind[iLay][iInLay] & 0x0fffffff)>>0);
		      trac.AddClusterMark(iLay,clmark[iLay][iInLay]);
		      clIndice[iLay] = clind[iLay][iInLay];
		    }
		  }

		  //fit with Kalman filter using RefitAtBase()
  		  AliITStrackU ot(trac, kFALSE);
		  ot.ResetCovariance(10.);
		  ot.ResetClusters();
              
		  // Propagate inside the innermost layer with a cluster 
		  if(ot.Propagate(ot.GetX()-0.1*ot.GetX())) {
		   if(RefitAtBase(AliITSRecoParam::GetrInsideITSscreen(),&ot,clIndice)){ //fit from layer 1 to layer N	
		      AliITStrackU otrack2(ot,kFALSE);
		      otrack2.ResetCovariance(10.); 
		      otrack2.ResetClusters();
		      //fit from last layer to layer 1
		      if(RefitAtBase(fSegmentation->GetRadius(0)-0.1,&otrack2,clIndice)) {
			new(arrSA[nFoundTracks]) AliITStrackU(otrack2);
			++nFoundTracks;
		      }
                              
		    }       
		  }

		}//end loop layer 8
              }//end loop layer 7
            }//end loop layer 6
          }//end loop layer 5
        }//end loop layer 4        
      }//end loop layer 3
    }//end loop layer 2 
  }//end loop layer 1




  if(fListOfUTracks->GetEntries()==0) return 0;

  Int_t lowchi2 = FindTrackLowChiSquare();
  AliITStrackU* otrack = (AliITStrackU*)fListOfUTracks->At(lowchi2);
   
  if(otrack==0) {
    return 0;
  }

  //  CookLabel(otrack,0.); //MI change - to see fake ratio
  //  printf("  label from CookLabel: %d  \n   ",otrack->GetLabel());

  Int_t label=FindLabel(otrack);
  otrack->SetLabel(label);
  //  printf("  label from FindLabel: %d  \n   ",otrack->GetLabel());

  Int_t indexc[fNLayers];
  for(Int_t i=0;i<fNLayers;i++) indexc[i]=0;
  for(Int_t nind=0;nind<otrack->GetNumberOfClusters();nind++){
    indexc[nind] = otrack->GetClusterIndex(nind);
  }      

  
  //remove clusters of found track
  for(Int_t nlay=0;nlay<fNLayers;nlay++){
    for(Int_t cln=0;cln<otrack->GetNumberOfMarked(nlay);cln++){
      Int_t index = otrack->GetClusterMark(nlay,cln);
      fCluLayer[nlay]->RemoveAt(index);
      RemoveClusterCoord(nlay,index);
      fCluLayer[nlay]->Compress();
    }    
  }
  
  return otrack;

}

//_______________________________________________________
void AliITStrackerUpgrade::StoreTrack(AliITStrackV2 *t,AliESDEvent *event, Bool_t pureSA) /*const*/ 
{
  //
  // Add new track to the ESD
  //
  AliESDtrack outtrack;
  outtrack.UpdateTrackParams(t,AliESDtrack::kITSin);
  if(pureSA) outtrack.SetStatus(AliESDtrack::kITSpureSA);
  outtrack.SetStatus(AliESDtrack::kITSout);
  outtrack.SetStatus(AliESDtrack::kITSrefit);
  Double_t sdedx[4]={0.,0.,0.,0.};
  for(Int_t i=0; i<4; i++) sdedx[i]=t->GetSampledEdx(i);
  outtrack.SetITSdEdxSamples(sdedx);
  event->AddTrack(&outtrack);
  return;
}


//_______________________________________________________
Int_t AliITStrackerUpgrade::SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, AliITStrackU* trs,Double_t /*zvertex*/,Int_t pflag){
  //function used to to find the clusters associated to the track

  AliDebug(2,"Starting...");
  if(ForceSkippingOfLayer(layer)) {
    AliDebug(2,Form("Forcing skipping of layer %i. Exiting",layer));
    return 0;
  }

  Int_t nc=0;
  Double_t r=fSegmentation->GetRadius(layer);
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

    AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(layer,index);

    Double_t lambda = arr->GetLambda();
    if (TMath::Abs(lambda-lamExpect)>lambdawindow) continue;

    Double_t phi = arr->GetPhi();
    Double_t deltaPhi = phi-phiExpect;
    if(deltaPhi>TMath::Pi()) deltaPhi-=2*TMath::Pi();
    else if(deltaPhi<-TMath::Pi()) deltaPhi+=2*TMath::Pi();
    if (TMath::Abs(deltaPhi)>phiwindow) continue;

    if(trs->GetNumberOfClustersU()==trs->GetMaxNumberOfClusters()) return 0;
    if(trs->GetNumberOfMarked(layer)==trs->GetMaxNMarkedPerLayer()) return 0;
    Int_t orind = arr->GetOrInd();
    trs->AddClusterU(layer,orind);
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
Bool_t AliITStrackerUpgrade::SetFirstPoint(Int_t lay, Int_t clu, Double_t* primaryVertex){
  // Sets the first point (seed) for tracking

  AliITSRecPointU* cl = (AliITSRecPointU*)fCluLayer[lay]->At(clu);
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
void AliITStrackerUpgrade::UpdatePoints(){
  //update of points for the estimation of the curvature  

  fPoint2[0]=fPoint3[0];
  fPoint2[1]=fPoint3[1];
  fPoint3[0]=fPointc[0];
  fPoint3[1]=fPointc[1];
}
//___________________________________________________________________
Int_t AliITStrackerUpgrade::FindEquation(Float_t x1, Float_t y1, Float_t x2, Float_t y2, Float_t x3, Float_t y3,Float_t& a, Float_t& b, Float_t& c){

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
Int_t AliITStrackerUpgrade::FindIntersection(Float_t a1, Float_t b1, Float_t c1, Float_t c2,Float_t& x1,Float_t& y1, Float_t& x2, Float_t& y2){
 
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
Double_t AliITStrackerUpgrade::Curvature(Double_t x1,Double_t y1,Double_t 
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
Double_t AliITStrackerUpgrade::ChoosePoint(Double_t p1, Double_t p2, Double_t pp){

  //Returns the point closest to pp

  Double_t diff1 = p1-pp;
  Double_t diff2 = p2-pp;
  
  if(TMath::Abs(diff1)<TMath::Abs(diff2)) fPhiEstimate=p1;
  else fPhiEstimate=p2;  
  return fPhiEstimate;
  
}
//_________________________________________________________________
Int_t AliITStrackerUpgrade::FindTrackLowChiSquare() const {
  // returns track with lowest chi square  
  Int_t dim=fListOfUTracks->GetEntries();
  if(dim<=1) return 0;
  AliITStrackV2* trk = (AliITStrackV2*)fListOfUTracks->At(0);
  Double_t minChi2=trk->GetChi2();
  Int_t index=0;
  for(Int_t i=1;i<dim;i++){
    trk = (AliITStrackV2*)fListOfUTracks->At(i);
    Double_t chi2=trk->GetChi2();
    if(chi2<minChi2){
      minChi2=chi2;
      index=i;
    }
  }
  return index;
}

//__________________________________________________________
Int_t AliITStrackerUpgrade::FindLabel(AliITStrackV2* track) const {
  //
  // max number of labels in AliITSRecPointU is 24;
  enum {kMaxLab=24};

  Int_t labl[fNLayers][kMaxLab];
  Int_t cnts[fNLayers][kMaxLab];
  for(Int_t j=0;j<fNLayers;j++){
    for(Int_t k=0;k<kMaxLab;k++){
      labl[j][k]=-2;
      cnts[j][k]=1;
    }
  }
  Int_t iNotLabel=0;
  for(Int_t i=0;i<track->GetNumberOfClusters(); i++) {
    Int_t indexc = track->GetClusterIndex(i);
    AliITSRecPointU* cl = (AliITSRecPointU*)GetCluster(indexc);
    //  AliDebug(2,Form("Layer %d; cluster index %i; MClabels=(%d,%d,%d) ",
    //		    cl->GetLayer(),indexc,cl->GetLabel(0),cl->GetLabel(1),cl->GetLabel(2)));
    Int_t iLayer=cl->GetLayer();
    for(Int_t k=0;k<kMaxLab;k++){
      labl[iLayer][k]=cl->GetTrackID(k);
      if(labl[iLayer][k]<0) iNotLabel++;
    }
  }
  if(iNotLabel==kMaxLab*track->GetNumberOfClusters()) return -2; // No MC label found at all

  for(Int_t j1=0;j1<fNLayers; j1++) {
    for(Int_t j2=0; j2<j1;  j2++){
      for(Int_t k1=0; k1<kMaxLab; k1++){
        for(Int_t k2=0; k2<kMaxLab; k2++){
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
  for(Int_t j=0;j<fNLayers;j++){
    for(Int_t k=0;k<kMaxLab;k++){
      if(cnts[j][k]>cntMax && labl[j][k]>=0){
        cntMax=cnts[j][k];
        label=labl[j][k];
      }
    }
  }

  Int_t lflag=0;
  for(Int_t i=0;i<fNLayers;i++) {
    Bool_t labelOkInLayer = kFALSE;
    for(Int_t k=0;k<kMaxLab;k++){
      if (labl[i][k]==label) {
	labelOkInLayer=kTRUE;
	break;
      }
    }
    if(labelOkInLayer) lflag++;
  }
  if(lflag<track->GetNumberOfClusters()) {
    if (label==0) {
      label = -1; // to avoid label 0 which can not be identified as as fake with -0
    } else {
      label = -label;
    }
  }
  if (label<0) {
    for(Int_t i=0;i<track->GetNumberOfClusters(); i++) {
      Int_t indexc = track->GetClusterIndex(i);
      AliITSRecPointU* cl = (AliITSRecPointU*)GetCluster(indexc);
      AliDebug(2,Form("Layer %d; cluster index %i; MClabels=(%d,%d,%d...) ",
		      cl->GetLayer(),indexc,cl->GetTrackID(0),cl->GetTrackID(1),cl->GetTrackID(2)));
    }
  }


  return label;
}

//_____________________________________________________________________________
void AliITStrackerUpgrade::SetCalculatedWindowSizes(Int_t n, Float_t phimin, Float_t phimax, Float_t lambdamin, Float_t lambdamax){
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
void AliITStrackerUpgrade::SetFixedWindowSizes(Int_t n, Double_t *phi, Double_t *lam){
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

    Double_t factor=AliITSReconstructor::GetRecoParam()->GetFactorSAWindowSizes(); 
  
    for(Int_t k=0;k<fNloop;k++){
      fPhiWin[k]=phid[k]*factor;
      fLambdaWin[k]=lambdad[k]*factor;
    }
  
  }

}
//_______________________________________________________________________
void AliITStrackerUpgrade::GetCoorAngles(AliITSRecPointU* cl,Double_t &phi,Double_t &lambda, Double_t &x, Double_t &y,Double_t &z,Double_t* vertex){
  //Returns values of phi (azimuthal) and lambda angles for a given cluster
  
  Double_t xz[2];
  xz[0]= cl->GetDetLocalX(); 
  xz[1]= cl->GetDetLocalZ() ; 
  Int_t ilayer;
  ilayer = cl->GetLayer();
  fSegmentation->DetToGlobal(ilayer,cl->GetModule(),xz[0], xz[1],x,y,z);

  if(x!=0 && y!=0)  
    phi=TMath::ATan2(y-vertex[1],x-vertex[0]);
  lambda=TMath::ATan2(z-vertex[2],TMath::Sqrt((x-vertex[0])*(x-vertex[0])+(y-vertex[1])*(y-vertex[1])));
}

//________________________________________________________________________
void AliITStrackerUpgrade::GetCoorErrors(AliITSRecPointU* cl,Float_t &sx,Float_t &sy, Float_t &sz){

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



  //  Float_t covm[6];
  //  cl->GetGlobalCov(covm);
  //  sx=12*1e-04;//TMath::Sqrt(covm[0]);
  //  sy=12*1e-04;//TMath::Sqrt(covm[3]);
  //  sz=40*1e-04;//TMath::Sqrt(covm[5]);
  /*
    sx = TMath::Sqrt(sp*sp*cl->GetSigmaY2());
    sy = TMath::Sqrt(cp*cp*cl->GetSigmaY2());
    sz = TMath::Sqrt(cl->GetSigmaZ2());
  */
}
//_____________________________________________________________
void AliITStrackerUpgrade::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads ITS clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fNLayers; i++) fLayers[i]->ResetClusters();
}
//______________________________________________________________________
Bool_t AliITStrackerUpgrade::RefitAtBase(Double_t xx,AliITStrackU *track,
					 const Int_t *clusters)
{
  //--------------------------------------------------------------------
  // Simplified version for ITS upgrade studies -- does not use module info
  //--------------------------------------------------------------------
  Int_t index[fNLayers];
  Int_t k;
  for (k=0; k<fNLayers; k++) index[k]=-1;
  //
  for (k=0; k<fNLayers; k++) {
    index[k]=clusters[k];
  }

  ULong_t trStatus=0;
  if(track->GetESDtrack()) trStatus=track->GetStatus();
  Int_t innermostlayer=0;
  if(trStatus&AliESDtrack::kTPCin)  {
    innermostlayer=fNLayers-1;
    Double_t drphi = TMath::Abs(track->GetD(0.,0.));
    for(innermostlayer=0; innermostlayer<fNLayers; innermostlayer++) {
      if( (drphi < (fSegmentation->GetRadius(innermostlayer)+1.)) ||
          index[innermostlayer] >= 0 ) break;
    }
    AliDebug(2,Form(" drphi  %f  innermost %d",drphi,innermostlayer));
  }

  Int_t from, to, step;
  if (xx > track->GetX()) {
    from=innermostlayer; to=fNLayers;
    step=+1;
  } else {
    from=fNLayers-1; to=innermostlayer-1;
    step=-1;
  }
  TString dir = (step>0 ? "outward" : "inward");

  for (Int_t ilayer = from; ilayer != to; ilayer += step) {
    Double_t r=0.;
    r=fSegmentation->GetRadius(ilayer);

    if (step<0 && xx>r){
      break;
    }


    Double_t oldGlobXYZ[3];

    if (!track->GetXYZ(oldGlobXYZ)) {
      return kFALSE;
    }
    // continue if we are already beyond this layer
    Double_t oldGlobR = TMath::Sqrt(oldGlobXYZ[0]*oldGlobXYZ[0]+oldGlobXYZ[1]*oldGlobXYZ[1]);
    if(step>0 && oldGlobR > r){
      continue; // going outward
    }
    if(step<0 && oldGlobR < r) {
      continue; // going inward
    }
    Double_t phi,z;
    if (!track->GetPhiZat(r,phi,z)){
      return kFALSE;
    } // only for ITS-SA tracks refit
    track->SetCheckInvariant(kFALSE);
      
    if (phi<0) phi+=TMath::TwoPi();//from 0 to 360 (rad) 
    else if (phi>=TMath::TwoPi()) phi-=TMath::TwoPi();//

    Int_t trkModule = fSegmentation->GetModule(phi);

    Double_t alpha = fSegmentation->GetAlpha(trkModule);
    if (!track->Propagate(alpha,r)) {
      return kFALSE;
    }

    const AliITSRecPointU *clAcc=0;
    Double_t maxchi2=1000.*AliITSReconstructor::GetRecoParam()->GetMaxChi2();

    Int_t idx=index[ilayer];
    if (idx>=0) { // cluster in this layer   
      Int_t cli = idx&0x0fffffff;
      AliITSRecPointU *cl=(AliITSRecPointU *)fLayers[ilayer]->GetCluster(cli);

      if (cl) {                                                                  
	Int_t cllayer = (idx & 0xf0000000) >> 28;;                               
       // edge effect correction 
        if(trkModule!=cl->GetModule()){
	if (!track->Propagate(fSegmentation->GetAlpha(cl->GetModule()),fSegmentation->GetRadius(ilayer))) return kFALSE;
        }
	Double_t chi2=GetPredictedChi2MI(track,cl,cllayer);                      
	if (chi2<maxchi2) {                                                      
	  clAcc=cl;                                                              
	  maxchi2=chi2;                                                          
	} else {                                                                 
	  return kFALSE;                                                         
	}                                                                        
      }                                                                          
    }              

    if (clAcc) {
      if (!UpdateMI(track,clAcc,maxchi2,idx)){
	return kFALSE;
      }
      track->SetSampledEdx(clAcc->GetQ(),ilayer-2);
    }


    // Correct for material of the current layer
    // cross material
    // add time if going outward
    if(!CorrectForLayerMaterial(track,ilayer,oldGlobXYZ,dir)){
      return kFALSE;
    }
    track->SetCheckInvariant(kFALSE);
  } // end loop on layers

  if (!track->PropagateTo(xx,0.,0.)){
    return kFALSE;
  } 
 
  return kTRUE;
}

//_____________________________________________________________________
Int_t AliITStrackerUpgrade::UpdateMI(AliITStrackU* track, const AliITSRecPointU* cl,Double_t chi2,Int_t index) const
{
  //
  // Update ITS track
  //
  Int_t layer = (index & 0xf0000000) >> 28;
  track->SetClIndex(layer, index);
  if (TMath::Abs(cl->GetQ())<1.e-13) return 0;  // ingore the "virtual" clusters


  // Take into account the mis-alignment (bring track to cluster plane)
  Double_t xTrOrig=track->GetX();
  Double_t trxyz[3]; track->GetXYZ(trxyz);
  Float_t xclu1 = cl->GetX();//upgrade clusters tracking coordinate
  //Float_t yclu1 = cl->GetY();
  //Float_t zclu1 = cl->GetZ();
  if (!track->Propagate(xTrOrig+xclu1)){
    return 0;
  }
  AliCluster c(*cl);
  c.SetSigmaY2(track->GetSigmaY(layer)*track->GetSigmaY(layer));
  c.SetSigmaZ2(track->GetSigmaZ(layer)*track->GetSigmaZ(layer));
  c.SetSigmaYZ(track->GetSigmaYZ(layer));

  // Bring the track back to detector plane in ideal geometry  
  Int_t updated = track->UpdateMI(&c,chi2,index);
  if (!track->Propagate(xTrOrig)) {
    return 0;
  }
  if(!updated){
    AliDebug(2,"update failed");
  }
  return updated;
}
//____________________________________________________________________-
Int_t AliITStrackerUpgrade::CorrectForLayerMaterial(AliITStrackU *t,
						    Int_t layerindex,
						    Double_t oldGlobXYZ[3],
						    TString direction) {
  //-------------------------------------------------------------------
  // Propagate beyond layer and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  // Define budget mode:
  // 0: material from AliITSRecoParam (hard coded)
  // 1: material from TGeo in stepsof X cm (on the fly)
  //    X = AliITSRecoParam::GetStepSizeTGeo()
  // 2: material from lut
  // 3: material from TGeo in one step (same for all hypotheses)
  Int_t mode;
  switch(fUseTGeo) {
  case 0:
    mode=0;
    break;
  case 1:
    mode=1;
    break;
  case 2:
    mode=2;
    break;
  case 3:
    if(fTrackingPhase.Contains("Clusters2Tracks"))
      { mode=3; } else { mode=1; }
    break;
  case 4:
    if(fTrackingPhase.Contains("Clusters2Tracks"))
      { mode=3; } else { mode=2; }
    break;
  default:
    mode=0;
    break;
  }
  if(fTrackingPhase.Contains("Default")) {
    mode=0;
  }
  Float_t  dir = (direction.Contains("inward") ? 1. : -1.);
  //Double_t r = fSegmentation->GetRadius(layerindex);
  //Double_t deltar=(layerindex<2 ? 0.10*r : 0.05*r);
  Double_t deltar=0.05;
  Double_t rToGo=TMath::Sqrt(t->GetX()*t->GetX()+t->GetY()*t->GetY())-deltar*dir;
  Double_t xToGo;
  if (!t->GetLocalXat(rToGo,xToGo)) {
    return 0;
  }
  Int_t index=fNLayers*fCurrentEsdTrack+layerindex;

  Double_t xOverX0=0.0,x0=0.0,lengthTimesMeanDensity=0.0;
  Int_t nsteps=1;
  // back before material (no correction)
  Double_t rOld,xOld;
  rOld=TMath::Sqrt(oldGlobXYZ[0]*oldGlobXYZ[0]+oldGlobXYZ[1]*oldGlobXYZ[1]);
  if (!t->GetLocalXat(rOld,xOld)) return 0;
  if (!t->Propagate(xOld)) return 0;

  switch(mode) {
  case 1:
    x0=21.82;
    xOverX0 = fSegmentation->GetThickness(layerindex)/x0;

    lengthTimesMeanDensity = xOverX0*x0;
    lengthTimesMeanDensity *= dir;
    // Bring the track beyond the material
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)){
      return 0;
    }
    break;
  case 0:
    nsteps = (Int_t)(TMath::Abs(xOld-xToGo)/AliITSReconstructor::GetRecoParam()->GetStepSizeTGeo())+1;
    if (!t->PropagateToTGeo(xToGo,nsteps)) return 0; // cross the material and apply correction
    break;
  case 2:
    if(fxOverX0Layer[layerindex]<0) BuildMaterialLUT("Layers");
    xOverX0 = fxOverX0Layer[layerindex];
    lengthTimesMeanDensity = fxTimesRhoLayer[layerindex];
    lengthTimesMeanDensity *= dir;
    // Bring the track beyond the material
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 3:
    if(!fxOverX0LayerTrks || index<0 || index>=fNLayers*fNtracks) Error("CorrectForLayerMaterial","Incorrect usage of UseTGeo option!\n");
    if(fxOverX0LayerTrks[index]<0) {
      nsteps = (Int_t)(TMath::Abs(xOld-xToGo)/AliITSReconstructor::GetRecoParam()->GetStepSizeTGeo())+1;
      if (!t->PropagateToTGeo(xToGo,nsteps,xOverX0,lengthTimesMeanDensity)) return 0;
      Double_t angle=TMath::Sqrt((1.+t->GetTgl()*t->GetTgl())/
                                 ((1.-t->GetSnp())*(1.+t->GetSnp())));
      fxOverX0LayerTrks[index] = TMath::Abs(xOverX0)/angle;
      fxTimesRhoLayerTrks[index] = TMath::Abs(lengthTimesMeanDensity)/angle;
      return 1;
    }
    xOverX0 = fxOverX0LayerTrks[index];
    lengthTimesMeanDensity = fxTimesRhoLayerTrks[index];
    lengthTimesMeanDensity *= dir;
    // Bring the track beyond the material
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)){
      return 0;
    }
    break;
  }


  return 1;
}


//_____________________________________________________________________________
Double_t AliITStrackerUpgrade::GetPredictedChi2MI(AliITStrackU* track, const AliITSRecPointU *cluster,Int_t layer)
{
  //
  // Compute predicted chi2
  //
  Float_t erry,errz,covyz;
  Float_t theta = track->GetTgl();
  Float_t phi   = track->GetSnp();
  phi = TMath::Abs(phi)*TMath::Sqrt(1./((1.-phi)*(1.+phi)));
  GetError(layer,cluster,theta,phi,track->GetExpQ(),erry,errz,covyz);
  AliDebug(3,Form(" chi2: tr-cl   %f  %f   tr X %f cl X %f",track->GetY()-cluster->GetY(),track->GetZ()-cluster->GetZ(),track->GetX(),cluster->GetX()));
  // Take into account the mis-alignment (bring track to cluster plane)
  Double_t xTrOrig=track->GetX();
  if (!track->Propagate(xTrOrig+cluster->GetX())) return 1000.;
  AliDebug(3,Form(" chi2: tr-cl   %f  %f   tr X %f cl X %f",track->GetY()-cluster->GetY(),track->GetZ()-cluster->GetZ(),track->GetX(),cluster->GetX()));
  Double_t chi2 = track->GetPredictedChi2MI(cluster->GetY(),cluster->GetZ(),erry,errz,covyz);
  // Bring the track back to detector plane in ideal geometry
  // [mis-alignment will be accounted for in UpdateMI()]
  if (!track->Propagate(xTrOrig)) return 1000.;
  Float_t ny=0.,nz=0.;
  //GetNTeor(layer,cluster,theta,phi,ny,nz);
  Double_t delta = cluster->GetNy()+cluster->GetNz()-nz-ny;
  if (delta>1){
    chi2+=0.5*TMath::Min(delta/2,2.);
    chi2+=2.*cluster->GetDeltaProbability();
  }
  //
  track->SetNy(layer,ny);
  track->SetNz(layer,nz);
  track->SetSigmaY(layer,erry);
  track->SetSigmaZ(layer, errz);
  track->SetSigmaYZ(layer,covyz);

  track->SetNormQ(layer,cluster->GetQ()/TMath::Sqrt((1.+ track->GetTgl()*track->GetTgl())/((1.-track->GetSnp())*(1.+track->GetSnp()))));
  return chi2;
}
//________________________________________________________________
Int_t AliITStrackerUpgrade::GetError(Int_t /*layer*/,
				     const AliITSRecPointU *cl,
				     Float_t /*tgl*/,Float_t /*tgphitr*/,Float_t /*expQ*/,
				     Float_t &erry,Float_t &errz,Float_t &covyz,
				     Bool_t /*addMisalErr*/)
{
  //
  // Calculate cluster position error for the upgrade
  //
  Int_t retval=0;
  covyz=0.;
  retval = GetErrorOrigRecPoint(cl,erry,errz,covyz);



  return retval;
                                                                                   
}
//____________________________________________________________________-

Int_t AliITStrackerUpgrade::GetErrorOrigRecPoint(const AliITSRecPointU *cl,
						 Float_t &erry,Float_t &errz,Float_t &covyz)
{
  //
  // Calculate cluster position error (just take error from AliITSRecPoint)
  //
  erry   = TMath::Sqrt(cl->GetSigmaY2());
  errz   = TMath::Sqrt(cl->GetSigmaZ2());
  covyz  = cl->GetSigmaYZ();
  return 1;
}
//__________________________
//void AliITStrackerUpgrade::GetNTeor(Int_t layer,const AliITSRecPointU* /*cl*/,
/*				    Float_t tgl,Float_t tgphitr,
				    Float_t &ny,Float_t &nz)
{
  //
  // Get "mean shape" (original parametrization from AliITStrackerMI)
  //
  tgl = TMath::Abs(tgl);
  tgphitr = TMath::Abs(tgphitr);

  // SPD
  if (layer==0) {
    ny = 1.+tgphitr*3.2;
    nz = 1.+tgl*0.34;
    return;
  }
  if (layer==1) {
    ny = 1.+tgphitr*3.2;
    nz = 1.+tgl*0.28;
    return;
  }
  // SSD
  if (layer==4 || layer==5) {
    ny = 2.02+tgphitr*1.95;
    nz = 2.02+tgphitr*2.35;
    return;
  }
  // SDD
  ny  = 6.6-2.7*tgphitr;
  nz  = 2.8-3.11*tgphitr+0.45*tgl;
  return;
  }*/

//_________________________________________________________________
Int_t AliITStrackerUpgrade::PropagateBack(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This functions propagates reconstructed ITS tracks back
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  fTrackingPhase="PropagateBack";
  Int_t nentr=event->GetNumberOfTracks();
  Info("PropagateBack", "Number of ESD tracks: %d\n", nentr);
  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);

    if ((esd->GetStatus()&AliESDtrack::kITSin)==0) continue;
    if (esd->GetStatus()&AliESDtrack::kITSout) continue;

    /*
      AliITStrackMI *t=0;
      try {
      t=new AliITStrackMI(*esd);
      } catch (const Char_t *msg) {
      //Warning("PropagateBack",msg);
      delete t;
      continue;
      }
    */
    AliITStrackU *t = new AliITStrackU(*esd);
    t->SetExpQ(TMath::Max(0.8*t->GetESDtrack()->GetTPCsignal(),30.));

    ResetTrackToFollow(*t);

    /*
    // propagate to vertex [SR, GSI 17.02.2003]
    // Start Time measurement [SR, GSI 17.02.2003], corrected by I.Belikov
    if (CorrectForPipeMaterial(&fTrackToFollow,"inward")) {
    if (fTrackToFollow.PropagateToVertex(event->GetVertex()))
    fTrackToFollow.StartTimeIntegral();
    // from vertex to outside pipe
    CorrectForPipeMaterial(&fTrackToFollow,"outward");
    }*/
    // Start time integral and add distance from current position to vertex
    Double_t xyzTrk[3],xyzVtx[3]={GetX(),GetY(),GetZ()};
    fTrackToFollow.GetXYZ(xyzTrk);
    Double_t dst2 = 0.;
    for (Int_t icoord=0; icoord<3; icoord++) {
      Double_t di = xyzTrk[icoord] - xyzVtx[icoord];
      dst2 += di*di;
    }
    fTrackToFollow.StartTimeIntegral();
    fTrackToFollow.AddTimeStep(TMath::Sqrt(dst2));

    fTrackToFollow.ResetCovariance(10.); 
    fTrackToFollow.ResetClusters();
    //
    Int_t inx[fNLayers];
    for (Int_t k=0; k<fNLayers; k++) inx[k]=-1;
    Int_t nclusters = t->GetNumberOfClusters();
    for(Int_t ncl=0;ncl<nclusters;ncl++){
      Int_t index = t-> GetClIndex(ncl);
      Int_t lay = (index & 0xf0000000) >> 28;
      inx[lay]=index;
    }
    //   
    if (RefitAtBase(fSegmentation->GetRadius(fNLayers-1),&fTrackToFollow,inx)) {
      //fTrackToFollow.SetLabel(t->GetLabel());//
      //fTrackToFollow.CookdEdx();
      //CookLabel(&fTrackToFollow,0.); //For comparison only
      fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSout);
      //UseClusters(&fTrackToFollow);
      ntrk++;
    }
    delete t;
  }

  fTrackingPhase="Default";

  return 0;
}
//_________________________________________________________________________
AliCluster *AliITStrackerUpgrade::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t l=(index & 0xf0000000) >> 28;
  Int_t c=(index & 0x0fffffff) >> 0;

  AliCluster *cl = fLayers[l]->GetCluster(c);

  //  AliDebug(2,Form("index %i; cluster index %i; layer %i; Pointer to cluster:%p", index,c,l,cl));

  return cl;
}
//______________________________________________________________________________
Int_t AliITStrackerUpgrade::CorrectForPipeMaterial(AliITStrackU *t, TString direction) {
  //-------------------------------------------------------------------
  // Propagate beyond beam pipe and correct for material
  // (material budget in different ways according to fUseTGeo value)
  // Add time if going outward (PropagateTo or PropagateToTGeo)
  //-------------------------------------------------------------------

  // Define budget mode:
  // 0: material from AliITSRecoParam (hard coded)
  // 1: material from TGeo in one step (on the fly)
  // 2: material from lut
  // 3: material from TGeo in one step (same for all hypotheses)
  Int_t mode;
  switch(fUseTGeo) {
  case 0:
    mode=0;
    break;
  case 1:
    mode=1;
    break;
  case 2:
    mode=2;
    break;
  case 3:
    if(fTrackingPhase.Contains("Clusters2Tracks"))
      { mode=3; } else { mode=1; }
    break;
  case 4:
    if(fTrackingPhase.Contains("Clusters2Tracks"))
      { mode=3; } else { mode=2; }
    break;
  default:
    mode=0;
    break;
  }
  if(fTrackingPhase.Contains("Default")) mode=0;
  Int_t index=fCurrentEsdTrack;

  Float_t  dir = (direction.Contains("inward") ? 1. : -1.);
  Double_t rToGo=(dir>0 ? AliITSRecoParam::GetrInsidePipe() : AliITSRecoParam::GetrOutsidePipe());
  Double_t xToGo;
  if (!t->GetLocalXat(rToGo,xToGo)) return 0;

  Double_t xOverX0,x0,lengthTimesMeanDensity;

  switch(mode) {
  case 0:
    xOverX0 = AliITSRecoParam::GetdPipe();
    x0 = AliITSRecoParam::GetX0Be();
    lengthTimesMeanDensity = xOverX0*x0;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 1:
    if (!t->PropagateToTGeo(xToGo,1)) return 0;    
    break;
  case 2:
    if(fxOverX0Pipe<0) BuildMaterialLUT("Pipe");
    xOverX0 = fxOverX0Pipe;
    lengthTimesMeanDensity = fxTimesRhoPipe;
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  case 3:
    if(!fxOverX0PipeTrks || index<0 || index>=fNtracks) Error("CorrectForPipeMaterial","Incorrect usage of UseTGeo option!\n");
    if(fxOverX0PipeTrks[index]<0) {
      if (!t->PropagateToTGeo(xToGo,1,xOverX0,lengthTimesMeanDensity)) return 0;
      Double_t angle=TMath::Sqrt((1.+t->GetTgl()*t->GetTgl())/
				 ((1.-t->GetSnp())*(1.+t->GetSnp())));
      fxOverX0PipeTrks[index] = TMath::Abs(xOverX0)/angle;
      fxTimesRhoPipeTrks[index] = TMath::Abs(lengthTimesMeanDensity)/angle;
      return 1;
    }
    xOverX0 = fxOverX0PipeTrks[index];
    lengthTimesMeanDensity = fxTimesRhoPipeTrks[index];
    lengthTimesMeanDensity *= dir;
    if (!t->PropagateTo(xToGo,xOverX0,lengthTimesMeanDensity/xOverX0)) return 0;
    break;
  }

  return 1;
}
//__________________________________________________________________________
Int_t AliITStrackerUpgrade::RefitInward(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This functions refits ITS tracks using the
  // "inward propagated" TPC tracks
  // The clusters must be loaded !
  //--------------------------------------------------------------------
  fTrackingPhase="RefitInward";
  Int_t nentr=event->GetNumberOfTracks();
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);

    // if ((esd->GetStatus()&AliESDtrack::kITSout) == 0) continue;
    // if (esd->GetStatus()&AliESDtrack::kITSrefit) continue;
    //if (esd->GetStatus()&AliESDtrack::kTPCout)
    //  if ((esd->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;
    AliITStrackU *t=0;
    try {
      t=new AliITStrackU(*esd);
    } catch (const Char_t *msg) {
      //Warning("RefitInward",msg);
      delete t;
      continue;
    }

    ResetTrackToFollow(*t);
    // fTrackToFollow.ResetClusters();
    //  if ((esd->GetStatus()&AliESDtrack::kTPCin)==0)
    //  fTrackToFollow.ResetCovariance(10.);
    //Refitting...
    //The beam pipe
    if (CorrectForPipeMaterial(&fTrackToFollow,"inward")) {
      fTrackToFollow.UpdateESDtrack(AliESDtrack::kITSrefit);
      AliESDtrack  *esdTrack =fTrackToFollow.GetESDtrack();
      Double_t r[3]={0.,0.,0.};
      Double_t maxD=3.;
      esdTrack->RelateToVertex(event->GetVertex(),GetBz(r),maxD);
    }
    delete t;
  }

  fTrackingPhase="Default";
  return 0;
}



