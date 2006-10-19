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
//  Stand alone tracker class                     //
//  Origin:  Elisabetta Crescio                   //
//  e-mail:  crescio@to.infn.it                   //
//  tracks are saved as AliITStrackV2 objects     //
////////////////////////////////////////////////////

#include <stdlib.h>

#include <TArrayI.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliITSVertexer.h"
#include "AliITSclusterTable.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliITStrackSA.h"
#include "AliITStrackerSA.h"
#include "AliRun.h"

ClassImp(AliITStrackerSA)

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA():AliITStrackerMI(),
fPhiEstimate(0),
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
fGeom(0),
fListOfTracks(0),
fITSclusters(0),
fSixPoints(0),
fCluLayer(0),
fCluCoord(0){
  // Default constructor
  Init();
 
}
//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom):AliITStrackerMI(geom),
fPhiEstimate(0),
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
fGeom(0),
fListOfTracks(0),
fITSclusters(0),
fSixPoints(0),
fCluLayer(0),
fCluCoord(0) 
{
  // Standard constructor (Vertex is known and passed to this obj.)
  Init();
  fVert = 0;
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom, AliESDVertex *vert):AliITStrackerMI(geom),
fPhiEstimate(0),
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
fGeom(0),
fListOfTracks(0),
fITSclusters(0),
fSixPoints(0),
fCluLayer(0),
fCluCoord(0) 
{
  // Standard constructor (Vertex is known and passed to this obj.)
  Init();
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom, AliITSVertexer *vertexer):AliITStrackerMI(geom),
fPhiEstimate(0),
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
fGeom(0),
fListOfTracks(0),
fITSclusters(0),
fSixPoints(0),
fCluLayer(0),
fCluCoord(0)  
{
  // Standard constructor (Vertex is unknown - vertexer is passed to this obj)
  Init();
  fVertexer = vertexer;
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(const AliITStrackerSA& tracker):AliITStrackerMI(),
fPhiEstimate(tracker.fPhiEstimate),
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
fGeom(tracker.fGeom),
fListOfTracks(tracker.fListOfTracks),
fITSclusters(tracker.fITSclusters),
fSixPoints(tracker.fSixPoints),
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
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
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
  if(fCluLayer){
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      if(fCluLayer[i]){
	fCluLayer[i]->Delete();
	delete fCluLayer[i];
      }
    }
    delete [] fCluLayer;
  }
  if(fCluCoord){
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      if(fCluCoord[i]){
	fCluCoord[i]->Delete();
	delete fCluCoord[i];
      }
    }
    delete [] fCluCoord;
  }
  
  fGeom = 0; 
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
    fGeom = 0;
    SetWindowSizes();
    fITSclusters = 0;
    SetSixPoints();
    fListOfTracks=new TObjArray(0,0);
    fCluLayer = 0;
    fCluCoord = 0;
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
    fListOfTracks->Delete();
}

 

//______________________________________________________________________
Int_t AliITStrackerSA::FindTracks(AliESD* event){

// Track finder using the ESD object


  //controllare numero cluster sui layer1 e 2 (morti?)
  //non trova tracce...controllare..

  if(!fITSclusters){
    Fatal("FindTracks","ITS cluster tree is not accessed - Abort!!!\n Please use method SetClusterTree to pass the pointer to the tree\n");
    return -1;
  }
  
   
  //Reads event and mark clusters of traks already found, with flag kITSin
   Int_t nentr=event->GetNumberOfTracks();
   while (nentr--) {
     AliESDtrack *track=event->GetTrack(nentr);
     if (track->GetStatus()&AliESDtrack::kITSin==AliESDtrack::kITSin){
       Int_t idx[12];
       Int_t ncl = track->GetITSclusters(idx);
       for(Int_t k=0;k<ncl;k++){
       	 AliITSRecPoint* cll = (AliITSRecPoint*)GetCluster(idx[k]);
	 cll->SetBit(kSAflag);
       }
     }
   }

   Double_t primaryVertex[3];
   event->GetVertex()->GetXYZ(primaryVertex);
   //Creates TClonesArray with clusters for each layer. The clusters already used
   //by AliITStrackerMI are not considered
   
   Int_t nclusters[6]={0,0,0,0,0,0};
   Int_t dmar[6]={0,0,0,0,0,0};
   fCluLayer = new TClonesArray*[fGeom->GetNlayers()];
   fCluCoord = new TClonesArray*[fGeom->GetNlayers()];

   for(Int_t i=0;i<fGeom->GetNlayers();i++){
     AliITSlayer &layer=fgLayers[i];
     for(Int_t cli=0;cli<layer.GetNumberOfClusters();cli++){
       AliITSRecPoint* cls = (AliITSRecPoint*)layer.GetCluster(cli);
       if(cls->TestBit(kSAflag)==kTRUE) continue; //clusters used by TPC prol.
       if(cls->GetQ()==0) continue; //fake clusters dead zones
       nclusters[i]++;
     }
     dmar[i]=0;
     fCluLayer[i] = new TClonesArray("AliITSRecPoint",nclusters[i]);
     fCluCoord[i] = new TClonesArray("AliITSclusterTable",nclusters[i]);
   }

   
   Int_t * firstmod = new Int_t[fGeom->GetNlayers()];
   for(Int_t i=0;i<fGeom->GetNlayers();i++){
     firstmod[i]=fGeom->GetModuleIndex(i+1,1,1);
   }
   
   for(Int_t ilay=0;ilay<fGeom->GetNlayers();ilay++){
     TClonesArray &clulay = *fCluLayer[ilay];
     TClonesArray &clucoo = *fCluCoord[ilay];
     AliITSlayer &layer=fgLayers[ilay];
     for(Int_t cli=0;cli<layer.GetNumberOfClusters();cli++){
       AliITSRecPoint* cls = (AliITSRecPoint*)layer.GetCluster(cli);
       if(cls->TestBit(kSAflag)==kTRUE) continue;
       if(cls->GetQ()==0) continue;
       Double_t phi=0;Double_t lambda=0;
       Float_t x=0;Float_t y=0;Float_t z=0;
       Float_t sx=0;Float_t sy=0;Float_t sz=0;
       Int_t module = cls->GetDetectorIndex()+firstmod[ilay];
       GetCoorAngles(cls,module,phi,lambda,x,y,z,primaryVertex);
       GetCoorErrors(cls,module,sx,sy,sz);
       new (clulay[dmar[ilay]]) AliITSRecPoint(*cls);
       new (clucoo[dmar[ilay]]) AliITSclusterTable(x,y,z,sx,sy,sz,phi,lambda,cli);
       dmar[ilay]++;
     }
   }
   
 
  //Get primary vertex

   Int_t ntrack=0;
   //loop on the different windows
   for(Int_t nloop=0;nloop<fNloop;nloop++){
     for(Int_t ncl=0;ncl<fCluLayer[0]->GetEntries();ncl++){ //loop starting from layer 0
       
       ResetForFinding();
       Int_t pflag=0;
       
       AliITSRecPoint* cl = (AliITSRecPoint*)fCluLayer[0]->At(ncl);

       if(!cl) continue;


       if (cl->GetQ()<=0) continue;
       
       AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(0,ncl); 
       fPhic = arr->GetPhi();
       fLambdac = arr->GetLambda();
       if (TMath::Abs(fLambdac)>0.26*TMath::Pi()) continue;
       fPhiEstimate = fPhic;
       AliITStrackSA* trs = new AliITStrackSA(); 
       fPoint1[0]=primaryVertex[0];
       fPoint1[1]=primaryVertex[1];


       fPoint2[0]=arr->GetX();
       fPoint2[1]=arr->GetY();
       Int_t * nn = new Int_t[fGeom->GetNlayers()];//counter for clusters on each layer
       for(Int_t i=0;i<fGeom->GetNlayers();i++){ nn[i]=0;}
       nn[0] = SearchClusters(0,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag);
       
       nn[1] = SearchClusters(1,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag);
      if(nn[1]>0){
        pflag=1;
        fPoint3[0] = fPointc[0];
        fPoint3[1] = fPointc[1];
      }
      nn[2] = SearchClusters(2,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag);
      if(nn[1]==0 && nn[2]==0) pflag=0;
      if(nn[2]!=0 && nn[1]!=0){ pflag=1; UpdatePoints();}
      if(nn[2]!=0 && nn[1]==0){
        pflag=1;
        fPoint3[0]=fPointc[0];
        fPoint3[1]=fPointc[1];
      }

      nn[3] = SearchClusters(3,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag);
      pflag=1;
      if(nn[3]!=0) UpdatePoints();
      nn[4] = SearchClusters(4,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag); 
      pflag=1;
      if(nn[4]!=0) UpdatePoints();
      nn[5] = SearchClusters(5,fPhiWin[nloop],fLambdaWin[nloop],trs,primaryVertex[2],pflag); 
          

      Int_t layOK=0;
      Int_t numberofpoints;
      if(fSixPoints) numberofpoints=6;  //check of the candidate track
      else numberofpoints=5;           //if track is good (with the required number        
      for(Int_t nnp=0;nnp<fGeom->GetNlayers();nnp++){    //of points) it is written on file
        if(nn[nnp]!=0) layOK+=1;
      }
      if(layOK>=numberofpoints){

        AliITStrackV2* tr2 = FitTrack(trs,primaryVertex);
	
        if(tr2==0) continue;

	AliESDtrack outtrack;
	outtrack.UpdateTrackParams(tr2,AliESDtrack::kITSin);
	event->AddTrack(&outtrack);
	ntrack++;
	
      }
	

      delete trs;
      delete[] nn;
      
     }//end loop on clusters of layer1
    
     //end loop2
   }
   
   //if 5/6 points are required, second loop starting 
   //from second layer, to find tracks with point of 
   //layer 1 missing
   
   if(!fSixPoints){
     //   counter for clusters on each layer  
     Int_t * nn = new Int_t[fGeom->GetNlayers()-1];      
     for(Int_t nloop=0;nloop<fNloop;nloop++){
       Int_t ncl2=fCluLayer[1]->GetEntries();
       while(ncl2--){ //loop starting from layer 2
	 ResetForFinding();
	 Int_t pflag=0;
	 AliITSRecPoint* cl = (AliITSRecPoint*)fCluLayer[1]->At(ncl2);
	 
	 if(!cl) continue;
	 AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(1,ncl2);
	 fPhic = arr->GetPhi();
	 fLambdac = arr->GetLambda();
	 fPhiEstimate = fPhic;
         
	 AliITStrackSA* trs = new AliITStrackSA(); 
	 fPoint1[0]=primaryVertex[0];
	 fPoint1[1]=primaryVertex[1];
	 
	 fPoint2[0]=arr->GetX();
	 fPoint2[1]=arr->GetY();
	 for(Int_t kk=0;kk<fGeom->GetNlayers()-1;kk++)nn[kk] = 0;
	 nn[0] = SearchClusters(1,fPhiWin[nloop],fLambdaWin[nloop],
				   trs,primaryVertex[2],pflag);
	 nn[1] = SearchClusters(2,fPhiWin[nloop],fLambdaWin[nloop],
				   trs,primaryVertex[2],pflag);
	 if(nn[1]!=0){
	   pflag=1;
	   fPoint3[0]=fPointc[0];
	   fPoint3[1]=fPointc[1];
	 }
	 nn[2]= SearchClusters(3,fPhiWin[nloop],fLambdaWin[nloop],
				   trs,primaryVertex[2],pflag);
	 if(nn[2]!=0){
	   pflag=1;
	   UpdatePoints();
	 }
	 nn[3]= SearchClusters(4,fPhiWin[nloop],fLambdaWin[nloop],
				   trs,primaryVertex[2],pflag);
	 if(nn[3]!=0){
	   pflag=1;
	   UpdatePoints();
	 }
	 nn[4]=SearchClusters(5,fPhiWin[nloop],fLambdaWin[nloop],
				   trs,primaryVertex[2],pflag);

	 Int_t fl=0;
	 for(Int_t nnp=0;nnp<fGeom->GetNlayers()-1;nnp++){
	   if(nn[nnp]!=0) fl+=1;
	 }
	 if(fl>=5){  // 5/6       
	   AliITStrackV2* tr2 = FitTrack(trs,primaryVertex);
	   if(tr2==0){
	     continue;
	   }
	   
	   AliESDtrack outtrack;
	   outtrack.UpdateTrackParams(tr2,AliESDtrack::kITSin);
	   event->AddTrack(&outtrack);
	   ntrack++;

	 }   
	
       
        delete trs;
      }//end loop on clusters of layer2
     }
   
    delete [] nn;
  }  //end opt="5/6"  
   
  delete [] firstmod;
  Info("FindTracks","Number of found tracks: %d",event->GetNumberOfTracks());
  return 0;

}
 




//________________________________________________________________________

AliITStrackV2* AliITStrackerSA::FitTrack(AliITStrackSA* tr,Double_t *primaryVertex){
  //fit of the found track

  
  Int_t * firstmod = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    firstmod[i]=fGeom->GetModuleIndex(i+1,1,1);
  }  

  Int_t nclusters = tr->GetNumberOfClustersSA();
  TObjArray** listlayer = new TObjArray*[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    listlayer[i] = new TObjArray(0,0);
  }

  TArrayI clind0(20);
  TArrayI clind1(20);
  TArrayI clind2(20);
  TArrayI clind3(20);
  TArrayI clind4(20);
  TArrayI clind5(20);

  TArrayI mark0(20);
  TArrayI mark1(20);
  TArrayI mark2(20);
  TArrayI mark3(20);
  TArrayI mark4(20);
  TArrayI mark5(20);


  Int_t * nnn = new Int_t[fGeom->GetNlayers()];
  Int_t * kkk = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++) {nnn[i]=0;kkk[i]=0;}
  
  for(Int_t ncl=0;ncl<nclusters;ncl++){
    Int_t index = tr->GetClusterIndexSA(ncl); 
    AliITSRecPoint* cl = (AliITSRecPoint*)GetCluster(index);
    if(cl->TestBit(kSAflag)==kTRUE) cl->ResetBit(kSAflag);
    Int_t lay = (index & 0xf0000000) >> 28;
    if(lay==0) { listlayer[0]->AddLast(cl); clind0[nnn[0]]=index;nnn[0]++;}
    if(lay==1) { listlayer[1]->AddLast(cl); clind1[nnn[1]]=index;nnn[1]++;}
    if(lay==2) { listlayer[2]->AddLast(cl); clind2[nnn[2]]=index;nnn[2]++;}
    if(lay==3) { listlayer[3]->AddLast(cl); clind3[nnn[3]]=index;nnn[3]++;}
    if(lay==4) { listlayer[4]->AddLast(cl); clind4[nnn[4]]=index;nnn[4]++;}
    if(lay==5) { listlayer[5]->AddLast(cl); clind5[nnn[5]]=index;nnn[5]++;}    
  }
  delete [] nnn;

  for(Int_t nlay=0;nlay<fGeom->GetNlayers();nlay++){
    for(Int_t ncl=0;ncl<tr->GetNumberOfMarked(nlay);ncl++){
      Int_t mark = tr->GetClusterMark(nlay,ncl);
      if(nlay==0) { mark0[kkk[0]]=mark;kkk[0]++;}
      if(nlay==1) { mark1[kkk[1]]=mark;kkk[1]++;}
      if(nlay==2) { mark2[kkk[2]]=mark;kkk[2]++;}
      if(nlay==3) { mark3[kkk[3]]=mark;kkk[3]++;}
      if(nlay==4) { mark4[kkk[4]]=mark;kkk[4]++;}
      if(nlay==5) { mark5[kkk[5]]=mark;kkk[5]++;}

    }
  }

  delete [] kkk;

 
  Int_t * end = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    if(listlayer[i]->GetEntries()==0) end[i]=1;
    else end[i]=listlayer[i]->GetEntries();
  }

  TClonesArray* listSA = new TClonesArray("AliITStrackSA");
  TClonesArray &tri = *listSA;
  Int_t nlist=0;

  if(end[0]==0) end[0]=1; //for tracks with cluster on layer 0 missing
  for(Int_t l1=0;l1<end[0];l1++){//loop on layer 1
    AliITSRecPoint* cl0 = (AliITSRecPoint*)listlayer[0]->At(l1); 
    Double_t x1,y1,z1,sx1,sy1,sz1;
    Double_t x2,y2,z2,sx2,sy2,sz2;
    AliITSRecPoint* p1=0;
    AliITSRecPoint* p2=0;
    Int_t index1=clind0[l1];
    Int_t index2=0;
    Int_t mrk1 = mark0[l1];
    Int_t mrk2 = 0;
    Int_t lay1=0;
    Int_t lay2=1;
    Int_t module1=-1;
    for(Int_t l2=0;l2<end[1];l2++){//loop on layer 2
      AliITSRecPoint* cl1 = (AliITSRecPoint*)listlayer[1]->At(l2); 
      index2=clind1[l2];
      mrk2 = mark1[l2];
      for(Int_t l3=0;l3<end[2];l3++){  //loop on layer 3
        AliITSRecPoint* cl2 = (AliITSRecPoint*)listlayer[2]->At(l3);

        if(cl0==0 && cl1!=0) {
          p2 = cl2;index1=clind2[l3];mrk1=mark2[l3];lay1=2;
          p1=cl1;
          module1 = p1->GetDetectorIndex()+firstmod[1]; 
        }
        if(cl0!=0 && cl1==0){
          p1=cl0;
          p2=cl2;index2=clind2[l3];mrk2=mark2[l3];lay2=2;
	  module1 = p1->GetDetectorIndex()+firstmod[0];
        }
        if(cl0!=0 && cl1!=0){
          p1=cl0;
          p2=cl1;
	  module1 = p1->GetDetectorIndex()+firstmod[0];
        }
	
	Int_t cln1=mrk1;
	Int_t cln2=mrk2;
	AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(lay1,cln1);
	AliITSclusterTable* arr1= (AliITSclusterTable*)GetClusterCoord(lay2,cln2);
	x1 = arr->GetX();
	x2 = arr1->GetX();
	y1 = arr->GetY();
	y2 = arr1->GetY();
	z1 = arr->GetZ();
	z2 = arr1->GetZ();
	sx1 = arr->GetSx();
	sx2 = arr1->GetSx();
	sy1 = arr->GetSy();
	sy2 = arr1->GetSy();
	sz1 = arr->GetSz();
	sz2 = arr1->GetSz();
	
	Int_t layer,ladder,detector;
	fGeom->GetModuleId(module1,layer,ladder,detector);
        Float_t yclu1 = p1->GetY();
        Float_t zclu1 = p1->GetZ();
        Double_t cv=Curvature(primaryVertex[0],primaryVertex[1],x1,y1,x2,y2);
        
	Double_t tgl2 = (z2-z1)/TMath::Sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	Double_t phi2 = TMath::ATan2((y2-y1),(x2-x1));

        for(Int_t l4=0;l4<end[3];l4++){ //loop on layer 4   
          AliITSRecPoint* cl3 = (AliITSRecPoint*)listlayer[3]->At(l4);
          for(Int_t l5=0;l5<end[4];l5++){ //loop on layer 5
            AliITSRecPoint* cl4 = (AliITSRecPoint*)listlayer[4]->At(l5);
            for(Int_t l6=0;l6<end[5];l6++){ //loop on layer 6  
              AliITSRecPoint* cl5 = (AliITSRecPoint*)listlayer[5]->At(l6);
              AliITStrackSA* trac = new AliITStrackSA(fGeom,layer,ladder,detector,yclu1,zclu1,phi2,tgl2,cv,1);
                              
              if(cl5!=0) {
		trac->AddClusterV2(5,(clind5[l6] & 0x0fffffff)>>0);
		trac->AddClusterMark(5,mark5[l6]);
	      }
              if(cl4!=0){
		trac->AddClusterV2(4,(clind4[l5] & 0x0fffffff)>>0);
		trac->AddClusterMark(4,mark4[l5]);
	      }
              if(cl3!=0){
		trac->AddClusterV2(3,(clind3[l4] & 0x0fffffff)>>0);
		trac->AddClusterMark(3,mark3[l4]);
	      }
              if(cl2!=0){
		trac->AddClusterV2(2,(clind2[l3] & 0x0fffffff)>>0);
		trac->AddClusterMark(2,mark2[l3]);
	      }
              if(cl1!=0){
		trac->AddClusterV2(1,(clind1[l2] & 0x0fffffff)>>0);
		trac->AddClusterMark(1,mark1[l2]);
	      }
              if(cl0!=0){
		trac->AddClusterV2(0,(clind0[l1] & 0x0fffffff)>>0);
		trac->AddClusterMark(0,mark0[l1]);
	      }
              //fit with Kalman filter using AliITStrackerMI::RefitAt()
          

              AliITStrackMI* ot = new AliITStrackSA(*trac);
              
              ot->ResetCovariance(10.);
              ot->ResetClusters();
              
              if(RefitAt(49.,ot,trac)){ //fit from layer 1 to layer 6
                AliITStrackMI *otrack2 = new AliITStrackMI(*ot);
                otrack2->ResetCovariance(10.); 
                otrack2->ResetClusters();
                //fit from layer 6 to layer 1
                if(RefitAt(3.7,otrack2,ot)) {
		  fListOfTracks->AddLast(otrack2);
		  new (tri[nlist]) AliITStrackSA(*trac);
		  nlist++;
		} else {
		  delete otrack2;
		}
                              
              }       
          
              delete ot;
              delete trac;
            }//end loop layer 6
          }//end loop layer 5
        }//end loop layer 4
        
      }//end loop layer 3
    }//end loop layer 2 
  }//end loop layer 1

  delete [] end;

 

  Int_t dim=fListOfTracks->GetEntries();
  if(dim==0){
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      delete listlayer[i];
    }
    delete [] listlayer;
    listSA->Delete();
    delete listSA;
    delete [] firstmod;
    return 0;
  }

  Int_t lowchi2 = FindTrackLowChiSquare(fListOfTracks,dim);
  AliITStrackV2* otrack =(AliITStrackV2*)fListOfTracks->At(lowchi2);
  AliITStrackSA* trsa = (AliITStrackSA*)listSA->At(lowchi2);
 
  if(otrack==0) {
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      delete listlayer[i];
    }
    delete [] listlayer; 
    listSA->Delete();
    delete listSA;
    delete [] firstmod;
    return 0;
  }
  Int_t * indexc = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++) indexc[i]=0;
  for(Int_t nind=0;nind<otrack->GetNumberOfClusters();nind++){
    indexc[nind] = otrack->GetClusterIndex(nind);
  }      
  AliITSRecPoint* cl0 = (AliITSRecPoint*)GetCluster(indexc[0]);
  AliITSRecPoint* cl1 = (AliITSRecPoint*)GetCluster(indexc[1]);     
  AliITSRecPoint* cl2 = (AliITSRecPoint*)GetCluster(indexc[2]);     
  AliITSRecPoint* cl3 = (AliITSRecPoint*)GetCluster(indexc[3]);
  AliITSRecPoint* cl4 = (AliITSRecPoint*)GetCluster(indexc[4]);
  Int_t labl[3]={-1,-1,-1};
  if(otrack->GetNumberOfClusters()==fGeom->GetNlayers()){
    AliITSRecPoint* cl5 = (AliITSRecPoint*)GetCluster(indexc[5]);
    labl[0]=cl5->GetLabel(0);
    labl[1]=cl5->GetLabel(1);
    labl[2]=cl5->GetLabel(2);
  }
  delete [] indexc;
  if(otrack->GetNumberOfClusters()==(fGeom->GetNlayers()-1)){
    labl[0]=-1;
    labl[1]=-1;
    labl[2]=-1;
  }
  Int_t numberofpoints;
  if(fSixPoints) numberofpoints=6;
  else numberofpoints=5;
  CookLabel(otrack,0.); //MI change - to see fake ratio
  Int_t label =  Label(cl0->GetLabel(0),cl1->GetLabel(0), 
                       cl2->GetLabel(0),cl3->GetLabel(0),
                       cl4->GetLabel(0),labl[0],
                       cl0->GetLabel(1),cl1->GetLabel(1),
                       cl2->GetLabel(1),cl3->GetLabel(1),
                       cl4->GetLabel(1),labl[1],
                       cl0->GetLabel(2),cl1->GetLabel(2),
                       cl2->GetLabel(2),cl3->GetLabel(2),
                       cl4->GetLabel(2),labl[2],numberofpoints);
  
  otrack->SetLabel(label);  
  //remove clusters of found track

  for(Int_t nlay=0;nlay<fGeom->GetNlayers();nlay++){
    for(Int_t cln=0;cln<trsa->GetNumberOfMarked(nlay);cln++){
      Int_t index = trsa->GetClusterMark(nlay,cln);
      fCluLayer[nlay]->RemoveAt(index);
      RemoveClusterCoord(nlay,index);
      fCluLayer[nlay]->Compress();
    }    
  }
  listSA->Delete();
  delete listSA;

  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    delete listlayer[i];
  }
  delete [] listlayer; 
  delete [] firstmod;
  return otrack;

}



//_______________________________________________________
Int_t AliITStrackerSA::SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, AliITStrackSA* trs,Double_t /*zvertex*/,Int_t pflag){
  //function used to to find the clusters associated to the track
  Int_t nc=0;
  AliITSlayer &lay = fgLayers[layer];
  Double_t r=lay.GetR();
   if(pflag==1){      
    Float_t cx1,cx2,cy1,cy2;
    FindEquation(fPoint1[0],fPoint1[1],fPoint2[0],fPoint2[1],fPoint3[0],fPoint3[1],fCoef1,fCoef2,fCoef3);
    if (FindIntersection(fCoef1,fCoef2,fCoef3,-r*r,cx1,cy1,cx2,cy2)==0)
       return 0;
    Double_t fi1=TMath::ATan2(cy1,cx1);
    Double_t fi2=TMath::ATan2(cy2,cx2);
    fPhiEstimate=ChoosePoint(fi1,fi2,fPhic);
  }

 
  Int_t ncl = fCluLayer[layer]->GetEntries();
  for (Int_t index=0; index<ncl; index++) {
    AliITSRecPoint *c = (AliITSRecPoint*)fCluLayer[layer]->At(index);
    if (!c) continue;
    if (c->GetQ()<=0) continue;
    
     AliITSclusterTable* arr = (AliITSclusterTable*)GetClusterCoord(layer,index);
     Double_t phi = arr->GetPhi();
     if (TMath::Abs(phi-fPhiEstimate)>phiwindow) continue;

     Double_t lambda = arr->GetLambda();
     if (TMath::Abs(lambda-fLambdac)>lambdawindow) continue;

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
Int_t AliITStrackerSA::FindTrackLowChiSquare(TObjArray* tracklist, Int_t dim) const {
  // returns track with lowes chi square  
  if(dim==1){
    //AliITStrackV2* trk = (AliITStrackV2*)tracklist->At(0);
    //return trk;
    return 0;
  }
  //if(dim==0) return 0;
  Double_t * chi2 = new Double_t[dim];
  Int_t * index = new Int_t[dim];
  for(Int_t i=0;i<dim;i++){
    AliITStrackV2* trk = (AliITStrackV2*)tracklist->At(i);
    chi2[i]=trk->GetChi2();
    index[i]=i;
  }

  Int_t w=0;Double_t value;
  Int_t lp;
  while(w<dim){
    for(Int_t j=w+1;j<dim;j++){
      if(chi2[w]<chi2[j]){
        value=chi2[w];
        chi2[w]=chi2[j];
        chi2[j]=value;
        lp=index[w];
        index[w]=index[j];
        index[j]=lp;
      }
    }
    w++;
  }

  delete [] chi2;
  delete [] index;
  return index[dim-1];
}

//__________________________________________________________
Int_t AliITStrackerSA::FindLabel(Int_t l1, Int_t l2, Int_t l3, Int_t l4, Int_t l5, Int_t l6){

  //function used to determine the track label
  
  Int_t lb[6] = {l1,l2,l3,l4,l5,l6};
  Int_t aa[6]={1,1,1,1,1,1};
  Int_t ff=0; 
  Int_t ll=0;
  Int_t k=0;Int_t w=0;Int_t num=6;
  if(lb[5]==-1) num=5;
  
  while(k<num){
  
    for(Int_t i=k+1;i<num;i++){
    
      if(lb[k]==lb[i] && aa[k]!=0){
      
        aa[k]+=1;
        aa[i]=0;
      }
    }
  k++;
  }

  while(w<num){
  
    for(Int_t j=0;j<6;j++){
      if(aa[w]<aa[j]){
      ff=aa[w];
      aa[w]=aa[j];
      aa[j]=ff;
      ll=lb[w];
      lb[w]=lb[j];
      lb[j]=ll;
     }
    }
  w++;
  }
  if(num==6)  return lb[5];
  else return lb[4];
}

//_____________________________________________________________________________
Int_t AliITStrackerSA::Label(Int_t gl1, Int_t gl2, Int_t gl3, Int_t gl4, Int_t gl5, Int_t gl6,Int_t gl7, Int_t gl8, Int_t gl9, Int_t gl10,Int_t gl11,
Int_t gl12, Int_t gl13, Int_t gl14,Int_t gl15, Int_t gl16, Int_t gl17, Int_t gl18, Int_t numberofpoints){

 
  //function used to assign label to the found track. If track is fake, the label is negative

  Int_t lb0[6] = {gl1,gl2,gl3,gl4,gl5,gl6};
  Int_t lb1[6] = {gl7,gl8,gl9,gl10,gl11,gl12};
  Int_t lb2[6] = {gl13,gl14,gl15,gl16,gl17,gl18};
  Int_t ll=FindLabel(lb0[0],lb0[1],lb0[2],lb0[3],lb0[4],lb0[5]);
  Int_t lflag=0;Int_t num=6;
  if(lb0[5]==-1 && lb1[5]==-1 && lb2[5]==-1) num=5;

  for(Int_t i=0;i<num;i++){
    if(lb0[i]==ll || lb1[i]==ll || lb2[i]==ll) lflag+=1;
  }

  if(lflag>=numberofpoints) return ll;
  else return -ll;

  
}

//_____________________________________________________________________________
void AliITStrackerSA::SetWindowSizes(Int_t n, Double_t *phi, Double_t *lam){
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
   
    for(Int_t k=0;k<fNloop;k++){
      fPhiWin[k]=phid[k];
      fLambdaWin[k]=lambdad[k];
    }
  
  }

}
//_______________________________________________________________________
void AliITStrackerSA::GetCoorAngles(AliITSRecPoint* cl,Int_t module,Double_t &phi,Double_t &lambda, Float_t &x, Float_t &y,Float_t &z,Double_t* vertex){
  //Returns values of phi (azimuthal) and lambda angles for a given cluster
  
  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; fGeom->GetModuleId(module,lay,lad,det);
  Float_t tx,ty,tz;  fGeom->GetTrans(lay,lad,det,tx,ty,tz);     

  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi1=TMath::Pi()/2+alpha;
  if (lay==1) phi1+=TMath::Pi();

  Float_t cp=TMath::Cos(phi1), sp=TMath::Sin(phi1);
  Float_t r=tx*cp+ty*sp;

  x= r*cp - cl->GetY()*sp;
  y= r*sp + cl->GetY()*cp;
  z=cl->GetZ();
  
  phi=TMath::ATan2(y,x);
  lambda=TMath::ATan2(z-vertex[2],TMath::Sqrt((x-vertex[0])*(x-vertex[0])+(y-vertex[1])*(y-vertex[1])));
}

//________________________________________________________________________
void AliITStrackerSA::GetCoorErrors(AliITSRecPoint* cl, Int_t module,Float_t &sx,Float_t &sy, Float_t &sz){

  //returns x,y,z of cluster in global coordinates

  Double_t rot[9];     fGeom->GetRotMatrix(module,rot);
  Int_t lay,lad,det; fGeom->GetModuleId(module,lay,lad,det);
 
  Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
  Double_t phi=TMath::Pi()/2+alpha;
  if (lay==1) phi+=TMath::Pi();

  Float_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);

  sx = TMath::Sqrt(sp*sp*cl->GetSigmaY2());
  sy = TMath::Sqrt(cp*cp*cl->GetSigmaY2());
  sz = TMath::Sqrt(cl->GetSigmaZ2());

}
