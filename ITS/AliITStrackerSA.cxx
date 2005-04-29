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
#include <TMath.h>
#include <TObjArray.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliITSRiemannFit.h"
#include "AliITSVertexer.h"
#include "AliITSclusterTable.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliITStrackSA.h"
#include "AliITStrackerSA.h"
#include "AliRun.h"

ClassImp(AliITStrackerSA)

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA():AliITStrackerMI(){
  // Default constructor
  Init();
 
}
//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom):AliITStrackerMI(geom) 
{
  // Standard constructor (Vertex is known and passed to this obj.)
  Init();
  fVert = 0;
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom, AliESDVertex *vert):AliITStrackerMI(geom) 
{
  // Standard constructor (Vertex is known and passed to this obj.)
  Init();
  fVert = vert;
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITSgeom *geom, AliITSVertexer *vertexer):AliITStrackerMI(geom) 
{
  // Standard constructor (Vertex is unknown - vertexer is passed to this obj)
  Init();
  fVertexer = vertexer;
  fGeom = geom;
 
}

//____________________________________________________________________________
AliITStrackerSA::AliITStrackerSA(AliITStrackerSA& tracker):AliITStrackerMI(){
  // Copy constructor
  fPhiEstimate = tracker.fPhiEstimate;
  for(Int_t i=0;i<2;i++){
    fPoint1[i]=tracker.fPoint1[i];
    fPoint2[i]=tracker.fPoint2[i];
    fPoint3[i]=tracker.fPoint3[i];
    fPointc[i]=tracker.fPointc[i];
  }
  fLambdac = tracker.fLambdac;
  fPhic = tracker.fPhic;
  fCoef1 = tracker.fCoef1;
  fCoef2 = tracker.fCoef2;
  fCoef3 = tracker.fCoef3;
  fNloop = tracker.fNloop;
  fPhiWin = tracker.fPhiWin;
  fLambdaWin = tracker.fLambdaWin;
  if(tracker.fVertexer && tracker.fVert){
    fVert = new AliESDVertex(*tracker.fVert);
  }
  else {
    fVert = tracker.fVert;
  }
  fVertexer = tracker.fVertexer;
  fGeom = tracker.fGeom;
  fTable = tracker.fTable;
  fListOfTracks = tracker.fListOfTracks;
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
  fGeom = 0;
  if(fPhiWin)delete []fPhiWin;
  if(fLambdaWin)delete []fLambdaWin;
  fTable =0;
  fListOfTracks->Delete();
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
    fTable = 0;
    fITSclusters = 0;
    SetSixPoints();
    fListOfTracks=new TObjArray(0,0);
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
//____________________________________________________________________________
void AliITStrackerSA::FindTracks(TTree *out,Int_t evnumber){

  /**************************************************************************
   * This function finds primary tracks.
   *                                                                        * 
   *                                                                        *
   * Example: to execute function with only the ITS (no combined tracking   *
   *          with TPC+ITS) and requiring 5/6 points to define a good track *
   *          call SetSixPoinbts(kFALSE) in advance and then                *
   *          use: FindTracks(treein,treeout,evnumber)                      *
   *          to execute combined tracking, before using FindTracks, use    *
   *          UseFoundTracksV2                                              *
   *************************************************************************/

  if(!fITSclusters){
    Fatal("FindTracks","ITS cluster tree is not accessed - Abort!!!\n Please use method SetClusterTree to pass the pointer to the tree\n");
    exit(1);
  }
  //Get primary vertex
  if(fVertexer){
    if(fVert)delete fVert;
    fVert = fVertexer->FindVertexForCurrentEvent(evnumber);
  }
  else {
    gAlice->GetEvent(evnumber);
    if(!fVert){
      Fatal("FindTracks","Vertex is missing\n");
      return;
    }
  }
  Double_t primaryVertex[3];
  Double_t errorsprimvert[3];
  fVert->GetXYZ(primaryVertex);
  fVert->GetSigmaXYZ(errorsprimvert);
  if(errorsprimvert[0]==0 || errorsprimvert[1]==0){
    Warning("FindTracks","Set errors on vertex positions x and y at 0.0001");
    errorsprimvert[0]=0.0001;
    errorsprimvert[1]=0.0001;
  }
  fVert->PrintStatus();


  //Fill array with cluster indices for each module
  if(!fTable){
    fTable = new AliITSclusterTable(fGeom,this,primaryVertex);
    fTable->FillArray(fITSclusters);
    fTable->FillArrayCoorAngles(); 
  }

   
//Fill tree for found tracks
  AliITStrackV2* outrack=0;
  TBranch* branch=out->Branch("tracks","AliITStrackV2",&outrack,32000,0);
  if (!branch) out->Branch("tracks","AliITStrackV2",&outrack,32000,3);
  else branch->SetAddress(&outrack);

  
  Int_t * firstmod = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    firstmod[i]=fGeom->GetModuleIndex(i+1,1,1);
  }
  // firstmod [i] number of the first module in the ITS layer i.
    
  AliITSlayer &layer=fgLayers[0];   // first layer
  Int_t ntrack=0;
  Int_t dim=layer.GetNumberOfClusters();
  //loop on the different windows
  for(Int_t nloop=0;nloop<fNloop;nloop++){
    for(Int_t ncl=0;ncl<dim;ncl++){ 
      //loop starting from layer 0
      ResetForFinding();
      Int_t pflag=0;
      AliITSclusterV2* cl = layer.GetCluster(ncl);
      if(cl->IsUsed()==1) continue;
      if(cl->TestBit(kSAflag)==kTRUE) continue;
      
      fPhic = fTable->GetPhiCluster(0,ncl);
      fLambdac = fTable->GetLambdaCluster(0,ncl);
      fPhiEstimate = fPhic;
      AliITStrackSA* trs = new AliITStrackSA();      
      fPoint1[0]=primaryVertex[0];
      fPoint1[1]=primaryVertex[1];
      fPoint2[0]=fTable->GetXCluster(0,ncl);
      fPoint2[1]=fTable->GetYCluster(0,ncl);
      
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
        AliITStrackV2* tr2 = FitTrack(trs,primaryVertex,errorsprimvert);
        if(tr2==0){
          Int_t nct = trs->GetNumberOfClustersSA();
          while(nct--){
            Int_t index = trs->GetClusterIndexSA(nct);
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
	    if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
          }
	  delete [] nn;
	  delete trs;
	  continue;
        }
        outrack=tr2;  
        out->Fill();
        ntrack++;
        Int_t nct = tr2->GetNumberOfClusters();

        while(nct--){
          Int_t index = tr2->GetClusterIndex(nct);     
          AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);      
      	  kl->SetBit(kSAflag);

        }
      } 
      else{
        Int_t nct = trs->GetNumberOfClustersSA();
        while(nct--){
          Int_t index = trs->GetClusterIndexSA(nct);
          AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
       	  if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
        }
      }
      delete [] nn;
      delete trs;

    }//end loop on clusters of layer1

  }//end loop2

  //if 5/6 points are required, second loop starting 
  //from second layer, to find tracks with point of 
  //layer 1 missing
   
  if(!fSixPoints){
    //   counter for clusters on each layer  
    Int_t * nn = new Int_t[fGeom->GetNlayers()-1];      
    for(Int_t nloop=0;nloop<fNloop;nloop++){
      AliITSlayer &layer2=fgLayers[1]; //loop on layer 2
      Int_t ncl2=layer2.GetNumberOfClusters();
      while(ncl2--){ //loop starting from layer 2
        ResetForFinding();
        Int_t pflag=0;
        AliITSclusterV2* cl = layer2.GetCluster(ncl2);
        if(cl->IsUsed()==1) continue;
	if(cl->TestBit(kSAflag)==kTRUE) continue;
	fPhic = fTable->GetPhiCluster(1,ncl2);
	fLambdac = fTable->GetLambdaCluster(1,ncl2);
        fPhiEstimate = fPhic;
        AliITStrackSA* trs = new AliITStrackSA(); 
     
        fPoint1[0]=primaryVertex[0];
        fPoint1[1]=primaryVertex[1];
        fPoint2[0]=fTable->GetXCluster(1,ncl2);;
        fPoint2[1]=fTable->GetYCluster(1,ncl2);;
   
        for(Int_t kk=0;kk<fGeom->GetNlayers()-1;kk++)nn[kk] = 0;
        for(Int_t kk=0;kk<fGeom->GetNlayers()-1;kk++){
          nn[kk] = SearchClusters(kk+1,fPhiWin[nloop],fLambdaWin[nloop],
                   trs,primaryVertex[2],pflag);
          if(nn[kk]==0)break;
          if(kk>0){
            UpdatePoints();
            pflag = 1;
          }
        }
        Int_t fl=0;
        for(Int_t nnp=0;nnp<fGeom->GetNlayers()-1;nnp++){
          if(nn[nnp]!=0) fl+=1;
        }
        if(fl>=5){  // 5/6       
          AliITStrackV2* tr2 = FitTrack(trs,primaryVertex,errorsprimvert);
          if(tr2==0){
            Int_t nct = trs->GetNumberOfClustersSA();
            while(nct--){
              Int_t index = trs->GetClusterIndexSA(nct);
              AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
	      if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
            }
	    delete trs;
            continue;
          }
          outrack=tr2;
          out->Fill();
          Int_t nct = tr2->GetNumberOfClusters();
          while(nct--){
            Int_t index = tr2->GetClusterIndex(nct);     
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
            if(kl==0) continue;
	    kl->SetBit(kSAflag);
          }
        }       
        else{
          Int_t nct = trs->GetNumberOfClustersSA();
          while(nct--){
            Int_t index = trs->GetClusterIndexSA(nct);
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
            if(kl==0) continue;
	    if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
	  }
        }
        delete trs;
      }//end loop on clusters of layer2
    }
    delete [] nn;
  }  // if(!fSixPoints....  

  delete [] firstmod;
  delete fTable; fTable=0;
}


//______________________________________________________________________
Int_t AliITStrackerSA::FindTracks(AliESD* event){

  // Track finder using the ESD object


  if(!fITSclusters){
    Fatal("FindTracks","ITS cluster tree is not accessed - Abort!!!\n Please use method SetClusterTree to pass the pointer to the tree\n");
    return -1;
  }
 
  //Get primary vertex
  Double_t errorsprimvert[3];
  Double_t primaryVertex[3];
  event->GetVertex()->GetXYZ(primaryVertex);
  event->GetVertex()->GetSigmaXYZ(errorsprimvert);

  if(errorsprimvert[0]==0 || errorsprimvert[1]==0){
    //    Warning("FindTracks","Set errors on vertex positions x and y at 0.005");
    errorsprimvert[0]=0.005;
    errorsprimvert[1]=0.005;
  }

  //Fill array with cluster indices for each module
  if(!fTable){
    fTable = new AliITSclusterTable(fGeom,this,primaryVertex);
    fTable->FillArray(fITSclusters);
    fTable->FillArrayCoorAngles();
  }

  Int_t * firstmod = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    firstmod[i]=fGeom->GetModuleIndex(i+1,1,1);
  }
  // firstmod [i] number of the first module in the ITS layer i.
  
    
  AliITSlayer &layer=fgLayers[0];   
  Int_t ntrack=0;
   Int_t dim=layer.GetNumberOfClusters();
  //loop on the different windows
  for(Int_t nloop=0;nloop<fNloop;nloop++){
    for(Int_t ncl=0;ncl<dim;ncl++){ //loop starting from layer 0

      ResetForFinding();
      Int_t pflag=0;
      AliITSclusterV2* cl = (AliITSclusterV2*)layer.GetCluster(ncl);
      if(cl==0) continue;
      if(cl->IsUsed()==1) continue;
      if(cl->TestBit(kSAflag)==kTRUE) continue;
      if (cl->GetQ()<=0) continue;
      
      fPhic = fTable->GetPhiCluster(0,ncl);
      fLambdac = fTable->GetLambdaCluster(0,ncl);

      if (TMath::Abs(fLambdac)>0.26*TMath::Pi()) continue;

      fPhiEstimate = fPhic;
      AliITStrackSA* trs = new AliITStrackSA();      
      fPoint1[0]=primaryVertex[0];
      fPoint1[1]=primaryVertex[1];
      fPoint2[0]=fTable->GetXCluster(0,ncl);
      fPoint2[1]=fTable->GetYCluster(0,ncl);
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
        AliITStrackV2* tr2 = FitTrack(trs,primaryVertex,errorsprimvert);
        if(tr2==0){
          Int_t nct = trs->GetNumberOfClustersSA();
          while(nct--){
            Int_t index = trs->GetClusterIndexSA(nct);
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
     	    if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
	    
          }
	  continue;
        }
	
	AliESDtrack outtrack;
	outtrack.UpdateTrackParams(tr2,AliESDtrack::kITSin);
	event->AddTrack(&outtrack);
	ntrack++;
	Int_t nct = tr2->GetNumberOfClusters();
	while(nct--){
          Int_t index = tr2->GetClusterIndex(nct);     
          AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);      
	  kl->SetBit(kSAflag);
          
	} 
      }
      else{
	Int_t nct = trs->GetNumberOfClustersSA();
	while(nct--){
	  Int_t index = trs->GetClusterIndexSA(nct);
	  AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
	  if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
	  
	}
      }
      delete trs;
      delete[] nn;
      
    }//end loop on clusters of layer1
    
  }//end loop2



  //if 5/6 points are required, second loop starting 
  //from second layer, to find tracks with point of 
  //layer 1 missing
   
  if(!fSixPoints){
    //   counter for clusters on each layer  
    Int_t * nn = new Int_t[fGeom->GetNlayers()-1];      
    for(Int_t nloop=0;nloop<fNloop;nloop++){
      AliITSlayer &layer2=fgLayers[1]; 
      Int_t ncl2=layer2.GetNumberOfClusters();
      while(ncl2--){ //loop starting from layer 2
        ResetForFinding();
        Int_t pflag=0;
        AliITSclusterV2* cl = layer2.GetCluster(ncl2);
        if(cl->IsUsed()==1) continue;
	if(cl->TestBit(kSAflag)==kTRUE) continue;
	fPhic = fTable->GetPhiCluster(1,ncl2);
	fLambdac = fTable->GetLambdaCluster(1,ncl2);
	fPhiEstimate = fPhic;
        AliITStrackSA* trs = new AliITStrackSA(); 
        fPoint1[0]=primaryVertex[0];
        fPoint1[1]=primaryVertex[1];
        fPoint2[0]=fTable->GetXCluster(1,ncl2);
        fPoint2[1]=fTable->GetYCluster(1,ncl2);
   
        for(Int_t kk=0;kk<fGeom->GetNlayers()-1;kk++)nn[kk] = 0;
        for(Int_t kk=0;kk<fGeom->GetNlayers()-1;kk++){
          nn[kk] = SearchClusters(kk+1,fPhiWin[nloop],fLambdaWin[nloop],
                   trs,primaryVertex[2],pflag);
          if(nn[kk]==0)break;
          if(kk>0){
            UpdatePoints();
            pflag = 1;
          }
        }
        Int_t fl=0;
        for(Int_t nnp=0;nnp<fGeom->GetNlayers()-1;nnp++){
          if(nn[nnp]!=0) fl+=1;
        }
        if(fl>=5){  // 5/6       
          AliITStrackV2* tr2 = FitTrack(trs,primaryVertex,errorsprimvert);
          if(tr2==0){
	    Int_t nct = trs->GetNumberOfClustersSA();
            while(nct--){
              Int_t index = trs->GetClusterIndexSA(nct);
              AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
  	      if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
	      
            }

            continue;
          }

	  AliESDtrack outtrack;
	  outtrack.UpdateTrackParams(tr2,AliESDtrack::kITSin);
	  event->AddTrack(&outtrack);
	  ntrack++;
          Int_t nct = tr2->GetNumberOfClusters();
          while(nct--){
            Int_t index = tr2->GetClusterIndex(nct);     
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
            if(kl==0) continue;
	    kl->SetBit(kSAflag);
	  }
        }       
        else{
          Int_t nct = trs->GetNumberOfClustersSA();
          while(nct--){
            Int_t index = trs->GetClusterIndexSA(nct);
            AliITSclusterV2* kl = (AliITSclusterV2*)GetCluster(index);
            if(kl==0) continue;
	    if(kl->TestBit(kSAflag)==kTRUE) kl->ResetBit(kSAflag);
	    
          }
        }
        delete trs;
      }//end loop on clusters of layer2
    }
    delete [] nn;
  }  //end opt="5/6"  

  delete [] firstmod;
  delete fTable;fTable=0;   
  Info("FindTracks","Number of found tracks: %d",event->GetNumberOfTracks());
  return 0;

}


//________________________________________________________________________

AliITStrackV2* AliITStrackerSA::FitTrack(AliITStrackSA* tr,Double_t *primaryVertex,Double_t *errorsprimvert){
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

  Int_t * nnn = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++)nnn[i]=0;
  
  for(Int_t ncl=0;ncl<nclusters;ncl++){
    Int_t index = tr->GetClusterIndexSA(ncl);   
    AliITSclusterV2* cl = (AliITSclusterV2*)GetCluster(index);

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

 
  Int_t * end = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    if(listlayer[i]->GetEntries()==0) end[i]=1;
    else end[i]=listlayer[i]->GetEntries();
  }

  for(Int_t l1=0;l1<end[0];l1++){//loop on layer 1
    AliITSclusterV2* cl0 = (AliITSclusterV2*)listlayer[0]->At(l1); 
    TVector3** recp = new TVector3*[3];
    TVector3** errs = new TVector3*[3];
    recp[0] = new TVector3(primaryVertex[0],primaryVertex[1],primaryVertex[2]);
    errs[0] = new TVector3(errorsprimvert[0],errorsprimvert[1],errorsprimvert[2]);
    Double_t x1,y1,z1,sx1,sy1,sz1;
    Double_t x2,y2,z2,sx2,sy2,sz2;
    AliITSclusterV2* p1=0;
    AliITSclusterV2* p2=0;
    Int_t index1=clind0[l1];
    Int_t index2=0;
    for(Int_t l2=0;l2<end[1];l2++){//loop on layer 2
      AliITSclusterV2* cl1 = (AliITSclusterV2*)listlayer[1]->At(l2); 
      index2=clind1[l2];
      for(Int_t l3=0;l3<end[2];l3++){  //loop on layer 3
        AliITSclusterV2* cl2 = (AliITSclusterV2*)listlayer[2]->At(l3);

        if(cl0==0 && cl1!=0) {
          p2 = cl2;index1=clind2[l3];
          p1=cl1;
                
        }
        if(cl0!=0 && cl1==0){
          p1=cl0;
          p2=cl2;index2=clind2[l3];
        }
        if(cl0!=0 && cl1!=0){
          p1=cl0;
          p2=cl1;
        }
	Int_t lay1=(index1 & 0xf0000000) >> 28;
	Int_t cln1=(index1 & 0x0fffffff) >> 00;
	Int_t lay2=(index2 & 0xf0000000) >> 28;
	Int_t cln2=(index2 & 0x0fffffff) >> 00;
	x1 = fTable->GetXCluster(lay1,cln1);
	x2 = fTable->GetXCluster(lay2,cln2);
	y1 = fTable->GetYCluster(lay1,cln1);
	y2 = fTable->GetYCluster(lay2,cln2);
	z1 = fTable->GetZCluster(lay1,cln1);
	z2 = fTable->GetZCluster(lay2,cln2);
	sx1 = fTable->GetXClusterError(lay1,cln1);
	sx2 = fTable->GetXClusterError(lay2,cln2);
	sy1 = fTable->GetYClusterError(lay1,cln1);
	sy2 = fTable->GetYClusterError(lay2,cln2);
	sz1 = fTable->GetZClusterError(lay1,cln1);
	sz2 = fTable->GetZClusterError(lay2,cln2);
        Double_t phi1 = fTable->GetPhiCluster(lay1,cln1);
        Int_t module1 = p1->GetDetectorIndex()+firstmod[0];
	recp[1] = new TVector3(x1,y1,z1);
        errs[1] = new TVector3(sx1,sy1,sz1);
        recp[2] = new TVector3(x2,y2,z2);
        errs[2] = new TVector3(sx2,sy2,sz2);
        
        //fit on the Riemann sphere
        Float_t seed1,seed2,seed3;
        AliITSRiemannFit fit;
        Int_t rf = fit.FitHelix(3,recp,errs,seed1,seed2,seed3); //this gives phi,tgl,curvature to start Kalman Filter
        if(rf==0) {
	  for(Int_t i=1;i<3;i++){
	    delete recp[i];
	    delete errs[i];
	  }     
	  continue;
	}  
        Double_t phi=seed1;
        Double_t tgl=seed2;
        
        if(phi1>0){ 
          if(seed1>-TMath::Pi() && seed1<-0.5*TMath::Pi()){
            phi=seed1+1.5*TMath::Pi();
            tgl=seed2; 
          }
          if(seed1>-0.5*TMath::Pi() && seed1<0.5*TMath::Pi()){
            phi=seed1+0.5*TMath::Pi();
            tgl=(-1)*seed2; 
          }
          if(seed1>0.5*TMath::Pi() && seed1<TMath::Pi()){
            phi=seed1-0.5*TMath::Pi();
            tgl=seed2; 
          }
        }      
        if(phi1<0){
          if(seed1>-TMath::Pi() && seed1<-0.5*TMath::Pi()){
            phi=seed1+0.5*TMath::Pi();
            tgl=(-1)*seed2; 
          }
          if(seed1>-0.5*TMath::Pi() && seed1<0.5*TMath::Pi()){
            phi=seed1-0.5*TMath::Pi();
            tgl=seed2; 
          }
          if(seed1>0.5*TMath::Pi() && seed1<TMath::Pi()){
            phi=seed1-1.5*TMath::Pi();
            tgl=(-1)*seed2; 
            }
          }
        
	Int_t layer,ladder,detector;
	fGeom->GetModuleId(module1,layer,ladder,detector);
        Float_t yclu1 = p1->GetY();
        Float_t zclu1 = p1->GetZ();
        Double_t cv=Curvature(primaryVertex[0],primaryVertex[1],x1,y1,x2,y2);
              
        for(Int_t l4=0;l4<end[3];l4++){ //loop on layer 4   
          AliITSclusterV2* cl3 = (AliITSclusterV2*)listlayer[3]->At(l4);
          for(Int_t l5=0;l5<end[4];l5++){ //loop on layer 5
            AliITSclusterV2* cl4 = (AliITSclusterV2*)listlayer[4]->At(l5);
            for(Int_t l6=0;l6<end[5];l6++){ //loop on layer 6  
              AliITSclusterV2* cl5 = (AliITSclusterV2*)listlayer[5]->At(l6);
              AliITStrackSA* trac = new AliITStrackSA(layer,ladder,detector,yclu1,zclu1,phi,tgl,cv,1);
                              
              if(cl5!=0) trac->AddClusterV2(5,(clind5[l6] & 0x0fffffff)>>0);
              if(cl4!=0) trac->AddClusterV2(4,(clind4[l5] & 0x0fffffff)>>0);
              if(cl3!=0) trac->AddClusterV2(3,(clind3[l4] & 0x0fffffff)>>0);
              if(cl2!=0) trac->AddClusterV2(2,(clind2[l3] & 0x0fffffff)>>0);
              if(cl1!=0) trac->AddClusterV2(1,(clind1[l2] & 0x0fffffff)>>0);
              if(cl0!=0) trac->AddClusterV2(0,(clind0[l1] & 0x0fffffff)>>0);
            
              //fit with Kalman filter using AliITStrackerMI::RefitAt()
          
              AliITStrackMI* ot = new AliITStrackSA(*trac);
              
              ot->ResetCovariance();
              ot->ResetClusters();
              
              if(RefitAt(49.,ot,trac)){ //fit from layer 1 to layer 6

                AliITStrackMI *otrack2 = new AliITStrackMI(*ot);
                otrack2->ResetCovariance(); 
                otrack2->ResetClusters();
                //fit from layer 6 to layer 1
                if(RefitAt(3.7,otrack2,ot)) {
		  fListOfTracks->AddLast(otrack2);
		} else {
		  delete otrack2;
		}
                              
              }       
          
              delete ot;
              delete trac;
            }//end loop layer 6
          }//end loop layer 5
        }//end loop layer 4
        
        for(Int_t i=1;i<3;i++){
          delete recp[i];
          delete errs[i];
        }
      }//end loop layer 3
    }//end loop layer 2 
    delete recp[0];
    delete errs[0];
    delete[] recp;
    delete[] errs;    
  }//end loop layer 1

  delete [] end;

  Int_t dim=fListOfTracks->GetEntries();
  if(dim==0){
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      delete listlayer[i];
    }
    delete [] listlayer;
    delete [] firstmod;
    return 0;
  }

  AliITStrackV2* otrack =(AliITStrackV2*)FindTrackLowChiSquare(fListOfTracks,dim);

  if(otrack==0) {
    for(Int_t i=0;i<fGeom->GetNlayers();i++){
      delete listlayer[i];
    }
    delete [] listlayer; 
    delete [] firstmod;
    return 0;
  }
  Int_t * indexc = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++) indexc[i]=0;
  for(Int_t nind=0;nind<otrack->GetNumberOfClusters();nind++){
    indexc[nind] = otrack->GetClusterIndex(nind);
  }      
  AliITSclusterV2* cl0 = (AliITSclusterV2*)GetCluster(indexc[0]);
  AliITSclusterV2* cl1 = (AliITSclusterV2*)GetCluster(indexc[1]);     
  AliITSclusterV2* cl2 = (AliITSclusterV2*)GetCluster(indexc[2]);     
  AliITSclusterV2* cl3 = (AliITSclusterV2*)GetCluster(indexc[3]);
  AliITSclusterV2* cl4 = (AliITSclusterV2*)GetCluster(indexc[4]);
  Int_t labl[3]={-1,-1,-1};
  if(otrack->GetNumberOfClusters()==fGeom->GetNlayers()){
    AliITSclusterV2* cl5 = (AliITSclusterV2*)GetCluster(indexc[5]);
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
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    delete listlayer[i];
  }
  delete [] listlayer; 
  delete [] firstmod;
  return otrack;

}

//_______________________________________________________________________
void AliITStrackerSA::UseFoundTracksV2(Int_t evnum,TTree* treev2){
  // Marks as used clusters belonging to tracks found with V2 TPC+ITS tracking
  //(or AliITStrackV2 tracks found with function FindTracks of this class)
  

  //Get primary vertex
  if(fVertexer){
    if(fVert)delete fVert;
    fVert = fVertexer->FindVertexForCurrentEvent(evnum);
  }
  else {
    gAlice->GetEvent(evnum);
    if(!fVert){
      Fatal("FindTracks","Vertex is missing\n");
      return;
    }
  }
  Double_t primaryVertex[3];
  fVert->GetXYZ(primaryVertex);

  if(!fTable){
    fTable = new AliITSclusterTable(fGeom,this,primaryVertex);
    fTable->FillArray(fITSclusters);
    fTable->FillArrayCoorAngles();
  }

  TBranch* bra = (TBranch*)treev2->GetBranch("tracks");
  if(!bra) Warning("UseFoundTracksV2","No branch for track tree");
  AliITStrackV2* ttrrt = new AliITStrackV2;
  bra->SetAddress(&ttrrt);

  for(Int_t nj=0;nj<treev2->GetEntries();nj++){
    treev2->GetEvent(nj);
    Int_t ncl = ttrrt->GetNumberOfClusters();
    for(Int_t k=0;k<ncl;k++){
      Int_t index = ttrrt->GetClusterIndex(k);
      AliITSclusterV2* clui = (AliITSclusterV2*)GetCluster(index);
      if(clui->IsUsed()==0) clui->Use();
      
    }
  }
  delete ttrrt;
  
}

//_______________________________________________________________________
void AliITStrackerSA::UseFoundTracksV2(AliESD *event){
  // Marks as used clusters belonging to tracks found with V2 TPC+ITS tracking

  //Get primary vertex

  Double_t primaryVertex[3];
  event->GetVertex()->GetXYZ(primaryVertex);

  if(!fTable){
    fTable = new AliITSclusterTable(fGeom,this,primaryVertex);
    fTable->FillArray(fITSclusters);
    fTable->FillArrayCoorAngles(); 
  }

  Int_t ntracks = event->GetNumberOfTracks();
  while (ntracks--) {
    AliESDtrack *esd=event->GetTrack(ntracks);
    if ((esd->GetStatus()&
	 AliESDtrack::kITSin|AliESDtrack::kTPCin)==0) continue; 
    UInt_t idx[6];
    Int_t ncl = esd->GetITSclusters(idx);
    for(Int_t  clu=0; clu<ncl; clu++){
      AliITSclusterV2* cl = (AliITSclusterV2*)GetCluster(idx[clu]);
      if(cl->IsUsed()==0) cl->Use();
    }
  }
  
  Info("UseFoundTracksV2","Clusters of tracks prolonged from TPC deleted");
  

}

/*
//_______________________________________________________
Int_t AliITStrackerSA::SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, AliITStrackSA* trs,Double_t zvertex,Int_t pflag){
  //function used to to find the clusters associated to the track
  Int_t nc=0;
  AliITSlayer &lay = fgLayers[layer];
  Int_t * firstmod = new Int_t[fGeom->GetNlayers()];
  for(Int_t i=0;i<fGeom->GetNlayers();i++){
    firstmod[i]=fGeom->GetModuleIndex(i+1,1,1);
  }
  if(pflag==1){
      
    Float_t cx1,cx2,cy1,cy2;
    FindEquation(fPoint1[0],fPoint1[1],fPoint2[0],fPoint2[1],fPoint3[0],fPoint3[1],fCoef1,fCoef2,fCoef3);
    Int_t fun = FindIntersection(fCoef1,fCoef2,fCoef3,-(lay.GetR()*lay.GetR()),cx1,cy1,cx2,cy2);
    if(fun==0) {
      delete[] firstmod;
      return 0;
    }

    Double_t fi1 =TMath::ATan2(cy1,cx1);
    Double_t fi2 =TMath::ATan2(cy2,cx2);
    fPhiEstimate = ChoosePoint(fi1,fi2,fPhic);
  }

  Double_t zed = TMath::Tan(fLambdac)*lay.GetR()+zvertex;
  Double_t zed1  =  TMath::Tan(fLambdac+lambdawindow)*lay.GetR()+zvertex;
  Double_t zed2  =  TMath::Tan(fLambdac-lambdawindow)*lay.GetR()+zvertex;
  
  Double_t fi = fPhiEstimate;
  Int_t nmod  = lay.FindDetectorIndex(fi,zed);
  if (nmod < 0) {
    delete[] firstmod;
    return 0;
  }
  nmod += firstmod[layer];
 
  Int_t nm[8]={0,0,0,0,0,0,0,0};
  nm[0] = lay.FindDetectorIndex(fi+phiwindow,zed);
  nm[1] = lay.FindDetectorIndex(fi-phiwindow,zed);
  nm[2] = lay.FindDetectorIndex(fi,zed1);
  nm[3] = lay.FindDetectorIndex(fi,zed2);
  nm[4] = lay.FindDetectorIndex(fi+phiwindow,zed1);
  nm[5] = lay.FindDetectorIndex(fi-phiwindow,zed1);
  nm[6] = lay.FindDetectorIndex(fi+phiwindow,zed2);
  nm[7] = lay.FindDetectorIndex(fi-phiwindow,zed2);


  Int_t nn=0;
  TArrayI* array =(TArrayI*)fTable->GetListOfClusters(nmod);
  TArrayI* listc = new TArrayI(array->GetSize());
  for(Int_t i=0;i<array->GetSize();i++){
    Int_t in=(Int_t)array->At(i);
    listc->AddAt(in,nn);
    nn++;
  }
    
  Int_t k=0;
  Int_t val;
  while(k<8){
    for(Int_t h=k+1;h<8;h++){
      if(nm[k]>nm[h]){
	val=nm[k];
	nm[k]=nm[h];
	nm[h]=val;
      }
     
   }
    k++;
  }
 
  Int_t value=-5;
 
  for(Int_t ii=0;ii<8;ii++){
    if(nm[ii]!=value && nm[ii]!=nmod && nm[ii]>=0){
      TArrayI* ar =(TArrayI*)fTable->GetListOfClusters(nm[ii]+firstmod[layer]);
      listc->Set(listc->GetSize()+ar->GetSize());
      for(Int_t j=0;j<ar->GetSize();j++){
	Int_t in=(Int_t)ar->At(j);
	listc->AddAt(in,nn);
	nn++;
	value=nm[ii];
      }
    }
  }
 
  
  for(Int_t i=0;i<listc->GetSize();i++){
    Int_t index = (Int_t)listc->At(i);
    AliITSclusterV2* cllay = lay.GetCluster(index);
    if(cllay==0) continue;
    if(cllay->IsUsed()==1) continue;
    if(cllay->TestBit(kSAflag)==kTRUE) continue;
    Double_t phi   = fTable->GetPhiCluster(layer,index);
    Double_t lambda= fTable->GetLambdaCluster(layer,index);
 
    if(TMath::Abs(fLambdac-lambda)<lambdawindow && 
       TMath::Abs(fPhiEstimate-phi)<phiwindow){
      nc+=1;
      fLambdac = lambda;
      if(trs->GetNumberOfClustersSA()==15){
	delete[] firstmod;
        delete listc;
        return 0;
      }
      trs->AddClusterSA(layer,index);
      cllay->SetBit(kSAflag);
      fPhiEstimate=phi;
      fPointc[0]=fTable->GetXCluster(layer,index);
      fPointc[1]=fTable->GetYCluster(layer,index);
    }

  }
  delete listc;
  delete [] firstmod;
  return nc;

}
*/


//_______________________________________________________
Int_t AliITStrackerSA::SearchClusters(Int_t layer,Double_t phiwindow,Double_t lambdawindow, AliITStrackSA* trs,Double_t zvertex,Int_t pflag){
  //function used to to find the clusters associated to the track
  Int_t nc=0;
  AliITSlayer &lay = fgLayers[layer];
  Double_t r=lay.GetR(),tgl=TMath::Tan(fLambdac);

  if(pflag==1){      
    Float_t cx1,cx2,cy1,cy2;
    FindEquation(fPoint1[0],fPoint1[1],fPoint2[0],fPoint2[1],fPoint3[0],fPoint3[1],fCoef1,fCoef2,fCoef3);
    if (FindIntersection(fCoef1,fCoef2,fCoef3,-r*r,cx1,cy1,cx2,cy2)==0)
       return 0;
    Double_t fi1=TMath::ATan2(cy1,cx1);
    Double_t fi2=TMath::ATan2(cy2,cx2);
    fPhiEstimate=ChoosePoint(fi1,fi2,fPhic);
  }

  Double_t dz=r*lambdawindow*TMath::Sqrt(1+tgl*tgl) + 0.3*TMath::Abs(tgl);
  Double_t zmax=r*tgl + zvertex + dz;
  Double_t zmin=r*tgl + zvertex - dz;

  Int_t indmin = lay.FindClusterIndex(zmin);
  Int_t indmax = lay.FindClusterIndex(zmax);
  for (Int_t index=indmin; index<indmax; index++) {
     AliITSclusterV2 *c=lay.GetCluster(index);
     if (c->IsUsed()) continue;
     if (c->GetQ()<=0) continue;
     if (c->TestBit(kSAflag)==kTRUE) continue;

     Double_t phi   =fTable->GetPhiCluster(layer,index);
     if (TMath::Abs(phi-fPhiEstimate)>phiwindow) continue;

     Double_t lambda=fTable->GetLambdaCluster(layer,index);
     if (TMath::Abs(lambda-fLambdac)>lambdawindow) continue;

     if(trs->GetNumberOfClustersSA()==15) return 0;

     trs->AddClusterSA(layer,index);
     nc++;
     fLambdac=lambda;
     fPhiEstimate=phi;
     fPointc[0]=fTable->GetXCluster(layer,index);
     fPointc[1]=fTable->GetYCluster(layer,index);

     c->SetBit(kSAflag);
  }
  return nc;
}


//________________________________________________________________
void AliITStrackerSA::UpdatePoints(){
  //update of points for the estimation of the curvature  

  //fPoint1[0]=fPoint2[0]; 
  //fPoint1[1]=fPoint2[1];
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
AliITStrackV2* AliITStrackerSA::FindTrackLowChiSquare(TObjArray* tracklist, Int_t dim) const {
  // returns track with lowes chi square  
  if(dim==1){
    AliITStrackV2* trk = (AliITStrackV2*)tracklist->At(0);
    return trk;
  }
  if(dim==0) return 0;
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

  AliITStrackV2* trk = (AliITStrackV2*)tracklist->At(index[dim-1]);
  delete [] chi2;
  delete [] index;
  return trk;
  
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

