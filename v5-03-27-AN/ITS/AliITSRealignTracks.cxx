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

//---------------------------------
//Class to perform the realignment if the Inner Tracking System 
//with an iterative approach based on track to cluster residuals
// minimization. A chi2 function of the residuals is minimized with
//respect to alignment parameters. The class allows both single module
//realignment and set of modules realignment. Tracks are fitted with
//AliTrackFitter* fitters. AliTrackFitterKalman is more suited for
//straight line (e.g. cosmic in the absence of magnetic field) but can't
//work with helixes. AliTrackFitterRieman is suited for helixes. 
//The minimization is performed by AliTrackResiduals* classes: default
//one is AliTrackResidualsFast (analytic minimization by inversion). 
//For numerical minimization using MINUIT, use AliTrackResidualChi2.
//Methods are present to defined both the set of modules where the tracks
//are fittef and the set of modules which  are to be realigned
//The main method is  AlignVolumesITS
//
//Class by: A. Rossi, andrea,rossi@ts.infn.it

#include <TFile.h>
#include <TStopwatch.h>
#include <TNtuple.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "AliITSRealignTracks.h"
#include "AliAlignObjParams.h"
#include "AliAlignObj.h"
#include "AliGeomManager.h"
#include "AliTrackFitter.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackResidualsFast.h"
#include "AliTrackResidualsChi2.h"
#include "AliTrackResidualsLinear.h"

#include "AliLog.h"
#include <TSystem.h>
#include <TGeoManager.h>

class AliAlignmentTracks;
class TGeoMatrix;
class TArray;


/* $Id$ */


ClassImp(AliITSRealignTracks)

const Int_t kreferSect=2;


AliITSRealignTracks::AliITSRealignTracks(TString minimizer,Int_t fit,Bool_t covUsed,TString fileintro,TString geometryfile,TString misalignmentFile,TString startingfile):
  AliAlignmentTracks(),
  fSurveyObjs(0),
  fgeomfilename(),
  fmintracks(),
  fUpdateCov(kFALSE),
  fVarySigmaY(kFALSE),
  fCorrModules(0),
  fLimitCorr(0),
  fsigmaY(),
  fDraw(kFALSE),  
  fAlignDrawObjs(0), 
  fCanvPar(0), 
  fCanvGr(0), 
  fgrIterMeanX(0), 
  fgrIterRMSX(0),  
  fgrIterMeanY(0), 
  fgrIterRMSY(0),  
  fgrIterMeanZ(0), 
  fgrIterRMSZ(0),  
  fgrIterMeanPsi(0), 
  fgrIterRMSPsi(0),  
  fgrIterMeanTheta(0), 
  fgrIterRMSTheta(0),  
  fgrIterMeanPhi(0), 
  fgrIterRMSPhi(0)  
{

  // minimizer="fast"->AliTrackResidualFast minimizer
  //           "minuit"->AliTrackResidualChi2 minimizer
  //           "minuitnorot"->AliTrackResidualChi2 minimizer without rotations degrees of freedom
  //           "linear"->AliTrackResidualLinear minimizer    
  //fit=0-> Riemann Fitter, fit=1->Kalman
  //fileintro=file into which the Tree with the space points is stored
  //geometryfile=file containing the geometry  
  
  
  SetPointsFilename(fileintro.Data());
  SetGeomFilename(geometryfile);
  InitAlignObjs();
  if(!InitSurveyObjs(kFALSE))AliWarning("Unable to set Survey AlignObjs!");
  
  if(startingfile=="")printf("Starting from default geometry \n");
  else ReadAlignObjs(startingfile.Data(),"ITSAlignObjs");
  
  if(misalignmentFile=="")printf("NO FAKE MISALIGNMENT INTRODUCED \n");
  else {
    Bool_t misal=Misalign(misalignmentFile,"ITSAlignObjs");
    if(!misal)AliWarning("Incorrect fake misalignment filename!");;
  }
  
  if(!gGeoManager) AliGeomManager::LoadGeometry(fgeomfilename.Data());
  if(covUsed)SetCovIsUsed(kTRUE);
  if(!SelectFitter(fit))AliWarning("Incorrect fitter assignment!");
  if(!SelectMinimizer(minimizer))AliWarning("Incorrect minimizer assignment!");
  fsigmaY=1.;
  fmintracks=1;
  BuildIndex();
  
}

AliITSRealignTracks::AliITSRealignTracks(const AliITSRealignTracks &realignTracks):
  AliAlignmentTracks(),
  fSurveyObjs(new AliAlignObj**(*realignTracks.fSurveyObjs)),
  fgeomfilename(realignTracks.fgeomfilename),
  fmintracks(realignTracks.fmintracks),
  fUpdateCov(realignTracks.fUpdateCov),
  fVarySigmaY(realignTracks.fVarySigmaY),
  fCorrModules(new Double_t *(*realignTracks.fCorrModules)),
  fLimitCorr(realignTracks.fLimitCorr),
  fsigmaY(realignTracks.fsigmaY),
  fDraw(kFALSE),  
  fAlignDrawObjs(realignTracks.fAlignDrawObjs), 
  fCanvPar(realignTracks.fCanvPar), 
  fCanvGr(realignTracks.fCanvGr), 
  fgrIterMeanX(realignTracks.fgrIterMeanX), 
  fgrIterRMSX(realignTracks.fgrIterRMSX),  
  fgrIterMeanY(realignTracks.fgrIterMeanY), 
  fgrIterRMSY(realignTracks.fgrIterRMSY),  
  fgrIterMeanZ(realignTracks.fgrIterMeanZ), 
  fgrIterRMSZ(realignTracks.fgrIterRMSZ),  
  fgrIterMeanPsi(realignTracks.fgrIterMeanPsi), 
  fgrIterRMSPsi(realignTracks.fgrIterRMSPsi),  
  fgrIterMeanTheta(realignTracks.fgrIterMeanTheta), 
  fgrIterRMSTheta(realignTracks.fgrIterRMSTheta),  
  fgrIterMeanPhi(realignTracks.fgrIterMeanPhi), 
  fgrIterRMSPhi(realignTracks.fgrIterRMSPhi)  

{//Copy Constructor
  AliWarning("Can't copy AliAlignmentTracks Data member!");
}

AliITSRealignTracks& AliITSRealignTracks::operator=(const AliITSRealignTracks &obj){
  ////////////////////////
  // Assignment operator
  ////////////////////////
  this->~AliITSRealignTracks();
  new(this) AliITSRealignTracks(obj); 
  return *this;
}

AliITSRealignTracks::~AliITSRealignTracks(){
  //destructor
  
  if(fSurveyObjs)   DeleteSurveyObjs();
  if(fAlignDrawObjs) DeleteDrawHists();
  //delete [] fSurveyObjs; 
  
} 


//_____________________________
Bool_t AliITSRealignTracks::SelectFitter(Int_t fit,Int_t minTrackPoint){
  //Method to select the fitter: 0 for AliTrackFitterRieman (use this for helixes)
  //                             1 for AliTrackFitterKalman
  //minTrackPoint defines the minimum number of points (not rejected by the fit itself) 
  //a track should have to fit it
 
  if(fit==1){
     AliTrackFitterKalman *fitter= new AliTrackFitterKalman();
     fitter->SetMinNPoints(minTrackPoint);
     SetTrackFitter(fitter);    
  }
 
  else if(fit==0){
    AliTrackFitterRieman *fitter=new AliTrackFitterRieman();
    fitter->SetMinNPoints(minTrackPoint);
    SetTrackFitter(fitter);
  }
  else return kFALSE;
 
  return kTRUE;
}


Bool_t AliITSRealignTracks::SelectMinimizer(TString minimizer,Int_t minpoints,const Bool_t *coord){
  //Method to select the minimizer: "minuit" for AliTrackFitterChi2 (numerical minimization by MINUIT)
  //                                "fast" for AliTrackResidualsFast
  //                                "linear" for AliTrackResidualsLinear
  //                                "minuitnorot" for AliTrackFitterChi2 by 
  // coord[6] allows to fix the degrees of freedom in the minimization (e.g. look only for tranlsations, 
  //       or only for rotations). The coord are: dTx,dTy,dTz,dPsi,dTheta,dPhi where "d" stands for 
  //       "differential" and "T" for translation. If coord[i] is set to kTRUE then the i coord is fixed
  //       When a coordinate is fixed the value returnd for it is 0
  // minnpoints fix the minimum number of residuals to perform the minimization: it's not safe
  //     to align a module with a small number of track passing through it since the results could be
  //     not reliable. For single modules the number of residuals and the number of tracks passing 
  //     through it and accepted bu the fit procedure coincide. This is not the case for sets of modules,
  //     since a given track can pass through two or more modules in the set (e.g a cosmic track can give
  //     two cluster on a layer)

  AliTrackResiduals *res;
  if(minimizer=="minuit"){
    res = new AliTrackResidualsChi2();
    if(coord){
      for(Int_t j=0;j<6;j++){
	if(coord[j])res->FixParameter(j);
      }
    }
  }
  else if(minimizer=="minuitnorot"){
    res = new AliTrackResidualsChi2();
    res->FixParameter(3);
    res->FixParameter(4);
    res->FixParameter(5);
  }
  
  else if(minimizer=="fast"){
    res = new AliTrackResidualsFast();
    if(coord){
      for(Int_t j=0;j<6;j++){
	if(coord[j])res->FixParameter(j);
      }
    }
  } 
  else if(minimizer=="linear"){
    res = new AliTrackResidualsLinear();
  } 
  
  else {
    printf("Trying to set a non existing minimizer! \n");
    return kFALSE;
  }
  
  res->SetMinNPoints(minpoints);
  SetMinimizer(res);
  
  return kTRUE;
}

//____________________________________
void AliITSRealignTracks::SetVarySigmaY(Bool_t varysigmay,Double_t sigmaYfixed){
  //SigmaY is the value of the error along the track direction assigned 
  //to the AliTrackPoint constructed from the extrapolation of the fit of a track
  //to the module one is realigning. This error simulate the uncertainty on the
  //position of the cluster in space due to the misalingment of the module (i.e. 
  //you don't the real position and orientation of the plane of the desired extrapolation 
  //but you rely on the fit and consider the point to lie along the track)
  
  
  fVarySigmaY=varysigmay;
  if(!varysigmay){
    if(sigmaYfixed>0.)fsigmaY=sigmaYfixed;
    else {
      printf("Negative assignment to sigmaY! set it to default value (1 cm) \n");
      fsigmaY=1.;
    }
  }
}

//_______________________________________
void AliITSRealignTracks::RealignITSVolIndependent(Int_t iter1,Int_t iterations,Int_t minNtracks,Int_t layer,Int_t minTrackPoint){
  

  //iter1=#iterations inside AliAlignmentTracks::AliAlignVolumesITS method 
  //iterations=#iterations on the all procedure
  //layer=0->all ITS, otherways the usual notation is considered (1=SPD1,2=SPD2,3=SDD1,4=SDD2,5=SSD1,6=SSD2)
  //minNtracks=minimun number of tracks passing through a module in order to try to realign the module itsself
  //           if minNtracks<0, minimun number of tracks is |minNtracks|*minNumPoint[layer]/fact (see the code below): this allows a different
  //           choice of the number of tracks required on different layers and to vary these numbers once tuned the relative proportions.  
  //minTrackPoint=minimun number of "good" points required to a track (THE POINT ON THE MODULE THAT IS GOING TO BE REALIGNED 
  //IS NEVER CONSIDERED->max number can be required is 11 for cosmics tracks) for the track being considered in the minimization
  

  fTrackFitter->SetMinNPoints(minTrackPoint);
  TArrayI volIDs2(2200); 
  volIDs2.Reset(0);
  TArrayI volIDs(1);
  TString command;
  TArrayI volIDsFit;
  
  Int_t iLayer,iLayerToAlign;

  Int_t minNumPoint[6]={100,100,100,100,50,50}; 
  Double_t fact=10; 
  Int_t j=0;
  
  Int_t size=0;
  Int_t layerNum,modNum,lastVolid=0;
  TNtuple *ntVolumeAlign=new TNtuple("ntVolumeAlign","NTuple with volume tried to be realigned","layerNum:modNum:volumeIDnum"); 
 
  TStopwatch *timer=new TStopwatch();
  timer->Start(); 
  BuildIndex();
  j=0;
  UShort_t volid;
  
  for(Int_t iter=0;iter<iterations;iter++){
   
    //Starting Independent Modules Realignment
    for(iLayerToAlign=(Int_t)AliGeomManager::kSPD1;iLayerToAlign<=(Int_t)AliGeomManager::kSSD2;iLayerToAlign++){
      if(layer!=0&&iLayerToAlign!=layer)continue;
      j=0;
      size=0;
      for(Int_t k=(Int_t)AliGeomManager::kSPD1;k<=(Int_t)AliGeomManager::kSSD2;k++){
	size+=AliGeomManager::LayerSize(k);
	printf("size: %d \n",size);
      }
       
      for (Int_t iModule=0;iModule<AliGeomManager::LayerSize(iLayerToAlign);iModule++){      
	j=0;
	if(minNtracks<0){	
	  if(GetLastIndex(iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer,iModule)<minNumPoint[iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer]*(-1*minNtracks/fact))continue;	}	
	else if(GetLastIndex(iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer,iModule)<minNtracks)continue;
	
	UShort_t volidAl = AliGeomManager::LayerToVolUID(iLayerToAlign,iModule);
	
	volIDsFit.Reset(0);
	volIDsFit.Set(size-1);
	for (iLayer=AliGeomManager::kSPD1;iLayer<AliGeomManager::kTPC1;iLayer++){
	  for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){	    
	    volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	    if(AliGeomManager::LayerToVolUID(iLayer,iModule2)==volidAl)continue;
	    volIDsFit.AddAt(volid,j);
	    j++;
	  }
	}
	volIDs.AddAt((Int_t)volidAl,0);
	if(iter==iterations-1){
	  volIDs2.AddAt(volidAl,lastVolid);
	  lastVolid++;
	}
	volIDs2.AddAt(volidAl,lastVolid);
	AlignVolumesITS(&volIDs,&volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iter1);  
      }
    }
    
    
    if((iter+1)%5==0||iter==0||iter==1||iter==2||iter==3||iter==iterations-1){
      command="RealignObj";
      command+=(iter+1);
      command.Append(".root");
      WriteRealignObjArray(command.Data(),AliGeomManager::kSPD1,AliGeomManager::kSSD2);
    }
  }
  
  
  if(j==0){printf("j=0 \n");return;}
  for(Int_t k=0;k<volIDs2.GetSize();k++){
    if(volIDs2.At(k)==0)break;
    layerNum=AliGeomManager::VolUIDToLayer(volIDs2.At(k),modNum);
    ntVolumeAlign->Fill(layerNum,modNum,volIDs2.At(k));
  }
  printf("End of selecting modules cycle: %d modules selected \n",j);
  TFile *f=new TFile("RealignVolNt.root","RECREATE");
  f->cd();
  ntVolumeAlign->Write();
  f->Close();
  
  timer->Stop();
  timer->Print();
  return;
}


void AliITSRealignTracks::RealignITStracks(TString minimizer,Int_t fit=0,Int_t iter1=1,Int_t iterations=5,Int_t minNtracks=-10,Int_t layer=0,Int_t minTrackPoint=6,Bool_t covUsed=kFALSE,TString misalignmentFile="",TString startingfile="",Int_t doGlobal=1)
{
  
  // minimizer="fast"->AliTrackResidualFast minimizer
  //           "minuit"->AliTrackResidualChi2 minimizer
  //           "minuitnorot"->AliTrackResidualChi2 minimizer without rotations degrees of freedom
  //           "linear"->AliTrackResidualLinear minimizer    
  //fit=0-> Riemann Fitter, fit=1->Kalman
  //iter1=#iterations inside AliAlignmentTracks::AliAlignVolumesITS method 
  //iterations=#iterations on the all procedure
  //layer=0->all ITS, otherways the usual notation is considered (1=SPD1,2=SPD2,3=SDD1,4=SDD2,5=SSD1,6=SSD2)
  //minNtracks=minimun number of tracks passing through a module in order to try to realign the module itsself
  //           if minNtracks<0, minimun number of tracks is |minNtracks|*minNumPoint[layer]/fact (see the code below): this allows a different
  //           choice of the number of tracks required on different layers and to vary these numbers once tuned the relative proportions.  
  //minTrackPoint=minimun number of "good" points required to a track (THE POINT ON THE MODULE THAT IS GOING TO BE REALIGNED 
  //              IS NEVER CONSIDERED->max number that can be required is 11 for cosmics tracks) for the track being considered in the minimization
  //doGlobal : do global realignment, 0=no, 1= yes, 2=only global 


  TArrayI volIDs2(2200); 
  volIDs2.Reset(0);
  TArrayI volIDs(1);
  TString command;
  TArrayI volIDsFit;
  
  Int_t iLayer,iLayerToAlign;

  Int_t minNumPoint[6]={100,100,100,100,50,50}; 
  Double_t fact=10; 
  Int_t count=0;
 
  Int_t size=0;
  Int_t layerNum,modNum,lastVolid=0;
  TNtuple *ntVolumeAlign=new TNtuple("ntVolumeAlign","NTuple with volume tried to be realigned","layerNum:modNum:volumeIDnum"); 

 
  if(!SelectFitter(fit))AliWarning("Incorrect fitter assignment!");
  if(!SelectMinimizer(minimizer))AliWarning("Incorrect minimizer assignment!");
  if(misalignmentFile=="")printf("NO FAKE MISALIGNMENT INTRODUCED \n");
  else {
    Bool_t misal=Misalign(misalignmentFile,"ITSAlignObjs");
    if(!misal)return;
  }

 
  TStopwatch *timer=new TStopwatch();
  timer->Start(); 
  BuildIndex();
  count=0;
  UShort_t volid;

  if(startingfile=="")printf("Starting from default geometry \n");
  else {
    printf("Starting from AlignObjs file: %s",startingfile.Data());
    ReadAlignObjs(startingfile.Data(),"ITSAlignObjs");
  }
  
  for(Int_t iter=0;iter<iterations;iter++){
    if(covUsed)SetCovIsUsed(kTRUE);
    
   
    //START HIERARCHY REALIGNMENT

    if(layer==0&&(doGlobal==1||doGlobal==2)){
      for(Int_t siter=0;siter<5;siter++){
	fTrackFitter->SetMinNPoints(2);
	SetCovUpdate(kFALSE);
	AlignSPDHalfBarrelToSectorRef(kreferSect,3);
	//	AlignSPDBarrel(1);
	//	if(siter==0)SetCovUpdate(kFALSE);
	//	AlignSPDHalfBarrel(0,3);
	//	SetCovUpdate(kTRUE);
	AlignSPDHalfBarrelToHalfBarrel(1,3);
	//	AlignSPDHalfBarrelToSectorRef(kreferSect,3);
	for(Int_t sector=0;sector<10;sector++){
	  SetMinNtracks(100);
	  if(sector==kreferSect)continue;
	  AlignSPDSectorWithSectors(sector,1);
	}


	for(Int_t lay=1;lay<=6;lay++){
	  if(!AlignLayerToSector(lay,kreferSect,3))AlignLayerToSPDHalfBarrel(lay,0,3);
	}
	AlignSPDHalfBarrel(0,3);
	
	Int_t layers[6]={2,2,1,0,0,0};
	fTrackFitter->SetMinNPoints(4);
	AlignLayersToLayers(layers,1);
      
	fTrackFitter->SetMinNPoints(6);
	layers[2]=2;
	layers[3]=1;//{2,2,2,1,0,0};
	AlignLayersToLayers(layers,1);

	fTrackFitter->SetMinNPoints(6);
	layers[3]=2;
	layers[4]=1;//{2,2,2,2,1,0};
	AlignLayersToLayers(layers,1);

	fTrackFitter->SetMinNPoints(6);
	layers[4]=2;
	layers[5]=1;//{2,2,2,2,2,1};
	AlignLayersToLayers(layers,1);
	
	
	for(Int_t sector=0;sector<10;sector++){
	  AlignSPDSectorToOuterLayers(sector,1);
	}
	WriteRealignObjArray("AfterGlobal.root",AliGeomManager::kSPD1,AliGeomManager::kSSD2);
      }
    }
        
    if(doGlobal==2)return;    

    if(covUsed)SetCovUpdate(kTRUE);
    SetMinNtracks(1);


    // STARTS INDEPENDENT MOULES REALIGNMENT

    fTrackFitter->SetMinNPoints(minTrackPoint);
    for(iLayerToAlign=(Int_t)AliGeomManager::kSPD1;iLayerToAlign<=(Int_t)AliGeomManager::kSSD2;iLayerToAlign++){
      if(layer!=0&&iLayerToAlign!=layer)continue;
      count=0;
      size=0;
      for(Int_t k=(Int_t)AliGeomManager::kSPD1;k<=(Int_t)AliGeomManager::kSSD2;k++){
	size+=AliGeomManager::LayerSize(k);
	printf("size: %d \n",size);
      }
      
      for (Int_t iModule=0;iModule<AliGeomManager::LayerSize(iLayerToAlign);iModule++){      
	count=0;
	if(minNtracks<0){	
	  if(GetLastIndex(iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer,iModule)<minNumPoint[iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer]*(-1*minNtracks/fact))continue;
	}	
	else if(GetLastIndex(iLayerToAlign-(Int_t)AliGeomManager::kFirstLayer,iModule)<minNtracks)continue;
	
	UShort_t volidAl = AliGeomManager::LayerToVolUID(iLayerToAlign,iModule);
	
	volIDsFit.Reset(0);
	volIDsFit.Set(size-1);
	for (iLayer=(Int_t)AliGeomManager::kSPD1;iLayer<(Int_t)AliGeomManager::kTPC1;iLayer++){
	  for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){
	    
	  
	    volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	  
	    if(AliGeomManager::LayerToVolUID(iLayer,iModule2)==volidAl)continue;
	    volIDsFit.AddAt(volid,count);
	    count++;
	  }
	}
    
	volIDs.AddAt((Int_t)volidAl,0);
	if(iter==iterations-1){
	  volIDs2.AddAt(volidAl,lastVolid);
	  lastVolid++;
	}
	volIDs2.AddAt(volidAl,lastVolid);	
	AlignVolumesITS(&volIDs,&volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iter1);  	
      }
    }   
    
    if((iter+1)%2==0||(iter+1)%5==0||iter==0||iter==1||iter==2||iter==3||iter==iterations-1){
      command="RealignObj";
      command+=(iter+1);
      command.Append(".root");
      WriteRealignObjArray(command.Data(),AliGeomManager::kSPD1,AliGeomManager::kSSD2);
    }
    
  }
  
  if(count==0){printf("count=0 \n");return;}
  for(Int_t k=0;k<volIDs2.GetSize();k++){
    if(volIDs2.At(k)==0)break;
    layerNum=AliGeomManager::VolUIDToLayer(volIDs2.At(k),modNum);
    ntVolumeAlign->Fill(layerNum,modNum,volIDs2.At(k));
  }
  printf("End of selecting modules cycle: %d modules selected \n",count);
  TFile *f=new TFile("RealignVolNt.root","RECREATE");
  f->cd();
  ntVolumeAlign->Write();
  f->Close();
  
  timer->Stop();
  timer->Print();
  return;
  
}


//______________________________________________________________________________
void AliITSRealignTracks::InitAlignObjs()
{
  // Initialize the alignment objects array
  TMatrixDSym c(6);
  Double_t cov[21];
  for(Int_t i=0;i<21;i++)cov[i]=0.;
  for(Int_t i=0;i<3;i++)cov[i*(i+1)/2+i]=0.05*0.05;//Set Default Error to 500 micron for Translations 
  for(Int_t i=3;i<6;i++)cov[i*(i+1)/2+i]=0.001*0.001*180*180/3.14/3.14;//and 1 mrad for rotations (global ref. sysytem->~40 micron for SPD1,~450 micron for SSD2)

  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fAlignObjs = new AliAlignObj**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fAlignObjs[iLayer] = new AliAlignObj*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer+ AliGeomManager::kFirstLayer,iModule);
      fAlignObjs[iLayer][iModule] = new AliAlignObjParams(AliGeomManager::SymName(volid),volid,0,0,0,0,0,0,kTRUE);
      fAlignObjs[iLayer][iModule]->SetCorrMatrix(cov);
      fAlignObjs[iLayer][iModule]->SetUniqueID(0);
    }
  }
}

//______________________________________________________________________________
void AliITSRealignTracks::ResetCorrModules(){
  // Initialize and reset to 0 the array with the information on correlation
  if(!fCorrModules){
    Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
    fCorrModules = new Double_t*[nLayers];
    for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
      fCorrModules[iLayer] = new Double_t[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
	fCorrModules[iLayer][iModule]=0.;
      }
    }
  }
  else{
    for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
	fCorrModules[iLayer][iModule]=0.;
      }
    }
  }
}

//______________________________________________________________________________
Bool_t AliITSRealignTracks::InitSurveyObjs(Bool_t infinite,Double_t factor,TString filename,TString arrayName){
  //Initialize the Survey Objects. There is the possibility to set them equal to external objects
  //   stored in file filename. Otherwuse they are set equal to 0 and with default values for the variances
  //infinite: set the cov matrix to extremly large values, so that the results of a minimization
  //   are never rejected by the comparison with the survey
  //factor: multiplication factor for the variances of the cov. matrix of the survey obj.

  if(fSurveyObjs)DeleteSurveyObjs();
  Bool_t fromfile=kFALSE;
  TFile *surveyObj;
  TClonesArray *clnarray;
  if(!filename.IsNull()){
    //Initialize from file
    if(gSystem->AccessPathName(filename.Data(),kFileExists)){
      printf("Wrong Survey AlignObjs File Name \n");
      return kFALSE;
    } 
    if(arrayName.IsNull()){
      printf("Null Survey Object Name! \n");
      return kFALSE;
    }
   
    fromfile=kTRUE;
  }

  // Initialize the alignment objects array with default values
  Double_t v=1.*factor;
  if(infinite)v*=100000.;
  TMatrixDSym c(6);
  Double_t cov[21];
  for(Int_t i=0;i<21;i++)cov[i]=0.;
  for(Int_t i=0;i<3;i++)cov[i*(i+1)/2+i]=0.1*0.1*v;//Set Default Error to 1 mm  for Translation 
  for(Int_t i=3;i<6;i++)cov[i*(i+1)/2+i]=0.03*0.03*180.*180./3.14/3.14*v;//and 30 mrad (~1.7 degrees)for rotations (global ref. sysytem)
  
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fSurveyObjs = new AliAlignObj**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fSurveyObjs[iLayer] = new AliAlignObj*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      fSurveyObjs[iLayer][iModule] = 0x0;
    }
  }
  
  if(fromfile){
    surveyObj=TFile::Open(filename.Data());
    if (!surveyObj || !surveyObj->IsOpen()) {
      AliError(Form("Could not open SurveyObjs file: %s !",filename.Data()));
      return kFALSE;
    }
    printf("Getting TClonesArray \n");
    clnarray=(TClonesArray*)surveyObj->Get(arrayName);
    Int_t size=clnarray->GetSize();
    UShort_t volid;
    for(Int_t ivol=0;ivol<size;ivol++){
      AliAlignObjParams *a=(AliAlignObjParams*)clnarray->At(ivol);
      volid=a->GetVolUID();
      Int_t iModule;
      AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
      if(iLayer<=0)continue;
      if(a->GetUniqueID()==0)continue;
      printf("Updating survey for volume: %d ,layer: %d module: %d from file\n",volid,iLayer,iModule);
      fSurveyObjs[iLayer-AliGeomManager::kFirstLayer][iModule] = new AliAlignObjParams(*a);
      fSurveyObjs[iLayer-AliGeomManager::kFirstLayer][iModule]->SetUniqueID(a->GetUniqueID());
      fSurveyObjs[iLayer-AliGeomManager::kFirstLayer][iModule]->Print("");
    }
    delete clnarray;
    surveyObj->Close();
  }
 
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer+ AliGeomManager::kFirstLayer,iModule);
      if(!fSurveyObjs[iLayer][iModule]){
	printf("Updating survey for volume: %d ,layer: %d module: %d with default values \n",volid,iLayer,iModule);
	fSurveyObjs[iLayer][iModule] = new AliAlignObjParams(AliGeomManager::SymName(volid),volid,0,0,0,0,0,0,kTRUE);
	fSurveyObjs[iLayer][iModule]->SetCorrMatrix(cov);
	fSurveyObjs[iLayer][iModule]->SetUniqueID(0);
      }
      
    }
  }
  
 
  return kTRUE;
}


//______________________________________________________________________________
Int_t AliITSRealignTracks::CheckWithSurvey(Double_t factor,const TArrayI *volids){
  
  // Check the parameters of the alignment objects in volids (or of all objects if volids is null) 
  // are into the boundaries set by the cov. matrix of the survey objs
  // Returns the number of objects out of boudaries
  AliAlignObj *alignObj;	
  Int_t outofsurv=0;
  UShort_t volid;
  Double_t surveycov[21],transl[3],rot[3],survtransl[3],survrot[3];
  if(volids==0x0){
    for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++){
	volid=AliGeomManager::LayerToVolUIDSafe(iLayer+AliGeomManager::kFirstLayer,iModule);
	alignObj=GetAlignObj(volid);
	alignObj->GetPars(transl,rot);
	fSurveyObjs[iLayer][iModule]->GetCovMatrix(surveycov);
	fSurveyObjs[iLayer][iModule]->GetPars(survtransl,survrot);
	if(TMath::Sqrt(TMath::Abs(surveycov[0]))*factor<TMath::Abs(transl[0]-survtransl[0])||TMath::Sqrt(TMath::Abs(surveycov[2]))*factor<TMath::Abs(transl[1]-survtransl[1])||TMath::Sqrt(TMath::Abs(surveycov[5]))*factor<TMath::Abs(transl[2]-survtransl[2])||TMath::Sqrt(TMath::Abs(surveycov[9]))*factor<TMath::Abs(rot[0]-survrot[0])||TMath::Sqrt(TMath::Abs(surveycov[14]))*factor<TMath::Abs(rot[1]-survrot[1])||TMath::Sqrt(TMath::Abs(surveycov[20]))*factor<TMath::Abs(rot[2]-survrot[2])){
	  printf("Results for module %d out of Survey: reinitializing it from survey \n",volid);
	  //	  *alignObj = *alignObjSurv;
	  alignObj->SetPars(survtransl[0],survtransl[1],survtransl[2],survrot[0],survrot[1],survrot[2]);
	  alignObj->SetUniqueID(0);
	  if(fUpdateCov)alignObj->SetCorrMatrix(surveycov);
	  outofsurv++;
	}
      }
    }
  }
  else{
    Int_t iLayer;
    Int_t iModule;
    for(Int_t j=0;j<volids->GetSize();j++){
      volid=volids->At(j);
      alignObj=GetAlignObj(volid);
      alignObj->GetPars(transl,rot);
      iLayer=(Int_t)AliGeomManager::VolUIDToLayerSafe(volid,iModule)-(Int_t)AliGeomManager::kFirstLayer;
      fSurveyObjs[iLayer][iModule]->GetCovMatrix(surveycov);
      fSurveyObjs[iLayer][iModule]->GetPars(survtransl,survrot);
      if(TMath::Sqrt(TMath::Abs(surveycov[0]))*factor<TMath::Abs(transl[0]-survtransl[0])||TMath::Sqrt(TMath::Abs(surveycov[2]))*factor<TMath::Abs(transl[1]-survtransl[1])||TMath::Sqrt(TMath::Abs(surveycov[5]))*factor<TMath::Abs(transl[2]-survtransl[2])||TMath::Sqrt(TMath::Abs(surveycov[9]))*factor<TMath::Abs(rot[0]-survrot[0])||TMath::Sqrt(TMath::Abs(surveycov[14]))*factor<TMath::Abs(rot[1]-survrot[1])||TMath::Sqrt(TMath::Abs(surveycov[20]))*factor<TMath::Abs(rot[2]-survrot[2])){
	printf("Results for module %d out of Survey: reinitializing it from survey \n",volid);
	//	  *alignObj = *alignObjSurv;
	alignObj->SetPars(survtransl[0],survtransl[1],survtransl[2],survrot[0],survrot[1],survrot[2]);
	alignObj->SetUniqueID(0);
	if(fUpdateCov)alignObj->SetCorrMatrix(surveycov);
	outofsurv++;
      }
    }
  }  
  return outofsurv;
}  

//___________________________________________________________________

void AliITSRealignTracks::ResetAlignObjs(Bool_t all,TArrayI *volids)
{
  // Reset the alignment objects in volids or all if all=kTRUE
  if(all){
    for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++)
	fAlignObjs[iLayer][iModule]->SetPars(0,0,0,0,0,0);
    }
  }
  else{
    Int_t layer;
    Int_t mod;
    for(Int_t j=0;j<volids->GetSize();j++){
      layer=(Int_t)AliGeomManager::VolUIDToLayer(volids->At(j),mod)-(Int_t)AliGeomManager::kFirstLayer;
      fAlignObjs[layer][mod]->SetPars(0,0,0,0,0,0);
    }
  }
}

//______________________________________________-
void AliITSRealignTracks::DeleteSurveyObjs()
{//destructor for the survey objs. array

  if(!fSurveyObjs)return;
  // Delete the alignment objects array
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++){
      if (fSurveyObjs[iLayer][iModule])	delete fSurveyObjs[iLayer][iModule];
    }
    
    if(fSurveyObjs[iLayer])delete [] fSurveyObjs[iLayer];
  }
    
  delete [] fSurveyObjs;
  fSurveyObjs = 0;
}


//______________________________________________________________________________
Bool_t AliITSRealignTracks::ReadAlignObjs(const char *alignObjFileName, const char* arrayName){

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

//_________________________________________
Bool_t AliITSRealignTracks::FirstAlignmentLayers(const Bool_t *layers,Int_t minNtracks,Int_t iterations,Bool_t fitall,const TArrayI *volidsSet){

  //Align all modules in the set of layers independently according to a sequence based on the number of tracks passing through a given module
  
  BuildIndex();
  TString name="DrawFirstAlignment_Layers";
  UShort_t voluid;
  Int_t **lastIndex;
  Int_t laymax = 0;
  Int_t modmax = 0;
  Int_t maxntr=0,nMod=0,modAligned=0,size=0;
  for(Int_t i=0;i<6;i++){
    if(layers[i]==1){
      size+=AliGeomManager::LayerSize(i+AliGeomManager::kFirstLayer);
      name+=i+1;
    }
  }

  TArrayI *volFit=new TArrayI(size);
  TArrayI *volFit2=new TArrayI(size-1);
  TArrayI *sequence=new TArrayI(size);
  TArrayI *volIn=new TArrayI(1);
  
  // Initialize the index arrays
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  lastIndex = new Int_t*[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    lastIndex[iLayer] = new Int_t[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      lastIndex[iLayer][iModule] =  fLastIndex[iLayer][iModule];
      if(iLayer<=(AliGeomManager::kSSD2-AliGeomManager::kFirstLayer)&&layers[iLayer]==1){
	volFit->AddAt(AliGeomManager::LayerToVolUID(iLayer+AliGeomManager::kFirstLayer,iModule),maxntr);
	maxntr++;
      }
    }
  }
  Int_t found=0;
  maxntr=minNtracks+1;
  while (maxntr>minNtracks){
    maxntr=minNtracks;
    for (Int_t iLayer = 0; iLayer <= (AliGeomManager::kSSD2 - AliGeomManager::kFirstLayer); iLayer++) {
      if(layers[iLayer]==0)continue;
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
	if(lastIndex[iLayer][iModule]>maxntr){
	  maxntr=lastIndex[iLayer][iModule];
	  laymax=iLayer;
	  modmax=iModule;
	}
      }
    }
    if(maxntr>minNtracks){
      voluid=AliGeomManager::LayerToVolUID(laymax+AliGeomManager::kFirstLayer,modmax);
      sequence->AddAt(voluid,nMod);
      lastIndex[laymax][modmax]=0;
      nMod++;
    }
  }

  sequence->Set(nMod);

 Int_t ilayer,imod;
  for(Int_t iter=0;iter<iterations;iter++){
    if(iter>0&&fDraw)UpdateDraw(sequence,iter,iter);
    modAligned=0;
    for(Int_t k=0;k<nMod;k++){
      TArrayI *volFit3;
      voluid=sequence->At(k);
      ilayer=AliGeomManager::VolUIDToLayer(voluid,imod);
      volIn->AddAt(voluid,0);
      found=0;
      if(!fitall){
	for(Int_t j=0;j<nMod;j++){
	  if(j==k){
	    found=1;
	    continue;
	  }
	  else volFit2->AddAt(sequence->At(j),j-found);
	}
	volFit2->Set(nMod-1);
      }
      else{
	for(Int_t j=0;j<volFit->GetSize();j++){
	  if(volFit->At(j)!=volIn->At(0))volFit2->AddAt(volFit->At(j),j-found);
	  else found=1;
	}
      }
      
      if(volidsSet){
	volFit3=IntersectVolArray(volidsSet,volFit2);
      }
      else volFit3=new TArrayI(*volFit2);
      
      
      if(AlignVolumesITS(volIn,volFit3,AliGeomManager::kSPD1,AliGeomManager::kTPC1,2))modAligned++;
      delete volFit3;
    
    }
  }
  Int_t noutofsurv=CheckWithSurvey(2.,sequence);
  printf("%d modules into the sequence \n %d modules re-aligned \n %d modules moved far away from survey (-> reset) \n",nMod,modAligned,noutofsurv);
  name.Append("_iter");
  name+=iterations;
  name.Append(".root");
  if(fDraw)WriteHists(name.Data());
  delete volFit;
  delete volFit2;
  delete sequence;
  for(Int_t m=0;m<nLayers;m++){
    delete [] lastIndex[m];
  }
  delete [] lastIndex;
  return kTRUE;
  
}

//__________________________________________
Bool_t AliITSRealignTracks::FirstAlignmentSPD(Int_t minNtracks,Int_t iterations,Bool_t fitall,const TArrayI *volidsSet){

  //OBSOLETE METHOD: perform a stand-alone realignment of the SPD modules
  //                 based on a sequence constructed accordingly to the number of tracks
  //                 passing through each module
  
  BuildIndex();
   
  UShort_t voluid;
  Int_t **lastIndex;
  Int_t laymax = 0;
  Int_t modmax = 0;
  Int_t maxntr=0,nMod=0,modAligned=0;
  TArrayI *volFit=new TArrayI(AliGeomManager::LayerSize(1)+AliGeomManager::LayerSize(2));
  TArrayI *volFit2=new TArrayI(AliGeomManager::LayerSize(1)+AliGeomManager::LayerSize(2)-1);
  TArrayI *sequence=new TArrayI(AliGeomManager::LayerSize(1)+AliGeomManager::LayerSize(2));
  TArrayI *volIn=new TArrayI(1);
 
  // Initialize the index arrays
  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  lastIndex = new Int_t*[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    lastIndex[iLayer] = new Int_t[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      lastIndex[iLayer][iModule] =  fLastIndex[iLayer][iModule];
      if(iLayer+AliGeomManager::kFirstLayer<=AliGeomManager::kSPD2){
	volFit->AddAt(AliGeomManager::LayerToVolUID(iLayer+AliGeomManager::kFirstLayer,iModule),maxntr);
	maxntr++;
      }
    }
  }
  Int_t found=0;
  maxntr=minNtracks+1;
  while (maxntr>minNtracks){
    maxntr=minNtracks;
    for (Int_t iLayer = 0; iLayer <= (AliGeomManager::kSPD2 - AliGeomManager::kFirstLayer); iLayer++) {
      for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
	if(lastIndex[iLayer][iModule]>maxntr){
	  laymax=iLayer;
	  modmax=iModule;
	  maxntr=lastIndex[iLayer][iModule];
	}
      }
    }
    if(maxntr>minNtracks){
      voluid=AliGeomManager::LayerToVolUID(laymax+AliGeomManager::kFirstLayer,modmax);
      sequence->AddAt(voluid,nMod);
      lastIndex[laymax][modmax]=0;
      nMod++;
      volIn->AddAt(voluid,0);
    }
  }
  sequence->Set(nMod);
  
  Int_t ilayer,imod;
  for(Int_t iter=0;iter<iterations;iter++){ 
    modAligned=0;
    for(Int_t k=0;k<nMod;k++){
      TArrayI *volFit3;
      voluid=sequence->At(k);
      ilayer=AliGeomManager::VolUIDToLayer(voluid,imod);
      volIn->AddAt(voluid,0);
      found=0;
      if(!fitall){
	for(Int_t j=0;j<nMod;j++){
	  if(j==k){
	    found=1;
	    continue;
	  }
	  else volFit2->AddAt(sequence->At(j),j-found);
	}
	volFit2->Set(nMod-1);
      }
      else{
	for(Int_t j=0;j<volFit->GetSize();j++){
	  if(volFit->At(j)!=volIn->At(0))volFit2->AddAt(volFit->At(j),j-found);
	  else found=1;
	}
      }
      
      if(volidsSet){
	volFit3=IntersectVolArray(volidsSet,volFit2);
      }
      else volFit3=new TArrayI(*volFit2);
      
      
      if(AlignVolumesITS(volIn,volFit3,AliGeomManager::kSPD1,AliGeomManager::kSDD1,2))modAligned++;
      delete volFit3;
      //      if(volidsSet)delete volFit3;
    }
  }
  Int_t noutofsurv=CheckWithSurvey(2.,sequence);
  printf("%d modules into the sequence \n %d modules re-aligned \n %d modules moved far away from survey (-> reset) \n",nMod,modAligned,noutofsurv);
  delete volFit;
  delete volFit2;
  delete sequence;
  for(Int_t m=0;m<nLayers;m++){
    delete [] lastIndex[m];
  }
  delete [] lastIndex;

  return kTRUE;
}


//__________________________________
Bool_t AliITSRealignTracks::SPDmodulesAlignToSSD(Int_t minNtracks,Int_t iterations){
  //Align each SPD module with at least minNtracks passing through it with respect to SSD
  //The selection based on the minimum number of tracks is a fast one:
  // the number considere here doesn't coincide with the tracks effectively used then in the
  // minimization, it's just the total number of tracks in the sample passing through the module
  // The procedure is iterated "iterations" times
  Int_t volSSD[6]={0,0,0,0,1,1};
  TArrayI *volOuter=GetLayersVolUID(volSSD);
  TArrayI *voluid=new TArrayI(1);
  for(Int_t iter=0;iter<iterations;iter++){
    //SPD1
    for(Int_t imod=0;imod<AliGeomManager::LayerSize(AliGeomManager::kSPD1);imod++){    if(GetLastIndex(AliGeomManager::kSPD1-AliGeomManager::kFirstLayer,imod)<minNtracks){
      printf("Not enough tracks for module: lay %d mod %d \n",1,imod );
      continue;
    }
    voluid->AddAt(AliGeomManager::LayerToVolUID(AliGeomManager::kSPD1,imod),0);
    AlignVolumesITS(voluid,volOuter,AliGeomManager::kSSD1,AliGeomManager::kSSD2,2);  
    }
    //SPD2
    for(Int_t imod=0;imod<AliGeomManager::LayerSize(AliGeomManager::kSPD2);imod++){ 
      if(GetLastIndex(AliGeomManager::kSPD2-AliGeomManager::kFirstLayer,imod)<minNtracks){
	printf("Not enough tracks for module: lay %d mod %d \n",2,imod );
	continue;
      }
      voluid->AddAt(AliGeomManager::LayerToVolUID(AliGeomManager::kSPD2,imod),0);
      AlignVolumesITS(voluid,volOuter,AliGeomManager::kSSD1,AliGeomManager::kSSD2,2);  
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSRealignTracks::AlignVolumesITS(const TArrayI *volids, const TArrayI *volidsfit,
				     AliGeomManager::ELayerID layerRangeMin,
				     AliGeomManager::ELayerID layerRangeMax,
				     Int_t iterations){
  
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
  Bool_t correlated=kFALSE;
  Double_t surveycov[21],transl[3],rot[3],survtransl[3],survrot[3];
  Double_t frac;

  TGeoHMatrix hM;
  Double_t phiglob,smearing,rotorig[9],normplanevect[3]={0.,0.,0.},normplanevect2[3]={0.,0.,0.};  
  Double_t *deltarot;
  TMatrixDSym covmatrx(6);
  AliAlignObj *modAlign;

  // Load only the tracks with at least one
  // space point in the set of volume (volids)
  BuildIndex();
  AliTrackPointArray **points;
  Bool_t failed=kFALSE;
  Int_t pointsdim=0,skipped=0;
  // Start the iterations
  while (iterations > 0){
    normplanevect2[0]=0.;
    normplanevect2[1]=0.;
    normplanevect2[2]=0.;
    if(fLimitCorr>0.){
      ResetCorrModules();
      skipped=0;
    }
    Int_t nArrays = LoadPoints(volids, points,pointsdim);
    
    if (nArrays < fmintracks||nArrays<=0){
      failed=kTRUE;
      printf("Not enough tracks to try minimization: volUID %d and following in volids \n", volids->At(0));
      UnloadPoints(pointsdim, points);
      break;
    }
    frac=1./(Double_t)nArrays;
    AliTrackResiduals *minimizer = CreateMinimizer();
    minimizer->SetNTracks(nArrays);
    minimizer->InitAlignObj();
    AliTrackFitter *fitter = CreateFitter();

    //Here prepare to set the plane for GetPCArot
                                                       //    if(volids->GetSize()==1){//TEMPORARY: to be improved
    AliGeomManager::GetOrigRotation(volids->At(0),rotorig);
    if((Int_t)AliGeomManager::VolUIDToLayer(volids->At(0))==1){//TEMPORARY: to be improved  
      normplanevect[0]=-rotorig[1];
      normplanevect[1]=-rotorig[4];
      normplanevect[2]=0.;
    }
    else{
	  normplanevect[0]=rotorig[1];
	  normplanevect[1]=rotorig[4];
	  normplanevect[2]=0.;
    }
    
    phiglob=TMath::ATan2(normplanevect[1],normplanevect[0]);
    
    modAlign=GetAlignObj(volids->At(0));
    modAlign->GetMatrix(hM);
    deltarot=hM.GetRotationMatrix();
    for(Int_t j=0;j<3;j++){
      for(Int_t i=0;i<3;i++){
	normplanevect2[j]+=deltarot[j*3+i]*normplanevect[i];
      }
      // printf("Here the difference: norm1[%d]=%f  norm2[%d]=%f \n",j,normplanevect[j],j,normplanevect2[j]);
    }
    
    if(fVarySigmaY){
      if(modAlign->GetUniqueID()==0)smearing=fsigmaY;
      else{
	modAlign->GetCovMatrix(covmatrx);
	smearing=5.*5.*(covmatrx(0,0)+covmatrx(1,1)+covmatrx(2,2)+10.*10.*covmatrx(3,3)+10.*10.*covmatrx(4,4)+10.*10.*covmatrx(5,5))/6.; 
	//This is a sort of average: the trace with the variances of the angles 
	//weighted with 10 cm divided per 6 and the result multiplied per 25 
	// (the sqrt would be 5 times a sort of "mean sigma" )
	//	 
	  }
    }
    else smearing=fsigmaY;
    printf("This is the sigmaY value: %f \n",smearing);
    // the plane will be set into the loop on tracks
    

    for (Int_t iArray = 0; iArray < nArrays; iArray++) {
      if (!points[iArray]) continue;
      points[iArray]->Sort(kTRUE);
      fitter->SetTrackPointArray(points[iArray], kFALSE);
      // printf("Here normplane vect: %f \n",normplanevect2[1]); //TO BE REPLACED BY      fitter->SetNormPlaneVect(normplanevect2);
      if (fitter->Fit(volids,volidsfit,layerRangeMin,layerRangeMax) == kFALSE) continue;
      
       if(fLimitCorr>0.){
	correlated=kFALSE;
	AliTrackPoint p;
	Int_t layer,module;
	TArrayI *volparray=new TArrayI(points[iArray]->GetNPoints());
	for(Int_t point=0;point<points[iArray]->GetNPoints();point++){
	  points[iArray]->GetPoint(p,point);
	  volparray->AddAt(p.GetVolumeID(),point);
	}
	TArrayI	*volpArray=ExcludeVolidsFromVolidsArray(volids,volparray);
	for(Int_t point=0;point<volpArray->GetSize();point++){
	  layer=(Int_t)AliGeomManager::VolUIDToLayerSafe(volpArray->At(point),module);  
	  if(fCorrModules[layer-AliGeomManager::kFirstLayer][module]>fLimitCorr){
	    correlated=kTRUE;
	    //	    printf("volid %d, iarray = %d : skipping %d for Volume: %d \n",volids->At(0),iArray,skipped,volpArray->At(point));
	    skipped++;
	    break;
	  }
	}
	if(!correlated){
	  for(Int_t point=0;point<volpArray->GetSize();point++){
	    layer=(Int_t)AliGeomManager::VolUIDToLayerSafe(volpArray->At(point),module);  
	    //printf("Number of common tracks: %d \n",fCorrModules[layer-AliGeomManager::kFirstLayer][module]);
	    fCorrModules[layer-AliGeomManager::kFirstLayer][module]+=frac;
	    delete volparray;
	    delete volpArray;
	  }
	}
       	else { 
	  delete volparray;
	  delete volpArray;
	  continue;
	}
       }
       
      AliTrackPointArray *pVolId,*pTrack;
      fitter->GetTrackResiduals(pVolId,pTrack);
      minimizer->AddTrackPointArrays(pVolId,pTrack);
    }
    
    printf("Number of tracks considered: %d \n",nArrays);
    frac=(Double_t)skipped/(Double_t)nArrays;
    printf("Number of tracks skipped cause of correlation: %d (fraction: %f )\n",skipped,frac);
    
    Int_t ntracks=minimizer->GetNFilledTracks();
    frac=(Double_t)ntracks/(Double_t)nArrays;
    printf("Number of tracks into the minimizer: %d (fraction: %f )\n",ntracks,frac);
    if(ntracks<=fmintracks){
      printf("Not enough good tracks found: could not find parameter for volume %d (and following in volids)\n",volids->At(0));
      UnloadPoints(pointsdim, points);
      failed=kTRUE;
      break;
    }
    
    failed=(!minimizer->Minimize());
    
    // Update the alignment object(s)
    if (fDoUpdate) for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
      UShort_t volid = (*volids)[iVolId];
      if(!failed){
	Int_t iModule;
	AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
	AliAlignObj *alignObj = fAlignObjs[iLayer-AliGeomManager::kFirstLayer][iModule];  

	//Check the last minimization is not too large
	minimizer->GetAlignObj()->GetPars(transl,rot);     
	fSurveyObjs[iLayer-AliGeomManager::kFirstLayer][iModule]->GetCovMatrix(surveycov);
	fSurveyObjs[iLayer-AliGeomManager::kFirstLayer][iModule]->GetPars(survtransl,survrot);
	if(TMath::Sqrt(TMath::Abs(surveycov[0]))*2<TMath::Abs(transl[0]-survtransl[0])||TMath::Sqrt(TMath::Abs(surveycov[2]))*2<TMath::Abs(transl[1]-survtransl[1])||TMath::Sqrt(TMath::Abs(surveycov[5]))*2<TMath::Abs(transl[2]-survtransl[2])||TMath::Sqrt(TMath::Abs(surveycov[9]))*2<TMath::Abs(rot[0]-survrot[0])||TMath::Sqrt(TMath::Abs(surveycov[14]))*2<TMath::Abs(rot[1]-survrot[1])||TMath::Sqrt(TMath::Abs(surveycov[20]))*2<TMath::Abs(rot[2]-survrot[2])){
	  printf("Results for module %d too large: can't update them \n",volid);
	  alignObj->SetUniqueID(2);
	  if(iterations==1){
	    failed=kTRUE;
	  }
	}
	else{
	  if(fUpdateCov){
	    *alignObj *= *minimizer->GetAlignObj();
	    alignObj->SetUniqueID(1);
	  }
	  else{
	    alignObj->GetCovMatrix(covmatrx);
	    *alignObj *= *minimizer->GetAlignObj();
	    alignObj->SetCorrMatrix(covmatrx);
	    alignObj->SetUniqueID(1);
	  }
	  
	  /*alignObj->GetPars(transl,rot);
	  
	  if(TMath::Sqrt(TMath::Abs(surveycov[0]))*20<TMath::Abs(transl[0])||TMath::Sqrt(TMath::Abs(surveycov[2]))*20<TMath::Abs(transl[1])||TMath::Sqrt(TMath::Abs(surveycov[5]))*20<TMath::Abs(transl[2])||TMath::Sqrt(TMath::Abs(surveycov[9]))*20<TMath::Abs(rot[0])||TMath::Sqrt(TMath::Abs(surveycov[14]))*20<TMath::Abs(rot[1])||TMath::Sqrt(TMath::Abs(surveycov[20]))*20<TMath::Abs(rot[2])){
	  printf("Results for module %d out of Survey: reinitializing it from survey \n",volid);
	    //	  *alignObj = *alignObjSurv;
	    alignObj->SetPars(0.,0.,0.,0.,0.,0.);
	    alignObj->SetUniqueID(0);
	    if(fUpdateCov)alignObj->SetCorrMatrix(surveycov);
	    if(iterations==1){
	    failed=kTRUE;
	    }
	    }*/
	}
	if(iterations==1)alignObj->Print("");
      }
      else {
	printf("Minimization failed: cannot update AlignObj for volume: %d \n",volid);
      }
    }
    UnloadPoints(pointsdim,points);
    if(failed)break;
    minimizer->InitAlignObj();
    iterations--;
  }
  
  printf("\n \n");

  return (!failed);
}



//______________________________________________
Bool_t AliITSRealignTracks::AlignSPDBarrel(Int_t iterations){
  //Align the SPD barrel "iterations" times
  
  Int_t size=0,size2=0;
  Int_t layers[6]={1,1,0,0,0,0};
  for(Int_t k=1;k<=2;k++){
    size+=AliGeomManager::LayerSize(k);
  }
  for(Int_t k=3;k<=6;k++){
    size2+=AliGeomManager::LayerSize(k);
    printf("size: %d \n",size2);
  }
  
  printf("Aligning SPDBarrel: nmodules: %d \n",size);  
  printf("Fitting modules: %d \n",size2);

  TArrayI *volIDs=GetLayersVolUID(layers);
  layers[0]=0;
  layers[1]=0;
  layers[2]=1;
  layers[3]=1;
  layers[4]=1;
  layers[5]=1;
  TArrayI *volIDsFit=GetLayersVolUID(layers);   

  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSDD1,AliGeomManager::kTPC1,iterations);
  
  return kTRUE; 
}

//______________________
Bool_t AliITSRealignTracks::AlignSPDHalfBarrel(Int_t method,Int_t iterations){
  //Align a SPD Half barrel "iterations" times
  //method 0 : align SPDHalfBarrel Up without using the points on SPD Half Barrel down in the fits (only outer layers)
  //method 1 : align SPDHalfBarrel Down without using the points on SPD Half Barrel up in the fits (only outer layers)
  //method 10 : align SPDHalfBarrel Up using also the points on SPD Half Barrel down in the fits (and points on outer layers)
  //method 11 : align SPDHalfBarrel Down using also the points on SPD Half Barrel up in the fits (and points on outer layers)

  Int_t size=0,size2=0;
  Int_t layers[6]={0,0,1,1,1,1};
  Int_t sectorsUp[10]={1,1,1,1,1,0,0,0,0,0}; 
  Int_t sectorsDown[10]={0,0,0,0,0,1,1,1,1,1}; 
  
  TString updownstr;
  if(method==0)updownstr="UpNoDown";
  else if (method==1)updownstr="DownNoUp";
  else if (method==10)updownstr="UpWithDown";
  else if (method==11)updownstr="DownWithUp";
  else {
    AliWarning("Wrong AlignSPDHalfBarrel method selected ");
    return kFALSE;
  }
  
  for(Int_t i=1;i<=2;i++){
    size+=AliGeomManager::LayerSize(i);
  }
  
  for(Int_t i=3;i<=6;i++){
    size2+=AliGeomManager::LayerSize(i);
  }
  
  size=size/2;
  if(method==10||method==11)size2+=size;
  
  printf("Aligning  SPDHalfBarrel %s: nmodules: %d \n",updownstr.Data(),size);
  printf("Fitting modules: %d \n",size2);
  TArrayI *volIDsFit2;
  TArrayI *volids = NULL;
  TArrayI *volIDsFit=GetLayersVolUID(layers);
  if(method==0||method==10)volids=GetSPDSectorsVolids(sectorsUp);
  if(method==1||method==11)volids=GetSPDSectorsVolids(sectorsDown);

  if(method==10)volIDsFit2=JoinVolArrays(GetSPDSectorsVolids(sectorsDown),volIDsFit);
  else if(method==11)volIDsFit2=JoinVolArrays(GetSPDSectorsVolids(sectorsUp),volIDsFit);
  else volIDsFit2=volIDsFit;
  
  AlignVolumesITS(volids,volIDsFit2,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iterations);
  
  return kTRUE; 
}


//______________________________________________________
Bool_t AliITSRealignTracks::AlignLayer(Int_t layer,Int_t iterations){
  //Align the layer "layer" iterations times

  Int_t size=0,size2=0;
  Int_t layers[6]={0,0,0,0,0,0};
  layers[layer-1]=1;
  TString layerstr[6]={"SPD1","SPD2","SDD1","SDD2","SSD1","SSD2"};
  for(Int_t k=1;k<=6;k++){
    if(k!=layer)size2+=AliGeomManager::LayerSize(k);
  }
  size=AliGeomManager::LayerSize(layer);
  
  printf("Aligning layer %s, nmodules %d ,fitted modules %d \n",layerstr[layer-1].Data(),size,size2);
  
  
  TArrayI *volIDs=GetLayersVolUID(layers);
  layers[0]=1;
  layers[1]=1;
  layers[2]=1;
  layers[3]=1;
  layers[4]=1;
  layers[5]=1;
  layers[layer]=0;
  TArrayI *volIDsFit=GetLayersVolUID(layers);   
  
  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSDD1,AliGeomManager::kSSD2,iterations);
  
  return kTRUE; 
}

//___________________________________________

Bool_t AliITSRealignTracks::AlignLayersToLayers(const Int_t *layer,Int_t iterations){

  //Align the set of layers A with respect to the set of layers B iterations time.
  //The two sets A and B are defined into *layer==layer[6] the following way:
  //   layer[i]=0 the layer is skipped both in the fits than in the minimization
  //   layer[i]=1 the layer is skipped in the fits and considered in the minimization
  //   layer[i]=2 the layer is considered in the fits and skipped in the minimization
  //   layer[i]=3 the layer is considered both in the fits and in the minimization
  
  UShort_t volid;
  Int_t size=0,size2=0,j=0,k=0;
  Int_t iLayer;
  TString layerstr[6]={"SPD1","SPD2","SDD1","SDD2","SSD1","SSD2"};
  TString command="",str;
  for(Int_t i=1;i<=6;i++){
    if(layer[i-1]==1||layer[i-1]==3){
      size+=AliGeomManager::LayerSize(i);
      command.Append(" ");
      command.Append(layerstr[i-1]);
    }
    if(layer[i-1]==2||layer[i-1]==3){
      size2+=AliGeomManager::LayerSize(i);
      str.Append(" ");
      str.Append(layerstr[i-1]);
    }
  }
  
  printf("Aligning layers %s To layers %s, nmodules %d ,fitted modules %d \n",command.Data(),str.Data(),size,size2);
  
  
  TArrayI volIDs(size);
  TArrayI volIDsFit(size2);   
  
  for (iLayer=(Int_t)AliGeomManager::kSPD1;iLayer<(Int_t)AliGeomManager::kTPC1;iLayer++){
    if(layer[iLayer-AliGeomManager::kFirstLayer]==0)continue;
    if(layer[iLayer-AliGeomManager::kFirstLayer]==1){
      for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){
	volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	volIDs.AddAt(volid,j);
	j++;
      }
    }
    else if(layer[iLayer-AliGeomManager::kFirstLayer]==2){
      for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){
	volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	volIDsFit.AddAt(volid,k);
	k++;
      }
    }
    else if(layer[iLayer-AliGeomManager::kFirstLayer]==3){
      for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){
	volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	volIDs.AddAt(volid,j);
	j++;
	volIDsFit.AddAt(volid,k);
	k++;
      }
    }
  }
  
  AlignVolumesITS(&volIDs,&volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iterations);
   
  return kTRUE; 
}

//______________________________________________

Bool_t AliITSRealignTracks::AlignSPDSectorToOuterLayers(Int_t sector,Int_t iterations){
  //Align the SPD sector "sector" with respect to outer layers iterations times

  
  Int_t layers[6]={0,0,1,1,1,1};
  Bool_t spd=kFALSE;
  Int_t sectorsIN[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t sectorsFit[10]={1,1,1,1,1,1,1,1,1,1};
  
  if(sector<0){
    sector=-sector;
    spd=kTRUE;
  }
  sectorsIN[sector]=1;
  sectorsFit[sector]=0;
  TArrayI *volIDs=GetSPDSectorsVolids(sectorsIN);
  TArrayI *volIDsFit;
  if(spd){
    volIDsFit=JoinVolArrays(GetSPDSectorsVolids(sectorsFit),GetLayersVolUID(layers));
  }
  else volIDsFit=GetLayersVolUID(layers);
  
  printf("Aligning SPD sector %d: nmodules: %d \n",sector,volIDs->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit->GetSize());
  
  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSPD2,iterations);
  
  return kTRUE; 
}

//______________________________________________
Bool_t AliITSRealignTracks::AlignSPDSectorWithSectors(Int_t sector,Int_t iterations){
  //Align the SPD sector "sector" with respect to the other SPD sectors iterations times

  Int_t sectorsIN[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t sectorsFit[10]={1,1,1,1,1,1,1,1,1,1};

  sectorsIN[sector]=1;
  sectorsFit[sector]=0;
  TArrayI *volIDs=GetSPDSectorsVolids(sectorsIN);
  TArrayI *volIDsFit=GetSPDSectorsVolids(sectorsFit);;   
  
  printf("Aligning SPD sector %d: nmodules: %d \n",sector,volIDs->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit->GetSize());
 

  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSDD1,iterations);
  
  return kTRUE; 
}



//___________________________________________________
Bool_t AliITSRealignTracks::AlignSPDSectorsWithSectors(const Int_t *sectorsIN,const Int_t *sectorsFit,Int_t iterations){
  //Align SPD sectors defined in "sectorsIN" with respect to 
  //SPD sectors defined in "sectorsFit" iterations time

  TArrayI *volIDs=GetSPDSectorsVolids(sectorsIN);
  TArrayI *volIDsFit=GetSPDSectorsVolids(sectorsFit);;   
  
  printf("Aligning SPD sectors: modules: %d \n",volIDs->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit->GetSize());
  
  
  
  return AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSDD1,iterations);; 
}

//___________________________________________________
Bool_t AliITSRealignTracks::AlignSPDStaves(const Int_t *staves,const Int_t *sectorsIN,const Int_t *sectorsFit,Int_t iterations){
  //Align SPD staves defined by staves and sectorsIN with respect to sectorsFit volumes iterations times

  TArrayI *volIDs=GetSPDStavesVolids(sectorsIN,staves);
  TArrayI *volIDsFit=GetSPDSectorsVolids(sectorsFit);   
  
  if(volIDs->GetSize()==0){
    printf("EMPTY ARRAY !! \n");
    return kFALSE;
  }
  printf("Aligning SPD staves: modules: %d \n",volIDs->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit->GetSize());

  TArrayI *volIDsFit2=ExcludeVolidsFromVolidsArray(volIDs,volIDsFit);
  return  AlignVolumesITS(volIDs,volIDsFit2,AliGeomManager::kSPD1,AliGeomManager::kSSD1,iterations); 
}


//___________________________________________

Bool_t AliITSRealignTracks::AlignLayerToSPDHalfBarrel(Int_t layer,Int_t updown,Int_t iterations){
  //Align the layer "layer" with respect to SPD Half Barrel Up (updowon=0) 
  //or Down (updown=1) iterations times
  


  Int_t sectorsDown[10]={0,0,0,0,0,1,1,1,1,1};
  Int_t sectorsUp[10]={1,1,1,1,1,0,0,0,0,0};
  TString layerstr[6]={"SPD1","SPD2","SDD1","SDD2","SSD1","SSD2"};  
  TArrayI *volIDsFit;
  Int_t layers[6]={0,0,0,0,0,0};
  layers[layer-1]=1;
  Int_t size=AliGeomManager::LayerSize(layer);
  TArrayI *volIDs=GetLayersVolUID(layers);

  if(updown==0){
    volIDsFit=GetSPDSectorsVolids(sectorsUp);   
    printf("Aligning layer %s, nmodules %d ,to half barrel Up \n",layerstr[layer-1].Data(),size);
  }
  else if(updown==1){
    volIDsFit=GetSPDSectorsVolids(sectorsDown);
    printf("Aligning layer %s, nmodules %d ,to half barrel Down \n",layerstr[layer-1].Data(),size);
  }
  else {
    printf("Wrong Half Barrel selection! \n");
    return kFALSE;
  }
 
  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iterations);
  
  return kTRUE; 
}

//___________________________________________

Bool_t AliITSRealignTracks::AlignLayerToSector(Int_t layer,Int_t sector,Int_t iterations){
  //Align the layer "layer" with respect to SPD sector "sector" iterations times

  if(sector>9){
    printf("Wrong Sector selection! \n");
    return kFALSE;
  }
  Int_t sectors[10]={0,0,0,0,0,0,0,0,0,0};
  sectors[sector]=1;
  TString layerstr[6]={"SPD1","SPD2","SDD1","SDD2","SSD1","SSD2"};  
  TArrayI *volIDsFit;
  Int_t layers[6]={0,0,0,0,0,0};
  layers[layer-1]=1;
  TArrayI *volIDs=GetLayersVolUID(layers);
  Int_t size=AliGeomManager::LayerSize(layer); 
  
 
  volIDsFit=GetSPDSectorsVolids(sectors);   
  printf("Aligning layer %s, nmodules %d ,to half barrel Up \n",layerstr[layer-1].Data(),size);
  
  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iterations);
  
  return kTRUE; 
}

//_______________________________________________

Bool_t AliITSRealignTracks::AlignSPDHalfBarrelToHalfBarrel(Int_t updown,Int_t iterations){
  //Align the SPD Half Barrel Up[Down] with respect to HB Down[Up] iterations time if
  //updown=0[1]

  
  Int_t sectorsDown[10]={0,0,0,0,0,1,1,1,1,1};
  Int_t sectorsUp[10]={1,1,1,1,1,0,0,0,0,0};
  
  TArrayI *volIDsUp=GetSPDSectorsVolids(sectorsUp);
  TArrayI *volIDsDown=GetSPDSectorsVolids(sectorsDown);   
  
  if(updown==0){
    printf("Aligning SPD HalfBarrel up to half Barrel down : nmodules: %d \n",volIDsUp->GetSize());  
    printf("Fitting modules: %d \n",volIDsDown->GetSize());
    AlignVolumesITS(volIDsUp,volIDsDown,AliGeomManager::kSPD1,AliGeomManager::kSPD2,iterations);
  }
  else if(updown==1){
    printf("Aligning SPD HalfBarrel down to half Barrel Up : nmodules: %d \n",volIDsDown->GetSize());  
    printf("Fitting modules: %d \n",volIDsUp->GetSize()); 
    AlignVolumesITS(volIDsDown,volIDsUp,AliGeomManager::kSPD1,AliGeomManager::kSPD2,iterations);
  }
  else {
    printf("Wrong Half Barrel selection! \n");
    return kFALSE;
  }
  
  return kTRUE; 
}


//_______________________
Bool_t AliITSRealignTracks::AlignSPDHalfBarrelToSectorRef(Int_t sector,Int_t iterations){
  //Align the SPD Half Barrel Down with respect to sector "sector" iterations times

  Int_t sectorsIN[10]={0,0,0,0,0,1,1,1,1,1};
  Int_t sectorsFit[10]={0,0,0,0,0,0,0,0,0,0};

  sectorsFit[sector]=1;

  TArrayI *volIDs=GetSPDSectorsVolids(sectorsIN);
  TArrayI *volIDsFit=GetSPDSectorsVolids(sectorsFit);   
  
  printf("Aligning SPD HalfBarrel to sector 0 %d: nmodules: %d \n",sector,volIDs->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit->GetSize());
 

  AlignVolumesITS(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSPD2,iterations);

 
  return kTRUE; 
}
//_________________________________________
Bool_t AliITSRealignTracks::AlignSPD1SectorRef(Int_t sector,Int_t iterations){
  //OBSOLETE METHOD: Align the SPD1 modules of sector "sector" with respect 
  // to the other SPD volumes iterations times

  Int_t sectorsIN[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t sectorsFit[10]={1,1,1,1,1,1,1,1,1,1};
  sectorsIN[sector]=1;
  sectorsFit[sector]=0;
  TArrayI *volIDs=GetSPDSectorsVolids(sectorsIN);
  TArrayI *volIDsFit=GetSPDSectorsVolids(sectorsFit);   
  Int_t size=volIDs->GetSize();
  Int_t size2=volIDsFit->GetSize();
  UShort_t volID;
  Int_t k=0;

  TArrayI *volIDsSPD1=new TArrayI(size-8);
  TArrayI *volIDsFit2=new TArrayI(size2+8);
  
  for(Int_t j=0;j<size;j++){
    volID=volIDs->At(j);
    if(AliGeomManager::VolUIDToLayer(volID)==AliGeomManager::kSPD1){
      volIDsSPD1->AddAt(volID,size2+k);
      k++;
    }
    else volIDsFit2->AddAt(volID,j-k);
  }
  
  
  for(Int_t j=0;j<size2;j++){
    volID=volIDsFit->At(j);
    volIDsFit2->AddAt(volID,size-k+j);
  }
  
  printf("Aligning SPD Sector %d: nmodules: %d \n",sector,volIDsSPD1->GetSize());  
  printf("Fitting modules: %d \n",volIDsFit2->GetSize());
  
  AlignVolumesITS(volIDsSPD1,volIDsFit2,AliGeomManager::kSPD1,AliGeomManager::kSPD2,iterations);
  
  return kTRUE; 
}

//_____________________________________________

AliAlignObjParams* AliITSRealignTracks::MediateAlignObj(const TArrayI *volIDs,Int_t lastVolid){
  //TEMPORARY METHOD: perform an average of the values of the parameters of the AlignObjs 
  // defined by the array volIDs up to lastVolid position in this array
  //The aim of such a method is to look for collective movement of a given set of modules

  UShort_t volid;

  TGeoHMatrix hm;
  Double_t *rot,*transl;
  Double_t rotSum[9],translSum[3]={0.,0.,0.};
  for(Int_t k=0;k<8;k++)rotSum[k]=0.;


  for(Int_t ivol=0;ivol<lastVolid;ivol++){
    volid=volIDs->At(ivol);
  
    GetAlignObj(volIDs->At(ivol))->GetMatrix(hm); 
   
    rot=hm.GetRotationMatrix();
    transl=hm.GetTranslation();
   
    for(Int_t j=0;j<9;j++)rotSum[j]+=rot[j];
    for(Int_t jt=0;jt<3;jt++)translSum[jt]+=transl[jt];
  }
  if(lastVolid!=0){
    for(Int_t j=0;j<9;j++)rotSum[j]=rotSum[j]/lastVolid;
    for(Int_t jt=0;jt<3;jt++)translSum[jt]=translSum[jt]/lastVolid;
  }
  else printf("Try to mediate results for zero modules \n");
 
  hm.SetRotation(rotSum);
  hm.SetTranslation(translSum);


  
  AliAlignObjParams *alignObj=new AliAlignObjParams("average", 0,hm, kTRUE);
  return alignObj;
  
}


//________________________________________________
TArrayI* AliITSRealignTracks::GetSPDStavesVolids(const Int_t *sectors,const Int_t* staves){

  
  // This method gets the volID Array for the chosen staves into the 
  // chosen sectors. You have to pass an array (10 dim) with a 1 for each 
  // selected sector and an array (6 dim) with a 1 for each chosen stave.    
  // The staves are numbered in this way: 0,1 for SPD1 and 2,3,4,5 for SPD2
  // i.e. sectors[10] = {1,1,0,0,0,0,0,0,1,0} -> Sector 0, 1, 9 selected.
  // staves[6]={0,1,1,0,0,1} -> Staves 1 on SPD1 and 0 and 3 on SPD2 selected
  
  Int_t nSect=0,nStaves=0;
  Int_t last=0;

 
  for(Int_t co=0;co<10;co++){ //counts the number of sectors chosen
    if(sectors[co]==1) nSect++;
  }

  for(Int_t co=0;co<6;co++){ //counts the number of sectors chosen
    if(staves[co]==1) nStaves++;
  }
  
  if(nSect<1||nStaves<1){ //if no sector chosen -> exit
    Printf("Error! No Sector/s or staves Selected!");
    return 0x0;
  }

  TArrayI *volIDs = new TArrayI(nSect*nStaves*4);
  TString stave="/Stave",str,symn,laystr;
  
  TArrayI *sectvol=GetSPDSectorsVolids(sectors); 
  //SPD1 
  laystr="SPD0";
  for(Int_t k=0;k<2;k++){
    if(staves[k]==1){
      str=stave;
      str+=k;
      for(Int_t i=0;i<sectvol->GetSize();i++){
	symn=AliGeomManager::SymName(sectvol->At(i));
	if(symn.Contains(str)&&symn.Contains(laystr)){
	  //	  printf("Adding: %s \n",symn.Data());
	  volIDs->AddAt(sectvol->At(i),last);
	  last++;
	}
      }
    }
  }
  //SPD1 
  laystr="SPD1";
  for(Int_t k=2;k<6;k++){
    if(staves[k]==1){
      str=stave;
      str+=k-2;
      for(Int_t i=0;i<sectvol->GetSize();i++){
	symn=AliGeomManager::SymName(sectvol->At(i));
	if(symn.Contains(str)&&symn.Contains(laystr)){
	  volIDs->AddAt(sectvol->At(i),last);
	  printf("Adding: %s \n",symn.Data());
	  last++;
	}
      }
    }
  }

  volIDs->Set(last);
  return volIDs;
} 

//________________________________________________
TArrayI* AliITSRealignTracks::GetSPDSectorsVolids(const Int_t *sectors) 
{
  //
  // This method gets the volID Array for the chosen sectors.
  // You have to pass an array with a 1 for each selected sector.
  // i.e. sectors[10] = {1,1,0,0,0,0,0,0,1,0} -> Sector 0, 1, 9 selected.
  //

  Int_t nSect=0;
  Int_t iModule=0;

 
  for(Int_t co=0;co<10;co++){ //counts the number of sectors chosen
    if(sectors[co]==1) nSect++;
  }
  
  if(nSect<1){ //if no sector chosen -> exit
    Printf("Error! No Sector/s Selected!");
    return 0x0;
  }

  TArrayI *volIDs = new TArrayI(nSect*24);
  
    if(sectors[0]==1){ //--->cSect = 0 <---
      volIDs->AddAt(2048,iModule); iModule++;
      volIDs->AddAt(2049,iModule); iModule++;
      volIDs->AddAt(2050,iModule); iModule++;
      volIDs->AddAt(2051,iModule); iModule++;
      volIDs->AddAt(2052,iModule); iModule++;
      volIDs->AddAt(2053,iModule); iModule++;
      volIDs->AddAt(2054,iModule); iModule++;
      volIDs->AddAt(2055,iModule); iModule++;
      volIDs->AddAt(4096,iModule); iModule++;
      volIDs->AddAt(4097,iModule); iModule++;
      volIDs->AddAt(4098,iModule); iModule++;
      volIDs->AddAt(4099,iModule); iModule++;
      volIDs->AddAt(4100,iModule); iModule++;
      volIDs->AddAt(4101,iModule); iModule++;
      volIDs->AddAt(4102,iModule); iModule++;
      volIDs->AddAt(4103,iModule); iModule++;
      volIDs->AddAt(4104,iModule); iModule++;
      volIDs->AddAt(4105,iModule); iModule++;
      volIDs->AddAt(4106,iModule); iModule++;
      volIDs->AddAt(4107,iModule); iModule++;
      volIDs->AddAt(4108,iModule); iModule++;
      volIDs->AddAt(4109,iModule); iModule++;
      volIDs->AddAt(4110,iModule); iModule++;
      volIDs->AddAt(4111,iModule); iModule++;
    }
    if(sectors[1]==1){ //--->cSect = 1 <//---
      volIDs->AddAt(2056,iModule); iModule++;
      volIDs->AddAt(2057,iModule); iModule++;
      volIDs->AddAt(2058,iModule); iModule++;
      volIDs->AddAt(2059,iModule); iModule++;
      volIDs->AddAt(2060,iModule); iModule++;
      volIDs->AddAt(2061,iModule); iModule++;
      volIDs->AddAt(2062,iModule); iModule++;
      volIDs->AddAt(2063,iModule); iModule++;
      volIDs->AddAt(4112,iModule); iModule++;
      volIDs->AddAt(4113,iModule); iModule++;
      volIDs->AddAt(4114,iModule); iModule++;
      volIDs->AddAt(4115,iModule); iModule++;
      volIDs->AddAt(4116,iModule); iModule++;
      volIDs->AddAt(4117,iModule); iModule++;
      volIDs->AddAt(4118,iModule); iModule++;
      volIDs->AddAt(4119,iModule); iModule++;
      volIDs->AddAt(4120,iModule); iModule++;
      volIDs->AddAt(4121,iModule); iModule++;
      volIDs->AddAt(4122,iModule); iModule++;
      volIDs->AddAt(4123,iModule); iModule++;
      volIDs->AddAt(4124,iModule); iModule++;
      volIDs->AddAt(4125,iModule); iModule++;
      volIDs->AddAt(4126,iModule); iModule++;
      volIDs->AddAt(4127,iModule); iModule++;
    }
    if(sectors[2]==1){//--->cSect = 2 <//---
      volIDs->AddAt(2064,iModule); iModule++;
      volIDs->AddAt(2065,iModule); iModule++;
      volIDs->AddAt(2066,iModule); iModule++;
      volIDs->AddAt(2067,iModule); iModule++;
      volIDs->AddAt(2068,iModule); iModule++;
      volIDs->AddAt(2069,iModule); iModule++;
      volIDs->AddAt(2070,iModule); iModule++;
      volIDs->AddAt(2071,iModule); iModule++;
      volIDs->AddAt(4128,iModule); iModule++;
      volIDs->AddAt(4129,iModule); iModule++;
      volIDs->AddAt(4130,iModule); iModule++;
      volIDs->AddAt(4131,iModule); iModule++;
      volIDs->AddAt(4132,iModule); iModule++;
      volIDs->AddAt(4133,iModule); iModule++;
      volIDs->AddAt(4134,iModule); iModule++;
      volIDs->AddAt(4135,iModule); iModule++;
      volIDs->AddAt(4136,iModule); iModule++;
      volIDs->AddAt(4137,iModule); iModule++;
      volIDs->AddAt(4138,iModule); iModule++;
      volIDs->AddAt(4139,iModule); iModule++;
      volIDs->AddAt(4140,iModule); iModule++;
      volIDs->AddAt(4141,iModule); iModule++;
      volIDs->AddAt(4142,iModule); iModule++;
      volIDs->AddAt(4143,iModule); iModule++;
    }
    if(sectors[3]==1){//--->cSect = 3 <//---
      volIDs->AddAt(2072,iModule); iModule++;
      volIDs->AddAt(2073,iModule); iModule++;
      volIDs->AddAt(2074,iModule); iModule++;
      volIDs->AddAt(2075,iModule); iModule++;
      volIDs->AddAt(2076,iModule); iModule++;
      volIDs->AddAt(2077,iModule); iModule++;
      volIDs->AddAt(2078,iModule); iModule++;
      volIDs->AddAt(2079,iModule); iModule++;
      volIDs->AddAt(4144,iModule); iModule++;
      volIDs->AddAt(4145,iModule); iModule++;
      volIDs->AddAt(4146,iModule); iModule++;
      volIDs->AddAt(4147,iModule); iModule++;
      volIDs->AddAt(4148,iModule); iModule++;
      volIDs->AddAt(4149,iModule); iModule++;
      volIDs->AddAt(4150,iModule); iModule++;
      volIDs->AddAt(4151,iModule); iModule++;
      volIDs->AddAt(4152,iModule); iModule++;
      volIDs->AddAt(4153,iModule); iModule++;
      volIDs->AddAt(4154,iModule); iModule++;
      volIDs->AddAt(4155,iModule); iModule++;
      volIDs->AddAt(4156,iModule); iModule++;
      volIDs->AddAt(4157,iModule); iModule++;
      volIDs->AddAt(4158,iModule); iModule++;
      volIDs->AddAt(4159,iModule); iModule++;
    }
    if(sectors[4]==1){//--->cSect = 4 <//---
      volIDs->AddAt(2080,iModule); iModule++;
      volIDs->AddAt(2081,iModule); iModule++;
      volIDs->AddAt(2082,iModule); iModule++;
      volIDs->AddAt(2083,iModule); iModule++;
      volIDs->AddAt(2084,iModule); iModule++;
      volIDs->AddAt(2085,iModule); iModule++;
      volIDs->AddAt(2086,iModule); iModule++;
      volIDs->AddAt(2087,iModule); iModule++;
      volIDs->AddAt(4160,iModule); iModule++;
      volIDs->AddAt(4161,iModule); iModule++;
      volIDs->AddAt(4162,iModule); iModule++;
      volIDs->AddAt(4163,iModule); iModule++;
      volIDs->AddAt(4164,iModule); iModule++;
      volIDs->AddAt(4165,iModule); iModule++;
      volIDs->AddAt(4166,iModule); iModule++;
      volIDs->AddAt(4167,iModule); iModule++;
      volIDs->AddAt(4168,iModule); iModule++;
      volIDs->AddAt(4169,iModule); iModule++;
      volIDs->AddAt(4170,iModule); iModule++;
      volIDs->AddAt(4171,iModule); iModule++;
      volIDs->AddAt(4172,iModule); iModule++;
      volIDs->AddAt(4173,iModule); iModule++;
      volIDs->AddAt(4174,iModule); iModule++;
      volIDs->AddAt(4175,iModule); iModule++;
    }
    if(sectors[5]==1){//--->cSect = 5 <//---
      volIDs->AddAt(2088,iModule); iModule++;
      volIDs->AddAt(2089,iModule); iModule++;
      volIDs->AddAt(2090,iModule); iModule++;
      volIDs->AddAt(2091,iModule); iModule++;
      volIDs->AddAt(2092,iModule); iModule++;
      volIDs->AddAt(2093,iModule); iModule++;
      volIDs->AddAt(2094,iModule); iModule++;
      volIDs->AddAt(2095,iModule); iModule++;
      volIDs->AddAt(4176,iModule); iModule++;
      volIDs->AddAt(4177,iModule); iModule++;
      volIDs->AddAt(4178,iModule); iModule++;
      volIDs->AddAt(4179,iModule); iModule++;
      volIDs->AddAt(4180,iModule); iModule++;
      volIDs->AddAt(4181,iModule); iModule++;
      volIDs->AddAt(4182,iModule); iModule++;
      volIDs->AddAt(4183,iModule); iModule++;
      volIDs->AddAt(4184,iModule); iModule++;
      volIDs->AddAt(4185,iModule); iModule++;
      volIDs->AddAt(4186,iModule); iModule++;
      volIDs->AddAt(4187,iModule); iModule++;
      volIDs->AddAt(4188,iModule); iModule++;
      volIDs->AddAt(4189,iModule); iModule++;
      volIDs->AddAt(4190,iModule); iModule++;
      volIDs->AddAt(4191,iModule); iModule++;
    }
    if(sectors[6]==1){//--->cSect = 6 <//---
      volIDs->AddAt(2096,iModule); iModule++;
      volIDs->AddAt(2097,iModule); iModule++;
      volIDs->AddAt(2098,iModule); iModule++;
      volIDs->AddAt(2099,iModule); iModule++;
      volIDs->AddAt(2100,iModule); iModule++;
      volIDs->AddAt(2101,iModule); iModule++;
      volIDs->AddAt(2102,iModule); iModule++;
      volIDs->AddAt(2103,iModule); iModule++;
      volIDs->AddAt(4192,iModule); iModule++;
      volIDs->AddAt(4193,iModule); iModule++;
      volIDs->AddAt(4194,iModule); iModule++;
      volIDs->AddAt(4195,iModule); iModule++;
      volIDs->AddAt(4196,iModule); iModule++;
      volIDs->AddAt(4197,iModule); iModule++;
      volIDs->AddAt(4198,iModule); iModule++;
      volIDs->AddAt(4199,iModule); iModule++;
      volIDs->AddAt(4200,iModule); iModule++;
      volIDs->AddAt(4201,iModule); iModule++;
      volIDs->AddAt(4202,iModule); iModule++;
      volIDs->AddAt(4203,iModule); iModule++;
      volIDs->AddAt(4204,iModule); iModule++;
      volIDs->AddAt(4205,iModule); iModule++;
      volIDs->AddAt(4206,iModule); iModule++;
      volIDs->AddAt(4207,iModule); iModule++;
    }
     if(sectors[7]==1){ //--->cSect = 7 <//---
       volIDs->AddAt(2104,iModule); iModule++;
       volIDs->AddAt(2105,iModule); iModule++;
       volIDs->AddAt(2106,iModule); iModule++;
       volIDs->AddAt(2107,iModule); iModule++;
       volIDs->AddAt(2108,iModule); iModule++;
       volIDs->AddAt(2109,iModule); iModule++;
       volIDs->AddAt(2110,iModule); iModule++;
       volIDs->AddAt(2111,iModule); iModule++;
       volIDs->AddAt(4208,iModule); iModule++;
       volIDs->AddAt(4209,iModule); iModule++;
       volIDs->AddAt(4210,iModule); iModule++;
       volIDs->AddAt(4211,iModule); iModule++;
       volIDs->AddAt(4212,iModule); iModule++;
       volIDs->AddAt(4213,iModule); iModule++;
       volIDs->AddAt(4214,iModule); iModule++;
       volIDs->AddAt(4215,iModule); iModule++;
       volIDs->AddAt(4216,iModule); iModule++;
       volIDs->AddAt(4217,iModule); iModule++;
       volIDs->AddAt(4218,iModule); iModule++;
       volIDs->AddAt(4219,iModule); iModule++;
       volIDs->AddAt(4220,iModule); iModule++;
       volIDs->AddAt(4221,iModule); iModule++;
       volIDs->AddAt(4222,iModule); iModule++;
       volIDs->AddAt(4223,iModule); iModule++;
     }
     if(sectors[8]==1){//--->cSect = 8 <//---
       volIDs->AddAt(2112,iModule); iModule++;
       volIDs->AddAt(2113,iModule); iModule++;
       volIDs->AddAt(2114,iModule); iModule++;
       volIDs->AddAt(2115,iModule); iModule++;
       volIDs->AddAt(2116,iModule); iModule++;
       volIDs->AddAt(2117,iModule); iModule++;
       volIDs->AddAt(2118,iModule); iModule++;
       volIDs->AddAt(2119,iModule); iModule++;
       volIDs->AddAt(4224,iModule); iModule++;
       volIDs->AddAt(4225,iModule); iModule++;
       volIDs->AddAt(4226,iModule); iModule++;
       volIDs->AddAt(4227,iModule); iModule++;
       volIDs->AddAt(4228,iModule); iModule++;
       volIDs->AddAt(4229,iModule); iModule++;
       volIDs->AddAt(4230,iModule); iModule++;
       volIDs->AddAt(4231,iModule); iModule++;
       volIDs->AddAt(4232,iModule); iModule++;
       volIDs->AddAt(4233,iModule); iModule++;
       volIDs->AddAt(4234,iModule); iModule++;
       volIDs->AddAt(4235,iModule); iModule++;
       volIDs->AddAt(4236,iModule); iModule++;
       volIDs->AddAt(4237,iModule); iModule++;
       volIDs->AddAt(4238,iModule); iModule++;
       volIDs->AddAt(4239,iModule); iModule++;
     }
     if(sectors[9]==1){//--->cSect = 9 <//---
       volIDs->AddAt(2120,iModule); iModule++;
       volIDs->AddAt(2121,iModule); iModule++;
       volIDs->AddAt(2122,iModule); iModule++;
       volIDs->AddAt(2123,iModule); iModule++;
       volIDs->AddAt(2124,iModule); iModule++;
       volIDs->AddAt(2125,iModule); iModule++;
       volIDs->AddAt(2126,iModule); iModule++;
       volIDs->AddAt(2127,iModule); iModule++;
       volIDs->AddAt(4240,iModule); iModule++;
       volIDs->AddAt(4241,iModule); iModule++;
       volIDs->AddAt(4242,iModule); iModule++;
       volIDs->AddAt(4243,iModule); iModule++;
       volIDs->AddAt(4244,iModule); iModule++;
       volIDs->AddAt(4245,iModule); iModule++;
       volIDs->AddAt(4246,iModule); iModule++;
       volIDs->AddAt(4247,iModule); iModule++;
       volIDs->AddAt(4248,iModule); iModule++;
       volIDs->AddAt(4249,iModule); iModule++;
       volIDs->AddAt(4250,iModule); iModule++;
       volIDs->AddAt(4251,iModule); iModule++;
       volIDs->AddAt(4252,iModule); iModule++;
       volIDs->AddAt(4253,iModule); iModule++;
       volIDs->AddAt(4254,iModule); iModule++;
       volIDs->AddAt(4255,iModule); iModule++;
     }

  return volIDs;
}

//___________________________________
TArrayI* AliITSRealignTracks::GetLayersVolUID(const Int_t *layer){

  //return a TArrayI with the volUIDs of the modules into the set of layers
  //defined by layer[6]
  
  TArrayI *out=new TArrayI(2198);
  Int_t last=0;
  UShort_t voluid;
  for(Int_t i=0;i<6;i++){
    if(layer[i]==1){
      for(Int_t mod=0;mod<AliGeomManager::LayerSize(i+AliGeomManager::kFirstLayer);mod++){
	voluid=AliGeomManager::LayerToVolUID(i+AliGeomManager::kFirstLayer,mod);
	out->AddAt(voluid,last);
	//	printf("voluid %d at position %d \n",out->At(last),last);
	last++;
      }
    }
  }  
  out->Set(last);
  return out;
}

//_________________
TArrayI* AliITSRealignTracks::SelectLayerInVolids(const TArrayI *volidsIN,AliGeomManager::ELayerID layer){
  //Select between the modules specified by their volUIDs in volidsIN only those
  // of a given layer "layer"

  Int_t size=volidsIN->GetSize();
  Int_t count=0;
  for(Int_t j=0;j<size;j++){
    if(AliGeomManager::VolUIDToLayer(volidsIN->At(j))==layer)count++;
  }
  TArrayI *volidsOUT=new TArrayI(count);
  count=0;
  for(Int_t j=0;j<size;j++){
    if(AliGeomManager::VolUIDToLayer(volidsIN->At(j))==layer){
      volidsOUT->AddAt(volidsIN->At(j),count);
      count++;
    }
  }
  return volidsOUT;
}

//______________________________________________

TArrayI* AliITSRealignTracks::IntersectVolArray(const TArrayI *vol1,const TArrayI *vol2){

  //Perform the intersection between the array vol1 and vol2
  
  Int_t size1=vol1->GetSize();
  Int_t size2=vol2->GetSize();
  Int_t last=0,volid;
  Bool_t found;
  TArrayI *volidOut=new TArrayI(size1+size2);  
  
  for(Int_t k=0;k<size1;k++){
    found=kFALSE;
    volid=vol1->At(k);
    for(Int_t j=0;j<size2;j++){
      if(vol2->At(j)==volid)found=kTRUE;
    }
    if(found){
      volidOut->AddAt(volid,last);
      last++;
    }
  }
  volidOut->Set(last);
  return volidOut;
}
//_________________________________________

TArrayI* AliITSRealignTracks::JoinVolArrays(const TArrayI *vol1,const TArrayI *vol2){
  //!BE CAREFUL: If an index is repeated into vol1 or into vol2 will be repeated also in the final array
  
  Int_t size1=vol1->GetSize();
  Int_t size2=vol2->GetSize();
  Int_t count=0;
  UShort_t volid;
  Bool_t found;
  TArrayI *volidOut=new TArrayI(size1+size2);  
  
  for(Int_t k=0;k<size1;k++){
    volid=vol1->At(k);
    volidOut->AddAt(volid,k);
  }
 
  for(Int_t k=0;k<size2;k++){
    found=kFALSE;
    volid=vol2->At(k);
    for(Int_t j=0;j<size1;j++){
      if(volidOut->At(j)==volid)found=kTRUE;
    }
    if(!found){
      volidOut->AddAt(volid,size1+count);
      count++;
    }
  }
  volidOut->Set(size1+count);
  return volidOut;
}

//______________________________________

TArrayI* AliITSRealignTracks::ExcludeVolidsFromVolidsArray(const TArrayI *volidsToExclude,const TArrayI *volStart){
  //Excludes the modules defined by their volUID in the array volidsToExclude from the array volStart

  Int_t size1=volidsToExclude->GetSize();
  Int_t size2=volStart->GetSize();
  Int_t last=0;
  UShort_t volid;
  Bool_t found;
  TArrayI *volidOut=new TArrayI(size2);  

  for(Int_t k=0;k<size2;k++){
    found=kFALSE;
    volid=volStart->At(k);
    for(Int_t j=0;j<size1;j++){
      if(volidsToExclude->At(j)==volid){
	found=kTRUE;
	break;
      }
    }
    if(!found){
      volidOut->AddAt(volid,last);
      last++;
    }
  }
  volidOut->Set(last);
  return volidOut;
}


//________________________________________

TArrayI* AliITSRealignTracks::GetLayerVolumes(const Int_t *layer){
  //returns a TArrayI with the volUIDs of the modules of the layers
  //specified into *layer
  
  TArrayI *out=new TArrayI(2198);
  Int_t last=0;
  UShort_t voluid;
  for(Int_t i=0;i<6;i++){
    if(layer[i]==1){
      for(Int_t mod=0;mod<AliGeomManager::LayerSize(i+AliGeomManager::kFirstLayer);mod++){
	voluid=AliGeomManager::LayerToVolUID(i+AliGeomManager::kFirstLayer,mod);
	out->AddAt(voluid,last);
	//	printf("voluid %d at position %d \n",out->At(last),last);
	last++;
      }
    }
  }  
  out->Set(last);
  return out;
}



//______________________________
TArrayI* AliITSRealignTracks::GetAlignedVolumes(char *filename){
  //Open the file "filename" which is expected to contain
  //a TClonesArray named "ITSAlignObjs" with stored a set of AlignObjs
  //returns an array with the volumes UID of the modules considered realigned 
  
  if(gSystem->AccessPathName(filename)){
    printf("Wrong Realignment file name \n");
    return 0x0;
  }
  TFile *f=TFile::Open(filename,"READ");
  TClonesArray *array=(TClonesArray*)f->Get("ITSAlignObjs");
  AliAlignObjParams *a;
  Int_t last=0;
  TArrayI *volidOut=new TArrayI(2200);
  for(Int_t j=0;j<array->GetSize();j++){
    a=(AliAlignObjParams*)array->At(j);
    if(a->GetUniqueID()==0)continue;
    
    else {
      volidOut->AddAt(a->GetVolUID(),last);
      last++;								      
    }
  }
  volidOut->Set(last);
  f->Close();
  return volidOut;
}


//________________________________________
void AliITSRealignTracks::SetDraw(Bool_t draw,Bool_t refresh){
  //TEPMORARY METHOD: method to switch on/off the drawing of histograms
  // if refresh=kTRUE deletes the old histos and constructs new ones
  
  if(refresh){
    // WriteHists();
    if(fAlignDrawObjs)DeleteDrawHists();
    InitDrawHists();
  }
  fDraw=draw;
  return;
}

void AliITSRealignTracks::DeleteDrawHists(){
  //Delete the pointers to the histograms

  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      delete fAlignDrawObjs[iLayer][iModule];
    }
    if(fAlignDrawObjs[iLayer])delete [] fAlignDrawObjs[iLayer];
  }
  
  delete [] fAlignDrawObjs;
  fAlignDrawObjs = 0;

  
  delete fCanvPar;
  delete fCanvGr; 
  delete fgrIterMeanX;
  delete fgrIterRMSX; 
  delete fgrIterMeanY;
  delete fgrIterRMSY;
  delete fgrIterMeanZ;
  delete fgrIterRMSZ;
  delete fgrIterMeanPsi;
  delete fgrIterRMSPsi;
  delete fgrIterMeanTheta;
  delete fgrIterRMSTheta;
  delete fgrIterMeanPhi;
  delete fgrIterRMSPhi;  

} 

void AliITSRealignTracks::InitDrawHists(){
  //Initialize the histograms to monitor the results

  Int_t nLayers = AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer;
  fAlignDrawObjs = new AliAlignObj**[nLayers];
  for (Int_t iLayer = 0; iLayer < (AliGeomManager::kLastLayer - AliGeomManager::kFirstLayer); iLayer++) {
    fAlignDrawObjs[iLayer] = new AliAlignObj*[AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer)];
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer + AliGeomManager::kFirstLayer); iModule++) {
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer+ AliGeomManager::kFirstLayer,iModule);
      fAlignDrawObjs[iLayer][iModule] = new AliAlignObjParams(AliGeomManager::SymName(volid),volid,0,0,0,0,0,0,kTRUE);
      fAlignDrawObjs[iLayer][iModule]->SetUniqueID(1);
    }
  }


  TH1F *hX=new TH1F("hX","hX",1000,-10000.,10000.);
  TH1F *hY=new TH1F("hY","hY",1000,-10000.,10000.);
  TH1F *hZ=new TH1F("hZ","hZ",1000,-10000.,10000.);
  TH1F *hPsi=new TH1F("hPsi","hPsi",1000,-5000.,5000.);
  TH1F *hTheta=new TH1F("hTheta","hTheta",1000,-5000.,5000.);
  TH1F *hPhi=new TH1F("hPhi","hPhi",1000,-5000.,5000.);

  fCanvPar=new TCanvas("fCanvPar","Parameters trend during iterations: Convergence \n");
  fCanvPar->Divide(3,2);
  fCanvPar->cd(1);
  hX->Draw();
  hX->SetXTitle("#mum");
  fCanvPar->cd(2);
  hY->Draw();
  hY->SetXTitle("#mum");
  fCanvPar->cd(3);
  hZ->SetXTitle("#mum");
  hZ->Draw();
  fCanvPar->cd(4);
  hPsi->SetXTitle("mdeg");
  hPsi->Draw();
  fCanvPar->cd(5);
  hTheta->SetXTitle("mdeg");
  hTheta->Draw();
  fCanvPar->cd(6);
  hPhi->SetXTitle("mdeg");
  hPhi->Draw();
  fCanvPar->Update();

  
  fCanvGr=new TCanvas("fCanvGr","Parameters trend during iterations: Convergence \n");
  fCanvGr->Divide(3,2);

  fCanvGr->cd(1);
  fgrIterMeanX=new TGraph(1);
  fgrIterRMSX=new TGraph(1);
  fgrIterRMSX->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSX->SetName("fgrIterRMSX");
  fgrIterRMSX->SetLineColor(2);
  fgrIterMeanX->SetName("fgrIterMeanX");
  fgrIterMeanX->SetTitle("Convergence of #deltaX \n");
  fgrIterMeanX->GetXaxis()->SetTitle("#mum");
  fgrIterRMSX->Draw("acp");
  fgrIterMeanX->Draw("cp");

  fCanvGr->cd(2);
  fgrIterMeanY=new TGraph(1);
  fgrIterRMSY=new TGraph(1);
  fgrIterRMSY->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSY->SetName("fgrIterRMSY");
  fgrIterRMSY->SetLineColor(2);
  fgrIterMeanY->SetName("fgrIterMeanY");
  fgrIterMeanY->SetTitle("Convergence of #deltaY \n");
  fgrIterMeanY->GetXaxis()->SetTitle("#mum");
  fgrIterRMSY->Draw("acp");
  fgrIterMeanY->Draw("cp");

  fCanvGr->cd(3);
  fgrIterMeanZ=new TGraph(1);
  fgrIterRMSZ=new TGraph(1);
  fgrIterRMSZ->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSZ->SetName("fgrIterRMSZ");
  fgrIterRMSZ->SetLineColor(2);
  fgrIterMeanZ->SetName("fgrIterMeanZ");
  fgrIterMeanZ->SetTitle("Convergence of #deltaZ \n");
  fgrIterMeanZ->GetXaxis()->SetTitle("#mum");
  fgrIterRMSZ->Draw("acp");
  fgrIterMeanZ->Draw("cp");

  fCanvGr->cd(4);
  fgrIterMeanPsi=new TGraph(1);
  fgrIterRMSPsi=new TGraph(1);
  fgrIterRMSPsi->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSPsi->SetName("fgrIterRMSPsi");
  fgrIterRMSPsi->SetLineColor(2);
  fgrIterMeanPsi->SetName("fgrIterMeanPsi");
  fgrIterMeanPsi->SetTitle("Convergence of #deltaPsi \n");
  fgrIterMeanPsi->GetXaxis()->SetTitle("mdeg");
  fgrIterRMSPsi->Draw("acp");
  fgrIterMeanPsi->Draw("cp");

  fCanvGr->cd(5);
  fgrIterMeanTheta=new TGraph(1);
  fgrIterRMSTheta=new TGraph(1);
  fgrIterRMSTheta->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSTheta->SetName("fgrIterRMSTheta");
  fgrIterRMSTheta->SetLineColor(2);
  fgrIterMeanTheta->SetName("fgrIterMeanTheta");
  fgrIterMeanTheta->SetTitle("Convergence of #deltaTheta \n");
  fgrIterMeanTheta->GetXaxis()->SetTitle("mdeg");
  fgrIterRMSTheta->Draw("acp");
  fgrIterMeanTheta->Draw("cp");

  fCanvGr->cd(6);
  fgrIterMeanPhi=new TGraph(1);
  fgrIterRMSPhi=new TGraph(1);
  fgrIterRMSPhi->GetYaxis()->SetRangeUser(-1000.,1000.);
  fgrIterRMSPhi->SetName("fgrIterRMSPhi");
  fgrIterRMSPhi->SetLineColor(2);
  fgrIterMeanPhi->SetName("fgrIterMeanPhi");
  fgrIterMeanPhi->SetTitle("Convergence of #deltaPhi \n");
  fgrIterMeanPhi->GetXaxis()->SetTitle("mdeg");
  fgrIterRMSPhi->Draw("acp");
  fgrIterMeanPhi->Draw("cp");


  
}

void AliITSRealignTracks::UpdateDraw(TArrayI *volids,Int_t iter,Int_t color){
  //Updates the histograms to monitor the results. Only the histograms
  //of the volumes specified in *volids will be updated.
  // iter is just a flag for the names of the histo
  // color specifies the color of the lines of the histograms for this update
  
  TString name="hX_";
  name+=iter;
  name.Append("iter");
  TH1F *hX=new TH1F("hX",name.Data(),1000,-10000.,10000.);

  name="hY_";
  name+=iter;
  name.Append("iter");
  TH1F *hY=new TH1F("hY",name.Data(),1000,-10000.,10000.);

  name="hZ_";
  name+=iter;
  name.Append("iter");
  TH1F *hZ=new TH1F("hZ",name.Data(),1000,-10000.,10000.);

  name="hPsi_";
  name+=iter;
  name.Append("iter");
  TH1F *hPsi=new TH1F("hPsi",name.Data(),1000,-5000.,5000.);
  
  name="hTheta_";
  name+=iter;
  name.Append("iter");
  TH1F *hTheta=new TH1F("hTheta",name.Data(),1000,-5000.,5000.);

  name="hPhi_";
  name+=iter;
  name.Append("iter");
  TH1F *hPhi=new TH1F("hPhi",name.Data(),1000,-5000.,5000.);
  
  Int_t layer,mod;
  Double_t transl[3],rot[3],transldr[3],rotdr[3];
  
  for(Int_t i=0;i<volids->GetSize();i++){
    layer=AliGeomManager::VolUIDToLayer(volids->At(i),mod); 
    fAlignObjs[layer-AliGeomManager::kFirstLayer][mod]->GetPars(transl,rot);
    fAlignDrawObjs[layer-AliGeomManager::kFirstLayer][mod]->GetPars(transldr,rotdr);
    
    hX->Fill(10000.*(transl[0]-transldr[0]));
    hY->Fill(10000.*(transl[1]-transldr[1]));
    hZ->Fill(10000.*(transl[2]-transldr[2]));
    hPsi->Fill(1000.*(rot[0]-rotdr[0]));
    hTheta->Fill(1000.*(rot[1]-rotdr[1]));
    hPhi->Fill(1000.*(rot[1]-rotdr[2]));
    //Update the pars of the draw object
    fAlignDrawObjs[layer-AliGeomManager::kFirstLayer][mod]->SetPars(transl[0],transl[1],transl[2],rot[0],rot[1],rot[2]);
  }

  hX->SetLineColor(color);
  hY->SetLineColor(color);
  hZ->SetLineColor(color);
  hPsi->SetLineColor(color);
  hTheta->SetLineColor(color);
  hPhi->SetLineColor(color);
  
  
  fCanvPar->cd(1);
  hX->Draw("Same");
  fCanvPar->cd(2);
  hY->Draw("Same");
  fCanvPar->cd(3);
  hZ->Draw("Same");
  fCanvPar->cd(4);
  hPsi->Draw("Same");
  fCanvPar->cd(5);
  hTheta->Draw("Same");
  fCanvPar->cd(6);
  hPhi->Draw("Same");
  gPad->Modified();
  fCanvPar->Update();
  fCanvPar->Modified();
  
  fgrIterMeanX->SetPoint(fgrIterMeanX->GetN()+1,iter,hX->GetMean());
  fgrIterRMSX->SetPoint(fgrIterRMSX->GetN()+1,iter,hX->GetRMS());
  fgrIterMeanY->SetPoint(fgrIterMeanY->GetN()+1,iter,hY->GetMean());
  fgrIterRMSY->SetPoint(fgrIterRMSY->GetN()+1,iter,hY->GetRMS());
  fgrIterMeanZ->SetPoint(fgrIterMeanZ->GetN()+1,iter,hZ->GetMean());
  fgrIterRMSZ->SetPoint(fgrIterRMSZ->GetN()+1,iter,hZ->GetRMS());
  fgrIterMeanPsi->SetPoint(fgrIterMeanPsi->GetN()+1,iter,hPsi->GetMean());
  fgrIterRMSPsi->SetPoint(fgrIterRMSPsi->GetN()+1,iter,hPsi->GetRMS());
  fgrIterMeanTheta->SetPoint(fgrIterMeanTheta->GetN()+1,iter,hTheta->GetMean());
  fgrIterRMSTheta->SetPoint(fgrIterRMSTheta->GetN()+1,iter,hTheta->GetRMS());
  fgrIterMeanPhi->SetPoint(fgrIterMeanPhi->GetN()+1,iter,hPhi->GetMean());
  fgrIterRMSPhi->SetPoint(fgrIterRMSPhi->GetN()+1,iter,hPhi->GetRMS());

  gPad->Modified();
  fCanvGr->Update();
  fCanvGr->Update();
}

void AliITSRealignTracks::WriteHists(const char *outfile){
  //Writes the histograms for the monitoring of the results
  // in a file named "outfile"

  TFile *f=new TFile(outfile,"RECREATE");
  f->cd();
  fCanvPar->Write();
  fCanvGr->Write(); 
  fgrIterMeanX->Write(); 
  fgrIterRMSX->Write();  
  fgrIterMeanY->Write(); 
  fgrIterRMSY->Write();  
  fgrIterMeanZ->Write(); 
  fgrIterRMSZ->Write();  
  fgrIterMeanPsi->Write(); 
  fgrIterRMSPsi->Write();  
  fgrIterMeanTheta->Write(); 
  fgrIterRMSTheta->Write();  
  fgrIterMeanPhi->Write();  
  fgrIterRMSPhi->Write();

  f->Close();
  return;
  
}
