#include <TArray.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TArray.h>
#include <TNtuple.h>
#include <TGeoManager.h>
#include "AliGeomManager.h"
#include "AliAlignmentTracks.h"
#include "AliTrackFitter.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackResidualsFast.h"
#include "AliTrackResidualsChi2.h"
#include "AliTrackResidualsLinear.h"


void RealignmentITSIndependVol(TString minimizer="fast",const Int_t fit=0,const Int_t iter1=1,const Int_t iterations=5,const Int_t minNtracks=-10,const Int_t layer=0,const Int_t minTrackPoint=6,TString fileintro="AliTrackPoints.root",TString geometryfile="geometry.root"){
  
  // minimizer="fast"->AliTrackResidualFast minimizer
  //           "minuit"->AliTrackResidualChi2 minimizer
  //           "minuitnorot"->AliTrackResidualChi2 minimizer without rotations degrees of freedom
  //           "linear"->AliTrackResidualLinear minimizer    
  //fit=0-> Riemann Fitter, fit=1->Kalman
  //iter1=#iterations inside AliAlignmentTracks::AliAlignVolumes method 
  //iterations=#iterations of the entire procedure
  //layer=0->all ITS, otherways the usual notation is considered (1=SPD1,2=SPD2,3=SDD1,4=SDD2,5=SSD1,6=SSD2)
  //minNtracks=minimun number of tracks passing through a module in order to try to realign the module itsself
  //           if minNtracks<0, minimun number of tracks is |minNtracks|*minNumPoint[layer]/fact (see the code below): this allows a different
  //           choice of the number of tracks required on different layers and to vary these numbers once tuned the relative proportions.  
  //minTrackPoint=minimun number of "good" points required to a track (THE POINT ON THE MODULE THAT IS GOING TO BE REALIGNED 
  //              IS NEVER CONSIDERED->max number can be required is 11 for cosmics tracks) for the track being considered in the minimization
  //fileintro=file into which the Tree with the space points is stored
  //geometryfile=file containing the geometry


  TArrayI volIDs2(2200); 
  volIDs2.Reset(0);
  TArrayI volIDs(1);
  TString outname;
  Int_t layerNum,modNum,iLayer,iLayerToAlign,j=0,size=0,lastVolid=0;
  Int_t minNumPoint[6]={200,200,200,200,100,100}; 
  Double_t fact=10; 
  TNtuple *ntVolumeAlign=new TNtuple("ntVolumeAlign","NTuple with volume tried to be realigned","layerNum:modNum:volumeIDnum"); 


  AliAlignmentTracks *AliAlTrack=new AliAlignmentTracks();
 
  AliAlTrack->SetPointsFilename(fileintro.Data());

  AliTrackFitter *fitter;
  if(fit==1)fitter= new AliTrackFitterKalman();
  else fitter=new AliTrackFitterRieman();

  fitter->SetMinNPoints(minTrackPoint);

  AliAlTrack->SetTrackFitter(fitter);


  AliTrackResiduals *res;
  if(minimizer=="minuit"){
    res = new AliTrackResidualsChi2();
  }
  else if(minimizer=="minuitnorot"){
    res = new AliTrackResidualsChi2();
    res->FixParameter(3);
    res->FixParameter(4);
    res->FixParameter(5);
  }

  else if(minimizer=="fast"){
    res = new AliTrackResidualsFast();
  } 
  else if(minimizer=="linear"){
    res = new AliTrackResidualsLinear();
  } 

  else {
    printf("Trying to set a non existing minimizer! \n");
    return;
  }

  res->SetMinNPoints(1);
  AliAlTrack->SetMinimizer(res);
  
  if(!gGeoManager) TGeoManager::Import(geometryfile.Data());
   
  
  TStopwatch *timer=new TStopwatch();
  timer->Start(); 
  AliAlTrack->BuildIndex();
  j=0;
  UShort_t volid;
  for(Int_t iter=0;iter<iterations;iter++){
    for(iLayerToAlign=1;iLayerToAlign<=6;iLayerToAlign++){
      if(layer!=0&&(Int_t)iLayerToAlign!=layer)continue;
      j=0;
      size=0;
      for(Int_t k=1;k<=6;k++){
	size+=(Int_t)AliGeomManager::LayerSize(k);
	printf("size: %d \n",size);
      }
  
      
      for (Int_t iModule=0;iModule<AliGeomManager::LayerSize(iLayerToAlign);iModule++){      
	j=0;
	if(minNtracks<0){	
	  if(AliAlTrack->GetLastIndex(iLayerToAlign-1,iModule)<minNumPoint[iLayerToAlign-1]*(-1*minNtracks/fact))continue;
	}	
	else if(AliAlTrack->GetLastIndex(iLayerToAlign-1,iModule)<minNtracks)continue;

	UShort_t volidAl = AliGeomManager::LayerToVolUID(iLayerToAlign,iModule);
	
	TArrayI volIDsFit(size-1);
	for (iLayer=(Int_t)AliGeomManager::kSPD1;iLayer<(Int_t)AliGeomManager::kTPC1;iLayer++){
	  for (Int_t iModule2=0;iModule2<AliGeomManager::LayerSize(iLayer);iModule2++){
	    volid = AliGeomManager::LayerToVolUID(iLayer,iModule2);
	    if(volid==volidAl)continue;
	    volIDsFit.AddAt(volid,j);
	    j++;
	  }
	}
	volIDs.AddAt((Int_t)volidAl,0);
	if(iter==iterations-1){
	  volIDs2.AddAt(volidAl,lastVolid);
	  lastVolid++;
	}
	AliAlTrack->AlignVolumes(&volIDs,&volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,iter1);  
      }
    }
    
    if((iter+1)%5==0||iter==0||iter==1||iter==2||iter==3||iter==iterations-1){
      outname="RealignObj";
      outname+=(iter+1);
      outname.Append(".root");
      AliAlTrack->WriteRealignObjArray(outname.Data(),AliGeomManager::kSPD1,AliGeomManager::kSSD2);
    }
  }
  
  if(lastVolid==0){printf("No modules could be realigned \n");return;}
  printf("End of selecting modules cycle: %d modules selected \n",lastVolid);
  for(Int_t k=0;k<volIDs2.GetSize();k++){
    if(volIDs2.At(k)==0)break;
    layerNum=AliGeomManager::VolUIDToLayer(volIDs2.At(k),modNum);
    ntVolumeAlign->Fill(layerNum,modNum,volIDs2.At(k));
  }
  TFile *f=new TFile("RealignVolNt.root","RECREATE");
  f->cd();
  ntVolumeAlign->Write();
  f->Close();

  timer->Stop();
  timer->Print();
  return;
}


