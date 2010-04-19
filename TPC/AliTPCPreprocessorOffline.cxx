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



/*
   .x ~/rootlogon.C
   gSystem->Load("libANALYSIS");
   gSystem->Load("libTPCcalib");
   .L $ALICE_ROOT/TPC/AliTPCPreprocessorOffline.cxx+

   AliTPCPreprocessorOffline proces;
   proces.CalibTimeGain(
 
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMap.h"
#include "TGraphErrors.h"
#include "AliExternalTrackParam.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2D.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"

#include "AliESDfriend.h"


#include "AliTPCcalibTime.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTPCcalibBase.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "AliRelAlignerKalman.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibTimeGain.h"
#include "AliSplineFit.h"
#include "AliTPCPreprocessorOffline.h"
#endif


ClassImp(AliTPCPreprocessorOffline)

AliTPCPreprocessorOffline::AliTPCPreprocessorOffline():
  TNamed("PreprocessorOffline","PreprocessorOffline"),
  kMinEntries(500),                      // minimal number of entries for fit
  startRun(0),                         // start Run - used to make fast selection in THnSparse
  endRun(0),                           // end   Run - used to make fast selection in THnSparse
  startTime(0),                        // startTime - used to make fast selection in THnSparse
  endTime(0),                          // endTime   - used to make fast selection in THnSparse
  ocdbStorage(""),                   // path to the OCDB storage
  fVdriftArray(new TObjArray),
  fTimeDrift(0),
  fGraphMIP(0),                // graph time dependence of MIP
  fGraphCosmic(0),             // graph time dependence at Plateu
  fFitMIP(0),                  // fit of dependence - MIP
  fFitCosmic(0),               // fit of dependence - Plateu
  fGainArray(new TObjArray),               // array to be stored in the OCDB
  fGainMIP(0),          // calibration component for MIP
  fGainCosmic(0)       // calibration component for cosmic
{
  //
}

AliTPCPreprocessorOffline::~AliTPCPreprocessorOffline() {
  //
  // Destructor
  //
}




void AliTPCPreprocessorOffline::GetRunRange(AliTPCcalibTime* fTimeDrift){
  //
  // find the fist and last run
  //
  TObjArray *hisArray =fTimeDrift->GetHistoDrift();
  {for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
    THnSparse* addHist=(THnSparse*)hisArray->UncheckedAt(i);
    if (addHist->GetEntries()<kMinEntries) continue;
    if (!addHist) continue;
    TH1D* histo    =addHist->Projection(3);
    TH1D* histoTime=addHist->Projection(0);
    printf("%s\t%f\t%d\t%d\n",histo->GetName(), histo->GetEntries(),histo->FindFirstBinAbove(0),histo->FindLastBinAbove(0));

    if (startRun<=0){ 
      startRun=histo->FindFirstBinAbove(0);
      endRun  =histo->FindLastBinAbove(0);
    }else{
      startRun=TMath::Min(histo->FindFirstBinAbove(0),startRun);
      endRun  =TMath::Max(histo->FindLastBinAbove(0),endRun);
    }
    if (startTime==0){ 
      startTime=histoTime->FindFirstBinAbove(0);
      endTime  =histoTime->FindLastBinAbove(0);
    }else{
      startTime=TMath::Min(histoTime->FindFirstBinAbove(0),startTime);
      endTime  =TMath::Max(histoTime->FindLastBinAbove(0),endTime);
    }
    delete histo;
    delete histoTime;
  }}
  if (startRun<0) startRun=0;
  if (endRun<0) endRun=100000000;
  printf("Run range  :\t%d-%d\n", startRun, endRun);
  printf("Time range :\t%d-%d\n", startTime, endTime);

}



void AliTPCPreprocessorOffline::CalibTimeVdrift(Char_t* file, Int_t ustartRun, Int_t uendRun, TString pocdbStorage){
  //
  //
  //
  const Int_t    kMinEntries=500;     // minimal number of entries
  if (pocdbStorage.Length()>0) ocdbStorage=pocdbStorage;
  else
  ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  // 1. Initialization and run range setting
  TFile fcalib(file);
  fTimeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");
  startRun=ustartRun;
  endRun=ustartRun; 
  TObjArray *hisArray =fTimeDrift->GetHistoDrift();  
  GetRunRange(fTimeDrift);
  for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
    THnSparse* addHist=(THnSparse*)hisArray->At(i);
    if (!addHist) continue;
    if (startTime<endTime) addHist->GetAxis(0)->SetRange(startTime-1,endTime+1);
    if (startRun<endRun) addHist->GetAxis(3)->SetRange(startRun-1,endRun+1);
  }
  //
  //
  // 2. extraction of the information
  //
  fVdriftArray = new TObjArray();
  AddAlignmentGraphs(fVdriftArray,fTimeDrift);
  AddHistoGraphs(fVdriftArray,fTimeDrift,kMinEntries);
  AddLaserGraphs(fVdriftArray,fTimeDrift);
  //
  // 3. Append QA plots
  //
  MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. update of OCDB
  //
  //
  
  UpdateOCDBDrift(ustartRun,uendRun,ocdbStorage);
}

void AliTPCPreprocessorOffline::UpdateOCDBDrift( Int_t ustartRun, Int_t uendRun,  const char* storagePath ){
  //
  // Update OCDB 
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Calibration of the time dependence of the drift velocity");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/TimeDrift", ustartRun, uendRun);
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
  gStorage->Put(fVdriftArray, (*id1), metaData);
}



void AliTPCPreprocessorOffline::UpdateDriftParam(AliTPCParam *param, TObjArray *arr, Int_t startRun){
  //
  //  update the OCDB entry for the nominal time0
  //
  //
  //  AliTPCParam * param = AliTPCcalibDB::Instance()->GetParameters();
  AliTPCParam *paramNew = (AliTPCParam *)param->Clone();
  TGraphErrors *grT =  (TGraphErrors *)arr->FindObject("ALIGN_ITSM_TPC_T0");
  Double_t deltaTcm = TMath::Median(grT->GetN(),grT->GetY());
  Double_t deltaT   = deltaTcm/param->GetDriftV();
  paramNew->SetL1Delay(param->GetL1Delay()-deltaT);
  paramNew->Update();

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-02"); //root version
  metaData->SetComment("Updated calibration of nominal time 0");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Parameters", startRun, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(param, (*id1), metaData);

}


void AliTPCPreprocessorOffline::PrintArray(TObjArray *array){
  //
  //
  //
  Int_t entries = array->GetEntries();
  for (Int_t i=0; i<entries; i++){
    if (!array->At(i)) continue;
    printf("%d\t %s\n", i,  array->At(i)->GetName());
  }
}



TGraphErrors* AliTPCPreprocessorOffline::FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs){
  // 2 filters:
  //    1. filter graph - error cut errSigmaCut
  //    2. filter graph - medianCutAbs around median
  //
  // errSigmaCut   - cut on error
  // medianCutAbs  - cut on value around median
  Double_t dummy=0;               //   
  //
  // 1. filter graph - error cut errSigmaCut
  //              
  TGraphErrors *graphF; 
  graphF = AliTPCcalibDButil::FilterGraphMedianErr(graph,errSigmaCut,dummy);
  delete graph;
  if (!graphF) return 0;
  graph = AliTPCcalibDButil::FilterGraphMedianErr(graphF,errSigmaCut,dummy);
  delete graphF;
  if (!graph) return 0;
  //
  // filter graph - kMedianCutAbs around median
  // 
  graphF=FilterGraphMedianAbs(graph, medianCutAbs,dummy);
  delete graph;
  if (!graphF) return 0;
  graph=FilterGraphMedianAbs(graphF, medianCutAbs,dummy);
  delete graphF;
  if (!graph) return 0;
  return graph;
}



TGraphErrors* AliTPCPreprocessorOffline::FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY){
  //
  // filter outlyer measurement
  // Only points around median +- cut filtered 
  //
  if (!graph) return  0;
  Int_t kMinPoints=2;
  Int_t npoints0 = graph->GetN();
  Int_t npoints=0;
  Float_t  rmsY=0;
  Double_t *outx=new Double_t[npoints0];
  Double_t *outy=new Double_t[npoints0];
  Double_t *errx=new Double_t[npoints0];
  Double_t *erry=new Double_t[npoints0];
  //
  //
  if (npoints0<kMinPoints) return 0;
  for (Int_t iter=0; iter<3; iter++){
    npoints=0;
    for (Int_t ipoint=0; ipoint<npoints0; ipoint++){
      if (graph->GetY()[ipoint]==0) continue;
      if (iter>0 &&TMath::Abs(graph->GetY()[ipoint]-medianY)>cut) continue;  
      outx[npoints]  = graph->GetX()[ipoint];
      outy[npoints]  = graph->GetY()[ipoint];
      errx[npoints]  = graph->GetErrorX(ipoint);
      erry[npoints]  = graph->GetErrorY(ipoint);
      npoints++;
    }
    if (npoints<=1) break;
    medianY  =TMath::Median(npoints,outy);
    rmsY   =TMath::RMS(npoints,outy);
  }
  TGraphErrors *graphOut=0;
  if (npoints>1) graphOut= new TGraphErrors(npoints,outx,outy,errx,erry); 
  return graphOut;
}


void AliTPCPreprocessorOffline::AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift, Int_t minEntries){
  //
  // Add graphs corresponding to the alignment
  //
  const Double_t kErrSigmaCut=5;      // error sigma cut - for filtering
  const Double_t kMedianCutAbs=0.03;  // error sigma cut - for filtering
  //
  TObjArray * array=timeDrift->GetHistoDrift();
  if (array){
    THnSparse* hist=NULL;
    // 2.a) cosmics with different triggers
    for (Int_t i=0; i<array->GetEntriesFast();i++){
      hist=(THnSparseF*)array->UncheckedAt(i);
      if(!hist) continue;
      if (hist->GetEntries()<minEntries) continue;
      //hist->Print();
      TString name=hist->GetName();
      Int_t dim[4]={0,1,2,3};
      THnSparse* newHist=hist->Projection(4,dim);
      newHist->SetName(name);
      TGraphErrors* graph=AliTPCcalibBase::FitSlices(newHist,2,0,400,100,0.05,0.95, kTRUE);
      printf("name=%s graph=%i, N=%i\n", name.Data(), graph==0, graph->GetN());
      Int_t pos=name.Index("_");
      name=name(pos,name.Capacity()-pos);
      TString graphName=graph->ClassName();
      graphName+=name;
      graphName.ToUpper();
      //
      graph = FilterGraphDrift(graph, kErrSigmaCut, kMedianCutAbs);
      if (!graph) {
	printf("Graph =%s filtered out\n", name.Data());
	continue;
      }
      //
      graph->SetMarkerStyle(i%8+20);
      graph->SetMarkerColor(i%7);
      graph->GetXaxis()->SetTitle("Time");
      graph->GetYaxis()->SetTitle("v_{dcor}");
      graph->SetName(graphName);
      graph->SetTitle(graphName);
      printf("Graph %d\t=\t%s\n", i, graphName.Data());
      vdriftArray->Add(graph);
    }
  }
}




void AliTPCPreprocessorOffline::AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift){
  //
  // Add graphs corresponding to alignment to the object array
  //
  TObjArray *arrayITS=0;
  TObjArray *arrayTOF=0;
  TObjArray *arrayTRD=0;
  TMatrixD *mstatITS=0;
  TMatrixD *mstatTOF=0;
  TMatrixD *mstatTRD=0;
  //
  arrayITS=timeDrift->GetAlignITSTPC();
  arrayTRD=timeDrift->GetAlignTRDTPC();
  arrayTOF=timeDrift->GetAlignTOFTPC();

  if (arrayITS->GetEntries()>0) mstatITS= AliTPCcalibDButil::MakeStatRelKalman(arrayITS,0.9,50,0.025);
  if (arrayTOF->GetEntries()>0) mstatTOF= AliTPCcalibDButil::MakeStatRelKalman(arrayTOF,0.9,1000,0.025);
  if (arrayTRD->GetEntries()>0) mstatTRD= AliTPCcalibDButil::MakeStatRelKalman(arrayTRD,0.9,50,0.025);
  //
  TObjArray * arrayITSP= AliTPCcalibDButil::SmoothRelKalman(arrayITS,*mstatITS, 0, 5.);
  TObjArray * arrayITSM= AliTPCcalibDButil::SmoothRelKalman(arrayITS,*mstatITS, 1, 5.);
  TObjArray * arrayITSB= AliTPCcalibDButil::SmoothRelKalman(arrayITSP,arrayITSM);
  TObjArray * arrayTOFP= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,*mstatTOF, 0, 5.);
  TObjArray * arrayTOFM= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,*mstatTOF, 1, 5.);
  TObjArray * arrayTOFB= AliTPCcalibDButil::SmoothRelKalman(arrayTOFP,arrayTOFM);

  TObjArray * arrayTRDP= 0x0;
  TObjArray * arrayTRDM= 0x0;
  TObjArray * arrayTRDB= 0x0;
  arrayTRDP= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,*mstatTRD, 0, 5.);
  arrayTRDM= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,*mstatTRD, 1, 5.);
  arrayTRDB= AliTPCcalibDButil::SmoothRelKalman(arrayTRDP,arrayTRDM);
  //
  //
  Int_t entries=TMath::Max(arrayITS->GetEntriesFast(),arrayTOF->GetEntriesFast());
  TObjArray *arrays[12]={arrayITS, arrayITSP, arrayITSM, arrayITSB,
			 arrayTRD, arrayTRDP, arrayTRDM, arrayTRDB,
			 arrayTOF, arrayTOFP, arrayTOFM, arrayTOFB};
  TString   grnames[12]={"ALIGN_ITS", "ALIGN_ITSP", "ALIGN_ITSM", "ALIGN_ITSB",
			 "ALIGN_TRD", "ALIGN_TRDP", "ALIGN_TRDM","ALIGN_TRDB",
			 "ALIGN_TOF", "ALIGN_TOFP", "ALIGN_TOFM","ALIGN_TOFB"};
  TString   grpar[9]={"DELTAPSI", "DELTATHETA", "DELTAPHI",
		      "DELTAX", "DELTAY", "DELTAZ",
		      "DRIFTVD", "T0", "VDGY"};

  
  TVectorD vX(entries);
  TVectorD vY(entries);
  TVectorD vEx(entries);
  TVectorD vEy(entries);
  TObjArray *arr=0;
  for (Int_t iarray=0; iarray<12; iarray++){
    arr = arrays[iarray];
    if (arr==0) continue;
    for (Int_t ipar=0; ipar<9; ipar++){      
      Int_t counter=0;
      for (Int_t itime=0; itime<arr->GetEntriesFast(); itime++){
	AliRelAlignerKalman * kalman = (AliRelAlignerKalman *) arr->UncheckedAt(itime);
	if (!kalman) continue;
	vX[counter]=kalman->GetTimeStamp();
	vY[counter]=(*(kalman->GetState()))[ipar];
	if (ipar==6) vY[counter]=1./(*(kalman->GetState()))[ipar]-1;
	vEx[counter]=0;
	vEy[counter]=TMath::Sqrt((*(kalman->GetStateCov()))(ipar,ipar));
	counter++;
      }
    
      TGraphErrors * graph=new TGraphErrors(counter, vX.GetMatrixArray(),
					  vY.GetMatrixArray(),
					  vEx.GetMatrixArray(),
					  vEy.GetMatrixArray());
      TString grName=grnames[iarray];
      grName+="_TPC_";
      grName+=grpar[ipar];
      graph->SetName(grName.Data());
      vdriftArray->AddLast(graph);
    }
  }  
}




void AliTPCPreprocessorOffline::AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift){
  //
  // add graphs for laser
  //
  const Double_t delayL0L1 = 0.071;  //this is hack for 1/2 weeks
  THnSparse *hisN=0;
  TGraphErrors *grLaser[6]={0,0,0,0,0,0};
  hisN = timeDrift->GetHistVdriftLaserA(0);
  if (timeDrift->GetHistVdriftLaserA(0)){
    grLaser[0]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(0),0,2,5,delayL0L1);
    grLaser[0]->SetName("GRAPH_MEAN_DELAY_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[0]);
  }    
  if (timeDrift->GetHistVdriftLaserA(1)){
    grLaser[1]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(1),0,2,5);
    grLaser[1]->SetName("GRAPH_MEAN_DRIFT_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[1]);
  }    
  if (timeDrift->GetHistVdriftLaserA(2)){
    grLaser[2]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserA(2),0,2,5);
    grLaser[2]->SetName("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
    vdriftArray->AddLast(grLaser[2]);
  }    
  if (timeDrift->GetHistVdriftLaserC(0)){
    grLaser[3]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(0),0,2,5,delayL0L1);
    grLaser[3]->SetName("GRAPH_MEAN_DELAY_LASER_ALL_C");
    vdriftArray->AddLast(grLaser[3]);
  }    
  if (timeDrift->GetHistVdriftLaserC(1)){
    grLaser[4]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(1),0,2,5);
    grLaser[4]->SetName("GRAPH_MEAN_DRIFT_LASER_ALL_C");
    vdriftArray->AddLast(grLaser[4]);
  }    
  if (timeDrift->GetHistVdriftLaserC(2)){
    grLaser[5]=MakeGraphFilter0(timeDrift->GetHistVdriftLaserC(2),0,2,5);
    grLaser[5]->SetName("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");    
    vdriftArray->AddLast(grLaser[5]);
  }    
  for (Int_t i=0; i<6;i++){
    if (grLaser[i]) {
      SetDefaultGraphDrift(grLaser[i], 1,(i+20));
      grLaser[i]->GetYaxis()->SetTitle("Laser Correction");
    }
  }
}
 
 
TGraphErrors * AliTPCPreprocessorOffline::MakeGraphFilter0(THnSparse *hisN, Int_t itime, Int_t ival, Int_t minEntries, Double_t offset){
  //
  // Make graph with mean values and rms
  //
  hisN->GetAxis(itime)->SetRange(0,100000000);
  hisN->GetAxis(ival)->SetRange(0,100000000);
  TH1 * hisT      = hisN->Projection(itime);
  TH1 * hisV      = hisN->Projection(ival);
  //
  Int_t firstBinA = hisT->FindFirstBinAbove(2);
  Int_t lastBinA  = hisT->FindLastBinAbove(2);    
  Int_t firstBinV = hisV->FindFirstBinAbove(0);
  Int_t lastBinV  = hisV->FindLastBinAbove(0);    
  hisN->GetAxis(itime)->SetRange(firstBinA,lastBinA);
  hisN->GetAxis(ival)->SetRange(firstBinV,lastBinV);
  Int_t entries=0;
  for (Int_t ibin=firstBinA; ibin<lastBinA; ibin++){
    Double_t cont = hisT->GetBinContent(ibin);
    if (cont<minEntries) continue;
    entries++;
  }
  TVectorD vecTime(entries);
  TVectorD vecMean0(entries);
  TVectorD vecRMS0(entries);
  TVectorD vecMean1(entries);
  TVectorD vecRMS1(entries);
  entries=0;
  {for (Int_t ibin=firstBinA; ibin<lastBinA; ibin++){
      Double_t cont = hisT->GetBinContent(ibin);
      if (cont<minEntries) continue;
      hisN->GetAxis(itime)->SetRange(ibin-1,ibin+1);
      Double_t time = hisT->GetBinCenter(ibin);
      TH1 * his = hisN->Projection(ival);
      Double_t nentries0= his->GetBinContent(his->FindBin(0));
      if (cont-nentries0<minEntries) continue;
      //
      his->SetBinContent(his->FindBin(0),0);
      vecTime[entries]=time;
      vecMean0[entries]=his->GetMean()+offset;
      vecMean1[entries]=his->GetMeanError();
      vecRMS0[entries] =his->GetRMS();
      vecRMS1[entries] =his->GetRMSError();
      delete his;  
      entries++;
    }}
  delete hisT;
  delete hisV;
  TGraphErrors * graph =  new TGraphErrors(entries,vecTime.GetMatrixArray(), vecMean0.GetMatrixArray(),					   0, vecMean1.GetMatrixArray());
  return graph;
}








void AliTPCPreprocessorOffline::SetDefaultGraphDrift(TGraph *graph, Int_t color, Int_t style){
  //
  //
  //
  graph->GetXaxis()->SetTimeDisplay(kTRUE);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
  graph->SetMaximum( 0.025);
  graph->SetMinimum(-0.025);
  graph->GetXaxis()->SetTitle("Time");
  graph->GetYaxis()->SetTitle("v_{dcorr}");
  //
  graph->GetYaxis()->SetLabelSize(0.03);
  graph->GetXaxis()->SetLabelSize(0.03);
  //
  graph->GetXaxis()->SetNdivisions(10,5,0);
  graph->GetYaxis()->SetNdivisions(10,5,0);
  //
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetYaxis()->SetLabelOffset(0.005);
  //
  graph->GetXaxis()->SetTitleOffset(1.3);
  graph->GetYaxis()->SetTitleOffset(1.2);
  //
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetMarkerStyle(style);
}

void AliTPCPreprocessorOffline::SetPadStyle(TPad *pad, Float_t mx0, Float_t mx1, Float_t my0, Float_t my1){
  //
  //
  pad->SetTicks(1,1);
  pad->SetMargin(mx0,mx1,my0,my1);
}


void AliTPCPreprocessorOffline::MakeDefaultPlots(TObjArray * arr, TObjArray *picArray){
  //
  //
  //
  // margins
  Float_t mx0=0.12, mx1=0.1, my0=0.15, my1=0.1;
  //
  TGraphErrors* laserA       =(TGraphErrors*)arr->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
  TGraphErrors* laserC       =(TGraphErrors*)arr->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
  TGraphErrors* cosmic       =(TGraphErrors*)arr->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  TGraphErrors* cross        =(TGraphErrors*)arr->FindObject("TGRAPHERRORS_VDRIFT_CROSS_ALL");
  TGraphErrors* itstpcP       =(TGraphErrors*)arr->FindObject("ALIGN_ITSP_TPC_DRIFTVD");
  TGraphErrors* itstpcM       =(TGraphErrors*)arr->FindObject("ALIGN_ITSM_TPC_DRIFTVD");
  TGraphErrors* itstpcB       =(TGraphErrors*)arr->FindObject("ALIGN_ITSB_TPC_DRIFTVD");
  //
  if (laserA)  SetDefaultGraphDrift(laserA,2,25);
  if (laserC)  SetDefaultGraphDrift(laserC,4,26);
  if (cosmic)  SetDefaultGraphDrift(cosmic,3,27);
  if (cross)   SetDefaultGraphDrift(cross,4,28);
  if (itstpcP) SetDefaultGraphDrift(itstpcP,2,29);
  if (itstpcM) SetDefaultGraphDrift(itstpcM,4,30);
  if (itstpcB) SetDefaultGraphDrift(itstpcB,1,31);
  //
  //
  TPad *pad=0;
  //
  // Laser-Laser
  //
  if (laserA&&laserC){
    pad = new TCanvas("TPCLaserVDrift","TPCLaserVDrift");
    laserA->Draw("alp");
    SetPadStyle(pad,mx0,mx1,my0,my1);
    laserA->Draw("apl");
    laserC->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(laserA,"Laser A side");
    legend->AddEntry(laserC,"Laser C side");
    legend->Draw();    
    picArray->AddLast(pad);
  }

  if (itstpcP&&itstpcM){
    pad = new TCanvas("ITSTPC","ITSTPC");
    itstpcP->Draw("alp");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    gPad->Clear();
    itstpcM->Draw("apl");
    itstpcP->Draw("p");
    itstpcB->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(itstpcP,"ITS-TPC smooth plus");
    legend->AddEntry(itstpcM,"ITS-TPC smooth minus");
    legend->AddEntry(itstpcB,"ITS-TPC smooth ");
    legend->Draw();    
    picArray->AddLast(pad);
  }

  if (itstpcB&&laserA){
    pad = new TCanvas("ITSTPC_LASER","ITSTPC_LASER");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    laserA->Draw("alp");
    itstpcP->Draw("p");
    itstpcM->Draw("p");
    itstpcB->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(laserA,"TPC laser");
    legend->AddEntry(itstpcP,"ITS-TPC smooth plus");   
    legend->AddEntry(itstpcM,"ITS-TPC smooth minus");   
    legend->AddEntry(itstpcB,"ITS-TPC smooth ");
    legend->Draw();
    picArray->AddLast(pad);
  }

  if (itstpcP&&cross){ 
    pad = new TCanvas("ITSTPC_CROSS","ITSTPC_CROSS");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    pad->Clear();
    cross->Draw("ap");
    itstpcP->Draw("p");
    //
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");

    legend->AddEntry(cross,"TPC cross tracks");
    legend->AddEntry(itstpcB,"ITS-TPC smooth");
    legend->Draw();        
    picArray->AddLast(pad);
  }
  if (itstpcP&&cosmic){ 
    pad = new TCanvas("ITSTPC_COSMIC","ITSTPC_COSMIC");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    itstpcP->Draw("alp");
    pad->Clear();
    cosmic->Draw("ap");
    itstpcP->Draw("p");
    //
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");

    legend->AddEntry(cosmic,"TPC cross tracks0 up-down");
    legend->AddEntry(itstpcB,"ITS-TPC smooth");
    legend->Draw();        
    picArray->AddLast(pad);
  }
}




void AliTPCPreprocessorOffline::CalibTimeGain(Char_t* fileName, Int_t startRunNumber, Int_t endRunNumber,  TString  ocdbStorage){
  //
  // Update OCDB gain
  //
  ReadGainGlobal(fileName);
  AnalyzeGain(startRunNumber,endRunNumber, 1000,1.43);
  MakeQAPlot(1.43);  
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  UpdateOCDBGain( startRunNumber, endRunNumber, ocdbStorage.Data());
}




void AliTPCPreprocessorOffline::ReadGainGlobal(Char_t* fileName){
  //
  // read calibration entries from file
  // 
  TFile fcalib(fileName);
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  if (array){
    fGainMIP    = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
    fGainCosmic = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGainCosmic");
  }else{
    fGainMIP    = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGain");
    fGainCosmic = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGainCosmic");
  }
  TH1 * hisT=0;
  Int_t firstBinA =0, lastBinA=0;

  if (fGainCosmic){ 
    hisT= fGainCosmic->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    fGainCosmic->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

  if (fGainMIP){ 
    hisT= fGainMIP->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    fGainMIP->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

}



Bool_t AliTPCPreprocessorOffline::AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit,  Float_t FPtoMIPratio){
  //
  //
  //
  fGainMIP->GetHistGainTime()->GetAxis(5)->SetRangeUser(startRunNumber, endRunNumber);
  // 1.) try to create MIP spline
  fGainMIP->GetHistGainTime()->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
  fGainMIP->GetHistGainTime()->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions
  //
  fGraphMIP = AliTPCcalibBase::FitSlices(fGainMIP->GetHistGainTime(),0,1,minEntriesGaussFit,10,0.1,0.7);
  if (fGraphMIP->GetN()==0) fGraphMIP = 0x0;
  if (fGraphMIP) fFitMIP = AliTPCcalibTimeGain::MakeSplineFit(fGraphMIP);
  if (fGraphMIP) fGraphMIP->SetName("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");// set proper names according to naming convention
  fGainArray->AddAt(fFitMIP,0);
  

  // 2.) try to create Cosmic spline
  fGainCosmic->GetHistGainTime()->GetAxis(2)->SetRangeUser(0.51,1.49); // only cosmics
  fGainCosmic->GetHistGainTime()->GetAxis(4)->SetRangeUser(20,100);    // only Fermi-Plateau muons
  //
  fGraphCosmic = AliTPCcalibBase::FitSlices(fGainCosmic->GetHistGainTime(),0,1,minEntriesGaussFit,10);
  if (fGraphCosmic->GetN()==0) fGraphCosmic = 0x0;
  //
  if (fGraphCosmic) {
    for(Int_t i=0; i < fGraphCosmic->GetN(); i++) {
      fGraphCosmic->GetY()[i] /= FPtoMIPratio;
      fGraphCosmic->GetEY()[i] /= FPtoMIPratio;
    }
  }
  //
  if (fGraphCosmic) fFitCosmic = AliTPCcalibTimeGain::MakeSplineFit(fGraphCosmic);
  if (fGraphCosmic) fGraphCosmic->SetName("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL"); // set proper names according to naming convention
  fGainArray->AddAt(fFitCosmic,1);
  // with naming convention and backward compatibility
  fGainArray->AddAt(fGraphMIP,2);
  fGainArray->AddAt(fGraphCosmic,3);
  cout << "fGraphCosmic: " << fGraphCosmic << " fGraphMIP " << fGraphMIP << endl;
  return kTRUE;

}



void AliTPCPreprocessorOffline::UpdateOCDBGain(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
  //
  // Update OCDB entry
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-24-00"); //root version
  metaData->SetComment("Calibration of the time dependence of the gain due to pressure and temperature changes.");
  AliCDBId id1("TPC/Calib/TimeGain", startRunNumber, endRunNumber);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
  gStorage->Put(fGainArray, id1, metaData);    
}

void AliTPCPreprocessorOffline::MakeQAPlot(Float_t  FPtoMIPratio) {
  //
  // Make QA plot to visualize results
  //
  //
  //
  if (fGraphCosmic) {
    TCanvas * canvasCosmic = new TCanvas("gain Cosmic", "time dependent gain QA histogram cosmic");
    canvasCosmic->cd();
    TH2D * gainHistoCosmic = fGainCosmic->GetHistGainTime()->Projection(0,1);
    gainHistoCosmic->SetDirectory(0);
    gainHistoCosmic->SetName("GainHistoCosmic");
    gainHistoCosmic->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoCosmic->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoCosmic->Draw("colz");
    fGraphCosmic->SetMarkerStyle(25);
    fGraphCosmic->Draw("lp");
    fGraphCosmic->SetMarkerStyle(25);
    TGraph * grfFitCosmic = fFitCosmic->MakeGraph(fGraphCosmic->GetX()[0],fGraphCosmic->GetX()[fGraphCosmic->GetN()-1],50000,0);
    if (grfFitCosmic) {
      for(Int_t i=0; i < grfFitCosmic->GetN(); i++) {
 	grfFitCosmic->GetY()[i] *= FPtoMIPratio;	
      }
      for(Int_t i=0; i < fGraphCosmic->GetN(); i++) {
 	fGraphCosmic->GetY()[i] *= FPtoMIPratio;	
      }
    }
    fGraphCosmic->Draw("lp");
    grfFitCosmic->SetLineColor(2);
    grfFitCosmic->Draw("lu");
    fGainArray->AddLast(gainHistoCosmic);
    fGainArray->AddLast(canvasCosmic->Clone());
    delete canvasCosmic;    
  }
  if (fFitMIP) {
    TCanvas * canvasMIP = new TCanvas("gain MIP", "time dependent gain QA histogram MIP");
    canvasMIP->cd();
    TH2D * gainHistoMIP    = fGainMIP->GetHistGainTime()->Projection(0,1);
    gainHistoMIP->SetName("GainHistoCosmic");
    gainHistoMIP->SetDirectory(0);
    gainHistoMIP->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoMIP->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoMIP->Draw("colz");
    fGraphMIP->SetMarkerStyle(25);
    fGraphMIP->Draw("lp");
    TGraph * grfFitMIP = fFitMIP->MakeGraph(fGraphMIP->GetX()[0],fGraphMIP->GetX()[fGraphMIP->GetN()-1],50000,0);
    grfFitMIP->SetLineColor(2);
    grfFitMIP->Draw("lu");    
    fGainArray->AddLast(gainHistoMIP);
    fGainArray->AddLast(canvasMIP->Clone());
    delete canvasMIP;    
  }  
}



