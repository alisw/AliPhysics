/* small macro to generate and update OCDB entries for a given run:

// How to use it...
//
// 1. go to the directory with the "CalibObjectsTrain1.root" file and start your aliroot
//
// 2. copy and paste the following commands
//
gSystem->Load("libSTEER");
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");
gSystem->AddIncludePath("-I$ALICE_ROOT/STEER");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
.L $ALICE_ROOT/TPC/CalibMacros/CalibTimeVdrift.C+

//
// 3. setup your OCDB output path e.g:
//
ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
//
// 4. simple execution
//
CalibTimeVdriftGlobal()
//
// 5. try to visualize new entry
//
Int_t run=90000; //
//.x ConfigOCDB.C(run);

ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB")
AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/TimeDrift",ocdbStorage.Data());

AliCDBEntry* entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeDrift",run)
TObjArray * arr = (TObjArray*)entry->GetObject();

TObjArray *picArray = new TObjArray;
MakeDefaultPlots(arr,picArray);

//
//  7. Tricky part - hopefully not neccessary
//     change the default time 0
//     BUG in OCDB - we can not change the SpecificStorage once we read given
//                   entry

   .x ConfigOCDB.C(

*/

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

TString ocdbStorage="dummy";
//
//
//
void AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift, Int_t minEntries);
void AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift);
void AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift);
//

TGraphErrors* FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY);
TGraphErrors* FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs);




void CalibTimeVdriftGlobal(Char_t* file="CalibObjectsTrain1.root"){

  const Int_t    kMinEntries=500;     // minimal number of entries
  //
  // 1. Initialization and run range setting
  TFile fcalib(file);
  AliTPCcalibTime* timeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");
  Int_t startRun=0, endRun=0;
  Int_t startTime=0, endTime=0;
  //Int_t startPT=-1, endPT=-1;
  //
  // find the fist and last run
  //
  TObjArray *hisArray =timeDrift->GetHistoDrift();
  {for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
    THnSparse* addHist=(THnSparse*)hisArray->UncheckedAt(i);
    if (addHist->GetEntries()<kMinEntries) continue;
    if (!addHist) continue;
    TH1D* histo    =addHist->Projection(3);
    TH1D* histoTime=addHist->Projection(0);
    if (startRun==0){ 
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
  {for (Int_t i=0; i<hisArray->GetEntriesFast(); i++){
      THnSparse* addHist=(THnSparse*)hisArray->At(i);
      if (!addHist) continue;
      addHist->GetAxis(0)->SetRange(startTime-1,endTime+1);
      addHist->GetAxis(3)->SetRange(startRun-1,endRun+1);
    }
  }
  printf("Run range  :\t%d-%d\n", startRun, endRun);
  printf("Time range :\t%d-%d\n", startTime, endTime);
  //
  // 2. extraction of the information
  //
  TObjArray * vdriftArray = new TObjArray();
  AddAlignmentGraphs(vdriftArray,timeDrift);
  AddHistoGraphs(vdriftArray,timeDrift,kMinEntries);
  AddLaserGraphs(vdriftArray,timeDrift);
  //
  //
  // 3. update of OCDB
  //
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Dag Toppe Larsen");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Calibration of the time dependence of the drift velocity due to pressure and temperature changes");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/TimeDrift", startRun, AliCDBRunRange::Infinity());
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(vdriftArray, (*id1), metaData);
 
}


void UpdateDriftParam(AliTPCParam *param, TObjArray *arr, Int_t startRun){
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

void CalibTimeVdrift(Int_t runNumber, AliTPCcalibTime* vdrift, Int_t end=false){
  TObjArray* vdriftArray = new TObjArray();
  TObjArray * array=vdrift->GetHistoDrift();
  if(!array) return;
  TIterator* iterator = array->MakeIterator();
  iterator->Reset();
  THnSparse* hist=NULL;
  while((hist=(THnSparseF*)iterator->Next())){
    if(!hist) continue;
    hist->Print();
    if(end) hist->GetAxis(3)->SetRangeUser(runNumber-0.5,       end+0.5);
    else    hist->GetAxis(3)->SetRangeUser(runNumber-0.5, runNumber+0.5);
    TString name=hist->GetName();
    Int_t dim[4]={0,1,2,3};
    THnSparse* newHist=hist->Projection(4,dim);
    newHist->SetName(name);
    vdriftArray->Add(newHist);
    TGraphErrors* graph=AliTPCcalibBase::FitSlices(newHist,2,0,400,100,0.05,0.95, kTRUE);
    printf("name=%s graph=%i\n", name.Data(), graph==0);
    if(!graph || !graph->GetN()) continue;
    printf("name=%s graph=%i, N=%i\n", name.Data(), graph==0, graph->GetN());
    Int_t pos=name.Index("_");
    name=name(pos,name.Capacity()-pos);
    TString graphName=graph->ClassName();
    graphName+=name;
    graphName.ToUpper();
    graph->SetName(graphName);
    printf("name=%s\n", graphName.Data());
    vdriftArray->Add(graph);
//      AliSplineFit* fit=new AliSplineFit();
//      fit->SetGraph(graph);
//      fit->SetMinPoints(graph->GetN()+1);
//      fit->InitKnots(graph,2,0,0.001);
//      fit->SplineFit(0);
//      TString fiName=fit->ClassName();
//      fiName+=type;
//      fiName+=trigger;
//      fiName.ToUpper();
//      fit->SetName(fiName.Data());
//      printf("name=%s\n", fiName.Data());
//      vdriftArray->Add(fit);
  }
  THnSparse* laserHist=NULL;
  TGraphErrors* laserGraph=NULL;
  TString laserName="";

  laserHist=vdrift->GetHistVdriftLaserA(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DRIFT_LASER_ALL_A";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  vdriftArray->Add(laserHist);
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,400,100,0.05,0.95, kTRUE);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DRIFT_LASER_ALL_A";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  laserHist=vdrift->GetHistVdriftLaserC(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DRIFT_LASER_ALL_C";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  vdriftArray->Add(laserHist);
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,400,100,0.05,0.95, kTRUE);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DRIFT_LASER_ALL_C";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Dag Toppe Larsen");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Calibration of the time dependence of the drift velocity due to pressure and temperature changes");
  AliCDBId* id1=NULL;
  if(end) id1=new AliCDBId("TPC/Calib/TimeDrift", runNumber, end);
  else    id1=new AliCDBId("TPC/Calib/TimeDrift", runNumber, runNumber);
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(vdriftArray, (*id1), metaData);
  printf("done runNumber=%i, end=%i\n", runNumber, end);
}

void CalibrateAll(Char_t* file="CalibObjectsTrain1.root"){
  TFile fcalib(file);
  AliTPCcalibTime* timeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");
  THnSparse* addHist=timeDrift->GetHistoDrift("all");
  TH1D* histo=addHist->Projection(3);
  Int_t start=histo->FindFirstBinAbove(0);
  Int_t end  =histo->FindLastBinAbove(0);
  CalibTimeVdrift(start-1, timeDrift, end-1);
  for(Int_t i=start;i<end;i++){
    if(histo->GetBinContent(i)){
      printf("start=%i i=%i end=%i i=%i center=%f", start, i, end, i, histo->GetBinCenter(i));
      CalibTimeVdrift(i, timeDrift);
    }
  }
}

void MakePlot(Char_t* file="CalibObjects.root", Char_t* trigger="", /*Int_t useError=true,*/ Int_t separateFits=false, Int_t pointLimit=0, Int_t showHisto=true, Int_t otherColour=false, Int_t spline=true, Double_t yMin=0.0, Double_t yMax=0.0){
  TFile fcalib(file);
  AliTPCcalibTime* timeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");

//Get graphs
  TObjArray* grArray=new TObjArray();
  TObjArray* hiArray=new TObjArray();
  Int_t counter=0;
  Int_t pos=counter+1;
  Int_t colour=2;
  TObjArray * addArray=timeDrift->GetHistoDrift();
  if(!addArray) return;
  TIterator* iterator = addArray->MakeIterator();
  iterator->Reset();
  THnSparse* addHist=NULL;

  THnSparse* laserHist=NULL;
  TString laserName="";

  laserHist=timeDrift->GetHistVdriftLaserA(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_VDRIFT_LASER_ALL_A";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  addArray->Add(laserHist);

  laserHist=timeDrift->GetHistVdriftLaserC(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_VDRIFT_LASER_ALL_C";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  addArray->Add(laserHist);

  while((addHist=(THnSparseF*)iterator->Next())){
    if (!addHist) continue;
    TString histName=addHist->GetName();
    if(!histName.Contains(trigger) && !histName.Contains("COSMICS_ALL")) continue;
    addHist->Print();
    TGraph* graph=AliTPCcalibBase::FitSlices(addHist,2,0,400,100,0.05,0.95, kTRUE);
    TH2D* histo=addHist->Projection(2,0);
    if(!graph || !histo || !graph->GetN()) continue;
    if(spline){
      AliSplineFit fit;
      fit.SetGraph(graph);
      fit.SetMinPoints(graph->GetN()+1);
      fit.InitKnots(graph,2,0,0.001);
      fit.SplineFit(0);
      graph=fit.MakeGraph(graph->GetX()[0],graph->GetX()[graph->GetN()-1],50000,0);
    }
    graph->SetNameTitle(addHist->GetName(),"Relative drift velocity change");
    graph->GetXaxis()->SetTitle("time");
    graph->GetXaxis()->SetTimeFormat();
    graph->GetXaxis()->SetTimeDisplay(1);
    graph->GetYaxis()->SetTitleOffset(1.2);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitle("relative drift velocity change");
    if(yMin) graph->SetMinimum(yMin);
    if(yMax) graph->SetMaximum(yMax);
    printf("bins=%i name=%s couter=%i pos=%i colour=%i Ymax=%f Ymin=%f Xmax=%f Xmin=%f\n", graph->GetXaxis()->GetNbins(), addHist->GetName(), counter, pos, colour, graph->GetYaxis()->GetXmax(), graph->GetYaxis()->GetXmin(), graph->GetXaxis()->GetXmax(), graph->GetXaxis()->GetXmin());
    if(graph->GetXaxis()->GetNbins()<=pointLimit) continue;
//    if(!useError&&!spline) for(Int_t j=0; j<graph->GetXaxis()->GetNbins(); j++) graph->SetPointError(j,0,0);
//    if( *((Int_t*)(graph->GetName()))==*((Int_t*)("all")) ){
    TString name=graph->GetName();
    if(name.Contains("COSMICS_ALL")){
      graph->SetLineColor(0+1);
      graph->SetMarkerColor(0+1);
      grArray->AddAtAndExpand(graph,0);
      histo->SetLineColor(1);
      histo->SetMarkerColor(1);
      hiArray->AddAtAndExpand(histo,0);
      pos=counter;
    }
    else{
      graph->SetLineColor(colour);
      graph->SetMarkerColor(colour);
      grArray->AddAtAndExpand(graph,pos);
      histo->SetLineColor(colour);
      histo->SetMarkerColor(colour);
      hiArray->AddAtAndExpand(histo,pos);
      colour++;
      if(colour==10) colour++;
    }
    pos++;
    counter++;
    graph=0;
    addHist=0;
  }

  if(!separateFits){
    TCanvas* plot = new TCanvas("plot", "Relative drift velocity change", 1500, 900);
    //Draw graphs
    TGraph* grTemp=0;
    TH2D* hiTemp=0;
    for(Int_t i=0; i<counter; i++){
      if(!(*grArray)[i] || !(*hiArray)[i]) continue;
      grTemp=(TGraph*)(*grArray)[i];
      hiTemp=(TH2D*)(*hiArray)[i];
      if(i==0) grTemp->Draw("alp");
      else grTemp->Draw("lp");
      if(showHisto) hiTemp->Draw("same");
    }
    //Add legends
     TLegend* leg=NULL;
     if(pointLimit){
       leg = new TLegend(0.4,0.75,1,0.95,NULL,"brNDC");
       leg->SetTextSize(0.03);
     }
     else{
       leg = new TLegend(0.4,0.65,1,1,NULL,"brNDC");
       leg->SetTextSize(0.02);
     }
    leg->SetLineColor(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    TLegendEntry* entry=0;
    for(Int_t i=0; i<counter; i++){
      if(!(*grArray)[i]) continue;
      TString name=(*grArray)[i]->GetName();
      Int_t pos=name.Index("_")+1+12;
      name=name(pos,name.Capacity()-pos);
      entry = leg->AddEntry((*grArray)[i],name,"l");
      leg->Draw();
    }
    plot->SaveAs("relative_drift_velocity_change.pdf");
    plot->SaveAs("relative_drift_velocity_change.png");
  }
  else{
    //Draw graphs
    TGraph* grTemp=0;
    TH2D* hiTemp=0;
    for(Int_t i=0; i<counter; i++){
      if(!(*grArray)[i] || !(*hiArray)[i]) continue;
      grTemp=(TGraph*)(*grArray)[i];
      hiTemp=(TH2D*)(*hiArray)[i];
      TString name=grTemp->GetName();
      Int_t pos=name.Index("_")+1+12;
      name=name(pos,name.Capacity()-pos);
      grTemp->SetTitle(name);
      if(otherColour){
        grTemp->SetLineColor(2);
        grTemp->SetMarkerColor(2);
        hiTemp->SetLineColor(1);
        hiTemp->SetMarkerColor(1);
      }
      TCanvas* plot = new TCanvas(name, name, 1500, 900);
      grTemp->Draw("alp");
      if(showHisto) hiTemp->Draw("same");
      TString name2=name;
      name+=".pdf";
      name2+=".png";
      plot->SaveAs(name);
      plot->SaveAs(name2);
    }
  }
}

void MakeQaPlot(Char_t* file="CalibObjects.root"){
  TFile fcalib(file);
  AliTPCcalibTime* timeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");
  TCanvas* plotQa = new TCanvas("plotQa", "Relative drift velocity change QA", 1600, 600);
  plotQa->Divide(5,2);
  for(Int_t i=0; i<10; i++){
    plotQa->cd(i+1);
    timeDrift->GetCosmiMatchingHisto(i)->Draw();
  }
  plotQa->cd(0);
  plotQa->SaveAs("relative_drift_velocity_change_qa.pdf");
}

void MakeHistPlot(Char_t* file="CalibObjects.root", Char_t* trigger="COSMICS_ALL"){
  TFile fcalib(file);
  AliTPCcalibTime* timeDrift=(AliTPCcalibTime*)fcalib.Get("calibTime");
  THnSparse* addHist=timeDrift->GetHistoDrift(trigger);
  TH2D* histo=addHist->Projection(2,0);

  histo->SetTitle("Relative drift velocity change");
  histo->GetXaxis()->SetTitle("time");
  histo->GetXaxis()->SetTimeFormat();
  histo->GetXaxis()->SetTimeDisplay(1);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->GetYaxis()->SetTitle("relative drift velocity change");
  histo->SetLineColor(1);
  histo->SetMarkerColor(1);
  TGraphErrors* graph=AliTPCcalibBase::FitSlices(addHist,2,0,400,100,0.05,0.95, kTRUE);
  if(!graph) return;
  graph->SetNameTitle(addHist->GetName(),"Relative drift velocity change");
  graph->GetXaxis()->SetTitle("time");
  graph->GetXaxis()->SetTimeFormat();
  graph->GetXaxis()->SetTimeDisplay(1);
  graph->GetYaxis()->SetTitleOffset(1.3);  
  graph->GetYaxis()->SetTitle("relative drift velocity change");
  graph->SetMaximum(0.003);
  graph->SetMinimum(-0.003);
  graph->SetLineColor(2);
  graph->SetMarkerColor(2);
  TCanvas* histo_fit = new TCanvas("histo_fit", "Relative drift velocity change histogram/fit", 800, 600);
  histo->Draw();
  graph->Draw("lp");
  histo_fit->SaveAs("relative_drift_velocity_change_histo_fit.pdf");
  histo_fit->SaveAs("relative_drift_velocity_change_histo_fit.png");
}



void PrintArray(TObjArray *array){
  //
  //
  //
  Int_t entries = array->GetEntries();
  for (Int_t i=0; i<entries; i++){
    printf("%s\n", array->At(i)->GetName());
  }
}


TMultiGraph * MakeMultiGraph(TObjArray *array, const char  *mask, Int_t color, Int_t style){
  //
  // Make a multi graph with graphs selected
  //
  if (!array) return 0;
  TMultiGraph *mgraph=new TMultiGraph("MultiGraph", "MultiGraph");
  TObjArray * maskArray=0;
  if (mask){
    TString grmask(mask);
    maskArray = grmask.Tokenize("*");
  }
  Double_t ymin=0,ymax=-1;
  for (Int_t i=0; i<array->GetEntries();i++){
    TGraphErrors *graph= (TGraphErrors*)array->At(i);
    if (!graph) continue;
    if (maskArray){
      Bool_t isOK=kTRUE;
      TString str(graph->GetName());
      for (Int_t imask=0; imask<maskArray->GetEntries();imask++)
	if (!str.Contains(maskArray->At(imask)->GetName())==0) isOK=kFALSE;
      if (!isOK) continue;
    }
    if (ymax<ymin){
      ymin=graph->GetMinimum();
      ymax=graph->GetMaximum();
    }
    ymin = TMath::Min(ymin,graph->GetMinimum());
    ymax = TMath::Max(ymax,graph->GetMaximum());
  }
  ymin=ymin-(ymax-ymin)*0.2;
  ymax=ymax+(ymax-ymin)*0.5;
  Bool_t first=kTRUE;
  for (Int_t i=0; i<array->GetEntries();i++){
    TGraphErrors *graph= (TGraphErrors*)array->At(i);
    if (!graph) continue;
    if (maskArray){
      Bool_t isOK=kTRUE;
      TString str(graph->GetName());
      for (Int_t imask=0; imask<maskArray->GetEntries();imask++)
	if (str.Contains(maskArray->At(imask)->GetName())==0) isOK=kFALSE;
      if (!isOK) continue;
    }
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(style);
    graph->SetLineColor(color);
    graph->SetMaximum(ymax);
    graph->SetMinimum(ymin);
    if (first)  {first=kFALSE;}
    mgraph->Add(graph);
  }
  return mgraph;
}


TGraphErrors* FilterGraphDrift(TGraphErrors * graph, Float_t errSigmaCut, Float_t medianCutAbs){
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



TGraphErrors* FilterGraphMedianAbs(TGraphErrors * graph, Float_t cut,Double_t &medianY){
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


void AddHistoGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift, Int_t minEntries){
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




void AddAlignmentGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift){
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
  mstatITS= AliTPCcalibDButil::MakeStatRelKalman(arrayITS,0.9,50,0.025);
  mstatTOF= AliTPCcalibDButil::MakeStatRelKalman(arrayTOF,0.9,1000,0.025);
  mstatTRD= AliTPCcalibDButil::MakeStatRelKalman(arrayTRD,0.9,50,0.025);
  //
  TObjArray * arrayITSP= AliTPCcalibDButil::SmoothRelKalman(arrayITS,*mstatITS, 0, 5.);
  TObjArray * arrayITSM= AliTPCcalibDButil::SmoothRelKalman(arrayITS,*mstatITS, 1, 5.);
  TObjArray * arrayITSB= AliTPCcalibDButil::SmoothRelKalman(arrayITSP,arrayITSM);
  TObjArray * arrayTOFP= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,*mstatTOF, 0, 5.);
  TObjArray * arrayTOFM= AliTPCcalibDButil::SmoothRelKalman(arrayTOF,*mstatTOF, 1, 5.);
  TObjArray * arrayTOFB= AliTPCcalibDButil::SmoothRelKalman(arrayTOFP,arrayTOFM);

  TObjArray * arrayTRDP= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,*mstatTRD, 0, 5.);
  TObjArray * arrayTRDM= AliTPCcalibDButil::SmoothRelKalman(arrayTRD,*mstatTRD, 1, 5.);
  TObjArray * arrayTRDB= AliTPCcalibDButil::SmoothRelKalman(arrayTRDP,arrayTRDM);

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



void AddLaserGraphs(  TObjArray * vdriftArray, AliTPCcalibTime *timeDrift){
  //
  // Add laser graphs
  //
  THnSparse* laserHist=NULL;
  TGraphErrors* laserGraph=NULL;
  TString laserName="";
  laserHist=timeDrift->GetHistVdriftLaserA(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DRIFT_LASER_ALL_A";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,300,10,0.05,0.95, kTRUE);
  laserGraph=FilterGraphDrift(laserGraph, 5,0.04);
  //
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DRIFT_LASER_ALL_A";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  laserHist=NULL;
  laserGraph=NULL;
  laserName="";

  laserHist=timeDrift->GetHistVdriftLaserC(1);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DRIFT_LASER_ALL_C";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  //vdriftArray->Add(laserHist);// really add the whole histogram ??
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,300,10,0.05,0.95, kTRUE);
  laserGraph = FilterGraphDrift(laserGraph, 5,0.04);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DRIFT_LASER_ALL_C";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  // 2.c) additional informaion from laser for A and C side : DELAY
  laserHist=NULL;
  laserGraph=NULL;
  laserName="";
  
  laserHist=timeDrift->GetHistVdriftLaserA(0);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DELAY_LASER_ALL_A";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,100,10,0.05,0.95, kTRUE);
  
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DELAY_LASER_ALL_A";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  laserHist=NULL;
  laserGraph=NULL;
  laserName="";

  laserHist=timeDrift->GetHistVdriftLaserC(0);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_DELAY_LASER_ALL_C";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  //vdriftArray->Add(laserHist);// really add the whole histogram ??
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,100,10,0.05,0.95, kTRUE);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_DELAY_LASER_ALL_C";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  // 2.d) additional informaion from laser for A and C side : GlobalY Gradient
  laserHist=NULL;
  laserGraph=NULL;
  laserName="";

  laserHist=timeDrift->GetHistVdriftLaserA(2);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_GLOBALYGRADIENT_LASER_ALL_A";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  //vdriftArray->Add(laserHist);  // really add the whole histogram ??
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,100,10,0.05,0.95, kTRUE);
  laserGraph=FilterGraphDrift(laserGraph, 5,0.04);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_GLOBALYGRADIENT_LASER_ALL_A";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }

  laserHist=NULL;
  laserGraph=NULL;
  laserName="";

  laserHist=timeDrift->GetHistVdriftLaserC(2);
  laserName=laserHist->ClassName();
  laserName+="_MEAN_GLOBALYGRADIENT_LASER_ALL_C";
  laserName.ToUpper();
  laserHist->SetName(laserName);
  laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,100,10,0.05,0.95, kTRUE);
  laserGraph=FilterGraphDrift(laserGraph, 5,0.04);
  if(laserGraph && laserGraph->GetN()){
    laserName=laserGraph->GetName();
    laserName+="_MEAN_GLOBALYGRADIENT_LASER_ALL_C";
    laserName.ToUpper();
    laserGraph->SetName(laserName);
    vdriftArray->Add(laserGraph);
  }
}

void SetDefaultGraphDrift(TGraph *graph, Int_t color, Int_t style){
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

void SetPadStyle(TPad *pad, Float_t mx0, Float_t mx1, Float_t my0, Float_t my1){
  //
  //
  pad->SetTicks(1,1);
  pad->SetMargin(mx0,mx1,my0,my1);
}


void MakeDefaultPlots(TObjArray * arr, TObjArray *picArray){
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
  //
  if (laserA)  SetDefaultGraphDrift(laserA,2,25);
  if (laserC)  SetDefaultGraphDrift(laserC,4,26);
  if (cosmic)  SetDefaultGraphDrift(cosmic,3,27);
  if (cross)   SetDefaultGraphDrift(cross,4,28);
  if (itstpcP) SetDefaultGraphDrift(itstpcP,2,29);
  if (itstpcM) SetDefaultGraphDrift(itstpcM,4,30);
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
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(itstpcP,"ITS-TPC smooth forward");
    legend->AddEntry(itstpcM,"ITS-TPC smooth back");
    legend->Draw();    
    picArray->AddLast(pad);
  }

  if (itstpcP&&laserA){
    pad = new TCanvas("ITSTPC_LASER","ITSTPC_LASER");
    SetPadStyle(pad,mx0,mx1,my0,my1);    
    laserA->Draw("alp");
    itstpcM->Draw("p");
    TLegend *legend = new TLegend(mx0+0.01,1-my1-0.2, 0.5, 1-my1-0.01, "Drift velocity correction");
    legend->AddEntry(laserA,"TPC laser");
    legend->AddEntry(itstpcM,"ITS-TPC smooth back");
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
    legend->AddEntry(itstpcP,"ITS-TPC smooth back");
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
    legend->AddEntry(itstpcP,"ITS-TPC smooth back");
    legend->Draw();        
    picArray->AddLast(pad);
  }
}


TObjArray * SelectEntries(TObjArray* array, char * mask){
  //
  //
  //
  TObjArray * arraySelect = new TObjArray;
  TObjArray * maskArray=0;
  if (mask){
    TString grmask(mask);
    maskArray = grmask.Tokenize("*");
  }
  for (Int_t i=0; i<array->GetEntries();i++){
    TObject *graph= array->At(i);
    if (!graph) continue;
    if (maskArray){
      Bool_t isOK=kTRUE;
      TString str(graph->GetName());
      for (Int_t imask=0; imask<maskArray->GetEntries();imask++)
	if (str.Contains(maskArray->At(imask)->GetName())==0) isOK=kFALSE;
      if (!isOK) continue;
      if (isOK) {
	printf("%s\n",graph->GetName());
	arraySelect->AddLast(graph->Clone());
      }
    }
  }
  return arraySelect;
}


void Substract(TGraphErrors *refgr, TObjArray *arraySelect){
  //
  //
  /*
  TObjArray *arraySelect = SelectEntries(arr,"COSMIC");
  TGraphErrors *refgr= (TGraphErrors*)arr->FindObject("ALIGN_ITSP_TPC_DRIFTVD");
  Substract(refgr,arraySelect);
  */
  Float_t mx0=0.12, mx1=0.1, my0=0.15, my1=0.5;
  {
  for (Int_t igr=0; igr<arraySelect->GetEntries(); igr++){
    TGraphErrors *gr = (TGraphErrors *)arraySelect->At(igr);
    Int_t ngr = gr->GetN();
     for (Int_t i=0; i<ngr;i++){
      gr->GetY()[i]-=refgr->Eval(gr->GetX()[i]);
      gr->GetY()[i]*=250/0.264;
    }
    gr->GetXaxis()->SetRangeUser(refgr->GetXaxis()->GetXmin(),refgr->GetXaxis()->GetXmax());
    SetDefaultGraphDrift(gr,((igr)),20+igr);
    gr->SetMaximum(8.);
    gr->SetMinimum(4.);
    gr->GetYaxis()->SetTitle("#Delta_{t}(Time Bin)");
  }
  //
  TPad * pad = new TCanvas("Delays","Delays",1000,800);
  SetPadStyle(pad,mx0,mx1,my0,my1);    
  TLegend *legend = new TLegend(mx0+0.01,1-my1, 1-mx1, 1-0.01, "Time Offset");
  {
    arraySelect->At(0)->Draw("ap");
    for (Int_t igr=0; igr<arraySelect->GetEntries(); igr++){
      TGraphErrors *gr = (TGraphErrors *)arraySelect->At(igr);
      if (gr->GetN()<5) continue;
      if (gr) gr->Draw("p");
      legend->AddEntry(gr,gr->GetName());
    }
  }
  legend->Draw();
  }
}
