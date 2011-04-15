// Macro to generate and update OCDB entries for a given run:
// this is a TObjArray which has at 0 the MIP-Spline and at 1 the Fermi-Plateau-Spline ...
// Responsible: marian.ivanov@cern.ch
// Responsible: A.Kalweit@gsi.de

/* How to use it locally:
//
// Load libraries
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");
.L $ALICE_ROOT/TPC/CalibMacros/CalibTimeGain.C+

//Make calibration
CalibTimeGain("CalibObjectsTrain1.root",0,120000);
TBrowser b;
b.Add(gainArray);

*/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"

#include "AliTPCcalibTimeGain.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif




TGraphErrors * graphMIP    = 0;       // graph time dependence of MIP
TGraphErrors * graphCosmic = 0;       // graph time dependence at Plateu
AliSplineFit * fitMIP = 0;            // fit of dependence - MIP
AliSplineFit * fitCosmic = 0;         // fit of dependence - Plateu
TObjArray    * gainArray = new TObjArray(4); // array to be stored in the OCDB
AliTPCcalibTimeGain * gainMIP =0;     // calibration component for MIP
AliTPCcalibTimeGain * gainCosmic =0;  // calibration component for cosmic


void UpdateOCDBGain(Int_t  startRunNumber, Int_t endRunNumber, const char* storagePath);
void ReadGainGlobal(Char_t* fileName="CalibObjectsTrain1.root");
void MakeQAPlot(Float_t  FPtoMIPratio);
Bool_t AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit = 500, Float_t FPtoMIPratio = 1.43); 



void CalibTimeGain(Char_t* fileName="CalibObjectsTrain1.root", Int_t startRunNumber=0, Int_t endRunNumber=AliCDBRunRange::Infinity(),  TString  ocdbStorage=""){
  //
  // Update OCDB gain
  //
  ReadGainGlobal(fileName);
  AnalyzeGain(startRunNumber,endRunNumber, 1000,1.43);
  MakeQAPlot(1.43);  
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  UpdateOCDBGain( startRunNumber, endRunNumber, ocdbStorage.Data());
}




void ReadGainGlobal(Char_t* fileName){
  //
  // read calibration entries from file
  // 
  TFile fcalib(fileName);
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  if (array){
    gainMIP    = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
    gainCosmic = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGainCosmic");
  }else{
    gainMIP    = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGain");
    gainCosmic = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGainCosmic");
  }
  TH1 * hisT=0;
  Int_t firstBinA =0, lastBinA=0;

  if (gainCosmic){ 
    hisT= gainCosmic->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    gainCosmic->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

  if (gainMIP){ 
    hisT= gainCosmic->GetHistGainTime()->Projection(1);
    firstBinA = hisT->FindFirstBinAbove(2);
    lastBinA  = hisT->FindLastBinAbove(2);    
    gainMIP->GetHistGainTime()->GetAxis(1)->SetRange(firstBinA,lastBinA);
    delete hisT;
  }

}



Bool_t AnalyzeGain(Int_t startRunNumber, Int_t endRunNumber, Int_t minEntriesGaussFit,  Float_t FPtoMIPratio){
  //
  //
  //
  gainMIP->GetHistGainTime()->GetAxis(5)->SetRangeUser(startRunNumber, endRunNumber);
  // 1.) try to create MIP spline
  gainMIP->GetHistGainTime()->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
  gainMIP->GetHistGainTime()->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions
  //
  graphMIP = AliTPCcalibBase::FitSlices(gainMIP->GetHistGainTime(),0,1,minEntriesGaussFit,10,0.1,0.7);
  if (graphMIP->GetN()==0) graphMIP = 0x0;
  if (graphMIP) fitMIP = AliTPCcalibTimeGain::MakeSplineFit(graphMIP);
  if (graphMIP) graphMIP->SetName("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");// set proper names according to naming convention
  gainArray->AddAt(fitMIP,0);
  

  // 2.) try to create Cosmic spline
  gainCosmic->GetHistGainTime()->GetAxis(2)->SetRangeUser(0.51,1.49); // only cosmics
  gainCosmic->GetHistGainTime()->GetAxis(4)->SetRangeUser(20,100);    // only Fermi-Plateau muons
  //
  graphCosmic = AliTPCcalibBase::FitSlices(gainCosmic->GetHistGainTime(),0,1,minEntriesGaussFit,10);
  if (graphCosmic->GetN()==0) graphCosmic = 0x0;
  //
  if (graphCosmic) {
    for(Int_t i=0; i < graphCosmic->GetN(); i++) {
      graphCosmic->GetY()[i] /= FPtoMIPratio;
      graphCosmic->GetEY()[i] /= FPtoMIPratio;
    }
  }
  //
  if (graphCosmic) fitCosmic = AliTPCcalibTimeGain::MakeSplineFit(graphCosmic);
  if (graphCosmic) graphCosmic->SetName("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL"); // set proper names according to naming convention
  gainArray->AddAt(fitCosmic,1);
  // with naming convention and backward compatibility
  gainArray->AddAt(graphMIP,2);
  gainArray->AddAt(graphCosmic,3);
  cout << "graphCosmic: " << graphCosmic << " graphMIP " << graphMIP << endl;
  return kTRUE;

}



void UpdateOCDBGain(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
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
  gStorage->Put(gainArray, id1, metaData);    
}

void MakeQAPlot(Float_t  FPtoMIPratio) {
  //
  // Make QA plot to visualize results
  //
  //
  //
  if (graphCosmic) {
    TCanvas * canvasCosmic = new TCanvas("gain Cosmic", "time dependent gain QA histogram cosmic");
    canvasCosmic->cd();
    TH2D * gainHistoCosmic = gainCosmic->GetHistGainTime()->Projection(0,1);
    gainHistoCosmic->SetDirectory(0);
    gainHistoCosmic->SetName("GainHistoCosmic");
    gainHistoCosmic->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoCosmic->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoCosmic->Draw("colz");
    graphCosmic->SetMarkerStyle(25);
    graphCosmic->Draw("lp");
    graphCosmic->SetMarkerStyle(25);
    TGraph * grfitCosmic = fitCosmic->MakeGraph(graphCosmic->GetX()[0],graphCosmic->GetX()[graphCosmic->GetN()-1],50000,0);
    if (grfitCosmic) {
      for(Int_t i=0; i < grfitCosmic->GetN(); i++) {
 	grfitCosmic->GetY()[i] *= FPtoMIPratio;	
      }
      for(Int_t i=0; i < graphCosmic->GetN(); i++) {
 	graphCosmic->GetY()[i] *= FPtoMIPratio;	
      }
    }
    graphCosmic->Draw("lp");
    grfitCosmic->SetLineColor(2);
    grfitCosmic->Draw("lu");
    gainArray->AddLast(gainHistoCosmic);
    gainArray->AddLast(canvasCosmic->Clone());
    delete canvasCosmic;    
  }
  if (fitMIP) {
    TCanvas * canvasMIP = new TCanvas("gain MIP", "time dependent gain QA histogram MIP");
    canvasMIP->cd();
    TH2D * gainHistoMIP    = gainMIP->GetHistGainTime()->Projection(0,1);
    gainHistoMIP->SetName("GainHistoCosmic");
    gainHistoMIP->SetDirectory(0);
    gainHistoMIP->GetXaxis()->SetTimeDisplay(kTRUE);
    gainHistoMIP->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    gainHistoMIP->Draw("colz");
    graphMIP->SetMarkerStyle(25);
    graphMIP->Draw("lp");
    TGraph * grfitMIP = fitMIP->MakeGraph(graphMIP->GetX()[0],graphMIP->GetX()[graphMIP->GetN()-1],50000,0);
    grfitMIP->SetLineColor(2);
    grfitMIP->Draw("lu");    
    gainArray->AddLast(gainHistoMIP);
    gainArray->AddLast(canvasMIP->Clone());
    delete canvasMIP;    
  }  
}


