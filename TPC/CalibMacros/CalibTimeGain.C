
// small macro to generate and update OCDB entries for a given run:
//
// this is a TObjArray which has at 0 the MIP-Spline and at 1 the Fermi-Plateau-Spline ...

/* How to use it...

gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");

.L $ALICE_ROOT/TPC/CalibMacros/CalibTimeGain.C+

TFile fcalib("CalibObjectsTrain1.root");
//TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
//AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)array->FindObject("calibTimeGain");
AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGain")

ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
CalibTimeGain(gain, ocdbStorage.Data(),70000,100000,kTRUE)

*/

#include "TObjArray.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliTPCcalibTimeGain.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"


Bool_t CalibTimeGain(AliTPCcalibTimeGain * gain, Char_t * storagePath, Int_t startRunNumber, Int_t endRunNumber, Bool_t updateOCDB = kFALSE, Int_t minEntriesGaussFit = 500, Bool_t makeQAplot=kTRUE){

  TObjArray * splineArray = new TObjArray(4);
  gain->GetHistGainTime()->GetAxis(5)->SetRangeUser(startRunNumber, endRunNumber);

  // 1.) try to create MIP spline
  gain->GetHistGainTime()->GetAxis(2)->SetRangeUser(1.51,2.49); // only beam data
  gain->GetHistGainTime()->GetAxis(4)->SetRangeUser(0.39,0.51); // only MIP pions
  //
  TGraphErrors * graphMIP = AliTPCcalibBase::FitSlices(gain->GetHistGainTime(),0,1,minEntriesGaussFit,10);
  if (graphMIP->GetN()==0) graphMIP = 0x0;
  AliSplineFit * fitMIP = 0;
  if (graphMIP) fitMIP = AliTPCcalibTimeGain::MakeSplineFit(graphMIP);
  if (graphMIP) graphMIP->SetName("TGRAPHERRORS_MEAN_GAIN_BEAM_ALL");// set proper names according to naming convention
  (*splineArray)[0] = fitMIP;
  
  // 2.) try to create Cosmic spline
  gain->GetHistGainTime()->GetAxis(2)->SetRangeUser(0.51,1.49); // only cosmics
  gain->GetHistGainTime()->GetAxis(4)->SetRangeUser(20,100);    // only Fermi-Plateau muons
  //
  TGraphErrors * graphCosmic = AliTPCcalibBase::FitSlices(gain->GetHistGainTime(),0,1,minEntriesGaussFit,10);
  if (graphCosmic->GetN()==0) graphCosmic = 0x0;
  AliSplineFit * fitCosmic = 0;
  if (graphCosmic) fitCosmic = AliTPCcalibTimeGain::MakeSplineFit(graphCosmic);
  if (graphCosmic) graphCosmic->SetName("TGRAPHERRORS_MEAN_GAIN_COSMIC_ALL"); // set proper names according to naming convention
  (*splineArray)[1] = fitCosmic;

  // with naming convention and backward compatibility
  (*splineArray)[2] = graphMIP;
  (*splineArray)[3] = graphCosmic;
  cout << "graphCosmic: " << graphCosmic << " graphMIP " << graphMIP << endl;
  //
  // store OCDB entry
  //
  if (!fitCosmic && !fitMIP) return kFALSE;
  if (!graphCosmic && !graphMIP) return kFALSE;
  //
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-24-00"); //root version
  metaData->SetComment("Calibration of the time dependence of the gain due to pressure and temperature changes.");
  AliCDBId id1("TPC/Calib/TimeGain", startRunNumber, endRunNumber);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
  if (updateOCDB) gStorage->Put(splineArray, id1, metaData);

  if (makeQAplot) {
    TCanvas * canvQA = new TCanvas("canvQA", "time dependent gain QA histogram");
    canvQA->cd();
    TGraphErrors * gr = gain->GetGraphGainVsTime(0,minEntriesGaussFit);
    TH2D * GainTime = gain->GetHistGainTime()->Projection(0,1);
    GainTime->GetXaxis()->SetTimeDisplay(kTRUE);
    GainTime->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    GainTime->Draw("colz");
    //
    if (graphCosmic) {
      graphCosmic->SetMarkerStyle(25);
      graphCosmic->Draw("lp");
      TGraph * grfitCosmic = fitCosmic->MakeGraph(graphCosmic->GetX()[0],graphCosmic->GetX()[graphCosmic->GetN()-1],50000,0);
      grfitCosmic->SetLineColor(2);
      grfitCosmic->Draw("lu");
    }
    if (fitMIP) {
      graphMIP->SetMarkerStyle(25);
      graphMIP->Draw("lp");
      TGraph * grfitMIP = fitMIP->MakeGraph(graphMIP->GetX()[0],graphMIP->GetX()[graphMIP->GetN()-1],50000,0);
      grfitMIP->SetLineColor(2);
      grfitMIP->Draw("lu");

    }
 
  }

  return kTRUE;

}


/* Make dummy entry - pseudo code 

 TObjArray * splineArray = new TObjArray(2);

 AliSplineFit * fitMIP = 0;
 AliSplineFit * fitCosmic = 0;
 (*splineArray)[0] = fitMIP;
 (*splineArray)[1] = fitCosmic;
 
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Alexander Kalweit");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-23-02"); //root version
  metaData->SetComment("Calibration of the time dependence of the gain due to pressure and temperature changes");
  AliCDBId id1("TPC/Calib/TimeGain", 0, AliCDBRunRange::Infinity());
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage("local:///d/alice05/akalweit/projects/TimeGainCalibration");
  gStorage->Put(splineArray, id1, metaData);


*/




