
// small macro to generate and update OCDB entries for a given run:
//
// this is a TObjArray which has at 0 the MIP-Spline and at 1 the Fermi-Plateau-Spline ...

/* How to use it...

gSystem->Load("libSTEER");
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");
gSystem->AddIncludePath("-I$ALICE_ROOT/STEER");
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");

.L $ALICE_ROOT/TPC/CalibMacros/CalibTimeVdrift.C+

TFile fcalib("CalibObjects.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTime * timeDrift = ( AliTPCcalibTime *)array->FindObject("calibTime");

CalibTimeVdrift(55079, gain)
  
*/

#include "TObjArray.h"
#include "TGraphErrors.h"
#include "AliExternalTrackParam.h"
#include "TFile.h"
#include "TGraph.h"

#include "AliTPCcalibTime.h"
#include "AliSplineFit.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"


void CalibTimeVdrift(Int_t runNumber, AliTPCcalibTime * vdrift){

  TObjArray * splineArray = new TObjArray(6);

  // 0. Vdrift, best value
  (*splineArray)[0] = vdrift->GetFitDrift(); //So far only cosmics available

  // 1. Vdrift from cosmics
  (*splineArray)[1] = (*splineArray)[0]; //Only cosmics so far
  
  // 2. Vdrift from beam
  (*splineArray)[2] = NULL; //Not yet implemnted

  // 3. Vdrift from ITS-TPC
  (*splineArray)[3] = NULL; //Not yet implemnted

  // 4. Vdrift from goophie
  (*splineArray)[4] = NULL; //Not yet implemnted

  // 5. Vdrift from laser
  (*splineArray)[4] = NULL; //Not yet implemnted
    
  //
  // store OCDB entry
  //

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Dag Toppe Larsen");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-23-02"); //root version
  metaData->SetComment("Calibration of the time dependence of the drift velocity due to pressure and temperature changes");
  AliCDBId id1("TPC/Calib/TimeGain", runNumber, runNumber);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
  gStorage->Put(splineArray, id1, metaData);


}



void MakePlot(){
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibTime * timeDrift = ( AliTPCcalibTime *)array->FindObject("calibTime");
  
  //
  //
  //
  TGraph * grAll   = timeDrift->GetGraphDrift("all");
  TGraph * grD0SCO = timeDrift->GetGraphDrift(" D0SCO ");
  TGraph * grD0ASL = timeDrift->GetGraphDrift(" D0ASL ");
  TGraph * grDEMTY = timeDrift->GetGraphDrift(" DEMPTY ");
  TGraph * grtrdbytof = timeDrift->GetGraphDrift(" trdbytof ");
  TGraph * grtoftrd = timeDrift->GetGraphDrift(" toftrd ");
  //
  grAll->SetLineColor(2);
  grD0SCO->SetLineColor(3);
  grD0ASL->SetLineColor(4);
  grDEMTY->SetLineColor(5);
  grtoftrd->SetLineColor(6);
  //
  grAll->SetMaximum(0.02);
  grAll->SetMinimum(-0.02);
  grAll->Draw("alp");
  grD0SCO->Draw("lp");
  grD0ASL->Draw("lp");
  grDEMTY->Draw("lp");
  grtoftrd->Draw("lp");
}
