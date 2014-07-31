/*
g++ makeTRDdEdxOADB.C  -g -O3 -Wall -Werror -I$ROOTSYS/include -I$ALICE_BUILD/include -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TRD/Cal  -L$ROOTSYS/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit -lEG -lGeom -lVMC -lProof -lProofPlayer -lXMLParser -lXMLIO -lSpectrum -lTreePlayer -lMLP -lGui -L$ALICE_BUILD/lib/tgt_linuxx8664gcc -lCDB -lSTEER -lRAWDatarec -lESD -lSTEERBase -lANALYSIS -lRAWDatabase  -lANALYSISalice -lAOD -lTPCrec -lTPCbase -lTRDbase -lTRDrec -lSTAT   -o makeTRDdEdxOADB

#generate OADB
makeTRDdEdxOADB 1

#read and print OADB content
makeTRDdEdxOADB 0
 */

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "TASImage.h"
#include "TAxis.h"
#include "TColor.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphPolar.h"
#include "TGrid.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLinearFitter.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TUUID.h"
#include "TVector3.h"
#include "TVectorD.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAODInputHandler.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCentralitySelectionTask.h"
#include "AliCTPRawStream.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliGeomManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMathBase.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPID.h"
#include "AliRawReaderRoot.h"
#include "AliReconstruction.h"
//#include "AliTPCcalibDB.h"
#include "AliTPCclusterMI.h"
//#include "AliTPCParam.h"
//#include "AliTPCRecoParam.h"
//#include "AliTPCseed.h"
//#include "AliTPCTransform.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"
#include "AliTrackPointArray.h"
#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDrawStream.h"
#include "AliTRDseedV1.h"
#include "AliTriggerAnalysis.h"
#include "AliSysInfo.h"


#include "AliTRDdEdxParams.h"
#include "AliOADBContainer.h"
#include "AliLog.h"

AliTRDdEdxParams * getDefault()
{
  //
  //version 31/07/2014 by Xianguo Lu based on pPb small sample only for nch=6 ncls/nch>=18, no binning in eta, no binning on multiplicity
  //

  AliTRDdEdxParams *defobj=new AliTRDdEdxParams("default","default");

  const Float_t meanpars[]={1.521e+00 , 7.294e-01 , 8.066e+00 , 2.146e-02 , 3.137e+01 , 7.160e-15 , 2.838e+00 , 1.100e+01};
  const Float_t respars[]={-4.122e-01 , 1.019e+00 , -1.319e-01};

  for(Int_t ipar=0; ipar< MAXNPAR; ipar++){
    if(ipar == AliPID::kProton){
      const Float_t tmpproton[]={2.026e+00 , -1.462e-04 , 1.202e+00 , 1.162e-01 , 2.092e+00 , -3.018e-02 , 3.665e+00 , 3.487e-07}; 
      defobj->SetMeanParameter(ipar, sizeof(tmpproton)/sizeof(Float_t), tmpproton);
    }
    else if(ipar == AliPID::kPion || ipar ==AliPID::kElectron){       
      const Float_t tmppe[]={1.508e+00 , 7.356e-01 , 8.002e+00 , 1.932e-02 , 3.058e+01 , 5.114e-18 , 3.561e+00 , 1.408e+01};  
      defobj->SetMeanParameter(ipar, sizeof(tmppe)/sizeof(Float_t), tmppe);
    }
    else{
      defobj->SetMeanParameter(ipar, sizeof(meanpars)/sizeof(Float_t), meanpars);
    }

    defobj->SetSigmaParameter(ipar,  sizeof(respars)/sizeof(Float_t), respars);
  }

  defobj->Print();

  return defobj;
}

void makeTRDdEdxOADB(const Bool_t kMC=0)
{
  //
  //make OADB object
  //currently only default
  //non-default has to be tested after finer inspection on dependence on eta, multiplicity, etc. 
  //

  AliTRDdEdxParams * defobj = getDefault();

  TString containerName = "TRDdEdxParamsContainer";
  AliOADBContainer* cont = new AliOADBContainer(containerName.Data());

  //current convention as other detetor is only to provide for data
  //no MC is seen in current (31/07/2014)  AliROOT except for HMPID

  const TString filePathNamePackage=Form("$ALICE_ROOT/OADB/COMMON/PID/%s/TRDdEdxParams.root", kMC?"MC":"data");
  Int_t statusCont = cont->InitFromFile(filePathNamePackage.Data(), cont->GetName());
  if (statusCont) {
    printf("No OADBContainer for the current settings found - creating a new one...\n");
  }

  //non-default params
  /*
  const Double_t runLow = 1;
  const Double_t runUp =  10;
  const Int_t index = cont->GetIndexForRun((runUp + runLow) / 2.0);
  if (index < 0) {
    printf("Creating new object for run range %d - %d...\n", runLow, runUp);
    cont->AppendObject(hh, runLow, runUp);
  }
  else {
    printf("Updating existing object for run range %d - %d...\n", runLow, runUp);
    cont->UpdateObject(index, hh, runLow, runUp);
  }
  */

  //default params
  cont->CleanDefaultList(); // Remove old default objects at first
  cont->AddDefaultObject(defobj);


  TFile* f = new TFile(filePathNamePackage.Data(), "update");
  f->Delete(cont->GetName());
  cont->Write(0, TObject::kOverwrite);
  f->Purge();
  f->Close();

  printf("makeTRDdEdxOADB done\n");
}

void readTRDdEdxOADB(const Bool_t kMC=0)
{
  //
  //only read in default
  //non-default has to be uncommented if needed
  //

  const TString containerName = "TRDdEdxParamsContainer";
  AliOADBContainer cont(containerName.Data()); 
  
  const TString filePathNamePackage=Form("$ALICE_ROOT/OADB/COMMON/PID/%s/TRDdEdxParams.root", kMC?"MC":"data");

  const Int_t statusCont = cont.InitFromFile(filePathNamePackage.Data(), cont.GetName());
  if (statusCont){
    printf("Failed initializing settings from OADB\n"); exit(1);
  }

  /*
  Int_t fRun = 5;
  AliTRDdEdxParams *jj=(AliTRDdEdxParams*)(cont.GetObject(fRun, defaultObj.Data()));
  jj->Print();
  */

  TString defaultObj = "default";
  AliTRDdEdxParams *jj=(AliTRDdEdxParams*)(cont.GetObject(1000, defaultObj.Data()));
  jj->Print();

}

int main(int argc, char * argv[])
{
  for(Int_t ii=0; ii<argc; ii++){
    printf("%d: %s\n", ii, argv[ii]);
  }

  if(argc!=2){
    printf("argc!=2\n");
    return 1; 
  }

  if(atoi(argv[1])){
    makeTRDdEdxOADB();
  }
  else{
    readTRDdEdxOADB();
  }

  return 0;
}
