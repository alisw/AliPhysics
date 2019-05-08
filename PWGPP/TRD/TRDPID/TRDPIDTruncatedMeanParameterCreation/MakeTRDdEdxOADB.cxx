//
// macro to write TRD dE/dx parameters to OADB container
//

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
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
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMathBase.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliRawReaderRoot.h"
#include "AliReconstruction.h"
#include "AliTPCclusterMI.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"
#include "AliTrackPointArray.h"
#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDrawStream.h"
#include "AliTRDseedV1.h"
#include "AliSysInfo.h"

#include "AliCentralitySelectionTask.h"
#include "AliPhysicsSelectionTask.h"
#include "AliTriggerAnalysis.h"

#include "AliTRDdEdxParams.h"
//#include "AliTRDdEdxParams.cxx"
#include "AliOADBContainer.h"
#include "AliLog.h"

AliTRDdEdxParams * getParameter( Float_t*  tmpresParsEtaCorr,  Float_t* tmpresParsEtaUncorr,  Float_t* tmpParsEtaCorr,  Float_t* tmpParsEtaUncorr) {
    //
    //version 31/07/2014 by Xianguo Lu based on pPb small sample only for nch=6 ncls/nch>=18, no binning in eta, no binning on multiplicity
    
    // modified 05/05/2015 by Lucas Altenkämper based on pPb sample for nch=6, 5, 4 and ncls/nch>=17, no binning in eta, no binning on multiplicity   
    
    // comment 20/04/2016 by FHerrmann - we decided to use just nch >=4, but to have the possibility to return to the old scheme, we wrote the same parameters in all nch below
    
    AliTRDdEdxParams *cont = new AliTRDdEdxParams("parameter", "parameter");

    const Float_t meanParsDummy[8]  = {0};
    const Float_t resParsDummy[2]   = {0};
    
    Float_t resParsEtaCorr[2];
    Float_t resParsEtaUncorr[2];
    Float_t ParsEtaCorr[8], ParsEtaUncorr[8];

    // temporary solution
    for (int i=0; i<2; i++) { // resoultuion parameter w/wo corrections (all corrections or no correction)
      resParsEtaCorr[i]=tmpresParsEtaCorr[i];
      resParsEtaUncorr[i]=tmpresParsEtaUncorr[i];
    }
    
    for (int i=0; i<8; i++) { // signal parameter w/wo corrections
      ParsEtaCorr[i]=tmpParsEtaCorr[i];
      ParsEtaUncorr[i]=tmpParsEtaUncorr[i];
    }

    //const Float_t resParsEtaCorr[2]     = {3.014e-02, 1.390e+00	};
    //const Float_t resParsEtaUncorr[2]   = {4.715e-02, 1.368e+00	};
    
    for (Int_t icls = 0; icls < 2; icls++) { // charge?? or cluster
      for (Int_t ich = 0; ich < 4; ich++) { //6,5,4,x chamber *x4
	for (Int_t ipar = 0; ipar < 10; ipar++) {// proton, electron/pion, other x10
    
                TVectorF tmp;
                // !!! has to be consistent with
                // Int_t GetIter(const Int_t itype, const Int_t nch, const Int_t ncls) const
                // in AliTRDdEdxParams

                Int_t nch   = 1;
                Int_t ncls  = 1;
                
                if (icls == 0 && ich == 0) {
                    nch     = 6;
                    ncls    = 180;   // just to satisfy ncls/nch>=17 - actual value not of interest (Florian)
                    
		    //   const Float_t tmpParsEtaCorr[8]     = ParsEtaCorr; //{4.369e-01, 1.959e+00, 7.726e+00, 3.070e-01, 2.942e+00, 2.959e-02, 2.269e+00, 8.138e-01};
		    //  const Float_t tmpParsEtaUncorr[8]    = ParsEtaUnCorr; //{6.129e-01, 1.816e+00, 7.767e+00, 2.242e-01, 3.836e+00, 1.001e-02, 2.326e+00, 1.060e+00};

                    if (ipar == AliPID::kProton) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else if (ipar == AliPID::kPion || ipar == AliPID::kElectron) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    }
                } else if (icls == 0 && ich == 1) {
                    nch     = 5;
                    ncls    = 180;

		    //  const Float_t tmpParsEtaCorr[]      = {4.369e-01, 1.959e+00, 7.726e+00, 3.070e-01, 2.942e+00, 2.959e-02, 2.269e+00, 8.138e-01};
		    // const Float_t tmpParsEtaUncorr[]    = {6.129e-01, 1.816e+00, 7.767e+00, 2.242e-01, 3.836e+00, 1.001e-02, 2.326e+00, 1.060e+00};
                    
                    if (ipar == AliPID::kProton) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else if (ipar == AliPID::kPion || ipar == AliPID::kElectron) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    }
                } else if (icls == 0 && ich == 2) {
                    nch     = 4;
                    ncls    = 180;
                    
		    //cnst Float_t tmpParsEtaCorr[]      = {4.369e-01, 1.959e+00, 7.726e+00, 3.070e-01, 2.942e+00, 2.959e-02, 2.269e+00, 8.138e-01};
		    //cnst Float_t tmpParsEtaUncorr[]    = {6.129e-01, 1.816e+00, 7.767e+00, 2.242e-01, 3.836e+00, 1.001e-02, 2.326e+00, 1.060e+00};

                    if (ipar == AliPID::kProton) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else if (ipar == AliPID::kPion || ipar == AliPID::kElectron) {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    } else {
                        // eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaCorr)/sizeof(Float_t),    ParsEtaCorr,     kTRUE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaCorr)/sizeof(Float_t),    resParsEtaCorr,     kTRUE);
                        
                        // NOT eta-corrected
                        cont->SetMeanParameter(   ipar, nch, ncls, sizeof(ParsEtaUncorr)/sizeof(Float_t),  ParsEtaUncorr,   kFALSE);
                        cont->SetSigmaParameter(  ipar, nch, ncls, sizeof(resParsEtaUncorr)/sizeof(Float_t),  resParsEtaUncorr,   kFALSE);
                    }
                } else {
                    nch     = 1;
                    ncls    = 1;
                    
                    //bad nch and ncls, 0-size vector
                    cont->SetMeanParameter( ipar,   nch, ncls, 0, meanParsDummy,  kTRUE);
                    cont->SetSigmaParameter(ipar,   nch, ncls, 0, resParsDummy,   kTRUE);
                    
                    cont->SetMeanParameter( ipar,   nch, ncls, 0, meanParsDummy,  kFALSE);
                    cont->SetSigmaParameter(ipar,   nch, ncls, 0, resParsDummy,   kFALSE);
                }
            }
        }
    }

    cont->Print();

    return cont;
}


void makeTRDdEdxOADB(const Bool_t kMC=0) {
  //
  // make OADB object
  //
  
  
  // Read in ResolutionParameter
  TString ParaFileNameCorr = ("output/LHC18d_meanParam_pid3_LHC18d_from4to6Tracklets_etaTPCtglCorr_NclsCorr.txt");
  TString ParaFileNameUnCorr = ("output/LHC18d_meanParam_pid3_LHC18d_from4to6Tracklets.txt"); // "output/LHC15o-compareFullAndLIR/FULL/LHC15o_Pass1FullData_resParam_Iter0_pid3_17ncls_Test__NEwPara__ClusterMeanOverChamber.txt");
  //    TString ParaFileNameCorr = ("output/LHC15o-compareFullAndLIR/FULL/LHC15o_Pass1FullData_meanParam_pid3_6nch_17ncls_Test__NEwPara__ClusterMeanOverChamber_CentCorr_etaTPCtglCorr_NclsCorr.txt"); //
  TString ResFileNameUnCorr = ("output/LHC18d_resParam_Iter0_pid3_LHC18d_from4to6Tracklets.txt");
  //    TString ParaFileNameUnCorr = ("output/LHC15o-compareFullAndLIR/FULL/LHC15o_Pass1FullData_meanParam_pid3_6nch_17ncls_Test__NEwPara__ClusterMeanOverChamber.txt");
  TString ResFileNameCorr = ("output/LHC18d_resParam_Iter0_pid3_LHC18d_from4to6Tracklets_etaTPCtglCorr_NclsCorr.txt");


  ifstream ResFileCorr(ResFileNameCorr);
  ifstream ResFileUnCorr(ResFileNameUnCorr);
  if (!ResFileCorr) cout << "Could not open file " << ResFileNameCorr << endl;
  if (!ResFileUnCorr) cout << "Could not open file " << ResFileNameUnCorr << endl;

  string line;
  TString dummy;

  Float_t ResolutionCorr[2]; 
  Float_t ResolutionUnCorr[2];
   
  while (getline(ResFileCorr, line)) {
    stringstream lineInfo(stringstream::in | stringstream::out);
    size_t found=string::npos;
    for (int i=0; i<2; i++) {
      if (found!=line.find(Form("par[%i]", i))){
	lineInfo << line;
	lineInfo >> dummy >> dummy >> ResolutionCorr[i];
      }
    }
  }
  while (getline(ResFileUnCorr, line)) {
    stringstream lineInfo(stringstream::in | stringstream::out);
    size_t found=string::npos;
    for (int i=0; i<2; i++) {
      if (found!=line.find(Form("par[%i]", i))){
	lineInfo << line;
	lineInfo >> dummy >> dummy >> ResolutionUnCorr[i];
      }
    }
  }
  for (int i=0; i<2; i++) {
    cout << "Resolution Parameter Corr " << i << "\t" << ResolutionCorr[i] << "\t Resolution Parameter UnCorr" << i << "\t" << ResolutionUnCorr[i] <<  endl;
  }

  // Read in Parametrization
  Float_t ParaCorr[8];
  Float_t ParaUnCorr[8];
  ifstream ParaFileCorr(ParaFileNameCorr);
  ifstream ParaFileUnCorr(ParaFileNameUnCorr);
  if (!ParaFileCorr) cout << "Could not open file " << ParaFileNameCorr << endl;
  if (!ParaFileUnCorr) cout << "Could not open file " << ParaFileNameUnCorr << endl;

  cout << "Corrected Parameter:" << endl;
  while (getline(ParaFileCorr, line)) {
    stringstream lineInfo(stringstream::in | stringstream::out);
    size_t found=string::npos;
    for (int i=0;i<8; i++) {
      if (found!=line.find(Form("par[%i]", i))) {
	lineInfo << line;
	lineInfo >> dummy >> dummy >> ParaCorr[i];
      }
    }
  }

  cout << "Uncorrected Parameter:" << endl;
  while (getline(ParaFileUnCorr, line)) {
    stringstream lineInfo(stringstream::in | stringstream::out);
    size_t found=string::npos;
    for (int i=0;i<8; i++) {
      if (found!=line.find(Form("par[%i]", i))) {
	lineInfo << line;
	lineInfo >> dummy >> dummy >> ParaUnCorr[i];
      }
    }
  }

  for (int i=0; i<8; i++) {
    cout << "Parameter Corr " << i << "\t" << ParaCorr[i] << "\t" << "Parameter UnCorr " << i << "\t" << ParaUnCorr[i] << endl; 
  }



  AliTRDdEdxParams * params = getParameter(ResolutionCorr, ResolutionUnCorr, ParaCorr, ParaUnCorr);
    
  TString containerName = "TRDdEdxParamsContainer";
  AliOADBContainer cont(containerName.Data());

  // current convention as other detetor is only to provide for data
  // no MC is seen in current (31/07/2014)  AliROOT except for HMPID
  
  // const TString filePathNamePackage=Form("$ALICE_ROOT/OADB/COMMON/PID/%s/TRDdEdxParams.root", kMC?"MC":"data");
  const TString filePathNamePackage=Form("/home/ikp/alice/AliPhysics/OADB/COMMON/PID/%s/TRDdEdxParams.root", kMC?"MC":"data");
  //const TString filePathNamePackage = Form("TRDdEdxParams.root");
  Int_t statusCont = cont.InitFromFile(filePathNamePackage.Data(), cont.GetName());
  if (statusCont) {
    printf("No OADBContainer for the current settings found - creating a new one...\n");
  }

  // run range

  Int_t fRunRange[2]={286350, 285978}; // LHC18d
  // Int_t fRunRange[2]= {244824,246994}; // LHC15o
  //  Int_t fRunRange[2] ={195344, 197388};// LHC13 bc


  // cont.CleanDefaultList();               // Remove old default objects at first
  // cont.AddDefaultObject(params);
  cout << cont.GetOADBPath() << endl;
  Int_t CheckForExisitingData=-1;
  CheckForExisitingData=cont.GetIndexForRun(fRunRange[0]); //fRunRange[0]);
  cout << "ExistingData" << endl;
  cout << CheckForExisitingData << endl;
  if (CheckForExisitingData!=-1) {
    cout << "UPDATE" << endl;
    cont.UpdateObject(cont.GetIndexForRun(fRunRange[0]), params, fRunRange[0], fRunRange[1]);
  }
  else {
    cout << "APPEND" << endl;
    cont.AppendObject(params, fRunRange[0], fRunRange[1]);
  }
  
  cout << "Store" << endl;
  cont.List();
  
  
  TFile* f = new TFile(filePathNamePackage.Data(), "RECREATE");
  if(f){
    f->Delete("*");
    f->Close();
  }
  cont.WriteToFile(filePathNamePackage.Data());
  f->Print();
  // f->Delete(cont.GetName());
  //cont.Write(0,TObject::kOverwrite); // TObject::kOverwrite);
  f->Purge();
  f ->Close();

  printf("makeTRDdEdxOADB done\n");
}


void readTRDdEdxOADB(const Bool_t kMC=0) {
  //
  //only read in default
  //non-default has to be uncommented if needed
  //

  const TString containerName = "TRDdEdxParamsContainer";
  AliOADBContainer cont(containerName.Data());
  
  const TString filePathNamePackage=Form("/home/ikp/alice/AliPhysics/OADB/COMMON/PID/%s/TRDdEdxParams.root", kMC?"MC":"data");
  //const TString filePathNamePackage=Form("TRDdEdxParams.root");
  
  const Int_t statusCont = cont.InitFromFile(filePathNamePackage.Data(), cont.GetName());
  if (statusCont) {
    printf("Failed initializing settings from OADB\n");
    exit(1);
  }
  
  Int_t run = 244340;
  cout << "Reading OADB entries for run " << run << endl;
  
  TString defaultObj = "default";
  AliTRDdEdxParams *jj = (AliTRDdEdxParams*)(cont.GetObject(run, defaultObj.Data()));             // cont.GetObject() returns TObject
  jj->AliTRDdEdxParams::Print();
}


int main(int argc, char * argv[]) {
  for (Int_t ii = 0; ii < argc; ii++) {
    printf("%d: %s\n", ii, argv[ii]);
  }
  
  if (argc != 2) {
    printf("argc!=2\n");
    return 1;
  }

  if (atoi(argv[1])) {
    cout << "write" << endl;
    makeTRDdEdxOADB();
  
  }
  else {
    cout << "read" << endl;
    readTRDdEdxOADB();
  }

  return 0;
}
