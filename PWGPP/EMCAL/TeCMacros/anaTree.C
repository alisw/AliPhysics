#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <AliOADBContainer.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TString.h>
#include "createTree.C"
#endif

Bool_t writeDetailed  = kFALSE;       // switch on detailed information storing
Int_t debugInfo       = 1;            // enable different levels of debuggin information
Bool_t useWeight      = kTRUE;
Int_t nBinsT          = 1000;

void anaTree(
              const char *ifile     ="treefile.root",
              const char *ofile     ="outhist.root",
              Bool_t applyDurLimit  = kFALSE,                   // bolean to switch on duration limit
              Float_t durMin        = 0,                        // minimum length of a run in minutes
              Float_t durMax        = 10000,                    // maximum length of a run in minutes
              Bool_t appBC          = kFALSE,                   // boolean to switch on bad channel
              Int_t referenceRun    = -1                        // define reference run to which all are being calibrated
            )
{
  // Load EMCAL geometry for run 2
  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  const Int_t kSM=g->GetNumberOfSuperModules();
  const Int_t kNcells=g->GetNCells();
  const Int_t gain = 1;

  // Return info about optional settings
  if (applyDurLimit){
    cout << "INFO: only runs with a length of " << durMin << " to "  << durMax << " minutes will be considered in the analysis" << endl;
  }
  if (appBC){
    cout << "INFO: will be using bad channel map" << endl;
  }
  if (referenceRun != -1){
    cout << "INFO: will be storing ratios from run: " << referenceRun << " as reference numbers" << endl;
  }

  // initialize info from tree created by $ALICE_PHYSICS/PWGPP/EMCAL/TeCMacros/createTree.C
  TCalInfo *info = 0;
  TFile *in = TFile::Open(ifile,"read");
  TTree *tt = (TTree*)in->Get("tcal");
  tt->SetBranchAddress("event",&info);
  tt->Branch("event", &info, 32000, 99);
  Int_t Nev=tt->GetEntries();

  Int_t minRunNo = 1e6;
  Int_t maxRunNo = -1;
  for (Int_t i=0;i<Nev;++i) {
    tt->GetEvent(i);
    if (info->fRunNo < minRunNo) minRunNo = info->fRunNo;
    if (info->fRunNo > maxRunNo) maxRunNo = info->fRunNo;
  }

  TProfile *gLedVsT[20];
  TProfile *gLedMonVsT[20];
  TH2F* hLedVsLength[20];
  TH2F* hLedMonVsLength[20];
  TProfile *gRatVsT[20];
  TH1D* hRefRunCellIdVsRat = new TH1D("", "Led/LedMon Ref run; Cell ID", kNcells+1, -0.5, kNcells+1-0.5);
  hRefRunCellIdVsRat->SetName("ReferenceRunRatios");
  TH1D* hRefSMVsT = new TH1D("", "T Ref run; SM", 20, -0.5, 20-0.5);
  hRefSMVsT->SetName("ReferenceRunTemperatures");
  TH1F* hAverageT[20];
  TH1F* hAverageTSorted[20];

  const char* opt = ""; //"S" for spread
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i] = new TProfile("","Led info;T;",nBinsT,15,40,opt);
    gLedVsT[i]->SetName(Form("ledsm%d",i));
    gLedMonVsT[i] = new TProfile("","Led info;T;",nBinsT,15,40,opt);
    gLedMonVsT[i]->SetName(Form("ledmonsm%d",i));
    gRatVsT[i] = new TProfile("","Led/LedMon;T",nBinsT,15,40,opt);
    gRatVsT[i]->SetName(Form("ledovermonsm%d",i));
    hLedVsLength[i]     = new TH2F ("","Led info; t [h];",500,0,20, 20000, 0, 20000);
    hLedVsLength[i]->SetName(Form("ledVsLength%d",i));
    hLedMonVsLength[i]  = new TH2F ("","Led Mon info; t [h];",500,0,20, 20000, 0, 40000);
    hLedMonVsLength[i]->SetName(Form("ledMonVsLength%d",i));
    hAverageT[i]     = new TH1F ("",Form("T SM %d ; run ID; T",i),Nev,0.5,Nev+0.5);
    hAverageT[i]->SetName(Form("TAverageSM%dvsRunId",i));
    hAverageTSorted[i]     = new TH1F ("",Form("T SM %d ; run number; T",i),maxRunNo-minRunNo+1,minRunNo-0.5,maxRunNo+0.5);
    hAverageTSorted[i]->SetName(Form("TAverageSM%dvsRunNumber",i));
  }

  TH1F* hSensorsT[160];
  TH1F* hSensorsTSorted[160];
  for (Int_t sens = 0; sens< 160; sens++){
    hSensorsT[sens]     = new TH1F ("",Form("T sensor %d ; run ID; T",sens),Nev,0.5,Nev+0.5);
    hSensorsT[sens]->SetName(Form("Tsensor%dvsRunId",sens));
    hSensorsTSorted[sens]     = new TH1F ("",Form("T sensor %d ; run number; T",sens),maxRunNo-minRunNo+1,minRunNo-0.5,maxRunNo+0.5);
    hSensorsTSorted[sens]->SetName(Form("Tsensor%dvsRunNumber",sens));
  }

  TProfile *gLedCellVsT[kNcells+1];
  TProfile *gLedCellRMSDiffMeanVsT[kNcells+1];
  TProfile *gLedMonCellVsT[kNcells+1];
  TProfile *gLedMonCellRMSDiffMeanVsT[kNcells+1];
  TProfile *gRatCellVsT[kNcells+1];
  TProfile *gRatECellVsT[kNcells+1];

  cout << "Initializing cell histos" << endl;
  for (Int_t j=0;j<kNcells+1;++j) {
    if (j%500 == 0) cout << "-->next 500: " << j << endl;
    if (writeDetailed){
      gLedCellVsT[j] = new TProfile("",Form("Led info cell ID%i ;T;",j),nBinsT,15,40,opt);
      gLedCellVsT[j]->SetName(Form("ledCell%d",j));
      gLedMonCellVsT[j] = new TProfile("",Form("LedMon info cell ID%i ;T;",j),nBinsT,15,40,opt);
      gLedMonCellVsT[j]->SetName(Form("ledMonCell%d",j));
      gLedCellRMSDiffMeanVsT[j] = new TProfile("",Form("Led rms/mean info cell ID%i ;T;",j),nBinsT,15,40,opt);
      gLedCellRMSDiffMeanVsT[j]->SetName(Form("ledCellRMSDiffMean%d",j));
      gLedMonCellRMSDiffMeanVsT[j] = new TProfile("",Form("LedMon rms/mean info cell ID%i ;T;",j),nBinsT,15,40,opt);
      gLedMonCellRMSDiffMeanVsT[j]->SetName(Form("ledMonRMSDiffMeanCell%d",j));
    }
    gRatCellVsT[j] = new TProfile("",Form("Led/LedMon cell ID%i ;T;",j),nBinsT,15,40,opt);
    gRatCellVsT[j]->SetName(Form("ledovermonCell%d",j));
    gRatECellVsT[j] = new TProfile("",Form("Led/LedMon Error cell ID%i ;T;",j),nBinsT,15,40,opt);
    gRatECellVsT[j]->SetName(Form("ledovermonECell%d",j));
  }
  cout << "-> done initializing histos" << endl;

  Bool_t isRefRun   = kFALSE;
  Bool_t hadRefRun  = kFALSE;

  for (Int_t i=0;i<Nev;++i) {
    tt->GetEvent(i);
    Float_t deltaTime = ((Float_t)info->fLastTime-(Float_t)info->fFirstTime)/60.; // run duration in minutes
    cout << info->fRunNo << "\t" << info->fFirstTime << "\t"<<info->fLastTime << "\t"<< deltaTime<<  endl;

    if (referenceRun != -1 && referenceRun == info->fRunNo){
      isRefRun  = kTRUE;
      hadRefRun = kTRUE;
    } else {
      isRefRun  = kFALSE;
    }

    if (applyDurLimit){
      if (deltaTime < durMin || deltaTime > durMax){
        cout << "INFO: skipped run due to mismatch in run length" << endl;
        continue;
      }
    }

    TClonesArray &sms = info->fSMs;
    for (Int_t sm=0;sm<sms.GetEntries();++sm) {
      TCalSM *smInfot = static_cast<TCalSM*>(sms.At(sm));
      hAverageT[sm]->SetBinContent(i+1,smInfot->fAvgT);
      hAverageTSorted[sm]->SetBinContent(hAverageTSorted[sm]->FindBin(info->fRunNo),smInfot->fAvgT);
      hSensorsT[sm*8]->SetBinContent(i+1,smInfot->fT1);
      hSensorsT[sm*8+1]->SetBinContent(i+1,smInfot->fT2);
      hSensorsT[sm*8+2]->SetBinContent(i+1,smInfot->fT3);
      hSensorsT[sm*8+3]->SetBinContent(i+1,smInfot->fT4);
      hSensorsT[sm*8+4]->SetBinContent(i+1,smInfot->fT5);
      hSensorsT[sm*8+5]->SetBinContent(i+1,smInfot->fT6);
      hSensorsT[sm*8+6]->SetBinContent(i+1,smInfot->fT7);
      hSensorsT[sm*8+7]->SetBinContent(i+1,smInfot->fT8);
      hSensorsTSorted[sm*8]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT1);
      hSensorsTSorted[sm*8+1]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT2);
      hSensorsTSorted[sm*8+2]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT3);
      hSensorsTSorted[sm*8+3]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT4);
      hSensorsTSorted[sm*8+4]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT5);
      hSensorsTSorted[sm*8+5]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT6);
      hSensorsTSorted[sm*8+6]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT7);
      hSensorsTSorted[sm*8+7]->SetBinContent(hSensorsTSorted[sm*8]->FindBin(info->fRunNo),smInfot->fT8);
    }


    TClonesArray &cells = info->fCells;
    for (Int_t j=0;j<cells.GetEntries();++j) {
      TCalCell *cell = static_cast<TCalCell*>(cells.At(j));
      Int_t cellID  = cell->fId;
      Int_t sm      = cell->fSM;
      Int_t badcell = cell->fBad;
      Double_t ledM = cell->fLedM;
      Double_t ledR = cell->fLedR;
      Double_t monM = cell->fMonM;
      Double_t monR = cell->fMonR;
      Double_t locT = cell->fLocT;
      Double_t smT  = cell->fSMT;

      if (appBC && badcell > 0){
        if (debugInfo > 1) cout << "found bad cell for " << cellID <<  "\t " << info->fRunNo << endl;
        continue;
      }

      Double_t T = smT; // use local or SM T, change here
      if ((T<5)||(T>45))
        continue;

      if (isRefRun){
        if (hRefSMVsT->GetBinContent(sm+1) < 5) hRefSMVsT->SetBinContent(sm+1,smT);
      }

      if ((ledM<=0)||(ledR<=0))
        continue;
      if ((monM<=0)||(monR<=0))
        continue;
      Double_t w        = 1./(ledR*ledR);
      if (!useWeight) w = 1;
      gLedVsT[sm]->Fill(T,ledM,w);
      if (writeDetailed){
        gLedCellVsT[cellID]->Fill(smT,ledM,w);
        gLedCellRMSDiffMeanVsT[cellID]->Fill(smT,ledR/ledM,1);
      }
      Double_t w3       = 1./(monR*monR);
      if (!useWeight) w3 = 1;
      gLedMonVsT[sm]->Fill(T,monM,w3);
      if (writeDetailed){
        gLedMonCellVsT[cellID]->Fill(smT,monM,w3);
        gLedMonCellRMSDiffMeanVsT[cellID]->Fill(smT,monR/monM,1);
      }
      Double_t ratErr = ledM/monM * TMath::Sqrt((ledR*ledR)/(ledM*ledM)+(monR*monR)/(monM*monM));
      Double_t w2=1./(ratErr*ratErr);
      if (!useWeight) w2 = 1;
      gRatVsT[sm]->Fill(T,ledM/monM,w2);
      gRatCellVsT[cellID]->Fill(smT,ledM/monM,w2);
      gRatECellVsT[cellID]->Fill(smT,ratErr,w2);
      hLedVsLength[sm]->Fill(deltaTime/60,ledM,w);
      hLedMonVsLength[sm]->Fill(deltaTime/60,monM,w);
      if (isRefRun){
        hRefRunCellIdVsRat->SetBinContent(cellID+1,ledM/monM);
        hRefRunCellIdVsRat->SetBinError(cellID+1,ratErr);
      }
    }
  }

  TFile *out = TFile::Open(ofile,"recreate");
  if (!out)
    out = TFile::Open("dummyfile.root","update");
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i]->Write();
    gLedMonVsT[i]->Write();
    gRatVsT[i]->Write();
    hLedVsLength[i]->Write();
    hLedMonVsLength[i]->Write();
    hAverageT[i]->Write();
    hAverageTSorted[i]->Write();
  }
  if (hadRefRun){
    hRefRunCellIdVsRat->Write();
    hRefSMVsT->Write();
  }
  for (Int_t sens = 0; sens < 160; sens++){
    hSensorsT[sens]->Write();
    hSensorsTSorted[sens]->Write();
  }

  Int_t smcurr          = -1;
  Int_t cellIDFirstInSM = 0;
  for (Int_t j=0;j<kNcells+1;++j) {
    Int_t iTower = -1, iIphi = -1,  iIeta = -1, sm=-1;
    g->GetCellIndex(j,sm,iTower,iIphi,iIeta);
    if (smcurr != sm){
      smcurr = sm;
      out->mkdir(Form("cellsInSM%d",sm));
      out->cd(Form("cellsInSM%d",sm));
      if (debugInfo) cout << "SM: \t" << sm-1 << "\t"<< cellIDFirstInSM << "\t" << j << "\t" << j -cellIDFirstInSM<< endl;
      cellIDFirstInSM = j;
    }

    if (writeDetailed) {
      if (gLedCellVsT[j]->GetEntries() > 0)
        gLedCellVsT[j]->Write();
      if (gLedMonCellVsT[j]->GetEntries() > 0)
        gLedMonCellVsT[j]->Write();
      if (gLedCellRMSDiffMeanVsT[j]->GetEntries() > 0)
        gLedCellRMSDiffMeanVsT[j]->Write();
      if (gLedMonCellRMSDiffMeanVsT[j]->GetEntries() > 0)
        gLedMonCellRMSDiffMeanVsT[j]->Write();
    }
    if (gRatCellVsT[j]->GetEntries() > 0)
      gRatCellVsT[j]->Write();
    if (gRatECellVsT[j]->GetEntries() > 0)
      gRatECellVsT[j]->Write();
  }
}
