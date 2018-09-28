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
#include <TRandom3.h>
#include "createTree.C"
#endif

Bool_t writeDetailed  = kFALSE;       // switch on detailed information storing
Int_t debugInfo       = 1;            // enable different levels of debuggin information
Bool_t useWeight      = kTRUE;
Bool_t runWithRand    = kTRUE;

void anaTree(
              const char *ifile     ="treefile.root",
              const char *ofile     ="outhist.root",
              Bool_t applyDurLimit  = kFALSE,                   // bolean to switch on duration limit
              Float_t durMin        = 0,                        // minimum length of a run in minutes
              Float_t durMax        = 10000,                    // maximum length of a run in minutes
              Bool_t appBC          = kFALSE,                   // boolean to switch on bad channel
              Int_t referenceRun    = 286313,                   // define reference run to which all are being calibrated
              Int_t nBinsT          = 1000
            )
{
  // Load EMCAL geometry for reference run
  AliEMCALGeometry *g= AliEMCALGeometry::GetInstanceFromRunNumber(referenceRun);
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
  cout << minRunNo << "\t" << maxRunNo << endl;

  TProfile *gLedVsT[kSM];
  TProfile *gLedMonVsT[kSM];
  TH2F* hLedVsLength[kSM];
  TH2F* hLedMonVsLength[kSM];
  TProfile *gRatVsT[kSM];
  TProfile* gCellIdVsRat = new TProfile("", "Led/LedMon run; Cell ID", kNcells+1, -0.5, kNcells+1-0.5);
  gCellIdVsRat->SetName("LedDiffLedMonVsCellID");
  TH2D* gCellIdVsLed = new TH2D("", "Led run; Cell ID", kNcells+1, -0.5, kNcells+1-0.5, 1000, 0, 1000);
  gCellIdVsLed->SetName("LedVsCellID");
  TH2D* gCellIdVsLedMon = new TH2D("", "LedMon run; Cell ID", kNcells+1, -0.5, kNcells+1-0.5, 1000, 0, 1000);
  gCellIdVsLedMon->SetName("LedMonVsCellID");
  TH1F* hAverageT[kSM];
  TH1F* hAverageTSorted[kSM];
  TH2F* hAverageTPerSM      = new TH2F("","T per SM; SM; T", kSM, -0.5, kSM-0.5, nBinsT, 15, 40);
  hAverageTPerSM->SetName("MeanSMTemperature");
  hAverageTPerSM->Sumw2();
  const char* opt = "S"; //"S" for spread
  for (Int_t i=0;i<kSM;++i) {
    gLedVsT[i] = new TProfile("","Led info;T;",nBinsT,17, 27);
    gLedVsT[i]->SetName(Form("ledsm%d",i));
    gLedMonVsT[i] = new TProfile("","Led info;T;",nBinsT,17, 27);
    gLedMonVsT[i]->SetName(Form("ledmonsm%d",i));
    gRatVsT[i] = new TProfile("","Led/LedMon;T",nBinsT,17, 27);
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
    gRatCellVsT[j] = new TProfile("",Form("Led/LedMon cell ID%i ;T;",j),nBinsT,17, 27,opt);
    gRatCellVsT[j]->SetName(Form("ledovermonCell%d",j));
    gRatECellVsT[j] = new TProfile("",Form("Led/LedMon Error cell ID%i ;T;",j),nBinsT,17, 27,opt);
    gRatECellVsT[j]->SetName(Form("ledovermonECell%d",j));
  }
  cout << "-> done initializing histos" << endl;

  Bool_t isRefRun   = kFALSE;
  Bool_t hadRefRun  = kFALSE;

  cout << "-> there are " << Nev << " contained in this tree, starting to analyse them" << endl;

  for (Int_t i=0;i<Nev;++i) {
    tt->GetEvent(i);
    UInt_t deltaTimeS = ((UInt_t)info->fLastTime-(UInt_t)info->fFirstTime)/10;         // run duration in seconds
    Float_t deltaTime = ((Float_t)info->fLastTime-(Float_t)info->fFirstTime)/60.;   // run duration in minutes
    if (i%50 == 0)
      cout << "starting with run " << i << "/" << Nev << endl;
    cout << info->fRunNo << "\t" << info->fFirstTime << "\t"<<info->fLastTime << "\t"<< deltaTime<<  endl;

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
      hAverageTPerSM->Fill(sm,smInfot->fAvgT,deltaTime/60);
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

    TRandom3 ledRandom(1289790);
    TRandom3 ledMonRandom(909088);

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

      if (runWithRand){
        for ( UInt_t s = 0; s < deltaTimeS; s++){
          Double_t ledcurrent     = 0;
          Double_t ledMoncurrent  = 0;
          if (!((ledM<=0)||(ledR<=0)) ){
            ledcurrent              = ledRandom.Gaus(ledM,ledR);
            if (writeDetailed){
              gLedCellVsT[cellID]->Fill(smT,ledcurrent);
              if (s==0) gLedCellRMSDiffMeanVsT[cellID]->Fill(smT,ledR/ledM,1);

            }
            gLedVsT[sm]->Fill(T,ledcurrent);
            gCellIdVsLed->Fill(cellID,ledcurrent);
          }
          if (!((monM<=0)||(monR<=0)) ){
            ledMoncurrent         = ledMonRandom.Gaus(monM,monR);
            if (writeDetailed){
              gLedMonCellVsT[cellID]->Fill(smT,ledMoncurrent);
              if (s==0) gLedMonCellRMSDiffMeanVsT[cellID]->Fill(smT,monR/monM,1);
            }
            gLedMonVsT[sm]->Fill(T,ledMoncurrent);
            gCellIdVsLedMon->Fill(cellID,ledMoncurrent);
          }
          if ( (ledcurrent<=0) || (ledMoncurrent<=0) )
            continue;
          gCellIdVsRat->Fill(cellID,ledcurrent/ledMoncurrent);
          gRatVsT[sm]->Fill(T,ledcurrent/ledMoncurrent);
          gRatCellVsT[cellID]->Fill(smT,ledcurrent/ledMoncurrent);
          gRatECellVsT[cellID]->Fill(smT,1);
          hLedVsLength[sm]->Fill(deltaTime/60,ledcurrent);
          hLedMonVsLength[sm]->Fill(deltaTime/60,ledMoncurrent);
        }
      } else {
        Double_t w        = 1;
        if (ledR != 0 && useWeight)
          w               = 1./(ledR*ledR);
        Double_t w3       = 1.;
        if (monR != 0 && useWeight)
          w3              = 1./(monR*monR);;

        if (!((ledM<=0)||(ledR<=0)) ){
          if (writeDetailed){
            gLedCellVsT[cellID]->Fill(smT,ledM,w);
            gLedCellRMSDiffMeanVsT[cellID]->Fill(smT,ledR/ledM,1);
          }
          gLedVsT[sm]->Fill(T,ledM,w);
        }

        if (!((monM<=0)||(monR<=0)) ){
          if (writeDetailed){
            gLedMonCellVsT[cellID]->Fill(smT,monM,w3);
            gLedMonCellRMSDiffMeanVsT[cellID]->Fill(smT,monR/monM,1);
          }
          gLedMonVsT[sm]->Fill(T,monM,w3);
        }
        if ( (ledM<=0)||(ledR<=0) || (monM<=0)||(monR<=0) )
          continue;

        Double_t ratErr = ledM/monM * TMath::Sqrt((ledR*ledR)/(ledM*ledM)+(monR*monR)/(monM*monM));
        Double_t w2=1./(ratErr*ratErr);
        if (!useWeight) w2 = 1;
        gRatVsT[sm]->Fill(T,ledM/monM,w2);
        gRatCellVsT[cellID]->Fill(smT,ledM/monM,w2);
        gRatECellVsT[cellID]->Fill(smT,ratErr,w2);
        hLedVsLength[sm]->Fill(deltaTime/60,ledM,w);
        hLedMonVsLength[sm]->Fill(deltaTime/60,monM,w);
      }
    }
  }

  TFile *out = TFile::Open(ofile,"recreate");
  if (!out)
    out = TFile::Open("dummyfile.root","update");
  for (Int_t i=0;i<kSM;++i) {
    gLedVsT[i]->Write();
    gLedMonVsT[i]->Write();
    gRatVsT[i]->Write();
    hLedVsLength[i]->Write();
    hLedMonVsLength[i]->Write();
    hAverageT[i]->Write();
    hAverageTSorted[i]->Write();
  }
  hAverageTPerSM->Write();
  gCellIdVsRat->Write();
  gCellIdVsLed->Write();
  gCellIdVsLedMon->Write();
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
