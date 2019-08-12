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

class TCalRun : public TObject {
public:
  TCalRun() : fLedM(-1), fLedR(-1), fMonM(-1), fMonR(-1), fSMT(-1), fLength(0), fRunNo(-1), fBad(0) {;}
  virtual ~TCalRun() {;}
  Double32_t fLedM;       //[0,0,16] led mean
  Double32_t fLedR;       //[0,0,16] led rms
  Double32_t fMonM;       //[0,0,16] mon mean
  Double32_t fMonR;       //[0,0,16] mon rms
  Double32_t fSMT;        //[0,0,16] sm T
  UInt_t     fLength;     // run length in seconds
  Int_t      fRunNo;      // run number
  Short_t    fBad;        // bad cell
  ClassDef(TCalRun, 1); // CalRun class
};

class TCalCellInfo : public TObject {
public:
  TCalCellInfo() : fCellId(-1),  fRuns("TCalRun") {}
  virtual ~TCalCellInfo() {;}
  Int_t         fCellId;        // cell IDs
  TClonesArray  fRuns;
  ClassDef(TCalCellInfo, 1); // CalCellInfo class
};

Bool_t writeDetailed  = kFALSE;       // switch on detailed information storing
Int_t debugInfo       = 1;            // enable different levels of debuggin information
Int_t testCells       = 100;

void anaTreeV3(
                const char *ifile     ="treefile.root",
                const char *ofile     ="outhist.root",
                Bool_t applyDurLimit  = kFALSE,                   // bolean to switch on duration limit
                Float_t durMin        = 0,                        // minimum length of a run in minutes
                Float_t durMax        = 10000,                    // maximum length of a run in minutes
                Bool_t appBC          = kFALSE,                   // boolean to switch on bad channel
                Int_t referenceRun    = 286313,                   // define reference run to which all are being calibrated
                Int_t nBinsT          = 1000,
                TString iterFile      = "",
                Int_t singleSM        = -1,
                Int_t nBins2D         = 500
              )
{
  // Load EMCAL geometry for reference run
  AliEMCALGeometry *g   = AliEMCALGeometry::GetInstanceFromRunNumber(referenceRun);
  const Int_t kSM       = g->GetNumberOfSuperModules();
  const Int_t kNcells   = g->GetNCells();

  // Return info about optional settings
  if (applyDurLimit){
    cout << "INFO: only runs with a length of " << durMin << " to "  << durMax << " minutes will be considered in the analysis" << endl;
  }
  if (appBC){
    cout << "INFO: will be using bad channel map" << endl;
  }
  Bool_t enable2DWriting        = kFALSE;
  TH1D* referenceLEDDiffLEDMon  = NULL;
  if (iterFile.CompareTo("") != 0){
    enable2DWriting             = kTRUE;
    TFile* fileIter             = new TFile(iterFile);
    referenceLEDDiffLEDMon      = (TH1D*)fileIter->Get("LedDiffLedMonVsCellID");
  }
  Bool_t runSingleSM            = kFALSE;
  if (singleSM != -1 && (singleSM >= 0 && singleSM < 21))
      runSingleSM               = kTRUE;

  // initialize info from tree created by $ALICE_PHYSICS/PWGPP/EMCAL/TeCMacros/createTree.C
  TCalInfo *info          = 0;
  TFile *in               = TFile::Open(ifile,"read");
  TTree *treeRuns         = (TTree*)in->Get("tcal");
  treeRuns->SetBranchAddress("event",&info);
  treeRuns->Branch("event", &info, 32000, 99);
  Int_t Nev               = treeRuns->GetEntries();

  TCalCellInfo *cellinfo  = NULL;
  TTree *treeCells        = (TTree*)in->Get("tcalcell");
  treeCells->SetBranchAddress("cells",&cellinfo);
  treeCells->Branch("cells", &cellinfo, 32000, 99);
  Int_t nCellsTree        = treeCells->GetEntries();

  Int_t minRunNo = 1e6;
  Int_t maxRunNo = -1;
  for (Int_t i=0;i<Nev;++i) {
    treeRuns->GetEvent(i);
    if (info->fRunNo < minRunNo) minRunNo = info->fRunNo;
    if (info->fRunNo > maxRunNo) maxRunNo = info->fRunNo;
  }
  cout << minRunNo << "\t" << maxRunNo << endl;

  TProfile *gLedVsT[kSM];
  TProfile *gLedMonVsT[kSM];
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
    if (runSingleSM && i != singleSM )
        continue;
    gLedVsT[i] = new TProfile("","Led info;T;",nBinsT,17, 27);
    gLedVsT[i]->SetName(Form("ledsm%d",i));
    gLedMonVsT[i] = new TProfile("","Led info;T;",nBinsT,17, 27);
    gLedMonVsT[i]->SetName(Form("ledmonsm%d",i));
    gRatVsT[i] = new TProfile("","Led/LedMon;T",nBinsT,17, 27);
    gRatVsT[i]->SetName(Form("ledovermonsm%d",i));
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

  Bool_t isRefRun   = kFALSE;
  Bool_t hadRefRun  = kFALSE;

  cout << "-> there are " << Nev << " contained in this tree, starting to analyse them" << endl;

  for (Int_t i=0;i<Nev;++i) {
    treeRuns->GetEvent(i);
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
      if (runSingleSM && sm != singleSM )
        continue;
      hAverageT[sm]->SetBinContent(i+1,smInfot->fAvgT);
      hAverageTSorted[sm]->SetBinContent(hAverageTSorted[sm]->FindBin(info->fRunNo),smInfot->fAvgT);
    }
    for (Int_t sm=0;sm<sms.GetEntries();++sm) {
      TCalSM *smInfot = static_cast<TCalSM*>(sms.At(sm));
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
  }

  TFile *out =  NULL;
  if (runSingleSM)
    out       = TFile::Open(ofile,"update");
  else
    out       = TFile::Open(ofile,"recreate");
  for (Int_t i=0;i<kSM;++i) {
    if (runSingleSM && i != singleSM )
      continue;
    hAverageT[i]->Write();
    hAverageTSorted[i]->Write();
  }
  hAverageTPerSM->Write(hAverageTPerSM->GetName(),TObject::kOverwrite);
  for (Int_t sens = 0; sens < 160; sens++){
    hSensorsT[sens]->Write(hSensorsT[sens]->GetName(),TObject::kOverwrite);
    hSensorsTSorted[sens]->Write(hSensorsTSorted[sens]->GetName(),TObject::kOverwrite);
  }


  Int_t smcurr          = -1;
  Int_t cellIDFirstInSM = 0;

  for (Int_t j = 0; j < nCellsTree; j++){
    treeCells->GetEvent(j);

    if (cellinfo->fCellId > kNcells)
      continue;
    if (cellinfo->fCellId != j)
      cout << "WARNING: cellID and counter are out of sync" << endl;

    if (cellinfo->fCellId%50 == 0)
      cout << "next 50 cells: " << cellinfo->fCellId << "/" << kNcells<< endl;

    Int_t iTower = -1, iIphi = -1,  iIeta = -1, sm=-1;
    g->GetCellIndex(cellinfo->fCellId,sm,iTower,iIphi,iIeta);

    TRandom3 ledRandom(1289790);
    TRandom3 ledMonRandom(909088);

    TClonesArray &runs  = cellinfo->fRuns;

    if (runSingleSM && sm != singleSM )
      continue;

    if (smcurr != sm){
      smcurr = sm;
      out->mkdir(Form("cellsInSM%d",sm));
      out->cd(Form("cellsInSM%d",sm));
      if (debugInfo) cout << "SM: \t" << sm-1 << "\t"<< cellIDFirstInSM << "\t" << j << "\t" << j -cellIDFirstInSM<< endl;
      cellIDFirstInSM = j;
    }

    vector<Double_t> vecCellInfo[4];//            = new vector<Double_t>[4];     // T x led x ledM x RunNo
    Int_t nEntriesCell                       = 0;

    TProfile *gRatCellVsT = new TProfile("",Form("Led/LedMon cell ID%i ;T;",j),nBinsT,17, 27,opt);
    gRatCellVsT->SetName(Form("ledovermonCell%d",j));

    for (Int_t r=0;r<runs.GetEntries();++r) {
      TCalRun *run = static_cast<TCalRun*>(runs.At(r));

//       if (r%50 == 0)
//         cout << "starting with run " << r << "/" << Nev << endl;

      Int_t badcell     = run->fBad;
      Double_t ledM     = run->fLedM;
      Double_t ledR     = run->fLedR;
      Double_t monM     = run->fMonM;
      Double_t monR     = run->fMonR;
      Double_t smT      = run->fSMT;
      ULong_t deltaTimeS= (ULong_t)run->fLength/10;
      Int_t runNo       = run->fRunNo;

      if (runSingleSM && sm != singleSM )
          continue;

      if (appBC && badcell > 0){
        if (debugInfo > 1) cout << "found bad cell for " << cellinfo->fCellId <<  "\t " << runNo << endl;
        continue;
      }
      if ((smT<5)||(smT>45))
        continue;

      for ( UInt_t s = 0; s < deltaTimeS; s++){
        Double_t ledcurrent     = 0;
        Double_t ledMoncurrent  = 0;
        if (!((ledM<=0)||(ledR<=0)) ){
          ledcurrent              = ledRandom.Gaus(ledM,ledR);
          gLedVsT[sm]->Fill(smT,ledcurrent);
          gCellIdVsLed->Fill(cellinfo->fCellId,ledcurrent);
        }
        if (!((monM<=0)||(monR<=0)) ){
          ledMoncurrent         = ledMonRandom.Gaus(monM,monR);
          gLedMonVsT[sm]->Fill(smT,ledMoncurrent);
          gCellIdVsLedMon->Fill(cellinfo->fCellId,ledMoncurrent);
        }
        if ( (ledcurrent<=0) || (ledMoncurrent<=0) )
          continue;
        gCellIdVsRat->Fill(cellinfo->fCellId,ledcurrent/ledMoncurrent);
        gRatVsT[sm]->Fill(smT,ledcurrent/ledMoncurrent);
        gRatCellVsT->Fill(smT,ledcurrent/ledMoncurrent);

        vecCellInfo[0].push_back(smT);
        vecCellInfo[1].push_back(ledcurrent);
        vecCellInfo[2].push_back(ledMoncurrent);
        vecCellInfo[3].push_back(runNo);
        nEntriesCell++;
      }
    }
    if (nEntriesCell > 0){
      TGraph* graphRatVsT      = new TGraph(nEntriesCell);
      graphRatVsT->SetName(Form("graphRatvsT_Cell%d", cellinfo->fCellId));
      TGraph* graphLedVsT      = new TGraph(nEntriesCell);
      graphLedVsT->SetName(Form("graphLedvsT_Cell%d", cellinfo->fCellId));
      TGraph* graphLedMonVsT   = new TGraph(nEntriesCell);
      graphLedMonVsT->SetName(Form("graphLedMonvsT_Cell%d", cellinfo->fCellId));
      TGraph* graphRatVsRun    = new TGraph(nEntriesCell);
      graphRatVsRun->SetName(Form("graphRatvsRun_Cell%d", cellinfo->fCellId));
      for (Int_t k = 0; k< nEntriesCell; k++){
        graphRatVsT->SetPoint(k, vecCellInfo[0].at(k), vecCellInfo[1].at(k)/vecCellInfo[2].at(k));
        graphLedVsT->SetPoint(k, vecCellInfo[0].at(k), vecCellInfo[1].at(k));
        graphLedMonVsT->SetPoint(k, vecCellInfo[0].at(k), vecCellInfo[2].at(k));
        graphRatVsRun->SetPoint(k, vecCellInfo[3].at(k), vecCellInfo[1].at(k)/vecCellInfo[2].at(k));
      }
      graphRatVsT->Sort();
      graphLedVsT->Sort();
      graphLedMonVsT->Sort();
      graphRatVsRun->Sort();

      if (gRatCellVsT->GetEntries() > 0)
        gRatCellVsT->Write();
      if (graphRatVsT){
        graphRatVsT->Write();
      }
      if (graphRatVsRun){
        graphRatVsRun->Write();
      }
      if (graphLedVsT){
        graphLedVsT->Write();
      }
      if (graphLedMonVsT){
        graphLedMonVsT->Write();
      }
      delete graphLedVsT;
      delete graphRatVsT;
      delete graphRatVsRun;
      delete graphLedMonVsT;
    }
    delete gRatCellVsT;
  }

  for (Int_t i=0;i<kSM;++i) {
    if (runSingleSM && i != singleSM )
        continue;
    gLedVsT[i]->Write();
    gLedMonVsT[i]->Write();
    gRatVsT[i]->Write();
    hAverageT[i]->Write();
    hAverageTSorted[i]->Write();
  }
  if (runSingleSM){
    if (out->Get("LedDiffLedMonVsCellID")){
      TProfile* tempgCellIdVsRat  = (TProfile*)out->Get("LedDiffLedMonVsCellID");
      TH2D* tempgCellIdVsLed      = (TH2D*)out->Get("LedVsCellID");
      TH2D* tempgCellIdVsLedMon   = (TH2D*)out->Get("LedMonVsCellID");
      gCellIdVsRat->Add(gCellIdVsRat);
      gCellIdVsLed->Add(tempgCellIdVsLed);
      gCellIdVsLedMon->Add(tempgCellIdVsLedMon);
      gCellIdVsRat->Write(gCellIdVsRat->GetName(),TObject::kOverwrite);
      gCellIdVsLed->Write(gCellIdVsLed->GetName(),TObject::kOverwrite);
      gCellIdVsLedMon->Write(gCellIdVsLedMon->GetName(),TObject::kOverwrite);
    } else {
      gCellIdVsRat->Write();
      gCellIdVsLed->Write();
      gCellIdVsLedMon->Write();
    }
  } else {
    hAverageTPerSM->Write();
    gCellIdVsRat->Write();
    gCellIdVsLed->Write();
    gCellIdVsLedMon->Write();
  }
}
