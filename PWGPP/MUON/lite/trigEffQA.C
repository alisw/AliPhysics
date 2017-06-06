
#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TList.h"
#include "TObjString.h"
#include "TString.h"
#include "TGrid.h"
#include "TArrayD.h"
#include "TArrayI.h"
#include "TMap.h"
#include "TGridResult.h"
#include "TF1.h"
#include "TPad.h"
#include "TLatex.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliCDBStorage.h"
#include "AliMTRChEffAnalysis.h"
#include "AliMUONTriggerUtilities.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONVDigit.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONCalibrationData.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliCounterCollection.h"
#include "AliTriggerConfiguration.h"
#endif

const Int_t kNch = 4;
const Double_t kZero = 1.e-7; // Avoid problems when comparing to 0.

//_____________________________________________________________________________
void SetMyStyle()
{
  /// Set graphic style
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(10);

  gStyle->SetTitleXSize(0.03);
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYSize(0.03);
  gStyle->SetTitleYOffset(1.9);

  gStyle->SetMarkerSize(0.7);
  gStyle->SetHistLineWidth(2);

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.03);

  gROOT->ForceStyle();
}

//_____________________________________________________________________________
Bool_t IsRunNum ( TString stringToken )
{
  return ( stringToken.IsDigit() && stringToken.Length()>=6 && stringToken.Length()<=9 );
}


//_____________________________________________________________________________
void SetRunAxisRange ( TH1* histo, TString runAxis = "X" )
{
  /// Set axis range
  runAxis.ToUpper();
  TAxis* axis = histo->GetXaxis();
  if ( runAxis == "Y" ) axis = histo->GetYaxis();
  else if ( runAxis == "Z" ) axis = histo->GetZaxis();
  histo->LabelsOption("v",runAxis.Data());
  axis->SetLabelSize(0.02);
  for ( Int_t ibin=1; ibin<=axis->GetNbins(); ibin++ ) {
    TString binLabel = axis->GetBinLabel(ibin);
    if ( ! binLabel.IsNull()) continue;
    axis->SetRange(1, ibin-1);
    return;
  }
}

//_____________________________________________________________________________
Int_t GetRunNumber(TString filePath)
{
  /// Get run number from file path
  TObjArray* array = filePath.Tokenize("/");
  array->SetOwner();
  TString auxString = "";
  Int_t runNum = -1;
  for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
    auxString = array->At(ientry)->GetName();
    if ( IsRunNum(auxString) ) {
      runNum = auxString.Atoi();
      break;
    }
  }
  delete array;

  if ( runNum < 0 ) {
    array = auxString.Tokenize("_");
    array->SetOwner();
    auxString = array->Last()->GetName();
    auxString.ReplaceAll(".root","");
    if ( IsRunNum(auxString) ) runNum = auxString.Atoi();
    delete array;
  }

  return runNum;
}

//_____________________________________________________________________________
Double_t* GetProdErr(Double_t* effErr, Int_t exclude, Int_t nFactors = kNch)
{
  /// Error of product
  Double_t prod = 1.;
  Double_t relErr = 0., relProdErrSquare = 0.;
  for ( Int_t iprod=0; iprod<nFactors; iprod++ ) {
    if ( iprod == exclude ) continue;
    prod *= effErr[iprod];
    relErr = ( effErr[iprod] > kZero ) ? effErr[iprod+nFactors]/effErr[iprod] : 0.;
    relProdErrSquare += relErr*relErr;
    //printf("%f +- %f  ", effErr[iprod], effErr[iprod+nFactors]); // alBER TO CUT
  }
  Double_t* prodErr = new Double_t[2];
  prodErr[0] = prod;
  prodErr[1] = prod*TMath::Sqrt(relProdErrSquare);
  //printf("-> %f %f\n", prodErr[0], prodErr[1]); // REMEMBER TO CUT
  return prodErr;
}


//_____________________________________________________________________________
Double_t* GetConditionalEffErr(Double_t* effErr1, Double_t* effErr2, Double_t* effErrBoth, Int_t exclude = -1)
{
  /// Error on conditional efficiency
  Double_t* effErr = new Double_t[2*kNch];
  for ( Int_t ich=0; ich<kNch; ich++ ) {
    if ( ich == exclude ) {
      effErr[ich] = ( effErr1[ich] < 1. ) ? ( effErr2[ich] - effErrBoth[ich] ) / ( 1. - effErr1[ich] ) : 0.;
      effErr[ich+kNch] = 0;
      if ( effErr1[ich] < 1. ) {
        Double_t err2 = effErr2[ich+kNch] / ( 1. - effErr1[ich] );
        Double_t errBoth = effErrBoth[ich+kNch] / ( 1. - effErr1[ich] );
        Double_t err1 = effErr1[ich+kNch] * effErr[ich] / ( 1. - effErr1[ich] );
        effErr[ich+kNch] = TMath::Sqrt(err2*err2 + errBoth*errBoth + err1*err1);
      }
    }
    else {
      effErr[ich] = ( effErr1[ich] > kZero ) ? effErrBoth[ich]/effErr1[ich] : 0.;
      Double_t relErr1 = ( effErr1[ich] > kZero ) ? effErr1[ich+kNch]/effErr1[ich] : 0.;
      Double_t relErrBoth = ( effErrBoth[ich] > kZero ) ? effErrBoth[ich+kNch]/effErrBoth[ich] : 0.;
      effErr[ich+kNch] = effErr[ich] * TMath::Sqrt(relErr1*relErr1 + relErrBoth*relErrBoth);
    }
    //printf("%f  %f  %f -> %f\n", effErr1[ich], effErr2[ich], effErrBoth[ich], effErr[ich]); // REMEMBER TO CUT
  } // loop on chambers
  return effErr;
}


//_____________________________________________________________________________
Double_t* GetBinomial(Double_t* effErr1, Double_t* effErr2 = 0x0, Double_t* effErrBoth = 0x0)
{
  /// Binomial error
  Double_t effProd[4];
  Double_t defaultEffErr[2] = {1.,0.};
  Double_t* auxBinomial = 0x0;
  Double_t* currEffErr44 = 0x0;
  Double_t* effErrBinomial = new Double_t[2];
  effErrBinomial[0] = 0.;
  effErrBinomial[1] = 0.;

  for ( Int_t ich = -1; ich<kNch; ich++ ) {
    Double_t* currEffErr = GetProdErr(effErr1, ich);
    if ( ich >= 0 ) {
      currEffErr[0] = currEffErr[0] - currEffErr44[0];
      currEffErr[1] = TMath::Sqrt(currEffErr[1]*currEffErr[1] + currEffErr44[1]*currEffErr44[1]);
    }
    if ( effErr2 ) {
      Double_t* auxEffErr = GetConditionalEffErr(effErr1, effErr2, effErrBoth, ich);
      auxBinomial = GetBinomial(auxEffErr);
      delete [] auxEffErr;
    }
    for ( Int_t ival=0; ival<2; ival++ ) {
      effProd[2*ival] = currEffErr[ival];
      effProd[2*ival+1] = ( effErr2 ) ? auxBinomial[ival] : defaultEffErr[ival];
    }
    if ( ich < 0 ) currEffErr44 = currEffErr;
    else delete [] currEffErr;
    delete [] auxBinomial;

    Double_t* effErr = GetProdErr(effProd, -1, 2);
    //printf("%f * %f = %f\n", effProd[0], effProd[1], effErr[0]); // REMEMBER TO CUT
    effErrBinomial[0] += effErr[0];
    effErrBinomial[1] += effErr[1]*effErr[1];
    delete [] effErr;
  } // loop on chambers

  delete [] currEffErr44;

  effErrBinomial[1] = TMath::Sqrt(effErrBinomial[1]);

  return effErrBinomial;
}

//_____________________________________________________________________________
Bool_t CheckOCDBFile ( TString cdbDir, Int_t runNum )
{
  /// Check if (run-by-run) CDB object is there
  /// This is needed for example for the scalers
  /// when the default storage is cvmfs.
  /// Indeed, the cvmfs OCDB has sometimes synchro problem
  /// and the latest files are not copied
  if ( runNum >= 0 ) AliCDBManager::Instance()->SetRun(runNum);
  TList* list = AliCDBManager::Instance()->GetAll(cdbDir.Data());
  if ( list->GetEntries() == 0 ) {
    printf("Warning: no entry found in %s for run %i\n",cdbDir.Data(),runNum);
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
TList* GetOCDBList ( TString ocdbDirs )
{
  /// Get list of CDB objetcs
  TString storageType = AliCDBManager::Instance()->GetDefaultStorage()->GetType();
  Bool_t isGrid = storageType.Contains("alien");
  TString baseFolder = AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder();

  TList* outList = new TList();
  outList->SetOwner();
  TObjArray* dirNameList = ocdbDirs.Tokenize(",");
  for ( Int_t idir=0; idir<dirNameList->GetEntries(); idir++ ) {
    TString fullPath = Form("%s/%s",baseFolder.Data(),dirNameList->At(idir)->GetName());
    if ( isGrid ) {
      TGridResult *res = gGrid->Ls(fullPath.Data());
      if (!res) return 0x0;
      for ( Int_t ires=0; ires<res->GetEntries(); ires++ ) {
        TString currFile = static_cast<TMap*>(res->At(ires))->GetValue("name")->GetName();
        outList->Add(new TObjString(currFile));
      }
      delete res;
    }
    else {
      TString fileListStr = gSystem->GetFromPipe(Form("ls %s",fullPath.Data()));
      TObjArray* fileList = fileListStr.Tokenize("\n");
      for ( Int_t ires=0; ires<fileList->GetEntries(); ires++ ) {
        TString currFile = fileList->At(ires)->GetName();
        outList->Add(new TObjString(currFile));
      }
      delete fileList;
    }
  }
  delete dirNameList;
  return outList;
}

//_____________________________________________________________________________
Bool_t SetAndCheckOCDB ( TString defaultStorage )
{
  /// Set the default storage and check if it is ok
  if ( defaultStorage.IsNull() ) {
    printf("Default storage not specified. Nothing done\n");
    return kFALSE;
  }

  if ( AliCDBManager::Instance()->IsDefaultStorageSet() ) {
    printf("Default storage already set: nothing done\n");
    return kTRUE;
  }

  if ( defaultStorage.Contains("alien://") || defaultStorage.Contains("raw://") ) {
    if ( ! gGrid ) TGrid::Connect("alien://");
    if ( ! gGrid ) {
      printf("Error: Problem connetting to grid: default storage not set\n");
      return kFALSE;
    }
  }

  AliCDBManager::Instance()->SetDefaultStorage(defaultStorage.Data());

  if ( defaultStorage.Contains("raw://") ) return kTRUE;

  Bool_t isOk = kTRUE;

  if ( AliCDBManager::Instance()->IsDefaultStorageSet() ) {
    TString searchDir = "MUON/Calib/MappingData";
    TString fullPath = Form("%s/%s",defaultStorage.Data(),searchDir.Data());
    TList* ocdbList = GetOCDBList(searchDir);
    if ( ocdbList->GetEntries() == 0 ) {
      printf("No entries in %s\n",fullPath.Data());
      isOk = kFALSE;
    }
    else {
      TString checkFile = Form("%s/%s",fullPath.Data(),ocdbList->At(0)->GetName());
      checkFile.ReplaceAll("local://","");
      checkFile.ReplaceAll("folder=","");
      checkFile.ReplaceAll("Folder=","");
      TFile* file = TFile::Open(checkFile.Data());
      if ( ! file ) {
        printf("Cannot access test file: %s\n", checkFile.Data());
        isOk = kFALSE;
      }
      delete file;
    }
    delete ocdbList;
  }
  else {
    printf("Tried to set the default storage, but something went wrong.\n");
    isOk = kFALSE;
  }

  if ( ! isOk ) printf("Please check path %s\n",defaultStorage.Data());

  return  isOk;
}

//_____________________________________________________________________________
Bool_t IsOCDBChanged ( Int_t currRun, Int_t previousRun, TList* fileList )
{
  /// Check if the OCDB object is changed w.r.t. the previous run
  if ( ! fileList ) return kTRUE;
  for ( Int_t ifile=0; ifile<fileList->GetEntries(); ifile++ ) {
    TString filename = static_cast<TObjString*>(fileList->At(ifile))->GetString();
    filename.ReplaceAll("Run","");
    TObjArray* array = filename.Tokenize("_");
    Int_t firstRun = static_cast<TObjString*>(array->At(0))->GetString().Atoi();
    Int_t lastRun = static_cast<TObjString*>(array->At(1))->GetString().Atoi();
    delete array;
    Bool_t isCurrRunInside = ( currRun >= firstRun && currRun <= lastRun );
    Bool_t isPreviousRunInside = ( previousRun >= firstRun && previousRun <= lastRun );
    if ( isCurrRunInside != isPreviousRunInside ) return kTRUE;
  }
  return kFALSE;
}


//_____________________________________________________________________________
void TrigEffTrending ( TString fileNameList, TList& outList)
{

  TString physSel = "PhysSelPass,PhysSelReject";
  TString trigClass = "ANY";
  TString centrClass = "-5_105";
  Int_t matchTrig = AliTrigChEffOutput::kMatchApt;
  for ( Int_t ieff=0; ieff<2; ieff++ ) {
    Int_t trackSel = ( ieff == 0 ) ? AliTrigChEffOutput::kNoSelectTrack : AliTrigChEffOutput::kSelectTrack;
    Int_t effType = ( ieff == 0 ) ? AliTrigChEffOutput::kEffFromTrig : AliTrigChEffOutput::kEffFromTrack;

    AliMTRChEffAnalysis an;
    an.SetEffConditions(physSel,trigClass,centrClass,trackSel,matchTrig,effType);
    an.InitFromLocal(fileNameList);

    Int_t nCanvases = gROOT->GetListOfCanvases()->GetEntries();

    an.DrawEffTrend(AliTrigChEffOutput::kHchamberEff,-1,3.,0.9,1.05);
    an.DrawEffTrend(AliTrigChEffOutput::kHslatEff,-1,3.,0.8,1.05);

    TString baseName = ( ieff == 0 ) ? "Trigger" : "Tracker";
    for ( Int_t ican=nCanvases;   ican<gROOT->GetListOfCanvases()->GetEntries(); ican++ ) {
      TCanvas* can = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->At(ican));
      can->SetName(Form("%s_from%s",baseName.Data(),can->GetName()));
      can->SetTitle(Form("%s_from%s",baseName.Data(),can->GetTitle()));
      can->cd();
      TPad* pad = new TPad(Form("%s_text",can->GetName()),"",0.35,0.97,0.65,1.);
      pad->SetFillStyle(4000);
      pad->Draw();
      pad->cd();
      TLatex tex;
      tex.SetTextSize(0.75);
      tex.DrawLatexNDC(0.1,0.1,Form("From %s track",baseName.Data()));
      outList.Add(can);
    } // loop on canvases
  } // loop on efficiency types
}

//_____________________________________________________________________________
void MaskTrending ( TObjArray runNumArray, TString defaultStorage, TList& outList )
{
  /// Get the masks vs. run number

  if ( ! SetAndCheckOCDB(defaultStorage) ) return;

  TObjArray maskedList(8);
  TObjArray auxList(8);
  auxList.SetOwner();
  TString histoName = "", histoTitle = "";
  for(Int_t icath=0; icath<2; icath++){
    TString cathName = ( icath==0 ) ? "bendPlane" : "nonBendPlane";
    for(Int_t ich=0; ich<kNch; ich++){
      histoName = Form("%sMaskCh%i", cathName.Data(), 11+ich);
      histoTitle = Form("Chamber %i - %s: fraction of masked channels", 11+ich, cathName.Data());
      TH2* histo = new TH2D(histoName.Data(), histoTitle.Data(),1,0.,1., 234, 0.5, 234. + 0.5);
      histo->GetYaxis()->SetTitle("Board Id");
//      histo->SetOption("COLZ");
      Int_t imask = 2*ich + icath;
      maskedList.AddAt(histo, imask);
      auxList.AddAt(histo->Clone(Form("%s_aux",histoName.Data())), imask);
    } // loop on chambers
  } // loop on cathodes

  TArrayS xyPatternAll[2];
  for(Int_t icath=0; icath<2; icath++){
    xyPatternAll[icath].Set(kNch);
    xyPatternAll[icath].Reset(0xFFFF);
  }

  TList* ocdbFileList = 0x0;
  Int_t previousRun = -1;
  AliMUONDigitMaker* digitMaker = 0x0;
  AliMUONDigitStoreV2R digitStore;

  AliMUONCalibrationData* calibData = 0x0;
  AliMUONTriggerUtilities* trigUtilities = 0x0;
  for ( Int_t irun=0; irun<runNumArray.GetEntries(); irun++ ) {
    TString runNumString = runNumArray.At(irun)->GetName();
    Int_t runNumber = runNumString.Atoi();

    if ( IsOCDBChanged(runNumber, previousRun, ocdbFileList) ) {
      AliCDBManager::Instance()->SetRun(runNumber);

      if ( ! digitMaker ) {
        digitMaker = new AliMUONDigitMaker(kFALSE);
        // Create a store with all digits in trigger
        for ( Int_t iboard=1; iboard<=234; iboard++ ) {
          digitMaker->TriggerDigits(iboard, xyPatternAll, digitStore, kFALSE);
        }
      }

      if ( ! ocdbFileList ) ocdbFileList = GetOCDBList("MUON/Calib/GlobalTriggerCrateConfig,MUON/Calib/RegionalTriggerConfig,MUON/Calib/LocalTriggerBoardMasks");

      delete calibData;
      calibData = new AliMUONCalibrationData (runNumber);
      delete trigUtilities;
      trigUtilities = new AliMUONTriggerUtilities (calibData);
    }

    previousRun = runNumber;

    TIter next(digitStore.CreateIterator());
    AliMUONVDigit* dig = 0x0;
    while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) ) {
      Int_t icath = dig->Cathode();
      Int_t detElemId = dig->DetElemId();
      Int_t ich = detElemId/100-11;
      Int_t iboard = dig->ManuId();
      Int_t imask = 2*ich + icath;
      static_cast<TH2*>(auxList.At(imask))->Fill(runNumString.Data(),iboard,1.);
      static_cast<TH2*>(maskedList.At(imask))->Fill(runNumString.Data(),iboard,(Double_t)trigUtilities->IsMasked(*dig));
    }
  } // loop on runs
  delete calibData;
  delete trigUtilities;
  delete digitMaker;

  TString canName = "";
  for ( Int_t imask=0; imask<maskedList.GetEntries(); imask++ ) {
    TH2* histo = static_cast<TH2*>(maskedList.At(imask));
    histo->Divide(static_cast<TH2*>(auxList.At(imask)));
    SetRunAxisRange(histo);

    canName = Form("%sCan", histo->GetName());
    TCanvas* can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
    can->SetRightMargin(0.14);
    histo->SetStats(kFALSE);
    if ( histo->GetMinimum() == 0. && histo->GetMaximum() == 0. ) histo->SetMaximum(0.1);
    histo->DrawCopy("COLZ");
    outList.Add(can);
  }
}

//_____________________________________________________________________________
Bool_t CheckPattern ( TString trigName, TObjArray* keepArray, TObjArray* rejectArray )
{
  /// Check pattern
  for ( Int_t ipat=0; ipat<rejectArray->GetEntries(); ++ipat ) {
    if ( trigName.Contains(rejectArray->At(ipat)->GetName() ) ) return kFALSE;
  } // loop on reject pattern

  for ( Int_t ipat=0; ipat<keepArray->GetEntries(); ++ipat ) {
    if ( trigName.Contains(keepArray->At(ipat)->GetName() ) ) return kTRUE;
  } // loop on keep pattern

  return ( keepArray->GetEntries() == 0 ) ? kTRUE : kFALSE;
}

//_____________________________________________________________________________
TObjArray* BuildListOfTrigger ( const TObjArray* triggerArray, TString keepPattern = "", TString rejectPattern="OTHER,TRUE,PHI,ANY,EMC,-ACE-,-ABCE-,WU,MUP,SPI,SHM" )
{
  /// Build list of trigger classes
  TObjArray* selectedList = new TObjArray();
  selectedList->SetOwner();
  TObjArray* rejectArray = rejectPattern.Tokenize(",");
  TObjArray* keepArray = keepPattern.Tokenize(",");

  for ( Int_t iTrig = 0; iTrig < triggerArray->GetEntries(); iTrig++ ){
    TString currTrigName = ((TObjString*)triggerArray->At(iTrig))->GetName();
    if ( CheckPattern(currTrigName, keepArray, rejectArray) ) selectedList->AddLast(new TObjString(currTrigName.Data()));
  }

  delete rejectArray;
  delete keepArray;

  return selectedList;

}

//_____________________________________________________________________________
TString FindCorrespondingTrigger ( TString checkTrigger, TObjArray* triggerArray )
{
  /// Find trigger from pattern
  TString foundName = "";
  for ( Int_t iTrig = 0; iTrig < triggerArray->GetEntries(); iTrig++ ){
    TString currTrigName = ((TObjString*)triggerArray->At(iTrig))->GetName();
    TObjArray* array = currTrigName.Tokenize("-");
    TString collisionType = array->At(1)->GetName();
    delete array;
    collisionType.Append("-");
    collisionType.Prepend("-");
    if ( checkTrigger.Contains(collisionType.Data()) ) {
      foundName = currTrigName;
      break;
    }
  }

  return foundName;
}

//_____________________________________________________________________________
void ScalerTrending ( TObjArray runNumArray, TString mergedFileName, TString defaultStorage, TList& outList )
{
  /// Get the scalers vs. run number
  if ( ! SetAndCheckOCDB(defaultStorage) ) return;

  //trigger count from ESDs
  TFile *file = TFile::Open(mergedFileName.Data());
  AliCounterCollection* ccol = (AliCounterCollection*)((TDirectoryFile*)file->FindObjectAny("MUON_QA"))->FindObjectAny("eventCounters");

  //Build the trigger list for trigger with muon only in readout and min. bias triggers
  TString triggerListName = ccol->GetKeyWords("trigger");

  TObjArray selectedTriggerArray, selectedL0TriggerArray;
  selectedTriggerArray.SetOwner();
  selectedL0TriggerArray.SetOwner();

  const Int_t nScaler = 3;
  TString sScaler[nScaler] = {"L0B","L2A","L0BRATE"};
  enum eScaler {kL0B = 0, kL2A=1, kL0BRATE=2};
  Float_t maxScaler[nScaler] = {1e8,1e7,1e6};
  TObjArray hFromQA;
  TObjArray hFromScalers;
  TObjArray hOutput;

  TString cdbDir = "GRP/CTP/Config";


  TString sHistName, sHistNameFull, sTitleName;
  Int_t nRuns = runNumArray.GetEntries();

  //
  //Fill histos for Scalers and QA
  //
  //loop on run list
  for ( Int_t iRun = 0; iRun < runNumArray.GetEntries(); iRun++ ) {

    TString sRunNr = ((TObjString*)runNumArray.At(iRun))->GetString();
    Int_t runNr = sRunNr.Atoi();
    if ( ! CheckOCDBFile(cdbDir, runNr)) continue;
    AliAnalysisTriggerScalers triggerScaler(runNr);
    AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(triggerScaler.GetOCDBObject(cdbDir.Data(),runNr));
    const TObjArray& trClasses = tc->GetClasses();

    Int_t ibin = iRun+1;

    for ( Int_t itype=0; itype<2; itype++ ) {
      TObjArray* currSelectedList = ( itype == 0 ) ? &selectedTriggerArray : &selectedL0TriggerArray;
      TString matchTrig = ( itype == 0 ) ? "" : "C0TVX";
      TObjArray* selectedTrigArrayForRun = BuildListOfTrigger(&trClasses, matchTrig);

      //loop on trigger list
      for ( Int_t iTrig = 0; iTrig < selectedTrigArrayForRun->GetEntries(); iTrig++ ) {

        TString currTrigName = selectedTrigArrayForRun->At(iTrig)->GetName();
        if ( itype == 0 && ! triggerListName.Contains(currTrigName.Data()) ) continue;
        if ( ! currSelectedList->FindObject(currTrigName.Data()) ) currSelectedList->Add(new TObjString(currTrigName));

        //loop on scaler list
        for ( Int_t iScaler = 0; iScaler < nScaler; iScaler++ ) {

          if ( itype == 1 && iScaler != kL0B ) continue;

          //from Scalers
          TGraph* graph = triggerScaler.PlotTrigger(currTrigName.Data(),sScaler[iScaler].Data(),kFALSE);

          sHistName = Form("%s_%s",currTrigName.Data(),sScaler[iScaler].Data());
          sHistNameFull = Form("Scalers_%s",sHistName.Data());

          TH1* hist = (TH1*) hFromScalers.FindObject(sHistNameFull);
          if ( ! hist ) {
            hist = new TH1D(sHistNameFull,sHistName,nRuns,1.,1.+(Double_t)nRuns);
            hist->LabelsOption("v");
            hist->GetXaxis()->SetLabelSize(0.02);
            hist->SetDirectory(0);
            hist->SetMinimum(1);
            hist->SetMaximum(maxScaler[0]);
            hFromScalers.AddLast(hist);
            hOutput.AddLast(hist);
            if ( iScaler == kL2A ) {
              sHistNameFull = "QA_" + sHistName;
              hFromQA.AddLast(hist->Clone(sHistNameFull.Data()));
            }
          }
          Double_t *tab = (Double_t*) graph->GetY();
          if ( tab ) hist->SetBinContent(ibin,tab[0]);
          hist->GetXaxis()->SetBinLabel(ibin,sRunNr.Data());
          delete graph;

          //from QA
          if ( iScaler != kL2A ) continue;
          TH1* histCounters = static_cast<TH1*>(ccol->Get("run",Form("run:%s/trigger:%s",sRunNr.Data(),currTrigName.Data())));
          sHistNameFull = sHistNameFull = "QA_" + sHistName;
          hist = (TH1*) hFromQA.FindObject(sHistNameFull);
          if ( histCounters ) hist->SetBinContent(ibin,histCounters->GetSumOfWeights());
          hist->GetXaxis()->SetBinLabel(ibin,sRunNr.Data());
          delete histCounters;
        }//end loop on scaler list
      }//end loop on trigger list
    } // end loop on type
  }//end loop on run list


  if ( selectedTriggerArray.GetEntries() == 0 ) {
    printf("No trigger selected from trigger list %s\n",triggerListName.Data());
    return;
  }
  printf("Nr of triggers selected %i\n",selectedTriggerArray.GetEntries());

  printf("Nr of T0 triggers selected %i\n",selectedL0TriggerArray.GetEntries());

  //Set options for QA and Scalers histos

  for ( Int_t itype=0; itype<2; itype++ ) {
    TObjArray* currList = ( itype == 0 ) ? &hFromScalers : &hFromQA;
    for ( Int_t ihisto=0; ihisto<currList->GetEntriesFast(); ihisto++ ) {
      TH1* histo = static_cast<TH1*> ( currList->At(ihisto) );
      if (!histo) continue;
      // Write run number to each bin
      for ( Int_t iRun = 0; iRun < runNumArray.GetEntries(); iRun++ ) {
        TString sRunNr = ((TObjString*)runNumArray.At(iRun))->GetString();
        Int_t ibin = iRun+1;
        TString binLabel = histo->GetXaxis()->GetBinLabel(ibin);
        if ( ! binLabel.IsNull() ) continue;
        histo->GetXaxis()->SetBinLabel(ibin,sRunNr.Data());
      }
      histo->SetStats(kFALSE);
    }
  }


  //Loop on histos from scalers and QA and create resulting histos from scalers
  const Int_t nHisto = 3;
  TString sHisto[nHisto] = {"L0BoverL0BC0TVX","L2AoverL0B","L2AQAoverSCALERS"};
  TString sTitleHisto[nHisto] = {"L0B trigger / L0BC0TVX","L2A / L0B","L2A from QA / L2A from SCALERS"};
  //  TString sHisto[nHisto] = {"L2AoverL0B","L2AQAoverSCALERS"};

  //loop on trigger list
  for ( Int_t iTrig = 0; iTrig < selectedTriggerArray.GetEntries(); iTrig++ ) {

    sHistNameFull = Form("Scalers_%s_L0B",((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
    TH1* histo1 = static_cast<TH1*> ( hFromScalers.FindObject(sHistNameFull) );
    if (!histo1) continue;


    //C0TVX
    TString sTrig = ( (TObjString*) selectedTriggerArray.At(iTrig) )->GetName();
    TString sL0Trig = FindCorrespondingTrigger(sTrig, &selectedL0TriggerArray);

    sHistNameFull = Form("Scalers_%s_L0B",sL0Trig.Data());

    TH1* histo0 = static_cast<TH1*> ( hFromScalers.FindObject(sHistNameFull) );
    if ( histo0 ) {
      sHistNameFull = Form("%s_%s",sHisto[0].Data(),((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
      TH1* histo10 = (TH1*) histo1->Clone(sHistNameFull);
      histo10->SetTitle(sTitleHisto[0].Data());
      histo10->Sumw2();
      histo10->Divide(histo0);
      histo10->SetMaximum(10);
      histo10->SetMinimum(1e-5);
      //outList.Add(histo10);
      hOutput.AddLast(histo10);
      //outList.Add(histo0);
      //outList.Add(histo1);
    }

    //DEADTIME
    sHistNameFull = Form("Scalers_%s_L2A",((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
    TH1* histo2 = static_cast<TH1*> ( hFromScalers.FindObject(sHistNameFull) );
    if (!histo2) continue;

    sHistNameFull = Form("%s_%s",sHisto[1].Data(),((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
    TH1* histo3 = (TH1*) histo2->Clone(sHistNameFull);
    histo3->SetTitle(sTitleHisto[1]);
    histo3->Sumw2();
    histo3->Divide(histo1);
    histo3->SetMaximum(1.2);
    histo3->SetMinimum(1e-5);
    //outList.Add(histo3);
    hOutput.AddLast(histo3);

    //QA over Scalers
    sHistNameFull = Form("QA_%s_L2A",((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
    TH1* histo4 = static_cast<TH1*> ( hFromQA.FindObject(sHistNameFull) );
    if (!histo4) continue;

    sHistNameFull = Form("%s_%s",sHisto[2].Data(),((TObjString*) selectedTriggerArray.At(iTrig))->GetName());
    TH1* histo5 = (TH1*) histo4->Clone(sHistNameFull);
    histo5->SetTitle(sTitleHisto[2]);
    histo5->Sumw2();
    histo5->Divide(histo2);
    histo5->SetMaximum(1.2);
    histo5->SetMinimum(5e-1);
    //outList.Add(histo5);
    hOutput.AddLast(histo5);
  }

  // Plot all on canvases (only canvases will be saved)
  const Int_t nCanvases = nScaler + nHisto;
  TString sCanvases[nCanvases];
  for (Int_t iScaler = 0; iScaler < nScaler; iScaler++) sCanvases[iScaler] = sScaler[iScaler];
  for (Int_t iHisto = 0; iHisto < nHisto; iHisto++) sCanvases[nScaler+iHisto] = sHisto[iHisto];

  //loop on canvases
  for ( Int_t iCan = 0; iCan < nCanvases; iCan++) {
    TCanvas* canvas = new TCanvas(sCanvases[iCan],sCanvases[iCan],200,10,600,600);
    TLegend* leg  = new TLegend(0.72,0.7,0.9,0.85);
    leg->SetBorderSize(1);
    if ( iCan != 4 ) canvas->SetLogy();
    TString optDraw = "e";

    //loop on trigger list
    Int_t icolor = 1;
    for ( Int_t iTrig = 0; iTrig < selectedTriggerArray.GetEntries(); iTrig++ ) {

      if ( iCan < nScaler ) sHistNameFull = Form("Scalers_%s_%s",selectedTriggerArray.At(iTrig)->GetName(),sCanvases[iCan].Data());
      else sHistNameFull = Form("%s_%s",sCanvases[iCan].Data(),selectedTriggerArray.At(iTrig)->GetName());
      TH1* histo1 = static_cast<TH1*> ( hOutput.FindObject(sHistNameFull) );
      if (!histo1) continue;

      if ( icolor == 10 ) icolor++;
      histo1->SetLineColor(icolor++);
      histo1->Draw(optDraw);
      optDraw = "esame";

      leg->AddEntry(histo1,selectedTriggerArray.At(iTrig)->GetName(),"l");
    }

    leg->Draw();
    outList.Add(canvas);
  }

  file->Close();
}

//_____________________________________________________________________________
void trigEffQA(TString fileListName, TString outFilename = "", TString defaultStorage = "raw://", Bool_t doScalers = kFALSE, TString trackerQAmergedOut="QAresults_merged.root")
{
  /// Main function
  ifstream inFile(fileListName.Data());
  TObjArray fileNameArray, runNumArray;
  fileNameArray.SetOwner();
  runNumArray.SetOwner();
  TString currString = "";
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      currString.ReadLine(inFile); // Read line
      if ( ! currString.Contains(".root") ||
          currString.BeginsWith("#") ) continue;
      fileNameArray.AddLast(new TObjString(currString.Data()));
      Int_t runNum = GetRunNumber(currString);
      runNumArray.AddLast(new TObjString(Form("%i",runNum)));
    }
    inFile.close();
  }
  else {
    printf("Fatal: cannot open input file %s\n",fileListName.Data());
    return;
  }

  runNumArray.Sort();

  TList outList;
  TrigEffTrending(fileListName, outList);
  if ( SetAndCheckOCDB(defaultStorage) ) {
    MaskTrending(runNumArray, defaultStorage, outList);
    if ( doScalers ) {
      if ( gSystem->AccessPathName(trackerQAmergedOut.Data()) ) {
        printf("Warning: cannot perform scaler trending:\n merged QA from tracker\n  %s\n  does not exist\n",trackerQAmergedOut.Data());
      }
      else {
        ScalerTrending(runNumArray, trackerQAmergedOut, defaultStorage, outList);
      }
    }
  }

  if ( outFilename.IsNull() ) return;

  TString outCanName = outFilename;
  outCanName.ReplaceAll(".root",".pdf");
  for ( Int_t ican=0; ican<outList.GetEntries(); ican++ ) {
    TString canName = outCanName;
    if ( ican == 0 ) canName.Append("("); // open pdf file
    else if ( ican == outList.GetEntries()-1 ) canName.Append(")"); // close pdf file
    static_cast<TCanvas*>(outList.At(ican))->Print(canName.Data());
  }
  // There is a bug when creating a pdf
  // So create a ps and then convert via epstopdf
  if ( outCanName.Contains(".ps") || outCanName.Contains(".eps") ) {
    gSystem->Exec(Form("epstopdf %s", outCanName.Data()));
    gSystem->Exec(Form("rm %s", outCanName.Data()));
  }

  TFile* outFile = new TFile(outFilename.Data(), "recreate");
  outList.Write();
  outFile->Close();
}
