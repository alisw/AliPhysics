
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

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliCDBStorage.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"
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
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadTopMargin(0.08);
  
  gROOT->ForceStyle();
}


//_____________________________________________________________________________
void SetRunAxisRange ( TAxis* axis )
{
  /// Set axis range
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
    if ( auxString.IsDigit() && auxString.Length()>=6 && auxString.Length()<=9 ) {
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
    runNum = auxString.Atoi();
    delete array;
  }
  
  return runNum;
}

//_____________________________________________________________________________
Bool_t ChangeFilenames ( TObjArray &fileNameArray )
{
  /// Use custom output
  /// We used to perform the QA on the MTR chamber efficiency
  /// but since it is obtained form tracks matching with the tracker
  /// it is not that good for QA since we are dependent on the tracker status.
  /// In recent versions, the task also calculates the efficiency from all trigger tracks
  /// (including ghosts). Analyse this output instead.
  for ( Int_t ifile=0; ifile<fileNameArray.GetEntries(); ifile++ ) {
    TObjString* currObjString = static_cast<TObjString*>(fileNameArray.At(ifile));
    TString currFile = currObjString->GetString();
    TString dirName = gSystem->DirName(currFile.Data());
    TString fileName = gSystem->BaseName(currFile.Data());
    Int_t runNum = GetRunNumber(fileName);
    TString newFilename = Form("%s/terminateRuns/%i/trigChEff_ANY_Apt_allTrig.root",dirName.Data(),runNum);
    if ( ! gSystem->AccessPathName(newFilename.Data()) ) {
      printf("New output not found. Use the standard efficiency instead\n");
      return kFALSE;
    }
    currObjString->SetString(newFilename);
  }
  return kTRUE;
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
      delete auxEffErr;
    }
    for ( Int_t ival=0; ival<2; ival++ ) {
      effProd[2*ival] = currEffErr[ival];
      effProd[2*ival+1] = ( effErr2 ) ? auxBinomial[ival] : defaultEffErr[ival];
    }
    if ( ich < 0 ) currEffErr44 = currEffErr;
    else delete currEffErr;
    delete auxBinomial;
    
    Double_t* effErr = GetProdErr(effProd, -1, 2);
    //printf("%f * %f = %f\n", effProd[0], effProd[1], effErr[0]); // REMEMBER TO CUT
    effErrBinomial[0] += effErr[0];
    effErrBinomial[1] += effErr[1]*effErr[1];
    delete effErr;
  } // loop on chambers
  
  delete currEffErr44;
  
  effErrBinomial[1] = TMath::Sqrt(effErrBinomial[1]);
  
  return effErrBinomial;
}


//_____________________________________________________________________________
TH1* GetHisto(TString histoName, TFile* file, TList* histoList)
{
  /// Get histogram
  TH1* histo = 0x0;
  if ( histoList )
    histo = (TH1*)histoList->FindObject(histoName.Data());
  else
    histo = (TH1*)file->FindObjectAny(histoName.Data());
  
  return histo;
}

//_____________________________________________________________________________
Int_t GetEffIndex ( Int_t iel, Int_t icount, Int_t ich = -1 )
{
  /// Get efficienct histogram index
  if ( iel == 0 ) return icount;
  return 3 + 4*3*(iel-1) + 3*ich + icount;
}

//_____________________________________________________________________________
TList* GetOCDBList ( )
{
  /// Get list of CDB objetcs
  TString storageType = AliCDBManager::Instance()->GetDefaultStorage()->GetType();
  Bool_t isGrid = storageType.Contains("alien");
  TString baseFolder = AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder();
  
  TList* outList = new TList();
  outList->SetOwner();
  TString dirName[3] = {"GlobalTriggerCrateConfig", "RegionalTriggerConfig", "LocalTriggerBoardMasks"};
  for ( Int_t idir=0; idir<3; idir++ ) {
    if ( isGrid ) {
      TGridResult *res = gGrid->Ls(Form("%s/MUON/Calib/%s", baseFolder.Data(), dirName[idir].Data()));
      if (!res) return 0x0;
      for ( Int_t ires=0; ires<res->GetEntries(); ires++ ) {
        TString currFile = static_cast<TMap*>(res->At(ires))->GetValue("name")->GetName();
        outList->Add(new TObjString(currFile));
      }
    }
    else {
      TString fileListStr = gSystem->GetFromPipe(Form("ls %s/MUON/Calib/%s", baseFolder.Data(), dirName[idir].Data()));
      TObjArray* fileList = fileListStr.Tokenize("\n");
      for ( Int_t ires=0; ires<fileList->GetEntries(); ires++ ) {
        TString currFile = fileList->At(ires)->GetName();
        outList->Add(new TObjString(currFile));
      }
      delete fileList;
    }
  }
  return outList;
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
void TrigEffTrending(TObjArray runNumArray, TObjArray fileNameArray, TList& outCanList, TList& outList)
{
  /// Get the efficiency vs. run number
  TString elementName[3] = { "Chamber", "RPC", "Board" };
  TString countTypeName[4] = { "allTracks", "bendPlane", "nonBendPlane", "bothPlanes" };

  TString filename = "", effName = "", effTitle = "";

  SetMyStyle();
  Double_t effValues[3][2*kNch];
  const Int_t kNgraphs = kNch*3*2+3;
  TObjArray effList(kNgraphs);
  const Int_t kNeffVsRun = kNgraphs+1;
  TObjArray effVsRunList(kNeffVsRun);
  
  effName = "totalEffEvolution";
  effTitle = "Multinomial probability of firing at least 3/4 chambers";
  TH1D* totalEff = new TH1D(effName.Data(), effTitle.Data(), 0, 0., 1.);
  effVsRunList.AddAt(totalEff, kNeffVsRun-1);

  TString runNumString = "";
  for ( Int_t irun=0; irun<runNumArray.GetEntries(); irun++ ) {
    runNumString = runNumArray.At(irun)->GetName();

    // Search corresponding file (for sorting)
    for ( Int_t ifile=0; ifile<fileNameArray.GetEntries(); ifile++ ) {
      filename = fileNameArray.At(ifile)->GetName();
      if ( filename.Contains(runNumString.Data()) ) break;
    }

    if ( filename.Contains("alien://") && ! gGrid ) gGrid->Connect("alien://");
    
    //
    // First get the list of efficiency graphs
    //
    
    // Chamber efficiency
    TFile* file = TFile::Open(filename.Data());
    if ( ! file ) {
      printf("Warning: cannot find %s\n", filename.Data());
      continue;
    }
    
    TList* trigEffList = (TList*)file->FindObjectAny("triggerChamberEff");
    if ( ! trigEffList ) printf("Warning: histo list not found in %s. Check directly in file\n", filename.Data());
    
    TH1* histoDen = 0x0;
    for ( Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts; icount++ ) {
      effName = countTypeName[icount] + "Count" + elementName[0];
      if ( icount == 0 ) {
        histoDen = GetHisto(effName, file, trigEffList);
        continue;
      }
      
      TH1* histoNum = GetHisto(effName, file, trigEffList);
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(histoNum, histoDen);
      effName.ReplaceAll("Count","Eff");
      graph->SetName(effName.Data());
      effList.AddAt(graph, GetEffIndex(0, icount-1));
    }
    file->Close();
  
    if ( ! histoDen ) {
      printf("Error: cannot find histograms in file %s. Skip to next\n", filename.Data());
      continue;
    }
  
    // RPC/board efficiency
    AliMUONTriggerChamberEfficiency trigChEff(filename);
    for ( Int_t iel=1; iel<3; iel++ ) {
      for ( Int_t ich=0; ich<kNch; ich++ ) {
        for ( Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts-1; icount++ ) {
          TObject* obj = trigChEff.GetEffObject(2-iel, icount, ich);
          effList.AddAt(obj, GetEffIndex(iel, icount, ich));
        }
      }
    }
    
    // Fill efficiency vs run
    for ( Int_t iel=0; iel<3; iel++ ) {
      for ( Int_t ich=0; ich<kNch; ich++ ) {
        for ( Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts-1; icount++ ) {
          TGraphAsymmErrors* graph = static_cast<TGraphAsymmErrors*>(effList.At(GetEffIndex(iel, icount, ich)));
          Int_t nPoints = ( iel == 0 ) ? 1 : graph->GetN();
          for ( Int_t ipoint=0; ipoint<nPoints; ipoint++ ) {
            Int_t currPoint = ( iel == 0 ) ? ich : ipoint;
            Double_t xpt, ypt;
            graph->GetPoint(currPoint, xpt, ypt);
            effValues[icount][ich] = ypt;
            Int_t ihisto = GetEffIndex(iel,icount,ich);
            TH2* effHisto = static_cast<TH2*>(effVsRunList.At(ihisto));
            if ( ! effHisto ) {
              effName = Form("effEvolution%s%s", countTypeName[icount+1].Data(), elementName[iel].Data());
              effTitle = Form("Trigger chamber efficiency vs run");
              if ( iel>0 ) {
                effName += Form("Ch%i", 11+ich);
                effTitle += Form(" for chamber %i", 11+ich);
              }
              effHisto = new TH2D(effName.Data(), effTitle.Data(), 0, 0., 1., graph->GetN(), xpt-0.5, xpt-0.5+(Double_t)graph->GetN());
              effVsRunList.AddAt(effHisto, ihisto);
            }
            Int_t currBin = effHisto->Fill(runNumString.Data(), xpt, ypt);
            Double_t err = 0.5*(graph->GetErrorYlow(ipoint) + graph->GetErrorYhigh(ipoint));
            Int_t binx, biny, binz;
            effHisto->GetBinXYZ(currBin, binx, biny, binz);
            effHisto->SetBinError(binx, biny, err);
            effValues[icount][ich+kNch] = err;
          } // loop on points
        } // loop on counts
      } // loop on chambers
      if ( iel > 0 ) continue;
      Double_t* binomialEff = GetBinomial(effValues[0], effValues[1], effValues[2]);
      Int_t currBin = totalEff->Fill(runNumString, binomialEff[0]);
      // CAVEAT!!!!
      // Well, error calculation of the binomial efficiency is a mess...
      // Sometimes it happens that efficiency + error > 1.
      // In that case reduce the error.
      totalEff->SetBinError(currBin, TMath::Min(binomialEff[1], 1.-binomialEff[0]));
      delete binomialEff;
    } // loop on detection elements
  } // loop on runs
  
  // Set correct range (do not show last empty bins)
  for ( Int_t ihisto=0; ihisto<effVsRunList.GetEntries(); ihisto++ ) {
    TH1* histo = static_cast<TH1*>(effVsRunList.At(ihisto));
    SetRunAxisRange(histo->GetXaxis());
    outList.Add(histo);
    //histo->GetXaxis()->SetLabelSize(0.03);
  }

  TString canName = "totalEff";
  TCanvas* can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
  TH1* totEff = (TH1*)effVsRunList.At(kNeffVsRun-1);
  totEff->GetYaxis()->SetRangeUser(0.9,1.05);
  totEff->GetYaxis()->SetTitle("Probability to satisfy trigger conditions (3/4)");
  totEff->SetStats(kFALSE);
  totEff->DrawCopy();
  outCanList.Add(can);

  Int_t color[3] = {kBlack, kRed, kBlue};
  Int_t markStyle[3] = {20, 24, 26};
  TLegend* leg = 0x0;

  for ( Int_t ich=0; ich<kNch; ich++ ) {
    canName = Form("trigEffCh%i", 11+ich);
    can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
    can->SetGridy();
    leg = new TLegend(0.6, 0.2, 0.9, 0.4);
    leg->SetBorderSize(1);
    //can->Divide(2,2);
    TString drawOpt = "e";
    for(Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts-1; icount++) {
      //can->cd(icount+1);
      TH2* histo = static_cast<TH2*>(effVsRunList.At(GetEffIndex(0, icount)));
      TH1* chEff = histo->ProjectionX(Form("effEvolutionCh%i",11+ich), ich+1, ich+1);
      chEff->SetTitle(Form("%s for chamber %i", histo->GetTitle(), 11+ich));
      chEff->GetYaxis()->SetRangeUser(0.9,1.);
      chEff->SetStats(kFALSE);
      chEff->GetYaxis()->SetTitle("Trigger chamber efficiency");
      TH1* copyEff = chEff->DrawCopy(drawOpt.Data());
      copyEff->SetLineColor(color[icount]);
      copyEff->SetMarkerColor(color[icount]);
      copyEff->SetMarkerStyle(markStyle[icount]);
      leg->AddEntry(copyEff, countTypeName[icount+1].Data(), "lp");
      drawOpt = "esame";
    } // loop on counts
    leg->Draw("same");
    outCanList.Add(can);
  } // loop on chambers
  
  for ( Int_t iel=1; iel<3; iel++ ) {
    for ( Int_t ich=0; ich<kNch; ich++ ) {
      Int_t icount = AliMUONTriggerEfficiencyCells::kBothPlanesEff; // Just plot the efficiency for both
//      for ( Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts-1; icount++ ) {
        canName = Form("trigEff%sCh%i", elementName[iel].Data(), 11+ich);
        can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
        can->SetRightMargin(0.14);
        TH2* histo = static_cast<TH2*>(effVsRunList.At(GetEffIndex(iel, icount,ich)));
        histo->SetStats(kFALSE);
        histo->GetYaxis()->SetTitle(elementName[iel].Data());
        histo->DrawCopy("COLZ");
//      } // loop on counts
      outCanList.Add(can);
    } // loop on chambers
  } // loop on detection element type
}

//_____________________________________________________________________________
void MaskTrending ( TObjArray runNumArray, TString defaultStorage, TList& outCanList, TList& outList )
{
  /// Get the masks vs. run number
  if ( defaultStorage.Contains("alien://") || defaultStorage.Contains("raw://") ) {
    if ( ! gGrid ) TGrid::Connect("alien://");
    if ( ! gGrid ) {
      printf("Error: Problem connetting to grid: nothing done");
      return;
    }
  }
  
  
  TObjArray maskedList(8);
  TObjArray auxList(8);
  auxList.SetOwner();
  TString histoName = "", histoTitle = "";
  for(Int_t icath=0; icath<2; icath++){
    TString cathName = ( icath==0 ) ? "bendPlane" : "nonBendPlane";
    for(Int_t ich=0; ich<kNch; ich++){
      histoName = Form("%sMaskCh%i", cathName.Data(), 11+ich);
      histoTitle = Form("Chamber %i - %s: fraction of masked channels", 11+ich, cathName.Data());
      TH2* histo = new TH2D(histoName.Data(), histoTitle.Data(),0,0.,1., 234, 0.5, 234. + 0.5);
      histo->GetYaxis()->SetTitle("Board Id");
      histo->SetOption("COLZ");
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
  
  if ( ! AliCDBManager::Instance()->IsDefaultStorageSet() ) AliCDBManager::Instance()->SetDefaultStorage(defaultStorage.Data());
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
      
      if ( ! ocdbFileList ) ocdbFileList = GetOCDBList();
  
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
    SetRunAxisRange(histo->GetXaxis());
    outList.Add(histo);
    
    canName = Form("%sCan", histo->GetName());
    TCanvas* can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
    can->SetRightMargin(0.14);
    histo->SetStats(kFALSE);
    histo->DrawCopy("COLZ");
    outCanList.Add(can);
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
void ScalerTrending ( TObjArray runNumArray, TString mergedFileName, TString defaultStorage, TList& outCanList, TList& outList )
{
  /// Get the scalers vs. run number
  if ( defaultStorage.Contains("alien://") || defaultStorage.Contains("raw://") ) {
    if ( ! gGrid ) TGrid::Connect("alien://");
    if ( ! gGrid ) {
      printf("Error: Problem connetting to grid: nothing done");
      return;
    }
  }
  
  if ( ! AliCDBManager::Instance()->IsDefaultStorageSet() ) AliCDBManager::Instance()->SetDefaultStorage(defaultStorage.Data());

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
  
  
  TString sHistName, sHistNameFull, sTitleName;
  Int_t nRuns = runNumArray.GetEntries();

  //
  //Fill histos for Scalers and QA
  //
  //loop on run list
  for ( Int_t iRun = 0; iRun < runNumArray.GetEntries(); iRun++ ) {
    
    TString sRunNr = ((TObjString*)runNumArray.At(iRun))->GetString();
    Int_t runNr = sRunNr.Atoi();
    AliAnalysisTriggerScalers triggerScaler(runNr);
    AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(triggerScaler.GetOCDBObject("GRP/CTP/Config",runNr));
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
          TGraph* graph = triggerScaler.PlotTrigger(currTrigName.Data(),sScaler[iScaler].Data());
        
          sHistName = Form("%s_%s",currTrigName.Data(),sScaler[iScaler].Data());
          sHistNameFull = Form("Scalers_%s",sHistName.Data());
          
          TH1* hist = (TH1*) hFromScalers.FindObject(sHistNameFull);
          if ( ! hist ) {
            hist = new TH1D(sHistNameFull,sHistName,nRuns,1.,1.+(Double_t)nRuns);
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
    outCanList.Add(canvas);
  }
  
  file->Close();
}

//_____________________________________________________________________________
void trigEffQA(TString fileListName, TString outFilename = "", TString defaultStorage = "raw://", Bool_t doScalers = kFALSE)
{
  /// Main function
  ifstream inFile(fileListName.Data());
  TObjArray fileNameArray, runNumArray;
  fileNameArray.SetOwner();
  runNumArray.SetOwner();
  TString currString = "";
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      currString.ReadLine(inFile,kTRUE); // Read line
      if ( currString.IsNull() || ! currString.Contains(".root") ||
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
  
  // Instead of using the efficiency stored in the QA output
  // search for the new efficiency produced with trigger tracks only
  TObjArray tmpArray = fileNameArray;
  TObjArray* finalFileNameArray =  ChangeFilenames(tmpArray) ? &tmpArray : &fileNameArray;
  
  TList outCanList, outList;
  TrigEffTrending(runNumArray, *finalFileNameArray, outCanList, outList);
  if ( ! defaultStorage.IsNull() ) MaskTrending(runNumArray, defaultStorage, outCanList, outList);
  if ( ! defaultStorage.IsNull() && doScalers ) ScalerTrending(runNumArray, "QAresults_Merged.root", defaultStorage, outCanList, outList);
  
  if ( outFilename.IsNull() ) return;
  
  TString outCanName = outFilename;
  outCanName.ReplaceAll(".root",".pdf");
  for ( Int_t ican=0; ican<outCanList.GetEntries(); ican++ ) {
    TString canName = outCanName;
    if ( ican == 0 ) canName.Append("("); // open pdf file
    else if ( ican == outCanList.GetEntries()-1 ) canName.Append(")"); // close pdf file
    static_cast<TCanvas*>(outCanList.At(ican))->Print(canName.Data());
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
