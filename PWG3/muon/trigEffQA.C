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
#include "TLine.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TGrid.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerChamberEfficiency.h"
#endif

const Int_t kNch = 4;
const Double_t kZero = 1.e-7; // Avoid problems when comparing to 0.
void SetMyStyle();
Int_t GetRunNumber(TString);
//Double_t GetBinomial(Double_t*, Double_t* eff2 = 0x0, Double_t* effBoth = 0x0);
Double_t* GetBinomial(Double_t*, Double_t* effErr2 = 0x0, Double_t* effErrBoth = 0x0);
Double_t* GetProdErr(Double_t*, Int_t, Int_t nFactors=kNch);
Double_t* GetConditionalEffErr(Double_t*, Double_t*, Double_t*, Int_t exclude=-1);
TH1* GetHisto(TString, TFile*, TList*);

void trigEffQA(TString fileListName, TString outFilename = "")
{
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
    printf("Fatal: cannot open input file\n");
    return;
  }

  runNumArray.Sort();

  //TString elementName[2] = {"Board", "RPC"};
  TString countTypeName[3] = {"bendPlane", "nonBendPlane","bothPlanes"};

  TObjArray chEffList(13);

  TCanvas* can = 0x0;
  TString filename = "", histoName = "", histoTitle = "";

  histoName = "totalEffEvolution";
  histoTitle = "Multinomial probability of firing at least 3/4 chambers";
  TH1* totalEff = new TH1D(histoName.Data(), histoTitle.Data(), 0, 0., 1.);
  chEffList.AddAtAndExpand(totalEff, 12);

  SetMyStyle();
  Double_t effValues[3][2*kNch];

  TString runNumString = "";
  for ( Int_t irun=0; irun<runNumArray.GetEntries(); irun++ ) {
    runNumString = runNumArray.At(irun)->GetName();

    // Search corresponding file (for sorting)
    for ( Int_t ifile=0; ifile<fileNameArray.GetEntries(); ifile++ ) {
      filename = fileNameArray.At(ifile)->GetName();
      if ( filename.Contains(runNumString.Data()) ) break;
    }

    if ( filename.Contains("alien://") && ! gGrid )
      gGrid->Connect("alien://");

    TFile* file = TFile::Open(filename.Data());
    if ( ! file ) {
      printf("Warning: cannot find %s\n", filename.Data());
      continue;
    }
    TList* trigEffList = (TList*)file->FindObjectAny("triggerChamberEff");
    if ( ! trigEffList )
      printf("Warning: histo list not found in %s. Check directly in file\n", filename.Data());

    TH1* histoDen = GetHisto("allTracksCountChamber", file, trigEffList);
    if ( ! histoDen ) {
      printf("Error: cannot find histograms in file %s. Skip to next\n", filename.Data());
      continue;
    }
    for(Int_t icount=0; icount<AliMUONTriggerEfficiencyCells::kNcounts-1; icount++) {
      histoName = countTypeName[icount] + "CountChamber";
      TH1* histoNum = GetHisto(histoName, file, trigEffList);
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(histoNum, histoDen);
      Double_t xpt, ypt;
      for ( Int_t ich=0; ich<kNch; ich++ ) {
        graph->GetPoint(ich, xpt, ypt);
        effValues[icount][ich] = ypt;
        Int_t ihisto = kNch*icount+ich;
        TH1* chEff = (TH1*)chEffList.At(ihisto);
        if ( ! chEff ) {
          histoName = Form("effEvolution%sCh%i", countTypeName[icount].Data(), 11+ich);
          histoTitle = Form("Efficiency of trigger chamber %i", 11+ich);
          chEff = new TH1D(histoName.Data(), histoTitle.Data(), 0, 0., 1.);
          chEffList.AddAtAndExpand(chEff, ihisto);
        }
        Int_t currBin = chEff->Fill(runNumString.Data(), ypt);
        Double_t err = 0.5*(graph->GetErrorYlow(ich) + graph->GetErrorYhigh(ich));
        chEff->SetBinError(currBin, err);
        effValues[icount][ich+kNch] = err;
      } // loop on chambers
    } // loop on counts
    Double_t* binomialEff = GetBinomial(effValues[0], effValues[1], effValues[2]);
    Int_t currBin = totalEff->Fill(runNumString, binomialEff[0]);
    // CAVEAT!!!!
    // Well, error calculation of the binomial efficiency is a mess...
    // Sometimes it happens that efficiency + error > 1.
    // In that case reduce the error.
    totalEff->SetBinError(currBin, TMath::Min(binomialEff[1], 1.-binomialEff[0]));
    delete binomialEff;
  } // loop on runs

  TString outCanName = outFilename;
  outCanName.ReplaceAll(".root",".pdf");
  TString canName = "";
  canName = "totalEff";
  can = new TCanvas(canName.Data(), canName.Data(), 200, 10, 600, 600);
  if ( ! outCanName.IsNull() )
    can->Print(Form("%s[", outCanName.Data())); // Open pdf file

  // Set correct range (do not show last empty bins)
  for ( Int_t ihisto=0; ihisto<chEffList.GetEntries(); ihisto++ ) {
    TH1* histo = (TH1*)chEffList.At(ihisto);
    histo->GetXaxis()->SetRange(1,fileNameArray.GetEntries());
  }

  TH1* totEff = (TH1*)chEffList.At(12);
  totEff->GetYaxis()->SetRangeUser(0.,1.1);
  totEff->GetYaxis()->SetTitle("Probability to satisfy trigger conditions (3/4)");
  totEff->SetStats(kFALSE);
  totEff->DrawCopy();
  if ( ! outCanName.IsNull() )
    can->Print(outCanName.Data());

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
      TH1* chEff = (TH1*)chEffList.At(4*icount+ich);
      chEff->GetYaxis()->SetRangeUser(0.48,1.05);
      chEff->SetStats(kFALSE);
      chEff->GetYaxis()->SetTitle("Trigger chamber efficiency");
      TH1* copyEff = chEff->DrawCopy(drawOpt.Data());
      copyEff->SetLineColor(color[icount]);
      copyEff->SetMarkerColor(color[icount]);
      copyEff->SetMarkerStyle(markStyle[icount]);
      leg->AddEntry(copyEff, countTypeName[icount].Data(), "lp");
      drawOpt = "esame";
    } // loop on counts
    leg->Draw("same");
    if ( ! outCanName.IsNull() )
      can->Print(outCanName.Data());
  } // loop on chambers
  if ( ! outCanName.IsNull() ) {
    can->Print(Form("%s]", outCanName.Data())); // close pdf file
    // There is a bug when creating a pdf
    // So create a ps and then convert via epstopdf
    if ( outCanName.Contains(".ps") || outCanName.Contains(".eps") ) {
      gSystem->Exec(Form("epstopdf %s", outCanName.Data()));
      gSystem->Exec(Form("rm %s", outCanName.Data()));
    }
  }

  TFile* outFile = 0x0;
  if ( ! outFilename.IsNull() ) {
    outFile = new TFile(outFilename.Data(), "recreate");
    chEffList.Write();
    outFile->Close();
  }
}


//_____________________________________________________________________________
void SetMyStyle()
{
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
Int_t GetRunNumber(TString filePath)
{
  TObjArray* array = filePath.Tokenize("/");
  array->SetOwner();
  TString auxString = "";
  Int_t runNum = -1;
  for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
    auxString = array->At(ientry)->GetName();
    if ( auxString.Contains("000") ) {
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


// //_____________________________________________________________________________
// Double_t GetBinomial(Double_t* eff1, Double_t* eff2, Double_t* effBoth)
// {
//   Double_t binomialEff = 0.;
//   Double_t eff44 = 1.;
//   Double_t auxEff[4];
//   Double_t auxBinomial = 1.;

//   // All chamber efficient
//   for ( Int_t ich=0; ich<4; ich++ ) {
//     eff44 *= eff1[ich];
//   }
//   if ( eff2 ) {
//     for ( Int_t ich=0; ich<4; ich++ ) {
//       auxEff[ich] = ( eff1[ich] > kZero ) ? effBoth[ich]/eff1[ich] : 0.;
//     }
//     auxBinomial = GetBinomial(auxEff);
//     eff44 *= auxBinomial;
//   }

//   binomialEff += eff44;

//   // 1 chamber inefficient
//   for ( Int_t ich=0; ich<4; ich++ ) {
//     Double_t eff34 = 1.;
//     for ( Int_t jch=0; jch<4; jch++ ) {
//       eff34 *= ( ich == jch ) ? ( 1. - eff1[jch] ) : eff1[jch];
//     }
//     if ( eff2 ) {
//       for ( Int_t jch=0; jch<4; jch++ ) {
//         if ( ich == jch ) {
//           auxEff[jch] = ( eff1[jch] < 1. ) ? ( eff2[jch] - effBoth[jch] ) / ( 1. - eff1[jch]) : 0.;
//         }
//         else {
//           auxEff[jch] = ( eff1[ich] > kZero ) ? effBoth[ich]/eff1[ich] : 0.;
//         }
//         auxBinomial = GetBinomial(auxEff);
//         eff34 *= auxBinomial;
//       }
//     }

//     binomialEff += eff34;
//   }

//   return binomialEff;
// }


//_____________________________________________________________________________
Double_t* GetBinomial(Double_t* effErr1, Double_t* effErr2, Double_t* effErrBoth)
{
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
Double_t* GetProdErr(Double_t* effErr, Int_t exclude, Int_t nFactors)
{
  Double_t prod = 1.;
  Double_t relErr = 0., relProdErrSquare = 0.;
  for ( Int_t iprod=0; iprod<nFactors; iprod++ ) {
    if ( iprod == exclude ) continue;
    prod *= effErr[iprod];
    relErr = ( effErr[iprod] > kZero ) ? effErr[iprod+nFactors]/effErr[iprod] : 0.;
    relProdErrSquare += relErr*relErr;
    //printf("%f +- %f  ", effErr[iprod], effErr[iprod+nFactors]); // REMEMBER TO CUT
  }
  Double_t* prodErr = new Double_t[2];
  prodErr[0] = prod;
  prodErr[1] = prod*TMath::Sqrt(relProdErrSquare);
  //printf("-> %f %f\n", prodErr[0], prodErr[1]); // REMEMBER TO CUT
  return prodErr;
}


//_____________________________________________________________________________
Double_t* GetConditionalEffErr(Double_t* effErr1, Double_t* effErr2, Double_t* effErrBoth, Int_t exclude)
{
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
TH1* GetHisto(TString histoName, TFile* file, TList* histoList)
{
  TH1* histo = 0x0;
  if ( histoList )
    histo = (TH1*)histoList->FindObject(histoName.Data());
  else 
    histo = (TH1*)file->FindObjectAny(histoName.Data());
  
  return histo;
}
