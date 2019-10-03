/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliMTRParameterizedResponse.h"

// ROOT includes
#include <Riostream.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TSystem.h"

#include "AliLog.h"
#include "AliMergeableCollection.h"

/// \cond CLASSIMP
ClassImp(AliMTRParameterizedResponse) // Class implementation in ROOT context
/// \endcond

//________________________________________________________________________
AliMTRParameterizedResponse::AliMTRParameterizedResponse() :
  TObject(),
  fResponseBoard(0x0),
  fResponseEta(0x0),
  fEtaBinning(0x0),
  fPtBinning(0x0)
{
  /// Default Ctor.

}


//________________________________________________________________________
AliMTRParameterizedResponse::~AliMTRParameterizedResponse()
{
  //
  /// Destructor
  //
  delete fResponseBoard;
  delete fResponseEta;
  delete fEtaBinning;
  delete fPtBinning;
}

//________________________________________________________________________
AliMTRParameterizedResponse::AliMTRParameterizedResponse ( const AliMTRParameterizedResponse& obj ) :
  TObject(obj),
  fResponseBoard(( obj.fResponseBoard ) ? static_cast<TObjArray*>(obj.fResponseBoard->Clone() ) : 0x0),
  fResponseEta(( obj.fResponseEta ) ? static_cast<TObjArray*>(obj.fResponseEta->Clone() ) : 0x0),
  fEtaBinning(( obj.fEtaBinning ) ? static_cast<TAxis*>(obj.fEtaBinning->Clone() ) : 0x0),
  fPtBinning(( obj.fPtBinning ) ? static_cast<TAxis*>(obj.fPtBinning->Clone() ) : 0x0)
{
  /// Copy constructor
}


//________________________________________________________________________
AliMTRParameterizedResponse& AliMTRParameterizedResponse::operator = ( const AliMTRParameterizedResponse& obj)
{
  /// Assignment operator
  if ( this != &obj ) {
    TObject::operator=(obj);
    delete fResponseBoard;
    fResponseBoard = ( obj.fResponseBoard ) ? static_cast<TObjArray*>(obj.fResponseBoard->Clone() ) : 0x0;
    delete fResponseEta;
    fResponseEta = ( obj.fResponseEta ) ? static_cast<TObjArray*>(obj.fResponseEta->Clone() ) : 0x0;
    delete fEtaBinning;
    fEtaBinning = ( obj.fEtaBinning ) ? static_cast<TAxis*>(obj.fEtaBinning->Clone() ) : 0x0;
    delete fPtBinning;
    fPtBinning = ( obj.fPtBinning ) ? static_cast<TAxis*>(obj.fPtBinning->Clone() ) : 0x0;
  }
  return *this;
}

//________________________________________________________________________
void AliMTRParameterizedResponse::InitFunctionParams ( TF1* func, Bool_t fitLowPtIncrease, Int_t itype ) const
{
  /// Initialize function parameters for fit
  Double_t defParams[8] = {0.5,1.,0.3,1.,0.2,0.1,0.35,1.};
//  Double_t defParams[7] = {0.5,1.,0.3,100.,2.5,0.5,0.5};
  Double_t fixLowPtIncrease[] = {0.,0.,1.,-1.};

  for ( Int_t ipar=0; ipar<func->GetNpar(); ipar++ ) func->ReleaseParameter(ipar);

  func->SetParameters(defParams);
  if ( itype == kLptOverApt ) func->SetParameter(7,1.);
  if ( ! fitLowPtIncrease ) {
    for ( Int_t ipar=3; ipar<7; ipar++ ) func->FixParameter(ipar,fixLowPtIncrease[ipar-3]);
  }

  func->SetParLimits(7,0.,1.);
}


//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::FitResponses ( Bool_t buildDataAptOverAllFromGraph, Bool_t fitLowPtIncrease )
{
  /// Fit Response

  Double_t minPt = 0., maxPt = 100.;
  TF1* fitFuncData = new TF1("fitFuncData",this,&AliMTRParameterizedResponse::FitFunctionErf,minPt,maxPt,8);
  fitFuncData->SetNpx(1000);
  TF1* fitFuncMC = new TF1("fitFuncMC",this,&AliMTRParameterizedResponse::FitFunctionErf,minPt,maxPt,8);
  fitFuncMC->SetNpx(1000);
  TString parNames[8] = {"Norm","p_{T}^{cut}","#sigma","Norm_{low}","p_{T,low}^{cut}","#sigma_{low}","p_{T}^{limit}","maxEff"};

  TF1* fitFuncMCApt = new TF1("fitFuncApt",this,&AliMTRParameterizedResponse::FitFunctionErfApt,minPt,maxPt,3);
  fitFuncMCApt->SetNpx(1000);
  TString parNamesApt[3] = {"Norm","p_{T}^{cut}","#sigma"};

  for ( Int_t ipar=0; ipar<8; ipar++ ) {
    fitFuncData->SetParName(ipar,parNames[ipar].Data());
    fitFuncMC->SetParName(ipar,parNames[ipar].Data());
  }
  for ( Int_t ipar=0; ipar<3; ipar++ ) {
    fitFuncMCApt->SetParName(ipar,parNames[ipar].Data());
  }


  TString funcName = "fitResponse";
  TString fitOpt = "QNWR";

  TGraphAsymmErrors* graph = 0x0;
  for ( Int_t iresp=0; iresp<2; iresp++ ) {

    Int_t checkMCAptAll = kTRUE;
    TObjArray* mcAptAll = GetResponse(kAptOverAll,kTRUE,iresp,kFALSE);
//    TObjArray* mcAptAll = GetResponse(kHptOverLpt,kTRUE,iresp); // Just for tests
    TObjArray* dataAptAll = 0x0;
    if ( mcAptAll ) {
      dataAptAll = static_cast<TObjArray*>(mcAptAll->Clone());
      AddResponse(kAptOverAll,kFALSE,iresp,dataAptAll);
      for ( Int_t igraph=0; igraph<dataAptAll->GetEntriesFast(); igraph++ ) {
        graph = static_cast<TGraphAsymmErrors*>(dataAptAll->UncheckedAt(igraph));
        TString graphName = graph->GetName();
        graphName.ReplaceAll("MC","SimulatedData");
        graph->SetName(graphName.Data());
      }
    }
    else checkMCAptAll = kFALSE;

    TGraphAsymmErrors* graphData = 0x0, *graphMC = 0x0;
    Bool_t isMCAptAllDone = kFALSE;
    for ( Int_t itype=kLptOverApt; itype<=kHptOverLpt; itype++ ) {
      TObjArray* dataLptApt = GetResponse(itype,kFALSE,iresp,kFALSE);
      if ( ! dataLptApt ) continue;
      TObjArray* mcLptApt = GetResponse(itype,kTRUE,iresp,kFALSE);

      Double_t minFit = ( itype == kLptOverApt ) ? 0.6 : 1.;
      Double_t maxFit = ( itype == kLptOverApt ) ? 4. : 10.;

      for ( Int_t igraph=0; igraph<dataLptApt->GetEntriesFast(); igraph++ ) {
        graph = static_cast<TGraphAsymmErrors*>(dataLptApt->UncheckedAt(igraph));
        InitFunctionParams(fitFuncData, fitLowPtIncrease, itype);
        fitFuncData->SetParLimits(0,0.,0.5);
        graph->Fit(fitFuncData,fitOpt.Data(),"",minFit,maxFit);
        TF1* clonedFitFuncData = static_cast<TF1*>(fitFuncData->Clone(funcName.Data()));
        graph->GetListOfFunctions()->Add(clonedFitFuncData);
        graphData = graph;

        if ( ! mcLptApt ) continue;
        graph = static_cast<TGraphAsymmErrors*>(mcLptApt->UncheckedAt(igraph));
        InitFunctionParams(fitFuncMC, fitLowPtIncrease, itype);
        if ( mcAptAll && itype == kLptOverApt ) {
          fitFuncMC->FixParameter(0,fitFuncData->GetParameter(0));
          fitFuncMC->FixParameter(7,fitFuncData->GetParameter(7));
        }
        graph->Fit(fitFuncMC,fitOpt.Data(),"",minFit,maxFit);
        TF1* clonedFitFuncMC = static_cast<TF1*>(fitFuncMC->Clone(funcName.Data()));
        graph->GetListOfFunctions()->Add(clonedFitFuncMC);
        graphMC = graph;

        if ( ! checkMCAptAll ) continue;
        isMCAptAllDone = kTRUE;

        graph = static_cast<TGraphAsymmErrors*>(mcAptAll->UncheckedAt(igraph));
        fitFuncMCApt->SetParameters(0.5,0.5,0.2);
        graph->Fit(fitFuncMCApt,fitOpt.Data(),"",0.,4.);
        graph->GetListOfFunctions()->Add(fitFuncMCApt->Clone(funcName.Data()));

        // We need to build a "data response" for the Apt/All case.
        // We have 2 different ways to do it.

        // First way:
        // Re-weight the Apt/All MC graph using the ratio Data/MC from Lpt/Apt
        // This can be done in two ways:
        // a) directly from the Lpt/Apt graph.
        // b) scaling the pt_cut and sigma parameters of the fit results in Lpt/Apt
        // In both cases, one need to find the pt_cut value in Lpt/Apt and Apt/All
        // While the former is rather strightforward, the latter is complicated
        // since the Apt/All shape is weird at forward rapidities (acceptance effect)
        // In order to be consistent, we estimate the pt_cut directly from the graphs.
        // In particular, we define it as the point on the graph for which the efficiency
        // is the average between the maximum value and the minimum one.

        // So, first let's search fot the pt_cut
        Double_t xpt, ypt;
        Double_t ptCuts[2] = {1., 0.5};
        for ( Int_t igr=0; igr<2; igr++ ) {
          TGraphAsymmErrors* auxGr = ( igr == 0 ) ? graphMC : graph;
          Double_t minEff = 10., maxEff = -1.;
          Int_t iptMin = -1, iptMax = -1;
          // First get minimum and maximum
          for ( Int_t ipt=0; ipt<auxGr->GetN(); ipt++ ) {
            auxGr->GetPoint(ipt,xpt,ypt);
            if ( ypt == 0. ) continue;
            if ( ypt < minEff ) {
              minEff = ypt;
              iptMin = ipt;
            }
            if ( ypt > maxEff ) {
              maxEff = ypt;
              iptMax = ipt;
            }
          }
          // Then search the halfway pt
          Double_t halfEff = 0.5*(minEff+maxEff);
          for ( Int_t ipt=iptMin; ipt<=iptMax; ipt++ ) {
            auxGr->GetPoint(ipt,xpt,ypt);
            if ( ypt > halfEff ) {
              ptCuts[igr] = xpt;
              break;
            }
          }
          AliDebug(1, Form("%s half efficiency reached at %g",auxGr->GetName(),ptCuts[igr]));
          AliInfo(Form("%s half efficiency reached at %g",auxGr->GetName(),ptCuts[igr]));
        }


        graph = static_cast<TGraphAsymmErrors*>(dataAptAll->UncheckedAt(igraph));

        // We now build the ratio data/MC from Lpt/Apt
        TGraph graphRatio(graphData->GetN());
        Double_t xptMC, yptMC;
        if ( buildDataAptOverAllFromGraph ) {
          // We can do it directly using the graphs
          for ( Int_t ipt=0; ipt<graphRatio.GetN(); ipt++ ) {
            graphData->GetPoint(ipt,xpt,ypt);
            graphMC->GetPoint(ipt,xptMC,yptMC);
            Double_t val = ( yptMC > 0. ) ? ypt/yptMC : 1.;
            graphRatio.SetPoint(ipt,xpt,val);
          }
        }
        else {
          // Or we can do it using the fit functions
          Double_t scaleFactor = ( ptCuts[0] > 0 ) ? ptCuts[1]/ptCuts[0] : 1.;
          // We scale the pt_cut and sigma parameters
          // according to the ratio of the values found above
          for ( Int_t ifunc=0; ifunc<2; ifunc++ ) {
            TF1* func = ( ifunc == 0 ) ? fitFuncData : fitFuncMC;
            for ( Int_t ipar=1; ipar<=2; ipar++ ) {
              func->SetParameter(ipar,func->GetParameter(ipar)*scaleFactor);
            }
          }
        }

        // Finally we apply the ratio data/MC to the Apt/All MC
        // Here we use the ratio of the pt_cut values as a stretch factor
        Double_t stretchFactor = ( ptCuts[1] > 0 ) ? ptCuts[0]/ptCuts[1] : 1.;
        Double_t val = 1.;
        for ( Int_t ipt=0; ipt<graph->GetN(); ipt++ ) {
          graph->GetPoint(ipt,xptMC,yptMC);
          if ( buildDataAptOverAllFromGraph ) {
            graphRatio.GetPoint(ipt,xpt,ypt);
            Double_t stretchX = TMath::Min(stretchFactor*xptMC,fPtBinning->GetBinUpEdge(fPtBinning->GetNbins()));
            val = yptMC*graphRatio.Eval(stretchX);
          }
          else {
            ypt = fitFuncMC->Eval(xptMC);
            val = ( ypt == 0. ) ? 1. : fitFuncData->Eval(xptMC) / ypt;
            val *= yptMC;
          }
          if ( val < 0.) val = 0.;
          else if ( val > 1. ) val = 1.;
          graph->SetPoint(ipt,xptMC,val);
        }

        // Second way:
        // assume that the relative difference between the pt_cut and sigma parameters
        // fitted from Lpt/Apt is the same also for the Apt/All case
        for ( Int_t ipar=1; ipar<=2; ipar++ ) {
          Double_t parMC = clonedFitFuncMC->GetParameter(ipar);
          Double_t parRatio = ( parMC == 0. ) ? 1. : clonedFitFuncData->GetParameter(ipar)/parMC;
          fitFuncMCApt->SetParameter(ipar,fitFuncMCApt->GetParameter(ipar)*parRatio);
        }
        graph->GetListOfFunctions()->Add(fitFuncMCApt->Clone(funcName.Data()));
      } // loop on graphs

      if ( isMCAptAllDone ) checkMCAptAll = kFALSE;
    } // loop on Lpt/Apt and Hpt/Lpt
  } // loop on responses per board or per eta

  return kTRUE;
}

//________________________________________________________________________
Double_t AliMTRParameterizedResponse::FitFunctionErfApt ( Double_t* xVal, Double_t *par )
{
  /// Fit function

  Double_t xx = xVal[0];
  Double_t sqrtTwo = TMath::Sqrt(2.);
  Double_t yVal = par[0]*(1.+TMath::Erf((xx-par[1])/par[2]/sqrtTwo));

  return yVal;
}

//________________________________________________________________________
Double_t AliMTRParameterizedResponse::FitFunctionErf ( Double_t* xVal, Double_t *par )
{
  /// Fit function

  Double_t xx = xVal[0];
  Double_t currX = TMath::Max(xx,par[6]);
  Double_t sqrtTwo = TMath::Sqrt(2.);
  Double_t yVal = par[7]+par[0]*(TMath::Erf((currX-par[1])/par[2]/sqrtTwo)-1.);
  if ( xx < par[6] ) yVal += par[3]*(TMath::Erf((-xx-par[4])/par[5]/sqrtTwo) - TMath::Erf((-par[6]-par[4])/par[5]/sqrtTwo));
//  Double_t zz = xx-par[6];
//  if ( xx < par[6] ) yVal += par[3]*zz+par[4]*zz*zz+par[5]*zz*zz*zz;

  return yVal;
}


//________________________________________________________________________
TString AliMTRParameterizedResponse::GetStdName ( Int_t itype, Bool_t isMC, Bool_t perBoard ) const
{
  TString str = "";
  TString ratioName[3] = { "AptOverAll", "LptOverApt", "HptOverLpt" };
  TString dataType = isMC ? "MC" : "Data";
  TString varType = perBoard ? "Board" : "Eta";

  str = ratioName[itype] + dataType + "Per" + varType;

  return str;
}

////________________________________________________________________________
//TH1* AliMTRParameterizedResponse::GetHisto ( Int_t type, Int_t isMC, Bool_t perBoard, Double_t val )
//{
//  /// Get response histogram
//  TObjArray* respArr = ( iresp==0 ) ? fResponseBoard : fResponseEta;
//  if ( ! respArr ) return 0x0;
//  Int_t itype = type;
//  if ( isMC ) itype += 3;
//  TObjArray* arr = static_cast<TObjArray*>(respArr->UncheckedAt(itype));
//  if ( ! arr ) return 0x0;
//  Int_t ibin = 0;
//  if ( perBoard ) ibin = TMath::Nint(val);
//  else ibin = fEtaBinning->FindBin(val);
//
//  return static_cast<TH1*>(arr->UncheckedAt(ibin-1));
//}

//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::SetAptOverAllMC ( TH2* histoApt, TH2* histoAll )
{
  /// Set ratio Apt/All
  return SetRatio(histoApt,histoAll,kAptOverAll,kTRUE);
}


//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::SetHptOverLpt ( TH2* histoHpt, TH2* histoLpt,
                                                    TH2* histoMCHpt, TH2* histoMCLpt )
{
  /// Set ratio Hpt/Apt
  SetRatio(histoHpt,histoLpt,kHptOverLpt,kFALSE);
  if ( histoMCHpt && histoMCLpt ) SetRatio(histoMCHpt,histoMCLpt,kHptOverLpt,kTRUE);
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::SetLptOverApt ( TH2* histoLpt, TH2* histoApt,
                                                    TH2* histoMCLpt, TH2* histoMCApt )
{
  /// Set ratio Lpt/Apt
  SetRatio(histoLpt,histoApt,kLptOverApt,kFALSE);
  if ( histoMCLpt && histoMCApt ) SetRatio(histoMCLpt,histoMCApt,kLptOverApt,kTRUE);
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::SetRatio ( TH2* histoNum, TH2* histoDen, Int_t itype, Bool_t isMC, Int_t rebin )
{
  /// Set ratio Lpt/Apt

  Bool_t perBoard = ( histoNum->GetYaxis()->GetNbins() == 234 ) ;

  if ( perBoard && itype == kAptOverAll ) {
    AliWarning("Cannot build Apt/All per board. You can do it per eta bins instead");
    return kFALSE;
  }

//  TH2* histoRatio = static_cast<TH2*>(histoNum->Clone("tmpRatio"));
//  if ( histoRatio->GetSumw2N() == 0 ) histoRatio->Sumw2();
//  histoRatio->Divide(histoDen);

  Bool_t isCloned = kFALSE;

  if ( ! perBoard && rebin > 1 ) {
    histoNum = static_cast<TH2*>(histoNum->Clone("tmpNum2D"));
    histoDen = static_cast<TH2*>(histoDen->Clone("tmpDen2D"));
    isCloned = kTRUE;
    histoNum->RebinY(rebin);
    histoDen->RebinY(rebin);
  }

  TAxis* yAxis = histoNum->GetYaxis();

  if ( ! perBoard && ! fEtaBinning ) fEtaBinning = static_cast<TAxis*>(yAxis->Clone());
  if ( ! fPtBinning ) fPtBinning = static_cast<TAxis*>(histoNum->GetXaxis()->Clone());

  TString stdName = GetStdName(itype,isMC,perBoard);

  for ( Int_t ibin=1; ibin<=yAxis->GetNbins(); ibin++ ) {
    TH1* hNum = histoNum->ProjectionX("tmp_hNum",ibin,ibin,"e");
    TH1* hDen = histoDen->ProjectionX("tmp_hDen",ibin,ibin,"e");
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(hNum,hDen,"cpe0");
    TString title = Form("%g < %s < %g",yAxis->GetBinLowEdge(ibin),yAxis->GetTitle(),yAxis->GetBinUpEdge(ibin));
    graph->SetName(Form("graph%s_%i",stdName.Data(),ibin));
    graph->SetTitle(title.Data());
//    TH1F* hCopy = new TH1F();
//    hNum->Copy(*hCopy);
//    hCopy->SetName(Form("histo%s_%i",stdName.Data(),ibin));
//    hCopy->SetTitle();
//    hCopy->SetDirectory(0);
//    hCopy->SetStats(0);
//    gr->SetHistogram(hCopy);
    delete hNum;
    delete hDen;

    // Smooth response:
    Double_t xpt, xpt1, xpt2, ypt, ypt1, ypt2, yptMod;
    Bool_t posDer = kFALSE;
    for ( Int_t ipt=1; ipt<graph->GetN(); ipt++ ) {
      graph->GetPoint(ipt,xpt,ypt);
      graph->GetPoint(ipt-1,xpt1,ypt1);
      Double_t der = (ypt-ypt1)/(xpt-xpt1);
      if ( der > 0. ) posDer = kTRUE;
      if ( posDer && ypt == 0. ) {
        yptMod = ypt1;
        for ( Int_t jpt=ipt+1; jpt<graph->GetN(); jpt++ ) {
          graph->GetPoint(jpt,xpt2,ypt2);
          if ( ypt2>0. ) {
            yptMod += (xpt - xpt1) * ( ypt2 - ypt1 ) / ( xpt2 - xpt1 );
            break;
          }
        }
        graph->SetPoint(ipt,xpt,yptMod);
        graph->SetPointEYlow(ipt,yptMod);
        graph->SetPointEYhigh(ipt,1.- yptMod);
        AliWarning(Form("%s: change point %i at pt = %g from %g to %g",graph->GetName(),ipt,xpt,ypt,yptMod));
      }
    }

//    TH1* projection = histoRatio->ProjectionX(Form("histo%s%sPer%s_%i",dataType.Data(),ratioName[type].Data(),varType.Data(),ibin),ibin,ibin,"e");
//    projection->SetDirectory(0);
    if ( ibin == 1 ) AddResponse(itype,isMC,perBoard,new TObjArray(yAxis->GetNbins()));
    GetResponse(itype,isMC,perBoard)->AddAt(graph,ibin-1);
  }

  if ( isCloned ) {
    delete histoNum;
    delete histoDen;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::SetFromMTRResponseTaskOutput ( Bool_t perBoard,
                                                                  const char* filenameData,
                                                                  const char* filenameMC,
                                                                  const char* filenameMCApt,
                                                                  Int_t rebin )
//const char* outputNameData, const char* outputNameMC, const char* identifierData, const char* identifierMC)
{
  /// Set histos from MTR response task output
  /// Filenames expected in the form: filename:outputName
  /// If no outputname is found, the default will be used.

  TString filenames[3] = {filenameData,filenameMC,filenameMCApt};
//  TString outNames[3] = {outputNameData,outputNameMC};
//  TString identifiers[2] = {identifierData,identifierMC};

  if ( perBoard && ! filenames[2].IsNull() ) AliWarning("Cannot build Apt/AllPt per local board (no board information on denominator). Please do it vs. eta instead");


  TString matchNames[4] = {"MatchNo", "MatchAllPt", "MatchLowPt", "MatchHighPt"};
  TString varType = perBoard ? "PerBoard" : "PerEta";

  TH2* histos[2] = {0x0,0x0};

  for ( Int_t ifile=0; ifile<3; ifile++ ) {
    if ( ifile > 0 && filenames[ifile].IsNull() ) break;
    gSystem->ExpandPathName(filenames[ifile]);
    TObjArray* arr = filenames[ifile].Tokenize(":");
    TString filename = arr->UncheckedAt(0)->GetName();
    TString outName = ( arr->GetEntriesFast() == 2 ) ? arr->UncheckedAt(1)->GetName() : "MTRResponseOut";
    delete arr;
    TFile* file = TFile::Open(filename.Data());
    if ( ! file ) return kFALSE;

    AliMergeableCollection* mcol = static_cast<AliMergeableCollection*>(file->FindObjectAny(outName.Data()));
    if ( ! mcol ) {
      AliError(Form("Cannot find %s in %s",outName.Data(),filename.Data()));
      file->ls();
      return kFALSE;
    }

    Bool_t isMC = ( ifile > 0 );
    Int_t firstType = ( ifile == 2 ) ? kAptOverAll : kLptOverApt;
    Int_t lastType = ( ifile == 2 ) ? kAptOverAll : kHptOverLpt;
    for ( Int_t itype=firstType; itype<=lastType; itype++ ) {
      TString identifier = ( itype == kHptOverLpt && ifile == 0 ) ? "MuTrig" : "MB";
      Int_t ihisto = 0;
      for ( Int_t imatch=itype; imatch<=itype+1; imatch++ ) {
        TString histoName = Form("histo%s%s",matchNames[imatch].Data(),varType.Data());
        AliDebug(1,Form("Getting /%s/%s from file %s\n",identifier.Data(),histoName.Data(),filename.Data()));
        histos[ihisto++] = static_cast<TH2*>(mcol->GetObject(identifier,histoName));
      }
      if ( histos[0] && histos[1] ) SetRatio(histos[1],histos[0],itype,isMC,rebin);
    }

//    delete mcol;
    file->Close();
  }
  return kTRUE;
}

//________________________________________________________________________
void AliMTRParameterizedResponse::ZoomPad()
{
  if ( gPad->GetEvent() != kButton1Double ) return;
  TVirtualPad* pad = gPad;
  Int_t px = pad->GetEventX();
  Int_t py = pad->GetEventY();
  TCanvas* can = new TCanvas("zoom","zoom",px,py,600,600);
  for ( Int_t iobj=0; iobj<pad->GetListOfPrimitives()->GetEntries(); iobj++ ) {
    TObject* obj = pad->GetListOfPrimitives()->At(iobj);
    obj = obj->Clone(Form("%s_zoom",obj->GetName()));
    TString drawOpt = obj->GetOption();
    if ( drawOpt.IsNull() ) {
      if ( obj->InheritsFrom(TGraph::Class()) ) {
        drawOpt = "p";
        if ( iobj == 1 ) drawOpt.Append("a");
      }
      else if ( obj->InheritsFrom(TH1::Class()) ) {
        drawOpt = "e";
        if ( iobj == 1 ) drawOpt.Append("same");
      }
    }
    obj->Draw(drawOpt.Data());
  }
  can->Modified();
  can->Update();
}

//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::CompareResponses ( Int_t itype, Bool_t perBoard ) const
{
  /// Compare data and MC responses
  TCanvas* can = 0x0;
  TGraphAsymmErrors* graph = 0x0;
  for ( Int_t imc=0; imc<2; imc++ ) {
    TObjArray* resp = GetResponse(itype,imc,perBoard);
//    if ( itype == kAptOverAll && imc == 1 ) resp = GetResponse(kHptOverLpt,0,perBoard); // Just for tests
    if ( ! resp ) return kFALSE;

    if ( imc == 0 ) {
      TString canName = GetStdName(itype,imc,perBoard);
      canName.Prepend("compareResp_");
      can = new TCanvas(canName.Data(),canName.Data(),0,0,1200,800);
      can->DivideSquare(resp->GetEntriesFast());
    }

    TIter next(resp);
    Int_t ipad = 1;
    while ( (graph=static_cast<TGraphAsymmErrors*>(next())) ) {
      can->cd(ipad++);
      if ( gPad->GetListOfExecs()->GetEntries() == 0 ) gPad->AddExec("ZoomPad","AliMTRParameterizedResponse::ZoomPad()");
      TGraphAsymmErrors* clonedGraph = static_cast<TGraphAsymmErrors*>(graph->Clone());
      Int_t icolor = imc+1;
      clonedGraph->SetLineColor(icolor);
      clonedGraph->SetMarkerColor(icolor);
      TF1* func = static_cast<TF1*>(clonedGraph->GetListOfFunctions()->At(0));
      if ( func ) func->SetLineColor(icolor);
      clonedGraph->GetYaxis()->SetRangeUser(0.,1.02);
      clonedGraph->GetXaxis()->SetRangeUser(0.,itype==kHptOverLpt ?12.:4.);
      clonedGraph->Draw(imc==0?"ap":"p");
      if ( func ) {
        TPaveText* pt = new TPaveText(0.45,0.12+0.3*imc,0.88,0.12+0.28+0.3*imc,"nbNDC");
        pt->SetBorderSize(1);
        pt->SetTextColor(icolor);
        pt->SetFillColor(0);
        pt->AddText(imc==0?"Data":"MC");
        for ( Int_t selpar=0; selpar<4; selpar++ ) {
          Int_t ipar = ( selpar < 3 ) ? selpar : 7;
          pt->AddText(Form("%s: %g #pm %g",func->GetParName(ipar),func->GetParameter(ipar),func->GetParError(ipar)));
        }
        pt->Draw();
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::ShowResponses ( Int_t itype, Bool_t isMC, Bool_t perBoard ) const
{
  /// Show responses
  TObjArray* resp = GetResponse(itype,isMC,perBoard);
  if ( ! resp ) return kFALSE;

  TString canName = GetStdName(itype,isMC,perBoard);
  canName.Prepend("showResp_");
  TCanvas* can = new TCanvas(canName.Data(),canName.Data(),0,0,1200,800);
  can->DivideSquare(resp->GetEntriesFast());

  TIter next(resp);
  TGraphAsymmErrors* graph = 0x0;
  Int_t ipad = 1;
  while ( (graph=static_cast<TGraphAsymmErrors*>(next())) ) {
    can->cd(ipad++);
    TGraphAsymmErrors* clonedGraph = static_cast<TGraphAsymmErrors*>(graph->Clone());
    clonedGraph->GetYaxis()->SetRangeUser(0.,1.02);
    clonedGraph->GetXaxis()->SetRangeUser(0.,itype==kHptOverLpt ?12.:4.);
    clonedGraph->Draw("ap");
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRParameterizedResponse::AddResponse ( Int_t itype, Bool_t isMC, Bool_t perBoard, TObjArray* resp )
{
  /// Add response array


  TObjArray* respArr = perBoard ? fResponseBoard : fResponseEta;

  if ( ! respArr ) {
    respArr = new TObjArray(6);
    respArr->SetOwner();
    if ( perBoard ) fResponseBoard = respArr;
    else fResponseEta = respArr;
  }

  Int_t idx = itype;
  if ( isMC ) idx += 3;

  TObjArray* arr = static_cast<TObjArray*>(respArr->UncheckedAt(idx));
  if ( arr ) {
    AliWarning(Form("Replacing existing response array for %s",GetStdName(itype,isMC,perBoard).Data()));
    delete arr;
  }

  resp->SetOwner();
  respArr->AddAt(resp,idx);

  return kTRUE;
}

//________________________________________________________________________
TObjArray* AliMTRParameterizedResponse::GetResponse ( Int_t itype, Bool_t isMC, Bool_t perBoard, Bool_t warn ) const
{
  /// Get response array
  TObjArray* respArr = perBoard ? fResponseBoard : fResponseEta;

  TString stdName = GetStdName(itype,isMC,perBoard);

  if ( ! respArr ) {
    if ( warn ) AliError(Form("Weights per %s are not set",perBoard?"board":"eta"));
    return 0x0;
  }

  Int_t idx = itype;
  if ( isMC ) idx += 3;

  TObjArray* arr = static_cast<TObjArray*>(respArr->UncheckedAt(idx));
  if ( ! arr ) {
    if ( warn ) AliError(Form("Weights are not set for %s",stdName.Data()));
    return 0x0;
  }
  return arr;
}

//________________________________________________________________________
Double_t AliMTRParameterizedResponse::GetWeight ( Double_t pt, Int_t ibin, Int_t itype, Bool_t isMC, Bool_t useFit, Bool_t perBoard ) const
{
  /// Get weight

  TObjArray* arr = GetResponse(itype,isMC,perBoard);
  if ( ! arr ) return 0.;

  TGraphAsymmErrors* graph = static_cast<TGraphAsymmErrors*>(arr->UncheckedAt(ibin));

  if ( useFit ) {
    TF1* fitFunc = static_cast<TF1*>(graph->GetListOfFunctions()->At(0));
    if ( fitFunc ) return fitFunc->Eval(pt);
    else AliWarning("Fit function not found for %s. Use histogram instead");
  }

  Int_t ipoint = fPtBinning->FindBin(pt) - 1;
  if ( ipoint < 0 ) {
    AliWarning("Underflow: set first bin instead");
    ipoint = 0;
  }
  else if ( ipoint >= graph->GetN() ) {
    AliWarning("Overflow: set last bin instead");
    ipoint = graph->GetN() - 1;
  }

  Double_t  xpt, ypt;
  graph->GetPoint(ipoint,xpt,ypt);

  AliDebug(3,Form("GraphName: %s pt: %g  xpt %g  ypt %g",graph->GetName(), pt, xpt, ypt));

  return ypt;
}

//________________________________________________________________________
Double_t AliMTRParameterizedResponse::WeightPerBoard ( Double_t pt, Int_t iboard, Int_t itype, Bool_t isMC, Bool_t useFit ) const
{
  /// Get weight per board

  if ( iboard == 0 || iboard > 234 ) return 0.;
  return GetWeight(pt,iboard-1,itype,isMC,useFit,kTRUE);
}


//________________________________________________________________________
Double_t AliMTRParameterizedResponse::WeightPerEta ( Double_t pt, Double_t eta, Int_t itype, Bool_t isMC, Bool_t useFit ) const
{
  /// Get weight per eta
  if ( ! fEtaBinning ) {
    AliError("Eta binning not initialised");
    return 0.;
  }

  Int_t ibin = fEtaBinning->FindBin(eta);
  if ( ibin == 0 || ibin > fEtaBinning->GetNbins() ) return 0.;
  return GetWeight(pt,ibin-1,itype,isMC,useFit,kFALSE);
}
