///
/// \file MakeNDataNMCComparisonPerSMClusterEbin.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Compare cluster shapes per SM or TCard channel
///
/// Compare different shape parameter histograms in several data and/or seveal MC productions
/// depending on the SM or T-Card channel number, the sum of all or pre-defined groups and
/// per cluster energy bis
/// Treatment of the output of the class AliAnaClusterShapeCorrelStudies
/// Input are TH3 histograms where x=energy, y= SM number or T-Card channel and z is one cluster parameter
/// see the different options in the lines below.
/// Most of the job is done in CompareTH3DataAndMCProd.C, here we define the data/MC samples, 
/// where are the files and what energy ranges to be used.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TGraphErrors.h>

#include "CompareTH3DataAndMCProd.C"

#endif

//------------------------------------------------------------------------------
/// Main method
///
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakeNDataNMCComparisonPerSMClusterEbin
( 
 Bool_t  debug    = kFALSE
)
{  
  // Energy bins and productions
  // Depending on the production more or less energy bins considered
  
  const Int_t nEBins = 5;
  TArrayD binE; binE.Set(nEBins+1);  
  Double_t energyBins[] = { 8, 10, 12, 14, 16, 18, 20 };
  for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
  
  const Int_t nProd = 2; 
  TString prod   [nProd];
  TString prodLeg[nProd]; 
  
  if ( debug ) printf("N prod %d\n",nProd);
  
  // Put here the path to the output files per Data or MC production
  
  Int_t nproditer = 0;
  prod   [nproditer++] = Form("data/module/TCardChannel2/LHC11cd_EMC7");
  prodLeg[nproditer-1] = Form("Data, pp @ 7 TeV, kEMC7, LHC11c+d");

  prod   [nproditer++] = Form("data/LHC16l_EG1");
  prodLeg[nproditer-1] = Form("Data, pp @ 13 TeV, kEMCEGA-EG1, LHC16l");
  
  TString titleData = "pp_7_13_TeV";
  TString daLeg     = "pp @ 7 & 13 TeV";
  
  TString titleMC   = "";
  TString mcLeg     = "";
  
  //
  Int_t   firstMC   = 2; // All data
  //
  
  if ( debug ) 
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      printf("iprod %d,\n \t file: %s,\n \t legend: %s\n",
             iprod,prod[iprod].Data(),prodLeg[iprod].Data());
    }
  } 

  Int_t nHisto = 26;
  TString histoName[] = 
  {
    "SMM02NoCut","SMM02","SMNCell",
    "SMM20LowM02NoCut","SMM20LowM02","SMM20HighM02NoCut","SMM20HighM02",
    "SMEMaxEClusterRat","SMNCellModuleMax","SMNCellModuleOut",
    "SMECellModuleMax","SMECellModuleOut","SMNCellModuleMaxOutRat","SMECellModuleMaxRat",
    "SMECellModuleMaxOutRat","SMECellModuleMaxTot","SMECellModuleMaxTotRat","SMECellModuleMaxTotRatClus",
    "SMNCellModuleOutModDiff","SMNCellModuleOutModSame","SMECellModuleOutModDiff","SMECellModuleOutModSame",
    "TCardChannelNCellModMax","TCardChannelNCell","TCardChannelM02","TCardChannelM02NoCut"
  };
  
  TString clusterization = "";//"_V1_Ecell100_Eseed500",
  TString tm        = "_TMDep";//"TMFix", 
  TString pid       = "_Neutral";//"_Electron", "_Hadron","Neutral"
  Bool_t  plotRatio = kFALSE;
  Bool_t  saveHisto = kFALSE;
  Int_t   firstP    = 0;
  Int_t   lastP     = 9;
  TString outputFileName = "";
  TString opt = "";
  
  for(Int_t ihisto = 0; ihisto < 3; ihisto++)
  {
    if(histoName[ihisto].Contains("TCard")) lastP = 15;
    printf("histogram %d %s\n",ihisto,histoName[ihisto].Data());
    
    CompareTH3DataAndMCProd
    (outputFileName, 
     histoName[ihisto], nProd, prod, prodLeg,
     clusterization,tm,pid,
     titleMC  , mcLeg,     
     titleData, daLeg,     
     firstP, lastP, 
     binE, firstMC, 
     plotRatio, saveHisto, 
     opt, debug);
  } // histogram loop
  
}


