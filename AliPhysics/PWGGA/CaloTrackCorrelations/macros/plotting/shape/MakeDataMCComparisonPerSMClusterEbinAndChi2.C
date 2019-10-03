///
/// \file MakeDataMCComparisonPerSMClusterEbinAndChi2.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Calculate Chi2 of shape parameter data vs different MC mimic 
///
/// Compare different shape parameter histograms in data and MC productions
/// depending on the SM per cluster energy bis.
/// Treatment of the output of the class AliAnaClusterShapeCorrelStudies
/// Input are TH3 histograms where x=energy, y= SM number and z is one cluster parameter
/// see the different options in the lines below.
/// Most of the job is done in CompareTH3DataAndMCProd.C and CalculateParamChi2MCxTalkDataPerSM.C,
/// Here we define the data/MC samples, where are the files and what energy ranges to be used
/// and call the corresponding methods.
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
#include "CalculateParamChi2MCxTalkDataPerSM.C"

#endif

//------------------------------------------------------------------------------
/// Main method
///
/// \param titleMC   : Simplified string acronym of input MC
/// \param titleData : Simplified string acronym of input data
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakeDataMCComparisonPerSMClusterEbinAndChi2
( 
 TString titleMC  = "JJDecLow", // "JJDecHigh",  "JJDecLow"
 TString titleData= "LHC11cd_EMC7",
 Bool_t  debug    = kFALSE
)
{  
  // Energy bins 
  const Int_t nEBins = 5;//5;
  TArrayD binE; binE.Set(nEBins+1);
  Double_t binEErr[] = { 1, 1, 1, 1, 1, 1, 1};

  TString daLeg = "pp@7 TeV, LHC11c+d EMC7";
  if(titleData == "LHC11cd_EMC7")
    daLeg = "pp@7 TeV, LHC11c+d EMC7";  
  if(titleData == "LHC11cd_INT7")
    daLeg = "pp@7 TeV, LHC11c+d INT7";
  
  TString mcDir = "JJ_Dec_GJ";
  TString mcLeg = "MC: #gammaJ+JJ({p^{EMCal}_{T,#gamma}>3.5 GeV/#it{c})";
  Int_t nProdParam = 29;
  TArrayD prodParam; 

  if ( titleMC.Contains("Low") )
  {
    // JJDecLow
    Double_t energyBins[] = { 8, 10, 12, 14, 16, 18, 20 };
    for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
    mcDir = "JJ_Dec_GJ";
    mcLeg = "#gammaJ+JJ(p^{EMCal}_{T,#gamma}>3.5 GeV/#it{c})";
    
    // Production per induced fraction bins
    Double_t prodParamArr[] = 
    {0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.050,1.100,1.125,
      1.150,1.175,1.200,1.225,1.250,1.275,1.300,1.325,1.350,1.375,1.400,1.450,1.50,1.60};
    nProdParam = 29; // make sure it is the same number as the hand made array below
    prodParam.Set(nProdParam-1);
    for(Int_t iprod = 0; iprod < nProdParam-1; iprod++) 
      prodParam[iprod] = prodParamArr[iprod];
  }
  else if ( titleMC.Contains("High") )
  {
    // JJDecHigh
    Double_t energyBins[] = { 18, 20, 25, 30, 35, 40, 50 };
    for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
    mcDir = "JJ_Dec_High_GJ";
    mcLeg = "#gammaJ+JJ(p^{EMCal}_{T,#gamma}>7 GeV/#it{c})}";
    // Production per induced fraction bins
    Double_t prodParamArr[] = 
    {0,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,
      0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
      1.000,1.050,1.100,1.150,1.200,1.250,1.300,1.350,1.400,1.450};
    nProdParam = 31; // make sure it is the same number as the hand made array below
    prodParam.Set(nProdParam-1);
    for(Int_t iprod = 0; iprod < nProdParam-1; iprod++) 
      prodParam[iprod] = prodParamArr[iprod];
  }
  
  const Int_t nProd = nProdParam; 
  TString prod   [nProd];
  TString prodLeg[nProd]; 
  
  if ( debug ) printf("N prod %d\n",nProd);
  
  // Here we have the files for each production done for a different fraction
  // of induced energy mu1
  Int_t nproditer = 0;
  if ( titleMC.Contains("Low") )
  {
    prod   [nproditer++] = "data/LHC11cd_EMC7";
    prodLeg[nproditer-1] = "Data";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/default/sumw2on";
    prodLeg[nproditer-1] = "MC, Default";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_79/sumw2on"; // 0.1%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.1%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_78/sumw2on"; // 0.2%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.2%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_71/sumw2on"; // 0.3%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.3%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_77/sumw2on"; // 0.4%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.4%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_70/sumw2on"; // 0.5%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.5%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_76/sumw2on"; // 0.6%    
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.6%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_75/sumw2on"; // 0.7%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.7%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_69/sumw2on"; // 0.8%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.8%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_74/sumw2on"; // 0.9%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.9%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_68/sumw2on"; // 1.0%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.0%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_89/sumw2on"; // 1.05%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.05%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_73/sumw2on"; // 1.1%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.1%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_95/sumw2on"; // 1.125%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.125%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_85/sumw2on"; // 1.15%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.15%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_90/sumw2on"; // 1.175%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.175%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_67/sumw2on"; // 1.2%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.2%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_91/sumw2on"; // 1.225%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.225%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_86/sumw2on"; // 1.25%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.25%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_92/sumw2on"; // 1.275%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.275%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_66/sumw2on"; // 1.3%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.3%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_93/sumw2on"; // 1.325%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.325%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_87/sumw2on"; // 1.35%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.35%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_94/sumw2on"; // 1.375%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.375%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_65/sumw2on"; // 1.4%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.4%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_88/sumw2on"; // 1.45%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.45%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_64/sumw2on"; // 1.5%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.5%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_72/sumw2on"; // 1.6%
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.6%";
  }
  else if ( titleMC.Contains("High") )
  {
    prod   [nproditer++] = "data/higherE/LHC11cd_EMC7";
    prodLeg[nproditer-1] = "Data";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.000PerCent.root" ;
    prodLeg[nproditer-1] = "MC, Default";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.050PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.05%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.100PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.1%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.150PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.15%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.200PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.2%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.250PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.25%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.300PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.3%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.350PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.35%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.400PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.4%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.450PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.45%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.500PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.5%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.550PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.55%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.600PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.6%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.650PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.65%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.700PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.7%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.750PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.75%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.800PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.8%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.850PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.85%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.900PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.9%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_0.950PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.95%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.000PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.0%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.050PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.05%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.100PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.1%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.150PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.15%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.200PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.2%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.250PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.25%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.300PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.3%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.350PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.35%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.400PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.4%";
    prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_High_GJ/param/ScaledMerged_V1_Ecell100_Eseed500_Cell_IndEFrac_1.450PerCent.root" ;
    prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.45%";
  }
    
  if ( debug ) 
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      printf("iprod %d,\n \t file: %s,\n \t legend: %s\n",
             iprod,prod[iprod].Data(),prodLeg[iprod].Data());
    }
  } 

  Int_t nHisto = 7;
  TString histoName[] = 
  {
    "SMM02NoCut","SMM02",
    "SMM20LowM02NoCut","SMM20LowM02",
    "SMM20HighM02NoCut","SMM20HighM02",
    "SMNCell"
  }; // This was done with an old version of the code with few TH3 defined, 
     // if exercise repited now, more histogranms could be used.
  
  // Different options selecing the file or histogram name
  TString clusterization = "_V1_Ecell100_Eseed500";
  TString tm        = "_TMDep";//"TMFix", 
  TString pid       = "_Neutral";//"_Electron", "_Hadron","Neutral"
  
  Bool_t  plotRatio = kFALSE;
  Bool_t  saveHisto = kTRUE;
  Int_t   firstSM   = 0; 
  Int_t   lastSM    = 9; 
  Int_t   firstMC   = 1;
  TString outputFileName = "";
  
  for(Int_t ihisto = 0; ihisto < nHisto; ihisto++)
  {
    printf("histogram %d %s\n",ihisto,histoName[ihisto].Data());
    CompareTH3DataAndMCProd
    (outputFileName, 
     histoName[ihisto], nProd, prod, prodLeg,
     clusterization,tm,pid,
     titleMC  , mcLeg,     
     titleData, daLeg,     
     firstSM  , lastSM, 
     binE     , firstMC, 
     plotRatio, saveHisto,
     "ForChi2", debug);
  
    if(saveHisto)
    {
      printf("inputFile %s\n",outputFileName.Data());
      CalculateParamChi2MCxTalkDataPerSM
      (histoName[ihisto], 
       outputFileName, 
       0, // dataProd index
       firstSM, lastSM, 
       binE, prodParam,
       1, // extra rebin
       debug); 
    }    
  } // histo loop

}


