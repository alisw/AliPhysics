///
/// \file MakeDataMCComparisonPerSMClusterEbin.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Compare cluster shapes per SM or TCard channel
///
/// Compare different shape parameter histograms in data and MC productions
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
/// \param titleMC   : Simplified string acronym of input MC
/// \param titleData : Simplified string acronym of input data
/// \param debug     : Bool to activate printfs.
//------------------------------------------------------------------------------
void MakeDataMCComparisonPerSMClusterEbin
( 
 TString titleMC  = "JJDecLow", // "JJDecHigh", "MB", "JJDecLow"
 TString titleData= "LHC11cd_EMC7",//"LHC11cd_EMC7", "LHC11cd_INT7"
 Bool_t  debug    = kFALSE
)
{  
  // Energy bins and productions
  // Depending on the production more or less energy bins considered
  
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
  if ( titleMC.Contains("Low") )
  {
    // JJDecLow
    Double_t energyBins[] = { 8, 10, 12, 14, 16, 18, 20 };
    for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
    mcDir = "JJ_Dec_GJ";
    mcLeg = "#gammaJ+JJ(p^{EMCal}_{T,#gamma}>3.5 GeV/#it{c})";
  }
  
  if ( titleMC.Contains("High") )
  {
    // JJDecHigh
    Double_t energyBins[] = { 18, 20, 25, 30, 35, 40, 50 };
    for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
    mcDir = "JJ_Dec_High_GJ";
    mcLeg = "#gammaJ+JJ(p^{EMCal}_{T,#gamma}>7 GeV/#it{c})}";
  }
  
  if ( titleMC.Contains("MB") )
  {
    // MB  
    Double_t energyBins[] = { 5, 7, 12};
    for(Int_t ie = 0; ie < nEBins+1; ie++) binE[ie] = energyBins[ie];
    mcDir = "MB";
    mcLeg = "Min. Bias";
  }
  
  Int_t nprodtmp = 3;
  
  Bool_t bIndFracCase = kFALSE;
  if ( bIndFracCase )  nprodtmp = 7;
  
  const Int_t nProd = nprodtmp; 
  TString prod   [nProd];
  TString prodLeg[nProd]; 
  
  if ( debug ) printf("N prod %d\n",nProd);
  
  Int_t nproditer = 0;
  if ( !bIndFracCase )
  {
    prod   [nproditer++] = Form("data/module/TCardChannel3/%s",titleData.Data());
    prodLeg[nproditer-1] = "Data";
    prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic0_Scaled2_v2",mcDir.Data()); 
    prodLeg[nproditer-1] = "MC, default";
    prod   [nproditer++] =  Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_EcellCut_Scaled2_v2",mcDir.Data());
    prodLeg[nproditer-1] = "MC, mimic";//", #it{E}_{ind}>100 MeV";
    
//    prod   [nproditer++] = Form("data/module/TCardChannel2/ClV2/%s",titleData.Data());
//    prodLeg[nproditer-1] = "Data, Clus. V2";
//    prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic0_Scaled2/ClV2",mcDir.Data()) ;
//    prodLeg[nproditer-1] = "MC, default, Clus. V2";
//    prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_EcellCut_Scaled2/ClV2",mcDir.Data());
//    prodLeg[nproditer-1] = "MC, mimic, Clust. V2";//", #it{E}_{ind}>100 MeV, Clus. V2";

    //prod   [nproditer++] = Form("data/module/TCardChannel_NLoff/%s",titleData.Data());
    //prodLeg[nproditer-1] = "Data, NL off";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Default2" ,mcDir.Data());
    //prodLeg[nproditer-1] = "MC, default";

    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic0_NLoff",mcDir.Data()); 
    //prodLeg[nproditer-1] = "MC, default, NL off";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c",mcDir.Data());
    //prodLeg[nproditer-1] = "MC, mimic";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_v1",mcDir.Data()) 
    //prodLeg[nproditer-1] = "MC, mimic";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_EcellCut",mcDir.Data()) 
    //prodLeg[nproditer-1] = "MC, mimic, #it{E}_{ind} > 100 MeV";
 
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_MimicEmCorr_v1",mcDir.Data());
    //prodLeg[nproditer-1] = "MC, mimic, EMCal Corr.";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_NLoff",mcDir.Data());
    //prodLeg[nproditer-1] = "MC, mimic, NL off";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_NLoff_EcellCut",mcDir.Data()) ;
    //prodLeg[nproditer-1] = "MC, mimic, #it{E}_{ind}>100 MeV, NL off";
    //prod   [nproditer++] = Form("simu/module//pp_7TeV_%s/TCardChannel_Mimic10c_NLoff_10MeV"   ,mcDir.Data()) ;
    //prodLeg[nproditer-1] = "MC, mimic, #it{E}_{ind}>100 MeV, #it{E}_{leak}< 10 MeV, NL off";
  }
  else
  {
    printf("INDUCED Fraction cases!!!! %d\n",nproditer);
    if ( titleMC.Contains("Low") )
    {
      prod   [nproditer++] = "data/LHC11cd_EMC7";
      prodLeg[nproditer-1] = "Data";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/default/sumw2on";
      prodLeg[nproditer-1] = "MC, Default";
      //      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_79/sumw2on"; // 0.1%
      //      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.1%";
      //      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_78/sumw2on"; // 0.2%
      //      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.2%";
      //      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_71/sumw2on"; // 0.3%
      //      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.3%";
      //      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_77/sumw2on"; // 0.4%
      //      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.4%";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_70/sumw2on"; // 0.5%
      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.5%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_76/sumw2on"; // 0.6%    
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.6%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_75/sumw2on"; // 0.7%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.7%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_69/sumw2on"; // 0.8%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.8%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_74/sumw2on"; // 0.9%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=0.9%";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_68/sumw2on"; // 1.0%
      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.0%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_89/sumw2on"; // 1.05%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.05%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_73/sumw2on"; // 1.1%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.1%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_95/sumw2on"; // 1.125%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.125%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_85/sumw2on"; // 1.15%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.15%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_90/sumw2on"; // 1.175%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.175%";
      //      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_67/sumw2on"; // 1.2%
      //      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.2%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_91/sumw2on"; // 1.225%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.225%";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_86/sumw2on"; // 1.25%
      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.25%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_92/sumw2on"; // 1.275%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.275%";
      //            prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_66/sumw2on"; // 1.3%
      //            prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.3%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_93/sumw2on"; // 1.325%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.325%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_87/sumw2on"; // 1.35%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.35%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_94/sumw2on"; // 1.375%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.375%";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_65/sumw2on"; // 1.4%
      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.4%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_88/sumw2on"; // 1.45%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.45%";
      prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_64/sumw2on"; // 1.5%
      prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.5%";
      //          prod   [nproditer++] = "simu/mimic/pp_7TeV_JJ_Dec_GJ/TCard_C_72/sumw2on"; // 1.6%
      //          prodLeg[nproditer-1] = "MC, mimic #mu_{1}=1.6%";
    }
  }
  
  if ( debug ) 
  {
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      printf("iprod %d,\n \t file: %s,\n \t legend: %s\n",
             iprod,prod[iprod].Data(),prodLeg[iprod].Data());
    }
  } 

  TString histoName[] = 
  {
    "SMM02NoCut","SMM02","SMM20LowM02NoCut","SMM20LowM02","SMM20HighM02NoCut","SMM20HighM02","SMNCell",
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
  Int_t   firstMC   = 1;
  TString outputFileName = "";
  TString opt = "";
  if(bIndFracCase)
  {
    if ( clusterization =="" ) clusterization = "_V1_Ecell100_Eseed500";
    opt = "_IndFracCases";
  }
  
  for(Int_t ihisto = 0; ihisto < 26; ihisto++)
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


