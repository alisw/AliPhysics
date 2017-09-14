#ifndef AliUtilTOFParams_h
#define AliUtilTOFParams_h

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Set of parameters and utilities for the Pi/Ka/Pr analysis with TOF    //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <TROOT.h>
#include "TMath.h"
#endif

namespace AliUtilTOFParams {
  
  enum {nMultBin = 12, kPtBins = 59, kSpecies = 3, kFModes = 4};
  //Indexes for PID
  enum {kExpSpecies = 6, ke = 0, kmu, kpi, kK, kp, kd};
  enum {kpos = 0, kneg, kCharges = 2};
  //Indexes to store event
  //Track info
  enum fTrkMaskIndex {kNegTrk, kIsMismatch, kT0_0, kT0_1, kT0_2, kIsTOFout, kIsTOFTime, kIsTRDout, kPassGoldenChi2, kLimitfTrkMask};//Track information bitmask fTrkMask - kT0_0 (T0 TOF) kT0_1 (T0 T0A) kT0_2 (T0C)
  //Track cuts
  enum fTrkCutMaskIndex {kTPCSetL, kTPCSetT, kTPCChi2SetL, kTPCChi2SetT, kDCAzSetL, kDCAzSetT, kPrimSetL, kPrimSetL01, kPrimSetT01, kPrimSetT, kGeoCutSet1, kGeoCutSet2, kLimitfTrkCutMask};//Track cut information bitmask fTrkCutMask
  //Track PID
  enum fTPCPIDMaskIndex {kIsTPCElectron, kIsTPCMuon, kIsTPCPion, kIsTPCKaon, kIsTPCProton, kIsTPCDeuteron, kLimitfTPCPIDMask};//TPC PID information bitmask fTPCPIDMask
  
  //DCA binning
  const Int_t kDCAXYBins = 65536 - 3;
  const Int_t kDCAZBins = 65536 - 3;
  //DCA range for binning but also histograms
  const Double_t fDCAXYRange = 3.;
  const Double_t fDCAZRange = 3.;
  //DCA binning in histograms 
  const Int_t fDCAXYbins = 2000;
  
  //Dimensions of the screen
  const Double_t screendim[2] = {1366, 768};
  
  
  //Useful keywords
  const Int_t nMultEstimators = 2;
  const TString multlabel[nMultEstimators] = {"N_{Ch}", "%V0M"};
  const TString multest[nMultEstimators] = {"RefMult", "V0Mult"};
  const TString mEstimator[nMultEstimators] = {"Reference Multiplicity", "V0M Percentile"};
  const Int_t nCharges = 2;
  const TString pC[nCharges] = {"Pos", "Neg"};
  const TString pc[nCharges] = {"pos", "neg"};
  const TString pCharge[nCharges] = {"Positive", "Negative"};
  const TString pSign[nCharges] = {"+", "-"};
  const Int_t nSpecies = 4;
  const TString pS[nSpecies] = {"Pi", "K", "P", "D"};
  const TString ps[nSpecies] = {"pi", "k", "p", "d"};
  const TString pSpecies[nSpecies] = {"Pion", "Kaon", "Proton", "Deuteron"};
  const TString pspecies[nSpecies] = {"pion", "kaon", "proton", "deuteron"};
  const TString speciesRoot[nSpecies*2] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}", "D", "#bar{D}"};
  const TString speciesRootNoSign[nSpecies] = {"#pi", "K", "p", "D"};
  const TString pS_all[kExpSpecies] = {"El", "Mu", "Pi", "K", "P", "D"};
  const TString ps_all[kExpSpecies] = {"el", "mu", "pi", "k", "p", "d"};
  const TString pSpecies_all[kExpSpecies] = {"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron"};
  const TString pspecies_all[kExpSpecies] = {"electron", "muon", "pion", "kaon", "proton", "deuteron"};
  const TString speciesRoot_all[kExpSpecies*2] = {"e^{+}", "e^{-}", "#mu^{+}", "#bar{#mu}", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}", "D", "#bar{D}"};
  const TString speciesRootNoSign_all[kExpSpecies] = {"e", "#mu", "#pi", "K", "p", "d"};
  
  const TString gevoverc = "GeV/#it{c}";
  const TString ptstringOnly = "#it{p}_{T}";
  const TString ptstring = ptstringOnly + " (" + gevoverc + ")";
  const TString pstringOnly = "#it{p}";
  const TString pstring = pstringOnly + " (" + gevoverc + ")";
  const TString etastring = "#eta";
  const TString phistring = "#phi";
  const TString spectrastring = "1/N_{ev} d^{2}N/dp_{T}dy (Gev/c)^{-1}";
  const TString spectrastringdNdeta = "1/N_{ev} d^{2}N/dp_{T}d#eta (Gev/c)^{-1}";
  const TString spectrastring_2pipt = "1/N_{ev} 1/(2#pi p_{T}) d^{2}N/dp_{T}dy (Gev/c)^{-2}";
  const TString spectrastringdNdeta_2pipt = "1/N_{ev} 1/(2#pi p_{T}) d^{2}N/dp_{T}d#eta (Gev/c)^{-2}";
  const TString nsigmastring = "(T-T_{0}-T_{exp})/#sigma";
  const TString nsigmastringSpecies[kExpSpecies] = {"(T-T_{0}-T_{exp e})/#sigma", "(T-T_{0}-T_{exp #mu})/#sigma", "(T-T_{0}-T_{exp #pi})/#sigma", "(T-T_{0}-T_{exp K})/#sigma", "(T-T_{0}-T_{exp p})/#sigma", "(T-T_{0}-T_{exp d})/#sigma"};
  const TString tofsignalstring = "T-T_{0}-T_{exp} (ps)";
  const TString tofsignalstringSpecies[kExpSpecies] = {"(T-T_{0}-T_{exp e})", "(T-T_{0}-T_{exp #mu})", "(T-T_{0}-T_{exp #pi})", "(T-T_{0}-T_{exp K})", "(T-T_{0}-T_{exp p})", "(T-T_{0}-T_{exp d})"};
  const TString collenergystring = "#sqrt{#it{s}}";
  const TString HIcollenergystring = "#sqrt{#it{s}_{NN}}";
  const TString fitmodes[kFModes] = {"TFF", "RooFit", "CD", "Functions"};
  const TString DCAxystring = "DCA_{xy} (cm)";
  
  const TString systemString[2] = {"PbPb", "pp"};
  
  const Double_t CSPEED = TMath::C() * 1.e2 / 1.e12; /* cm/ps */
  
  //Pt binning
  const Double_t fBinPt[kPtBins+1] = {0.01, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};
  
  const Double_t PtRange_Hi[kCharges*nSpecies] = {.3, 5., .4, 3., .6, 5., 1.5, 5.};
  const Double_t PtRange_pp[kCharges*nSpecies] = {.3, 3., .4, 3., .6, 3., 1.5, 5.};
  extern const Double_t *PtRange[2];// = {PtRange_Hi, PtRange_pp};
  const Double_t MultBin[nMultBin] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100};
  const TString MultBinString[nMultBin] = { "0to5", "5to10", "10to20", "20to30", "30to40", "40to50", "50to60", "60to70", "70to80", "80to90", "90to100", "MB"};
  const TString MultBinStringInSquare[nMultBin] = { "[0-5]%", "[5-10]%", "[10-20]%", "[20-30]%", "[30-40]%", "[40-50]%", "[50-60]%", "[60-70]%", "[70-80]%", "[80-90]%", "[90-100]%", "MB"};
  const TString MultBinStringDash[nMultBin] = { "0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%", "MB"};
  // const Int_t multcolor[nMultBin] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, -1};//Blue Central
  const Int_t multcolor[nMultBin] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1};//Red Central
  const Int_t multdraw[nMultBin] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  
  //Cuts applied used for cut variation
  const Int_t nCuts = 5;
  const TString Cuts[nCuts] = {"TPCRows", "TrkChi2", "DCAz", "DCAxy", "GeoCut"};
  //WARNING -> in order to work correctly please make sure that the looser cuts are in the first position!! This is fundamental!!
  //==> TPC crossed rows
  const UInt_t kTPCrows = 0;
  const UInt_t nCutTPCRows = 3;
  const Double_t CutValueTPCRows[nCutTPCRows] = {60, 70, 80};
  //==> Track chi2
  const UInt_t kTrkChi2 = 1;
  const UInt_t nCutMaxChi2 = 3;
  const Double_t CutValueMaxChi2[nCutMaxChi2] = {5, 4, 3};
  //==> Track DCAz 
  const UInt_t kDCAz = 2;
  const UInt_t nCutDCAz = 3;
  const Double_t CutValueMaxDCAz[nCutDCAz] = {3, 2, 1};
  //==> Track DCAxy
  const UInt_t kDCAxy = 3;
  const UInt_t nCutDCAxy = 5;
  const Double_t CutValueMaxDCAxy[nCutDCAxy] = {10, 1., 0.9, 1.1, 0.1};
  //==> Geo cut
  const UInt_t kGeo = 4;
  const UInt_t nCutGeo = 3;
  //SetCutGeoNcrNcl(0., 0., 0.0, 0.0, 0.0)
  const Double_t GeoSetLoose[5] = {0.0, 0.0, 0.0, 0.0, 0.0};//Loose
  //SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0)
  const Double_t GeoSetStd[5] = {2., 130., 0.0, 0.0, 0.0};//Standard
  //SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
  const Double_t GeoSetTight[5] = {3., 130., 1.5, 0.85, 0.7};//Tigth
  const Double_t CutValueGeo[5*nCutGeo] = {
    GeoSetLoose[0], GeoSetLoose[1], GeoSetLoose[2], GeoSetLoose[3], GeoSetLoose[4], //Loose
    GeoSetStd[0],     GeoSetStd[1],   GeoSetStd[2],   GeoSetStd[3],   GeoSetStd[4], //Standard
    GeoSetTight[0], GeoSetTight[1], GeoSetTight[2], GeoSetTight[3], GeoSetTight[4], //Tigth
  };
  
  const UInt_t CutStdIndexInMask[nCuts] = {kTPCSetL, kTPCChi2SetL, kDCAzSetL, kPrimSetL, kGeoCutSet1};//Position in the mask of the standard cuts
  const UInt_t CutStdIndex[nCuts] = {1, 1, 1, 1, 1};//Position in the array of the standard cuts
  const UInt_t CutIndex[nCuts] = {nCutTPCRows, nCutMaxChi2, nCutDCAz, nCutDCAxy, nCutGeo};//Number of cuts values for each cut type
  extern const Double_t *CutValues[nCuts];// = {CutValueTPCRows, CutValueMaxChi2, CutValueMaxDCAz, CutValueMaxDCAxy, CutValueGeo};//Values of the cuts for each cut type
  
  const Int_t nCutVars = nCutTPCRows + nCutMaxChi2 + nCutDCAz + nCutDCAxy + nCutGeo - 5;//Number of cut present in the tree mask
  const TString CutVarsName[nCutVars + 1] = {"TPCRows_60", "TPCRows_80", "TrkChi2_5", "TrkChi2_3", "DCAz_3", "DCAz_1", "DCAxy_10", "DCAxy_11", "DCAxy_09", "DCAxy_0", "GeoCut_0", "GeoCut_1", ""};
  const TString CutVarsTitle[nCutVars + 1] = {"TPCRows 60", "TPCRows 80", "Chi2 5", "Chi2 3", "DCAz 3", "DCAz 1", "DCAxy x10", "DCAxy x0.9", "DCAxy x1.1", "DCAxy x0.1", "GeoCut 0", "GeoCut 1", "Std"};
  const TString primfunct = "0.0105+0.0350/pt^1.1";//Standard pt dependence of the DCAxy cut
  
  //Golden Chis cut
  const Double_t primchi2 = 36;
  
  const UInt_t nMCs_Hi = 2;//Number of different MC productions for PbPb
  const UInt_t nMCs_pp = 2;//Number of different MC productions for pp
  const TString MCProduction_Hi[nMCs_Hi] = {"HIJING", "DPMJET"};
  const TString MCProduction_pp[nMCs_pp] = {"PYTHIA8", "PYTHIA6"};
  const UInt_t nMCs[2]= {nMCs_Hi, nMCs_pp};
  extern const TString *MCProduction[2];// = {MCProduction_Hi, MCProduction_pp};
  
  //
  //General Mask methods
  //
  //**UChar_t**
  inline void SetMaskBit(UChar_t &mask, Int_t bit, Bool_t value){
    //     cout<< "Changing "<<bit<< " to "<<value<<endl; 
    if(value) mask = mask | 1<<bit;
    else mask = mask & ~(1<<bit);
  };
  inline Bool_t GetMaskBit(UChar_t mask, Int_t bit){
    if(mask & 1<<bit) return kTRUE; 
    return kFALSE;
  };
  inline void ResetMask(UChar_t &mask){
    mask = mask&0;
  };
  inline void PrintMaskBit(UChar_t mask){
    for(UInt_t c = 0; c < 8*sizeof(mask); c++) std::cout<<c<<" ";
    std::cout<<std::endl;
    for(UInt_t c = 0; c < 8*sizeof(mask); c++){
      if(c > 9) std::cout<<" ";
      std::cout<<GetMaskBit(mask, c)<<" ";
    }
    std::cout<<std::endl;
  }
  
  //**UShort_t**
  inline void SetMaskBit(UShort_t &mask, Int_t bit, Bool_t value){
    //     cout<< "Changing "<<bit<< " to "<<value<<endl; 
    if(value) mask = mask | 1<<bit;
    else mask = mask & ~(1<<bit);
  };
  inline Bool_t GetMaskBit(UShort_t mask, Int_t bit){
    if(mask & 1<<bit) return kTRUE; 
    return kFALSE;
  };
  inline void ResetMask(UShort_t &mask){
    mask = mask&0;
  };
  inline void PrintMaskBit(UShort_t mask){
    for(UInt_t c = 0; c < 8*sizeof(mask); c++) std::cout<<c<<" ";
    std::cout<<std::endl;
    for(UInt_t c = 0; c < 8*sizeof(mask); c++){
      if(c > 9) std::cout<<" ";
      std::cout<<GetMaskBit(mask, c)<<" ";
    }
    std::cout<<std::endl;
  }
  
}
#endif
