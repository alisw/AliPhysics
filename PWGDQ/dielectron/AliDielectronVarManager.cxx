/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron Variables Manager class                     //
//                                                                       //
/*

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliDielectronVarManager.h"

ClassImp(AliDielectronVarManager)

const char* AliDielectronVarManager::fgkParticleNames[AliDielectronVarManager::kNMaxValues] = {
  "Px",
  "Py",
  "Pz",
  "Pt",
  "P",
  "Xv",
  "Yv",
  "Zv",
  "OneOverPt",
  "Phi",
  "Theta",
  "Eta",
  "Y",
  "E",
  "M",
  "Charge",
  "NclsITS",
  "ITSchi2PerCluster",
  "NclsTPC",
  "NclsSTPC",
  "NclsSFracTPC",
  "NclsTPCiter1",
  "NFclsTPC",
  "NFclsTPCrobust",
  "NFclsTPCrobustFraction",
  "TPCsignalN",
  "TPCsignalNfrac",
  "TPCchi2PerCluster",
  "TPCclsDiff",
  "TrackStatus",
    
  "NclsTRD",
  "TRDntracklets",
  "TRDpidQuality",
  "TRDpidProb_Electrons",
  "TRDpidProb_Pions",
  "TRDphi",
  "TRDpidEffLeg",

  "ImpactParXY",
  "ImpactParZ",
  "TrackLength",

  "PdgCode",
  "PdgCodeMother",
  "PdgCodeGrandMother",
  "NumberOfDaughters",
  "HaveSameMother",
  "IsJpsiPrimary",
  "NumberOfJPsisIncl",
  "NumberOfJPsisPrompt",
  "NumberOfJPsisNPrompt",
  
  "ITS_signal",
  "SSD1_signal",
  "SSD2_signal",
  "SDD1_signal",
  "SDD2_signal",
  "ITS_clusterMap",
  "ITSLayerFirstCls",
  "ITS_nSigma_Electrons",
  "ITS_nSigma_Pions",
  "ITS_nSigma_Muons",
  "ITS_nSigma_Kaons",
  "ITS_nSigma_Protons",

  "P_InnerParam",
  "P_OuterParam",
  "Y_signed_InnerParam",
  "TPC_signal",
  "TOF_signal",
  "TOF_beta",
  "TOF_PIDbit",
  
  "TPC_nSigma_Electrons",
  "TPC_nSigma_Pions",
  "TPC_nSigma_Muons",
  "TPC_nSigma_Kaons",
  "TPC_nSigma_Protons",

  "TOF_nSigma_Electrons",
  "TOF_nSigma_Pions",
  "TOF_nSigma_Muons",
  "TOF_nSigma_Kaons",
  "TOF_nSigma_Protons",

  "EMCAL_nSigma_Electrons",

  "KinkIndex0",
  //
  "Chi2NDF",
  "DecayLength",
  "R",
  "OpeningAngle",
  "ThetaHE",
  "PhiHE",
  "Cos2PhiHE",
  "ThetaCS",
  "PhiCS",
  "Cos2PhiCS",
  "DeltaPhiV0ArpH2",        
  "DeltaPhiV0CrpH2",        
  "DeltaPhiV0ACrpH2",       
  "V0ArpH2FlowV2",         
  "V0CrpH2FlowV2",      
  "V0ACrpH2FlowV2",
  "LegDistance",
  "LegDistanceXY",
  "DeltaEta",
  "DeltaPhi",
  "Merr",
  "DCA",
  "PairType",
  "PseudoProperTime",
  "PseudoProperTimeResolution",
  "PseudoProperTimePull",
  "TRDpidEffPair",
  //
  "X",
  "Y",
  "Z",
  "XRes",
  "YRes",
  "ZRes",

  "v0ArpH2",
  "v0CrpH2",
  "v0ACrpH2",
  "v0ATPCDiffH2",
  "v0CTPCDiffH2",
  "v0Av0CDiffH2",

  "MultV0A", 
  "MultV0C",
  "MultV0",
  "AdcV0A",
  "AdcV0C",
  "AdcV0",  
  "VZERO_ch0",  "VZERO_ch1",  "VZERO_ch2",  "VZERO_ch3",  "VZERO_ch4",  "VZERO_ch5",  "VZERO_ch6",  "VZERO_ch7",  "VZERO_ch8",  "VZERO_ch9",
  "VZERO_ch10", "VZERO_ch11", "VZERO_ch12", "VZERO_ch13", "VZERO_ch14", "VZERO_ch15", "VZERO_ch16", "VZERO_ch17", "VZERO_ch18", "VZERO_ch19",
  "VZERO_ch20", "VZERO_ch21", "VZERO_ch22", "VZERO_ch23", "VZERO_ch24", "VZERO_ch25", "VZERO_ch26", "VZERO_ch27", "VZERO_ch28", "VZERO_ch29",
  "VZERO_ch30", "VZERO_ch31", "VZERO_ch32", "VZERO_ch33", "VZERO_ch34", "VZERO_ch35", "VZERO_ch36", "VZERO_ch37", "VZERO_ch38", "VZERO_ch39",
  "VZERO_ch40", "VZERO_ch41", "VZERO_ch42", "VZERO_ch43", "VZERO_ch44", "VZERO_ch45", "VZERO_ch46", "VZERO_ch47", "VZERO_ch48", "VZERO_ch49",
  "VZERO_ch50", "VZERO_ch51", "VZERO_ch52", "VZERO_ch53", "VZERO_ch54", "VZERO_ch55", "VZERO_ch56", "VZERO_ch57", "VZERO_ch58", "VZERO_ch59",
  "VZERO_ch60", "VZERO_ch61", "VZERO_ch62", "VZERO_ch63",
  "V0AxH2", 
  "V0AyH2",
  "V0ArpH2",
  "V0CxH2",
  "V0CyH2",
  "V0CrpH2",
  "V0ACxH2",
  "V0ACyH2",
  "V0ACrpH2",
  "V0ArpResH2",
  "V0CrpResH2",
  "V0ACrpResH2",
  "V0XaXcH2",
  "V0XaYaH2",
  "V0XaYcH2",
  "V0YaXcH2",
  "V0YaYcH2",
  "V0XcYcH2",
  "V0ATPCDiffH2",
  "V0CTPCDiffH2",
  "V0AV0CDiffH2",
  "TPCxH2",
  "TPCyH2",
  "TPCrpH2",
  "TPCsub1xH2",
  "TPCsub1yH2",
  "TPCsub1rpH2",
  "TPCsub2xH2",
  "TPCsub2yH2",
  "TPCsub2rpH2",
  "TPCsub12DiffH2",
  "TPCsub12DiffH2Sin",

  "TPCxH2uc",
  "TPCyH2uc",
  "TPCrpH2uc",
  "TPCsub1xH2uc",
  "TPCsub1yH2uc",
  "TPCsub1rpH2uc",
  "TPCsub2xH2uc",
  "TPCsub2yH2uc",
  "TPCsub2rpH2uc",
  "TPCsub12DiffH2uc",

  "NTrk",
  "Tracks",
  "Nacc",
  "NaccTrcklts",
  "NaccTrcklts0916",
  
  "NaccTrckltsEsd05",
  "NaccTrckltsEsd10",
  "NaccTrckltsEsd16",
  "NaccTrckltsEsd05Corr",
  "NaccTrckltsEsd10Corr",
  "NaccTrckltsEsd16Corr",
  "NaccItsTpcEsd05",
  "NaccItsTpcEsd10",
  "NaccItsTpcEsd16",
  "NaccItsTpcEsd05Corr",
  "NaccItsTpcEsd10Corr",
  "NaccItsTpcEsd16Corr",
  "NaccItsPureEsd05",
  "NaccItsPureEsd10",
  "NaccItsPureEsd16",
  "NaccItsPureEsd05Corr",
  "NaccItsPureEsd10Corr",
  "NaccItsPureEsd16Corr",  
    
  "Nch",                    // Number of charged MC tracks in |eta|<1.6
  "Nch05",                  // Number of charged MC tracks in |eta|<0.5
  "Nch10",                  // Number of charged MC tracks in |eta|<1.0
  "Centrality",
  "Nevents",
  "RunNumber",
  "MixingBin"
};

AliPIDResponse* AliDielectronVarManager::fgPIDResponse      = 0x0;
AliVEvent*      AliDielectronVarManager::fgEvent            = 0x0;
AliEventplane*  AliDielectronVarManager::fgTPCEventPlane    = 0x0;
AliKFVertex*    AliDielectronVarManager::fgKFVertex         = 0x0;
TProfile*       AliDielectronVarManager::fgMultEstimatorAvg[4][9] = {{0x0}};
TH3D*           AliDielectronVarManager::fgTRDpidEff[10][4] = {{0x0}};
Double_t        AliDielectronVarManager::fgTRDpidEffCentRanges[10][4] = {{0.0}};
Double_t        AliDielectronVarManager::fgData[AliDielectronVarManager::kNMaxValues] = {};
//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager() :
  TNamed("AliDielectronVarManager","AliDielectronVarManager")
{
  //
  // Default constructor
  //
  for(Int_t i=0; i<4; ++i)
    for(Int_t j=0; j<9; ++j)
      fgMultEstimatorAvg[i][j] = 0x0;
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      fgTRDpidEff[i][j] = 0x0;
}

//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager(const char* name, const char* title) :
  TNamed(name,title)
{
  //
  // Named constructor
  //
  for(Int_t i=0; i<4; ++i)
    for(Int_t j=0; j<9; ++j)
      fgMultEstimatorAvg[i][j] = 0x0;
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      fgTRDpidEff[i][j] = 0x0;  
}

//________________________________________________________________
AliDielectronVarManager::~AliDielectronVarManager()
{
  //
  // Default destructor
  //
  for(Int_t i=0; i<4; ++i)
    for(Int_t j=0; j<9; ++j)
      if(fgMultEstimatorAvg[i][j]) delete fgMultEstimatorAvg[i][j];
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      if(fgTRDpidEff[i][j]) delete fgTRDpidEff[i][j];    
}

