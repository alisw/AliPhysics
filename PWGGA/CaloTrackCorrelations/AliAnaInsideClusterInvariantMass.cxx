/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
//
// Split clusters with some criteria and calculate invariant mass
// to identify them as pi0 or conversion
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)  
//_________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TList.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TH2F.h>

// --- Analysis system --- 
#include "AliAnaInsideClusterInvariantMass.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliEMCALGeoParams.h"

// --- Detectors --- 
//#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnaInsideClusterInvariantMass)
  
//__________________________________________________________________
AliAnaInsideClusterInvariantMass::AliAnaInsideClusterInvariantMass() : 
  AliAnaCaloTrackCorrBaseClass(),  
  fCalorimeter(""), 
  fM02MaxCut(0),    fM02MinCut(0),       
  fMinNCells(0),    fMinBadDist(0),
  fFillAngleHisto(kFALSE),
  fFillTMResidualHisto(kFALSE),
  fFillSSExtraHisto(kFALSE),
  fFillMCFractionHisto(kFALSE),
  fhMassM02CutNLocMax1(0),    fhMassM02CutNLocMax2(0),    fhMassM02CutNLocMaxN(0),
  fhAsymM02CutNLocMax1(0),    fhAsymM02CutNLocMax2(0),    fhAsymM02CutNLocMaxN(0),
  fhMassSplitECutNLocMax1(0), fhMassSplitECutNLocMax2(0), fhMassSplitECutNLocMaxN(0),
  fhMCGenFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenFracAfterCutsNLocMaxNMCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0(0),
  fhSplitEFractionAfterCutsNLocMax1(0),
  fhSplitEFractionAfterCutsNLocMax2(0),
  fhSplitEFractionAfterCutsNLocMaxN(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 7; i++)
  {
    fhMassAfterCutsNLocMax1[i] = 0;
    fhMassAfterCutsNLocMax2[i] = 0;
    fhMassAfterCutsNLocMaxN[i] = 0;
    
    for(Int_t j = 0; j < 2; j++)
    {
      fhMassNLocMax1[i][j]  = 0;
      fhMassNLocMax2[i][j]  = 0;
      fhMassNLocMaxN[i][j]  = 0;
      fhNLocMax[i][j]       = 0;
      fhNLocMaxM02Cut[i][j] = 0;
      fhM02NLocMax1[i][j]   = 0;
      fhM02NLocMax2[i][j]   = 0;
      fhM02NLocMaxN[i][j]   = 0;
      fhNCellNLocMax1[i][j] = 0;
      fhNCellNLocMax2[i][j] = 0;
      fhNCellNLocMaxN[i][j] = 0;
      fhM02Pi0LocMax1[i][j] = 0;
      fhM02EtaLocMax1[i][j] = 0;
      fhM02ConLocMax1[i][j] = 0;
      fhM02Pi0LocMax2[i][j] = 0;
      fhM02EtaLocMax2[i][j] = 0;
      fhM02ConLocMax2[i][j] = 0;
      fhM02Pi0LocMaxN[i][j] = 0;
      fhM02EtaLocMaxN[i][j] = 0;
      fhM02ConLocMaxN[i][j] = 0;
      
      fhMassPi0LocMax1[i][j] = 0;
      fhMassEtaLocMax1[i][j] = 0;
      fhMassConLocMax1[i][j] = 0;
      fhMassPi0LocMax2[i][j] = 0;
      fhMassEtaLocMax2[i][j] = 0;
      fhMassConLocMax2[i][j] = 0;
      fhMassPi0LocMaxN[i][j] = 0;
      fhMassEtaLocMaxN[i][j] = 0;
      fhMassConLocMaxN[i][j] = 0;
      
      
      fhAsyPi0LocMax1[i][j] = 0;
      fhAsyEtaLocMax1[i][j] = 0;
      fhAsyConLocMax1[i][j] = 0;
      fhAsyPi0LocMax2[i][j] = 0;
      fhAsyEtaLocMax2[i][j] = 0;
      fhAsyConLocMax2[i][j] = 0;
      fhAsyPi0LocMaxN[i][j] = 0;
      fhAsyEtaLocMaxN[i][j] = 0;
      fhAsyConLocMaxN[i][j] = 0;      
      
      fhMassM02NLocMax1[i][j]= 0;
      fhMassM02NLocMax2[i][j]= 0;
      fhMassM02NLocMaxN[i][j]= 0;   
      fhMassDispEtaNLocMax1[i][j]= 0;
      fhMassDispEtaNLocMax2[i][j]= 0;
      fhMassDispEtaNLocMaxN[i][j]= 0;      
      fhMassDispPhiNLocMax1[i][j]= 0;
      fhMassDispPhiNLocMax2[i][j]= 0;
      fhMassDispPhiNLocMaxN[i][j]= 0;      
      fhMassDispAsyNLocMax1[i][j]= 0;
      fhMassDispAsyNLocMax2[i][j]= 0;
      fhMassDispAsyNLocMaxN[i][j]= 0;      
      
      fhSplitEFractionNLocMax1[i][j]=0;
      fhSplitEFractionNLocMax2[i][j]=0;
      fhSplitEFractionNLocMaxN[i][j]=0;
      
      fhMCGenFracNLocMax1[i][j]= 0;
      fhMCGenFracNLocMax2[i][j]= 0;
      fhMCGenFracNLocMaxN[i][j]= 0;
      
      fhMCGenSplitEFracNLocMax1[i][j]= 0;
      fhMCGenSplitEFracNLocMax2[i][j]= 0;
      fhMCGenSplitEFracNLocMaxN[i][j]= 0;    
      
      fhMCGenEFracvsSplitEFracNLocMax1[i][j]= 0;
      fhMCGenEFracvsSplitEFracNLocMax2[i][j]= 0;
      fhMCGenEFracvsSplitEFracNLocMaxN[i][j]= 0;    
      
      fhMCGenEvsSplitENLocMax1[i][j]= 0;
      fhMCGenEvsSplitENLocMax2[i][j]= 0;
      fhMCGenEvsSplitENLocMaxN[i][j]= 0;     
      
      fhAsymNLocMax1 [i][j] = 0;
      fhAsymNLocMax2 [i][j] = 0;
      fhAsymNLocMaxN [i][j] = 0;
    }
    
    for(Int_t jj = 0; jj < 4; jj++)
    {
      fhM02MCGenFracNLocMax1Ebin[i][jj] = 0;
      fhM02MCGenFracNLocMax2Ebin[i][jj] = 0;
      fhM02MCGenFracNLocMaxNEbin[i][jj] = 0;
      
      fhMassMCGenFracNLocMax1Ebin[i][jj]= 0;
      fhMassMCGenFracNLocMax2Ebin[i][jj]= 0;
      fhMassMCGenFracNLocMaxNEbin[i][jj]= 0;
      
      fhMCGenFracNLocMaxEbin[i][jj]       = 0;
      fhMCGenFracNLocMaxEbinMatched[i][jj]= 0;
      
      fhMassSplitEFractionNLocMax1Ebin[i][jj] = 0;
      fhMassSplitEFractionNLocMax2Ebin[i][jj] = 0;
      fhMassSplitEFractionNLocMaxNEbin[i][jj] = 0;
    }
    
    fhTrackMatchedDEtaLocMax1[i] = 0; 
    fhTrackMatchedDPhiLocMax1[i] = 0;
    fhTrackMatchedDEtaLocMax2[i] = 0;
    fhTrackMatchedDPhiLocMax2[i] = 0; 
    fhTrackMatchedDEtaLocMaxN[i] = 0; 
    fhTrackMatchedDPhiLocMaxN[i] = 0; 
    
  }
   
  for(Int_t i = 0; i < 2; i++)
  {
    fhAnglePairLocMax1    [i] = 0;
    fhAnglePairLocMax2    [i] = 0;
    fhAnglePairLocMaxN    [i] = 0;
    fhAnglePairMassLocMax1[i] = 0;
    fhAnglePairMassLocMax2[i] = 0;
    fhAnglePairMassLocMaxN[i] = 0;
    fhSplitEFractionvsAsyNLocMax1[i] = 0;
    fhSplitEFractionvsAsyNLocMax2[i] = 0; 
    fhSplitEFractionvsAsyNLocMaxN[i] = 0;    
  }
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhMassM02NLocMax1Ebin[i] = 0 ;
    fhMassM02NLocMax2Ebin[i] = 0 ;
    fhMassM02NLocMaxNEbin[i] = 0 ;

    fhMassAsyNLocMax1Ebin[i] = 0 ;
    fhMassAsyNLocMax2Ebin[i] = 0 ;
    fhMassAsyNLocMaxNEbin[i] = 0 ;

    fhAsyMCGenRecoNLocMax1EbinPi0[i] = 0 ;
    fhAsyMCGenRecoNLocMax2EbinPi0[i] = 0 ;
    fhAsyMCGenRecoNLocMaxNEbinPi0[i] = 0 ;
    
    fhMassDispEtaNLocMax1Ebin[i] = 0 ;
    fhMassDispEtaNLocMax2Ebin[i] = 0 ;
    fhMassDispEtaNLocMaxNEbin[i] = 0 ;
    
    fhMassDispPhiNLocMax1Ebin[i] = 0 ;
    fhMassDispPhiNLocMax2Ebin[i] = 0 ;
    fhMassDispPhiNLocMaxNEbin[i] = 0 ;    
    
    fhMassDispAsyNLocMax1Ebin[i] = 0 ;
    fhMassDispAsyNLocMax2Ebin[i] = 0 ;
    fhMassDispAsyNLocMaxNEbin[i] = 0 ;    

    fhMCAsymM02NLocMax1MCPi0Ebin[i] = 0 ;
    fhMCAsymM02NLocMax2MCPi0Ebin[i] = 0 ;
    fhMCAsymM02NLocMaxNMCPi0Ebin[i] = 0 ;
  }
  
  InitParameters();
  
}

//_______________________________________________________________
TObjString *  AliAnaInsideClusterInvariantMass::GetAnalysisCuts()
{	
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaInsideClusterInvariantMass ---\n") ;
  parList+=onePar ;	
  
  snprintf(onePar,buffersize,"Calorimeter: %s\n",        fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fLocMaxCutE =%2.2f \n",    GetCaloUtils()->GetLocalMaximaCutE()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fLocMaxCutEDiff =%2.2f \n",GetCaloUtils()->GetLocalMaximaCutEDiff()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"%2.2f< M02 < %2.2f \n",    fM02MinCut, fM02MaxCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinNCells =%d \n",        fMinNCells) ;
  parList+=onePar ;    
  snprintf(onePar,buffersize,"fMinBadDist =%1.1f \n",    fMinBadDist) ;
  parList+=onePar ;  

  return new TObjString(parList) ;
  
}

//________________________________________________________________
TList * AliAnaInsideClusterInvariantMass::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("InsideClusterHistos") ;
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();           Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();           Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax  = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin  = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t mbins    = GetHistogramRanges()->GetHistoMassBins();         Float_t mmax   = GetHistogramRanges()->GetHistoMassMax();         Float_t mmin   = GetHistogramRanges()->GetHistoMassMin();
  Int_t ncbins   = GetHistogramRanges()->GetHistoNClusterCellBins(); Int_t   ncmax  = GetHistogramRanges()->GetHistoNClusterCellMax(); Int_t   ncmin  = GetHistogramRanges()->GetHistoNClusterCellMin(); 

  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();          
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();          
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();          
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();          
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();  
  
  TString ptype[] ={"","#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron"}; 
  TString pname[] ={"","Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron"};
  
  Int_t n = 1;
  
  if(IsDataMC()) n = 7;
  
  Int_t nMaxBins = 10;
  
  TString sMatched[] = {"","Matched"};
  
  
  fhMassSplitECutNLocMax1  = new TH2F("hMassSplitECutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, (E1+E2)/E cut",
                                   nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassSplitECutNLocMax1->SetYTitle("M (GeV/c^{2})");
  fhMassSplitECutNLocMax1->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassSplitECutNLocMax1) ;   
  
  fhMassSplitECutNLocMax2  = new TH2F("hMassSplitECutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, (E1+E2)/E cut",
                                   nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassSplitECutNLocMax2->SetYTitle("M (GeV/c^{2})");
  fhMassSplitECutNLocMax2->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassSplitECutNLocMax2) ;   
  
  fhMassSplitECutNLocMaxN  = new TH2F("hMassSplitECutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, (E1+E2)/E cut",
                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassSplitECutNLocMaxN->SetYTitle("M (GeV/c^{2})");
  fhMassSplitECutNLocMaxN->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassSplitECutNLocMaxN) ;   

  fhMassM02CutNLocMax1  = new TH2F("hMassM02CutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, M02 cut",
                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassM02CutNLocMax1->SetYTitle("M (GeV/c^{2})");
  fhMassM02CutNLocMax1->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassM02CutNLocMax1) ;   
  
  fhMassM02CutNLocMax2  = new TH2F("hMassM02CutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, M02 cut",
                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassM02CutNLocMax2->SetYTitle("M (GeV/c^{2})");
  fhMassM02CutNLocMax2->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassM02CutNLocMax2) ;   
  
  fhMassM02CutNLocMaxN  = new TH2F("hMassM02CutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, M02 cut",
                             nptbins,ptmin,ptmax,mbins,mmin,mmax); 
  fhMassM02CutNLocMaxN->SetYTitle("M (GeV/c^{2})");
  fhMassM02CutNLocMaxN->SetXTitle("E (GeV)");
  outputContainer->Add(fhMassM02CutNLocMaxN) ;
  
  fhAsymM02CutNLocMax1  = new TH2F("hAsymM02CutNLocMax1","Asymmetry of NLM=1  vs cluster Energy, M02Cut", nptbins,ptmin,ptmax,200,-1,1);
  fhAsymM02CutNLocMax1->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
  fhAsymM02CutNLocMax1->SetXTitle("E (GeV)");
  outputContainer->Add(fhAsymM02CutNLocMax1) ;
  
  fhAsymM02CutNLocMax2  = new TH2F("hAsymM02CutNLocMax2","Asymmetry of NLM=2  vs cluster Energy, M02Cut", nptbins,ptmin,ptmax,200,-1,1);
  fhAsymM02CutNLocMax2->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
  fhAsymM02CutNLocMax2->SetXTitle("E (GeV)");
  outputContainer->Add(fhAsymM02CutNLocMax2) ;

  fhAsymM02CutNLocMaxN  = new TH2F("hAsymM02CutNLocMaxN","Asymmetry of NLM>2  vs cluster Energy, M02Cut", nptbins,ptmin,ptmax,200,-1,1);
  fhAsymM02CutNLocMaxN->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
  fhAsymM02CutNLocMaxN->SetXTitle("E (GeV)");
  outputContainer->Add(fhAsymM02CutNLocMaxN) ;
  
  fhSplitEFractionAfterCutsNLocMax1     = new TH2F("hSplitEFractionAfterCutsNLocMax1",
                                                   "(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 1, M02 and Asy cut on",
                                                nptbins,ptmin,ptmax,120,0,1.2);
  fhSplitEFractionAfterCutsNLocMax1   ->SetXTitle("E_{cluster} (GeV)");
  fhSplitEFractionAfterCutsNLocMax1   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionAfterCutsNLocMax1) ;

  fhSplitEFractionAfterCutsNLocMax2     = new TH2F("hSplitEFractionAfterCutsNLocMax2",
                                                   "(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 2, M02 and Asy cut on",
                                                nptbins,ptmin,ptmax,120,0,1.2);
  fhSplitEFractionAfterCutsNLocMax2   ->SetXTitle("E_{cluster} (GeV)");
  fhSplitEFractionAfterCutsNLocMax2   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionAfterCutsNLocMax2) ;
  
  fhSplitEFractionAfterCutsNLocMaxN    = new TH2F("hSplitEFractionAfterCutsNLocMaxN",
                                                  "(E1+E2)/E_{cluster} vs E_{cluster} for N max  > 2, M02 and Asy cut on",
                                               nptbins,ptmin,ptmax,120,0,1.2);
  fhSplitEFractionAfterCutsNLocMaxN   ->SetXTitle("E_{cluster} (GeV)");
  fhSplitEFractionAfterCutsNLocMaxN   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionAfterCutsNLocMaxN) ;
  
  if(IsDataMC() && fFillMCFractionHisto)
  {
    
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0     = new TH2F("hMCGenSplitEFracAfterCutsNLocMax1MCPi0",
                                                           "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 1 MC Pi0, after M02 and Asym cut",
                                                           nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMax1MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMax1MCPi0) ;
    
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0    = new TH2F("hMCGenSplitEFracAfterCutsNLocMax2MCPi0",
                                                          "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 2 MC Pi0, after M02 and Asym cut",
                                                          nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0  ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMax2MCPi0  ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMax2MCPi0) ;
    
    
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0    = new TH2F("hMCGenSplitEFracAfterCutsNLocMaxNMCPi0",
                                                          "E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  > 2 MC Pi0, after M02 and Asym cut",
                                                          nptbins,ptmin,ptmax,200,0,2);
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0  ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
    fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0  ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0) ;
    
    fhMCGenFracAfterCutsNLocMax1MCPi0     = new TH2F("hMCGenFracAfterCutsNLocMax1MCPi0",
                                                                  "E_{gen} / E_{reco} vs E_{reco} for N max  = 1 MC Pi0, after M02 and Asym cut",
                                                                  nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMax1MCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMax1MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMax1MCPi0) ;
    
    fhMCGenFracAfterCutsNLocMax2MCPi0    = new TH2F("hMCGenFracAfterCutsNLocMax2MCPi0",
                                                                 " E_{gen} / E_{reco} vs E_{reco} for N max  = 2 MC Pi0, after M02 and Asym cut",
                                                                 nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMax2MCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMax2MCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMax2MCPi0) ;
    
    
    fhMCGenFracAfterCutsNLocMaxNMCPi0   = new TH2F("hMCGenFracAfterCutsNLocMaxNMCPi0",
                                                                " E_{gen} / E_{reco}  vs E_{reco} for N max  > 2 MC Pi0, after M02 and Asym cut",
                                                                nptbins,ptmin,ptmax,200,0,2);
    fhMCGenFracAfterCutsNLocMaxNMCPi0   ->SetYTitle("E_{gen} / E_{reco}");
    fhMCGenFracAfterCutsNLocMaxNMCPi0   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhMCGenFracAfterCutsNLocMaxNMCPi0) ;
    
  }
  
  for(Int_t i = 0; i < n; i++)
  {
    for(Int_t j = 0; j < 2; j++)
    {
      
      fhMassNLocMax1[i][j]  = new TH2F(Form("hMassNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=1 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax1[i][j]) ;   
      
      fhMassNLocMax2[i][j]  = new TH2F(Form("hMassNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM=2 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax2[i][j]) ;   
      
      fhMassNLocMaxN[i][j]  = new TH2F(Form("hMassNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of splitted cluster with NLM>2 vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMaxN[i][j]) ;
     
      if(j==0)
      {
        fhMassAfterCutsNLocMax1[i]     = new TH2F(Form("hMassAfterCutsNLocMax1%s",pname[i].Data()),
                                                 Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 1, m02 and asy cut",
                                                      GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                                 nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMax1[i]   ->SetYTitle("Mass (MeV/c^{2})");
        fhMassAfterCutsNLocMax1[i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMax1[i]) ;
        
        fhMassAfterCutsNLocMax2[i]     = new TH2F(Form("hMassAfterCutsNLocMax2%s",pname[i].Data()),
                                                    Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 2, asy cut",
                                                         GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                                    nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMax2[i]   ->SetYTitle("Mass (MeV/c^{2})");
        fhMassAfterCutsNLocMax2[i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMax2[i]) ;
        
        
        fhMassAfterCutsNLocMaxN[i]     = new TH2F(Form("hMassAfterCutsNLocMaxN%s",pname[i].Data()),
                                                    Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max > 2, asy cut",
                                                         GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                                    nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassAfterCutsNLocMaxN[i]   ->SetYTitle("Mass (MeV/c^{2})");
        fhMassAfterCutsNLocMaxN[i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassAfterCutsNLocMaxN[i]) ;
      }
      
      fhMassM02NLocMax1[i][j]  = new TH2F(Form("hMassM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=1, #lambda_{0}^{2}, E > 8 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1[i][j]) ;   
      
      fhMassM02NLocMax2[i][j]  = new TH2F(Form("hMassM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=2, #lambda_{0}^{2}, E > 8 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax2[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax2[i][j]) ;   
      
      fhMassM02NLocMaxN[i][j]  = new TH2F(Form("hMassM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM>2, vs #lambda_{0}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMaxN[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMaxN[i][j]) ;   
      
      
      fhAsymNLocMax1[i][j]  = new TH2F(Form("hAsymNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                    Form("Asymmetry of NLM=1  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                    nptbins,ptmin,ptmax,200,-1,1); 
      fhAsymNLocMax1[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMax1[i][j]) ;   
      
      fhAsymNLocMax2[i][j]  = new TH2F(Form("hAsymNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                    Form("Asymmetry of NLM=2  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                    nptbins,ptmin,ptmax,200,-1,1); 
      fhAsymNLocMax2[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMax2[i][j]) ;   
      
      fhAsymNLocMaxN[i][j]  = new TH2F(Form("hAsymNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                    Form("Asymmetry of NLM>2  vs cluster Energy, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                    nptbins,ptmin,ptmax,200,-1,1); 
      fhAsymNLocMaxN[i][j]->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
      fhAsymNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsymNLocMaxN[i][j]) ;   
      
      
      if(fFillSSExtraHisto)
      {
        fhMassDispEtaNLocMax1[i][j]  = new TH2F(Form("hMassDispEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=1, #sigma_{#eta #eta}^{2}, E > 8 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax1[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax1[i][j]) ;   
        
        fhMassDispEtaNLocMax2[i][j]  = new TH2F(Form("hMassDispEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=2 #sigma_{#eta #eta}^{2}, E > 8 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax2[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax2[i][j]) ;   
        
        fhMassDispEtaNLocMaxN[i][j]  = new TH2F(Form("hMassDispEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM>2, #sigma_{#eta #eta}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMaxN[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMaxN[i][j]) ;   
        
        fhMassDispPhiNLocMax1[i][j]  = new TH2F(Form("hMassDispPhiNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, E > 8 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax1[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax1[i][j]) ;   
        
        fhMassDispPhiNLocMax2[i][j]  = new TH2F(Form("hMassDispPhiNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, E > 8 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax2[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax2[i][j]) ;   
        
        fhMassDispPhiNLocMaxN[i][j]  = new TH2F(Form("hMassDispPhiNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of N>2 local maxima cells vs #sigma_{#phi #phi}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMaxN[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMaxN[i][j]) ;   
        
        fhMassDispAsyNLocMax1[i][j]  = new TH2F(Form("hMassDispAsyNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 8 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax1[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax1[i][j]) ;   
        
        fhMassDispAsyNLocMax2[i][j]  = new TH2F(Form("hMassDispAsyNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 8 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax2[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax2[i][j]) ;   
        
        fhMassDispAsyNLocMaxN[i][j]  = new TH2F(Form("hMassDispAsyNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of N>2 local maxima cells vsA = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMaxN[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMaxN[i][j]) ;   
      }
      
      fhNLocMax[i][j]     = new TH2F(Form("hNLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in cluster %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMax[i][j]   ->SetYTitle("N maxima");
      fhNLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMax[i][j]) ; 
            
      fhNLocMaxM02Cut[i][j] = new TH2F(Form("hNLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Number of local maxima in cluster %s for %2.2f < M02 < %2.2f",ptype[i].Data(),fM02MinCut,fM02MaxCut),
                                       nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMaxM02Cut[i][j]->SetYTitle("N maxima");
      fhNLocMaxM02Cut[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMaxM02Cut[i][j]) ; 
      
      
      fhM02NLocMax1[i][j]     = new TH2F(Form("hM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02NLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMax1[i][j]) ; 
      
      fhM02NLocMax2[i][j]     = new TH2F(Form("hM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02NLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMax2[i][j]) ; 
      
      fhM02NLocMaxN[i][j]    = new TH2F(Form("hM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02NLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02NLocMaxN[i][j]) ; 
      
      
      fhSplitEFractionNLocMax1[i][j]     = new TH2F(Form("hSplitEFractionNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,120,0,1.2); 
      fhSplitEFractionNLocMax1[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMax1[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMax1[i][j]) ; 
      
      fhSplitEFractionNLocMax2[i][j]     = new TH2F(Form("hSplitEFractionNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins,ptmin,ptmax,120,0,1.2); 
      fhSplitEFractionNLocMax2[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMax2[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMax2[i][j]) ; 
      
      fhSplitEFractionNLocMaxN[i][j]    = new TH2F(Form("hSplitEFractionNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                        Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,120,0,1.2); 
      fhSplitEFractionNLocMaxN[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionNLocMaxN[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionNLocMaxN[i][j]) ; 
      
      
      if(i > 0 && fFillMCFractionHisto) // skip first entry in array, general case not filled
      {
        fhMCGenFracNLocMax1[i][j]     = new TH2F(Form("hMCGenFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMax1[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax1[i][j]) ; 
        
        fhMCGenFracNLocMax2[i][j]     = new TH2F(Form("hMCGenFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMax2[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMax2[i][j]) ; 
        
        
        fhMCGenFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenFracNLocMaxN[i][j]   ->SetYTitle("E_{gen} / E_{reco}");
        fhMCGenFracNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenFracNLocMaxN[i][j]) ; 
      
        fhMCGenSplitEFracNLocMax1[i][j]     = new TH2F(Form("hMCGenSplitEFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMax1[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax1[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax1[i][j]) ; 
        
        fhMCGenSplitEFracNLocMax2[i][j]     = new TH2F(Form("hMCGenSplitEFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                 Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                 nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMax2[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMax2[i][j]) ; 
        
        
        fhMCGenSplitEFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenSplitEFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E_{gen} / (E_{1 split}+E_{2 split}) vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                nptbins,ptmin,ptmax,200,0,2); 
        fhMCGenSplitEFracNLocMaxN[i][j]   ->SetYTitle("E_{gen} / (E_{1 split}+E_{2 split})");
        fhMCGenSplitEFracNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhMCGenSplitEFracNLocMaxN[i][j]) ; 
       
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]     = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                       Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                       200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMax1[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMax1[i][j]) ; 
        
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]     = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                       Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                       200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMax2[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMax2[i][j]) ; 
        
        
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]    = new TH2F(Form("hMCGenEFracvsSplitEFracNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                      Form("(E_{1 split}+E_{2 split})/E_{reco} vs E_{gen} / E_{reco} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                      200,0,2,200,0,2); 
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]   ->SetYTitle("(E_{1 split}+E_{2 split})/E_{reco}");
        fhMCGenEFracvsSplitEFracNLocMaxN[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
        outputContainer->Add(fhMCGenEFracvsSplitEFracNLocMaxN[i][j]) ; 
        
        
        fhMCGenEvsSplitENLocMax1[i][j]     = new TH2F(Form("hMCGenEvsSplitENLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                             Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMax1[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMax1[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMax1[i][j]) ; 
        
        fhMCGenEvsSplitENLocMax2[i][j]     = new TH2F(Form("hMCGenEvsSplitENLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                             Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                             nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMax2[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMax2[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMax2[i][j]) ; 
        
        
        fhMCGenEvsSplitENLocMaxN[i][j]    = new TH2F(Form("hMCGenEvsSplitENLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                            Form("E_{1 split}+E_{2 split} vs E_{gen} for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                            nptbins,ptmin,ptmax,nptbins,ptmin,ptmax); 
        fhMCGenEvsSplitENLocMaxN[i][j]   ->SetYTitle("E_{1 split}+E_{2 split} (GeV)");
        fhMCGenEvsSplitENLocMaxN[i][j]   ->SetXTitle("E_{gen} (GeV)");
        outputContainer->Add(fhMCGenEvsSplitENLocMaxN[i][j]) ; 
        
      }
      
      if(fFillSSExtraHisto)
      {
        fhNCellNLocMax1[i][j]  = new TH2F(Form("hNCellNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMax1[i][j] ->SetYTitle("N cells");
        fhNCellNLocMax1[i][j] ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax1[i][j]) ; 
        
        fhNCellNLocMax2[i][j]     = new TH2F(Form("hNCellNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E for N max  = 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMax2[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMax2[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMax2[i][j]) ; 
        
        
        fhNCellNLocMaxN[i][j]     = new TH2F(Form("hNCellNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                             Form("#lambda_{0}^{2} vs E for N max  > 2 %s %s",ptype[i].Data(),sMatched[j].Data()),
                                             nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
        fhNCellNLocMaxN[i][j]   ->SetYTitle("N cells");
        fhNCellNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhNCellNLocMaxN[i][j]) ;
      }
      
      fhM02Pi0LocMax1[i][j]     = new TH2F(Form("hM02Pi0LocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMax1[i][j]) ; 
      
      fhM02EtaLocMax1[i][j]     = new TH2F(Form("hM02EtaLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMax1[i][j]) ; 
      
      fhM02ConLocMax1[i][j]    = new TH2F(Form("hM02ConLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMax1[i][j]) ; 
      
      fhM02Pi0LocMax2[i][j]     = new TH2F(Form("hM02Pi0LocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMax2[i][j]) ; 
      
      fhM02EtaLocMax2[i][j]     = new TH2F(Form("hM02EtaLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMax2[i][j]) ;
      
      fhM02ConLocMax2[i][j]    = new TH2F(Form("hM02ConLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMax2[i][j]) ; 
      
      fhM02Pi0LocMaxN[i][j]     = new TH2F(Form("hM02Pi0LocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMaxN[i][j]) ; 
      
      fhM02EtaLocMaxN[i][j]     = new TH2F(Form("hM02EtaLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max > 2", 
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMaxN[i][j]) ; 
      
      fhM02ConLocMaxN[i][j]    = new TH2F(Form("hM02ConLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMaxN[i][j]) ;
            
      
      fhMassPi0LocMax1[i][j]     = new TH2F(Form("hMassPi0LocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0LocMax1[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassPi0LocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0LocMax1[i][j]) ; 

      
      fhMassEtaLocMax1[i][j]     = new TH2F(Form("hMassEtaLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaLocMax1[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassEtaLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaLocMax1[i][j]) ; 
      
      fhMassConLocMax1[i][j]    = new TH2F(Form("hMassConLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConLocMax1[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassConLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConLocMax1[i][j]) ; 
      
      fhMassPi0LocMax2[i][j]     = new TH2F(Form("hMassPi0LocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0LocMax2[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassPi0LocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0LocMax2[i][j]) ; 
      
      
      fhMassEtaLocMax2[i][j]     = new TH2F(Form("hMassEtaLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaLocMax2[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassEtaLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaLocMax2[i][j]) ; 
      
      fhMassConLocMax2[i][j]    = new TH2F(Form("hMassConLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConLocMax2[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassConLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConLocMax2[i][j]) ; 
      
      fhMassPi0LocMaxN[i][j]     = new TH2F(Form("hMassPi0LocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0LocMaxN[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassPi0LocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0LocMaxN[i][j]) ; 
      
      fhMassEtaLocMaxN[i][j]     = new TH2F(Form("hMassEtaLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max > 2", 
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaLocMaxN[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassEtaLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaLocMaxN[i][j]) ; 
      
      fhMassConLocMaxN[i][j]    = new TH2F(Form("hMassConLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConLocMaxN[i][j]   ->SetYTitle("Mass (MeV/c^{2})");
      fhMassConLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConLocMaxN[i][j]) ; 
      
      
      fhAsyPi0LocMax1[i][j]     = new TH2F(Form("hAsyPi0LocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 1",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0LocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0LocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0LocMax1[i][j]) ; 
      
      fhAsyEtaLocMax1[i][j]     = new TH2F(Form("hAsyEtaLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaLocMax1[i][j]) ; 
      
      fhAsyConLocMax1[i][j]    = new TH2F(Form("hAsyConLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConLocMax1[i][j]) ; 
      
      fhAsyPi0LocMax2[i][j]     = new TH2F(Form("hAsyPi0LocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0LocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0LocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0LocMax2[i][j]) ; 
      
      fhAsyEtaLocMax2[i][j]     = new TH2F(Form("hAsyEtaLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaLocMax2[i][j]) ; 
      
      fhAsyConLocMax2[i][j]    = new TH2F(Form("hAsyConLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConLocMax2[i][j]) ; 
      
      fhAsyPi0LocMaxN[i][j]     = new TH2F(Form("hAsyPi0LocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max > 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0LocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0LocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0LocMaxN[i][j]) ; 
      
      fhAsyEtaLocMaxN[i][j]     = new TH2F(Form("hAsyEtaLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max > 2", 
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaLocMaxN[i][j]) ; 
      
      fhAsyConLocMaxN[i][j]    = new TH2F(Form("hAsyConLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConLocMaxN[i][j]) ; 
      
    } // matched, not matched
    
      for(Int_t j = 0; j < 4; j++)
      {  
        
        fhMassSplitEFractionNLocMax1Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax1%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 highest energy cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax); 
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax1Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax1Ebin[i][j]) ;   
        
        fhMassSplitEFractionNLocMax2Ebin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMax2%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of 2 local maxima cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax); 
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMax2Ebin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMax2Ebin[i][j]) ;   
        
        fhMassSplitEFractionNLocMaxNEbin[i][j]  = new TH2F(Form("hMassSplitEFractionNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                           Form("Invariant mass of N>2 local maxima cells vs (E1+E2)/Ecluster, %s, E bin %d",ptype[i].Data(),j),
                                                           120,0,1.2,mbins,mmin,mmax); 
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassSplitEFractionNLocMaxNEbin[i][j]->SetXTitle("(E_{split1}+E_{split2})/E_{cluster}");
        outputContainer->Add(fhMassSplitEFractionNLocMaxNEbin[i][j]) ;   
        
        if(i>0 && fFillMCFractionHisto) // skip first entry in array, general case not filled
        {
          fhMCGenFracNLocMaxEbin[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%d",pname[i].Data(),j),
                                                   Form("NLM vs E, %s, E bin %d",ptype[i].Data(),j),
                                                   200,0,2,nMaxBins,0,nMaxBins); 
          fhMCGenFracNLocMaxEbin[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbin[i][j]) ;           
          
          fhMCGenFracNLocMaxEbinMatched[i][j]  = new TH2F(Form("hMCGenFracNLocMax%sEbin%dMatched",pname[i].Data(),j),
                                                          Form("NLM vs E, %s, E bin %d, matched to a track",ptype[i].Data(),j),
                                                          200,0,2,nMaxBins,0,nMaxBins); 
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetYTitle("NLM");
          fhMCGenFracNLocMaxEbinMatched[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMCGenFracNLocMaxEbinMatched[i][j]) ;   
          
          fhMassMCGenFracNLocMax1Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 highest energy cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax); 
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax1Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax1Ebin[i][j]) ;   
          
          fhMassMCGenFracNLocMax2Ebin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of 2 local maxima cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax); 
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMax2Ebin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMax2Ebin[i][j]) ;   
          
          fhMassMCGenFracNLocMaxNEbin[i][j]  = new TH2F(Form("hMassMCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                        Form("Invariant mass of N>2 local maxima cells vs E, %s, E bin %d",ptype[i].Data(),j),
                                                        200,0,2,mbins,mmin,mmax); 
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetYTitle("M (GeV/c^{2})");
          fhMassMCGenFracNLocMaxNEbin[i][j]->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhMassMCGenFracNLocMaxNEbin[i][j]) ;   
          
          fhM02MCGenFracNLocMax1Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax1%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 1 %s, E bin %d",ptype[i].Data(), j),
                                                          200,0,2,ssbins,ssmin,ssmax); 
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax1Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax1Ebin[i][j]) ; 
          
          fhM02MCGenFracNLocMax2Ebin[i][j]     = new TH2F(Form("hM02MCGenFracNLocMax2%sEbin%d",pname[i].Data(),j),
                                                          Form("#lambda_{0}^{2} vs E for N max  = 2 %s, E bin %d",ptype[i].Data(),j),
                                                          200,0,2,ssbins,ssmin,ssmax); 
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMax2Ebin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMax2Ebin[i][j]) ; 
          
          fhM02MCGenFracNLocMaxNEbin[i][j]    = new TH2F(Form("hM02MCGenFracNLocMaxN%sEbin%d",pname[i].Data(),j),
                                                         Form("#lambda_{0}^{2} vs E for N max  > 2 %s, E bin %d",ptype[i].Data(),j),
                                                         200,0,2,ssbins,ssmin,ssmax); 
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetYTitle("#lambda_{0}^{2}");
          fhM02MCGenFracNLocMaxNEbin[i][j]   ->SetXTitle("E_{gen} / E_{reco}");
          outputContainer->Add(fhM02MCGenFracNLocMaxNEbin[i][j]) ; 
        }
      }
  } // MC particle list
 
  for(Int_t i = 0; i < 4; i++)
  {  
    fhMassM02NLocMax1Ebin[i]  = new TH2F(Form("hMassM02NLocMax1Ebin%d",i),
                                        Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=1, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMax1Ebin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMax1Ebin[i]) ;   
    
    fhMassM02NLocMax2Ebin[i]  = new TH2F(Form("hMassM02NLocMax2Ebin%d",i),
                                        Form("Invariant mass of split clusters vs #lambda_{0}^{2}, NLM=2, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMax2Ebin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMax2Ebin[i]) ;   
    
    fhMassM02NLocMaxNEbin[i]  = new TH2F(Form("hMassM02NLocMaxNEbin%d",i),
                                        Form("Invariant mass of split clusters vs vs #lambda_{0}^{2}, NLM>2, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMaxNEbin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMaxNEbin[i]) ; 
    
    
    fhMassAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassAsyNLocMax1Ebin%d",i),
                                         Form("Invariant mass of split clusters vs split asymmetry, NLM=1, E bin %d",i),
                                         200,-1,1,mbins,mmin,mmax);
    fhMassAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassAsyNLocMax1Ebin[i]->SetXTitle("asymmetry");
    outputContainer->Add(fhMassAsyNLocMax1Ebin[i]) ;
    
    fhMassAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassAsyNLocMax2Ebin%d",i),
                                         Form("Invariant mass of split clusters vs split asymmetry, NLM=2, E bin %d",i),
                                         200,-1,1,mbins,mmin,mmax);
    fhMassAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassAsyNLocMax2Ebin[i]->SetXTitle("asymmetry");
    outputContainer->Add(fhMassAsyNLocMax2Ebin[i]) ;
    
    fhMassAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassAsyNLocMaxNEbin%d",i),
                                         Form("Invariant mass of split clusters vs split asymmetry, NLM>2, E bin %d",i),
                                         200,-1,1,mbins,mmin,mmax);
    fhMassAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassAsyNLocMaxNEbin[i]->SetXTitle("asymmetry");
    outputContainer->Add(fhMassAsyNLocMaxNEbin[i]) ;

    
    if(IsDataMC())
    {
      fhMCAsymM02NLocMax1MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax1MCPi0Ebin%d",i),
                                                  Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=1, E bin %d",i),
                                                  ssbins,ssmin,ssmax,100,0,1);
      fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
      fhMCAsymM02NLocMax1MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMCAsymM02NLocMax1MCPi0Ebin[i]) ;
      
      fhMCAsymM02NLocMax2MCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMax2MCPi0Ebin%d",i),
                                                  Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM=2, E bin %d",i),
                                                  ssbins,ssmin,ssmax,100,0,1);
      fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetYTitle("Decay asymmetry");
      fhMCAsymM02NLocMax2MCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMCAsymM02NLocMax2MCPi0Ebin[i]) ;
      
      fhMCAsymM02NLocMaxNMCPi0Ebin[i]  = new TH2F(Form("hMCAsymM02NLocMaxNMCPi0Ebin%d",i),
                                                  Form("Asymmetry of MC #pi^{0} vs #lambda_{0}^{2}, NLM>2, E bin %d",i),
                                                  ssbins,ssmin,ssmax,100,0,1);
      fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetYTitle("Decay asymmetry");
      fhMCAsymM02NLocMaxNMCPi0Ebin[i]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMCAsymM02NLocMaxNMCPi0Ebin[i]) ;    
      
      
      fhAsyMCGenRecoNLocMax1EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax1Ebin%dPi0",i),
                                                Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=1, E bin %d",i),
                                                200,-1,1,200,-1,1);
      fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
      fhAsyMCGenRecoNLocMax1EbinPi0[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhAsyMCGenRecoNLocMax1EbinPi0[i]) ;
      
      fhAsyMCGenRecoNLocMax2EbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMax2Ebin%dPi0",i),
                                                Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM=2, E bin %d",i),
                                                200,-1,1,200,-1,1);
      fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetYTitle("M (GeV/c^{2})");
      fhAsyMCGenRecoNLocMax2EbinPi0[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhAsyMCGenRecoNLocMax2EbinPi0[i]) ;
      
      fhAsyMCGenRecoNLocMaxNEbinPi0[i]  = new TH2F(Form("hAsyMCGenRecoNLocMaxNEbin%dPi0",i),
                                                Form("Generated vs reconstructed asymmetry of split clusters from pi0, NLM>2, E bin %d",i),
                                                200,-1,1,200,-1,1);
      fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetYTitle("M (GeV/c^{2})");
      fhAsyMCGenRecoNLocMaxNEbinPi0[i]->SetXTitle("asymmetry");
      outputContainer->Add(fhAsyMCGenRecoNLocMaxNEbinPi0[i]) ;
    }
    
    if(fFillSSExtraHisto)
    {
      fhMassDispEtaNLocMax1Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax1Ebin%d",i),
                                               Form("Invariant mass of 2 highest energy cells #sigma_{#eta #eta}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispEtaNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispEtaNLocMax1Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
      outputContainer->Add(fhMassDispEtaNLocMax1Ebin[i]) ;   
      
      fhMassDispEtaNLocMax2Ebin[i]  = new TH2F(Form("hMassDispEtaNLocMax2Ebin%d",i),
                                               Form("Invariant mass of 2 local maxima cells #sigma_{#eta #eta}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispEtaNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispEtaNLocMax2Ebin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
      outputContainer->Add(fhMassDispEtaNLocMax2Ebin[i]) ;   
      
      fhMassDispEtaNLocMaxNEbin[i]  = new TH2F(Form("hMassDispEtaNLocMaxNEbin%d",i),
                                               Form("Invariant mass of N>2 local maxima cells vs #sigma_{#eta #eta}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispEtaNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispEtaNLocMaxNEbin[i]->SetXTitle("#sigma_{#eta #eta}^{2}");
      outputContainer->Add(fhMassDispEtaNLocMaxNEbin[i]) ;   
      
      fhMassDispPhiNLocMax1Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax1Ebin%d",i),
                                               Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispPhiNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispPhiNLocMax1Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
      outputContainer->Add(fhMassDispPhiNLocMax1Ebin[i]) ;   
      
      fhMassDispPhiNLocMax2Ebin[i]  = new TH2F(Form("hMassDispPhiNLocMax2Ebin%d",i),
                                               Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispPhiNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispPhiNLocMax2Ebin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
      outputContainer->Add(fhMassDispPhiNLocMax2Ebin[i]) ;   
      
      fhMassDispPhiNLocMaxNEbin[i]  = new TH2F(Form("hMassDispPhiNLocMaxNEbin%d",i),
                                               Form("Invariant mass of N>2 local maxima cells vs #sigma_{#phi #phi}^{2}, E bin %d",i),
                                               ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassDispPhiNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispPhiNLocMaxNEbin[i]->SetXTitle("#sigma_{#phi #phi}^{2}");
      outputContainer->Add(fhMassDispPhiNLocMaxNEbin[i]) ;   
      
      fhMassDispAsyNLocMax1Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax1Ebin%d",i),
                                               Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                               200,-1,1,mbins,mmin,mmax); 
      fhMassDispAsyNLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispAsyNLocMax1Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
      outputContainer->Add(fhMassDispAsyNLocMax1Ebin[i]) ;   
      
      fhMassDispAsyNLocMax2Ebin[i]  = new TH2F(Form("hMassDispAsyNLocMax2Ebin%d",i),
                                               Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                               200,-1,1,mbins,mmin,mmax); 
      fhMassDispAsyNLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispAsyNLocMax2Ebin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
      outputContainer->Add(fhMassDispAsyNLocMax2Ebin[i]) ;   
      
      fhMassDispAsyNLocMaxNEbin[i]  = new TH2F(Form("hMassDispAsyNLocMaxNEbin%d",i),
                                               Form("Invariant mass of N>2 local maxima cells vs A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E bin %d",i),
                                               200,-1,1,mbins,mmin,mmax); 
      fhMassDispAsyNLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
      fhMassDispAsyNLocMaxNEbin[i]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
      outputContainer->Add(fhMassDispAsyNLocMaxNEbin[i]) ;   
    }
  }  
  
  if(fFillTMResidualHisto)
  {
    for(Int_t i = 0; i < n; i++)
    {  
      
      fhTrackMatchedDEtaLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDEtaLocMax1%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaLocMax1[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDPhiLocMax1%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiLocMax1[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaLocMax1[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiLocMax1[i]) ;
      
      fhTrackMatchedDEtaLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDEtaLocMax2%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaLocMax2[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDPhiLocMax2%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiLocMax2[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaLocMax2[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiLocMax2[i]) ;
      
      fhTrackMatchedDEtaLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDEtaLocMaxN%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaLocMaxN[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDPhiLocMaxN%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiLocMaxN[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaLocMaxN[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiLocMaxN[i]) ;    
    }
  }
  
  if(fFillAngleHisto)
  {
    for(Int_t j = 0; j < 2; j++)
    {  
      
      fhAnglePairLocMax1[j]  = new TH2F(Form("hAnglePairLocMax1%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 highest energy cells vs pair Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairLocMax1[j]->SetYTitle("#alpha (rad)");
      fhAnglePairLocMax1[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairLocMax1[j]) ;   
      
      fhAnglePairLocMax2[j]  = new TH2F(Form("hAnglePairLocMax2%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairLocMax2[j]->SetYTitle("#alpha (rad)");
      fhAnglePairLocMax2[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairLocMax2[j]) ;   
      
      fhAnglePairLocMaxN[j]  = new TH2F(Form("hAnglePairLocMaxN%s",sMatched[j].Data()),
                                        Form("Opening angle of N>2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairLocMaxN[j]->SetYTitle("#alpha (rad)");
      fhAnglePairLocMaxN[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairLocMaxN[j]) ;   
      
      fhAnglePairMassLocMax1[j]  = new TH2F(Form("hAnglePairMassLocMax1%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 highest energy cells vs Mass for E > 8 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassLocMax1[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassLocMax1[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassLocMax1[j]) ;   
      
      fhAnglePairMassLocMax2[j]  = new TH2F(Form("hAnglePairMassLocMax2%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 local maxima cells vs Mass for E > 8 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassLocMax2[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassLocMax2[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassLocMax2[j]) ;   
      
      fhAnglePairMassLocMaxN[j]  = new TH2F(Form("hAnglePairMassLocMaxN%s",sMatched[j].Data()),
                                            Form("Opening angle of N>2 local maxima cells vs Mass for E > 8 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassLocMaxN[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassLocMaxN[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassLocMaxN[j]) ;  
      
    }
  }
  
  for(Int_t j = 0; j < 2; j++)
  {  
  fhSplitEFractionvsAsyNLocMax1[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax1%s",sMatched[j].Data()),
                                                Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 1, E>8, %s",sMatched[j].Data()),
                                                100,-1,1,120,0,1.2); 
  fhSplitEFractionvsAsyNLocMax1[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
  fhSplitEFractionvsAsyNLocMax1[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionvsAsyNLocMax1[j]) ; 
  
  fhSplitEFractionvsAsyNLocMax2[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax2%s",sMatched[j].Data()),
                                                Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 2,E>8, %s",sMatched[j].Data()),
                                                100,-1,1,120,0,1.2); 
  fhSplitEFractionvsAsyNLocMax2[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
  fhSplitEFractionvsAsyNLocMax2[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionvsAsyNLocMax2[j]) ; 
  
  fhSplitEFractionvsAsyNLocMaxN[j]    = new TH2F(Form("hSplitEFractionvsAsyNLocMaxN%s",sMatched[j].Data()),
                                               Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  > 2, E>8, %s",sMatched[j].Data()),
                                               100,-1,1,120,0,1.2); 
  fhSplitEFractionvsAsyNLocMaxN[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
  fhSplitEFractionvsAsyNLocMaxN[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
  outputContainer->Add(fhSplitEFractionvsAsyNLocMaxN[j]) ; 
  }
   
  
  return outputContainer ;
  
}

//___________________________________________
void AliAnaInsideClusterInvariantMass::Init()
{ 
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD())
  {
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD())
  {
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kMC )
  {
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use pure MC data!!\n");
    abort();
    
  }
  
}

//_____________________________________________________
void AliAnaInsideClusterInvariantMass::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0InsideClusterInvariantMass_");
  
  fCalorimeter = "EMCAL" ;

  fM02MinCut   = 0.26 ;
  fM02MaxCut   = 10 ;
  
  fMinNCells   = 4 ;
  fMinBadDist  = 2 ;
    
}


//__________________________________________________________________
void  AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl       = 0x0; 
  AliVCaloCells* cells = 0x0;

  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS")
  {
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL")
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  const Float_t ecut = 8.; // Fixed cut for some histograms
  
  if(!pl || !cells) 
  {
    Info("MakeAnalysisFillHistograms","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }  
  
	if(fCalorimeter == "PHOS") return; // Not implemented for PHOS yet

  for(Int_t icluster = 0; icluster < pl->GetEntriesFast(); icluster++)
  {
    AliVCluster * cluster = (AliVCluster*) (pl->At(icluster));	

    // Study clusters with large shape parameter
    Float_t en = cluster->E();
    Float_t l0 = cluster->GetM02();
    Int_t   nc = cluster->GetNCells();
    Float_t bd = cluster->GetDistanceToBadChannel() ; 

    
    //If too small or big E or low number of cells, or close to a bad channel skip it
    if( en < GetMinEnergy() || en > GetMaxEnergy() || nc < fMinNCells || bd < fMinBadDist) continue ; 
    
    //printf("en %2.2f, GetMinEnergy() %2.2f, GetMaxEnergy() %2.2f, nc %d, fMinNCells %d,  bd %2.2f, fMinBadDist %2.2f\n",
    //       en,GetMinEnergy(), GetMaxEnergy(), nc, fMinNCells, bd, fMinBadDist);
    
    // Get more Shower Shape parameters
    Float_t ll0  = 0., ll1  = 0.;
    Float_t disp= 0., dispEta = 0., dispPhi    = 0.; 
    Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;  
   
    GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), GetReader()->GetInputEvent()->GetEMCALCells(), cluster,
                                                                                 ll0, ll1, disp, dispEta, dispPhi, sEta, sPhi, sEtaPhi);
    
    Float_t dispAsy = -1;
    if(dispEta+dispPhi >0 ) dispAsy = (dispPhi-dispEta) / (dispPhi+dispEta);
    
    Int_t    nMax = 0;
    Double_t mass = 0., angle = 0.;
    Double_t e1   = 0., e2    = 0.;
    Int_t pidTag = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(cluster,cells,GetCaloUtils(),
                                                                               GetVertex(0), nMax, mass, angle,e1,e2);    
    if (nMax <= 0) 
    {
      if(GetDebug() > 0 )
        printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found! It did not pass CaloPID selection criteria \n");
      
      return;
    }
    
    Float_t splitFrac = (e1+e2)/en;
    Float_t asym = -10;
    if(e1+e2>0) asym = (e1-e2)/(e1+e2);
    
    Bool_t  matched   = IsTrackMatched(cluster,GetReader()->GetInputEvent());
    
    fhNLocMax[0][matched]->Fill(en,nMax);
    
    if     ( nMax == 1  ) 
    { 
      fhM02NLocMax1[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax1[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut) fhSplitEFractionvsAsyNLocMax1[matched]->Fill(asym,splitFrac) ; 
      if(fFillSSExtraHisto) fhNCellNLocMax1[0][matched]->Fill(en,nc) ; 
    }
    else if( nMax == 2  ) 
    { 
      fhM02NLocMax2[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax2[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut) fhSplitEFractionvsAsyNLocMax2[matched]->Fill(asym,splitFrac) ; 
      if(fFillSSExtraHisto) fhNCellNLocMax2[0][matched]->Fill(en,nc) ; }
    else if( nMax >= 3  ) 
    { 
      fhM02NLocMaxN[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMaxN[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut) fhSplitEFractionvsAsyNLocMaxN[matched]->Fill(asym,splitFrac) ; 
      if(fFillSSExtraHisto) fhNCellNLocMaxN[0][matched]->Fill(en,nc) ; 
    }
    else printf("N max smaller than 1 -> %d \n",nMax);
    
    
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
    {
      dR = 2000., dZ = 2000.;
      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }    
    //printf("Pi0EbE: dPhi %f, dEta %f\n",dR,dZ);
    
    if(TMath::Abs(dR) < 999 && fFillTMResidualHisto)
    {
      if     ( nMax == 1  ) { fhTrackMatchedDEtaLocMax1[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMax1[0]->Fill(en,dR); }
      else if( nMax == 2  ) { fhTrackMatchedDEtaLocMax2[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMax2[0]->Fill(en,dR); }
      else if( nMax >= 3  ) { fhTrackMatchedDEtaLocMaxN[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMaxN[0]->Fill(en,dR); }
    }
    
    // Play with the MC stack if available
    // Check origin of the candidates
    Int_t mcindex   = -1;
    Float_t eprim   =  0;
    Float_t asymGen = -2; 
    Int_t mcLabel   = cluster->GetLabel();
    if(IsDataMC())
    {
      Int_t tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader());
            
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = kmcEta;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPhoton;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
                GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcConversion;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))   mcindex = kmcElectron;
      else                                                                                mcindex = kmcHadron;

      fhNLocMax[mcindex][matched]->Fill(en,nMax);
            
      if     (nMax == 1 ) { fhM02NLocMax1[mcindex][matched]->Fill(en,l0) ; fhSplitEFractionNLocMax1[mcindex][matched]->Fill(en,splitFrac) ; if(fFillSSExtraHisto) fhNCellNLocMax1[mcindex][matched]->Fill(en,nc) ; }
      else if(nMax == 2 ) { fhM02NLocMax2[mcindex][matched]->Fill(en,l0) ; fhSplitEFractionNLocMax2[mcindex][matched]->Fill(en,splitFrac) ; if(fFillSSExtraHisto) fhNCellNLocMax2[mcindex][matched]->Fill(en,nc) ; }
      else if(nMax >= 3 ) { fhM02NLocMaxN[mcindex][matched]->Fill(en,l0) ; fhSplitEFractionNLocMaxN[mcindex][matched]->Fill(en,splitFrac) ; if(fFillSSExtraHisto) fhNCellNLocMaxN[mcindex][matched]->Fill(en,nc) ; }
      
      if(TMath::Abs(dR) < 999 && fFillTMResidualHisto)
      {
        if     ( nMax == 1  ) { fhTrackMatchedDEtaLocMax1[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMax1[mcindex]->Fill(en,dR); }
        else if( nMax == 2  ) { fhTrackMatchedDEtaLocMax2[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMax2[mcindex]->Fill(en,dR); }
        else if( nMax >= 3  ) { fhTrackMatchedDEtaLocMaxN[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMaxN[mcindex]->Fill(en,dR); }
      }
      
      Bool_t ok = kFALSE;
      TLorentzVector primary = GetMCAnalysisUtils()->GetMother(mcLabel,GetReader(),ok);
      eprim = primary.E();
      
      if(mcindex == kmcPi0 || mcindex == kmcEta)
      {
        if(mcindex == kmcPi0) 
        {
          asymGen = TMath::Abs(GetMCAnalysisUtils()->GetMCDecayAsymmetryForPDG(mcLabel,111,GetReader(),ok));
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,111,GetReader(),ok); 
          if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
        }
        else 
        {
          asymGen = TMath::Abs(GetMCAnalysisUtils()->GetMCDecayAsymmetryForPDG(mcLabel,221,GetReader(),ok));
          TLorentzVector grandmom = GetMCAnalysisUtils()->GetMotherWithPDG(mcLabel,221,GetReader(),ok); 
          if(grandmom.E() > 0 && ok) eprim =  grandmom.E();
        }
      }
    } 
    
    Float_t efrac      = eprim/en;
    Float_t efracSplit = 0;
    if(e1+e2 > 0) efracSplit = eprim/(e1+e2);

    //printf("e1 %2.2f, e2 %2.2f, eprim %2.2f, ereco %2.2f, esplit/ereco %2.2f, egen/ereco %2.2f, egen/esplit %2.2f\n",
    //       e1,e2,eprim,en,splitFrac,efrac,efracSplit);
    
    Int_t ebin = -1;
    if(en > 8  && en <= 12) ebin = 0; 
    if(en > 12 && en <= 16) ebin = 1;
    if(en > 16 && en <= 20) ebin = 2;
    if(en > 20)             ebin = 3; 
    
    if(ebin >= 0 && IsDataMC() && fFillMCFractionHisto)
    {
      if( !matched ) fhMCGenFracNLocMaxEbin       [mcindex][ebin]->Fill(efrac,nMax);
      else           fhMCGenFracNLocMaxEbinMatched[mcindex][ebin]->Fill(efrac,nMax);
    }
    
    if     (nMax==1) 
    { 
      if( en > ecut ) 
      {      
        fhMassM02NLocMax1    [0][matched]->Fill(l0     ,  mass ); 
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMax1[0][matched]->Fill(dispEta,  mass ); 
          fhMassDispPhiNLocMax1[0][matched]->Fill(dispPhi,  mass ); 
          fhMassDispAsyNLocMax1[0][matched]->Fill(dispAsy,  mass );
        }
        
        if(IsDataMC()) 
        {
          fhMassM02NLocMax1          [mcindex][matched]->Fill(l0     ,  mass  ); 
          if(fFillMCFractionHisto)
          {
            fhMCGenFracNLocMax1        [mcindex][matched]->Fill(en     ,  efrac ); 
            fhMCGenSplitEFracNLocMax1  [mcindex][matched]->Fill(en     ,  efracSplit ); 
            fhMCGenEvsSplitENLocMax1   [mcindex][matched]->Fill(eprim  ,  e1+e2); 
            fhMCGenEFracvsSplitEFracNLocMax1[mcindex][matched]->Fill(efrac,splitFrac ); 
          }
          
          if(!matched && ebin >= 0)
          {
            if(fFillMCFractionHisto)
            {
              fhM02MCGenFracNLocMax1Ebin [mcindex][ebin]->Fill(efrac  ,  l0    ); 
              fhMassMCGenFracNLocMax1Ebin[mcindex][ebin]->Fill(efrac  ,  mass  ); 
            }
            fhMCAsymM02NLocMax1MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
            fhAsyMCGenRecoNLocMax1EbinPi0[ebin]->Fill(asym,  asymGen );
          }
          
          if(fFillSSExtraHisto)
          {
            fhMassDispEtaNLocMax1[mcindex][matched]->Fill(dispEta,  mass ); 
            fhMassDispPhiNLocMax1[mcindex][matched]->Fill(dispPhi,  mass ); 
            fhMassDispAsyNLocMax1[mcindex][matched]->Fill(dispAsy,  mass ); 
          }
        }
      }
      
      if(!matched && ebin >= 0)
      {
        fhMassSplitEFractionNLocMax1Ebin[0][ebin]->Fill(splitFrac,  mass);
        if(IsDataMC())fhMassSplitEFractionNLocMax1Ebin[mcindex][ebin]->Fill(splitFrac,  mass);

        fhMassM02NLocMax1Ebin    [ebin]->Fill(l0  ,  mass );
        fhMassAsyNLocMax1Ebin    [ebin]->Fill(asym,  mass );
        
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMax1Ebin[ebin]->Fill(dispEta,  mass );
          fhMassDispPhiNLocMax1Ebin[ebin]->Fill(dispPhi,  mass );
          fhMassDispAsyNLocMax1Ebin[ebin]->Fill(dispAsy,  mass );
        }
      }
    }  
    else if(nMax==2) 
    {
      if( en > ecut ) 
      {      
        fhMassM02NLocMax2    [0][matched]->Fill(l0     ,  mass ); 
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMax2[0][matched]->Fill(dispEta,  mass ); 
          fhMassDispPhiNLocMax2[0][matched]->Fill(dispPhi,  mass ); 
          fhMassDispAsyNLocMax2[0][matched]->Fill(dispAsy,  mass ); 
        }
        
        if(IsDataMC()) 
        {
          fhMassM02NLocMax2        [mcindex][matched]->Fill(l0     ,  mass ); 
          if(fFillMCFractionHisto)
          {
            fhMCGenFracNLocMax2      [mcindex][matched]->Fill(en     ,  efrac ); 
            fhMCGenSplitEFracNLocMax2[mcindex][matched]->Fill(en     ,  efracSplit ); 
            fhMCGenEvsSplitENLocMax2 [mcindex][matched]->Fill(eprim  ,  e1+e2); 
            fhMCGenEFracvsSplitEFracNLocMax2[mcindex][matched]->Fill(efrac,splitFrac ); 
          }
          
          if(!matched && ebin >= 0)
          {
            if(fFillMCFractionHisto)
            {
              fhM02MCGenFracNLocMax2Ebin [mcindex][ebin]->Fill(efrac  ,  l0    ); 
              fhMassMCGenFracNLocMax2Ebin[mcindex][ebin]->Fill(efrac  ,  mass  ); 
            }
            fhMCAsymM02NLocMax2MCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
            fhAsyMCGenRecoNLocMax2EbinPi0[ebin]->Fill(asym,  asymGen );
          }
          if(fFillSSExtraHisto)
          {
            fhMassDispEtaNLocMax2[mcindex][matched]->Fill(dispEta,  mass ); 
            fhMassDispPhiNLocMax2[mcindex][matched]->Fill(dispPhi,  mass ); 
            fhMassDispAsyNLocMax2[mcindex][matched]->Fill(dispAsy,  mass ); 
          }
        }
      }
      
      if(!matched && ebin >= 0)
      {
        fhMassSplitEFractionNLocMax2Ebin[0][ebin]->Fill(splitFrac,  mass);
        if(IsDataMC())fhMassSplitEFractionNLocMax2Ebin[mcindex][ebin]->Fill(splitFrac,  mass);

        fhMassM02NLocMax2Ebin    [ebin]->Fill(l0  ,  mass );
        fhMassAsyNLocMax2Ebin    [ebin]->Fill(asym,  mass );
        
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMax2Ebin[ebin]->Fill(dispEta,  mass );
          fhMassDispPhiNLocMax2Ebin[ebin]->Fill(dispPhi,  mass );
          fhMassDispAsyNLocMax2Ebin[ebin]->Fill(dispAsy,  mass );
        }
      }   
    }
    else if(nMax > 2 ) 
    {
      if( en > ecut ) 
      {      
        fhMassM02NLocMaxN    [0][matched]->Fill(l0     ,  mass ); 
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMaxN[0][matched]->Fill(dispEta,  mass ); 
          fhMassDispPhiNLocMaxN[0][matched]->Fill(dispPhi,  mass ); 
          fhMassDispAsyNLocMaxN[0][matched]->Fill(dispAsy,  mass ); 
        }
        
        if(IsDataMC()) 
        {
          fhMassM02NLocMaxN        [mcindex][matched]->Fill(l0     ,  mass ); 
          if(fFillMCFractionHisto)
          {
            fhMCGenFracNLocMaxN      [mcindex][matched]->Fill(en     ,  efrac ); 
            fhMCGenSplitEFracNLocMaxN[mcindex][matched]->Fill(en     ,  efracSplit ); 
            fhMCGenEvsSplitENLocMaxN [mcindex][matched]->Fill(eprim  ,  e1+e2); 
            fhMCGenEFracvsSplitEFracNLocMaxN[mcindex][matched]->Fill(efrac,  splitFrac ); 
          }
          
          if(!matched && ebin >= 0)
          {
            if(fFillMCFractionHisto)
            {
              fhM02MCGenFracNLocMaxNEbin [mcindex][ebin]->Fill(efrac  ,  l0     ); 
              fhMassMCGenFracNLocMaxNEbin[mcindex][ebin]->Fill(efrac  ,  mass   ); 
            }
            fhMCAsymM02NLocMaxNMCPi0Ebin [ebin]->Fill(l0  ,  asymGen );
            fhAsyMCGenRecoNLocMaxNEbinPi0[ebin]->Fill(asym,  asymGen );
          }
          if(fFillSSExtraHisto)
          {
            fhMassDispEtaNLocMaxN[mcindex][matched]->Fill(dispEta,  mass ); 
            fhMassDispPhiNLocMaxN[mcindex][matched]->Fill(dispPhi,  mass ); 
            fhMassDispAsyNLocMaxN[mcindex][matched]->Fill(dispAsy,  mass ); 
          }
        }
      }
      
      if(!matched && ebin >= 0)
      {
        fhMassSplitEFractionNLocMaxNEbin[0][ebin]->Fill(splitFrac,  mass);
        if(IsDataMC())fhMassSplitEFractionNLocMaxNEbin[mcindex][ebin]->Fill(splitFrac,  mass);

        fhMassM02NLocMaxNEbin    [ebin]->Fill(l0  ,  mass );
        fhMassAsyNLocMaxNEbin    [ebin]->Fill(asym,  mass );
        
        if(fFillSSExtraHisto)
        {
          fhMassDispEtaNLocMaxNEbin[ebin]->Fill(dispEta,  mass );
          fhMassDispPhiNLocMaxNEbin[ebin]->Fill(dispPhi,  mass );
          fhMassDispAsyNLocMaxNEbin[ebin]->Fill(dispAsy,  mass );
        }
      }
    }
    
    //---------------------------------------------------------------------
    // From here only if M02 is large but not too large, fill histograms 
    //---------------------------------------------------------------------
    
    if( l0 < fM02MinCut || l0 > fM02MaxCut ) continue ;
    
    Bool_t m02OK = GetCaloPID()->IsInPi0M02Range(en,l0,nMax);
    Bool_t asyOK = GetCaloPID()->IsInPi0SplitAsymmetryRange(en,asym,nMax);
    
    fhNLocMaxM02Cut[0][matched]->Fill(en,nMax);
    if(IsDataMC()) fhNLocMaxM02Cut[mcindex][matched]->Fill(en,nMax);
        
    if     (nMax==1) 
    { 
      fhMassNLocMax1[0][matched]->Fill(en,mass ); 
      fhAsymNLocMax1[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms
      if(splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched)
      {
        fhMassSplitECutNLocMax1->Fill(en,mass );
        if(m02OK)
        {
          fhMassM02CutNLocMax1->Fill(en,mass);
          fhAsymM02CutNLocMax1->Fill(en,asym );
          if(asyOK) fhMassAfterCutsNLocMax1[0]->Fill(en,mass);
        } // m02
      } // split frac
      
      if(m02OK && asyOK && !matched)
      {
        fhSplitEFractionAfterCutsNLocMax1->Fill(en,splitFrac);
        if(IsDataMC() && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          fhMCGenFracAfterCutsNLocMax1MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax1MCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto) 
      {
        fhAnglePairLocMax1[matched]->Fill(en,angle);
      if( en > ecut ) 
        fhAnglePairMassLocMax1[matched]->Fill(mass,angle);
      }
      
      if(asyOK && m02OK)
      {
      }
      
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMax1[0][matched]->Fill(en,l0); fhMassConLocMax1[0][matched]->Fill(en,mass);  fhAsyConLocMax1[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMax1[0][matched]->Fill(en,l0); fhMassPi0LocMax1[0][matched]->Fill(en,mass);  fhAsyPi0LocMax1[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kEta)    { fhM02EtaLocMax1[0][matched]->Fill(en,l0); fhMassEtaLocMax1[0][matched]->Fill(en,mass);  fhAsyEtaLocMax1[0][matched]->Fill(en,asym); }
    }
    else if(nMax==2) 
    {
      fhMassNLocMax2[0][matched]->Fill(en,mass );
      fhAsymNLocMax2[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms
      if(splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched)
      {
        fhMassSplitECutNLocMax2->Fill(en,mass );
        if(m02OK)
        {
          fhMassM02CutNLocMax2->Fill(en,mass);
          fhAsymM02CutNLocMax2->Fill(en,asym );
          if(asyOK) fhMassAfterCutsNLocMax2[0]->Fill(en,mass);
        } // m02
      } // split frac
      
      if(m02OK && asyOK && !matched)
      {
        fhSplitEFractionAfterCutsNLocMax2->Fill(en,splitFrac);
        if(IsDataMC()  && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          fhMCGenFracAfterCutsNLocMax2MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax2MCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto) 
      {
        fhAnglePairLocMax2[matched]->Fill(en,angle);
        if( en > ecut ) 
          fhAnglePairMassLocMax2[matched]->Fill(mass,angle);
      }
            
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMax2[0][matched]->Fill(en,l0); fhMassConLocMax2[0][matched]->Fill(en,mass);  fhAsyConLocMax2[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMax2[0][matched]->Fill(en,l0); fhMassPi0LocMax2[0][matched]->Fill(en,mass);  fhAsyPi0LocMax2[0][matched]->Fill(en,asym); }        
      else if(pidTag==AliCaloPID::kEta)    { fhM02EtaLocMax2[0][matched]->Fill(en,l0); fhMassEtaLocMax2[0][matched]->Fill(en,mass);  fhAsyEtaLocMax2[0][matched]->Fill(en,asym); }
    }
    else if(nMax >2) 
    {
      fhMassNLocMaxN[0][matched]->Fill(en,mass);
      fhAsymNLocMaxN[0][matched]->Fill(en,asym);
      
      // Effect of cuts in mass histograms
      if(splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched)
      {
        fhMassSplitECutNLocMaxN->Fill(en,mass );
        if(m02OK)
        {
          fhMassM02CutNLocMaxN->Fill(en,mass);
          fhAsymM02CutNLocMaxN->Fill(en,asym );
          if(asyOK) fhMassAfterCutsNLocMaxN[0]->Fill(en,mass);
        } // m02
      } // split frac
      
      if(m02OK && asyOK && !matched)
      {
        fhSplitEFractionAfterCutsNLocMaxN->Fill(en,splitFrac);
        if(IsDataMC() && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          fhMCGenFracAfterCutsNLocMaxNMCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto)
      {
        fhAnglePairLocMaxN[matched]->Fill(en,angle);
        if( en > ecut ) 
          fhAnglePairMassLocMaxN[matched]->Fill(mass,angle);
      }
            
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMaxN[0][matched]->Fill(en,l0); fhMassConLocMaxN[0][matched]->Fill(en,mass);  fhAsyConLocMaxN[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMaxN[0][matched]->Fill(en,l0); fhMassPi0LocMaxN[0][matched]->Fill(en,mass);  fhAsyPi0LocMaxN[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kEta)    { fhM02EtaLocMaxN[0][matched]->Fill(en,l0); fhMassEtaLocMaxN[0][matched]->Fill(en,mass);  fhAsyEtaLocMaxN[0][matched]->Fill(en,asym); } 
    }
    
    
    if(IsDataMC())
    {
      if     (nMax==1) 
      { 
        fhMassNLocMax1[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax1[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK && splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched) fhMassAfterCutsNLocMax1[mcindex]->Fill(en,mass);

        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMax1[mcindex][matched]->Fill(en,l0); fhMassConLocMax1[mcindex][matched]->Fill(en,mass); fhAsyConLocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMax1[mcindex][matched]->Fill(en,l0); fhMassPi0LocMax1[mcindex][matched]->Fill(en,mass); fhAsyPi0LocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaLocMax1[mcindex][matched]->Fill(en,l0); fhMassEtaLocMax1[mcindex][matched]->Fill(en,mass); fhAsyEtaLocMax1[mcindex][matched]->Fill(en,asym); } 
      }  
      else if(nMax==2) 
      {
        fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax2[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK && splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched) fhMassAfterCutsNLocMax2[mcindex]->Fill(en,mass);
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMax2[mcindex][matched]->Fill(en,l0); fhMassConLocMax2[mcindex][matched]->Fill(en,mass); fhAsyConLocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMax2[mcindex][matched]->Fill(en,l0); fhMassPi0LocMax2[mcindex][matched]->Fill(en,mass); fhAsyPi0LocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaLocMax2[mcindex][matched]->Fill(en,l0); fhMassEtaLocMax2[mcindex][matched]->Fill(en,mass); fhAsyEtaLocMax2[mcindex][matched]->Fill(en,asym); } 
        
      }
      else if(nMax >2) 
      {
        fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMaxN[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK && splitFrac > GetCaloPID()->GetSplitEnergyFractionMinimum() && !matched) fhMassAfterCutsNLocMaxN[mcindex]->Fill(en,mass);
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConLocMaxN[mcindex][matched]->Fill(en,l0); fhMassConLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyConLocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0LocMaxN[mcindex][matched]->Fill(en,l0); fhMassPi0LocMaxN[mcindex][matched]->Fill(en,mass); fhAsyPi0LocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaLocMaxN[mcindex][matched]->Fill(en,l0); fhMassEtaLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyEtaLocMaxN[mcindex][matched]->Fill(en,asym); } 
      }
      
    }//Work with MC truth first
    
  }//loop
  
  if(GetDebug() > 1) printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - END \n");  

}

//______________________________________________________________________
void AliAnaInsideClusterInvariantMass::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print("");
  printf("Calorimeter     =     %s\n",  fCalorimeter.Data()) ;
  printf("Loc. Max. E > %2.2f\n",       GetCaloUtils()->GetLocalMaximaCutE());
  printf("Loc. Max. E Diff > %2.2f\n",  GetCaloUtils()->GetLocalMaximaCutEDiff());
  printf("Min. N Cells =%d \n",         fMinNCells) ;
  printf("Min. Dist. to Bad =%1.1f \n", fMinBadDist) ;
  printf("%2.2f < lambda_0^2 <%2.2f \n",fM02MinCut,fM02MaxCut);
 
  printf("    \n") ;
  
} 



