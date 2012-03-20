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
  fMinNCells(0),    fMinBadDist(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 7; i++)
  {
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
      
      fhMassM02NLocMax1[i][j]= 0;
      fhMassM02NLocMax2[i][j]= 0;
      fhMassM02NLocMaxN[i][j]= 0;      
      
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
  }
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhMassM02NLocMax1Ebin[i] = 0 ;
    fhMassM02NLocMax2Ebin[i] = 0 ;
    fhMassM02NLocMaxNEbin[i] = 0 ;
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
  
  for(Int_t i = 0; i < n; i++)
  {  
    
    for(Int_t j = 0; j < 2; j++)
    {  
      
      fhMassNLocMax1[i][j]  = new TH2F(Form("hMassNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of 2 highest energy cells vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax1[i][j]) ;   
      
      fhMassNLocMax2[i][j]  = new TH2F(Form("hMassNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of 2 local maxima cells vs E,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax2[i][j]) ;   
      
      fhMassNLocMaxN[i][j]  = new TH2F(Form("hMassNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of N>2 local maxima cells vs E, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMaxN[i][j]) ;   
      
      fhMassM02NLocMax1[i][j]  = new TH2F(Form("hMassM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of 2 highest energy cells #lambda_{0}^{2}, E > 7 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1[i][j]) ;   
      
      fhMassM02NLocMax2[i][j]  = new TH2F(Form("hMassM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of 2 local maxima cells #lambda_{0}^{2}, E > 7 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax2[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax2[i][j]) ;   
      
      fhMassM02NLocMaxN[i][j]  = new TH2F(Form("hMassM02NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of N>2 local maxima cells vs #lambda_{0}^{2}, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMaxN[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMaxN[i][j]) ;   
      
      
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
      
    } // matched, not matched
    
    
  } // MC particle list
 
  for(Int_t i = 0; i < 4; i++)
  {  
    fhMassM02NLocMax1Ebin[i]  = new TH2F(Form("hMassM02NLocMax1Ebin%d",i),
                                        Form("Invariant mass of 2 highest energy cells #lambda_{0}^{2}, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMax1Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMax1Ebin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMax1Ebin[i]) ;   
    
    fhMassM02NLocMax2Ebin[i]  = new TH2F(Form("hMassM02NLocMax2Ebin%d",i),
                                        Form("Invariant mass of 2 local maxima cells #lambda_{0}^{2}, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMax2Ebin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMax2Ebin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMax2Ebin[i]) ;   
    
    fhMassM02NLocMaxNEbin[i]  = new TH2F(Form("hMassM02NLocMaxNEbin%d",i),
                                        Form("Invariant mass of N>2 local maxima cells vs #lambda_{0}^{2}, E bin %d",i),
                                        ssbins,ssmin,ssmax,mbins,mmin,mmax); 
    fhMassM02NLocMaxNEbin[i]->SetYTitle("M (GeV/c^{2})");
    fhMassM02NLocMaxNEbin[i]->SetXTitle("#lambda_{0}^{2}");
    outputContainer->Add(fhMassM02NLocMaxNEbin[i]) ;   
  }  
  
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
                                          Form("Opening angle of 2 highest energy cells vs Mass for E > 7 GeV, %s",sMatched[j].Data()),
                                          mbins,mmin,mmax,200,0,0.2); 
    fhAnglePairMassLocMax1[j]->SetXTitle("M (GeV/c^{2})");
    fhAnglePairMassLocMax1[j]->SetYTitle("#alpha (rad)");
    outputContainer->Add(fhAnglePairMassLocMax1[j]) ;   
    
    fhAnglePairMassLocMax2[j]  = new TH2F(Form("hAnglePairMassLocMax2%s",sMatched[j].Data()),
                                          Form("Opening angle of 2 local maxima cells vs Mass for E > 7 GeV, %s",sMatched[j].Data()),
                                          mbins,mmin,mmax,200,0,0.2); 
    fhAnglePairMassLocMax2[j]->SetXTitle("M (GeV/c^{2})");
    fhAnglePairMassLocMax2[j]->SetYTitle("#alpha (rad)");
    outputContainer->Add(fhAnglePairMassLocMax2[j]) ;   
    
    fhAnglePairMassLocMaxN[j]  = new TH2F(Form("hAnglePairMassLocMaxN%s",sMatched[j].Data()),
                                          Form("Opening angle of N>2 local maxima cells vs Mass for E > 7 GeV, %s",sMatched[j].Data()),
                                          mbins,mmin,mmax,200,0,0.2); 
    fhAnglePairMassLocMaxN[j]->SetXTitle("M (GeV/c^{2})");
    fhAnglePairMassLocMaxN[j]->SetYTitle("#alpha (rad)");
    outputContainer->Add(fhAnglePairMassLocMaxN[j]) ;  
    
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
    
   
    Int_t    nMax = 0;
    Double_t mass = 0;
    Double_t angle= 0;
    Int_t pidTag = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(cluster,cells,GetCaloUtils(),
                                                                               GetVertex(0), nMax, mass, angle);
                                                                    
    if (nMax <= 0) 
    {
      if(GetDebug() > 0 )
        printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found! It did not pass CaloPID selection criteria \n");
      
      return;
    }
    
    Bool_t   matched   = IsTrackMatched(cluster,GetReader()->GetInputEvent());
    
    fhNLocMax[0][matched]->Fill(en,nMax);
    
    if     ( nMax == 1  ) { fhM02NLocMax1[0][matched]->Fill(en,l0) ; fhNCellNLocMax1[0][matched]->Fill(en,nc) ; }
    else if( nMax == 2  ) { fhM02NLocMax2[0][matched]->Fill(en,l0) ; fhNCellNLocMax2[0][matched]->Fill(en,nc) ; }
    else if( nMax >= 3  ) { fhM02NLocMaxN[0][matched]->Fill(en,l0) ; fhNCellNLocMaxN[0][matched]->Fill(en,nc) ; }
    else printf("N max smaller than 1 -> %d \n",nMax);
    
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
    if(cluster->IsEMCAL() && GetCaloUtils()->IsRecalculationOfClusterTrackMatchingOn())
    {
      dR = 2000., dZ = 2000.;
      GetCaloUtils()->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }    
    //printf("Pi0EbE: dPhi %f, dEta %f\n",dR,dZ);
    
    if(TMath::Abs(dR) < 999)
    {
      if     ( nMax == 1  ) { fhTrackMatchedDEtaLocMax1[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMax1[0]->Fill(en,dR); }
      else if( nMax == 2  ) { fhTrackMatchedDEtaLocMax2[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMax2[0]->Fill(en,dR); }
      else if( nMax >= 3  ) { fhTrackMatchedDEtaLocMaxN[0]->Fill(en,dZ); fhTrackMatchedDPhiLocMaxN[0]->Fill(en,dR); }
    }

    // Play with the MC stack if available
    // Check origin of the candidates
    Int_t mcindex = -1;
    if(IsDataMC())
    {
      Int_t tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader(), 0);
            
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = kmcEta;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPhoton;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
                GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcConversion;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))   mcindex = kmcElectron;
      else                                                                                mcindex = kmcHadron;

      fhNLocMax[mcindex][matched]->Fill(en,nMax);
            
      if     (nMax == 1 ) { fhM02NLocMax1[mcindex][matched]->Fill(en,l0) ; fhNCellNLocMax1[mcindex][matched]->Fill(en,nc) ; }
      else if(nMax == 2 ) { fhM02NLocMax2[mcindex][matched]->Fill(en,l0) ; fhNCellNLocMax2[mcindex][matched]->Fill(en,nc) ; }
      else if(nMax >= 3 ) { fhM02NLocMaxN[mcindex][matched]->Fill(en,l0) ; fhNCellNLocMaxN[mcindex][matched]->Fill(en,nc) ; }
      
      if(TMath::Abs(dR) < 999)
      {
        if     ( nMax == 1  ) { fhTrackMatchedDEtaLocMax1[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMax1[mcindex]->Fill(en,dR); }
        else if( nMax == 2  ) { fhTrackMatchedDEtaLocMax2[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMax2[mcindex]->Fill(en,dR); }
        else if( nMax >= 3  ) { fhTrackMatchedDEtaLocMaxN[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiLocMaxN[mcindex]->Fill(en,dR); }
      }
      
    }  

    if     (nMax==1) 
    { 
      if( en > 7 ) 
      {      
        fhMassM02NLocMax1[0]  [matched]->Fill(l0,  mass ); 
        if(IsDataMC()) 
          fhMassM02NLocMax1[mcindex][matched]->Fill(l0,  mass ); 
      }
      
      if(!matched)
      {
        if(en > 6  && en <= 10) fhMassM02NLocMax1Ebin[0]->Fill(l0,  mass ); 
        if(en > 10 && en <= 15) fhMassM02NLocMax1Ebin[1]->Fill(l0,  mass ); 
        if(en > 15 && en <= 20) fhMassM02NLocMax1Ebin[2]->Fill(l0,  mass ); 
        if(en > 20)             fhMassM02NLocMax1Ebin[3]->Fill(l0,  mass ); 
      }
    }  
    else if(nMax==2) 
    {
      if( en > 7 )
      {
        fhMassM02NLocMax2[0][matched]  ->Fill(l0,  mass ); 
        if(IsDataMC()) 
          fhMassM02NLocMax2[mcindex][matched]->Fill(l0,  mass ); 
      }
      
      if(!matched)
      {
        if(en > 6  && en <= 10) fhMassM02NLocMax2Ebin[0]->Fill(l0,  mass ); 
        if(en > 10 && en <= 15) fhMassM02NLocMax2Ebin[1]->Fill(l0,  mass ); 
        if(en > 15 && en <= 20) fhMassM02NLocMax2Ebin[2]->Fill(l0,  mass ); 
        if(en > 20)             fhMassM02NLocMax2Ebin[3]->Fill(l0,  mass );       
      }
    }
    else if(nMax >2) 
    {
      if( en > 7 ) 
      {      
        fhMassM02NLocMaxN[0]  [matched]->Fill(l0  ,mass );
        if(IsDataMC()) 
          fhMassM02NLocMaxN[mcindex][matched]->Fill(l0,  mass ); 
      }
      
      if(!matched)
      {
        if(en > 6  && en <= 10) fhMassM02NLocMaxNEbin[0]->Fill(l0,  mass ); 
        if(en > 10 && en <= 15) fhMassM02NLocMaxNEbin[1]->Fill(l0,  mass ); 
        if(en > 15 && en <= 20) fhMassM02NLocMaxNEbin[2]->Fill(l0,  mass ); 
        if(en > 20)             fhMassM02NLocMaxNEbin[3]->Fill(l0,  mass );
      }
    }
    
    //---------------------------------------------------------------------
    // From here only if M02 is large but not too large, fill histograms 
    //---------------------------------------------------------------------
    
    if( l0 < fM02MinCut || l0 > fM02MaxCut ) continue ;
    
    fhNLocMaxM02Cut[0][matched]->Fill(en,nMax);
    if(IsDataMC()) fhNLocMaxM02Cut[mcindex][matched]->Fill(en,nMax);
    
    if     (nMax==1) 
    { 
      fhAnglePairLocMax1[matched]->Fill(en,angle);
      fhMassNLocMax1[0][matched] ->Fill(en,mass ); 
      
      if( en > 7 ) 
      {      
        fhAnglePairMassLocMax1[matched]->Fill(mass,angle);
      }
      
      if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMax1[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kPi0   ) fhM02Pi0LocMax1[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kEta)    fhM02EtaLocMax1[0][matched]->Fill(en,l0);
    }  
    else if(nMax==2) 
    {
      fhAnglePairLocMax2[matched]->Fill(en,angle);
      fhMassNLocMax2[0] [matched]->Fill(en,mass );
      
      if( en > 7 )
      {
        fhAnglePairMassLocMax2[matched]->Fill(mass,angle);        
      }
      
      if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMax2[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kPi0   ) fhM02Pi0LocMax2[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kEta)    fhM02EtaLocMax2[0][matched]->Fill(en,l0);        
    }
    else if(nMax >2) 
    {
      fhAnglePairLocMaxN[matched]->Fill(en,angle);
      fhMassNLocMaxN[0] [matched]->Fill(en,mass );
      
      if( en > 7 ) 
      {      
        fhAnglePairMassLocMaxN[matched]->Fill(mass,angle);
      }
      
      if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMaxN[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kPi0   ) fhM02Pi0LocMaxN[0][matched]->Fill(en,l0);
      else if(pidTag==AliCaloPID::kEta   ) fhM02EtaLocMaxN[0][matched]->Fill(en,l0);
    }
    
    
    if(IsDataMC())
    {
      if     (nMax==1) 
      { 
        fhMassNLocMax1[mcindex][matched]->Fill(en,mass); 
        if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMax1[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kPi0) fhM02Pi0LocMax1[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kEta) fhM02EtaLocMax1[mcindex][matched]->Fill(en,l0);
      }  
      else if(nMax==2) 
      {
        fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
        if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMax2[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kPi0   ) fhM02Pi0LocMax2[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kEta   ) fhM02EtaLocMax2[mcindex][matched]->Fill(en,l0);        
      }
      else if(nMax >2) 
      {
        fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
        if     (pidTag==AliCaloPID::kPhoton) fhM02ConLocMaxN[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kPi0   ) fhM02Pi0LocMaxN[mcindex][matched]->Fill(en,l0);
        else if(pidTag==AliCaloPID::kEta) fhM02EtaLocMaxN[mcindex][matched]->Fill(en,l0);
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



