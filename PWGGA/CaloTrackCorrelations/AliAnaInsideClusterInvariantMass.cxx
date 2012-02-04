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
#include <TH3F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>

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
  fM02Cut(0),       fMinNCells(0),  
  fMassEtaMin(0),   fMassEtaMax(0),
  fMassPi0Min(0),   fMassPi0Max(0),
  fMassConMin(0),   fMassConMax(0)
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
      fhNLocMaxEMax[i][j]   = 0;
      fhNLocMaxEFrac[i][j]  = 0;
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
  snprintf(onePar,buffersize,"fM02Cut =%2.2f \n",        fM02Cut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinNCells =%d \n",        fMinNCells) ;
  parList+=onePar ;  
  snprintf(onePar,buffersize,"pi0 : %2.1f < m <%2.1f\n", fMassPi0Min,fMassPi0Max);
  parList+=onePar ;
  snprintf(onePar,buffersize,"eta : %2.1f < m <%2.1f\n", fMassEtaMin,fMassEtaMax);
  parList+=onePar ;
  snprintf(onePar,buffersize,"conv: %2.1f < m <%2.1f\n", fMassConMin,fMassConMax);
  parList+=onePar ;

  return new TObjString(parList) ;
  
}


//_____________________________________________________________________________________
TLorentzVector AliAnaInsideClusterInvariantMass::GetCellMomentum(const Int_t absId,
                                                                 Float_t en,
                                                                 AliVCaloCells * cells)
{

  // Cell momentum calculation for invariant mass
  
  Double_t cellpos[] = {0, 0, 0};
  GetEMCALGeometry()->GetGlobal(absId, cellpos);
  
  if(GetVertex(0)){//calculate direction from vertex
    cellpos[0]-=GetVertex(0)[0];
    cellpos[1]-=GetVertex(0)[1];
    cellpos[2]-=GetVertex(0)[2];  
  }
  
  Double_t r = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]+cellpos[2]*cellpos[2] ) ;
  
  //If not calculated before, get the energy
  if(en <=0 )
  {
    en = cells->GetCellAmplitude(absId);
    GetCaloUtils()->RecalibrateCellAmplitude(en,fCalorimeter,absId);  
  }
  
  TLorentzVector cellMom ;   
  cellMom.SetPxPyPzE( en*cellpos[0]/r,  en*cellpos[1]/r, en*cellpos[2]/r,  en) ;   

  return cellMom;
  
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
                                       Form("Invariant mass of 2 highest energy cells %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax1[i][j]) ;   
      
      fhMassNLocMax2[i][j]  = new TH2F(Form("hMassNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of 2 local maxima cells %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMax2[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMax2[i][j]) ;   
      
      fhMassNLocMaxN[i][j]  = new TH2F(Form("hMassNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Invariant mass of N>2 local maxima cells %s %s",ptype[i].Data(),sMatched[j].Data()),
                                       nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassNLocMaxN[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassNLocMaxN[i][j]) ;   
      
      
      fhNLocMax[i][j]     = new TH2F(Form("hNLocMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                     Form("Number of local maxima in cluster %s %s",ptype[i].Data(),sMatched[j].Data()),
                                     nptbins,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMax[i][j]   ->SetYTitle("N maxima");
      fhNLocMax[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhNLocMax[i][j]) ; 
      
      fhNLocMaxEMax[i][j]     = new TH2F(Form("hNLocMaxEMax%s%s",pname[i].Data(),sMatched[j].Data()),
                                         Form("Number of local maxima in cluster vs energy of maxima %s %s",ptype[i].Data(),sMatched[j].Data()),
                                         nptbins*10,ptmin,ptmax,nMaxBins,0,nMaxBins); 
      fhNLocMaxEMax[i][j]   ->SetYTitle("N maxima");
      fhNLocMaxEMax[i][j]   ->SetXTitle("E of maxima (GeV)");
      outputContainer->Add(fhNLocMaxEMax[i][j]) ; 
      
      fhNLocMaxEFrac[i][j]     = new TH2F(Form("hNLocMaxEFrac%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Number of local maxima in cluster vs fraction of cluster energy of maxima %s %s",ptype[i].Data(),sMatched[j].Data()),
                                          100,0,1,nMaxBins,0,nMaxBins); 
      fhNLocMaxEFrac[i][j]   ->SetYTitle("N maxima");
      fhNLocMaxEFrac[i][j]   ->SetXTitle("E maxima / E cluster");
      outputContainer->Add(fhNLocMaxEFrac[i][j]) ; 
      
      fhNLocMaxM02Cut[i][j] = new TH2F(Form("hNLocMaxM02Cut%s%s",pname[i].Data(),sMatched[j].Data()),
                                       Form("Number of local maxima in cluster %s for M02 > %2.2f",ptype[i].Data(),fM02Cut),
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
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 1",fMassPi0Min,fMassPi0Max,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMax1[i][j]) ; 
      
      fhM02EtaLocMax1[i][j]     = new TH2F(Form("hM02EtaLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",fMassEtaMin,fMassEtaMax,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMax1[i][j]) ; 
      
      fhM02ConLocMax1[i][j]    = new TH2F(Form("hM02ConLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 1",fMassConMin,fMassConMax,ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMax1[i][j]) ; 
      
      fhM02Pi0LocMax2[i][j]     = new TH2F(Form("hM02Pi0LocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max = 2",fMassPi0Min,fMassPi0Max,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMax2[i][j]) ; 
      
      fhM02EtaLocMax2[i][j]     = new TH2F(Form("hM02EtaLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",fMassEtaMin,fMassEtaMax,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMax2[i][j]) ; 
      
      fhM02ConLocMax2[i][j]    = new TH2F(Form("hM02ConLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max = 2",fMassConMin,fMassConMax,ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMax2[i][j]) ; 
      
      fhM02Pi0LocMaxN[i][j]     = new TH2F(Form("hM02Pi0LocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2} %s, for N Local max > 2",fMassPi0Min,fMassPi0Max,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0LocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0LocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0LocMaxN[i][j]) ; 
      
      fhM02EtaLocMaxN[i][j]     = new TH2F(Form("hM02EtaLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] MeV/c^{2}, %s, for N Local max > 2", fMassEtaMin,fMassEtaMax,ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaLocMaxN[i][j]) ; 
      
      fhM02ConLocMaxN[i][j]    = new TH2F(Form("hM02ConLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",fMassConMin,fMassConMax,ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConLocMaxN[i][j]) ; 
      
    } // matched, not matched
    
  } // MC particle list
 
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
                                          Form("Opening angle of 2 highest energy cells vs Mass for E > 5 GeV, %s",sMatched[j].Data()),
                                          mbins,mmin,mmax,200,0,0.2); 
    fhAnglePairMassLocMax1[j]->SetXTitle("M (GeV/c^{2})");
    fhAnglePairMassLocMax1[j]->SetYTitle("#alpha (rad)");
    outputContainer->Add(fhAnglePairMassLocMax1[j]) ;   
    
    fhAnglePairMassLocMax2[j]  = new TH2F(Form("hAnglePairMassLocMax2%s",sMatched[j].Data()),
                                          Form("Opening angle of 2 local maxima cells vs Mass for E > 5 GeV, %s",sMatched[j].Data()),
                                          mbins,mmin,mmax,200,0,0.2); 
    fhAnglePairMassLocMax2[j]->SetXTitle("M (GeV/c^{2})");
    fhAnglePairMassLocMax2[j]->SetYTitle("#alpha (rad)");
    outputContainer->Add(fhAnglePairMassLocMax2[j]) ;   
    
    fhAnglePairMassLocMaxN[j]  = new TH2F(Form("hAnglePairMassLocMaxN%s",sMatched[j].Data()),
                                          Form("Opening angle of N>2 local maxima cells vs Mass for E > 5 GeV, %s",sMatched[j].Data()),
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
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kMC ){
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

  fM02Cut      = 0.26 ;
  fMinNCells   = 4 ;
  
  fMassEtaMin  = 0.4;
  fMassEtaMax  = 0.6;
  
  fMassPi0Min  = 0.08;
  fMassPi0Max  = 0.20;
  
  fMassConMin  = 0.0;
  fMassConMax  = 0.05;
  
}


//__________________________________________________________________
void  AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl       = 0x0; 
  AliVCaloCells* cells = 0x0;

  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS"){
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL"){
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl || !cells) {
    Info("MakeAnalysisFillHistograms","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }  
  
	if(fCalorimeter == "PHOS") return; // Not implemented for PHOS yet

  for(Int_t icluster = 0; icluster < pl->GetEntriesFast(); icluster++){
    AliVCluster * cluster = (AliVCluster*) (pl->At(icluster));	

    // Study clusters with large shape parameter
    Float_t en = cluster->E();
    Float_t l0 = cluster->GetM02();
    Int_t   nc = cluster->GetNCells();

    //If too small or big E or low number of cells, skip it
    if( en < GetMinEnergy() || en > GetMaxEnergy() || nc < fMinNCells) continue ; 
  
    Int_t    absId1    = -1; Int_t absId2 = -1;
    Int_t   *absIdList = new Int_t  [nc]; 
    Float_t *maxEList  = new Float_t[nc]; 
    Int_t    nMax      = GetCaloUtils()->GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList) ;
    Bool_t   matched   = IsTrackMatched(cluster,GetReader()->GetInputEvent());
    
    if (nMax <= 0) 
    {
      printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found!\n");
      
      /*
      for(Int_t iDigit  = 0; iDigit < cluster->GetNCells(); iDigit++ ) {
        Float_t ec = cells->GetCellAmplitude(cluster->GetCellsAbsId()[iDigit]);
        GetCaloUtils()->RecalibrateCellAmplitude(ec,fCalorimeter,cluster->GetCellsAbsId()[iDigit]);
        printf("iDigit %d, absId %d, Ecell %f\n",iDigit,cluster->GetCellsAbsId()[iDigit], ec);
      }
      */
      
      delete [] absIdList ;
      delete [] maxEList  ;
      return;
    }
    
    fhNLocMax[0][matched]->Fill(en,nMax);
    for(Int_t imax = 0; imax < nMax; imax++)
    {
      fhNLocMaxEMax [0][matched]->Fill(maxEList[imax]   ,nMax);
      fhNLocMaxEFrac[0][matched]->Fill(maxEList[imax]/en,nMax);
    }
    
    
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
    if(IsDataMC()){
      
      Int_t tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader(), 0);
            
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = kmcEta;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPhoton;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
                GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcConversion;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))   mcindex = kmcElectron;
      else                                                                                mcindex = kmcHadron;
      
      //GetMCAnalysisUtils()->PrintMCTag(tag);
      //printf("\t MC index Assigned %d \n",mcindex);
      
      fhNLocMax[mcindex][matched]->Fill(en,nMax);
      for(Int_t imax = 0; imax < nMax; imax++)
      {
        fhNLocMaxEMax [mcindex][matched]->Fill(maxEList[imax]   ,nMax);
        fhNLocMaxEFrac[mcindex][matched]->Fill(maxEList[imax]/en,nMax);
      }
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
    
    //---------------------------------------------------------------------
    // From here only if M02 is large, fill histograms or split the cluster
    //---------------------------------------------------------------------

    if( l0 < fM02Cut ) 
    {
      delete [] absIdList ;
      delete [] maxEList  ;
      continue;    
    }
        
    fhNLocMaxM02Cut[0][matched]->Fill(en,nMax);
    if(IsDataMC()) fhNLocMaxM02Cut[mcindex][matched]->Fill(en,nMax);
    
    //---------------------------------------------------------------------
    // Get the 2 max indeces and do inv mass
    //---------------------------------------------------------------------

    if     ( nMax == 2 ) {
      absId1 = absIdList[0];
      absId2 = absIdList[1];
    }
    else if( nMax == 1 )
    {
      
      absId1 = absIdList[0];

      //Find second highest energy cell

      Float_t enmax = 0 ;
      for(Int_t iDigit = 0 ; iDigit < cluster->GetNCells() ; iDigit++){
        Int_t absId = cluster->GetCellsAbsId()[iDigit];
        if( absId == absId1 ) continue ; 
        Float_t endig = cells->GetCellAmplitude(absId);
        GetCaloUtils()->RecalibrateCellAmplitude(endig,fCalorimeter,absId); 
        if(endig > enmax) {
          enmax  = endig ;
          absId2 = absId ;
        }
      }// cell loop
    }// 1 maxima 
    else {  // n max > 2
      // loop on maxima, find 2 highest
      
      // First max
      Float_t enmax = 0 ;
      for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++){
        Float_t endig = maxEList[iDigit];
        if(endig > enmax) {
          enmax  = endig ;
          absId1 = absIdList[iDigit];
        }
      }// first maxima loop

      // Second max 
      Float_t enmax2 = 0;
      for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++){
        if(absIdList[iDigit]==absId1) continue;
        Float_t endig = maxEList[iDigit];
        if(endig > enmax2) {
          enmax2  = endig ;
          absId2 = absIdList[iDigit];
        }
      }// second maxima loop

    } // n local maxima > 2
    
    //---------------------------------------------------------------------
    // Split the cluster energy in 2, around the highest 2 local maxima
    //---------------------------------------------------------------------

    //Float_t en1 = 0, en2 = 0;
    //SplitEnergy(absId1,absId2,cluster, cells, en1, en2, nMax /*absIdList, maxEList,*/);

    AliAODCaloCluster *cluster1 = new AliAODCaloCluster(0, 0,NULL,0.,NULL,NULL,1,0);
    AliAODCaloCluster *cluster2 = new AliAODCaloCluster(1, 0,NULL,0.,NULL,NULL,1,0);
    
    SplitEnergy(absId1,absId2,cluster, cells, cluster1, cluster2, nMax /*absIdList, maxEList,*/);

    //---------------------------------------------------------------------
    // Get mass of pair of clusters
    //---------------------------------------------------------------------

    // First set position of cluster as local maxima position, 
    // assign splitted energy to calculate momentum
    
    //TLorentzVector cellMom1 = GetCellMomentum(absId1, en1, cells);
    //TLorentzVector cellMom2 = GetCellMomentum(absId2, en2, cells);

    TLorentzVector cellMom1; 
    TLorentzVector cellMom2;  
    
    cluster1->GetMomentum(cellMom1,GetVertex(0));
    cluster2->GetMomentum(cellMom2,GetVertex(0));
    
    Float_t  mass  = (cellMom1+cellMom2).M();
    Double_t angle = cellMom2.Angle(cellMom1.Vect());

    if     (nMax==1) 
    { 
      Int_t icol1 = -1, irow1 = -1, iRCU1 = -1;
      Int_t icol2 = -1, irow2 = -1, iRCU2 = -1;
      Int_t sm1 = GetCaloUtils()->GetModuleNumberCellIndexes(absId1, fCalorimeter, icol1, irow1, iRCU1) ;
      Int_t sm2 = GetCaloUtils()->GetModuleNumberCellIndexes(absId2, fCalorimeter, icol2, irow2, iRCU2) ;
      Float_t en1 = cluster1->E();
      Float_t en2 = cluster2->E();

     printf("Angle %f, E %f, E1 %f, E2 %f, ecell1 %f, ecell2 %f, Asym %f, Fcell1 %f, Fcell 2 %f, SM1 %d, SM2 %d, icol1 %d, icol2 %d, irow1 %d, irow2 %d \n  ",
             angle,en,en1,en2,
             cells->GetCellAmplitude(absId1),cells->GetCellAmplitude(absId2), 
             TMath::Abs(en1-en2)/(en1+en2),cells->GetCellAmplitude(absId1)/en, cells->GetCellAmplitude(absId2)/en,
             sm1,sm2,icol1,icol2,irow1,irow2);
      
      fhAnglePairLocMax1[matched]->Fill(en,angle);
      if( en > 5 )       
        fhAnglePairMassLocMax1[matched]->Fill(mass,angle);
      
      fhMassNLocMax1[0][matched]->Fill(en,mass); 
      if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMax1[0][matched]->Fill(en,l0);
      else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMax1[0][matched]->Fill(en,l0);
      else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMax1[0][matched]->Fill(en,l0);
    }  
    else if(nMax==2) 
    {
      fhAnglePairLocMax2[matched]->Fill(en,angle);
      if( en > 5 )       
        fhAnglePairMassLocMax2[matched]->Fill(mass,angle);
      
      fhMassNLocMax2[0][matched]->Fill(en,mass);
      if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMax2[0][matched]->Fill(en,l0);
      else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMax2[0][matched]->Fill(en,l0);
      else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMax2[0][matched]->Fill(en,l0);        
    }
    else if(nMax >2) 
    {
      fhAnglePairLocMaxN[matched]->Fill(en,angle);
      if( en > 5 )       
        fhAnglePairMassLocMaxN[matched]->Fill(mass,angle);
      
      fhMassNLocMaxN[0][matched]->Fill(en,mass);
      if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMaxN[0][matched]->Fill(en,l0);
      else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMaxN[0][matched]->Fill(en,l0);
      else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMaxN[0][matched]->Fill(en,l0);
    }
    
    if(IsDataMC()){
            
      if     (nMax==1) 
      { 
        fhMassNLocMax1[mcindex][matched]->Fill(en,mass); 
        if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMax1[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMax1[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMax1[mcindex][matched]->Fill(en,l0);
      }  
      else if(nMax==2) {
        fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
        if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMax2[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMax2[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMax2[mcindex][matched]->Fill(en,l0);        
      }
      else if(nMax >2) {
        fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
        if     (mass < fMassConMax && mass > fMassConMin) fhM02ConLocMaxN[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassPi0Max && mass > fMassPi0Min) fhM02Pi0LocMaxN[mcindex][matched]->Fill(en,l0);
        else if(mass < fMassEtaMax && mass > fMassEtaMin) fhM02EtaLocMaxN[mcindex][matched]->Fill(en,l0);
      }
      
    }//Work with MC truth first
    
    delete cluster1 ;
    delete cluster2 ;
    delete [] absIdList ;
    delete [] maxEList  ;

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
  printf("lambda_0^2 >  %2.1f \n",      fM02Cut);
  printf("pi0 : %2.2f<m<%2.2f \n",      fMassPi0Min,fMassPi0Max);
  printf("eta : %2.2f<m<%2.2f \n",      fMassEtaMin,fMassEtaMax);
  printf("conv: %2.2f<m<%2.2f \n",      fMassConMin,fMassConMax);

  printf("    \n") ;
  
} 



//________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::SplitEnergy(const Int_t absId1, const Int_t absId2,
                                                   AliVCluster* cluster, 
                                                   AliVCaloCells* cells,
                                                   //Float_t & e1, Float_t & e2,
                                                   AliAODCaloCluster* cluster1,
                                                   AliAODCaloCluster* cluster2,
                                                   const Int_t nMax)
{

  // Split energy of cluster between the 2 local maxima, sum energy on 3x3, and if the 2 
  // maxima are too close and have common cells, split the energy between the 2
  
/*
  TH2F* hClusterMap    = new TH2F("hClusterMap","Cluster Map",48,0,48,24,0,24);
  TH2F* hClusterLocMax = new TH2F("hClusterLocMax","Cluster Local Maxima",48,0,48,24,0,24);
  TH2F* hCluster1      = new TH2F("hCluster1","Cluster 1",48,0,48,24,0,24);
  TH2F* hCluster2      = new TH2F("hCluster2","Cluster 2",48,0,48,24,0,24);
*/
  
  const Int_t ncells  = cluster->GetNCells();  
  Int_t absIdList[ncells]; 
  
  Float_t e1 = 0,  e2   = 0 ;
  Int_t icol = -1, irow = -1, iRCU = -1, sm = -1;
  
  //printf("Split Local Max: 1) %d - 2) %d\n",absId1,absId2);
  Float_t eCluster = 0;
  for(Int_t iDigit  = 0; iDigit < ncells; iDigit++ ) {
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit];
    
    //printf("iDigit %d, absId %d, Ecell %f\n",iDigit,absIdList[iDigit], cells->GetCellAmplitude(absIdList[iDigit]));

    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList[iDigit], fCalorimeter, icol, irow, iRCU) ;
    Float_t ec = cells->GetCellAmplitude(absIdList[iDigit]);
    GetCaloUtils()->RecalibrateCellAmplitude(ec,fCalorimeter, absIdList[iDigit]);
    eCluster+=ec;
    
/*    hClusterMap->Fill(icol,irow,ec); */
     
  }
    
  // Init counters and variables
  Int_t ncells1 = 1 ;
  UShort_t absIdList1[9] ;  
  Double_t fracList1 [9] ;  
  absIdList1[0] = absId1 ;
  fracList1 [0] = 1. ;
  
  Float_t ecell1 = cells->GetCellAmplitude(absId1);
  GetCaloUtils()->RecalibrateCellAmplitude(ecell1, fCalorimeter, absId1);
  e1 =  ecell1;  
  
  Int_t ncells2 = 1 ;
  UShort_t absIdList2[9] ;  
  Double_t fracList2 [9] ; 
  absIdList2[0] = absId2 ;
  fracList2 [0] = 1. ;

  Float_t ecell2 = cells->GetCellAmplitude(absId2);
  GetCaloUtils()->RecalibrateCellAmplitude(ecell2, fCalorimeter, absId2);
  e2 =  ecell2;  
  
  /*
  Int_t icol1 = -1, irow1 = -1, icol2 = -1, irow2 = -1;
  sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId1, fCalorimeter, icol1, irow1, iRCU) ;
  hClusterLocMax->Fill(icol1,irow1,ecell1);
  sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId2, fCalorimeter, icol2, irow2, iRCU) ;
  hClusterLocMax->Fill(icol2,irow2,ecell2);
  */
  
  // Very rough way to share the cluster energy
  Float_t eRemain = (eCluster-ecell1-ecell2)/2;
  Float_t shareFraction1 = ecell1/eCluster+eRemain/eCluster;
  Float_t shareFraction2 = ecell2/eCluster+eRemain/eCluster;
 
  for(Int_t iDigit = 0; iDigit < ncells; iDigit++){
    Int_t absId = absIdList[iDigit];
    
    if(absId==absId1 || absId==absId2 || absId < 0) continue;
    
    Float_t ecell = cells->GetCellAmplitude(absId);
    GetCaloUtils()->RecalibrateCellAmplitude(ecell, fCalorimeter, absId);
    
    if(GetCaloUtils()->AreNeighbours(fCalorimeter, absId1,absId ))
    { 
       absIdList1[ncells1]= absId;
    
      if(GetCaloUtils()->AreNeighbours(fCalorimeter, absId2,absId ))
      { 
        fracList1[ncells1] = shareFraction1; 
        e1 += ecell*shareFraction1;
      }
      else 
      {
        fracList1[ncells1] = 1.; 
        e1 += ecell;
      }
      
      ncells1++;
      
     } // neigbour to cell1
    
    if(GetCaloUtils()->AreNeighbours(fCalorimeter, absId2,absId ))
    { 
      absIdList2[ncells2]= absId;
     
      if(GetCaloUtils()->AreNeighbours(fCalorimeter, absId1,absId ))
      { 
        fracList2[ncells2] = shareFraction2; 
        e2 += ecell*shareFraction2;
      }
      else
      { 
        fracList2[ncells2] = 1.; 
        e2 += ecell;
      }

      ncells2++;

    } // neigbour to cell2

  }
  
   if(GetDebug() > 1) printf("AliAnaInsideClusterInvariantMass::SplitEnergy() - n Local Max %d, Cluster energy  = %f, Ecell1 = %f, Ecell2 = %f, Enew1 = %f, Enew2 = %f, Remain %f, \n ncells %d, ncells1 %d, ncells2 %d, f1 %f, f2  %f, sum f12 = %f \n",
         nMax, eCluster,ecell1,ecell2,e1,e2,eCluster-e1-e2,ncells,ncells1,ncells2,shareFraction1,shareFraction2,shareFraction1+shareFraction2);

  cluster1->SetE(e1);
  cluster2->SetE(e2);  
  
  cluster1->SetNCells(ncells1);
  cluster2->SetNCells(ncells2);  
  
  cluster1->SetCellsAbsId(absIdList1);
  cluster2->SetCellsAbsId(absIdList2);
  
  cluster1->SetCellsAmplitudeFraction(fracList1);
  cluster2->SetCellsAmplitudeFraction(fracList2);
  
  GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterPosition(GetEMCALGeometry(), cells, cluster1);
  GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterPosition(GetEMCALGeometry(), cells, cluster2);
  
  
/*  
  printf("Cells of cluster1: ");
  for(Int_t iDigit  = 0; iDigit < ncells1; iDigit++ ) 
  {
    printf(" %d ",absIdList1[iDigit]);
    
    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList1[iDigit], fCalorimeter, icol, irow, iRCU) ;
    
    if( GetCaloUtils()->AreNeighbours(fCalorimeter, absId2,absIdList1[iDigit]) )
      hCluster1->Fill(icol,irow,cells->GetCellAmplitude(absIdList1[iDigit])*shareFraction1);
    else 
      hCluster1->Fill(icol,irow,cells->GetCellAmplitude(absIdList1[iDigit]));
  }
  
  printf(" \n ");
  printf("Cells of cluster2: ");

  for(Int_t iDigit  = 0; iDigit < ncells2; iDigit++ ) 
  {
    printf(" %d ",absIdList2[iDigit]);

    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList2[iDigit], fCalorimeter, icol, irow, iRCU) ;
    if( GetCaloUtils()->AreNeighbours(fCalorimeter, absId1,absIdList2[iDigit]) )
      hCluster2->Fill(icol,irow,cells->GetCellAmplitude(absIdList2[iDigit])*shareFraction2);
    else
      hCluster2->Fill(icol,irow,cells->GetCellAmplitude(absIdList2[iDigit]));

  }
  printf(" \n ");
   
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadLeftMargin(0.1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(000000);
   
   TCanvas  * c= new TCanvas("canvas", "canvas", 8000, 4000) ;
   c->Divide(2,2);  
   c->cd(1);
   gPad->SetGridy();
   gPad->SetGridx();
   hClusterMap->Draw("colz");
   c->cd(2);
   gPad->SetGridy();
   gPad->SetGridx();
   hClusterLocMax ->Draw("colz");
   c->cd(3);
   gPad->SetGridy();
   gPad->SetGridx();
   hCluster1      ->Draw("colz");
   c->cd(4);
   gPad->SetGridy();
   gPad->SetGridx();
   hCluster2      ->Draw("colz");
   
   c->Print(Form("Event%d_nMax%d_NCell1_%d_NCell2_%d.eps",GetEventNumber(),nMax,ncells1,ncells2));
   
   delete c;
   delete hClusterMap;
   delete hClusterLocMax;
   delete hCluster1;
   delete hCluster2;
*/
}


//________________________________________________________________________________________
//void AliAnaInsideClusterInvariantMass::SplitEnergy(Int_t absId1, Int_t absId2,
//                                                   AliVCluster* cluster, 
//                                                   AliVCaloCells* cells,
//                                                   Float_t & e1, Float_t & e2,
//                                                   const Int_t nMax, Int_t *listMax, Float_t *eListMax,)
//{
//  
//  // Split energy of cluster between the 2 local maxima.
//  const Int_t ncells  = cluster->GetNCells();  
//  Int_t     absIdList[ncells]; 
//  //Int_t icol = -1, irow = -1, iRCU = -1, sm = -1;
///*
//  TH2F* hClusterMap    = new TH2F("hClusterMap","Cluster Map",48,0,48,24,0,24);
//  TH2F* hClusterLocMax = new TH2F("hClusterLocMax","Cluster Local Maxima",48,0,48,24,0,24);
//  TH2F* hClusterLocMax2= new TH2F("hClusterLocMax2","Cluster Highest Local Maxima",48,0,48,24,0,24);
//  TH2F* hCluster1      = new TH2F("hCluster1","Cluster 1",48,0,48,24,0,24);
//  TH2F* hCluster2      = new TH2F("hCluster2","Cluster 2",48,0,48,24,0,24);
//  TH2F* hClusterNo     = new TH2F("hClusterNo","Cluster No",48,0,48,24,0,24);
// */
//  Float_t ec = 0;
//  //printf("Split Local Max: 1) %d - 2) %d\n",absId1,absId2);
//  for(Int_t iDigit  = 0; iDigit < ncells; iDigit++ ) {
//    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit];
//    
//    //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList[iDigit], fCalorimeter, icol, irow, iRCU) ;
//    //ec = cells->GetCellAmplitude(absIdList[iDigit]);
//    //GetCaloUtils()->RecalibrateCellAmplitude(ec,fCalorimeter, absIdList[iDigit]);
//    //hClusterMap->Fill(icol,irow,ec);
//    
//    //printf("iDigit %d, absId %d, Ecell %f\n",iDigit,absIdList[iDigit], cells->GetCellAmplitude(absIdList[iDigit]));
//
//  }
//   /* 
//  printf("N Local Maxima %d \n",nMax);
//  for(Int_t imax = 0; imax < nMax; imax++)
//  {
//    sm = GetCaloUtils()->GetModuleNumberCellIndexes(listMax[imax], fCalorimeter, icol, irow, iRCU) ;
//    printf("LocalMaxima absId %d, Ecell %f\n",absIdList[imax], cells->GetCellAmplitude(listMax[imax]));
//    hClusterLocMax->Fill(icol,irow,eListMax[imax]);
//  }
//  */
//  
//  //Int_t icol1 = -1, irow1 = -1, icol2 = -1, irow2 = -1;
//  //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId1, fCalorimeter, icol1, irow1, iRCU) ;
//  Float_t ec1 = cells->GetCellAmplitude(absId1);
//  GetCaloUtils()->RecalibrateCellAmplitude(ec1,fCalorimeter, absId1);
//  //hClusterLocMax2->Fill(icol1,irow1,ec1);
//  
//  //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId2, fCalorimeter, icol2, irow2, iRCU) ;
//  Float_t ec2 = cells->GetCellAmplitude(absId2);
//  GetCaloUtils()->RecalibrateCellAmplitude(ec2,fCalorimeter, absId2);
//  //hClusterLocMax2->Fill(icol2,irow2,ec2);
//
//  Int_t absIdtmp = 0;
//  if(ec2>ec1){
//    absIdtmp = absId2;
//    absId2   = absId1;
//    absId1   = absIdtmp;
//  }
//  
//  // SubCluster 1
//
//  // Init counters and variables
//  Int_t ncells1 = 1;
//  Int_t absIdList1[ncells];  
//  absIdList1[0] = absId1;
//  //printf("First local max : absId1 %d %d \n",absId1, absIdList1[0]);  
//  for(Int_t iDigit1 = 1; iDigit1 < ncells; iDigit1++) absIdList1[iDigit1] = -1;
//  
//  Float_t ecell1 = cells->GetCellAmplitude(absId1);
//  GetCaloUtils()->RecalibrateCellAmplitude(ecell1,fCalorimeter, absId1);
//  e1 =  ecell1;  
//  
//  //Int_t icolNew = -1, irowNew = -1, iRCUNew = -1;
//  //Int_t jcol = -1, jrow = -1, jRCU = -1;
//
//  Bool_t added = kTRUE;
//  Int_t cellj = 0;
//  while (cellj < ncells1) 
//  {
//    added = kFALSE;
//    Int_t absId1New = absIdList1[cellj];
//    //printf("\t absId %d added \n",absId1New);
//    
//    Float_t e1New = cells->GetCellAmplitude(absId1New);
//    GetCaloUtils()->RecalibrateCellAmplitude(e1New,fCalorimeter, absId1New);
//
//    //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId1New, fCalorimeter, icolNew, irowNew, iRCUNew) ;
//    
//    for(Int_t iDigit = 0; iDigit < ncells ; iDigit++)
//    {
//      Int_t absId = absIdList[iDigit] ;
//      if(absId!=absId1New && absId!=absId2 && absId>=0)
//      {
//        Float_t en = cells->GetCellAmplitude(absId);
//        GetCaloUtils()->RecalibrateCellAmplitude(en,fCalorimeter, absId);
//        //printf("\t \t iDig %d, absId %d, absIdNew %d, en %f, enNew %f\n",iDigit,absId, absId1New,en, e1New);
//        //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId, fCalorimeter, jcol, jrow, jRCU) ;
//        //printf("\t \t \t (col,row) New  (%d,%d), check (%d,%d) \n",icolNew, irowNew,jcol,jrow);
//        if(GetCaloUtils()->AreNeighbours(fCalorimeter, absId1New,absId )){ 
//          //printf("\t \t \t neighbours\n");
//          if(e1New > en-GetCaloUtils()->GetLocalMaximaCutEDiff()){
//            absIdList1[ncells1++] = absId; 
//            
//            if((absId1New==absId1 && GetCaloUtils()->AreNeighbours(fCalorimeter, absId1,absId ) && GetCaloUtils()->AreNeighbours( absId2,absId ))) {
//              e1+=en/2;
//            }
//            else {
//              absIdList [iDigit]    = -1; 
//              e1+=en;
//            }
//          } // Decreasing energy with respect reference
//        } // Neighbours
//      } //Not local maxima or already removed
//    } // cell loop
//    cellj++;
//  }// while cells added to list of cells for cl1
//  
//  // SubCluster 2
//  
//  // Init counters and variables
//  Int_t ncells2 = 1;
//  Int_t absIdList2[ncells];  
//  absIdList2[0] = absId2;
//  //printf("Second local max : absId2 %d %d \n",absId2, absIdList2[0]);  
//  for(Int_t iDigit2 = 1; iDigit2 < ncells; iDigit2++) absIdList2[iDigit2] = -1;
//  
//  Float_t ecell2 = cells->GetCellAmplitude(absId2);
//  GetCaloUtils()->RecalibrateCellAmplitude(ecell2,fCalorimeter, absId2);
//  e2 =  ecell2;  
//    
//  added = kTRUE;
//  cellj = 0;
//  while (cellj < ncells2) 
//  {
//    added = kFALSE;
//    Int_t absId2New = absIdList2[cellj];
//    //printf("\t absId %d added \n",absId2New);
//    
//    Float_t e2New = cells->GetCellAmplitude(absId2New);
//    GetCaloUtils()->RecalibrateCellAmplitude(e2New,fCalorimeter,absId2New);
//    //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId2New, fCalorimeter, icolNew, irowNew, iRCU) ;
//
//    for(Int_t iDigit = 0; iDigit < ncells ; iDigit++)
//    {
//      Int_t absId = absIdList[iDigit] ;
//      if(absId!=absId2New && absId>=0)
//      {
//        Float_t en = cells->GetCellAmplitude(absId);
//        GetCaloUtils()->RecalibrateCellAmplitude(en,fCalorimeter, absId);
//        //printf("\t \t iDig %d, absId %d, absIdNew %d, en %f, enNew %f\n",iDigit,absId, absId2New,en, e2New);
//        //sm = GetCaloUtils()->GetModuleNumberCellIndexes(absId, fCalorimeter, jcol, jrow, jRCU) ;
//        //printf("\t \t \t (col,row) New  (%d,%d), check (%d,%d) \n",icolNew, irowNew,jcol,jrow);        
//        if(GetCaloUtils()->AreNeighbours( fCalorimeter, absId2New,absId )){ 
//          //printf("\t \t \t neighbours\n");
//          if(e2New > en-GetCaloUtils()->GetLocalMaximaCutEDiff()){
//            absIdList2[ncells2++] = absId; 
//            absIdList [iDigit]    = -1; 
//            if(absId2New==absId2 && GetCaloUtils()->AreNeighbours(fCalorimeter, absId1,absId ) && GetCaloUtils()->AreNeighbours( fCalorimeter,absId2,absId )){
//              e2+=en/2;
//            }
//            else {
//              e2+=en;
//            }
//          } // Decreasing energy with respect reference
//        } // Neighbours
//      } //Not local maxima or already removed
//    } // cell loop
//    cellj++;
//  }// while cells added to list of cells for cl2  
// /* 
//  for(Int_t iDigit  = 0; iDigit < ncells1; iDigit++ ) {
//    
//    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList1[iDigit], fCalorimeter, icol, irow, iRCU) ;
//    
//    hCluster1->Fill(icol,irow,cells->GetCellAmplitude(absIdList1[iDigit]));
//    
//    
//  }
//  
//  for(Int_t iDigit  = 0; iDigit < ncells2; iDigit++ ) {
//    
//    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList2[iDigit], fCalorimeter, icol, irow, iRCU) ;
//    
//    hCluster2->Fill(icol,irow,cells->GetCellAmplitude(absIdList2[iDigit]));
//    
//    
//  }
//  
//  
//  for(Int_t iDigit  = 0; iDigit < ncells; iDigit++ ) {
//    if(absIdList[iDigit] < 0 || absIdList[iDigit]==absId1 || absIdList[iDigit]==absId2) continue;
//    sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList[iDigit], fCalorimeter, icol, irow, iRCU) ;
//    hClusterNo->Fill(icol,irow,cells->GetCellAmplitude(absIdList[iDigit]));
//  }
//  
//  
//  printf("Cluster energy  = %f, Ecell1 = %f, Ecell2 = %f, Enew1 = %f, Enew2 = %f, ncells %d, ncells1 %d, ncells2 %d \n",
//         cluster->E(),ecell1,ecell2,e1,e2,ncells,ncells1,ncells2);
//  if(ncells!=(ncells1+ncells2)) printf("\t Not all cells!\n");
//  
//  gStyle->SetPadRightMargin(0.15);
//  gStyle->SetPadLeftMargin(0.1);
//  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(000000);
//
//  TCanvas  * c= new TCanvas("canvas", "canvas", 8000, 4000) ;
//  c->Divide(3,2);  
//  c->cd(1);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hClusterMap->Draw("colz");
//  c->cd(2);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hClusterLocMax ->Draw("colz");
//  c->cd(3);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hClusterLocMax2->Draw("colz");
//  c->cd(4);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hCluster1      ->Draw("colz");
//  c->cd(5);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hCluster2      ->Draw("colz");
//  c->cd(6);
//  gPad->SetGridy();
//  gPad->SetGridx();
//  hClusterNo     ->Draw("colz");
//
//  c->Print(Form("Event%d_nMax%d_NCell1_%d_NCell2_%d_Left%d.eps",GetEventNumber(),nMax,ncells1,ncells2,ncells-ncells1-ncells2));
//  
//  delete c;
//  delete hClusterMap;
//  delete hClusterLocMax;
//  delete hClusterLocMax2;
//  delete hCluster1;
//  delete hCluster2;
//  delete hClusterNo;
//*/
//}
//

