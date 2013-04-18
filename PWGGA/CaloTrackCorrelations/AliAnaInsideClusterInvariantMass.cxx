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
  fFillTMHisto(kFALSE),
  fFillTMResidualHisto(kFALSE),
  fFillSSExtraHisto(kFALSE),
  fFillMCFractionHisto(kFALSE),
  fFillSSWeightHisto(kFALSE),
  fFillEbinHisto(0),
  fSSWeightN(0),              fSSECellCutN(0),
  fWSimu(0),
  fhMassM02CutNLocMax1(0),    fhMassM02CutNLocMax2(0),    fhMassM02CutNLocMaxN(0),
  fhAsymM02CutNLocMax1(0),    fhAsymM02CutNLocMax2(0),    fhAsymM02CutNLocMaxN(0),
  fhMassSplitECutNLocMax1(0), fhMassSplitECutNLocMax2(0), fhMassSplitECutNLocMaxN(0),
  fhMCGenFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenFracAfterCutsNLocMaxNMCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax1MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMax2MCPi0(0),
  fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0(0),
  fhEventPlanePi0NLocMax1(0), fhEventPlaneEtaNLocMax1(0),
  fhEventPlanePi0NLocMax2(0), fhEventPlaneEtaNLocMax2(0),
  fhEventPlanePi0NLocMaxN(0), fhEventPlaneEtaNLocMaxN(0),
  fhClusterEtaPhiNLocMax1(0), fhClusterEtaPhiNLocMax2(0),  fhClusterEtaPhiNLocMaxN(0),
  fhPi0EtaPhiNLocMax1(0),     fhPi0EtaPhiNLocMax2(0),      fhPi0EtaPhiNLocMaxN(0),
  fhEtaEtaPhiNLocMax1(0),     fhEtaEtaPhiNLocMax2(0),      fhEtaEtaPhiNLocMaxN(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 8; i++)
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
      fhM02Pi0NLocMax1[i][j] = 0;
      fhM02EtaNLocMax1[i][j] = 0;
      fhM02ConNLocMax1[i][j] = 0;
      fhM02Pi0NLocMax2[i][j] = 0;
      fhM02EtaNLocMax2[i][j] = 0;
      fhM02ConNLocMax2[i][j] = 0;
      fhM02Pi0NLocMaxN[i][j] = 0;
      fhM02EtaNLocMaxN[i][j] = 0;
      fhM02ConNLocMaxN[i][j] = 0;
      
      fhMassPi0NLocMax1[i][j] = 0;
      fhMassEtaNLocMax1[i][j] = 0;
      fhMassConNLocMax1[i][j] = 0;
      fhMassPi0NLocMax2[i][j] = 0;
      fhMassEtaNLocMax2[i][j] = 0;
      fhMassConNLocMax2[i][j] = 0;
      fhMassPi0NLocMaxN[i][j] = 0;
      fhMassEtaNLocMaxN[i][j] = 0;
      fhMassConNLocMaxN[i][j] = 0;
      
      
      fhAsyPi0NLocMax1[i][j] = 0;
      fhAsyEtaNLocMax1[i][j] = 0;
      fhAsyConNLocMax1[i][j] = 0;
      fhAsyPi0NLocMax2[i][j] = 0;
      fhAsyEtaNLocMax2[i][j] = 0;
      fhAsyConNLocMax2[i][j] = 0;
      fhAsyPi0NLocMaxN[i][j] = 0;
      fhAsyEtaNLocMaxN[i][j] = 0;
      fhAsyConNLocMaxN[i][j] = 0;      
      
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
      
      fhMassAfterCutsNLocMax1[i][j] = 0;
      fhMassAfterCutsNLocMax2[i][j] = 0;
      fhMassAfterCutsNLocMaxN[i][j] = 0;
      
      fhSplitEFractionAfterCutsNLocMax1[i][j] = 0 ;
      fhSplitEFractionAfterCutsNLocMax2[i][j] = 0 ;
      fhSplitEFractionAfterCutsNLocMaxN[i][j] = 0 ;
      
      fhCentralityPi0NLocMax1[i][j] = 0 ;
      fhCentralityEtaNLocMax1[i][j] = 0 ;

      fhCentralityPi0NLocMax2[i][j] = 0 ;
      fhCentralityEtaNLocMax2[i][j] = 0 ;

      fhCentralityPi0NLocMaxN[i][j] = 0 ;
      fhCentralityEtaNLocMaxN[i][j] = 0 ;      
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
    
    fhTrackMatchedDEtaNLocMax1[i] = 0; 
    fhTrackMatchedDPhiNLocMax1[i] = 0;
    fhTrackMatchedDEtaNLocMax2[i] = 0;
    fhTrackMatchedDPhiNLocMax2[i] = 0; 
    fhTrackMatchedDEtaNLocMaxN[i] = 0; 
    fhTrackMatchedDPhiNLocMaxN[i] = 0; 
    
  }
   
  for(Int_t i = 0; i < 2; i++)
  {
    fhAnglePairNLocMax1    [i] = 0;
    fhAnglePairNLocMax2    [i] = 0;
    fhAnglePairNLocMaxN    [i] = 0;
    fhAnglePairMassNLocMax1[i] = 0;
    fhAnglePairMassNLocMax2[i] = 0;
    fhAnglePairMassNLocMaxN[i] = 0;
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
  
  
  for(Int_t nlm = 0; nlm < 3; nlm++)
  {
    fhPi0CellE       [nlm] = 0 ;
    fhPi0CellEFrac   [nlm] = 0 ;
    fhPi0CellLogEFrac[nlm] = 0 ;
    
    fhPi0CellEMaxEMax2Frac   [nlm] = 0 ;
    fhPi0CellEMaxClusterFrac [nlm] = 0 ;
    fhPi0CellEMax2ClusterFrac[nlm] = 0 ;

    fhPi0CellEMaxFrac [nlm] = 0 ;
    fhPi0CellEMax2Frac[nlm] = 0 ;
    
    for(Int_t i = 0; i < 10; i++)
    {
      fhM02WeightPi0  [nlm][i] = 0;
      fhM02ECellCutPi0[nlm][i] = 0;
    }
  }
  
  InitParameters();

}

//__________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::FillSSWeightHistograms(AliVCluster *clus,  const Int_t nlm,
                                                              const Int_t absId1, const Int_t absId2)
{
  // Calculate weights and fill histograms
  
  if(!fFillSSWeightHisto) return;
  
  AliVCaloCells* cells = 0;
  if(fCalorimeter == "EMCAL") cells = GetEMCALCells();
  else                        cells = GetPHOSCells();
  
  // First recalculate energy in case non linearity was applied
  Float_t  energy = 0;
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    energy    += amp;
      
  } // energy loop
  
  if(energy <=0 )
  {
    printf("AliAnaInsideClusterInvatiantMass::WeightHistograms()- Wrong calculated energy %f\n",energy);
    return;
  }
  
  //Get amplitude of  main local maxima, recalibrate if needed
  Float_t amp1 = cells->GetCellAmplitude(absId1);
  GetCaloUtils()->RecalibrateCellAmplitude(amp1,fCalorimeter, absId1);
  Float_t amp2 = cells->GetCellAmplitude(absId2);
  GetCaloUtils()->RecalibrateCellAmplitude(amp2,fCalorimeter, absId2);

  if(amp1 < amp2)        printf("Bad local maxima E ordering : id1 E %f, id2 E %f\n ",amp1,amp2);
  if(amp1==0 || amp2==0) printf("Null E local maxima : id1 E %f, id2 E %f\n "        ,amp1,amp2);
  
  if(amp1>0)fhPi0CellEMaxEMax2Frac   [nlm]->Fill(energy,amp2/amp1);
  fhPi0CellEMaxClusterFrac [nlm]->Fill(energy,amp1/energy);
  fhPi0CellEMax2ClusterFrac[nlm]->Fill(energy,amp2/energy);
  
  //Get the ratio and log ratio to all cells in cluster
  for (Int_t ipos = 0; ipos < clus->GetNCells(); ipos++)
  {
    Int_t id       = clus->GetCellsAbsId()[ipos];
    
    //Recalibrate cell energy if needed
    Float_t amp = cells->GetCellAmplitude(id);
    GetCaloUtils()->RecalibrateCellAmplitude(amp,fCalorimeter, id);
    
    if(amp > 0)fhPi0CellE       [nlm]->Fill(energy,amp);
    fhPi0CellEFrac   [nlm]->Fill(energy,amp/energy);
    fhPi0CellLogEFrac[nlm]->Fill(energy,TMath::Log(amp/energy));
    
    if     (id!=absId1 && id!=absId2)
    {
      if(amp1>0)fhPi0CellEMaxFrac [nlm]->Fill(energy,amp/amp1);
      if(amp2>0)fhPi0CellEMax2Frac[nlm]->Fill(energy,amp/amp2);
    }

  }
  
  //Recalculate shower shape for different W0
  if(fCalorimeter=="EMCAL")
  {
    Float_t l0org = clus->GetM02();
    Float_t l1org = clus->GetM20();
    Float_t dorg  = clus->GetDispersion();
    Float_t w0org =  GetCaloUtils()->GetEMCALRecoUtils()->GetW0();
    
    for(Int_t iw = 0; iw < fSSWeightN; iw++)
    {
      GetCaloUtils()->GetEMCALRecoUtils()->SetW0(fSSWeight[iw]);
      //GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParameters(GetEMCALGeometry(), cells, clus);
      
      Float_t l0   = 0., l1   = 0.;
      Float_t disp = 0., dEta = 0., dPhi    = 0.;
      Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
      
      RecalculateClusterShowerShapeParametersWithCellCut(GetEMCALGeometry(), cells, clus,l0,l1,disp,
                                                         dEta, dPhi, sEta, sPhi, sEtaPhi,0);

      
      fhM02WeightPi0[nlm][iw]->Fill(energy,clus->GetM02());
      
    } // w0 loop
    
    // Set the original values back
    clus->SetM02(l0org);
    clus->SetM20(l1org);
    clus->SetDispersion(dorg);
    GetCaloUtils()->GetEMCALRecoUtils()->SetW0(w0org);

    for(Int_t iec = 0; iec < fSSECellCutN; iec++)
    {
      Float_t l0   = 0., l1   = 0.;
      Float_t disp = 0., dEta = 0., dPhi    = 0.;
      Float_t sEta = 0., sPhi = 0., sEtaPhi = 0.;
      
      RecalculateClusterShowerShapeParametersWithCellCut(GetEMCALGeometry(), cells, clus,l0,l1,disp,
                                                         dEta, dPhi, sEta, sPhi, sEtaPhi,fSSECellCut[iec]);
      
      //printf("E %f, l0 org %f, l0 new %f, slope %f\n",clus->E(),l0org,l0,fSSECellCut[iec]);
      fhM02ECellCutPi0[nlm][iec]->Fill(energy,l0);
      
    } // w0 loop
  
  }// EMCAL
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
  snprintf(onePar,buffersize,"fNLocMaxCutE =%2.2f \n",    GetCaloUtils()->GetLocalMaximaCutE()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNLocMaxCutEDiff =%2.2f \n",GetCaloUtils()->GetLocalMaximaCutEDiff()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"%2.2f< M02 < %2.2f \n",    fM02MinCut, fM02MaxCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinNCells =%d \n",        fMinNCells) ;
  parList+=onePar ;    
  snprintf(onePar,buffersize,"fMinBadDist =%1.1f \n",    fMinBadDist) ;
  parList+=onePar ;  
  if(fFillSSWeightHisto)
  {
    snprintf(onePar,buffersize," N w %d - N e cut %d \n",fSSWeightN,fSSECellCutN);
    parList+=onePar ;
  }
  
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
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();          Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();          Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();          Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();          Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();

  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();          
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();          
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();          
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();          
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();  
  
  TString ptype[] ={"","#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron","#pi^{0} (#gamma->e^{#pm})"}; 
  TString pname[] ={"","Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron","Pi0Conv"};
  
  Int_t n = 1;
  
  if(IsDataMC()) n = 8;
  
  Int_t nMaxBins = 10;
  
  TString sMatched[] = {"","Matched"};
  
  Int_t nMatched = 2;
  if(!fFillTMHisto) nMatched = 1;
  
  for(Int_t i = 0; i < n; i++)
  {
    for(Int_t j = 0; j < nMatched; j++)
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
     
      if(i==0 && j==0)
      {
        fhMassSplitECutNLocMax1  = new TH2F("hMassSplitECutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMax1->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMax1->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMax1) ;
        
        fhMassSplitECutNLocMax2  = new TH2F("hMassSplitECutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMax2->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMax2->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMax2) ;
        
        fhMassSplitECutNLocMaxN  = new TH2F("hMassSplitECutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, (E1+E2)/E cut, no TM",
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassSplitECutNLocMaxN->SetYTitle("M (GeV/c^{2})");
        fhMassSplitECutNLocMaxN->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassSplitECutNLocMaxN) ;
        
        fhMassM02CutNLocMax1  = new TH2F("hMassM02CutNLocMax1","Invariant mass of splitted cluster with NLM=1 vs E, M02 cut, no TM",
                                         nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassM02CutNLocMax1->SetYTitle("M (GeV/c^{2})");
        fhMassM02CutNLocMax1->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassM02CutNLocMax1) ;
        
        fhMassM02CutNLocMax2  = new TH2F("hMassM02CutNLocMax2","Invariant mass of splitted cluster with NLM=2 vs E, M02 cut, no TM",
                                         nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassM02CutNLocMax2->SetYTitle("M (GeV/c^{2})");
        fhMassM02CutNLocMax2->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassM02CutNLocMax2) ;
        
        fhMassM02CutNLocMaxN  = new TH2F("hMassM02CutNLocMaxN","Invariant mass of splitted cluster with NLM>2 vs E, M02 cut, no TM",
                                         nptbins,ptmin,ptmax,mbins,mmin,mmax);
        fhMassM02CutNLocMaxN->SetYTitle("M (GeV/c^{2})");
        fhMassM02CutNLocMaxN->SetXTitle("E (GeV)");
        outputContainer->Add(fhMassM02CutNLocMaxN) ;
        
      }
      
      fhMassAfterCutsNLocMax1[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max = 1, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMax1[i][j]) ;
      
      fhMassAfterCutsNLocMax2[i][j]     = new TH2F(Form("hMassAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max = 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMax2[i][j]) ;
      
      
      fhMassAfterCutsNLocMaxN[i][j]     = new TH2F(Form("hMassAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                   Form("Mass vs E, %s %s, for N Local max > 2, M02 and asy cut",ptype[i].Data(),sMatched[j].Data()),
                                                   nptbins,ptmin,ptmax,mbins,mmin,mmax);
      fhMassAfterCutsNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassAfterCutsNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassAfterCutsNLocMaxN[i][j]) ;
            
      fhSplitEFractionAfterCutsNLocMax1[i][j]     = new TH2F(Form("hSplitEFractionAfterCutsNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                             Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 1, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                             nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionAfterCutsNLocMax1[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionAfterCutsNLocMax1[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionAfterCutsNLocMax1[i][j]) ;
      
      fhSplitEFractionAfterCutsNLocMax2[i][j]     = new TH2F(Form("hSplitEFractionAfterCutsNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                             Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  = 2, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                             nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionAfterCutsNLocMax2[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionAfterCutsNLocMax2[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionAfterCutsNLocMax2[i][j]) ;
      
      fhSplitEFractionAfterCutsNLocMaxN[i][j]    = new TH2F(Form("hSplitEFractionAfterCutsNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                            Form("(E1+E2)/E_{cluster} vs E_{cluster} for N max  > 2, M02 and Asy cut on, %s %s",ptype[i].Data(),sMatched[j].Data()),
                                                            nptbins,ptmin,ptmax,120,0,1.2);
      fhSplitEFractionAfterCutsNLocMaxN[i][j]   ->SetXTitle("E_{cluster} (GeV)");
      fhSplitEFractionAfterCutsNLocMaxN[i][j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
      outputContainer->Add(fhSplitEFractionAfterCutsNLocMaxN[i][j]) ;
      
      
      fhMassM02NLocMax1[i][j]  = new TH2F(Form("hMassM02NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=1, #lambda_{0}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                          ssbins,ssmin,ssmax,mbins,mmin,mmax); 
      fhMassM02NLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
      fhMassM02NLocMax1[i][j]->SetXTitle("#lambda_{0}^{2}");
      outputContainer->Add(fhMassM02NLocMax1[i][j]) ;   
      
      fhMassM02NLocMax2[i][j]  = new TH2F(Form("hMassM02NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Invariant mass of splitted cluster with NLM=2, #lambda_{0}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
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
      
      if(i==0 && j==0)
      {
        fhAsymM02CutNLocMax1  = new TH2F("hAsymM02CutNLocMax1","Asymmetry of NLM=1  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
        fhAsymM02CutNLocMax1->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
        fhAsymM02CutNLocMax1->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsymM02CutNLocMax1) ;
        
        fhAsymM02CutNLocMax2  = new TH2F("hAsymM02CutNLocMax2","Asymmetry of NLM=2  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
        fhAsymM02CutNLocMax2->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
        fhAsymM02CutNLocMax2->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsymM02CutNLocMax2) ;
        
        fhAsymM02CutNLocMaxN  = new TH2F("hAsymM02CutNLocMaxN","Asymmetry of NLM>2  vs cluster Energy, M02Cut, no TM", nptbins,ptmin,ptmax,200,-1,1);
        fhAsymM02CutNLocMaxN->SetYTitle("(E_{1}-E_{2})/(E_{1}+E_{2})");
        fhAsymM02CutNLocMaxN->SetXTitle("E (GeV)");
        outputContainer->Add(fhAsymM02CutNLocMaxN) ;
      }
      
      if(fFillSSExtraHisto)
      {
        fhMassDispEtaNLocMax1[i][j]  = new TH2F(Form("hMassDispEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=1, #sigma_{#eta #eta}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispEtaNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispEtaNLocMax1[i][j]->SetXTitle("#sigma_{#eta #eta}^{2}");
        outputContainer->Add(fhMassDispEtaNLocMax1[i][j]) ;   
        
        fhMassDispEtaNLocMax2[i][j]  = new TH2F(Form("hMassDispEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of splitted cluster with NLM=2 #sigma_{#eta #eta}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
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
                                                Form("Invariant mass of 2 highest energy cells #sigma_{#phi #phi}^{2}, E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                ssbins,ssmin,ssmax,mbins,mmin,mmax); 
        fhMassDispPhiNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispPhiNLocMax1[i][j]->SetXTitle("#sigma_{#phi #phi}^{2}");
        outputContainer->Add(fhMassDispPhiNLocMax1[i][j]) ;   
        
        fhMassDispPhiNLocMax2[i][j]  = new TH2F(Form("hMassDispPhiNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells #sigma_{#phi #phi}^{2}, E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
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
                                                Form("Invariant mass of 2 highest energy cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 12 GeV,%s %s",ptype[i].Data(),sMatched[j].Data()),
                                                200,-1,1,mbins,mmin,mmax); 
        fhMassDispAsyNLocMax1[i][j]->SetYTitle("M (GeV/c^{2})");
        fhMassDispAsyNLocMax1[i][j]->SetXTitle("A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2})");
        outputContainer->Add(fhMassDispAsyNLocMax1[i][j]) ;   
        
        fhMassDispAsyNLocMax2[i][j]  = new TH2F(Form("hMassDispAsyNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("Invariant mass of 2 local maxima cells A = (#sigma_{#phi #phi}^{2} - #sigma_{#eta #eta}^{2}) / (#sigma_{#phi #phi}^{2} + #sigma_{#eta #eta}^{2}), E > 12 GeV, %s %s",ptype[i].Data(),sMatched[j].Data()),
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
      
      
      // E vs centrality
      
      fhCentralityPi0NLocMax1[i][j]  = new TH2F(Form("hCentralityPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMax1[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMax1[i][j]) ;
      
      fhCentralityPi0NLocMax2[i][j]  = new TH2F(Form("hCentralityPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=2, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMax2[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMax2[i][j]) ;
      
      fhCentralityPi0NLocMaxN[i][j]  = new TH2F(Form("hCentralityPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM>1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityPi0NLocMaxN[i][j]->SetYTitle("Centrality");
      fhCentralityPi0NLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityPi0NLocMaxN[i][j]) ;
      
      fhCentralityEtaNLocMax1[i][j]  = new TH2F(Form("hCentralityEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax1[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMax1[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax1[i][j]) ;
      
      fhCentralityEtaNLocMax2[i][j]  = new TH2F(Form("hCentralityEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM=2, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMax2[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMax2[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMax2[i][j]) ;
      
      fhCentralityEtaNLocMaxN[i][j]  = new TH2F(Form("hCentralityEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                                Form("E vs Centrality, selected pi0 cluster with NLM>1, %s",ptype[i].Data()),
                                                nptbins,ptmin,ptmax,100,0,100);
      fhCentralityEtaNLocMaxN[i][j]->SetYTitle("Centrality");
      fhCentralityEtaNLocMaxN[i][j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhCentralityEtaNLocMaxN[i][j]) ;
      
      
      fhM02Pi0NLocMax1[i][j]     = new TH2F(Form("hM02Pi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02Pi0NLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax1[i][j]) ;
      
      fhM02EtaNLocMax1[i][j]     = new TH2F(Form("hM02EtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
      fhM02EtaNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMax1[i][j]) ;
      
      fhM02ConNLocMax1[i][j]    = new TH2F(Form("hM02ConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMax1[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMax1[i][j]) ; 
      
      fhM02Pi0NLocMax2[i][j]     = new TH2F(Form("hM02Pi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMax2[i][j]) ; 
      
      fhM02EtaNLocMax2[i][j]     = new TH2F(Form("hM02EtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMax2[i][j]) ;
      
      fhM02ConNLocMax2[i][j]    = new TH2F(Form("hM02ConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMax2[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMax2[i][j]) ; 
      
      fhM02Pi0NLocMaxN[i][j]     = new TH2F(Form("hM02Pi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02Pi0NLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02Pi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02Pi0NLocMaxN[i][j]) ; 
      
      fhM02EtaNLocMaxN[i][j]     = new TH2F(Form("hM02EtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02EtaNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02EtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02EtaNLocMaxN[i][j]) ; 
      
      fhM02ConNLocMaxN[i][j]    = new TH2F(Form("hM02ConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("#lambda_{0}^{2} vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
      fhM02ConNLocMaxN[i][j]   ->SetYTitle("#lambda_{0}^{2}");
      fhM02ConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhM02ConNLocMaxN[i][j]) ;
            
      
      fhMassPi0NLocMax1[i][j]     = new TH2F(Form("hMassPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax1[i][j]) ; 

      
      fhMassEtaNLocMax1[i][j]     = new TH2F(Form("hMassEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMax1[i][j]) ; 
      
      fhMassConNLocMax1[i][j]    = new TH2F(Form("hMassConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMax1[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMax1[i][j]) ; 
      
      fhMassPi0NLocMax2[i][j]     = new TH2F(Form("hMassPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMax2[i][j]) ; 
      
      
      fhMassEtaNLocMax2[i][j]     = new TH2F(Form("hMassEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMax2[i][j]) ; 
      
      fhMassConNLocMax2[i][j]    = new TH2F(Form("hMassConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMax2[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMax2[i][j]) ; 
      
      fhMassPi0NLocMaxN[i][j]     = new TH2F(Form("hMassPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassPi0NLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassPi0NLocMaxN[i][j]) ; 
      
      fhMassEtaNLocMaxN[i][j]     = new TH2F(Form("hMassEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Mass vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2", 
                                                GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassEtaNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassEtaNLocMaxN[i][j]) ; 
      
      fhMassConNLocMaxN[i][j]    = new TH2F(Form("hMassConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                          Form("Mass vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                               GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                          nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhMassConNLocMaxN[i][j]   ->SetYTitle("Mass (GeV/c^{2})");
      fhMassConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhMassConNLocMaxN[i][j]) ; 
      
      
      fhAsyPi0NLocMax1[i][j]     = new TH2F(Form("hAsyPi0NLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 1",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax1[i][j]) ; 
      
      fhAsyEtaNLocMax1[i][j]     = new TH2F(Form("hAsyEtaNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMax1[i][j]) ; 
      
      fhAsyConNLocMax1[i][j]    = new TH2F(Form("hAsyConNLocMax1%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 1",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMax1[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMax1[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMax1[i][j]) ; 
      
      fhAsyPi0NLocMax2[i][j]     = new TH2F(Form("hAsyPi0NLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max = 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMax2[i][j]) ; 
      
      fhAsyEtaNLocMax2[i][j]     = new TH2F(Form("hAsyEtaNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMax2[i][j]) ; 
      
      fhAsyConNLocMax2[i][j]    = new TH2F(Form("hAsyConNLocMax2%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max = 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMax2[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMax2[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMax2[i][j]) ; 
      
      fhAsyPi0NLocMaxN[i][j]     = new TH2F(Form("hAsyPi0NLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2} %s, for N Local max > 2",
                                                 GetCaloPID()->GetPi0MinMass(),GetCaloPID()->GetPi0MaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyPi0NLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyPi0NLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyPi0NLocMaxN[i][j]) ; 
      
      fhAsyEtaNLocMaxN[i][j]     = new TH2F(Form("hAsyEtaNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                            Form("Asymmetry vs E for mass range [%2.2f-%2.2f] GeV/c^{2}, %s, for N Local max > 2",
                                                 GetCaloPID()->GetEtaMinMass(),GetCaloPID()->GetEtaMaxMass(),ptype[i].Data()),
                                            nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyEtaNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyEtaNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyEtaNLocMaxN[i][j]) ; 
      
      fhAsyConNLocMaxN[i][j]    = new TH2F(Form("hAsyConNLocMaxN%s%s",pname[i].Data(),sMatched[j].Data()),
                                           Form("Asymmetry vs E for mass range [%2.2f-%2.2f], %s, for N Local max > 2",
                                                GetCaloPID()->GetPhotonMinMass(),GetCaloPID()->GetPhotonMaxMass(),ptype[i].Data()),
                                           nptbins,ptmin,ptmax,mbins,mmin,mmax); 
      fhAsyConNLocMaxN[i][j]   ->SetYTitle("Asymmetry");
      fhAsyConNLocMaxN[i][j]   ->SetXTitle("E (GeV)");
      outputContainer->Add(fhAsyConNLocMaxN[i][j]) ; 
      
    } // matched, not matched
    
    if(fFillEbinHisto)
    {
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
    }
  } // MC particle list
 
  // E vs Event plane angle
  
  fhEventPlanePi0NLocMax1  = new TH2F("hEventPlanePi0NLocMax1","E vs Event Plane Angle, selected pi0 cluster with NLM=1",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlanePi0NLocMax1->SetYTitle("Event Plane Angle (rad)");
  fhEventPlanePi0NLocMax1->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlanePi0NLocMax1) ;
  
  fhEventPlanePi0NLocMax2  = new TH2F("hEventPlanePi0NLocMax2","E vs Event Plane Angle, selected pi0 cluster with NLM=2",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlanePi0NLocMax2->SetYTitle("Event Plane Angle (rad)");
  fhEventPlanePi0NLocMax2->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlanePi0NLocMax2) ;
  
  fhEventPlanePi0NLocMaxN  = new TH2F("hEventPlanePi0NLocMaxN","E vs Event Plane Angle, selected pi0 cluster with NLM>1",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlanePi0NLocMaxN->SetYTitle("Event Plane Angle (rad)");
  fhEventPlanePi0NLocMaxN->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlanePi0NLocMaxN) ;
  
  fhEventPlaneEtaNLocMax1  = new TH2F("hEventPlaneEtaNLocMax1","E vs Event Plane Angle, selected pi0 cluster with NLM=1",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlaneEtaNLocMax1->SetYTitle("Event Plane Angle (rad)");
  fhEventPlaneEtaNLocMax1->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlaneEtaNLocMax1) ;
  
  fhEventPlaneEtaNLocMax2  = new TH2F("hEventPlaneEtaNLocMax2","E vs Event Plane Angle, selected pi0 cluster with NLM=2",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlaneEtaNLocMax2->SetYTitle("Event Plane Angle (rad)");
  fhEventPlaneEtaNLocMax2->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlaneEtaNLocMax2) ;
  
  fhEventPlaneEtaNLocMaxN  = new TH2F("hEventPlaneEtaNLocMaxN","E vs Event Plane Angle, selected pi0 cluster with NLM>1",
                                      nptbins,ptmin,ptmax,100,0,TMath::Pi());
  fhEventPlaneEtaNLocMaxN->SetYTitle("Event Plane Angle (rad)");
  fhEventPlaneEtaNLocMaxN->SetXTitle("E (GeV)");
  outputContainer->Add(fhEventPlaneEtaNLocMaxN) ;
  
  if(fFillEbinHisto)
  {
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
  }
    
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
  
  if(fFillTMResidualHisto && fFillTMHisto)
  {
    for(Int_t i = 0; i < n; i++)
    {  
      
      fhTrackMatchedDEtaNLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax1%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMax1[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax1[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax1%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 1 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMax1[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax1[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax1[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMax1[i]) ;
      
      fhTrackMatchedDEtaNLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMax2%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMax2[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMax2[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMax2%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, 2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMax2[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMax2[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMax2[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMax2[i]) ;
      
      fhTrackMatchedDEtaNLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDEtaNLocMaxN%s",pname[i].Data()),
       Form("d#eta of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax); 
      fhTrackMatchedDEtaNLocMaxN[i]->SetYTitle("d#eta");
      fhTrackMatchedDEtaNLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      fhTrackMatchedDPhiNLocMaxN[i]  = new TH2F
      (Form("hTrackMatchedDPhiNLocMaxN%s",pname[i].Data()),
       Form("d#phi of cluster-track vs cluster energy, N>2 Local Maxima, %s",ptype[i].Data()),
       nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax); 
      fhTrackMatchedDPhiNLocMaxN[i]->SetYTitle("d#phi (rad)");
      fhTrackMatchedDPhiNLocMaxN[i]->SetXTitle("E_{cluster} (GeV)");
      
      outputContainer->Add(fhTrackMatchedDEtaNLocMaxN[i]) ; 
      outputContainer->Add(fhTrackMatchedDPhiNLocMaxN[i]) ;    
    }
  }
  
  if(fFillAngleHisto)
  {
    for(Int_t j = 0; j < nMatched; j++)
    {  
      
      fhAnglePairNLocMax1[j]  = new TH2F(Form("hAnglePairNLocMax1%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 highest energy cells vs pair Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMax1[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMax1[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMax1[j]) ;   
      
      fhAnglePairNLocMax2[j]  = new TH2F(Form("hAnglePairNLocMax2%s",sMatched[j].Data()),
                                        Form("Opening angle of 2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMax2[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMax2[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMax2[j]) ;   
      
      fhAnglePairNLocMaxN[j]  = new TH2F(Form("hAnglePairNLocMaxN%s",sMatched[j].Data()),
                                        Form("Opening angle of N>2 local maxima cells vs Energy, %s",sMatched[j].Data()),
                                        nptbins,ptmin,ptmax,200,0,0.2); 
      fhAnglePairNLocMaxN[j]->SetYTitle("#alpha (rad)");
      fhAnglePairNLocMaxN[j]->SetXTitle("E (GeV)");
      outputContainer->Add(fhAnglePairNLocMaxN[j]) ;   
      
      fhAnglePairMassNLocMax1[j]  = new TH2F(Form("hAnglePairMassNLocMax1%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 highest energy cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMax1[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMax1[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMax1[j]) ;   
      
      fhAnglePairMassNLocMax2[j]  = new TH2F(Form("hAnglePairMassNLocMax2%s",sMatched[j].Data()),
                                            Form("Opening angle of 2 local maxima cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMax2[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMax2[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMax2[j]) ;   
      
      fhAnglePairMassNLocMaxN[j]  = new TH2F(Form("hAnglePairMassNLocMaxN%s",sMatched[j].Data()),
                                            Form("Opening angle of N>2 local maxima cells vs Mass for E > 12 GeV, %s",sMatched[j].Data()),
                                            mbins,mmin,mmax,200,0,0.2); 
      fhAnglePairMassNLocMaxN[j]->SetXTitle("M (GeV/c^{2})");
      fhAnglePairMassNLocMaxN[j]->SetYTitle("#alpha (rad)");
      outputContainer->Add(fhAnglePairMassNLocMaxN[j]) ;  
      
    }
  }
  
  for(Int_t j = 0; j < nMatched; j++)
  {
    fhSplitEFractionvsAsyNLocMax1[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax1%s",sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 1, E>12, %s",sMatched[j].Data()),
                                                    100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMax1[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMax1[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMax1[j]) ; 
    
    fhSplitEFractionvsAsyNLocMax2[j]     = new TH2F(Form("hSplitEFractionvsAsyNLocMax2%s",sMatched[j].Data()),
                                                    Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  = 2,E>12, %s",sMatched[j].Data()),
                                                    100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMax2[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMax2[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMax2[j]) ; 
    
    fhSplitEFractionvsAsyNLocMaxN[j]    = new TH2F(Form("hSplitEFractionvsAsyNLocMaxN%s",sMatched[j].Data()),
                                                   Form("(E1+E2)/E_{cluster} vs (E_{split1}-E_{split2})/(E_{split1}+E_{split2}) for N max  > 2, E>12, %s",sMatched[j].Data()),
                                                   100,-1,1,120,0,1.2); 
    fhSplitEFractionvsAsyNLocMaxN[j]   ->SetXTitle("(E_{split1}-E_{split2})/(E_{split1}+E_{split2})");
    fhSplitEFractionvsAsyNLocMaxN[j]   ->SetYTitle("(E_{split1}+E_{split2})/E_{cluster}");
    outputContainer->Add(fhSplitEFractionvsAsyNLocMaxN[j]) ; 
  }
  
  
  fhClusterEtaPhiNLocMax1  = new TH2F
  ("hClusterEtaPhiNLocMax1","Neutral Clusters with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMax1->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMax1->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMax1) ;

  fhClusterEtaPhiNLocMax2  = new TH2F
  ("hClusterEtaPhiNLocMax2","Neutral Clusters with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMax2->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMax2->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMax2) ;
  
  fhClusterEtaPhiNLocMaxN  = new TH2F
  ("hClusterEtaPhiNLocMaxN","Neutral Clusters with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhClusterEtaPhiNLocMaxN->SetYTitle("#phi (rad)");
  fhClusterEtaPhiNLocMaxN->SetXTitle("#eta");
  outputContainer->Add(fhClusterEtaPhiNLocMaxN) ;
  
  fhPi0EtaPhiNLocMax1  = new TH2F
  ("hPi0EtaPhiNLocMax1","Selected #pi^{0}'s with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMax1->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMax1->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMax1) ;
  
  fhPi0EtaPhiNLocMax2  = new TH2F
  ("hPi0EtaPhiNLocMax2","Selected #pi^{0}'s with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMax2->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMax2->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMax2) ;
  
  fhPi0EtaPhiNLocMaxN  = new TH2F
  ("hPi0EtaPhiNLocMaxN","Selected #pi^{0}'s with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhPi0EtaPhiNLocMaxN->SetYTitle("#phi (rad)");
  fhPi0EtaPhiNLocMaxN->SetXTitle("#eta");
  outputContainer->Add(fhPi0EtaPhiNLocMaxN) ;

  fhEtaEtaPhiNLocMax1  = new TH2F
  ("hEtaEtaPhiNLocMax1","Selected #eta's with E > 8 GeV, NLM = 1: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaEtaPhiNLocMax1->SetYTitle("#phi (rad)");
  fhEtaEtaPhiNLocMax1->SetXTitle("#eta");
  outputContainer->Add(fhEtaEtaPhiNLocMax1) ;
  
  fhEtaEtaPhiNLocMax2  = new TH2F
  ("hEtaEtaPhiNLocMax2","Selected #eta's with E > 8 GeV, NLM = 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaEtaPhiNLocMax2->SetYTitle("#phi (rad)");
  fhEtaEtaPhiNLocMax2->SetXTitle("#eta");
  outputContainer->Add(fhEtaEtaPhiNLocMax2) ;
  
  fhEtaEtaPhiNLocMaxN  = new TH2F
  ("hEtaEtaPhiNLocMaxN","Selected #eta's with E > 8 GeV, NLM > 2: #eta vs #phi",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaEtaPhiNLocMaxN->SetYTitle("#phi (rad)");
  fhEtaEtaPhiNLocMaxN->SetXTitle("#eta");
  outputContainer->Add(fhEtaEtaPhiNLocMaxN) ;

  
  if(fFillSSWeightHisto)
  {
    TString snlm[] = {"1","2","N"};
    for(Int_t nlm = 0; nlm < 3; nlm++)
    {
      fhPi0CellE[nlm]  = new TH2F(Form("hPi0CellENLocMax%s",snlm[nlm].Data()),
                                  Form("Selected #pi^{0}'s, NLM = %s: cluster E vs cell E",snlm[nlm].Data()),
                                  nptbins,ptmin,ptmax, nptbins,ptmin,ptmax);
      fhPi0CellE[nlm]->SetYTitle("E_{cell}");
      fhPi0CellE[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellE[nlm]) ;
      
      fhPi0CellEFrac[nlm]  = new TH2F(Form("hPi0CellEFracNLocMax%s",snlm[nlm].Data()),
                                      Form("Selected #pi^{0}'s, NLM = %s: cluster E vs cell E / cluster E",snlm[nlm].Data()),
                                      nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEFrac[nlm]->SetYTitle("E_{cell} / E_{cluster}");
      fhPi0CellEFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEFrac[nlm]) ;
      
      fhPi0CellLogEFrac[nlm]  = new TH2F(Form("hPi0CellLogEFracNLocMax%s",snlm[nlm].Data()),
                                      Form("Selected #pi^{0}'s, NLM = %s: cluster E vs Log(cell E / cluster E)",snlm[nlm].Data()),
                                      nptbins,ptmin,ptmax, 100,-10,0);
      fhPi0CellLogEFrac[nlm]->SetYTitle("Log(E_{cell} / E_{cluster})");
      fhPi0CellLogEFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellLogEFrac[nlm]) ;

      
      fhPi0CellEMaxEMax2Frac[nlm]  = new TH2F(Form("hPi0CellEMaxEMax2FracNLocMax%s",snlm[nlm].Data()),
                                      Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / 1st loc. max.  E",snlm[nlm].Data()),
                                      nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxEMax2Frac[nlm]->SetYTitle("E_{Loc Max 2} / E_{Loc Max 1}");
      fhPi0CellEMaxEMax2Frac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxEMax2Frac[nlm]) ;
      
      fhPi0CellEMaxClusterFrac[nlm]  = new TH2F(Form("hPi0CellEMaxClusterFracNLocMax%s",snlm[nlm].Data()),
                                              Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 1st loc. max. E / E cluster",snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxClusterFrac[nlm]->SetYTitle("E_{Loc Max 1} / E_{cluster}");
      fhPi0CellEMaxClusterFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxClusterFrac[nlm]) ;

      fhPi0CellEMax2ClusterFrac[nlm]  = new TH2F(Form("hPi0CellEMax2ClusterFracNLocMax%s",snlm[nlm].Data()),
                                                Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / E cluster",snlm[nlm].Data()),
                                                nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMax2ClusterFrac[nlm]->SetYTitle("E_{Loc Max 2} / E_{cluster}");
      fhPi0CellEMax2ClusterFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMax2ClusterFrac[nlm]) ;
      
      fhPi0CellEMaxFrac[nlm]  = new TH2F(Form("hPi0CellEMaxFracNLocMax%s",snlm[nlm].Data()),
                                                Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 1st loc. max. E / E cell i",snlm[nlm].Data()),
                                                nptbins,ptmin,ptmax, 100,0,1);
      fhPi0CellEMaxFrac[nlm]->SetYTitle("E_{Loc Max 1} / E_{cell i}");
      fhPi0CellEMaxFrac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMaxFrac[nlm]) ;
      
      fhPi0CellEMax2Frac[nlm]  = new TH2F(Form("hPi0CellEMax2FracNLocMax%s",snlm[nlm].Data()),
                                                 Form("Selected #pi^{0}'s, NLM = %s: cluster E vs 2nd loc. max. E / E cell i",snlm[nlm].Data()),
                                                 nptbins,ptmin,ptmax, 200,0,2);
      fhPi0CellEMax2Frac[nlm]->SetYTitle("E_{Loc Max 2} / E_{cell i}");
      fhPi0CellEMax2Frac[nlm]->SetXTitle("E_{cluster}");
      outputContainer->Add(fhPi0CellEMax2Frac[nlm]) ;

      
      for(Int_t i = 0; i < fSSWeightN; i++)
      {
        fhM02WeightPi0[nlm][i]     = new TH2F(Form("hM02Pi0NLocMax%s_W%d",snlm[nlm].Data(),i),
                                              Form("#lambda_{0}^{2} vs E, with W0 = %2.2f, for N Local max = %s", fSSWeight[i], snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02WeightPi0[nlm][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02WeightPi0[nlm][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02WeightPi0[nlm][i]) ;
      }
      
      for(Int_t i = 0; i < fSSECellCutN; i++)
      {
        fhM02ECellCutPi0[nlm][i]     = new TH2F(Form("hM02Pi0NLocMax%s_Ecell%d",snlm[nlm].Data(),i),
                                              Form("#lambda_{0}^{2} vs E, with Ecell > %2.2f, for N Local max = %s", fSSECellCut[i], snlm[nlm].Data()),
                                              nptbins,ptmin,ptmax,ssbins,ssmin,ssmax);
        fhM02ECellCutPi0[nlm][i]   ->SetYTitle("#lambda_{0}^{2}");
        fhM02ECellCutPi0[nlm][i]   ->SetXTitle("E (GeV)");
        outputContainer->Add(fhM02ECellCutPi0[nlm][i]) ;
      }

    }
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
  
  fSSWeightN   = 5;
  fSSWeight[0] = 4.6;  fSSWeight[1] = 4.7; fSSWeight[2] = 4.8; fSSWeight[3] = 4.9; fSSWeight[4] = 5.0;
  fSSWeight[5] = 5.1;  fSSWeight[6] = 5.2; fSSWeight[7] = 5.3; fSSWeight[8] = 5.4; fSSWeight[9] = 5.5;
  
  fSSECellCutN   = 10;
  fSSECellCut[0] = 0.16;  fSSECellCut[1] = 0.18; fSSECellCut[2] = 0.2; fSSECellCut[3] = 0.22; fSSECellCut[4] = 0.24;
  fSSECellCut[5] = 0.26;  fSSECellCut[6] = 0.28; fSSECellCut[7] = 0.3; fSSECellCut[8] = 0.32; fSSECellCut[9] = 0.34;

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

    Bool_t  matched   = IsTrackMatched(cluster,GetReader()->GetInputEvent());
    if(!fFillTMHisto && matched) continue ;
    
    // Study clusters with large shape parameter
    Float_t en = cluster->E();
    Float_t l0 = cluster->GetM02();
    Int_t   nc = cluster->GetNCells();
    Float_t bd = cluster->GetDistanceToBadChannel() ; 
    
    //If too small or big E or low number of cells, or close to a bad channel skip it
    if( en < GetMinEnergy() || en > GetMaxEnergy() || nc < fMinNCells || bd < fMinBadDist) continue ;
    
    TLorentzVector lv;
    cluster->GetMomentum(lv, GetVertex(0));
    Float_t eta = lv.Eta();
    Float_t phi = lv.Phi();
    if(phi<0) phi=+TMath::TwoPi();
    
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
    Int_t    absId1 = -1; Int_t absId2 = -1;

    Int_t pidTag = GetCaloPID()->GetIdentifiedParticleTypeFromClusterSplitting(cluster,cells,GetCaloUtils(),
                                                                               GetVertex(0), nMax, mass, angle,
                                                                               e1,e2,absId1,absId2);
    if (nMax <= 0) 
    {
      if(GetDebug() > 0 )
        printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found! It did not pass CaloPID selection criteria \n");
      
      return;
    }
    
    Float_t splitFrac = (e1+e2)/en;
    Float_t asym = -10;
    if(e1+e2>0) asym = (e1-e2)/(e1+e2);
        
    fhNLocMax[0][matched]->Fill(en,nMax);
    
    Int_t inlm = -1;
    if     (nMax == 1) inlm = 0;
    else if(nMax == 2) inlm = 1;
    else if(nMax > 2 ) inlm = 2;
    
    if     ( nMax == 1  ) 
    { 
      fhM02NLocMax1[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax1[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut)
      {
        fhSplitEFractionvsAsyNLocMax1[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMax1->Fill(eta,phi);
      }
      if(fFillSSExtraHisto) fhNCellNLocMax1[0][matched]->Fill(en,nc) ; 
    }
    else if( nMax == 2  ) 
    { 
      fhM02NLocMax2[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMax2[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut)
      {
        fhSplitEFractionvsAsyNLocMax2[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMax2->Fill(eta,phi);
      }
      if(fFillSSExtraHisto) fhNCellNLocMax2[0][matched]->Fill(en,nc) ; }
    else if( nMax >= 3  ) 
    { 
      fhM02NLocMaxN[0][matched]->Fill(en,l0) ; 
      fhSplitEFractionNLocMaxN[0][matched]->Fill(en,splitFrac) ; 
      if(en > ecut)
      {
        fhSplitEFractionvsAsyNLocMaxN[matched]->Fill(asym,splitFrac) ;
        if(!matched)fhClusterEtaPhiNLocMaxN->Fill(eta,phi);
      }
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
      if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1[0]->Fill(en,dR); }
      else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2[0]->Fill(en,dR); }
      else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxN[0]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxN[0]->Fill(en,dR); }
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
            
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0) &&
               !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = kmcPi0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = kmcPi0Conv;
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
        if     ( nMax == 1  ) { fhTrackMatchedDEtaNLocMax1[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax1[mcindex]->Fill(en,dR); }
        else if( nMax == 2  ) { fhTrackMatchedDEtaNLocMax2[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMax2[mcindex]->Fill(en,dR); }
        else if( nMax >= 3  ) { fhTrackMatchedDEtaNLocMaxN[mcindex]->Fill(en,dZ); fhTrackMatchedDPhiNLocMaxN[mcindex]->Fill(en,dR); }
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
    
    if(fFillEbinHisto && ebin >= 0 && IsDataMC() && fFillMCFractionHisto)
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
          
          if(!matched && ebin >= 0 && fFillEbinHisto)
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
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
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
          
          if(!matched && ebin >= 0 && fFillEbinHisto )
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
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
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
          
          if(!matched && ebin >= 0 && fFillEbinHisto)
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
      
      if(!matched && ebin >= 0 && fFillEbinHisto)
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
    
    Float_t cent = GetEventCentrality();
    Float_t evp  = GetEventPlaneAngle();
    
    fhNLocMaxM02Cut[0][matched]->Fill(en,nMax);
    if(IsDataMC()) fhNLocMaxM02Cut[mcindex][matched]->Fill(en,nMax);
    
    Float_t splitFracMin = GetCaloPID()->GetSplitEnergyFractionMinimum(inlm) ;
    
    if     (nMax==1) 
    { 
      fhMassNLocMax1[0][matched]->Fill(en,mass ); 
      fhAsymNLocMax1[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms

      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMax1->Fill(en,mass);
          fhAsymM02CutNLocMax1->Fill(en,asym );
          if(splitFrac > splitFracMin) fhMassSplitECutNLocMax1->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMax1[0][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin) fhMassAfterCutsNLocMax1[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          
          fhMCGenFracAfterCutsNLocMax1MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax1MCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto) 
      {
        fhAnglePairNLocMax1[matched]->Fill(en,angle);
      if( en > ecut ) 
        fhAnglePairMassNLocMax1[matched]->Fill(mass,angle);
      }
            
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax1[0][matched]->Fill(en,l0); fhMassConNLocMax1[0][matched]->Fill(en,mass);  fhAsyConNLocMax1[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMax1[0][matched]->Fill(en,l0); fhMassPi0NLocMax1[0][matched]->Fill(en,mass);  fhAsyPi0NLocMax1[0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMax1[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMax1->Fill(en,evp) ;
          if(en > ecut)fhPi0EtaPhiNLocMax1->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 0, absId1, absId2);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMax1[0][matched]->Fill(en,l0); fhMassEtaNLocMax1[0][matched]->Fill(en,mass);  fhAsyEtaNLocMax1[0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMax1[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMax1->Fill(en,evp) ;
          if(en > ecut)fhEtaEtaPhiNLocMax1->Fill(eta,phi);
        }
      }
      
    }
    else if(nMax==2) 
    {
      fhMassNLocMax2[0][matched]->Fill(en,mass );
      fhAsymNLocMax2[0][matched]->Fill(en,asym );
      
      // Effect of cuts in mass histograms
      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMax2->Fill(en,mass);
          fhAsymM02CutNLocMax2->Fill(en,asym );
          if(splitFrac > splitFracMin) fhMassSplitECutNLocMax2->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMax2[0][matched]->Fill(en,splitFrac);
        if(splitFrac >splitFracMin) fhMassAfterCutsNLocMax2[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          
          fhMCGenFracAfterCutsNLocMax2MCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMax2MCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto) 
      {
        fhAnglePairNLocMax2[matched]->Fill(en,angle);
        if( en > ecut ) 
          fhAnglePairMassNLocMax2[matched]->Fill(mass,angle);
      }
      
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax2[0][matched]->Fill(en,l0); fhMassConNLocMax2[0][matched]->Fill(en,mass);  fhAsyConNLocMax2[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMax2[0][matched]->Fill(en,l0); fhMassPi0NLocMax2[0][matched]->Fill(en,mass);  fhAsyPi0NLocMax2[0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMax2[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMax2->Fill(en,evp) ;
          if(en > ecut)fhPi0EtaPhiNLocMax2->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 1, absId1, absId2);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMax2[0][matched]->Fill(en,l0); fhMassEtaNLocMax2[0][matched]->Fill(en,mass);  fhAsyEtaNLocMax2[0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMax2[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMax2->Fill(en,evp) ;
          if(en > ecut)fhEtaEtaPhiNLocMax2->Fill(eta,phi);
        }
      }
      
    }
    else if(nMax >2) 
    {
      fhMassNLocMaxN[0][matched]->Fill(en,mass);
      fhAsymNLocMaxN[0][matched]->Fill(en,asym);
      
      // Effect of cuts in mass histograms
      if(!matched)
      {
        if(m02OK)
        {
          fhMassM02CutNLocMaxN->Fill(en,mass);
          fhAsymM02CutNLocMaxN->Fill(en,asym );
          if(splitFrac > splitFracMin)fhMassSplitECutNLocMaxN->Fill(en,mass );
        } // m02
      } // split frac
      
      if(m02OK && asyOK)
      {
        fhSplitEFractionAfterCutsNLocMaxN[0][matched]->Fill(en,splitFrac);
        if(splitFrac > splitFracMin) fhMassAfterCutsNLocMaxN[0][matched]->Fill(en,mass);
        
        if(!matched && IsDataMC() && fFillMCFractionHisto && mcindex==kmcPi0)
        {
          
          fhMCGenFracAfterCutsNLocMaxNMCPi0      ->Fill(en   ,  efrac     );
          fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0->Fill(en   ,  efracSplit);
        }
      }
      
      if(fFillAngleHisto)
      {
        fhAnglePairNLocMaxN[matched]->Fill(en,angle);
        if( en > ecut ) 
          fhAnglePairMassNLocMaxN[matched]->Fill(mass,angle);
      }
            
      if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMaxN[0][matched]->Fill(en,l0); fhMassConNLocMaxN[0][matched]->Fill(en,mass);  fhAsyConNLocMaxN[0][matched]->Fill(en,asym); }
      else if(pidTag==AliCaloPID::kPi0   )
      {
        fhM02Pi0NLocMaxN[0][matched]->Fill(en,l0); fhMassPi0NLocMaxN[0][matched]->Fill(en,mass);  fhAsyPi0NLocMaxN[0][matched]->Fill(en,asym);
        fhCentralityPi0NLocMaxN[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlanePi0NLocMaxN->Fill(en,evp) ;
          if(en > ecut)fhPi0EtaPhiNLocMaxN->Fill(eta,phi);
          FillSSWeightHistograms(cluster, 2,  absId1, absId2);
        }
      }
      else if(pidTag==AliCaloPID::kEta)
      {
        fhM02EtaNLocMaxN[0][matched]->Fill(en,l0); fhMassEtaNLocMaxN[0][matched]->Fill(en,mass);  fhAsyEtaNLocMaxN[0][matched]->Fill(en,asym);
        fhCentralityEtaNLocMaxN[0][matched]->Fill(en,cent) ;
        if(!matched)
        {
          fhEventPlaneEtaNLocMaxN->Fill(en,evp) ;
          if(en > ecut)fhEtaEtaPhiNLocMaxN->Fill(eta,phi);
        }
      }
      
    }
    
    
    if(IsDataMC())
    {
      if     (nMax==1) 
      { 
        fhMassNLocMax1[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax1[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMax1[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac > splitFracMin)
            fhMassAfterCutsNLocMax1[mcindex][matched]->Fill(en,mass);
        }

        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax1[mcindex][matched]->Fill(en,l0); fhMassConNLocMax1[mcindex][matched]->Fill(en,mass); fhAsyConNLocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMax1[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMax1[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMax1[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMax1[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMax1[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMax1[mcindex][matched]->Fill(en,asym); }
        
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMax1[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMax1[mcindex][matched]->Fill(en,cent) ;
      }
      else if(nMax==2) 
      {
        fhMassNLocMax2[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMax2[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMax2[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac >splitFracMin)
            fhMassAfterCutsNLocMax2[mcindex][matched]->Fill(en,mass);
        }
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMax2[mcindex][matched]->Fill(en,l0); fhMassConNLocMax2[mcindex][matched]->Fill(en,mass); fhAsyConNLocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMax2[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMax2[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMax2[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMax2[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMax2[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMax2[mcindex][matched]->Fill(en,asym); } 
       
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMax2[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMax2[mcindex][matched]->Fill(en,cent) ;
        
      }
      else if(nMax >2) 
      {
        fhMassNLocMaxN[mcindex][matched]->Fill(en,mass);
        fhAsymNLocMaxN[mcindex][matched]->Fill(en,asym);
        
        if(asyOK && m02OK)
        {
          fhSplitEFractionAfterCutsNLocMaxN[mcindex][matched]->Fill(en,splitFrac);
          if(splitFrac > splitFracMin )
            fhMassAfterCutsNLocMaxN[mcindex][matched]->Fill(en,mass);
        }
        
        if     (pidTag==AliCaloPID::kPhoton) { fhM02ConNLocMaxN[mcindex][matched]->Fill(en,l0); fhMassConNLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyConNLocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kPi0   ) { fhM02Pi0NLocMaxN[mcindex][matched]->Fill(en,l0); fhMassPi0NLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyPi0NLocMaxN[mcindex][matched]->Fill(en,asym); }
        else if(pidTag==AliCaloPID::kEta   ) { fhM02EtaNLocMaxN[mcindex][matched]->Fill(en,l0); fhMassEtaNLocMaxN[mcindex][matched]->Fill(en,mass); fhAsyEtaNLocMaxN[mcindex][matched]->Fill(en,asym); }
        
        if     (pidTag==AliCaloPID::kPi0) fhCentralityPi0NLocMaxN[mcindex][matched]->Fill(en,cent) ;
        else if(pidTag==AliCaloPID::kEta) fhCentralityEtaNLocMaxN[mcindex][matched]->Fill(en,cent) ;
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
  if(GetCaloUtils()) printf("Loc. Max. E > %2.2f\n",       GetCaloUtils()->GetLocalMaximaCutE());
  if(GetCaloUtils()) printf("Loc. Max. E Diff > %2.2f\n",  GetCaloUtils()->GetLocalMaximaCutEDiff());
  printf("Min. N Cells =%d \n",         fMinNCells) ;
  printf("Min. Dist. to Bad =%1.1f \n", fMinBadDist) ;
  printf("%2.2f < lambda_0^2 <%2.2f \n",fM02MinCut,fM02MaxCut);
  if(fFillSSWeightHisto) printf(" N w %d - N e cut %d \n",fSSWeightN,fSSECellCutN);

  printf("    \n") ;
  
} 


//___________________________________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::RecalculateClusterShowerShapeParametersWithCellCut(const AliEMCALGeometry * geom,
                                                                                          AliVCaloCells* cells,
                                                                                          AliVCluster * cluster,
                                                                                          Float_t & l0,   Float_t & l1,
                                                                                          Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                                                          Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi,
                                                                                          Float_t eCellMin)
{
  // Calculates new center of gravity in the local EMCAL-module coordinates
  // and tranfers into global ALICE coordinates
  // Calculates Dispersion and main axis
  
  if(!cluster)
  {
    AliInfo("Cluster pointer null!");
    return;
  }
  
  Double_t eCell       = 0.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  
  Int_t    iSupMod = -1;
  Int_t    iTower  = -1;
  Int_t    iIphi   = -1;
  Int_t    iIeta   = -1;
  Int_t    iphi    = -1;
  Int_t    ieta    = -1;
  Double_t etai    = -1.;
  Double_t phii    = -1.;
  
  Int_t    nstat   = 0 ;
  Float_t  wtot    = 0.;
  Double_t w       = 0.;
  Double_t etaMean = 0.;
  Double_t phiMean = 0.;
    
  //Loop on cells, calculate the cluster energy, in case a cut on cell energy is added
  // and to check if the cluster is between 2 SM in eta
  Int_t   iSM0   = -1;
  Bool_t  shared = kFALSE;
  Float_t energy = 0;

  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Check if there are cells of different SM
    if     (iDigit == 0   ) iSM0 = iSupMod;
    else if(iSupMod!= iSM0) shared = kTRUE;
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    //if(eCell > eCellMin)
      energy += eCell;
    
  }//cell loop
  
  //Loop on cells, get weighted parameters
  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if(energy > 0)// && eCell > eCellMin)
    {
      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);

      //correct weight, ONLY in simulation
      //w *= (1 - fWSimu * w );
      w *= (1 - eCellMin * w );

      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      
      if(w > 0.0)
      {
        wtot += w ;
        nstat++;
        //Shower shape
        sEta     += w * etai * etai ;
        etaMean  += w * etai ;
        sPhi     += w * phii * phii ;
        phiMean  += w * phii ;
        sEtaPhi  += w * etai * phii ;
      }
    }
    //else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, energy));
    
  }//cell loop
  
  //Normalize to the weight
  if (wtot > 0)
  {
    etaMean /= wtot ;
    phiMean /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));
  
  //Calculate dispersion
  for(Int_t iDigit=0; iDigit < cluster->GetNCells(); iDigit++)
  {
    //Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    //Get the cell energy, if recalibration is on, apply factors
    fraction  = cluster->GetCellAmplitudeFraction(iDigit);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    if (GetCaloUtils()->GetEMCALRecoUtils()->IsRecalibrationOn())
    {
      recalFactor = GetCaloUtils()->GetEMCALRecoUtils()->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cluster->GetCellAbsId(iDigit))*fraction*recalFactor;
    
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if(shared && iSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    if(energy > 0)// && eCell > eCellMin)
    {
      w  = GetCaloUtils()->GetEMCALRecoUtils()->GetCellWeight(eCell,energy);
      
      //correct weight, ONLY in simulation
      //w *= (1 - fWSimu * w );
      w *= (1 - eCellMin * w );

      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      if(w > 0.0)
      {
        disp +=  w *((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
        dEta +=  w * (etai-etaMean)*(etai-etaMean) ;
        dPhi +=  w * (phii-phiMean)*(phii-phiMean) ;
      }
    }
    //else if(energy == 0 || (eCellMin <0.01 && eCell == 0)) AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, energy));
  }// cell loop
  
  //Normalize to the weigth and set shower shape parameters
  if (wtot > 0 && nstat > 1)
  {
    disp    /= wtot ;
    dEta    /= wtot ;
    dPhi    /= wtot ;
    sEta    /= wtot ;
    sPhi    /= wtot ;
    sEtaPhi /= wtot ;
    
    sEta    -= etaMean * etaMean ;
    sPhi    -= phiMean * phiMean ;
    sEtaPhi -= etaMean * phiMean ;
    
    l0 = (0.5 * (sEta + sPhi) + TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    l1 = (0.5 * (sEta + sPhi) - TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
  }
  else
  {
    l0   = 0. ;
    l1   = 0. ;
    dEta = 0. ; dPhi = 0. ; disp    = 0. ;
    sEta = 0. ; sPhi = 0. ; sEtaPhi = 0. ;
  }
  
}


