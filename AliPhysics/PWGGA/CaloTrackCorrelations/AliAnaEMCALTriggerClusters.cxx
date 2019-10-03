/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// --- ROOT system ---
#include <TH2F.h>
#include <TClonesArray.h>
#include <TObjString.h>

// --- Analysis system ---
#include "AliAnaEMCALTriggerClusters.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"

// --- Detectors ---
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaEMCALTriggerClusters)
/// \endcond

//______________________________________________________
/// Default constructor. Initialize parameters.
//______________________________________________________
AliAnaEMCALTriggerClusters::AliAnaEMCALTriggerClusters() :
AliAnaCaloTrackCorrBaseClass(),
fRejectTrackMatch(0),         fNCellsCut(0),
fMinM02(0),                   fMaxM02(0),
fMomentum(),
// Histograms
fhE(0),                       fhESelected(0),
fhEtaPhi(0),                  fhEtaPhiSelected(0),
fhEtaPhiEMCALBC0(0),          fhEtaPhiEMCALBC1(0),          fhEtaPhiEMCALBCN(0),
fhTimeTriggerEMCALBCCluster(0),
fhTimeTriggerEMCALBCUMCluster(0),
fhEtaPhiTriggerEMCALBCClusterOverTh(0),
fhEtaPhiTriggerEMCALBCUMClusterOverTh(0),
fhEtaPhiTriggerEMCALBCClusterBelowTh1(0),
fhEtaPhiTriggerEMCALBCUMClusterBelowTh1(0),
fhEtaPhiTriggerEMCALBCClusterBelowTh2(0),
fhEtaPhiTriggerEMCALBCUMClusterBelowTh2(0),
fhEtaPhiTriggerEMCALBCExotic(0),             fhTimeTriggerEMCALBCExotic(0),
fhEtaPhiTriggerEMCALBCUMExotic(0),           fhTimeTriggerEMCALBCUMExotic(0),
fhEtaPhiTriggerEMCALBCBad(0),                fhTimeTriggerEMCALBCBad(0),
fhEtaPhiTriggerEMCALBCUMBad(0),              fhTimeTriggerEMCALBCUMBad(0),
fhEtaPhiTriggerEMCALBCBadExotic(0),          fhTimeTriggerEMCALBCBadExotic(0),
fhEtaPhiTriggerEMCALBCUMBadExotic(0),        fhTimeTriggerEMCALBCUMBadExotic(0),
fhEtaPhiTriggerEMCALBCExoticCluster(0),      fhTimeTriggerEMCALBCExoticCluster(0),
fhEtaPhiTriggerEMCALBCUMExoticCluster(0),    fhTimeTriggerEMCALBCUMExoticCluster(0),
fhEtaPhiTriggerEMCALBCBadCluster(0),         fhTimeTriggerEMCALBCBadCluster(0),
fhEtaPhiTriggerEMCALBCUMBadCluster(0),       fhTimeTriggerEMCALBCUMBadCluster(0),
fhEtaPhiTriggerEMCALBCBadExoticCluster(0),   fhTimeTriggerEMCALBCBadExoticCluster(0),
fhEtaPhiTriggerEMCALBCUMBadExoticCluster(0), fhTimeTriggerEMCALBCUMBadExoticCluster(0),
fhTimeTriggerEMCALBCBadMaxCell(0),           fhTimeTriggerEMCALBCUMBadMaxCell(0),
fhTimeTriggerEMCALBCBadMaxCellExotic(0),     fhTimeTriggerEMCALBCUMBadMaxCellExotic(0),
fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster (0), fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster(0),
fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster(0),fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster(0),
fhEtaPhiTriggerEMCALBCUMReMatchBothCluster(0),      fhTimeTriggerEMCALBCUMReMatchBothCluster(0),
fhTimeTriggerEMCALBC0UMReMatchOpenTime(0),
fhTimeTriggerEMCALBC0UMReMatchCheckNeigh(0),
fhTimeTriggerEMCALBC0UMReMatchBoth(0),
fhEtaPhiNoTrigger(0),                        fhTimeNoTrigger(0),
fhEtaPhiSelectedEMCALBC0(0),    fhEtaPhiSelectedEMCALBC1(0),   fhEtaPhiSelectedEMCALBCN(0),
fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime(0),
fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh(0),
fhTimeSelectedTriggerEMCALBC0UMReMatchBoth(0)
{
  for(Int_t i = 0; i < 11; i++)
  {
    fhEtaPhiTriggerEMCALBC             [i] = 0 ;
    fhTimeTriggerEMCALBC               [i] = 0 ;
    fhEtaPhiTriggerEMCALBCUM           [i] = 0 ;
    fhTimeTriggerEMCALBCUM             [i] = 0 ;
    
    fhEtaPhiSelectedTriggerEMCALBC       [i] = 0 ;
    fhTimeSelectedTriggerEMCALBC         [i] = 0 ;
    fhEtaPhiSelectedTriggerEMCALBCUM     [i] = 0 ;
    fhTimeSelectedTriggerEMCALBCUM       [i] = 0 ;
    
    fhTimeSelectedTriggerEMCALBCPileUpSPD[i] = 0 ;
    fhTimeTriggerEMCALBCPileUpSPD      [i] = 0 ;
    
    fhEtaPhiTriggerEMCALBCCluster      [i] = 0 ;
    fhEtaPhiTriggerEMCALBCUMCluster    [i] = 0 ;    
  }
  
  InitParameters();
}

//_____________________________________________________________
/// Fill Bad events histo, study bad/exotic trigger BC.
//_____________________________________________________________
void AliAnaEMCALTriggerClusters::FillBadTriggerEventHistogram()
{
  Int_t  idTrig = GetReader()->GetTriggerClusterIndex();
  Bool_t exotic = GetReader()->IsExoticEvent();
  Bool_t bad    = GetReader()->IsBadCellTriggerEvent();
  
  Bool_t ok = kFALSE;
  if(( bad || exotic )  && idTrig >= 0 && !GetReader()->AreBadTriggerEventsRemoved()) ok = kTRUE;
  
  if(!ok) return;
  
  //    printf("Index %d, Id %d,  bad %d, exo %d\n",
  //           GetReader()->GetTriggerClusterIndex(),
  //           GetReader()->GetTriggerClusterId(),
  //           GetReader()->IsBadCellTriggerEvent(),
  //           GetReader()->IsExoticEvent() );
  
  TClonesArray * clusterList = 0;
  TString  clusterListName   = GetReader()->GetEMCALClusterListName();
  if     (GetReader()->GetInputEvent()->FindListObject(clusterListName))
    clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetInputEvent() ->FindListObject(clusterListName));
  else if(GetReader()->GetOutputEvent())
    clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetOutputEvent()->FindListObject(clusterListName));
  
  AliVCluster  *  badClusTrig = 0;
  if(clusterList) badClusTrig = (AliVCluster*) clusterList->At(idTrig);
  else            badClusTrig = GetReader()->GetInputEvent()->GetCaloCluster(idTrig);
  
  if(!badClusTrig)
  {
    AliWarning(Form("No cluster (bad-exotic trigger) found with requested index %d \n",idTrig));
    return;
  }
  
  badClusTrig->GetMomentum(fMomentum,GetVertex(0));
  
  Float_t etaclusterBad = fMomentum.Eta();
  Float_t phiclusterBad = fMomentum.Phi();
  if( phiclusterBad < 0 ) phiclusterBad+=TMath::TwoPi();
  Float_t tofclusterBad = badClusTrig->GetTOF()*1.e9;
  Float_t eclusterBad   = badClusTrig->E();
  
  if( bad && exotic )
  {
    if(GetReader()->IsTriggerMatched())
    {
      fhEtaPhiTriggerEMCALBCBadExoticCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCBadExoticCluster  ->Fill(eclusterBad,   tofclusterBad, GetEventWeight());
    }
    else
    {
      fhEtaPhiTriggerEMCALBCUMBadExoticCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCUMBadExoticCluster  ->Fill(eclusterBad,   tofclusterBad, GetEventWeight());
    }
  }
  else if( bad && !exotic )
  {
    if(GetReader()->IsTriggerMatched())
    {
      fhEtaPhiTriggerEMCALBCBadCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCBadCluster  ->Fill(eclusterBad,   tofclusterBad, GetEventWeight());
    }
    else
    {
      fhEtaPhiTriggerEMCALBCUMBadCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCUMBadCluster  ->Fill(eclusterBad,   tofclusterBad, GetEventWeight());
    }
  }// Bad cluster trigger
  else if( !bad && exotic )
  {
    if(GetReader()->IsTriggerMatched())
    {
      fhEtaPhiTriggerEMCALBCExoticCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCExoticCluster  ->Fill(eclusterBad  , tofclusterBad, GetEventWeight());
    }
    else
    {
      fhEtaPhiTriggerEMCALBCUMExoticCluster->Fill(etaclusterBad, phiclusterBad, GetEventWeight());
      fhTimeTriggerEMCALBCUMExoticCluster  ->Fill(eclusterBad  , tofclusterBad, GetEventWeight());
    }
  }
}

//____________________________________________________________________________________________________________________________
/// Fill trigger related histograms.
//____________________________________________________________________________________________________________________________
void  AliAnaEMCALTriggerClusters::FillRawClusterTriggerBCHistograms(Int_t idcalo,       Float_t ecluster,  Float_t tofcluster,
                                                                    Float_t etacluster, Float_t phicluster)

{
  Float_t tofclusterUS = TMath::Abs(tofcluster);
  
  if(ecluster > 2)
  {
    if      (tofclusterUS < 25) fhEtaPhiEMCALBC0->Fill(etacluster, phicluster, GetEventWeight());
    else if (tofclusterUS < 75) fhEtaPhiEMCALBC1->Fill(etacluster, phicluster, GetEventWeight());
    else                        fhEtaPhiEMCALBCN->Fill(etacluster, phicluster, GetEventWeight());
  }
  
  Int_t  bc     = GetReader()->GetTriggerClusterBC();
  Int_t  id     = GetReader()->GetTriggerClusterId();
  Bool_t badMax = GetReader()->IsBadMaxCellTriggerEvent();
  
  Int_t histoBC = bc+5;
  if(GetReader()->AreBadTriggerEventsRemoved()) histoBC=0; // histograms created only for one BC since the others where rejected
  
  if(id==-2)
  {
    //printf("AliAnaEMCALTriggerClusters::ClusterSelected() - No trigger found bc=%d\n",bc);
    fhEtaPhiNoTrigger->Fill(etacluster, phicluster, GetEventWeight());
    fhTimeNoTrigger  ->Fill(ecluster  , tofcluster, GetEventWeight());
  }
  else if(TMath::Abs(bc) < 6)
  {
    if(!GetReader()->IsBadCellTriggerEvent() && !GetReader()->IsExoticEvent() )
    {
      if(GetReader()->IsTriggerMatched())
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBC       [histoBC]->Fill(etacluster, phicluster, GetEventWeight());
          
        fhTimeTriggerEMCALBC             [histoBC]->Fill(ecluster  , tofcluster, GetEventWeight());
          
        if(GetReader()->IsPileUpFromSPD())
            fhTimeTriggerEMCALBCPileUpSPD[histoBC]->Fill(ecluster  , tofcluster, GetEventWeight());
        
        if(idcalo ==  GetReader()->GetTriggerClusterId())
        {
          fhEtaPhiTriggerEMCALBCCluster[histoBC]->Fill(etacluster, phicluster, GetEventWeight());
          fhTimeTriggerEMCALBCCluster           ->Fill(ecluster  , tofcluster, GetEventWeight());
          
          if(bc==0)
          {
            Float_t threshold = GetReader()->GetEventTriggerL1Threshold() ;
            if(GetReader()->IsEventEMCALL0()) threshold = GetReader()->GetEventTriggerL0Threshold() ;
            
            if(ecluster > threshold)
              fhEtaPhiTriggerEMCALBCClusterOverTh  ->Fill(etacluster, phicluster, GetEventWeight());
            else if(ecluster > threshold-1)
              fhEtaPhiTriggerEMCALBCClusterBelowTh1->Fill(etacluster, phicluster, GetEventWeight());
            else
              fhEtaPhiTriggerEMCALBCClusterBelowTh2->Fill(etacluster, phicluster, GetEventWeight());
          }
        }
      }
      else
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCUM[histoBC]->Fill(etacluster, phicluster, GetEventWeight());
          
        fhTimeTriggerEMCALBCUM      [histoBC]->Fill(ecluster  , tofcluster, GetEventWeight());
        
        if(bc==0)
        {
          if(GetReader()->IsTriggerMatchedOpenCuts(0)) fhTimeTriggerEMCALBC0UMReMatchOpenTime   ->Fill(ecluster, tofcluster, GetEventWeight());
          if(GetReader()->IsTriggerMatchedOpenCuts(1)) fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ->Fill(ecluster, tofcluster, GetEventWeight());
          if(GetReader()->IsTriggerMatchedOpenCuts(2)) fhTimeTriggerEMCALBC0UMReMatchBoth       ->Fill(ecluster, tofcluster, GetEventWeight());
        }
        
        if(idcalo ==  GetReader()->GetTriggerClusterId())
        {
          fhEtaPhiTriggerEMCALBCUMCluster[histoBC]->Fill(etacluster, phicluster, GetEventWeight());
          fhTimeTriggerEMCALBCUMCluster           ->Fill(ecluster  , tofcluster, GetEventWeight());
          if(bc==0)
          {
            Float_t threshold = GetReader()->GetEventTriggerL1Threshold() ;
            if(GetReader()->IsEventEMCALL0()) threshold = GetReader()->GetEventTriggerL0Threshold() ;
            
            if(ecluster > threshold)
              fhEtaPhiTriggerEMCALBCUMClusterOverTh  ->Fill(etacluster, phicluster, GetEventWeight());
            else if(ecluster > threshold-1)
              fhEtaPhiTriggerEMCALBCUMClusterBelowTh1->Fill(etacluster, phicluster, GetEventWeight());
            else
              fhEtaPhiTriggerEMCALBCUMClusterBelowTh2->Fill(etacluster, phicluster, GetEventWeight());
            
            if(GetReader()->IsTriggerMatchedOpenCuts(0))
            {
              fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster->Fill(etacluster, phicluster, GetEventWeight());
              fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster  ->Fill(ecluster  , tofcluster, GetEventWeight());
            }
            if(GetReader()->IsTriggerMatchedOpenCuts(1))
            {
              fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster->Fill(etacluster, phicluster, GetEventWeight());
              fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster  ->Fill(ecluster  , tofcluster, GetEventWeight());
            }
            if(GetReader()->IsTriggerMatchedOpenCuts(2))
            {
              fhEtaPhiTriggerEMCALBCUMReMatchBothCluster->Fill(etacluster, phicluster, GetEventWeight());
              fhTimeTriggerEMCALBCUMReMatchBothCluster  ->Fill(ecluster  , tofcluster, GetEventWeight());
            }
            
          }
        }
      }
    }// neither bad nor exotic
    else if(GetReader()->IsBadCellTriggerEvent() && GetReader()->IsExoticEvent())
    {
      if(GetReader()->IsTriggerMatched())
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCBadExotic   ->Fill(etacluster, phicluster, GetEventWeight());
        
        fhTimeTriggerEMCALBCBadExotic         ->Fill(ecluster  , tofcluster, GetEventWeight());
        
        if(badMax)
          fhTimeTriggerEMCALBCBadMaxCellExotic->Fill(ecluster  , tofcluster, GetEventWeight());
      }
      else
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCUMBadExotic     ->Fill(etacluster, phicluster, GetEventWeight());
          
        fhTimeTriggerEMCALBCUMBadExotic           ->Fill(ecluster  , tofcluster, GetEventWeight());
          
        if(badMax)
            fhTimeTriggerEMCALBCUMBadMaxCellExotic->Fill(ecluster  , tofcluster, GetEventWeight());
        
      }
    }// Bad and exotic cluster trigger
    else if(GetReader()->IsBadCellTriggerEvent() )
    {
      if(GetReader()->IsTriggerMatched())
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCBad     ->Fill(etacluster, phicluster, GetEventWeight());
          
        fhTimeTriggerEMCALBCBad           ->Fill(ecluster  , tofcluster, GetEventWeight());
          
        if(badMax)
            fhTimeTriggerEMCALBCBadMaxCell->Fill(ecluster  , tofcluster, GetEventWeight());
      }
      else
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCUMBad     ->Fill(etacluster, phicluster, GetEventWeight());
          
        fhTimeTriggerEMCALBCUMBad           ->Fill(ecluster  , tofcluster, GetEventWeight());
          
        if(badMax)
            fhTimeTriggerEMCALBCUMBadMaxCell->Fill(ecluster  , tofcluster, GetEventWeight());
      }
    }// Bad cluster trigger
    else if(GetReader()->IsExoticEvent() )
    {
      if(GetReader()->IsTriggerMatched())
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCExotic->Fill(etacluster, phicluster, GetEventWeight());
        
        fhTimeTriggerEMCALBCExotic      ->Fill(ecluster  , tofcluster, GetEventWeight());
      }
      else
      {
        if(ecluster > 2)
            fhEtaPhiTriggerEMCALBCUMExotic->Fill(etacluster, phicluster, GetEventWeight());
        
        fhTimeTriggerEMCALBCUMExotic      ->Fill(ecluster  , tofcluster, GetEventWeight());
      }
    }
  }
  else if(TMath::Abs(bc) >= 6) AliWarning(Form("Trigger BC not expected = %d\n",bc));
}

//_________________________________________________________
/// Save parameters used for analysis.
//_________________________________________________________
TObjString *  AliAnaEMCALTriggerClusters::GetAnalysisCuts()
{
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaEMCALTriggerClusters ---:") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRejectTrackMatch: %d;",fRejectTrackMatch) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinM02: %2.2f, fMaxM02: %2.2f;",fMinM02,fMaxM02) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fNCellsCut: %d;",fNCellsCut) ;
  parList+=onePar ;
  
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  return new TObjString(parList) ;
}

//___________________________________________________________
/// Create histograms to be saved in output file and
/// store them in outputContainer.
//___________________________________________________________
TList *  AliAnaEMCALTriggerClusters::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName("EMCALTriggerClusters") ;
	
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();   Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();    Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins();  Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();   Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();  Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();   Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  Int_t ntimebins= GetHistogramRanges()->GetHistoTimeBins(); Float_t timemax = GetHistogramRanges()->GetHistoTimeMax(); Float_t timemin = GetHistogramRanges()->GetHistoTimeMin();
  
  Int_t nTrigBC  = 1;
  Int_t iBCShift = 0;
  if(!GetReader()->AreBadTriggerEventsRemoved())
  {
    nTrigBC = 11;
    iBCShift = 5;
  }
  
  fhE = new TH1F("hE","raw cluster #it{E}",nptbins,ptmin,ptmax);
  fhE->SetYTitle("d#it{N}/d#it{E} ");
  fhE->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhE) ;
  
  fhESelected = new TH1F("hESelected","selected cluster #it{E}",nptbins,ptmin,ptmax);
  fhESelected->SetYTitle("d#it{N}/d#it{E} ");
  fhESelected->SetXTitle("#it{E} (GeV)");
  outputContainer->Add(fhESelected) ;

  fhEtaPhi  = new TH2F
  ("hEtaPhi","cluster,#it{E} > 0.5 GeV, #eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhi->SetYTitle("#phi (rad)");
  fhEtaPhi->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhi) ;

  fhEtaPhiSelected  = new TH2F
  ("hEtaPhiSelected","selected cluster,#it{E} > 0.5 GeV, #eta vs #phi",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiSelected->SetYTitle("#phi (rad)");
  fhEtaPhiSelected->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiSelected) ;
  
  fhEtaPhiEMCALBC0  = new TH2F
  ("hEtaPhiEMCALBC0","cluster,#it{E} > 2 GeV, #eta vs #phi, for clusters with |time| < 25 ns, EMCAL-BC=0",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiEMCALBC0->SetYTitle("#phi (rad)");
  fhEtaPhiEMCALBC0->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiEMCALBC0) ;
  
  fhEtaPhiEMCALBC1  = new TH2F
  ("hEtaPhiEMCALBC1","cluster,#it{E} > 2 GeV, #eta vs #phi, for clusters with 25 < |time| < 75 ns, EMCAL-BC=1",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiEMCALBC1->SetYTitle("#phi (rad)");
  fhEtaPhiEMCALBC1->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiEMCALBC1) ;
  
  fhEtaPhiEMCALBCN  = new TH2F
  ("hEtaPhiEMCALBCN","cluster,#it{E} > 2 GeV, #eta vs #phi, for clusters with |time| > 75 ns, EMCAL-BC>1",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiEMCALBCN->SetYTitle("#phi (rad)");
  fhEtaPhiEMCALBCN->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiEMCALBCN) ;
  
  for(Int_t i = 0; i < nTrigBC; i++)
  {
    fhEtaPhiTriggerEMCALBC[i] = new TH2F
    (Form("hEtaPhiTriggerEMCALBC%d",i-iBCShift),
     Form("cluster #it{E} > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBC[i]->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBC[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBC[i]) ;
    
    fhTimeTriggerEMCALBC[i] = new TH2F
    (Form("hTimeTriggerEMCALBC%d",i-iBCShift),
     Form("cluster #it{time} vs #it{E} of clusters, Trigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBC[i]->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBC[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBC[i]);
    
    fhTimeTriggerEMCALBCPileUpSPD[i] = new TH2F
    (Form("hTimeTriggerEMCALBC%dPileUpSPD",i-iBCShift),
     Form("cluster #it{time} vs #it{E} of clusters, Trigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCPileUpSPD[i]->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCPileUpSPD[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCPileUpSPD[i]);
    
    fhEtaPhiTriggerEMCALBCUM[i] = new TH2F
    (Form("hEtaPhiTriggerEMCALBC%d_UnMatch",i-iBCShift),
     Form("cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUM[i]->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUM[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUM[i]) ;
    
    fhTimeTriggerEMCALBCUM[i] = new TH2F
    (Form("hTimeTriggerEMCALBC%d_UnMatch",i-iBCShift),
     Form("cluster #it{time} vs #it{E} of clusters, unmatched trigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUM[i]->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUM[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUM[i]);
    
    fhEtaPhiTriggerEMCALBCCluster[i] = new TH2F
    (Form("hEtaPhiTriggerEMCALBC%d_OnlyTrigger",i-iBCShift),
     Form("trigger cluster, #eta vs #phi, Trigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCCluster[i]->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCCluster[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCCluster[i]) ;
    
    fhEtaPhiTriggerEMCALBCUMCluster[i] = new TH2F
    (Form("hEtaPhiTriggerEMCALBC%d_OnlyTrigger_UnMatch",i-iBCShift),
     Form("trigger cluster, #eta vs #phi, unmatched trigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMCluster[i]->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMCluster[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMCluster[i]) ;
  }
  
  fhTimeTriggerEMCALBCCluster = new TH2F("hTimeTriggerEMCALBC_OnlyTrigger",
                                         "trigger cluster #it{time} vs #it{E} of clusters",
                                         nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBCCluster->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBCCluster->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBCCluster);
  
  fhTimeTriggerEMCALBCUMCluster = new TH2F("hTimeTriggerEMCALBC_OnlyTrigger_UnMatch",
                                           "trigger cluster #it{time} vs #it{E} of clusters, unmatched trigger",
                                           nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBCUMCluster->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBCUMCluster->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBCUMCluster);
  
  fhEtaPhiTriggerEMCALBCClusterOverTh = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_OverThreshold",
   "trigger cluster #it{E} > trigger threshold, #eta vs #phi, Trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCClusterOverTh->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCClusterOverTh->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCClusterOverTh) ;
  
  fhEtaPhiTriggerEMCALBCUMClusterOverTh = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_OverThreshold_UnMatch",
   "trigger cluster #it{E} > trigger threshold, #eta vs #phi, unmatched trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMClusterOverTh->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMClusterOverTh->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMClusterOverTh) ;
  
  fhEtaPhiTriggerEMCALBCClusterBelowTh1 = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_BelowThreshold1",
   "trigger cluster thresh-1 < #it{E} < thres, #eta vs #phi, Trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCClusterBelowTh1->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCClusterBelowTh1->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCClusterBelowTh1) ;
  
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh1 = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_BelowThreshold1_UnMatch",
   "trigger cluster thresh-1 < #it{E} < thres, #eta vs #phi, unmatched trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh1->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh1->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMClusterBelowTh1) ;
  
  fhEtaPhiTriggerEMCALBCClusterBelowTh2 = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_BelowThreshold2",
   "trigger cluster thresh-2 < #it{E} < thres, #eta vs #phi, Trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCClusterBelowTh2->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCClusterBelowTh2->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCClusterBelowTh2) ;
  
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh2 = new TH2F
  ("hEtaPhiTriggerEMCALBC0_OnlyTrigger_BelowThreshold2_UnMatch",
   "trigger cluster thresh-2 < #it{E} < thres, #eta vs #phi, unmatched trigger EMCAL-BC=0",
   netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh2->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMClusterBelowTh2->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMClusterBelowTh2) ;
  
  if(!GetReader()->AreBadTriggerEventsRemoved())
  {
    fhEtaPhiTriggerEMCALBCExotic = new TH2F
    ("hEtaPhiTriggerExotic",
     "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCExotic->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCExotic->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCExotic) ;
    
    fhTimeTriggerEMCALBCExotic = new TH2F
    ("hTimeTriggerExotic",
     "cluster #it{time} vs #it{E} of clusters, Trigger Exotic ",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCExotic);
    
    fhEtaPhiTriggerEMCALBCUMExotic = new TH2F
    ("hEtaPhiTriggerExotic_UnMatch",
     "cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMExotic->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMExotic->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMExotic) ;
    
    fhTimeTriggerEMCALBCUMExotic = new TH2F
    ("hTimeTriggerExotic_UnMatch",
     "cluster #it{time} vs #it{E} of clusters, unmatched trigger Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMExotic);
    
    fhEtaPhiTriggerEMCALBCExoticCluster = new TH2F
    ("hEtaPhiTriggerExotic_OnlyTrigger",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCExoticCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCExoticCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCExoticCluster) ;
    
    fhTimeTriggerEMCALBCExoticCluster = new TH2F
    ("hTimeTriggerExotic_OnlyTrigger",
     "trigger cluster #it{time} vs #it{E} of clusters, Trigger Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCExoticCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCExoticCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCExoticCluster);
    
    fhEtaPhiTriggerEMCALBCUMExoticCluster = new TH2F
    ("hEtaPhiTriggerExotic_OnlyTrigger_UnMatch",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMExoticCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMExoticCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMExoticCluster) ;
    
    fhTimeTriggerEMCALBCUMExoticCluster = new TH2F
    ("hTimeTriggerExotic_OnlyTrigger_UnMatch",
     "trigger cluster #it{time} vs #it{E} of clusters, unmatched trigger Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMExoticCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMExoticCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMExoticCluster);
    
    fhEtaPhiTriggerEMCALBCBad = new TH2F
    ("hEtaPhiTriggerBad",
     "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Bad",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCBad->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCBad->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCBad) ;
    
    fhTimeTriggerEMCALBCBad = new TH2F
    ("hTimeTriggerBad",
     "cluster #it{time} vs #it{E} of clusters, Trigger Bad ",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBad->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBad->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBad);
    
    fhEtaPhiTriggerEMCALBCUMBad = new TH2F
    ("hEtaPhiTriggerBad_UnMatch",
     "cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Bad",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMBad->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMBad->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMBad) ;
    
    fhTimeTriggerEMCALBCUMBad = new TH2F
    ("hTimeTriggerBad_UnMatch",
     "cluster #it{time} vs #it{E} of clusters, unmatched trigger Bad",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBad->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBad->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBad);
    
    fhEtaPhiTriggerEMCALBCBadCluster = new TH2F
    ("hEtaPhiTriggerBad_OnlyTrigger",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Bad",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCBadCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCBadCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCBadCluster) ;
    
    fhTimeTriggerEMCALBCBadCluster = new TH2F
    ("hTimeTriggerBad_OnlyTrigger",
     "trigger cluster #it{time} vs #it{E} of clusters, Trigger Bad",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBadCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBadCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBadCluster);
    
    fhEtaPhiTriggerEMCALBCUMBadCluster = new TH2F
    ("hEtaPhiTriggerBad_OnlyTrigger_UnMatch",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Bad",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMBadCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMBadCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMBadCluster) ;
    
    fhTimeTriggerEMCALBCUMBadCluster = new TH2F
    ("hTimeTriggerBad_OnlyTrigger_UnMatch",
     "trigger cluster time vs #it{E} of clusters, unmatched trigger Bad",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBadCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBadCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBadCluster);
    
    fhEtaPhiTriggerEMCALBCBadExotic = new TH2F
    ("hEtaPhiTriggerBadExotic",
     "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Bad&Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCBadExotic->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCBadExotic->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCBadExotic) ;
    
    fhTimeTriggerEMCALBCBadExotic = new TH2F
    ("hTimeTriggerBadExotic",
     "cluster #it{time} vs #it{E} of clusters, Trigger Bad&Exotic ",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBadExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBadExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBadExotic);
    
    fhEtaPhiTriggerEMCALBCUMBadExotic = new TH2F
    ("hEtaPhiTriggerBadExotic_UnMatch",
     "cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Bad&Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMBadExotic->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMBadExotic->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMBadExotic) ;
    
    fhTimeTriggerEMCALBCUMBadExotic = new TH2F
    ("hTimeTriggerBadExotic_UnMatch",
     "cluster #it{time} vs #it{E} of clusters, unmatched trigger Bad&Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBadExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBadExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBadExotic);
    
    fhEtaPhiTriggerEMCALBCBadExoticCluster = new TH2F
    ("hEtaPhiTriggerBadExotic_OnlyTrigger",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, Trigger Bad&Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCBadExoticCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCBadExoticCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCBadExoticCluster) ;
    
    fhTimeTriggerEMCALBCBadExoticCluster = new TH2F
    ("hTimeTriggerBadExotic_OnlyTrigger",
     "trigger cluster #it{time} vs #it{E} of clusters, Trigger Bad&Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBadExoticCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBadExoticCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBadExoticCluster);
    
    fhEtaPhiTriggerEMCALBCUMBadExoticCluster = new TH2F
    ("hEtaPhiTriggerBadExotic_OnlyTrigger_UnMatch",
     "trigger cluster #it{E} > 2 GeV, #eta vs #phi, unmatched trigger Bad&Exotic",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiTriggerEMCALBCUMBadExoticCluster->SetYTitle("#phi (rad)");
    fhEtaPhiTriggerEMCALBCUMBadExoticCluster->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiTriggerEMCALBCUMBadExoticCluster) ;
    
    fhTimeTriggerEMCALBCUMBadExoticCluster = new TH2F
    ("hTimeTriggerBadExotic_OnlyTrigger_UnMatch",
     "trigger cluster #it{time} vs #it{E} of clusters, unmatched trigger Bad&Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBadExoticCluster->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBadExoticCluster->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBadExoticCluster);
    
    fhTimeTriggerEMCALBCBadMaxCell = new TH2F
    ("hTimeTriggerBadMaxCell",
     "cluster #it{time} vs #it{E} of clusters, Trigger BadMaxCell",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBadMaxCell->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBadMaxCell->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBadMaxCell);
    
    fhTimeTriggerEMCALBCUMBadMaxCell = new TH2F
    ("hTimeTriggerBadMaxCell_UnMatch",
     "cluster #it{time} vs #it{E} of clusters, unmatched trigger BadMaxCell",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBadMaxCell->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBadMaxCell->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBadMaxCell);
    
    
    fhTimeTriggerEMCALBCBadMaxCellExotic = new TH2F
    ("hTimeTriggerBadMaxCellExotic",
     "cluster #it{time} vs #it{E} of clusters, Trigger BadMaxCell&Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCBadMaxCellExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCBadMaxCellExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCBadMaxCellExotic);
    
    fhTimeTriggerEMCALBCUMBadMaxCellExotic = new TH2F
    ("hTimeTriggerBadMaxCellExotic_UnMatch",
     "cluster #it{time} vs #it{E} of clusters, unmatched trigger BadMaxCell&Exotic",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeTriggerEMCALBCUMBadMaxCellExotic->SetXTitle("#it{E} (GeV)");
    fhTimeTriggerEMCALBCUMBadMaxCellExotic->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeTriggerEMCALBCUMBadMaxCellExotic);
    
    fhTimeNoTrigger = new TH2F
    ("hTimeNoTrigger",
     "events with no foundable trigger, time vs e of clusters",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeNoTrigger->SetXTitle("#it{E} (GeV)");
    fhTimeNoTrigger->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeNoTrigger);
    
    fhEtaPhiNoTrigger = new TH2F
    ("hEtaPhiNoTrigger",
     "events with no foundable trigger, eta vs phi of clusters",
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiNoTrigger->SetYTitle("#phi (rad)");
    fhEtaPhiNoTrigger->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiNoTrigger) ;
  }
  
  fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster = new TH2F("hEtaPhiTriggerEMCALBC0_OnlyTrigger_UnMatch_ReMatch_OpenTime",
                                                            "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=0, un match, rematch open time",
                                                            netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster) ;
  
  fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster = new TH2F("hTimeTrigger_OnlyTrigger_UnMatch_ReMatch_OpenTime",
                                                          "cluster #it{time} vs #it{E} of clusters, no match, rematch open time",
                                                          nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster);
  
  
  fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster = new TH2F("hEtaPhiTriggerEMCALBC0_OnlyTrigger_UnMatch_ReMatch_CheckNeighbours",
                                                              "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=0, un match, rematch with neighbour patches",
                                                              netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster) ;
  
  fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster = new TH2F("hTimeTrigger_OnlyTrigger_UnMatch_ReMatch_CheckNeighbours",
                                                            "cluster #it{time} vs #it{E} of clusters, no match, rematch with neigbour parches",
                                                            nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster);
  
  fhEtaPhiTriggerEMCALBCUMReMatchBothCluster = new TH2F("hEtaPhiTriggerEMCALBC0_OnlyTrigger_UnMatch_ReMatch_Both",
                                                        "cluster #it{E} > 2 GeV, #eta vs #phi, Trigger EMCAL-BC=0, un match, rematch open time and neighbour",
                                                        netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiTriggerEMCALBCUMReMatchBothCluster->SetYTitle("#phi (rad)");
  fhEtaPhiTriggerEMCALBCUMReMatchBothCluster->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiTriggerEMCALBCUMReMatchBothCluster) ;
  
  fhTimeTriggerEMCALBCUMReMatchBothCluster = new TH2F("hTimeTrigger_OnlyTrigger_UnMatch_ReMatch_Both",
                                                      "cluster #it{time} vs #it{E} of clusters, no match, rematch open time and neigbour",
                                                      nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBCUMReMatchBothCluster->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBCUMReMatchBothCluster->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBCUMReMatchBothCluster);
  
  fhTimeTriggerEMCALBC0UMReMatchOpenTime = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_OpenTime",
                                                    "cluster #it{time} vs #it{E} of clusters, no match, rematch open time",
                                                    nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBC0UMReMatchOpenTime->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchOpenTime);
  
  
  fhTimeTriggerEMCALBC0UMReMatchCheckNeigh = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_CheckNeighbours",
                                                      "cluster #it{time} vs #it{E} of clusters, no match, rematch with neigbour parches",
                                                      nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBC0UMReMatchCheckNeigh->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchCheckNeigh);
  
  fhTimeTriggerEMCALBC0UMReMatchBoth = new TH2F("hTimeTriggerBC0_UnMatch_ReMatch_Both",
                                                "cluster #it{time} vs #it{E} of clusters, no match, rematch open time and neigbour",
                                                nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeTriggerEMCALBC0UMReMatchBoth->SetXTitle("#it{E} (GeV)");
  fhTimeTriggerEMCALBC0UMReMatchBoth->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeTriggerEMCALBC0UMReMatchBoth);
  
  fhEtaPhiSelectedEMCALBC0  = new TH2F
  ("hEtaPhiSelectedEMCALBC0","Selected, #it{E} > 2 GeV, #eta vs #phi, for clusters with |time| < 25 ns, EMCAL-BC=0",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiSelectedEMCALBC0->SetYTitle("#phi (rad)");
  fhEtaPhiSelectedEMCALBC0->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiSelectedEMCALBC0) ;
  
  fhEtaPhiSelectedEMCALBC1  = new TH2F
  ("hEtaPhiSelectedEMCALBC1","Selected, #it{E} > 2 GeV, #eta vs #phi, for clusters with 25 < |time| < 75 ns, EMCAL-BC=1",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiSelectedEMCALBC1->SetYTitle("#phi (rad)");
  fhEtaPhiSelectedEMCALBC1->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiSelectedEMCALBC1) ;
  
  fhEtaPhiSelectedEMCALBCN  = new TH2F
  ("hEtaPhiSelectedEMCALBCN","Selected, #it{E} > 2 GeV, #eta vs #phi, for clusters with |time| > 75 ns, EMCAL-BC>1",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhEtaPhiSelectedEMCALBCN->SetYTitle("#phi (rad)");
  fhEtaPhiSelectedEMCALBCN->SetXTitle("#eta");
  outputContainer->Add(fhEtaPhiSelectedEMCALBCN) ;
  
  for(Int_t i = 0; i < nTrigBC; i++)
  {
    fhEtaPhiSelectedTriggerEMCALBC[i] = new TH2F
    (Form("hEtaPhiSelectedTriggerEMCALBC%d",i-iBCShift),
     Form("photon #it{E} > 2 GeV, #eta vs #phi, SelectedTrigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiSelectedTriggerEMCALBC[i]->SetYTitle("#phi (rad)");
    fhEtaPhiSelectedTriggerEMCALBC[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiSelectedTriggerEMCALBC[i]) ;
    
    fhTimeSelectedTriggerEMCALBC[i] = new TH2F
    (Form("hTimeSelectedTriggerEMCALBC%d",i-iBCShift),
     Form("photon #it{time} vs #it{E} of clusters, SelectedTrigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeSelectedTriggerEMCALBC[i]->SetXTitle("#it{E} (GeV)");
    fhTimeSelectedTriggerEMCALBC[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeSelectedTriggerEMCALBC[i]);
    
    fhTimeSelectedTriggerEMCALBCPileUpSPD[i] = new TH2F
    (Form("hTimeSelectedTriggerEMCALBC%dPileUpSPD",i-iBCShift),
     Form("photon #it{time} vs #it{E}, SelectedTrigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeSelectedTriggerEMCALBCPileUpSPD[i]->SetXTitle("#it{E} (GeV)");
    fhTimeSelectedTriggerEMCALBCPileUpSPD[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeSelectedTriggerEMCALBCPileUpSPD[i]);
    
    fhEtaPhiSelectedTriggerEMCALBCUM[i] = new TH2F
    (Form("hEtaPhiSelectedTriggerEMCALBC%d_UnMatch",i-iBCShift),
     Form("photon #it{E} > 2 GeV, #eta vs #phi, unmatched trigger EMCAL-BC=%d",i-iBCShift),
     netabins,etamin,etamax,nphibins,phimin,phimax);
    fhEtaPhiSelectedTriggerEMCALBCUM[i]->SetYTitle("#phi (rad)");
    fhEtaPhiSelectedTriggerEMCALBCUM[i]->SetXTitle("#eta");
    outputContainer->Add(fhEtaPhiSelectedTriggerEMCALBCUM[i]) ;
    
    fhTimeSelectedTriggerEMCALBCUM[i] = new TH2F
    (Form("hTimeSelectedTriggerEMCALBC%d_UnMatch",i-iBCShift),
     Form("photon #it{time} vs #it{E}, unmatched trigger EMCAL-BC=%d",i-iBCShift),
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
    fhTimeSelectedTriggerEMCALBCUM[i]->SetXTitle("#it{E} (GeV)");
    fhTimeSelectedTriggerEMCALBCUM[i]->SetYTitle("#it{time} (ns)");
    outputContainer->Add(fhTimeSelectedTriggerEMCALBCUM[i]);
  }
  
  fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime = new TH2F("hTimeSelectedTriggerBC0_UnMatch_ReMatch_OpenTime",
                                                          "cluster #it{time} vs #it{E} of photons, no match, rematch open time",
                                                          nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime->SetXTitle("#it{E} (GeV)");
  fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime);
  
  
  fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh = new TH2F("hTimeSelectedTriggerBC0_UnMatch_ReMatch_CheckNeighbours",
                                                            "cluster #it{time} vs #it{E} of photons, no match, rematch with neigbour parches",
                                                            nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh->SetXTitle("#it{E} (GeV)");
  fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh);
  
  fhTimeSelectedTriggerEMCALBC0UMReMatchBoth = new TH2F("hTimeSelectedTriggerBC0_UnMatch_ReMatch_Both",
                                                      "cluster #it{time} vs #it{E} of photons, no match, rematch open time and neigbour",
                                                      nptbins,ptmin,ptmax, ntimebins,timemin,timemax);
  fhTimeSelectedTriggerEMCALBC0UMReMatchBoth->SetXTitle("#it{E} (GeV)");
  fhTimeSelectedTriggerEMCALBC0UMReMatchBoth->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeSelectedTriggerEMCALBC0UMReMatchBoth);
  
  return outputContainer ;
}

//_____________________________________
/// Init. Check that EMCal data is recovered
/// and it is no MC data, if not abort.
//_____________________________________
void AliAnaEMCALTriggerClusters::Init()
{
  if(!GetReader()->IsEMCALSwitchedOn() || GetReader()->GetDataType() == AliCaloTrackReader::kMC)
  {
    AliFatal("You want to use EMCAL real data in analysis but it is not read!! \n!!Check the configuration file!!\n");
  }
}

//_______________________________________________
/// Initialize the parameters of the analysis with default values.
//_______________________________________________
void AliAnaEMCALTriggerClusters::InitParameters()
{
  AddToHistogramsName("AnaEMCALTriggerCluster_");
	
  fRejectTrackMatch = kTRUE;
  fMinM02           = 0.1;
  fMaxM02           = 0.3;
  fNCellsCut        = 2;
}

//____________________________________________________________
/// Main method. Fill some histograms with cluster kinematics
/// or time dependent on the trigger.
//____________________________________________________________
void  AliAnaEMCALTriggerClusters::MakeAnalysisFillHistograms()
{
  TObjArray * pl = GetEMCALClusters();
  
  if(!pl)
  {
    AliWarning("TObjArray with clusters is NULL!");
    return;
  }
  
  FillBadTriggerEventHistogram();

  // Loop on raw clusters before filtering in the reader and fill control histogram
  
  Int_t nCaloClusters = pl->GetEntriesFast();
  Int_t idTrig        = GetReader()->GetTriggerClusterIndex();

  AliDebug(1,Form("Input cluster entries %d", nCaloClusters));
  
  // Loop on clusters
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));
    //printf("calo %d, %f\n",icalo,calo->E());
    
    calo->GetMomentum(fMomentum,GetVertex(0)) ;
    
    Float_t tofcluster = calo->GetTOF()*1.e9;
    Float_t ecluster   = fMomentum.E();
    Float_t etacluster = fMomentum.Eta();
    Float_t phicluster = fMomentum.Phi();
    if(phicluster < 0) phicluster+=TMath::TwoPi();
    
    FillRawClusterTriggerBCHistograms(calo->GetID(),ecluster,tofcluster,etacluster,phicluster);

    // Select clusters
    
    if(idTrig < 0) continue;
    
    fhE->Fill(ecluster, GetEventWeight());
    if(ecluster > 0.5) fhEtaPhi->Fill(etacluster, phicluster, GetEventWeight());
    
    //.......................................
    //If too small or big energy, skip it
    if(ecluster < GetMinEnergy() || ecluster > GetMaxEnergy() ) continue ;
    
    //.......................................
    if(calo->GetNCells() <= fNCellsCut) continue;
    
    //.......................................
    //Check acceptance selection
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(fMomentum.Eta(),fMomentum.Phi(),kEMCAL) ;
      if(! in ) continue ;
    }
    
    //.......................................
    //Skip matched clusters with tracks
    if(fRejectTrackMatch && IsTrackMatched(calo,GetReader()->GetInputEvent())) continue;
    
    //.......................................
    //Skip matched clusters with Large shower shape
    if(calo->GetM02() < fMinM02 || calo->GetM02() > fMaxM02) continue;
    
    fhESelected ->Fill(ecluster, GetEventWeight());
    if(ecluster > 0.5) fhEtaPhiSelected->Fill(etacluster, phicluster, GetEventWeight());
    
    Float_t  tofUS = TMath::Abs(tofcluster);
    
    if(calo->E() > 2)
    {
      if      (tofUS < 25) fhEtaPhiSelectedEMCALBC0->Fill(etacluster, phicluster, GetEventWeight());
      else if (tofUS < 75) fhEtaPhiSelectedEMCALBC1->Fill(etacluster, phicluster, GetEventWeight());
      else                 fhEtaPhiSelectedEMCALBCN->Fill(etacluster, phicluster, GetEventWeight());
    }
    
    Int_t bc = GetReader()->GetTriggerClusterBC();
    Int_t histoBC = bc-5;
    if(GetReader()->AreBadTriggerEventsRemoved()) histoBC = 0 ; // histograms created only for one BC since the others where rejected
    
    if(TMath::Abs(bc) < 6 && !GetReader()->IsBadCellTriggerEvent() && !GetReader()->IsExoticEvent())
    {
      if(GetReader()->IsTriggerMatched())
      {
        if(calo->E() > 2) fhEtaPhiSelectedTriggerEMCALBC[histoBC]->Fill(etacluster, phicluster, GetEventWeight());
        fhTimeSelectedTriggerEMCALBC[histoBC]->Fill(ecluster, tofcluster, GetEventWeight());
        if(GetReader()->IsPileUpFromSPD()) fhTimeSelectedTriggerEMCALBCPileUpSPD[histoBC]->Fill(ecluster, tofcluster, GetEventWeight());
      }
      else
      {
        if(calo->E() > 2) fhEtaPhiSelectedTriggerEMCALBCUM[histoBC]->Fill(etacluster, phicluster, GetEventWeight());
        fhTimeSelectedTriggerEMCALBCUM[histoBC]->Fill(calo->E(), tofcluster, GetEventWeight());
        
        if(bc==0)
        {
          if(GetReader()->IsTriggerMatchedOpenCuts(0)) fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime   ->Fill(ecluster, tofcluster, GetEventWeight());
          if(GetReader()->IsTriggerMatchedOpenCuts(1)) fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh ->Fill(ecluster, tofcluster, GetEventWeight());
          if(GetReader()->IsTriggerMatchedOpenCuts(2)) fhTimeSelectedTriggerEMCALBC0UMReMatchBoth       ->Fill(ecluster, tofcluster, GetEventWeight());
        }
      }
    }
    else if(TMath::Abs(bc) >= 6)
      AliWarning(Form("Trigger BC not expected = %d",bc));
    
  }// cluster loop
  
  AliDebug(1,"End fill histograms");
}

//________________________________________________________________
/// Print some relevant parameters set for the analysis.
//________________________________________________________________
void AliAnaEMCALTriggerClusters::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  //AliAnaCaloTrackCorrBaseClass::Print(" ");
  printf("Reject clusters with a track matched = %d\n", fRejectTrackMatch);
  printf("M02 Cut: %2.2f < m02  < %2.2f\n"            , fMinM02,fMaxM02);
  printf("Number of cells in cluster is > %d \n"      , fNCellsCut);
}
