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

//------------------------------------------------------------------------//
//  Fill histograms with basic QA information for EMCAL offline trigger   //
//  Author: Nicolas Arbor (LPSC-Grenoble), Rachid Guernane (LPSC-Grenoble)//
//          Gustavo Conesa Balbastre  (LPSC-Grenoble)                     //
//                                                                        //
//  $Id$ //
//------------------------------------------------------------------------//

#include <TList.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TStreamerInfo.h>
#include <TFile.h>

#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"
#include "AliCentrality.h"

#include "AliVEvent.h"
//#include "AliVVZERO.h"
#include "AliVCaloTrigger.h"

//#include "AliESDEvent.h"
#include "AliESDVZERO.h"
//#include "AliESDCaloTrigger.h"

#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskEMCALTriggerQA.h"

ClassImp(AliAnalysisTaskEMCALTriggerQA)

//______________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA() :
AliAnalysisTaskSE(),
fOutputList(0),            fRecoUtils(0x0),
fGeoSet(0),                fGeometry(0),         fGeoName(""),
fOADBSet(kFALSE),          fAccessOADB(kTRUE),   fOADBFilePath(""),
fBitEGA(0),                fBitEJE(0),
fEtaPhiEnMin(10.),
fSTUTotal(0),              fTRUTotal(0),
fV0Trigger(0),             fV0A(0),              fV0C(0),
fEventMB   (0),            fEventL0   (0),
fEventL1G  (0),            fEventL1G2 (0),
fEventL1J  (0),            fEventL1J2 (0),
fEventCen  (0),            fEventSem  (0),
//Histograms
fhNEvents(0),              fhFORAmp(0),
fhFORAmpL1G(0),            fhFORAmpL1G2(0),
fhFORAmpL1J(0),            fhFORAmpL1J2(0),
fhL0Amp(0),                fhL0AmpL1G(0),        fhL0AmpL1J(0),
fhL1Amp(0),                fhL1GAmp(0),          fhL1G2Amp(0),
fhL1JAmp(0),               fhL1J2Amp(0),         fhL1FOREnergy(0),
fhL0Patch(0),              fhL1GPatch(0),        fhL1G2Patch(0),
fhL1GPatchNotFake(0),      fhL1GPatchFake(0),    fhL1GPatchNotAllFake(0),
fhL1GPatchAllFake(0),      fhL1GPatchNotAllFakeMax(0),
fhL1GPatchAllFakeMax(0),   fhL1GPatchNotAllFakeMaxE(0),
fhL1GPatchAllFakeMaxE(0),  fhL1GPatchNotAllFakeE(0),
fhL1GPatchAllFakeE(0),     fhL1GPatchFakeE(0),
fhL1GPatchNotFakeE(0),     fhNPatchFake(0),      fhNPatchNotFake(0),
fhL1JPatch(0),             fhL1J2Patch(0),
fhFEESTU(0),               fhTRUSTU(0),          fhV0STU(0),
fhGPMaxVV0TT(0),           fhJPMaxVV0TT(0),
fhFORMeanAmp(0),           fhL0MeanAmp(0),       fhL1MeanAmp(0),
fhL1GPatchMax(0),          fhL1G2PatchMax(0),
fhL1JPatchMax(0),          fhL1J2PatchMax(0),
//Histogram settings
fNBinsSTUSignal  (2000),   fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000),   fMaxTRUSignal  (200000),
fNBinsV0Signal   (5000),   fMaxV0Signal   (100000),
fNBinsSTUFEERatio(1000),   fMaxSTUFEERatio(100),
fNBinsSTUTRURatio(1000),   fMaxSTUTRURatio(100),
fNBinsClusterE   (500),    fMaxClusterE   (200),
fMapCell     (),
fMapCellL1G  (),           fMapCellL1G2 (),
fMapCellL1J  (),           fMapCellL1J2 (),
fMapTrigL0   (),           fMapTrigL1   (),
fMapTrigL0L1G(),           fMapTrigL0L1J(),
fMapTrigL1G  (),           fMapTrigL1G2 (),
fMapTrigL1J  (),           fMapTrigL1J2 ()
{
  // Constructor
	
  InitHistogramArrays();
  
  DefineOutput(1, TList::Class());
  
}

//______________________________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA(const char *name) :
AliAnalysisTaskSE(name),
fOutputList(0),            fRecoUtils(0x0),
fGeoSet(0),                fGeometry(0),         fGeoName(""),
fOADBSet(kFALSE),          fAccessOADB(kTRUE),   fOADBFilePath(""),
fBitEGA(0),                fBitEJE(0),
fEtaPhiEnMin(10.),
fSTUTotal(0),              fTRUTotal(0),
fV0Trigger(0),             fV0A(0),              fV0C(0),
fEventMB   (0),            fEventL0   (0),
fEventL1G  (0),            fEventL1G2 (0),
fEventL1J  (0),            fEventL1J2 (0),
fEventCen  (0),            fEventSem  (0),
//Histograms
fhNEvents(0),              fhFORAmp(0),
fhFORAmpL1G(0),            fhFORAmpL1G2(0),
fhFORAmpL1J(0),            fhFORAmpL1J2(0),
fhL0Amp(0),                fhL0AmpL1G(0),        fhL0AmpL1J(0),
fhL1Amp(0),                fhL1GAmp(0),          fhL1G2Amp(0),
fhL1JAmp(0),               fhL1J2Amp(0),         fhL1FOREnergy(0),
fhL0Patch(0),              fhL1GPatch(0),        fhL1G2Patch(0),
fhL1GPatchNotFake(0),      fhL1GPatchFake(0),    fhL1GPatchNotAllFake(0),
fhL1GPatchAllFake(0),      fhL1GPatchNotAllFakeMax(0),
fhL1GPatchAllFakeMax(0),   fhL1GPatchNotAllFakeMaxE(0),
fhL1GPatchAllFakeMaxE(0),  fhL1GPatchNotAllFakeE(0),
fhL1GPatchAllFakeE(0),     fhL1GPatchFakeE(0),
fhL1GPatchNotFakeE(0),     fhNPatchFake(0),      fhNPatchNotFake(0),
fhL1JPatch(0),             fhL1J2Patch(0),
fhFEESTU(0),               fhTRUSTU(0),          fhV0STU(0),
fhGPMaxVV0TT(0),           fhJPMaxVV0TT(0),
fhFORMeanAmp(0),           fhL0MeanAmp(0),       fhL1MeanAmp(0),
fhL1GPatchMax(0),          fhL1G2PatchMax(0),
fhL1JPatchMax(0),          fhL1J2PatchMax(0),
//Histogram settings
fNBinsSTUSignal  (2000),   fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000),   fMaxTRUSignal  (200000),
fNBinsV0Signal   (5000),   fMaxV0Signal   (100000),
fNBinsSTUFEERatio(1000),   fMaxSTUFEERatio(100),
fNBinsSTUTRURatio(1000),   fMaxSTUTRURatio(100),
fNBinsClusterE   (500),    fMaxClusterE   (200),
fMapCell     (),
fMapCellL1G  (),           fMapCellL1G2 (),
fMapCellL1J  (),           fMapCellL1J2 (),
fMapTrigL0   (),           fMapTrigL1   (),
fMapTrigL0L1G(),           fMapTrigL0L1J(),
fMapTrigL1G  (),           fMapTrigL1G2 (),
fMapTrigL1J  (),           fMapTrigL1J2 ()
{
  // Constructor
	
  InitHistogramArrays();
	
  DefineOutput(1, TList::Class());
  
}

//______________________________________________
void AliAnalysisTaskEMCALTriggerQA::AccessOADB()
{
  // Set the AODB  bad channels at least once
  
  //Set it only once
  if(fOADBSet) return ;
  
  if(fOADBFilePath == "") fOADBFilePath = "$ALICE_ROOT/OADB/EMCAL" ;
  
  Int_t   runnumber = InputEvent()->GetRunNumber() ;
  
  if(DebugLevel() > 0)
    printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Get AODB parameters from EMCAL in %s for run %d\n",fOADBFilePath.Data(),runnumber);
  
  Int_t nSM = fGeometry->GetNumberOfSuperModules();
  
  // Bad map
  if(fRecoUtils->IsBadChannelsRemovalSwitchedOn())
  {
    AliOADBContainer *contBC=new AliOADBContainer("");
    contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fOADBFilePath.Data()),"AliEMCALBadChannels");
    
    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runnumber);
    
    if(arrayBC)
    {
      if(DebugLevel() > 0)
        printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Remove EMCAL bad cells \n");
      
      for (Int_t i=0; i<nSM; ++i)
      {
        TH2I *hbm = fRecoUtils->GetEMCALChannelStatusMap(i);
        
        if (hbm)
          delete hbm;
        
        hbm=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
        
        if (!hbm)
        {
          AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
          continue;
        }
        
        hbm->SetDirectory(0);
        fRecoUtils->SetEMCALChannelStatusMap(i,hbm);
        
      } // loop
    } else if(DebugLevel() > 0)
      printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT remove EMCAL bad channels\n"); // run array
  }  // Remove bad
  
  fOADBSet = kTRUE;
  
}

//________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillCellMaps()
{
  // Cells analysis
  // Fill FEE energy per channel array
  Int_t posX    = -1, posY = -1;
  Int_t nSupMod = -1, ieta = -1, iphi = -1, nModule = -1, nIphi = -1, nIeta = -1;
  Short_t absId = -1;
  Int_t nCells  =  0;
  
  AliVCaloCells& cells= *(InputEvent()->GetEMCALCells());
  
  if(cells.IsEMCAL())
  {
    for (Int_t icell = 0; icell <  cells.GetNumberOfCells(); icell++)
    {
      nCells ++;
      
      Double_t amp =0., time = 0., efrac = 0;
      Int_t mclabel = -1;
      
      cells.GetCell(icell, absId, amp, time,mclabel,efrac);
      
      fGeometry->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
      fGeometry->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta);
      
      posX = (nSupMod % 2) ? ieta + AliEMCALGeoParams::fgkEMCALCols : ieta;
      posY = iphi + AliEMCALGeoParams::fgkEMCALRows * int(nSupMod / 2);
      
      Int_t indexX = Int_t(posX/2);
      Int_t indexY = Int_t(posY/2);
      
      if(indexX > fgkFALTROCols || indexY > fgkFALTRORows )
      {
        if(DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Wrong Position (x,y) = (%d,%d)\n",posX,posY);
        continue;
      }
      
      // here it is the amplitude for each cell
      fMapCell[indexY][indexX] += amp;
      
      if(fEventL1G)
      {
        fMapCellL1G[indexY][indexX] += amp;
        //printf("L1G cell[%i,%i] amp=%f\n",indexY,indexX,fMapCellL1G[indexY][indexX]);
      }
      
      if(fEventL1G2)
      {
        fMapCellL1G2[indexY][indexX] += amp;
        //printf("L1G2 cell[%i,%i] amp=%f\n",indexY,indexX,fMapCellL1G2[indexY][indexX]);
      }
			
      if(fEventL1J)  fMapCellL1J [indexY][indexX] += amp;
      if(fEventL1J2) fMapCellL1J2[indexY][indexX] += amp;
			
      //printf("cell[%i,%i] amp=%f\n",indexY,indexX,fMapCell[indexY][indexX]);
			
    }
  }

}

//______________________________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillTriggerPatchMaps(TString triggerclasses)
{
  // Trigger analysis, fill L0, L1 arrays
  
  AliVCaloTrigger& trg= * (InputEvent()->GetCaloTrigger("EMCAL"));
  
  fV0Trigger = trg.GetL1V0(0)+trg.GetL1V0(1); // used elsewhere

  Int_t posX    = -1, posY = -1;

  Int_t    nL0Patch = 0 ;
  //Int_t    nL1Patch = 0 ;
  fSTUTotal = 0;
  fTRUTotal = 0;
  
  trg.Reset();
  // loop on FASTOR
	
  while (trg.Next())
  {
    trg.GetPosition(posX,posY);
    
    if (posX > -1 && posY > -1)
    {
      //L0 analysis
      Int_t nTimes = 0;
      trg.GetNL0Times(nTimes);
      Int_t l0Times[10];
      trg.GetL0Times(l0Times);
			
      Float_t ampL0 = 0.;
      trg.GetAmplitude(ampL0);
      if (ampL0 > 0) fMapTrigL0[posY][posX] = ampL0;
			
      if(triggerclasses.Contains("CEMC7EGA-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EGA") || triggerclasses.Contains("CPBI2EG1")) fMapTrigL0L1G[posY][posX] += ampL0;
      if(triggerclasses.Contains("CEMC7EJE-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EJE") || triggerclasses.Contains("CPBI2EJ1")) fMapTrigL0L1J[posY][posX] += ampL0;
      fTRUTotal += ampL0;
      
      int l0fired = 0;
      for (int itime = 0; itime < nTimes; itime++)
      {
        if (l0Times[itime] > 7 && l0Times[itime] < 10) l0fired = 1;
        // time bin too open? restrict to time bin 8-9?
      }
			
      if (l0fired)
      {
        nL0Patch += nTimes;
        fhL0Patch->Fill(posX,posY);
      }
			
      //L1 analysis
      Int_t bit = 0;
      trg.GetTriggerBits(bit);
      
      Int_t ts = 0;
      trg.GetL1TimeSum(ts);
      if (ts > 0) fMapTrigL1[posY][posX] = ts;
      fSTUTotal += ts;
      // cout << "ts =" <<ts<<endl;
			
      //L1
      Bool_t isEGA1 = ((bit >>  fBitEGA   ) & 0x1) && fEventL1G  ;
      Bool_t isEGA2 = ((bit >> (fBitEGA+1)) & 0x1) && fEventL1G2 ;
      Bool_t isEJE1 = ((bit >>  fBitEJE   ) & 0x1) && fEventL1J  ;
      Bool_t isEJE2 = ((bit >> (fBitEJE+1)) & 0x1) && fEventL1J2 ;
      
      //if(isEGA1 || isEGA2 || isEJE1 || isEJE2) nL1Patch++;
      //if(isEJE1 || isEJE2) printf("Jet STU patch %d, time sum %d, posX %d , posY %d\n",nL1Patch,ts,posX, posY);
      
      //   L1-Gamma
      
      if (isEGA1)
      {
        fhL1GPatch ->Fill(posX,posY);
        if (ts > 0)  fMapTrigL1G [posY][posX] = ts;
      }
      
      if (isEGA2)
      {
        fhL1G2Patch->Fill(posX,posY);
        if (ts > 0) fMapTrigL1G2[posY][posX] = ts;
      }
      
      //    L1-Jet
      if (isEJE1)
      {
        fhL1JPatch ->Fill(posX,posY);
        if (ts > 0) fMapTrigL1J [posY][posX] = ts;
        
      }
      
      if (isEJE2)
      {
        fhL1J2Patch->Fill(posX,posY);
        if (ts > 0) fMapTrigL1J2[posY][posX] = ts;
      }
      
    }
  }
	
  if (!nL0Patch)
  {
    fEventL0 = kFALSE;
    if (!triggerclasses.Contains("CPBI2")) fEventL1G = fEventL1G2 = fEventL1J = fEventL1J2 = kFALSE; // pp running
  }
	
  if(fTRUTotal > fMaxTRUSignal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::FillTriggerPatchMaps() - Large fTRUTotal %f\n",fTRUTotal);
  if(fSTUTotal > fMaxSTUSignal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::FillTriggerPatchMaps() - Large fSTUTotal %d\n",fSTUTotal);
	
}

//___________________________________________________
void AliAnalysisTaskEMCALTriggerQA::ClusterAnalysis()
{
  // Loop on clusters and fill corresponding histograms
  
  //Get Vertex
  Double_t v[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(v);
  
  //clusters distribution
  TRefArray* caloClus = new TRefArray();
  InputEvent()->GetEMCALClusters(caloClus);
  
  Int_t nCaloClusters = caloClus->GetEntriesFast();
  
  Float_t emax   = 0;
  Float_t etamax = 0;
  Float_t phimax = 0;
  Float_t ietamax=-1;
  Float_t iphimax=-1;
  
  Float_t e      = 0;
  Float_t eta    = 0;
  Float_t phi    = 0;
  
  Int_t nSupMod = -1, ieta = -1, iphi = -1;

  TLorentzVector mom;
  
  //Get vertex for momentum calculation
  Double_t vertex[] = {0.0,0.0,0.0};
  //InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  Float_t centrality = -1;
  if(InputEvent()->GetCentrality()) centrality = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
  
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (caloClus->At(icalo));
    
    if(!clus->IsEMCAL()) continue;
    
    if(!fRecoUtils->IsGoodCluster(clus,fGeometry,InputEvent()->GetEMCALCells(),InputEvent()->GetBunchCrossNumber()))
    continue;
    
    if(clus->GetNCells() < 2) continue ; // Avoid 1 cell clusters, noisy, exotic.
    
    clus->GetMomentum(mom, vertex);
    
    Bool_t shared = kFALSE;
    Int_t  idAbs  = -1, iphi0 = -1, ieta0 = -1;
    fRecoUtils->GetMaxEnergyCell(fGeometry, InputEvent()->GetEMCALCells(),clus,
                                 idAbs,nSupMod,ieta0,iphi0,shared);
    //Change index to be continuous over SM
    ieta = (nSupMod % 2) ? ieta0 + AliEMCALGeoParams::fgkEMCALCols : ieta0;
    iphi = iphi0 + AliEMCALGeoParams::fgkEMCALRows * int(nSupMod / 2);
    ieta/=2;
    iphi/=2;
    
    if(ieta > fgkFALTROCols || iphi > fgkFALTRORows )
    printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Wrong Position (x,y) = (%d,%d)\n",ieta,iphi);
    
    e   = clus->E();
    eta = mom.Eta();
    phi = mom.Phi();
    
    if(e > emax)
    {
      emax    = e;
      etamax  = eta;
      phimax  = phi;
      ietamax = ieta;
      iphimax = iphi;
    }
    
    // Fill cluster histograms depending on the event trigger selection
    if( fEventMB  ) FillClusterHistograms(kMBTrig         ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventCen ) FillClusterHistograms(kCentralTrig    ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventSem ) FillClusterHistograms(kSemiCentralTrig,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    
    if( fEventL0  ) FillClusterHistograms(kL0Trig         ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1G ) FillClusterHistograms(kL1GammaTrig    ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1G2) FillClusterHistograms(kL1GammaTrig2   ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1J ) FillClusterHistograms(kL1JetTrig      ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1J2) FillClusterHistograms(kL1JetTrig2     ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);

    if( fEventL1G2 && !fEventL1G) FillClusterHistograms(kL1Gamma2OnlyGammaTrig,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1J2 && !fEventL1J) FillClusterHistograms(kL1Jet2OnlyJetTrig    ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);

    if( fEventL1G && !fEventL1J )
    FillClusterHistograms       (kL1GammaOnlyTrig,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    if( fEventL1J && !fEventL1G )
    FillClusterHistograms       (kL1JetOnlyTrig  ,kFALSE,e,eta,phi,ieta,iphi,centrality,fV0A+fV0C);
    
    if( fEventMB  && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMBPure[0]  ->Fill(e);
    if( fEventCen && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMBPure[1]  ->Fill(e);
    if( fEventSem && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMBPure[2]  ->Fill(e);
    
  }
  
  // Maximum energy cluster per event histograms
  
  if( fEventMB  ) FillClusterHistograms(kMBTrig         ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventCen ) FillClusterHistograms(kCentralTrig    ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventSem ) FillClusterHistograms(kSemiCentralTrig,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  
  if( fEventL0  ) FillClusterHistograms(kL0Trig         ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1G ) FillClusterHistograms(kL1GammaTrig    ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1G2) FillClusterHistograms(kL1GammaTrig2   ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1J ) FillClusterHistograms(kL1JetTrig      ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1J2) FillClusterHistograms(kL1JetTrig2     ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);

  if( fEventL1G2 && !fEventL1G) FillClusterHistograms(kL1Gamma2OnlyGammaTrig,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1J2 && !fEventL1J) FillClusterHistograms(kL1Jet2OnlyJetTrig    ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);

  if( fEventL1G && !fEventL1J )
  FillClusterHistograms         (kL1GammaOnlyTrig,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  if( fEventL1J && !fEventL1G )
  FillClusterHistograms         (kL1JetOnlyTrig  ,kTRUE,emax,etamax,phimax,ietamax,iphimax,centrality,fV0A+fV0C);
  
  if( fEventMB  && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMaxMBPure[0]  ->Fill(emax);
  if( fEventCen && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMaxMBPure[1]  ->Fill(emax);
  if( fEventSem && !fEventL1G && !fEventL1J && !fEventL0  ) fhClusMaxMBPure[2]  ->Fill(emax);

}

//______________________________________________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillClusterHistograms(Int_t triggerNumber, Bool_t max,
                                                          Float_t e,  Float_t eta, Float_t phi,
                                                          Float_t ieta, Float_t iphi,
                                                          Float_t centrality, Float_t fV0AC)
{
  //Fill normal cluster related histograms depending on the trigger type selection
  // (10 options, MB, L0, L1 Gamma ... defined in enum triggerType)
	
  if(!max)
  {
    fhClus   [triggerNumber]->Fill(e);
    fhClusCen[triggerNumber]->Fill(e,centrality);
    fhClusV0 [triggerNumber]->Fill(e,fV0AC);
    fhClusEta[triggerNumber]->Fill(e,eta);
    fhClusPhi[triggerNumber]->Fill(e,phi);
		
    if(e > fEtaPhiEnMin)
    {
      fhClusEtaPhiHigh       [triggerNumber]->Fill( eta, phi);
      fhClusEtaPhiHighCellMax[triggerNumber]->Fill(ieta,iphi);
    }
    else
    {
      fhClusEtaPhiLow       [triggerNumber]->Fill( eta, phi);
      fhClusEtaPhiLowCellMax[triggerNumber]->Fill(ieta,iphi);
    }
  }
  else
  {
    fhClusMax   [triggerNumber]->Fill(e);
    fhClusCenMax[triggerNumber]->Fill(e,centrality);
    fhClusV0Max [triggerNumber]->Fill(e,fV0AC);
    fhClusEtaMax[triggerNumber]->Fill(e,eta);
    fhClusPhiMax[triggerNumber]->Fill(e,phi);
		
    if(e > fEtaPhiEnMin)
    {
      fhClusEtaPhiHighCluMax       [triggerNumber]->Fill( eta, phi);
      fhClusEtaPhiHighCellMaxCluMax[triggerNumber]->Fill(ieta,iphi);
    }
    else
    {
      fhClusEtaPhiLowCluMax       [triggerNumber]->Fill( eta, phi);
      fhClusEtaPhiLowCellMaxCluMax[triggerNumber]->Fill(ieta,iphi);
    }
  }
}

//_____________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillCorrelationHistograms()
{
  //FEE-TRU-STU correlation checks
  
  Double_t ampFOR[30] = {0.}, ampL0[30] = {0.}, ampL1[30] = {0.};
  for (Int_t i = 0; i < fgkFALTRORows; i++)
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++)
    {
      //method to get TRU number
      Int_t idFOR = -1;
      fGeometry->GetAbsFastORIndexFromPositionInEMCAL(j,i,idFOR);
      Int_t iTRU  = -1;
      Int_t iADC  = -1;
      fGeometry->GetTRUFromAbsFastORIndex(idFOR,iTRU,iADC);
      
      //printf("i %d, j %d, iTRU %d, iADC %d, idFOR %d; cell %f, L0 %f, L1 %f\n",
      //       i,j,iTRU,iADC,idFOR, fMapCell  [i][j],fMapTrigL0[i][j],fMapTrigL1[i][j]);
      
      if (iTRU >= 0)
      {
        ampFOR[iTRU] += fMapCell  [i][j];
        ampL0[iTRU]  += fMapTrigL0[i][j];
        ampL1[iTRU]  += fMapTrigL1[i][j];
      }
    }
  }
	
  // FEE vs STU and TRU vs STU ratios
  for (Int_t i = 0; i < 30; i++)
  {
    if (ampFOR[i] != 0 && ampL1[i] != 0)
    {
      fhFEESTU->Fill(ampL1[i]/ampFOR[i],i);
      if(ampL1[i]/ampFOR[i] > fMaxSTUFEERatio  && DebugLevel() > 0 ) printf("AliAnalysisTaskEMCALTriggerQA::FillCorrelationHistograms() - Large STU/FEE ratio %f\n",ampL1[i]/ampFOR[i]);
    }
    
    if (ampL0[i]  != 0 && ampL1[i] != 0)
    {
      fhTRUSTU->Fill(ampL1[i]/ampL0[i] ,i);
      if(ampL1[i]/ampL0[i] > fMaxSTUTRURatio  && DebugLevel() > 0 ) printf("AliAnalysisTaskEMCALTriggerQA::FillCorrelationHistograms() - Large STU/TRU ratio %f\n",ampL1[i]/ampL0[i]);
    }
  }
}

//_____________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillEventCounterHistogram()
{
  // Fill a TH1 histogram, each bin corresponds to a even trigger type
  
  fhNEvents->Fill(0.5); // All physics events
  
  if( fEventMB )
  {
    fhNEvents->Fill(1.5);
    if( !fEventL1G && !fEventL1J && !fEventL1G2 && !fEventL1J2 && !fEventL0 ) fhNEvents->Fill(12.5);
  }
  else
  {
    if( fEventL0  ) fhNEvents->Fill(13.5);
    if( fEventL1G ) fhNEvents->Fill(14.5);
    if( fEventL1J ) fhNEvents->Fill(15.5);
  }
	
  if( fEventCen)   fhNEvents->Fill(2.5);
  if( fEventSem)   fhNEvents->Fill(3.5);
  
  if( fEventL0 )
  {
    fhNEvents->Fill(4.5);
  }
  
  if( fEventL1G )
  {
    fhNEvents->Fill(5.5);
    if(!fEventL1J)  fhNEvents->Fill(9.5);
    if(fEventCen || fEventSem) fhNEvents->Fill(16.5);
  }
	
  if( fEventL1G2 )
  {
    fhNEvents->Fill(6.5);
    if( !fEventL1G ) fhNEvents->Fill(18.5);
    
  }
  
  if( fEventL1J )
  {
    fhNEvents->Fill(7.5);
    if(!fEventL1G)  fhNEvents->Fill(10.5);
    if(fEventCen || fEventSem) fhNEvents->Fill(17.5);
  }
  
  if( fEventL1J2 )
  {
    fhNEvents->Fill(8.5);
    if( !fEventL1J ) fhNEvents->Fill(19.5);
  }
  
  if(fEventL1J && fEventL1G) fhNEvents->Fill(11.5);
  
}

//______________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillL1GammaPatchHistograms()
{
  // L1 Gamma
  
  // Study fakes  - Make it more understandable!!!
  
  Int_t    areAllFakes        = 2;
  Int_t    numberpatchNotFake = 0;
  Int_t    numberpatchFake    = 0;

  Int_t    threshold = 10;// 10 GeV !it's not GeV it's ADC !!
  Bool_t   enoughE   = kFALSE;
  Double_t patchMax  = 0;
  Int_t    colMax    = -1;
  Int_t    rowMax    = -1;
  Int_t    shiftCol  = -1;
  Int_t    shiftRow  = -1;
  
  // loop on patchs
  for (Int_t posx = 0; posx < fgkFALTROCols; posx++)
  {
    for (Int_t posy = 0; posy < fgkFALTRORows; posy++)
    {
      Double_t patchEnergy = 0;
      
      if(fMapTrigL1G[posy][posx] > 0)
      {
        for(Int_t irow = 0; irow < 2; irow++)
        {
          for(Int_t icol = 0; icol < 2; icol++)
          {
            // loop on cells
            shiftCol = posx+icol;
            shiftRow = posy+irow;
            
            //	    printf("cell[%i,%i]=%f\n",posy+icol,posx+irow, fMapCellL1G[posy+icol][posx+irow]);
            if(shiftRow < fgkFALTRORows && shiftCol < fgkFALTROCols)
              patchEnergy += fMapCellL1G[shiftRow][shiftCol] ;
            
            if( fMapCellL1G[shiftRow][shiftCol] >threshold/2) enoughE = kTRUE;
          }
        }
        
        if (patchEnergy > patchMax)
        {
          patchMax = patchEnergy;
          colMax = posx;
          rowMax = posy;
        }
        
        if(patchEnergy>threshold || (patchEnergy>threshold-3 && enoughE))
        {
          numberpatchNotFake++;
          fhL1GPatchNotFake ->Fill(posx,posy);
          fhL1GPatchNotFakeE->Fill(patchEnergy);
          areAllFakes = 1;
        }
        else
        {
          numberpatchFake++;
          areAllFakes = 0;
          fhL1GPatchFake ->Fill(posx,posy);
          fhL1GPatchFakeE->Fill(patchEnergy);
        }
      }
    }
  }
  
  fhNPatchNotFake->Fill(areAllFakes,numberpatchNotFake);
  fhNPatchFake   ->Fill(areAllFakes,numberpatchFake);
  
  if(areAllFakes == 0)
  {
    // loop on patchs
    for (Int_t col = 0; col < fgkFALTROCols; col++)
    {
      for (Int_t row = 0; row < fgkFALTRORows; row++)
      {
        if(fMapTrigL1G[row][col] > 0)
        {
          //	cout <<"checking fMapTrigL1G[row][col]"<<fMapTrigL1G[row][col]<<endl;
          fhL1GPatchAllFake->Fill(col,row);
          
          Double_t patchEnergy=0;
          for(Int_t irow = 0; irow < 2; irow++)
          {
            for(Int_t icol = 0; icol < 2; icol++)
            {
              shiftCol = col+icol;
              shiftRow = row+irow;
              
              if(shiftRow < fgkFALTRORows && shiftCol < fgkFALTROCols)
                patchEnergy += fMapCellL1G[shiftRow][shiftCol] ;

            }
          }
          fhL1GPatchAllFakeE->Fill(patchEnergy);
        }
      }
    }
    //  cout << "row max"<<rowMax<<" colmax"<<colMax<< " fMapTrigL1G[rowMax][colMax]"<< fMapTrigL1G[rowMax][colMax]<<endl;
    
    if(fMapTrigL1G[rowMax][colMax] > 0)
    {
      //  printf("\npatch max [%i,%i] = %f\n",rowMax,colMax,patchMax);
      fhL1GPatchAllFakeMax ->Fill(colMax,rowMax);
      fhL1GPatchAllFakeMaxE->Fill(patchMax);
    }
  }
  else
  {
    // loop on patches
    for (Int_t col = 0; col < fgkFALTROCols; col++)
    {
      for (Int_t row = 0; row < fgkFALTRORows; row++)
      {
        if(fMapTrigL1G[row][col] > 0)
        {
          fhL1GPatchNotAllFake->Fill(col,row);
          
          Double_t patchEnergy = 0;
          for(Int_t irow = 0; irow < 2; irow++)
          {
            for(Int_t icol = 0; icol < 2; icol++)
            {
              shiftCol = col+icol;
              shiftRow = row+irow;
              
              if(shiftRow < fgkFALTRORows && shiftCol < fgkFALTROCols)
                patchEnergy += fMapCellL1G[shiftRow][shiftCol] ;
            }
          }
          
          fhL1GPatchNotAllFakeE->Fill(patchEnergy);
          
        }
      }
    }
    
    if(fMapTrigL1G[rowMax][colMax] > 0 )
    {
      fhL1GPatchNotAllFakeMax ->Fill(colMax,rowMax);
      fhL1GPatchNotAllFakeMaxE->Fill(patchMax);
    }
  }
  
  fhGPMaxVV0TT->Fill(fV0Trigger, patchMax);
  if( fEventL1G )  fhL1GPatchMax ->Fill(colMax,rowMax);
  if( fEventL1G2 ) fhL1G2PatchMax->Fill(colMax,rowMax);
  
}

//____________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillL1JetPatchHistograms()
{
  // L1 Jet
  
  Double_t patchMax =  0;
  Int_t    colMax   = -1;
  Int_t    rowMax   = -1;
  Int_t    col, row =  0;
  
  for (Int_t i = 0; i < 9; i++)
  {
    for (Int_t j = 0; j < 12; j++)
    {
      Int_t patchJ = 0;
      col = i;
      row = j;
      
      for (Int_t k = 0; k < 16; k++)
      {
        for (Int_t l = 0; l < 16; l++)
        {
          patchJ += Int_t(fMapTrigL1[4*j + l][4*i + k]);
        }
      }
      
      if (patchJ > patchMax)
      {
        patchMax = patchJ;
        colMax = 4*col;
        rowMax = 4*row;
      }
    }
  }
  
  fhJPMaxVV0TT->Fill(fV0Trigger, patchMax);
  if( fEventL1J )  fhL1JPatchMax ->Fill(colMax,rowMax);
  if( fEventL1J2 ) fhL1J2PatchMax->Fill(colMax,rowMax);
  
}

//______________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillMapHistograms()
{
  //Matrix with signal per channel
  
  for (Int_t i = 0; i < fgkFALTRORows; i++)
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++) //check x,y direction for reading FOR ((0,0) = top left);
    {
      fhFORAmp    ->Fill( j, i, fMapCell     [i][j]);
      fhFORAmpL1G ->Fill( j, i, fMapCellL1G  [i][j]);
      fhFORAmpL1G2->Fill( j, i, fMapCellL1G2 [i][j]);
      fhFORAmpL1J ->Fill( j, i, fMapCellL1J  [i][j]);
      fhFORAmpL1J2->Fill( j, i, fMapCellL1J2 [i][j]);
      fhL0Amp     ->Fill( j, i, fMapTrigL0   [i][j]);
      fhL0AmpL1G  ->Fill( j, i, fMapTrigL0L1G[i][j]);
      fhL0AmpL1J  ->Fill( j, i, fMapTrigL0L1J[i][j]);
      fhL1Amp     ->Fill( j, i, fMapTrigL1   [i][j]);
      
      fhL1FOREnergy->Fill(i+fgkFALTRORows*j,  fMapTrigL1   [i][j]);
      fhL1GAmp    ->Fill( j, i, fMapTrigL1G  [i][j]);
      fhL1G2Amp   ->Fill( j, i, fMapTrigL1G2 [i][j]);
      fhL1JAmp    ->Fill( j, i, fMapTrigL1J  [i][j]);
      fhL1J2Amp   ->Fill( j, i, fMapTrigL1J2 [i][j]);
      fhFORMeanAmp->Fill( j, i, fMapCell     [i][j]);
      fhL0MeanAmp ->Fill( j, i, fMapTrigL0   [i][j]);
      fhL1MeanAmp ->Fill( j, i, fMapTrigL1   [i][j]);
    }
  }
}

//____________________________________________________
void AliAnalysisTaskEMCALTriggerQA::FillV0Histograms()
{
  //V0 analysis, only for ESDs
  
  AliESDVZERO* eventV0 = dynamic_cast<AliESDVZERO*> (InputEvent()->GetVZEROData());
  
  if(eventV0)
  {
    for (Int_t i = 0; i < 32; i++)
    {
      fV0C += eventV0->GetAdcV0C(i);
      fV0A += eventV0->GetAdcV0A(i);
    }
    
    if (fSTUTotal != 0)
    {
      fhV0STU->Fill(fV0A+fV0C,fSTUTotal);
      if( fV0A+fV0C > fMaxV0Signal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large fV0A+fV0C %f\n",fV0A+fV0C);
    }
    
    if( fEventL1G )  fhV0[kL1GammaTrig]    ->Fill(fV0A+fV0C);
    if( fEventL1G2 ) fhV0[kL1GammaTrig2]   ->Fill(fV0A+fV0C);
    if( fEventL1J )  fhV0[kL1JetTrig]      ->Fill(fV0A+fV0C);
    if( fEventL1J2 ) fhV0[kL1JetTrig2]     ->Fill(fV0A+fV0C);
    if( fEventMB  )  fhV0[kMBTrig]         ->Fill(fV0A+fV0C);
    if( fEventL0  )  fhV0[kL0Trig]         ->Fill(fV0A+fV0C);
    if( fEventCen )  fhV0[kCentralTrig]    ->Fill(fV0A+fV0C);
    if( fEventSem )  fhV0[kSemiCentralTrig]->Fill(fV0A+fV0C);
    if( fEventL1G  && !fEventL1J) fhV0[kL1GammaOnlyTrig]      ->Fill(fV0A+fV0C);
    if( fEventL1J  && !fEventL1G) fhV0[kL1JetOnlyTrig]        ->Fill(fV0A+fV0C);
    if( fEventL1G2 && !fEventL1G) fhV0[kL1Gamma2OnlyGammaTrig]->Fill(fV0A+fV0C);
    if( fEventL1J2 && !fEventL1J) fhV0[kL1Jet2OnlyJetTrig]    ->Fill(fV0A+fV0C);
    //if(nL0Patch!=0 || nL1Patch!=0) printf("total TRU %f, total STU %f, V0C+V0A %f; nL0 %d, nL1 %d \n",
    //       fTRUTotal,fSTUTotal,fV0A+fV0C,nL0Patch,nL1Patch);
  }
}

//________________________________________
void AliAnalysisTaskEMCALTriggerQA::Init()
{
  //Init analysis parameters not set before
	
  if(!fRecoUtils)
  {
    fRecoUtils    = new AliEMCALRecoUtils ;
    fRecoUtils->SwitchOnBadChannelsRemoval();
  }
  
}

//_____________________________________________________
void AliAnalysisTaskEMCALTriggerQA::InitCellPatchMaps()
{
  //Init to 0 and per event different cell/patch maps, depending on trigger type

  for (Int_t i = 0; i < fgkFALTRORows; i++)
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++)
    {
      fMapTrigL0   [i][j] = 0.;
      fMapTrigL0L1G[i][j] = 0.;
      fMapTrigL0L1J[i][j] = 0.;
      fMapTrigL1G  [i][j] = 0.;
      fMapTrigL1G2 [i][j] = 0.;
      fMapTrigL1J  [i][j] = 0.;
      fMapTrigL1J2 [i][j] = 0.;
      fMapTrigL1   [i][j] = 0.;
      fMapCell     [i][j] = 0.;
      fMapCellL1G  [i][j] = 0.;
      fMapCellL1G2 [i][j] = 0.;
      fMapCellL1J  [i][j] = 0.;
      fMapCellL1J2 [i][j] = 0.;
    }
  }
}

//________________________________________________
void AliAnalysisTaskEMCALTriggerQA::InitGeometry()
{
  // Init geometry and set the geometry matrix, for the first event, skip the rest
  // Also set once the run dependent calibrations
  
  if(fGeoSet) return;
  
  // Init the trigger bit once, correct depending on version
  fBitEGA = 4;
  fBitEJE = 5;
  
  TFile* file = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
	
  const TList *clist = file->GetStreamerInfoCache();
  
  if(clist)
  {
    TStreamerInfo *cinfo = (TStreamerInfo*)clist->FindObject("AliESDCaloTrigger");
    Int_t verid = 5; // newer ESD header version
    if(!cinfo)
    {
      cinfo = (TStreamerInfo*)clist->FindObject("AliAODCaloTrigger");
      verid = 2; // newer AOD header version
    }
    if(cinfo)
    {
      Int_t classversionid = cinfo->GetClassVersion();
      //printf("********* Header class version %d *********** \n",classversionid);
      
      if (classversionid >= verid)
      {
        fBitEGA = 6;
        fBitEJE = 8;
      }
    }  else printf("AliAnalysisTaskEMCALTriggerQA - Streamer info for trigger class not available, bit not changed\n");
  } else printf("AliAnalysisTaskEMCALTriggerQA - Streamer list not available!, bit not changed\n");
  
  Int_t runnumber = InputEvent()->GetRunNumber() ;
  
  if (!fGeometry)
  {
    if(fGeoName=="")
    {
      if     (runnumber < 140000) fGeoName = "EMCAL_FIRSTYEARV1";
      else if(runnumber < 171000) fGeoName = "EMCAL_COMPLETEV1";
      else                        fGeoName = "EMCAL_COMPLETE12SMV1";
      if(DebugLevel() > 0)
        printf("AliAnalysisTaskEMCALTriggerQA::InitGeometry() - Set EMCAL geometry name to <%s> for run %d\n",fGeoName.Data(),runnumber);
    }
    
    fGeometry = AliEMCALGeometry::GetInstance(fGeoName);
  }
  
  fGeoSet = kTRUE;
	
}

//_______________________________________________________
void AliAnalysisTaskEMCALTriggerQA::InitHistogramArrays()
{
  //Histograms array initialization
  
  for (Int_t i = 0; i < fgkTriggerCombi; i++)
  {
    fhV0     [i] = 0;
    fhClus   [i] = 0;              fhClusMax   [i] = 0;
    fhClusCen[i] = 0;              fhClusCenMax[i] = 0;
    fhClusV0 [i] = 0;              fhClusV0Max [i] = 0;
    fhClusEta[i] = 0;              fhClusEtaMax[i] = 0;
    fhClusPhi[i] = 0;              fhClusPhiMax[i] = 0;
    
    fhClusEtaPhiHigh       [i] = 0;  fhClusEtaPhiHighCluMax       [i] = 0;
    fhClusEtaPhiHighCellMax[i] = 0;  fhClusEtaPhiHighCellMaxCluMax[i] = 0;
    fhClusEtaPhiLow        [i] = 0;  fhClusEtaPhiLowCluMax        [i] = 0;
    fhClusEtaPhiLowCellMax [i] = 0;  fhClusEtaPhiLowCellMaxCluMax [i] = 0;
    
    if(i<3){ fhClusMBPure[i] = 0;   fhClusMaxMBPure[i] = 0; }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::SetTriggerEventBit( TString triggerclasses)
{
  // Check what trigger is the event, set the corresponding bit
  
  // Init trigger event bit
  fEventMB   = kFALSE;
  fEventL0   = kFALSE;
  fEventL1G  = kFALSE;
  fEventL1G2 = kFALSE;
  fEventL1J  = kFALSE;
  fEventL1J2 = kFALSE;
  fEventCen  = kFALSE;
  fEventSem  = kFALSE;
  
  //Min bias event trigger?
  if((triggerclasses.Contains("CINT") || triggerclasses.Contains("CPBI2_B1") ) &&
     (triggerclasses.Contains("-B-")  || triggerclasses.Contains("-I-"))       &&
		 triggerclasses.Contains("-NOPF-ALLNOTRD") )   fEventMB  = kTRUE;
  
  // EMC triggered event? Which type?
  if( triggerclasses.Contains("-B-") || triggerclasses.Contains("-S-") || triggerclasses.Contains("-I-") )
  {
    if( triggerclasses.Contains("CEMC") &&
       !triggerclasses.Contains("EGA" ) &&
       !triggerclasses.Contains("EJE" )	&&
       !triggerclasses.Contains("EG1" ) &&
       !triggerclasses.Contains("EJ1" )	&&
       !triggerclasses.Contains("EG2" ) &&
       !triggerclasses.Contains("EJ2" )    ) fEventL0  = kTRUE;
    
    if( triggerclasses.Contains("EGA" ) || triggerclasses.Contains("EG1" )   ) fEventL1G  = kTRUE;
    if( triggerclasses.Contains("EG2" ) )                                      fEventL1G2 = kTRUE;
    
    if( triggerclasses.Contains("EJE" ) || triggerclasses.Contains("EJ1" )   ) fEventL1J  = kTRUE;
    if( triggerclasses.Contains("EJ2" ) )                                      fEventL1J2 = kTRUE;
  }
  
  // Semi/Central PbPb trigger
  if     (triggerclasses.Contains("CCENT_R2-B-NOPF-ALLNOTRD")) fEventCen = kTRUE;
  else if(triggerclasses.Contains("CSEMI_R1-B-NOPF-ALLNOTRD")) fEventSem = kTRUE;

  //  printf("MB : %d; L0 : %d; L1-Gam : %d; L1-Jet : %d; Central : %d; SemiCentral : %d; Trigger Names : %s \n ",
	//       fEventMB,fEventL0,fEventL1G,fEventL1J,fEventCen,fEventSem,triggerclasses.Data());
  
}

//___________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserCreateOutputObjects()
{
  // Init histograms and geometry
  
  fOutputList  = new TList;
  fOutputList ->SetOwner(kTRUE);
  
  fhNEvents    = new TH1F("hNEvents","Number of selected events",20,0,20);
  fhNEvents   ->SetYTitle("N events");
  fhNEvents   ->GetXaxis()->SetBinLabel(1 ,"All");
  fhNEvents   ->GetXaxis()->SetBinLabel(2 ,"MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(3 ,"Central Pb");
  fhNEvents   ->GetXaxis()->SetBinLabel(4 ,"SemiCentral Pb");
  fhNEvents   ->GetXaxis()->SetBinLabel(5 ,"L0");
  fhNEvents   ->GetXaxis()->SetBinLabel(6 ,"L1-G1");
  fhNEvents   ->GetXaxis()->SetBinLabel(7 ,"L1-G2");
  fhNEvents   ->GetXaxis()->SetBinLabel(8 ,"L1-J1");
  fhNEvents   ->GetXaxis()->SetBinLabel(9 ,"L1-J2");
  fhNEvents   ->GetXaxis()->SetBinLabel(10 ,"L1-G1 & !L1-J1");
  fhNEvents   ->GetXaxis()->SetBinLabel(11 ,"L1-J1 & !L1-G1");
  fhNEvents   ->GetXaxis()->SetBinLabel(12 ,"L1-J1 & L1-G1");
  fhNEvents   ->GetXaxis()->SetBinLabel(13 ,"MB & !L1 & !L0");
  fhNEvents   ->GetXaxis()->SetBinLabel(14,"L0 & !MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(15,"L1-G1 & !MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(16,"L1-J1 & !MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(17,"L1-G1 & (Cen | Semi)");
  fhNEvents   ->GetXaxis()->SetBinLabel(18,"L1-J1 & (Cen | Semi)");
  fhNEvents   ->GetXaxis()->SetBinLabel(19,"L1-G2 & !L1-G1");
  fhNEvents   ->GetXaxis()->SetBinLabel(20,"L1-J2 & !L1-J1");

  fhFORAmp     = new TH2F("hFORAmp", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmp    ->SetXTitle("Index #eta (columnns)");
  fhFORAmp    ->SetYTitle("Index #phi (rows)");
  fhFORAmp    ->SetZTitle("Amplitude");
  
  fhFORAmpL1G  = new TH2F("hFORAmpL1G", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1G trigger condition",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1G ->SetXTitle("Index #eta (columnns)");
  fhFORAmpL1G ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1G ->SetZTitle("Amplitude");
	
  fhFORAmpL1G2  = new TH2F("hFORAmpL1G2", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1G2 trigger condition",
                           fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1G2 ->SetXTitle("Index #eta (columnns)");
  fhFORAmpL1G2 ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1G2 ->SetZTitle("Amplitude");
  
  fhFORAmpL1J  = new TH2F("hFORAmpL1J", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1J trigger condition",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1J ->SetXTitle("Index #eta (columnns)");
  fhFORAmpL1J ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1J ->SetZTitle("Amplitude");
	
  fhFORAmpL1J2  = new TH2F("hFORAmpL1J2", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1J2 trigger condition",
                           fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1J2 ->SetXTitle("Index #eta (columnns)");
  fhFORAmpL1J2 ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1J2 ->SetZTitle("Amplitude");
  
  
  fhL0Amp      = new TH2F("hL0Amp","FALTRO signal per Row and Column",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Amp     ->SetXTitle("Index #eta (columnns)");
  fhL0Amp     ->SetYTitle("Index #phi (rows)");
  fhL0Amp     ->SetZTitle("Amplitude");
  
  fhL0AmpL1G   = new TH2F("hL0AmpL1G","FALTRO signal per Row and Column, with L1G trigger condition",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0AmpL1G  ->SetXTitle("Index #eta (columnns)");
  fhL0AmpL1G  ->SetYTitle("Index #phi (rows)");
  fhL0AmpL1G  ->SetZTitle("Amplitude");
	
  
  fhL0AmpL1J   = new TH2F("hL0AmpL1J","FALTRO signal per Row and Column, with L1j trigger condition",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0AmpL1J  ->SetXTitle("Index #eta (columnns)");
  fhL0AmpL1J  ->SetYTitle("Index #phi (rows)");
  fhL0AmpL1J  ->SetZTitle("Amplitude");
	
  
  fhL1Amp      = new TH2F("hL1Amp","STU signal per Row and Column",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1Amp     ->SetXTitle("Index #eta (columnns)");
  fhL1Amp     ->SetYTitle("Index #phi (rows)");
  fhL1Amp     ->SetZTitle("Amplitude");
  
  fhL1GAmp     = new TH2F("hL1GAmp","STU signal per Row and Column for L1 Gamma",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GAmp    ->SetXTitle("Index #eta (columnns)");
  fhL1GAmp    ->SetYTitle("Index #phi (rows)");
  fhL1GAmp    ->SetZTitle("Amplitude");
	
  fhL1G2Amp     = new TH2F("hL1G2Amp","STU signal per Row and Column for L1 Gamma2",
                           fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1G2Amp    ->SetXTitle("Index #eta (columnns)");
  fhL1G2Amp    ->SetYTitle("Index #phi (rows)");
  fhL1G2Amp    ->SetZTitle("Amplitude");
  
  fhL1JAmp     = new TH2F("hL1JAmp","STU signal per Row and Column for L1 Jet",
                          fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JAmp    ->SetXTitle("Index #eta (columnns)");
  fhL1JAmp    ->SetYTitle("Index #phi (rows)");
  fhL1JAmp    ->SetZTitle("Amplitude");
	
  fhL1J2Amp     = new TH2F("hL1J2Amp","STU signal per Row and Column for L1 Jet2",
                           fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1J2Amp    ->SetXTitle("Index #eta (columnns)");
  fhL1J2Amp    ->SetYTitle("Index #phi (rows)");
  fhL1J2Amp    ->SetZTitle("Amplitude");
	
  fhL1FOREnergy     = new TH2F("hL1FOREnergy","FOR index vs FOR energy",
                               fgkFALTROCols*fgkFALTRORows,0,fgkFALTROCols*fgkFALTRORows,200,0,200);
  fhL1FOREnergy    ->SetXTitle("Index FOR");
  fhL1FOREnergy    ->SetYTitle("Energy (ADC)");
  
  fhL0Patch    = new TH2F("hL0Patch","FOR with associated L0 Patch",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Patch   ->SetXTitle("Index #eta (columnns)");
  fhL0Patch   ->SetYTitle("Index #phi (rows)");
  fhL0Patch   ->SetZTitle("counts");
  
  fhL1GPatch   = new TH2F("hL1GPatch","FOR with associated L1 Gamma Patch",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatch  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatch  ->SetYTitle("Index #phi (rows)");
  fhL1GPatch  ->SetZTitle("counts");
	
  fhL1G2Patch   = new TH2F("hL1G2Patch","FOR with associated L1 Gamma2 Patch",
                           fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1G2Patch  ->SetXTitle("Index #eta (columnns)");
  fhL1G2Patch  ->SetYTitle("Index #phi (rows)");
  fhL1G2Patch  ->SetZTitle("counts");
  
  fhL1GPatchNotFake   = new TH2F("hL1GPatchNotFake","FOR with L1 Gamma Patch associated to energetic cells",
                                 fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchNotFake  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchNotFake  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchNotFake  ->SetZTitle("counts");
	
  fhL1GPatchFake   = new TH2F("hL1GPatchFake","FOR without L1 Gamma Patch associated to energetic cells",
                              fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchFake  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchFake  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchFake  ->SetZTitle("counts");
	
	
  fhL1GPatchNotAllFake   = new TH2F("hL1GPatchNotAllFake","FOR with one L1 Gamma Patch associated to an energetic cell",
                                    fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchNotAllFake  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchNotAllFake  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchNotAllFake  ->SetZTitle("counts");
	
  fhL1GPatchAllFake   = new TH2F("hL1GPatchAllFake","FOR without any L1 Gamma Patch associated to an energetic cell",
                                 fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchAllFake  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchAllFake  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchAllFake  ->SetZTitle("counts");
	
  fhL1GPatchAllFakeMax   = new TH2F("hL1GPatchAllFakeMax","FOR with L1 Gamma Patch Max not associated to an energetic cell",
                                    fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchAllFakeMax  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchAllFakeMax  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchAllFakeMax  ->SetZTitle("counts");
	
  fhL1GPatchNotAllFakeMax   = new TH2F("hL1GPatchNotAllFakeMax","FOR with one L1 Gamma Patch Max associated to an energetic cell",
                                       fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchNotAllFakeMax  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchNotAllFakeMax  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchNotAllFakeMax  ->SetZTitle("counts");
	
	fhL1GPatchNotAllFakeMaxE   = new TH1F("hL1GPatchNotAllFakeMaxE","Energy distribution of FOR in events with L1 Gamma Patch Max associated to an energetic cell",
                                        fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchNotAllFakeMaxE ->SetXTitle("Energy (GeV)");
  
	
  fhL1GPatchAllFakeMaxE   = new TH1F("hL1GPatchAllFakeMaxE","Energy distribution of FOR in events with L1 Gamma Patch Max not associated to an energetic cell",
                                     fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchAllFakeMaxE ->SetXTitle("Energy (GeV)");
  
  fhL1GPatchNotAllFakeE   = new TH1F("hL1GPatchNotAllFakeE","Energy distribution of FOR in events with L1 Gamma Patch not associated to an energetic cell",
                                     fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchNotAllFakeE ->SetXTitle("Energy (GeV)");
	
  fhL1GPatchAllFakeE   = new TH1F("hL1GPatchAllFakeE","Energy distribution of FOR in events with L1 Gamma Patch  associated to an energetic cell",
                                  fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchAllFakeE ->SetXTitle("Energy (GeV)");
	
	
  fhL1GPatchFakeE   = new TH1F("hL1GPatchFakeE","Energy distribution of FOR with L1 Gamma Patch not associated to an energetic cell",
                               fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchFakeE ->SetXTitle("Energy (GeV)");
	
  fhL1GPatchNotFakeE   = new TH1F("hL1GPatchNotFakeE","Energy distribution of FOR with L1 Gamma Patch  associated to an energetic cell",
                                  fNBinsClusterE,0,fMaxClusterE);
	fhL1GPatchNotFakeE ->SetXTitle("Energy (GeV)");
	
  fhNPatchFake   = new TH2F("hNPatchFake","number of fake patchs vs. all patchs are fake",
                            3,0,3, 2880,0,2880);
  fhNPatchFake  ->SetYTitle("number of fake patchs");
  fhNPatchFake  ->SetXTitle("all fake event");
  fhNPatchFake  ->SetZTitle("counts");
	
	
  fhNPatchNotFake   = new TH2F("hNPatchNotFake","number of Not fake patchs vs. all patchs are fake",
                               3, 0, 3, 2000,0,2000);
  fhNPatchNotFake  ->SetYTitle("number of Not fake patchs");
  fhNPatchNotFake  ->SetXTitle("all fake event");
  fhNPatchNotFake  ->SetZTitle("counts");
	
	
  fhL1JPatch   = new TH2F("hL1JPatch","FOR with associated L1 Jet Patch",
                          fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JPatch  ->SetXTitle("Index #eta (columnns)");
  fhL1JPatch  ->SetYTitle("Index #phi (rows)");
  fhL1JPatch  ->SetZTitle("counts");
	
  fhL1J2Patch   = new TH2F("hL1J2Patch","FOR with associated L1 Jet2 Patch",
                           fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1J2Patch  ->SetXTitle("Index #eta (columnns)");
  fhL1J2Patch  ->SetYTitle("Index #phi (rows)");
  fhL1J2Patch  ->SetZTitle("counts");
  
  fhV0STU      = new TH2I("hV0STU","Total signal STU vs V0C+V0S",
                          fNBinsV0Signal,0,fMaxV0Signal,fNBinsSTUSignal,0,fMaxSTUSignal);
  fhV0STU     ->SetXTitle("Signal V0C+V0A");
  fhV0STU     ->SetYTitle("Total signal STU");
  fhV0STU     ->SetZTitle("counts");
  
  
  fhFEESTU     = new TH2F("hFEESTU","STU / FEE vs channel", fNBinsSTUFEERatio,0,fMaxSTUFEERatio,30,0,30);
  fhFEESTU    ->SetXTitle("STU/FEE signal");
  fhFEESTU    ->SetYTitle("channel");
  fhFEESTU    ->SetZTitle("counts");
  
  fhTRUSTU     = new TH2F("hTRUSTU","STU / TRU vs channel", fNBinsSTUTRURatio,0,fMaxSTUTRURatio,30,0,30);
  fhTRUSTU    ->SetXTitle("STU/TRU signal");
  fhTRUSTU    ->SetYTitle("channel");
  fhTRUSTU    ->SetZTitle("counts");
  
  fhGPMaxVV0TT = new TH2F("hGPMaxVV0TT","Maximum patch of L1-Gamma vs V0 signal in STU",fNBinsV0Signal,0,fMaxV0Signal,1000,0,1000);
  fhGPMaxVV0TT ->SetXTitle("V0 from STU");
  fhGPMaxVV0TT ->SetYTitle("Patch Max");
  
  fhJPMaxVV0TT = new TH2F("hJPMaxVV0TT","Maximum patch of L1-Jet   vs V0 signal in STU",fNBinsV0Signal,0,fMaxV0Signal,1000,0,1000);
  fhJPMaxVV0TT ->SetXTitle("V0 from STU");
  fhJPMaxVV0TT ->SetYTitle("Patch Max");
  
  fhFORMeanAmp = new TProfile2D("hFORMeanAmp", "Mean FastOR(FEE) signal per Row and Column", fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORMeanAmp->SetXTitle("Index #eta");
  fhFORMeanAmp->SetYTitle("Index #phi");
  
  fhL0MeanAmp = new TProfile2D("hL0MeanAmp", "Mean FastOR(TRU) signal per Row and Column", fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0MeanAmp->SetXTitle("Index #eta");
  fhL0MeanAmp->SetYTitle("Index #phi");
  
  fhL1MeanAmp = new TProfile2D("hL1MeanAmp", "Mean FastOR(STU) signal per Row and Column", fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1MeanAmp->SetXTitle("Index #eta");
  fhL1MeanAmp->SetYTitle("Index #phi");
  
  fhL1GPatchMax   = new TH2F("hL1GPatchMax","FOR of max amplitude patch with associated L1 Gamma Patch",
                             fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchMax  ->SetZTitle("counts");
	
  fhL1G2PatchMax   = new TH2F("hL1G2PatchMax","FOR of max amplitude patch with associated L1 Gamma2 Patch",
                              fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1G2PatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1G2PatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1G2PatchMax  ->SetZTitle("counts");
  
  fhL1JPatchMax   = new TH2F("hL1JPatchMax","FOR of max amplitude patch with associated L1 Jet Patch",
                             fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JPatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1JPatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1JPatchMax  ->SetZTitle("counts");
	
  fhL1J2PatchMax   = new TH2F("hL1JPatchMax","FOR of max amplitude patch with associated L1 Jet2 Patch",
                              fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1J2PatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1J2PatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1J2PatchMax  ->SetZTitle("counts");
  
  fOutputList->Add(fhNEvents);
  fOutputList->Add(fhV0STU);
  fOutputList->Add(fhFORAmp);
  fOutputList->Add(fhFORAmpL1G);
  fOutputList->Add(fhFORAmpL1G2);
  fOutputList->Add(fhFORAmpL1J);
  fOutputList->Add(fhFORAmpL1J2);
  fOutputList->Add(fhL0Amp);
  fOutputList->Add(fhL0AmpL1G);
  fOutputList->Add(fhL0AmpL1J);
  fOutputList->Add(fhL1Amp);
  fOutputList->Add(fhL1GAmp);
  fOutputList->Add(fhL1G2Amp);
  fOutputList->Add(fhL1JAmp);
  fOutputList->Add(fhL1J2Amp);
  fOutputList->Add(fhL1FOREnergy);
  fOutputList->Add(fhL0Patch);
  fOutputList->Add(fhL1GPatch);
  fOutputList->Add(fhL1G2Patch);
  fOutputList->Add(fhL1GPatchNotFake);
  fOutputList->Add(fhL1GPatchFake);
  fOutputList->Add(fhL1GPatchAllFake);
  fOutputList->Add(fhL1GPatchNotAllFake);
  fOutputList->Add(fhL1GPatchNotAllFakeMax);
  fOutputList->Add(fhL1GPatchAllFakeMax);
  fOutputList->Add(fhL1GPatchNotAllFakeMaxE);
  fOutputList->Add(fhL1GPatchAllFakeMaxE);
  fOutputList->Add(fhL1GPatchNotAllFakeE);
  fOutputList->Add(fhL1GPatchAllFakeE);
  fOutputList->Add(fhL1GPatchFakeE);
  fOutputList->Add(fhL1GPatchNotFakeE);
  fOutputList->Add(fhNPatchFake);
  fOutputList->Add(fhNPatchNotFake);
	
  fOutputList->Add(fhL1JPatch);
  fOutputList->Add(fhL1J2Patch);
  fOutputList->Add(fhFEESTU);
  fOutputList->Add(fhTRUSTU);
  
  fOutputList->Add(fhGPMaxVV0TT);
  fOutputList->Add(fhJPMaxVV0TT);
  
  fOutputList->Add(fhFORMeanAmp);
  fOutputList->Add(fhL0MeanAmp);
  fOutputList->Add(fhL1MeanAmp);
	
  fOutputList->Add(fhL1GPatchMax);
  fOutputList->Add(fhL1G2PatchMax);
  fOutputList->Add(fhL1JPatchMax);
  fOutputList->Add(fhL1J2PatchMax);
  
  // Cluster histograms, E
  TString hName  [] = {"MB","L0","L1G1","L1G2","L1J1","L1J2","L1G1NoL1J1","L1J1NoLG1","L1G2NoL1G1","L1J2NoL1J1","Central","SemiCentral"};
  TString hTitle [] = {"MB trigger","L0 trigger","L1 Gamma1 trigger","L1 Gamma2 trigger","L1 Jet1 trigger","L1 Jet2 trigger",
                       "L1 Gamma1 trigger and not L1 Jet1"  ,"L1 Jet1 trigger and not L1 Gamma1",
                       "L1 Gamma2 trigger and not L1 Gamma1","L1 Jet2 trigger and not L1 Jet1",
                        "Central trigger","SemiCentral trigger"};
	
  for(Int_t i=0; i < 3; i++)
  {
    Int_t j = i+5;
    if(i==0)j=0;
    
    fhClusMBPure[i]  = new TH1F(Form("hClus%sPure",hName[j].Data()),
                                Form("clusters E distribution for %s, no other EMCAL trigger on",hTitle[j].Data()),
                                fNBinsClusterE,0,fMaxClusterE);
    fhClusMBPure[i] ->SetXTitle("Energy (GeV)");
    fOutputList->Add(fhClusMBPure[i]);
		
    fhClusMaxMBPure[i]  = new TH1F(Form("hClusMax%sPure",hName[j].Data()),
                                   Form("maximum energy cluster per event for %s, no other EMCAL trigger on",hTitle[j].Data()),
                                   fNBinsClusterE,0,fMaxClusterE);
    fhClusMaxMBPure[i] ->SetXTitle("Energy (GeV)");
    fOutputList->Add(fhClusMaxMBPure[i]);
  }
  
  for(Int_t i=0; i < fgkTriggerCombi; i++)
  {
    fhV0[i] = new TH1F(Form("hV0%s",hName[i].Data()),
                       Form("V0 distribution for %s",hTitle[i].Data()),
                       fNBinsV0Signal,0,fMaxV0Signal);
    fhV0[i]->SetXTitle("V0");
    fOutputList->Add(fhV0[i] );
		
    fhClus[i]    = new TH1F(Form("hClus%s",hName[i].Data()),
                            Form("clusters E distribution for %s",hTitle[i].Data()),
                            fNBinsClusterE,0,fMaxClusterE);
    fhClus[i]   ->SetXTitle("Energy (GeV)");
    fOutputList->Add(fhClus[i]);
    
    fhClusMax[i] = new TH1F(Form("hClusMax%s",hName[i].Data()),
                            Form("maximum energy cluster per event for %s",hTitle[i].Data()),
                            fNBinsClusterE,0,fMaxClusterE);
    fhClusMax[i]->SetXTitle("Energy (GeV)");
    fOutputList->Add(fhClusMax[i]);
    
    // Cluster histograms, E vs Cen
    
    fhClusCen[i]    = new TH2F(Form("hClusCen%s",hName[i].Data()),
                               Form("clusters E distribution vs centrality for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
    fhClusCen[i]   ->SetXTitle("Energy (GeV)");
    fhClusCen[i]   ->SetYTitle("Centrality");
    fOutputList->Add(fhClusCen[i]);
    
    fhClusCenMax[i] = new TH2F(Form("hClusCenMax%s",hName[i].Data()),
                               Form("maximum energy cluster per event vs centrality for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
    fhClusCenMax[i]->SetXTitle("Energy (GeV)");
    fhClusCenMax[i]->SetYTitle("Centrality");
    fOutputList->Add(fhClusCenMax[i]);
    
    // Cluster histograms, E vs V0
    
    fhClusV0[i]    = new TH2F(Form("hClusV0%s",hName[i].Data()),
                              Form("clusters E distribution vs V0 for %s",hTitle[i].Data()),
                              fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
    fhClusV0[i]   ->SetXTitle("Energy (GeV)");
    fhClusV0[i]   ->SetYTitle("V0");
    fOutputList->Add(fhClusV0[i]);
    
    fhClusV0Max[i] = new TH2F(Form("hClusV0Max%s",hName[i].Data()),
                              Form("maximum energy cluster per event vs V0 for %s",hTitle[i].Data()),
                              fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
    fhClusV0Max[i]->SetXTitle("Energy (GeV)");
    fhClusV0Max[i]->SetYTitle("V0");
    fOutputList->Add(fhClusV0Max[i]);
    
    // Cluster histograms, E vs Pseudorapidity
    Float_t etamin =-0.8;
    Float_t etamax = 0.8;
    Int_t neta     = 160;
    fhClusEta[i]    = new TH2F(Form("hClusEta%s",hName[i].Data()),
                               Form("clusters distribution vs #eta for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
    fhClusEta[i]   ->SetXTitle("Energy (GeV)");
    fhClusEta[i]   ->SetYTitle("#eta");
    fOutputList->Add(fhClusEta[i]);
    
    fhClusEtaMax[i] = new TH2F(Form("hClusEtaMax%s",hName[i].Data()),
                               Form("maximum energy cluster per event vs #eta for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
    fhClusEtaMax[i]->SetXTitle("Energy (GeV)");
    fhClusEtaMax[i]->SetYTitle("#eta");
    fOutputList->Add(fhClusEtaMax[i]);
    
    // Cluster histograms, E vs Azimuthal angle
    Float_t phimin = 80. *TMath::DegToRad();
    Float_t phimax = 190.*TMath::DegToRad();
    Int_t   nphi   = 110;
    
    fhClusPhi[i]    = new TH2F(Form("hClusPhi%s",hName[i].Data()),
                               Form("clusters distribution vs #phi for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
    fhClusPhi[i]   ->SetXTitle("Energy (GeV)");
    fhClusPhi[i]   ->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusPhi[i]);
    
    fhClusPhiMax[i] = new TH2F(Form("hClusPhiMax%s",hName[i].Data()),
                               Form("maximum energy cluster per event vs #phi for %s",hTitle[i].Data()),
                               fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
    fhClusPhiMax[i]->SetXTitle("Energy (GeV)");
    fhClusPhiMax[i]->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusPhiMax[i]);
    
    // Cluster histograms, Pseudorapidity vs Azimuthal angle
    
    fhClusEtaPhiHigh[i]    = new TH2F(Form("hClusEtaPhiHigh%s",hName[i].Data()),
                                      Form("clusters distribution #eta vs #phi for %s, E > 10 GeV",hTitle[i].Data()),
                                      neta, etamin, etamax,nphi, phimin, phimax);
    fhClusEtaPhiHigh[i]   ->SetXTitle("#eta");
    fhClusEtaPhiHigh[i]   ->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiHigh[i]);
    
    fhClusEtaPhiHighCluMax[i] = new TH2F(Form("hClusEtaPhiHighCluMax%s",hName[i].Data()),
                                         Form("maximum energy cluster per event #eta  vs #phi for %s, E > 10 GeV",hTitle[i].Data()),
                                         neta, etamin, etamax,nphi, phimin, phimax);
    fhClusEtaPhiHighCluMax[i]->SetXTitle("#eta");
    fhClusEtaPhiHighCluMax[i]->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiHighCluMax[i]);
    
    fhClusEtaPhiLow[i]    = new TH2F(Form("hClusEtaPhiLow%s",hName[i].Data()),
                                     Form("clusters distribution #eta vs #phi for %s, E < 10 GeV",hTitle[i].Data()),
                                     neta, etamin, etamax,nphi, phimin, phimax);
    fhClusEtaPhiLow[i]   ->SetXTitle("#eta");
    fhClusEtaPhiLow[i]   ->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiLow[i]);
    
    fhClusEtaPhiLowCluMax[i] = new TH2F(Form("hClusEtaPhiLowCluMax%s",hName[i].Data()),
                                        Form("maximum energy cluster per event #eta  vs #phi for %s, E < 10 GeV",hTitle[i].Data()),
                                        neta, etamin, etamax,nphi, phimin, phimax);
    fhClusEtaPhiLowCluMax[i]->SetXTitle("#eta");
    fhClusEtaPhiLowCluMax[i]->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiLowCluMax[i]);
    
    fhClusEtaPhiHighCellMax[i]    = new TH2F(Form("hClusEtaPhiHighCellMax%s",hName[i].Data()),
                                             Form("Cluster hit map in calorimeter (max cell), column vs row for %s, E > 10 GeV",hTitle[i].Data()),
                                             fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
    fhClusEtaPhiHighCellMax[i]   ->SetXTitle("Index #eta (columnns)");
    fhClusEtaPhiHighCellMax[i]   ->SetYTitle("Index #phi (rows)");
    fOutputList->Add(fhClusEtaPhiHighCellMax[i]);
    
    fhClusEtaPhiHighCellMaxCluMax[i] = new TH2F(Form("hClusEtaPhiHighCellMaxCluMax%s",hName[i].Data()),
                                                Form("Max E cluster hit map in calorimeter (max cell), column vs row  for %s, E > 10 GeV",
                                                     hTitle[i].Data()),fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
    fhClusEtaPhiHighCellMaxCluMax[i]->SetXTitle("Index #eta (columnns)");
    fhClusEtaPhiHighCellMaxCluMax[i]->SetYTitle("Index #phi (rows)");
    fOutputList->Add(fhClusEtaPhiHighCellMaxCluMax[i]);
    
    fhClusEtaPhiLowCellMax[i]    = new TH2F(Form("hClusEtaPhiLowCellMax%s",hName[i].Data()),
                                            Form("Cluster hit map in calorimeter (max cell), column vs row for %s, E < 10 GeV",hTitle[i].Data()),
                                            fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
    fhClusEtaPhiLowCellMax[i]   ->SetXTitle("Index #eta (columnns)");
    fhClusEtaPhiLowCellMax[i]   ->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiLowCellMax[i]);
    
    fhClusEtaPhiLowCellMaxCluMax[i] = new TH2F(Form("hClusEtaPhiLowCellMaxCluMax%s",hName[i].Data()),
                                               Form("Max E cluster hit map in calorimeter (max cell), column vs row  for %s, E < 10 GeV",
                                                    hTitle[i].Data()),fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
    fhClusEtaPhiLowCellMaxCluMax[i]->SetXTitle("Index #eta (columnns)");
    fhClusEtaPhiLowCellMaxCluMax[i]->SetYTitle("#phi (rad)");
    fOutputList->Add(fhClusEtaPhiLowCellMaxCluMax[i]);
  }
  
  PostData(1, fOutputList);
}

//______________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserExec(Option_t *)
{
  // Main loop
	
  AliVEvent* event = InputEvent();
  
  if (!event)
  {
    AliError("No Event, exit");
    return;
  }
  
  InitGeometry(); // only once, must be done before OADB, geo OADB accessed here
  
  if(fAccessOADB) AccessOADB(); // only once
  
  InitCellPatchMaps();   //init to 0 map for cells and patches
  
  //trigger configuration
  TString triggerclasses = event->GetFiredTriggerClasses();
  
  Int_t eventType = ((AliVHeader*)event->GetHeader())->GetEventType();
  //std::cout << "trigger = " << triggerclasses << std::endl;

  // physics events eventType=7, select only those
  if(triggerclasses=="" || eventType != 7) return;
  
  SetTriggerEventBit(triggerclasses);
		
  FillEventCounterHistogram();
  
  FillCellMaps();
  
  FillTriggerPatchMaps(triggerclasses);
  
  FillMapHistograms();
  
  FillV0Histograms();
  
  FillL1GammaPatchHistograms();
  
  FillL1JetPatchHistograms();

  // FEE vs STU vs TRU
  FillCorrelationHistograms();
	
  ClusterAnalysis();
  
  PostData(1, fOutputList);
  
}
