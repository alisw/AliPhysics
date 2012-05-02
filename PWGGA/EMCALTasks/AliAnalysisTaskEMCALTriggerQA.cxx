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
//  $Id$                                                                //
//------------------------------------------------------------------------//

#include <TList.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile2D.h> 

#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

#include "AliAnalysisTaskEMCALTriggerQA.h"

ClassImp(AliAnalysisTaskEMCALTriggerQA)

//______________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA() : 
AliAnalysisTaskSE(), 
fOutputList(0),            fRecoUtils(0x0),
fGeometry(0),              fGeoName(""),         
fhNEvents(0),              fhFORAmp(0),
fhFORAmpL1G(0),            fhFORAmpL1J(0),
fhL0Amp(0),                fhL0AmpL1G(0),        fhL0AmpL1J(0),
fhL1Amp(0),                fhL1GAmp(0),          fhL1JAmp(0),
fhL0Patch(0),              fhL1GPatch(0),        fhL1JPatch(0),
fhFEESTU(0),               fhTRUSTU(0),          fhV0STU(0),
fhGPMaxVV0TT(0),           fhJPMaxVV0TT(0),
fhFORMeanAmp(0),           fhL0MeanAmp(0),       fhL1MeanAmp(0),
fhV0MB(0),                 fhV0L1G(0),           fhV0L1J(0),
fhL1GPatchMax(0),          fhL1JPatchMax(0),

fhClusMB(0),               fhClusMBPure(0),      fhClusL0(0),
fhClusL1G(0),              fhClusL1J(0),
fhClusL1GOnly(0),          fhClusL1JOnly(0),
fhClusMaxMB(0),            fhClusMaxMBPure(0),   fhClusMaxL0(0),
fhClusMaxL1G(0),           fhClusMaxL1J(0),
fhClusMaxL1GOnly(0),       fhClusMaxL1JOnly(0),

fhClusCenMB(0),            fhClusCenL0(0),
fhClusCenL1G(0),           fhClusCenL1J(0),
fhClusCenL1GOnly(0),       fhClusCenL1JOnly(0),
fhClusCenMaxMB(0),         fhClusCenMaxL0(0),
fhClusCenMaxL1G(0),        fhClusCenMaxL1J(0),
fhClusCenMaxL1GOnly(0),    fhClusCenMaxL1JOnly(0),

fhClusV0MB(0),             fhClusV0L0(0),
fhClusV0L1G(0),            fhClusV0L1J(0),
fhClusV0L1GOnly(0),        fhClusV0L1JOnly(0),
fhClusV0MaxMB(0),          fhClusV0MaxL0(0),
fhClusV0MaxL1G(0),         fhClusV0MaxL1J(0),
fhClusV0MaxL1GOnly(0),     fhClusV0MaxL1JOnly(0),

fhClusEtaMB(0),            fhClusEtaL0(0),
fhClusEtaL1G(0),           fhClusEtaL1J(0),
fhClusEtaL1GOnly(0),       fhClusEtaL1JOnly(0),
fhClusEtaMaxMB(0),         fhClusEtaMaxL0(0),
fhClusEtaMaxL1G(0),        fhClusEtaMaxL1J(0),
fhClusEtaMaxL1GOnly(0),    fhClusEtaMaxL1JOnly(0),

fhClusPhiMB(0),            fhClusPhiL0(0),
fhClusPhiL1G(0),           fhClusPhiL1J(0),
fhClusPhiL1GOnly(0),       fhClusPhiL1JOnly(0),
fhClusPhiMaxMB(0),         fhClusPhiMaxL0(0),
fhClusPhiMaxL1G(0),        fhClusPhiMaxL1J(0),
fhClusPhiMaxL1GOnly(0),    fhClusPhiMaxL1JOnly(0),

fhClusEtaPhiHighMB(0),            fhClusEtaPhiHighL0(0),
fhClusEtaPhiHighL1G(0),           fhClusEtaPhiHighL1J(0),
fhClusEtaPhiHighL1GOnly(0),       fhClusEtaPhiHighL1JOnly(0),
fhClusEtaPhiHighCluMaxMB(0),      fhClusEtaPhiHighCluMaxL0(0),
fhClusEtaPhiHighCluMaxL1G(0),     fhClusEtaPhiHighCluMaxL1J(0),
fhClusEtaPhiHighCluMaxL1GOnly(0), fhClusEtaPhiHighCluMaxL1JOnly(0),

fhClusEtaPhiHighCellMaxMB(0),            fhClusEtaPhiHighCellMaxL0(0),
fhClusEtaPhiHighCellMaxL1G(0),           fhClusEtaPhiHighCellMaxL1J(0),
fhClusEtaPhiHighCellMaxL1GOnly(0),       fhClusEtaPhiHighCellMaxL1JOnly(0),
fhClusEtaPhiHighCellMaxCluMaxMB(0),      fhClusEtaPhiHighCellMaxCluMaxL0(0),
fhClusEtaPhiHighCellMaxCluMaxL1G(0),     fhClusEtaPhiHighCellMaxCluMaxL1J(0),
fhClusEtaPhiHighCellMaxCluMaxL1GOnly(0), fhClusEtaPhiHighCellMaxCluMaxL1JOnly(0),

fhClusEtaPhiLowMB(0),            fhClusEtaPhiLowL0(0),
fhClusEtaPhiLowL1G(0),           fhClusEtaPhiLowL1J(0),
fhClusEtaPhiLowL1GOnly(0),       fhClusEtaPhiLowL1JOnly(0),
fhClusEtaPhiLowCluMaxMB(0),      fhClusEtaPhiLowCluMaxL0(0),
fhClusEtaPhiLowCluMaxL1G(0),     fhClusEtaPhiLowCluMaxL1J(0),
fhClusEtaPhiLowCluMaxL1GOnly(0), fhClusEtaPhiLowCluMaxL1JOnly(0),

fhClusEtaPhiLowCellMaxMB(0),            fhClusEtaPhiLowCellMaxL0(0),
fhClusEtaPhiLowCellMaxL1G(0),           fhClusEtaPhiLowCellMaxL1J(0),
fhClusEtaPhiLowCellMaxL1GOnly(0),       fhClusEtaPhiLowCellMaxL1JOnly(0),
fhClusEtaPhiLowCellMaxCluMaxMB(0),      fhClusEtaPhiLowCellMaxCluMaxL0(0),
fhClusEtaPhiLowCellMaxCluMaxL1G(0),     fhClusEtaPhiLowCellMaxCluMaxL1J(0),
fhClusEtaPhiLowCellMaxCluMaxL1GOnly(0), fhClusEtaPhiLowCellMaxCluMaxL1JOnly(0),

//Histogram settings
fNBinsSTUSignal  (2000),   fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000),   fMaxTRUSignal  (200000),
fNBinsV0Signal   (5000),   fMaxV0Signal   (100000),
fNBinsSTUFEERatio(1000),   fMaxSTUFEERatio(100),
fNBinsSTUTRURatio(1000),   fMaxSTUTRURatio(100),
fNBinsClusterE   (500),    fMaxClusterE   (200)

{
  // Constructor
  
  fGeoName   = "EMCAL_COMPLETEV1"; 
  fRecoUtils = new AliEMCALRecoUtils;
  
  DefineOutput(1, TList::Class());
  
}		      

//______________________________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA(const char *name) : 
AliAnalysisTaskSE(name), 
fOutputList(0),            fRecoUtils(0x0),
fGeometry(0),              fGeoName(""),         
fhNEvents(0),              fhFORAmp(0),
fhFORAmpL1G(0),            fhFORAmpL1J(0),
fhL0Amp(0),                fhL0AmpL1G(0),        fhL0AmpL1J(0),
fhL1Amp(0),                fhL1GAmp(0),          fhL1JAmp(0),
fhL0Patch(0),              fhL1GPatch(0),        fhL1JPatch(0),
fhFEESTU(0),               fhTRUSTU(0),          fhV0STU(0),
fhGPMaxVV0TT(0),           fhJPMaxVV0TT(0),
fhFORMeanAmp(0),           fhL0MeanAmp(0),       fhL1MeanAmp(0),
fhV0MB(0),                 fhV0L1G(0),           fhV0L1J(0),
fhL1GPatchMax(0),          fhL1JPatchMax(0),

fhClusMB(0),               fhClusMBPure(0),      fhClusL0(0),
fhClusL1G(0),              fhClusL1J(0),
fhClusL1GOnly(0),          fhClusL1JOnly(0),
fhClusMaxMB(0),            fhClusMaxMBPure(0),   fhClusMaxL0(0),
fhClusMaxL1G(0),           fhClusMaxL1J(0),
fhClusMaxL1GOnly(0),       fhClusMaxL1JOnly(0),

fhClusCenMB(0),            fhClusCenL0(0),
fhClusCenL1G(0),           fhClusCenL1J(0),
fhClusCenL1GOnly(0),       fhClusCenL1JOnly(0),
fhClusCenMaxMB(0),         fhClusCenMaxL0(0),
fhClusCenMaxL1G(0),        fhClusCenMaxL1J(0),
fhClusCenMaxL1GOnly(0),    fhClusCenMaxL1JOnly(0),

fhClusV0MB(0),             fhClusV0L0(0),
fhClusV0L1G(0),            fhClusV0L1J(0),
fhClusV0L1GOnly(0),        fhClusV0L1JOnly(0),
fhClusV0MaxMB(0),          fhClusV0MaxL0(0),
fhClusV0MaxL1G(0),         fhClusV0MaxL1J(0),
fhClusV0MaxL1GOnly(0),     fhClusV0MaxL1JOnly(0),

fhClusEtaMB(0),            fhClusEtaL0(0),
fhClusEtaL1G(0),           fhClusEtaL1J(0),
fhClusEtaL1GOnly(0),       fhClusEtaL1JOnly(0),
fhClusEtaMaxMB(0),         fhClusEtaMaxL0(0),
fhClusEtaMaxL1G(0),        fhClusEtaMaxL1J(0),
fhClusEtaMaxL1GOnly(0),    fhClusEtaMaxL1JOnly(0),

fhClusPhiMB(0),            fhClusPhiL0(0),
fhClusPhiL1G(0),           fhClusPhiL1J(0),
fhClusPhiL1GOnly(0),       fhClusPhiL1JOnly(0),
fhClusPhiMaxMB(0),         fhClusPhiMaxL0(0),
fhClusPhiMaxL1G(0),        fhClusPhiMaxL1J(0),
fhClusPhiMaxL1GOnly(0),    fhClusPhiMaxL1JOnly(0),

fhClusEtaPhiHighMB(0),            fhClusEtaPhiHighL0(0),
fhClusEtaPhiHighL1G(0),           fhClusEtaPhiHighL1J(0),
fhClusEtaPhiHighL1GOnly(0),       fhClusEtaPhiHighL1JOnly(0),
fhClusEtaPhiHighCluMaxMB(0),      fhClusEtaPhiHighCluMaxL0(0),
fhClusEtaPhiHighCluMaxL1G(0),     fhClusEtaPhiHighCluMaxL1J(0),
fhClusEtaPhiHighCluMaxL1GOnly(0), fhClusEtaPhiHighCluMaxL1JOnly(0),

fhClusEtaPhiHighCellMaxMB(0),            fhClusEtaPhiHighCellMaxL0(0),
fhClusEtaPhiHighCellMaxL1G(0),           fhClusEtaPhiHighCellMaxL1J(0),
fhClusEtaPhiHighCellMaxL1GOnly(0),       fhClusEtaPhiHighCellMaxL1JOnly(0),
fhClusEtaPhiHighCellMaxCluMaxMB(0),      fhClusEtaPhiHighCellMaxCluMaxL0(0),
fhClusEtaPhiHighCellMaxCluMaxL1G(0),     fhClusEtaPhiHighCellMaxCluMaxL1J(0),
fhClusEtaPhiHighCellMaxCluMaxL1GOnly(0), fhClusEtaPhiHighCellMaxCluMaxL1JOnly(0),

fhClusEtaPhiLowMB(0),            fhClusEtaPhiLowL0(0),
fhClusEtaPhiLowL1G(0),           fhClusEtaPhiLowL1J(0),
fhClusEtaPhiLowL1GOnly(0),       fhClusEtaPhiLowL1JOnly(0),
fhClusEtaPhiLowCluMaxMB(0),      fhClusEtaPhiLowCluMaxL0(0),
fhClusEtaPhiLowCluMaxL1G(0),     fhClusEtaPhiLowCluMaxL1J(0),
fhClusEtaPhiLowCluMaxL1GOnly(0), fhClusEtaPhiLowCluMaxL1JOnly(0),

fhClusEtaPhiLowCellMaxMB(0),            fhClusEtaPhiLowCellMaxL0(0),
fhClusEtaPhiLowCellMaxL1G(0),           fhClusEtaPhiLowCellMaxL1J(0),
fhClusEtaPhiLowCellMaxL1GOnly(0),       fhClusEtaPhiLowCellMaxL1JOnly(0),
fhClusEtaPhiLowCellMaxCluMaxMB(0),      fhClusEtaPhiLowCellMaxCluMaxL0(0),
fhClusEtaPhiLowCellMaxCluMaxL1G(0),     fhClusEtaPhiLowCellMaxCluMaxL1J(0),
fhClusEtaPhiLowCellMaxCluMaxL1GOnly(0), fhClusEtaPhiLowCellMaxCluMaxL1JOnly(0),

//Histogram settings
fNBinsSTUSignal  (2000),   fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000),   fMaxTRUSignal  (200000),
fNBinsV0Signal   (5000),   fMaxV0Signal   (100000),
fNBinsSTUFEERatio(1000),   fMaxSTUFEERatio(100),
fNBinsSTUTRURatio(1000),   fMaxSTUTRURatio(100),
fNBinsClusterE   (500),    fMaxClusterE   (200)

{
  // Constructor
  
  fGeoName   = "EMCAL_COMPLETEV1"; 
  fRecoUtils = new AliEMCALRecoUtils;
  
  DefineOutput(1, TList::Class());
  
}


//___________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserCreateOutputObjects() 
{
  // Init histograms and geometry 
  
  fGeometry = AliEMCALGeometry::GetInstance(fGeoName);
  
  fOutputList  = new TList;
  fOutputList ->SetOwner(kTRUE);
  
  fhNEvents    = new TH1F("hNEvents","Number of selected events",12,0,12);
  fhNEvents   ->SetYTitle("N events");
  fhNEvents   ->GetXaxis()->SetBinLabel(1 ,"All");
  fhNEvents   ->GetXaxis()->SetBinLabel(2 ,"MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(3 ,"L0");
  fhNEvents   ->GetXaxis()->SetBinLabel(4 ,"L1-G");
  fhNEvents   ->GetXaxis()->SetBinLabel(5 ,"L1-J");
  fhNEvents   ->GetXaxis()->SetBinLabel(6 ,"L1-G & !L1-J");
  fhNEvents   ->GetXaxis()->SetBinLabel(7 ,"L1-J & !L1-G");
  fhNEvents   ->GetXaxis()->SetBinLabel(8 ,"L1-J & L1-G");  
  fhNEvents   ->GetXaxis()->SetBinLabel(9 ,"MB & !L1 & !L0");
  fhNEvents   ->GetXaxis()->SetBinLabel(10,"L0 & !MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(11,"L1-G & !MB");
  fhNEvents   ->GetXaxis()->SetBinLabel(12,"L1-J & !MB");

  
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
  
  fhFORAmpL1J  = new TH2F("hFORAmpL1J", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1J trigger condition",
                          fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1J ->SetXTitle("Index #eta (columnns)");
  fhFORAmpL1J ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1J ->SetZTitle("Amplitude");
  
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
  
  fhL1JAmp     = new TH2F("hL1JAmp","STU signal per Row and Column for L1 Jet",
                          fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JAmp    ->SetXTitle("Index #eta (columnns)");
  fhL1JAmp    ->SetYTitle("Index #phi (rows)");
  fhL1JAmp    ->SetZTitle("Amplitude");
  
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
  
  fhL1JPatch   = new TH2F("hL1JPatch","FOR with associated L1 Jet Patch",
                          fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JPatch  ->SetXTitle("Index #eta (columnns)");
  fhL1JPatch  ->SetYTitle("Index #phi (rows)");
  fhL1JPatch  ->SetZTitle("counts");
  
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
  
  fhV0MB = new TH1F("hV0MB","V0 distribution for MB triggered event",fNBinsV0Signal,0,fMaxV0Signal);
  fhV0MB->SetXTitle("V0");
  
  fhV0L1G = new TH1F("hV0L1G","V0 distribution for L1G triggered event",fNBinsV0Signal,0,fMaxV0Signal);
  fhV0L1G->SetXTitle("V0");
  
  fhV0L1J = new TH1F("hV0L1J","V0 distribution for L1J triggered event",fNBinsV0Signal,0,fMaxV0Signal);
  fhV0L1J->SetXTitle("V0");
  
  fhL1GPatchMax   = new TH2F("hL1GPatchMax","FOR of max amplitude patch with associated L1 Gamma Patch",
                             fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1GPatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1GPatchMax  ->SetZTitle("counts");
  
  fhL1JPatchMax   = new TH2F("hL1JPatchMax","FOR of max amplitude patch with associated L1 Jet Patch",
                             fgkFALTROCols/4,0,fgkFALTROCols,fgkFALTRORows/4,0,fgkFALTRORows);
  fhL1JPatchMax  ->SetXTitle("Index #eta (columnns)");
  fhL1JPatchMax  ->SetYTitle("Index #phi (rows)");
  fhL1JPatchMax  ->SetZTitle("counts");
  
  fOutputList->Add(fhNEvents);
  fOutputList->Add(fhV0STU);
  fOutputList->Add(fhFORAmp);
  fOutputList->Add(fhFORAmpL1G);
  fOutputList->Add(fhFORAmpL1J);
  fOutputList->Add(fhL0Amp);
  fOutputList->Add(fhL0AmpL1G);
  fOutputList->Add(fhL0AmpL1J);
  fOutputList->Add(fhL1Amp);
  fOutputList->Add(fhL1GAmp);
  fOutputList->Add(fhL1JAmp);
  fOutputList->Add(fhL0Patch);
  fOutputList->Add(fhL1GPatch);
  fOutputList->Add(fhL1JPatch);
  fOutputList->Add(fhFEESTU);
  fOutputList->Add(fhTRUSTU);
  
  fOutputList->Add(fhGPMaxVV0TT);
  fOutputList->Add(fhJPMaxVV0TT);
  
  fOutputList->Add(fhFORMeanAmp);
  fOutputList->Add(fhL0MeanAmp);
  fOutputList->Add(fhL1MeanAmp);
  
  fOutputList->Add(fhV0MB );
  fOutputList->Add(fhV0L1G);
  fOutputList->Add(fhV0L1J);
  
  fOutputList->Add(fhL1GPatchMax);
  fOutputList->Add(fhL1JPatchMax);
  
  // Cluster histograms, E
  
  fhClusMB     = new TH1F("hClusMB","clusters E distribution for MB trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusMB    ->SetXTitle("Energy (GeV)");
  
  fhClusMBPure  = new TH1F("hClusMBPure","clusters E distribution for MB trigger, no other EMCAL trigger on",fNBinsClusterE,0,fMaxClusterE);
  fhClusMBPure ->SetXTitle("Energy (GeV)");
  
  fhClusL0     = new TH1F("hClusL0","clusters E distribution for L0 trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusL0    ->SetXTitle("Energy (GeV)");
  
  fhClusL1G    = new TH1F("hClusL1G","clusters E distribution for L1G trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusL1G   ->SetXTitle("Energy (GeV)");
  
  fhClusL1J    = new TH1F("hClusL1J","clusters E distribution for L1J trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusL1J   ->SetXTitle("Energy (GeV)");
  
  fhClusL1GOnly = new TH1F("hClusL1GOnly","clusters E distribution for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE);
  fhClusL1GOnly->SetXTitle("Energy (GeV)");
  
  fhClusL1JOnly = new TH1F("hClusL1JOnly","clusters E distribution for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE);
  fhClusL1JOnly->SetXTitle("Energy (GeV)");
  
  fhClusMaxMB  = new TH1F("hClusMaxMB","maximum energy cluster per event for MB trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxMB ->SetXTitle("Energy (GeV)");
  
  fhClusMaxMBPure  = new TH1F("hClusMaxMBPure","maximum energy cluster per event for MB trigger, no other EMCAL trigger on",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxMBPure ->SetXTitle("Energy (GeV)");
  
  fhClusMaxL0  = new TH1F("hClusMaxL0","maximum energy cluster per event for L0 trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxL0 ->SetXTitle("Energy (GeV)");
  
  fhClusMaxL1G = new TH1F("hClusMaxL1G","maximum energy cluster per event for L1G trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxL1G->SetXTitle("Energy (GeV)");
  
  fhClusMaxL1J = new TH1F("hClusMaxL1J","maximum energy cluster per event for L1J trigger",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxL1J->SetXTitle("Energy (GeV)");
  
  fhClusMaxL1GOnly = new TH1F("hClusMaxL1GOnly","maximum energy cluster per event for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxL1GOnly->SetXTitle("Energy (GeV)");
  
  fhClusMaxL1JOnly = new TH1F("hClusMaxL1JOnly","maximum energy cluster per event for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE);
  fhClusMaxL1JOnly->SetXTitle("Energy (GeV)");
  
  
  fOutputList->Add(fhClusMB);
  fOutputList->Add(fhClusMBPure);
  fOutputList->Add(fhClusL0);
  fOutputList->Add(fhClusL1G);
  fOutputList->Add(fhClusL1J);
  fOutputList->Add(fhClusL1GOnly);
  fOutputList->Add(fhClusL1JOnly);
  
  fOutputList->Add(fhClusMaxMB);
  fOutputList->Add(fhClusMaxMBPure);
  fOutputList->Add(fhClusMaxL0);
  fOutputList->Add(fhClusMaxL1G);
  fOutputList->Add(fhClusMaxL1J);
  fOutputList->Add(fhClusMaxL1GOnly);
  fOutputList->Add(fhClusMaxL1JOnly);
  
  // Cluster histograms, E vs Cen

  fhClusCenMB     = new TH2F("hClusCenMB","clusters E distribution vs centrality for MB trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMB    ->SetXTitle("Energy (GeV)");
  fhClusCenMB    ->SetYTitle("Centrality");
  
  fhClusCenL0     = new TH2F("hClusCenL0","clusters E distribution vs centrality for L0 trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenL0    ->SetXTitle("Energy (GeV)");
  fhClusCenL0    ->SetYTitle("Centrality");

  fhClusCenL1G    = new TH2F("hClusCenL1G","clusters E distribution vs centrality for L1G trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenL1G   ->SetXTitle("Energy (GeV)");
  
  fhClusCenL1J    = new TH2F("hClusCenL1J","clusters E distribution vs centrality for L1J trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenL1J   ->SetXTitle("Energy (GeV)");
  fhClusCenL1J   ->SetYTitle("Centrality");

  fhClusCenL1GOnly = new TH2F("hClusCenL1GOnly","clusters E distribution vs centrality for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenL1GOnly->SetXTitle("Energy (GeV)");
  fhClusCenL1GOnly->SetYTitle("Centrality");

  fhClusCenL1JOnly = new TH2F("hClusCenL1JOnly","clusters E distribution vs centrality for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenL1JOnly->SetXTitle("Energy (GeV)");
  fhClusCenL1JOnly->SetYTitle("Centrality");

  fhClusCenMaxMB  = new TH2F("hClusCenMaxMB","maximum energy cluster per event vs centrality for MB trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxMB ->SetXTitle("Energy (GeV)");
  fhClusCenMaxMB ->SetYTitle("Centrality");

  fhClusCenMaxL0  = new TH2F("hClusCenMaxL0","maximum energy cluster per event vs centrality for L0 trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxL0 ->SetXTitle("Energy (GeV)");
  fhClusCenMaxL0 ->SetYTitle("Centrality");

  fhClusCenMaxL1G = new TH2F("hClusCenMaxL1G","maximum energy cluster per event vs centrality for L1G trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxL1G->SetXTitle("Energy (GeV)");
  fhClusCenMaxL1G->SetYTitle("Centrality");

  fhClusCenMaxL1J = new TH2F("hClusCenMaxL1J","maximum energy cluster per event vs centrality for L1J trigger",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxL1J->SetXTitle("Energy (GeV)");
  fhClusCenMaxL1J->SetYTitle("Centrality");

  fhClusCenMaxL1GOnly = new TH2F("hClusCenMaxL1GOnly","maximum energy cluster per event vs centrality for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxL1GOnly->SetXTitle("Energy (GeV)");
  fhClusCenMaxL1GOnly->SetYTitle("Centrality");

  fhClusCenMaxL1JOnly = new TH2F("hClusCenMaxL1JOnly","maximum energy cluster per event vs centrality for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,100, 0, 100);
  fhClusCenMaxL1JOnly->SetXTitle("Energy (GeV)");
  fhClusCenMaxL1JOnly->SetYTitle("Centrality");

  
  fOutputList->Add(fhClusCenMB);
  fOutputList->Add(fhClusCenL0);
  fOutputList->Add(fhClusCenL1G);
  fOutputList->Add(fhClusCenL1J);
  fOutputList->Add(fhClusCenL1GOnly);
  fOutputList->Add(fhClusCenL1JOnly);
  
  fOutputList->Add(fhClusCenMaxMB);
  fOutputList->Add(fhClusCenMaxL0);
  fOutputList->Add(fhClusCenMaxL1G);
  fOutputList->Add(fhClusCenMaxL1J);
  fOutputList->Add(fhClusCenMaxL1GOnly);
  fOutputList->Add(fhClusCenMaxL1JOnly);
  
  
  // Cluster histograms, E vs V0
  
  fhClusV0MB     = new TH2F("hClusV0MB","clusters E distribution vs V0 for MB trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MB    ->SetXTitle("Energy (GeV)");
  fhClusV0MB    ->SetYTitle("V0");
  
  fhClusV0L0     = new TH2F("hClusV0L0","clusters E distribution vs V0 for L0 trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0L0    ->SetXTitle("Energy (GeV)");
  fhClusV0L0    ->SetYTitle("V0");
  
  fhClusV0L1G    = new TH2F("hClusV0L1G","clusters E distribution vs V0 for L1G trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0L1G   ->SetXTitle("Energy (GeV)");
  
  fhClusV0L1J    = new TH2F("hClusV0L1J","clusters E distribution vs V0 for L1J trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0L1J   ->SetXTitle("Energy (GeV)");
  fhClusV0L1J   ->SetYTitle("V0");
  
  fhClusV0L1GOnly = new TH2F("hClusV0L1GOnly","clusters E distribution vs V0 for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0L1GOnly->SetXTitle("Energy (GeV)");
  fhClusV0L1GOnly->SetYTitle("V0");
  
  fhClusV0L1JOnly = new TH2F("hClusV0L1JOnly","clusters E distribution vs V0 for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0L1JOnly->SetXTitle("Energy (GeV)");
  fhClusV0L1JOnly->SetYTitle("V0");
  
  fhClusV0MaxMB  = new TH2F("hClusV0MaxMB","maximum energy cluster per event vs V0 for MB trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxMB ->SetXTitle("Energy (GeV)");
  fhClusV0MaxMB ->SetYTitle("V0");
  
  fhClusV0MaxL0  = new TH2F("hClusV0MaxL0","maximum energy cluster per event vs V0 for L0 trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxL0 ->SetXTitle("Energy (GeV)");
  fhClusV0MaxL0 ->SetYTitle("V0");
  
  fhClusV0MaxL1G = new TH2F("hClusV0MaxL1G","maximum energy cluster per event vs V0 for L1G trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxL1G->SetXTitle("Energy (GeV)");
  fhClusV0MaxL1G->SetYTitle("V0");
  
  fhClusV0MaxL1J = new TH2F("hClusV0MaxL1J","maximum energy cluster per event vs V0 for L1J trigger",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxL1J->SetXTitle("Energy (GeV)");
  fhClusV0MaxL1J->SetYTitle("V0");
  
  fhClusV0MaxL1GOnly = new TH2F("hClusV0MaxL1GOnly","maximum energy cluster per event vs V0 for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxL1GOnly->SetXTitle("Energy (GeV)");
  fhClusV0MaxL1GOnly->SetYTitle("V0");
  
  fhClusV0MaxL1JOnly = new TH2F("hClusV0MaxL1JOnly","maximum energy cluster per event vs V0 for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,fNBinsV0Signal,0,fMaxV0Signal);
  fhClusV0MaxL1JOnly->SetXTitle("Energy (GeV)");
  fhClusV0MaxL1JOnly->SetYTitle("V0");
  
  
  fOutputList->Add(fhClusV0MB);
  fOutputList->Add(fhClusV0L0);
  fOutputList->Add(fhClusV0L1G);
  fOutputList->Add(fhClusV0L1J);
  fOutputList->Add(fhClusV0L1GOnly);
  fOutputList->Add(fhClusV0L1JOnly);
  
  fOutputList->Add(fhClusV0MaxMB);
  fOutputList->Add(fhClusV0MaxL0);
  fOutputList->Add(fhClusV0MaxL1G);
  fOutputList->Add(fhClusV0MaxL1J);
  fOutputList->Add(fhClusV0MaxL1GOnly);
  fOutputList->Add(fhClusV0MaxL1JOnly);

  // Cluster histograms, E vs Pseudorapidity
  Float_t etamin =-0.8;
  Float_t etamax = 0.8;
  Int_t neta     = 160;
  fhClusEtaMB     = new TH2F("hClusEtaMB","clusters distribution vs #eta for MB trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMB    ->SetXTitle("Energy (GeV)");
  fhClusEtaMB    ->SetYTitle("#eta");
  
  fhClusEtaL0     = new TH2F("hClusEtaL0","clusters distribution vs #eta for L0 trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaL0    ->SetXTitle("Energy (GeV)");
  fhClusEtaL0    ->SetYTitle("#eta");
  
  fhClusEtaL1G    = new TH2F("hClusEtaL1G","clusters distribution vs #eta for L1G trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaL1G   ->SetXTitle("Energy (GeV)");
  
  fhClusEtaL1J    = new TH2F("hClusEtaL1J","clusters distribution vs #eta for L1J trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaL1J   ->SetXTitle("Energy (GeV)");
  fhClusEtaL1J   ->SetYTitle("#eta");
  
  fhClusEtaL1GOnly = new TH2F("hClusEtaL1GOnly","clusters distribution vs #eta for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaL1GOnly->SetXTitle("Energy (GeV)");
  fhClusEtaL1GOnly->SetYTitle("#eta");
  
  fhClusEtaL1JOnly = new TH2F("hClusEtaL1JOnly","clusters distribution vs #eta for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaL1JOnly->SetXTitle("Energy (GeV)");
  fhClusEtaL1JOnly->SetYTitle("#eta");
  
  fhClusEtaMaxMB  = new TH2F("hClusEtaMaxMB","maximum energy cluster per event vs #eta for MB trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxMB ->SetXTitle("Energy (GeV)");
  fhClusEtaMaxMB ->SetYTitle("#eta");
  
  fhClusEtaMaxL0  = new TH2F("hClusEtaMaxL0","maximum energy cluster per event vs #eta for L0 trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxL0 ->SetXTitle("Energy (GeV)");
  fhClusEtaMaxL0 ->SetYTitle("#eta");
  
  fhClusEtaMaxL1G = new TH2F("hClusEtaMaxL1G","maximum energy cluster per event vs #eta for L1G trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxL1G->SetXTitle("Energy (GeV)");
  fhClusEtaMaxL1G->SetYTitle("#eta");
  
  fhClusEtaMaxL1J = new TH2F("hClusEtaMaxL1J","maximum energy cluster per event vs #eta for L1J trigger",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxL1J->SetXTitle("Energy (GeV)");
  fhClusEtaMaxL1J->SetYTitle("#eta");
  
  fhClusEtaMaxL1GOnly = new TH2F("hClusEtaMaxL1GOnly","maximum energy cluster per event vs #eta for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxL1GOnly->SetXTitle("Energy (GeV)");
  fhClusEtaMaxL1GOnly->SetYTitle("#eta");
  
  fhClusEtaMaxL1JOnly = new TH2F("hClusEtaMaxL1JOnly","maximum energy cluster per event vs #eta for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,neta, etamin, etamax);
  fhClusEtaMaxL1JOnly->SetXTitle("Energy (GeV)");
  fhClusEtaMaxL1JOnly->SetYTitle("#eta");
  
  
  fOutputList->Add(fhClusEtaMB);
  fOutputList->Add(fhClusEtaL0);
  fOutputList->Add(fhClusEtaL1G);
  fOutputList->Add(fhClusEtaL1J);
  fOutputList->Add(fhClusEtaL1GOnly);
  fOutputList->Add(fhClusEtaL1JOnly);
  
  fOutputList->Add(fhClusEtaMaxMB);
  fOutputList->Add(fhClusEtaMaxL0);
  fOutputList->Add(fhClusEtaMaxL1G);
  fOutputList->Add(fhClusEtaMaxL1J);
  fOutputList->Add(fhClusEtaMaxL1GOnly);
  fOutputList->Add(fhClusEtaMaxL1JOnly);
  
 
  // Cluster histograms, E vs Azimuthal angle
  Float_t phimin = 80. *TMath::DegToRad();
  Float_t phimax = 190.*TMath::DegToRad();
  Int_t   nphi   = 110;
  
  fhClusPhiMB     = new TH2F("hClusPhiMB","clusters distribution vs #phi for MB trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMB    ->SetXTitle("Energy (GeV)");
  fhClusPhiMB    ->SetYTitle("#phi (rad)");
  
  fhClusPhiL0     = new TH2F("hClusPhiL0","clusters distribution vs #phi for L0 trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiL0    ->SetXTitle("Energy (GeV)");
  fhClusPhiL0    ->SetYTitle("#phi (rad)");
  
  fhClusPhiL1G    = new TH2F("hClusPhiL1G","clusters distribution vs #phi for L1G trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiL1G   ->SetXTitle("Energy (GeV)");
  
  fhClusPhiL1J    = new TH2F("hClusPhiL1J","clusters distribution vs #phi for L1J trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiL1J   ->SetXTitle("Energy (GeV)");
  fhClusPhiL1J   ->SetYTitle("#phi (rad)");
  
  fhClusPhiL1GOnly = new TH2F("hClusPhiL1GOnly","clusters distribution vs #phi for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiL1GOnly->SetXTitle("Energy (GeV)");
  fhClusPhiL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusPhiL1JOnly = new TH2F("hClusPhiL1JOnly","clusters distribution vs #phi for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiL1JOnly->SetXTitle("Energy (GeV)");
  fhClusPhiL1JOnly->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxMB  = new TH2F("hClusPhiMaxMB","maximum energy cluster per event vs #phi for MB trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxMB ->SetXTitle("Energy (GeV)");
  fhClusPhiMaxMB ->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxL0  = new TH2F("hClusPhiMaxL0","maximum energy cluster per event vs #phi for L0 trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxL0 ->SetXTitle("Energy (GeV)");
  fhClusPhiMaxL0 ->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxL1G = new TH2F("hClusPhiMaxL1G","maximum energy cluster per event vs #phi for L1G trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxL1G->SetXTitle("Energy (GeV)");
  fhClusPhiMaxL1G->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxL1J = new TH2F("hClusPhiMaxL1J","maximum energy cluster per event vs #phi for L1J trigger",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxL1J->SetXTitle("Energy (GeV)");
  fhClusPhiMaxL1J->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxL1GOnly = new TH2F("hClusPhiMaxL1GOnly","maximum energy cluster per event vs #phi for L1G trigger and not L1J",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxL1GOnly->SetXTitle("Energy (GeV)");
  fhClusPhiMaxL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusPhiMaxL1JOnly = new TH2F("hClusPhiMaxL1JOnly","maximum energy cluster per event vs #phi for L1J trigger and not L1G",fNBinsClusterE,0,fMaxClusterE,nphi, phimin, phimax);
  fhClusPhiMaxL1JOnly->SetXTitle("Energy (GeV)");
  fhClusPhiMaxL1JOnly->SetYTitle("#phi (rad)");
  
  
  fOutputList->Add(fhClusPhiMB);
  fOutputList->Add(fhClusPhiL0);
  fOutputList->Add(fhClusPhiL1G);
  fOutputList->Add(fhClusPhiL1J);
  fOutputList->Add(fhClusPhiL1GOnly);
  fOutputList->Add(fhClusPhiL1JOnly);
  
  fOutputList->Add(fhClusPhiMaxMB);
  fOutputList->Add(fhClusPhiMaxL0);
  fOutputList->Add(fhClusPhiMaxL1G);
  fOutputList->Add(fhClusPhiMaxL1J);
  fOutputList->Add(fhClusPhiMaxL1GOnly);
  fOutputList->Add(fhClusPhiMaxL1JOnly);
  
  // Cluster histograms, Pseudorapidity vs Azimuthal angle
  
  fhClusEtaPhiHighMB     = new TH2F("hClusEtaPhiHighMB","clusters distribution #eta vs #phi for MB trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighMB    ->SetXTitle("#eta");
  fhClusEtaPhiHighMB    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighL0     = new TH2F("hClusEtaPhiHighL0","clusters distribution #eta  vs #phi for L0 trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighL0    ->SetXTitle("#eta");
  fhClusEtaPhiHighL0    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighL1G    = new TH2F("hClusEtaPhiHighL1G","clusters distribution #eta  vs #phi for L1G trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighL1G   ->SetXTitle("#eta");
  
  fhClusEtaPhiHighL1J    = new TH2F("hClusEtaPhiHighL1J","clusters distribution #eta  vs #phi for L1J trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighL1J   ->SetXTitle("#eta");
  fhClusEtaPhiHighL1J   ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighL1GOnly = new TH2F("hClusEtaPhiHighL1GOnly","clusters distribution #eta  vs #phi for L1G trigger and not L1J, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighL1GOnly->SetXTitle("#eta");
  fhClusEtaPhiHighL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighL1JOnly = new TH2F("hClusEtaPhiHighL1JOnly","clusters distribution #eta  vs #phi for L1J trigger and not L1G, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighL1JOnly->SetXTitle("#eta");
  fhClusEtaPhiHighL1JOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxMB  = new TH2F("hClusEtaPhiHighCluMaxMB","maximum energy cluster per event #eta  vs #phi for MB trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxMB ->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxMB ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxL0  = new TH2F("hClusEtaPhiHighCluMaxL0","maximum energy cluster per event #eta  vs #phi for L0 trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxL0 ->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxL0 ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxL1G = new TH2F("hClusEtaPhiHighCluMaxL1G","maximum energy cluster per event #eta  vs #phi for L1G trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxL1G->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxL1G->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxL1J = new TH2F("hClusEtaPhiHighCluMaxL1J","maximum energy cluster per event #eta  vs #phi for L1J trigger, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxL1J->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxL1J->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxL1GOnly = new TH2F("hClusEtaPhiHighCluMaxL1GOnly","maximum energy cluster per event #eta  vs #phi for L1G trigger and not L1J, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxL1GOnly->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiHighCluMaxL1JOnly = new TH2F("hClusEtaPhiHighCluMaxL1JOnly","maximum energy cluster per event #eta  vs #phi for L1J trigger and not L1G, E > 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiHighCluMaxL1JOnly->SetXTitle("#eta");
  fhClusEtaPhiHighCluMaxL1JOnly->SetYTitle("#phi (rad)");
  
  
  fOutputList->Add(fhClusEtaPhiHighMB);
  fOutputList->Add(fhClusEtaPhiHighL0);
  fOutputList->Add(fhClusEtaPhiHighL1G);
  fOutputList->Add(fhClusEtaPhiHighL1J);
  fOutputList->Add(fhClusEtaPhiHighL1GOnly);
  fOutputList->Add(fhClusEtaPhiHighL1JOnly);
  
  fOutputList->Add(fhClusEtaPhiHighCluMaxMB);
  fOutputList->Add(fhClusEtaPhiHighCluMaxL0);
  fOutputList->Add(fhClusEtaPhiHighCluMaxL1G);
  fOutputList->Add(fhClusEtaPhiHighCluMaxL1J);
  fOutputList->Add(fhClusEtaPhiHighCluMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiHighCluMaxL1JOnly);
  
  fhClusEtaPhiLowMB     = new TH2F("hClusEtaPhiLowMB","clusters distribution #eta vs #phi for MB trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowMB    ->SetXTitle("#eta");
  fhClusEtaPhiLowMB    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowL0     = new TH2F("hClusEtaPhiLowL0","clusters distribution #eta  vs #phi for L0 trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowL0    ->SetXTitle("#eta");
  fhClusEtaPhiLowL0    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowL1G    = new TH2F("hClusEtaPhiLowL1G","clusters distribution #eta  vs #phi for L1G trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowL1G   ->SetXTitle("#eta");
  
  fhClusEtaPhiLowL1J    = new TH2F("hClusEtaPhiLowL1J","clusters distribution #eta  vs #phi for L1J trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowL1J   ->SetXTitle("#eta");
  fhClusEtaPhiLowL1J   ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowL1GOnly = new TH2F("hClusEtaPhiLowL1GOnly","clusters distribution #eta  vs #phi for L1G trigger and not L1J, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowL1GOnly->SetXTitle("#eta");
  fhClusEtaPhiLowL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowL1JOnly = new TH2F("hClusEtaPhiLowL1JOnly","clusters distribution #eta  vs #phi for L1J trigger and not L1G, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowL1JOnly->SetXTitle("#eta");
  fhClusEtaPhiLowL1JOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxMB  = new TH2F("hClusEtaPhiLowCluMaxMB","maximum energy cluster per event #eta  vs #phi for MB trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxMB ->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxMB ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxL0  = new TH2F("hClusEtaPhiLowCluMaxL0","maximum energy cluster per event #eta  vs #phi for L0 trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxL0 ->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxL0 ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxL1G = new TH2F("hClusEtaPhiLowCluMaxL1G","maximum energy cluster per event #eta  vs #phi for L1G trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxL1G->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxL1G->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxL1J = new TH2F("hClusEtaPhiLowCluMaxL1J","maximum energy cluster per event #eta  vs #phi for L1J trigger, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxL1J->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxL1J->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxL1GOnly = new TH2F("hClusEtaPhiLowCluMaxL1GOnly","maximum energy cluster per event #eta  vs #phi for L1G trigger and not L1J, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxL1GOnly->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCluMaxL1JOnly = new TH2F("hClusEtaPhiLowCluMaxL1JOnly","maximum energy cluster per event #eta  vs #phi for L1J trigger and not L1G, E < 10 GeV",neta, etamin, etamax,nphi, phimin, phimax);
  fhClusEtaPhiLowCluMaxL1JOnly->SetXTitle("#eta");
  fhClusEtaPhiLowCluMaxL1JOnly->SetYTitle("#phi (rad)");
  
  
  fOutputList->Add(fhClusEtaPhiLowMB);
  fOutputList->Add(fhClusEtaPhiLowL0);
  fOutputList->Add(fhClusEtaPhiLowL1G);
  fOutputList->Add(fhClusEtaPhiLowL1J);
  fOutputList->Add(fhClusEtaPhiLowL1GOnly);
  fOutputList->Add(fhClusEtaPhiLowL1JOnly);
  
  fOutputList->Add(fhClusEtaPhiLowCluMaxMB);
  fOutputList->Add(fhClusEtaPhiLowCluMaxL0);
  fOutputList->Add(fhClusEtaPhiLowCluMaxL1G);
  fOutputList->Add(fhClusEtaPhiLowCluMaxL1J);
  fOutputList->Add(fhClusEtaPhiLowCluMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiLowCluMaxL1JOnly);
  
  
  fhClusEtaPhiHighCellMaxMB     = new TH2F("hClusEtaPhiHighCellMaxMB","Cluster hit map in calorimeter (max cell), column vs row for MB trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxMB    ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxMB    ->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxL0     = new TH2F("hClusEtaPhiHighCellMaxL0","Cluster hit map in calorimeter (max cell), column vs row for L0 trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxL0    ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxL0    ->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxL1G    = new TH2F("hClusEtaPhiHighCellMaxL1G","Cluster hit map in calorimeter (max cell), column vs row for L1G trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxL1G   ->SetXTitle("Index #eta (columnns)");
  
  fhClusEtaPhiHighCellMaxL1J    = new TH2F("hClusEtaPhiHighCellMaxL1J","Cluster hit map in calorimeter (max cell), column vs row for L1J trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxL1J   ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxL1J   ->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxL1GOnly = new TH2F("hClusEtaPhiHighCellMaxL1GOnly","Cluster hit map in calorimeter (max cell), column vs row for L1G trigger and not L1J, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxL1GOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxL1GOnly->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxL1JOnly = new TH2F("hClusEtaPhiHighCellMaxL1JOnly","Cluster hit map in calorimeter (max cell), column vs row for L1J trigger and not L1G, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxL1JOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxL1JOnly->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxMB  = new TH2F("hClusEtaPhiHighCellMaxCluMaxMB","Max E cluster hit map in calorimeter (max cell), column vs row  for MB trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxMB ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxMB ->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxL0  = new TH2F("hClusEtaPhiHighCellMaxCluMaxL0","Max E cluster hit map in calorimeter (max cell), column vs row  for L0 trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxL0 ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxL0 ->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxL1G = new TH2F("hClusEtaPhiHighCellMaxCluMaxL1G","Max E cluster hit map in calorimeter (max cell), column vs row  for L1G trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxL1G->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxL1G->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxL1J = new TH2F("hClusEtaPhiHighCellMaxCluMaxL1J","Max E cluster hit map in calorimeter (max cell), column vs row  for L1J trigger, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxL1J->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxL1J->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxL1GOnly = new TH2F("hClusEtaPhiHighCellMaxCluMaxL1GOnly","Max E cluster hit map in calorimeter (max cell), column vs row  for L1G trigger and not L1J, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxL1GOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxL1GOnly->SetYTitle("Index #phi (rows)");
  
  fhClusEtaPhiHighCellMaxCluMaxL1JOnly = new TH2F("hClusEtaPhiHighCellMaxCluMaxL1JOnly","Max E cluster hit map in calorimeter (max cell), column vs row  for L1J trigger and not L1G, E > 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiHighCellMaxCluMaxL1JOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiHighCellMaxCluMaxL1JOnly->SetYTitle("Index #phi (rows)");
  
  
  fOutputList->Add(fhClusEtaPhiHighCellMaxMB);
  fOutputList->Add(fhClusEtaPhiHighCellMaxL0);
  fOutputList->Add(fhClusEtaPhiHighCellMaxL1G);
  fOutputList->Add(fhClusEtaPhiHighCellMaxL1J);
  fOutputList->Add(fhClusEtaPhiHighCellMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiHighCellMaxL1JOnly);
  
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxMB);
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxL0);
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxL1G);
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxL1J);
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiHighCellMaxCluMaxL1JOnly);
  
  fhClusEtaPhiLowCellMaxMB     = new TH2F("hClusEtaPhiLowCellMaxMB","Cluster hit map in calorimeter (max cell), column vs row for MB trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxMB    ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxMB    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxL0     = new TH2F("hClusEtaPhiLowCellMaxL0","Cluster hit map in calorimeter (max cell), column vs row for L0 trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxL0    ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxL0    ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxL1G    = new TH2F("hClusEtaPhiLowCellMaxL1G","Cluster hit map in calorimeter (max cell), column vs row for L1G trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxL1G   ->SetXTitle("Index #eta (columnns)");
  
  fhClusEtaPhiLowCellMaxL1J    = new TH2F("hClusEtaPhiLowCellMaxL1J","Cluster hit map in calorimeter (max cell), column vs row for L1J trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxL1J   ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxL1J   ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxL1GOnly = new TH2F("hClusEtaPhiLowCellMaxL1GOnly","Cluster hit map in calorimeter (max cell), column vs row for L1G trigger and not L1J, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxL1GOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxL1JOnly = new TH2F("hClusEtaPhiLowCellMaxL1JOnly","Cluster hit map in calorimeter (max cell), column vs row for L1J trigger and not L1G, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxL1JOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxL1JOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxMB  = new TH2F("hClusEtaPhiLowCellMaxCluMaxMB","Max E cluster hit map in calorimeter (max cell), column vs row  for MB trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxMB ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxMB ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxL0  = new TH2F("hClusEtaPhiLowCellMaxCluMaxL0","Max E cluster hit map in calorimeter (max cell), column vs row  for L0 trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxL0 ->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxL0 ->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxL1G = new TH2F("hClusEtaPhiLowCellMaxCluMaxL1G","Max E cluster hit map in calorimeter (max cell), column vs row  for L1G trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxL1G->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxL1G->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxL1J = new TH2F("hClusEtaPhiLowCellMaxCluMaxL1J","Max E cluster hit map in calorimeter (max cell), column vs row  for L1J trigger, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxL1J->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxL1J->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxL1GOnly = new TH2F("hClusEtaPhiLowCellMaxCluMaxL1GOnly","Max E cluster hit map in calorimeter (max cell), column vs row  for L1G trigger and not L1J, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxL1GOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxL1GOnly->SetYTitle("#phi (rad)");
  
  fhClusEtaPhiLowCellMaxCluMaxL1JOnly = new TH2F("hClusEtaPhiLowCellMaxCluMaxL1JOnly","Max E cluster hit map in calorimeter (max cell), column vs row  for L1J trigger and not L1G, E < 10 GeV",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhClusEtaPhiLowCellMaxCluMaxL1JOnly->SetXTitle("Index #eta (columnns)");
  fhClusEtaPhiLowCellMaxCluMaxL1JOnly->SetYTitle("#phi (rad)");
  
  
  fOutputList->Add(fhClusEtaPhiLowCellMaxMB);
  fOutputList->Add(fhClusEtaPhiLowCellMaxL0);
  fOutputList->Add(fhClusEtaPhiLowCellMaxL1G);
  fOutputList->Add(fhClusEtaPhiLowCellMaxL1J);
  fOutputList->Add(fhClusEtaPhiLowCellMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiLowCellMaxL1JOnly);
  
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxMB);
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxL0);
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxL1G);
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxL1J);
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxL1GOnly);
  fOutputList->Add(fhClusEtaPhiLowCellMaxCluMaxL1JOnly);  
  
  PostData(1, fOutputList);  
  
}
//______________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserExec(Option_t *) 
{
  // Main loop
  
  AliVEvent* event = InputEvent();
  
  //Remove next lines when AODs ready
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
  
  if (!esdEvent) 
  {
    AliError("Work only with ESDs, not available, exit");
    return;
  }
  
  //trigger configuration
  TString triggerclasses = esdEvent->GetFiredTriggerClasses();
  
  Int_t eventType = ((AliVHeader*)InputEvent()->GetHeader())->GetEventType();
  // physics events eventType=7, select only those
  
  if(triggerclasses=="" || eventType != 7) return;
  
  // Check trigger
  Bool_t bMB  = kFALSE;
  Bool_t bL0  = kFALSE;
  Bool_t bL1G = kFALSE;
  Bool_t bL1J = kFALSE;
  
  if(triggerclasses.Contains("CINT7-B-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CINT7-I-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CINT1-I-NOPF-ALLNOTRD") || 
     triggerclasses.Contains("CINT1-B-NOPF-ALLNOTRD") ||
     triggerclasses.Contains("CPBI2_B1-B-NOPF-ALLNOTRD") )   bMB  = kTRUE;
  
  if(triggerclasses.Contains("CEMC7-B-NOPF-ALLNOTRD") || 
     triggerclasses.Contains("CEMC1-B-NOPF-ALLNOTRD")    )   bL0  = kTRUE;
  
  if(triggerclasses.Contains("CEMC7EGA-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CPBI2EGA")                 )   bL1G = kTRUE;
  
  if(triggerclasses.Contains("CEMC7EJE-B-NOPF-CENTNOTRD") ||
     triggerclasses.Contains("CPBI2EJE")                 )   bL1J = kTRUE;
  
  // Fill event histo
  fhNEvents->Fill(0.5); // All physics events
  
  if( bMB ) 
  { 
    fhNEvents->Fill(1.5);
    if( !bL1G && !bL1J && !bL0 ) fhNEvents->Fill(8.5);
  }
  else 
  {
    if( bL0  ) fhNEvents->Fill( 9.5);
    if( bL1G ) fhNEvents->Fill(10.5);
    if( bL1J ) fhNEvents->Fill(11.5);
  }

  if( bL0  )  fhNEvents->Fill(2.5);
  
  if( bL1G ) 
  {
    fhNEvents->Fill(3.5);
    if(!bL1J)  fhNEvents->Fill(5.5);
  }
  
  if( bL1J ) 
  {
    fhNEvents->Fill(4.5);
    if(!bL1G)  fhNEvents->Fill(6.5);
  }
  
  if(bL1J && bL1G) fhNEvents->Fill(7.5);

    
  //std::cout << "trigger = " << triggerclasses << std::endl;
  
  //map for cells and patches
  
  Double_t emcalCell     [fgkFALTRORows][fgkFALTROCols], emcalCellL1G  [fgkFALTRORows][fgkFALTROCols];
  Double_t emcalCellL1J  [fgkFALTRORows][fgkFALTROCols], emcalTrigL0   [fgkFALTRORows][fgkFALTROCols];  
  Double_t emcalTrigL0L1G[fgkFALTRORows][fgkFALTROCols], emcalTrigL0L1J[fgkFALTRORows][fgkFALTROCols]; 
  Double_t emcalTrigL1G  [fgkFALTRORows][fgkFALTROCols], emcalTrigL1J  [fgkFALTRORows][fgkFALTROCols], emcalTrigL1  [fgkFALTRORows][fgkFALTROCols];
  Double_t emcalPatchL0  [fgkFALTRORows][fgkFALTROCols], emcalPatchL1G [fgkFALTRORows][fgkFALTROCols], emcalPatchL1J[fgkFALTRORows][fgkFALTROCols];
  
  for (Int_t i = 0; i < fgkFALTRORows; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++) 
    {   
      emcalTrigL0[i][j]   = 0.;
      emcalTrigL0L1G[i][j]= 0.;
      emcalTrigL0L1J[i][j]= 0.;
      emcalTrigL1G[i][j]  = 0.;
      emcalTrigL1J[i][j]  = 0.;
      emcalTrigL1[i][j]   = 0.;
      emcalCell[i][j]     = 0.;
      emcalCellL1G[i][j]  = 0.;
      emcalCellL1J[i][j]  = 0.;
      emcalPatchL0[i][j]  = 0.;
      emcalPatchL1G[i][j] = 0.;
      emcalPatchL1J[i][j] = 0.;
    }
  }
  
  // ---------------------------------
  // Cells analysis
  // Fill FEE energy per channel array
  // ---------------------------------
  
  Int_t posX    = -1, posY = -1;
  Int_t nSupMod = -1, ieta = -1, iphi = -1, nModule = -1, nIphi = -1, nIeta = -1;
  Short_t absId = -1;
  Int_t nCells  =  0;
  
  AliVCaloCells& cells= *(event->GetEMCALCells());
  
  if (cells.IsEMCAL()) 
  {
    for (Int_t icell = 0; icell <  cells.GetNumberOfCells(); icell++) 
    {
      nCells ++;
      
      Double_t amp =0., time = 0., efrac = 0;
      Short_t mclabel = -1;
      
      cells.GetCell(icell, absId, amp, time,mclabel,efrac);	
      
      fGeometry->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
      fGeometry->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta); 
      
      posX = (nSupMod % 2) ? ieta + AliEMCALGeoParams::fgkEMCALCols : ieta;				
      posY = iphi + AliEMCALGeoParams::fgkEMCALRows * int(nSupMod / 2);
      
      if(int(posX/2) > fgkFALTROCols || int(posY/2) > fgkFALTRORows ) {
        if(DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Wrong Position (x,y) = (%d,%d)\n",posX,posY);
        continue;
      }
      
      emcalCell[int(posY/2)][int(posX/2)] += amp; 
      
      if(triggerclasses.Contains("CEMC7EGA-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EGA")) emcalCellL1G[int(posY/2)][int(posX/2)] += amp;
      if(triggerclasses.Contains("CEMC7EJE-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EJE")) emcalCellL1J[int(posY/2)][int(posX/2)] += amp;
      
    }
  }
  
  //-------------------------------------  
  // Trigger analysis, fill L0, L1 arrays
  //------------------------------------- 
  
  AliESDCaloTrigger& trg= * (esdEvent->GetCaloTrigger("EMCAL"));
  
  Int_t    nL0Patch = 0 ;
  Int_t    nL1Patch = 0 ;
  Double_t totSTU   = 0.;
  Double_t totTRU   = 0.;
  
  trg.Reset();
  while (trg.Next())
  {
    trg.GetPosition(posX,posY);
    
    
    if (posX > -1 && posY > -1) 
    {
      //L0 analysis  
      Int_t nTimes = 0;
      trg.GetNL0Times(nTimes);
      
      Float_t ampL0 = 0.;
      trg.GetAmplitude(ampL0);
      if (ampL0 > 0) emcalTrigL0[posY][posX] = ampL0;
      if(triggerclasses.Contains("CEMC7EGA-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EGA")) emcalTrigL0L1G[posY][posX] += ampL0;
      if(triggerclasses.Contains("CEMC7EJE-B-NOPF-CENTNOTRD") || triggerclasses.Contains("CPBI2EJE")) emcalTrigL0L1J[posY][posX] += ampL0;
      totTRU += ampL0;
      
      if (nTimes) 
	    {
	      nL0Patch += nTimes;
	      emcalPatchL0[posY][posX] = 1.;
	      fhL0Patch->Fill(posX,posY);
	    }
      
      //L1 analysis
      Int_t bit = 0;
      trg.GetTriggerBits(bit);
      
      Int_t ts = 0;
      trg.GetL1TimeSum(ts);
      if (ts > 0) emcalTrigL1[posY][posX] = ts;
      totSTU += ts;
      
      //L1-Gamma
      if (bit >> 6 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1G[posY][posX] += 1.;
        fhL1GPatch->Fill(posX,posY);
        
        if (ts > 0) emcalTrigL1G[posY][posX] = ts;
        
        //printf("Gamma STU patch %d, time sum %d, posX %d , posY %d\n",nL1Patch,ts,posX, posY);
      }
      
      //L1-Jet
      if (bit >> 8 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1J[posY][posX] += 1.;
        fhL1JPatch->Fill(posX,posY);
        
        if (ts > 0) emcalTrigL1J[posY][posX] = ts;
        
        //printf("Jet STU patch %d, time sum %d, posX %d , posY %d\n",nL1Patch,ts,posX, posY);
        
      }
      
    }
  }
  
  if(totTRU > fMaxTRUSignal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large totTRU %f\n",totTRU);
  if(totSTU > fMaxSTUSignal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large totSTU %f\n",totSTU);
 
  //V0 analysis 
  AliESDVZERO* eventV0 = esdEvent->GetVZEROData(); 
	
  Float_t v0C = 0, v0A = 0, v0TT = trg.GetL1V0(0)+trg.GetL1V0(1);
	
  if (eventV0) 
  {
    for (Int_t i = 0; i < 32; i++)
    {
      v0C += eventV0->GetAdcV0C(i);
      v0A += eventV0->GetAdcV0A(i);
    }
  }
  
  if (totSTU != 0) {
    fhV0STU->Fill(v0A+v0C,totSTU);
    if( v0A+v0C > fMaxV0Signal && DebugLevel() > 0) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large v0A+v0C %f\n",v0A+v0C);
  }
  
  if( bL1G ) fhV0L1G->Fill(v0A+v0C);                              ;
  if( bL1J ) fhV0L1J->Fill(v0A+v0C);
  if( bMB  ) fhV0MB ->Fill(v0A+v0C);
  
  //if(nL0Patch!=0 || nL1Patch!=0) printf("total TRU %f, total STU %f, V0C+V0A %f; nL0 %d, nL1 %d \n",
  //       totTRU,totSTU,v0A+v0C,nL0Patch,nL1Patch);
  
  //Maximum amplitude patch analysis
  Int_t patchMax = 0;
  Int_t colMax = -1;
  Int_t rowMax = -1;
  Int_t col = 0;
  Int_t row = 0;
  
  for (Int_t i = 0; i < 47; i++)
  {
    for (Int_t j = 0; j < 59; j++)
    {				
      Int_t patchG = 0;
      col = i;
      row = j;
      
      for (Int_t k = 0; k < 2; k++) 
	    {
	      for (Int_t l = 0; l < 2; l++) 
        {
          patchG += int(emcalTrigL1[j + l][i + k]);
        }
	    }
      
      if (patchG > patchMax) 
      {
        patchMax = patchG;
        colMax = col;
        rowMax = row;
      }
    }
  }
  
  fhGPMaxVV0TT->Fill(v0TT, patchMax);
  if( bL1G ) fhL1GPatchMax->Fill(colMax,rowMax);
  
  patchMax = 0;
  colMax = -1;
  rowMax = -1;
  
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
          patchJ += int(emcalTrigL1[4*j + l][4*i + k]);
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
  
  fhJPMaxVV0TT->Fill(v0TT, patchMax);
  if( bL1J ) fhL1JPatchMax->Fill(colMax,rowMax);
  
  //Matrix with signal per channel
  for (Int_t i = 0; i < fgkFALTRORows; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++) //check x,y direction for reading FOR ((0,0) = top left);
    {
      fhFORAmp    ->Fill( j, i, emcalCell     [i][j]);
      fhFORAmpL1G ->Fill( j, i, emcalCellL1G  [i][j]);
      fhFORAmpL1J ->Fill( j, i, emcalCellL1J  [i][j]);
      fhL0Amp     ->Fill( j, i, emcalTrigL0   [i][j]);
      fhL0AmpL1G  ->Fill( j, i, emcalTrigL0L1G[i][j]);
      fhL0AmpL1J  ->Fill( j, i, emcalTrigL0L1J[i][j]);
      fhL1Amp     ->Fill( j, i, emcalTrigL1   [i][j]);
      fhL1GAmp    ->Fill( j, i, emcalTrigL1G  [i][j]);
      fhL1JAmp    ->Fill( j, i, emcalTrigL1J  [i][j]);
      fhFORMeanAmp->Fill( j, i, emcalCell     [i][j]);
      fhL0MeanAmp ->Fill( j, i, emcalTrigL0   [i][j]);
      fhL1MeanAmp ->Fill( j, i, emcalTrigL1   [i][j]);
    }
  }
  
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
      //       i,j,iTRU,iADC,idFOR, emcalCell  [i][j],emcalTrigL0[i][j],emcalTrigL1[i][j]);
      
      if (iTRU >= 0)
      {
        ampFOR[iTRU] += emcalCell  [i][j];
        ampL0[iTRU]  += emcalTrigL0[i][j];
        ampL1[iTRU]  += emcalTrigL1[i][j];
      }
    }
  }
  
  // FEE vs STU and TRU vs STU ratios  
  for (Int_t i = 0; i < 30; i++)
  {
    
    if (ampFOR[i] != 0 && ampL1[i] != 0) { 
      fhFEESTU->Fill(ampL1[i]/ampFOR[i],i);
      if(ampL1[i]/ampFOR[i] > fMaxSTUFEERatio  && DebugLevel() > 0 ) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large STU/FEE ratio %f\n",ampL1[i]/ampFOR[i]);
    }
    
    if (ampL0[i]  != 0 && ampL1[i] != 0) {
      fhTRUSTU->Fill(ampL1[i]/ampL0[i] ,i);
      if(ampL1[i]/ampL0[i] > fMaxSTUTRURatio  && DebugLevel() > 0 ) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large STU/TRU ratio %f\n",ampL1[i]/ampL0[i]);
    }
    
  }
  
  //Get Vertex
  Double_t v[3] = {0,0,0};
  esdEvent->GetVertex()->GetXYZ(v);
  
  //clusters distribution
  TRefArray* caloClus = new TRefArray();
  esdEvent->GetEMCALClusters(caloClus);
  
  Int_t nCaloClusters = caloClus->GetEntriesFast();
  
  Float_t emax   = 0;
  Float_t etamax = 0;
  Float_t phimax = 0;
  Float_t ietamax=-1;
  Float_t iphimax=-1;

  Float_t e      = 0;
  Float_t eta    = 0;
  Float_t phi    = 0;

  //Energy threshold to fill Eta vs Phi histograms
  Float_t etaphiEnMin = 10.;
  
  TLorentzVector mom;
  
  //Get vertex for momentum calculation
  Double_t vertex[] = {0.0,0.0,0.0};
  //InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  Float_t centrality = -1;
  if(InputEvent()->GetCentrality()) centrality = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
    
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
    AliESDCaloCluster *clus = (AliESDCaloCluster*) (caloClus->At(icalo));
		
    if(!clus->IsEMCAL()) continue;
    
    if(!fRecoUtils->IsGoodCluster(clus,fGeometry,InputEvent()->GetEMCALCells(),InputEvent()->GetBunchCrossNumber()))
    { 
      continue;
    }
    
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

    if(ieta > fgkFALTROCols || iphi > fgkFALTRORows ) {
      printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Wrong Position (x,y) = (%d,%d)\n",ieta,iphi);
    }
           
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
    
    if( bMB  ) 
    { 
      fhClusMB      ->Fill(e); 
      fhClusCenMB   ->Fill(e,centrality); 
      fhClusV0MB    ->Fill(e,v0A+v0C);
      fhClusEtaMB   ->Fill(e,eta);
      fhClusPhiMB   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighMB        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxMB ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowMB         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxMB  ->Fill(ieta,iphi);
      }
    }
    
    if( bL0  )    { 
      fhClusL0      ->Fill(e); 
      fhClusCenL0   ->Fill(e,centrality); 
      fhClusV0L0    ->Fill(e,v0A+v0C);
      fhClusEtaL0   ->Fill(e,eta);
      fhClusPhiL0   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighL0        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxL0 ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowL0         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxL0  ->Fill(ieta,iphi);
      }
    }
    
    if( bL1G )    { 
      fhClusL1G      ->Fill(e); 
      fhClusCenL1G   ->Fill(e,centrality); 
      fhClusV0L1G    ->Fill(e,v0A+v0C);
      fhClusEtaL1G   ->Fill(e,eta);
      fhClusPhiL1G   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighL1G        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxL1G ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowL1G         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxL1G  ->Fill(ieta,iphi);
      }
    }
    
    if( bL1J )    { 
      fhClusL1J      ->Fill(e); 
      fhClusCenL1J   ->Fill(e,centrality); 
      fhClusV0L1J    ->Fill(e,v0A+v0C);
      fhClusEtaL1J   ->Fill(e,eta);
      fhClusPhiL1J   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighL1J        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxL1J ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowL1J         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxL1J  ->Fill(ieta,iphi);
      }
    }   
    
    if( bL1G && !bL1J )    
    { 
      fhClusL1GOnly      ->Fill(e); 
      fhClusCenL1GOnly   ->Fill(e,centrality); 
      fhClusV0L1GOnly    ->Fill(e,v0A+v0C);
      fhClusEtaL1GOnly   ->Fill(e,eta);
      fhClusPhiL1GOnly   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighL1GOnly        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxL1GOnly ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowL1GOnly         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxL1GOnly  ->Fill(ieta,iphi);
      }
    }    
    
    if( bL1J && !bL1G ) 
    { 
      fhClusL1JOnly      ->Fill(e); 
      fhClusCenL1JOnly   ->Fill(e,centrality); 
      fhClusV0L1JOnly    ->Fill(e,v0A+v0C);
      fhClusEtaL1JOnly   ->Fill(e,eta);
      fhClusPhiL1JOnly   ->Fill(e,phi);
      if(e > etaphiEnMin) 
      {
        fhClusEtaPhiHighL1JOnly        ->Fill( eta, phi);
        fhClusEtaPhiHighCellMaxL1JOnly ->Fill(ieta,iphi);
      }
      else {
        fhClusEtaPhiLowL1JOnly         ->Fill( eta, phi);
        fhClusEtaPhiLowCellMaxL1JOnly  ->Fill(ieta,iphi);
      }
    } 
    
    if( bMB && !bL1G && !bL1J && !bL0  ) fhClusMBPure  ->Fill(e);

  }
  
  // Maximum energy cluster per event histograms

  if( bMB  ) 
  { 
    fhClusMaxMB      ->Fill(emax); 
    fhClusCenMaxMB   ->Fill(emax,centrality); 
    fhClusV0MaxMB    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxMB   ->Fill(emax,etamax);
    fhClusPhiMaxMB   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxMB        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxMB ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxMB         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxMB  ->Fill(ietamax,iphimax);
    }
  }
  
  if( bL0  )    { 
    fhClusMaxL0      ->Fill(emax); 
    fhClusCenMaxL0   ->Fill(emax,centrality); 
    fhClusV0MaxL0    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxL0   ->Fill(emax,etamax);
    fhClusPhiMaxL0   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxL0        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxL0 ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxL0         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxL0  ->Fill(ietamax,iphimax);
    }
  }
  
  if( bL1G )    { 
    fhClusMaxL1G      ->Fill(emax); 
    fhClusCenMaxL1G   ->Fill(emax,centrality); 
    fhClusV0MaxL1G    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxL1G   ->Fill(emax,etamax);
    fhClusPhiMaxL1G   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxL1G        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxL1G ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxL1G         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxL1G  ->Fill(ietamax,iphimax);
    }
  }
  
  if( bL1J )    { 
    fhClusMaxL1J      ->Fill(emax); 
    fhClusCenMaxL1J   ->Fill(emax,centrality); 
    fhClusV0MaxL1J    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxL1J   ->Fill(emax,etamax);
    fhClusPhiMaxL1J   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxL1J        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxL1J ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxL1J         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxL1J  ->Fill(ietamax,iphimax);
    }
  }   
  
  if( bL1G && !bL1J )    
  { 
    fhClusMaxL1GOnly      ->Fill(emax); 
    fhClusCenMaxL1GOnly   ->Fill(emax,centrality);
    fhClusV0MaxL1GOnly    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxL1GOnly   ->Fill(emax,etamax);
    fhClusPhiMaxL1GOnly   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxL1GOnly        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxL1GOnly ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxL1GOnly         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxL1GOnly  ->Fill(ietamax,iphimax);
    }
  }    
  
  if( bL1J && !bL1G ) 
  { 
    fhClusMaxL1JOnly      ->Fill(emax); 
    fhClusCenMaxL1JOnly   ->Fill(emax,centrality); 
    fhClusV0MaxL1JOnly    ->Fill(emax,v0A+v0C);
    fhClusEtaMaxL1JOnly   ->Fill(emax,etamax);
    fhClusPhiMaxL1JOnly   ->Fill(emax,phimax);
    if(emax > etaphiEnMin) 
    {
      fhClusEtaPhiHighCluMaxL1JOnly        ->Fill( etamax, phimax);
      fhClusEtaPhiHighCellMaxCluMaxL1JOnly ->Fill(ietamax,iphimax);
    }
    else {
      fhClusEtaPhiLowCluMaxL1JOnly         ->Fill( etamax, phimax);
      fhClusEtaPhiLowCellMaxCluMaxL1JOnly  ->Fill(ietamax,iphimax);
    }
  } 
    
  if( bMB && !bL1G && !bL1J && !bL0 ) fhClusMaxMBPure  ->Fill(emax);
  
  PostData(1, fOutputList);  
  
}
