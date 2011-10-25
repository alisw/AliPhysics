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
//  Author: Nicolas Arbor (LPSC-Grenoble)                                  //
//          Gustavo Conesa Balbastre  (LPSC-Grenoble)                     //
//                                                                        //
//------------------------------------------------------------------------//


#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALGeometry.h"

#include "AliAnalysisTaskEMCALTriggerQA.h"

ClassImp(AliAnalysisTaskEMCALTriggerQA)

//______________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA() : 
AliAnalysisTaskSE(), 
fOutputList(0),
fGeometry(0),  fGeoName("EMCAL_COMPLETEV1"),
fhNEvents(0),
fhFORAmp(0),
fhFORAmpL1G(0),
fhFORAmpL1J(0),
fhL0Amp(0),
fhL0AmpL1G(0),
fhL0AmpL1J(0),
fhL1Amp(0),
fhL1GAmp(0),
fhL1JAmp(0),
fhL0Patch(0),
fhL1GPatch(0),
fhL1JPatch(0),
fhFEESTU(0),
fhTRUSTU(0),
fhV0STU(0),
fhFullTRUSTU(0),
fhSTUChecks(0),
fNBinsSTUSignal  (2000), fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000), fMaxTRUSignal  (200000),
fNBinsV0Signal   (2000), fMaxV0Signal   (20000),
fNBinsSTUFEERatio(2000), fMaxSTUFEERatio(20000),
fNBinsSTUTRURatio(2000), fMaxSTUTRURatio(200)

{
  // Constructor
  
}		      

//______________________________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA(const char *name) : 
AliAnalysisTaskSE(name), 
fOutputList(0),
fGeometry(0), fGeoName("EMCAL_COMPLETEV1"),
fhNEvents(0),
fhFORAmp(0),
fhFORAmpL1G(0),
fhFORAmpL1J(0),
fhL0Amp(0),
fhL0AmpL1G(0),
fhL0AmpL1J(0),
fhL1Amp(0),
fhL1GAmp(0),
fhL1JAmp(0),
fhL0Patch(0),
fhL1GPatch(0),
fhL1JPatch(0),
fhFEESTU(0),
fhTRUSTU(0),
fhV0STU(0),
fhFullTRUSTU(0),
fhSTUChecks(0),
fNBinsSTUSignal  (2000), fMaxSTUSignal  (200000),
fNBinsTRUSignal  (2000), fMaxTRUSignal  (200000),
fNBinsV0Signal   (2000), fMaxV0Signal   (20000),
fNBinsSTUFEERatio(2000), fMaxSTUFEERatio(20000),
fNBinsSTUTRURatio(2000), fMaxSTUTRURatio(200)

{
  // Constructor
  
  DefineOutput(1, TList::Class());
  
}


//___________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserCreateOutputObjects() 
{
  // Init histograms and geometry 
  
  fGeometry = AliEMCALGeometry::GetInstance(fGeoName);
  
  fOutputList  = new TList;
  fOutputList ->SetOwner(kTRUE);
  
  fhNEvents    = new TH1F("hNEvents","Number of selected events",1,0,1);
  fhNEvents   ->SetYTitle("N events");
  
  fhFORAmp     = new TH2F("hFORAmp", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmp    ->SetXTitle("Index #eta (collumns)");
  fhFORAmp    ->SetYTitle("Index #phi (rows)");
  fhFORAmp    ->SetZTitle("Amplitude");

  fhFORAmpL1G  = new TH2F("hFORAmpL1G", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1G trigger condition",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1G ->SetXTitle("Index #eta (collumns)");
  fhFORAmpL1G ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1G ->SetZTitle("Amplitude");

  fhFORAmpL1J  = new TH2F("hFORAmpL1J", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1J trigger condition",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORAmpL1J ->SetXTitle("Index #eta (collumns)");
  fhFORAmpL1J ->SetYTitle("Index #phi (rows)");
  fhFORAmpL1J ->SetZTitle("Amplitude");
  
  fhL0Amp      = new TH2F("hL0Amp","FALTRO signal per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Amp     ->SetXTitle("Index #eta (collumns)");
  fhL0Amp     ->SetYTitle("Index #phi (rows)");
  fhL0Amp     ->SetZTitle("Amplitude");

  fhL0AmpL1G   = new TH2F("hL0AmpL1G","FALTRO signal per Row and Column, with L1G trigger condition",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0AmpL1G  ->SetXTitle("Index #eta (collumns)");
  fhL0AmpL1G  ->SetYTitle("Index #phi (rows)");
  fhL0AmpL1G  ->SetZTitle("Amplitude");

  fhL0AmpL1J   = new TH2F("hL0AmpL1J","FALTRO signal per Row and Column, with L1j trigger condition",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0AmpL1J  ->SetXTitle("Index #eta (collumns)");
  fhL0AmpL1J  ->SetYTitle("Index #phi (rows)");
  fhL0AmpL1J  ->SetZTitle("Amplitude");
  
  fhL1Amp      = new TH2F("hL1Amp","STU signal per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1Amp     ->SetXTitle("Index #eta (collumns)");
  fhL1Amp     ->SetYTitle("Index #phi (rows)");
  fhL1Amp     ->SetZTitle("Amplitude");
  
  fhL1GAmp     = new TH2F("hL1GAmp","STU signal per Row and Column for L1 Gamma",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GAmp    ->SetXTitle("Index #eta (collumns)");
  fhL1GAmp    ->SetYTitle("Index #phi (rows)");
  fhL1GAmp    ->SetZTitle("Amplitude");
  
  fhL1JAmp     = new TH2F("hL1JAmp","STU signal per Row and Column for L1 Jet",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1JAmp    ->SetXTitle("Index #eta (collumns)");
  fhL1JAmp    ->SetYTitle("Index #phi (rows)");
  fhL1JAmp    ->SetZTitle("Amplitude");
  
  fhL0Patch    = new TH2F("hL0Patch","FOR with associated L0 Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Patch   ->SetXTitle("Index #eta (collumns)");
  fhL0Patch   ->SetYTitle("Index #phi (rows)");
  fhL0Patch   ->SetZTitle("counts");
  
  fhL1GPatch   = new TH2F("hL1GPatch","FOR with associated L1 Gamma Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatch  ->SetXTitle("Index #eta (collumns)");
  fhL1GPatch  ->SetYTitle("Index #phi (rows)");
  fhL1GPatch  ->SetZTitle("counts");
  
  fhL1JPatch   = new TH2F("hL1JPatch","FOR with associated L1 Jet Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1JPatch  ->SetXTitle("Index #eta (collumns)");
  fhL1JPatch  ->SetYTitle("Index #phi (rows)");
  fhL1JPatch  ->SetZTitle("counts");
  
  fhFullTRUSTU = new TH2I("hFullTRUSTU","Total signal STU vs TRU",fNBinsTRUSignal,0,fMaxTRUSignal,fNBinsSTUSignal,0,fMaxSTUSignal);
  fhFullTRUSTU->SetXTitle("Total signal TRU");
  fhFullTRUSTU->SetYTitle("Total signal STU");
  fhFullTRUSTU->SetZTitle("counts");
  
  fhV0STU      = new TH2I("hV0STU","Total signal STU vs V0C+V0S",fNBinsV0Signal,0,fMaxV0Signal,fNBinsSTUSignal,0,fMaxSTUSignal);
  fhV0STU     ->SetXTitle("Signal V0C+V0A");
  fhV0STU     ->SetYTitle("Total signal STU");
  fhV0STU     ->SetZTitle("counts");
  
  fhSTUChecks  = new TH2I("hSTUChecks","Check FEE/STU link",2,0,2,15,0,15);
  fhSTUChecks ->SetXTitle("Index #eta");
  fhSTUChecks ->SetYTitle("Index #phi");
  
  fhFEESTU     = new TH2F("hFEESTU","STU / FEE vs channel", fNBinsSTUFEERatio,0,fMaxSTUFEERatio,30,0,30);
  fhFEESTU    ->SetXTitle("STU/FEE signal");
  fhFEESTU    ->SetYTitle("channel");
  fhFEESTU    ->SetZTitle("counts");
  
  fhTRUSTU     = new TH2F("hTRUSTU","STU / TRU vs channel", fNBinsSTUTRURatio,0,fMaxSTUTRURatio,30,0,30);
  fhTRUSTU    ->SetXTitle("STU/TRU signal");
  fhTRUSTU    ->SetYTitle("channel");
  fhTRUSTU    ->SetZTitle("counts");
  
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
  fOutputList->Add(fhFullTRUSTU);
  fOutputList->Add(fhSTUChecks);
  fOutputList->Add(fhFEESTU);
  fOutputList->Add(fhTRUSTU);
  
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
  
  fhNEvents->Fill(0);

  //trigger configuration
  TString triggerclasses = esdEvent->GetFiredTriggerClasses();
  std::cout << "trigger = " << triggerclasses << std::endl;
  
  //map for cells and patches
  
  Double_t emcalCell   [fgkFALTRORows][fgkFALTROCols], emcalCellL1G   [fgkFALTRORows][fgkFALTROCols], emcalCellL1J   [fgkFALTRORows][fgkFALTROCols], emcalTrigL0  [fgkFALTRORows][fgkFALTROCols],  emcalTrigL0L1G  [fgkFALTRORows][fgkFALTROCols],  emcalTrigL0L1J  [fgkFALTRORows][fgkFALTROCols], emcalTrigL1G  [fgkFALTRORows][fgkFALTROCols], emcalTrigL1J  [fgkFALTRORows][fgkFALTROCols], emcalTrigL1[fgkFALTRORows][fgkFALTROCols];
  Double_t emcalPatchL0[fgkFALTRORows][fgkFALTROCols], emcalPatchL1G[fgkFALTRORows][fgkFALTROCols], emcalPatchL1J[fgkFALTRORows][fgkFALTROCols];
  
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
      
      Double_t amp =0., time = 0.;
      
      cells.GetCell(icell, absId, amp, time);	
      
      fGeometry->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
      fGeometry->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta); 
      
      posX = (nSupMod % 2) ? ieta + AliEMCALGeoParams::fgkEMCALCols : ieta;				
      posY = iphi + AliEMCALGeoParams::fgkEMCALRows * int(nSupMod / 2);
      
      if(int(posX/2) > fgkFALTROCols || int(posY/2) > fgkFALTRORows ) {
        printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Wrong Position (x,y) = (%d,%d)\n",posX,posY);
        continue;
      }
      
      emcalCell[int(posY/2)][int(posX/2)] += amp; 

      if(triggerclasses.Contains("CEMCTEGA-B-NOPF-CENTNOTRD")) emcalCellL1G[int(posY/2)][int(posX/2)] += amp;
      if(triggerclasses.Contains("CEMCTEJE-B-NOPF-CENTNOTRD")) emcalCellL1J[int(posY/2)][int(posX/2)] += amp;
      
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
      emcalTrigL0[posY][posX] += ampL0;
      if(triggerclasses.Contains("CEMCTEGA-B-NOPF-CENTNOTRD")) emcalTrigL0L1G[posY][posX] += ampL0;
      if(triggerclasses.Contains("CEMCTEJE-B-NOPF-CENTNOTRD")) emcalTrigL0L1J[posY][posX] += ampL0;
      totTRU += ampL0;
      
      if (nTimes) 
	    {
	      nL0Patch += nTimes;
	      emcalPatchL0[posY][posX] = 1.;
	      fhL0Patch->Fill(posX,59-posY);//59 is due to FOR reference
	    }
      
      //L1 analysis
      Int_t bit = 0;
      trg.GetTriggerBits(bit);
      
      Int_t ts = 0;
      trg.GetL1TimeSum(ts);
      emcalTrigL1 [posY][posX] += ts;
      totSTU += ts;
      
      //L1-Gamma
      if (bit >> 4 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1G[posY][posX] = 1.;
        fhL1GPatch->Fill(posX,59-posY);
        
        emcalTrigL1G[posY][posX] += ts;
        
        //printf("Gamma STU patch %d, time sum %d, posX %d , posY %d\n",nL1Patch,ts,posX, posY);
      }
      
      //L1-Jet
      if (bit >> 5 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1J[posY][posX] = 1.;
        fhL1JPatch->Fill(posX,59-posY);
        
        emcalTrigL1J[posY][posX] += ts;
        
        //printf("Jet STU patch %d, time sum %d, posX %d , posY %d\n",nL1Patch,ts,posX, posY);
        
      }
      
    }
  }
  
  if(totTRU > fMaxTRUSignal)printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large totTRU %f\n",totTRU);
  if(totSTU > fMaxSTUSignal)printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large totSTU %f\n",totSTU);
  
  if (totTRU != 0) fhFullTRUSTU->Fill(totTRU,totSTU);
  
  //V0 analysis 
  AliESDVZERO* eventV0 = esdEvent->GetVZEROData(); 
	
  Float_t v0C = 0, v0A = 0;
	
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
    if( v0A+v0C > fMaxV0Signal) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large v0A+v0C %f\n",v0A+v0C);
  }
  
  //if(nL0Patch!=0 || nL1Patch!=0) printf("total TRU %f, total STU %f, V0C+V0A %f; nL0 %d, nL1 %d \n",
  //       totTRU,totSTU,v0A+v0C,nL0Patch,nL1Patch);
  
  //Matrix with signal per channel
  for (Int_t i = 0; i < fgkFALTRORows; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++) //check x,y direction for reading FOR ((0,0) = top left);
    {
      fhFORAmp->Fill( j, fgkFALTRORows-i-1, emcalCell   [i][j]);
      fhFORAmpL1G->Fill( j, fgkFALTRORows-i-1, emcalCellL1G   [i][j]);
      fhFORAmpL1J->Fill( j, fgkFALTRORows-i-1, emcalCellL1J   [i][j]);
      fhL0Amp ->Fill( j, fgkFALTRORows-i-1, emcalTrigL0 [i][j]);
      fhL0AmpL1G ->Fill( j, fgkFALTRORows-i-1, emcalTrigL0L1G [i][j]);
      fhL0AmpL1J ->Fill( j, fgkFALTRORows-i-1, emcalTrigL0L1J [i][j]);
      fhL1Amp ->Fill( j, fgkFALTRORows-i-1, emcalTrigL1 [i][j]);
      fhL1GAmp->Fill( j, fgkFALTRORows-i-1, emcalTrigL1G[i][j]);
      fhL1JAmp->Fill( j, fgkFALTRORows-i-1, emcalTrigL1J[i][j]);
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
      if(ampL1[i]/ampFOR[i] > fMaxSTUFEERatio) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large STU/FEE ratio %f\n",ampL1[i]/ampFOR[i]);
    }
    
    if (ampL0[i]  != 0 && ampL1[i] != 0) {
      fhTRUSTU->Fill(ampL1[i]/ampL0[i] ,i);
      if(ampL1[i]/ampL0[i] > fMaxSTUTRURatio) printf("AliAnalysisTaskEMCALTriggerQA::UserExec() - Large STU/TRU ratio %f\n",ampL1[i]/ampL0[i]);
    }
    
  }
  
  PostData(1, fOutputList);  
  
}

//_______________________________________________________
void AliAnalysisTaskEMCALTriggerQA::Terminate(Option_t *)
{
  // Terminate analysis
  // Do some plots
  
  Int_t checkSTU[30];
  for (Int_t i = 0; i < 30; i++)
  {
    checkSTU[i] = 1 ; // Init array
    TH1F* hTRUSTUChannel = (TH1F*) fhTRUSTU->ProjectionY(Form("hTRUSTUChannel%d",i),i,i);
    //printf("AliAnalysisTaskEMCALTriggerQA::Terminate() - Channel %d TRUSTU Entries %d, Integral(10,20) %f ",
    //      i, hTRUSTUChannel->GetEntries(), hTRUSTUChannel->Integral(10,20));
    Int_t binMin = hTRUSTUChannel->FindBin(10);
    Int_t binMax = hTRUSTUChannel->FindBin(20);
    if     (hTRUSTUChannel->GetEntries() > 0 && 
            hTRUSTUChannel->Integral(binMin,binMax)/hTRUSTUChannel->GetEntries() < 0.9) 
      checkSTU[i] = 0;
    else if(hTRUSTUChannel->GetEntries() <= 0)
      checkSTU[i] = 0;
    
    //printf(" - %d\n",checkSTU[i]);
    delete hTRUSTUChannel;
  }
  
  for (Int_t i = 0; i < 30; i++)
  {
    if (i<15) fhSTUChecks->Fill(0.,i,checkSTU[i]);
    else      fhSTUChecks->Fill(1.,i,checkSTU[i]);
  }
  
}
