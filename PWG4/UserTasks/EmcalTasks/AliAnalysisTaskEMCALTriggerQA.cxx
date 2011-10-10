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

//_________________________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA() : 
AliAnalysisTaskSE(), 
fOutputList(0),
fGeometry(0),  fGeoName("EMCAL_COMPLETEV1"),
fhNEvents(0),
fhFORPos(0),
fhL0Pos(0),
fhL1Pos(0),
fhL0Patch(0),
fhL1GPatch(0),
fhL1JPatch(0),
fhV0STU(0),
fhFullTRUSTU(0),
fhSTUChecks(0),
fNBinsSTUSignal(1000),fMaxSTUSignal(100000),
fNBinsTRUSignal(1000),fMaxTRUSignal(100000),
fNBinsV0Signal (1000),fMaxV0Signal (10000)

{
  // Constructor
  for (int i = 0; i < 30; i++)
  {
    fhFEESTU[i] = 0;
    fhTRUSTU[i] = 0;
  }
  
}		      

//________________________________________________________________________
AliAnalysisTaskEMCALTriggerQA::AliAnalysisTaskEMCALTriggerQA(const char *name) : 
AliAnalysisTaskSE(name), 
fOutputList(0),
fGeometry(0), fGeoName("EMCAL_COMPLETEV1"),
fhNEvents(0),
fhFORPos(0),
fhL0Pos(0),
fhL1Pos(0),
fhL0Patch(0),
fhL1GPatch(0),
fhL1JPatch(0),
fhV0STU(0),
fhFullTRUSTU(0),
fhSTUChecks(0),
fNBinsSTUSignal(1000),fMaxSTUSignal(100000),
fNBinsTRUSignal(1000),fMaxTRUSignal(100000),
fNBinsV0Signal (1000),fMaxV0Signal (10000)

{
  // Constructor
  for (int i = 0; i < 30; i++)
  {
    fhFEESTU[i] = 0;
    fhTRUSTU[i] = 0;
  }
  
  DefineOutput(1, TList::Class());
  
}


//________________________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserCreateOutputObjects() 
{
  // Init histograms and geometry 
  
  fGeometry = AliEMCALGeometry::GetInstance(fGeoName);
  
  fOutputList = new TList;
  fOutputList->SetOwner(kTRUE);
  
  fhNEvents = new TH1F("hNEvents","Number of selected events",1,0,1);
  fhNEvents->SetYTitle("N events");
  
  fhFORPos  = new TH2F("hFORPos", "FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhFORPos->SetXTitle("Index #eta (collumns)");
  fhFORPos->SetYTitle("Index #phi (rows)");
  
  fhL0Pos   = new TH2F("hL0Pos","FALTRO signal per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Pos->SetXTitle("Index #eta (collumns)");
  fhL0Pos->SetYTitle("Index #phi (rows)");
  
  fhL1Pos   = new TH2F("hL1Pos","STU signal per Row and Column",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1Pos->SetXTitle("Index #eta (collumns)");
  fhL1Pos->SetYTitle("Index #phi (rows)");
  
  fhL0Patch   = new TH2F("fhL0Patch","FOR with associated L0 Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL0Patch->SetXTitle("Index #eta (collumns)");
  fhL0Patch->SetYTitle("Index #phi (rows)");
  
  fhL1GPatch   = new TH2F("fhL1GPatch","FOR with associated L1 Gamma Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1GPatch->SetXTitle("Index #eta (collumns)");
  fhL1GPatch->SetYTitle("Index #phi (rows)");
  
  fhL1JPatch   = new TH2F("fhL1JPatch","FOR with associated L1 Jet Patch",fgkFALTROCols,0,fgkFALTROCols,fgkFALTRORows,0,fgkFALTRORows);
  fhL1JPatch->SetXTitle("Index #eta (collumns)");
  fhL1JPatch->SetYTitle("Index #phi (rows)");
  
  fhFullTRUSTU  = new TH2I("hFullTRUSTU","Total signal STU vs TRU",fNBinsTRUSignal,0,fMaxTRUSignal,fNBinsSTUSignal,0,fMaxSTUSignal);
  fhFullTRUSTU->SetXTitle("Total signal TRU");
  fhFullTRUSTU->SetYTitle("Total signal STU");
  
  fhV0STU   = new TH2I("hV0STU","Total signal STU vs V0C+V0S",fNBinsV0Signal,0,fMaxV0Signal,fNBinsSTUSignal,0,fMaxSTUSignal);
  fhV0STU->SetXTitle("Signal V0C+V0A");
  fhV0STU->SetYTitle("Total signal STU");
  
  fhSTUChecks = new TH2I("hSTUChecks","Check FEE/STU link",2,0,2,15,0,15);
  fhSTUChecks->SetXTitle("Index #eta");
  fhSTUChecks->SetYTitle("Index #phi");
  
  for (int i = 0; i < 30; i++)
  {
    fhFEESTU[i] = new TH1F(Form("hFEESTU%d",i),Form("STU / FEE, channel %d",i),1000,0,100);
    fhFEESTU[i]->SetXTitle("STU/FEE signal");
    
    fhTRUSTU[i] = new TH1F(Form("hTRUSTU%d",i),Form("STU / TRU, channel %d",i),1000,0,100);
    fhTRUSTU[i]->SetXTitle("STU/TRU signal");
  }
  
  fOutputList->Add(fhNEvents);
  fOutputList->Add(fhV0STU);
  fOutputList->Add(fhFORPos);
  fOutputList->Add(fhL0Pos);
  fOutputList->Add(fhL1Pos);
  fOutputList->Add(fhL0Patch);
  fOutputList->Add(fhL1GPatch);
  fOutputList->Add(fhL1JPatch);
  fOutputList->Add(fhFullTRUSTU);
  fOutputList->Add(fhSTUChecks);
  
  
  for (int i = 0; i < 30; i++)
  {
    fOutputList->Add(fhFEESTU[i]);
    fOutputList->Add(fhTRUSTU[i]);
  }
  
  PostData(1, fOutputList);  
  
}
//________________________________________________________________________
void AliAnalysisTaskEMCALTriggerQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  AliVEvent* event = InputEvent();
  
  //Remove next lines when AODs ready
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
  
  if (!esdEvent) 
  {
    AliError("Work only with ESDs, not available, exit");
    return;
  }
  
  fhNEvents->Fill(0);
  
  
  //map for cells and patches
  
  Double_t emcalCell   [fgkFALTRORows][fgkFALTROCols], emcalTrigL0 [fgkFALTRORows][fgkFALTROCols], emcalTrigL1[fgkFALTRORows][fgkFALTROCols];
  Double_t emcalPatchL0[fgkFALTRORows][fgkFALTROCols], emcalPatchL1G[fgkFALTRORows][fgkFALTROCols], emcalPatchL1J[fgkFALTRORows][fgkFALTROCols];
  
  for (Int_t i = 0; i < fgkFALTRORows; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols; j++) 
    {   
      emcalTrigL0[i][j]  = 0.;
      emcalTrigL1[i][j]  = 0.;
      emcalCell[i][j]    = 0.;
      emcalPatchL0[i][j] = 0.;
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
      if(int(posX/2) > fgkFALTROCols) printf("Wrong X, posX %d\n",posX);
      if(int(posY/2) > fgkFALTRORows) printf("Wrong Y, posY %d\n",posX);
      
      emcalCell[int(posY/2)][int(posX/2)] += amp;     
    }
  }
  
  // -----------------------------------  
  //Trigger analysis, fill L0, L1 arrays
  // ----------------------------------- 
  
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
      if (nTimes) 
	    {
	      nL0Patch += nTimes;
	      Float_t ampL0 = 0.;
	      trg.GetAmplitude(ampL0);
	      emcalTrigL0[posY][posX] = ampL0;
	      emcalPatchL0[posY][posX] = 1.;
	      totTRU += ampL0;
	      fhL0Patch->Fill(posX-1,posY-1);//-1 is used to compare in good way patch L0 and patch L1
	    }
      
      //L1 analysis
      Int_t bit = 0;
      trg.GetTriggerBits(bit);
      
      //L1-Gamma
      if (bit >> 4 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1G[posY][posX] = 1.;
        fhL1GPatch->Fill(posX,posY);
        
        Int_t ts = 0;
        trg.GetL1TimeSum(ts);
        emcalTrigL1[posY][posX] += ts;
        totSTU += ts;
      }
      
      //L1-Jet
      if (bit >> 5 & 0x1) 
      {
        nL1Patch ++;
        emcalPatchL1J[posY][posX] = 1.;
        fhL1JPatch->Fill(posX,posY);
        
        Int_t ts = 0;
        trg.GetL1TimeSum(ts);
        emcalTrigL1[posY][posX] += ts;
        totSTU += ts;
      }
      
    }
  }
  
  if(totTRU > fMaxTRUSignal)printf("large totTRU %f\n",totTRU);
  if(totSTU > fMaxSTUSignal)printf("large totSTU %f\n",totSTU);
  
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
    if( v0A+v0C > fMaxV0Signal) printf("large v0A+v0C %f\n",v0A+v0C);
  }
  
  //Matrix with signal per channel
  for (Int_t i = 0; i < fgkFALTRORows-1; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols-1; j++) //check x,y direction for reading FOR ((0,0) = top left);
    {
      fhFORPos->Fill(fgkFALTROCols-j, i, emcalCell  [i][j]);
      fhL0Pos ->Fill(fgkFALTROCols-j, i, emcalTrigL0[i][j]);
      fhL1Pos ->Fill(fgkFALTROCols-j, i, emcalTrigL1[i][j]);
    }
  }
  
  //TRU checks
  Double_t ampFOR[30] = {0.}, ampL0[30] = {0.}, ampL1[30] = {0.};
  for (Int_t i = 0; i < fgkFALTRORows-1; i++) 
  {
    for (Int_t j = 0; j < fgkFALTROCols-1; j++) //A-side
    {
      
      //method to get TRU number
      Int_t FORid = -1;
      fGeometry->GetAbsFastORIndexFromPositionInEMCAL(j,i,FORid);
      Int_t iTRU = -1;
      Int_t iADC = -1;
      fGeometry->GetTRUFromAbsFastORIndex(FORid,iTRU,iADC);	
      
      if (iTRU >= 0)
      {
        ampFOR[iTRU] += emcalCell  [i][j];
        ampL0[iTRU]  += emcalTrigL0[i][j];
        ampL1[iTRU]  += emcalTrigL1[i][j];
      }
    }
  }
  
  
  //   Double_t ampFOR[30] = {0.}, ampL0[30] = {0.}, ampL1[30] = {0.};
  //   for (Int_t i = 0; i < fgkFALTRORows-1; i++) 
  //   {
  //     for (Int_t j = 0; j < fgkFALTROCols-1; j++) //A-side
  //     {
  //       if (i < 3)
  //       {
  //         if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[0] += emcalCell  [i][j];
  // 	    ampL0[0]  += emcalTrigL0[i][j];
  // 	    ampL1[0]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[15] += emcalCell  [i][j];
  // 	    ampL0[15]  += emcalTrigL0[i][j];
  // 	    ampL1[15]  += emcalTrigL1[i][j];
  // 	  }
	
  //       }
  //       else if (i > 3 && i < 8)
  //       {
  //         if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[1] += emcalCell  [i][j];
  // 	    ampL0[1]  += emcalTrigL0[i][j];
  // 	    ampL1[1]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[16] += emcalCell  [i][j];
  // 	    ampL0[16]  += emcalTrigL0[i][j];
  // 	    ampL1[16]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 8 && i < 12)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[2] += emcalCell  [i][j];
  // 	    ampL0[2]  += emcalTrigL0[i][j];
  // 	    ampL1[2]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[17] += emcalCell  [i][j];
  // 	    ampL0[17]  += emcalTrigL0[i][j];
  // 	    ampL1[17]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 12 && i < 16)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[3] += emcalCell  [i][j];
  // 	    ampL0[3]  += emcalTrigL0[i][j];
  // 	    ampL1[3]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[18] += emcalCell  [i][j];
  // 	    ampL0[18]  += emcalTrigL0[i][j];
  // 	    ampL1[18]  += emcalTrigL1[i][j];
  // 	  }
	
  //       }
  //       else if (i > 16 && i < 20)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[4] += emcalCell  [i][j];
  // 	    ampL0[4]  += emcalTrigL0[i][j];
  // 	    ampL1[4]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[19] += emcalCell  [i][j];
  // 	    ampL0[19]  += emcalTrigL0[i][j];
  // 	    ampL1[19]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 20 && i < 24)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[5] += emcalCell  [i][j];
  // 	    ampL0[5]  += emcalTrigL0[i][j];
  // 	    ampL1[5]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[20] += emcalCell  [i][j];
  // 	    ampL0[20]  += emcalTrigL0[i][j];
  // 	    ampL1[20]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //        else if (i > 24 && i < 28)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[6] += emcalCell  [i][j];
  // 	    ampL0[6]  += emcalTrigL0[i][j];
  // 	    ampL1[6]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[21] += emcalCell  [i][j];
  // 	    ampL0[21]  += emcalTrigL0[i][j];
  // 	    ampL1[21]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 28 && i < 32)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[7] += emcalCell  [i][j];
  // 	    ampL0[7]  += emcalTrigL0[i][j];
  // 	    ampL1[7]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[22] += emcalCell  [i][j];
  // 	    ampL0[22]  += emcalTrigL0[i][j];
  // 	    ampL1[22]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 32 && i < 36)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[8] += emcalCell  [i][j];
  // 	    ampL0[8]  += emcalTrigL0[i][j];
  // 	    ampL1[8]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[23] += emcalCell  [i][j];
  // 	    ampL0[23]  += emcalTrigL0[i][j];
  // 	    ampL1[23]  += emcalTrigL1[i][j];
  // 	  }
	
  //       }
  //       else if (i > 36 && i < 40)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[9] += emcalCell  [i][j];
  // 	    ampL0[9]  += emcalTrigL0[i][j];
  // 	    ampL1[9]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[24] += emcalCell  [i][j];
  // 	    ampL0[24]  += emcalTrigL0[i][j];
  // 	    ampL1[24]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 40 && i < 44)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[10] += emcalCell  [i][j];
  // 	    ampL0[10]  += emcalTrigL0[i][j];
  // 	    ampL1[10]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[25] += emcalCell  [i][j];
  // 	    ampL0[25]  += emcalTrigL0[i][j];
  // 	    ampL1[25]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 44 && i < 48)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[11] += emcalCell  [i][j];
  // 	    ampL0[11]  += emcalTrigL0[i][j];
  // 	    ampL1[11]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[26] += emcalCell  [i][j];
  // 	    ampL0[26]  += emcalTrigL0[i][j];
  // 	    ampL1[26]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //        else if (i > 48 && i < 52)
  //       {
  // 	if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[12] += emcalCell  [i][j];
  // 	    ampL0[12]  += emcalTrigL0[i][j];
  // 	    ampL1[12]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[27] += emcalCell  [i][j];
  // 	    ampL0[27]  += emcalTrigL0[i][j];
  // 	    ampL1[27]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 52 && i < 56)
  //       {
  //         if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[13] += emcalCell  [i][j];
  // 	    ampL0[13]  += emcalTrigL0[i][j];
  // 	    ampL1[13]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[28] += emcalCell  [i][j];
  // 	    ampL0[28]  += emcalTrigL0[i][j];
  // 	    ampL1[28]  += emcalTrigL1[i][j];
  // 	  }
  //       }
  //       else if (i > 56 && i < 60)
  //       {
  //         if(j-23 <= 0) //A-side
  // 	  {
  // 	    ampFOR[14] += emcalCell  [i][j];
  // 	    ampL0[14]  += emcalTrigL0[i][j];
  // 	    ampL1[14]  += emcalTrigL1[i][j];
  // 	  }
  // 	else //C-side
  // 	  {
  // 	    ampFOR[29] += emcalCell  [i][j];
  // 	    ampL0[29]  += emcalTrigL0[i][j];
  // 	    ampL1[29]  += emcalTrigL1[i][j];
  // 	  }	
  //       }
  
  //     }
  //   }
  
  for (Int_t i = 0; i < 30; i++)
  {
    if (ampFOR[i] != 0 && ampL1[i] != 0) fhFEESTU[i]->Fill(ampL1[i]/ampFOR[i]);
    if (ampL0[i] != 0 && ampL1[i] != 0) fhTRUSTU[i]->Fill(ampL1[i]/ampL0[i]);
  }
  
  //  Int_t TRUCheck[30] = {1};
  //   Int_t STUCheck[30] = {1};
  //   for (Int_t i = 0; i < 30; i++)
  //   {
  //     if (fhTRUSTU[i]->GetEntries()>0) if(fhTRUSTU[i]->Integral(10,20)/fhTRUSTU[i]->GetEntries() < 0.9) STUCheck[i] = 0;
  //     if (fhTRUSTU[i]->GetEntries()==0) STUCheck[i] = 0;
  //   }
  
  //   for (Int_t i = 0; i < 30; i++)
  //     {
  //       if (i<15) fhSTUChecks->Fill(0.,i,STUCheck[i]);
  //       else fhSTUChecks->Fill(1.,i,STUCheck[i]);
  //     }
  
  PostData(1, fOutputList);  
  
}

