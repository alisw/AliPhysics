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

/**************************************************************************
 * Analysis Class Implimentation for the MC Truth and responce Matrix 
 *    Auther:    Satyajit Jena, IIT Bombay |  sjena@cern.ch
 * 
 *                     Mon Nov 22 19:54:27 CET 2010
 * 
 * You may need some type of tunning to histograms and filling schemes,
 * Please contact me if you want something to change this code 
 **************************************************************************/

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliESD.h"
#include "AliESDEvent.h"

#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliPID.h"

#include "AliESDInputHandler.h"
#include "AliPMDAnalysisMCTaskSE.h"

ClassImp(AliPMDAnalysisMCTaskSE)

//________________________________________________________________________
AliPMDAnalysisMCTaskSE::AliPMDAnalysisMCTaskSE(const char *name) 
  : AliAnalysisTaskSE(name),
  fPhysList(0),
//  fhCounter(0),
//  fhVtx(0),
//  fhVty(0),
//  fhVtz(0),
  fhResponseAll(0),
  fhTrueAll(0),
  fhMeasuredAll(0)
    
{
  for (Int_t i = 0; i < 10; i++) {
    fhResponse[i] = NULL;
    fhTrue[i] = NULL;
    fhMeasured[i] = NULL;
  }
  
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
void AliPMDAnalysisMCTaskSE::UserCreateOutputObjects()
{
  fPhysList   = new TList();
  fPhysList->SetOwner();
    
  //___________ Counter ________________
  // fCounter = new TH1D("hCounter","EVENT COUNTERS", 100, 0, 100);;
  // fPhysList->Add(fCounter);

  Int_t   kMultBin = 100;  
  Float_t gMultMin = -0.5; 
  Float_t gMultMax = 99.5;
  

  fhResponseAll = new TH2F("hResponseAll","Responce Matrix For All",
			   kMultBin,gMultMin,gMultMax,kMultBin,gMultMin,gMultMax);
  
  fhResponseAll->GetXaxis()->SetTitle("Measured Multiplicity #gamma_{like}");
  fhResponseAll->GetXaxis()->SetTitle("Ture Multiplicity #gamma");

  fhTrueAll = new TH2F("hTrueAll","True Multiplicity of #gamma",kMultBin,gMultMin,gMultMax);
  fhTrueAll->GetXaxis()->SetTitle("True Multiplicity #gamma");

  fhMeasuredAll = new TH2F("hMeasuredAll","Measured Multiplicity of #gamma",kMultBin,gMultMin,gMultMax);
  fhMeasuredAll->GetXaxis()->SetTitle("Measured Multiplicity #gamma");
 
  fPhysList->Add(fhResponseAll);
  fPhysList->Add(fhTrueAll);
  fPhysList->Add(fhMeasuredAll);
  
  for(Int_t i = 0; i < 10; i++) {
    Float_t gbin = 2.3 + i*0.2;
    
    char a[40];
    sprintf(a,"hResponseEtaBin%2.1fto%2.1f",gbin, gbin+0.2);
    char b[40];
    sprintf(b,"hTrueEtaBin%2.1fto%2.1f",gbin, gbin+0.2);
    char c[40];
    sprintf(c,"hMeasuredEtaBin%2.1fto%2.1f",gbin, gbin+0.2);

    char aa[40];
    sprintf(aa," Responce Matrix from \eta - %2.1f to %2.1f",gbin, gbin+0.2);
    char bb[40];
    sprintf(bb," True Multiplicity of #gamma from \eta %2.1f to %2.1f",gbin, gbin+0.2);
    char cc[40];
    sprintf(cc," Measured multiplicity of gamma [#gamma_{like}] from \eta %2.1f to %2.1f",gbin, gbin+0.2);
    

    fhResponse[i] = new TH2F(a,aa,kMultBin,gMultMin,gMultMax,kMultBin,gMultMin,gMultMax);
    fhResponse[i]->GetXaxis()->SetTitle("Measured Multiplicity #gamma_{like}");
    fhResponse[i]->GetYaxis()->SetTitle("Ture Multiplicity #gamma");


    fhTrue[i] = new TH1F(b,bb,kMultBin,gMultMin,gMultMax);
    fhResponse[i]->GetXaxis()->SetTitle("True Multiplicity #gamma");

    fhMeasured[i] = new TH1F(c,cc,kMultBin,gMultMin,gMultMax);
    fhMeasured[i]->GetXaxis()->SetTitle("True Multiplicity #gamma");

    fPhysList->Add(fhResponse[i]);
    fPhysList->Add(fhTrue[i]);
    fPhysList->Add(fhMeasured[i]);

  }

  PostData(1,fPhysList);
}

//________________________________________________________________________
void AliPMDAnalysisMCTaskSE::UserExec(Option_t *) 
{
  
  fCntr++;
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  AliStack *stack = 0;
  AliMCEventHandler *mcH = 0;
  AliMCEvent* mcEvent = 0;
  
  if(MCEvent()){
    stack = MCEvent()->Stack();
    mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    mcEvent = mcH->MCEvent();
    AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
  }
  if(!stack) return;
  if(!esd) return;
  if(!mcEvent)return;
  
  EventByEvent(esd,mcEvent);
  
}

//----------------------------
void AliPMDAnalysisMCTaskSE::EventByEvent(AliESDEvent* esd, AliMCEvent* mcEvent)
{ 
  Int_t trnn = mcEvent->GetNumberOfTracks();

  Float_t McEta = 0.;
 
  Int_t kcntr[0] ={0};
  Int_t kcntrr = 0;

  AliStack *stack = mcEvent->Stack();
  Int_t sln = 0;
  for (Int_t iTracks = 0; iTracks < trnn; iTracks++) {
    AliMCParticle* track = (AliMCParticle*)mcEvent->GetTrack(iTracks);
    // TParticle* part = track->Particle();

    TParticle* part = stack->Particle(iTracks);
    if (!part) continue;
    Int_t pid = part->GetPdgCode();
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    TParticlePDG* pdgPart = pdgDB->GetParticle(pid);
    Float_t charge = (pdgPart ? pdgPart->Charge() : 0);
    
    McEta = part->Eta();

    if(McEta < 2.0 || McEta > 4.2) continue;

    kcntrr++;
    if (McEta >= 2.1 && McEta < 2.3 ) kcntr[0]++;
    else if (McEta >= 2.3 && McEta < 2.5 ) kcntr[1]++;
    else if (McEta >= 2.5 && McEta < 2.7 ) kcntr[2]++;
    else if (McEta >= 2.7 && McEta < 2.9 ) kcntr[3]++;
    else if (McEta >= 2.9 && McEta < 3.1 ) kcntr[4]++;
    else if (McEta >= 3.1 && McEta < 3.3 ) kcntr[5]++;
    else if (McEta >= 3.3 && McEta < 3.5 ) kcntr[6]++;
    else if (McEta >= 3.5 && McEta < 3.7 ) kcntr[7]++;
    else if (McEta >= 3.7 && McEta < 3.9 ) kcntr[8]++;  
    else if (McEta >= 3.9 && McEta < 4.1 ) kcntr[9]++;
    
  }
  
  // ESD
  
  Int_t npmdcl   = esd->GetNumberOfPmdTracks();
  Float_t rpxpy  = 0.;
  Float_t theta  = 0.; 
  Float_t pEta   = 0.;
 
  Int_t cntr[10] = {0};
  Int_t cntrr = 0;

  while (npmdcl--)
    {
      
      AliESDPmdTrack *pmdtr = esd->GetPmdTrack(npmdcl);
      Int_t det     = pmdtr->GetDetector();
      if( det == 1) continue;
      Int_t   ncell = pmdtr->GetClusterCells();
      if(ncell < 2) continue;
      Float_t adc   = pmdtr->GetClusterADC();
      if( adc < 216 ) continue;

      Float_t clsX  = pmdtr->GetClusterX();
      Float_t clsY  = pmdtr->GetClusterY();
      Float_t clsZ  = pmdtr->GetClusterZ();
      
      rpxpy  =  TMath::Sqrt(clsX*clsX + clsY*clsY);
      theta  =  TMath::ATan2(rpxpy,clsZ);
      pEta   = -TMath::Log(TMath::Tan(0.5*theta));
      cntrr++;
      if (pEta >= 2.1 && pEta < 2.3 ) cntr[0]++;
      else if (pEta >= 2.3 && pEta < 2.5 ) cntr[1]++;
      else if (pEta >= 2.5 && pEta < 2.7 ) cntr[2]++;
      else if (pEta >= 2.7 && pEta < 2.9 ) cntr[3]++;
      else if (pEta >= 2.9 && pEta < 3.1 ) cntr[4]++;
      else if (pEta >= 3.1 && pEta < 3.3 ) cntr[5]++;
      else if (pEta >= 3.3 && pEta < 3.5 ) cntr[6]++;
      else if (pEta >= 3.5 && pEta < 3.7 ) cntr[7]++;
      else if (pEta >= 3.7 && pEta < 3.9 ) cntr[8]++;  
      else if (pEta >= 3.9 && pEta < 4.1 ) cntr[9]++;
    }
  
    fhResponseAll->Fill(cntrr,kcntrr);
    fhTrueAll->Fill(kcntrr);
    fhMeasuredAll->Fill(cntrr);

  for(int j = 0; j < 10; j++ ) {
    fhResponse[i]->Fill(cntr[i],kcntr[i]);
    fhTrue[i]->Fill(kcntr[i]);
    fhMeasured[i]->Fill(cntr[i]);
  }
}


//________________________________________________________________________
void AliPMDAnalysisMCTaskSE::Terminate(Option_t *) 
{

  printf("Event number %d\n", fCntr);

     
}
