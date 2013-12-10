/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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


//=========================================================================//
//                 AliEbyE Fluctuaion Analysis for PID                     //
//                   Deepika Jena | drathee@cern.ch                        //
//                   Satyajit Jena | sjena@cern.ch                         //
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

#include "THnSparse.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "TDatabasePDG.h"

#include "AliEbyEPidTaskFastGen.h"

ClassImp(AliEbyEPidTaskFastGen)

//-----------------------------------------------------------------------
AliEbyEPidTaskFastGen::AliEbyEPidTaskFastGen( const char *name )
: AliAnalysisTaskSE( name ),
  fThnList(0),
  fVxMax(3.),
  fVyMax(3.),
  fVzMax(15.),
  fPtLowerLimit(0.2),
  fPtHigherLimit(5.),
  fEtaLowerLimit(-1.),
  fEtaHigherLimit(1.),
  fHistoCorrelationMC(NULL)
{ 
 
  DefineOutput(1, TList::Class()); 
}

AliEbyEPidTaskFastGen::~AliEbyEPidTaskFastGen()
{
  if (fThnList) delete fThnList;
}


//---------------------------------------------------------------------------------
void AliEbyEPidTaskFastGen::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);
  Int_t fgSparseDataBins[kNSparseData]   = {50000,   500,     440,  30000,   15000,  15000,   10000,  5000,   5000,  5000,  2500,  2500,  2500,  1250, 1250};
  Double_t fgSparseDataMin[kNSparseData] = {0.,       0.,      0.,     0.,      0.,     0.,      0.,    0.,     0.,    0.,    0.,    0.,    0.,   0.,    0.};
  Double_t fgSparseDataMax[kNSparseData] = {50000.,  100.,    440., 30000.,  15000., 15000.,  10000., 5000.,  5000., 5000., 2500., 2500., 2500., 1250., 1250.};
  
  const Char_t *fgkSparseDataTitle[] = {
    "Nnumber of Partiles",
    "Impact Parameters",
    "N_{p}",
    "N_{ch}",
    "N_{+}",
    "N_{-}",
    "N_{#pi}",
    "N_{#pi^{+}}",
    "N_{#pi^{-}}",
    "N_{K}",
    "N_{K^{+}}",
    "N_{K^{-}}",
    "N_{pr}",
    "N_{p}",
    "N_{#bar{p}}"
  };
  
  fHistoCorrelationMC = new THnSparseI("hHistoPR", "", kNSparseData, fgSparseDataBins, fgSparseDataMin, fgSparseDataMax);
  for (Int_t iaxis = 0; iaxis < kNSparseData; iaxis++)
    fHistoCorrelationMC->GetAxis(iaxis)->SetTitle(fgkSparseDataTitle[iaxis]);
  fThnList->Add(fHistoCorrelationMC);
  
  PostData(1, fThnList);
}


//----------------------------------------------------------------------------------
void AliEbyEPidTaskFastGen::UserExec( Option_t * ){
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  AliStack *stack = mcEvent->Stack();
 
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  if(!genHeader){
    printf("  Event generator header not available!!!\n");
    return;
  }
  
  Double_t imp = -1;
  Double_t np = -1;
  Double_t nt = -1;
  
  if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
    imp = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
    nt  = ((AliGenHijingEventHeader*) genHeader)->TargetParticipants();
    np  = ((AliGenHijingEventHeader*) genHeader)->ProjectileParticipants();
  }  
  

  Int_t count[4][2];
  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 2; j++) count[i][j] = 0; 
  }
  
  Int_t nParticles = stack->GetNtrack();
  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
    TParticle* part = stack->Particle(iParticle);
    if (!part) {
      Printf(" No Particle Available ");
      continue;
    }
    
    //  TDatabasePDG* pdgDB   = TDatabasePDG::Instance();
    //  TParticlePDG* pdgPart = pdgDB->GetParticle(part->GetPdgCode());
    
    Float_t pt = part->Pt();
    if(pt < fPtLowerLimit || pt > fPtHigherLimit ) continue;
    Float_t gEta = part->Eta();
    if(gEta < fEtaLowerLimit || gEta > fEtaHigherLimit ) continue;
    
    
    //  Float_t gCharge = (pdgPart ? pdgPart->Charge() : 0);
    
    Int_t pid = part->GetPdgCode();
    if(pid == 211)         { count[1][0]++;  count[0][0]++; } 
    else if(pid ==  -211)  { count[1][1]++;  count[0][1]++; } 
    else if(pid ==   321)  { count[2][0]++;  count[0][0]++; } 
    else if(pid ==  -321)  { count[2][1]++;  count[0][1]++; } 
    else if(pid ==  2212)  { count[3][0]++;  count[0][0]++; } 
    else if(pid == -2212)  { count[3][1]++;  count[0][1]++; } 
  }
  
  Double_t vsparseMC[kNSparseData];
  vsparseMC[0]  = nParticles;
  vsparseMC[1]  = imp;
  vsparseMC[2]  = np + nt;
  vsparseMC[3]  = count[0][0] + count[0][1];
  vsparseMC[4]  = count[0][0];
  vsparseMC[5]  = count[0][1];

  vsparseMC[6]  = count[1][0] + count[1][1];
  vsparseMC[7]  = count[1][0];
  vsparseMC[8]  = count[1][1];

  vsparseMC[9]  = count[2][0] + count[2][1];
  vsparseMC[10]  = count[2][0];
  vsparseMC[11]  = count[2][1];

  vsparseMC[12] = count[3][0] + count[3][1];
  vsparseMC[13] = count[3][0];
  vsparseMC[14] = count[3][1];
  

   /* Printf(" %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f",
  	 vsparseMC[0], vsparseMC[1], vsparseMC[2], 
	 vsparseMC[3], vsparseMC[4], vsparseMC[5], 
	 vsparseMC[6], vsparseMC[7], vsparseMC[8], 
	 vsparseMC[9], vsparseMC[10], vsparseMC[11],
	 vsparseMC[12], vsparseMC[13], vsparseMC[13]);
   */
  fHistoCorrelationMC->Fill(vsparseMC);
  PostData(1, fThnList);
  
}

void AliEbyEPidTaskFastGen::Terminate( Option_t * ){

  Info("AliEbyEPiDTask"," Task Successfully finished");
  
}


