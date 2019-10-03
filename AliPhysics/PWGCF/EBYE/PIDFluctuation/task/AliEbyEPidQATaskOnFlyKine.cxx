/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Offline.                                                       *
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
//             AliEbyE OnFLy QA Tasks for Charge and PID                   //
//                         For Testing Only                                //
//                   Satyajit Jena | sjena@cern.ch                         //
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"

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

#include "AliEbyEPidQATaskOnFlyKine.h"

ClassImp(AliEbyEPidQATaskOnFlyKine)

//-----------------------------------------------------------------------
AliEbyEPidQATaskOnFlyKine::AliEbyEPidQATaskOnFlyKine( const char *name )
: AliAnalysisTaskSE( name ),
  fThnList(0),fCentrality(-1), fEtaCut(0.5),fPtCut(5.),fVtxZ(30.),
  fHistImpNpart(NULL),  fHistImpMult(NULL),   
  fHistNpartMult(NULL),  fHistStat(NULL)  {
  for (Int_t i = 0; i < 4; ++i) {
    fHistPt[i][0]     = NULL;
    fHistPt[i][1]     = NULL;
    fHistEtaY[i]      = NULL;
    fHistPhi[i]       = NULL;
    fHistPhiPt[i]       = NULL;
    fHistMult[i][0]   = NULL;
    fHistMult[i][1]   = NULL;
    fHistMultTot[i][0]= NULL;
    fHistMultTot[i][1]= NULL;
  }
  DefineOutput(1, TList::Class()); 
}

AliEbyEPidQATaskOnFlyKine::~AliEbyEPidQATaskOnFlyKine() {
  if (fThnList) delete fThnList;
}

//---------------------------------------------------------------------------------
void AliEbyEPidQATaskOnFlyKine::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

  fHistImpNpart  = new TH2F("HistImpNpart","Impart Parameter ~ N_{part};N_{part};#it{b}",420,-0.5,419.5,20,-0.5,19.5);
  fHistImpMult   = new TH2F("HistImpMult","Impart Parameter ~ N_{tot};#it{b};N_{part}",20,-0.5,19.5, 9000,0,90000);
  fHistNpartMult = new TH2F("HistNpartMult","N_{part} ~ N_{tot};N_{part};Multiplicity",420,-0.5,419.5, 9000,0,90000);
  fHistStat = new TH1F("HistEventStat","Event Statistics;#it{b};N_{evt}",20,-0.5,19.5);
  
 fThnList->Add(fHistImpNpart);
 fThnList->Add(fHistImpMult);
 fThnList->Add(fHistNpartMult);  
 fThnList->Add(fHistStat);

  TString charge[] = {"Munus","Plus"};
 
  Int_t nBin[] = {30000,20000,15000,6000};
  Double_t gBinUp[] = {29999.5,19999.5,14999.5,5999.5};

  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    fHistPhi[iPid] = new TH2F(Form("HistPhi%s",fgkPidName[iPid]), Form("Impact Parameter vs #phi : (%s+%s);#it{b} (fm);#phi",fgkPidLatex[iPid][0],fgkPidLatex[iPid][1]),20,-0.5,19.5,80,-1,7);

    fHistPhiPt[iPid] = new TH2F(Form("HistPhiPt%s",fgkPidName[iPid]), Form("#phi vs p_{T}: (%s+%s);#it{b} (fm);#phi",fgkPidLatex[iPid][0],fgkPidLatex[iPid][1]),80,-1,7,1000,0,20);

    TString title = (iPid == 0) ? Form("Impact Parameter vs #eta : (%s+%s);b(fm);#eta",fgkPidLatex[iPid][0],fgkPidLatex[iPid][1]) : Form("Impact Parameter vs #it{y} : (%s+%s);#it{b} (fm);#it{y}",fgkPidLatex[iPid][0],fgkPidLatex[iPid][1]);
    fHistEtaY[iPid]  = new TH2F(Form("HistEtaY%s",fgkPidName[iPid]),title.Data(),20,-0.5,19.5,400,-20,20);

    fThnList->Add(fHistPhi[iPid]);
    fThnList->Add(fHistPhiPt[iPid]);
    fThnList->Add(fHistEtaY[iPid]);
    
    for (Int_t iSg = 0; iSg < 2; ++iSg) {
      fHistPt[iPid][iSg] = new TH2F(Form("HistPt%s%s",fgkPidName[iPid],charge[iSg].Data()), Form("Impact Parameter vs p_{T} : %s;#it{b} (fm);p_{T}",fgkPidLatex[iPid][iSg]), 20,-0.5,19.5,1000,0,20);
      fHistMult[iPid][iSg] = new TH2F(Form("HistMultiplicity%s%s",fgkPidName[iPid],charge[iSg].Data()), Form("Impact Parameter vs Multplicity : %s;#it{b} (fm);Multiplicity",fgkPidLatex[iPid][iSg]), 20,-0.5,19.5,nBin[iPid],0.5,gBinUp[iPid]);
      fHistMultTot[iPid][iSg] = new TH2F(Form("HistTotalMultiplicity%s%s",fgkPidName[iPid],charge[iSg].Data()), Form("Impact Parameter vs Multiplicity : %s;#it{b} (fm);Multiplicity",fgkPidLatex[iPid][iSg]), 20,-0.5,19.5,nBin[iPid],0.5,gBinUp[iPid]);
      fThnList->Add(fHistPt[iPid][iSg]);
      fThnList->Add(fHistMult[iPid][iSg]);
      fThnList->Add(fHistMultTot[iPid][iSg]);
    }
   
  }
  PostData(1, fThnList);
}


//----------------------------------------------------------------------------------
void AliEbyEPidQATaskOnFlyKine::UserExec( Option_t * ){
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
 
  const AliVVertex *vtxMC = mcEvent->GetPrimaryVertex();
  if (vtxMC->GetZ() > fVtxZ) return;
 
  AliStack *stack = mcEvent->Stack();
 
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  if(!genHeader){
    printf("  Event generator header not available!!!\n");
    return;
  }
  
  Double_t imp   = 21;
  Double_t npart = 0;

  Double_t np = -1;
  Double_t nt = -1;  

  if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
    imp = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
    nt  = ((AliGenHijingEventHeader*) genHeader)->TargetParticipants();
    np  = ((AliGenHijingEventHeader*) genHeader)->ProjectileParticipants();
  }  

 

  if (imp > -1 && imp < 20) fCentrality = Int_t(imp);
  else return;
  npart = nt + np;
  if (npart < 0 ||  npart > 419) return; 
  
  fHistStat->Fill(fCentrality);
  fHistImpNpart->Fill(npart,fCentrality);
  Int_t nParticles = stack->GetNtrack();
  
  fHistImpMult->Fill(fCentrality,nParticles);
  fHistNpartMult->Fill(npart,nParticles);

  Float_t multFull[4][2]; Float_t mult[4][2];
  for (Int_t i = 0; i < 4; ++i) {
    mult[i][0] = 0; mult[i][1] = 0;
    multFull[i][0] = 0; multFull[i][1] = 0;
  }

  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
    TParticle* part = stack->Particle(iParticle);
    if (!part) {
      Printf(" No Particle Available ");
      continue;
    }

    if(!part->GetPDG()) continue;
    if (part->GetPDG()->Charge() == 0.) continue;

    Int_t iCharge = part->GetPDG()->Charge() < 0 ? 0 : 1;
    
    Int_t iPid = -1;
    Int_t pid = part->GetPdgCode();
    if(pid == 211)         { iPid = 1; } 
    else if(pid ==  -211)  { iPid = 1; } 
    else if(pid ==   321)  { iPid = 2; } 
    else if(pid ==  -321)  { iPid = 2; } 
    else if(pid ==  2212)  { iPid = 3; } 
    else if(pid == -2212)  { iPid = 3; } 
    else    iPid = -1;   
    
    multFull[0][iCharge] += 1;
    if (iPid > 0 && iPid < 4)
      multFull[iPid][iCharge] += 1;

    Float_t pt = part->Pt();    
    if (pt > fPtCut) continue;
    
    mult[0][iCharge]     += 1;
    if (iPid < 1 || iPid > 3)   continue;   
    mult[iPid][iCharge] += 1;

    Float_t gEta = part->Eta();
    Float_t gY   = part->Y();

    fHistEtaY[0]->Fill(fCentrality,gEta);
    fHistEtaY[iPid]->Fill(fCentrality,gY);
    
    if (TMath::Abs(part->Eta()) > fEtaCut)  continue;
    Float_t phi = part->Phi();
    
    fHistPhiPt[0]->Fill(phi,pt);
    fHistPhiPt[iPid]->Fill(phi,pt);

    fHistPt[0][iCharge]->Fill(fCentrality,pt);
    fHistPt[iPid][iCharge]->Fill(fCentrality,pt);

    fHistPhi[0]->Fill(fCentrality,phi);
    fHistPhi[iPid]->Fill(fCentrality,phi);

  }
  
  for (Int_t i = 0; i < 4; ++i) {
    fHistMult[i][0]->Fill(fCentrality,mult[i][0]);
    fHistMult[i][1]->Fill(fCentrality,mult[i][1]);

    fHistMultTot[i][0]->Fill(fCentrality,multFull[i][0]);
    fHistMultTot[i][1]->Fill(fCentrality,multFull[i][1]);
    
  }
  
  PostData(1, fThnList);
}

void AliEbyEPidQATaskOnFlyKine::Terminate( Option_t * ){
  Info("AliEbyEPidQATaskOnFlyKine"," Task Successfully finished");
}
//________________________________________________________________________
const Char_t* AliEbyEPidQATaskOnFlyKine::fgkPidName[4]      = {"Nch","Npi","Nka","Npr"};
//________________________________________________________________________
const Char_t* AliEbyEPidQATaskOnFlyKine::fgkPidLatex[4][2]  = {{"N_{-}","N_{+}"}, {"N_{#pi^{-}}","N_{#pi^{+}}"},{"N_{K^{-}}","N_{K^{+}}"}, {"N_{#bar{p}}","N_{p}"}};
//________________________________________________________________________

