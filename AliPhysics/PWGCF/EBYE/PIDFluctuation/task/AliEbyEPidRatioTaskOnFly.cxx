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

#include "AliEbyEPidRatioTaskOnFly.h"

ClassImp(AliEbyEPidRatioTaskOnFly)

//-----------------------------------------------------------------------
AliEbyEPidRatioTaskOnFly::AliEbyEPidRatioTaskOnFly( const char *name )
: AliAnalysisTaskSE( name ),
  fThnList(0),
  fPtLowerLimit(0.2),
  fPtHigherLimit(5.),
  fEtaLowerLimit(-1.),
  fEtaHigherLimit(1.),
  fCentrality(-1),
  fOrder(8),
  fRedFactp(NULL)
  
 
{ 
 
  DefineOutput(1, TList::Class()); 
}

AliEbyEPidRatioTaskOnFly::~AliEbyEPidRatioTaskOnFly()
{
  if (fThnList) delete fThnList;
  
  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    if (fRedFactp[ii]) delete[] fRedFactp[ii];
  if (fRedFactp) delete[] fRedFactp;
  
}




//---------------------------------------------------------------------------------
void AliEbyEPidRatioTaskOnFly::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

  const Char_t *name = "Mc";
  const Char_t *title = Form(" #eta [%2.1f-%2.1f]",fEtaLowerLimit,fEtaHigherLimit);


  TString sName(name);
  TString sTitle(title);
 
    
  fRedFactp = new Double_t*[fOrder+1];
  for (Int_t ii = 0 ; ii <= fOrder; ++ii)
    fRedFactp[ii] = new Double_t[2];
  
  //TList *list[4];
  fThnList->Add(new TList);
  TList *list =  static_cast<TList*>(fThnList->Last());
  list->SetName(Form("f%s",name));
  list->SetOwner(kTRUE);
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    TString sNetTitle(Form("%s - %s",fgkPidLatex[iPid][1],fgkPidLatex[iPid][0]));

    list->Add(new TProfile(Form("fProfTot%sPlus%s",fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   20,-0.5,19.5));

    list->Add(new TProfile(Form("fProfTot%sMinus%s",fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   20,-0.5,19.5));
    
    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("fProf%s%sNet%dM",fgkPidName[iPid],name, idx), 
			     Form("(%s)^{%d} : %s;Centrality(100);(%s)^{%d}",sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     20,-0.5,19.5));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile(Form("fProf%s%sNetF%02d%02d",fgkPidName[iPid], name, ii, kk),
			       Form("f_{%02d%02d} : %s;Centrality(100);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
			       20,-0.5,19.5));
      }
    }
  
  }  
			     

  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile(Form("fProf%sNu%02d",name,iPhy),Form("Physics Variable for index %d | %s ; Centrality;",iPhy,name),20,-0.5,19.5));
  }
  
 
  PostData(1, fThnList);
}


//----------------------------------------------------------------------------------
void AliEbyEPidRatioTaskOnFly::UserExec( Option_t * ){
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
 
  for (Int_t ii = 0 ; ii < 4; ++ii) 
    for (Int_t kk = 0 ; kk < 2; ++kk) 
    fNp[ii][kk] = 0;
  

  const AliVVertex *vtxMC = mcEvent->GetPrimaryVertex();
  if (vtxMC->GetZ() > 10.) return;
 
  AliStack *stack = mcEvent->Stack();
 
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  if(!genHeader){
    printf("  Event generator header not available!!!\n");
    return;
  }
  
  Double_t imp = 21;
  
  if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
    imp = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
  }  
  if (imp < 20) fCentrality = Int_t(imp);
  else return;

  Int_t nParticles = stack->GetNtrack();
  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
    TParticle* part = stack->Particle(iParticle);
    if (!part) {
      Printf(" No Particle Available ");
      continue;
    }
    
    Float_t pt = part->Pt();
    if(pt < fPtLowerLimit || pt > fPtHigherLimit ) continue;
    Float_t gEta = part->Eta();
    if(gEta < fEtaLowerLimit || gEta > fEtaHigherLimit ) continue;
    
    //  Float_t gCharge = (pdgPart ? pdgPart->Charge() : 0);
    
    Int_t pid = part->GetPdgCode();
    if(pid == 211)         { fNp[1][1]++;  fNp[0][1]++; } 
    else if(pid ==  -211)  { fNp[1][0]++;  fNp[0][0]++; } 
    else if(pid ==   321)  { fNp[2][1]++;  fNp[0][1]++; } 
    else if(pid ==  -321)  { fNp[2][0]++;  fNp[0][0]++; } 
    else if(pid ==  2212)  { fNp[3][1]++;  fNp[0][1]++; } 
    else if(pid == -2212)  { fNp[3][0]++;  fNp[0][0]++; } 
    
  }
    
  /*Printf("%6d %6d %6d %6d %6d %6d %6d %6d %6d", fCentrality,  
	 fNp[0][0], fNp[0][1], 
	 fNp[1][0], fNp[1][1], 
	 fNp[2][0], fNp[2][1], 
	 fNp[3][0], fNp[3][1]);*/


  FillHistSetCent();
  PostData(1, fThnList);
}

void AliEbyEPidRatioTaskOnFly::Terminate( Option_t * ){
  Info("AliEbyEPidRatioTaskOnFly"," Task Successfully finished");
}
//________________________________________________________________________
const Char_t* AliEbyEPidRatioTaskOnFly::fgkPidName[4]      = {"Nch","Npi","Nka","Npr"};
//________________________________________________________________________
const Char_t* AliEbyEPidRatioTaskOnFly::fgkPidLatex[4][2]  = {{"N_{-}","N_{+}"}, {"N_{#pi^{-}}","N_{#pi^{+}}"},{"N_{K^{-}}","N_{K^{+}}"}, {"N_{#bar{p}}","N_{p}"}};
//________________________________________________________________________
const Char_t* AliEbyEPidRatioTaskOnFly::fgkPidTitles[4][2] = {{"Negative","Positive"},{"Anti-Pions","Pions"},{"Anti-Kaons","Kaons"}, {"Anti-Protons","Protons"}};
//________________________________________________________________________
void AliEbyEPidRatioTaskOnFly::FillHistSetCent()  {

  const Char_t *name = "Mc";
    
  TList *list = static_cast<TList*>(fThnList->FindObject(Form("f%s",name)));
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t deltaNp = fNp[iPid][1]-fNp[iPid][0];  
    Double_t delta = 1.;
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sPlus%s", fgkPidName[iPid], name))))->Fill(fCentrality, fNp[iPid][1]);
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sMinus%s", fgkPidName[iPid], name))))->Fill(fCentrality, fNp[iPid][0]);
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;
      (static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNet%dM", fgkPidName[iPid], name, idxOrder))))->Fill(fCentrality, delta);
    }
    
    for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = 1.;
      fRedFactp[idxOrder][1]  = 1.;
    }
    
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(fNp[iPid][0]-(idxOrder-1));
      fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(fNp[iPid][1]-(idxOrder-1));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {  
      for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   
	(static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNetF%02d%02d", fgkPidName[iPid], name, ii, kk))))->Fill(fCentrality, fik);
      }
    }
  }
 
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 fNp[0][1],  fNp[0][0], 
  //	 fNp[1][1],  fNp[1][0], 
  ///	 fNp[2][1],  fNp[2][0], 
  //	 fNp[3][1],  fNp[3][0]);
  //

   Int_t a[6][4]; Int_t b[22];
   for (Int_t iPid = 0; iPid < 4; ++iPid) {
     a[0][iPid] = fNp[iPid][1]+fNp[iPid][0];       // 0  n+ + n-
     a[1][iPid] = fNp[iPid][1];                        // 1  n+
     a[2][iPid] = fNp[iPid][0];                        // 2  n-
     a[3][iPid] = fNp[iPid][1]*fNp[iPid][0];       // 3  n+ . n-
     a[4][iPid] = fNp[iPid][1]*(fNp[iPid][1]-1);   // 4  n+ (n+ - 1)
     a[5][iPid] = fNp[iPid][0]*(fNp[iPid][0]-1);   // 5  n- (n- - 1)
     
     // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
     //	   a[0][iPid], a[1][iPid], a[2][iPid]);

  }
  
  b[0]  = a[0][0]*a[0][2];       // 24 N   K
  b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
  b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
  b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
  b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
  b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
  
  b[6]  = a[0][0]*a[0][3];       // 30 N   P
  b[7]  = a[0][2]*a[0][3];       // 31 K   P
  b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
  b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
  b[10] = a[2][2]*a[1][3];       // 34 k-  p+
  b[11] = a[2][2]*a[2][3];       // 35 k-  p-
  
  b[12] = a[0][0]*a[0][1];       // 36 N  Pi
  b[13] = a[0][3]*a[0][1];       // 37 P  Pi
  b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
  b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
  b[16] = a[2][3]*a[1][1];       // 40 p- pi+
  b[17] = a[2][3]*a[2][1];       // 41 p- pi-
  
  b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
  b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
  b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
  b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
  // TList *list_nu = static_cast<TList*>(fOutList->FindObject(Form("f%s_nu",name)));
  Int_t k = 0;
  for (Int_t j = 0; j < 4; j++) {
    for (Int_t i = 0; i < 6; i++) {
      (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,k))))->Fill(fCentrality,a[i][j]); 
      k++;
    }
  }

  for (Int_t j = 0; j < 22; j++) {
    (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,j+24))))->Fill(fCentrality,b[j]); 
  }
  
  /*
  Printf("%6d %6d %6d %6d %6d %6d %6d %6d %6d", fCentrality,  
	 fNp[0][0], fNp[0][1], 
	 fNp[1][0], fNp[1][1], 
	 fNp[2][0], fNp[2][1], 
	 fNp[3][0], fNp[3][1]);
  */

  return;
}

