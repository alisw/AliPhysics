/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Class AliAnalysisTaskSEmcCorr
// AliAnalysisTaskSE for studying HF-(hadron,electrons) and hadron-hadron correlations
//     at MC level
// Author: Andrea Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////


#include <TH1F.h>
#include <TH2F.h>
#include <TAxis.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include <TArrayI.h>
#include <TString.h>
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNormalizationCounter.h"
//#include "/Users/administrator/ALICE/CHARM/HFCJ/TestsForProduction/AliAnalysisTaskSEmcCorr.h"
#include "AliAnalysisTaskSEmcCorr.h"
#include "AliGenEventHeader.h"

class TCanvas;
class TTree;
class TChain;
class AliAnalysisTaskSE;


ClassImp(AliAnalysisTaskSEmcCorr)

AliAnalysisTaskSEmcCorr::AliAnalysisTaskSEmcCorr() 
: AliAnalysisTaskSE(),
  fReadMC(kTRUE),
  fFastSimul(kFALSE),
  fCheckDecay(kTRUE),
  fDoHFCorrelations(kTRUE),
  fOnlyKine(kTRUE),
  fDoHadronHadron(kFALSE),
  fNentries(0),
  fhNprongsD0(0),
  fhNprongsD0chargedOnly(0),
  fhNprongsD0chargedRef(0),
  fYrangeTrig(0.5),
  fEtarangeEleTrig(0.8),
  fminDpt(1.),
  fArrayDaugh(0x0),
  flastdaugh(0),
  fArrayAssoc(0x0),
  fLastAss(0),
  fminPtass(0.3),
  fmaxEtaAss(99.),
  fhMCcorrel(0x0),
  fhMCtrigPart(0x0),
  fhMChadroncorrel(0x0),
  fhMChadrontrigPart(0x0),
  fhDzeroDecay(0x0),  
  fhDplusDecay(0x0),
  fhLambdaCDecay(0x0),
  fhAllBDecay(0x0),
  fnpionPlus(0),                         
  fnpionMinus(0),                       
  fnpionZero(0),                         
  fnkaonPlus(0),                        
  fnkaonMinus(0),                       
  fnkaonZeroShort(0),                        
  fnProton(0),                        
  fnAntiProton(0),                       
  fnElePlus(0),                        
  fnEleMinus(0),                        
  fnMuonPlus(0),
  fnMuonMinus(0),
  fnNeutrino(0),
  fnJpsi(0),                            
  fnPhysPrim(0),
  fpxtotdaugh(0),
  fpytotdaugh(0),
  fpztotdaugh(0),
  fGeneratorString(0)
{// standard constructor
  
}



//________________________________________________________________________
AliAnalysisTaskSEmcCorr::AliAnalysisTaskSEmcCorr(const char *name) 
  : AliAnalysisTaskSE(name),
    fReadMC(kTRUE),
    fFastSimul(kFALSE),
    fCheckDecay(kTRUE),
    fDoHFCorrelations(kTRUE),
    fOnlyKine(kTRUE),
    fDoHadronHadron(kFALSE),
    fNentries(0),
    fhNprongsD0(0),
    fhNprongsD0chargedOnly(0),
    fhNprongsD0chargedRef(0),
    fYrangeTrig(0.5),
    fEtarangeEleTrig(0.8),
    fminDpt(1.),
    fArrayDaugh(0x0),
    flastdaugh(0),
    fArrayAssoc(0x0),
    fLastAss(0),
    fminPtass(0.3),
    fmaxEtaAss(99.),
    fhMCcorrel(0x0),
    fhMCtrigPart(0x0),
    fhMChadroncorrel(0x0),
    fhMChadrontrigPart(0x0),
    fhDzeroDecay(0x0),
    fhDplusDecay(0x0),
    fhLambdaCDecay(0x0),
    fhAllBDecay(0x0),
    fnpionPlus(0),                         
    fnpionMinus(0),                       
    fnpionZero(0),                         
    fnkaonPlus(0),                        
    fnkaonMinus(0),                       
    fnkaonZeroShort(0),                        
    fnProton(0),                        
    fnAntiProton(0),                       
    fnElePlus(0),                        
    fnEleMinus(0),                        
    fnMuonPlus(0),
    fnMuonMinus(0),
    fnNeutrino(0),
    fnJpsi(0),                            
    fnPhysPrim(0),
    fpxtotdaugh(0),
    fpytotdaugh(0),
    fpztotdaugh(0),
    fGeneratorString(0)
{ // default constructor
  
  DefineOutput(1, TH1F::Class());
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TH1F::Class());
  DefineOutput(5, THnSparseF::Class());
  DefineOutput(6, THnSparseF::Class());
  DefineOutput(7, THnSparseF::Class());
  DefineOutput(8, THnSparseF::Class());
  DefineOutput(9,TH1F::Class());
  DefineOutput(10,TH1F::Class());
  DefineOutput(11,TH1F::Class());
  DefineOutput(12,TH1F::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEmcCorr::~AliAnalysisTaskSEmcCorr(){
  // destructor
  delete fNentries;
  delete fhNprongsD0;
  delete fhNprongsD0chargedOnly;
  delete fhNprongsD0chargedRef;
  delete fhMCcorrel;
  delete fhMCtrigPart;
  delete fhMChadroncorrel;
  delete fhMChadrontrigPart;
  delete fArrayDaugh;
  delete fArrayAssoc;
  delete fhDzeroDecay;  
  delete fhDplusDecay;
  delete fhLambdaCDecay;
  delete fhAllBDecay;
}

//________________________________________________________________________
void AliAnalysisTaskSEmcCorr::Init()
{
  // Initialization
  fArrayDaugh=new TArrayI(20);

  fArrayAssoc=new TArrayI(200);

}
//_______________________________________________
void AliAnalysisTaskSEmcCorr::GetFinalStateParticles(AliAODMCParticle *part,TClonesArray *clarr){
  
  Int_t fdaugh=part->GetDaughterLabel(0);
  Int_t ldaugh=part->GetDaughterLabel(1);  
  if(fdaugh<0||ldaugh<0||(fdaugh>ldaugh))return;
  for(Int_t j=fdaugh;j<=ldaugh;j++){
    AliAODMCParticle *pd=(AliAODMCParticle*)clarr->At(j);
    Int_t pdg=pd->GetPdgCode();
    fnPhysPrim++;// increment it now and decrement it after all if, in case the daughter is not a phys prim
    fpxtotdaugh+=pd->Px(); // increment it now and decrement it after all if, in case the daughter is not a phys prim
    fpytotdaugh+=pd->Py();
    fpztotdaugh+=pd->Pz();
    
    if(pdg==211)fnpionPlus++;                         
    else if(pdg==-211)fnpionMinus++;                           
    else if(pdg==111)fnpionZero++;                         
    else if(pdg==321)fnkaonPlus++;                        
    else if(pdg==-321)fnkaonMinus++;                       
    else if(pdg==310)fnkaonZeroShort++;                        
    else if(pdg==2212)fnProton++;                        
    else if(pdg==-2212)fnAntiProton++;                        
    else if(pdg==-11)fnElePlus++;                        
    else if(pdg==11)fnEleMinus++;                        
    else if(pdg==-13)fnMuonPlus++;                        
    else if(pdg==13)fnMuonMinus++;                        
    else if(pdg==12||pdg==14||pdg==-12||pdg==-14)fnNeutrino++;
    else if(pdg==443)fnJpsi++;                                 
    else if(pdg!=130&&pdg!=22){
      fnPhysPrim--;
      fpxtotdaugh-=pd->Px();
      fpytotdaugh-=pd->Py();
      fpztotdaugh-=pd->Pz();
      GetFinalStateParticles(pd,clarr);
    }
  }
}

//_______________________________________________
void AliAnalysisTaskSEmcCorr::CheckDzeroChannels(AliAODMCParticle *part,TClonesArray *clarr){
  Int_t pdg=part->GetPdgCode();
  if(!(TMath::Abs(pdg)==421))return;
  fhDzeroDecay->Fill(0);
  fnpionPlus=0;                         
  fnpionMinus=0;                       
  fnpionZero=0;                         
  fnkaonPlus=0;                        
  fnkaonMinus=0;                       
  fnkaonZeroShort=0;                        
  fnProton=0;                        
  fnAntiProton=0;
  fnElePlus=0;                        
  fnEleMinus=0;                        
  fnMuonPlus=0;
  fnMuonMinus=0;
  fnNeutrino=0;
  fnJpsi=0;                            
  fnPhysPrim=0;
  fpxtotdaugh=0.;
  fpytotdaugh=0.; 
  fpztotdaugh=0.;

  GetFinalStateParticles(part,clarr);
  
  Double_t px=part->Px(); 
  Double_t py=part->Py(); 
  Double_t pz=part->Pz(); 
  // check exclusive channels in the list
  if(TMath::Abs(px-fpxtotdaugh)<1.e-6&&TMath::Abs(py-fpytotdaugh)<1.e-6&&TMath::Abs(pz-fpztotdaugh)<1.e-6){
    if(fnPhysPrim==2&&((pdg==421&&fnpionPlus==1&&fnkaonMinus==1)||(pdg==-421&&fnpionMinus==1&&fnkaonPlus==1))){
      fhDzeroDecay->Fill(6);
    }
    if(fnPhysPrim==3&&((pdg==421&&fnpionPlus==1&&fnkaonMinus==1&&fnpionZero==1)||(pdg==-421&&fnpionMinus==1&&fnkaonPlus==1&&fnpionZero==1))){
      fhDzeroDecay->Fill(7);
    }
    if(fnPhysPrim==4&&((pdg==421&&fnpionPlus==2&&fnkaonMinus==1&&fnpionMinus==1)||(pdg==-421&&fnpionMinus==2&&fnkaonPlus==1&&fnpionPlus==1))){
      fhDzeroDecay->Fill(8);
    }
    if(fnPhysPrim==3&&((pdg==421&&fnElePlus==1&&fnkaonMinus==1&&fnNeutrino==1)||(pdg==-421&&fnEleMinus==1&&fnkaonPlus==1&&fnNeutrino==1))){
      fhDzeroDecay->Fill(9);
    }
    if(fnPhysPrim==3&&((pdg==421&&fnMuonPlus==1&&fnkaonMinus==1&&fnNeutrino==1)||(pdg==-421&&fnMuonMinus==1&&fnkaonPlus==1&&fnNeutrino==1))){
      fhDzeroDecay->Fill(10);
    }
  }
  
  //now check decay to electron + X
  if(fnElePlus>=1||fnEleMinus>=1){
    fhDzeroDecay->Fill(11);
  }
  // no check number of prong
  Int_t nprong=fnPhysPrim;
  nprong-=fnpionZero;
  nprong-=fnNeutrino;
  nprong-=fnkaonZeroShort;
  if(nprong==0){
    fhDzeroDecay->Fill(1);
  }
  else if(nprong==1){// should never happen for D0
    fhDzeroDecay->Fill(2);
  }
  else if(nprong==2){// should never happen for D0
    fhDzeroDecay->Fill(3);
  }
  else if(nprong==3){
    fhDzeroDecay->Fill(4);
  }
  else if(nprong==4){
    fhDzeroDecay->Fill(5);
  }
  
  return;
}


//_______________________________________________
void AliAnalysisTaskSEmcCorr::CheckDplusChannels(AliAODMCParticle *part,TClonesArray *clarr){
  Int_t pdg=part->GetPdgCode();
  if(!(TMath::Abs(pdg)==411))return;
  fhDplusDecay->Fill(0);
  fnpionPlus=0;                         
  fnpionMinus=0;                       
  fnpionZero=0;                         
  fnkaonPlus=0;                        
  fnkaonMinus=0;                       
  fnkaonZeroShort=0;                        
  fnProton=0;                        
  fnAntiProton=0;
  fnElePlus=0;                        
  fnEleMinus=0;                        
  fnMuonPlus=0;
  fnMuonMinus=0;
  fnNeutrino=0;
  fnJpsi=0;                            
  fnPhysPrim=0;
  fpxtotdaugh=0.;
  fpytotdaugh=0.; 
  fpztotdaugh=0.;
  GetFinalStateParticles(part,clarr);
  
  Double_t px=part->Px(); 
  Double_t py=part->Py(); 
  Double_t pz=part->Pz(); 
  // check exclusive channels in the list
  if(TMath::Abs(px-fpxtotdaugh)<1.e-6&&TMath::Abs(py-fpytotdaugh)<1.e-6&&TMath::Abs(pz-fpztotdaugh)<1.e-6){
    if(fnPhysPrim==3&&((pdg==411&&fnpionPlus==2&&fnkaonMinus==1)||(pdg==-411&&fnpionMinus==2&&fnkaonPlus==1))){
      fhDplusDecay->Fill(6);
    }
    if(fnPhysPrim==3&&((pdg==411&&fnkaonZeroShort==1&&fnpionPlus==1&&fnpionZero==1)||(pdg==-411&&fnkaonZeroShort==1&&fnpionMinus==1&&fnpionZero==1))){
      fhDplusDecay->Fill(7);
    }
    if(fnPhysPrim==4&&((pdg==411&&fnpionPlus==2&&fnkaonMinus==1&&fnpionZero==1)||(pdg==-411&&fnpionMinus==2&&fnkaonPlus==1&&fnpionZero==1))){
      fhDzeroDecay->Fill(8);
    }
    if(fnPhysPrim==3&&((pdg==411&&fnElePlus==1&&fnkaonZeroShort==1&&fnNeutrino==1)||(pdg==-411&&fnEleMinus==1&&fnkaonZeroShort==1&&fnNeutrino==1))){
      fhDplusDecay->Fill(9);
    }
    if(fnPhysPrim==3&&((pdg==411&&fnMuonPlus==1&&fnkaonZeroShort==1&&fnNeutrino==1)||(pdg==-411&&fnMuonMinus==1&&fnkaonZeroShort==1&&fnNeutrino==1))){
      fhDplusDecay->Fill(10);
    }
  }
  
  //now check decay to electron + X
  if(fnElePlus>=1||fnEleMinus>=1){
    fhDplusDecay->Fill(11);
  }
  // no check number of prong
  Int_t nprong=fnPhysPrim;
  nprong-=fnpionZero;
  nprong-=fnNeutrino;
  nprong-=fnkaonZeroShort;
  if(nprong==0){
    fhDplusDecay->Fill(1);
  }
  else if(nprong==1){// should never happen for D0
    fhDplusDecay->Fill(2);
  }
  else if(nprong==2){// should never happen for D0
    fhDplusDecay->Fill(3);
  }
  else if(nprong==3){
    fhDplusDecay->Fill(4);
  }
  else if(nprong==4){
    fhDplusDecay->Fill(5);
  }
  
  return;
}
//________________________________________________________________________
// void AliAnalysisTaskSEmcCorr::CheckLambdaChannels(AliAODMCParticle *part,TClonesArray *clarr){
//   // not implemented yet
// }

//________________________________________________________________________
// void AliAnalysisTaskSEmcCorr::CheckBallChannels(AliAODMCParticle *part,TClonesArray *clarr){
//   // not implemented yet
// }


//________________________________________________________________________
void AliAnalysisTaskSEmcCorr::UserCreateOutputObjects(){
 
  if(!fArrayDaugh)  fArrayDaugh=new TArrayI(20);

  if(!fArrayAssoc)fArrayAssoc=new TArrayI(200);
  

  fNentries=new TH1F("nentriesChFr", "Analyzed sample properties", 21,0.5,21.5);
  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nEvSel");
  fNentries->GetXaxis()->SetBinLabel(3,"nEvPile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(4,"nEvGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"nEvRejVtxZ");
  fNentries->GetXaxis()->SetBinLabel(6,"nD0");
  fNentries->GetXaxis()->SetBinLabel(7,"nEventsSelGenName");

  fhNprongsD0=new TH1F("fhNprongsD0","fhNprongsD0",20,-0.5,19.5);
  fhNprongsD0chargedOnly=new TH1F("fhNprongsD0chargedOnly","fhNprongsD0chargedOnly",20,-0.5,19.5);
  fhNprongsD0chargedRef=new TH1F("fhNprongsD0chargedRef","fhNprongsD0chargedRef",20,-0.5,19.5);


  Int_t nbinsCorrMC[8]={15,20,20,20,50,63,20,11};
  Double_t binlowCorrMC[8]={-7.5,0.,-1.,0.,-5.,-0.5*TMath::Pi(),-2,-1.5};
  Double_t binupCorrMC[8]={7.5,20.,1.,2.,5.,1.5*TMath::Pi(),2.,9.5};
  fhMCcorrel=new THnSparseF("fhMCcorrel","fhMCcorrel;pdg;ptTrig;etaTrig;ptAss;etaAss;deltaPhi;deltaEta;pdgAss;",8,nbinsCorrMC,binlowCorrMC,binupCorrMC);

  Int_t nbinsTrigMC[3]={15,20,20};
  Double_t binlowTrigMC[3]={-7.5,0.,-1.};
  Double_t binupTrigMC[3]={7.5,20.,1.};
  fhMCtrigPart=new THnSparseF("fhMCtrigPart","fhMCcorrel;pdg;ptTrig;etaTrig;",3,nbinsTrigMC,binlowTrigMC,binupTrigMC);

  fhDzeroDecay=new TH1F("fhDzeroDecay","Dzero Decay channels",14,-0.5,13.5);
  fhDzeroDecay->GetXaxis()->SetBinLabel(1,"All D0");
  fhDzeroDecay->GetXaxis()->SetBinLabel(2,"0 prong");
  fhDzeroDecay->GetXaxis()->SetBinLabel(3,"1 prong");
  fhDzeroDecay->GetXaxis()->SetBinLabel(4,"2 prong");
  fhDzeroDecay->GetXaxis()->SetBinLabel(5,"3 prong");
  fhDzeroDecay->GetXaxis()->SetBinLabel(6,"4 prong");
  fhDzeroDecay->GetXaxis()->SetBinLabel(7,"K pi (3.88)");
  fhDzeroDecay->GetXaxis()->SetBinLabel(8,"K pi pi0 (13.9)");
  fhDzeroDecay->GetXaxis()->SetBinLabel(9,"K pi pi pi (8.09)");
  fhDzeroDecay->GetXaxis()->SetBinLabel(10,"e K nu (3.55)");
  fhDzeroDecay->GetXaxis()->SetBinLabel(11,"mu K nu (3.31)");
  fhDzeroDecay->GetXaxis()->SetBinLabel(12,"e X (6.49)");


  fhDplusDecay=new TH1F("fhDplusDecay","Dplus Decay channels",14,-0.5,13.5);
  fhDplusDecay->GetXaxis()->SetBinLabel(1,"All D+");
  fhDplusDecay->GetXaxis()->SetBinLabel(2,"0 prong");
  fhDplusDecay->GetXaxis()->SetBinLabel(3,"1 prong");
  fhDplusDecay->GetXaxis()->SetBinLabel(4,"2 prong");
  fhDplusDecay->GetXaxis()->SetBinLabel(5,"3 prong");
  fhDplusDecay->GetXaxis()->SetBinLabel(6,"4 prong");
  fhDplusDecay->GetXaxis()->SetBinLabel(7,"K pi pi (9.4)");
  fhDplusDecay->GetXaxis()->SetBinLabel(8,"K0s pi pi0 (6.9)");
  fhDplusDecay->GetXaxis()->SetBinLabel(9,"K pi pi pi0 (6.08)");
  fhDplusDecay->GetXaxis()->SetBinLabel(10,"e K0bar nu (8.83)");
  fhDplusDecay->GetXaxis()->SetBinLabel(11,"mu K0bar nu (9.4)");
  fhDplusDecay->GetXaxis()->SetBinLabel(12,"e X (16.07");

  fhLambdaCDecay=new TH1F("fhLambdaCDecay","LambdaC Decay channels",10,-0.5,9.5);
  fhLambdaCDecay->GetXaxis()->SetBinLabel(1,"All LambdaC");
  fhLambdaCDecay->GetXaxis()->SetBinLabel(2,"p K0s");
  fhLambdaCDecay->GetXaxis()->SetBinLabel(3,"p K pi");

  fhAllBDecay=new TH1F("fhAllBDecay","All bhadron (also LambdaB) Decay channels",10,-0.5,9.5);
  fhAllBDecay->GetXaxis()->SetBinLabel(1,"All b-hadrons");
  fhAllBDecay->GetXaxis()->SetBinLabel(2,"D0bar X (59.6)");
  fhAllBDecay->GetXaxis()->SetBinLabel(3,"D0bar D0 X (5.1) ");
  fhAllBDecay->GetXaxis()->SetBinLabel(4,"D+- D0 X (2.7)");
  fhAllBDecay->GetXaxis()->SetBinLabel(5,"D- X (22.7)");
  fhAllBDecay->GetXaxis()->SetBinLabel(6,"J/Psi X (1.16)");


  
  Int_t nbinsCorrMChadron[8]={6,20,20,20,50,63,20,7};
  Double_t binlowCorrMChadron[8]={-1.5,0.,-1.,0.,-5.,-0.5*TMath::Pi(),-2,-1.5};
  Double_t binupCorrMChadron[8]={4.5,20.,1.,2.,5.,1.5*TMath::Pi(),2.,5.5};
  fhMChadroncorrel=new THnSparseF("fhMChadroncorrel","fhMChadroncorrel;pdg;ptTrig;etaTrig;ptAss;etaAss;deltaPhi;deltaEta;pdgAss;",8,nbinsCorrMChadron,binlowCorrMChadron,binupCorrMChadron);

  Int_t nbinsTrigMChadron[3]={6,20,20};
  Double_t binlowTrigMChadron[3]={-1.5,0.,-1.};
  Double_t binupTrigMChadron[3]={4.5,20.,1.};
  fhMChadrontrigPart=new THnSparseF("fhMChadrontrigPart","fhMChadroncorrel;pdg;ptTrig;etaTrig;",3,nbinsTrigMChadron,binlowTrigMChadron,binupTrigMChadron);
  

  PostData(1,fNentries);
  PostData(2,fhNprongsD0);
  PostData(3,fhNprongsD0chargedOnly);
  PostData(4,fhNprongsD0chargedRef);
  PostData(5,fhMCcorrel);
  PostData(6,fhMCtrigPart);
  PostData(7,fhMChadroncorrel);
  PostData(8,fhMChadrontrigPart);
  PostData(9,fhDzeroDecay);
  PostData(10,fhDplusDecay);
  PostData(11,fhLambdaCDecay);
  PostData(12,fhAllBDecay);


}



//________________________________________________________________________
void AliAnalysisTaskSEmcCorr::UserExec(Option_t */*option*/){



  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!aod) {
    Printf("ERROR: aod not available");
    return;
  }



  fNentries->Fill(1);

  // AOD primary vertex
  
  
  if(!fOnlyKine){
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    TString primTitle = vtx1->GetTitle();
    if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) { 
      
      fNentries->Fill(3);
      
    }
    else {
      PostData(1,fNentries);
      return;
      
    }
  }
  

  TClonesArray *arrayMC=0x0;
  AliAODMCHeader *aodmcHeader=0x0;
  Double_t vtxTrue[3];
 
  
  if(fReadMC){
    //    aod->GetList()->Print();
    // load MC particles
    if(fFastSimul){
      arrayMC = 
	(TClonesArray*)aod->GetList()->At(23);//FindObject("mcparticles");
    }
    else{
      arrayMC = 
	(TClonesArray*)aod->GetList()->FindObject("mcparticles");
    }
    if(!arrayMC) {
      Printf("AliAnalysisTaskSEmcCorr::UserExec: MC particles branch not found!\n");
      fNentries->Fill(2);
      PostData(1,fNentries);
      return;
    }
    // load MC header
    aodmcHeader = 
      (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!aodmcHeader) {
      Printf("AliAnalysisTaskSEmcCorr::UserExec: MC header branch not found!\n");
      fNentries->Fill(3);
      PostData(1,fNentries);
      return;
    }
    
    if(!(fGeneratorString.IsNull())){
      TString strGenTitle=((AliGenEventHeader*)aodmcHeader->GetCocktailHeader(0))->GetTitle();
      // Printf("Gen Title: %s",strGenTitle.Data());
      if(!fGeneratorString.Contains(strGenTitle.Data())){
	//	Printf("Event rejected from generator title");
	return;
      }
    }
    fNentries->Fill(4);
    //    printf("N MC entries: %d \n",arrayMC->GetEntriesFast());
    // MC primary vertex
    aodmcHeader->GetVertex(vtxTrue);
    
    
    // FILL HISTOS FOR D0 mesons, c quarks and D0 from B properties
    if(fDoHFCorrelations)      SelectAssociatedParticles(arrayMC);
    for(Int_t jp=0;jp<arrayMC->GetEntriesFast();jp++){
      AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->At(jp);
      Int_t pdg=TMath::Abs(part->GetPdgCode());
      if(fDoHFCorrelations){	
	if(pdg==421||pdg==411||pdg==413){
	  
	  if(TMath::Abs(part->Y())>fYrangeTrig)continue;
	  if(TMath::Abs(part->Pt())<fminDpt)continue;
	  FillSkipParticleArray(part,arrayMC,jp);
	  FillCorrelationPlots(part,arrayMC);
	}
	else if(pdg==11){	  
	  if(TMath::Abs(part->Eta())>fEtarangeEleTrig)continue;
	  if(TMath::Abs(part->Pt())<fminDpt)continue;
	  FillSkipParticleArray(part,arrayMC,jp);
	  FillCorrelationPlots(part,arrayMC);	
	}       
      }
      if(fCheckDecay){
	CheckDzeroChannels(part,arrayMC);
	CheckDplusChannels(part,arrayMC);
	//	  CheckLambdaChannels(part,arrayMC);
	//	  CheckBallChannels(part,arrayMC);
	  
	// NOW STUDY OF N-PRONGS (OLD PART OF THE TASK)
	if(pdg==421){
	  Bool_t semilept=kFALSE;
	  Int_t nprongsD0=0;
	  Int_t nprongsD0charged=0;
	  
	  Int_t lab1=part->GetDaughterLabel(0);
	  Int_t lab2=part->GetDaughterLabel(1);
	  if(lab1>=0&&lab2>=0){
	    for(Int_t jd=lab1;jd<=lab2;jd++){
	      AliAODMCParticle *d=(AliAODMCParticle*)arrayMC->At(jd);
	      Int_t pdgd=TMath::Abs(d->GetPdgCode());
	      if(pdgd>=11&&pdgd<=16){
		semilept=kTRUE;	   
		break;
	      }
	      if(pdgd==211||pdgd==111||pdgd==321||pdgd==311||pdgd==310||pdgd==130){
		nprongsD0++;
	      }
	      if(pdgd==211||pdgd==321){
		nprongsD0charged++;
	      }	    
	      else{	    
		Int_t lab1d=d->GetDaughterLabel(0);
		Int_t lab2d=d->GetDaughterLabel(1);
		if(lab1d>=0&&lab2d>=0){
		  for(Int_t jdd=lab1d;jdd<=lab2d;jdd++){
		    AliAODMCParticle *dd=(AliAODMCParticle*)arrayMC->At(jdd);
		    Int_t pdgdd=TMath::Abs(dd->GetPdgCode());
		    if(pdgdd==211||pdgdd==111||pdgdd==321||pdgdd==311||pdgdd==310||pdgdd==130){
		      nprongsD0++;
		    }
		    if(pdgd==211||pdgd==321){
		      nprongsD0charged++;
		    } 
		  }
		}
	      }
	    }
	  }         
	  if(!semilept){
	    fhNprongsD0->Fill(nprongsD0);
	    fhNprongsD0chargedOnly->Fill(nprongsD0charged);
	    fhNprongsD0chargedRef->Fill(nprongsD0charged);
	  }
	}
      }
    }
    
    if(fDoHadronHadron){
      for(Int_t jp=0;jp<fLastAss;jp++){
	
	AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->At(fArrayAssoc->At(jp));
	Int_t pdg=TMath::Abs(part->GetPdgCode());      
	if(pdg==11||pdg==13||pdg==211||pdg==2212||pdg==321){
	  if(TMath::Abs(part->Eta())>fEtarangeEleTrig)continue;
	  if(TMath::Abs(part->Pt())<fminDpt)continue;
	  FillSkipParticleArray(part,arrayMC,fArrayAssoc->At(jp));
	  FillCorrelationPlotsHadrons(part,arrayMC);
	}
      }
      
    }
    
  }
  

  PostData(1,fNentries);
  PostData(2,fhNprongsD0);
  PostData(3,fhNprongsD0chargedOnly);
  PostData(4,fhNprongsD0chargedRef);
  PostData(5,fhMCcorrel);
  PostData(6,fhMCtrigPart);
  PostData(7,fhMChadroncorrel);
  PostData(8,fhMChadrontrigPart);
  PostData(9,fhDzeroDecay);
  PostData(10,fhDplusDecay);
  PostData(11,fhLambdaCDecay);
  PostData(12,fhAllBDecay);
}



void AliAnalysisTaskSEmcCorr::FillSkipParticleArray(AliAODMCParticle *part,TClonesArray *arrayMC,Int_t partPos){
  fArrayDaugh->Reset(0);
  fArrayDaugh->AddAt(partPos,0);
  flastdaugh=1;
  
  // Loop on daughters and nephews
  if(part->GetNDaughters()>0){
    Int_t pos1=part->GetDaughterLabel(0);
    Int_t pos2=part->GetDaughterLabel(1);
    if(pos2<0)pos2=pos1;
    if(pos1>=0){
      for(Int_t jd=pos1;jd<=pos2;jd++){ 
	AliAODMCParticle *pd=(AliAODMCParticle*)arrayMC->At(jd);
	Int_t pdgd=TMath::Abs(pd->GetPdgCode());
	if(pdgd==11||pdgd==13||pdgd==211||pdgd==321||pdgd==2212){
	  fArrayDaugh->AddAt(jd,flastdaugh);
	  flastdaugh++;
	}
	else if(pdgd!=12&&pdgd!=14&&pdgd!=16&&pdgd!=310&&pdgd!=111){// CHECK NEPHEWS (FOR RESONANT CHANNELS)
	  Int_t nephw=pd->GetNDaughters();
	  if(nephw>0){
	    Int_t posN1=pd->GetDaughterLabel(0);
	    Int_t posN2=pd->GetDaughterLabel(1);
	    if(posN2<0)posN2=posN1;
	    if(posN1>=0){
	      for(Int_t jN=posN1;jN<=posN2;jN++){ 
		AliAODMCParticle *pN=(AliAODMCParticle*)arrayMC->At(jN);
		Int_t pdgN=TMath::Abs(pN->GetPdgCode());
		if(pdgN==11||pdgN==13||pdgN==211||pdgN==321||pdgN==2212){
		  fArrayDaugh->AddAt(jN,flastdaugh);
		  flastdaugh++;
		}
		else if(pdgN!=12&&pdgN!=14&&pdgN!=16&&pdgN!=310&&pdgN!=111){// CHECK one step more (FOR RESONANT CHANNELS)	 
		  
		  Int_t gnephw=pN->GetNDaughters();
		  if(gnephw>0){
		    Int_t posGN1=pN->GetDaughterLabel(0);
		    Int_t posGN2=pN->GetDaughterLabel(1);
		    if(posGN2<0)posGN2=posGN1;
		    if(posGN1>=0){
		      for(Int_t jGN=posGN1;jGN<=posGN2;jGN++){ 
			AliAODMCParticle *pGN=(AliAODMCParticle*)arrayMC->At(jGN);
			Int_t pdgGN=TMath::Abs(pGN->GetPdgCode());
			if(pdgGN==11||pdgGN==13||pdgGN==211||pdgGN==321||pdgGN==2212){
			  fArrayDaugh->AddAt(jGN,flastdaugh);
			  flastdaugh++;			 
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } 
  }
  //  printf("particles to be skipped %d\n",flastdaugh++);
}

void AliAnalysisTaskSEmcCorr::SelectAssociatedParticles(TClonesArray *arrayMC){
  fLastAss=0;
  fArrayAssoc->Reset(0);
  Double_t pt,eta;
  for(Int_t jAss=0;jAss<arrayMC->GetEntriesFast();jAss++){
    AliAODMCParticle *partAss=(AliAODMCParticle*)arrayMC->At(jAss);
    Int_t pdgAss=TMath::Abs(partAss->GetPdgCode());
    if(pdgAss==11||pdgAss==13||pdgAss==211||pdgAss==321||pdgAss==2212){// check type: keep only e,mu,pi,K+,protons
      pt=partAss->Pt();
      eta=TMath::Abs(partAss->Eta());
      if(pt<fminPtass)continue;
      if(eta>fmaxEtaAss)continue;
      if(!partAss->IsPhysicalPrimary()){
	if(pdgAss!=11)	continue;// keep ele
      }
      fArrayAssoc->AddAt(jAss,fLastAss);
      fLastAss++;
      if(fLastAss==fArrayAssoc->GetSize()-1){
	fArrayAssoc->Set(fArrayAssoc->GetSize()+200);
      }
    }
  }  
  //  printf("n ass part: %d\n",fLastAss++);
}

void AliAnalysisTaskSEmcCorr::FillCorrelationPlots(AliAODMCParticle *part,TClonesArray *arrayMC){
  

  Bool_t hasToSkip=kFALSE;  
  Int_t index=0,softpi=-1;
  Int_t pdgTrig=TMath::Abs(part->GetPdgCode());
  if(pdgTrig==421){
    pdgTrig=1;
  }
  else if(pdgTrig==411){
    pdgTrig=2;
  }
  else if(pdgTrig==413){
    pdgTrig=3;
  }
  else if(pdgTrig==11){
    
    pdgTrig=4;
  }
  else pdgTrig=0;
  // check if it is prompt or from B and, for electrons the origin
  Int_t labmum=part->GetMother();
  if(labmum>=0){
    AliAODMCParticle *moth=(AliAODMCParticle*)arrayMC->At(labmum);
    Int_t pdgmum=TMath::Abs(moth->GetPdgCode());
    if(pdgTrig==4){// for electrons go to the parent hadron
      //      labmum=part->GetMother();
      //      moth=(AliAODMCParticle*)arrayMC->At(labmum);
      //      pdgmum=TMath::Abs(moth->GetPdgCode());
      if(!(pdgmum==5||pdgmum==4||(400<pdgmum&&pdgmum<600)||(4000<pdgmum&&pdgmum<6000))){// NON HF electron
	if(pdgmum==22){// GAMMA CONV
	  pdgTrig=5;
	}
	else if(pdgmum==111||pdgmum==113||pdgmum==221||pdgmum==223){ //DALITZ DECAY
	  pdgTrig=6;
	}
	else pdgTrig=7;
      }
    }
    if(pdgmum==413){ // D from Dstar: check soft pion  and go to the grandmother
      if(pdgTrig==1){
	for(Int_t isp=moth->GetDaughterLabel(0);isp<=moth->GetDaughterLabel(1);isp++){
	  AliAODMCParticle *sfp=(AliAODMCParticle*)arrayMC->At(isp);
	  Int_t pdgsp=TMath::Abs(sfp->GetPdgCode());
	  if(pdgsp==211)softpi=isp;
	}
      }
      labmum=moth->GetMother();
      if(labmum>=0){
	moth=(AliAODMCParticle*)arrayMC->At(labmum);
	pdgmum=TMath::Abs(moth->GetPdgCode());
      }
      else pdgmum=-1;
    }
    else if(pdgmum==423){// D*0 -> go to the mother
      labmum=moth->GetMother();
      if(labmum>=0){
	moth=(AliAODMCParticle*)arrayMC->At(labmum);
	pdgmum=TMath::Abs(moth->GetPdgCode());
      }
      else pdgmum=-1;
    }
    if(pdgmum==5||(500<pdgmum&&pdgmum<600)||(5000<pdgmum&&pdgmum<6000)){
      pdgTrig*=-1;
    }
    else if(pdgmum==4){
      labmum=moth->GetMother();
      if(labmum>=0){
	moth=(AliAODMCParticle*)arrayMC->At(labmum);
	pdgmum=TMath::Abs(moth->GetPdgCode());
	if(pdgmum==5||(500<pdgmum&&pdgmum<600)||(5000<pdgmum&&pdgmum<6000)){
	  //	  Printf("c quark from b quark");
	  pdgTrig*=-1;
	}   
      }
    }
    else {
      //      printf("Charm not from charm or beauty \n");
    }
  }
  Double_t phiTrig=part->Phi();
  Double_t etaTrig=part->Eta();
  Double_t deltaPhi;
  Double_t point[8]={ static_cast<Double_t>(pdgTrig),part->Pt(),etaTrig,0,0,0,0,0};//pdg,ptTrig,etaTrig,ptAss,etaAss,deltaPhi,deltaEta,pdgAss
  fhMCtrigPart->Fill(point);
  for(Int_t jAss=0;jAss<fLastAss;jAss++){
    index=fArrayAssoc->At(jAss);
    hasToSkip=kFALSE;  
      for(Int_t jd=0;jd<flastdaugh;jd++){
	if(index==fArrayDaugh->At(jd)){
	  hasToSkip=kTRUE;  
	  break;
	}
      }
      if(hasToSkip)continue;
      AliAODMCParticle *partAss=(AliAODMCParticle*)arrayMC->At(index);
      point[3]=partAss->Pt();
      point[4]=partAss->Eta();
      deltaPhi=partAss->Phi()-phiTrig;
      if(deltaPhi<-0.5*TMath::Pi())deltaPhi+=2.*TMath::Pi();
      else if(deltaPhi>1.5*TMath::Pi())deltaPhi-=2.*TMath::Pi();
      point[5]=deltaPhi;
      point[6]=partAss->Eta()-etaTrig;      
      UInt_t pdgassmum=TMath::Abs(partAss->GetPdgCode());
      if(index==softpi){
	point[7]=-1;
      }
      else if(pdgassmum==11){
	if(partAss->IsPhysicalPrimary())point[7]=1;
	Int_t elemum=partAss->GetMother();
	if(elemum>0){
	  AliAODMCParticle *partEleMum=(AliAODMCParticle*)arrayMC->At(elemum);
	  UInt_t pdgelemum=TMath::Abs(partEleMum->GetPdgCode());
	  if(pdgelemum==111||pdgelemum==113||pdgelemum==221||pdgelemum==223){// DALITZ DECAY
	    point[7]=6;
	  }
	  else if(pdgelemum==22){// GAMMA CONV
	    point[7]=7;
	  }
	  else if(pdgelemum==4||pdgelemum==5||(pdgelemum>400&&pdgelemum<600)||(pdgelemum>4000&&pdgelemum<6000)){// HF e
	    point[7]=8;
	  }
	  else point[7]=9;
	}
      }
      else if(pdgassmum==13){
	point[7]=2;
      }
      else if(pdgassmum==211){
	point[7]=3;
      }
      else if(pdgassmum==321){
	point[7]=4;
      }
      else if(pdgassmum==2212){
	point[7]=5;
      }
      else point[7]=0;

      fhMCcorrel->Fill(point);
  }

}



//______________________________________________________________
void AliAnalysisTaskSEmcCorr::FillCorrelationPlotsHadrons(AliAODMCParticle *part,TClonesArray *arrayMC){
  

  Bool_t hasToSkip=kFALSE;  
  Int_t index=0;
  Int_t pdgTrig=TMath::Abs(part->GetPdgCode());
  if(pdgTrig==11){
    pdgTrig=0;
  }
  else if(pdgTrig==13){
    pdgTrig=1;
  }
  else if(pdgTrig==211){
    pdgTrig=2;
  }
  else if(pdgTrig==321){
    pdgTrig=3;
  }
  else if(pdgTrig==2212){
    pdgTrig=4;
  }
  else pdgTrig=-1;
  
  
  Double_t phiTrig=part->Phi();
  Double_t etaTrig=part->Eta();
  Double_t deltaPhi;
  Double_t point[8]={ static_cast<Double_t>(pdgTrig),part->Pt(),etaTrig,0,0,0,0,0};//pdg,ptTrig,etaTrig,ptAss,etaAss,deltaPhi,deltaEta,pdgAss
  fhMChadrontrigPart->Fill(point);
  for(Int_t jAss=0;jAss<fLastAss;jAss++){
    index=fArrayAssoc->At(jAss);
    hasToSkip=kFALSE;  
    for(Int_t jd=0;jd<flastdaugh;jd++){
      if(index==fArrayDaugh->At(jd)){
	hasToSkip=kTRUE;  
	break;
      }
    }
    if(hasToSkip)continue;
    AliAODMCParticle *partAss=(AliAODMCParticle*)arrayMC->At(index);
    point[3]=partAss->Pt();
    if(point[3]>point[1])continue;
    point[4]=partAss->Eta();
    deltaPhi=partAss->Phi()-phiTrig;
    if(deltaPhi<-0.5*TMath::Pi())deltaPhi+=2.*TMath::Pi();
    else if(deltaPhi>1.5*TMath::Pi())deltaPhi-=2.*TMath::Pi();
    point[5]=deltaPhi;
    point[6]=partAss->Eta()-etaTrig;
    point[7]=TMath::Abs(partAss->GetPdgCode());
    if(point[7]==11){
      point[7]=0;
    }
    else if(point[7]==13){
      point[7]=1;
    }
    else if(point[7]==211){
      point[7]=2;
    }
    else if(point[7]==321){
      point[7]=3;
    }
    else if(point[7]==2212){
      point[7]=4;
    }
    else point[7]=-1;
    
    fhMChadroncorrel->Fill(point);
  }
 
}

//_______________________________________________________________
void AliAnalysisTaskSEmcCorr::Terminate(const Option_t*){
  //TERMINATE METHOD: NOTHING TO DO



}

