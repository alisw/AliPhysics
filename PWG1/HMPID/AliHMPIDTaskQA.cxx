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

//==============================================================================
// AliHMPIDTaskQA - Class representing a quality check tool of HMPID 
// A set of histograms is created.
//==============================================================================

#ifndef AliHMPIDTaskQA_CXX
#define AliHMPIDTaskQA_CXX


#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliHMPIDTaskQA.h"

ClassImp(AliHMPIDTaskQA)

//__________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA() :
  fESD(0x0),fMC(0x0),fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpPesdPhmp(0x0),fHmpCkovPesd(0x0),fHmpCkovPhmp(0x0),
  fHmpMipCharge3cm(0x0),fHmpMipTrkDist(0x0),
  fHmpTrkFlags(0x0),
  fN1(6),
  fN2(8),
  fPionEff(0x0),
  fKaonEff(0x0),
  fProtEff(0x0),
  fPionTot(0x0),
  fKaonTot(0x0),
  fProtTot(0x0),
  fPionNot(0x0),
  fKaonNot(0x0),
  fProtNot(0x0),
  fPionCon(0x0),
  fKaonCon(0x0),
  fProtCon(0x0)
{
  //
  //Default ctor
  //
  for (Int_t i=0; i<7; i++){
    fHmpPhotons[i]    = 0x0;
    fHmpPhotP[i]      = 0x0;
    fHmpPhotSin2th[i] = 0x0;
    fHmpMipTrkDistPosX[i] = fHmpMipTrkDistNegX[i] = 0x0;
    fHmpMipTrkDistPosY[i] = fHmpMipTrkDistNegY[i] = 0x0;
    fHmpMipCharge[i] = 0x0;
  }
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fMC(0x0), fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpPesdPhmp(0x0),fHmpCkovPesd(0x0),fHmpCkovPhmp(0x0),
  fHmpMipCharge3cm(0x0),fHmpMipTrkDist(0x0),
  fHmpTrkFlags(0x0),
  fN1(6),
  fN2(8),
  fPionEff(0x0),
  fKaonEff(0x0),
  fProtEff(0x0),
  fPionTot(0x0),
  fKaonTot(0x0),
  fProtTot(0x0),
  fPionNot(0x0),
  fKaonNot(0x0),
  fProtNot(0x0),
  fPionCon(0x0),
  fKaonCon(0x0),
  fProtCon(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<7; i++){
    fHmpPhotons[i] = 0x0;
    fHmpPhotP[i]   = 0x0;
    fHmpPhotSin2th[i] = 0x0;
    fHmpMipTrkDistPosX[i] = fHmpMipTrkDistNegX[i] = 0x0;
    fHmpMipTrkDistPosY[i] = fHmpMipTrkDistNegY[i] = 0x0;
    fHmpMipCharge[i] = 0x0;
  }

  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliHMPIDTaskQA& AliHMPIDTaskQA::operator=(const AliHMPIDTaskQA& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fESD             = c.fESD;
    fMC              = c.fMC;
    fUseMC           = c.fUseMC;
    fHmpHistList     = c.fHmpHistList;
    fHmpPesdPhmp     = c.fHmpPesdPhmp;
    fHmpCkovPesd     = c.fHmpCkovPesd;
    fHmpCkovPhmp     = c.fHmpCkovPhmp;
    fHmpMipCharge3cm = c.fHmpMipCharge3cm;
    fHmpMipTrkDist   = c.fHmpMipTrkDist;
    fHmpTrkFlags     = c.fHmpTrkFlags;
    fN1              = c.fN1;
    fN2              = c.fN2;
    fPionEff         = c.fPionEff;
    fKaonEff         = c.fKaonEff;
    fProtEff         = c.fProtEff;
    fPionTot         = c.fPionTot;
    fKaonTot         = c.fKaonTot;
    fProtTot         = c.fProtTot;
    fPionNot         = c.fPionNot;
    fKaonNot         = c.fKaonNot;
    fProtNot         = c.fProtNot;
    fPionCon         = c.fPionCon;
    fKaonCon         = c.fKaonCon;
    fProtCon         = c.fProtCon;
    for (Int_t i=0; i<7; i++){
      fHmpPhotons[i] = c.fHmpPhotons[i];
      fHmpPhotP[i]   = c.fHmpPhotP[i];
      fHmpPhotSin2th[i] = c.fHmpPhotSin2th[i];
      fHmpMipTrkDistPosX[i] = c.fHmpMipTrkDistPosX[i];
      fHmpMipTrkDistNegX[i] = c.fHmpMipTrkDistNegX[i];
      fHmpMipTrkDistPosY[i] = c.fHmpMipTrkDistPosY[i];
      fHmpMipTrkDistNegY[i] = c.fHmpMipTrkDistNegY[i];
      fHmpMipCharge[i] = c.fHmpMipCharge[i];
    }
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const AliHMPIDTaskQA& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fMC(c.fMC),fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpPesdPhmp(c.fHmpPesdPhmp),fHmpCkovPesd(c.fHmpCkovPesd),fHmpCkovPhmp(c.fHmpCkovPhmp),
  fHmpMipCharge3cm(c.fHmpMipCharge3cm),fHmpMipTrkDist(c.fHmpMipTrkDist),
  fHmpTrkFlags(c.fHmpTrkFlags),
  fN1(c.fN1),
  fN2(c.fN2),
  fPionEff(c.fPionEff),
  fKaonEff(c.fKaonEff),
  fProtEff(c.fProtEff),
  fPionTot(c.fPionTot),
  fKaonTot(c.fKaonTot),
  fProtTot(c.fProtTot),
  fPionNot(c.fPionNot),
  fKaonNot(c.fKaonNot),
  fProtNot(c.fProtNot),
  fPionCon(c.fPionCon),
  fKaonCon(c.fKaonCon),
  fProtCon(c.fProtCon)
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<7; i++){
    fHmpPhotons[i] = c.fHmpPhotons[i];
    fHmpPhotP[i]   = c.fHmpPhotP[i];
    fHmpPhotSin2th[i] = c.fHmpPhotSin2th[i];
    fHmpMipTrkDistPosX[i] = c.fHmpMipTrkDistPosX[i];
    fHmpMipTrkDistNegX[i] = c.fHmpMipTrkDistNegX[i];
    fHmpMipTrkDistPosY[i] = c.fHmpMipTrkDistPosY[i];
    fHmpMipTrkDistNegY[i] = c.fHmpMipTrkDistNegY[i];
    fHmpMipCharge[i] = c.fHmpMipCharge[i];
  }
}
 
//___________________________________________________________________________
AliHMPIDTaskQA::~AliHMPIDTaskQA() {
  //
  //destructor
  //
  Info("~AliHMPIDTaskQA","Calling Destructor");
  if (fHmpHistList /*&& !AliAnalysisManager::GetAnalysisManager()->IsProofMode()*/) delete fHmpHistList;
}

//___________________________________________________________________________
void AliHMPIDTaskQA::ConnectInputData(Option_t *option)
{
  AliAnalysisTaskSE::ConnectInputData(option);

  // Connect ESD here
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
  } else
    fESD = esdH->GetEvent();

  if (fUseMC){
    // Connect MC
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcH) {
      AliDebug(2,Form("ERROR: Could not get MCEventHandler"));
      fUseMC = kFALSE;
    } else
      fMC = mcH->MCEvent();
      if (!fMC) AliDebug(2,Form("ERROR: Could not get MCEvent"));
  }
}

//___________________________________________________________________________
void AliHMPIDTaskQA::UserExec(Option_t *)
{
  Double_t priors[10]={1.,1.,1.,1.,1.,0.,0.,0.,0.,0.}; //{0.01,0.01,0.83,0.10,0.5};
  Double_t probs[5];
  AliPID *pPid = new AliPID();
  pPid->SetPriors(priors);
  AliESDtrack *track=0;
  TParticle *pPart=0;
  AliStack* pStack = 0;
  Int_t label = -1;
  if (fUseMC){
    pStack = fMC->Stack();
  }

  //
  // Main loop function, executed on Event basis
  //
  
  if( !fESD->GetPrimaryVertex()) return;
  if( fESD->GetPrimaryVertex()->GetNContributors() < 1 ) return;
  if( TMath::Abs(fESD->GetPrimaryVertex()->GetZ()) > 10.0 /* cm */) return;
  
  
  for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {

    track = fESD->GetTrack(iTrack);
    if(!track) continue;
    if(! track->IsOn(AliESDtrack::kITSrefit) || !track->IsOn(AliESDtrack::kITSrefit) ) continue;
     
    Double_t rin[3], rout[3];
    track->GetInnerXYZ(rin);
    track->GetOuterXYZ(rout);
    Double_t ktol = 0.001;
    Double_t thetaCh = track->GetHMPIDsignal();

    if(Equal(thetaCh,-20.,ktol))      fHmpTrkFlags->Fill(0);
    else if(Equal(thetaCh,-9.,ktol))  fHmpTrkFlags->Fill(1);
    else if(Equal(thetaCh,-5.,ktol))  fHmpTrkFlags->Fill(2);
    else if(Equal(thetaCh,-11.,ktol)) fHmpTrkFlags->Fill(3);
    else if(Equal(thetaCh,-22.,ktol)) fHmpTrkFlags->Fill(4);
    else fHmpTrkFlags->Fill(4);

    if(Equal(thetaCh,-20.,ktol)) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;

    Int_t q, nph;
    Float_t x, y;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);

    if(Equal(x,0.,ktol) && Equal(y,0.,ktol) && Equal(xpc,0.,ktol) && Equal(ypc,0.,ktol)) continue;

    Float_t sign = track->GetSign();
    Int_t ch = track->GetHMPIDcluIdx()/1000000; 
    if( ch < 0 || ch > 6 ) continue; //Chamber numbering is [0,6]
    
    if (sign > 0.) fHmpMipTrkDistPosX[ch]->Fill(xpc - x), fHmpMipTrkDistPosY[ch]->Fill(ypc - y);
    if (sign < 0.) fHmpMipTrkDistNegX[ch]->Fill(xpc - x), fHmpMipTrkDistNegY[ch]->Fill(ypc - y);
    fHmpMipCharge[ch]->Fill(q);
    Double_t dist = TMath::Sqrt((xpc - x)*(xpc - x) + (ypc - y)*(ypc - y));
    if (dist <= 3.0) fHmpMipCharge3cm->Fill(q);
    fHmpMipTrkDist->Fill(dist);

    track->GetHMPIDpid(probs);
    pPid->SetProbabilities(probs);
    if (fUseMC){
      if ((label = track->GetLabel()) < 0) continue;
      pPart = pStack->Particle(label);
    }

    if(thetaCh > 0 ){
      Double_t pHmp[3] = {0}, pHmp3 = 0;
      if (track->GetOuterHmpPxPyPz(pHmp)) pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);
      if (pHmp3) fHmpPesdPhmp->Fill(track->P(),pHmp3);
      fHmpCkovPesd->Fill(track->P(),thetaCh);
      if (pHmp3) fHmpCkovPhmp->Fill(pHmp3,thetaCh);
      fHmpPhotons[ch]->Fill(nph);
      fHmpPhotP[ch]->Fill(pHmp3,nph);
      Double_t sin2 = TMath::Power(TMath::Sin(th),2);
      fHmpPhotSin2th[ch]->Fill(sin2,nph);

      if (fUseMC && dist<0.5 && TMath::Abs(th)<0.13){
        if (!pStack->IsPhysicalPrimary(label)) continue;
        if (track->Pt()<1. || track->Pt()>5.) continue;
        Int_t pdgCode = TMath::Abs(pPart->GetPdgCode());
        Int_t ptBin=(Int_t) (2*(track->Pt()-1));
        if (pdgCode!=2212) fProtCon->Fill(ptBin);
        if (pdgCode==2212){
          fProtTot->Fill(ptBin);
          fProtEff->Fill(ptBin,pPid->GetProbability(AliPID::kProton));
          fPionNot->Fill(ptBin,pPid->GetProbability(AliPID::kPion));
          fKaonNot->Fill(ptBin,pPid->GetProbability(AliPID::kKaon));
        }
        if (pdgCode!=211) fPionCon->Fill(ptBin);
        if (pdgCode!=321) fKaonCon->Fill(ptBin);
        if (pdgCode==211){
          if (ptBin < 6){
            Float_t weight=pPid->GetProbability(AliPID::kElectron)+
                           pPid->GetProbability(AliPID::kMuon)+
                           pPid->GetProbability(AliPID::kPion);
            fPionTot->Fill(ptBin);
            fPionEff->Fill(ptBin,weight);
            fKaonNot->Fill(ptBin,pPid->GetProbability(AliPID::kKaon));
          }
          fProtNot->Fill(ptBin,pPid->GetProbability(AliPID::kProton));
        }
        if (pdgCode==321){
          if (ptBin < 6){
            fKaonTot->Fill(ptBin);
            fKaonEff->Fill(ptBin,pPid->GetProbability(AliPID::kKaon));
            fPionNot->Fill(ptBin,(pPid->GetProbability(AliPID::kPion)));
          }
          fProtNot->Fill(ptBin,(pPid->GetProbability(AliPID::kProton)));
        }
      }
    }//there is signal
  }//track loop
  delete pPid;

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
}

//___________________________________________________________________________
void AliHMPIDTaskQA::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate"," ");
  AliAnalysisTaskSE::Terminate();

}

//___________________________________________________________________________
void AliHMPIDTaskQA::UserCreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //

  //slot #1
//   OpenFile(1);
   fHmpHistList = new TList();
   fHmpHistList->SetOwner();

   fHmpPesdPhmp = new TH2F("fHmpPesdPhmp","HMPID: ESD p - running p;HMP p (GeV/c);ESD p (Gev/c)",100,0,10,100,0,10);
   fHmpHistList->Add(fHmpPesdPhmp);

   fHmpCkovPesd = new TH2F("fHmpCkovPesd","HMPID: ThetaCherenkov vs P;p_esd (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPesd);

   fHmpCkovPhmp = new TH2F("fHmpCkovPhmp","HMPID: ThetaCherenkov vs P;p_hmp (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPhmp);

   for (Int_t i=0; i<7; i++){
     TString title=Form("MIP-Track distance in local X, Chamber %d;distance (cm);Entries",i);
     fHmpMipTrkDistPosX[i] = new TH1F(Form("fHmpMipTrkDistPosX%d",i),title.Data(),800,-20,20);
     fHmpHistList->Add(fHmpMipTrkDistPosX[i]);
     fHmpMipTrkDistNegX[i] = new TH1F(Form("fHmpMipTrkDistNegX%d",i),title.Data(),800,-20,20);
     fHmpHistList->Add(fHmpMipTrkDistNegX[i]);

     title=Form("MIP-Track distance in local Y, Chamber %d;distance (cm);Entries",i);
     fHmpMipTrkDistPosY[i] = new TH1F(Form("fHmpMipTrkDistPosY%d",i),title.Data(),800,-20,20);
     fHmpHistList->Add(fHmpMipTrkDistPosY[i]);
     fHmpMipTrkDistNegY[i] = new TH1F(Form("fHmpMipTrkDistNegY%d",i),title.Data(),800,-20,20);
     fHmpHistList->Add(fHmpMipTrkDistNegY[i]);

     title=Form("Mip charge distribution, Chamber %d;Charge;Entries",i);
     fHmpMipCharge[i] = new TH1F(Form("fHmpMipCharge%d",i),title.Data(),5001,-0.5,5000.5);
     fHmpHistList->Add(fHmpMipCharge[i]);

     title=Form("Photons, Chamber %d;N of Photons;Entries",i);
     fHmpPhotons[i] = new TH1F(Form("fHmpPhotons%d",i),title.Data(),100,0,100);
     fHmpHistList->Add(fHmpPhotons[i]);

     title=Form("Photons versus HMP momentum, Chamber %d;P (GeV);N of Photons",i);
     fHmpPhotP[i] = new TH2F(Form("fHmpPhotP%d",i),title.Data(),200,0.,20.,100,0,100);
     fHmpHistList->Add(fHmpPhotP[i]);

     title=Form("Photons versus sin(th)^2 (uncorrected), Chamber %d;P (GeV);N of Photons",i);
     fHmpPhotSin2th[i] = new TH2F(Form("fHmpPhotSin2th%d",i),title.Data(),100,0.,1.,100,0,100);
     fHmpHistList->Add(fHmpPhotSin2th[i]);
   }

   fHmpMipCharge3cm = new TH1F("fHmpMipCharge3cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge3cm);

   fHmpMipTrkDist = new TH1F("fHmpMipTrkDist","Mip-track distance in all the chambers",100,0.,10.);
   fHmpHistList->Add(fHmpMipTrkDist);

   fHmpTrkFlags = new TH1F("fHmpTrkFlags","HMPID track flags",6,0,6);
   TString summary[6] =  {"NotPerformed","MipDistCut", "MipQdcCut", "NoPhotAccept", "kNoRad", "other"};
   for(Int_t ibin = 0; ibin < 6; ibin++) fHmpTrkFlags->GetXaxis()->SetBinLabel(ibin+1,Form("%i  %s",ibin+1,summary[ibin].Data()));
   fHmpHistList->Add(fHmpTrkFlags);

   fPionEff = new TH1F("PionEff","Identified pions",fN1,0,fN1);
   fKaonEff = new TH1F("KaonEff","Identified kaons",fN1,0,fN1);
   fProtEff = new TH1F("ProtEff","Identified protons",fN2,0,fN2);
   fPionTot = new TH1I("PionTot","Total MC pions",fN1,0,fN1);
   fKaonTot = new TH1I("KaonTot","Total MC kaons",fN1,0,fN1);
   fProtTot = new TH1I("ProtTot","Total MC protons",fN2,0,fN2);
   fPionNot = new TH1F("PionNot","Misidentified pions",fN1,0,fN1);
   fKaonNot = new TH1F("KaonNot","Misidentified kaons",fN1,0,fN1);
   fProtNot = new TH1F("ProtNot","Misidentified protons",fN2,0,fN2);
   fPionCon = new TH1I("PionCon","Total not MC pions",fN1,0,fN1);
   fKaonCon = new TH1I("KaonCon","Total not MC kaons",fN1,0,fN1);
   fProtCon = new TH1I("ProtCon","Total not MC protons",fN2,0,fN2);

   fHmpHistList->Add(fPionEff); fHmpHistList->Add(fKaonEff); fHmpHistList->Add(fProtEff);
   fHmpHistList->Add(fPionTot); fHmpHistList->Add(fKaonTot); fHmpHistList->Add(fProtTot);
   fHmpHistList->Add(fPionNot); fHmpHistList->Add(fKaonNot); fHmpHistList->Add(fProtNot);
   fHmpHistList->Add(fPionCon); fHmpHistList->Add(fKaonCon); fHmpHistList->Add(fProtCon);

   PostData(1,fHmpHistList);
}

//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDTaskQA::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
   
#endif
