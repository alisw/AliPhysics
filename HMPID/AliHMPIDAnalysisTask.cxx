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
// AliHMPIDAnalysysTask - Class representing a basic analysis tool of HMPID data at
// level of ESD.
// A set of histograms is created.
//==============================================================================
//
// By means of AliHMPIDAnalysisTask.C macro it is possible to use this class
// to perform the analysis on local data, on data on alien using local machine
// and on CAF.

#ifndef AliHMPIDAnalysisTASK_CXX
#define AliHMPIDAnalysisTASK_CXX


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
#include "AliHMPIDAnalysisTask.h"

ClassImp(AliHMPIDAnalysisTask)

//__________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask() :
  fESD(0x0),fHmpHistList(0x0),
  fNevts(0),
  fTrigNevts(0),
  fTrigger(0),
  fHmpPesdPhmp(0x0),fHmpCkovPesd(0x0),fHmpCkovPhmp(0x0),
  fHmpMipTrkDist(0x0),fHmpMipTrkDistX(0x0),fHmpMipTrkDistY(0x0),fHmpMipCharge3cm(0x0),fHmpMipCharge1cm(0x0),fHmpNumPhots(0x0),
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
  fProtCon(0x0),
  fThetavsPiFromK(0x0),
  fThetapivsPesd(0x0),
  fThetaKvsPesd(0x0),
  fThetaPvsPesd(0x0)
{
  //
  //Default ctor
  //
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0), fHmpHistList(0x0), fNevts(0),
  fTrigNevts(0),
  fTrigger(0),
  fHmpPesdPhmp(0x0), fHmpCkovPesd(0x0), fHmpCkovPhmp(0x0),
  fHmpMipTrkDist(0x0), fHmpMipTrkDistX(0x0), fHmpMipTrkDistY(0x0),
  fHmpMipCharge3cm(0x0), fHmpMipCharge1cm(0x0),
  fHmpNumPhots(0x0), fHmpTrkFlags(0x0),
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
  fProtCon(0x0),
  fThetavsPiFromK(0x0),
  fThetapivsPesd(0x0),
  fThetaKvsPesd(0x0),
  fThetaPvsPesd(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //

  DefineOutput(0,TList::Class());
}

//___________________________________________________________________________
AliHMPIDAnalysisTask& AliHMPIDAnalysisTask::operator=(const AliHMPIDAnalysisTask& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fESD             = c.fESD;
    fHmpHistList     = c.fHmpHistList;
    fNevts           = c.fNevts;
    fTrigNevts       = c.fTrigNevts;
    fTrigger         = c.fTrigger;
    fHmpPesdPhmp     = c.fHmpPesdPhmp;
    fHmpCkovPesd     = c.fHmpCkovPesd;
    fHmpCkovPhmp     = c.fHmpCkovPhmp;
    fHmpMipTrkDist   = c.fHmpMipTrkDist;
    fHmpMipTrkDistX  = c.fHmpMipTrkDistX;
    fHmpMipTrkDistY  = c.fHmpMipTrkDistY;
    fHmpMipCharge3cm = c.fHmpMipCharge3cm;
    fHmpMipCharge1cm = c.fHmpMipCharge1cm;
    fHmpNumPhots     = c.fHmpNumPhots;
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
    fThetavsPiFromK  = c.fThetavsPiFromK;
    fThetapivsPesd   = c.fThetapivsPesd;
    fThetaKvsPesd    = c.fThetaKvsPesd;
    fThetaPvsPesd    = c.fThetaPvsPesd;
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fHmpHistList(c.fHmpHistList), fNevts(c.fNevts),
  fTrigNevts(c.fTrigNevts),
  fTrigger(c.fTrigger),
  fHmpPesdPhmp(c.fHmpPesdPhmp),fHmpCkovPesd(c.fHmpCkovPesd),fHmpCkovPhmp(c.fHmpCkovPhmp),
  fHmpMipTrkDist(c.fHmpMipTrkDist),fHmpMipTrkDistX(c.fHmpMipTrkDistX),fHmpMipTrkDistY(c.fHmpMipTrkDistY),
  fHmpMipCharge3cm(c.fHmpMipCharge3cm), fHmpMipCharge1cm(c.fHmpMipCharge1cm),
  fHmpNumPhots(c.fHmpNumPhots), fHmpTrkFlags(c.fHmpTrkFlags),
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
  fProtCon(c.fProtCon),
  fThetavsPiFromK(c.fThetavsPiFromK),
  fThetapivsPesd(c.fThetapivsPesd),
  fThetaKvsPesd(c.fThetaKvsPesd),
  fThetaPvsPesd(c.fThetaPvsPesd)
{
  //
  // Copy Constructor
  //
}
 
//___________________________________________________________________________
AliHMPIDAnalysisTask::~AliHMPIDAnalysisTask() {
  //
  //destructor
  //
  Info("~AliHMPIDAnalysisTask","Calling Destructor");
  if (fHmpHistList) {fHmpHistList->Clear(); delete fHmpHistList;}
}

//___________________________________________________________________________
void AliHMPIDAnalysisTask::ConnectInputData(Option_t *)
{
  // Connect ESD here

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
    } else
      fESD = esdH->GetEvent();

  // Connect MC

  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcH) {
      AliDebug(2,Form("ERROR: Could not get MCEventHandler"));
    } else
      fMC = mcH->MCEvent();
}

//_________________________________________________
void AliHMPIDAnalysisTask::Exec(Option_t *)
{
  AliStack* pStack = fMC->Stack();
  Int_t label;
  Double_t priors[5]={1.,1.,1.,1.,1.}; //{0.01,0.01,0.83,0.10,0.5};
  Double_t probs[5];
  AliPID *pPid = new AliPID();
  pPid->SetPriors(priors);
  Double_t n = 1.293;
  Double_t dGeVMass[] = {0.000511,0.105658,0.13957018,0.493677,0.938272};

  AliESDtrack *track=0;
  TParticle *pPart=0;
  //
  // Main loop function, executed on Event basis
  //
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    track = fESD->GetTrack(iTrack);
    if(!track) continue;
    Double_t dPtr = track->P();
    Double_t rin[3], rout[3];
    track->GetInnerXYZ(rin);
    track->GetOuterXYZ(rout);
    Double_t ktol = 0.001;

    if(Equal(track->GetHMPIDsignal(),-20.,ktol))      fHmpTrkFlags->Fill(0);
    else if(Equal(track->GetHMPIDsignal(),-9.,ktol))  fHmpTrkFlags->Fill(1);
    else if(Equal(track->GetHMPIDsignal(),-5.,ktol))   fHmpTrkFlags->Fill(2);
    else if(Equal(track->GetHMPIDsignal(),-11.,ktol))  fHmpTrkFlags->Fill(3);
    else if(Equal(track->GetHMPIDsignal(),-22.,ktol)) fHmpTrkFlags->Fill(4);
    else fHmpTrkFlags->Fill(4);

    if(Equal(track->GetHMPIDsignal(),-20.,ktol)) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;

    //Int_t ch = track->GetHMPIDcluIdx()/1000000;
    Int_t q, nph;
    Float_t x, y;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);

    if(Equal(x,0.,ktol) && Equal(y,0.,ktol) && Equal(xpc,0.,ktol) && Equal(ypc,0.,ktol)) continue;

    Double_t dist = TMath::Sqrt( (xpc-x)*(xpc-x) + (ypc - y)*(ypc - y));
    fHmpMipTrkDist->Fill(dist);
    fHmpMipTrkDistX->Fill(xpc - x);
    fHmpMipTrkDistY->Fill(ypc - y);
    Double_t pHmp[3] = {0}, pHmp3 = 0;
    if (track->GetOuterHmpPxPyPz(pHmp)) pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);
    if (dist <= 3.0) fHmpMipCharge3cm->Fill(q);

    Float_t thetaTheorPion = 999., thetaTheorKaon = 999., thetaTheorProt = 999.,
            thetaHmpTheorPion = 999., thetaHmpTheorKaon = 999., thetaHmpTheorProt = 999.;
    if(dPtr != 0){
      thetaTheorPion = TMath::ACos(TMath::Sqrt(dPtr*dPtr + dGeVMass[2]*dGeVMass[2])/(n*dPtr));
      thetaTheorKaon = TMath::ACos(TMath::Sqrt(dPtr*dPtr + dGeVMass[3]*dGeVMass[3])/(n*dPtr));
      thetaTheorProt = TMath::ACos(TMath::Sqrt(dPtr*dPtr + dGeVMass[4]*dGeVMass[4])/(n*dPtr));
    }
    if(pHmp3 != 0){
      thetaHmpTheorPion = TMath::ACos(TMath::Sqrt(pHmp3*pHmp3 + dGeVMass[2]*dGeVMass[2])/(n*pHmp3));
      thetaHmpTheorKaon = TMath::ACos(TMath::Sqrt(pHmp3*pHmp3 + dGeVMass[3]*dGeVMass[3])/(n*pHmp3));
      thetaHmpTheorProt = TMath::ACos(TMath::Sqrt(pHmp3*pHmp3 + dGeVMass[4]*dGeVMass[4])/(n*pHmp3));
    }

    track->GetHMPIDpid(probs);
    pPid->SetProbabilities(probs);
    if ((label = track->GetLabel()) < 0) continue;
    pPart = pStack->Particle(label);


    if(track->GetHMPIDsignal() > 0 )
    {
      fHmpPesdPhmp->Fill(track->P(),pHmp3);
      if (dist<=1.0) fHmpMipCharge1cm->Fill(q);
      fHmpNumPhots->Fill(nph);
      fHmpCkovPesd->Fill(track->P(),track->GetHMPIDsignal());
      fHmpCkovPhmp->Fill(pHmp3,track->GetHMPIDsignal());

      if (!pStack->IsPhysicalPrimary(label)) continue;
      Int_t pdgCode = TMath::Abs(pPart->GetPdgCode());
      if (pdgCode==211){
        fThetapivsPesd->Fill(track->P(),track->GetHMPIDsignal());
        Int_t mot=pPart->GetFirstMother();
        if (mot > -1){
          TParticle *pMot=pStack->Particle(mot);
          TString str=pMot->GetName();
          if (str.Contains("K")) fThetavsPiFromK->Fill(pHmp3,track->GetHMPIDsignal());
        }
      }
      if (pdgCode==321) fThetaKvsPesd->Fill(track->P(),track->GetHMPIDsignal());
      if (pdgCode==2212) fThetaPvsPesd->Fill(track->P(),track->GetHMPIDsignal());

      if (track->Pt()<1. || track->Pt()>5.) continue;
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
    }//there is signal
  }//track loop

  delete pPid;

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(0,fHmpHistList);
}

//___________________________________________________________________________
void AliHMPIDAnalysisTask::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  fHmpHistList = dynamic_cast<TList*> (GetOutputData(0));

  fPionEff = dynamic_cast<TH1F*> (fHmpHistList->FindObject("PionEff"));
  fPionTot = dynamic_cast<TH1I*> (fHmpHistList->FindObject("PionTot"));
  fPionNot = dynamic_cast<TH1F*> (fHmpHistList->FindObject("PionNot"));
  fPionCon = dynamic_cast<TH1I*> (fHmpHistList->FindObject("PionCon"));
  fKaonEff = dynamic_cast<TH1F*> (fHmpHistList->FindObject("KaonEff"));
  fKaonTot = dynamic_cast<TH1I*> (fHmpHistList->FindObject("KaonTot"));
  fKaonNot = dynamic_cast<TH1F*> (fHmpHistList->FindObject("KaonNot"));
  fKaonCon = dynamic_cast<TH1I*> (fHmpHistList->FindObject("KaonCon"));
  fProtEff = dynamic_cast<TH1F*> (fHmpHistList->FindObject("ProtEff"));
  fProtTot = dynamic_cast<TH1I*> (fHmpHistList->FindObject("ProtTot"));
  fProtNot = dynamic_cast<TH1F*> (fHmpHistList->FindObject("ProtNot"));
  fProtCon = dynamic_cast<TH1I*> (fHmpHistList->FindObject("ProtCon"));

  Float_t *pionEff=fPionEff->GetArray();
  Int_t   *pionTot=fPionTot->GetArray();
  Float_t *pionNot=fPionNot->GetArray();
  Int_t   *pionCon=fPionCon->GetArray();
  Float_t *kaonEff=fKaonEff->GetArray();
  Int_t   *kaonTot=fKaonTot->GetArray();
  Float_t *kaonNot=fKaonNot->GetArray();
  Int_t   *kaonCon=fKaonCon->GetArray();
  Float_t *protEff=fProtEff->GetArray();
  Int_t   *protTot=fProtTot->GetArray();
  Float_t *protNot=fProtNot->GetArray();
  Int_t   *protCon=fProtCon->GetArray();

  TGraphErrors *effPi = new TGraphErrors(fN1);
  TGraphErrors *effKa = new TGraphErrors(fN1);
  TGraphErrors *effPr = new TGraphErrors(fN2);
  TGraphErrors *conPi = new TGraphErrors(fN1);
  TGraphErrors *conKa = new TGraphErrors(fN1);
  TGraphErrors *conPr = new TGraphErrors(fN2);

  Float_t pt[8]={1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75};
  Float_t eff=0, effErr=0, con=0, conErr=0;
  for (Int_t i=0; i< fN2; i++){
    eff = protEff[i+1]/TMath::Max(protTot[i+1],1);
    effErr = TMath::Sqrt(eff*(1.-eff)/TMath::Max(protTot[i+1],1));
    con = protNot[i+1]/TMath::Max(protCon[i+1],1);
    conErr = TMath::Sqrt(con*(1.-con)/protCon[i+1]);
    effPr->SetPoint(i,pt[i],eff);
    effPr->SetPointError(i,0,effErr);
    conPr->SetPoint(i,pt[i],con);
    conPr->SetPointError(i,0,conErr);

    if (i>=fN1) continue;
    eff = pionEff[i+1]/pionTot[i+1];
    effErr = TMath::Sqrt(eff*(1.-eff)/pionTot[i+1]);
    con=pionNot[i+1]/pionCon[i+1];
    conErr = TMath::Sqrt(con*(1.-con)/pionCon[i+1]);
    effPi->SetPoint(i,pt[i],(Float_t)pionEff[i+1]/(Float_t)pionTot[i+1]);
    effPi->SetPointError(i,0,effErr);
    conPi->SetPoint(i,pt[i],(Float_t)pionNot[i+1]/(Float_t)pionCon[i+1]);
    conPi->SetPointError(i,0,conErr);

    eff = kaonEff[i+1]/TMath::Max(kaonTot[i+1],1);
    effErr = TMath::Sqrt(eff*(1.-eff)/kaonTot[i+1]);
    con = kaonNot[i+1]/TMath::Max(kaonCon[i+1],1);
    conErr = TMath::Sqrt(con*(1.-con)/kaonCon[i+1]);
    effKa->SetPoint(i,pt[i],(Float_t)kaonEff[i+1]/TMath::Max(kaonTot[i+1],1));
    effKa->SetPointError(i,0,effErr);
    conKa->SetPoint(i,pt[i],(Float_t)kaonNot[i+1]/TMath::Max(kaonCon[i+1],1));
    conKa->SetPointError(i,0,conErr);
  }

  TCanvas *pCan=new TCanvas("Hmp","Efficiency and contamination",500,900);
  pCan->Divide(1,3);

  pCan->cd(1);
  effPi->SetTitle("Pions");
  effPi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  effPi->GetYaxis()->SetRangeUser(0.,1.);
  effPi->SetMarkerStyle(20);
  effPi->Draw("ALP");
  conPi->SetMarkerStyle(21);
  conPi->Draw("sameLP");

  pCan->cd(2);
  effKa->SetTitle("Kaons");
  effKa->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  effKa->GetYaxis()->SetRangeUser(0.,1.);
  effKa->SetMarkerStyle(20);
  effKa->Draw("ALP");
  conKa->SetMarkerStyle(21);
  conKa->Draw("sameLP");

  pCan->cd(3);
  effPr->SetTitle("Protons");
  effPr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  effPr->GetYaxis()->SetRangeUser(0.,1.);
  effPr->SetMarkerStyle(20);
  effPr->Draw("ALP");
  conPr->SetMarkerStyle(21);
  conPr->Draw("sameLP");

  TFile *outFile = new TFile("HmpidGraphs.root","recreate");
  pCan->Write();
  outFile->Close();

  AliAnalysisTaskSE::Terminate();

}

//___________________________________________________________________________
void AliHMPIDAnalysisTask::CreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //

  //slot #1
   OpenFile(0);
   fHmpHistList = new TList();

   fHmpPesdPhmp = new TH2F("fHmpPesdPhmp","HMPID: ESD p - running p;HMP p (GeV/c);ESD p (Gev/c)",100,0,10,100,0,10);
   fHmpHistList->Add(fHmpPesdPhmp);

   fHmpCkovPesd = new TH2F("fHmpCkovPesd","HMPID: ThetaCherenkov vs P;p_esd (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPesd);

   fHmpCkovPhmp = new TH2F("fHmpCkovPhmp","HMPID: ThetaCherenkov vs P;p_hmp (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPhmp);

   fHmpMipTrkDist = new TH1F("fHmpMipTrkDist","HMPID MIP-Track distance;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDist);

   fHmpMipTrkDistX = new TH1F("fHmpMipTrkDistX","HMPID MIP-Track distance in local X;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistX);

   fHmpMipTrkDistY = new TH1F("fHmpMipTrkDistY","HMPID MIP-Track distance in local Y;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistY);

   fHmpMipCharge3cm = new TH1F("fHmpMipCharge3cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge3cm);

   fHmpMipCharge1cm = new TH1F("fHmpMipCharge1cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge1cm);

   fHmpNumPhots = new TH1F("fHmpNumPhots","HMPID Number of photon clusters on ring;#photon clus.;Entries",51,-0.5,50.5);
   fHmpHistList->Add(fHmpNumPhots);

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

   fThetavsPiFromK= new TH2F("ThetavsPiFromK","Theta vs p of pions from K;p_esd (GeV/c);#Theta_C",100,0,10,110,0,1.1);
   fHmpHistList->Add(fThetavsPiFromK);

   fThetapivsPesd = new TH2F("ThetapivsPesd","Theta of pions vs p of esd;p_esd (GeV/c);#Theta_C",100,0,10,110,0,1.1);
   fHmpHistList->Add(fThetapivsPesd);

   fThetaKvsPesd  = new TH2F("ThetaKvsPesd","Theta of kaons vs p of esd;p_esd (GeV/c);#Theta_C",100,0,10,110,0,1.1);
   fHmpHistList->Add(fThetaKvsPesd);

   fThetaPvsPesd  = new TH2F("ThetaPvsPesd","Theta of protons vs p of esd;p_esd (GeV/c);#Theta_C",100,0,10,110,0,1.1);
   fHmpHistList->Add(fThetaPvsPesd);

}
//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDAnalysisTask::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
   
#endif
