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
  fHmpMipTrkDistX(0x0),fHmpMipTrkDistY(0x0),fHmpMipCharge3cm(0x0),
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
  fTree(0x0)
{
  //
  //Default ctor
  //
  for (Int_t i=0; i<28; i++) fVar[i]=0;
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fMC(0x0), fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpMipTrkDistX(0x0), fHmpMipTrkDistY(0x0),
  fHmpMipCharge3cm(0x0),
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
  fTree(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<28; i++) fVar[i]=0;

  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
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
    fHmpMipTrkDistX  = c.fHmpMipTrkDistX;
    fHmpMipTrkDistY  = c.fHmpMipTrkDistY;
    fHmpMipCharge3cm = c.fHmpMipCharge3cm;
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
    fTree            = c.fTree;
    for (Int_t i=0; i<28; i++) fVar[i]=c.fVar[i];
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const AliHMPIDTaskQA& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fMC(c.fMC),fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpMipTrkDistX(c.fHmpMipTrkDistX),fHmpMipTrkDistY(c.fHmpMipTrkDistY),
  fHmpMipCharge3cm(c.fHmpMipCharge3cm),
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
  fProtCon(c.fProtCon),
  fTree(c.fTree)
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<28; i++) fVar[i]=c.fVar[i];
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
void AliHMPIDTaskQA::ConnectInputData(Option_t *)
{
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
  Double_t priors[5]={1.,1.,1.,1.,1.}; //{0.01,0.01,0.83,0.10,0.5};
  Double_t probs[5];
  AliPID *pPid = new AliPID();
  pPid->SetPriors(priors);
  Double_t n = 1.293;
  Double_t dGeVMass[] = {0.000511,0.105658,0.13957018,0.493677,0.938272};
  AliESDtrack *track=0;
  TParticle *pPart=0;
  AliStack* pStack = 0;
  Int_t label;
  if (fUseMC){
    pStack = fMC->Stack();
  }

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
    else if(Equal(track->GetHMPIDsignal(),-5.,ktol))  fHmpTrkFlags->Fill(2);
    else if(Equal(track->GetHMPIDsignal(),-11.,ktol)) fHmpTrkFlags->Fill(3);
    else if(Equal(track->GetHMPIDsignal(),-22.,ktol)) fHmpTrkFlags->Fill(4);
    else fHmpTrkFlags->Fill(4);

    if(Equal(track->GetHMPIDsignal(),-20.,ktol)) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;

    Int_t q, nph;
    Float_t x, y;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);

    if(Equal(x,0.,ktol) && Equal(y,0.,ktol) && Equal(xpc,0.,ktol) && Equal(ypc,0.,ktol)) continue;

    Double_t dist = TMath::Sqrt( (xpc-x)*(xpc-x) + (ypc - y)*(ypc - y));
    if(q>100.) fHmpMipTrkDistX->Fill(xpc - x);
    if(q>100.) fHmpMipTrkDistY->Fill(ypc - y);
    Double_t pHmp[3] = {0}, pHmp3 = 0;
    if (track->GetOuterHmpPxPyPz(pHmp)) pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);
    if (dist <= 3.0) fHmpMipCharge3cm->Fill(q);

    Float_t b[2];
    Float_t bCov[3];
    track->GetImpactParameters(b,bCov);    

    track->GetHMPIDpid(probs);
    pPid->SetProbabilities(probs);
    if (fUseMC){
      if ((label = track->GetLabel()) < 0) continue;
      pPart = pStack->Particle(label);
    }

    if(track->GetHMPIDsignal() > 0 ){

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

    fVar[0] = track->GetHMPIDcluIdx()/1000000;
    fVar[1] = pHmp3;
    fVar[2] = dPtr;
    fVar[3] = xpc;
    fVar[4] = ypc;
    fVar[5] = x;
    fVar[6] = y;
    fVar[7] = (Float_t)track->GetHMPIDsignal();
    fVar[8] = q;
    fVar[9] = th;
    fVar[10] = ph;
    fVar[11] = (Float_t)track->GetSign();
    fVar[12] = (Float_t)nph;
    fVar[13] = (Float_t)track->GetNcls(1);
    fVar[14] = (Float_t)probs[0];
    fVar[15] = (Float_t)probs[1];
    fVar[16] = (Float_t)probs[2];
    fVar[17] = (Float_t)probs[3];
    fVar[18] = (Float_t)probs[4];
    fVar[19] = (Float_t)track->GetTOFsignal();
    fVar[20] = (Float_t)track->GetKinkIndex(0);
    fVar[21] = (Float_t)track->Xv();
    fVar[22] = (Float_t)track->Yv();
    fVar[23] = (Float_t)track->Zv();
    fVar[24] = (Float_t)track->GetTPCchi2();
    fVar[25] = b[0];
    fVar[26] = b[1];
    fVar[27] = track->GetHMPIDcluIdx()%1000000/1000;
    fTree->Fill();
  }//track loop
  delete pPid;

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
  PostData(2,fTree);
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
   OpenFile(1);
   fHmpHistList = new TList();

   fHmpMipTrkDistX = new TH1F("fHmpMipTrkDistX","HMPID MIP-Track distance in local X;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistX);

   fHmpMipTrkDistY = new TH1F("fHmpMipTrkDistY","HMPID MIP-Track distance in local Y;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistY);

   fHmpMipCharge3cm = new TH1F("fHmpMipCharge3cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge3cm);

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

   OpenFile(2);
   fTree = new TTree("Tree","Tree with data");
   fTree->Branch("Chamber",&fVar[0]);
   fTree->Branch("pHmp3",&fVar[1]);
   fTree->Branch("P",&fVar[2]);
   fTree->Branch("Xpc",&fVar[3]);
   fTree->Branch("Ypc",&fVar[4]);
   fTree->Branch("X",&fVar[5]);
   fTree->Branch("Y",&fVar[6]);
   fTree->Branch("HMPIDsignal",&fVar[7]);
   fTree->Branch("Charge",&fVar[8]);
   fTree->Branch("Theta",&fVar[9]);
   fTree->Branch("Phi",&fVar[10]);
   fTree->Branch("Sign",&fVar[11]);
   fTree->Branch("NumPhotons",&fVar[12]);
   fTree->Branch("NumTPCclust",&fVar[13]);
   fTree->Branch("Prob0",&fVar[14]);
   fTree->Branch("Prob1",&fVar[15]);
   fTree->Branch("Prob2",&fVar[16]);
   fTree->Branch("Prob3",&fVar[17]);
   fTree->Branch("Prob4",&fVar[18]);
   fTree->Branch("TOFsignal",&fVar[19]);
   fTree->Branch("KinkIndex",&fVar[20]);
   fTree->Branch("Xv",&fVar[21]);
   fTree->Branch("Yv",&fVar[22]);
   fTree->Branch("Zv",&fVar[23]);
   fTree->Branch("TPCchi2",&fVar[24]);
   fTree->Branch("b0",&fVar[25]);
   fTree->Branch("b1",&fVar[26]);
   fTree->Branch("ClustSize",&fVar[27]);
}

//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDTaskQA::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
   
#endif
