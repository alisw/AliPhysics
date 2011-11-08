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
//
// QA class of Heavy Flavor quark and fragmeted/decayed particles
// -Check kinematics of Heavy Quarks/hadrons, and decayed leptons
//    pT, rapidity
//    decay lepton kinematics w/wo acceptance
//    heavy hadron decay length, electron pT fraction carried from decay
// -Check yield of Heavy Quarks/hadrons
//    Number of produced heavy quark
//    Number of produced hadron of given pdg code
//
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include <TH2F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <AliAODMCParticle.h>
#include <AliStack.h>

#include "AliHFEmcQA.h"
#include "AliHFEtools.h"
#include "AliHFEcollection.h"

ClassImp(AliHFEmcQA)

//_______________________________________________________________________________________________
    AliHFEmcQA::AliHFEmcQA() :
         fMCEvent(NULL)
        ,fMCHeader(NULL)
        ,fMCArray(NULL)
        ,fQAhistos(NULL)
        ,fMCQACollection(NULL)
        ,fNparents(0)
        ,fCentrality(0)
        ,fIsPbPb(kFALSE)
        ,fIsppMultiBin(kFALSE)
{
        // Default constructor
  for(Int_t mom = 0; mom < 9; mom++){
    fhD[mom] = NULL;
    fhDLogbin[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 50; mom++){
    fHeavyQuark[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 2; mom++){
    fIsHeavy[mom] = 0;
  }
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
}

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA(const AliHFEmcQA&p):
TObject(p)
        ,fMCEvent(NULL)
        ,fMCHeader(NULL)
        ,fMCArray(NULL)
        ,fQAhistos(p.fQAhistos)
        ,fMCQACollection(p.fMCQACollection)
        ,fNparents(p.fNparents)
        ,fCentrality(0)
        ,fIsPbPb(kFALSE)
        ,fIsppMultiBin(kFALSE)
{
        // Copy constructor
  for(Int_t mom = 0; mom < 9; mom++){
    fhD[mom] = NULL;
    fhDLogbin[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 50; mom++){
    fHeavyQuark[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 2; mom++){
    fIsHeavy[mom] = 0;
  }
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
}

//_______________________________________________________________________________________________
AliHFEmcQA&
AliHFEmcQA::operator=(const AliHFEmcQA &)
{
  // Assignment operator

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEmcQA::~AliHFEmcQA()
{
        // Destructor

        AliInfo("Analysis Done.");
}

//_______________________________________________________________________________________________
void AliHFEmcQA::PostAnalyze() const
{
        //
        // Post analysis
        //
}

//_______________________________________________________________________________________________
void AliHFEmcQA::SetBackgroundWeightFactor(Double_t *elecBackgroundFactor, Double_t *binLimit)
{
   //
   // copy background weighting factors into data member
   //
  
  memcpy(fElecBackgroundFactor,elecBackgroundFactor,sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memcpy(fBinLimit,binLimit,sizeof(Double_t) * (kBgPtBins+1));
}

//__________________________________________
void AliHFEmcQA::CreatDefaultHistograms(TList * const qaList)
{      
  //
  // make default histograms
  //
  
  if(!qaList) return;

  fQAhistos = qaList;
  fQAhistos->SetName("MCqa");

  CreateHistograms(AliHFEmcQA::kCharm,0,"mcqa_");               // create histograms for charm
  CreateHistograms(AliHFEmcQA::kBeauty,0,"mcqa_");              // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kOthers,0,"mcqa_");              // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kCharm,1,"mcqa_barrel_");        // create histograms for charm 
  CreateHistograms(AliHFEmcQA::kBeauty,1,"mcqa_barrel_");       // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kOthers,1,"mcqa_barrel_");       // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kCharm,2,"mcqa_unitY_");         // create histograms for charm 
  CreateHistograms(AliHFEmcQA::kBeauty,2,"mcqa_unitY_");        // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kOthers,2,"mcqa_unitY_");        // create histograms for beauty
 
// check D spectra
  TString kDspecies[9];
  kDspecies[0]="411";
  kDspecies[1]="421";
  kDspecies[2]="431";
  kDspecies[3]="4122";
  kDspecies[4]="4132";
  kDspecies[5]="4232";
  kDspecies[6]="4332";
  kDspecies[7]="413";
  kDspecies[8]="423";

  const Double_t kPtbound[2] = {0.1, 20.}; //bin taken for considering inclusive e analysis binning
  Int_t iBin[2];
  iBin[0] = 44; // bins in pt
  iBin[1] = 23; // bins in pt
  Double_t* binEdges[1];
  binEdges[0] =  AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);

  // bin size is chosen to consider ALICE D measurement
  Double_t xbins[15]={0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,12,16,20,30,40,50};
  Double_t ybins[10]={-7.5,-1.0,-0.9,-0.8,-0.5,0.5,0.8,0.9,1.0,7.5};
  TString hname;
  for (Int_t iDmeson=0; iDmeson<9; iDmeson++){
     hname = "Dmeson"+kDspecies[iDmeson];
     fhD[iDmeson] = new TH2F(hname,hname+";p_{T} (GeV/c)",14,xbins,9,ybins);
     hname = "DmesonLogbin"+kDspecies[iDmeson];
     fhDLogbin[iDmeson] = new TH1F(hname,hname+";p_{T} (GeV/c)",iBin[0],binEdges[0]);
     if(fQAhistos) fQAhistos->Add(fhD[iDmeson]);
     if(fQAhistos) fQAhistos->Add(fhDLogbin[iDmeson]);
  }

  const Double_t kPtRange[24] = {0.,0.3,0.4,0.5,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,5.,6.,7.,20.,30.}; // to cope with Ana's bin

  Int_t kNcent;
  if(fIsPbPb) kNcent=11;
  else
  {
      if(fIsppMultiBin) kNcent=8;
      else kNcent = 1;
  }

  fMCQACollection = new AliHFEcollection("TaskMCQA", "MC QA histos for meason pt spectra");

  for(Int_t centbin=0; centbin<kNcent; centbin++)
  {
      fMCQACollection->CreateTH1Farray(Form("pionspectra_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("etaspectra_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("omegaspectra_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("phispectra_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("etapspectra_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("rhospectra_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[1],kPtRange);

      fMCQACollection->CreateTH1F(Form("pionspectraLog_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etaspectraLog_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("omegaspectraLog_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("phispectraLog_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etapspectraLog_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("rhospectraLog_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("kaonspectraLog_centrbin%i",centbin), "kaon yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("k0LspectraLog_centrbin%i",centbin), "k0L yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);

      fMCQACollection->CreateTH1F(Form("piondaughters_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etadaughters_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("omegadaughters_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("phidaughters_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etapdaughters_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("rhodaughters_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
  }

  fQAhistos->Add(fMCQACollection->GetList());

}
  
//__________________________________________
void AliHFEmcQA::CreateHistograms(const Int_t kquark, Int_t icut, TString hnopt) 
{
  // create histograms

  if (!(kquark == kCharm || kquark == kBeauty || kquark == kOthers)) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 

  TString kqTypeLabel[fgkqType];
  if (kquark == kCharm){
    kqTypeLabel[kQuark]="c";
    kqTypeLabel[kantiQuark]="cbar";
    kqTypeLabel[kHadron]="cHadron";
    kqTypeLabel[keHadron]="ceHadron";
    kqTypeLabel[kDeHadron]="nullHadron";
    kqTypeLabel[kElectron]="ce";
    kqTypeLabel[kElectron2nd]="nulle";
  } else if (kquark == kBeauty){
    kqTypeLabel[kQuark]="b";
    kqTypeLabel[kantiQuark]="bbar";
    kqTypeLabel[kHadron]="bHadron";
    kqTypeLabel[keHadron]="beHadron";
    kqTypeLabel[kDeHadron]="bDeHadron";
    kqTypeLabel[kElectron]="be";
    kqTypeLabel[kElectron2nd]="bce";
  } else if (kquark == kOthers){
    kqTypeLabel[kGamma-4]="gammae";
    kqTypeLabel[kPi0-4]="pi0e";
    kqTypeLabel[kElse-4]="elsee";
    kqTypeLabel[kMisID-4]="miside";
  }

  const Double_t kPtbound[2] = {0.1, 20.}; //bin taken for considering inclusive e analysis binning
  Int_t iBin[1];
  iBin[0] = 44; // bins in pt
  Double_t* binEdges[1];
  binEdges[0] =  AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);


  Double_t xcorrbin[501];
  for (int icorrbin = 0; icorrbin< 501; icorrbin++){
    xcorrbin[icorrbin]=icorrbin*0.1;
  }

  TString hname; 
  if(kquark == kOthers){
    for (Int_t iqType = 0; iqType < 4; iqType++ ){
       hname = hnopt+"Pt_"+kqTypeLabel[iqType];
       //fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",iBin[0],binEdges[0]); //mj to compare with FONLL
       fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",60,0.25,30.25);
       hname = hnopt+"Y_"+kqTypeLabel[iqType];
       fHist[iq][iqType][icut].fY = new TH1F(hname,hname,150,-7.5,7.5);
       hname = hnopt+"Eta_"+kqTypeLabel[iqType];
       fHist[iq][iqType][icut].fEta = new TH1F(hname,hname,150,-7.5,7.5);
       // Fill List
       if(fQAhistos) fHist[iq][iqType][icut].FillList(fQAhistos);
    }
    return;
  }
  for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){
     if (iqType < keHadron && icut > 0) continue; // don't duplicate histogram for quark and hadron
     hname = hnopt+"PdgCode_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fPdgCode = new TH1F(hname,hname,20001,-10000.5,10000.5);
     hname = hnopt+"Pt_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",iBin[0],binEdges[0]);
     //fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",60,0.,30.); // mj to compare with FONLL
     //fHist[iq][iqType][icut].fPt = new TH1F(hname,hname+";p_{T} (GeV/c)",60,0.25,30.25); // mj to compare with FONLL
     hname = hnopt+"Y_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fY = new TH1F(hname,hname,150,-7.5,7.5);
     hname = hnopt+"Eta_"+kqTypeLabel[iqType];
     fHist[iq][iqType][icut].fEta = new TH1F(hname,hname,150,-7.5,7.5);
     // Fill List
     if(fQAhistos) fHist[iq][iqType][icut].FillList(fQAhistos);
  }

  if (icut == 0){ 
    hname = hnopt+"Nq_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fNq = new TH1F(hname,hname,50,-0.5,49.5);
    hname = hnopt+"ProcessID_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fProcessID = new TH1F(hname,hname,21,-10.5,10.5);
  }
  hname = hnopt+"ePtRatio_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",200,0,20,100,0,1);
  hname = hnopt+"PtCorr_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fPtCorr = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",500,xcorrbin,44,binEdges[0]);
  //fHistComm[iq][icut].fPtCorr = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",200,0,20,200,0,20);
  hname = hnopt+"PtCorrDp_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fPtCorrDp= new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",500,xcorrbin,44,binEdges[0]);
  //fHistComm[iq][icut].fPtCorrDp= new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",200,0,20,200,0,20);
  hname = hnopt+"PtCorrD0_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fPtCorrD0 = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",500,xcorrbin,44,binEdges[0]);
  //fHistComm[iq][icut].fPtCorrD0 = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",200,0,20,200,0,20);
  hname = hnopt+"PtCorrDrest_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fPtCorrDrest = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",500,xcorrbin,44,binEdges[0]);
  //fHistComm[iq][icut].fPtCorrDrest = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",200,0,20,200,0,20);
  hname = hnopt+"DePtRatio_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fDePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,20,100,0,1);
  hname = hnopt+"eDistance_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].feDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,20,200,0,2);
  hname = hnopt+"DeDistance_"+kqTypeLabel[kQuark];
  fHistComm[iq][icut].fDeDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,20,200,0,2);
  if(fQAhistos) fHistComm[iq][icut].FillList(fQAhistos);
}

//__________________________________________
void AliHFEmcQA::Init()
{
  // called at begining every event
  
  for (Int_t i=0; i<2; i++){
     fIsHeavy[i] = 0;
  } 

  fNparents = 7;

  fParentSelect[0][0] =  411;  
  fParentSelect[0][1] =  421;
  fParentSelect[0][2] =  431;
  fParentSelect[0][3] = 4122;
  fParentSelect[0][4] = 4132;
  fParentSelect[0][5] = 4232;
  fParentSelect[0][6] = 4332;

  fParentSelect[1][0] =  511;
  fParentSelect[1][1] =  521;
  fParentSelect[1][2] =  531;
  fParentSelect[1][3] = 5122;
  fParentSelect[1][4] = 5132;
  fParentSelect[1][5] = 5232;
  fParentSelect[1][6] = 5332;

}

//__________________________________________
void AliHFEmcQA::GetMesonKine() 
{
  //
  // get meson pt spectra
  //

  AliVParticle *mctrack2 = NULL;
  AliMCParticle *mctrack0 = NULL;
  AliVParticle *mctrackdaugt= NULL;
  AliMCParticle *mctrackd= NULL;
  Int_t id1=0, id2=0;


  if(fCentrality>11) printf("warning centrality out of histogram array limits \n");


  for(Int_t imc = 0; imc <fMCEvent->GetNumberOfPrimaries(); imc++){
     if(!(mctrack2 = fMCEvent->GetTrack(imc))) continue;
     TParticle* mcpart0 = fMCEvent->Stack()->Particle(imc);
     if(!mcpart0) continue;
     mctrack0 = dynamic_cast<AliMCParticle *>(mctrack2);
     if(!mctrack0) continue;

     if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

     if(abs(mctrack0->PdgCode()) == 111) // pi0 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("pionspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("pionspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==2)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("piondaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 221) // eta 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("etaspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          } 
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==2||(id2-id1)==3)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("etadaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 223) // omega
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("omegaspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("omegaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==1||(id2-id1)==2)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("omegadaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 333) // phi 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("phispectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("phispectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          } 
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==1)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("phidaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 331) // eta prime
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("etapspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etapspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==2||(id2-id1)==3)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("etapdaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 113) // rho
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("rhospectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("rhospectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
          id1=mctrack0->GetFirstDaughter();
          id2=mctrack0->GetLastDaughter();
          if(!((id2-id1)==1)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("rhodaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 321) // kaon+-
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("kaonspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     else if(abs(mctrack0->PdgCode()) == 130) // k0L
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("k0LspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     }

}
//__________________________________________
void AliHFEmcQA::GetQuarkKine(TParticle *part, Int_t iTrack, const Int_t kquark) 
{
  // get heavy quark kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 

    if (iTrack < 0 || !part) { 
      AliDebug(1, "Stack label is negative or no mcparticle, return\n");
      return; 
    }

    AliMCParticle *mctrack = NULL;
    Int_t partPdgcode = TMath::Abs(part->GetPdgCode());

    // select heavy hadron or not fragmented heavy quark 
    if ( int(partPdgcode/100.)==kquark || int(partPdgcode/1000.)==kquark || (partPdgcode==kquark && (part->GetNDaughters()==0 && iTrack>5)) ){ 

      TParticle *partMother;
      Int_t iLabel;

      if (partPdgcode == kquark){ // in case of not fragmented heavy quark  
        partMother = part; 
        iLabel = iTrack;
      } else{ // in case of heavy hadron, start to search for mother heavy parton 
        iLabel = part->GetFirstMother(); 
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
        if (iLabel>-1) { partMother = mctrack->Particle(); }
        else {
          AliDebug(1, "Stack label is negative, return\n");
          return; 
        }
      }

      // heavy parton selection as a mother of heavy hadron 
      // if the heavy particle comes from string which is denoted as particle status 12|12|12...12|11,[PYTHIA p.60]
      // in this case, the mother of heavy particle can be one of the fragmented parton of the string
      // should I make a condition that partMother should be quark or diquark? -> not necessary
      if ( abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11) ){
      //if ( abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11 || partMother->GetStatusCode() == 12) ){

        if ( abs(partMother->GetPdgCode()) != kquark ){
          // search fragmented partons in the same string
          Bool_t isSameString = kTRUE; 
          for (Int_t i=1; i<fgkMaxIter; i++){
             iLabel = iLabel - 1;
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
             if (iLabel>-1) { partMother = mctrack->Particle(); }
             else {
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }
             if ( abs(partMother->GetPdgCode()) == kquark ) break;
             if ( partMother->GetStatusCode() != 12 ) isSameString = kFALSE;
             if (!isSameString) return; 
          }
        }
        AliDebug(1, "Can not find heavy parton of this heavy hadron in the string, return\n");
        if (abs(partMother->GetPdgCode()) != kquark) return; 

        if (fIsHeavy[iq] >= 50) return;  
        fHeavyQuark[fIsHeavy[iq]] = partMother;
        fIsHeavy[iq]++;

        // fill kinematics for heavy parton
        if (partMother->GetPdgCode() > 0) { // quark
          fHist[iq][kQuark][0].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kQuark][0].fPt->Fill(partMother->Pt());
          fHist[iq][kQuark][0].fY->Fill(AliHFEtools::GetRapidity(partMother));
          fHist[iq][kQuark][0].fEta->Fill(partMother->Eta());
        } else{ // antiquark
          fHist[iq][kantiQuark][0].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kantiQuark][0].fPt->Fill(partMother->Pt());
          fHist[iq][kantiQuark][0].fY->Fill(AliHFEtools::GetRapidity(partMother));
          fHist[iq][kantiQuark][0].fEta->Fill(partMother->Eta());
        }

      } // end of heavy parton slection loop 

    } // end of heavy hadron or quark selection

}

//__________________________________________
void AliHFEmcQA::EndOfEventAna(const Int_t kquark)
{
  // end of event analysis

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 


  // # of heavy quark per event
  AliDebug(1,Form("Number of heavy quark in this event = %d \n",fIsHeavy[iq]));
  fHistComm[iq][0].fNq->Fill(fIsHeavy[iq]);

  Int_t motherID[fgkMaxGener];
  Int_t motherType[fgkMaxGener];
  Int_t motherLabel[fgkMaxGener];
  Int_t ancestorPdg[fgkMaxGener];
  Int_t ancestorLabel[fgkMaxGener];

  for (Int_t i = 0; i < fgkMaxGener; i++){ // initialization
     motherID[i] = 0;
     motherType[i] = 0;
     motherLabel[i] = 0;
     ancestorPdg[i] = 0;
     ancestorLabel[i] = 0;
  }


  // check history of found heavy quarks
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){

     if(!fHeavyQuark[i]) return;

     ancestorLabel[0] = i;
     ancestorPdg[0] = fHeavyQuark[i]->GetPdgCode(); 
     ancestorLabel[1] = fHeavyQuark[i]->GetFirstMother(); 

     AliDebug(1,Form("pdg code= %d\n",ancestorPdg[0]));
     AliDebug(1,Form("ancestor label= %d\n",ancestorLabel[1]));

     Int_t ig = 1;
     while (ancestorLabel[ig] != -1){
          // in case there is mother, get mother's pdg code and grandmother's label
          IdentifyMother(ancestorLabel[ig], ancestorPdg[ig], ancestorLabel[ig+1]); 
          // if mother is still heavy, find again mother's ancestor
          if (ancestorPdg[ig-1] == ancestorPdg[ig]) {
            ig++;
            continue; // if it is from same heavy
          }
          // if the heavy's mother is not heavy, check the mother's label to know if it comes from inital or final parton shower
          if (IsFromInitialShower(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          if (IsFromFinalParton(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          // if it is not the above case, something is strange
          ReportStrangeness(motherID[i],motherType[i],motherLabel[i]);
          break;
     } 
     if (ancestorLabel[ig] == -1){ // from hard scattering
       HardScattering(kquark, motherID[i],motherType[i], motherLabel[i]);
     }

  } // end of found heavy quark loop


  // check process type
  Int_t processID = 0;
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){
     AliDebug(1,Form("Mother ID= %d type= %d label= %d\n",motherID[i],motherType[i],motherLabel[i]));
  }


  Int_t nheavypair = Int_t(fIsHeavy[iq]/2.); 
  for (Int_t ipair = 0; ipair < nheavypair; ipair++){

     Int_t id1 = ipair*2;
     Int_t id2 = ipair*2 + 1;

     if (motherType[id1] == 2 && motherType[id2] == 2){
       if (motherLabel[id1] == motherLabel[id2]) processID = kGluonSplitting; // gluon spliting
       else processID = -9;
     }
     else if (motherType[id1] == -1 && motherType[id2] == -1) {
       if (motherLabel[id1] == -1 && motherLabel[id2] == -1) {
         if (motherID[id1] == fgkGluon) processID = kPairCreationFromg; // gluon fusion
         else processID = kPairCreationFromq; // q-qbar pair creation
       }
       else processID = -8;
     }
     else if (motherType[id1] == -1 || motherType[id2] == -1) {
       if ((motherLabel[id1] == -1 || motherLabel[id2] == -1) && (motherLabel[id1]*motherLabel[id2] == -2 || motherLabel[id1]*motherLabel[id2] == -3)) {
         if(motherID[id1]*motherID[id2] == kquark*fgkGluon) processID = kFlavourExitation; // flavour exitation 
         else processID = kLightQuarkShower;
       }
       else processID = -7;
     }
     else if (motherType[id1] == -2 || motherType[id2] == -2) {
       if (motherLabel[id1] == motherLabel[id2]) processID = kInitialPartonShower; // initial parton shower
       else processID = -6;
       
     }
     else processID = -5;

     if (nheavypair >1) AliDebug(1,Form("Multi pair found : process ID = %d\n",processID));
     else fHistComm[iq][0].fProcessID->Fill(processID);
     AliDebug(1,Form("Process ID = %d\n",processID));
  } // end of # heavy quark pair loop

}

//__________________________________________
void AliHFEmcQA::GetHadronKine(TParticle* mcpart, const Int_t kquark)
{
    // decay electron kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return;
    }
    Int_t iq = kquark - kCharm;

    if(!mcpart){
      AliDebug(1, "no mc particle, return\n");
      return;
    }

    Int_t iLabel = mcpart->GetFirstMother();
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return;
    }

    TParticle *partCopy = mcpart;
    Int_t pdgcode = mcpart->GetPdgCode();
    Int_t pdgcodeCopy = pdgcode;

    AliMCParticle *mctrack = NULL;

    // if the mother is charmed hadron  
    Bool_t isDirectCharm = kFALSE;
    if ( int(abs(pdgcode)/100.) == kCharm || int(abs(pdgcode)/1000.) == kCharm ) {

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = mcpart->GetFirstMother();
             if (jLabel == -1){
               isDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return;
             }
             // if there is an ancester
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return; 
             TParticle* mother = mctrack->Particle();
             Int_t motherPDG = mother->GetPdgCode();
    
             for (Int_t j=0; j<fNparents; j++){
                if (abs(motherPDG)==fParentSelect[1][j]) return; // return if this hadron is originated from b
             }

             mcpart = mother;
          } // end of iteration 
    } // end of if
    if((isDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(pdgcodeCopy)==fParentSelect[iq][i]){

              // fill hadron kinematics
              fHist[iq][kHadron][0].fPdgCode->Fill(pdgcodeCopy);
              fHist[iq][kHadron][0].fPt->Fill(partCopy->Pt());
              fHist[iq][kHadron][0].fY->Fill(AliHFEtools::GetRapidity(partCopy));
              fHist[iq][kHadron][0].fEta->Fill(partCopy->Eta());

              if(iq==0) {
               fhD[i]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
               if(TMath::Abs(AliHFEtools::GetRapidity(partCopy))<0.5) fhDLogbin[i]->Fill(partCopy->Pt());
              }
            }
         }
	 // I also want to store D* info to compare with D* measurement 
	 if (abs(pdgcodeCopy)==413 && iq==0) { //D*+
               fhD[7]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
               if(TMath::Abs(AliHFEtools::GetRapidity(partCopy))<0.5) fhDLogbin[7]->Fill(partCopy->Pt());
	 }
	 if (abs(pdgcodeCopy)==423 && iq==0) { //D*0
               fhD[8]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
               if(TMath::Abs(AliHFEtools::GetRapidity(partCopy))<0.5) fhDLogbin[8]->Fill(partCopy->Pt());
	 }
    } // end of if
}

//__________________________________________
void AliHFEmcQA::GetDecayedKine(TParticle* mcpart, const Int_t kquark, Int_t kdecayed, Int_t icut) 
{
    // decay electron kinematics
    
    if (!(kquark == kCharm || kquark == kBeauty || kquark == kOthers)){
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 
    Bool_t isFinalOpenCharm = kFALSE;

    if(!mcpart){
      AliDebug(1, "no mcparticle, return\n");
      return;
    }

    if(kquark==kOthers){
      Int_t esource = -1;
      if ( abs(mcpart->GetPdgCode()) != kdecayed ) esource = kMisID-4;
      else esource =GetSource(mcpart)-4; // return for the cases kGamma=4, kPi0=5, kElse=6
      if(esource==0|| esource==1 || esource==2 || esource==3){
        fHist[iq][esource][icut].fPt->Fill(mcpart->Pt());
        fHist[iq][esource][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
        fHist[iq][esource][icut].fEta->Fill(mcpart->Eta());
        return; 
      }
      else {
        AliDebug(1, "e source is out of defined ranges, return\n");
        return;
      }
    }

    if ( abs(mcpart->GetPdgCode()) != kdecayed ) return;

    Int_t iLabel = mcpart->GetFirstMother(); 
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    AliMCParticle *mctrack = NULL;
    if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
    TParticle *partMother = mctrack->Particle();
    TParticle *partMotherCopy = partMother;
    Int_t maPdgcode = partMother->GetPdgCode();
    Int_t maPdgcodeCopy = maPdgcode;

    // get mc primary vertex
    /*
    TArrayF mcPrimVtx;
    if(fMCHeader) fMCHeader->PrimaryVertex(mcPrimVtx);

    // get electron production vertex   
    TLorentzVector ePoint;
    mcpart->ProductionVertex(ePoint);

    // calculated production vertex to primary vertex (in xy plane)
    Float_t decayLxy = TMath::Sqrt((mcPrimVtx[0]-ePoint[0])*(mcPrimVtx[0]-ePoint[0])+(mcPrimVtx[1]-ePoint[1])*(mcPrimVtx[1]-ePoint[1]));
    */ 
    Float_t decayLxy = 0;

    // if the mother is charmed hadron  
    Bool_t isMotherDirectCharm = kFALSE;
    if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) { 

         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcode)==fParentSelect[0][i]){
              isFinalOpenCharm = kTRUE;
            } 
         }  
         if (!isFinalOpenCharm) return ;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = partMother->GetFirstMother(); 
             if (jLabel == -1){
               isMotherDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }

             // if there is an ancester
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return; 
             TParticle* grandMa = mctrack->Particle();
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){

                  if (kquark == kCharm) return;
                  // fill electron kinematics
                  fHist[iq][kElectron2nd][icut].fPdgCode->Fill(mcpart->GetPdgCode());
                  fHist[iq][kElectron2nd][icut].fPt->Fill(mcpart->Pt());
                  fHist[iq][kElectron2nd][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                  fHist[iq][kElectron2nd][icut].fEta->Fill(mcpart->Eta());

                  // fill mother hadron kinematics
                  fHist[iq][kDeHadron][icut].fPdgCode->Fill(grandMaPDG); 
                  fHist[iq][kDeHadron][icut].fPt->Fill(grandMa->Pt());
                  fHist[iq][kDeHadron][icut].fY->Fill(AliHFEtools::GetRapidity(grandMa));
                  fHist[iq][kDeHadron][icut].fEta->Fill(grandMa->Eta());
                 
                  // ratio between pT of electron and pT of mother B hadron 
                  if(grandMa->Pt()) fHistComm[iq][icut].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());

                  // distance between electron production point and primary vertex
                  fHistComm[iq][icut].fDeDistance->Fill(grandMa->Pt(),decayLxy);
                  return;
                }
             } 

             partMother = grandMa;
          } // end of iteration 
    } // end of if
    if((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcodeCopy)==fParentSelect[iq][i]){

//mj weighting to consider measured spectra!!!
              Double_t mpt=partMotherCopy->Pt();
              Double_t wfactor=2*(703.681*mpt/TMath::Power((1+TMath::Power(mpt/1.73926,2)),2.34821))/(368.608*mpt/TMath::Power((1+TMath::Power(mpt/2.74868,2)),2.34225)); // 2* considering in pythia I used particle + antiparticle differently from the measurement
              //Double_t wfactor=(703.681*mpt/TMath::Power((1+TMath::Power(mpt/1.73926,2)),2.34821))/(368.608*mpt/TMath::Power((1+TMath::Power(mpt/2.74868,2)),2.34225));
              // fill electron kinematics
              if(iq==0){
                fHist[iq][kElectron2nd][icut].fPdgCode->Fill(mcpart->GetPdgCode(),wfactor);
                fHist[iq][kElectron2nd][icut].fPt->Fill(mcpart->Pt(),wfactor);
                fHist[iq][kElectron2nd][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart),wfactor);
                fHist[iq][kElectron2nd][icut].fEta->Fill(mcpart->Eta(),wfactor);  

                fHist[iq][kElectron][icut].fPdgCode->Fill(mcpart->GetPdgCode());
                fHist[iq][kElectron][icut].fPt->Fill(mcpart->Pt());
                fHist[iq][kElectron][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                fHist[iq][kElectron][icut].fEta->Fill(mcpart->Eta());  
              } 
              else{
                fHist[iq][kElectron][icut].fPdgCode->Fill(mcpart->GetPdgCode());
                fHist[iq][kElectron][icut].fPt->Fill(mcpart->Pt());
                fHist[iq][kElectron][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                fHist[iq][kElectron][icut].fEta->Fill(mcpart->Eta());  
              }
//--------------
              // fill mother hadron kinematics
              fHist[iq][keHadron][icut].fPdgCode->Fill(maPdgcodeCopy); 
              fHist[iq][keHadron][icut].fPt->Fill(partMotherCopy->Pt());
              fHist[iq][keHadron][icut].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
              fHist[iq][keHadron][icut].fEta->Fill(partMotherCopy->Eta());

              // ratio between pT of electron and pT of mother B or direct D hadron 
              if(partMotherCopy->Pt()) fHistComm[iq][icut].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());
              fHistComm[iq][icut].fPtCorr->Fill(partMotherCopy->Pt(),mcpart->Pt());
              if(TMath::Abs(partMotherCopy->GetPdgCode())==411) fHistComm[iq][icut].fPtCorrDp->Fill(partMotherCopy->Pt(),mcpart->Pt());
              else if(TMath::Abs(partMotherCopy->GetPdgCode())==421) fHistComm[iq][icut].fPtCorrD0->Fill(partMotherCopy->Pt(),mcpart->Pt());
              else fHistComm[iq][icut].fPtCorrDrest->Fill(partMotherCopy->Pt(),mcpart->Pt());

              // distance between electron production point and primary vertex
              fHistComm[iq][icut].feDistance->Fill(partMotherCopy->Pt(),decayLxy);
            }
         }
    } // end of if
}

//____________________________________________________________________
void  AliHFEmcQA::GetDecayedKine(AliAODMCParticle *mcpart, const Int_t kquark, Int_t kdecayed, Int_t icut)
{
  // decay electron kinematics

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return;
  }

  Int_t iq = kquark - kCharm;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return;
  }

  if ( abs(mcpart->GetPdgCode()) != kdecayed ) return;

  // mother
  Int_t iLabel = mcpart->GetMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return;
  }

  AliAODMCParticle *partMother = (AliAODMCParticle*)fMCArray->At(iLabel);
  AliAODMCParticle *partMotherCopy = partMother;
  Int_t maPdgcode = partMother->GetPdgCode();
  Int_t maPdgcodeCopy = maPdgcode;

  Bool_t isMotherDirectCharm = kFALSE;
  if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    } 
    if (!isFinalOpenCharm) return;

    for (Int_t i=1; i<fgkMaxIter; i++){

       Int_t jLabel = partMother->GetMother();
       if (jLabel == -1){
         isMotherDirectCharm = kTRUE;
         break; // if there is no ancester
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return;
       }

       // if there is an ancester
       AliAODMCParticle* grandMa = (AliAODMCParticle*)fMCArray->At(jLabel);
       Int_t grandMaPDG = grandMa->GetPdgCode();

       for (Int_t j=0; j<fNparents; j++){
          if (abs(grandMaPDG)==fParentSelect[1][j]){

            if (kquark == kCharm) return;
            // fill electron kinematics
            fHist[iq][kElectron2nd][icut].fPdgCode->Fill(mcpart->GetPdgCode());
            fHist[iq][kElectron2nd][icut].fPt->Fill(mcpart->Pt());
            fHist[iq][kElectron2nd][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
            fHist[iq][kElectron2nd][icut].fEta->Fill(mcpart->Eta());

            // fill mother hadron kinematics
            fHist[iq][kDeHadron][icut].fPdgCode->Fill(grandMaPDG);
            fHist[iq][kDeHadron][icut].fPt->Fill(grandMa->Pt());
            fHist[iq][kDeHadron][icut].fY->Fill(AliHFEtools::GetRapidity(grandMa));
            fHist[iq][kDeHadron][icut].fEta->Fill(grandMa->Eta());

            return;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  if ((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcodeCopy)==fParentSelect[iq][i]){

         // fill electron kinematics
         fHist[iq][kElectron][icut].fPdgCode->Fill(mcpart->GetPdgCode());
         fHist[iq][kElectron][icut].fPt->Fill(mcpart->Pt());
         fHist[iq][kElectron][icut].fY->Fill(AliHFEtools::GetRapidity(mcpart));
         fHist[iq][kElectron][icut].fEta->Fill(mcpart->Eta());

         // fill mother hadron kinematics
         fHist[iq][keHadron][icut].fPdgCode->Fill(maPdgcodeCopy);
         fHist[iq][keHadron][icut].fPt->Fill(partMotherCopy->Pt());
         fHist[iq][keHadron][icut].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
         fHist[iq][keHadron][icut].fEta->Fill(partMotherCopy->Eta());

       }
    }
  } // end of if

}

//__________________________________________
void AliHFEmcQA::IdentifyMother(Int_t motherlabel, Int_t &motherpdg, Int_t &grandmotherlabel)
{
       // find mother pdg code and label 

       if (motherlabel < 0) { 
         AliDebug(1, "Stack label is negative, return\n");
         return; 
       }
       AliMCParticle *mctrack = NULL;
       if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(motherlabel))))) return; 
       TParticle *heavysMother = mctrack->Particle();
       motherpdg = heavysMother->GetPdgCode();
       grandmotherlabel = heavysMother->GetFirstMother();
       AliDebug(1,Form("ancestor pdg code= %d\n",motherpdg));
}

//__________________________________________
void AliHFEmcQA::HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -1 means this heavy quark coming from hard vertex

       AliMCParticle *mctrack1 = NULL;
       AliMCParticle *mctrack2 = NULL;
       if(!(mctrack1 = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(4))))) return; 
       if(!(mctrack2 = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(5))))) return; 
       TParticle *afterinitialrad1  = mctrack1->Particle();
       TParticle *afterinitialrad2  = mctrack2->Particle();
           
       motherlabel = -1;

       if (abs(afterinitialrad1->GetPdgCode()) == fgkGluon && abs(afterinitialrad2->GetPdgCode()) == fgkGluon){
         AliDebug(1,"heavy from gluon gluon pair creation!\n");
         mothertype = -1;
         motherID = fgkGluon;
       }
       else if (abs(afterinitialrad1->GetPdgCode()) == kquark || abs(afterinitialrad2->GetPdgCode()) == kquark){ // one from Q and the other from g
         AliDebug(1,"heavy from flavor exitation!\n");
         mothertype = -1;
         motherID = kquark;
       }
       else if  (abs(afterinitialrad1->GetPdgCode()) == abs(afterinitialrad2->GetPdgCode())){
         AliDebug(1,"heavy from q-qbar pair creation!\n");
         mothertype = -1;
         motherID = 1;
       }
       else {
         AliDebug(1,"something strange!\n");
         mothertype = -999;
         motherlabel = -999;
         motherID = -999;
       }
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -2 means this heavy quark coming from initial state 

       AliMCParticle *mctrack = NULL;
       if (inputmotherlabel==2 || inputmotherlabel==3){ // mother exist before initial state radiation
         if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(inputmotherlabel))))) return kFALSE; 
         TParticle *heavysMother = mctrack->Particle();
         motherID = heavysMother->GetPdgCode(); 
         mothertype = -2; // there is mother before initial state radiation
         motherlabel = inputmotherlabel;
         AliDebug(1,"initial parton shower! \n");

         return kTRUE;
       }

       return kFALSE;
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype 2 means this heavy quark coming from final state 

       AliMCParticle *mctrack = NULL;
       if (inputmotherlabel > 5){ // mother exist after hard scattering
         if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(inputmotherlabel))))) return kFALSE; 
         TParticle *heavysMother = mctrack->Particle();
         motherID = heavysMother->GetPdgCode(); 
         mothertype = 2; // 
         motherlabel = inputmotherlabel;
         AliDebug(1,Form("heavy quark from %d after hard scattering! \n",motherID));

         return kTRUE;
       }
       return kFALSE;
}

//__________________________________________
void AliHFEmcQA::ReportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
      // mark strange behavior  

       mothertype = -888;
       motherlabel = -888;
       motherID = -888;
       AliDebug(1,"something strange!\n");
}

//__________________________________________
Int_t AliHFEmcQA::GetSource(const AliAODMCParticle * const mcpart)
{        
  // decay particle's origin 

  //if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;
       
  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "Stack label is negative or no mcparticle, return\n");
    return -1;
  }

  // mother
  Int_t iLabel = mcpart->GetMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  } 
       
  AliAODMCParticle *partMother = (AliAODMCParticle*)fMCArray->At(iLabel);
  Int_t maPdgcode = partMother->GetPdgCode();
  
  // if the mother is charmed hadron  
  if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {
    
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    }
    if (!isFinalOpenCharm) return -1;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<fgkMaxIter; i++){

       Int_t jLabel = partMother->GetMother();
       if (jLabel == -1){
         origin = kDirectCharm;
         return origin;
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return -1;
       }

       // if there is an ancester
       AliAODMCParticle* grandMa = (AliAODMCParticle*)fMCArray->At(jLabel);
       Int_t grandMaPDG = grandMa->GetPdgCode();

       for (Int_t j=0; j<fNparents; j++){
          if (abs(grandMaPDG)==fParentSelect[1][j]){
            origin = kBeautyCharm;
            return origin;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
    for (Int_t i=0; i<fNparents; i++){
       if (abs(maPdgcode)==fParentSelect[1][i]){
         origin = kDirectBeauty;
         return origin;
       }
    }
  } // end of if
  else if ( abs(maPdgcode) == 22 ) {
    origin = kGamma;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 111 ) {
    origin = kPi0;
    return origin;
  } // end of if

  return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetSource(const TParticle * const mcpart)
{
  // decay particle's origin 

  //if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  AliMCParticle *mctrack = NULL;
  if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return -1; 
  TParticle *partMother = mctrack->Particle();
  Int_t maPdgcode = partMother->GetPdgCode();

   // if the mother is charmed hadron  
   if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[0][i]){
          isFinalOpenCharm = kTRUE;
        }
     }
     if (!isFinalOpenCharm) return -1;

     // iterate until you find B hadron as a mother or become top ancester 
     for (Int_t i=1; i<fgkMaxIter; i++){

        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          origin = kDirectCharm;
          return origin;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return -1; 
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t j=0; j<fNparents; j++){
           if (abs(grandMaPDG)==fParentSelect[1][j]){
             origin = kBeautyCharm;
             return origin;
           }
        }

        partMother = grandMa;
     } // end of iteration 
   } // end of if
   else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[1][i]){
          origin = kDirectBeauty;
          return origin;
        }
     }
   } // end of if
   else if ( abs(maPdgcode) == 22 ) {
     origin = kGamma;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 111 ) {
     origin = kPi0;
     return origin;
   } // end of if
   else origin = kElse;

   return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetElecSource(TParticle * const mcpart)
{
  // decay particle's origin 

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return kMisID;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  AliMCParticle *mctrack = NULL;
  Int_t tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return -1; 
  TParticle *partMother = mctrack->Particle();
  TParticle *partMotherCopy = mctrack->Particle();
  Int_t maPdgcode = partMother->GetPdgCode();

   // if the mother is charmed hadron  
   if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[0][i]){
          isFinalOpenCharm = kTRUE;
        }
     }
     if (!isFinalOpenCharm) return -1;

     // iterate until you find B hadron as a mother or become top ancester 
     for (Int_t i=1; i<fgkMaxIter; i++){

        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          origin = kDirectCharm;
          return origin;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return -1; 
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t j=0; j<fNparents; j++){
           if (abs(grandMaPDG)==fParentSelect[1][j]){
             origin = kBeautyCharm;
             return origin;
           }
        }

        partMother = grandMa;
     } // end of iteration 
   } // end of if
   else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
     for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[1][i]){
          origin = kDirectBeauty;
          return origin;
        }
     }
   } // end of if
   else if ( abs(maPdgcode) == 22 ) { //conversion

     tmpMomLabel = partMotherCopy->GetFirstMother();
     if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tmpMomLabel))))) return -1;
     partMother = mctrack->Particle();
     maPdgcode = partMother->GetPdgCode();
     if ( abs(maPdgcode) == 111 ) {
       origin = kGammaPi0;
       return origin;
     } 
     else if ( abs(maPdgcode) == 221 ) {
       origin = kGammaEta;
       return origin;
     } 
     else if ( abs(maPdgcode) == 223 ) {
       origin = kGammaOmega;
       return origin;
     } 
     else if ( abs(maPdgcode) == 333 ) {
       origin = kGammaPhi;
       return origin;
     }
     else if ( abs(maPdgcode) == 331 ) {
       origin = kGammaEtaPrime;
       return origin; 
     }
     else if ( abs(maPdgcode) == 113 ) {
       origin = kGammaRho0;
       return origin;
     }
     else origin = kElse;
     //origin = kGamma; // finer category above
     return origin;

   } // end of if
   else if ( abs(maPdgcode) == 111 ) {
     origin = kPi0;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 221 ) {
     origin = kEta;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 223 ) {
     origin = kOmega;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 333 ) {
     origin = kPhi;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 331 ) {
     origin = kEtaPrime;
     return origin;
   } // end of if
   else if ( abs(maPdgcode) == 113 ) {
     origin = kRho0;
     return origin;
   } // end of if
   else origin = kElse;

   return origin;
}
//__________________________________________
Double_t AliHFEmcQA::GetWeightFactor(AliMCParticle *mctrack, const Int_t iBgLevel){
  //
  // Get weighting factor for the realistic background estimation, for three possible background yield levels, indicated by the argument "iLevel": the best estimate (0), the lower uncertainty level (1), and the upper uncertainty level (2)
  //
  AliMCParticle *mctrackmother = NULL;  
  Double_t weightElecBg = 0.;
  Double_t mesonPt = 0.;
  Double_t bgcategory = 0.;
  Int_t mArr = -1;  
  Int_t mesonID = GetElecSource(mctrack->Particle());
  if(mesonID==kGammaPi0 || mesonID==kPi0) mArr=0;                //pion
  else if(mesonID==kGammaEta || mesonID==kEta) mArr=1;           //eta
  else if(mesonID==kGammaOmega || mesonID==kOmega) mArr=2;       //omega
  else if(mesonID==kGammaPhi || mesonID==kPhi) mArr=3;           //phi
  else if(mesonID==kGammaEtaPrime || mesonID==kEtaPrime) mArr=4; //etaprime
  else if(mesonID==kGammaRho0 || mesonID==kRho0) mArr=5;         //rho

  if(!(mArr<0)){
     if(mesonID>=kGammaPi0) {  // conversion electron, be careful with the enum odering 
        Int_t glabel=TMath::Abs(mctrack->GetMother()); // gamma label
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's label
          if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
            mesonPt = mctrackmother->Pt(); //meson pt
            bgcategory = 1.;
          } 
        }
     }
     else{ // nonHFE except for the conversion electron
        Int_t glabel=TMath::Abs(mctrack->GetMother()); 
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          mesonPt=mctrackmother->Pt(); //meson pt
          bgcategory = -1.;
        }
     }

     if(fIsPbPb){
       if(fCentrality < 0)return 0.;
       weightElecBg=fElecBackgroundFactor[iBgLevel][fCentrality][mArr][kBgPtBins-1];                        
       for(int ii=0; ii<kBgPtBins; ii++){              
	 if((mesonPt > fBinLimit[ii]) && (mesonPt < fBinLimit[ii+1])){
	   weightElecBg = fElecBackgroundFactor[iBgLevel][fCentrality][mArr][ii];
	   break;
	 }
       }
     }
     else{
       weightElecBg=fElecBackgroundFactor[iBgLevel][0][mArr][kBgPtBins-1];                         
       for(int ii=0; ii<kBgPtBins; ii++){              
	 if((mesonPt > fBinLimit[ii]) && (mesonPt < fBinLimit[ii+1])){
	   weightElecBg = fElecBackgroundFactor[iBgLevel][0][mArr][ii];
	   break;
	 }
       }
     }    
  }
  return bgcategory*weightElecBg;
}

//__________________________________________
void AliHFEmcQA::AliHists::FillList(TList *l) const {
  //
  // Fill Histos into a list for output
  //
  if(fPdgCode) l->Add(fPdgCode);
  if(fPt) l->Add(fPt);
  if(fY) l->Add(fY);
  if(fEta) l->Add(fEta);
}

//__________________________________________
void AliHFEmcQA::AliHistsComm::FillList(TList *l) const { 
  //
  // Fill Histos into a list for output
  //
  if(fNq) l->Add(fNq);
  if(fProcessID) l->Add(fProcessID);
  if(fePtRatio) l->Add(fePtRatio);
  if(fPtCorr) l->Add(fPtCorr);
  if(fPtCorrDp) l->Add(fPtCorrDp);
  if(fPtCorrD0) l->Add(fPtCorrD0);
  if(fPtCorrDrest) l->Add(fPtCorrDrest);
  if(fDePtRatio) l->Add(fDePtRatio);
  if(feDistance) l->Add(feDistance);
  if(fDeDistance) l->Add(fDeDistance);
}

