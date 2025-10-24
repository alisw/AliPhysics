/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.   *
 *                                                                          *
 * Author: Ewa Glimos                                                       *
 * Based on AliAnalysisTaskGammaPureMC macro by:                            *
 *                 Baldo Sahlmueller, Friederike Bock                       *
 * Version 1.0                                                              *
 *                                                                          *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on eta prime and its decay products
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskEtaPrimeMCStudies.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include <algorithm>
#include <array>
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskEtaPrimeMCStudies)

//________________________________________________________________________
AliAnalysisTaskEtaPrimeMCStudies::AliAnalysisTaskEtaPrimeMCStudies(): 
    AliAnalysisTaskSE(),
    fOutputContainer(nullptr),
    fInputEvent(NULL),
    fMCEvent(NULL),
    fHistNEvents(nullptr),
    fHistXSection(nullptr),
    fHistPtHard(nullptr),
    fHistEtaPrime_Pt(nullptr),
    fHistEtaPrimeAll_Pt(nullptr),
    fHistEtaFromEtaPrime_Pt(nullptr),
    fHistPiPl_Pt(nullptr),
    fHistPiMi_Pt(nullptr),
    fHistGamma1_Pt(nullptr),
    fHistGamma2_Pt(nullptr),
    fHistPi0_Pt(nullptr),
    fHistEta_Pt(nullptr),
    fHistEtaPrime_InAcc_EMC(nullptr),
    fHistEtaPrime_InAcc_PCMECM(nullptr),
    fTreeEtaPrimeKinematics(nullptr),
    fBuffer_Gamma1_pt(0),
    fBuffer_Gamma1_eta(0),
    fBuffer_Gamma1_phi(0),
    fBuffer_Gamma1_y(0),
    fBuffer_Gamma1_E(0),
    fBuffer_Gamma2_pt(0),
    fBuffer_Gamma2_eta(0),
    fBuffer_Gamma2_phi(0),
    fBuffer_Gamma2_y(0),
    fBuffer_Gamma2_E(0),
    fBuffer_PiPl_pt(0),
    fBuffer_PiPl_E(0),
    fBuffer_PiPl_eta(0),
    fBuffer_PiPl_phi(0),
    fBuffer_PiPl_y(0),
    fBuffer_PiMi_pt(0),
    fBuffer_PiMi_E(0),
    fBuffer_PiMi_eta(0),
    fBuffer_PiMi_phi(0),
    fBuffer_PiMi_y(0),
    fBuffer_Eta_pt(0),
    fBuffer_Eta_E(0),
    fBuffer_Eta_eta(0),
    fBuffer_Eta_phi(0),
    fBuffer_Eta_y(0),
    fBuffer_EtaPrime_pt(0),
    fBuffer_EtaPrime_E(0),
    fBuffer_EtaPrime_eta(0),
    fBuffer_EtaPrime_phi(0),
    fBuffer_EtaPrime_y(0)
{
    
}


//________________________________________________________________________
AliAnalysisTaskEtaPrimeMCStudies::AliAnalysisTaskEtaPrimeMCStudies(const char *name):
    AliAnalysisTaskSE(name),
    fOutputContainer(nullptr),
    fInputEvent(NULL),
    fMCEvent(NULL),
    fHistNEvents(nullptr),
    fHistXSection(nullptr),
    fHistPtHard(nullptr),
    fHistEtaPrime_Pt(nullptr),
    fHistEtaPrimeAll_Pt(nullptr),
    fHistEtaFromEtaPrime_Pt(nullptr),
    fHistPiPl_Pt(nullptr),
    fHistPiMi_Pt(nullptr),
    fHistGamma1_Pt(nullptr),
    fHistGamma2_Pt(nullptr),
    fHistPi0_Pt(nullptr),
    fHistEta_Pt(nullptr),
    fHistEtaPrime_InAcc_EMC(nullptr),
    fHistEtaPrime_InAcc_PCMECM(nullptr),
    fTreeEtaPrimeKinematics(nullptr),
    fBuffer_Gamma1_pt(0),
    fBuffer_Gamma1_eta(0),
    fBuffer_Gamma1_phi(0),
    fBuffer_Gamma1_y(0),
    fBuffer_Gamma1_E(0),
    fBuffer_Gamma2_pt(0),
    fBuffer_Gamma2_eta(0),
    fBuffer_Gamma2_phi(0),
    fBuffer_Gamma2_y(0),
    fBuffer_Gamma2_E(0),
    fBuffer_PiPl_pt(0),
    fBuffer_PiPl_E(0),
    fBuffer_PiPl_eta(0),
    fBuffer_PiPl_phi(0),
    fBuffer_PiPl_y(0),
    fBuffer_PiMi_pt(0),
    fBuffer_PiMi_E(0),
    fBuffer_PiMi_eta(0),
    fBuffer_PiMi_phi(0),
    fBuffer_PiMi_y(0),
    fBuffer_Eta_pt(0),
    fBuffer_Eta_E(0),
    fBuffer_Eta_eta(0),
    fBuffer_Eta_phi(0),
    fBuffer_Eta_y(0),
    fBuffer_EtaPrime_pt(0),
    fBuffer_EtaPrime_E(0),
    fBuffer_EtaPrime_eta(0),
    fBuffer_EtaPrime_phi(0),
    fBuffer_EtaPrime_y(0)
{
    DefineOutput(1,TList::Class());
    DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEtaPrimeMCStudies::~AliAnalysisTaskEtaPrimeMCStudies()
{

}


//________________________________________________________________________
void AliAnalysisTaskEtaPrimeMCStudies::UserCreateOutputObjects(){

    if( fOutputContainer != nullptr){
        delete fOutputContainer;
        fOutputContainer        = nullptr;
    }
    if( fOutputContainer == nullptr){
        fOutputContainer        = new TList();
        fOutputContainer->SetOwner(kTRUE);
    }
    
    // general histograms
    fHistNEvents                = new TH1F("NEvents","NEvents", 3, -0.5, 2.5);
    fHistNEvents->Sumw2();
    fOutputContainer->Add(fHistNEvents);

    fHistXSection               = new TH1D("XSection","XSection", 1000000, 0, 1e4);
    fHistXSection->Sumw2();
    fOutputContainer->Add(fHistXSection);

    fHistPtHard                 = new TH1F("PtHard","PtHard", 400, 0, 200);
    fHistPtHard->Sumw2();
    fOutputContainer->Add(fHistPtHard);

    // eta prime histograms
    fHistEtaPrime_Pt            = new TH1F("EtaPrime_Pt","EtaPrime_Pt", 160, 0, 80);
    fHistEtaPrime_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEtaPrime_Pt->GetYaxis()->SetTitle("N_{#eta'}");
    fHistEtaPrime_Pt->Sumw2();
    fOutputContainer->Add(fHistEtaPrime_Pt);

    fHistEtaPrimeAll_Pt         = new TH1F("EtaPrimeAll_Pt","EtaPrimeAll_Pt", 160, 0, 80);
    fHistEtaPrimeAll_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEtaPrimeAll_Pt->GetYaxis()->SetTitle("N_{#eta'^{all}}");
    fHistEtaPrimeAll_Pt->Sumw2();
    fOutputContainer->Add(fHistEtaPrimeAll_Pt);

    fHistEtaFromEtaPrime_Pt     = new TH1F("EtaFromEtaPrime_Pt","EtaFromEtaPrime_Pt", 160, 0, 80);
    fHistEtaFromEtaPrime_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEtaFromEtaPrime_Pt->GetYaxis()->SetTitle("N_{#eta}");
    fHistEtaFromEtaPrime_Pt->Sumw2();
    fOutputContainer->Add(fHistEtaFromEtaPrime_Pt);

    fHistPiPl_Pt                = new TH1F("PiPl_Pt","PiPl_Pt", 160, 0, 80);
    fHistPiPl_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPiPl_Pt->GetYaxis()->SetTitle("N_{#pi^{+}}");
    fHistPiPl_Pt->Sumw2();
    fOutputContainer->Add(fHistPiPl_Pt);

    fHistPiMi_Pt                = new TH1F("PiMi_Pt","PiMi_Pt", 160, 0, 80);
    fHistPiMi_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPiMi_Pt->GetYaxis()->SetTitle("N_{#pi^{-}}");
    fHistPiMi_Pt->Sumw2();
    fOutputContainer->Add(fHistPiMi_Pt);

    fHistGamma1_Pt              = new TH1F("Gamma1_Pt","Gamma1_Pt", 160, 0, 80);
    fHistGamma1_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistGamma1_Pt->GetYaxis()->SetTitle("N_{#gamma_{1}}");
    fHistGamma1_Pt->Sumw2();
    fOutputContainer->Add(fHistGamma1_Pt);

    fHistGamma2_Pt              = new TH1F("Gamma2_Pt","Gamma2_Pt", 160, 0, 80);
    fHistGamma2_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistGamma2_Pt->GetYaxis()->SetTitle("N_{#gamma_{2}}");\
    fHistGamma2_Pt->Sumw2();
    fOutputContainer->Add(fHistGamma2_Pt);

    fHistPi0_Pt                 = new TH1F("Pi0_Pt","Pi0_Pt", 160, 0, 80);
    fHistPi0_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPi0_Pt->GetYaxis()->SetTitle("N_{#pi^{0}}");
    fHistPi0_Pt->Sumw2();
    fOutputContainer->Add(fHistPi0_Pt);

    fHistEta_Pt                 = new TH1F("Eta_Pt","Eta_Pt", 160, 0, 80);
    fHistEta_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEta_Pt->GetYaxis()->SetTitle("N_{#eta^{all}}");
    fHistEta_Pt->Sumw2();
    fOutputContainer->Add(fHistEta_Pt);

    fHistEtaPrime_InAcc_EMC     = new TH1F("EtaPrime_Pt_InEMCAcc","EtaPrime_Pt_InEMCAcc", 160, 0, 80);
    fHistEtaPrime_InAcc_EMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEtaPrime_InAcc_EMC->GetYaxis()->SetTitle("A_{EMC} * N_{#eta'}");
    fHistEtaPrime_InAcc_EMC->Sumw2();
    fOutputContainer->Add(fHistEtaPrime_InAcc_EMC);

    fHistEtaPrime_InAcc_PCMECM  = new TH1F("EtaPrime_Pt_InPCM-EMCAcc","EtaPrime_Pt_InPCM-EMCAcc", 160, 0, 80);
    fHistEtaPrime_InAcc_PCMECM->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistEtaPrime_InAcc_PCMECM->GetYaxis()->SetTitle("A_{PCM-EMC} * N_{#eta'}");
    fHistEtaPrime_InAcc_PCMECM->Sumw2();
    fOutputContainer->Add(fHistEtaPrime_InAcc_PCMECM);

    PostData(1, fOutputContainer);

    fTreeEtaPrimeKinematics     = new TTree("EtaPrimeKinematics","EtaPrimeKinematics");
    
    fTreeEtaPrimeKinematics->Branch( "Gamma1_pt",       &fBuffer_Gamma1_pt,         "Gamma1_pt/s");
    fTreeEtaPrimeKinematics->Branch( "Gamma1_eta",      &fBuffer_Gamma1_eta,        "Gamma1_eta/S");
    fTreeEtaPrimeKinematics->Branch( "Gamma1_phi",      &fBuffer_Gamma1_phi,        "Gamma1_phi/s");
    fTreeEtaPrimeKinematics->Branch( "Gamma1_y",        &fBuffer_Gamma1_y,          "Gamma1_y/S");
    fTreeEtaPrimeKinematics->Branch( "Gamma1_E",        &fBuffer_Gamma1_E,          "Gamma1_E/s");

    fTreeEtaPrimeKinematics->Branch( "Gamma2_pt",       &fBuffer_Gamma2_pt,         "Gamma2_pt/s");
    fTreeEtaPrimeKinematics->Branch( "Gamma2_eta",      &fBuffer_Gamma2_eta,        "Gamma2_eta/S");
    fTreeEtaPrimeKinematics->Branch( "Gamma2_phi",      &fBuffer_Gamma2_phi,        "Gamma2_phi/s");
    fTreeEtaPrimeKinematics->Branch( "Gamma2_y",        &fBuffer_Gamma2_y,          "Gamma2_y/S");
    fTreeEtaPrimeKinematics->Branch( "Gamma2_E",        &fBuffer_Gamma2_E,          "Gamma2_E/s");

    fTreeEtaPrimeKinematics->Branch( "PiPl_pt",         &fBuffer_PiPl_pt,           "PiPl_pt/s");
    fTreeEtaPrimeKinematics->Branch( "PiPl_E",          &fBuffer_PiPl_E,            "PiPl_E/s");
    fTreeEtaPrimeKinematics->Branch( "PiPl_eta",        &fBuffer_PiPl_eta,          "PiPl_eta/S");
    fTreeEtaPrimeKinematics->Branch( "PiPl_phi",        &fBuffer_PiPl_phi,          "PiPl_phi/s");
    fTreeEtaPrimeKinematics->Branch( "PiPl_y",          &fBuffer_PiPl_y,            "PiPl_y/S");

    fTreeEtaPrimeKinematics->Branch( "PiMi_pt",         &fBuffer_PiMi_pt,           "PiMi_pt/s");
    fTreeEtaPrimeKinematics->Branch( "PiMi_E",          &fBuffer_PiMi_E,            "PiMi_E/s");
    fTreeEtaPrimeKinematics->Branch( "PiMi_eta",        &fBuffer_PiMi_eta,          "PiMi_eta/S");
    fTreeEtaPrimeKinematics->Branch( "PiMi_phi",        &fBuffer_PiMi_phi,          "PiMi_phi/s");
    fTreeEtaPrimeKinematics->Branch( "PiMi_y",          &fBuffer_PiMi_y,            "PiMi_y/S");

    fTreeEtaPrimeKinematics->Branch( "Eta_pt",          &fBuffer_Eta_pt,            "Eta_pt/s");
    fTreeEtaPrimeKinematics->Branch( "Eta_E",           &fBuffer_Eta_E,             "Eta_E/s");
    fTreeEtaPrimeKinematics->Branch( "Eta_eta",         &fBuffer_Eta_eta,           "Eta_eta/S");
    fTreeEtaPrimeKinematics->Branch( "Eta_phi",         &fBuffer_Eta_phi,           "Eta_phi/s");
    fTreeEtaPrimeKinematics->Branch( "Eta_y",           &fBuffer_Eta_y,             "Eta_y/S");
    
    fTreeEtaPrimeKinematics->Branch( "EtaPrime_pt",     &fBuffer_EtaPrime_pt,       "EtaPrime_pt/s");
    fTreeEtaPrimeKinematics->Branch( "EtaPrime_E",      &fBuffer_EtaPrime_E,        "EtaPrime_E/s");
    fTreeEtaPrimeKinematics->Branch( "EtaPrime_eta",    &fBuffer_EtaPrime_eta,      "EtaPrime_eta/S");
    fTreeEtaPrimeKinematics->Branch( "EtaPrime_phi",    &fBuffer_EtaPrime_phi,      "EtaPrime_phi/s");
    fTreeEtaPrimeKinematics->Branch( "EtaPrime_y",      &fBuffer_EtaPrime_y,        "EtaPrime_y/S");

    OpenFile(2);
    PostData(2, fTreeEtaPrimeKinematics);
}

void AliAnalysisTaskEtaPrimeMCStudies::UserExec(Option_t *){
    
    fInputEvent     = InputEvent();
    fMCEvent        = MCEvent();
    if( fMCEvent == nullptr ) return;

    const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
    Double_t mcProdVtxZ   = primVtxMC->GetZ();


    if (TMath::Abs(mcProdVtxZ) < 10 ){
      fHistNEvents->Fill(0);
    } else {
      fHistNEvents->Fill(1);
    }

    AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
    AliGenPythiaEventHeader *pyH  = dynamic_cast<AliGenPythiaEventHeader*>(mcEH);
    if( !pyH ) return;

    Float_t         ntrials = 0;
    ntrials         = pyH->Trials();
    if(ntrials)     fHistNEvents->Fill(2,ntrials);

    Double_t    xSection = 0;
    Double_t    ptHard = 0;
    xSection    = pyH->GetXsection();
    ptHard      = pyH->GetPtHard();
    if( xSection )  fHistXSection->Fill(xSection);
    if( ptHard )    fHistPtHard->Fill(ptHard);

    ProcessMCParticles();

    PostData( 1, fOutputContainer );
}

void AliAnalysisTaskEtaPrimeMCStudies::ProcessMCParticles(){
    for(Long_t i=0; i<fMCEvent->GetNumberOfTracks(); i++){
        // fill primary histograms
        AliMCParticle* etaPrime     = nullptr;
        etaPrime                    = (AliMCParticle*)(fMCEvent->GetTrack(i));
        if( !etaPrime )                     continue;
        if( etaPrime->PdgCode() == 221 )    fHistEta_Pt->Fill(etaPrime->Pt());
        if( etaPrime->PdgCode() == 111 )    fHistPi0_Pt->Fill(etaPrime->Pt()); 
        if( etaPrime->PdgCode() != 331 )    continue;
        fHistEtaPrimeAll_Pt->Fill( etaPrime->Pt() );

        AliMCParticle* piPlus       = nullptr;
        AliMCParticle* piMinus      = nullptr;
        AliMCParticle* etaMeson     = nullptr;
        AliMCParticle* gamma1       = nullptr;
        AliMCParticle* gamma2       = nullptr;

        Bool_t  etaToGammaGamma     = kFALSE;

        for(Int_t index=etaPrime->GetDaughterFirst(); index<=etaPrime->GetDaughterLast(); index++){
            AliMCParticle* temp     = (AliMCParticle*)fMCEvent->GetTrack(index);
            switch( temp->PdgCode()) {
                case 211:
                    piPlus      = temp;
                    break;
                case -211:
                    piMinus     = temp;
                    break;
                case 221:
                    etaMeson    = temp;
                    gamma1  = (AliMCParticle*)(fMCEvent->GetTrack( temp->GetDaughterFirst() ));
                    gamma2  = (AliMCParticle*)(fMCEvent->GetTrack( temp->GetDaughterLast() ));
                    if( etaMeson->GetNDaughters() == 2 && gamma1->PdgCode() == 22 && gamma2->PdgCode() == 22) etaToGammaGamma = kTRUE;
                    break;
            }
        }

        if( !(etaMeson && piPlus && piMinus && etaToGammaGamma) ) continue; // skip if eta prime doesn't decay into pi+ pi- eta

        fHistEtaPrime_Pt->Fill(etaPrime->Pt());
        fBuffer_EtaPrime_pt     = static_cast<Short_t>  (etaPrime->Pt()     * 1000);
        fBuffer_EtaPrime_E      = static_cast<Short_t>  (etaPrime->E()      * 1000);
        fBuffer_EtaPrime_eta    = static_cast<UShort_t> (etaPrime->Eta()    * 1000);
        fBuffer_EtaPrime_phi    = static_cast<Short_t>  (etaPrime->Phi()    * 1000);
        fBuffer_EtaPrime_y      = static_cast<UShort_t> (etaPrime->Y()      * 1000);

        fBuffer_PiPl_pt     = static_cast<Short_t>  (piPlus->Pt()     * 1000);
        fBuffer_PiPl_E      = static_cast<Short_t>  (piPlus->E()      * 1000);
        fBuffer_PiPl_eta    = static_cast<UShort_t> (piPlus->Eta()    * 1000);
        fBuffer_PiPl_phi    = static_cast<Short_t>  (piPlus->Phi()    * 1000);
        fBuffer_PiPl_y      = static_cast<UShort_t> (piPlus->Y()      * 1000);
        fHistPiPl_Pt->Fill( piPlus->Pt() );

        fBuffer_PiMi_pt     = static_cast<Short_t>  (piMinus->Pt()     * 1000);
        fBuffer_PiMi_E      = static_cast<Short_t>  (piMinus->E()      * 1000);
        fBuffer_PiMi_eta    = static_cast<UShort_t> (piMinus->Eta()    * 1000);
        fBuffer_PiMi_phi    = static_cast<Short_t>  (piMinus->Phi()    * 1000);
        fBuffer_PiMi_y      = static_cast<UShort_t> (piMinus->Y()      * 1000);
        fHistPiMi_Pt->Fill( piMinus->Pt() );

        fBuffer_Eta_pt      = static_cast<Short_t>  (etaMeson->Pt()     * 1000);
        fBuffer_Eta_E       = static_cast<Short_t>  (etaMeson->E()      * 1000);
        fBuffer_Eta_eta     = static_cast<UShort_t> (etaMeson->Eta()    * 1000);
        fBuffer_Eta_phi     = static_cast<Short_t>  (etaMeson->Phi()    * 1000);
        fBuffer_Eta_y       = static_cast<UShort_t> (etaMeson->Y()      * 1000);
        fHistEtaFromEtaPrime_Pt->Fill( etaMeson->Pt() );

        fBuffer_Gamma1_pt       = static_cast<Short_t>  (gamma1->Pt()   * 1000);
        fBuffer_Gamma1_E        = static_cast<Short_t>  (gamma1->E()    * 1000);
        fBuffer_Gamma1_eta      = static_cast<UShort_t> (gamma1->Eta()  * 1000);
        fBuffer_Gamma1_phi      = static_cast<Short_t>  (gamma1->Phi()  * 1000);
        fBuffer_Gamma1_y        = static_cast<UShort_t> (gamma1->Y()    * 1000);
        fHistGamma1_Pt->Fill( gamma1->Pt() );

        fBuffer_Gamma2_pt       = static_cast<Short_t>  (gamma2->Pt()   * 1000);
        fBuffer_Gamma2_E        = static_cast<Short_t>  (gamma2->E()    * 1000);
        fBuffer_Gamma2_eta      = static_cast<UShort_t> (gamma2->Eta()  * 1000);
        fBuffer_Gamma2_phi      = static_cast<Short_t>  (gamma2->Phi()  * 1000);
        fBuffer_Gamma2_y        = static_cast<UShort_t> (gamma2->Y()    * 1000);
        fHistGamma2_Pt->Fill( gamma2->Pt() );

        fTreeEtaPrimeKinematics->Fill();

        if( !IsEtaSelected(etaMeson) || !IsPionSelected(piPlus) || !IsPionSelected(piMinus) )   continue;
        if( IsPhotonInEMCalAcceptance(gamma1) && IsPhotonInEMCalAcceptance(gamma2) )            
            fHistEtaPrime_InAcc_EMC->Fill(etaPrime->Pt());
        if( (IsPhotonInEMCalAcceptance(gamma1) && IsPhotonInPCMAcceptance(gamma2)) || (IsPhotonInPCMAcceptance(gamma1) && IsPhotonInEMCalAcceptance(gamma2)) ) 
            fHistEtaPrime_InAcc_PCMECM->Fill(etaPrime->Pt());
    }

    return;
}

//________________________________________________________________________
bool AliAnalysisTaskEtaPrimeMCStudies::IsPhotonInPCMAcceptance(AliMCParticle* part) const {
  const Double_t kBoundaryEta = 0.900001;
  if (part->Pt() > 0.050 && TMath::Abs(part->Eta()) < kBoundaryEta) return true;

  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskEtaPrimeMCStudies::IsPhotonInEMCalAcceptance(AliMCParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.6687;
  const Double_t kBoundaryEtaMax = 0.66465;
  const Double_t kBoundaryPhiMin = 1.39626;
  const Double_t kBoundaryPhiMax = 3.15;
  if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskEtaPrimeMCStudies::IsEtaSelected(AliMCParticle* part) const {   // based on the default cuts for eta prime analysis
    const Double_t  rapidityMin = -0.75;
    const Double_t  rapidityMax = 0.75;
    const Double_t  minPt       = 2.5;

    if( part->Y() < rapidityMin ) return false;
    if( part->Y() > rapidityMax ) return false;
    if( part->Pt() < minPt      ) return false;
    return true;
}

bool AliAnalysisTaskEtaPrimeMCStudies::IsPionSelected(AliMCParticle* part) const { // based on the default cuts for eta prime analysis
    const Double_t  etaMin  = -0.9;
    const Double_t  etaMax  = 0.9;

    if( part->Eta() > etaMax ) return false;
    if( part->Eta() < etaMin ) return false;
    return true;
}



void AliAnalysisTaskEtaPrimeMCStudies::Terminate(const Option_t *){
///Grid
}