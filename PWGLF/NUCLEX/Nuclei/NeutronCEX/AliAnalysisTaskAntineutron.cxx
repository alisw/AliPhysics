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

/* AliAnaysisTaskAntineutron
 *
 * Analysis for the detection of anti-neutrons
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TList.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDpid.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskAntineutron.h"

class AliAnalysisTaskAntineutron;    

using namespace std;           

ClassImp(AliAnalysisTaskAntineutron) 

AliAnalysisTaskAntineutron::AliAnalysisTaskAntineutron() : AliAnalysisTaskSE(),
    //Trees
    fTree_proton(0),
    fTree_Antip(0),
    fTree_pair(0),

    fESD(0),
    fOutputList(0),
    fMCevent(0),
    fSimulation(0),
    fPIDResponse(0),
    fESDpid(0),
    fHistPt(0),
    fHistNEvents(0),
    fITSmap(0),
    // MC Antineutron histograms
    fHistAntinEta(0),
    fHistAntinEk(0),
    fHistAntinEk_ITScuts(0),
    fHistAntipEk(0),
    fHistpEk(0),
    //MC CEX histograms
    fHistcexPx(0),
    fHistcexPy(0),
    fHistcexPz(0),
    fHistcexEk(0),
    fHistcexP(0), 
    fHistcexDP(0),
    fHistcexVrtx(0),
    fHistmomcexpdg(0),
    fHistcexEta_antin(0),
    fHistcex_ITScuts(0),
    fHistcex_angle(0),
    // MC CEX histograms normalized
    fHistncexPx(0),
    fHistncexPy(0),
    fHistncexPz(0),
    fHistncexEk(0),
    fHistncexP(0),
    fHistncex_ITScuts(0),
    
    // Data
    fHistRecoAngle(0),
    fHistRecoP(0),
    fHistRecoradio(0),
    fHistRecovtx(0),
    fHistDCApair(0),
    fHistRecoVtxSep(0),
    fHistRecosvtxXY(0),
    fHistRecosvtxZ(0)
{

}
//_____________________________________________________________________________
AliAnalysisTaskAntineutron::AliAnalysisTaskAntineutron(const char* name) : AliAnalysisTaskSE(name),
    //Trees
    fTree_proton(0),
    fTree_Antip(0),
    fTree_pair(0),

    fESD(0),
    fOutputList(0),
    fMCevent(0),
    fSimulation(0),
    fPIDResponse(0),
    fESDpid(0),
    fHistPt(0),
    fHistNEvents(0),
    fITSmap(0),
    // MC Antineutron histograms
    fHistAntinEta(0),
    fHistAntinEk(0),
    fHistAntinEk_ITScuts(0),
    fHistAntipEk(0),
    fHistpEk(0),
    // MC CEX histograms
    fHistcexPx(0), 
    fHistcexPy(0),
    fHistcexPz(0),
    fHistcexEk(0),
    fHistcexP(0), 
    fHistcexDP(0),
    fHistcexVrtx(0),
    fHistmomcexpdg(0),
    fHistcexEta_antin(0),
    fHistcex_ITScuts(0),
    fHistcex_angle(0),
    // MC CEX histograms normalized
    fHistncexPx(0),
    fHistncexPy(0),
    fHistncexPz(0),
    fHistncexEk(0),
    fHistncexP(0),
    fHistncex_ITScuts(0),
    
    // Data
    fHistRecoAngle(0),
    fHistRecoP(0),
    fHistRecoradio(0),
    fHistRecovtx(0),
    fHistDCApair(0),
    fHistRecoVtxSep(0),
    fHistRecosvtxXY(0),
    fHistRecosvtxZ(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                   
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        
    fSimulation = kTRUE;	        // Is a simulation?                                     
                                        
}
//_____________________________________________________________________________
AliAnalysisTaskAntineutron::~AliAnalysisTaskAntineutron()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskAntineutron::UserCreateOutputObjects()
{
    // Output objects

    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);      

    //Trees
    fTree_Antip = new TTree("fTree_Antip", "Datos MC Antiproton");
    fTree_Antip->Branch("fmom_pdg", &fmom_pdg, "fmom_pdg/I");
    fTree_Antip->Branch("fmom_process", &fmom_process, "fmom_process/I");
    fTree_Antip->Branch("fmom_E", &fmom_E, "fmom_E/D");
    fTree_Antip->Branch("fAntipP", &fAntipP, "fAntipP/D");
    fTree_Antip->Branch("fAntip_Px", &fAntip_Px, "fAntip_Px/D");
    fTree_Antip->Branch("fAntip_Py", &fAntip_Py, "fAntip_Py/D");
    fTree_Antip->Branch("fAntip_Pz", &fAntip_Pz, "fAntip_Pz/D");
    fTree_Antip->Branch("fAntipPt", &fAntipPt, "fAntipPt/D");
    fTree_Antip->Branch("fAntipE", &fAntipE, "fAntipE/D");
    fTree_Antip->Branch("fAntipY", &fAntipY, "fAntipY/D");
    fTree_Antip->Branch("fAntipEta", &fAntipEta, "fAntipEta/D");
    fTree_Antip->Branch("fAntipTPCsignal", &fAntipTPCsignal, "fAntipTPCsignal/D");
    fTree_Antip->Branch("fAntipITSsignal", &fAntipITSsignal, "fAntipITSsignal/D");
    fTree_Antip->Branch("fAntipTPCpoints", &fAntipTPCpoints, "fAntipTPCpoints/D");
    fTree_Antip->Branch("fAntipTOFsignal", &fAntipTOFsignal, "fAntipTOFsignal/D");
    fTree_Antip->Branch("fAntipncTPC", &fAntipncTPC, "fAntipncTPC/I");
    fTree_Antip->Branch("fAntipncITS", &fAntipncITS, "fAntipncITS/I");
    fTree_Antip->Branch("fAntipDCAxy", &fAntipDCAxy, "fDCAxy/D");
    fTree_Antip->Branch("fAntipDCAz", &fAntipDCAz, "fAntipDCAz/D");
    
    fTree_proton = new TTree("fTree_proton", "Datos MC proton");
    fTree_proton->Branch("fpmom_pdg", &fpmom_pdg, "fpmom_pdg/I");
    fTree_proton->Branch("fpmom_process", &fpmom_process, "fpmom_process/I");
    fTree_proton->Branch("fpmom_E", &fpmom_process, "fpmom_E/D");
    fTree_proton->Branch("fpP", &fpP, "fpP/D");
    fTree_proton->Branch("fp_Px", &fp_Px, "fp_Px/D");
    fTree_proton->Branch("fp_Py", &fp_Py, "fp_Py/D");
    fTree_proton->Branch("fp_Pz", &fp_Pz, "fp_Pz/D");
    fTree_proton->Branch("fpPt", &fpPt, "fpPt/D");
    fTree_proton->Branch("fpE", &fpE, "fpE/D");
    fTree_proton->Branch("fpY", &fpY, "fpY/D");
    fTree_proton->Branch("fpEta", &fpEta, "fpEta/D");
    fTree_proton->Branch("fpTPCsignal", &fpTPCsignal, "fpTPCsignal/D");
    fTree_proton->Branch("fpITSsignal", &fpITSsignal, "fpITSsignal/D");
    fTree_proton->Branch("fpTPCpoints", &fpTPCpoints, "fpTPCpoints/D");
    fTree_proton->Branch("fpTOFsignal", &fpTOFsignal, "fpTOFsignal/D");
    fTree_proton->Branch("fpncTPC", &fpncTPC, "fpncTPC/I");
    fTree_proton->Branch("fpncITS", &fpncITS, "fpncITS/I");
    fTree_proton->Branch("fpDCAxy", &fpDCAxy, "fpDCAxy/D");
    fTree_proton->Branch("fpDCAz", &fpDCAz, "fpDCAz/D");
    
    fTree_pair= new TTree("fTree_pair", "Datos track par");
    fTree_pair->Branch("fPairPSumMag", &fPairPSumMag, "fPairPSumMag/D");
    fTree_pair->Branch("fPairVtxSep", &fPairVtxSep, "fPairVtxSep/D");
    fTree_pair->Branch("fPairMinDCA", &fPairMinDCA, "fPairMinDCA/D");
    fTree_pair->Branch("fPairSecRadiu", &fPairSecRadius, "fPairSecRadiu/D");
    fTree_pair->Branch("fPairVtxDistance", &fPairVtxDistance, "fPairVtxDistance/D");
    fTree_pair->Branch("fSecVtxX", &fSecVtxX, "fSecVtxX/D");
    fTree_pair->Branch("fSecVtxY", &fSecVtxY, "fSecVtxY/D");
    fTree_pair->Branch("fSecVtxZ", &fSecVtxZ, "fSecVtxZ/D");
    fTree_pair->Branch("fPrimVtxX", &fPrimVtxX, "fPrimVtxX/D");
    fTree_pair->Branch("fPrimVtxY", &fPrimVtxY, "fPrimVtxY/D");
    fTree_pair->Branch("fPrimVtxZ", &fPrimVtxZ, "fPrimVtxZ/D");
    
    // Histograms
    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);
    fHistNEvents = new TH1F("fHistNEvents", "fHistNEvents", 100, 0, 100000000);
    // MC Antineutron histograms
    fHistAntinEta = new TH1F("fHistAntinEta", "fHistAntinEta", 200, -10, 10);
    fHistAntinEk = new TH1F("fHistAntinEk", "fHistAntinEk", 100, 0, 10);
    fHistAntinEk_ITScuts = new TH1F("fHistAntinEk_ITScuts", "fHistAntinEk_ITScuts", 100, 0, 10);
    fHistAntipEk = new TH1F("fHistAntipEk", "fHistAntipEk", 100, 0, 10);
    fHistpEk = new TH1F("fHistpEk", "fHistpEk", 100, 0, 10);
    // MC CEX histograms
    fHistcexPx = new TH2F("fHistcexPx", " ; Px_rec; Px_antin", 1000, -5, 5, 1000, -5, 5);
    fHistcexPy = new TH2F("fHistcexPy", " ; Py_rec; Py_antin", 1000, -5, 5, 1000, -5, 5);
    fHistcexPz = new TH2F("fHistcexPz", " ; Pz_rec; Pz_antin", 1000, -5, 5, 1000, -5, 5);
    fHistcexEk = new TH2F("fHistcexEk", " ; Ek_rec; Ek_antin", 100, 0, 10, 100, 0, 10);
    fHistcexP = new TH2F("fHistcexP", " ; P_rec; P_antin", 500, 0, 5, 500, 0, 5);
    fHistcexDP = new TH1F("fHistcexDP", "fHistcexDP", 100, 0, 10);
    fHistcexVrtx = new TH2F("fHistcexVrtx", " ; V_x; V_y", 8700, -43.5, 43.5, 8700, -43.5, 43.5);
    fHistmomcexpdg = new TH1F("fHistmomcexpdg", "fHistmomcexpdg", 6000, -3000, 3000);
    fHistcexEta_antin = new TH1F("fHistcexEta_antin", "fHistcexEta_antin", 200, -10, 10);
    fHistcex_ITScuts = new TH1F("fHistcex_ITScuts", "fHistcex_ITScuts", 200, -10, 10);
    fHistcex_angle = new TH1F("fHistcex_angle", "fHistcex_angle", 360, 0, 360);
    // MC CEX histograms normalized
    fHistncexPx = new TH1F("fHistncexPx", "fHistncexPx", 1000, -5, 5);
    fHistncexPy = new TH1F("fHistncexPy", "fHistncexPy", 1000, -5, 5);
    fHistncexPz = new TH1F("fHistncexPz", "fHistncexPz", 1000, -5, 5);
    fHistncexEk = new TH1F("fHistncexEk", "fHistncexEk", 1000, -5, 5);
    fHistncexP = new TH1F("fHistncexP", "fHistncexP", 500, 0, 5);
    fHistncex_ITScuts = new TH1F("fHistncex_ITScuts", "fHistncex_ITScuts", 200, -10, 10);
   
    //Data
    fHistRecoAngle = new TH1F("fHistRecoAngle", "fHistRecoAngle", 360, 0, 360);
    fHistRecoP = new TH1F("fHistRecoP", "fHistRecoP", 500, 0, 5);
    fHistRecoradio = new TH1F("fHistRecoradio", "fHistRecoradio",440, 0, 43.5);
    fHistRecovtx = new TH1F("fHistRecovtx", "fHistRecovtx",440, 0, 43.5);
    fHistDCApair = new TH1F("fHistDCApair", "fHistDCApair",500, 0, 50);
    fHistRecoVtxSep = new TH1F("fHistRecoVtxSep", "Distancia entre vértices individuales;#Delta r_{vtx} (cm);Entries", 100, 0, 5);
    fHistRecosvtxXY = new TH2F("fHistRecosvtxXY", "Vértice secundario estimado XY;x (cm);y (cm)", 8700, -43.5, 43.5, 8700, -43.5, 43.5);
    fHistRecosvtxZ = new TH1F("fHistRecosvtxZ", "Coordenada Z del vértice secundario;z (cm);Entries", 100, -50, 50);
    
    // Save in the output
    fOutputList->Add(fTree_proton);
    fOutputList->Add(fTree_Antip);
    fOutputList->Add(fTree_pair);
    fOutputList->Add(fHistPt);
    fOutputList->Add(fHistNEvents);
    // MC Antineutron histograms
    fOutputList->Add(fHistAntinEta);
    fOutputList->Add(fHistAntinEk);
    fOutputList->Add(fHistAntinEk_ITScuts);
    fOutputList->Add(fHistAntipEk);
    fOutputList->Add(fHistpEk);
    // MC CEX histograms
    fOutputList->Add(fHistcexPx); 
    fOutputList->Add(fHistcexPy);
    fOutputList->Add(fHistcexPz);
    fOutputList->Add(fHistcexEk);
    fOutputList->Add(fHistcexP);
    fOutputList->Add(fHistcexDP);
    fOutputList->Add(fHistcexVrtx);
    fOutputList->Add(fHistmomcexpdg);
    fOutputList->Add(fHistcexEta_antin);
    fOutputList->Add(fHistcex_ITScuts);
    fOutputList->Add(fHistcex_angle);
    // MC CEX histograms normalized
    fOutputList->Add(fHistncexPx);
    fOutputList->Add(fHistncexPy);
    fOutputList->Add(fHistncexPz);
    fOutputList->Add(fHistncexEk);
    fOutputList->Add(fHistncexP);
    fOutputList->Add(fHistncex_ITScuts);

    //data
    fOutputList->Add(fHistRecoAngle);
    fOutputList->Add(fHistRecoP);
    fOutputList->Add(fHistRecoradio);
    fOutputList->Add(fHistRecovtx);
    fOutputList->Add(fHistDCApair);
    fOutputList->Add(fHistRecoVtxSep);
    fOutputList->Add(fHistRecosvtxXY);
    fOutputList->Add(fHistRecosvtxZ);
    
    PostData(1, fOutputList);
}

//_____________________________________________________________________________
void AliAnalysisTaskAntineutron::UserExec(Option_t *)
{
    // this function is called once for each event
    
    fEventHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (fEventHandler == 0) {
        AliError("Could not get InputHandler");
        return;
    }
    
    // Monte Carlo
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliMCEventHandler *mch = dynamic_cast<AliMCEventHandler*> (man->GetMCtruthEventHandler());
    if (mch == 0){
        Printf("ERROR: Could not retrieve MC event handler");
        return;
    }
    fMCevent = mch->MCEvent();
    if(fMCevent == 0) {
        Printf("ERROR: Could not retrieve MC event");
        return;
    }
    Int_t nParticles = 0;
    AliStack* stack = fMCevent->Stack();
    if (!stack)
    {
        AliDebug(AliLog::kWarning, "stack not available");
        return;
    }
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(fESD == 0){
        Printf("ERROR: Could not retrieve ESD event");
        return;
    }
    
    //AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliESDInputHandler *esdh = dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler());
    if (esdh == 0){
        Printf("ERROR: Could not retrieve ESD event handler");
        return;
    }
    
// Select Antineutrons and cex events
       for (Int_t i = 0; i < fMCevent->GetNumberOfTracks(); ++i)
       {
           TParticle* iParticle = stack->Particle(i);
           if(!iParticle) continue;
           Int_t pdg = iParticle->GetPdgCode();
           
           // Select primary antineutrons
           if(stack->IsPhysicalPrimary(i) && pdg == -2112)
           {
               fHistAntinEta->Fill(iParticle->Eta());
               fHistAntinEk->Fill(iParticle->Ek());
               fHistAntipEk->Fill(iParticle->Px());
               fHistpEk->Fill(iParticle->Py());
               if (TMath::Abs(iParticle->Eta()) < 1.5 && TMath::Abs(iParticle->Vz())< 5.3) fHistAntinEk_ITScuts->Fill(iParticle->Ek());
           }
           
           // Select seconday antiprotons from material
           if(stack->IsSecondaryFromMaterial(i) && pdg == -2212)
           {
               Int_t imom = iParticle->GetFirstMother();
               if(imom < 0) return;
               TParticle* mom = stack->Particle(imom);
               Int_t pdgmom = mom->GetPdgCode();
               //if (pdgmom != -2112) continue;
               if (!stack->IsPhysicalPrimary(imom)) continue;
               
               Double_t PVx = mom->Vx();
               Double_t PVy = mom->Vy();
               Double_t PVz = mom->Vz();
               Double_t mPx = mom->Px();
               Double_t mPy = mom->Py();
               Double_t mPz = mom->Pz();
               Double_t mEk = mom->Ek();
               Double_t mP = TMath::Sqrt(TMath::Power(mPx,2)+TMath::Power(mPy,2)+TMath::Power(mPz,2));
               
               Double_t fAntipPx = iParticle->Px();
               Double_t fAntipPy = iParticle->Py();
               Double_t fAntipPz = iParticle->Pz();
               Double_t fAntipEk = iParticle->Ek();
               Double_t fAntipVx = iParticle->Vx();
               Double_t fAntipVy = iParticle->Vy();
               Double_t fAntipVz = iParticle->Vz();
               if (3.9<=TMath::Abs(fAntipVx) && TMath::Abs(fAntipVx)<=43.6 && 3.9<=TMath::Abs(fAntipVy) && TMath::Abs(fAntipVy)<=43.6 && TMath::Abs(fAntipVz)<=48.9)
               {
                   // Pion veto
                   Bool_t pion = kFALSE;
                   for (Int_t k = 0; k < fMCevent->GetNumberOfTracks(); ++k)
                   {
                       TParticle* kParticle = stack->Particle(k);
                       if(!kParticle) continue;
                       Int_t kpdg = kParticle->GetPdgCode();
                       
                       if(stack->IsSecondaryFromMaterial(k) && kpdg == 211)
                       {
                           Int_t kmom = kParticle->GetFirstMother();
                           if(kmom < 0) return;
                           TParticle* momk = stack->Particle(kmom);
                           Int_t kpdgmom = momk->GetPdgCode();
                           //if (kpdgmom != -2112) continue;
                           Double_t fpiVx = kParticle->Vx();
                           Double_t fpiVy = kParticle->Vy();
                           Double_t fpiVz = kParticle->Vz();
                           if(fAntipVx == fpiVx && fAntipVy == fpiVy && fAntipVz == fpiVz) pion = kTRUE;
                       }
                   }
                   
                   Bool_t pionm = kFALSE;
                   for (Int_t k = 0; k < fMCevent->GetNumberOfTracks(); ++k)
                   {
                       TParticle* kParticle = stack->Particle(k);
                       if(!kParticle) continue;
                       Int_t kpdg = kParticle->GetPdgCode();
                       
                       if(stack->IsSecondaryFromMaterial(k) && kpdg == -211)
                       {
                           Int_t kmom = kParticle->GetFirstMother();
                           if(kmom < 0) return;
                           TParticle* momk = stack->Particle(kmom);
                           Int_t kpdgmom = momk->GetPdgCode();
                           //if (kpdgmom != -2112) continue;
                           Double_t fpiVx = kParticle->Vx();
                           Double_t fpiVy = kParticle->Vy();
                           Double_t fpiVz = kParticle->Vz();
                           if(fAntipVx == fpiVx && fAntipVy == fpiVy && fAntipVz == fpiVz) pionm = kTRUE;
                       }
                   }
                   
                   if (pion == kFALSE && pionm == kFALSE)
                   {
                       // CEX
                       Double_t dplane = 10;
                       Double_t dplane_tmp = 0;
                       Double_t P = 0;
                       Double_t P_tmp = 0;
                       Double_t Px_p = 0;
                       Double_t Py_p = 0;
                       Double_t Pz_p = 0;
                       Double_t Ek = 0;
                       Double_t Ek_tmp = 0;
                       Double_t Ek_p = 0;
                       Int_t k_plane = 0;
                       Int_t k_Ek = 0;
                       Int_t k_P = 0;
                       
                       //Select secondary protons from material
                       for (Int_t k = 0; k < fMCevent->GetNumberOfTracks(); ++k)
                       {
                           TParticle* kParticle = stack->Particle(k);
                           if(!kParticle) continue;
                           Int_t kpdg = kParticle->GetPdgCode();
                           
                           if(stack->IsSecondaryFromMaterial(k) && kpdg == 2212)
                           {
                               Int_t kmom = kParticle->GetFirstMother();
                               if(kmom < 0) return;
                               TParticle* momk = stack->Particle(kmom);
                               Int_t kpdgmom = momk->GetPdgCode();
                               //if (kpdgmom != -2112) continue;
                               
                               Double_t fpPx = kParticle->Px();
                               Double_t fpPy = kParticle->Py();
                               Double_t fpPz = kParticle->Pz();
                               Double_t fpEk = kParticle->Ek();
                               Double_t fpVx = kParticle->Vx();
                               Double_t fpVy = kParticle->Vy();
                               Double_t fpVz = kParticle->Vz();
                               //CEX proton selection
                               if(fAntipVx == fpVx && fAntipVy == fpVy && fAntipVz == fpVz)
                               {
                                   dplane_tmp = (fpPy*fAntipPz - fpPz*fAntipPy)*(PVx-fAntipVx) + (fpPz*fAntipPx - fpPx*fAntipPz)*(PVy-fAntipVy)
                                   + (fpPx*fAntipPy - fpPy*fAntipPx)*(PVz-fAntipVz);
                                   if(TMath::Abs(dplane_tmp) < TMath::Abs(dplane))
                                   {
                                       k_plane = k;
                                       dplane = dplane_tmp;
                                   }
                                   
                                   Ek_tmp = fAntipEk + fpEk;
                                   if(TMath::Abs(Ek_tmp) > TMath::Abs(Ek))
                                   {
                                       k_Ek = k;
                                       Ek = Ek_tmp;
                                       Ek_p = fpEk;
                                   }
                                   
                                   P_tmp = TMath::Sqrt(TMath::Power((fpPx+fAntipPx),2)+TMath::Power((fpPy+fAntipPy),2)+TMath::Power((fpPz+fAntipPz),2));
                                   if(TMath::Abs(P_tmp) > TMath::Abs(P))
                                   {
                                       k_P = k;
                                       P = P_tmp;
                                       Px_p = fpPx;
                                       Py_p = fpPy;
                                       Pz_p = fpPz;
                                   }
                               }
                           }
                       }
                       if(k_plane == k_Ek && k_plane == k_P && k_plane != 0 )
                       {
                           fHistmomcexpdg->Fill(pdgmom);
                           //antin mother
                           if (pdgmom != -2112) continue;
                           if (TMath::Abs(mom->Eta()) < 1.5 && TMath::Abs(mom->Vz())< 5.3)
                           {
                               fHistcexPx->Fill(fAntipPx + Px_p, mPx);
                               fHistcexPy->Fill(fAntipPy + Py_p, mPy);
                               fHistcexPz->Fill(fAntipPz + Pz_p, mPz);
                               fHistcexP->Fill(P, mP);
                               fHistcexEk->Fill(fAntipEk + Ek_p, mEk);
                               fHistcexDP->Fill(TMath::Abs(dplane));
                               
                               //Angle
                               TVector3 pVecProton = TVector3(Px_p, Py_p, Pz_p);
                               TVector3 pVecAntiproton= TVector3(fAntipPx, fAntipPy, fAntipPz);
                               Double_t angleRad = pVecProton.Angle(pVecAntiproton);
                               Double_t angleDeg = angleRad * TMath::RadToDeg();
                               fHistcex_angle->Fill(angleDeg);
                               
                               if (TMath::Abs(mom->Eta()) < 0.9 && TMath::Abs(mom->Vz())< 5.3) fHistcex_ITScuts->Fill(P);
                               
                               fHistncexPx->Fill((fAntipPx + Px_p)/mPx);
                               fHistncexPy->Fill((fAntipPy + Py_p)/mPy);
                               fHistncexPz->Fill((fAntipPz + Pz_p)/mPz);
                               fHistncexP->Fill(P/mP);
                               fHistncexEk->Fill((fAntipEk + Ek_p)/mEk);
                               if (TMath::Abs(mom->Eta()) < 0.9 && TMath::Abs(mom->Vz())< 5.3) fHistncex_ITScuts->Fill(P/mP);
                               
                               fHistcexVrtx->Fill(fAntipVx,fAntipVy);
                               fHistcexEta_antin->Fill(mom->Eta());
                               
                               //Data
                               Int_t Antip = 0;
                               Int_t proton = 0;
                               AliESDtrack* trackAP = 0;
                               AliESDtrack* trackP = 0;
                               TVector3 ppvec;
                               TVector3 Antippvec;
                               Double_t AntipVx = 0;
                               Double_t pVx = 0;
                               Double_t AntipVy = 0;
                               Double_t pVy = 0;
                               Double_t AntipVz = 0;
                               Double_t pVz = 0;
                               
                               Int_t itracks = fESD->GetNumberOfTracks();
                               
                               //Track loop for antiprotons
                               Bool_t AntipL = 0;
                               for(Int_t j = 0; j < itracks; j++) {
                                   AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(j));
                                   if(!track) continue;
                                   fHistPt->Fill(track->Pt()); //test
                                   
                                   Int_t label = track->GetLabel();
                                   if (label < 0) continue;
                                   
                                   AliMCParticle* particle = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(label));
                                   Int_t mom_label = particle->GetMother();
                                   if(mom_label < 0) continue;
                                   AliMCParticle* mom_particle = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(mom_label));
                                   
                                   if (label == i){
                                       Antip = label;
                                       trackAP = track;
                                       
                                       fAntipP = track->P();
                                       fAntip_Px = track->Px();
                                       fAntip_Py = track->Py();
                                       fAntip_Pz = track->Pz();
                                       fAntipPt = track->Pt();
                                       fAntipE = track->E();
                                       fAntipY = track->Y();
                                       fAntipEta = track->Eta();
                                       fAntipTPCsignal = track->GetTPCsignal();
                                       fAntipTPCpoints = track->GetTPCPoints(4);
                                       fAntipTOFsignal = track->GetTOFsignal();
                                       fAntipITSsignal = track->GetITSsignal();
                                       fAntipncTPC = track->GetTPCNcls();
                                       fAntipncITS = track->GetITSNcls();
                                       Float_t dca_xy = 0;
                                       Float_t dca_z = 0;
                                       track->GetImpactParameters(dca_xy,dca_z);
                                       fAntipDCAxy = dca_xy;
                                       fAntipDCAz = dca_z;
                                       fmom_pdg = mom_particle->PdgCode();
                                       fmom_process = mom_particle->GetUniqueID();
                                       fmom_E= mom_particle->E();
                                       
                                       AntipVx = track->Xv();
                                       AntipVy = track->Yv();
                                       AntipVz = track->Zv();
                                       Antippvec.SetXYZ(fAntip_Px, fAntip_Py, fAntip_Pz);
                                       
                                       //fTree_Antip->Fill();
                                       
                                       fITSmap = track->GetITSClusterMap();
                                       if((fITSmap&kSPDL2) == kSPDL2 || (fITSmap&kSDDL1) == kSDDL1 || (fITSmap&kSDDL2) == kSDDL2 || (fITSmap&kSSDL1) == kSSDL1 || (fITSmap&kSSDL2) == kSSDL2)
                                       {
                                           if((fITSmap&kSPDL1) != kSPDL1 && (fITSmap&kSPDL2) != kSPDL2) AntipL = kTRUE;
                                       }
                                   }
                               }
                               
                               //Track loop for protons
                               Bool_t pL = 0;
                               for(Int_t j = 0; j < itracks; j++) {
                                   AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(j));
                                   if(!track) continue;
                                   
                                   Int_t label = track->GetLabel();
                                   if (label < 0) continue;
                                   
                                   AliMCParticle* particle = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(label));
                                   Int_t mom_label = particle->GetMother();
                                   if(mom_label < 0) continue;
                                   AliMCParticle* mom_particle = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(mom_label));
                                   
                                   if (label == k_plane){
                                       proton = label;
                                       trackP = track;
                                       
                                       fpP = track->P();
                                       fp_Px = track->Px();
                                       fp_Py = track->Py();
                                       fp_Pz = track->Pz();
                                       fpPt = track->Pt();
                                       fpE = track->E();
                                       fpY = track->Y();
                                       fpEta = track->Eta();
                                       fpTPCsignal = track->GetTPCsignal();
                                       fpTPCpoints = track->GetTPCPoints(4);
                                       fpTOFsignal = track->GetTOFsignal();
                                       fpITSsignal = track->GetITSsignal();
                                       fpncTPC = track->GetTPCNcls();
                                       fpncITS = track->GetITSNcls();
                                       Float_t pdca_xy = 0;
                                       Float_t pdca_z = 0;
                                       track->GetImpactParameters(pdca_xy,pdca_z);
                                       fpDCAxy = pdca_xy;
                                       fpDCAz = pdca_z;
                                       fpmom_pdg = mom_particle->PdgCode();
                                       fpmom_process = mom_particle->GetUniqueID();
                                       fpmom_E= mom_particle->E();
                                       
                                       pVx = track->Xv();
                                       pVy = track->Yv();
                                       pVz = track->Zv();
                                       ppvec.SetXYZ(fp_Px, fp_Py, fp_Pz);
                                       
                                       //fTree_proton->Fill();
                                       
                                       fITSmap = track->GetITSClusterMap();
                                       if((fITSmap&kSPDL2) == kSPDL2 || (fITSmap&kSDDL1) == kSDDL1 || (fITSmap&kSDDL2) == kSDDL2 || (fITSmap&kSSDL1) == kSSDL1 || (fITSmap&kSSDL2) == kSSDL2)
                                       {
                                           if((fITSmap&kSPDL1) != kSPDL1 && (fITSmap&kSPDL2) != kSPDL2) pL = kTRUE;
                                       }
                                   }
                               }
                               
                               if (Antip > 0 && proton > 0){
                                   Double_t angle = Antippvec.Angle(ppvec);
                                   Double_t angleD = angle * TMath::RadToDeg();
                                   fHistRecoAngle->Fill(angleD);
                                   //Cuts
                                   if (angleD > 50) continue;
                                   if (TMath::Abs(dplane) > 2) continue;
                                   if (AntipL == kTRUE && pL == kTRUE){
                                       
                                       fTree_Antip->Fill();
                                       fTree_proton->Fill();
                                       
                                       TVector3 pSum = Antippvec + ppvec;
                                       Double_t pSumMag = pSum.Mag();
                                       fHistRecoP->Fill(pSumMag);
                                       
                                       // Distance between individual vertices
                                       Double_t vtxSep = TMath::Sqrt(
                                                                     TMath::Power(pVx - AntipVx, 2) +
                                                                     TMath::Power(pVy - AntipVy, 2) +
                                                                     TMath::Power(pVz - AntipVz, 2)
                                                                     );
                                       fHistRecoVtxSep->Fill(vtxSep);
                                       
                                       //Radius scan to estimate secondary vertex
                                       Double_t b = fESD->GetMagneticField();
                                       AliExternalTrackParam paramP(*trackP);
                                       AliExternalTrackParam paramAP(*trackAP);
                                       
                                       Double_t r1[3], r2[3];
                                       Double_t secVertex[3] = {0};
                                       Double_t minDCA = 9999;
                                       
                                       for (Double_t x = 3.5; x <= 45.0; x += 0.1) {
                                           if (!paramP.GetXYZAt(x, b, r1)) continue;
                                           if (!paramAP.GetXYZAt(x, b, r2)) continue;
                                           
                                           Double_t dx = r1[0] - r2[0];
                                           Double_t dy = r1[1] - r2[1];
                                           Double_t dz = r1[2] - r2[2];
                                           Double_t dca = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
                                           
                                           if (dca < minDCA) {
                                               minDCA = dca;
                                               secVertex[0] = 0.5 * (r1[0] + r2[0]);
                                               secVertex[1] = 0.5 * (r1[1] + r2[1]);
                                               secVertex[2] = 0.5 * (r1[2] + r2[2]);
                                           }
                                       }
                                       
                                       fHistDCApair->Fill(minDCA);
                                       if (minDCA > 10) continue;
                                       
                                       //Radius of secondary vertex
                                       Double_t radius = TMath::Sqrt(secVertex[0]*secVertex[0] + secVertex[1]*secVertex[1]);
                                       fHistRecoradio->Fill(radius);
                                       
                                       //Distance to primary vertex
                                       Double_t primVertex[3];
                                       fESD->GetPrimaryVertex()->GetXYZ(primVertex);
                                       
                                       Double_t dx = secVertex[0] - primVertex[0];
                                       Double_t dy = secVertex[1] - primVertex[1];
                                       Double_t dz = secVertex[2] - primVertex[2];
                                       Double_t distToPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
                                       fHistRecovtx->Fill(distToPrimary);
                                       
                                       //Secondary vertex
                                       fHistRecosvtxXY->Fill(secVertex[0], secVertex[1]);
                                       fHistRecosvtxZ->Fill(secVertex[2]);
                                       
                                       //Tree fill
                                       fPairPSumMag     = pSumMag;
                                       fPairVtxSep      = vtxSep;
                                       fPairMinDCA      = minDCA;
                                       fPairSecRadius   = radius;
                                       fPairVtxDistance = distToPrimary;
                                       fSecVtxX = secVertex[0];
                                       fSecVtxY = secVertex[1];
                                       fSecVtxZ = secVertex[2];
                                       fPrimVtxX = primVertex[0];
                                       fPrimVtxY = primVertex[1];
                                       fPrimVtxZ = primVertex[2];
                                       
                                       fTree_pair->Fill();
                                   }
                               }
                           }
                       }
                   }
               }
           }
           
           PostData(1, fOutputList);
       }
}
//_____________________________________________________________________________
void AliAnalysisTaskAntineutron::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
