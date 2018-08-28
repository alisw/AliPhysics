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

/* AliAnaysisTaskDeuteronAbsorption
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include <iostream>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskDeuteronAbsorption.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliTRDCalDCS.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

class AliAnalysisTaskDeuteronAbsorption;    

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskDeuteronAbsorption) // classimp: necessary for root
  
AliAnalysisTaskDeuteronAbsorption::AliAnalysisTaskDeuteronAbsorption() : AliAnalysisTaskSE(), 
  fESD(0), fPIDResponse(0), fESDtrackCuts(0), fOutputList(0),
  fHistZv(0),
  fHist2PIDvP(0),
  fHist2PIDvka(0),
  fHist2PIDvDe(0),
  fHist3PIDvTr(0),
  fHist2PIDvTr(0),
  fHist2PIDv(0),  
  fHist2PIDf(0),  
  fHist2PIDDef(0),
  fHist3PIDvDe(0),
  fHist3PIDv(0),  
  fHist3PIDf(0),  
  fHist3PIDDef(0),
  fHistPhi(0), 
  fHistPhi2(0),
  fHistPhi2o(0),
  fHistPhi2n(0),
  fHistPhi2no(0),
  fHistmass(0),
  fHistmassDe(0),
  fHistmassP(0), 
  fHistmassDei(0),
  fHistmassDeo(0),
  fHistmassPri(0),
  fHistmassPro(0),
  fHistmassTr(0),
  fHistmassTri(0),
  fHistmassTro(0),
  fHistMatchAllDeuteronPos(0),
  fHistMatchTofDeuteronPos(0), 
  fHistMatchAllDeuteronNeg(0), 
  fHistMatchTofDeuteronNeg(0), 
  fHistMatchAllDeuteronPosMC(0),
  fHistMatchTofDeuteronPosMC(0),
  fHistMatchAllDeuteronNegMC(0),
  fHistMatchTofDeuteronNegMC(0),
  fHistMatchAllDeuteronPoso(0),
  fHistMatchTofDeuteronPoso(0),
  fHistMatchAllDeuteronNego(0),
  fHistMatchTofDeuteronNego(0),
  fHistMatchAllDeuteronPosMCo(0),
  fHistMatchTofDeuteronPosMCo(0),
  fHistMatchAllDeuteronNegMCo(0),
  fHistMatchTofDeuteronNegMCo(0),
  fHistMatchAllDeuteronPosTPCsigma(0),
  fHistMatchAllDeuteronPosTPCsigmao(0),
  fHistMatchAllDeuteronNegTPCsigma(0),
  fHistMatchAllDeuteronNegTPCsigmao(0),
  fHistMatchAllProtonPos(0),
  fHistMatchTofProtonPos(0),
  fHistMatchAllProtonNeg(0),
  fHistMatchTofProtonNeg(0),
  fHistMatchAllProtonPosMC(0),
  fHistMatchTofProtonPosMC(0),
  fHistMatchAllProtonNegMC(0),
  fHistMatchTofProtonNegMC(0),
  fHistMatchAllProtonPoso(0),
  fHistMatchTofProtonPoso(0),
  fHistMatchAllProtonNego(0),
  fHistMatchTofProtonNego(0),
  fHistMatchAllProtonPosMCo(0),
  fHistMatchTofProtonPosMCo(0),
  fHistMatchAllProtonNegMCo(0),
  fHistMatchTofProtonNegMCo(0),
  fHistMatchAllTritonPos(0),
  fHistMatchAllTritonNeg(0),
  fHistMatchTofTritonPos(0),
  fHistMatchTofTritonNeg(0),
  fHistMatchAllTritonPosMC(0),
  fHistMatchAllTritonNegMC(0),
  fHistMatchTofTritonPosMC(0),
  fHistMatchTofTritonNegMC(0),
  hptRecoDeut(0),
  hptRecoAntiDeut(0),
  hptMatchDeut(0),
  hptMatchAntiDeut(0),
  hptGoodMatchDeut(0),
  hptGoodMatchAntiDeut(0),
  hptRecoProt(0),
  hptRecoAntiProt(0),
  hptMatchProt(0),
  hptMatchAntiProt(0),
  hptGoodMatchProt(0),
  hptGoodMatchAntiProt(0)


{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}


//_____________________________________________________________________________
AliAnalysisTaskDeuteronAbsorption::AliAnalysisTaskDeuteronAbsorption(const char* name) : AliAnalysisTaskSE(name),
  fESD(0), fPIDResponse(0), fESDtrackCuts(0), fOutputList(0),
  fHistZv(0),
  fHist2PIDvP(0),
  fHist2PIDvka(0),
  fHist2PIDvDe(0),
  fHist2PIDv(0),  
  fHist2PIDf(0),  
  fHist2PIDDef(0),
  fHist3PIDvDe(0),
  fHist3PIDvTr(0),
  fHist2PIDvTr(0),
  fHist3PIDv(0),  
  fHist3PIDf(0),  
  fHist3PIDDef(0),
  fHistPhi(0), 
  fHistPhi2(0),
  fHistPhi2o(0),
  fHistPhi2n(0),
  fHistPhi2no(0),
  fHistmass(0),
  fHistmassDe(0),
  fHistmassP(0), 
  fHistmassDei(0),
  fHistmassDeo(0),
  fHistmassPri(0),
  fHistmassPro(0),
  fHistmassTr(0),
  fHistmassTri(0),
  fHistmassTro(0),
  fHistMatchAllDeuteronPos(0),
  fHistMatchTofDeuteronPos(0), 
  fHistMatchAllDeuteronNeg(0), 
  fHistMatchTofDeuteronNeg(0), 
  fHistMatchAllDeuteronPosMC(0),
  fHistMatchTofDeuteronPosMC(0),
  fHistMatchAllDeuteronNegMC(0),
  fHistMatchTofDeuteronNegMC(0),
  fHistMatchAllDeuteronPoso(0),
  fHistMatchTofDeuteronPoso(0),
  fHistMatchAllDeuteronNego(0),
  fHistMatchTofDeuteronNego(0),
  fHistMatchAllDeuteronPosMCo(0),
  fHistMatchTofDeuteronPosMCo(0),
  fHistMatchAllDeuteronNegMCo(0),
  fHistMatchTofDeuteronNegMCo(0),
  fHistMatchAllDeuteronPosTPCsigma(0),
  fHistMatchAllDeuteronPosTPCsigmao(0),
  fHistMatchAllDeuteronNegTPCsigma(0),
  fHistMatchAllDeuteronNegTPCsigmao(0),
  fHistMatchAllProtonPos(0),
  fHistMatchTofProtonPos(0),
  fHistMatchAllProtonNeg(0),
  fHistMatchTofProtonNeg(0),
  fHistMatchAllProtonPosMC(0),
  fHistMatchTofProtonPosMC(0),
  fHistMatchAllProtonNegMC(0),
  fHistMatchTofProtonNegMC(0),
  fHistMatchAllProtonPoso(0),
  fHistMatchTofProtonPoso(0),
  fHistMatchAllProtonNego(0),
  fHistMatchTofProtonNego(0),
  fHistMatchAllProtonPosMCo(0),
  fHistMatchTofProtonPosMCo(0),
  fHistMatchAllProtonNegMCo(0),
  fHistMatchTofProtonNegMCo(0),
  fHistMatchAllTritonPos(0),
  fHistMatchAllTritonNeg(0),
  fHistMatchTofTritonPos(0),
  fHistMatchTofTritonNeg(0),
  fHistMatchAllTritonPosMC(0),
  fHistMatchAllTritonNegMC(0),
  fHistMatchTofTritonPosMC(0),
  fHistMatchTofTritonNegMC(0),
  hptRecoDeut(0),
  hptRecoAntiDeut(0),
  hptMatchDeut(0),
  hptMatchAntiDeut(0),
  hptGoodMatchDeut(0),
  hptGoodMatchAntiDeut(0),
  hptRecoProt(0),
  hptRecoAntiProt(0),
  hptMatchProt(0),
  hptMatchAntiProt(0),
  hptGoodMatchProt(0),
  hptGoodMatchAntiProt(0)

{


  //
  // constructor
  //
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                       // this chain is created by the analysis manager, so no need to worry about it, 
                                       // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                       // you can add more output objects by calling DefineOutput(2, classname::Class())
                                       // if you add more output objects, make sure to call PostData for all of them, and to
	
}	



//_____________________________________________________________________________
AliAnalysisTaskDeuteronAbsorption::~AliAnalysisTaskDeuteronAbsorption() {
  //
  // destructor
  //
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
  
}


//_____________________________________________________________________________
void AliAnalysisTaskDeuteronAbsorption::UserCreateOutputObjects() {
  //
  // create output objects
  //
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler) fPIDResponse = inputHandler ->GetPIDResponse();
  }
  //
  // histograms used in the analysis
  // to an output file
  //
  fOutputList = new TList();          // this is a list which will contain all of your histograms                                                                        
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  //
  fHistZv = new TH1F("fHistZv", "fHistZv", 200, -40, 40);       // histogram to monitor z-position of the primary vertex -- quality assurance
  //
  fHist2PIDvP = new TH2F("fHist2PIDvP", "proton dE/dx; P(Gev/c); dE/dx (arb. units)", 1000, -10 ,10 ,1000, 0, 1000);  // proton dE/dx quality assurance
  fHist2PIDvka = new TH2F("fHist2PIDvka", "kaon dE/dx; P(Gev/c); dE/dx (arb. units)", 1000, -10 ,10 ,1000, 0, 1000);  // kaon dE/dx quality assurance
  fHist2PIDvDe = new TH2F("fHist2PIDvDe", "deuteron dE/dx; P(Gev/c); dE/dx (arb. units)", 1000, -10 ,10 ,1000, 0, 1000);  // deuteron dE/dx quality assurance
  fHist2PIDvTr = new TH2F("fHist2PIDvDe", "triton dE/dx; P(Gev/c); dE/dx (arb. units)", 1000, -10 ,10 ,1000, 0, 1000);  // triton dE/dx quality assurance
  fHist2PIDv = new TH2F("fHist2PIDv", "all particles dE/dx; P(Gev/c); dE/dx (arb. units)", 1000, -10 ,10 ,1000, 0, 1000);  // all particles dE/dx quality assurance
  fHist2PIDf = new TH2F("fHist2PIDf", "all paritlces TOF; P(Gev/c); beta", 1000, -10 ,10 ,1000, 0, 1.5);  // all particles TOF quality assurance
  fHist2PIDDef = new TH2F("fHist2PIDDef", "deuteron TOF; P(Gev/c); beta", 1000, -10 ,10 ,1000, 0, 1.5);  // deuteron TOF quality assurance
  //
  fHist3PIDvDe = new TH3F("fHist3PIDvDe", "deuteron dE/dx; P(Gev/c); dE/dx (arb. units); phi", 1000, -10 ,10 ,1000, 0, 500, 20, 0, 2*TMath::Pi());  // dE/dx deuteron vs phi
  fHist3PIDvTr = new TH3F("fHist3PIDvTr", "triton dE/dx; P(Gev/c); dE/dx (arb. units); phi", 1000, -10 ,10 ,1000, 0, 500, 20, 0, 2*TMath::Pi());  // dE/dx triton vs phi
  fHist3PIDv = new TH3F("fHist3PIDv", "all particles dE/dx; P(Gev/c); dE/dx (arb. units); phi", 1000, -10 ,10 ,1000, 0, 500, 20, 0, 2*TMath::Pi());  // dE/dx all particles vs phi
  fHist3PIDf = new TH3F("fHist3PIDf", "all particles TOF; P(Gev/c); beta; phi", 1000, -10 ,10 ,1000, 0, 1.5, 20, 0, 2*TMath::Pi());  // TOF all particles vs phi
  fHist3PIDDef = new TH3F("fHist3PIDDef", "deuteron TOF; P(Gev/c); beta; phi", 1000, -10 ,10 ,1000, 0, 1.5, 20, 0, 2*TMath::Pi());  // TOF deuteron vs phi
  //
  fHistPhi = new TH1F("fHistPhi", "fHistPhi; Phi(rad); dN/d#Phi", 100, 0, 2*TMath::Pi());       // QA histogram for phi
  fHistPhi2 = new TH2F("fHistPhi2", "fHistPhi2; Phi(rad); pT (GeV/c)", 100, 0, 2*TMath::Pi(), 100, 0, 7);       // QA 2d histogra for phi and pt for positive
  fHistPhi2n = new TH2F("fHistPhi2n", "fHistPhi2; Phi(rad); pT (GeV/c)", 100, 0, 2*TMath::Pi(), 100, 0, 7);       // QA 2d histogra for phi and pt for negative

  fHistPhi2o = new TH2F("fHistPhi2o", "fHistPhi2o; Phi(rad); pT (GeV/c)", 100, 0, 2*TMath::Pi(), 100, 0, 7);       // QA 2d histogra for phi and pt for positive
  fHistPhi2no = new TH2F("fHistPhi2no", "fHistPhi2o; Phi(rad); pT (GeV/c)", 100, 0, 2*TMath::Pi(), 100, 0, 7);       // QA 2d histogra for phi and pt for negative

  //
  fHistmass = new TH2F("fHistmass", "mass distr; P(Gev/c); Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       // momentum vs TOF mass
  fHistmassDe = new TH2F("fHistmassDe", "deuteron mass distr; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       // deuterons: momentum vs TOF mass

  fHistmassTr = new TH2F("fHistmassTr", "triton mass distr; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       // deuterons: momentum vs TOF mass

  fHistmassP = new TH2F("fHistmassP", "proton mass distr; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       // protons: momentum vs TOF mass
  //
  fHistmassDei = new TH2F("fHistmassDei", "deuteron mass TRDin; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       //  deuterons: momentum vs TOF mass with TRDin
  fHistmassDeo = new TH2F("fHistmassDeo", "deuteron mass all; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       //  deuterons: momentum vs TOF mass without TRDin
  //
  fHistmassPri = new TH2F("fHistmassPri", "Proton mass TRDin; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       //  deuterons: momentum vs TOF mass with TRDin
  fHistmassPro = new TH2F("fHistmassPro", "Proton mass all; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 6.5);       //  deuterons: momentum vs TOF mass without TRDin

  //
  fHistmassTri = new TH2F("fHistmassTri", "Triton mass TRDin; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 10.5);       //  deuterons: momentum vs TOF mass with TRDin
  fHistmassTro = new TH2F("fHistmassTro", "Triton mass all; P(Gev/c);  Mass^2(GeV^2/c^4)", 1000, -5, 5, 1000,0, 10.5);       //  deuterons: momentum vs TOF mass without TRDin

  // TOF matching histograms
  //
  fHistMatchAllDeuteronPos = new TH2F("fHistMatchAllDeuteronPos", "fHistMatchAllDeuteronPos; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofDeuteronPos = new TH2F("fHistMatchTofDeuteronPos", "fHistMatchTofDeuteronPos; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with TOF matched track
  fHistMatchAllDeuteronNeg = new TH2F("fHistMatchAllDeuteronNeg", "fHistMatchAllDeuteronNeg; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofDeuteronNeg = new TH2F("fHistMatchTofDeuteronNeg", "fHistMatchTofDeuteronNeg; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with TOF matched track
  
  fHistMatchAllDeuteronPoso = new TH2F("fHistMatchAllDeuteronPoso", "fHistMatchAllDeuteronPoso; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofDeuteronPoso = new TH2F("fHistMatchTofDeuteronPoso", "fHistMatchTofDeuteronPoso; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with TOF matched track
  fHistMatchAllDeuteronNego = new TH2F("fHistMatchAllDeuteronNego", "fHistMatchAllDeuteronNego; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofDeuteronNego = new TH2F("fHistMatchTofDeuteronNego", "fHistMatchTofDeuteronNego; momentum p;  mass in TOF", 100, 0, 4.0,  500, 0, 6.5); // filled with TOF matched track


  fHistMatchAllDeuteronPosMC = new TH2F("fHistMatchAllDeuteronPosMC", "fHistMatchAllDeuteronPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks

  fHistMatchTofDeuteronPosMC = new TH2F("fHistMatchTofDeuteronPosMC", "fHistMatchTofDeuteronPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track
  fHistMatchAllDeuteronNegMC = new TH2F("fHistMatchAllDeuteronNegMC", "fHistMatchAllDeuteronNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofDeuteronNegMC = new TH2F("fHistMatchTofDeuteronNegMC", "fHistMatchTofDeuteronNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track
  
  fHistMatchAllDeuteronPosMCo = new TH2F("fHistMatchAllDeuteronPosMCo", "fHistMatchAllDeuteronPosMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchAllDeuteronNegMCo = new TH2F("fHistMatchAllDeuteronNegMCo", "fHistMatchAllDeuteronNegMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofDeuteronPosMCo = new TH2F("fHistMatchTofDeuteronPosMCo", "fHistMatchTofDeuteronPosMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofDeuteronNegMCo = new TH2F("fHistMatchTofDeuteronNegMCo", "fHistMatchTofDeuteronNegMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks

  //
  fHistMatchAllTritonPos = new TH3F("fHistMatchAllTritonPos", "fHistMatchAllTritonPos; momentum p; phi; mass in TOF", 100, 0, 4.0, 100, 0, 2*TMath::Pi(), 500, 0, 10.5); // filled with all TPC tracks
  fHistMatchTofTritonPos = new TH3F("fHistMatchTofTritonPos", "fHistMatchTofTritonPos; momentum p; phi; mass in TOF", 100, 0, 4.0, 100, 0, 2*TMath::Pi(), 500, 0, 10.5); // filled with TOF matched track
  fHistMatchAllTritonNeg = new TH3F("fHistMatchAllTritonNeg", "fHistMatchAllTritonNeg; momentum p; phi; mass in TOF", 100, 0, 4.0, 100, 0, 2*TMath::Pi(), 500, 0, 10.5); // filled with all TPC tracks
  fHistMatchTofTritonNeg = new TH3F("fHistMatchTofTritonNeg", "fHistMatchTofTritonNeg; momentum p; phi; mass in TOF", 100, 0, 4.0, 100, 0, 2*TMath::Pi(), 500, 0, 10.5); // filled with TOF matched track
  fHistMatchAllTritonPosMC = new TH2F("fHistMatchAllTritonPosMC", "fHistMatchAllTritonPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofTritonPosMC = new TH2F("fHistMatchTofTritonPosMC", "fHistMatchTofTritonPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track
  fHistMatchAllTritonNegMC = new TH2F("fHistMatchAllTritonNegMC", "fHistMatchAllTritonNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofTritonNegMC = new TH2F("fHistMatchTofTritonNegMC", "fHistMatchTofTritonNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track

   
  fHistMatchAllProtonPosMC = new TH2F("fHistMatchAllProtonPosMC", "fHistMatchAllProtonPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofProtonPosMC = new TH2F("fHistMatchTofProtonPosMC", "fHistMatchTofProtonPosMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track
  fHistMatchAllProtonNegMC = new TH2F("fHistMatchAllProtonNegMC", "fHistMatchAllProtonNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofProtonNegMC = new TH2F("fHistMatchTofProtonNegMC", "fHistMatchTofProtonNegMC; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track

  fHistMatchAllProtonPosMCo = new TH2F("fHistMatchAllProtonPosMCo", "fHistMatchAllProtonPosMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofProtonPosMCo = new TH2F("fHistMatchTofProtonPosMCo", "fHistMatchTofProtonPosMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track
  fHistMatchAllProtonNegMCo = new TH2F("fHistMatchAllProtonNegMCo", "fHistMatchAllProtonNegMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with all TPC tracks
  fHistMatchTofProtonNegMCo = new TH2F("fHistMatchTofProtonNegMCo", "fHistMatchTofProtonNegMCo; momentum p; phi", 100, 0, 4.0, 100, 0, 2*TMath::Pi()); // filled with TOF matched track


  fHistMatchAllProtonPos = new TH2F("fHistMatchAllProtonPos", "fHistMatchAllProtonPos; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofProtonPos = new TH2F("fHistMatchTofProtonPos", "fHistMatchTofProtonPos; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with TOF matched track
  fHistMatchAllProtonNeg = new TH2F("fHistMatchAllProtonNeg", "fHistMatchAllProtonNeg; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofProtonNeg = new TH2F("fHistMatchTofProtonNeg", "fHistMatchTofProtonNeg; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with TOF matched track

  fHistMatchAllProtonPoso = new TH2F("fHistMatchAllProtonPoso", "fHistMatchAllProtonPoso; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofProtonPoso = new TH2F("fHistMatchTofProtonPoso", "fHistMatchTofProtonPoso; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with TOF matched track
  fHistMatchAllProtonNego = new TH2F("fHistMatchAllProtonNego", "fHistMatchAllProtonNego; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with all TPC tracks
  fHistMatchTofProtonNego = new TH2F("fHistMatchTofProtonNego", "fHistMatchTofProtonNego; momentum p; mass in TOF", 100, 0, 4.0, 500, 0, 6.5); // filled with TOF matched track

  fHistMatchAllDeuteronPosTPCsigma =  new TH2F("fHistMatchAllDeuteronPosTPCsigma", "fHistMatchAllDeuteronPosTPCsigma; momentum p;  n#sigma^{d}_{TPC}", 100, 0, 4.0,  100, -5, 5); // filled with all TPC tracks
  fHistMatchAllDeuteronPosTPCsigmao =  new TH2F("fHistMatchAllDeuteronPosTPCsigmao", "fHistMatchAllDeuteronPosTPCsigmao; momentum p;  n#sigma^{d}_{TPC}", 100, 0, 4.0,  100, -5, 5); // filled with all TPC tracks
  fHistMatchAllDeuteronNegTPCsigma =  new TH2F("fHistMatchAllDeuteronNegTPCsigma", "fHistMatchAllDeuteronNegTPCsigma; momentum p;  n#sigma^{d}_{TPC}", 100, 0, 4.0,  100, -5, 5); // filled with all TPC tracks
  fHistMatchAllDeuteronNegTPCsigmao =  new TH2F("fHistMatchAllDeuteronNegTPCsigmao", "fHistMatchAllDeuteronNegTPCsigmao; momentum p;  n#sigma^{d}_{TPC}", 100, 0, 4.0,  100, -5, 5); // filled with all TPC tracks

  hptRecoDeut = new TH1F("hptMCRecoDeut","MC reco. deuterons;pt(GeV/c)",500,0,10);
  hptRecoAntiDeut  = new TH1F("hptMCRecoAntiDeut","MC reco. antideuterons;pt(GeV/c)",500,0,10);

  hptMatchDeut  = new TH1F("hptMCMatchDeut","MC tof matched deuterons;pt(GeV/c)",500,0,10);
  hptMatchAntiDeut  = new TH1F("hptMCMatchAntiDeut","MC tof matched antideuterons;pt(GeV/c)",500,0,10);

  hptGoodMatchDeut  = new TH1F("hptMCGoodMatchDeut","MC good tof matched deuterons;pt(GeV/c)",500,0,10);
  hptGoodMatchAntiDeut  = new TH1F("hptMCGoodMatchAntiDeut","MC good tof matched antideuterons;pt(GeV/c)",500,0,10);

  hptRecoProt = new TH1F("hptMCRecoProt","MC reco. protons;pt(GeV/c)",500,0,10);
  hptRecoAntiProt  = new TH1F("hptMCRecoAntiProt","MC reco. antiprotons;pt(GeV/c)",500,0,10);

  hptMatchProt  = new TH1F("hptMCMatchProt","MC tof matched protons;pt(GeV/c)",500,0,10);
  hptMatchAntiProt  = new TH1F("hptMCMatchAntiProt","MC tof matched antiprotons;pt(GeV/c)",500,0,10);


  hptGoodMatchProt  = new TH1F("hptMCGoodMatchProt","MC good tof matched Protons;pt(GeV/c)",500,0,10);
  hptGoodMatchAntiProt  = new TH1F("hptMCGoodMatchAntiProt","MC good tof matched antiprotons;pt(GeV/c)",500,0,10);


//
  // add all histograms to output list
  //
  fOutputList->Add(fHistZv);
  //
  fOutputList->Add(fHist2PIDvP);
  fOutputList->Add(fHist2PIDvka);  
  fOutputList->Add(fHist2PIDvDe);  
  fOutputList->Add(fHist2PIDvTr);
  fOutputList->Add(fHist2PIDv);  
  fOutputList->Add(fHist2PIDf);  
  fOutputList->Add(fHist2PIDDef);
  //
  fOutputList->Add(fHist3PIDvDe);
  fOutputList->Add(fHist3PIDvTr);
  fOutputList->Add(fHist3PIDv);  
  fOutputList->Add(fHist3PIDf);  
  fOutputList->Add(fHist3PIDDef);
  //
  fOutputList->Add(fHistPhi);  
  fOutputList->Add(fHistPhi2); 
  fOutputList->Add(fHistPhi2n);
  fOutputList->Add(fHistPhi2o);
  fOutputList->Add(fHistPhi2no);


  //
  fOutputList->Add(fHistmass); 
  fOutputList->Add(fHistmassDe); 
  fOutputList->Add(fHistmassTr);

  fOutputList->Add(fHistmassP);  
  //
  fOutputList->Add(fHistmassDei);
  fOutputList->Add(fHistmassDeo);
  //
  fOutputList->Add(fHistmassPri);
  fOutputList->Add(fHistmassPro);
  //
  //
  fOutputList->Add(fHistmassTri);
  fOutputList->Add(fHistmassTro);
  //
 
  fOutputList->Add(fHistMatchAllDeuteronPos);
  fOutputList->Add(fHistMatchTofDeuteronPos); 
  fOutputList->Add(fHistMatchAllDeuteronNeg); 
  fOutputList->Add(fHistMatchTofDeuteronNeg); 
  //
  fOutputList->Add(fHistMatchAllDeuteronPoso);
  fOutputList->Add(fHistMatchTofDeuteronPoso);
  fOutputList->Add(fHistMatchAllDeuteronNego);
  fOutputList->Add(fHistMatchTofDeuteronNego);

  fOutputList->Add(fHistMatchAllDeuteronPosMC);
  fOutputList->Add(fHistMatchTofDeuteronPosMC);
  fOutputList->Add(fHistMatchAllDeuteronNegMC);
  fOutputList->Add(fHistMatchTofDeuteronNegMC);
  
  fOutputList->Add(fHistMatchAllDeuteronPosMCo);
  fOutputList->Add(fHistMatchAllDeuteronNegMCo);
  fOutputList->Add(fHistMatchTofDeuteronPosMCo);
  fOutputList->Add(fHistMatchTofDeuteronNegMCo);

  
  fOutputList->Add(fHistMatchAllTritonPos);
  fOutputList->Add(fHistMatchTofTritonPos);
  fOutputList->Add(fHistMatchAllTritonNeg);
  fOutputList->Add(fHistMatchTofTritonNeg);
  //
  fOutputList->Add(fHistMatchAllTritonPosMC);
  fOutputList->Add(fHistMatchTofTritonPosMC);
  fOutputList->Add(fHistMatchAllTritonNegMC);
  fOutputList->Add(fHistMatchTofTritonNegMC);

  //
  fOutputList->Add(fHistMatchAllProtonPos);
  fOutputList->Add(fHistMatchTofProtonPos);
  fOutputList->Add(fHistMatchAllProtonNeg);
  fOutputList->Add(fHistMatchTofProtonNeg);

  fOutputList->Add(fHistMatchAllProtonPoso);
  fOutputList->Add(fHistMatchTofProtonPoso);
  fOutputList->Add(fHistMatchAllProtonNego);
  fOutputList->Add(fHistMatchTofProtonNego);

  //
  
  fOutputList->Add(fHistMatchAllProtonPosMC);
  fOutputList->Add(fHistMatchTofProtonPosMC);
  fOutputList->Add(fHistMatchAllProtonNegMC);
  fOutputList->Add(fHistMatchTofProtonNegMC);
  //
  fOutputList->Add(fHistMatchAllProtonPosMCo);
  fOutputList->Add(fHistMatchTofProtonPosMCo);
  fOutputList->Add(fHistMatchAllProtonNegMCo);
  fOutputList->Add(fHistMatchTofProtonNegMCo);
  // 
  
  fOutputList->Add(hptRecoDeut); 
  fOutputList->Add(hptRecoAntiDeut); 
  fOutputList->Add(hptMatchDeut);  
  fOutputList->Add(hptMatchAntiDeut); 
  fOutputList->Add(hptGoodMatchDeut); 
  fOutputList->Add(hptGoodMatchAntiDeut); 

   //
  fOutputList->Add(hptRecoProt);
  fOutputList->Add(hptRecoAntiProt);
  fOutputList->Add(hptMatchProt);
  fOutputList->Add(hptMatchAntiProt);
  fOutputList->Add(hptGoodMatchProt);
  fOutputList->Add(hptGoodMatchAntiProt);

  fOutputList->Add(fHistMatchAllDeuteronPosTPCsigma);
  fOutputList->Add(fHistMatchAllDeuteronPosTPCsigmao);
  fOutputList->Add(fHistMatchAllDeuteronNegTPCsigma);
  fOutputList->Add(fHistMatchAllDeuteronNegTPCsigmao);


//
  // create track cuts object
  //
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,kTRUE);
  // fESDtrackCuts->SetMaxDCAToVertexXY(3);
  // fESDtrackCuts->SetMaxDCAToVertexZ(2);
   fESDtrackCuts->SetEtaRange(-0.8,0.8);    
    
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
  
  
  double paramsPos[4][4]={
    {1.38984e+00,-2.10187e+01,5.81724e-02,1.91938e+01},
    {2.02372e+00,-2.44456e+00,8.99000e-01,9.22399e-01},
    {4.21954e+00,-2.56555e+01,4.17557e-02,2.40301e+01},
    {5.17499e+00,-2.69241e+00,6.97167e-01,1.25974e+00}
  };
  for (int iFunction = 0; iFunction < 4; ++iFunction) {  
    fTRDboundariesPos[iFunction] = new TF1(Form("f%i",iFunction),"[0]-exp([1]*pow(x,[2])+[3])",0.2,10);
    for (int iParam = 0; iParam < 4; ++iParam) {
      fTRDboundariesPos[iFunction]->SetParameter(iParam, paramsPos[iFunction][iParam]);
    }
  }

  double paramsNeg[4][4]={
    {2.81984e+00,-1.81497e-01, -2.03494e+00,2.64148e-01},
    {5.79322e+00,-5.44966e-02,-1.10803e+00,1.29737e+00},
    {5.60000e+00,-2.06000e-01,-1.97130e+00,2.67181e-01},
    {9.72180e+00,-4.35801e-02,-1.14550e+00,1.49160e+00}
  };
  for (int iFunction = 0; iFunction < 4; ++iFunction) {
    fTRDboundariesNeg[iFunction] = new TF1(Form("f%i",iFunction),"[0]-exp([1]*pow(x,[2])+[3])",0.2,10);
    for (int iParam = 0; iParam < 4; ++iParam) {
      fTRDboundariesNeg[iFunction]->SetParameter(iParam, paramsNeg[iFunction][iParam]);
    }
  }

  
  }


//_____________________________________________________________________________
 void AliAnalysisTaskDeuteronAbsorption::UserExec(Option_t *) {
  //
  // main loop over events
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());                                                        
  if(!fESD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
  Int_t iTracks = fESD->GetNumberOfTracks();           // see how many tracks there are in the event
  
  Bool_t isMC = kTRUE;
  AliMCEventHandler* eventHandlerMC = 0x0;
  AliMCEvent* mcEvent = 0x0;                    
  AliStack* fStack = 0x0;

  eventHandlerMC = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!eventHandlerMC) {
    isMC = kFALSE;
  } else {
    mcEvent = eventHandlerMC->MCEvent();
    if (!mcEvent) {
      isMC = kFALSE; 
    } else {
      fStack = mcEvent->Stack();
    }
  }
  //
  // check for a proper primary vertex and monitor
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) vertex = 0x0;
  }  
  if (!vertex) return;
  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0) return; // remove events with a vertex which is more than 10cm away



  //
  // track loop
  //
  for(Int_t i = 0; i < iTracks; i++) { // loop ove rall these tracks
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));         // get a track (type AliESDDTrack) from the event
    if(!track) continue;
    //
    if (!fESDtrackCuts->AcceptTrack(track)) continue; // check if track passes the cuts
    if (!track->GetInnerParam()) continue; // check if track is a proper TPC track
    Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
    if (track->GetTPCsignalN() < 50) continue;
 
    //
    // check phi distriubtion after cuts
    //
    fHistPhi->Fill(track->Phi());
    //
    // get TPC nsigma PID
    //
    Double_t sign = track->GetSign();
    Double_t kaSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Double_t prSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Double_t DeSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);  
    Double_t TrSignal = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);

    //Float_t deutExp = AliExternalTrackParam::BetheBlochAleph(ptot/(0.938*2),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
    //
    // fill dE/dx QA histograms
    //
    fHist2PIDv->Fill(ptot*sign, track->GetTPCsignal());  // all particles
    fHist3PIDv->Fill(ptot*sign, track->GetTPCsignal(), track->Phi()); // all particles
    //
    // after nsigma TPC cut
    //
    if(TMath::Abs(prSignal)< 3 ) fHist2PIDvP->Fill(ptot*sign, track->GetTPCsignal());	       
    if(TMath::Abs(kaSignal)< 3 ) fHist2PIDvka->Fill(ptot*sign, track->GetTPCsignal());
    if(TMath::Abs(DeSignal)< 3 ) {    
      fHist2PIDvDe->Fill(ptot*sign, track->GetTPCsignal());
      fHist3PIDvDe->Fill(ptot*sign, track->GetTPCsignal(), track->Phi());  
    }

     if(TMath::Abs(TrSignal)< 3 ) {
      fHist2PIDvTr->Fill(ptot*sign, track->GetTPCsignal());
      fHist3PIDvTr->Fill(ptot*sign, track->GetTPCsignal(), track->Phi());
    }

    //
    // Process TOF information
    //      
    ULong_t status = (ULong_t)track->GetStatus();
    Bool_t hasTOF     = kFALSE;
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout;
    if (hasTOFout) hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength();
    if (length < 350.) hasTOF = kFALSE;
    //
    Double_t p    = track->P();
    Double_t tof  = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(p);
    //
    Float_t  beta = 0;
    Float_t  gamma = 0;
    Float_t  mass  = -99;
    //
    if (hasTOF) {
      beta = length / (2.99792457999999984e-02 * tof);
      //
      fHist2PIDf->Fill(track->P()*sign, beta); // QA histograms tof beta all particles
      fHist3PIDf->Fill(track->P()*sign, beta, track->Phi()); // QA histograms
      //
      if((1 - beta*beta) > 0){
	gamma = 1/TMath::Sqrt(1 - beta*beta);
	mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.}
      }
      else{
	gamma = 0;
	mass = 0;
      }
    }
    fHistmass->Fill(track->P()*sign, mass*mass); // QA histogram for mass calculation
    //
    // fill histograms for matching study (deuterons)
    //
    if(TMath::Abs(DeSignal)< 3 && track->GetTPCsignal() > 130) {    
      //
      fHistmassDe->Fill(ptot*sign, mass*mass);
      fHist2PIDDef->Fill(track->P()*sign, beta);
      fHist3PIDDef->Fill(track->P()*sign, beta, track->Phi());
      //
      // fill the TOF matching histograms
    }
    
     //triton

    if(TMath::Abs(TrSignal)< 3 ) {
      //
      fHistmassTr->Fill(ptot*sign, mass*mass);
      //
      if (sign > 0) fHistMatchAllTritonPos->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi(), mass*mass);
      if (sign > 0 && hasTOF) fHistMatchTofTritonPos->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi(), mass*mass);
      if (sign < 0) fHistMatchAllTritonNeg->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi(), mass*mass);
      if (sign < 0 && hasTOF) fHistMatchTofTritonNeg->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi(), mass*mass);
    }


    
    // study using TRDin
    //
    ULong_t hasTRDin  = status&AliESDtrack::kTRDin; // 2D phi pt for TRD

    if (sign > 0 ) {
	// TRD in
      float pt = track->Pt();
      float phi = track->Phi();
      while (phi < 0) phi += TMath::TwoPi();
      while (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
      if (phi < fTRDboundariesPos[0]->Eval(pt) ||
	  (phi > fTRDboundariesPos[1]->Eval(pt) && phi < fTRDboundariesPos[2]->Eval(pt)) ||
          phi > fTRDboundariesPos[3]->Eval(pt)){
        fHistPhi2->Fill(phi,pt);

       fHistMatchAllDeuteronPosTPCsigma->Fill(track->GetInnerParam()->GetP(), DeSignal);        
       if(TMath::Abs(DeSignal)< 3  && track->GetTPCsignal() > 100 ) {

        fHistMatchAllDeuteronPos->Fill(track->GetInnerParam()->GetP(), mass*mass);
        }
        if(TMath::Abs(prSignal)< 3 ) {

        fHistMatchAllProtonPos->Fill(track->GetInnerParam()->GetP(), mass*mass);
        }
 
        }
        else {

       fHistMatchAllDeuteronPosTPCsigmao->Fill(track->GetInnerParam()->GetP(), DeSignal);        
       fHistPhi2o->Fill(phi,pt);

        if(TMath::Abs(DeSignal)< 3 && track->GetTPCsignal() > 100) {

       fHistMatchAllDeuteronPoso->Fill(track->GetInnerParam()->GetP(), mass*mass);
     } 

      if(TMath::Abs(prSignal)< 3 ) {

       fHistMatchAllProtonPoso->Fill(track->GetInnerParam()->GetP(), mass*mass);
         }

         }

        }

        if (sign < 0 ) {
        // TRD in
      float pt = track->Pt();
      float phi = track->Phi();
      while (phi < 0) phi += TMath::TwoPi();
      while (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
      if (phi < fTRDboundariesNeg[0]->Eval(pt) ||
          (phi > fTRDboundariesNeg[1]->Eval(pt) && phi < fTRDboundariesNeg[2]->Eval(pt)) ||
          phi > fTRDboundariesNeg[3]->Eval(pt)){
        fHistPhi2n->Fill(phi,pt);
        fHistMatchAllDeuteronNegTPCsigma->Fill(track->GetInnerParam()->GetP(), DeSignal);

         if(TMath::Abs(DeSignal)< 3 && track->GetTPCsignal() > 100 ) {

        fHistMatchAllDeuteronNeg->Fill(track->GetInnerParam()->GetP(),  mass*mass);
      }
        if(TMath::Abs(prSignal)< 3 ) {

        fHistMatchAllProtonNeg->Fill(track->GetInnerParam()->GetP(), mass*mass);
        }

       // }
        }        
         else {

       fHistPhi2no->Fill(phi,pt);
       fHistMatchAllDeuteronNegTPCsigmao->Fill(track->GetInnerParam()->GetP(), DeSignal);

         if(TMath::Abs(DeSignal)< 3 && track->GetTPCsignal() > 100) {

        fHistMatchAllDeuteronNego->Fill(track->GetInnerParam()->GetP(),  mass*mass);
       }
         if(TMath::Abs(prSignal)< 3 ) {

        fHistMatchAllProtonNego->Fill(track->GetInnerParam()->GetP(), mass*mass);
        }

        // }
         }


      }
    
    //
    //---------- info on MC true:
    //
    if (!isMC) continue;
    //
    AliVParticle *mcpart;
    Int_t label=-999, Pdg=-1, t_label=-1;
    Bool_t isPrimary=kFALSE;
    Double_t t_pt=-1;

    label = TMath::Abs(track->GetLabel());

    mcpart = (AliVParticle *) mcEvent->GetTrack(label);

    Pdg = mcpart->PdgCode();//e.g. deuteron is 1000010020, antideuteron is -1000010020
    t_label = mcpart->GetLabel();
    isPrimary = fStack->IsPhysicalPrimary(t_label);
    //Double_t t_rapidity = mcpart->Y();
    t_pt = mcpart->Pt();
    //---------- info on MC true (end)

    if(!isPrimary) continue;

    if(Pdg==1000010020)
      hptRecoDeut->Fill(t_pt);//t_pt distribution of *reconstruced* particles     
    else if(Pdg==-1000010020)
      hptRecoAntiDeut->Fill(t_pt);//t_pt distribution of *reconstruced* anti particles
    
    if(Pdg==1000010020) {
      if (hasTRDin) {
	fHistMatchAllDeuteronPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      } else {
	fHistMatchAllDeuteronPosMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      }
    }

    if(Pdg==-1000010020) {
      if (hasTRDin) {
	fHistMatchAllDeuteronNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      } else {
	fHistMatchAllDeuteronNegMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      }
    }

         //triton  
   
    if(Pdg==1000010030)
      fHistMatchAllTritonPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
    else if(Pdg==-1000010030)
      fHistMatchAllTritonNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      

    if(Pdg==2212)

      hptRecoProt->Fill(t_pt);//t_pt distribution of *reconstruced* particles
    else if(Pdg==-2212)
      hptRecoAntiProt->Fill(t_pt);//t_pt distribution of *reconstruced* anti particles
      
      if(Pdg==2212){
       if (hasTRDin) {
      fHistMatchAllProtonPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
       } else {
       fHistMatchAllProtonPosMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
        
       
         
     }
    }

     if(Pdg==-2212)  {
       if (hasTRDin) {

      fHistMatchAllProtonNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
       } else {

       fHistMatchAllProtonNegMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
         
       }
      }
          
    if(!hasTOF) continue;

    //
    //For checking the TOF GOOD matching:
    Int_t *toflabel = new Int_t[3];
    ((AliESDtrack *)track)->GetTOFLabel(toflabel);
    
    //2) using nsigmaTOF:
    Double_t exptimes[9];//e, mu, pi, K, p, d, t, 3He, 4He
    track->GetIntegratedTimes(exptimes);
    Double_t m_proton = AliPID::ParticleMass(4);

    for(Int_t iN=5;iN<9;iN++) {
      Double_t massOverZ = AliPID::ParticleMassZ(iN);
      if(p>1e-18) exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ*massOverZ/p/p+1)/(m_proton*m_proton/p/p+1);
      exptimes[iN] = TMath::Sqrt(exptimes[iN]);
    }
    
    //getting TOF res:
    Double_t tofres[9];
    for(Int_t i=0;i<9;i++) tofres[i] = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, exptimes[i], (AliPID::EParticleType) i);

    //finally you compute your NsigmaTOF:
    Double_t nsigmaTOF[9];
    for(Int_t i=0;i<9;i++) {
      nsigmaTOF[i] = -99999.9;
      if(tofres[i]>1e-18) nsigmaTOF[i] = (tof-exptimes[i])/tofres[i];
    }
    
    
    if(Pdg==1000010020)
      hptMatchDeut->Fill(t_pt);//t_pt distribution of *TOF matched* particles
    
    else if(Pdg==-1000010020)
      hptMatchAntiDeut->Fill(t_pt);
    

    if(Pdg==1000010020) {
      if (hasTRDin) {
	if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofDeuteronPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      } else{
        if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofDeuteronPosMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      }
    } 
     
    if(Pdg==-1000010020) {
      if (hasTRDin) {
	if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofDeuteronNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      } else {
	if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofDeuteronNegMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      }
    } 


    if(Pdg==1000010030) {
      if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofTritonPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
    } else if(Pdg==-1000010030){
      if (toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<5) fHistMatchTofTritonNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
    }

    
    if(Pdg==2212)
      
      hptMatchProt->Fill(t_pt);//t_pt distribution of *TOF matched* particles
    else if(Pdg==-2212)
      hptMatchAntiProt->Fill(t_pt);
    
    if(Pdg==2212){
     if (hasTRDin) {
      if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<5) fHistMatchTofProtonPosMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      } else {
     if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<5) fHistMatchTofProtonPosMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      }
     }  

     if(Pdg==-2212) {
     if (hasTRDin) {

      if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<5) fHistMatchTofProtonNegMC->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
      
      } else {
      
      if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<5) fHistMatchTofProtonNegMCo->Fill(track->GetInnerParam()->GetP(), track->GetInnerParam()->Phi());
     }
      
    }

    //1) || 2)
    if(toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<3) {
      if(Pdg==1000010020) hptGoodMatchDeut->Fill(t_pt);//t_pt distribution of *TOF GOOD matched* particles
    }
    if(toflabel[0]==label || TMath::Abs(nsigmaTOF[5])<3) { 
      if(Pdg==-1000010020) hptGoodMatchAntiDeut->Fill(t_pt);
    }
    

    //1) || 2) Proton
     
    if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<3) {
      if(Pdg==2212) hptGoodMatchProt->Fill(t_pt);//t_pt distribution of *TOF GOOD matched* particles
    }
    if(toflabel[0]==label || TMath::Abs(nsigmaTOF[4])<3) {
      
      if(Pdg==-2212) hptGoodMatchAntiProt->Fill(t_pt);
    }
    
    
  } // end the track loop

  //
  // post the data
  //
  PostData(1, fOutputList);          
} // end the UserExec


//_____________________________________________________________________________
void AliAnalysisTaskDeuteronAbsorption::Terminate(Option_t *) {



} 


