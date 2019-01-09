/*An analysis task that plots the TPC dE/dx spectra Protons, Kaons, Pions and Antiprotons*/
/*A. Kalweit and S. Marium*/
/*26.01.2018*/

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TLatex.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskTpcSpectra.h"
#include "AliMultSelection.h"

class AliAnalysisTaskTpcSpectra;

using namespace std;

ClassImp(AliAnalysisTaskTpcSpectra)

//____________________________________________________________________________________//

AliAnalysisTaskTpcSpectra::AliAnalysisTaskTpcSpectra() : AliAnalysisTaskTpcSpectra("stdtask")     
{
  //
  // default contstructor: do nothing
  //
}
  
//____________________________________________________________________________________//

AliAnalysisTaskTpcSpectra::AliAnalysisTaskTpcSpectra(const char* name): AliAnalysisTaskSE(name), fESD(0), fOutputList(0),
  fEventCut(0),
  fESDtrackCuts(0),
  fESDtrackCutsWoDCA(0),
  fPIDResponse(0),
  fMultSel(0),
  fHistZv(0),
  fHistdEdxData(0),
  fHistTof(0),
  fHistCentBefEvSel(0),
  fHistCent(0),
  fHistMuonGen(0),
  fHistMuonReco(0),
  fHistMuonRecoTOF(0),
  fHistPionGen(0),
  fHistPionReco(0),
  fHistPionRecoTOF(0),
  fHistKaonGen(0),
  fHistKaonReco(0),
  fHistKaonRecoTOF(0),
  fHistProtGen(0),
  fHistProtReco(0),
  fHistProtRecoTOF(0),
  fHistAntiMuonGen(0),
  fHistAntiMuonReco(0),
  fHistAntiMuonRecoTOF(0),  
  fHistAntiPionGen(0),
  fHistAntiPionReco(0),
  fHistAntiPionRecoTOF(0),  
  fHistAntiKaonGen(0),
  fHistAntiKaonReco(0),
  fHistAntiKaonRecoTOF(0),  
  fHistAntiProtGen(0),
  fHistAntiProtReco(0),
  fHistAntiProtRecoTOF(0),  
  fMuonDeDxCent(0),
  fPionDeDxCent(0),
  fKaonDeDxCent(0),
  fProtonDeDxCent(0),
  fAntiMuonDeDxCent(0),
  fAntiPionDeDxCent(0),
  fAntiKaonDeDxCent(0),
  fAntiProtonDeDxCent(0),
  fDCAPion(0),                            
  fDCAAntiPion(0),                        
  fDCAKaon(0),                            
  fDCAAntiKaon(0),                        
  fDCAProton(0),                          
  fDCAAntiProton(0),
  fDCAPionMC(0),                            
  fDCAAntiPionMC(0),                        
  fDCAKaonMC(0),                            
  fDCAAntiKaonMC(0),                        
  fDCAProtonMC(0),                          
  fDCAAntiProtonMC(0)      
{
  //
  // main constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________________________//

AliAnalysisTaskTpcSpectra::~AliAnalysisTaskTpcSpectra()
{
  //
  // destructor
  //
  if(fOutputList) {
    delete fOutputList;
  }
  
}

//____________________________________________________________________________________//

void AliAnalysisTaskTpcSpectra::UserCreateOutputObjects()
{
  //
  // create the output objects
  //
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fOutputList);
  //
  // pt-bins
  //
  const Int_t nPtBins = 62;
  const Double_t binsPt[nPtBins + 1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
					0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1,
					1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
					2.7, 2.8, 2.9, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,
					10, 11, 12, 13, 14, 15, 16, 18, 20};
  //
  // create QA histograms
  //
  fHistZv = new TH1F("fHistZv", "fHistZv; vertex z (cm)", 200, -40, 40);  
  fHistdEdxData = new TH2F("fHistdEdxData", "fHistdEdxData; p/z (GeV/c); TPC signal", 500, -10, 10, 250, 0, 500);
  fHistTof = new TH2F("fHistTof", "all paritlces TOF; P(Gev/c); beta", 1000, -10 ,10 ,1000, 0, 1.5);  // all particles TOF quality assurance
  fHistCentBefEvSel = new TH1F("fHistCentBefEvSel", "fHistCentBefEvSel", 110, -10.0, 100.0);
  fHistCent = new TH1F("fHistCent", "fHistCent", 110, -10.0, 100.0);
  //
  // create histograms for efficiency calculation
  //
  fHistMuonGen = new TH1F("fHistMuonGen", "Generated Muon #phi vs p_{t}", nPtBins, binsPt);
  fHistMuonReco = new TH1F("fHistMuonReco", "Reconstructed Muon #phi vs p_{t}", nPtBins, binsPt);
  fHistMuonRecoTOF = new TH1F("fHistMuonRecoTOF", "Reconstructed Muon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistPionGen = new TH1F("fHistPionGen", "Generated Pion #phi vs p_{t}", nPtBins, binsPt);
  fHistPionReco = new TH1F("fHistPionReco", "Reconstructed Pion #phi vs p_{t}", nPtBins, binsPt);
  fHistPionRecoTOF = new TH1F("fHistPionRecoTOF", "Reconstructed Pion #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistKaonGen = new TH1F("fHistKaonGen", "Generated Kaon #phi vs p_{t}", nPtBins, binsPt);
  fHistKaonReco = new TH1F("fHistKaonReco", "Reconstructed Kaon #phi vs p_{t}", nPtBins, binsPt);
  fHistKaonRecoTOF = new TH1F("fHistKaonRecoTOF", "Reconstructed Kaon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistProtGen = new TH1F("fHistProtGen", "Generated Proton #phi vs p_{t}", nPtBins, binsPt);
  fHistProtReco = new TH1F("fHistProtReco", "Reconstructed Proton #phi vs p_{t}", nPtBins, binsPt);
  fHistProtRecoTOF = new TH1F("fHistProtRecoTOF", "Reconstructed Proton #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiMuonGen = new TH1F("fHistAntiMuonGen", "Generated AntiMuon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiMuonReco = new TH1F("fHistAntiMuonReco", "Reconstructed AntiMuon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiMuonRecoTOF = new TH1F("fHistAntiMuonRecoTOF", "Reconstructed AntiMuon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiPionGen = new TH1F("fHistAntiPionGen", "Generated AntiPion #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiPionReco = new TH1F("fHistAntiPionReco", "Reconstructed AntiPion #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiPionRecoTOF = new TH1F("fHistAntiPionRecoTOF", "Reconstructed AntiPion #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiKaonGen = new TH1F("fHistAntiKaonGen", "Generated AntiKaon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiKaonReco = new TH1F("fHistAntiKaonReco", "Reconstructed AntiKaon #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiKaonRecoTOF = new TH1F("fHistAntiKaonRecoTOF", "Reconstructed AntiKaon #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  fHistAntiProtGen = new TH1F("fHistAntiProtGen", "Generated AntiProton #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiProtReco = new TH1F("fHistAntiProtReco", "Reconstructed AntiProton #phi vs p_{t}", nPtBins, binsPt);
  fHistAntiProtRecoTOF = new TH1F("fHistAntiProtRecoTOF", "Reconstructed AntiProton #phi vs p_{t} with TOF", nPtBins, binsPt);
  //
  // create raw spectra histograms
  //
  fMuonDeDxCent   = new TH3F("fMuonDeDxCent", "fMuonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10,  200, -50, 150);
  fPionDeDxCent   = new TH3F("fPionDeDxCent", "fPionDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10,  200, -50, 150);
  fKaonDeDxCent   = new TH3F("fKaonDeDxCent", "fKaonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20,  300, -10, 10, 200, -50, 150);
  fProtonDeDxCent = new TH3F("fProtonDeDxCent", "fProtonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -50, 150);  
  fMuonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fPionDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fKaonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fProtonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  //
  fAntiMuonDeDxCent   = new TH3F("fAntiMuonDeDxCent", "fAntiMuonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10,  200, -50, 150);
  fAntiPionDeDxCent   = new TH3F("fAntiPionDeDxCent", "fAntiPionDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10,  200, -50, 150);
  fAntiKaonDeDxCent   = new TH3F("fAntiKaonDeDxCent", "fAntiKaonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20,  300, -10, 10, 200, -50, 150);
  fAntiProtonDeDxCent = new TH3F("fAntiProtonDeDxCent", "fAntiProtonDeDxCent; pT; nSigmaDeDx; centrality", nPtBins, 0, 20, 300, -10, 10, 200, -50, 150);  
  fAntiMuonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiPionDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiKaonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  fAntiProtonDeDxCent->GetXaxis()->Set(nPtBins, binsPt);
  //
  // create feed-down correction histograms
  //
  fDCAPion = new TH3F("fDCAPion", "fDCAPion", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 );                            
  fDCAAntiPion = new TH3F("fDCAAntiPion", "fDCAAntiPion", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 );                        
  fDCAKaon = new TH3F("fDCAKaon", "fDCAKaon", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 );                            
  fDCAAntiKaon = new TH3F("fDCAAntiKaon", "fDCAAntiKaon", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 );                        
  fDCAProton = new TH3F("fDCAProton", "fDCAProton", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 );                          
  fDCAAntiProton = new TH3F("fDCAAntiProton", "fDCAAntiProton", nPtBins, 0, 20, 600, -3, 3, 200, -50, 150 ); 
  fDCAPion->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiPion->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAKaon->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiKaon->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAProton->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiProton->GetXaxis()->Set(nPtBins, binsPt); 
  //
  // create feed-down correction histograms for MC
  //
  fDCAPionMC = new TH3F("fDCAPionMC", "fDCAPionMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 );                            
  fDCAAntiPionMC = new TH3F("fDCAAntiPionMC", "fDCAAntiPionMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 );                        
  fDCAKaonMC = new TH3F("fDCAKaonMC", "fDCAKaonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 );                            
  fDCAAntiKaonMC = new TH3F("fDCAAntiKaonMC", "fDCAAntiKaonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 );                        
  fDCAProtonMC = new TH3F("fDCAProtonMC", "fDCAProtonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 );                          
  fDCAAntiProtonMC = new TH3F("fDCAAntiProtonMC", "fDCAAntiProtonMC;#it{p}_{T};DCA_{xy};Production mode", nPtBins, 0, 20, 600, -3, 3, 4, -0.5, 3.5 ); 
  fDCAPionMC->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiPionMC->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAKaonMC->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiKaonMC->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAProtonMC->GetXaxis()->Set(nPtBins, binsPt); 
  fDCAAntiProtonMC->GetXaxis()->Set(nPtBins, binsPt); 
  //add histograms to output
  fOutputList->Add(fHistZv);
  fOutputList->Add(fHistdEdxData);
  fOutputList->Add(fHistTof);
  fOutputList->Add(fHistCentBefEvSel);
  fOutputList->Add(fHistCent);
  //
  fOutputList->Add(fHistMuonGen);
  fOutputList->Add(fHistMuonReco);
  fOutputList->Add(fHistMuonRecoTOF);
  //
  fOutputList->Add(fHistPionGen);
  fOutputList->Add(fHistPionReco);
  fOutputList->Add(fHistPionRecoTOF);
  //
  fOutputList->Add(fHistKaonGen);
  fOutputList->Add(fHistKaonReco);
  fOutputList->Add(fHistKaonRecoTOF);
  //
  fOutputList->Add(fHistProtGen);
  fOutputList->Add(fHistProtReco);
  fOutputList->Add(fHistProtRecoTOF);
  //
  fOutputList->Add(fHistAntiMuonGen);
  fOutputList->Add(fHistAntiMuonReco);
  fOutputList->Add(fHistAntiMuonRecoTOF);
  //
  fOutputList->Add(fHistAntiPionGen);
  fOutputList->Add(fHistAntiPionReco);
  fOutputList->Add(fHistAntiPionRecoTOF);
  //
  fOutputList->Add(fHistAntiKaonGen);
  fOutputList->Add(fHistAntiKaonReco);
  fOutputList->Add(fHistAntiKaonRecoTOF);
  //
  fOutputList->Add(fHistAntiProtGen);
  fOutputList->Add(fHistAntiProtReco);
  fOutputList->Add(fHistAntiProtRecoTOF);  
  //
  fOutputList->Add(fMuonDeDxCent);
  fOutputList->Add(fPionDeDxCent);
  fOutputList->Add(fKaonDeDxCent);
  fOutputList->Add(fProtonDeDxCent);
  //
  fOutputList->Add(fAntiMuonDeDxCent);
  fOutputList->Add(fAntiPionDeDxCent);
  fOutputList->Add(fAntiKaonDeDxCent);
  fOutputList->Add(fAntiProtonDeDxCent);
  //
  fOutputList->Add(fDCAPion);
  fOutputList->Add(fDCAAntiPion);
  fOutputList->Add(fDCAKaon);
  fOutputList->Add(fDCAAntiKaon);
  fOutputList->Add(fDCAProton);
  fOutputList->Add(fDCAAntiProton);
  //
  fOutputList->Add(fDCAPionMC);
  fOutputList->Add(fDCAAntiPionMC);
  fOutputList->Add(fDCAKaonMC);
  fOutputList->Add(fDCAAntiKaonMC);
  fOutputList->Add(fDCAProtonMC);
  fOutputList->Add(fDCAAntiProtonMC);
  //introduce track cuts
  //
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  //
  fESDtrackCutsWoDCA = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  fESDtrackCutsWoDCA->SetEtaRange(-0.8,0.8);
  //
  //
  // 
  PostData(1, fOutputList);
 
}



//____________________________________________________________________________________//
void AliAnalysisTaskTpcSpectra::UserExec(Option_t*) {
  //
  // loop over events
  //
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler) fPIDResponse = inputHandler->GetPIDResponse();
  } 
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD) return;
  bool EvtAcc = fEventCut.AcceptEvent(fESD);
  if(!fEventCut.PassedCut(AliEventCuts::kTrigger)) return;
  //
  // centrality
  //
  Float_t centrality = -999;                                          
  fMultSel = (AliMultSelection*)fESD->FindListObject("MultSelection");
  if (!fMultSel) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    AliDebug(2, "Estimating centrality");
    centrality = fMultSel->GetMultiplicityPercentile("V0M", kFALSE); //Event selection is embedded in the Multiplicity estimator so that the Multiplicit                                                                      y percentiles are well defined and refer to the same sample
  }
  //Fill centrality histogram before event selection
  fHistCentBefEvSel->Fill(centrality);
  if(!EvtAcc) return;
  //
  // primary vertex
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
  //Fill centrality histogram after event selection
  fHistCent->Fill(centrality);
  //
  // MONTE CARLO -- GENERATED PARTICLES
  //
  Bool_t mcTrue = kFALSE;
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandler) mcTrue = kTRUE;
  //
  AliMCEvent* mcEvent = 0x0;  
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (mcEvent) {
    //
    // loop over generated particles
    //
    for (Int_t i = 0; i < mcEvent->GetNumberOfTracks(); i++) { // start loop generated particles
      TParticle* trackMC = ((AliMCParticle*)mcEvent->GetTrack(i))->Particle();
      if (!trackMC) continue;
      Int_t mcCode = trackMC->GetPdgCode();
      if (mcEvent->IsPhysicalPrimary(i) && TMath::Abs(trackMC->Y()) < 0.5) {
	switch (mcCode) { 
	case 13 : fHistMuonGen->Fill(trackMC->Pt()); break; //muon
	case 211 : fHistPionGen->Fill(trackMC->Pt()); break; //pion
	case 321 : fHistKaonGen->Fill(trackMC->Pt()); break; //kaon
	case 2212 : fHistProtGen->Fill(trackMC->Pt()); break; //proton
	case -13 : fHistAntiMuonGen->Fill(trackMC->Pt()); break; //Antimuon
	case -211 : fHistAntiPionGen->Fill(trackMC->Pt()); break; //Antipion  
	case -321 : fHistAntiKaonGen->Fill(trackMC->Pt()); break; //Antikaon
	case -2212 : fHistAntiProtGen->Fill(trackMC->Pt()); break; //antiproton
	}
      }
    } // end loop generated particles
  }
  //
  // RECONSTRUCTED TRACKS
  //
  Int_t jTracks = fESD->GetNumberOfTracks();
  Float_t dcaxy, dcaz; 
  Float_t Prod = 0;
  for (Int_t j = 0; j < jTracks; j++) { // start loop over tracks
    //
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(j));
    if (!track) continue;
    if (!fESDtrackCutsWoDCA->AcceptTrack(track)) continue; // check if track passes the cuts
    //Get the DCAxy and DCAz
    track->GetImpactParameters(dcaxy, dcaz);
    //
    // calculate rapidities and kinematics
    //
    Double_t pvec[3];
    track->GetPxPyPz(pvec);
    Double_t energyMuon = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kMuon)*AliPID::ParticleMass(AliPID::kMuon));
    Double_t energyPion = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion));
    Double_t energyKaon = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon));
    Double_t energyProton = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
    //
    Double_t rapMuon = 0.5*TMath::Log((energyMuon + pvec[2])/(energyMuon - pvec[2]));
    Double_t rapPion = 0.5*TMath::Log((energyPion + pvec[2])/(energyPion - pvec[2]));
    Double_t rapKaon = 0.5*TMath::Log((energyKaon + pvec[2])/(energyKaon - pvec[2]));
    Double_t rapProton = 0.5*TMath::Log((energyProton + pvec[2])/(energyProton - pvec[2]));
    if(track->Charge() > 0){
        if (TMath::Abs(rapPion) < 0.5 && GetCombinedSigma(track, AliPID::kPion) < 2) fDCAPion->Fill(track->Pt(), dcaxy, centrality);
        if (TMath::Abs(rapKaon) < 0.5 && GetCombinedSigma(track, AliPID::kKaon) < 2) fDCAKaon->Fill(track->Pt(), dcaxy, centrality);
        if (TMath::Abs(rapProton) < 0.5 && GetCombinedSigma(track, AliPID::kProton) < 2) fDCAProton->Fill(track->Pt(), dcaxy, centrality);
    }
    else if(track->Charge() < 0){  
        if (TMath::Abs(rapPion) < 0.5 && GetCombinedSigma(track, AliPID::kPion) < 2) fDCAAntiPion->Fill(track->Pt(), dcaxy, centrality);
        if (TMath::Abs(rapKaon) < 0.5 && GetCombinedSigma(track, AliPID::kKaon) < 2) fDCAAntiKaon->Fill(track->Pt(), dcaxy, centrality);
        if (TMath::Abs(rapProton) < 0.5 && GetCombinedSigma(track, AliPID::kProton) < 2) fDCAAntiProton->Fill(track->Pt(), dcaxy, centrality);
    }
    if (mcEvent) {
      TParticle* trackReco = ((AliMCParticle*)mcEvent->GetTrack(abs(track->GetLabel())))->Particle();
      Int_t recoCode = trackReco->GetPdgCode();
      if (TMath::Abs(trackReco->Y()) < 0.5) {
        Prod = 3;
        if (mcEvent->IsPhysicalPrimary(abs(track->GetLabel())))
          Prod = 0;
        else if (mcEvent->IsSecondaryFromWeakDecay(abs(track->GetLabel())))
          Prod = 1;
        else if (mcEvent->IsSecondaryFromMaterial(abs(track->GetLabel())))
          Prod = 2;
        //
        if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kPion)) {
          if (recoCode > 0)
            fDCAPionMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiPionMC->Fill(track->Pt(), dcaxy, Prod);
        } else if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kKaon)) {
          if (recoCode > 0)
            fDCAKaonMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiKaonMC->Fill(track->Pt(), dcaxy, Prod);
        } else if (TMath::Abs(recoCode) == AliPID::ParticleCode(AliPID::kProton)) {
          if (recoCode > 0)
            fDCAProtonMC->Fill(track->Pt(), dcaxy, Prod);
          else
            fDCAAntiProtonMC->Fill(track->Pt(), dcaxy, Prod);
        }
      }
    }
    if (!fESDtrackCuts->AcceptTrack(track)) continue; // check if track passes the cuts
    if (!track->GetInnerParam()) continue;            // check if track is a proper TPC track
    //
    Double_t ptot = track->GetInnerParam()->GetP();   // momentum for dEdx determination
    Double_t sign = track->GetSign();                 // charge
    //
    // fill the QA dE/Dx histograms
    //
    if (track->GetTPCsignalN() > 80) fHistdEdxData->Fill(ptot*sign,track->GetTPCsignal());    
    //
    // fill the raw spectra
    //
    if(track->Charge() > 0){
        if (TMath::Abs(rapMuon) < 0.5 ) fMuonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon), centrality);
        if (TMath::Abs(rapPion) < 0.5 ) fPionDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion), centrality);
        if (TMath::Abs(rapKaon) < 0.5 ) fKaonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon), centrality);
        if (TMath::Abs(rapProton) < 0.5 ) fProtonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton), centrality);
    }
    else if(track->Charge() < 0){  
        if (TMath::Abs(rapMuon) < 0.5 ) fAntiMuonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon), centrality);
        if (TMath::Abs(rapPion) < 0.5 ) fAntiPionDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion), centrality);
        if (TMath::Abs(rapKaon) < 0.5 ) fAntiKaonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon), centrality);
        if (TMath::Abs(rapProton) < 0.5 ) fAntiProtonDeDxCent->Fill(track->Pt(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton), centrality);
    }
    //
    //Check TOF status
    //
    UInt_t status = 0;
    status = track->GetStatus();
    //  
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout;
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && hasTOFtime) hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength();
    if (length < 350.) hasTOF = kFALSE;
    //
    // access mcTruth for efficiency calculation
    //
    if (track->GetTPCsignal() < 10.0) continue; // make sure that PID efficiency is taken into account
    if (mcEvent) {
      TParticle* trackReco = ((AliMCParticle*)mcEvent->GetTrack(abs(track->GetLabel())))->Particle();
      Int_t recoCode = trackReco->GetPdgCode();
      //Without TOF
      if (mcEvent->IsPhysicalPrimary(abs(track->GetLabel())) && abs(trackReco->Y()) < 0.5) {
        switch (recoCode) {
        case 13:
          fHistMuonReco->Fill(trackReco->Pt());
          break; //muon
        case 211:
          fHistPionReco->Fill(trackReco->Pt());
          break; //pion
        case 321:
          fHistKaonReco->Fill(trackReco->Pt());
          break; //kaon
        case 2212:
          fHistProtReco->Fill(trackReco->Pt());
          break; //proton
        case -13:
          fHistAntiMuonReco->Fill(trackReco->Pt());
          break; //Antimuon
        case -211:
          fHistAntiPionReco->Fill(trackReco->Pt());
          break; //Antipion
        case -321:
          fHistAntiKaonReco->Fill(trackReco->Pt());
          break; //Antikaon
        case -2212:
          fHistAntiProtReco->Fill(trackReco->Pt());
          break; //Antiproton
        }
      }
      //With TOF
      if (mcEvent->IsPhysicalPrimary(abs(track->GetLabel())) && abs(trackReco->Y()) < 0.5 && track->GetNumberOfTPCClusters() > 70 && hasTOF) {
	switch (recoCode) {
	case 13 : fHistMuonRecoTOF->Fill(trackReco->Pt()); break; //muon
	case 211 : fHistPionRecoTOF->Fill(trackReco->Pt()); break; //pion
	case 321 : fHistKaonRecoTOF->Fill(trackReco->Pt()); break; //kaon
	case 2212 : fHistProtRecoTOF->Fill(trackReco->Pt()); break; //proton
	case -13 : fHistAntiMuonRecoTOF->Fill(trackReco->Pt()); break; //Antimuon
	case -211 : fHistAntiPionRecoTOF->Fill(trackReco->Pt()); break; //Antipion
	case -321 : fHistAntiKaonRecoTOF->Fill(trackReco->Pt()); break; //Antikaon
	case -2212 : fHistAntiProtRecoTOF->Fill(trackReco->Pt()); break; //Antiproton
	}
      }
    } // end access to MC truth

  } // end track loop
      

  //
  // post the data and end the event loop
  //
  PostData(1, fOutputList);


}

//____________________________________________________________________________________//

void AliAnalysisTaskTpcSpectra::Terminate(Option_t*)
{
  
}

//____________________________________________________________________________________//
							
