/*An analysis task that plots the TPC dE/dx spectra for Protons, Kaons, Pions and Antiprotons for GEANT3/4*/
/*A. Kalweit and S. Marium*/
/*26.01.2018*/

#ifndef AliAnalysisTaskTpcSpectra_H
#define AliAnalysisTaskTpcSpectra_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "TMath.h"


class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliMultSelection;
class AliEventCuts;


class AliAnalysisTaskTpcSpectra : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskTpcSpectra();
  AliAnalysisTaskTpcSpectra(const char* name);
  virtual ~AliAnalysisTaskTpcSpectra();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

 private:
  AliESDEvent* fESD;                          //! input event
  TList* fOutputList;                         //! output list
  AliEventCuts fEventCut;                     //! input event selection
  AliESDtrackCuts* fESDtrackCuts;             //! input track cuts
  AliESDtrackCuts* fESDtrackCutsWoDCA;        //! input track cuts without the cut on the DCA
  AliPIDResponse*        fPIDResponse;        //! pid response object
  AliMultSelection*      fMultSel;            //! centrality and multiplicity selection
  //
  // histograms
  //
  TH1F* fHistZv;                            //! z-Vertex distribution
  TH2F* fHistdEdxData;                      //! PID histogram dEdx all particles
  TH2F* fHistTof;                           //! PID histofram TOF all particles
  TH1F* fHistCentBefEvSel;                  //! centrality histogram before event selection, after physics selection
  TH1F* fHistCent;                          //! centrality histogram
  //
  TH1F* fHistMuonGen;                       //! generated muon tpc only
  TH1F* fHistMuonReco;                      //! reconstructed muon tpc only
  TH1F* fHistMuonRecoTOF;                   //! reconstructed muon tof extrapolation
  //
  TH1F* fHistPionGen;                       //! generated pion tpc only
  TH1F* fHistPionReco;                      //! reconstructed pion tpc only
  TH1F* fHistPionRecoTOF;                   //! reconstructed pion tof extrapolation
  //
  TH1F* fHistKaonGen;                       //! generated kaon tpc only
  TH1F* fHistKaonReco;                      //! reconstructed kaon tpc only
  TH1F* fHistKaonRecoTOF;                   //! reconstructed kaon tof extrapolation
  //
  TH1F* fHistProtGen;                       //! generated proton tpc only
  TH1F* fHistProtReco;                      //! reconstructed proton tpc only
  TH1F* fHistProtRecoTOF;                   //! reconstructed proton tof extrapolation
  //
  TH1F* fHistAntiMuonGen;                   //! generated antimuon tpc only
  TH1F* fHistAntiMuonReco;                  //! reconstructed antimuon tpc only
  TH1F* fHistAntiMuonRecoTOF;               //! reconstructed antimuon tof extrapolation
  //
  TH1F* fHistAntiPionGen;                   //! generated antipion tpc only
  TH1F* fHistAntiPionReco;                  //! reconstructed antipion tpc only
  TH1F* fHistAntiPionRecoTOF;               //! reconstructed antipion tof extrapolation
  //
  TH1F* fHistAntiKaonGen;                   //! generated antikaon tpc only
  TH1F* fHistAntiKaonReco;                  //! reconstructed antikaon tpc only
  TH1F* fHistAntiKaonRecoTOF;               //! reconstructed antikaon tof extrapolation
  //
  TH1F* fHistAntiProtGen;                   //! generated antiproton tpc only
  TH1F* fHistAntiProtReco;                  //! reconstructed antiproton tpc only
  TH1F* fHistAntiProtRecoTOF;               //! reconstructed antiproton tof extrapolation
  //
  TH3F* fMuonDeDxCent;                       //! pT histogram for muons
  TH3F* fPionDeDxCent;                       //! pT histogram for pions
  TH3F* fKaonDeDxCent;                       //! pT histogram for Kaons
  TH3F* fProtonDeDxCent;                     //! pT histogram for Protons
  //
  TH3F* fAntiMuonDeDxCent;                   //! pT histogram for Antimuons
  TH3F* fAntiPionDeDxCent;                   //! pT histogram for Antipions
  TH3F* fAntiKaonDeDxCent;                   //! pT histogram for AntiKaons
  TH3F* fAntiProtonDeDxCent;                 //! pT histogram for AntiProtons
  //
  TH3F* fDCAPion;                            //! DCA,pT,Centrality histogram for pions
  TH3F* fDCAAntiPion;                        //! DCA,pT,Centrality histogram for pions
  TH3F* fDCAKaon;                            //! DCA,pT,Centrality histogram for Kaons
  TH3F* fDCAAntiKaon;                        //! DCA,pT,Centrality histogram for Kaons
  TH3F* fDCAProton;                          //! DCA,pT,Centrality histogram for Protons
  TH3F* fDCAAntiProton;                      //! DCA,pT,Centrality histogram for Protons
  //
  TH3F* fDCAPionMC;                            //! DCA,pT,Production mode histogram for pions
  TH3F* fDCAAntiPionMC;                        //! DCA,pT,Production mode histogram for pions
  TH3F* fDCAKaonMC;                            //! DCA,pT,Production mode histogram for Kaons
  TH3F* fDCAAntiKaonMC;                        //! DCA,pT,Production mode histogram for Kaons
  TH3F* fDCAProtonMC;                          //! DCA,pT,Production mode histogram for Protons
  TH3F* fDCAAntiProtonMC;                      //! DCA,pT,Production mode histogram for Protons
  //
  //Utility functions
  Float_t GetCombinedSigma(AliESDtrack* track, AliPID::EParticleType type) { return TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(track, type), 2) + TMath::Power(fPIDResponse->NumberOfSigmasTOF(track, type), 2)); }
  //
  AliAnalysisTaskTpcSpectra(const AliAnalysisTaskTpcSpectra&);
  AliAnalysisTaskTpcSpectra& operator=(const AliAnalysisTaskTpcSpectra&);
  //
  ClassDef(AliAnalysisTaskTpcSpectra, 1);
  
};

#endif
  
  
