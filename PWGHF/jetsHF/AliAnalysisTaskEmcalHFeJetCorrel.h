#ifndef AliAnalysisTaskEmcalHFeJetCorrel_H
#define AliAnalysisTaskEmcalHFeJetCorrel_H

class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

#include "AliRDHFJetsCuts.h"
#include "AliPIDResponse.h"


class AliAnalysisTaskEmcalHFeJetCorrel : public AliAnalysisTaskEmcalJet {

 public:

  AliAnalysisTaskEmcalHFeJetCorrel();
  AliAnalysisTaskEmcalHFeJetCorrel(const char *name);
  virtual ~AliAnalysisTaskEmcalHFeJetCorrel();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();

//==============================================================================
//Setters
  void SetupPIDresponse();
  void SetWindowClusterTrack(Double_t dphi,Double_t deta){fPhiTrackClusterDistance=dphi;fEtaTrackClusterDistance=deta;}
  void SetMinPtElectron(Double_t minpt){fminpt=minpt;}
  void SetEtaTracks(double etamin, double etamax){fEtaMin=etamin;fEtaMax=etamax;}
  void SetNsigmaTPCelectron(double nsigmamin, double nsigmamax){fsigmaTPCmin=nsigmamin;fsigmaTPCmax=nsigmamax;}
  void SetPhotonicMassCut(Double_t mass){fMassPhotonicCut=mass;}
  void SetMinNumberOfClusterCells(Double_t ncell){fminNcell=ncell;}
  void SetMaxM20(Double_t maxm20){fmaxM20=maxm20;}
  void SetMaxM02(Double_t maxm02){fmaxM02=maxm02;}
  void SetEoverPlimits(Double_t min,Double_t max){fMinEoverP=min;fMaxEoverP=max;}
  void SetMCParticles(TString MCParticlesName){fMCParticlesName = MCParticlesName;}
  void SetReadMC(Int_t readMC){fReadMC=readMC;}
  void SetQA(bool QA){kQA=QA;}
  void SetDoAnalysis(bool p){kAnalysis=p;}
  void SetGeneralSpectra(bool p){kGeneralSpectra=p;}
  void SetDetectorsTest(bool p){kDetector=p;}
  void SetCheckClusterMatching(bool CheckClusterMatching){kCheckClusterMatching=CheckClusterMatching;}
  void SetJetCuts(AliRDHFJetsCuts *cuts){delete fCutsHFjets; fCutsHFjets=new AliRDHFJetsCuts(*cuts);}
  void SetElectronCuts(AliRDHFJetsCuts *cuts){delete fCutsElectron; fCutsElectron=new AliRDHFJetsCuts(*cuts);}
  void SetFilterBitElectron(Int_t fbit){ffilterbit=fbit;}
  void FindTPCElectrons(double TPCMax=3, double TPCMin=-3);
  void JetStudy();
  //void MCJetStudy();
  void SetMCJetsBranch(char* branch) {fMCJetsBranch = branch;}
//==============================================================================
//Getters
  AliRDHFJetsCuts* GetJetCuts(){return fCutsHFjets;}  


private:
//==============================================================================
//Functions

//void FillJetRecoHisto(const AliEmcalJet *jet,Int_t partonnat,Double_t contribution,Double_t ptpart);
  void FillCorrelationHisto(const AliVParticle *part, const double Energy, const Int_t idlab, const AliEmcalJet *jet, Double_t nsigma,THnSparseF* &h);
  void FillMCCorrelationHisto(const AliVParticle *part, const double Ep, const Int_t idlab, const AliEmcalJet *jet,Double_t nsigma, THnSparseF* &h, double IsHeavy, double JetHeavy, double Contribution);
  void FillJetHisto(const AliEmcalJet *jet);
  void FillMCJetHisto(const AliEmcalJet *jet);
  void FillMassHisto(const AliVParticle *part, const AliEmcalJet *jet, double pgamma, double mass, THnSparseF* &h);
  Int_t TagJetMC();
  AliAODMCParticle* GetMCPartonOrigin(TClonesArray* &arrayMC, AliAODMCParticle* &p , Int_t &idx);
  AliAODMCParticle* IsMCJet(TClonesArray *arrayMC,const AliEmcalJet *jet, Double_t &contribution);
  void TestContainerEvent(AliParticleContainer* &Cont, AliClusterContainer* &fCaloClustersCont, TClonesArray* &fTracks_tender, TClonesArray* &fCaloClusters_tender, AliVEvent* &Event, THnSparseF* &fTHN, THnSparseF* &fTHN2, THnSparseF* &fTHN3);
  void ClusterHisto(AliClusterContainer* &Cont, THnSparseF* &fTHN);
  //void TestHisto(const AliEmcalJet *jet, AliParticleContainer *fTRacksCont);
  void RandomCones(double RParameter = 0.4, double Rho = 0, int TMax = 1);
  double PhiInterval(double Phi, double Low = -0.5*TMath::Pi(), double Up = 1.5*TMath::Pi());
  double RPhiEta(double PhiRef, double EtaRef, double Phi, double Eta);
//============================================================================================================ 
//Containers
//TList						 *fOutput; 			         //! output list
  AliJetContainer            *fJetsCont;                //!Jets
  AliJetContainer            *fMCJetsCont;              //!MCJets
  AliParticleContainer       *fTracksCont;              //!Tracks
  AliClusterContainer        *fCaloClustersCont;        //!Clusters  
//============================================================================================================ 

//============================================================================================================ 
//Flags
  Int_t 			fReadMC;                     // 0=no read mc, 1=is MC but analysis is data-like, 2=MC based analysis 
  bool 				kQA;
  bool 				kCheckClusterMatching;		 //Flag to check the cluster-track matching phase space
  Bool_t 			fCheckVzero;                 // flag to check whether electron propagated to EMCal are included in some v0 candidate
  bool				kAnalysis;					 //flag to make the analysis plots
  bool				kGeneralSpectra;			 //flag to make pT spectra for electrons and jets
  bool				kDetector;

//============================================================================================================ 

//============================================================================================================ 
//Variables

  int				fdebug;						 //Debugging variable
  TH1F 				*fNentries;                  //! histo for event counting and checks
  TH1F              *fNRejected;                 //Event Selection Reasons
  TArrayI 			*fcandEleTPC;                //! array with electron candidates from TPC
  Int_t				fLastEleTPC;                 // number of electron candidates from TPC
  AliRDHFJetsCuts 	*fCutsHFjets;        		 // specific algo jet cut object (incude the above?)
  AliRDHFJetsCuts 	*fCutsElectron;      		 // cut object for electrn track cuts
  AliPIDResponse 	*fpidResp;            		 // pid response object
  Int_t 			ffilterbit;                  // filter bit
  Double_t 			fPhiTrackClusterDistance;    // max track-to-cluster distance in phi
  Double_t 			fEtaTrackClusterDistance;    // min track-to-cluster distance in eta
  Double_t 			fsigmaTPCmin;                // nsigma cut in tpc
  Double_t 			fsigmaTPCmax;                // nsigma cut in tpc
  Double_t 			fminpt;                      // min track pt
  Double_t 			fEtaMin;                     // min track eta
  Double_t 			fEtaMax;                     // max track eta
  Double_t 			fMinEoverP;                  // min E over p (EMCal E over TPC p)
  Double_t 			fMaxEoverP;                  // max E over p (EMCal E over TPC p)
  Double_t 			fMassPhotonicCut;            // mass cut to reject photonic electrons
  Double_t 			fminNcell;                   // min number of cell per EMCAL cluster
  Double_t 			fmaxM20;                     // max M20 (shower shape parameter)
  Double_t 			fmaxM02;                     // max M02 (shower shape parameter)
  TString			fMCParticlesName;			 // MC Array Name
  char              *fMCJetsBranch;              // name of the AOD MC-jets branch 
//============================================================================================================ 

//============================================================================================================ 
//Histograms
  TH3F 				*fhTrackRejection;           //! histo storing track rejection reason during electron ID
  TH3F 				*fhEleRejection;             //! histo storing electron (MC id) rejection
  THnSparseF 		*fSparseHFSpectrum;          //! sparse histo for ele-jet pairs properties
  THnSparseF 		*fSparseHFTriggers;          //! sparse to count the trigger electrons
  THnSparseF 		*fSparseJet;              	 //!sparse histo for jet properties
  THnSparseF 		*fSparseMCJet;               //!sparse histo for MCjet properties
  THnSparseF 		*fSparseQATracks;
  TH3F				*fHistPtDEtaDPhiClusTrack; 	 //!cluster pt, delta eta, delta phi to matched track
  TH3F				*fHistPtDEtaDPhiTrackClus;   //!track pt, delta eta, delta phi to matched cluster
  THnSparseF 		*fhPhotonicEle;              //! histo with unlike sign photonic mass candidate properties 
  THnSparseF 		*fDetector;
  THnSparseF 		*fRandomCones;
//============================================================================================================ 

  AliAnalysisTaskEmcalHFeJetCorrel(const AliAnalysisTaskEmcalHFeJetCorrel&);				// copy constructo not implemented yet
  AliAnalysisTaskEmcalHFeJetCorrel& operator=(const AliAnalysisTaskEmcalHFeJetCorrel&); 	// assignment operator not implemented yet

  ClassDef(AliAnalysisTaskEmcalHFeJetCorrel, 2) // jet sample analysis task
};
#endif
