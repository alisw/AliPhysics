#ifndef ALIANALYSISTASKJETPP_H
#define ALIANALYSISTASKJETPP_H


class TH1I;
class TH1F;
class TH2F;
class TH2D;
class TH1D;
class TArrayD;
class THnSparse;
class TProfile;
class TList;
class AliEmcalList;
class TClonesArray;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVParticle;
class AliLog;
class AliAnalysisUtils;
class TRandom3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include <vector>
#include "AliAnalysisTaskEmcalJet.h"

using std::vector;



class AliAnalysisTaskJetPP : public AliAnalysisTaskEmcalJet {
   public:

   enum MyContainer {
      kContainerOne = 0, 
      kContainerTwo = 1,
      kContainerThree = 2,
      kContainerFour = 3
   };



   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskJetPP();
   AliAnalysisTaskJetPP(const char *name);
   virtual ~AliAnalysisTaskJetPP();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   
  // ######### SETTERS/GETTERS

   void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}  
   void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;} 
   void        SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;} 
   void 	      SetMC(Bool_t a) {fIsMC=a;}

   void        SetAcceptanceWindows(Double_t trackEta, Double_t signalJetRadius){ 
                 fTrackEtaWindow  = trackEta; 
                 fSignalJetRadius = signalJetRadius; 
                 fSignalJetEtaWindow = fTrackEtaWindow-fSignalJetRadius; 
              } 

   void        SetCentralityType(const char* type){
                  fCentralityType = type;    
              } 

   void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
   void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}


   Bool_t   RetrieveEventObjects();
   Bool_t   Run();
   Bool_t   FillHistograms();

   private:

   Double_t EstimateBgKT(AliJetContainer *jetCont);  // median p/A of kt jets
   Double_t EstimateLocalBg(AliJetContainer *jetCont, AliParticleContainer *trkCont);  // local bg in a perpendicular cone

   Double_t Convert(Double_t input);
   // ######### CHECK FUNCTIONS
   Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isprimary=0);  
   Bool_t      IsEventInAcceptance(AliVEvent* event);     
   Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet, Bool_t suppressGhost=1); 


   // ######### STANDARD FUNCTIONS
   void      ExecOnceLocal();                    

   // ########## USAGE TRIGGERS 

   Bool_t              fUseDefaultVertexCut;   // trigger if automatic vertex cut from helper class should be done 
   Bool_t              fUsePileUpCut;          // trigger if pileup cut should be done
   Bool_t              fIsMC;          // Use MC

   // ########## JET/DIJET/RC PROPERTIES
   Double_t            fSignalJetRadius;       // Radius for the signal jets
   // ########## CUTS 
   Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets 
   Double_t            fTrackEtaWindow;        //gc +- window in eta for tracks  
   Double_t            fMinTrackPt;            //gc Min track pt to be accepted  
   Double_t            fMinJetArea;            // Min jet area to be accepted
   TString             fCentralityType;        //gc Used centrality estimate (V0A, V0C, V0M, ...) 

   Double_t            fZVertexCut;            //zvertex cut

   // ########## GENERAL ////VARS
   AliAnalysisUtils*   fHelperClass;           //! gc Vertex selection helper
   Bool_t              fInitializedLocal;           //! gc trigger if tracks/jets are loaded  initiates calling   ExecOnce 


   TH1D*  fhJetPt; //! pt spectrum of AKT jets
   TH1D* fhKTJetPt; //! ptspectrum of KT jets 
   TH1D*  fhJetPtRho; //! pt spectrum of AKT jets without kt bgk
   TH1D*  fhJetPtConeRho; //! pt spectrum of AKT jets without local bgk

   TH1D*  fhAtimesRho; //! Jet background times area
   TH1D*  fhCuts; //! Histogram for pilup/vertex cuts 
   TH1D* fhJetConstituentPt; //! Jet constituent pt
   TH1D* fhTrackPt; //! Track pt

   TH1D* fhRho; //! KT bgk
   TH1D* fhConeRho; //! local KT bgk
   TH2D* fhJetAreaPt; //! Jet Area-pt distribution
   TH2D* fhJetEtaPt; //! Jet eta-pt distribution

   TH2D* fhAktJetEtaPhi; //! Jet et-phi distribution
   TH2D* fhKtJetEtaPhi; //! Jet et-phi distribution
   TH2D* fhTrackEtaPt; //! Track eta-pt distribution
   TH2D* fhGenTrackEtaPt; //! Generated track eta-pt distribution
   TH2D* fhJetPhiPt; //! Jet phi-pt distribution

   TH2D* fhTrackPhiPt; //! Track phi-pt distribution
   TH2D* fhTrackEtaPhi; //! Track eta-phi distribution
   TH2D* fhRemx; //! Response matrix
   TH1D* fhZVertex; //! Z vertex

   TH1D* fhZVertexBC; //! Z vertex before cut
   TH1D* fhXVertex; //! X vertex
   TH1D* fhYVertex; //! Y vertex
   TH1D* fhPrimGenTrkPt; //! Pt spectrum of MC particle

   TH1D* fhGenJetPt; //! Pt spectrum of MC jets
   TH1D* fhRecTrkPt; //! Pt spectrum of correctly reconstructed tracks
   TH1D* fhFakeTrkPt; //! Pt spectrum of fake tracks
   TH1D* fhMult; //! Multiplicity

   TH2D* fhTrackPhiCG; //! Global tracks
   TH2D* fhTrackPhiTPCG; //! complementary tracks
   TH2D* fhInvPtQVsPhi[2];//! Inverse pt vs phi
   TH2D* fhInvPtQVsEta[2];//! Inverse pt vs eta

   TH2D* fhInvPtQVsPhiASide[2];//! Inverse pt vs phi in away-side rapidity
   TH2D* fhInvPtQVsPhiCSide[2];//! Inverse pt vs pi in close-side rapidity
   TH2D* fhSigmaPtOverPtVsPt[2];//! Pt vs sigmaPt/Pt


   AliAnalysisTaskJetPP(const AliAnalysisTaskJetPP&);
   AliAnalysisTaskJetPP& operator=(const AliAnalysisTaskJetPP&);

   ClassDef(AliAnalysisTaskJetPP, 2); // Charged jet analysis for pA,,, increase the last parametre by 1 after a change !!!!!!!!!!!!

};
#endif
