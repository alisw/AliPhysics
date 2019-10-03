#ifndef AliAnalysisTaskEmcalHFCJQA_H
#define AliAnalysisTaskEmcalHFCJQA_H

class TH1;
class TH2;
class TH3;
class TH3F;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"
#include "THnSparse.h"
#include "AliRDHFJetsCuts.h"
#include "AliPIDResponse.h"


class AliAnalysisTaskEmcalHFCJQA : public AliAnalysisTaskEmcalJet {

 public:

AliAnalysisTaskEmcalHFCJQA();
AliAnalysisTaskEmcalHFCJQA(const char *name);
virtual ~AliAnalysisTaskEmcalHFCJQA();

void                        UserCreateOutputObjects();
void                        Terminate(Option_t *option);
void                        ExecOnce();
Bool_t                      FillHistograms()   ;
Bool_t                      Run()              ;
void                        CheckClusTrackMatching();

//==============================================================================
//Setters
void SetupPIDresponse();
void SetReadMC(Int_t readMC){fReadMC=readMC;}
void SetDebug(int p){fDebug = p;}
void SetFilterBit(Int_t fbit){ffilterbit=fbit;}
void SetJetCuts(AliRDHFJetsCuts *cuts){delete fCuts; fCuts=new AliRDHFJetsCuts(*cuts);}
void SetCutObject(AliRDHFJetsCuts *cuts){fCuts=cuts;}
//==============================================================================
//Getters
AliRDHFJetsCuts* GetJetCuts(){return fCuts;}  


private:
//==============================================================================
//Functions

Int_t TagJetMC();
AliAODMCParticle* GetMCPartonOrigin(TClonesArray* &arrayMC, AliAODMCParticle* &p , Int_t &idx);
AliAODMCParticle* IsMCJet(TClonesArray *arrayMC,const AliEmcalJet *jet, Double_t &contribution);
void FillJetRecoHisto(const AliEmcalJet *jet,Int_t partonnat,Double_t contribution,Double_t ptpart);
 Bool_t FillTrackHistosAndSelectTrack(AliAODTrack *aodtr, const AliESDVertex *primary, const Double_t magfield,const Bool_t isPico);
void PrintDebug(int N, TString Section="DDG", TString Sub="", int LEVEL=0);
void TriggersHistogram(TString TriggerClass);
void TriggersMaskHistogram(int kMask);
void TriggersBitHistogram(AliAODEvent* aod, bool ReadMC);
void EnergyTriggers();
void ClustersEnergyDistribution(bool isL0, bool EGA1,bool EGA2,bool EJE1,bool EJE2);
//============================================================================================================ 
//Containers
//TList						 *fOutput; 			         //! output list
AliJetContainer            *fJetsCont;                //!Jets
AliParticleContainer       *fTracksCont;              //!Tracks
AliClusterContainer        *fCaloClustersCont;        //!Clusters  
//============================================================================================================ 

//============================================================================================================ 
//Flags
Int_t 			fReadMC;                     // 0=no read mc, 1=is MC but analysis is data-like, 2=MC based analysis 
int				fDebug;

//============================================================================================================ 

//============================================================================================================ 
//Variables
Int_t 					ffilterbit;                  // filter bit
Bool_t 					fKeepTrackNegID;             //  flag for rejecting track with neg ID
AliPIDResponse* 		fpidResp;              		 // !pid response object
TString					fMCParticlesName;			 // MC Array Name
AliRDHFJetsCuts*		fCuts;
//============================================================================================================ 

//============================================================================================================ 
//Histograms
TH1F* 			fhEventCounter;
TH1F* 			fhTriggerCounter;
TH1F* 			fhTriggerMaskCounter;
TH1F* 			fhTriggerBitCounter;
TH2F*           fClustersEnergydistribution;
TH1F* 			fEventsThreshold;
THnSparseF*		fSparseRecoJets;
THnSparseF *fhSparseFilterMask;          			//! sparse histo with track information
THnSparseF *fhSparseFilterMaskPico;          			//! sparse histo with track information
THnSparseF *fhSparseFilterMaskTrackAcc;   			//! sparse with filter bits and track kine/geometrical properties
THnSparseF *fhSparseFilterMaskTrackAccPico;   			//! sparse with filter bits and track kine/geometrical properties
THnSparseF *fhSparseFilterMaskImpPar;    			//! sparse with kine/geometrical prop and imp par xy
THnSparseF *fhSparseFilterMaskImpParPico;    			//! sparse with kine/geometrical prop and imp par xy
TH3F *fhImpParResolITSsel;              			//! histo with imp par distribution as a function of pt and ITS clusters
TH3F *fhImpParResolITSselGoodTracks;    			//! histo with imp par distribution as a function of pt and ITS clusters for selected tracks

TH3F *fhnSigmaTPCTOFEle;                			//! sparse with TPC-TOF nsigma informations, centered on ele hypo
TH3F *fhnSigmaTPCTOFPion;              				//! sparse with TPC-TOF nsigma informations, centered on pion hypo
TH3F *fhnSigmaTPCTOFKaon;               			//! sparse with TPC-TOF nsigma informations, centered on kaon hypo
TH3F *fhnSigmaTPCTOFProton;             			//! sparse with TPC-TOF nsigma informations, centered on proton hypo

THnSparseF *fhSparseEoverPeleTPC;       //! sparse histo with TPC-EMCal PID electron information
THnSparseF *fhSparseShowShapeEleTPC;    //! sparse histo with TPC-EMCAL PID electron info, including shower shape & Ncells
THnSparseF *fhTrackEMCal;              //! sparse with EMCal cluster properties related to clusters matched to tracks

 TH1F* 			fhptSpectrum;      //!
 TH1F* 			fhptSpectrumPico;      //!
 Bool_t fCheckPico;                      //!
AliAnalysisTaskEmcalHFCJQA(const AliAnalysisTaskEmcalHFCJQA&);				// copy constructo not implemented yet
AliAnalysisTaskEmcalHFCJQA& operator=(const AliAnalysisTaskEmcalHFCJQA&); 	// assignment operator not implemented yet

ClassDef(AliAnalysisTaskEmcalHFCJQA, 3) // jet sample analysis task
};
#endif
