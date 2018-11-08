#ifndef ALIANALYSISTASKSIMPLETREEMAKER_H
#define ALIANALYSISTASKSIMPLETREEMAKER_H
class TH1F;
class TList;
class TH2D;
class TH3D;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronPID.h"
#include "AliAnalysisFilter.h"
#include "AliDielectronEventCuts.h"
#ifndef ALIANALYSISTASKSE_H
#endif

/*************** Tree Maker Class **********************
*Designed to create simple trees for use in training   *
*MVA methods                                           *
*                                                      *
* Created: 05.10.2016                                  *
* Authors: Aaron Capon      (aaron.capon@cern.ch)      *
*          Sebastian Lehner (sebastian.lehner@cern.ch) *
*                                                      *
*******************************************************/


class AliAnalysisTaskSimpleTreeMaker : public AliAnalysisTaskSE {

  public:
		AliAnalysisTaskSimpleTreeMaker(const char *name);
		AliAnalysisTaskSimpleTreeMaker();
		~AliAnalysisTaskSimpleTreeMaker();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   FinishTaskOutput();
		virtual void   Terminate(Option_t *);

		void SetupTrackCuts(AliDielectronCutGroup* finalTrackCuts);
		void SetupEventCuts(AliDielectronEventCuts* finalEventCuts);
	
		//PID calibration function to correct the width and mean of detector
		//response (I.e should be unit guassian)
		void SetCorrWidthMeanTPC(TH3D* width, TH3D* mean){
			fMeanTPC  = mean;
			fWidthTPC = width;
		};
		
		void SetCorrWidthMeanITS(TH3D* width, TH3D* mean){
			fMeanITS  = mean;
			fWidthITS = width;
		};
		
		void SetCorrWidthMeanTOF(TH3D* width, TH3D* mean){
			fMeanTOF  = mean;
			fWidthTOF = width;
		};

		void SetCentralityPercentileRange(Float_t min, Float_t max){
				fCentralityPercentileMin = min;
				fCentralityPercentileMax = max;
		}

		void SetPtRange(Float_t min, Float_t max){
				fPtMin = min;
				fPtMax = max;
		}

		void SetEtaRange(Float_t min, Float_t max){
				fEtaMin = min;
				fEtaMax = max;
		}

		//Set inclusive electron PID cuts
		void SetESigRangeITS(Float_t min, Float_t max){
				fPIDcutITS = kTRUE;
				fESigITSMin = min;
				fESigITSMax = max;
		}

		void SetESigRangeTPC(Float_t min, Float_t max){
				fESigTPCMin = min;
				fESigTPCMax = max;
		}

		void SetESigRangeTOF(Float_t min, Float_t max){
				fPIDcutTOF = kTRUE;
				fESigTOFMin = min;
				fESigTOFMax = max;
		}

		//Set pion PID exclusion cut
		void SetPSigRangeTPC(Float_t min, Float_t max){
				fPionPIDcutTPC = kTRUE;
				fPSigTPCMin = min;
				fPSigTPCMax = max;
		}

		void SetMC(Bool_t answer){ hasMC = answer; }

		//Track cut setters. StandardITSTPC2011 cuts used if nothing specified
		void SetTPCminClusters(Int_t number){
				fESDtrackCuts->SetMinNClustersTPC(number);
		}

		void SetTPCminCrossedRows(Int_t number){
				fESDtrackCuts->SetMinNCrossedRowsTPC(number);
		}

		void SetTPCRatio(Double_t number){
				fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(number);
		}

		void SetTPCChi2PerCluster(Float_t number){
				fESDtrackCuts->SetMaxChi2PerClusterTPC(number);
		}

		void SetITSclusterRequirements(AliESDtrackCuts::Detector detector, AliESDtrackCuts::ITSClusterRequirement requirement){
				fESDtrackCuts->SetClusterRequirementITS(detector, requirement);
		}

		void SetITSChi2PerCluster(Float_t number){
				fESDtrackCuts->SetMaxChi2PerClusterITS(number);
		}

		void SetMaxDCAxy(Float_t number){
				fESDtrackCuts->SetMaxDCAToVertexXY(number);
		}

		void SetMaxDCAPtDep(const char* dist){
				fESDtrackCuts->SetMaxDCAToVertexXYPtDep(dist);
		}

		void SetMaxDCAz(Double_t number){
				fESDtrackCuts->SetMaxDCAToVertexZ(number);
		 }

		void SetTPCconstrainedGlobalChi2(Double_t number){
				fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(number);
		}

		void SetGridPID(std::string string){
				fGridPID = std::stoi(string);
		}
		void GRIDanalysis(Bool_t answer){
				fIsGRIDanalysis = answer;
		}
		void createV0tree(Bool_t answer){
				fIsV0tree = answer;
		}
		void setSDDstatus(Bool_t answer){
				fHasSDD = answer;
		}
		Bool_t isV0daughterAccepted(AliVTrack* track);

		void setFilterBitSelection(Int_t filterBit){
				fFilterBit = filterBit;
		}
		void analyseMC(Bool_t answer){
			hasMC = answer;
		}

		void SetUseTPCcorr(Bool_t answer){
			fUseTPCcorr = answer;
		}

		void SetUseITScorr(Bool_t answer){
			fUseITScorr = answer;
		}
		
		void SetUseTOFcorr(Bool_t answer){
			fUseTOFcorr = answer;
		}

		Bool_t GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0);

		// Check if the generator is on the list of generators
		// If found, assign track with integer value correspding to generator
		//0 = gen purp, 1=Pythia CC_1, 2= Pythia BB_1, 3=Pythia B_1, 4=Jpsi2ee_1, 5=B2Jpsi2ee_1";
		Int_t CheckGenerator(Int_t trackID);

  private:
 
		Int_t IsEventAccepted(AliVEvent* event);
		//Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
	
		AliAnalysisTaskSimpleTreeMaker(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

		AliAnalysisTaskSimpleTreeMaker& operator=(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

		Bool_t hasMC;

		AliESDtrackCuts* fESDtrackCuts;

		AliPIDResponse* fPIDResponse; //! PID response object
		AliMCEvent* fMCevent; //!

		TTree* fTree;
	
		//Dielectron cut classes needed to source cuts from LMEE cut libraries
		//The desired cut library should be specified in the AddTask
		AliDielectronEventCuts* eventCuts;
		AliAnalysisFilter* eventFilter;
		
		AliDielectronVarCuts* varCuts;
		AliDielectronTrackCuts *trackCuts;
		AliDielectronPID *pidCuts;
		AliDielectronCutGroup* cuts;
		AliAnalysisFilter* trackFilter;

		//Class needed to use PID within the Dielectron Framework
		AliDielectronVarManager* varManager;

		// TTree branch variables
		// Event variables
		Float_t primaryVertex[3];
		Float_t multiplicityV0A;
		Float_t multiplicityV0C;
		Float_t multiplicityCL1;
		Int_t runNumber;
		Int_t event;
		// Reconstructed
		Float_t pt;
		Float_t eta;
		Float_t phi;
		Float_t nTPCclusters;
		Float_t nTPCcrossed;
		Float_t fTPCcrossOverFind;
		Float_t nTPCfindable;
		TBits tpcSharedMap;
		Float_t nTPCshared;
		Float_t chi2TPC;
		Float_t DCA[2];
		Int_t nITS;
		Float_t chi2ITS;
		Float_t fITSshared;
		Bool_t SPDfirst;
		Int_t charge;
		Float_t EnSigmaITS;
		Float_t EnSigmaITScorr;
		Float_t EnSigmaTPC;
		Float_t EnSigmaTPCcorr;
		Float_t EnSigmaTOF;
		Float_t EnSigmaTOFcorr;
		Float_t PnSigmaTPC;
		Float_t PnSigmaITS;
		Float_t PnSigmaTOF;
		Float_t KnSigmaITS;
		Float_t KnSigmaTPC;
		Float_t KnSigmaTOF;
		Float_t ITSsignal;
		Float_t TPCsignal;
		Float_t TOFsignal;
		Float_t goldenChi2;
		//MC 
		Float_t mcEta;
		Float_t mcPhi;
		Float_t mcPt;
		Float_t mcVert[3];
		Int_t iPdg;
		Int_t iPdgMother;
		Bool_t HasMother;
		Int_t motherLabel;
		Int_t isInj; // If track is injected 
		// Pdg and label for initial particle in decay chain
		Int_t iPdgFirstMother;
		Int_t gLabelFirstMother;
		Int_t gLabelMinFirstMother;
		Int_t gLabelMaxFirstMother;
		//V0 features
		Float_t pointingAngle;
		Float_t daughtersDCA;
		Float_t decayLength;
		Float_t v0mass;
		Float_t ptArm;
		Float_t alpha;

		TH1F* fQAhist; //!
		// Currently no cut on centrality
		Float_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
		Float_t fCentralityPercentileMax;// maximum centrality threshold (default = 100)

		Float_t fPtMin;// minimum pT threshold (default = 0)
		Float_t fPtMax;// maximum pT threshold (default = 10)
		Float_t fEtaMin;// minimum eta threshold (default = -0.8)
		Float_t fEtaMax;// maximum eta threshold (default = 0.8)

		//Values and flags for PID cuts in ITS and TOF
		Float_t fESigITSMin;
		Float_t fESigITSMax;
		Float_t fESigTPCMin; 
		Float_t fESigTPCMax; 
		Float_t fESigTOFMin;
		Float_t fESigTOFMax;
		
		Bool_t fPIDcutITS;
		Bool_t fPIDcutTOF;
		
		//Values and flag for pion PID cuts in TPC
		Bool_t fPionPIDcutTPC;
		Float_t fPSigTPCMin;
		Float_t fPSigTPCMax;
		
		Bool_t fHasSDD;

		Bool_t fIsV0tree;
		TH2F* fArmPlot; //!

		Bool_t fIsAOD;
		Int_t fFilterBit;
		
		//Grid PID
		Bool_t fIsGRIDanalysis;
		Int_t fGridPID;
				
		Bool_t fUseTPCcorr;
		TH3D* fWidthTPC;
		TH3D* fMeanTPC;

		Bool_t fUseITScorr;
		TH3D* fWidthITS;
		TH3D* fMeanITS;
		
		Bool_t fUseTOFcorr;
		TH3D* fWidthTOF;
		TH3D* fMeanTOF;

		// Temp variable (needed for testing MC issues)
		Int_t TOFstartMask;

		// Store list of generator hashes which can be checked against to determine
		// whether or not the track was injected
		std::vector<UInt_t> fGeneratorHashes;
		ClassDef(AliAnalysisTaskSimpleTreeMaker, 5); //

};



#endif

