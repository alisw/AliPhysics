#ifndef ALIANALYSISTASKSIMPLETREEMAKER_H
#define ALIANALYSISTASKSIMPLETREEMAKER_H
class TH1F;
class TList;
class TH2D;
class TH3D;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "AliAnalysisTaskSE.h"
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
        virtual ~AliAnalysisTaskSimpleTreeMaker(){} 

        virtual void   UserCreateOutputObjects();
        virtual void   UserExec(Option_t *option);
        virtual void   FinishTaskOutput();
        virtual void   Terminate(Option_t *);
        //~ 

        void SetCentralityPercentileRange(Double_t min, Double_t max){
            fCentralityPercentileMin = min;
            fCentralityPercentileMax = max;
        }

        void SetPtRange(Double_t min, Double_t max){
            fPtMin = min;
            fPtMax = max;
        }

        void SetEtaRange(Double_t min, Double_t max){
            fEtaMin = min;
            fEtaMax = max;
        }

        //Set inclusive electron PID cuts
        void SetESigRangeITS(Double_t min, Double_t max){
            fPIDcutITS = kTRUE;
            fESigITSMin = min;
            fESigITSMax = max;
        }

        void SetESigRangeTPC(Double_t min, Double_t max){
            fESigTPCMin = min;
            fESigTPCMax = max;
        }

        void SetESigRangeTOF(Double_t min, Double_t max){
            fPIDcutTOF = kTRUE;
            fESigTOFMin = min;
            fESigTOFMax = max;
        }

        //Set pion PID exclusion cut
        void SetPSigRangeTPC(Double_t min, Double_t max){
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
				void createGenTree(Bool_t answer){
					fIsEffTree = answer;
				}
        Bool_t isV0daughterAccepted(AliVTrack* track);

        void setFilterBitSelection(Int_t filterBit){
            fFilterBit = filterBit;
        }
				void analyseMC(Bool_t answer){
					hasMC = answer;
				}

        Bool_t GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0);

    private:
 
        Int_t IsEventAccepted(AliVEvent* event);
        //Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
      
        AliAnalysisTaskSimpleTreeMaker(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

        AliAnalysisTaskSimpleTreeMaker& operator=(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

				Bool_t hasMC;

        AliESDtrackCuts* fESDtrackCuts; // ESD track cuts object

        AliPIDResponse* fPIDResponse; //! PID response object

        TTree* fTree;

				//TTree branch variables
				//Event variables
				Double_t primaryVertex[3];
				Double_t multiplicityV0A;
				Double_t multiplicityV0C;
				Double_t multiplicityCL1;
				Int_t runNumber;
				Int_t event;
				//Reconstructed
				Double_t pt;
				Double_t eta;
				Double_t phi;
				Double_t nTPCclusters;
				Double_t nTPCcrossed;
				Double_t fTPCcrossOverFind;
				Double_t nTPCfindable;
				TBits tpcSharedMap;
				Double_t nTPCshared;
				Double_t chi2TPC;
				Double_t DCA[2];
				Int_t nITS;
				Double_t chi2ITS;
				Double_t fITSshared;
				Bool_t SPDfirst;
				Int_t charge;
				Double_t EnSigmaITS;
				Double_t EnSigmaTPC;
				Double_t EnSigmaTOF;
				Double_t PnSigmaTPC;
				Double_t PnSigmaITS;
				Double_t PnSigmaTOF;
				Double_t KnSigmaITS;
				Double_t KnSigmaTPC;
				Double_t KnSigmaTOF;
				Double_t ITSsignal;
				Double_t TPCsignal;
				Double_t TOFsignal;
				Double_t goldenChi2;
				//MC 
				Double_t mcEta;
				Double_t mcPhi;
				Double_t mcPt;
				Double_t mcVert[3];
				Int_t iPdg;
				Int_t iPdgMother;
				Bool_t HasMother;
				Int_t motherLabel;
				//V0 features
				Double_t pointingAngle;
				Double_t daughtersDCA;
				Double_t decayLength;
				Double_t v0mass;
				Double_t ptArm;
				Double_t alpha;

        TH1F* fQAhist;
        Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
        Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

        Double_t fPtMin;// minimum pT threshold (default = 0)
        Double_t fPtMax;// maximum pT threshold (default = 10)
        Double_t fEtaMin;// minimum eta threshold (default = -0.8)
        Double_t fEtaMax;// maximum eta threshold (default = 0.8)

        //Values and flags for PID cuts in ITS and TOF
        Double_t fESigITSMin;
        Double_t fESigITSMax;
				Double_t fESigTPCMin; 
        Double_t fESigTPCMax; 
        Double_t fESigTOFMin;
        Double_t fESigTOFMax;
        
        Bool_t fPIDcutITS;
        Bool_t fPIDcutTOF;
        
        //Values and flag for pion PID cuts in TPC
        Bool_t fPionPIDcutTPC;
        Double_t fPSigTPCMin;
        Double_t fPSigTPCMax;
        
        Bool_t fHasSDD;

        Bool_t fIsV0tree;
        TH2F* fArmPlot;

				//Efficiency calculation flags
				Bool_t fIsEffTree;

        Bool_t fIsAOD;
        Int_t fFilterBit;
        
        //Grid PID
        Bool_t fIsGRIDanalysis;
        Int_t fGridPID;
        
        ClassDef(AliAnalysisTaskSimpleTreeMaker, 4); //

};



#endif

