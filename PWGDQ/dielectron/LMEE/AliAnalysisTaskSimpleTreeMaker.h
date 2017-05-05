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

    void SetMC(Bool_t answer){ fIsMC = answer; }

    //Set to request centrality value for events
    void SetIsIonColl(Bool_t answer){ fIsIonColl = answer; }
    
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

    void setFilterBitSelection( Int_t filterBit){
        fFilterBit = filterBit;
    }
    void useAODs( Bool_t answer){
        fIsAOD = answer;
    }

  
    private:
 
    Int_t IsEventAccepted(AliVEvent *event);
    //Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  
    AliAnalysisTaskSimpleTreeMaker(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

    AliAnalysisTaskSimpleTreeMaker& operator=(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

    AliMCEvent* mcEvent;    //MC object

    AliESDtrackCuts *fESDtrackCuts; // ESD track cuts object

    AliPIDResponse *fPIDResponse; //! PID response object

    TTreeStream* fStream;
    TTree* fTree;

    TH1F* fQAhist;
    Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
    Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

    Bool_t fIsMC;
  
    Double_t fPtMin;// minimum pT threshold (default = 0)
    Double_t fPtMax;// maximum pT threshold (default = 10)
    Double_t fEtaMin;// minimum eta threshold (default = -0.8)
    Double_t fEtaMax;// maximum eta threshold (default = 0.8)

 
    Double_t fESigTPCMin; 
    Double_t fESigTPCMax; 

    //Values and flags for PID cuts in ITS and TOF
    Double_t fESigITSMin;
    Double_t fESigITSMax;
    Double_t fESigTOFMin;
    Double_t fESigTOFMax;
    
    Bool_t fPIDcutITS;
    Bool_t fPIDcutTOF;
    
    //Values and flag for pion PID cuts in TPC
    Double_t fPSigTPCMin;
    Double_t fPSigTPCMax;
    
    Bool_t fPionPIDcutTPC;
    Bool_t fIsIonColl;
    
    Bool_t fIsV0tree;
    Bool_t fHasSDD;
    TH2F* fArmPlot;

    Bool_t fIsAOD;
    Int_t fFilterBit;
    
    //Grid PID
    Bool_t fIsGRIDanalysis;
    Int_t fGridPID;
    
    ClassDef(AliAnalysisTaskSimpleTreeMaker, 2); //

};



#endif

