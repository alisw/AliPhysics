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
*          Elisa Meninno    (elisa.meninno@cern.ch)    *
*                                                      *
*******************************************************/


class AliAnalysisTaskSimpleTreeMaker : public AliAnalysisTaskSE {

  public:
    AliAnalysisTaskSimpleTreeMaker(const char *name, Bool_t ExtraDCA = kFALSE);
    AliAnalysisTaskSimpleTreeMaker();
    ~AliAnalysisTaskSimpleTreeMaker();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   FinishTaskOutput();
    virtual void   Terminate(Option_t *);

    void SetupTrackCuts(AliDielectronCutGroup* finalTrackCuts);
    void SetupEventCuts(AliDielectronEventCuts* finalEventCuts);
  
    // PID calibration function to correct the width and mean of detector
    // response (I.e should be unit guassian)
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

    // Kinematic cuts set here only applied to V0 TTrees
    // For standard TTree a dielectron cut library is needed
    void SetPtRange(Float_t min, Float_t max){
      fPtMin = min;
      fPtMax = max;
    }

    void SetEtaRange(Float_t min, Float_t max){
      fEtaMin = min;
      fEtaMax = max;
    }

    // PID cut used for V0 TTrees
    void SetESigRangeTPC(Float_t min, Float_t max){
      fESigTPCMin = min;
      fESigTPCMax = max;
    }

    // Kaon PID values not saved by default
    void writeKaonPIDtoTree(Bool_t answer){
      storeKaonPID = answer;
    }
    // Proton PID values not saved by default
    void writeProtonPIDtoTree(Bool_t answer){
      storeProtonPID = answer;
    }

    void SetMC(Bool_t answer){ hasMC = answer; }

    // Track cut setters. StandardITSTPC2011 cuts used if nothing specified
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
    void createDCAbranches(Bool_t answer){
      fExtraDCA = answer;
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

    void SetMaxPtPIDcorrection(Float_t answer){
      maxPtPIDcorrection = answer;
    }

    inline Bool_t GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0){
      // this is a copy of the AliDielectronVarManager

      if(track->TestBit(AliAODTrack::kIsDCA)){
        d0z0[0] = track->DCA();
        d0z0[1] = track->ZAtDCA();
        // the covariance matrix is not stored in case of AliAODTrack::kIsDCA
        return kTRUE;
      }

      Bool_t ok = kFALSE;
      if(event){
        AliExternalTrackParam etp;
        etp.CopyFromVTrack(track);

        Float_t xstart = etp.GetX();
        if(xstart>3.){
          d0z0[0] = -999.;
          d0z0[1] = -999.;
          return kFALSE;
        }

        const AliAODVertex* vtx =dynamic_cast<const AliAODVertex*>((event->GetPrimaryVertex()));
        Double_t fBzkG = event->GetMagneticField(); // z componenent of field in kG
        ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
      }

      if(!ok){
        d0z0[0] = -999.;
        d0z0[1] = -999.;
      }
      return ok;
    }

    // Check if the generator is on the list of generators
    // If found, assign track with integer value correspding to generator
    // 0 = gen purp, 1=Pythia CC_1, 2= Pythia BB_1, 3=Pythia B_1, 4=Jpsi2ee_1, 5=B2Jpsi2ee_1";
    Int_t CheckGenerator(Int_t trackID);

  private:

    Int_t IsEventAccepted(AliVEvent* event);

    AliAnalysisTaskSimpleTreeMaker(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

    AliAnalysisTaskSimpleTreeMaker& operator=(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

    Int_t eventNum; //!
    Bool_t hasMC;

    AliESDtrackCuts* fESDtrackCuts;

    AliPIDResponse* fPIDResponse; //! PID response object
    AliMCEvent* fMCevent; //!

    TTree* fTree;

    // Dielectron cut classes needed to source cuts from LMEE cut libraries
    // The desired cut library should be specified in the AddTask
    AliDielectronEventCuts* eventCuts;
    AliAnalysisFilter* eventFilter;

    AliDielectronVarCuts* varCuts;
    AliDielectronTrackCuts *trackCuts;
    AliDielectronPID *pidCuts;
    AliDielectronCutGroup* cuts;
    AliAnalysisFilter* trackFilter;

    // Class needed to use PID within the Dielectron Framework
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
    //DCA
    Float_t DCA[2];
    Float_t DCAsigma[2];
    Float_t DCA3[3];
    //Double_t DCA3sigma[3];
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
    Float_t maxPtPIDcorrection;
    Float_t PnSigmaITS;
    Float_t PnSigmaTPC;
    Float_t PnSigmaTOF;
    Bool_t storeKaonPID;
    Float_t KnSigmaITS;
    Float_t KnSigmaTPC;
    Float_t KnSigmaTOF;
    Bool_t storeProtonPID;
    Float_t PrnSigmaITS;
    Float_t PrnSigmaTPC;
    Float_t PrnSigmaTOF;
    Float_t ITSsignal;
    Float_t TPCsignal;
    Float_t TOFsignal;
    Float_t goldenChi2;
    // MC
    Float_t mcEta;
    Float_t mcPhi;
    Float_t mcPt;
    Float_t mcVert[3];
    Int_t iPdg;
    Int_t iPdgMother;
    Bool_t HasMother;
    Int_t motherLabel;
    Int_t isInj;
    Bool_t isPhysPrimary;
    // Pdg and label for initial particle in decay chain
    Int_t iPdgFirstMother;
    Int_t gLabelFirstMother;
    Int_t gLabelMinFirstMother;
    Int_t gLabelMaxFirstMother;
    // V0 features
    Float_t pointingAngle;
    Float_t daughtersDCA;
    Float_t decayLength;
    Float_t v0mass;
    Float_t ptArm;
    Float_t alpha;

    TH1F* fQAhist; //!
    // Currently no cut on centrality
    Float_t fCentralityPercentileMin;
    Float_t fCentralityPercentileMax;

    // Kinematic cuts used in the creation of the V0 TTrees
    Float_t fPtMin;
    Float_t fPtMax;
    Float_t fEtaMin;
    Float_t fEtaMax;

    // TPC PID cuts for V0 TTree
    Float_t fESigTPCMin;
    Float_t fESigTPCMax;

    Bool_t fHasSDD;
    Bool_t fExtraDCA;
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

    // Store list of generator hashes which can be checked to determine whether
    // or not the track was injected
    std::vector<UInt_t> fGeneratorHashes;

    ClassDef(AliAnalysisTaskSimpleTreeMaker, 7);

};



#endif

