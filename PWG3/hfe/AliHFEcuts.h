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
//
// Cut container class for the ALICE HFE group
// Serves also as interface to the correction Framework
// Provides a set of standard cuts
//
#ifndef ALIHFECUTS_H
#define ALIHFECUTS_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ALIHFEEXTRACUTS_H
#include "AliHFEextraCuts.h"
#endif

class AliCFManager;
class AliESDtrack;
class AliMCParticle;

class TObjArray;
class TList;

class AliHFEcuts : public TNamed{
  public:
    typedef enum{
      kStepRecNoCut = 0,
      kStepRecKineITSTPC = 1,
      kStepRecPrim = 2,
      kStepHFEcutsITS = 3,
      kStepHFEcutsTOF = 4,
      kStepHFEcutsTPC = 5,
      kStepHFEcutsTRD = 6,
      kNcutStepsRecTrack = 7
    } RecoCutStep_t;
    typedef enum{
      kStepHFEcutsDca = 0, 
      kNcutStepsDETrack = 1
    } DECutStep_t;
    typedef enum{
      kStepHFEcutsSecvtx = 0, 
      kNcutStepsSecvtxTrack = 1
    } SecvtxCutStep_t;
    typedef enum{
      kStepMCGenerated = 0,
      kStepMCGeneratedZOutNoPileUpCentralityFine = 1,
      kStepMCGeneratedEventCut = 2,
      kStepMCInAcceptance = 3,
      kNcutStepsMCTrack =  4
    } MCCutStep_t;
    typedef enum{
      kEventStepGenerated = 0,
      kEventStepRecNoCut = 1,
      kEventStepRecNoPileUp = 2,
      kEventStepRecCentralityOk = 3,
      kEventStepZRange = 4,
      kEventStepReconstructed = 5,
      kNcutStepsEvent = 6
    } EventCutStep_t;

    AliHFEcuts();
    AliHFEcuts(const Char_t *name, const Char_t *title);
    AliHFEcuts(const AliHFEcuts &c);
    AliHFEcuts &operator=(const AliHFEcuts &c);
    void Copy(TObject &o) const;
    Long64_t Merge(const TCollection *list);
    ~AliHFEcuts();
    
    void Initialize(AliCFManager *cfm);
    void Initialize();

    Bool_t CheckParticleCuts(UInt_t step, TObject *o);
    Bool_t CheckEventCuts(const char*namestep, TObject *o);
  
    TList *GetQAhistograms() const { return fHistQA; }
    
    void SetQAOn() {SetBit(kDebugMode, kTRUE); };
    void UnsetQA() {SetBit(kDebugMode, kFALSE); };
    Bool_t IsQAOn() const { return TestBit(kDebugMode); };
    void SetAOD() { SetBit(kAOD, kTRUE); }
    void SetESD() { SetBit(kAOD, kFALSE); }
    Bool_t IsAOD() const { return TestBit(kAOD); }
    Bool_t IsESD() const { return !TestBit(kAOD); }

    // Cut Names
    static const Char_t *MCCutName(UInt_t step){
      if(step >= kNcutStepsMCTrack) return fgkUndefined;
      return fgkMCCutName[step];
    };
    static const Char_t *RecoCutName(UInt_t step){
      if(step >= kNcutStepsRecTrack) return fgkUndefined;
      return fgkRecoCutName[step];
    }
    static const Char_t *DECutName(UInt_t step){
      if(step >= kNcutStepsDETrack) return fgkUndefined;
      return fgkDECutName[step];
    }
    static const Char_t *SecvtxCutName(UInt_t step){
      if(step >= kNcutStepsSecvtxTrack) return fgkUndefined;
      return fgkSecvtxCutName[step];
    }
    static const Char_t *EventCutName(UInt_t step){
      if(step >= kNcutStepsEvent) return fgkUndefined;
      return fgkEventCutName[step];
    }
   
    // Getters
    Bool_t IsRequireITSpixel() const { return TESTBIT(fRequirements, kITSPixel); };
    Bool_t IsRequireMaxImpactParam() const { return TESTBIT(fRequirements, kMaxImpactParam); };
    Bool_t IsRequirePrimary() const { return TESTBIT(fRequirements, kPrimary); };
    Bool_t IsRequireProdVertex() const { return TESTBIT(fRequirements, kProductionVertex); };
    Bool_t IsRequireSigmaToVertex() const { return TESTBIT(fRequirements, kSigmaToVertex); };
    Bool_t IsRequireDCAToVertex() const {return TESTBIT(fRequirements, kDCAToVertex); };
    Bool_t IsRequireKineMCCuts() const {return TESTBIT(fRequirements, kKineMCCuts); };
    Double_t GetVertexRange() const {return fVertexRangeZ; };
    Int_t GetMinTrackletsTRD() const { return fMinTrackletsTRD; }
    
    // Setters
    inline void SetCutITSpixel(UChar_t cut);
    void SetCheckITSLayerStatus(Bool_t checkITSLayerStatus) { fCheckITSLayerStatus = checkITSLayerStatus; }
    void SetMinNClustersTPC(UChar_t minClustersTPC) { fMinClustersTPC = minClustersTPC; }
    void SetMinNClustersITS(UChar_t minClustersITS) { fMinClustersITS = minClustersITS; }
    void SetMinNTrackletsTRD(UChar_t minNtrackletsTRD) { fMinTrackletsTRD = minNtrackletsTRD; }
    void SetMaxChi2perClusterITS(Double_t chi2) { fMaxChi2clusterITS = chi2; };
    void SetMaxChi2perClusterTPC(Double_t chi2) { fMaxChi2clusterTPC = chi2; };
    inline void SetMaxImpactParam(Double_t radial, Double_t z);
    inline void SetIPcutParam(Float_t p0, Float_t p1, Float_t p2, Float_t p3, Bool_t isipsigma);
    void SetMinRatioTPCclusters(Double_t minRatioTPC) { fMinClusterRatioTPC = minRatioTPC; };
    void SetPtRange(Double_t ptmin, Double_t ptmax){fPtRange[0] = ptmin; fPtRange[1] = ptmax;};
    inline void SetProductionVertex(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax);
    inline void SetSigmaToVertex(Double_t sig);
    inline void SetSigmaToVertexXY(Double_t sig);
    inline void SetSigmaToVertexZ(Double_t sig);
    void SetTPCmodes(AliHFEextraCuts::ETPCclusterDef_t clusterDef, AliHFEextraCuts::ETPCclrDef_t ratioDef) {
      fTPCclusterDef= clusterDef;
      fTPCratioDef = ratioDef;
    }
    void SetVertexRange(Double_t zrange){fVertexRangeZ = zrange;};
    void SetTOFPIDStep(Bool_t tofPidStep) {fTOFPIDStep = tofPidStep;};
    void SetTOFMISMATCHStep(Bool_t tofMismatchStep) {fTOFMISMATCHStep = tofMismatchStep;};
    void SetTPCPIDCleanUpStep(Bool_t tpcPIDCleanUpStep) {fTPCPIDCLEANUPStep = tpcPIDCleanUpStep;};
    void SetUseMixedVertex(Bool_t useMixedVertex) {fUseMixedVertex = useMixedVertex;};    
    void SetFractionOfSharedTPCClusters( Bool_t fractionOfSharedTPCClusters) {fFractionOfSharedTPCClusters = fractionOfSharedTPCClusters;};
    void SetMaxImpactParameterRpar(Bool_t maxImpactParameterRpar) { fMaxImpactParameterRpar = maxImpactParameterRpar; };
    
    inline void CreateStandardCuts();
    
    // Requirements
    void SetRequireDCAToVertex() { SETBIT(fRequirements, kDCAToVertex); CLRBIT(fRequirements, kSigmaToVertex); };
    void SetRequireIsPrimary() { SETBIT(fRequirements, kPrimary); };
    void SetRequireITSPixel() { SETBIT(fRequirements, kITSPixel); }
    void SetRequireProdVertex() { SETBIT(fRequirements, kProductionVertex); };
    void SetRequireSigmaToVertex() { SETBIT(fRequirements, kSigmaToVertex); CLRBIT(fRequirements, kDCAToVertex); };
    void UnsetVertexRequirement() { CLRBIT(fRequirements, kDCAToVertex); CLRBIT(fRequirements, kSigmaToVertex); }
    void SetRequireKineMCCuts() { SETBIT(fRequirements, kKineMCCuts); };

    void SetDebugLevel(Int_t level) { fDebugLevel = level; };
    Int_t GetDebugLevel() const { return fDebugLevel; };

  private:
    enum{
      kDebugMode = BIT(14),
      kAOD = BIT(15)
    };
    typedef enum{
      kPrimary = 0,
      kProductionVertex = 1,
      kSigmaToVertex = 2,
      kDCAToVertex = 3,
      kITSPixel = 4,
      kMaxImpactParam = 5,
      kKineMCCuts = 6
    } Require_t;
    void SetParticleGenCutList();
    void SetAcceptanceCutList();
    void SetRecKineITSTPCCutList();
    void SetRecPrimaryCutList();
    void SetHFElectronITSCuts();
    void SetHFElectronTOFCuts();
    void SetHFElectronTPCCuts();
    void SetHFElectronTRDCuts();
    void SetHFElectronDcaCuts();
    void SetEventCutList(Int_t istep);

    static const Char_t* fgkMCCutName[kNcutStepsMCTrack];     // Cut step names for MC single Track cuts
    static const Char_t* fgkRecoCutName[kNcutStepsRecTrack];  // Cut step names for Rec single Track cuts
    static const Char_t* fgkDECutName[kNcutStepsDETrack];     // Cut step names for impact parameter cuts
    static const Char_t* fgkSecvtxCutName[kNcutStepsSecvtxTrack];     // Cut step names for secondary vertexing cuts
    static const Char_t* fgkEventCutName[kNcutStepsEvent];    // Cut step names for Event cuts
    static const Char_t* fgkUndefined;                        // Name for undefined (overflow)
  
    ULong64_t fRequirements;  	  // Bitmap for requirements
    AliHFEextraCuts::ETPCclusterDef_t fTPCclusterDef;       // TPC cluster definition
    AliHFEextraCuts::ETPCclrDef_t fTPCratioDef;             // TPC cluster ratio Definition
    Double_t fDCAtoVtx[2];	      // DCA to Vertex
    Double_t fProdVtx[4];	        // Production Vertex
    Double_t fPtRange[2];	        // pt range
    UChar_t fMinClustersTPC;	    // Min.Number of TPC clusters
    UChar_t fMinClustersITS;	    // Min.Number of TPC clusters
    UChar_t fMinTrackletsTRD;	    // Min. Number of TRD tracklets
    UChar_t fCutITSPixel;	        // Cut on ITS pixel
    Bool_t  fCheckITSLayerStatus;       // Check ITS layer status
    Double_t fMaxChi2clusterITS;	// Max Chi2 per ITS cluster
    Double_t fMaxChi2clusterTPC;	// Max Chi2 per TPC cluster
    Double_t fMinClusterRatioTPC;	// Min. Ratio findable / found TPC clusters
    Double_t fSigmaToVtx[3];	    // Sigma To Vertex
    Double_t fVertexRangeZ;             // Vertex Range reconstructed
    Bool_t   fTOFPIDStep;               // TOF matching step efficiency
    Bool_t   fTOFMISMATCHStep;        // TOF mismatch step
    Bool_t   fTPCPIDCLEANUPStep;      // TPC PIC cleanup step
    Bool_t   fUseMixedVertex;         // Use primary vertex from track only as before
    Float_t  fIPCutParams[4];         // Parameters of impact parameter cut parametrization
    Bool_t   fIsIPSigmacut;           // if abs IP cut or IP sigma cut 
    Double_t fFractionOfSharedTPCClusters; // Fraction of shared TPC clusters
    Bool_t   fMaxImpactParameterRpar;      // Max impact parameter


    
    TList *fHistQA;		            //! QA Histograms
    TObjArray *fCutList;	        //! List of cut objects(Correction Framework Manager)

    Int_t fDebugLevel;            // Debug Level
    
  ClassDef(AliHFEcuts, 2)         // Container for HFE cuts
};

//__________________________________________________________________
void AliHFEcuts::SetProductionVertex(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax){
  // Set the production vertex constraint
  SetRequireProdVertex();
  fProdVtx[0] = xmin;
  fProdVtx[1] = xmax;
  fProdVtx[2] = ymin;
  fProdVtx[3] = ymax;
}

//__________________________________________________________________
void AliHFEcuts::SetSigmaToVertex(Double_t sig){
  SetRequireSigmaToVertex();
  fSigmaToVtx[0] = sig;
}

//__________________________________________________________________
void AliHFEcuts::SetSigmaToVertexXY(Double_t sig){
  SetRequireSigmaToVertex();
  fSigmaToVtx[1] = sig;
}

//__________________________________________________________________
void AliHFEcuts::SetSigmaToVertexZ(Double_t sig){
  SetRequireSigmaToVertex();
  fSigmaToVtx[2] = sig;
}

//__________________________________________________________________
void AliHFEcuts::SetMaxImpactParam(Double_t radial, Double_t z){
  SetRequireDCAToVertex();
  fDCAtoVtx[0] = radial;
  fDCAtoVtx[1] = z;
}

//__________________________________________________________________
void AliHFEcuts::SetIPcutParam(Float_t p0, Float_t p1, Float_t p2, Float_t p3, Bool_t isipsigma){
  // Set parameters for impact parameter cut parametrization
  fIPCutParams[0] = p0;
  fIPCutParams[1] = p1;
  fIPCutParams[2] = p2;
  fIPCutParams[3] = p3;
  fIsIPSigmacut = isipsigma;
}

//__________________________________________________________________
void AliHFEcuts::SetCutITSpixel(UChar_t cut){
  SetRequireITSPixel();
  fCutITSPixel = cut;
}

//__________________________________________________________________
void AliHFEcuts::CreateStandardCuts(){
  //
  // Standard Cuts defined by the HFE Group
  //
  SetRequireProdVertex();
  fProdVtx[0] = 0;
  fProdVtx[1] = 3;
  fProdVtx[2] = 0;
  fProdVtx[3] = 3;
  //SetRequireDCAToVertex();
  //fDCAtoVtx[0] = 0.5;
  //fDCAtoVtx[1] = 1.5;
  fMinClustersTPC = 80;
  fMinClustersITS = 4;
  fMinTrackletsTRD = 0;
  SetRequireITSPixel();
  fCutITSPixel = AliHFEextraCuts::kFirst;
  fMaxChi2clusterITS = -1.;
  fMaxChi2clusterTPC = 4.;
  fMinClusterRatioTPC = 0.6;
  fPtRange[0] = 0.1;
  fPtRange[1] = 20.;
  SetRequireKineMCCuts();
}
#endif
