#ifndef AliAnalysisTaskPi0Flow_cxx
#define AliAnalysisTaskPi0Flow_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task to fill histograms with PHOS ESD or AOD clusters and cells
// Authors : Dmitri Peressounko
// Date    : 28.05.2011
// Modified: 03.08.2012 Henrik Qvigstad
/* $Id$ */

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliESDEvent ;
class AliPHOSCalibData;
class AliESDtrack ;
class AliESDCaloCluster ;
class AliEPFlattener;

#include "TArrayD.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0Flow : public AliAnalysisTaskSE {
public:
    enum Period { kUndefinedPeriod, kLHC10h, kLHC11h };

public:
    AliAnalysisTaskPi0Flow(const char *name = "AliAnalysisTaskPi0Flow", Period period = kUndefinedPeriod);
    virtual ~AliAnalysisTaskPi0Flow();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    /* virtual void   Terminate(Option_t *); */

    void SetPeriod(Period period) { fPeriod = period;}
    
    void SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed);
    void SetEventMixingRPBinning(UInt_t nBins) { fNEMRPBins = nBins; }
    void SetMaxAbsVertexZ(Float_t z) { fMaxAbsVertexZ = z; }
    void SetManualV0EPCalc(Bool_t manCalc = true) {fManualV0EPCalc = manCalc;}
    
    void SetPHOSBadMap(Int_t mod,TH2I * badMapHist);
    //Where to read AODB object with EP calibration if not default
    void SetEPcalibFileName(const TString filename) {fEPcalibFileName = filename; }   

private:
    AliAnalysisTaskPi0Flow(const AliAnalysisTaskPi0Flow&); // not implemented
    AliAnalysisTaskPi0Flow& operator=(const AliAnalysisTaskPi0Flow&); // not implemented

    // Step 0:
    AliVEvent* GetEvent();
    AliStack* GetMCStack();

    // Step 1:
    void SetGeometry();
    void SetMisalignment();
    void SetV0Calibration(); //V0 calibration
    void SetESDTrackCuts(); // AliESDtrack cuts ( for esd data )
    void SetPHOSCalibData(); // phos re-calibration ( for esd data)
    void SetFlatteningData(); // phos flattening

    // Step 2:
    void SetVertex();
    Bool_t RejectEventVertex();

    // Step 4:
    void SetCentrality();
    Bool_t RejectCentrality();

    // Step 5:
    void EvalReactionPlane();
    void EvalV0ReactionPlane();

    // Step 7: QA PHOS cells
    void FillPHOSCellQAHists();

    // Step 8: Event Photons (PHOS Clusters) selection
    void SelectPhotonClusters();
    void FillSelectedClusterHistograms();

    // Step 9: Consider pi0 (photon/cluster) pairs.
    void ConsiderPi0s();

    // Step 10; Mixing
    void ConsiderPi0sMix();

    // Step 11: Update lists
    void UpdateLists();

    Bool_t AreNeibors(Int_t id1,Int_t id2) ;
    Double_t ApplyFlattening(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
    Double_t ApplyFlatteningV0A(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
    Double_t ApplyFlatteningV0C(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
    Int_t ConvertToInternalRunNumber(Int_t run) ;
    Double_t CoreEnergy(AliVCluster * clu, AliVCaloCells * cells);
    void EvalCoreLambdas(AliVCluster * clu, AliVCaloCells * cells, Double_t &m02, Double_t &m20) ; 
    Bool_t TestCoreLambda(Double_t pt,Double_t l1,Double_t l2) ;




    void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key

    TVector3 GetVertexVector(const AliVVertex* vertex);
    Int_t GetCentralityBin(Float_t centralityV0M);
    Int_t GetRPBin();

    void LogProgress(int step);
    void LogSelection(int step, int internalRunNumber);

    Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS


    void Reclusterize(AliVCluster * clu) ;
    Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
    Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
    Bool_t TestLambda2(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
    
    UInt_t GetNumberOfCentralityBins() { return fCentEdges.GetSize()-1; }
    TList* GetCaloPhotonsPHOSList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin);
    



private:
    // constants:
    static const Double_t kLogWeight= 4.5 ; // log weight for recalibration.
    static const Double_t kAlphaCut=0.7 ;
    static const Bool_t doESDReCalibration = kTRUE;
    static const Int_t kNCenBins = 9; // see EvalV0ReactionPlane()

    // cluster cut variables:
    static const Double_t kMinClusterEnergy = 0.3;
    static const Double_t kMinBCDistance = 2.5;  //distance to nearest bad channel
    static const Int_t kMinNCells = 3;
    static const Double_t kMinM02 = 0.2;

    // Binning, [vtx, centrality, reaction-plane]
    static const Int_t kNVtxZBins = 1;
    static const Double_t kCentCutoff = 90.; // Ignore Centrality over 90%
    TArrayD fCentEdges;  // Centrality Bin Lower edges
    TArrayI fCentNMixed; // Number of mixed events for each centrality bin
    UInt_t fNEMRPBins;
    

    // Behavior / cuts
    Period fPeriod;
    Float_t fMaxAbsVertexZ; // in cm
    Bool_t fManualV0EPCalc;


    TList * fOutputContainer;        //final histogram container

    TF1 *fNonLinCorr;          // Non-linearity correction
//TF1 * fRecent[5][12] ;//Recentering corrections
    TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

    // Run variables

    // Step 0: Event Objects
    // fEvent, fMCStack
    AliVEvent* fEvent; //! Current event
    AliESDEvent* fEventESD; //! Current event, if ESD.
    AliAODEvent* fEventAOD; //! Current event, if AOD.
    AliStack * fMCStack ;

    // Step 1: Run Number, Misalignment Matrix, and Calibration
    Int_t fRunNumber; // run number
    Int_t fInternalRunNumber ;    //Current internal run number
    AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
    TProfile *fMultV0;                  // object containing VZERO calibration information
    Float_t fV0Cpol,fV0Apol;            // loaded by OADB
    Float_t fMeanQ[kNCenBins][2][2];    // and recentering
    Float_t fWidthQ[kNCenBins][2][2];   // ...
    AliESDtrackCuts *fESDtrackCuts; // Track cut
    AliPHOSCalibData *fPHOSCalibData; // PHOS calibration object
    TString fEPcalibFileName; 
    AliEPFlattener * fTPCFlat ; //Object for flattening of TPC
    AliEPFlattener * fV0AFlat ; //Object for flattening of V0A
    AliEPFlattener * fV0CFlat ; //Object for flattening of V0C
    
    
    // Step 2: Vertex
    Double_t fVertex[3];
    TVector3 fVertexVector;
    Int_t fVtxBin;

    // Step 4: Centrality
    Float_t fCentralityV0M ; //!Centrality of the currecnt event
    Int_t fCentBin ;       //! Current centrality bin

    // Step 5: Reaction Plane
    Bool_t fHaveTPCRP ; //! Is TPC RP defined?
    Float_t fRP ;       //!Reaction plane calculated with full TPC
    Float_t fRPV0A ;    //!Reaction plain calculated with A-side TPC: eta>0.15
    Float_t fRPV0C ;    //!Reaction plain calculated with C-side TPC: eta<-0.15
    Int_t fEMRPBin;       //! Event Mixing Reaction Plane Bin

    // Step 8: Event Photons (PHOS Clusters) selection
    TObjArray * fCaloPhotonsPHOS ;      //PHOS photons in current event

    // Step 11: Update lists for mixing.
    TObjArray* fCaloPhotonsPHOSLists; //! array of TList, Containers for events with PHOS photons


    ClassDef(AliAnalysisTaskPi0Flow, 1); // PHOS analysis task
};

#endif
