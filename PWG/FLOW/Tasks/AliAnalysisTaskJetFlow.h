/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetFlow_H
#define AliAnalysisTaskJetFlow_H

// root includes
#include <TMath.h>
// aliroot includes
#include <AliAnalysisTaskSE.h>
// forward declarations
class TString;
class TObject;
class TList;
class TArrayD;
class TClonesArray;
class AliFlowTrackCuts;
class AliFlowEventCuts;
class AliFlowEvent;
class TH1;
class AliAnalysisTaskRhoVnModulation;

class AliAnalysisTaskJetFlow : public AliAnalysisTaskSE
{
    public:
        // enumerators
        enum dataType           {kESD, kAOD, kESDMC, kAODMC };  // data type
        // constructors, destructor
                                AliAnalysisTaskJetFlow();
                                AliAnalysisTaskJetFlow(
                                        const char* name,
                                        AliAnalysisTaskRhoVnModulation* rhoTask, 
                                        Bool_t VPart,           // use jets or tracks as pois
                                        Bool_t VZEROEP,         // do vzero ep method
                                        Bool_t GQC2,            // do gapped qc2 method
                                        Bool_t QC2,             // do qc2 method
                                        Bool_t QC4,             // do simple qc4 method FIXME not implemented yet
                                        Bool_t FlowPackageSP,   // call flow package vzero scalar product
                                        Bool_t FlowPackageQC    // call flow package nth order q-cumulants
                                        );
        virtual                 ~AliAnalysisTaskJetFlow();
        // virtual methods
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        // setters
        void                    SetJetRadius(Double_t r)                {fJetRadius     = r;}
        void                    SetLocalRhoName(TString n)              {fLocalRhoName  = n;}
        void                    SetDebugMode(Int_t d)                   {fDebug         = d;}
        void                    SetCCMinPt(Float_t m)                   {fCCMinPt       = m;}
        void                    SetCCMaxPt(Float_t m)                   {fCCMaxPt       = m;}
        void                    SetCCBinsInPt(Int_t b)                  {fCCBinsInPt    = b;}
        void                    SetMinMaxCentrality(Float_t min, Float_t max)   {fCentralityMin = min; fCentralityMax = max; }
        void                    SetMinimizeDiffBins(Bool_t b)           {fMinimizeDiffBins = b; }
        void                    SetPtBins(TArrayD* pt)                  {fPtBins = pt; }
        void                    SetDoMultWeight(Bool_t m)               {fDoMultWeight = m; }
        void                    SetDoPtWeight(Bool_t p)                 {fDoPtWeight = p; }

        AliFlowTrackCuts*       GetRPCuts() const                       {return fCutsRP_VZERO;}
        AliAnalysisTaskRhoVnModulation* GetMaster() const               {return fRhoVn;}
 

        // cuts
        Bool_t                  PassesCuts();
        // analysis details
        void                    DoVZEROFlowAnalysis();
        void                    DoGappedQC2Analysis();
        void                    DoQC2FlowAnalysis();
        void                    DoQC4FlowAnalysis();
        Bool_t                  DoFlowPackageFlowAnalysis();
        // q-cumulant helper calculations TODO move to AliAnlaysisTaskRhoVnModulation for consistency
        void                    QCnDifferentialFlowVectors(Double_t* repn, Double_t* impn, Double_t *mp, Double_t *reqn, Double_t *imqn, Double_t* mq, Int_t n);

    private:

        // analysis flags and task setup specifics
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        TString                 fJetsName;              // name of jet list
        Double_t                fJetRadius;             // jet radius
        TString                 fTracksName;            // name of track list
        TString                 fLocalRhoName;          // name of local rho
        TClonesArray*           fPois;                  //! array with pois
        TClonesArray*           fRPs;                   //! array with rps
        TObject*                fLocalRho;              //! local energy density
        TList*                  fOutputList;            //! output list
        dataType                fDataType;              //! data type
        Bool_t                  fVParticleAnalysis;     // do the analysis on vparticles instead of jets
        Bool_t                  fMinimizeDiffBins;      // minimize variables (for low statistics)
        Bool_t                  fDoVZEROFlowAnalysis;   // do vzero flow analysis
        Bool_t                  fDoGappedQC2Analysis;   // do gapped qc2 analysis
        Bool_t                  fDoQC2FlowAnalysis;     // do qc2 flow analysis
        Bool_t                  fDoQC4FlowAnalysis;     // do qc4 flow analysis
        Bool_t                  fDoQCFPAnalysis;        // do qc fp analysis
        Bool_t                  fDoSPFPAnalysis;        // do sp fp analyis
        Bool_t                  fDoMultWeight;          // weight events with multiplicity
        Bool_t                  fDoPtWeight;            // introduce pt weighting for rp's and poi's
        Bool_t                  fInitialized;           //! check if the analysis is initialized
        // members
        Bool_t                  fUsePtWeight;           // use pt weights for the qc analysis
        Float_t                 fCCMinPt;               // min pt for flow analysis(common constants)
        Float_t                 fCCMaxPt;               // max pt for flow analysis (common constants)
        Int_t                   fCCBinsInPt;            // bins in pt for flow analysis (common constants)
        Float_t                 fCentralityMin;         // minimium centrality
        Float_t                 fCentralityMax;         // maximum centrality
        TArrayD*                fPtBins;                // custom pt bins for flow analysis
        // cut objects
        AliFlowTrackCuts*       fCutsRP_VZERO;          //! rp cuts for fzero
        AliFlowTrackCuts*       fCutsNull;              //! empty cuts
        AliFlowEventCuts*       fCutsEvent;             //! event cuts
        // containers, setup
        AliFlowEvent*           fFlowEvent_TPC;         //! container for flow analysis
        AliFlowEvent*           fFlowEvent_VZERO;       //! container for flow analysis
        AliAnalysisTaskRhoVnModulation* fRhoVn;         // common cuts and settings master object, see class header
        // histograms
        TH1F*                   fHistAnalysisSummary;   //! analysis summary
        TH1F*                   fCentralitySelection;   //! centrality selection
        // for event plane flow analysis
        TH1F*                   fVZEROAEP;              //! VZEROA EP
        TH1F*                   fVZEROCEP;              //! VZEROC EP
        TProfile*               fv2VZEROA;              //! v2 from VZEROA
        TProfile*               fv2VZEROC;              //! v2 from VZEROC
        // for qc flow analysis
        TProfile*               fRefCumulants;          //! (weighted) reference cumulant
        TProfile*               fDiffCumlantsV2;        //! (weighted) differential cumulant
        TProfile*               fDiffCumlantsV3;        //! (weighted) differential cumulant
        TH1F*                   fQC2v2;                 //! final qc2 result
        TH1F*                   fQC2v3;                 //! final qc2 result
        // additional histograms
        TProfile*               fTempA;                 //! internal bookkeeping
        TProfile*               fTempC;                 //! internal bookkeeping

        AliAnalysisTaskJetFlow(const AliAnalysisTaskJetFlow&);                  // not implemented
        AliAnalysisTaskJetFlow& operator=(const AliAnalysisTaskJetFlow&);       // not implemented

        ClassDef(AliAnalysisTaskJetFlow, 4);
};

#endif
