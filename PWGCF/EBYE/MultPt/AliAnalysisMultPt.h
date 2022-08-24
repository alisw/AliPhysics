#ifndef AliAnalysisMultPt_H
#define AliAnalysisMultPt_H


#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliEventCuts.h"

class TH1;
class THn;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TTree;



class AliAnalysisMultPt : public AliAnalysisTaskSE
{
    public:
                        AliAnalysisMultPt();
                        AliAnalysisMultPt(const char *name);
        virtual         ~AliAnalysisMultPt();

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t* option);
        virtual void    Terminate(Option_t* option);

        
        void            SetMCRead(Bool_t flag) { fIsMC = flag; }
        void            SetPtLimits(Double_t ptmin, Double_t ptmax) { fPtmin = ptmin; fPtmax=ptmax; }
        void            SetEtaMinLimit(Double_t etalowlimit) { fEtaMin = etalowlimit; }
        void            SetEtaMaxLimit(Double_t etauplimit) { fEtaMax = etauplimit; }
        void            SetFilterBit(Int_t filterbit) { fBit = filterbit; }
        void            SetPVzMinLimit(Float_t pvzmin) {fPVzMin=pvzmin;}
        void            SetPVzMaxLimit(Float_t pvzmax) {fPVzMax=pvzmax;}
 
    private:
        void            BuildData();
        void            BuildMC();
        AliAODEvent     *fAOD;             //! input event
        AliMCEvent 	    *fMC;              //! MC input event
        TList           *fOutputList;      //! output list
        
        TH1D            *MultHist;         // Mult histogram
        TH2D            *MultPtHist;       // pT-Mult histogram
        TH2D            *MultPtHistRec;    // pT-Mult MC-Rec histogram
        TH2D            *MultPtHistGen;    // pT-Mult MC-Gen histogram
        TH2D            *MultHistRatio;    // Ratio of Gen mult vs rec mult histogram
        TH2F		    *hV0MVsPtGen;      // V0 Gen
        TH2F		    *hV0MVsPtRec;      // V0 Rec
        TH1D            *fptRatio;         // pt ratio
    
        Bool_t          fIsMC;             // MC flag
        Double_t        fPtmin;            // min PT
        Double_t        fPtmax;            // max PT
        Double_t        fEtaMin;           // min eta
        Double_t        fEtaMax;           // max eta
        Int_t           fBit;              // filterbit
        Float_t         fPVzMax;           // max PVz
        Float_t         fPVzMin;           // min PVz
        Double_t        fChi2DoF;          // limit for chi2
        UShort_t        fTPCNCrossedRows;  // cr
        Bool_t          fIsRunFBOnly;      // only filterbit cuts
        Int_t           fTPCNcls;          // limit for TPC Ncls
        Bool_t          fIsPileUpCuts;     // pile up cuts flag
        Int_t           fPileUpLevel;      // PLevel
        TString         fGenName;          // MC generator name
    
        AliEventCuts    fEventCuts;

        AliAnalysisMultPt(const AliAnalysisMultPt&); // not implemented
        AliAnalysisMultPt& operator=(const AliAnalysisMultPt&); // not implemented
        ClassDef(AliAnalysisMultPt, 2);
};

#endif
