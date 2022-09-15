#ifndef AliAnalysisMultPt_H
#define AliAnalysisMultPt_H


class TH1;
class THn;
class TH1F;
class TH2D;
class TH2F;
class TH3D;
class TList;
class TTree;

#include "AliAnalysisTaskSE.h"

#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliEventCuts.h"


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
    
        AliAODEvent     *fAOD;             // input event
        AliMCEvent      *fMC;              // MC input event
        TList           *fOutputList;      // output list
        
        TH1D            *MultHist;         // Mult histogram
        TH2D            *MultPtHist;       // pT-Mult histogram
        TH2D            *MultPtHistRec;    // pT-Mult MC-Rec histogram
        TH2D            *MultPtHistGen;    // pT-Mult MC-Gen histogram
        TH2D            *MultHistRatio;    // Ratio of Gen mult vs rec mult histogram
        TH2D            *hV0MVsPtGen;      // V0 Gen
        TH2D            *hV0MVsPtRec;      // V0 Rec
        
    
        Bool_t          fIsMC;             // MC flag
        Double_t        fPtmin;            // min PT
        Double_t        fPtmax;            // max PT
        Double_t        fEtaMin;           // min eta
        Double_t        fEtaMax;           // max eta
        Int_t           fBit;              // filterbit
        Float_t         fPVzMax;           // max PVz
        Float_t         fPVzMin;           // min PVz
    
        AliEventCuts    fEventCuts;

        AliAnalysisMultPt(const AliAnalysisMultPt&); // not implemented
        AliAnalysisMultPt& operator=(const AliAnalysisMultPt&); // not implemented
        ClassDef(AliAnalysisMultPt, 2);
};

#endif
