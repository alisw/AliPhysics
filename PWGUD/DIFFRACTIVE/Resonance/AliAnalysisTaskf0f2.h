/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

 // Short comment describing what this class does needed!

 //===========================================================
 // Dummy comment, should be replaced by a real one
 //===========================================================

#ifndef ALIANALYSISTASKF0F2_H
#define ALIANALYSISTASKF0F2_H

#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDResponse.h"

class AliAnalysisTaskf0f2RunTable {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
        AliAnalysisTaskf0f2RunTable();
        AliAnalysisTaskf0f2RunTable(Int_t runnumber);
        ~AliAnalysisTaskf0f2RunTable();

        Bool_t IsPP(){
            return fCollisionType==kPP; 
        }
        Bool_t IsPA(){
            return fCollisionType==kPA; 
        }
        Bool_t IsAA(){
            return fCollisionType==kAA; 
        }
    private:
        Int_t  fCollisionType; //! Is proton-proton collisions?
};

class AliAnalysisTaskf0f2 : public AliAnalysisTaskSE {
    public:

        enum {  kPion=0, kKaon, kProton, kElectron, kUnknown};
        enum {  kSD=0, kDD, kND, kCD, kAllProc};
        //PN = unlike sign, PP and NN are like signs
        enum {  kPN=0, kPP, kNN, kAllSign}; //P=Positive charge, N=Negative

        AliAnalysisTaskf0f2();
        AliAnalysisTaskf0f2
        ( 
              const char *name
            , const char *option
        );

        AliAnalysisTaskf0f2
        (
              const AliAnalysisTaskf0f2& ap
        );

        AliAnalysisTaskf0f2& operator = 
        (
              const AliAnalysisTaskf0f2& ap
        );

        ~AliAnalysisTaskf0f2();
        

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);

        void SetOption(char * option) {fOption = option;}
        Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk); 
    
        Bool_t  GoodTracksSelection();
        void FillTracks();
        Bool_t FindPtBin(Double_t pt, UInt_t & bin);
        Bool_t FindCentBin(Double_t cent, UInt_t & bin);
    private:
        
        TString                         fOption;
        TList*                          fOutput; //!

        AliTriggerAnalysis*             fTrigger; //!
        AliESDtrackCuts*                fTrackCuts; //!
        AliESDEvent*                    fEsd; //!
        AliAODEvent*                    fAod; //!
        Bool_t                          IsFirstEvent;
        AliAnalysisTaskf0f2RunTable*             fRunTable; //!
        
        Double_t                        fCent; 
        std::vector < UInt_t >          goodtrackindices; //!
        
        AliPIDResponse                 *fPIDResponse; //!
        TH1D                           *fEventNumbers; //!
        TH2D                           *fPhiEta; //!

        //Histograms below are main
        std::vector< std::vector< TH2D* > > fMass2D; //! signbins, centbins
    

    ClassDef(AliAnalysisTaskf0f2, 1);
};

#endif

