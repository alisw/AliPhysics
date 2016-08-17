/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

 // Short comment describing what this class does needed!

 //===========================================================
 // Dummy comment, should be replaced by a real one
 //===========================================================

#ifndef ALIRSNF0F2TASK_H
#define ALIRSNF0F2TASK_H

#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDResponse.h"

class AliRsnf0f2RunTable {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
        AliRsnf0f2RunTable();
        AliRsnf0f2RunTable(Int_t runnumber);
        ~AliRsnf0f2RunTable();

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

class AliRsnf0f2Task : public AliAnalysisTaskSE {
    public:

        enum {  kPion=0, kKaon, kProton, kElectron, kUnknown};
        enum {  kSD=0, kDD, kND, kCD, kAllProc};
        //PN = unlike sign, PP and NN are like signs
        enum {  kPN=0, kPP, kNN, kAllSign}; //P=Positive charge, N=Negative

        AliRsnf0f2Task();
        AliRsnf0f2Task
        ( 
              const char *name
            , const char *option
        );

        AliRsnf0f2Task
        (
              const AliRsnf0f2Task& ap
        );

        AliRsnf0f2Task& operator = 
        (
              const AliRsnf0f2Task& ap
        );

        ~AliRsnf0f2Task();
        

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
        AliRsnf0f2RunTable*             fRunTable; //!
        
        Double_t                        fCent; 
        std::vector < UInt_t >          goodtrackindices; //!
        
        AliPIDResponse                 *fPIDResponse; //!
        TH1D                           *fEventNumbers; //!
        TH2D                           *fPhiEta; //!

        //Histograms below are main
        std::vector< std::vector< TH2D* > > fMass2D; //! signbins, centbins
    

    ClassDef(AliRsnf0f2Task, 1);
};

#endif

