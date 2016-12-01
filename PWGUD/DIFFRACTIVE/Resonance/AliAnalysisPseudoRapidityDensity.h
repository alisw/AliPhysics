/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

 // Short comment describing what this class does needed!

 //===========================================================
 // Dummy comment, should be replaced by a real one
 //===========================================================

#ifndef ALIANALYSISPSEUDORAPIDITYDENSITY_H
#define ALIANALYSISPSEUDORAPIDITYDENSITY_H

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "THistManager.h"
#include <deque>

class AliVMultiplicity;
class TTree;
class AliAnalysisUtils;
class AliOADBMultSelection;
class AliStack;
class TRandom3;

class AliAnalysisPseudoRapidityDensityRunTable {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
        AliAnalysisPseudoRapidityDensityRunTable();
        AliAnalysisPseudoRapidityDensityRunTable(Int_t runnumber);
        ~AliAnalysisPseudoRapidityDensityRunTable();

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

class AliAnalysisPseudoRapidityDensity : public AliAnalysisTaskSE {
    public:
        typedef std::vector<Bool_t> Bool_1d;

        enum {  kECbegin=1,kDATA=1, kINEL, kNSD, kINELg0, kSD, kDD, kND , kECend};
        enum {  kTrigbegin=1, kMBOR=1, kMBAND, kMBORg0 , kTrigend};
        enum {  kParTypebegin=1,kParDATA=1,kPion,kKaon,kProton, kOPar, kParTypeend};
        enum {  kV0Typebegin=1,kK0s, kLambda, kAntiLambda, kV0Typeend};
        //PN = unlike sign, PP and NN are like signs

        AliAnalysisPseudoRapidityDensity();
        AliAnalysisPseudoRapidityDensity
        ( 
              const char *name
            , const char *option
        );

        AliAnalysisPseudoRapidityDensity
        (
              const AliAnalysisPseudoRapidityDensity& ap
        );

        AliAnalysisPseudoRapidityDensity& operator = 
        (
              const AliAnalysisPseudoRapidityDensity& ap
        );

        ~AliAnalysisPseudoRapidityDensity();
        

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);

        void SetOption(char * option) {fOption = option;}
        void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}
        Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk); 
    
        void FillTracklets(Bool_1d bevtc, Bool_1d btrigc);
        void StrangenessMeasure( Bool_1d btrigc);
        void MeasureDiffMass();
        void SetIsAA (Bool_t isaa) {IsAA = isaa;}
        TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
        TAxis AxisVar( TString name, std::vector<Double_t> bin );
        TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax
            , Double_t xmin0);
        TAxis AxisStr( TString name, std::vector<TString> bin );
        THnSparse * CreateTHnSparse(TString name, TString title
            , Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
        THnSparse * CreateTHnSparse(TString name, TString title
            , TString templ, Option_t * opt="");
        Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1.);
        Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1.);


    private:
        typedef std::vector<AliVTrack*> tracklist;
        typedef std::deque<tracklist>  eventpool;
        typedef std::vector<std::vector<eventpool> > mixingpool; 
        typedef std::vector<Int_t> Int_1d;
        
        TString                         fOption;
        TList*                          fOutput=nullptr; //!
       

        AliTriggerAnalysis*             fTrigger=nullptr; //!
        AliESDtrackCuts*                fTrackCuts=nullptr; //!
        AliVEvent*                      fEvt=nullptr; //!
        UInt_t                          fFilterBit;
        Bool_t                          IsFirstEvent=kTRUE;
        AliAnalysisPseudoRapidityDensityRunTable*   fRunTable=nullptr; //!
        
        Double_t                        fCent=-1;
        Double_t                        fZ=-30;
        std::vector < UInt_t >          goodtrackindices; //!
        
        AliPIDResponse                 *fPIDResponse=nullptr; //!
        AliPIDCombined                 *fPIDCombined=nullptr; //!
        AliAnalysisUtils               *fUtils=nullptr;//!
        AliOADBMultSelection           *fOadbMultSelection=nullptr; //!

        //Histograms below are main
        std::vector< std::vector< TH2D* > > fMass2D; //! signbins, centbins
        //Histograms for pT_pair amd pT

        TAxis                           binCent; //! centrality bin 
        TAxis                           binZ; //! bin of zvtx
        TAxis                           binEta; //! bin of pseudorapidity
        TAxis                           binEventClass; //! bin of event classes
        TAxis                           binTriggClass; //! bin of trigg classes
        TAxis                           binParType; //! bin of particle types
        TAxis                           binV0Type; //! bin of v0 particle types
        AliStack*                       stack =nullptr;//!
        Int_t                           ismc = false;
        Int_t                           zmcbin = -1;
        Double_t                        fptcut = 0.2;
        Double_t                        fetacut = 0.9;
        Bool_t                          IsAA=kFALSE;
        THistManager*                   fHistos=nullptr; //!
        AliVMultiplicity*               fMultiplicity=nullptr; //!
        //TTree*                        fTree=nullptr;//!
        TH2*                            fForward2d=nullptr;//!
        TRandom3*                       fRandom=nullptr; //!
        Double_t                        sdweightingfactor=1.;
        
    

    ClassDef(AliAnalysisPseudoRapidityDensity, 1);
};

#endif

