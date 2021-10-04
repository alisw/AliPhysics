
#ifndef AliAnalysisMultPt_H
#define AliAnalysisMultPt_H

#include "AliAnalysisTaskSE.h"
#define PI 3.1415927


//class AliAODEvent;
//class AliMCEvent;
//class TList;
//class TTree;
//class TH3D;



class AliAnalysisMultPt : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisMultPt();
                                AliAnalysisMultPt(const char *name);
        virtual                 ~AliAnalysisMultPt();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                    ProcessMCParticles();
 
    private:
        AliAODEvent*            fAOD;             //! input event
        TList*                  fOutputList;      //! output list
        TTree*                  mixer;
        TH2D*                   fHistMultPt;      //! Pt-Mult histogram
        TH2D*                   fHistMultPtMC;    //! Pt-Mult MC histogram
        TH1D*                   fHistMult;        //! Multiplicity histogram
        TH1D*                   fHistMultMC;      //! Multiplicity MC histogram
        TH1D*                   fHistPt;          //! Pt histogram
        TH1D*                   fHistPtMC;        //! Pt MC histogram
        TH1D*                   fHistFit;         //! Fit histogram
        TH1D*                   fHistFitMC;       //! Fit MC histogram
        TH1D*                   fHistProj;        //! Projection histogram
        TH1D*                   fHistProjMC;      //! Projection MC histogram
        TH1D*                   fHistRatio;       //! Ratio of pTs
        //TH2D*                   fHistRatio;       //! Ratio of pTs
        AliMCEvent*             fMCEvent;         //! corresponding MC event
       
    
       Bool_t fIsMC;

    
      /*
        Bool_t fIsMC;
        Double_t fChi2DoF;
        Int_t fTPCNcls;

        AliEventCuts fEventCuts;
        Float_t mNMC;
        Float_t mSpTMC;
    */
    
     
    
       int mMode;
      

   AliAnalysisMultPt(const AliAnalysisMultPt&); // not implemented
   AliAnalysisMultPt& operator=(const AliAnalysisMultPt&); // not implemented

        ClassDef(AliAnalysisMultPt, 1);    /// yaaaa1 ya 2   ???????????
};

#endif
