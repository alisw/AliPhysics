#ifndef ALIANALYSISTASKNUCLEIV2PPB_H
#define ALIANALYSISTASKNUCLEIV2PPB_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TString.h>
#include <THnSparse.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliPIDResponse.h"


class AliAnalysisTaskNucleiv2pPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNucleiv2pPb();
  AliAnalysisTaskNucleiv2pPb(const char *name);

  virtual ~AliAnalysisTaskNucleiv2pPb();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 


  Float_t GetPhi0Pi(Float_t phi);
 
  void SetParticle(Int_t ptc)                    {fptc       = ptc;      };
  void SetVzMax(Float_t Vzmax)                   {fVzmax     = Vzmax;    };
  void SetCentralityEstimator(Short_t centEst)   {fCenCalV0  = centEst;  };
  void SetHarmonic(Int_t  harmonic)              {fNHarm     = harmonic; };
  void SetRecPass(Short_t recPass)               {fRecPass   = recPass;  };
  void SetFilterBit(Short_t filterBit)           {fFilterBit = filterBit;};
  

 private:
    void Analyze(AliVEvent* aod);
    void OpenInfoCalbration(Int_t run);

    AliAODEvent *fAODevent;                         // 
    AliVEvent   *fevent;                            // 
    
    Int_t        fRun;                // run number - for calibration
    TH1D*        fMultV0;             // profile from V0 multiplicity
    
    TH1D*        fQxnmV0A;            // <Qx2> V0A
    TH1D*        fQynmV0A;            // <Qy2> V0A
    TH1D*        fQxnsV0A;            // sigma Qx2 V0A
    TH1D*        fQynsV0A;            // sigma Qy2 V0A
    TH1D*        fQxnmV0C;            // <Qx2> V0C
    TH1D*        fQynmV0C;            // <Qy2> V0C
    TH1D*        fQxnsV0C;            // sigma Qx2 V0C
    TH1D*        fQynsV0C;            // sigma Qy2 V0C
    
    Double_t     fNHarm;              // harmonic number: 2, 3
    Short_t      fRecPass;            // flag for reconstruction pass: 0->Fast, 1->NoSDD, 2->SDD
    Short_t      fCenCalV0;           // flag for centrality estimators used for V0 recentering: 0->V0A, 1->V0, 2->V0AEq, 3-> CL1
    Short_t      fFilterBit;          // flag for AOD filterbit

    Int_t fptc ;
    Float_t fVzmax;

    
    //output hist
    TList	*fListHist;	           // List of  histograms
 
    TH1F  *fHistEventMultiplicity;           // event multiplicity
    TH2F  *fHistTrackMultiplicity;           // track multiplicity
    
    TH2F  *fhBB;                             // ScatterPlot Total
    TH2F  *fhBBDeu;                          // ScatterPlot Total
    TH2F  *fhTOF;                            // ScatterPlot Total TOF
    TH1F  *fhMassTOF;                        // Mass Distribution TOF
  
    // Event Plane vs Centrality
    
    TH2D *EPVzAvsCentrality  ; 
    TH2D *EPVzCvsCentrality  ; 
  
    // For SP resolution
    
    TH2F *hQVzAQVzCvsCentrality;
    
    // For NUA correction
    
    TH2F *hQxVzAvsCentrality;
    TH2F *hQyVzAvsCentrality;
    TH2F *hQxVzCvsCentrality;
    TH2F *hQyVzCvsCentrality;
     
    Int_t eventtype;

    // TTree
    TTree *ftree;                //! Some Information on the tracks
    Double_t tCentrality      ;
    Double_t tType            ;
    Double_t tHasTOF          ;
    Double_t tpT              ;
    Double_t tMassTOF         ;
    Double_t tuqV0A           ;
    Double_t tuqV0C           ;
    Double_t tCharge          ;
    Double_t tCosdeltaphiV0A  ;
    Double_t tCosdeltaphiV0C  ;
    Double_t timpactXY        ;
    Double_t timpactZ         ;
    Double_t tpull            ;
    Double_t tphi             ;
    
    //---------------------------------------------------------------------------
    AliPIDResponse  *fPIDResponse;   //! pointer to PID response
  
    
    ClassDef(AliAnalysisTaskNucleiv2pPb, 1);    //Analysis task for high pt analysis
    
};

#endif
