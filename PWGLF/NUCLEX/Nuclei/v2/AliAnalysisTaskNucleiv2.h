#ifndef ALIANALYSISTASKNUCLEIV2_H
#define ALIANALYSISTASKNUCLEIV2_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//        AliAnalysisTaskNucleivn class for 2015 data
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;

#include <AliPIDResponse.h>
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliQnCorrectionsQnVector.h"


class AliAnalysisTaskNucleiv2 : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskNucleiv2();
  AliAnalysisTaskNucleiv2(const char *name);
  virtual ~AliAnalysisTaskNucleiv2() {}
  
  virtual void  UserCreateOutputObjects(); 
  virtual void  Initialize();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void SetIsPrimCut(Bool_t  isPrimCut)           {fisPrimCut   = isPrimCut; }; 
  void SetParticle(Int_t ptc)                    {fptc         = ptc;       };
  void SetVzMax(Float_t Vzmax)                   {fVzmax       = Vzmax;     };
  void SetCentralityEstimator( TString centEst)  {fCentrality  = centEst;   };
  void SetAnalysisType (const char* analysisType = "ESD")   { fAnalysisType = analysisType; }
  void SetYear(Int_t  year)                      {fYear  = year; }
  void SetHarmonic(Int_t  harmonic)              {fHarmonic  = harmonic; }
  
  Float_t GetPhi0Pi(Float_t phi);

  const AliQnCorrectionsQnVector *GetQnVectorFromList( const TList *list,
						       const char *subdetector,
						       const char *expectedstep,
						       const char *altstep);

  AliEventCuts fEventCuts;


 private:

  AliESDEvent *fESDevent;                         // 
  AliAODEvent *fAODevent;                         // 
  AliVEvent   *fevent;                            // 

  TString        fAnalysisType;                  //
  Bool_t         fisPrimCut;                     // Boolean : kTRUE = isprimarycut 
  Int_t          fptc;                           // Selected ptc 1 = d; 2 = t; 3 =3He 
  Float_t        fVzmax;                         // Selected vz max
  TString        fCentrality;                    //
  Bool_t         fApplyFlatten;                  //
  Int_t          fYear;                          // Year of data tacking
  Int_t          fHarmonic;                      // Flow Harmonic

  TList	*fListHist;	           // List of  histograms
 
  TH1F  *fHistEventMultiplicity;           // event multiplicity
  TH2F  *fHistTrackMultiplicity;           // track multiplicity
  TH2F  *fHistTrackMultiplicityCentral;    // track multiplicity
  TH2F  *fHistTrackMultiplicitySemiCentral;// track multiplicity
  TH2F  *fHistTrackMultiplicityMB;         // track multiplicity
  TH2F  *fHistTrackMultiplicityINT7;       // track multiplicity

  TH2F  *fhBB;                             // ScatterPlot Total
  TH2F  *fhBBDeu;                          // ScatterPlot Total
  TH2F  *fhTOF;                            // ScatterPlot Total TOF
  TH1F  *fhMassTOF;                        // Mass Distribution TOF
  
   // Event Plane vs Centrality

  TH2D *EPVzAvsCentrality  ; 
  TH2D *EPVzCvsCentrality  ; 
  TH2D *EPTPCvsCentrality  ; 
  TH2D *EPVzvsCentrality   ; 
  TH2D *EPTPCpvsCentrality ; 
  TH2D *EPTPCnvsCentrality ; 
  
  // EP TPC vs EP VZ for different centralities 

  TH2F *hEvPlaneTPCvsEvPVz010  ;                      
  TH2F *hEvPlaneTPCvsEvPVz1020 ; 
  TH2F *hEvPlaneTPCvsEvPVz2030 ;
  TH2F *hEvPlaneTPCvsEvPVz3040 ;                      
  TH2F *hEvPlaneTPCvsEvPVz4050 ;                      
  TH2F *hEvPlaneTPCvsEvPVz5060 ;   
  TH2F *hEvPlaneTPCvsEvPVz6080 ;   
  
  // For EP Resolution

  TH2F *hCos2DeltaTPCVzAvsCentrality;
  TH2F *hCos2DeltaTPCVzCvsCentrality;
  TH2F *hCos2DeltaVzAVzCvsCentrality;
  TH2F *hCos2DeltaVzMVzAvsCentrality;
  TH2F *hCos2DeltaVzMVzCvsCentrality;
  TH2F *hCos2DeltaVzATPCvsCentrality;
  TH2F *hCos2DeltaVzCTPCvsCentrality;
  TH2F *hCos2DeltaVzCVzAvsCentrality;
  TH2F *hCos2DeltaVzMTPCpvsCentrality;
  TH2F *hCos2DeltaVzMTPCnvsCentrality;
  TH2F *hCos2DeltaTPCpTPCnvsCentrality;

  // For SP resolution

  TH2F *  hQVzAQVzCvsCentrality;
  TH2F *  hQVzAQTPCvsCentrality;
  TH2F *  hQVzCQTPCvsCentrality;
  
// For NUA correction

  TH2F *hQxVzAvsCentrality;
  TH2F *hQyVzAvsCentrality;
  TH2F *hQxVzCvsCentrality;
  TH2F *hQyVzCvsCentrality;
  TH2F *hQxVzMvsCentrality;
  TH2F *hQyVzMvsCentrality;
  
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
  Double_t tCosdeltaphiTPC  ;
  Double_t tCosdeltaphiV0M  ;
  Double_t tCosdeltaphiV0A  ;
  Double_t tCosdeltaphiV0C  ;
  Double_t timpactXY        ;
  Double_t timpactZ         ;
  Double_t tpull            ;
  Double_t tphi             ;

  //---------------------------------------------------------------------------
  AliESDtrackCuts *fESDtrackCuts; 
  AliESDtrackCuts *fESDtrackCutsEP;
  AliPIDResponse  *fPIDResponse;   //! pointer to PID response
  //_______________________________________________________________________
 
  AliAnalysisTaskNucleiv2(const AliAnalysisTaskNucleiv2&);            // not implemented
  AliAnalysisTaskNucleiv2& operator=(const AliAnalysisTaskNucleiv2&); // not implemented
  
  ClassDef(AliAnalysisTaskNucleiv2, 1);
};

#endif
