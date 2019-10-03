#ifndef ALIANALYSISTASKNUCLEIV2SP_H
#define ALIANALYSISTASKNUCLEIV2SP_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskNucleiv2SP class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;


#include <AliPIDResponse.h>
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisTaskNucleiv2SP : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskNucleiv2SP(const char *name); 
  AliAnalysisTaskNucleiv2SP();
  virtual ~AliAnalysisTaskNucleiv2SP() {}
      
  Float_t GetEventPlaneForCandidate(AliVTrack* track0, const TVector2* q,AliEventplane *pl);
  // Float_t GetEventPlaneForCandidate(AliAODtrack* track0, const TVector2* q,AliEventplane *pl);
  Float_t GetPhi0Pi(Float_t phi);
  Bool_t  Flatten(Float_t cent);
    
  virtual void  UserCreateOutputObjects();
  virtual void  Initialize();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);

  void SetIsPrimCut(Bool_t  isPrimCut)           {fisPrimCut   = isPrimCut; }; 
  void SetParticle(Int_t ptc)                    {fptc         = ptc;       };
  void SetVzMax(Float_t Vzmax)                   {fVzmax       = Vzmax;     };
  void SetCentralityEstimator( TString centEst)  {fCentrality  = centEst;   };
  void SetAnalysisType (const char* analysisType = "ESD")   { fAnalysisType = analysisType; }
  void SetApplyFlatten(Bool_t  applyFlatten)    {fApplyFlatten  = applyFlatten; }
  void SetYear(Int_t  year)    {fYear  = year; }
  void SetHarmonic(Int_t  harmonic)    {fHarmonic  = harmonic; }
  
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

  TH2F *hEvPlaneTPCvsEvPVz05;                      
  TH2F *hEvPlaneTPCvsEvPVz075; 
  TH2F *hEvPlaneTPCvsEvPVz1530;
  TH2F *hEvPlaneTPCvsEvPVz3050;                      
  TH2F *hEvPlaneTPCvsEvPVz2040;                      
  TH2F *hEvPlaneTPCvsEvPVz4060;                      

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

  TH2F *hQVzAQVzCvsCentrality;

  // For NUA correction

  TH2F *hQxVzAvsCentrality;
  TH2F *hQyVzAvsCentrality;
  TH2F *hQxVzCvsCentrality;
  TH2F *hQyVzCvsCentrality;
  TH2F *hQxVzMvsCentrality;
  TH2F *hQyVzMvsCentrality;
  
  // Flattness

  TH2F *hqEPCvsCentrality; 
  TH2F *hqEPAvsCentrality;
  TH2F *hqEPvsCentrality;
 
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
  
  AliAnalysisTaskNucleiv2SP(const AliAnalysisTaskNucleiv2SP&);            // not implemented
  AliAnalysisTaskNucleiv2SP& operator=(const AliAnalysisTaskNucleiv2SP&); // not implemented

  ClassDef(AliAnalysisTaskNucleiv2SP, 1);
};


#endif
