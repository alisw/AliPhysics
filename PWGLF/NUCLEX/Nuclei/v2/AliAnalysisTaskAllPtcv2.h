#ifndef ALIANALYSISTASKAllPtcV2_H
#define ALIANALYSISTASKAllPtcV2_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//          AliAnalysisTaskAllPtcv2 data for 2015 data 
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class AliAnalysisUtils;

#include <AliPIDResponse.h>
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliEventCuts.h"
#include "AliQnCorrectionsQnVector.h"


class AliAnalysisTaskAllPtcv2 : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskAllPtcv2();
  AliAnalysisTaskAllPtcv2(const char *name); 
  virtual ~AliAnalysisTaskAllPtcv2() {}
  
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

  TList	*fListHist;	           //! List of  histograms
 
  TH1F  *fHistEventMultiplicity;           // event multiplicity
  TH2F  *fHistTrackMultiplicity;           // track multiplicity
  TH2F  *fHistTrackMultiplicityCentral;    // track multiplicity
  TH2F  *fHistTrackMultiplicitySemiCentral;// track multiplicity
  TH2F  *fHistTrackMultiplicityMB;         // track multiplicity
  TH2F  *fHistTrackMultiplicityINT7;       // track multiplicity

  TH2F  *fhBB;                             // ScatterPlot Total
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

  TH2F *hEvPlaneTPCvsEvPVz0010 ;                      
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

  TH2F *hQVzAQVzCvsCentrality;

  // For NUA correction

  TH2F *hQxVzAvsCentrality;
  TH2F *hQyVzAvsCentrality;
  TH2F *hQxVzCvsCentrality;
  TH2F *hQyVzCvsCentrality;
  TH2F *hQxVzMvsCentrality;
  TH2F *hQyVzMvsCentrality;
  
  //Histograms ptc
  TH2F *hCosdeltaphiTPCvsPtAll0010;      
  TH2F *hCosdeltaphiV0AvsPtAll0010;      
  TH2F *hCosdeltaphiV0CvsPtAll0010;   
  TH2F *huqV0AvsPtAll0010;
  TH2F *huqV0CvsPtAll0010;
 
  TH2F *hCosdeltaphiTPCvsPtAll1020;      
  TH2F *hCosdeltaphiV0AvsPtAll1020;      
  TH2F *hCosdeltaphiV0CvsPtAll1020;       
  TH2F *huqV0AvsPtAll1020;
  TH2F *huqV0CvsPtAll1020;
 
  TH2F *hCosdeltaphiTPCvsPtAll2030;      
  TH2F *hCosdeltaphiV0AvsPtAll2030;      
  TH2F *hCosdeltaphiV0CvsPtAll2030;      
  TH2F *huqV0AvsPtAll2030;
  TH2F *huqV0CvsPtAll2030;
 
  TH2F *hCosdeltaphiTPCvsPtAll3040;      
  TH2F *hCosdeltaphiV0AvsPtAll3040;      
  TH2F *hCosdeltaphiV0CvsPtAll3040;       
  TH2F *huqV0AvsPtAll3040;
  TH2F *huqV0CvsPtAll3040;
  
  TH2F *hCosdeltaphiTPCvsPtAll4050;      
  TH2F *hCosdeltaphiV0AvsPtAll4050;      
  TH2F *hCosdeltaphiV0CvsPtAll4050;      
  TH2F *huqV0AvsPtAll4050;
  TH2F *huqV0CvsPtAll4050;

  TH2F *hCosdeltaphiTPCvsPtAll5060;      
  TH2F *hCosdeltaphiV0AvsPtAll5060;      
  TH2F *hCosdeltaphiV0CvsPtAll5060;      
  TH2F *huqV0AvsPtAll5060;
  TH2F *huqV0CvsPtAll5060;

  TH2F *hCosdeltaphiTPCvsPtAll6080;      
  TH2F *hCosdeltaphiV0AvsPtAll6080;      
  TH2F *hCosdeltaphiV0CvsPtAll6080;      
  TH2F *huqV0AvsPtAll6080;
  TH2F *huqV0CvsPtAll6080;
  
  TH2F *hphivsPtAll0010;
  TH2F *hphivsPtAll1020;
  TH2F *hphivsPtAll2030;
  TH2F *hphivsPtAll3040;
  TH2F *hphivsPtAll4050;
  TH2F *hphivsPtAll5060;
  TH2F *hphivsPtAll6080;
  
  TH2F *huqV0AvsuqV0C0010;
  TH2F *huqV0AvsuqV0C1020;
  TH2F *huqV0AvsuqV0C2030;
  TH2F *huqV0AvsuqV0C3040;
  TH2F *huqV0AvsuqV0C4050;
  TH2F *huqV0AvsuqV0C5060;
  TH2F *huqV0AvsuqV0C6080;
  
  TH2F *huqV0AxuqV0CvsPtAll0010;
  TH2F *huqV0AxuqV0CvsPtAll1020;
  TH2F *huqV0AxuqV0CvsPtAll2030;
  TH2F *huqV0AxuqV0CvsPtAll3040;
  TH2F *huqV0AxuqV0CvsPtAll4050;
  TH2F *huqV0AxuqV0CvsPtAll5060;
  TH2F *huqV0AxuqV0CvsPtAll6080;

  //---------------------------------------------------------------------------
  AliESDtrackCuts *fESDtrackCuts; 
  AliESDtrackCuts *fESDtrackCutsEP;
  AliPIDResponse  *fPIDResponse;   //! pointer to PID response
  //_______________________________________________________________________
 
  AliAnalysisTaskAllPtcv2(const AliAnalysisTaskAllPtcv2&);            // not implemented
  AliAnalysisTaskAllPtcv2& operator=(const AliAnalysisTaskAllPtcv2&); // not implemented
  
  ClassDef(AliAnalysisTaskAllPtcv2, 1);
};

#endif
