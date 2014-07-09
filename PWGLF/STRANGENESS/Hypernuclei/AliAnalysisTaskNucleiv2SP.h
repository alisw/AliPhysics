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
 
  AliAnalysisTaskNucleiv2SP();
  AliAnalysisTaskNucleiv2SP(const char *name); //select ptc 1 = d; 2 = t ; 3 = 3He
  virtual ~AliAnalysisTaskNucleiv2SP() {}
  
  virtual void  UserCreateOutputObjects();
  virtual void  Initialize();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
   
  Float_t GetEventPlaneForCandidate(AliESDtrack* track0, const TVector2* q,AliEventplane *pl);
  Float_t GetPhi0Pi(Float_t phi);

  void SetIsPrimCut        (Bool_t  isPrimCut           = kFALSE) { fisPrimCut       = isPrimCut;         } 
  void SetParticle         (Float_t ptc                 = 1.    ) { fptc             = ptc;               }
  void SetMaxPull          (Float_t pull                = 3.    ) { fmaxpull         = pull;               }
  void SetMaxVz            (Float_t maxVz               = 10.   ) { fmaxVz           = maxVz;               }

 private:

  Bool_t         fisPrimCut;                     // Boolean : kTRUE = isprimarycut 
  Float_t        fptc;                           // Selected ptc 1 = d; 2 = t; 3 =3He 
  Float_t        fmaxpull;                       // Selected ptc 1 = d; 2 = t; 3 =3He 
  Float_t        fmaxVz;                       // Selected ptc 1 = d; 2 = t; 3 =3He 

  TList	*fListHist;	           //! List of  histograms
 
  TH1F  *fHistEventMultiplicity;           //! event multiplicity
  TH2F  *fHistTrackMultiplicity;           //! track multiplicity
  TH2F  *fHistTrackMultiplicityCentral;    //! track multiplicity
  TH2F  *fHistTrackMultiplicitySemiCentral;//! track multiplicity
  TH2F  *fHistTrackMultiplicityMB;         //! track multiplicity

  TH2F  *fhBB;                             //! ScatterPlot Total
  TH2F  *fhBBDeu;                          //! ScatterPlot Total
  TH2F  *fhPtDeu;                          //! correctet vs non correcter d pt
  TH2F  *fhTOF;                            //! ScatterPlot Total TOF
  TH1F  *fhMassTOF;                        //! Mass Distribution TOF
  
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

  // THnSparse
  
  THnSparse * fHistRealTracks; 

  //---------------------------------------------------------------------------
  AliESDtrackCuts *fESDtrackCuts; 
  AliESDtrackCuts *fESDtrackCutsEP;
  AliPIDResponse  *fPIDResponse;   //! pointer to PID response
  //_______________________________________________________________________
 
  AliAnalysisTaskNucleiv2SP(const AliAnalysisTaskNucleiv2SP&);            // not implemented
  AliAnalysisTaskNucleiv2SP& operator=(const AliAnalysisTaskNucleiv2SP&); // not implemented
  
  ClassDef(AliAnalysisTaskNucleiv2SP, 0);
};

#endif
