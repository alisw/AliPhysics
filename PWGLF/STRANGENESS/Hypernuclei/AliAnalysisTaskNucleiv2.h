#ifndef ALIANALYSISTASKNUCLEIV2_H
#define ALIANALYSISTASKNUCLEIV2_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskNucleiv2 class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TNtuple;
class AliESDcascade;
class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;

//class AliCascadeVertexer; 
#include <AliPIDResponse.h>
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskNucleiv2 : public AliAnalysisTaskSE {
 public:
  //  AliAnalysisTaskNucleiv2(const char *datatype);
  AliAnalysisTaskNucleiv2();
  AliAnalysisTaskNucleiv2(const char *name);
  virtual ~AliAnalysisTaskNucleiv2() {}
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void    SetCollidingSystems(Short_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void    SetAnalysisType(const char* analysisType)         {fAnalysisType = analysisType;}
  void    SetDataType(const char* dataType)                 {fDataType = dataType;}
  void    SetFillTree(Bool_t ifFill)                        {fFillNtuple = ifFill;}
  void    SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax); //select centrality
  Float_t GetEventPlaneForCandidate(AliESDtrack* track0, const TVector2* q,AliEventplane *pl, const TVector2* qsub1, const TVector2* qsub2);
  template <typename T> void           SetNullCuts(T* aod);
  void    PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
  Float_t GetPhi0Pi(Float_t phi);
  void    Initialize();
 private:
  
  TString  fAnalysisType;	      //! "ESD" or "AOD" analysis type	
  
  Short_t  fCollidingSystems;	      //! 0 = pp collisions or 1 = AA collisions
  TString  fDataType;		      //! "REAL" or "SIM" data type	
  Bool_t   fFillNtuple;		      //! fill or not the tree	
  
  Double_t fCentralityMin;            // lower bound of cenrality bin
  Double_t fCentralityMax;            // upper bound of centrality bin

  AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
  AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
  AliFlowEvent         *fFlowEvent; //! flow events 
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
  
  //From Flow Analysis
  TH2D *EPVzAvsCentrality  ; 
  TH2D *EPVzCvsCentrality  ; 
  TH2D *EPTPCvsCentrality  ; 
  TH2D *EPVzvsCentrality   ; 
  TH2D *EPTPCpvsCentrality ; 
  TH2D *EPTPCnvsCentrality ; 
  
  //------------------------------
  
  TProfile *fSubEventDPhiv205; 
  TProfile *fSubEventDPhiv2new05;
  
  TProfile *fSubEventDPhiv22040; 
  TProfile *fSubEventDPhiv2new2040;
  
  TProfile *fSubEventDPhiv24060; 
  TProfile *fSubEventDPhiv2new4060;
 
  TH2F *hCos2DeltaPhiVzAvsCentrality;
  TH2F *hCos2DeltaPhiVzCvsCentrality;
  TH2F *hCos2DeltaPhiVzMvsCentrality;
  TH2F *hCos2DeltaPhiTPCfvsCentrality;
  TH2F *hCos2DeltaPhiTPCpvsCentrality;
  TH2F *hCos2DeltaPhiTPCnvsCentrality;
 
  //---------------------------------------------------------------------------
  TH2F *hEvPlaneTPCvsEvPVz05;                      
  TH2F *hEvPlaneTPCvsEvPVz075; 
  TH2F *hEvPlaneTPCvsEvPVz1530;
  TH2F *hEvPlaneTPCvsEvPVz3050;                      
  TH2F *hEvPlaneTPCvsEvPVz2040;                      
  TH2F *hEvPlaneTPCvsEvPVz4060;                      
  
  
  // From D meson analysis

  TH2F *hCos2DeltaPhivsPt075;                      
  TH2F *hCos2DeltaPhiVZEROvsPt075;                 
  TH2F *hCos2DeltaPhivsPt1530;                      
  TH2F *hCos2DeltaPhiVZEROvsPt1530;                 
  TH2F *hCos2DeltaPhivsPt3050;                      
  TH2F *hCos2DeltaPhiVZEROvsPt3050;                 
 
  TH2F *hCos2DeltaPhivsPt05;                      
  TH2F *hCos2DeltaPhiVZEROvsPt05;                 
  TH2F *hCos2DeltaPhivsPt2040;                      
  TH2F *hCos2DeltaPhiVZEROvsPt2040;                 
  TH2F *hCos2DeltaPhivsPt4060;                      
  TH2F *hCos2DeltaPhiVZEROvsPt4060;                 
  
  
  AliESDtrackCuts * fESDtrackCuts; 
  AliESDtrackCuts * fESDtrackCutsEP; 
  AliPIDResponse *fPIDResponse;   //! pointer to PID response

  //_______________________________________________________________________

  TTree *fNtuple1;                //! Some Information on the tracks
    
  Double_t tCentrality[2] ;       // Centrality + event type
  Double_t tPulls[3];             // Pulls
 
  Double_t tMomentum[4];          //pxpypz of the tracks +  corrected pt for d
  Double_t tDCA[2];               //dcaXY and dcaZ of the track
  Int_t    tisTOF[2];             //isTOF, isOuterTPCwall   
  Double_t tTOFtrack[3];          //poutTPC,timeTOF,trackLenghtTOF;
  Int_t    tCharge;               //Charge of the Track (Pos or Neg)
  Double_t tPhi;                  //Phi 
  Double_t trpangleTPC;           //rpangleTPC
  Double_t trpangleVZERO[3];      //rpangleVZERO: V0M, V0A, V0C
  
  // MC releted quantities

  Double_t tPDGCode;              //PDG code ptc
  Double_t tPDGCodeMum;           //PDG code mother ptc
  Double_t tIsPrimaryTr;
  Double_t tIsSecondaryTr[2];     //from material ; from weak deacy 
  
  //_______________________________________________________________________
  
  TTree *fNtuple2;                  //! MC tree

  Double_t tCentralityMC;
  Double_t tVertexCoordMC[3];
  Double_t tMomentumMC[3];          //pxpypz of the tracks

  Double_t tPDGCodeMC      ;
  Double_t tPDGCodeMumMC   ;
  Double_t tIsPrimary      ;
  Double_t tIsSecondary[2] ;      //from material ; from weak deacy 
  Double_t tEtaMC          ;
  Double_t tPtMC           ;
  Double_t tYMC            ;

  //_______________________________________________________________________
 
  AliAnalysisTaskNucleiv2(const AliAnalysisTaskNucleiv2&);            // not implemented
  AliAnalysisTaskNucleiv2& operator=(const AliAnalysisTaskNucleiv2&); // not implemented
  
  ClassDef(AliAnalysisTaskNucleiv2, 0);
};

#endif
