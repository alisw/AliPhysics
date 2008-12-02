#ifndef ALIANALYSISTASKSTRANGE_H
#define ALIANALYSISTASKSTRANGE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskStrange class
//       This task is for single strange study from ESD/AOD
//          Origin: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr
//-----------------------------------------------------------------

class TString;
class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliESDEvent;
class AliESDVertex;
class AliAODEvent;

class AliAnalysisTaskStrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskStrange();
  AliAnalysisTaskStrange(const char *name,  const char *optCuts);
  virtual ~AliAnalysisTaskStrange() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetCollidingSystems(Int_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  
 private:
  TString      fAnalysisType;                   //  ESD or AOD
  Int_t        fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb  
  TString      fOption;                         //  Option for the selections
  TList       *fListHist;                       //! List of histograms
  TH1F        *fHistPrimaryVertexPosX;          //! Primary vertex position in X
  TH1F        *fHistPrimaryVertexPosY;          //! Primary vertex position in Y
  TH1F        *fHistPrimaryVertexPosZ;          //! Primary vertex position in Z
  TH1F        *fHistTrackMultiplicity;          //! V0 multiplicity distribution
  TH1F        *fHistV0Multiplicity;

  TH2F        *fHistDcaPosToPrimVertex;         //! Dca of V0 positive daughter to primary vertex
  TH2F        *fHistDcaNegToPrimVertex;         //! Dca of V0 negative daughter to primary vertex
  TH2F        *fHistDcaPosToPrimVertexZoom;     //! Zoom
  TH2F        *fHistDcaNegToPrimVertexZoom;     //! Zoom
  TH2F        *fHistRadiusV0;                   //! V0 radial distance distribution
  TH2F        *fHistDecayLengthV0;              //! V0 decay length distribution
  TH2F        *fHistDcaV0Daughters;             //! Dca between V0 daughters
  TH2F        *fHistChi2;                       //! V0 chi2 distribution
  TH2F        *fHistCosPointAngle;              //! Cosine of V0 pointing angle
  TH2F        *fHistCosPointAngleZoom;          //! Zoom

              // V0 offline distributions
  TH1F        *fHistV0MultiplicityOff;          //! V0 multiplicity distribution offline
  TH2F        *fHistPtVsYK0sOff;                //! Pt vs.Y with K0s assumption
  TH2F        *fHistPtVsYLambdaOff;             //! Pt vs.Y with Lambda assumption
  TH2F        *fHistPtVsYAntiLambdaOff;         //! Pt vs.Y with Anti-Lambda assumption
  TH1F        *fHistMassK0sOff;                 //! Invariant mass of K0s
  TH1F        *fHistMassLambdaOff;		//! Invariant mass of Lambda
  TH1F        *fHistMassAntiLambdaOff;		//! Invariant mass of Anti-Lambda
  TH2F        *fHistMassVsRadiusK0sOff;         //! Invariant mass vs. radius of K0s	 
  TH2F        *fHistMassVsRadiusLambdaOff;	//! Invariant mass vs. radius of Lambda	 
  TH2F        *fHistMassVsRadiusAntiLambdaOff;	//! Invariant mass vs. radius of Anti-Lambda
  TH2F        *fHistPtVsMassK0sOff;             //! Pt vs. invariant mass of K0s	    
  TH2F        *fHistPtVsMassLambdaOff;		//! Pt vs. invariant mass of Lambda	    
  TH2F        *fHistPtVsMassAntiLambdaOff;	//! Pt vs. invariant mass of Anti-Lambda
  TH2F        *fHistArmenterosPodolanskiOff;    //! Armenteros-Podolanski distribution       

              // V0 on-the-fly distributions
  TH1F        *fHistV0MultiplicityOn;          //! V0 multiplicity distribution on-the-fly
  TH2F        *fHistPtVsYK0sOn;                //! Pt vs.Y with K0s assumption		    
  TH2F        *fHistPtVsYLambdaOn;	       //! Pt vs.Y with Lambda assumption	    
  TH2F        *fHistPtVsYAntiLambdaOn;	       //! Pt vs.Y with Anti-Lambda assumption	    
  TH1F        *fHistMassK0sOn;                 //! Invariant mass of K0s		    
  TH1F        *fHistMassLambdaOn;	       //! Invariant mass of Lambda		     
  TH1F        *fHistMassAntiLambdaOn;	       //! Invariant mass of Anti-Lambda	    
  TH2F        *fHistMassVsRadiusK0sOn;	       //! Invariant mass vs. radius of K0s	    
  TH2F        *fHistMassVsRadiusLambdaOn;      //! Invariant mass vs. radius of Lambda	    
  TH2F        *fHistMassVsRadiusAntiLambdaOn;  //! Invariant mass vs. radius of Anti-Lambda
  TH2F        *fHistPtVsMassK0sOn;	       //! Pt vs. invariant mass of K0s	    
  TH2F        *fHistPtVsMassLambdaOn;	       //! Pt vs. invariant mass of Lambda	    
  TH2F        *fHistPtVsMassAntiLambdaOn;      //! Pt vs. invariant mass of Anti-Lambda    
  TH2F        *fHistArmenterosPodolanskiOn;    //! Armenteros-Podolanski distribution
   
  AliAnalysisTaskStrange(const AliAnalysisTaskStrange&);            // not implemented 
  AliAnalysisTaskStrange& operator=(const AliAnalysisTaskStrange&); // not implemented 

  ClassDef(AliAnalysisTaskStrange, 1); 
};

#endif
