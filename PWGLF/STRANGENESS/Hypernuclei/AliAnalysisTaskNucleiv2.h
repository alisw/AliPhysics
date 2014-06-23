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
class TNtuple;
class AliESDcascade;
//class AliCascadeVertexer; 
#include <AliPIDResponse.h>
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskNucleiv2 : public AliAnalysisTaskSE {
 public:
  //  AliAnalysisTaskNucleiv2(const char *datatype);
  AliAnalysisTaskNucleiv2();
  AliAnalysisTaskNucleiv2(const char *name,const char *datatype,Bool_t filltree);
  virtual ~AliAnalysisTaskNucleiv2() {}
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType    (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
  void SetDataType    (const char* dataType = "REAL") {fDataType = dataType;}
  void SetFillTree    (Bool_t ifFill = kFALSE) {fFillNtuple = ifFill;}
  //  void SetDataType    (const char* dataTypeS = datatype) {fDataType = dataTypeS;}
  // void SetDataType    (const char* dataType) {fDataType = dataType;}
 
  //Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t optMC);
  //  Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t isPbPb);
  //  Bool_t IsTrackAccepted(AliESDtrack * const track);
  Float_t GetPhi0Pi(Float_t phi);
  void    Initialize();
 private:
  
  TString fAnalysisType;	      //! "ESD" or "AOD" analysis type	
  
  Short_t fCollidingSystems;	      //! 0 = pp collisions or 1 = AA collisions
  TString fDataType;		      //! "REAL" or "SIM" data type	
  Bool_t  fFillNtuple;		      //! fill or not the tree	
 
  TList	*fListHistCascade;	           //! List of Cascade histograms
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

  TH2F *hRPangleTPCvsCentrality;           //RESOLUTION Histrograms
  TH2F *hPlaneResoTPCvsCentrality;
  TH2F *hRPangleVZEROvsCentrality;
  TH2F *hRPangleVZEROAvsCentrality;
  TH2F *hRPangleVZEROCvsCentrality;
  TH2F *hPlaneResoVZEROvsCentrality;
  TH2F *hPlaneResoVZEROAvsCentrality;
  TH2F *hPlaneResoVZEROCvsCentrality;
  TH2F *hCosPhivsPt;                         
  TH2F *hSinPhivsPt;                         
  TH2F *hPhivsPt;                         
  TH2F *hAntiCosPhivsPt;                     
  TH2F *hAntiSinPhivsPt;                     
  TH2F *hAntiPhivsPt;                     
  TH2F *hCosDeltaPhivsPt075;                      
  TH2F *hSinDeltaPhivsPt075;                      
  TH2F *hDeltaPhivsPt075;                      
  TH2F *hCosDeltaPhiVZEROvsPt075;                 	      
  TH2F *hSinDeltaPhiVZEROvsPt075;                 	      
  TH2F *hDeltaPhiVZEROvsPt075;                 
  TH2F *hCosDeltaPhivsPt1530;                     
  TH2F *hSinDeltaPhivsPt1530;                     
  TH2F *hDeltaPhivsPt1530;                     
  TH2F *hCosDeltaPhiVZEROvsPt1530;                
  TH2F *hSinDeltaPhiVZEROvsPt1530;                
  TH2F *hDeltaPhiVZEROvsPt1530;                
  TH2F *hCosDeltaPhivsPt3050;                     
  TH2F *hSinDeltaPhivsPt3050;                     
  TH2F *hDeltaPhivsPt3050;                     
  TH2F *hCosDeltaPhiVZEROvsPt3050;                
  TH2F *hSinDeltaPhiVZEROvsPt3050;                
  TH2F *hDeltaPhiVZEROvsPt3050;                
  TH2F *hCosDeltaPhivsPt1550;                     
  TH2F *hSinDeltaPhivsPt1550;                     
  TH2F *hDeltaPhivsPt1550;                     
  TH2F *hCosDeltaPhiVZEROvsPt1550;                
  TH2F *hSinDeltaPhiVZEROvsPt1550;                
  TH2F *hDeltaPhiVZEROvsPt1550;                
  TH2F *hAntiCosDeltaPhivsPt075;                  
  TH2F *hAntiSinDeltaPhivsPt075;                  
  TH2F *hAntiDeltaPhivsPt075;                  
  TH2F *hAntiCosDeltaPhiVZEROvsPt075;             
  TH2F *hAntiSinDeltaPhiVZEROvsPt075;             
  TH2F *hAntiDeltaPhiVZEROvsPt075;             
  TH2F *hAntiCosDeltaPhivsPt1530;                 
  TH2F *hAntiSinDeltaPhivsPt1530;                 
  TH2F *hAntiDeltaPhivsPt1530;                 
  TH2F *hAntiCosDeltaPhiVZEROvsPt1530;            
  TH2F *hAntiSinDeltaPhiVZEROvsPt1530;            
  TH2F *hAntiDeltaPhiVZEROvsPt1530;            
  TH2F *hAntiCosDeltaPhivsPt3050;                 
  TH2F *hAntiSinDeltaPhivsPt3050;                 
  TH2F *hAntiDeltaPhivsPt3050;                 
  TH2F *hAntiCosDeltaPhiVZEROvsPt3050;            
  TH2F *hAntiSinDeltaPhiVZEROvsPt3050;            
  TH2F *hAntiDeltaPhiVZEROvsPt3050;            
  TH2F *hAntiCosDeltaPhivsPt1550;                 
  TH2F *hAntiSinDeltaPhivsPt1550;                 
  TH2F *hAntiDeltaPhivsPt1550;                 
  TH2F *hAntiCosDeltaPhiVZEROvsPt1550;            
  TH2F *hAntiSinDeltaPhiVZEROvsPt1550;            
  TH2F *hAntiDeltaPhiVZEROvsPt1550;     

  AliESDtrackCuts * fESDtrackCuts; 
  AliPIDResponse *fPIDResponse;   //! pointer to PID response

  //_______________________________________________________________________

  TTree *fNtuple1;                //! Some Information on the tracks
    
  Double_t tEventNumber[7];       //ev number; run number; Bunch Cross - Orbit -  Period Number; #tracks; event type
  Double_t tCentrality ;
  Double_t tVertexCoord[3];
  
  //TRACKS
  Double_t tPIDITS[9];            //n-signas different particle species ITS
  Double_t tPIDTPC[9];            //n-signas different particle species TPC
  Double_t tPIDTOF[9];            //n-signas different particle species TOF
  Double_t tPulls[3];          //Pulls
 
  Double_t tMomentum[3];          //pxpypz of the tracks
  Double_t tTPCMomentum;          //momentum at theh TPC inner wall
  Double_t tdEdx;                 //dEdx of the tracks
  Double_t tEta;                  //eta of the tracks
  Double_t tDCA[2];               //dcaXY and dcaZ of the track
  Double_t tTracksTPC[2];         //chi^2 and #TPCcluster
  Double_t tITSclustermap;        //ITS cluster map
  Double_t tITSsample[4];         //ITS samples 
  Int_t    tisTOF[2];             //isTOF, isOuterTPCwall   
  Double_t tTOFtrack[3];          //poutTPC,timeTOF,trackLenghtTOF;
  Int_t    tCharge;               //Charge of the Track
  Double_t tPtCorr;               //Corrected Momentum 
  Double_t tPhi;                  //Phi 
  Double_t trpangleTPC;           //rpangleTPC
  Double_t trpangleVZERO[3];      //rpangleVZERO: V0M, V0A, V0C
  
 
  Double_t tPDGCode;              //PDG code ptc
  Double_t tPDGCodeMum;           //PDG code mother ptc
  Double_t tIsPrimaryTr;
  Double_t tIsSecondaryTr[2];     //from material ; from weak deacy 
  
  //_______________________________________________________________________
  
  TTree *fNtuple2;                  //! MC tree

  Double_t tEventNumberMC[6];       //ev number; run number; Bunch Cross - Orbit -  Period Number; #tracks
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
