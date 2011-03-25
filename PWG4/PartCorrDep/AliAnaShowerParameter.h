#ifndef ALIANASHOWERPARAMETER_H
#define ALIANASHOWERPARAMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaShowerParameter.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
//
// Class cloned from AliAnaPhoton, main aim is shower shape studies
// 
// 
//
//-- Author: Jocelyn Mlynarz (WSU) and Gustavo Conesa (LPSC)

// --- ROOT system ---
class TH3F;
class TH2F ;
class TH1F;
class TString ;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
//#include "AliStack.h"
//#include "TParticle.h"
class AliStack;
class TParticle;

class TList ;
class AliEMCALGeoUtils;
class AliAnaShowerParameter : public AliAnaPartCorrBaseClass {

public: 

  AliAnaShowerParameter() ; // default ctor
  AliAnaShowerParameter(const AliAnaShowerParameter & g) ; // cpy ctor
  AliAnaShowerParameter & operator = (const AliAnaShowerParameter & g) ;//cpy assignment
  virtual ~AliAnaShowerParameter() ; //virtual dtor
  
  TList *  GetCreateOutputObjects();

  void Init();

  TObjString* GetAnalysisCuts();

  void MakeAnalysisFillAOD()  ;
    
  void MakeAnalysisFillHistograms() ; 
  
  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
	
  void InitParameters();

  void SetTimeCut(Double_t min, Double_t max) {fTimeCutMin = min; fTimeCutMax = max;}
  Double_t GetTimeCutMin() const {return fTimeCutMin;}
  Double_t GetTimeCutMax() const {return fTimeCutMax;}	
  
  void SetNCellsCut(Double_t min, Double_t max) {fNCellsCutMin = min; fNCellsCutMax = max;}
  Double_t GetNCellsCutMin() const {return fNCellsCutMin;}
  Double_t GetNCellsCutMax() const {return fNCellsCutMax;}	
	
  private:
 
  TString fCalorimeter ; // Calorimeter where the gamma is searched;
  Float_t fNCellsCutMin ;
  Float_t fNCellsCutMax ;
  Float_t fLambdaCut ;
  Double_t fTimeCutMin  ;    // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;    // Remove clusters/cells with time larger than this value, in ns

  //Histograms   
  TH1F * fhNClusters  ; 
  TH2F * fhNCellCluster;
  TH3F * fhEtaPhiPtCluster   ; 
  TH2F * fhLambdaPtCluster  ; 

  //MC
  TH2F * fhLambdaPtPhoton ;    
  TH2F * fhLambdaPtPi0 ;    
  TH2F * fhLambdaPtPion ; 
  TH1D * fhPtTruthPi0 ;

   ClassDef(AliAnaShowerParameter,1)

} ;
 

#endif//AliAnaShowerParameter_H



