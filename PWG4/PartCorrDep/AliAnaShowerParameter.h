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
class AliAnaShowerParameter : public AliAnaPartCorrBaseClass {

public: 

  AliAnaShowerParameter() ; // default ctor
  virtual ~AliAnaShowerParameter() ; //virtual dtor
  
  TList *  GetCreateOutputObjects();

  void Init();

  TObjString* GetAnalysisCuts();
    
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
 
  TString fCalorimeter ;    // Calorimeter where the gamma is searched;
  Float_t fNCellsCutMin ;   // N cells cut min
  Float_t fNCellsCutMax ;   // N cells cut max
  Float_t fLambdaCut ;      // l0 cut 
  Double_t fTimeCutMin  ;   // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;   // Remove clusters/cells with time larger than this value, in ns

  //Histograms   
  TH1F * fhNClusters  ;     //! cluster
  TH2F * fhNCellCluster;    //! cells per cluster
  TH3F * fhEtaPhiPtCluster; //! eta vs phi vs pt
  TH2F * fhLambdaPtCluster; //! l0 vs pt

  //MC
  TH2F * fhLambdaPtPhoton ; //! l0 vs pt mc photon   
  TH2F * fhLambdaPtPi0 ;    //! l0 vs pt mc pi0 
  TH2F * fhLambdaPtPion ;   //! l0 vs pt mc pi charged
  TH1D * fhPtTruthPi0 ;     //! pi0 pt mc

  AliAnaShowerParameter(const AliAnaShowerParameter & g) ;               // cpy ctor
  AliAnaShowerParameter & operator = (const AliAnaShowerParameter & g) ; // cpy assignment
  
   ClassDef(AliAnaShowerParameter,1)

} ;
 

#endif//AliAnaShowerParameter_H



