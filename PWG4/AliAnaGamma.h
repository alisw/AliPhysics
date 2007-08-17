#ifndef ALIANAGAMMA_H
#define ALIANAGAMMA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for prompt gamma and correlation analysis
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TList.h> 
#include <TClonesArray.h> 
#include <TH2F.h>
#include<TObject.h>
#include <TTree.h>
#include <AliLog.h>

class AliGammaReader ;
class AliAnaGammaDirect ;
class AliAnaGammaCorrelation ;
class AliAnaGammaJetLeadCone ;
class AliNeutralMesonSelection ;

// --- AliRoot
class AliAnaGamma : public TObject {

public: 
  
  AliAnaGamma() ; // default ctor
  AliAnaGamma(const AliAnaGamma & g) ; // cpy ctor
  AliAnaGamma & operator = (const AliAnaGamma & g) ;//cpy assignment
  virtual ~AliAnaGamma() ; //virtual dtor

  enum anatype_t {kPrompt, kCorrelation};

  //General methods
  TList * GetOutputContainer()      const {return fOutputContainer ; }
  
  void Init();
  void InitParameters();

  Int_t GetAnalysisType(){  return fAnaType ; }
  void SetAnalysisType(Int_t ana ){  fAnaType = ana ; }

  void Print(const Option_t * opt) const;

  void MakeAnalysis(TClonesArray * plCalo, TClonesArray * plNe, TClonesArray * plCTS, TClonesArray *plParton)  ;  
  Bool_t ProcessEvent(Long64_t entry) ;
  //TTree * MakeTreeG(TString name) ;

  TString GetCalorimeter() {return fCalorimeter ; }
  void SetCalorimeter(TString calo) {if (calo == "PHOS" || calo == "EMCAL") fCalorimeter = calo ;
    else AliFatal("Wrong calorimeter name") ; }

  TObject * GetData() {return fData ; }
  TObject * GetKine() {return fKine ;}
  void SetData(TObject * data) {fData = data ; }
  void SetKine(TObject * kine) {fKine = kine ; }

  AliGammaReader * GetReader() {return fReader ; }
  void SetReader(AliGammaReader * reader) { fReader = reader ; }

  void SetGammaDirect(AliAnaGammaDirect * dg) { fGammaDirect = dg ; }
  void SetGammaCorrelation(AliAnaGammaCorrelation * gc) { fGammaCorrelation = gc ;}
  void SetNeutralMesonSelection(AliNeutralMesonSelection * nms) { fNeutralMesonSelection = nms ; }
  
 private:
  
  //General Data members
  TList  *fOutputContainer ; //! output data container
  Int_t      fAnaType; //Analysis type to be done
  TString fCalorimeter; //Prompt photon detector
  TObject * fData ; //! ESD
  TObject * fKine ; //! Stack
  AliGammaReader *      fReader ; //! Pointer to reader 
  AliAnaGammaDirect *   fGammaDirect ; //! Pointer to prompt gamma algorithm 
  AliAnaGammaCorrelation *   fGammaCorrelation ; //! Pointer to gamma correlation algorithm
  AliNeutralMesonSelection *  fNeutralMesonSelection ; //! Pointer to pair selection for pi0 identification.
  
  ClassDef(AliAnaGamma,0)
} ;
 

#endif //ALIANAGAMMA_H



