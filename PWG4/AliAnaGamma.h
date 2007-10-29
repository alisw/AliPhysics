#ifndef ALIANAGAMMA_H
#define ALIANAGAMMA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for gamma and correlation analysis
// It is called by the task class AliAnalysisGammaTask and it connects the input (ESD/AOD/MonteCarlo)
// got with AliGammaReader (produces TClonesArrays of TParticles), with the analysis classes 
// AliAnaGammaDirect, AliAnaGammaCorrelation ....
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TList.h> 
#include <TClonesArray.h> 
#include <TH2F.h>
#include<TObject.h>
#include <TTree.h>

// --- AliRoot system ---
#include <AliLog.h>
#include "AliAODCaloCluster.h"

class AliGammaReader ;
class AliAnaGammaDirect ;
class AliAnaGammaCorrelation ;
class AliAnaGammaSelection ;
class AliNeutralMesonSelection ;
class AliAODEvent;

// --- AliRoot
class AliAnaGamma : public TObject {

public: 
  
  AliAnaGamma() ; // default ctor
  AliAnaGamma(const AliAnaGamma & g) ; // cpy ctor
  AliAnaGamma & operator = (const AliAnaGamma & g) ;//cpy assignment
  virtual ~AliAnaGamma() ; //virtual dtor

  enum Anatype {kPrompt, kCorrelation};

  //Setter and getters
  TList * GetOutputContainer()      const {return fOutputContainer ; }

  Int_t GetAnalysisType() const {return fAnaType ; }
  void SetAnalysisType(Int_t ana ){  fAnaType = ana ; }

  TString GetCalorimeter() const {return fCalorimeter ; }
  void SetCalorimeter(TString calo) {if (calo == "PHOS" || calo == "EMCAL") fCalorimeter = calo ;
    else AliFatal("Wrong calorimeter name") ; }

  TObject * GetData() const {return fData ; }
  TObject * GetKine() const {return fKine ;}
  void SetData(TObject * data) {fData = data ; }
  void SetKine(TObject * kine) {fKine = kine ; }

  AliGammaReader * GetReader() const {return fReader ; }
  void SetReader(AliGammaReader * reader) { fReader = reader ; }

  AliAnaGammaDirect * GetGammaDirect() const { return fGammaDirect ; }
  AliAnaGammaCorrelation * GetGammaCorrelation() const { return fGammaCorrelation  ;}
  AliAnaGammaSelection * GetGammaSelection() const { return fGammaSelection ;}
  AliNeutralMesonSelection * GetNeutralMesonSelection() const { return fNeutralMesonSelection  ; }

  void SetGammaDirect(AliAnaGammaDirect * dg) { fGammaDirect = dg ; }
  void SetGammaCorrelation(AliAnaGammaCorrelation * gc) { fGammaCorrelation = gc ;}
  void SetGammaSelection(AliAnaGammaSelection * gs) { fGammaSelection = gs ;}
  void SetNeutralMesonSelection(AliNeutralMesonSelection * nms) { fNeutralMesonSelection = nms ; }

  //AOD stuff  
  void   AddCluster(AliAODCaloCluster p);
  void   ConnectAOD(AliAODEvent* aod);
  void   FillAODs(TClonesArray * plPHOS, TClonesArray * plEMCAL);

  //Others
  void Init();
  void InitParameters();

  void MakeAnalysis(TClonesArray * plCalo, TClonesArray * plNe, TClonesArray * plCTS, TClonesArray *plParton, TClonesArray * plPrimCalo)  ;  

  void Print(const Option_t * opt) const;

  Bool_t ProcessEvent(Long64_t entry) ;

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
  AliAnaGammaSelection *   fGammaSelection ; //! Pointer to gamma selection algorithm
  AliNeutralMesonSelection *  fNeutralMesonSelection ; //! Pointer to pair selection for pi0 identification.
  TClonesArray* fAODclusters;        //! reconstructed jets
  Int_t         fNAODclusters;       //! number of reconstructed jets

  ClassDef(AliAnaGamma,1)
} ;
 

#endif //ALIANAGAMMA_H



