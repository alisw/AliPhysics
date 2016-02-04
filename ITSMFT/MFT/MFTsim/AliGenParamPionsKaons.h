#ifndef AliGenParamPionsKaons_H
#define AliGenParamPionsKaons_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Parametric generator of primary pions and kaons
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliGenerator.h"
#include "TH1D.h"
#include "TH2D.h"

class AliGenParamPionsKaons : public AliGenerator {

public:
  
  AliGenParamPionsKaons();
  AliGenParamPionsKaons(Int_t nPart, Char_t *inputFile);

  virtual ~AliGenParamPionsKaons() {}
  virtual void Generate();
  virtual void Init();
  virtual void SetPionOnly() { fGeneratePion=kTRUE; fGenerateKaon=kFALSE; }
  virtual void SetKaonOnly() { fGeneratePion=kFALSE; fGenerateKaon=kTRUE; }
  virtual void LoadInputHistos(Char_t *inputFile);

private:

  AliGenParamPionsKaons(const AliGenParamPionsKaons&);
  AliGenParamPionsKaons &operator=(const AliGenParamPionsKaons&);

protected:
  
  Bool_t fGeneratePion;
  Bool_t fGenerateKaon;

  TH2D *fPtVsRapidityPrimaryPosPions,  *fPtVsRapidityPrimaryNegPions;
  TH2D *fPtVsRapidityPrimaryPosKaons,  *fPtVsRapidityPrimaryNegKaons;

  TH1D *fHistPdgCode;

  ClassDef(AliGenParamPionsKaons, 1)

};

//====================================================================================================================================================

#endif


