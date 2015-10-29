#ifndef ALIDIELECTRONCFDRAW_H
#define ALIDIELECTRONCFDRAW_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#             Class AliDielectronCF                         #
//#       Dielectron Correction Framework Draw helper         #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################



#include <TNamed.h>
#include <TVectorD.h>
#include <AliCFContainer.h>
#include "AliDielectronVarManager.h"

class TObjArray;
class TSeqCollection;
class AliCFEffGrid;
class TH1;

class AliCFContainer;

class AliDielectronCFdraw : public TNamed {
public:
  enum ECollectType { kSE=0, kME, kMEOS, kROT, kAll };
  
  AliDielectronCFdraw();
  AliDielectronCFdraw(const char* name, const char* title);
  AliDielectronCFdraw(AliCFContainer *cont);
  AliDielectronCFdraw(const char* filename);
  
  virtual ~AliDielectronCFdraw();
  
  void SetCFContainer(AliCFContainer * const container) {fCfContainer=container;}
  void SetCFContainers(const TSeqCollection *arr);

  void SetCFContainers(const char* filename);

  AliCFContainer* GetCFContainer() const {return fCfContainer;}

  void SetRangeUser(Int_t ivar, Double_t min, Double_t max, const char* slices="");
  void SetRangeUser(const char* varname, Double_t min, Double_t max, const char* slices="");
  void SetRangeUser(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, const char* slices="", Bool_t leg=kFALSE);

  void UnsetRangeUser(Int_t ivar, const char* slices="");
  void UnsetRangeUser(const char* varname, const char* slices="");
  void UnsetRangeUser(AliDielectronVarManager::ValueTypes type, const char* slices="", Bool_t leg=kFALSE);

  TString FindSteps(const char* search="");
  Int_t   FindStep(const char* search="");
  Int_t   FindVar(AliDielectronVarManager::ValueTypes type, Bool_t leg=kFALSE);

  virtual void Draw(const Option_t* varnames = "") { Draw(varnames,"");}
  virtual void Print(const Option_t*) const { if (fCfContainer) fCfContainer->Print(""); }
  //Draw Projections
  void Draw(const Option_t* varnames, const char* opt, const char* slices="");
  void Draw(Int_t var, const char* opt="", const char* slices="");
  void Draw(Int_t var0, Int_t var1, const char* opt="", const char* slices="");
  void Draw(Int_t var0, Int_t var1, Int_t var2, const char* opt="", const char* slices="");

  TObjArray* CollectHistosProj(const Option_t* varnames, const char* slices);
  TObjArray* CollectHistosProj(const Int_t vars[3], const char* slices);
  TObjArray* CollectMinvProj(Int_t slice, ECollectType collect=kAll, TString var="M");
  TH1* Project(const Int_t vars[3], Int_t slice);
  TH1* Project(const Option_t* var, Int_t slice);
  
  //Draw efficiencies
  void DrawEfficiency(const char* varnames, const char* numerators, Int_t denominator=0, const char* opt="sameleg2");
  void DrawEfficiency(Int_t var, const char* numerators, Int_t denominator=0, const char* opt="sameleg", Int_t type=0);
  void DrawEfficiency(Int_t var0, Int_t var1, const char* numerators, Int_t denominator=0, const char* opt="sameleg", Int_t type=0);
  void DrawEfficiency(Int_t var0, Int_t var1, Int_t var2, const char* numerators, Int_t denominator=0, const char* opt="sameleg", Int_t type=0);
  
  TObjArray* CollectHistosEff(const Int_t vars[3], const char* numerators, Int_t denominator, Int_t type=0);
  TH1* ProjectEff(const Int_t vars[3]);

  Double_t GetAverageEfficiency(Int_t numerator, Int_t denominator, Double_t &effErr);
  
  const TVectorD& GetData() const {return fVdata;}
  void Draw(const TObjArray *arr, const char* opt="");
private:
  AliCFContainer *fCfContainer;                     // CF container
  AliCFEffGrid   *fEffGrid;                         // Efficiency calculation

  TVectorD fVdata;                                  // vector with data, like mean efficiencies
  
  AliDielectronCFdraw(const AliDielectronCFdraw &c);
  AliDielectronCFdraw &operator=(const AliDielectronCFdraw &c);
  
  ClassDef(AliDielectronCFdraw,0)                   // CF draw helper class
};

//
// Inline functions
//

inline Int_t AliDielectronCFdraw::FindVar(AliDielectronVarManager::ValueTypes type, Bool_t leg)
{
  //
  // find variable number in CF container
  //
  return ( fCfContainer->GetVar(Form("%s%s", leg?"Leg1_":"", AliDielectronVarManager::GetValueName(type))) );
}
#endif

