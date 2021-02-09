#ifndef ALIDIELECTRONCUTGROUP_H
#define ALIDIELECTRONCUTGROUP_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#################################################################
//#                                                               #
//#             Class AliDielectronCutGroup                       #
//#              Dielectron Group of cuts                         #
//#                                                               #
//#  Authors:                                                     #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de                 #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de                 #
//#   Christoph Baumann   uni Ffm / cbaumann@ikf.uni-frankfurt.de #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch           #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch     #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch             #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de                   #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch          #
//#                                                               #
//#################################################################

#include <AliAnalysisCuts.h>
#include <TList.h>

class TCollection;

class AliDielectronCutGroup : public AliAnalysisCuts {
  
public:
  enum TruthValues {
    kCompAND = kTRUE,
    kCompOR = kFALSE
  };
  
  AliDielectronCutGroup(Bool_t compOperator=kCompOR);
  AliDielectronCutGroup(const char*name, const char* title, Bool_t compOperator=kCompOR);
  
  virtual ~AliDielectronCutGroup();
  
  Int_t GetNCuts() const { return fCutGroupList.GetEntries(); }
  
  //Analysis cuts interface
  //
  virtual void Init();
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TObject* track, Double_t* values);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  
  void AddCut(AliAnalysisCuts* fCut);
  void SetCompOperator(Bool_t compOperator);

  virtual void Print(const Option_t* option = "") const;
  const AliAnalysisCuts* GetCut(Int_t iCut) const;

  
private:
  TList  fCutGroupList;  //for saving the different cuts
  Bool_t fCompOperator;  //determines whether the cuts are AND/OR compared
  
  ClassDef(AliDielectronCutGroup,2) //Group of cuts
};

#endif
