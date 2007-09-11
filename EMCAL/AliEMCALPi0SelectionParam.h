#ifndef ALIEMCALPI0SELECTIONPARAM_H
#define ALIEMCALPI0SELECTIONPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Set of parameters for pi0 selection 
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---
#include <TNamed.h>
#include <TObjArray.h>

// unit is GeV
class  pi0SelectionParam :  public TObject{
 public:
  pi0SelectionParam();
  virtual ~pi0SelectionParam() {};
  virtual const char* GetName() const {return "Pi0Par";}

  double eOfRpMin;   // minimal energy of em.cluster (rec point)
  double eOfRpMax;   // maximal energy of em.cluster (rec point)
  double massGGMin;  // minimal mass of gamma,gamma
  double massGGMax;  // maximal mass of gamma,gamma
  double momPi0Min;  // minimal pi0 momentum
  double momPi0Max;  // maximal pi0 momentum

  ClassDef(pi0SelectionParam,1);
};

class AliEMCALPi0SelectionParam : public TNamed {
 public:
  AliEMCALPi0SelectionParam(); // default constractor
  AliEMCALPi0SelectionParam(const char* name, const Int_t nrow);
  virtual ~AliEMCALPi0SelectionParam();

  AliEMCALPi0SelectionParam & operator = (const AliEMCALPi0SelectionParam  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };
  // 
  void AddAt(pi0SelectionParam* r);
  pi0SelectionParam* GetTable(Int_t i) const;
  Int_t       GetSize()  const {return fTable->GetSize();}
  Int_t       GetNRows() const {return fCurrentInd;}

 // Menu
  void PrintTable();                 // *MENU*
  void PrintTable(const Int_t i);    // *MENU*
  void PrintRec(pi0SelectionParam *r);

  // Set of parameter(s)
  static AliEMCALPi0SelectionParam* Set1();
  //
 protected:
  TObjArray *fTable;
  Int_t fCurrentInd;

  ClassDef(AliEMCALPi0SelectionParam, 2) // Set of Parameters For Pi0 Selection     
};

#endif // ALIEMCALPI0SELECTIONPARAM_H
