#ifndef ALIEMCALPI0SELECTIONPARAM_H
#define ALIEMCALPI0SELECTIONPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliEMCALPi0SelectionParam.h 24500 2008-03-13 23:39:38Z jklay $ */

//_________________________________________________________________________
//  Set of parameters for pi0 selection 
// unit is GeV
// pi0SelectionParam -> AliEMCALPi0SelectionParRec
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---
#include <TNamed.h>
#include <TObjArray.h>

// pi0SelectionParam -> AliEMCALPi0SelectionParRec

class  AliEMCALPi0SelectionParRec :  public TObject{
  friend class AliEMCALPi0SelectionParam;
  friend class AliEMCALPi0Calibration;
 public:
  AliEMCALPi0SelectionParRec();
  virtual ~AliEMCALPi0SelectionParRec() {};
  virtual const char* GetName() const {return "Pi0Par";}

 protected:
  double fEOfRpMin;   // minimal energy of em.cluster (rec point)
  double fEOfRpMax;   // maximal energy of em.cluster (rec point)
  double fMassGGMin;  // minimal mass of gamma,gamma
  double fMassGGMax;  // maximal mass of gamma,gamma
  double fMomPi0Min;  // minimal pi0 momentum
  double fMomPi0Max;  // maximal pi0 momentum

  ClassDef(AliEMCALPi0SelectionParRec,1);
};


class AliEMCALPi0SelectionParam : public TNamed {
 public:
  AliEMCALPi0SelectionParam(); // default constractor
  AliEMCALPi0SelectionParam(const AliEMCALPi0SelectionParam& param);
  AliEMCALPi0SelectionParam(const char* name, const Int_t nrow);
  virtual ~AliEMCALPi0SelectionParam();

  AliEMCALPi0SelectionParam & operator = (const AliEMCALPi0SelectionParam  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };
  // 
  void AddAt(AliEMCALPi0SelectionParRec* r);
  AliEMCALPi0SelectionParRec* GetTable(Int_t i) const;
  Int_t       GetSize()  const {return fTable->GetSize();}
  Int_t       GetNRows() const {return fCurrentInd;}

 // Menu
  void PrintTable();                 // *MENU*
  void PrintTable(const Int_t i);    // *MENU*
  void PrintRec(AliEMCALPi0SelectionParRec *r);

  // Set of parameter(s)
  static AliEMCALPi0SelectionParam* Set1();
  //
 protected:
  TObjArray *fTable; // Table of AliEMCALPi0SelectionParRec
  Int_t fCurrentInd; // Current index

  ClassDef(AliEMCALPi0SelectionParam, 2) // Set of Parameters For Pi0 Selection     
};

#endif // ALIEMCALPI0SELECTIONPARAM_H
