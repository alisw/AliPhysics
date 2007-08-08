#ifndef ALIEMCALPI0SELECTIONPARAM_H
#define ALIEMCALPI0SELECTIONPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//    Set of parameters for pi0 selection 
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---
#include <TTable.h>

// unit is GeV
struct  pi0SelectionParam {
  double eOfRpMin;   // minimal energy of em.cluster (rec point)
  double eOfRpMax;   // maximal energy of em.cluster (rec point)
  double massGGMin;  // minimal mass of gamma,gamma
  double massGGMax;  // maximal mass of gamma,gamma
  double momPi0Min;  // minimal pi0 momentum
  double momPi0Max;  // maximal pi0 momentum
};

class AliEMCALPi0SelectionParam : public TTable {
 public:
  ClassDefTable(AliEMCALPi0SelectionParam , pi0SelectionParam)

 // Menu
  void PrintTable();                 // *MENU*
  void PrintTable(const Int_t i);    // *MENU*
  void PrintRec(pi0SelectionParam *r);

  // Set of parameter(s)
  static AliEMCALPi0SelectionParam* Set1();

  ClassDef(AliEMCALPi0SelectionParam,1) // Set of Parameters For Pi0 Selection     
};

#endif // ALIEMCALPI0SELECTIONPARAM_H
