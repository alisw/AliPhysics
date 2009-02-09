/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* History of cvs commits:
*
* $Log$
* Revision 1.3  2007/10/16 14:36:39  pavlinov
* fixed code violation (almost)
*
* Revision 1.2  2007/09/11 19:38:15  pavlinov
* added pi0 calibration, linearity, shower profile
* co: warning: `$Log' is obsolescent; use ` * $Log'.

* Revision 1.1  2007/08/08 15:58:01  hristov
* New calibration classes. They depend on TTable, so libTable.so is added to the list of Root libraries. (Aleksei)
*/

//_________________________________________________________________________
//    Set of parameters for pi0 selection 
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALPi0SelectionParam.h"

ClassImp(AliEMCALPi0SelectionParRec)
//_________________________________________________________________________
AliEMCALPi0SelectionParRec::AliEMCALPi0SelectionParRec() : 
fEOfRpMin(0.3), fEOfRpMax(30.), fMassGGMin(0.03), fMassGGMax(0.28), fMomPi0Min(1.8), fMomPi0Max(12.)
{
  // Default constructor 
}




ClassImp(AliEMCALPi0SelectionParam)
//_________________________________________________________________________
AliEMCALPi0SelectionParam::AliEMCALPi0SelectionParam() : TNamed("",""), fTable(0), fCurrentInd(0)
{
  // Default constructor 
}

//_________________________________________________________________________
AliEMCALPi0SelectionParam::AliEMCALPi0SelectionParam(const AliEMCALPi0SelectionParam& param) 
  : TNamed(param), fTable(param.fTable), fCurrentInd(param.fCurrentInd)
{
  // Copy constructor 
}

//_________________________________________________________________________
AliEMCALPi0SelectionParam::AliEMCALPi0SelectionParam(const char* name, const Int_t nrow) : TNamed(name,"table of cell information") , fTable(0), fCurrentInd(0)
{
  // Oct 16, 2007
  fTable = new TObjArray(nrow);
}

//_________________________________________________________________________
void AliEMCALPi0SelectionParam::AddAt(AliEMCALPi0SelectionParRec* r)
{
  // Oct 16, 2007
  (*fTable)[fCurrentInd] = new AliEMCALPi0SelectionParRec(*r);
  fCurrentInd++;
}

//_________________________________________________________________________
AliEMCALPi0SelectionParam::~AliEMCALPi0SelectionParam()
{
  // Oct 16, 2007
  if(fTable) {
    fTable->Delete();
    delete fTable;
  }
}

//_________________________________________________________________________
AliEMCALPi0SelectionParRec* AliEMCALPi0SelectionParam::GetTable(Int_t i) const
{
  // Oct 16, 2007
  return (AliEMCALPi0SelectionParRec*)fTable->At(i);
}

//_________________________________________________________________________
void AliEMCALPi0SelectionParam::PrintTable()
{
  // Oct 16, 2007
  printf(" Table : %s : nrows %i \n", GetName(), int(GetNRows()));
  for(int i=0; i<GetNRows(); i++) PrintTable(i);
}

//_________________________________________________________________________
void AliEMCALPi0SelectionParam::PrintTable(const Int_t i)
{
  // Oct 16, 2007
  if(i>=GetNRows()) return;
  printf("row %i \n", i);
  PrintRec(GetTable(i));
}

//_________________________________________________________________________
void AliEMCALPi0SelectionParam::PrintRec(AliEMCALPi0SelectionParRec* r)
{
  // Oct 16, 2007
  if(r==0) return;
  printf(" cluster  energy  window %7.2f -> %7.2f \n", r->fEOfRpMin, r->fEOfRpMax);
  printf(" gamma,gamma mass window %7.2f -> %7.2f \n", r->fMassGGMin, r->fMassGGMax);
  printf(" pi0   momentum   window %7.2f -> %7.2f \n", r->fMomPi0Min, r->fMomPi0Max);
}

//_________________________________________________________________________
// Set 1;
AliEMCALPi0SelectionParam* AliEMCALPi0SelectionParam::Set1()
{
  // Initial set of pars of Pi0 selection
  AliEMCALPi0SelectionParRec r;
  r.fEOfRpMin  = 0.3; 
  r.fEOfRpMax  = 30.;
  r.fMassGGMin = 0.03;  
  r.fMassGGMax = 0.28; 
  r.fMomPi0Min = 1.8;
  r.fMomPi0Max = 12.0;
  AliEMCALPi0SelectionParam *t = new AliEMCALPi0SelectionParam("Pi0Set1",1);
  t->AddAt(&r);
  return t;
}
