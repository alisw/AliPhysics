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

/* $Log$co: warning: `/* $Log' is obsolescent; use ` * $Log'.

 * Revision 1.1  2007/08/08 15:58:01  hristov
 * New calibration classes. They depend on TTable, so libTable.so is added to the list of Root libraries. (Aleksei)
 * */

//_________________________________________________________________________
//    Set of parameters for pi0 selection 
//
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include "AliEMCALPi0SelectionParam.h"

ClassImp(pi0SelectionParam)
pi0SelectionParam::pi0SelectionParam() : 
eOfRpMin(0.3), eOfRpMax(30.), massGGMin(0.03), massGGMax(0.28), momPi0Min(1.8), momPi0Max(12.)
{ }
//_________________________________________________________________________

ClassImp(AliEMCALPi0SelectionParam)
AliEMCALPi0SelectionParam::AliEMCALPi0SelectionParam() : TNamed("",""), fTable(0), fCurrentInd(0)
{
}

AliEMCALPi0SelectionParam::AliEMCALPi0SelectionParam(const char* name, const Int_t nrow) : TNamed(name,"table of cell information") , fTable(0), fCurrentInd(0)
{
  fTable = new TObjArray(nrow);
}

void AliEMCALPi0SelectionParam::AddAt(pi0SelectionParam* r)
{
  (*fTable)[fCurrentInd] = new pi0SelectionParam(*r);
  fCurrentInd++;
}

AliEMCALPi0SelectionParam::~AliEMCALPi0SelectionParam()
{
  if(fTable) {
    fTable->Delete();
    delete fTable;
  }
}

pi0SelectionParam* AliEMCALPi0SelectionParam::GetTable(Int_t i) const
{
  return (pi0SelectionParam*)fTable->At(i);
}

void AliEMCALPi0SelectionParam::PrintTable()
{
  printf(" Table : %s : nrows %i \n", GetName(), int(GetNRows()));
  for(int i=0; i<GetNRows(); i++) PrintTable(i);
}

void AliEMCALPi0SelectionParam::PrintTable(const Int_t i)
{
  if(i>=GetNRows()) return;
  printf("row %i \n", i);
  PrintRec(GetTable(i));
}

void AliEMCALPi0SelectionParam::PrintRec(pi0SelectionParam* r)
{
  if(r==0) return;
  printf(" cluster  energy  window %7.2f -> %7.2f \n", r->eOfRpMin, r->eOfRpMax);
  printf(" gamma,gamma mass window %7.2f -> %7.2f \n", r->massGGMin, r->massGGMax);
  printf(" pi0   momentum   window %7.2f -> %7.2f \n", r->momPi0Min, r->momPi0Max);
}
// Set 1;
AliEMCALPi0SelectionParam* AliEMCALPi0SelectionParam::Set1()
{
  pi0SelectionParam r;
  r.eOfRpMin  = 0.3; 
  r.eOfRpMax  = 30.;
  r.massGGMin = 0.03;  
  r.massGGMax = 0.28; 
  r.momPi0Min = 1.8;
  r.momPi0Max = 12.0;
  AliEMCALPi0SelectionParam *t = new AliEMCALPi0SelectionParam("Pi0Set1",1);
  t->AddAt(&r);
  return t;
}
