/*-----------------------------------------------------------------------------------------
Simple constraint on the subunits of the module ID (if ID>=0) or all modules w/o 
parents (ID=-1): the mean or median of the GLOBAL corrections of each parameter requested
in the pattern must be = 0. When added explicitly to the fit it requires addition of 
Lagrange multipliers which may require more powerfull matrix preconditioners. For this 
reason we usually ommit the constrain from explicit fit and apply it afterwards to obtained
parameters (with median constraint this is the only method possible) 

Author: ruben.shahoyan@cern.ch
------------------------------------------------------------------------------------------*/
#include "AliITSAlignMille2Constraint.h"
#include "AliITSAlignMille2Module.h"



ClassImp(AliITSAlignMille2Constraint)

//________________________________________________________________________________________________________
AliITSAlignMille2Constraint::AliITSAlignMille2Constraint() :
TNamed(),
fType(kTypeMean),
fVal(0),
fModuleID(0),
fApplied(0),
fPattern(0)
{}

//________________________________________________________________________________________________________
AliITSAlignMille2Constraint::AliITSAlignMille2Constraint(const Char_t* name,Int_t t,Int_t mdID,Double_t val,UInt_t pattern) :
TNamed(name,""),
fType(t),
fVal(val),
fModuleID(mdID),
fApplied(0),
fPattern(pattern)
{
}

//________________________________________________________________________________________________________
AliITSAlignMille2Constraint::AliITSAlignMille2Constraint(const AliITSAlignMille2Constraint& src) :
TNamed(src),
fType(src.fType),
fVal(src.fVal),
fModuleID(src.fModuleID),
fApplied(src.fApplied),
fPattern(src.fPattern)
{/* DUMMY */} 

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2Constraint::IncludesModPar(const AliITSAlignMille2Module* mod, Int_t par) const
{
  // is this module/parameter mentioned in the list?
  if (!IncludesParam(par)) return kFALSE;
  if (fModuleID==-1 && !mod->GetParent()) return kTRUE;
  return IncludesModule( mod->GetUniqueID() );
}


//________________________________________________________________________________________________________
void AliITSAlignMille2Constraint::Print(Option_t* ) const
{
  // print data
  printf("#%3d Constraint %s of type %d on module %d to value %+e\n",GetConstraintID(),GetName(),GetType(),GetModuleID(),GetValue());
  printf("Paremeters: ");
  for (int i=0;i<=8;i++) if (TestBit(0x1<<i)) printf("%d ",i); printf("\n");
  //
}

