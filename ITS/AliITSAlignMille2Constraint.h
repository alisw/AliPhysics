#ifndef ALIITSALIGNMILLE2CONSTRAINT_H
#define ALIITSALIGNMILLE2CONSTRAINT_H

#include <TNamed.h>
#include <TArrayS.h>
#include <TArrayD.h>

class AliITSAlignMille2Module;

/*-----------------------------------------------------------------------------------------
Simple constraint on the subunits of the module ID (if ID>=0) or all modules w/o 
parents (ID=-1): the mean or median of the GLOBAL corrections of each parameter requested
in the pattern must be = 0. When added explicitly to the fit it requires addition of 
Lagrange multipliers which may require more powerfull matrix preconditioners. For this 
reason we usually ommit the constrain from explicit fit and apply it afterwards to obtained
parameters (with median constraint this is the only method possible) 

Author: ruben.shahoyan@cern.ch
------------------------------------------------------------------------------------------*/
class AliITSAlignMille2Constraint : public TNamed
{
 public:
  enum {kTypeMean,kTypeMedian};
  enum {kDisabledBit=31};
  //
  AliITSAlignMille2Constraint();
  AliITSAlignMille2Constraint(const Char_t* name,Int_t t,Int_t mdID,Double_t val=0,UInt_t pattern=0x0ffff);
  virtual ~AliITSAlignMille2Constraint() {}
  //
  UInt_t       GetConstraintID()        const {return GetUniqueID();}
  Int_t        GetType()                const {return fType;}
  Int_t        GetModuleID()            const {return fModuleID;}
  Double_t     GetValue()               const {return fVal;}
  UInt_t       GetPattern()             const {return fPattern;}
  UInt_t       GetAppliedPattern()      const {return fApplied;}
  UInt_t       GetRemainingPattern()    const {return (~fApplied)&GetPattern();}
  Bool_t       IsApplied(Int_t par)     const {return fApplied & (0x1<<par);}
  Bool_t       IsApplied()              const {return (fApplied&0xffff)==GetPattern();}
  Bool_t       IncludesParam(int id)    const {return fPattern&BIT(id);}
  void         Print(Option_t* opt="")  const;
  //
  void         SetConstraintID(UInt_t id)            {SetUniqueID(id);}
  void         SetType(Int_t t)                      {fType = t;}
  void         SetPattern(UInt_t pat)                {fPattern = pat;}
  void         SetValue(Double_t val)                {fVal = val;}
  void         SetApplied(Int_t par)                 {fApplied |= par<0 ? 0x0ffff : (0x1<<par);}
  void         Disable()                             {fApplied |= 0x1<<kDisabledBit;}
  void         Enable()                              {fApplied &= ~(0x1<<kDisabledBit);}
  Bool_t       IsDisabled()             const        {return (fApplied>>kDisabledBit)&0x1;}
  //
  virtual Bool_t IncludesModule(Int_t id)            const {return fModuleID==id;}
  virtual Bool_t IncludesModPar(Int_t id,Int_t par)  const {return IncludesModule(id) && IncludesParam(par);}
  virtual Bool_t IncludesModPar(const AliITSAlignMille2Module* mod, Int_t par) const;
  //
 protected:
  AliITSAlignMille2Constraint(const AliITSAlignMille2Constraint& src);
  AliITSAlignMille2Constraint& operator=(const AliITSAlignMille2Constraint& ) {return *this;}
  //
 protected:
  Int_t             fType;              // constriant type: mean, median ...
  Double_t          fVal;               // constraint value
  Int_t             fModuleID;          // Id of the module involved, -1 for orphans
  UInt_t            fApplied;           // was it already applied?
  UInt_t            fPattern;           // pattern of params involved
  //
  ClassDef(AliITSAlignMille2Constraint,0)
};


#endif
