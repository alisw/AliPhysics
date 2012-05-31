#ifndef ALIITSALIGNMILLE2CONSTRARRAY_H
#define ALIITSALIGNMILLE2CONSTRARRAY_H

/*-------------------------------------------------------------------------------------
Set of gaussian constraints on LOCAL corrections of the module parameters defined 
for a set of modules. 

Author: ruben.shahoyan@cern.ch
--------------------------------------------------------------------------------------*/
#include "AliITSAlignMille2Constraint.h"
class AliITSAlignMille2Modle;

class AliITSAlignMille2ConstrArray : public AliITSAlignMille2Constraint
{
 public:
  enum {kTypeGaussian=10};
  //
  AliITSAlignMille2ConstrArray();
  AliITSAlignMille2ConstrArray(const Char_t* name,Double_t *parcf,Int_t npar,Double_t val,Double_t err);
  virtual ~AliITSAlignMille2ConstrArray() {}
  //
  Double_t     GetError()               const {return fError;}
  Int_t        GetNModules()            const {return fModuleIDs.GetSize();}
  Int_t        GetNCoeffs()             const {return fCoeffs.GetSize();}
  Int_t        GetModuleID(Int_t i)     const {return fModuleIDs[i];}
  Double_t     GetCoeff(Int_t i)        const {return fCoeffs[i];}
  void         Print(Option_t* opt="")  const;
  //
  void         AddModule(AliITSAlignMille2Module* mod,Bool_t needGeom = kTRUE);
  void         SetError(Double_t err)         {fError = err;}
  //
  virtual Bool_t IncludesModule(Int_t id)            const;
  virtual Bool_t IncludesModPar(Int_t id,Int_t par)  const;
  virtual Bool_t IncludesModPar(const AliITSAlignMille2Module* mod, Int_t par) const;
  //
 protected:
  AliITSAlignMille2ConstrArray(const AliITSAlignMille2ConstrArray& src);
  AliITSAlignMille2ConstrArray& operator=(const AliITSAlignMille2ConstrArray& ) {return *this;}
  //
 protected:
  TArrayS           fModuleIDs;         // module id' to apply this constraint
  TArrayS           fModulePatt;        // pattern of variables involved (in the frame of varied params)
  TArrayD           fCoeffs;            // weight for each param          
  Double_t          fError;             // constraint error
  
  //
  ClassDef(AliITSAlignMille2ConstrArray,0)
};


#endif
