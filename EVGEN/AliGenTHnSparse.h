#ifndef ALIGENTHNSPARSE_H
#define ALIGENTHNSPARSE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Particle generator according to 4 correlated variables : here
// z, ptot, r, theta. The input is a THnSparse object included in
// the root file (path and name to be set via the SetTHnSparse method).
// This class is similar to AliGenFunction.

#include "AliGenerator.h"
#include "THnSparse.h"

class AliGenTHnSparse : public AliGenerator
{
public:

  AliGenTHnSparse();
  AliGenTHnSparse(const AliGenTHnSparse& func);
  AliGenTHnSparse &operator=(const AliGenTHnSparse& func);
  virtual ~AliGenTHnSparse();
  virtual void Generate();
  virtual void Init();
  virtual void SetPart(Int_t part, Bool_t both) {fIpart=part; fBoth=both;}
  virtual void SetThnSparse(char *file_name, char *thn_name);
  
private:

  THnSparse *fHn; // Pointer to THnSparse object
  TFile *fFile;   // Pointer to input file
  Int_t fIpart;   // Particle type
  Bool_t fBoth;   // Particle and anti-particle type
      
  ClassDef(AliGenTHnSparse,2)
};

#endif
