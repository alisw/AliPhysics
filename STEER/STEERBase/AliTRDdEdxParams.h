/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//

#ifndef ALITRDDEDXPARAMS_H
#define ALITRDDEDXPARAMS_H

#include "TNamed.h"
#include "TVectorT.h"
#include "TString.h"

//maximum number of paticle types
#define MAXNPAR 10

class AliTRDdEdxParams: public TNamed
{
 public:
  AliTRDdEdxParams(const TString name="name", const TString title="title");
  void Print(Option_t* option = "") const;

  const TVectorF GetMeanParameter(const Int_t itype) const { return GetParameter(fMeanPar, itype);}
  const TVectorF GetSigmaParameter(const Int_t itype) const { return GetParameter(fSigmaPar, itype);}

  void SetMeanParameter(const Int_t itype, const Int_t npar, const Float_t vals[]){ SetParameter(fMeanPar, itype, npar, vals); }
  void SetSigmaParameter(const Int_t itype, const Int_t npar, const Float_t vals[]){ SetParameter(fSigmaPar, itype, npar, vals); } 

 private:
  const TVectorF GetParameter(const TVectorF par[], const Int_t itype) const;
  void SetParameter(TVectorF par[], const Int_t itype, const Int_t npar, const Float_t vals[]);

  TVectorF fMeanPar[MAXNPAR];
  TVectorF fSigmaPar[MAXNPAR];

  void CheckType(const Int_t itype, const TString tag) const;

  ClassDef(AliTRDdEdxParams,1);
};

#endif
