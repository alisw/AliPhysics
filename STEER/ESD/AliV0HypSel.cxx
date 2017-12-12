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

//-------------------------------------------------------------------------
//    Implementation of the V0 selection hypothesis.
//    V0 is validated if it passes at least 1 hypotheses selection
//-------------------------------------------------------------------------


#include "AliV0HypSel.h"
#include "AliLog.h"

ClassImp(AliV0HypSel)

//________________________________________________________________________
AliV0HypSel::AliV0HypSel()
:  fM0(0.)
  ,fM1(0.)
  ,fMass(0.)
  ,fSigmaM(0.)
  ,fNSigma(0.)
  ,fMarginAdd(0.)
{}

//________________________________________________________________________
AliV0HypSel::AliV0HypSel(const AliV0HypSel& src)
  : TNamed(src)
  ,fM0(src.fM0)
  ,fM1(src.fM1)
  ,fMass(src.fMass)
  ,fSigmaM(src.fSigmaM)
  ,fNSigma(src.fNSigma)
  ,fMarginAdd(src.fMarginAdd)
{}

//________________________________________________________________________
AliV0HypSel::AliV0HypSel(const char *name, float m0,float m1, float mass, float sigma, float nsig, float margin) :
  TNamed(name,name)
  ,fM0(m0)
  ,fM1(m1)
  ,fMass(mass)
  ,fSigmaM(sigma)
  ,fNSigma(nsig)
  ,fMarginAdd(margin)
{
  Validate();
}

//________________________________________________________________________
void AliV0HypSel::Validate()
{
  // consistency check
  if (fM0<0.0005 || fM1<0.0005) {
    AliFatal("V0 decay product mass cannot be lighter than electron");
  }
  if (fMass<fM0+fM1) {
    AliFatal("V0 mass is less than sum of product masses");
  }
  if ( (fSigmaM<=0 || fNSigma<=0) && fMarginAdd<=1e-3) {
    AliFatal("No safety margin is provided");
  }
}

//________________________________________________________________________
void AliV0HypSel::Print(const Option_t *) const
{
  // print itself
  printf("%-15s | m0: %.4e m1: %.4e -> M: %.4e\nCut margin: %.1f*%.3e*(1+pT)+%.3e\n",GetName(),fM0,fM1,fMass,
	 fNSigma,fSigmaM,fMarginAdd);
}
