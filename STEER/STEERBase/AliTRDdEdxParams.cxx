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

#include "AliLog.h"
#include "AliTRDdEdxParams.h"

ClassImp(AliTRDdEdxParams);

AliTRDdEdxParams::AliTRDdEdxParams(const TString name, const TString title): TNamed(name,title)
{
  //
  //constructor
  //
}

void AliTRDdEdxParams::CheckType(const Int_t itype, const TString tag) const 
{
  //
  //check if itype is in range
  //

  if(itype<0 || itype>=MAXNPAR){
    AliError(Form("AliTRDdEdxParams::CheckType %s itype out of range %d\n", tag.Data(), itype));
  }
}

const TVectorF& AliTRDdEdxParams::GetParameter(const TVectorF par[], const Int_t itype)const
{
  //
  //return parameter for particle itype from par[]
  //

  CheckType(itype, "GetParameter");

  return par[itype];
}

void AliTRDdEdxParams::SetParameter(TVectorF par[], const Int_t itype, const Int_t npar, const Float_t vals[])
{
  //
  //set parameter, vals of dimension npar, for particle itype
  //

  CheckType(itype, "SetParameter");

  TVectorF p2(npar, vals);

  par[itype].ResizeTo(p2);
  par[itype] = p2;
}

void AliTRDdEdxParams::Print(Option_t* option) const
{
  //
  //print all members
  //

  TObject::Print(option);

  printf("\n======================= Mean ========================\n");
  for(Int_t ii=0; ii<MAXNPAR; ii++){
    printf("%d: Nrows() %d\n",ii, fMeanPar[ii].GetNrows());
    if(fMeanPar[ii].GetNrows()) fMeanPar[ii].Print();
  }

  printf("\n======================= Sigma ========================\n");

  for(Int_t ii=0; ii<MAXNPAR; ii++){
    printf("%d: Nrows() %d\n",ii, fSigmaPar[ii].GetNrows());
    if(fSigmaPar[ii].GetNrows()) fSigmaPar[ii].Print();
  }
  printf("AliTRDdEdxParams::Print done.\n\n");
}
