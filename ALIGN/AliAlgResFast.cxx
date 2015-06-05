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

#include "AliAlgResFast.h"
#include "AliAlgTrack.h"
#include "AliAlgPoint.h"
#include "AliAlgSens.h"
#include "AliLog.h"
#include <TString.h>
#include <TMath.h>
#include <stdio.h>

using namespace TMath;

ClassImp(AliAlgResFast)

//____________________________________
AliAlgResFast::AliAlgResFast() 
: fNPoints(0)
  ,fNMatSol(0)
  ,fNBook(0)
  ,fD0(0)
  ,fD1(0)
  ,fSig0(0)
  ,fSig1(0)
  ,fVolID(0)
  ,fLabel(0)
  ,fSolMat(0)
  ,fMatErr(0)
{
  // def c-tor
  for (int i=5;i--;) fTrCorr[i] = 0;
}

//________________________________________________
AliAlgResFast::~AliAlgResFast()
{
  // d-tor
  delete[] fD0;
  delete[] fD1;
  delete[] fSig0;
  delete[] fSig1;
  delete[] fVolID;
  delete[] fLabel;
  delete[] fSolMat;
  delete[] fMatErr;
}

//________________________________________________
void AliAlgResFast::Resize(Int_t np)
{
  // resize container
  if (np>fNBook) {
    delete[] fD0;
    delete[] fD1;
    delete[] fSig0;
    delete[] fSig1;
    delete[] fVolID;
    delete[] fLabel;
    delete[] fSolMat;
    delete[] fMatErr;
    //
    fNBook = 30+np;
    fD0      = new Float_t[fNBook];
    fD1      = new Float_t[fNBook];
    fSig0    = new Float_t[fNBook];
    fSig1    = new Float_t[fNBook];
    fVolID   = new Int_t[fNBook];
    fLabel   = new Int_t[fNBook];
    fSolMat  = new Float_t[fNBook*4]; // at most 4 material params per point
    fMatErr  = new Float_t[fNBook*4]; // at most 4 material params per point
    //
    memset(fD0   , 0,fNBook*sizeof(Float_t));
    memset(fD1   , 0,fNBook*sizeof(Float_t));
    memset(fSig0, 0,fNBook*sizeof(Float_t));
    memset(fSig1, 0,fNBook*sizeof(Float_t));
    memset(fVolID, 0,fNBook*sizeof(Int_t));
    memset(fLabel, 0,fNBook*sizeof(Int_t));
    memset(fSolMat, 0,4*fNBook*sizeof(Int_t));
    memset(fMatErr, 0,4*fNBook*sizeof(Int_t));
  }
  //
}

//____________________________________________
void AliAlgResFast::Clear(const Option_t *)
{
  // reset record
  fNPoints = 0;
  fNMatSol = 0;
  fTrCorr[4] = 0; // rest will be 100% overwritten
  //
}

//____________________________________________
void AliAlgResFast::Print(const Option_t */*opt*/) const
{
  // print info
}

//____________________________________________
void AliAlgResFast::SetResSigMeas(int ip, int ord, float res, float sig)
{
  // assign residual and error for measurement
  if (ord==0) {
    fD0[ip] = res;
    fSig0[ip] = sig;
  }
  else {
    fD1[ip] = res;
    fSig1[ip] = sig;
  }
}

//____________________________________________
void AliAlgResFast::SetMatCorr(int id, float res, float sig)
{
  // assign residual and error for material correction
  fSolMat[id] = res;
  fMatErr[id] = sig;
}

//____________________________________________
void AliAlgResFast::SetLabel(int ip, Int_t lab, Int_t vol)
{
  // set label/volid of measured volume
  fVolID[ip] = vol;
  fLabel[ip] = lab;
}
