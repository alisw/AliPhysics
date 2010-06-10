/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// class used to extract and store info of MC particles
//
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include "AliMuonInfoStoreMC.h"
#include "AliDimuInfoStoreRD.h"
#include "AliDimuInfoStoreMC.h"

ClassImp(AliDimuInfoStoreMC)

const TString AliDimuInfoStoreMC::fgkStdBranchName("DimuMC");
const Int_t   AliDimuInfoStoreMC::fgkSourcesN = 6;

//-----------------------------------------------------------------------------
AliDimuInfoStoreMC::AliDimuInfoStoreMC() :
AliDimuInfoStoreRD(),
fIsFull(kFALSE),
fLorentzP(),
fSource(-1)
{
  //
  // default constructor
  //
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreMC::AliDimuInfoStoreMC(AliMuonInfoStoreMC* const trk0, AliMuonInfoStoreMC* const trk1, Bool_t full) :
AliDimuInfoStoreRD(),
fIsFull(full),
fLorentzP(),
fSource(-1)
{
  //
  // default constructor
  //
  fMuonRef[0] = trk0;
  fMuonRef[1] = trk1;
  fLorentzP = trk0->LorentzP() + trk1->LorentzP();
  AliDimuInfoStoreRD::FillDimuInfo();
  if (fIsFull) this->FindDimuonSourceFull();
  else this->FindDimuonSourceFast();
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreMC::AliDimuInfoStoreMC(const AliDimuInfoStoreMC &src) :
AliDimuInfoStoreRD(src),
fIsFull(src.fIsFull),
fLorentzP(src.fLorentzP),
fSource(src.fSource)
{
  //
  // copy constructor
  //
}

//-----------------------------------------------------------------------------
AliDimuInfoStoreMC& AliDimuInfoStoreMC::operator=(const AliDimuInfoStoreMC &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;

  fIsFull   = src.fIsFull;
  fLorentzP = src.fLorentzP;
  fSource   = src.fSource;

  return *this;
}


//-----------------------------------------------------------------------------
AliDimuInfoStoreMC::~AliDimuInfoStoreMC()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliDimuInfoStoreMC::FindDimuonSourceFast()
{
  // find corr relation of two particles (fast for p-p)

  AliMuonInfoStoreMC *trk0 = (AliMuonInfoStoreMC*)fMuonRef[0].GetObject();
  Int_t src0 = trk0->Source();
  if (src0<0 || src0==4 || src0==3) {
    fSource=5; return;
  }

  AliMuonInfoStoreMC *trk1 = (AliMuonInfoStoreMC*)fMuonRef[1].GetObject();
  Int_t src1 = trk1->Source();
  if (src1<0 || src1==4 || src1==3) {
    fSource=5; return;
  }

  // Drell-Yan is expected very small at LHC, we ingore it
  Int_t np0 = trk0->ParentsN() - 1;
  if (np0<0) {
    fSource=5; return;
  }
  Int_t np1 = trk1->ParentsN() - 1;
  if (np1<0) {
    fSource=5; return;
  }

  if (trk0->IsMotherAResonance(np0) &&
      trk1->IsMotherAResonance(np1) &&
      trk0->ParentIndex(np0)==trk1->ParentIndex(np1)) {
    fSource=4; return;
  }

  if (src0==0 && src1==0) {
    if ((trk0->ParentIndex(0))==(trk1->ParentIndex(0)))
      fSource = 1;
    else
      fSource = 0;
    return;
  }

  if (src0==1 && src1==1) {
    if ((trk0->ParentIndex(0))==(trk1->ParentIndex(0)))
      fSource = 3;
    else
      fSource = 2;
    return;
  }

  fSource = 5;
  return;
}
