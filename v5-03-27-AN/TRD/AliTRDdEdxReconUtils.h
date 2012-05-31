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
//
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//
/*
grep " AliTRDdEdxReconUtils::" AliTRDdEdxReconUtils.cxx | grep "=" -v  | grep -v "[6]" | grep -v printf  |wc
grep "(" AliTRDdEdxReconUtils.h | grep ";" | grep -v grep | grep -v ClassDef | grep -v "{" | grep -v typedef | wc
*/


#ifndef ALITRDDEDXRECONUTILS_H
#define ALITRDDEDXRECONUTILS_H

#ifndef TVECTORD_H
#include "TVectorD.h"
#endif

#ifndef THNSPARSE_H
#include "THnBase.h"
#endif

#ifndef TTREESTREAM_H
#include "TTreeStream.h"
#endif 

class TH1D;
class TH2D;
class TObjArray;

class AliESDEvent;
class AliESDtrack;
class AliTRDcluster;
class AliTRDtrackV1;
class AliTRDseedV1;

class AliTRDdEdxReconUtils
{
 public:


  static Int_t ApplyCalib(const Int_t nc0, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj);
  static Double_t ToyCook(const Bool_t kinvq, Int_t &ncluster, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj=0x0);
  static Double_t CombineddEdx(const Bool_t kinvq, Int_t &concls, TVectorD *coarrayQ, TVectorD *coarrayX, const Int_t tpcncls, const TVectorD *tpcarrayQ, const TVectorD *tpcarrayX, const Int_t trdncls, const TVectorD *trdarrayQ, const TVectorD *trdarrayX);
  static Int_t GetArrayClusterQ(const Bool_t kinvq, TVectorD *arrayQ, TVectorD *arrayX, const AliTRDtrackV1 *trdtrack, Int_t timeBin0=-1, Int_t timeBin1=1000, Int_t tstep=1);
  static Int_t UpdateArrayX(const Int_t ncls, TVectorD* arrayX);

 private:
  static Double_t GetAngularCorrection(const AliTRDseedV1 *seed);
  static Double_t GetPadGain(const Int_t det, const Int_t icol, const Int_t irow);
  static Double_t GetRNDClusterQ(AliTRDcluster *cl);
  static Double_t GetClusterQ(const Bool_t kinvq, const AliTRDseedV1 * seed, const Int_t itb);

};

#endif
