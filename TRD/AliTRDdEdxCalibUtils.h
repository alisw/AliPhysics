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
grep " AliTRDdEdxCalibUtils::" AliTRDdEdxCalibUtils.cxx | grep "=" -v  | grep -v "[6]" | grep -v printf  |wc
grep "(" AliTRDdEdxCalibUtils.h | grep ";" | grep -v grep | grep -v ClassDef | grep -v "{" | grep -v typedef | wc
*/


#ifndef ALITRDDEDXCALIBUTILS_H
#define ALITRDDEDXCALIBUTILS_H

#ifndef TVECTORD_H
#include "TVectorD.h"
#endif

#ifndef THNSPARSE_H
#include "THnBase.h"
#endif

#ifndef TTREESTREAM_H
#include "TTreeStream.h"
#endif 

#ifndef ALITRDDEDXCALIBHISTARRAY_H
#include "AliTRDdEdxCalibHistArray.h"
#endif

class TH1D;
class TH2D;
class TObjArray;

class AliESDEvent;
class AliESDtrack;
class AliTRDcluster;
class AliTRDtrackV1;
class AliTRDseedV1;

class AliTRDdEdxCalibUtils
{
 public:

  static void SetObjArray(TObjArray * obj){fgObjArray = obj;}
  static TObjArray * GetObjArray();
  static TObjArray * GetObj(const Bool_t kinvq, const Double_t mag, const Int_t charge);
  static TObjArray* HistToObj(const THnBase *hh, Int_t run=-999, TList *lout=0x0, TTreeSRedirector *calibStream=0x0);
  static void DeleteObjArray();
  static Bool_t GenerateDefaultOCDB(const TString path="local://./");

  static AliTRDdEdxCalibHistArray * GetHistArray(){return fgHistArray;}
  static THnBase * GetHistAt(const Int_t iter);
  static void IniHistArray(TList *list, const Bool_t kNoInv);
  static Bool_t ReadHistArray(const TString filename, const TString listname);
  static void FillHist(const AliTRDtrackV1 *trdv1, const Bool_t kinvq, const Double_t mag, const Int_t charge, const Double_t scale) ;
  static void DeleteHistArray();

  static Double_t GetCalibTPCscale(const Int_t tpcncls, const Double_t tpcsig);
  static void Output(const TList *lin, Int_t run);
  
 private:
  static void FillHist(const Int_t ncls, const TVectorD *arrayQ, const TVectorD *arrayX, THnBase * hcalib, const Double_t scale);
  static void GetPHCountMeanRMS(const TH1D *hnor, TH1D *&hmean);

 
  static AliTRDdEdxCalibHistArray * fgHistArray;         //array containing 8 THnBase!
  static TObjArray * fgObjArray;                    //array containing 8 TObjArray!
 
};

#endif
