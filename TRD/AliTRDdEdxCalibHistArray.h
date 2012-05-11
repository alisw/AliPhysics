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

#ifndef ALITRDDEDXCALIBHISTARRAY_H
#define ALITRDDEDXCALIBHISTARRAY_H


#ifndef THNSPARSE_H
#include "THnBase.h"
#endif

class TObjArray;
class TCollection;

class AliTRDdEdxCalibHistArray: public TObjArray
{
 public:
  AliTRDdEdxCalibHistArray(const Bool_t kNoInv=kTRUE);
  AliTRDdEdxCalibHistArray(const AliTRDdEdxCalibHistArray &obj);
  AliTRDdEdxCalibHistArray & operator=(const AliTRDdEdxCalibHistArray &obj);
  virtual ~AliTRDdEdxCalibHistArray(){} //definition {} important for virtual
  virtual Long64_t Merge(const TCollection* list);

  static TString GetArrayName(){ return "TRDdEdxCalibHistArray"; }
  static TString GetNameAt(const Int_t iter){ return Form("TRDdEdxCalibHist%d", iter); }
  static Int_t GetIterator(const Bool_t kinvq, const Double_t mag, const Int_t charge){ return kinvq*4 + (mag>0)*2 + (charge>0); }


  THnBase * GetHist(const Bool_t kinvq, const Double_t mag, const Int_t charge){ return (THnBase*) At(GetIterator(kinvq, mag, charge)); }

 private:
  
  ClassDef(AliTRDdEdxCalibHistArray,1);
};

#endif
