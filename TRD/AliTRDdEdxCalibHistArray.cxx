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
// xx
// xx
// xx
// xx
// xx
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//  
//
#include "THnBase.h"
#include "THn.h"
#include "TCollection.h"

#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdEdxCalibHistArray.h"

ClassImp(AliTRDdEdxCalibHistArray);

AliTRDdEdxCalibHistArray::AliTRDdEdxCalibHistArray(const Bool_t kNoInv): 
  TObjArray(kNoInv ? 4: 8)
{
  //
  //constructor
  //
  SetName(GetArrayName());
  SetOwner(kTRUE);

  const Int_t    nbin[2]={AliTRDdEdxBaseUtils::NTRDtimebin(), 100};
  const Double_t xmin[2]={0,       0.1};
  const Double_t xmax[2]={nbin[0], 20};
  const TString aname[2]={"globalTimeBin", "trdqovertpc"};
  const TString atitle[2]={"det * AliTRDseedV1::kNtb + itb", "TRD-Cluster-Q / TPC-Signal"};

  for(Int_t iter=0; iter<GetSize(); iter++){
    THnBase *hi = new THnF(GetNameAt(iter), "", 2, nbin, xmin, xmax);
    for(Int_t iaxis=0; iaxis<2; iaxis++){
      TAxis *xi = hi->GetAxis(iaxis);
      xi->SetName(aname[iaxis]);
      xi->SetTitle(atitle[iaxis]);
      AliTRDdEdxBaseUtils::BinLogX(xi);
    }
    AddAt(hi, iter);
  }
}

AliTRDdEdxCalibHistArray::AliTRDdEdxCalibHistArray(const AliTRDdEdxCalibHistArray &obj): 
  TObjArray(obj)
{
  //
  //copy constructor
  //
}

AliTRDdEdxCalibHistArray & AliTRDdEdxCalibHistArray::operator=(const AliTRDdEdxCalibHistArray &obj)
{
  //
  //assignment operator
  //

  if(&obj == this) return *this;

  TObjArray::operator=(obj);

  return *this;
}

Long64_t AliTRDdEdxCalibHistArray::Merge(const TCollection* list) 
{
  //
  // Merge list of objects (needed by PROOF)
  //

  if(!list)
    return 0;
  
  if(list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
    {
      AliTRDdEdxCalibHistArray * entry = dynamic_cast<AliTRDdEdxCalibHistArray*>(obj);
      if (entry == 0) continue; 
      
      if(GetSize()!= entry->GetSize()){
        printf("AliTRDdEdxCalibHistArray::Merge GetSize()!= entry->GetSize() %d %d\n", GetSize(), entry->GetSize()); exit(1);
      }

      for(Int_t ii=0; ii<GetSize(); ii++){
        THnBase *h0 = (THnBase*) At(ii);
        THnBase *h1 = (THnBase*) entry->At(ii);
        h0->Add(h1);
      }
      
      count++;
    }
  
  return count;

}
