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
// Authors
//  Jan Wagner j.wagner@cern.ch
//
#include "AliHFEparamBag.h"
#include <TCollection.h>


ClassImp(AliHFEparamBag)


AliHFEparamBag::AliHFEparamBag():
    TObject(),
    useMC(),
    isAOD(),
    isBeauty(),
    appendix(),
    TPCcl(),
    TPCclPID(), 
    ITScl(),
    DCAxy(), 
    DCAz(), 
    TOFs(), 
    TOFmis(),
    itshitpixel(),
    spdcheck(),
    icent(),
    etami(),
    etama(),
    phimi(),
    phima(),
    assETAm(),
    assETAp(),
    assMinPt(),
    assITS(),
    assTPCcl(),
    assTPCPIDcl(),
    assDCAr(),
    assDCAz(),
    assITSpid(),
    assTOFs(),
    useCat1Tracks(),
    useCat2Tracks(),
    weightlevelback(),
    WhichWei(),
    etadalwei(),
    nonPhotonicElectronBeauty(),
    nonHFEsys(),
    ipCharge(),
    ipOpp(),
    mcQADebugTree(),
    RefMulti(),
    fName()
{
    // Default constructor
}

//__________________________________________________________________
AliHFEparamBag::AliHFEparamBag(const Char_t *name):
    TObject(),
    useMC(),
    isAOD(),
    isBeauty(),
    appendix(),
    TPCcl(),
    TPCclPID(), 
    ITScl(),
    DCAxy(), 
    DCAz(), 
    TOFs(), 
    TOFmis(),
    itshitpixel(),
    spdcheck(),
    icent(),
    etami(),
    etama(),
    phimi(),
    phima(),
    assETAm(),
    assETAp(),
    assMinPt(),
    assITS(),
    assTPCcl(),
    assTPCPIDcl(),
    assDCAr(),
    assDCAz(),
    assITSpid(),
    assTOFs(),
    useCat1Tracks(),
    useCat2Tracks(),
    weightlevelback(),
    WhichWei(),
    etadalwei(),
    nonPhotonicElectronBeauty(),
    nonHFEsys(),
    ipCharge(),
    ipOpp(),
    mcQADebugTree(),
    RefMulti(),
    fName(name)
{
}

//__________________________________________________________________
AliHFEparamBag::AliHFEparamBag(const AliHFEparamBag &ref):
    TObject(ref),
    useMC(ref.useMC),
    isAOD(ref.isAOD),
    isBeauty(ref.isBeauty),
    TPCcl(ref.TPCcl),
    TPCclPID(ref.TPCclPID), 
    ITScl(ref.ITScl),
    DCAxy(ref.DCAxy), 
    DCAz(ref.DCAz), 
    TOFs(ref.TOFs), 
    TOFmis(ref.TOFmis),
    itshitpixel(ref.itshitpixel),
    spdcheck(ref.spdcheck),
    icent(ref.icent),
    etami(ref.etami),
    etama(ref.etama),
    phimi(ref.phimi),
    phima(ref.phima),
    assETAm(ref.assETAm),
    assETAp(ref.assETAp),
    assMinPt(ref.assMinPt),
    assITS(ref.assITS),
    assTPCcl(ref.assTPCcl),
    assTPCPIDcl(ref.assTPCPIDcl),
    assDCAr(ref.assDCAr),
    assDCAz(ref.assDCAz),
    assITSpid(ref.assITSpid),
    assTOFs(ref.assTOFs),
    useCat1Tracks(ref.useCat1Tracks),
    useCat2Tracks(ref.useCat2Tracks),
    weightlevelback(ref.weightlevelback),
    WhichWei(ref.WhichWei),
    etadalwei(ref.etadalwei),
    nonPhotonicElectronBeauty(ref.nonPhotonicElectronBeauty),
    nonHFEsys(ref.nonHFEsys),
    ipCharge(ref.ipCharge),
    ipOpp(ref.ipOpp),
    mcQADebugTree(ref.mcQADebugTree),
    RefMulti(ref.RefMulti),
    fName(ref.fName)
{
    appendix = new TString(*ref.appendix);
    memcpy(tpcdEdxcutlow,ref.tpcdEdxcutlow,sizeof(Double_t)*12);
    memcpy(tpcdEdxcuthigh,ref.tpcdEdxcuthigh,sizeof(Double_t)*12);
    memcpy(assTPCSminus,ref.assTPCSminus,sizeof(Double_t)*12);
    memcpy(assTPCSplus,ref.assTPCSplus,sizeof(Double_t)*12);
    memcpy(ipPar,ref.ipPar,sizeof(Double_t)*3);
}

//__________________________________________________________________
AliHFEparamBag &AliHFEparamBag::operator=(const AliHFEparamBag &ref){
    //
    // Assignment operator
    //
    TObject::operator=(ref);
    fName = (ref.fName);
    return *this;
}

//__________________________________________________________________
AliHFEparamBag::~AliHFEparamBag(){
    //
    // Destructor
    //
}
const char *AliHFEparamBag::GetName() const
{
    // Return name 
    // if no name, return the class name.

    if (fName.Length() > 0) return fName.Data();
    return ClassName();
}


//__________________________________________________________________
Long64_t AliHFEparamBag::Merge(TCollection *coll){
  //
  // Merge paramBag 
  // idea: skip merging and keep only one cycle
  //
  if(!coll)
    return 0;
  if(coll->IsEmpty())
    return 1;

  TIter iter(coll);
  TObject *o = NULL;
  Long64_t count = 0;
  while((o = iter())){
    AliHFEparamBag *bag = dynamic_cast<AliHFEparamBag *>(o);
    if(!bag) continue;

    count++;
  }
  return count + 1;
}
