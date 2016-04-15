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
    fName("default"),
    useMC(kFALSE),
    isAOD(kFALSE),
    isBeauty(kFALSE),
    appendix(""),
    TPCcl(0),
    TPCclPID(0), 
    ITScl(0),
    DCAxy(0.0), 
    DCAz(0.0), 
    TOFs(0.0), 
    TOFmis(0),
    itshitpixel(0),
    spdcheck(0),
    icent(0),
    etami(0.0),
    etama(0.0),
    phimi(0.0),
    phima(0.0),
    assETAm(0.0),
    assETAp(0.0),
    assMinPt(0.0),
    assITS(0),
    assTPCcl(0),
    assTPCPIDcl(0),
    assDCAr(0.0),
    assDCAz(0.0),
    assITSpid(0.0),
    assTOFs(0.0),
    useCat1Tracks(kFALSE),
    useCat2Tracks(kFALSE),
    weightlevelback(-1),
    WhichWei(0),
    etadalwei(1),
    nonPhotonicElectronBeauty(kFALSE),
    nonHFEsys(kFALSE),
    ipCharge(kFALSE),
    ipOpp(kFALSE),
    mcQADebugTree(kFALSE),
    RefMulti(0.0),
    MultiSystem(0)
{
    memset(tpcdEdxcutlow,0,sizeof(Double_t)*12);
    memset(tpcdEdxcuthigh,0,sizeof(Double_t)*12);
    memset(assTPCSminus,0,sizeof(Double_t)*12);
    memset(assTPCSplus,0,sizeof(Double_t)*12);
    memset(ipPar,0,sizeof(Double_t)*3);
    // Default constructor
}

//__________________________________________________________________
AliHFEparamBag::AliHFEparamBag(const char *name):
    TObject(),
    fName(name),
    useMC(kFALSE),
    isAOD(kFALSE),
    isBeauty(kFALSE),
    appendix(""),
    TPCcl(0),
    TPCclPID(0), 
    ITScl(0),
    DCAxy(0.0), 
    DCAz(0.0), 
    TOFs(0.0), 
    TOFmis(0),
    itshitpixel(0),
    spdcheck(0),
    icent(0),
    etami(0.0),
    etama(0.0),
    phimi(0.0),
    phima(0.0),
    assETAm(0.0),
    assETAp(0.0),
    assMinPt(0.0),
    assITS(0),
    assTPCcl(0),
    assTPCPIDcl(0),
    assDCAr(0.0),
    assDCAz(0.0),
    assITSpid(0.0),
    assTOFs(0.0),
    useCat1Tracks(kFALSE),
    useCat2Tracks(kFALSE),
    weightlevelback(-1),
    WhichWei(0),
    etadalwei(1),
    nonPhotonicElectronBeauty(kFALSE),
    nonHFEsys(kFALSE),
    ipCharge(kFALSE),
    ipOpp(kFALSE),
    mcQADebugTree(kFALSE),
    RefMulti(0.0),
    MultiSystem(0)
{
    memset(tpcdEdxcutlow,0,sizeof(Double_t)*12);
    memset(tpcdEdxcuthigh,0,sizeof(Double_t)*12);
    memset(assTPCSminus,0,sizeof(Double_t)*12);
    memset(assTPCSplus,0,sizeof(Double_t)*12);
    memset(ipPar,0,sizeof(Double_t)*3);
}

//__________________________________________________________________
AliHFEparamBag::AliHFEparamBag(const AliHFEparamBag &ref):
    TObject(ref),
    fName(ref.fName),
    useMC(ref.useMC),
    isAOD(ref.isAOD),
    isBeauty(ref.isBeauty),
    appendix(ref.appendix),
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
    MultiSystem(ref.MultiSystem)
{
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
