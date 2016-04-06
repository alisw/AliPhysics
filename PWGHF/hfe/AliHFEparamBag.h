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
#ifndef ALIHFEPARAMBAG_H
#define ALIHFEPARAMBAG_H

#include <TObject.h>
#include <TString.h>

class TCollection;

class AliHFEparamBag : public TObject{ 

protected:
    TString fName;  // Variable Name

public:
    AliHFEparamBag();
    AliHFEparamBag(const char *name);
    AliHFEparamBag(const AliHFEparamBag &ref);
    AliHFEparamBag &operator=(const AliHFEparamBag &ref);
    ~AliHFEparamBag();

    Bool_t useMC;
    Bool_t isAOD;
    Bool_t isBeauty;
    TString appendix;
    Int_t TPCcl;
    Int_t TPCclPID; 
    Int_t ITScl;
    Double_t DCAxy; 
    Double_t DCAz; 
    Double_t tpcdEdxcutlow[12]; 
    Double_t tpcdEdxcuthigh[12];
    Double_t TOFs; 
    Int_t TOFmis;
    Int_t itshitpixel;
    Int_t spdcheck;
    Int_t icent;
    Double_t etami;
    Double_t etama;
    Double_t phimi;
    Double_t phima;
    Double_t assETAm;
    Double_t assETAp;
    Double_t assMinPt;
    Int_t assITS;
    Int_t assTPCcl;
    Int_t assTPCPIDcl;
    Double_t assDCAr;
    Double_t assDCAz;
    Double_t assTPCSminus[12];
    Double_t assTPCSplus[12];
    Double_t assITSpid;
    Double_t assTOFs;
    Bool_t useCat1Tracks;
    Bool_t useCat2Tracks;
    Int_t weightlevelback;
    Int_t WhichWei;
    Double_t etadalwei;
    Bool_t nonPhotonicElectronBeauty;
    Bool_t nonHFEsys;
    Bool_t ipCharge;
    Bool_t ipOpp;
    Double_t ipPar[3];
    Bool_t mcQADebugTree;
    Double_t RefMulti;
    Int_t MultiSystem;

    void SetName(const char *name) { fName = name; };
    const char *GetName() const;
    virtual Long64_t Merge(TCollection *coll);



    ClassDef(AliHFEparamBag, 4)
};
#endif
