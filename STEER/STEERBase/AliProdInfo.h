#ifndef ALIPRODINFO_H
#define ALIPRODINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------//
//        Base class for parsing alirootVersion                  //
//        TNamed object stored in ESD/AOD data                   //
//        and setup easy getters for end user                    //
//                                                               //
//   Origin: Pietro Antonioli, INFN-BO, pietro.antonioli@cern.ch //
//                                                               //
//---------------------------------------------------------------//
#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>
#include <TObjString.h>

class AliProdInfo : public TNamed {
public:
  AliProdInfo();
  AliProdInfo(const TString& name, const TString& title);
  AliProdInfo(TList *userInfo);
  virtual ~AliProdInfo();

  void Init(TList *userInfo);
  void List() const;
  TString GetLHCPeriod() const {return fPeriod;}
  TString GetAlirootVersion() const {return fAlirootVersion;}
  Int_t GetAlirootSvnVersion() const {return fAlirootSvnVersion;}
  TString GetRootVersion() const {return fRootVersion;}
  Int_t GetRootSvnVersion() const {return fRootSvnVersion;}
  Int_t GetRecoPass() const {return fRecoPass;}
  Bool_t IsMC() const {return fMcFlag;}

protected:
  void ParseProdInfo(TNamed *uList);

private:
  AliProdInfo(const AliProdInfo&);
  AliProdInfo &operator=(const AliProdInfo&);

  TString fPeriod;            // LHC period
  TString fAlirootVersion;    // aliroot version used producing data
  Int_t fAlirootSvnVersion;   // aliroot svn numbering
  TString fRootVersion;       // root version used producing data
  Int_t fRootSvnVersion;      // root svn numbering
  Bool_t fMcFlag;             // MC data: kTrue ; raw data: kFalse
  Int_t fRecoPass;            // Reconstruction pass
  ClassDef(AliProdInfo, 1);   // Combined PID using priors
};

#endif
