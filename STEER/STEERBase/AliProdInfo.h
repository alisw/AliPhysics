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
#include <AliOADBContainer.h>

class AliProdInfo : public TNamed {
public:
  enum {kParsedBit=BIT(14)};
  enum ETagType {
    kAliroot,
    kRoot,
    kOutDir,
    kPass,
    kProdType,
    kProdTag,
    kRunNumber,
    kAnchorRun,
    kAnchorProduction,
    kAnchorPassName,
    kNTags};

  AliProdInfo();
  AliProdInfo(const TString& name, const TString& title);
  AliProdInfo(TList *userInfo);
  virtual ~AliProdInfo();

  void Init(TList *userInfo);
  void List() const;
  //
  const TString& GetTag(ETagType tag) { return fTags[tag]; }
  //
  const TString& GetLHCPeriod() const {return fAnchorProduction;}
  const TString& GetAlirootVersion() const {return fAlirootVersion;}
  Int_t GetAlirootSvnVersion() const {return fAlirootSvnVersion;}
  const TString& GetRootVersion() const {return fRootVersion;}
  Int_t GetRootSvnVersion() const {return fRootSvnVersion;}
  Int_t GetRecoPass() const {return fRecoPass;}
  const TString& GetRecoPassName() const { return fRecoPassName; }
  const TString& GetAnchorProduction() const { return fAnchorProduction; }
  const TString& GetAnchorPassName() const { return fAnchorPassName; }

  void SetPassNameMapFileName(const TString& name) { fPassNameMapFileName = name; }
  const TString& GetPassNameMapFileName() const { return fPassNameMapFileName; }

  Bool_t IsMC() const {return fMcFlag;}
  //
  Bool_t HasLPMPass() const { return !fTags[kPass].IsNull() && fTags[kPass].IsDigit(); }
  //
  Bool_t IsParsed()                 const {return TestBit(kParsedBit);}
  void   SetParsed(Bool_t v=kTRUE)        {SetBit(kParsedBit);}
  //
protected:
  void ParseProdInfo(TNamed *uList);

private:
  AliProdInfo(const AliProdInfo&);
  AliProdInfo &operator=(const AliProdInfo&);
  //
  Bool_t  fMcFlag;                // MC data: kTrue ; raw data: kFalse
  Int_t   fAlirootSvnVersion;     // aliroot svn numbering
  Int_t   fRootSvnVersion;        // root svn numbering
  Int_t   fRecoPass;              // Reconstruction pass
  Int_t   fRunNumber;             // Run number
  Int_t   fAnchorRun;             // MC production anchor run number
  //
  TString fTags[kNTags];          // Array with tag values
  TString fProductionTag;         // production tag
  TString fAlirootVersion;        // aliroot version used producing data
  TString fRootVersion;           // root version used producing data
  TString fRecoPassName;          // Full name of the reconstruction pass, deduced from the output file structure
  TString fAnchorProduction;      // Anchor production LHC period
  TString fAnchorPassName;        // Full name of the anchored pass name

  TString fPassNameMapFileName;   // name of the OADB file containing the maps with production pass matching
  AliOADBContainer fPassNameMap;  //! container to map production to a pass name

  void UpdateMCInfoFromOADB();
  void CheckForSpecialPassName();
  void SetNumericRecoPass(const TString& name);
  //
  ClassDef(AliProdInfo, 5);     // Combined PID using priors
};

#endif
