#ifndef ALITRDEVENTINFO_H
#define ALITRDEVENTINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Event info for TRD performance train                                  //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

class AliESDHeader;
class AliESDRun;
class TH1D;
class AliTRDeventInfo : public TObject{
public:
  enum{
    kCentralityClasses = 5
   ,kLHCbunches = 3500
  };
  AliTRDeventInfo();
  AliTRDeventInfo(AliESDHeader *header, AliESDRun *run);
  AliTRDeventInfo(const AliTRDeventInfo &info);
  AliTRDeventInfo& operator=(const AliTRDeventInfo &info);
  virtual ~AliTRDeventInfo();
  virtual void  Delete(const Option_t *);

  AliESDHeader* GetEventHeader() const                 { return fHeader; }
  AliESDRun*    GetRunInfo() const                     { return fRun; }
  Int_t         GetCentrality() const                  { return fCentrality; }
  static Int_t  GetCentralityBin(Float_t cenPer);
  TString       GetFiredTriggerClasses();
  Int_t         GetMultiplicity() const                { return fMult; }
  static Int_t  GetMultiplicityBin(Int_t n);
  UShort_t      GetBunchFill() const;
  static void   GetListOfIsolatedBunches(TH1D *hbc, Int_t bunchSpacing=10);
  Bool_t        IsOwner() const                        { return TestBit(kOwner); }
  void          SetEventHeader(AliESDHeader *evHeader) { fHeader = evHeader; }
  void          SetRunInfo(AliESDRun *evRun)           { fRun = evRun; }
  void          SetCentrality(Float_t cent)            { fCentrality = cent>=0.?GetCentralityBin(cent):-1;}
  void          SetMultiplicity(Int_t n)               { fMult = n>=0?GetMultiplicityBin(n):-1;}
  void          SetOwner();

private:
  enum{
    kOwner = BIT(14)
  };
  static Int_t const   fgkMultBin[kCentralityClasses-1]; // multiplicity bins
  static Float_t const fgkCentBin[kCentralityClasses-1]; // centrality bins
  AliESDHeader* fHeader;      //! The ESD Header
  AliESDRun*    fRun;         //! The ESD Run Info
  Int_t         fCentrality;  //! Centrality class based on AliCentrality
  Int_t         fMult;        //! Centrality class based on AliMultiplicity

  ClassDef(AliTRDeventInfo, 2) // Event info  relevant for TRD analysis
};
#endif
