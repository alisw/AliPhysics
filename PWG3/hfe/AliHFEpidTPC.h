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
#ifndef ALIHFEPIDTPC_H
#define ALIHFEPIDTPC_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class TList;
class AliESDtrack;
class AliTPCpidESD;
class AliVParticle;

class AliHFEpidTPC : public AliHFEpidBase{
  typedef enum{
    kHistTPCelectron = 0,
    kHistTPCpion = 1,
    kHistTPCmuon = 2,
      kHistTPCkaon = 3,
      kHistTPCproton = 4,
      kHistTPCothers = 5,
      kHistTPCall = 6,
      kHistTPCprobEl = 7,
      kHistTPCprobPi = 8,
      kHistTPCprobMu = 9,
      kHistTPCprobKa = 10,
      kHistTPCprobPro = 11,
      kHistTPCprobOth = 12,
      kHistTPCprobAll = 13,
      kHistTPCsuppressPi = 14,
      kHistTPCsuppressMu = 15,
      kHistTPCsuppressKa = 16,
      kHistTPCsuppressPro = 17,
      kHistTPCenhanceElPi = 18,
      kHistTPCenhanceElMu = 19,
      kHistTPCenhanceElKa = 20,
      kHistTPCenhanceElPro = 21,
      kHistTPCElprobPi = 22,
      kHistTPCElprobMu = 23,
      kHistTPCElprobKa = 24,
      kHistTPCElprobPro = 25
  } QAHist_t;
  public:
    AliHFEpidTPC(const Char_t *name);
    AliHFEpidTPC(const AliHFEpidTPC &ref);
    AliHFEpidTPC &operator=(const AliHFEpidTPC &ref);
    virtual ~AliHFEpidTPC();
    
    virtual Bool_t InitializePID();
    virtual Int_t IsSelected(AliVParticle *track);
    virtual Bool_t HasQAhistos() const { return kTRUE; };

    void AddTPCdEdxLineCrossing(Int_t species, Double_t sigma);
    void SetTPCnSigma(Short_t nSigma) { fNsigmaTPC = nSigma; };
    Double_t Likelihood(const AliESDtrack *track, Int_t species, Float_t rsig = 2.);

    Double_t Suppression(const AliESDtrack *track, Int_t species);

  protected:
    void Copy(TObject &o) const;
    void AddQAhistograms(TList *qaList);
    void FillTPChistograms(const AliESDtrack *track);
 
  private:
    Double_t fLineCrossingSigma[AliPID::kSPECIES];          // with of the exclusion point
    UChar_t fLineCrossingsEnabled;                          // Bitmap showing which line crossing is set
    Short_t fNsigmaTPC;                                     // TPC sigma band
    AliPID *fPID;                                           //! PID Object
    AliTPCpidESD *fPIDtpcESD;                               //! TPC PID object
    TList *fQAList;                                         //! QA histograms

  ClassDef(AliHFEpidTPC, 1)   // TPC Electron ID class
};
#endif
