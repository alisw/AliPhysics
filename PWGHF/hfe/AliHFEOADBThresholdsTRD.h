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
// Container for TRD thresholds stored in the OADB
//
#ifndef ALIHFEOADBTHRESHOLDSTRD_H
#define ALIHFEOADBTHRESHOLDSTRD_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TSortedList;

class AliHFEOADBThresholdsTRD : public TNamed{
  public:
    AliHFEOADBThresholdsTRD();
    AliHFEOADBThresholdsTRD(const char *name);
    virtual ~AliHFEOADBThresholdsTRD();
    virtual void Print(Option_t *) const;

    Bool_t GetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params);
    void SetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params);

  private:
    class AliHFEthresholdParamsTRD : public TObject{
      public:
        AliHFEthresholdParamsTRD();
        AliHFEthresholdParamsTRD(Int_t nTracklets, Double_t eff, Double_t *params = NULL);
        AliHFEthresholdParamsTRD(const AliHFEthresholdParamsTRD &);
        AliHFEthresholdParamsTRD &operator=(const AliHFEthresholdParamsTRD &);
        virtual ~AliHFEthresholdParamsTRD() {}
        
        Int_t GetNTracklets() const { return fNTracklets; }
        Double_t GetElectronEfficiency() const { return fEfficiency; }
        const Double_t *GetThresholdParams() const { return fParams; }

        virtual Bool_t IsSortable() const { return kTRUE; }
        virtual Int_t Compare(const TObject *ref) const;
        virtual Bool_t IsEqual(const TObject *ref) const;

      private:
        Int_t fNTracklets;      // Number of tracklets
        Double_t fEfficiency;   // Efficiency level applied
        Double_t fParams[4];    // Cut parameterization

        ClassDef(AliHFEthresholdParamsTRD, 1);
    };

    static const Double_t fgkVerySmall;   // Comparison of efficiency values

    AliHFEOADBThresholdsTRD(const AliHFEOADBThresholdsTRD &ref);
    AliHFEOADBThresholdsTRD &operator=(const AliHFEOADBThresholdsTRD &ref);
    
    TSortedList *fEntries; // Container for Thresholds

    ClassDef(AliHFEOADBThresholdsTRD, 1);
};
#endif

