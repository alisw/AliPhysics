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
#ifndef ALITRDPIDPARAMS_H
#define ALITRDPIDPARAMS_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TSortedList;

class AliTRDPIDParams : public TNamed{
  public:
    AliTRDPIDParams();
    AliTRDPIDParams(const char *name);
    virtual ~AliTRDPIDParams();
    virtual void Print(Option_t *) const;

    Bool_t GetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params) const;
    void SetThresholdParameters(Int_t ntracklets, Double_t effMin, Double_t effMax, Double_t *params);

  private:
    class AliTRDPIDThresholds : public TObject{
      public:
        AliTRDPIDThresholds();
        AliTRDPIDThresholds(Int_t nTracklets, Double_t effMin, Double_t effMax, Double_t *params = NULL);
        AliTRDPIDThresholds(Int_t nTracklets, Double_t eff, Double_t *params = NULL);
        AliTRDPIDThresholds(const AliTRDPIDThresholds &);
        AliTRDPIDThresholds &operator=(const AliTRDPIDThresholds &);
        virtual ~AliTRDPIDThresholds() {}
        
        Int_t GetNTracklets() const { return fNTracklets; }
        Double_t GetElectronEfficiency(Int_t step = 0) const { if(step == 0) return fEfficiency[0]; else return fEfficiency[1]; }
        const Double_t *GetThresholdParams() const { return fParams; }

        virtual Bool_t IsSortable() const { return kTRUE; }
        virtual Int_t Compare(const TObject *ref) const;
        virtual Bool_t IsEqual(const TObject *ref) const;

      private:
        Int_t fNTracklets;          //
        Double_t fEfficiency[2];    //
        Double_t fParams[4];        //

        ClassDef(AliTRDPIDThresholds, 1);
    };

    static const Double_t kVerySmall;

    AliTRDPIDParams(const AliTRDPIDParams &);
    AliTRDPIDParams &operator=(const AliTRDPIDParams &);
    
    TSortedList *fEntries; //

    ClassDef(AliTRDPIDParams, 1);
};
#endif

