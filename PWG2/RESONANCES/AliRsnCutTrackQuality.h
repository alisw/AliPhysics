//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTTRACKQUALITY_H
#define ALIRSNCUTTRACKQUALITY_H

#include <TMath.h>
#include <TString.h>

#include "AliRsnCut.h"

class AliESDtrack;
class AliAODTrack;

class AliRsnCutTrackQuality : public AliRsnCut {
public:

   AliRsnCutTrackQuality(const char *name = "AliRsncutTrackQuality");
   AliRsnCutTrackQuality(const AliRsnCutTrackQuality& copy);
   AliRsnCutTrackQuality& operator=(const AliRsnCutTrackQuality& copy);
   virtual ~AliRsnCutTrackQuality() { }

   void      DisableAll();

   void      AddStatusFlag(ULong_t f, Bool_t on)       {if (on) fFlagsOn = fFlagsOn | f; else fFlagsOff = fFlagsOff | f;}
   void      SetStatusFlags(ULong_t f, Bool_t on)      {if (on) fFlagsOn = f; else fFlagsOff = f;}
   void      SetPtRange(Double_t a, Double_t b)        {fPt[0] = TMath::Min(a, b); fPt[1] = TMath::Max(a, b);}
   void      SetEtaRange(Double_t a, Double_t b)       {fEta[0] = TMath::Min(a, b); fEta[1] = TMath::Max(a, b);}

   void      SetDCARPtFormula(const char *formula)     {fDCARptFormula = formula; fDCARfixed = kFALSE;}
   void      SetDCARmax(Double_t value)                {fDCARmax = value; fDCARptFormula = ""; fDCARfixed = kTRUE;}
   void      SetDCAZPtFormula(const char *formula)     {fDCAZptFormula = formula; fDCAZfixed = kFALSE;}
   void      SetDCAZmax(Double_t value)                {fDCAZmax = value; fDCAZptFormula = ""; fDCAZfixed = kTRUE;}

   void      SetSPDminNClusters(Int_t value)           {fSPDminNClusters = value;}
   void      SetITSminNClusters(Int_t value)           {fITSminNClusters = value;}
   void      SetITSmaxChi2(Double_t value)             {fITSmaxChi2 = value;}

   void      SetTPCminNClusters(Int_t value)           {fTPCminNClusters = value;}
   void      SetTPCmaxChi2(Double_t value)             {fTPCmaxChi2 = value;}

   void      SetRejectKinkDaughters(Bool_t yn = kTRUE) {fRejectKinkDaughters = yn;}
   
   void      SetAODTestFilterBit(Int_t value)          {fAODTestFilterBit = value;}

   virtual Bool_t IsSelected(TObject *obj);
   virtual void   Print(const Option_t *option = "") const;

protected:

   Bool_t      CheckESD(AliESDtrack *track);
   Bool_t      CheckAOD(AliAODTrack *track);
   const char* Binary(UInt_t number);

   ULong_t    fFlagsOn;                // status flags which must be ON (used AliESDtrack ones, connected with '|')
   ULong_t    fFlagsOff;               // status flags which must be OFF (used AliESDtrack ones, connected with '|')
   Double_t   fPt[2];                  // pt range
   Double_t   fEta[2];                 // eta range
   Bool_t     fRejectKinkDaughters;    // switch to kTRUE if daughters of kinks must be rejected

   Bool_t     fDCARfixed;              // flag to switch between fixed and pt-dependent DCA cut
   TString    fDCARptFormula;          // expression to compute transverse DCA sigma w.r. to pt
   Double_t   fDCARmax;                // maximum value for transverse DCA

   Bool_t     fDCAZfixed;              // flag to switch between fixed and pt-dependent DCA cut
   TString    fDCAZptFormula;          // expression to compute longitudinal DCA sigma w.r. to pt
   Double_t   fDCAZmax;                // maximum value for longitudinal DCA

   Int_t      fSPDminNClusters;        // minimum number of required clusters in SPD
   Int_t      fITSminNClusters;        // minimum number of required clusters in ITS
   Double_t   fITSmaxChi2;             // maximum chi2 / number of clusters in ITS

   Int_t      fTPCminNClusters;        // minimum number of required clusters in TPC
   Double_t   fTPCmaxChi2;             // maximum chi2 / number of clusters in TPC
   Int_t      fAODTestFilterBit;       // test filter bit for AOD tracks

   ClassDef(AliRsnCutTrackQuality, 1)
};

//__________________________________________________________________________________________________
inline const char * AliRsnCutTrackQuality::Binary(UInt_t number)
{
//
// Convert an integer in binary
//

    static char b[15];
    b[0] = '\0';

    UInt_t z;
    for (z = 512; z > 0; z >>= 1)
        strcat(b, ((number & z) == z) ? "1" : "0");

    return b;
}

#endif
