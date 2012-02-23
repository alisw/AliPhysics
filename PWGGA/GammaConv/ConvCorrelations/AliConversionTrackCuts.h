#ifndef ALICONVERSIONTRACKCUTS_H
#define ALICONVERSIONTRACKCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: (this code is mostly copied from AliRsnTrackQuality) adapted by Svein Lindal 	*

class AliAODEvent;
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisCuts.h"

using namespace std;

class AliConversionTrackCuts : public AliAnalysisCuts {
	
public:

  Bool_t IsSelected(TObject * object, AliAODEvent * event) { return AcceptTrack(dynamic_cast<AliAODTrack*>(object), event); }
  Bool_t IsSelected(TList * list) { return kFALSE; }
  Bool_t AcceptTrack(AliAODTrack * track, AliAODEvent* event);

  AliConversionTrackCuts();
  AliConversionTrackCuts(TString name, TString title);
  

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
  void      SetDefaults2010();
  
  virtual void   Print(const Option_t *option = "") const;

protected :

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


  AliConversionTrackCuts(const AliConversionTrackCuts&); // not implemented
  AliConversionTrackCuts& operator=(const AliConversionTrackCuts&); // not implemented

  ClassDef(AliConversionTrackCuts,1)


};


//__________________________________________________________________________________________________
inline void AliConversionTrackCuts::SetDefaults2010()
{
//
// Default settings for cuts used in 2010
//

   AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   SetEtaRange(-0.8, 0.8);
   //SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   SetDCAZmax(2.0);
   SetSPDminNClusters(1);
   SetITSminNClusters(0);
   SetITSmaxChi2(1E+20);
   SetTPCminNClusters(70);
   SetTPCmaxChi2(10.0);
   SetRejectKinkDaughters();
   SetAODTestFilterBit(5);
}

#endif
