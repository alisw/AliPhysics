#ifndef ALICONVERSIONTRACKCUTS_H
#define ALICONVERSIONTRACKCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: (this code is mostly copied from AliRsnTrackQuality) adapted by Svein Lindal 	*

class TH2F;
class AliESDtrackCuts;
class TList;
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisCuts.h"

using namespace std;

class AliConversionTrackCuts : public AliAnalysisCuts {
	
public:
  
  enum CTCuts_t {
	kPreCut = -1,
	kCutNcls,
	kCutNclsFrac,
	kCutNDF,
	kCutKinc,
	kCutDCAZ,
	kCutDCAXY,
	kCutTPCRefit,
	kNCuts
  };

//  enum trackVal {
//	kPt = 0, 
//	kPhi,
//	kEta,
//	kNVar
//  };
  static const char * fgkCutNames[kNCuts];

  Bool_t IsSelected(TObject * object ); 
  Bool_t IsSelected(TList * /*list*/) { return kFALSE; }
  Bool_t AcceptTrack(AliAODTrack * track);
  Bool_t AcceptTrack(AliESDtrack * track);
  void FillDCAHist(Float_t dcaz, Float_t dcaxy, AliVTrack * track);
  AliConversionTrackCuts();
  AliConversionTrackCuts(TString name, TString title);
  ~AliConversionTrackCuts();

  void      AddStatusFlag(ULong_t f, Bool_t on)       {if (on) fFlagsOn = fFlagsOn | f; else fFlagsOff = fFlagsOff | f;}
  void      SetStatusFlags(ULong_t f, Bool_t on)      {if (on) fFlagsOn = f; else fFlagsOff = f;}
  
  void      SetDCARPtFormula(const char *formula)     {fDCARptFormula = formula; fDCARfixed = kFALSE;}
  void      SetDCARmax(Double_t value)                {fDCARmax = value; fDCARptFormula = ""; fDCARfixed = kTRUE;}
  void      SetDCAZPtFormula(const char *formula)     {fDCAZptFormula = formula; fDCAZfixed = kFALSE;}
  void      SetDCAZmax(Double_t value)                {fDCAZmax = value; fDCAZptFormula = ""; fDCAZfixed = kTRUE;}
  void      SetDCAXYmax(Double_t value)                {fDCAXYmax = value*value;}
  
  void      SetSPDminNClusters(Int_t value)           {fSPDminNClusters = value;}
  void      SetITSminNClusters(Int_t value)           {fITSminNClusters = value;}
  void      SetITSmaxChi2(Double_t value)             {fITSmaxChi2 = value;}
  
  void      SetTPCminNClusters(Int_t value)           {fTPCminNClusters = value;}
  void      SetTPCCFoundClusters(Double_t value)           {fTPCClusOverFindable = value;}
  void      SetTPCmaxChi2(Double_t value)             {fTPCmaxChi2 = value;}
  
  void      SetRejectKinkDaughters(Bool_t yn = kTRUE) {fRejectKinkDaughters = yn;}  
  void      SetAODTestFilterBit(Int_t value)          {fAODTestFilterBit = value;}

  void      SetRequireTPCRefit(Bool_t require)        { fRequireTPCRefit = require; }
  void      SetDefaults2010();
  
  TList * CreateHistograms();
  void FillHistograms(Int_t cutIndex, AliVTrack * track, Bool_t passed);
  virtual void   Print(const Option_t *option = "") const;

protected :

   ULong_t    fFlagsOn;                // status flags which must be ON (used AliESDtrack ones, connected with '|')
   ULong_t    fFlagsOff;               // status flags which must be OFF (used AliESDtrack ones, connected with '|')
   Bool_t     fRejectKinkDaughters;    // switch to kTRUE if daughters of kinks must be rejected

   Bool_t     fDCARfixed;              // flag to switch between fixed and pt-dependent DCA cut
   TString    fDCARptFormula;          // expression to compute transverse DCA sigma w.r. to pt
   Double_t   fDCARmax;                // maximum value for transverse DCA

   Bool_t     fDCAZfixed;              // flag to switch between fixed and pt-dependent DCA cut
   TString    fDCAZptFormula;          // expression to compute longitudinal DCA sigma w.r. to pt
   Double_t   fDCAZmax;                // maximum value for longitudinal DCA

  Double_t fDCAXYmax;                  // maximum xy value for dca

   Int_t      fSPDminNClusters;        // minimum number of required clusters in SPD
   Int_t      fITSminNClusters;        // minimum number of required clusters in ITS
   Double_t   fITSmaxChi2;             // maximum chi2 / number of clusters in ITS

   Int_t      fTPCminNClusters;        // minimum number of required clusters in TPC
   Double_t   fTPCClusOverFindable;        // minimum number of required clusters in TPC
   Double_t   fTPCmaxChi2;             // maximum chi2 / number of clusters in TPC
   Int_t      fAODTestFilterBit;       // test filter bit for AOD tracks
  Bool_t      fRequireTPCRefit;        // Require TPC refit

  AliESDtrackCuts * fESDCuts;

  TH2F * fhPhi;
  TH2F * fhPt;
  TH2F * fhPhiPt;
  TH2F * fhdcaxyPt;
  TH2F * fhdcazPt;
  TH2F * fhdca;
  TH2F * fhnclpt;
  TH2F * fhnclsfpt;
  
 // TAxis fCutAxis;
 // TAxisArr fCutVarAxes[kNCuts];
 // TAxisArr fVarAxes[kNVar];

  TList * fHistograms;

  AliConversionTrackCuts(const AliConversionTrackCuts&); // not implemented
  AliConversionTrackCuts& operator=(const AliConversionTrackCuts&); // not implemented

  ClassDef(AliConversionTrackCuts,2)


};


//__________________________________________________________________________________________________
inline void AliConversionTrackCuts::SetDefaults2010()
{
//
// Default settings for cuts used in 2010
//

   SetDCAZmax(3.0);
   SetDCAXYmax(2.5);
   SetTPCminNClusters(70);
   SetTPCmaxChi2(12.0);
   SetRejectKinkDaughters();
}

#endif
