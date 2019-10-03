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
class AliESDtrackCuts;

class AliRsnCutTrackQuality : public AliRsnCut {
public:

   AliRsnCutTrackQuality(const char *name = "AliRsncutTrackQuality");
   AliRsnCutTrackQuality(const AliRsnCutTrackQuality &copy);
   AliRsnCutTrackQuality &operator=(const AliRsnCutTrackQuality &copy);
   virtual ~AliRsnCutTrackQuality() { }

   void      DisableAll();

   void      AddStatusFlag(ULong_t f, Bool_t on)       {if (on) fFlagsOn = fFlagsOn | f; else fFlagsOff = fFlagsOff | f;}
   void      SetStatusFlags(ULong_t f, Bool_t on)      {if (on) fFlagsOn = f; else fFlagsOff = f;}
   void      SetPtRange(Double_t a, Double_t b);
   void      SetEtaRange(Double_t a, Double_t b);

   void      SetDCARPtFormula(const char *formula)     {fDCARptFormula = formula; fDCARmaxfixed = kFALSE;}
   void      SetDCARPtFormulaMin(const char *formula)  {fDCARptFormulaMin = formula; fDCARminfixed = kFALSE;}
   void      SetDCARmax(Double_t value)                {fDCARmax = value; fDCARptFormula = ""; fDCARmaxfixed = kTRUE;}
   void      SetDCARmin(Double_t value)                {fDCARmin = value; fDCARptFormulaMin = ""; fDCARminfixed = kTRUE;}
   void      SetDCAZPtFormula(const char *formula)     {fDCAZptFormula = formula; fDCAZfixed = kFALSE;}
   void      SetDCAZmax(Double_t value)                {fDCAZmax = value; fDCAZptFormula = ""; fDCAZfixed = kTRUE;}
   void      SetDCA2D(Double_t value)                  {fDCA2D = value;}

   void      SetSPDminNClusters(Int_t value)           {fSPDminNClusters = value;}
   void      SetITSminNClusters(Int_t value)           {fITSminNClusters = value;}
   void      SetITSmaxChi2(Double_t value)             {fITSmaxChi2 = value;}

   void      SetTPCminNClusters(Int_t value)           {fTPCminNClusters = value;}
   void      SetTPCmaxChi2(Double_t value)             {fTPCmaxChi2 = value;}
   void      SetMaxChi2TPCConstrainedGlobal(Float_t max) {fCutMaxChi2TPCConstrainedVsGlobal = max; }
   void      SetTrackMaxChi2(Double_t value)           {fTrackMaxChi2 = value;}
   void      SetMinNCrossedRowsTPC(Double_t min, Bool_t useTPCCrossedRows) {fTPCminNCrossedRows=min; fIsUseCrossedRowsCut=useTPCCrossedRows;}
   void      SetMinNCrossedRowsOverFindableClsTPC(Double_t min, Bool_t useTPCCrossedRows) {fTPCminCrossedRowsOverFindableCls=min; fIsUseCrossedRowsCut=useTPCCrossedRows;}
   void      SetMinLengthActiveVolumeTPC(Double_t min, Bool_t on = kFALSE) {fCutMinLengthActiveVolumeTPC=min; fIsUseLengthActiveVolumeTPCCut=on;}
   void      SetRejectKinkDaughters(Bool_t yn = kTRUE) {fRejectKinkDaughters = yn;}

   void      SetAODTestFilterBit(Int_t value)          {fAODTestFilterBit = value;}
   void      SetCheckOnlyFilterBit(Bool_t on=kTRUE)    {fCheckOnlyFilterBit=on;}
 
   void      SetDefaults2010(Bool_t useTPCCrossedRows = kTRUE, Bool_t useDefaultKinematicCuts=kTRUE);
   void      SetDefaults2011(Bool_t useTPCCrossedRows = kTRUE, Bool_t useDefaultKinematicCuts=kTRUE);
   void      SetDefaultsHighPt2011(Bool_t useTPCCrossedRows = kTRUE, Bool_t useDefaultKinematicCuts=kTRUE);
   void      SetDefaultsTPCOnly(Bool_t useDefaultKinematicCuts=kTRUE);
   void      SetESDtrackCuts(AliESDtrackCuts *esdTrackCuts) {fESDtrackCuts = esdTrackCuts;}
   AliESDtrackCuts  *GetESDtrackCuts() {return fESDtrackCuts;}
   Double_t   GetPtRange(Bool_t max) {return fPt[max];}
   Double_t   GetEtaRange(Bool_t max) {return fEta[max];}
   
   virtual Bool_t IsSelected(TObject *obj);
   virtual void   Print(const Option_t *option = "") const;

protected:

   Bool_t      CheckESD(AliESDtrack *track);
   Bool_t      CheckAOD(AliAODTrack *track);
   const char *Binary(UInt_t number);

   ULong_t    fFlagsOn;                // status flags which must be ON (used AliESDtrack ones, connected with '|')
   ULong_t    fFlagsOff;               // status flags which must be OFF (used AliESDtrack ones, connected with '|')
   Double_t   fPt[2];                  // pt range
   Double_t   fEta[2];                 // eta range
   Bool_t     fRejectKinkDaughters;    // switch to kTRUE if daughters of kinks must be rejected

   Bool_t     fDCARmaxfixed;           // flag to switch between fixed and pt-dependent DCA cut (maximum)
   Bool_t     fDCARminfixed;           // flag to switch between fixed and pt-dependent DCA cut (minimum)
   TString    fDCARptFormula;          // expression to compute transverse DCA sigma w.r. to pt (maximum)
   TString    fDCARptFormulaMin;       // expression to compute transverse DCA sigma w.r. to pt (minimum)
   Double_t   fDCARmax;                // maximum value for transverse DCA
   Double_t   fDCARmin;                // minimum value for transverse DCA

   Bool_t     fDCAZfixed;              // flag to switch between fixed and pt-dependent DCA cut
   TString    fDCAZptFormula;          // expression to compute longitudinal DCA sigma w.r. to pt
   Double_t   fDCAZmax;                // maximum value for longitudinal DCA
   Bool_t     fDCA2D;                  // use 2D DCA to vertex

   Int_t      fSPDminNClusters;        // minimum number of required clusters in SPD
   Int_t      fITSminNClusters;        // minimum number of required clusters in ITS
   Double_t   fITSmaxChi2;             // maximum chi2 / number of clusters in ITS

   Int_t      fTPCminNClusters;        // minimum number of required clusters in TPC
   Double_t   fTPCmaxChi2;             // maximum chi2 / number of clusters in TPC
   Float_t    fCutMaxChi2TPCConstrainedVsGlobal;  // max chi2 TPC track constrained with vtx vs. global track
   
   Double_t   fTrackMaxChi2;           // maximum track chi2/NDF
   
   Bool_t     fIsUseCrossedRowsCut;     //enable cut on minimum number of TPC crossed rows
   Float_t    fTPCminNCrossedRows;     // minimum number of TPC crossed rows
   Float_t    fTPCminCrossedRowsOverFindableCls;     // minimum number of crossed rows/findable clusters
   Bool_t     fIsUseLengthActiveVolumeTPCCut;     //enable cut on minimum track lenght in TPC active volume
   Float_t    fCutMinLengthActiveVolumeTPC; // mininum length (in cm) over which the track is sampled in the active volume of the TPC (outside boundaries)

   Int_t      fAODTestFilterBit;       // test filter bit for AOD tracks
   Bool_t     fCheckOnlyFilterBit;     // check only the filter bit
   AliESDtrackCuts *fESDtrackCuts;     // pointer to AliESDtrackCuts object

   ClassDef(AliRsnCutTrackQuality, 5)
};
#endif
