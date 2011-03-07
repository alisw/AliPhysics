//
// *** Class AliRsnCutESD2010 ***
//
// This class implements all cuts which have to be used for the 2010 runs
// for phi and generic resonance analysis.
// It contains an AliESDtrackCuts object for track quality selection
// and some criteria for particle identification with ITS, TPC and TOF.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTESD2010_H
#define ALIRSNCUTESD2010_H

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliRsnCut.h"

class AliESDpid;

class AliRsnCutESD2010 : public AliRsnCut {
public:

   AliRsnCutESD2010(const char *name = "cutESD2010", Bool_t isMC = kFALSE);
   AliRsnCutESD2010(const AliRsnCutESD2010& copy);
   AliRsnCutESD2010& operator=(const AliRsnCutESD2010& copy);
   virtual ~AliRsnCutESD2010() {;};

   AliESDpid*       GetESDpid()  {return &fESDpid;}
   AliESDtrackCuts* GetCutsTPC() {return &fESDtrackCutsTPC;}
   AliESDtrackCuts* GetCutsITS() {return &fESDtrackCutsITS;}
   void             CopyCutsTPC(const AliESDtrackCuts *cuts) {fESDtrackCutsTPC = (*cuts);}
   void             CopyCutsITS(const AliESDtrackCuts *cuts) {fESDtrackCutsITS = (*cuts);}
   void             CopyCutsTPC(AliESDtrackCuts cuts)        {fESDtrackCutsTPC = cuts;}
   void             CopyCutsITS(AliESDtrackCuts cuts)        {fESDtrackCutsITS = cuts;}
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

   void             SetMC(Bool_t yn = kTRUE);
   void             SetCheckITS(Bool_t yn = kTRUE) {fCheckITS = yn;}
   void             SetCheckTPC(Bool_t yn = kTRUE) {fCheckTPC = yn;}
   void             SetCheckTOF(Bool_t yn = kTRUE) {fCheckTOF = yn;}
   void             SetUseITSTPC(Bool_t yn = kTRUE) {fUseITSTPC = yn;}
   void             SetUseITSSA(Bool_t yn = kTRUE) {fUseITSSA = yn;}
   void             SetPID(AliPID::EParticleType t) {fPID = t;}

   void             SetMaxITSPIDmom(Double_t v) {fMaxITSPIDmom = v;}
   void             SetITSband(Double_t v)      {fMaxITSband = v;}

   void             SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
   void             SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
   void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
   {fTPCpar[0] = p0; fTPCpar[1] = p1; fTPCpar[2] = p2; fTPCpar[3] = p3; fTPCpar[4] = p4;}

   void             SetTOFrange(Double_t v1, Double_t v2) {fMinTOF = v1; fMaxTOF = v2;}

protected:

   Bool_t  OkQuality(AliESDtrack *d);  // check track quality parameters and DCA
   Bool_t  OkITSPID(AliESDtrack *d);   // check ITS PID
   Bool_t  OkTPCPID(AliESDtrack *d);   // check TPC PID
   Bool_t  OkTOFPID(AliESDtrack *d);   // check TOF PID
   Bool_t  IsITSTPC(AliESDtrack *d);   // check that the track is TPC+ITS
   Bool_t  IsITSSA(AliESDtrack *d);    // check that the track is ITS standalone
   Bool_t  MatchTOF(AliESDtrack *d);   // check that the track matches the TOF

   Bool_t                  fIsMC;             //  switch for MC analysis
   Bool_t                  fCheckITS;         //  switch for ITS dE/dx check
   Bool_t                  fCheckTPC;         //  switch for TPC dE/dx check
   Bool_t                  fCheckTOF;         //  switch for TOF time check
   Bool_t                  fUseITSTPC;        //  switch to use TPC global tracks
   Bool_t                  fUseITSSA;         //  switch to use ITS standalone tracks
   AliPID::EParticleType   fPID;              //  PID reference type used for checks

   Double_t                fMaxITSPIDmom;     //  maximum momentum where ITS PID is used for TPC+ITS tracks
   Double_t                fMaxITSband;       //  range for ITS de/dx band

   Double_t                fTPCpLimit;        //  limit to choose what band to apply
   Double_t                fTPCpar[5];        //  parameters for TPC bethe-Bloch
   Double_t                fMinTPCband;       //  range for TPC de/dx band - min
   Double_t                fMaxTPCband;       //  range for TPC de/dx band - max

   AliESDpid               fESDpid;           //  ESD PID object
   AliESDtrackCuts         fESDtrackCutsTPC;  //  ESD standard defined track cuts for TPC tracks
   AliESDtrackCuts         fESDtrackCutsITS;  //  ESD standard defined track cuts for ITS-SA tracks
   Double_t                fMinTOF;           //  range for TOF PID (min)
   Double_t                fMaxTOF;           //  range for TOF PID (max)

   ClassDef(AliRsnCutESD2010, 1)
};

inline Bool_t AliRsnCutESD2010::IsITSTPC(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for a global track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   return vtrack->IsOn(AliESDtrack::kTPCin);

   return kTRUE;
}

inline Bool_t AliRsnCutESD2010::IsITSSA(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTPCin     = vtrack->IsOn(AliESDtrack::kTPCin);
   Bool_t isITSrefit  = vtrack->IsOn(AliESDtrack::kITSrefit);
   Bool_t isITSpureSA = vtrack->IsOn(AliESDtrack::kITSpureSA);
   Bool_t isITSpid    = vtrack->IsOn(AliESDtrack::kITSpid);

   return ((!isTPCin) && isITSrefit && (!isITSpureSA) && isITSpid);

   return kTRUE;
}


inline Bool_t AliRsnCutESD2010::MatchTOF(AliESDtrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   // require a minimum length to have meaningful match
   if (vtrack->GetIntegratedLength() < 350.) return kFALSE;

   Bool_t isTOFout = vtrack->IsOn(AliESDtrack::kTOFout);
   Bool_t isTIME   = vtrack->IsOn(AliESDtrack::kTIME);

   return (isTOFout && isTIME);

   return kTRUE;
}

#endif
