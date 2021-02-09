#ifndef ALIRDHFCUTSDSTOK0SK_H
#define ALIRDHFCUTSDSTOK0SK_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
///
/// \class AliRDHFCutsDstoK0sK
/// \brief Class for cuts on AOD reconstructed Ds->K0S+K
///
/// \author J.Hamon, julien.hamon@cern.ch (IPHC)
///
/////////////////////////////////////////////////////////////


#include "AliRDHFCuts.h"

class AliAODEvent;
class AliAODRecoDecayHF;
class AliESDtrackCuts;


class AliRDHFCutsDstoK0sK : public AliRDHFCuts
{
public:


   AliRDHFCutsDstoK0sK(const char* name="CutsDstoK0sK");
   AliRDHFCutsDstoK0sK(const AliRDHFCutsDstoK0sK& source);
   AliRDHFCutsDstoK0sK& operator=(const AliRDHFCutsDstoK0sK& source);
   virtual ~AliRDHFCutsDstoK0sK();


   using   AliRDHFCuts::GetCutVarsForOpt;
   virtual void GetCutVarsForOpt(AliAODRecoDecayHF* obj, Float_t* vars, Int_t nvars, Int_t*pdgdaughters)
                           { return GetCutVarsForOpt(obj, vars, nvars, pdgdaughters, 0x0); }
   virtual void GetCutVarsForOpt(AliAODRecoDecayHF* obj, Float_t* vars, Int_t nvars, Int_t* pdgdaughters, AliAODEvent* aod);

   using   AliRDHFCuts::IsSelected;
   virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel)
                           { return IsSelected(obj, selectionLevel, 0); }
   virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod);

   using   AliRDHFCuts::IsSelectedPID;
   virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);

   using   AliRDHFCuts::PreSelect;
   Bool_t PreSelect(TObject* obj, AliAODv0 *v0, AliVTrack *bachelorTrack);

   using   AliRDHFCuts::IsInFiducialAcceptance;
   virtual Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y) const;


   Bool_t AreDtoK0sDaughtersSelected(AliAODRecoDecayHF *rd) const;
   Int_t  GetV0Type();


   void AddTrackCutsV0daughters(AliESDtrackCuts* v0daug)
                           { delete fV0daughtersCuts; fV0daughtersCuts = new AliESDtrackCuts(*v0daug); }
   virtual AliESDtrackCuts* GetTrackCutsV0daughters() const
                           { return fV0daughtersCuts; }

   Float_t GetMassCut(Int_t iPtBin=0) const
                           { return (GetCuts() ? fCutsRD[GetGlobalIndex(0, iPtBin)] : 1.e6);}
   Int_t GetExcludedCut()
                           { return fExcludedCut; }
   void SetExcludedCut(Int_t excludedCut)
                           { fExcludedCut = excludedCut; }
   Float_t GetV0PtCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(8,iPtBin)] : 0.);}
   Float_t GetMinV0PtCut() const {
     Float_t minPtCut=99999.;
     for(Int_t j=0; j<fnPtBins; j++){Float_t c=GetV0PtCut(j); if(c<minPtCut) minPtCut=c;}
     return minPtCut;
   }


protected:


   enum EPid {kConservative, kStrong};
   void SetPidOption(Int_t opt)        { fPidOption      = opt; }
   void SetMaxPtStrongPid(Float_t pid) { fMaxPtStrongPid = pid; }

   Int_t             fExcludedCut;           /// cut to be excluded (-1=none)
   Float_t           fV0Type;                /// V0 type -- should be defined as in AliRDHFCuts.h
   AliESDtrackCuts*  fV0daughtersCuts;       /// cuts for v0 daughters (AOD converted to ESD on the fly!)
   Int_t             fPidOption;             /// PID option
   Float_t           fMaxPtStrongPid;        /// Maximum pt of candidate to apply strong PID p dependent


   /// \cond CLASSIMP
   ClassDef(AliRDHFCutsDstoK0sK, 2);   /// class for cuts on AOD reconstructed D+ -> K0S + pi
   /// \endcond

};

#endif
