#ifndef ALIRDHFCUTSDPLUSTOK0SPI_H
#define ALIRDHFCUTSDPLUSTOK0SPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
///
/// \class AliRDHFCutsDplustoK0spi
/// \brief Class for cuts on AOD reconstructed D+->K0S+pi
///
/// \author J.Hamon, julien.hamon@cern.ch (IPHC)
///
/////////////////////////////////////////////////////////////


#include "AliRDHFCuts.h"

class AliAODEvent;
class AliAODRecoDecayHF;
class AliESDtrackCuts;


class AliRDHFCutsDplustoK0spi : public AliRDHFCuts
{
public:


   AliRDHFCutsDplustoK0spi(const char* name="CutsDplustoK0spi");
   AliRDHFCutsDplustoK0spi(const AliRDHFCutsDplustoK0spi& source);
   AliRDHFCutsDplustoK0spi& operator=(const AliRDHFCutsDplustoK0spi& source);
   virtual ~AliRDHFCutsDplustoK0spi();


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


   Int_t             fExcludedCut;           /// cut to be excluded (-1=none)
   Float_t           fV0Type;                /// V0 type -- should be defined as in AliRDHFCuts.h
   AliESDtrackCuts*  fV0daughtersCuts;       /// cuts for v0 daughters (AOD converted to ESD on the fly!)


   /// \cond CLASSIMP
   ClassDef(AliRDHFCutsDplustoK0spi, 1);   /// class for cuts on AOD reconstructed D+ -> K0S + pi
   /// \endcond

};

#endif
