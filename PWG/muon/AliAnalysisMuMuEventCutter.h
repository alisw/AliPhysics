#ifndef ALIANALYSISMUMUEVENTCUTTER_H
#define ALIANALYSISMUMUEVENTCUTTER_H

/**
 * 
 * \class AliAnalysisMuMuEventCutter
 *
 * \brief Various event cuts used in AliAnalysisTaskMuMu
 *
 * \author L. Aphecetche and J. Martin-Bianco
 *
 */

#include "TObject.h"
#include "TString.h"

class AliMuonEventCuts;
class TList;
class AliVEvent;
class AliVVertex;
class AliInputEventHandler;

class AliAnalysisMuMuEventCutter : public TObject
{
public:
  AliAnalysisMuMuEventCutter(TList* triggerClassesToConsider=0x0);
  virtual ~AliAnalysisMuMuEventCutter();
  
  Bool_t SelectTriggerClass(const TString& firedTriggerClasses, TString& acceptedClasses,
                            UInt_t L0, UInt_t L1, UInt_t L2) const;
  
  Bool_t IsTrue(const AliVEvent& /*event*/) const { return kTRUE; }
  void NameOfIsTrue(TString& name) const { name="ALL"; }

  Bool_t IsFalse(const AliVEvent& /*event*/) const { return kFALSE; }
  void NameOfIsFalse(TString& name) const { name="NONE"; }

  Bool_t IsPhysicsSelected(const AliInputEventHandler& eventHandler) const;
  void NameOfIsPhysicsSelected(TString& name) const { name="PSALL"; }

  Bool_t IsPhysicsSelectedVDM(const AliVEvent& event) const;
  void NameOfIsPhysicsSelectedVDM(TString& name) const { name="VDM"; }

  Bool_t IsAbsZBelowValue(const AliVEvent& event, const Double_t& z) const;
  void NameOfIsAbsZBelowValue(TString& name, const Double_t& z) const;
  
  Bool_t IsAbsZSPDBelowValue(const AliVEvent& event, const Double_t& z) const;
  void NameOfIsAbsZSPDBelowValue(TString& name, const Double_t& z) const;
  
  Bool_t IsSPDzVertexInRange(AliVEvent& event, const Double_t& zMin, const Double_t& zMax) const;
  void NameOfIsSPDzVertexInRange(TString& name, const Double_t& zMin, const Double_t& zMax) const;
  
  Bool_t IsSPDzQA(const AliVEvent& event/*, const AliVVertex& vertex2Test*/, const Double_t& zResCut, const Double_t& zDifCut) const;
  void NameOfIsSPDzQA(TString& name, const Double_t& zResCut, const Double_t& zDifCut) const;
  
  Bool_t HasSPDVertex(AliVEvent& event) const;
  void NameOfHasSPDVertex(TString& name) const { name = "HASSPD"; }
  
  Bool_t IsMeandNchdEtaInRange(AliVEvent& event, const Double_t& dNchdEtaMin, const Double_t& dNchdEtaMax) const;
  void NameOfIsMeandNchdEtaInRange(TString& name, const Double_t& dNchdEtaMin, const Double_t& dNchdEtaMax) const;
  
  Bool_t IsTZEROPileUp(const AliVEvent& event) const;
  void NameOfIsTZEROPileUp(TString& name) const { name="TZEROPILEUP"; }
  
  AliMuonEventCuts* MuonEventCuts() const;
  
//  enum EEventCut
//  {
//    kEventIR2PILEUP     = BIT( 6), /// events with pile-up (using AliAnalysisUtils::IsOutOfBunchPileUp)
//    
//    kEventGOODVERTEX    = BIT(10), /// events with a good vertex
//    kEventZPOS          = BIT(14), /// events with z > 0
//    kEventZNEG          = BIT(15), /// events with z < 0
//    
//    kEventTRKLETA1      = BIT(20), /// event with at least one tracklet in |eta| < fTrackletEtaCutValue[0]

//    kEvent0TVX          = BIT(31), /// events with 0TVX L0 input
//    kEventV0AND         = BIT(32), /// events with both 0V0C and 0V0A L0 inputs
//    kEvent0SM2          = BIT(33), /// events with 0SM2 L0 input
//    kEvent0MSL          = BIT(34), /// events with 0MSL input
//  };
  
private:
  
  /// not implemented on purpose
  AliAnalysisMuMuEventCutter& operator=(const AliAnalysisMuMuEventCutter& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuEventCutter(const AliAnalysisMuMuEventCutter& rhs);
  
  mutable AliMuonEventCuts* fMuonEventCuts; // common cuts for muon events (from Diego)
  
  ClassDef(AliAnalysisMuMuEventCutter,1) // default event cutters for AliAnalysisTaskMuMu
};

#endif
