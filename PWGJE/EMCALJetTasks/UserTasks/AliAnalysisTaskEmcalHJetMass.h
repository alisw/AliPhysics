#ifndef ALIANALYSISTASKEMCALHJETMASS_H
#define ALIANALYSISTASKEMCALHJETMASS_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class TArrayF;
class TRandom3;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJet;
class AliVParticle;

#include "AliAnalysisTaskEmcalJet.h"

namespace EmcalHJetMassAnalysis {
  class AliAnalysisTaskEmcalHJetMass : public AliAnalysisTaskEmcalJet {
  public:
    enum JetMassType {
      kRaw   = 0,  //mass form anti-kt 4-vector
      kDeriv = 1   //area based subtracted jet mass
    };

    enum TriggerTrackType {
      kInclusive       = 0,  //take all trigger tracks
      kSingleInclusive = 1   //take randomly trigger track within defined pt bin
    };


    AliAnalysisTaskEmcalHJetMass();
    AliAnalysisTaskEmcalHJetMass(const char *name);
    virtual ~AliAnalysisTaskEmcalHJetMass();

    void                                UserCreateOutputObjects();
    void                                Terminate(Option_t *option);

    //Setters
    void SetDoHJetAna(Bool_t b)                                        { fDoHJetAna         = b   ; }
    void SetDoNSHJetAna(Bool_t b)                                      { fDoNSHJetAna       = b   ; }
    void SetJetContainerBase(Int_t c)                                  { fContainerBase     = c   ; }
    void SetJetContainerUnsub(Int_t c)                                 { fContainerUnsub    = c   ; }
    void SetMinFractionShared(Double_t f, Bool_t useUnsubJet = kFALSE) { fMinFractionShared = f   ; fUseUnsubJet = useUnsubJet; }
    void SetJetMassType(JetMassType t)                                 { fJetMassType       = t   ; }
    void SetMaxDeltaPhi(Double_t dphi)                                 { fDPhiHJetMax       = dphi; }
    void SetTriggerTrackType(TriggerTrackType t)                       { fTriggerTrackType  = t   ; }
    void AddTriggerTrackPtCuts(Float_t min, Float_t max);
    void SelectConstituents(UInt_t constSel)                           { fEmbConstSel = constSel  ; }
    void SetMarkMCLabel(Int_t l)                                       { fMarkMCLabel = l         ; }
    void SetGapPhiLimits(Double_t min, Double_t max)                   { fGapPhiMin = min; fGapPhiMax = max; }

  protected:
    Bool_t                              RetrieveEventObjects();
    Bool_t                              Run();
    Bool_t                              FillHJetHistograms(const AliVParticle *vp, const AliEmcalJet *jet);
    Bool_t                              FillHJetHistogramsWithNS(const AliVParticle *vp, const AliEmcalJet *jet);
    AliEmcalJet                        *FindNearSideJet(const AliVParticle *vp);

    Double_t                            GetJetMass(const AliEmcalJet *jet) const;
    Double_t                            GetDeltaPhi(const AliVParticle *vp, const AliEmcalJet* jet) const;
    Double_t                            GetDeltaPhi(Double_t phi1,Double_t phi2) const; 
    AliVParticle                       *GetSingleInclusiveTT(AliParticleContainer *pCont, Double_t ptmin, Double_t ptmax) const;

    Bool_t                              fDoHJetAna;                  // do normal h-jet analysis
    Bool_t                              fDoNSHJetAna;                // do NS h-jet analysis
    Int_t                               fContainerBase;              // jets to be analyzed
    Int_t                               fContainerUnsub;             // unsubtracted jets
    Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
    Bool_t                              fUseUnsubJet;                // calc fraction of unsubtracted jet (relevant for constituent subtraction
    JetMassType                         fJetMassType;                // jet mass type to be used
    Double_t                            fDPhiHJetMax;                // maximum delta phi between hadron and jet
    TriggerTrackType                    fTriggerTrackType;           // method to select trigger track
    TArrayF                            *fPtTTMin;                    // minimum pt of trigger tracks
    TArrayF                            *fPtTTMax;                    // maximum pt of trigger tracks
    TRandom3                           *fRandom;                     //! Random number generator
    UInt_t                              fEmbConstSel;                // select embedded constituents using bit
    Int_t                               fMarkMCLabel;                // select embedded constituents using label
    Double_t                            fGapPhiMin;                  // min phi of acceptance gap
    Double_t                            fGapPhiMax;                  // max phi of acceptance gap

    TH1F            **fh1PtHadron;                        //!pt of hadrons
    TH1F            **fh1PtHadronMatch;                   //!pt of hadrons matched to MC
    TH1F            **fh1PhiHadron;                       //!phi of hadrons
    TH3F            **fh3PtHPtJDPhi;                      //!pt hadron vs pt jet vs delta phi
    TH3F            **fh3PtJet1VsMassVsHPtAllSel;         //!all jets after std selection pt vs mass vs track pt
    TH3F            **fh3PtJet1VsMassVsHPtAllSelMatch;    //!all jets after std selection pt vs mass vs track pt matched to MC
    TH3F            **fh3PtJet1VsMassVsHPtTagged;         //!tagged jets pt vs mass vs track pt
    TH3F            **fh3PtJet1VsMassVsHPtTaggedMatch;    //!tagged jets pt vs mass vs track pt matched to MC

    TH3F            **fh3PtJet1VsRatVsHPtAllSel;          //!all jets after std selection pt vs mass/pt vs track pt
    TH3F            **fh3PtJet1VsRatVsHPtAllSelMatch;     //!all jets after std selection pt vs mass/pt vs track pt matched to MC
    TH3F            **fh3PtJet1VsRatVsHPtTagged;          //!tagged jets pt vs mass/pt vs track pt
    TH3F            **fh3PtJet1VsRatVsHPtTaggedMatch;     //!tagged jets pt vs mas/pts vs track pt matched to MC

    THnSparse       **fhnAllSel;                          //!all jets after std selection pt vs mass vs track pt vs mass_NS
    THnSparse       **fhnAllSelMatch;                     //!all jets after std selection pt vs mass vs track pt vs mass_NS matched to MC
    THnSparse       **fhnTagged;                          //!tagged jets pt vs mass vs track pt vs mass_NS
    THnSparse       **fhnTaggedMatch;                     //!tagged jets pt vs mass vs track pt vs mass_NS matched to MC

  private:
    AliAnalysisTaskEmcalHJetMass(const AliAnalysisTaskEmcalHJetMass&);            // not implemented
    AliAnalysisTaskEmcalHJetMass &operator=(const AliAnalysisTaskEmcalHJetMass&); // not implemented

    ClassDef(AliAnalysisTaskEmcalHJetMass, 9)
      };
}
#endif

