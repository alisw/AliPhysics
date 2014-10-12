#ifndef ALIANALYSISTASKEMCALJETMASSSTRUCTURE_H
#define ALIANALYSISTASKEMCALJETMASSSTRUCTURE_H

class TH1;
class TH2;
class TH3;
class TProfile2D;
class THnSparse;
class TClonesArray;
class TArrayI;
class TRandom3;

class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJet;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetMassStructure : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassType {
    kRaw   = 0,  //mass form anti-kt 4-vector
    kDeriv = 1   //area based subtracted jet mass
  };

  AliAnalysisTaskEmcalJetMassStructure();
  AliAnalysisTaskEmcalJetMassStructure(const char *name);
  virtual ~AliAnalysisTaskEmcalJetMassStructure();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetJetMassType(JetMassType t)                            { fJetMassType       = t   ; }
  void SetFixedTrackEfficiency(Double_t eff)                    { fEfficiencyFixed   = eff ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Bool_t                              GetSortedArray(Int_t indexes[], Float_t dr[], AliEmcalJet *jet) const;
  Double_t                            GetJetMass(AliEmcalJet *jet);
  Double_t                            GetEfficiency(Double_t pt);
 
  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetMassType                         fJetMassType;                // jet mass type to be used

  TRandom3                           *fRandom;                     //! random number generator
  Double_t                            fEfficiencyFixed;            // fixed efficiency for all pT and all types of tracks
 
  TH3F                              **fh3PtDRMass;                 //! jet pT vs dr(jet axis, constituent) vs cumulative mass density
  TH3F                              **fh3PtDRRho;                  //! jet pT vs dr(jet axis, constituent) vs cumulative pt density
  TH3F                              **fh3PtDRMassCorr;             //! jet pT vs dr(jet axis, constituent) vs cumulative mass density corrected
  TH3F                              **fh3PtDRRhoCorr;              //! jet pT vs dr(jet axis, constituent) vs cumulative pt density corrected
  TH2F                              **fh2PtMass;                   //! jet pT vs mass
  TH2F                              **fh2PtMassCorr;               //! jet pT vs mass corrected
  THnSparse                          *fhnMassResponse;             //! response matrix
  THnSparse                          *fhnMassResponseCorr;         //! response matrix corrected

  TH3F                              **fh3JetPtDRTrackPt;           //! jet pt vs dr(jet axis, constituent) vs pT,track

 private:
  AliAnalysisTaskEmcalJetMassStructure(const AliAnalysisTaskEmcalJetMassStructure&);            // not implemented
  AliAnalysisTaskEmcalJetMassStructure &operator=(const AliAnalysisTaskEmcalJetMassStructure&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassStructure, 1)
};
#endif

