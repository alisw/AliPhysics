#ifndef ALIANALYSISTASKEMCALJETMASSSTRUCTURE_H
#define ALIANALYSISTASKEMCALJETMASSSTRUCTURE_H

class TH1;
class TH2;
class TH3;
class TProfile;
class THnSparse;
class TClonesArray;
class TArrayI;
class TRandom3;
class TList;

class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJet;
class AliEmcalJetByJetCorrection;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetMassStructure : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassType {
    kRaw   = 0,  //mass form anti-kt 4-vector
    kDeriv = 1   //area based subtracted jet mass
  };

  enum JetByJetCorrType {
    kNoCorr  = 0,
    kAnnulus = 1, //reproduce existing particles (over-corrects)
    kMeanPtR = 2
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
  void SetJetByJetCorrType(JetByJetCorrType t)                  { fCorrType          = t   ; }
  void SetJetByJetCorrObject(AliEmcalJetByJetCorrection *a)     { fEJetByJetCorr     = a   ; }
  void SetParticleArray(TString particles)                { fPartArrayN         = particles;}
  Int_t CalculateNMissingTracks(AliEmcalJet *jet1, AliEmcalJet *jPart);
  //Getters
  AliEmcalJetByJetCorrection *GetJetByJetCorrObject() const     { return fEJetByJetCorr    ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Bool_t                              GetSortedArray(Int_t indexes[], Float_t dr[], AliEmcalJet *jet) const;
  Double_t                            GetJetMass(AliEmcalJet *jet);
  Double_t                            GetEfficiency(Double_t pt);
 
  Int_t                               fContainerBase;              ///< jets to be analyzed
  Double_t                            fMinFractionShared;          ///< only fill histos for jets if shared fraction larger than X
  JetMassType                         fJetMassType;                ///< jet mass type to be used

  TH1F                               *fhtmppTRec;		   ///< Temporary stores the pT distribution of jet constituents (reco level)
  TH1F                               *fhtmppTGen;                  ///< Temporary stores the pT distribution of jet constituents (gen level)

  TRandom3                           *fRandom;                     //!<! random number generator
  Double_t                            fEfficiencyFixed;            ///< fixed efficiency for all pT and all types of tracks
                                                                     
  JetByJetCorrType                    fCorrType;                   ///< jet-by-jet correction method
  AliEmcalJetByJetCorrection         *fEJetByJetCorr;              ///< object to do jet-by-jet correction
 
  TH3F                              **fh3PtDRMass;                 //!<! jet pT vs dr(jet axis, constituent) vs cumulative mass density
  TH3F                              **fh3PtDRRho;                  //!<! jet pT vs dr(jet axis, constituent) vs cumulative pt density
  TH3F                              **fh3PtDRMassCorr;             //!<! jet pT vs dr(jet axis, constituent) vs cumulative mass density corrected
  TH3F                              **fh3PtDRRhoCorr;              //!<! jet pT vs dr(jet axis, constituent) vs cumulative pt density corrected
  TH2F                              **fh2PtMass;                   //!<! jet pT vs mass
  TH2F                              **fh2PtMassCorr;               //!<! jet pT vs mass corrected
  THnSparse                          *fhnMassResponse;             //!<! response matrix
  THnSparse                          *fhnMassResponseCorr;         //!<! response matrix corrected
  TH3F                              **fh3JetPtDRTrackPt;           //!<! jet pt vs dr(jet axis, constituent) vs pT,track
  THnSparse                          *fhnDeltaMass;                //!<! resolution on mass matrix
  THnSparse                          *fhnDeltaMassCorr;            //!<! resolution on mass matrix corrected
  TH2F                               *fhAllpTRec;
  //!<! histogram that stores the pT of the constituents (RECO level)
  TH2F                               *fhAllpTGen;
  //!<! histogram that stores the pT of the constituents (PART level)
  TH2F                               *fhAllpTCor;
  //!<! histogram that stores the pT of the constituents (CORR level)
  THnSparse                          *fhConstRecGen;           //!<! number of constituent correlation
  TList                              *fListOfOutputFromClass;      //!<! list of output from class AliEmcalJetByJetCorrection
  Bool_t                             fSwitchResolutionhn;          ///< switch on/off (default on) the 2 THnSparse for the mass resolution
  TString                       fPartArrayN;                 ///< Array of particles used for jet reconstruction at particle level (need to make it transient probably)
  private:
  AliAnalysisTaskEmcalJetMassStructure(const AliAnalysisTaskEmcalJetMassStructure&);            // not implemented
  AliAnalysisTaskEmcalJetMassStructure &operator=(const AliAnalysisTaskEmcalJetMassStructure&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassStructure, 6)
};
#endif

