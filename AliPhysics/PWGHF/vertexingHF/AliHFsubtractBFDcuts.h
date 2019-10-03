#ifndef ALIHFSUBTRACTBFDCUTS_H
#define ALIHFSUBTRACTBFDCUTS_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
/// \class AliHFsubtractBFDcuts                                           //
///                                                                       //
///  \brief Class for storing and handling D0 meson candidates properties //
///  for estimating the feed-down fraction using several sets of cuts     //
///     \author Andrea Rossi <andrea.rossi@cern.ch>                       //
///     \author Felix Reidt  <felix.reidt@cern.ch>                        //
///                                                                       //
////////////////////////////////////////////////////////////////////////////
#include <vector>

#include "TNamed.h"
#include "THnSparse.h"

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"

class TList;
class TClonesArray;
class AliAODMCHeader;
class AliAODEvent;
class AliAODVertex;
class AliNeutralTrackParam;

class AliHFsubtractBFDcuts : public TNamed{

public:
  AliHFsubtractBFDcuts();
  AliHFsubtractBFDcuts(const char* name,const char* title);
  ~AliHFsubtractBFDcuts();

  void InitHistos();

  void FillGenStep(AliAODMCParticle* dzeroMC,Double_t pt=-1,Double_t weight=1.,TClonesArray* mcArray=0x0, AliAODMCHeader* mcHeader=0x0);
  THnSparseF* GetSparseData() const {return fTHnData;}
  THnSparseF* GetSparseMC()   const {return fTHnMC;}
  THnSparseF* GetSparseMCgen()     const {return fTHnGenStep;}

  void SetHistoPtMCgen(THnSparseF* h){if(fTHnGenStep)delete fTHnGenStep;fTHnGenStep=(THnSparseF*)h->Clone();return;}
  void SetSparseData(THnSparseF* h){if(fTHnData)delete fTHnData;fTHnData=(THnSparseF*)h->Clone();return;}
  void SetSparseMC(THnSparseF* h){if(fTHnMC)delete fTHnMC;fTHnMC=(THnSparseF*)h->Clone();return;}

  void SetFillMC (Bool_t fillMC = kTRUE) {fIsMC = fillMC;}

  void FillSparses(AliAODRecoDecayHF2Prong* dzeroPart,Int_t isSelected,Double_t pt=-1,Double_t massD0=-1,Double_t massD0bar=-1,Double_t weight=1.,TClonesArray* mcArray=0x0, AliAODEvent* aodEvent=0x0, AliAODMCHeader* mcHeader=0x0);
  TList* GetDecayStrings() { return fDecayStrList; }
  TList* GetQAhists() { return fQAhists; }

private:
  AliHFsubtractBFDcuts(const AliHFsubtractBFDcuts& c);
  AliHFsubtractBFDcuts operator=(const AliHFsubtractBFDcuts& c);

  Bool_t GetCandidateLabel();
  Bool_t AnalyseDecay(Bool_t generateString, Bool_t mcOnly);  /// check in which decay process a particle was created
  void   CountProngs(Int_t labCurrMother, Int_t labCurrExcl, Bool_t generateString, Bool_t mConly);
  /// counting the prongs of labCurrMother, labCurrExcl is assumed to be a stable particle
  Bool_t IsStable(Int_t labProng);             /// Is that prong a stable particle?
  Bool_t IsInAcceptance(Int_t labProng) const; /// Is that prong within the fiducial acceptance
  AliAODVertex* RecBvtx(TObjArray* tracks) const; /// Reconstruct a secondary vertex with the supplied tracks
  Bool_t CheckBhypothesis(Int_t iAODtrack, Bool_t Bprong); /// Method to check Whether the current D0 candidate and the track originate from a B decay

  Bool_t      fIsMC;              /// flag for MC/Data
  Bool_t      fCheckAcceptance;   /// flag for checking whether the decay prongs are within acceptance
  Bool_t      fResolveResonances; /// flag resolve resonances in during the prong determination
  THnSparseF* fTHnGenStep;        //!<! histo with spectrum at generation level
  THnSparseF* fTHnData;           //!<! THnSparse for cut variables (data, with inv mass axis), first axis is always mass
  THnSparseF* fTHnMC;             //!<! THnSparse for cut variables (MC at PID level, w/o mass axis)y
  TList*      fQAhists;           //!<! List with QA histograms

  /// Event specific variables
  TClonesArray*             fMCarray;      //!<! TClonesArray holding the particles of the event to be processed
  TClonesArray*             fAODtracks;    //!<! TClonesArray holding the AliAODTracks of the event to be processed
  AliAODVertex*             fPriVtx;       //!<! Primary AOD vertex
  Double_t                  fBkG;          /// Magnetic field (z-direction) in units of kG
  AliAODRecoDecayHF2Prong*  fD0Cand;       /// Pointer to the D0 candidate from reconstruction
  AliNeutralTrackParam*     fD0CandParam;  /// Pointer to an AliNeutralTrackParam of the D0 candidata for DCA calculation
  Int_t                     fLabCand;      /// Label of the candidate D0 (charmed hadron in case of a chained decay)
  Double_t                  fPtCand;       /// pT of the candidate (from MC track, not following decay chain)
  Int_t                     fLabMother;    /// Label of the mother of the candidate D0 (or charmed hadron)
  UInt_t                    fNprongs;      /// Number of prongs, counting the first charmed hadron as one particle (simulation cuts can lead to loss of prongs!)
  UInt_t                    fNprongsInAcc; /// Number of prongs, counting only the particles within acceptance
  Bool_t                    fFoundElectron;/// Does the B meson decay contain an electron?
  Bool_t                    fDecayChain;   /// Chained decay of charmed hadrons
  Double_t                  fMotherPt;     /// Transverse momentum of the mother particle (B hadron in case of feed-down,
                                           /// the charmed hadron itsself in case of prompt production)
  Bool_t             fGenerateDecayList;   /// Generate the list containig strings with all PDG codes of the decay prongs
  std::vector<Int_t> fDecayProngs;         /// PDG codes of the daughters separated
  TList*             fDecayStrList;        //!<! List with all decay strings

  /// \cond CLASSIMP
  ClassDef(AliHFsubtractBFDcuts,10);
  /// \endcond
};

#endif
