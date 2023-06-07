#ifndef ALITWOPARTICLECORRELATIONSBASE_H
#define ALITWOPARTICLECORRELATIONSBASE_H

/// \file AliTwoParticleCorrelationsBase.h
/// \brief Base class for collecting data for two-particle correlation functions
///

#include <TNamed.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>

/// \class AliTwoParticleCorrelationsBase
/// \brief Establishes the protocol for collecting data for building
/// two-particle correlation functions
///
/// \author Víctor González <victor.gonzalez@cern.ch>, WSU
/// \date May 2021

class AliVEvent;
class AliVTrack;
class AliVParticle;

class AliTwoParticleCorrelationsBase : public TNamed
{
public:
  AliTwoParticleCorrelationsBase();
  AliTwoParticleCorrelationsBase(const char *name);
  virtual ~AliTwoParticleCorrelationsBase();
  /// \brief Use only singles, do not produce pairs data
  /// \param so flag to activate deactivate
  void SetSinglesOnly(Bool_t so) { fSinglesOnly = so; }
  /// \brief Use calibration weights
  /// \param uw flag to activate deactivate the use of weights
  void SetUseWeights(Bool_t uw) { fUseWeights = uw; }
  /// \brief Use simulated tracks
  /// \param us flag to activate deactivate the use of simulated tracks
  void SetUseSimulation(Bool_t us) { fUseSimulation = us; }
  /// \brief Get if simulation is being used
  /// \return kTRUE if track simulation is being used kFALSE otherwise
  Bool_t GetUseSimulation() { return fUseSimulation; }
  /// \brief Set if both tracks are of the same charge sign
  /// \param ss kTRUE if both tracks have the same charge kFALSE otherwise
  virtual void SetSameSign(Bool_t) {}
  /// \brief Produce data for all track combinations
  /// \param all kTRUE if all track combinations have to be produced kFALSE otherwise
  virtual void SetAllCombinations(Bool_t) {}
  /// \brief The requested charge for first track
  /// \param q charge, greater than one if positive negative otherwise
  virtual void SetRequestedCharge_1(Int_t) {}
  /// \brief The requested charge for second track
  /// \param q charge, greater than one if positive negative otherwise
  virtual void SetRequestedCharge_2(Int_t) {}

  virtual Bool_t ConfigureBinning(const char *configstring);
  void ConfigureResonances(const char *confstring);
  virtual TString GetBinningConfigurationString() const;
  TString GetResonancesConfigurationString() const;
  void SetSpeciesNames(std::vector<std::string> names);

  bool SetWeigths(std::vector<const TH3 *> h3);
  virtual bool SetPtAvg(std::vector<const TH2 *>) { return false; }
  bool SetEfficiencyCorrection(std::vector<const TH1 *> h1);
  bool SetPairEfficiencyCorrection(std::vector<std::vector<const THn *>> hn);
  bool SetSimultationPdfs(std::vector<const TObjArray *> pdf);
  /// \brief set the number of simulated events per passed event
  /// \param nevents the number of events to simulate per passed event
  void SetSimEventsPerEvent(Int_t nevents) { fSimEventsPerEvent = nevents; }
  /// \brief get the number of events being simulated per passed event
  /// \return the number of events being simulated per passed event
  int GetSimEventsPerEvent() { return fSimEventsPerEvent; }

  virtual void Initialize();
  /// Get the histograms list
  /// \return the histograms list
  TList *GetHistogramsList() { return fOutput; }
  virtual bool StartEvent(float centrality, float vertexZ);
  virtual bool IsTrackAccepted(AliVTrack *trk);
  virtual bool IsTrackAccepted(AliVParticle *part);
  virtual bool IsTrackAccepted(float pT, float eta, float phi);
  virtual bool ProcessTrack(int id, AliVTrack *trk);
  virtual bool ProcessTrack(int id, AliVParticle *part);
  virtual bool ProcessTrack(int pid, float pT, float eta, float phi) = 0;
  virtual void ProcessEventData() = 0;
  virtual void FinalizeProcess() = 0;

  /// \brief has the pre particle selection resonance rejection been configured
  /// \return true if the pre rejection has been configured otherwise false
  bool preRejectResonances() { return fPreRejectResonances; }
  /// \brief gets the index of the photon conversion pair rejection structures
  /// \return the index of the photon conversion (0)
  static int getPhotonConversionIndex() { return 0; }
  float checkIfResonance(Int_t ires, Bool_t fill, float pt1, float eta1, float phi1, float pt2, float eta2, float phi2);

protected:
  void FlagConversionsAndResonances();

protected:
  TList *fOutput; //!<! Output histograms list

  Float_t fVertexZ;    ///< the current event vertex \f$z\f$ coordinate
  Int_t fIxVertexZ;    ///< bin index for the current event vertex \f$z\f$ coordinate
  Float_t fCentrality; ///< current event centrality in percentage

  Bool_t fSinglesOnly;   ///< true if only collect singles data
  Bool_t fUseWeights;    ///< kTRUE if correction weights must be utilized
  Bool_t fUseSimulation; ///< kTRUE if particle production from stored profiles must be used

  /* the arrays with tracks information */
  int fArraySize;  ///< the size of the array to collect accepted track information
  int fNoOfTracks; ///< the number of stored tracks
  UInt_t* fFlags;    //!<! the array of track flags
  float* fPt;      //!<! the array of track \f$p_T\f$
  float* fEta;     //!<! the array of track \f$\eta\f$
  float* fPhi;     //!<! the array of track \f$\varphi\f$
  int* fPID;       //!<! the array of track Ids

  std::vector<float *> fCorrectionWeights;                ///< structure with the track correction weights
  std::vector<float *> fEfficiencyCorrection;             ///< structure with the track efficiency correction weights
  std::vector<std::vector<const THn *>> fPairsEfficiency; ///< pair efficiency correction for the different pair combinations
  std::vector<const TObjArray *> fTrackPdfs;              ///< the tracks density distributions
  std::vector<TH3F *> fTrackCurrentPdf;                   ///< current event tracks density distribution
  int fSimEventsPerEvent;                                 ///< the number of simulated events produced per real event

  int fNBins_vertexZ;    ///< the \f$z\f$ vertex component number of bins
  double fMin_vertexZ;   ///< the minimum \f$z\f$ vertex component value
  double fMax_vertexZ;   ///< the maximum \f$z\f$ vertex component value
  double fWidth_vertexZ; ///< the \f$z\f$ vertex component bin width

  double fNBinsPhiShift; ///< the number of bins the phi origin is shifted

  int fNBins_pt;        ///< the track \f$p_T\f$ number of bins
  double fMin_pt;       ///< the track minimum \f$p_T\f$ value
  double fMax_pt;       ///< the track maximum \f$p_T\f$ value
  double fWidth_pt;     ///< the track \f$p_T\f$ bin width
  int fNBins_phi;       ///< the track \f$\phi\f$ number of bins
  double fMin_phi;      ///< the track minimum \f$\phi\f$ value
  double fMax_phi;      ///< the track maximum \f$\phi\f$ value
  double fWidth_phi;    ///< the track \f$\phi\f$ bin width
  bool fExcludeTOFHole; ///< exclude azimuthally the TOF hole
  int fNBins_eta;       ///< the track \f$\eta\f$ number of bins
  double fMin_eta;      ///< the track minimum \f$\eta\f$ value
  double fMax_eta;      ///< the track maximum \f$\eta\f$ value
  double fWidth_eta;    ///< the track \f$\eta\f$ bin width
  int fNBins_etaPhi;    ///< the track combined \f$\eta, \phi\f$ number of bins
  int fNBins_etaPhiPt;  ///< the track combined \f$\eta, \phi, p_T\f$ number of bins
  int fNBins_zEtaPhiPt; ///< the combined event \f$z\f$ vertex component and track \f$\eta, \phi, p_T\f$ number of bins

  std::vector<double> fN1;       ///< weighted number of tracks per species for current event
  std::vector<double> fNnw1;     ///< not weighted number of tracks per species for current event
  std::vector<double> fSum1Pt;   ///< accumulated sum of weighted \f$p_T\f$ per species for current event
  std::vector<double> fSum1Ptnw; ///< accumulated sum of not weighted \f$p_T\f$ per species for current event

  /* histograms */
  std::vector<TH1F *> fhN1_vsPt;            ///< weighted single particle distribution vs \f$p_T\f$ per species
  std::vector<TH2F *> fhN1_vsEtaPhi;        ///< weighted single particle distribution vs \f$\eta,\;\phi\f$ per species
  std::vector<TH2F *> fhSum1Pt_vsEtaPhi;    ///< accumulated sum of weighted \f$p_T\f$ vs \f$\eta,\;\phi\f$ per species
  std::vector<TH3F *> fhN1_vsZEtaPhiPt;     ///< single particle distribution vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$ per species
  std::vector<TH3F *> fhSum1Pt_vsZEtaPhiPt; ///< accumulated sum of weighted \f$p_T\f$ vs \f$\mbox{vtx}_z,\; \eta,\;\phi,\;p_T\f$ per species
  /* versus centrality  profiles */
  std::vector<TProfile *> fhN1_vsC;       ///< weighted single particle distribution vs event centrality/multiplicity per species
  std::vector<TProfile *> fhSum1Pt_vsC;   ///< accumulated sum of weighted \f$p_T\f$ vs event centrality/multiplicity per species
  std::vector<TProfile *> fhN1nw_vsC;     ///< un-weighted single particle distribution vs event centrality/multiplicity per species
  std::vector<TProfile *> fhSum1Ptnw_vsC; ///< accumulated sum of un-weighted \f$p_T\f$ vs event centrality/multiplicity per species

 protected:
  std::vector<std::string> fSpeciesNames;                   ///< the name of the species to consider
  static Int_t fgkNoOfResonances;                           ///< the number of resonances conversions to consider
  Int_t                       fThresholdMult[16];           ///< the threshold multiplier, in 1/4 modulus units (i.e, four is one modulus, zero disable it)
  bool fPreRejectResonances;                                ///< flag if the resonances have to be prerejected before PID
  bool fPostRejectResonances;                               ///< flag if the resonances have to be rejected after PID once reached to the correlations engine, this one
 private:
  static Double_t             fgkMass[16];                  ///< the masses of resonances conversions to consider
  static Double_t             fgkChildMass[2][16];          ///< the masses of the resonances / conversions products
  static Double_t             fgkMassThreshold[16];         ///< the resonance / conversion mass threshold modulus
  TH2F                       *fhResonanceRoughMasses;       ///< the resonance approximate invariant mass histogram
  TH2F                       *fhResonanceMasses;            ///< the resonance invariant mass histogram
protected:
  TH2F                       *fhDiscardedResonanceMasses;   ///< the discarded resonance invariant mass histogram


private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliTwoParticleCorrelationsBase(const AliTwoParticleCorrelationsBase &);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliTwoParticleCorrelationsBase &operator=(const AliTwoParticleCorrelationsBase &);

  /// \cond CLASSIMP
  ClassDef(AliTwoParticleCorrelationsBase, 5);
  /// \endcond
};

inline Float_t GetSquaredInvMass(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  /* original reference */
  /* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */
  // encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
  // Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/

  Float_t tantheta1 = 1e10;

  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta1);
    tantheta1 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta2);
    tantheta2 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (TMath::Sqrt(e1squ * e2squ) - (pt1 * pt2 * (TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2)));

  return mass2;
}

inline Float_t GetSquaredInvMassCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  /* original reference */
  /* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */
  // encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
  // Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

  // calculate inv mass squared approximately

  Float_t tantheta1 = 1e10;

  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0 - eta1 + eta1 * eta1 / 2 - eta1 * eta1 * eta1 / 6 + eta1 * eta1 * eta1 * eta1 / 24;
    tantheta1 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = 1.0 - eta2 + eta2 * eta2 / 2 - eta2 * eta2 * eta2 / 6 + eta2 * eta2 * eta2 * eta2 / 24;
    tantheta2 = 2.0 * expTmp / (1.0 - expTmp * expTmp);
  }

  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

  // fold onto 0...pi
  Float_t deltaPhi = TMath::Abs(phi1 - phi2);
  while (deltaPhi > TMath::TwoPi())
    deltaPhi -= TMath::TwoPi();
  if (deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;

  Float_t cosDeltaPhi = 0;
  if (deltaPhi < TMath::Pi() / 3)
    cosDeltaPhi = 1.0 - deltaPhi * deltaPhi / 2 + deltaPhi * deltaPhi * deltaPhi * deltaPhi / 24;
  else if (deltaPhi < 2 * TMath::Pi() / 3)
    cosDeltaPhi = -(deltaPhi - TMath::Pi() / 2) + 1.0 / 6 * TMath::Power((deltaPhi - TMath::Pi() / 2), 3);
  else
    cosDeltaPhi = -1.0 + 1.0 / 2.0 * (deltaPhi - TMath::Pi()) * (deltaPhi - TMath::Pi()) - 1.0 / 24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);

  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (TMath::Sqrt(e1squ * e2squ) - (pt1 * pt2 * (cosDeltaPhi + 1.0 / tantheta1 / tantheta2)));

  return mass2;
}

inline Float_t AliTwoParticleCorrelationsBase::checkIfResonance(Int_t ires, Bool_t fill, float pt1, float eta1, float phi1, float pt2, float eta2, float phi2)
{
  /* inspired on */
  /* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */
  // encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
  // Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

  Bool_t itcouldbe = kFALSE;
  Float_t mass = GetSquaredInvMassCheap(pt1, eta1, phi1, pt2, eta2, phi2, fgkChildMass[0][ires], fgkChildMass[1][ires]);

  if (TMath::Abs(mass - fgkMass[ires] * fgkMass[ires]) < 5 * fgkMassThreshold[ires]) {
    if (fill)
      fhResonanceRoughMasses->Fill(ires, TMath::Sqrt(mass));
    mass = GetSquaredInvMass(pt1, eta1, phi1, pt2, eta2, phi2, fgkChildMass[0][ires], fgkChildMass[1][ires]);

    Float_t low = ((fgkMass[ires] != 0.0) ? (fgkMass[ires] - fThresholdMult[ires] * 0.5) * (fgkMass[ires] - fgkMassThreshold[ires] * 0.5) : 0.0);
    Float_t high = (fgkMass[ires] + fgkMassThreshold[ires] * 0.5) * (fgkMass[ires] + fgkMassThreshold[ires] * 0.5);

    if ((low < mass) && (mass < high)) {
      itcouldbe = kTRUE;
    } else if (fgkChildMass[0][ires] != fgkChildMass[1][ires]) {
      /* switch masses hypothesis */
      mass = GetSquaredInvMass(pt1, eta1, phi1, pt2, eta2, phi2, fgkChildMass[1][ires], fgkChildMass[0][ires]);

      Float_t low = ((fgkMass[ires] != 0.0) ? (fgkMass[ires] - fThresholdMult[ires] * 0.5) * (fgkMass[ires] - fgkMassThreshold[ires] * 0.5) : 0.0);
      Float_t high = (fgkMass[ires] + fgkMassThreshold[ires] * 0.5) * (fgkMass[ires] + fgkMassThreshold[ires] * 0.5);

      if ((low < mass) && (mass < high)) {
        itcouldbe = kTRUE;
      }
    }
  } else if (fgkChildMass[0][ires] != fgkChildMass[1][ires]) {
    /* switch masses hypothesis */
    mass = GetSquaredInvMassCheap(pt1, eta1, phi1, pt2, eta2, phi2, fgkChildMass[1][ires], fgkChildMass[0][ires]);

    if (TMath::Abs(mass - fgkMass[ires] * fgkMass[ires]) < 5 * fgkMassThreshold[ires]) {
      if (fill)
        fhResonanceRoughMasses->Fill(ires, TMath::Sqrt(mass));
      mass = GetSquaredInvMass(pt1, eta1, phi1, pt2, eta2, phi2, fgkChildMass[1][ires], fgkChildMass[0][ires]);

      Float_t low = ((fgkMass[ires] != 0.0) ? (fgkMass[ires] - fThresholdMult[ires] * 0.5) * (fgkMass[ires] - fgkMassThreshold[ires] * 0.5) : 0.0);
      Float_t high = (fgkMass[ires] + fgkMassThreshold[ires] * 0.5) * (fgkMass[ires] + fgkMassThreshold[ires] * 0.5);

      if ((low < mass) && (mass < high)) {
        itcouldbe = kTRUE;
      }
    }
  }
  if (itcouldbe)
    return mass;
  else
    return -1.0;
}

#endif // ALITWOPARTICLECORRELATIONSBASE_H
