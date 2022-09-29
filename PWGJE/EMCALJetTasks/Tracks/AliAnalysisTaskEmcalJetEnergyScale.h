/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIANALYSISTASKEMCALJETENERGYSCALE_H
#define ALIANALYSISTASKEMCALJETENERGYSCALE_H

#include <exception>
#include <iosfwd>
#include <string>
#include <TString.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliJetContainer.h"
#include "AliVCluster.h"

class THistManager;
class TRandom;

namespace PWGJE {

namespace EMCALJetTasks{

class AngularityHandler : public TObject{
public:
  class AngularityBin : public TObject {
  public:
    AngularityBin(): TObject(), fMin(-1.), fMax(-1.), fValue(-1.) {}
    AngularityBin(double min, double max, double value) : TObject(), fMin(min), fMax(max), fValue(value) {}
    virtual ~AngularityBin() {}

    double Min() const { return fMin; }
    double Max() const { return fMax; }
    double Value() const { return fValue; }

    void SetMin(double min) { fMin = min; }
    void SetMax(double max) { fMax = max; }
    void SetValue(double value) { fValue = value; }

    bool IsInRange(double entry) const { return entry >= fMin && entry < fMax; }
    bool IsOK() const { return fMin > -1. && fMax > -1.; }

    void PrintStream(std::ostream &stream) const;
  private:
    double fMin;          ///< Bin min
    double fMax;          ///< Bin max
    double fValue;        ///< Bin value

    ClassDef(AngularityBin, 1);
  };

  class BinNotFoundException : public std::exception {
  public:
    BinNotFoundException(double pt) : fMessage(), fPt(pt) {
      fMessage = Form("No bin found for pt %f", pt);
    }
    virtual ~BinNotFoundException() throw() {}

    const char *what() const throw() {
      return fMessage.c_str();
    }

    double getPt() const {return fPt; }
  private:
    std::string fMessage;   ///< Message
    double      fPt;        ///< Pt
  };

  AngularityHandler() : TObject(), fRadius(-1), fBins() { fBins.SetOwner(true); }
  AngularityHandler(double radius) : TObject(), fRadius(radius), fBins() { fBins.SetOwner(true); }
  virtual ~AngularityHandler() {}

  void InitFromFile(const char *filename);
  void SetBin(double min, double max, double value);
  void SetRadius(double radius) { fRadius = radius; }

  double GetValue(double pt) const;
  bool IsHigherAngularity(double pt, double angularity) const;


  virtual void Print(Option_t * option = "") const;
  void PrintStream(std::ostream &stream) const ;

private:
  AngularityBin FindBin(double pt) const;
  double fRadius;       ///< Jet raduius
  TList fBins;          ///< Pt bins

  ClassDef(AngularityHandler, 1);
};

class AliAnalysisTaskEmcalJetEnergyScale : public AliAnalysisTaskEmcalJet {
public:
  class AngularityException : public std::exception{
    public:
      AngularityException() {}
      virtual ~AngularityException() throw() {}

      virtual const char *what() const throw() {
        return "Angularity cannot be detrmined";
      }
  };


  enum EJetTypeOutliers_t {
    kOutlierPartJet,
    kOutlierDetJet
  };
  AliAnalysisTaskEmcalJetEnergyScale();
  AliAnalysisTaskEmcalJetEnergyScale(const char *name);
  virtual ~AliAnalysisTaskEmcalJetEnergyScale();


  AliJetContainer *GetPartLevelJetContainer() const { return GetJetContainer(fNameParticleJets); }
  AliJetContainer *GetDetLevelJetContainer() const { return GetJetContainer(fNameDetectorJets); }
  AliTrackContainer *GetTracks() const { if(fNameTracks.Length()) return GetTrackContainer(fNameTracks); return NULL;}
  AliClusterContainer *GetClusters() const { if(fNameClusters.Length()) return GetClusterContainer(fNameClusters); return NULL; }
  AliMCParticleContainer *GetMCParticles() const { if(fNameMCParticles.Length()) return GetMCParticleContainer(fNameMCParticles); return NULL; }
  const TString &GetNamePartLevelLets() const { return fNameParticleJets; }
  const TString &GetNameDetLevelJets() const { return fNameDetectorJets; }
  AngularityHandler *GetAngularityHandler() const {return fAngularityHandler; }
  bool HasAngularitySplitting() const { return fAngularityHandler != NULL; }

  void SetNameDetJetContainer(const char *name)  { fNameDetectorJets = name; }
  void SetNamePartJetContainer(const char *name) { fNameParticleJets = name; }
  void SetNameTracks(const char *name)           { fNameTracks = name;}
  void SetNameClusters(const char *name)         { fNameClusters = name; }
  void SetNameMCParticles(const char *name)      { fNameMCParticles = name; }
  void SetTriggerName(const char *name)          { fTriggerSelectionString = name; }
  void SetFractionResponseClosure(double fraction) { fFractionResponseClosure = fraction; }
  void SetFillHSparse(Bool_t doFill)             { fFillHSparse = doFill; }
  void SetEnergyScaleShift(Double_t scaleshift)  { fScaleShift = scaleshift; }
  void SetUseStandardOutlierRejection(bool doUse) { fUseStandardOutlierRejection = doUse; }
  void SetDebugMaxJetOutliers(bool doDebug)      { fDebugMaxJetOutliers = doDebug; }
  void SetJetTypeOutlierCut(EJetTypeOutliers_t jtype) { fJetTypeOutliers = jtype; }
  void SetRequireSameAcceptance(Bool_t doRequire) { fRequireSameAcceptance = doRequire; }
  void SetDoBkgSubtraction(bool doBkg = true)             { fDoBkgSub = doBkg; }
  void SetAngularitySpitting(bool doSplit = true);

  void ConfigurePtHard(MCProductionType_t mcprodtype, const TArrayI &pthardbinning, Bool_t doMCFilter, Double_t jetptcut);
  void ConfigureMinBias(MCProductionType_t mcprodtype);
  void ConfigureJetSelection(Double_t minJetPtPart, Double_t minJetPtDet, Double_t maxTrackPtPart, Double_t maxTrackPtDet, Double_t maxClusterPt, Double_t minAreaPerc);

  //*** Standard add functions. No random cones background subtraction. ***//
  // Base add function called by both of the two add functions following it.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScaleBase(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      trigger,
    const char *                      nametrackcont,
    const char *                      nameclustercont,
    const char *                      suffix
  );

  // Standard add function. Calls the base add fuction above.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScale(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      trigger,
    const char *                      suffix
  );

  // Add function with embedded MC background rejection handling. Calls the base add function above.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScale(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      trigger,
    const char *                      nametrackcont,
    const char *                      nameclustercont,
    const char *                      suffix
  );


  //*** Add functions for use with random cones background subtraction. ***//
  // Base add function called by both of the two add functions following it.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScaleBkgSubBase(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      nRho,
    const char *                      nRhoMC,
    const char *                      trigger,
    const char *                      nametrackcont,
    const char *                      nameclustercont,
    const char *                      suffix
  );

  // Standard add function. Calls the base add fuction above.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScaleBkgSub(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      nRho,
    const char *                      nRhoMC,
    const char *                      trigger,
    const char *                      suffix
  );

  // Add function with embedded MC background rejection handling. Calls the base add function above.
  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScaleBkgSub(
    AliJetContainer::EJetType_t       jetType,
    AliJetContainer::ERecoScheme_t    recoscheme,
    AliVCluster::VCluUserDefEnergy_t  energydef,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      namepartcont,
    const char *                      nRho,
    const char *                      nRhoMC,
    const char *                      trigger,
    const char *                      nametrackcont,
    const char *                      nameclustercont,
    const char *                      suffix
  );

protected:
  virtual void UserCreateOutputObjects();
  virtual void UserRunBeforeEventSelection();
  virtual Bool_t Run();
  virtual Bool_t CheckMCOutliers();
  bool IsSelectEmcalTriggers(const TString &triggerstring) const;
  Double_t MakeAngularity(const AliEmcalJet &jet, AliVCluster::VCluUserDefEnergy_t energydef) const;

private:
  THistManager                *fHistos;                       //!<! Histogram collection
  AngularityHandler           *fAngularityHandler;            ///< Handler for angularity splitting
  TString                     fNameTracks;                    ///< Name of the container for tracks
  TString                     fNameClusters;                  ///< Name of the container for EMCAL clusters
  TString                     fNameMCParticles;               ///< Name of the container for MC particles
  TString                     fNameDetectorJets;              ///< Name of the data jet container
  TString                     fNameParticleJets;              ///< Name of the MC jet container
  TString                     fTriggerSelectionString;        ///< Trigger selection string
  TString                     fNameTriggerDecisionContainer;  ///< Global trigger decision container
  Double_t                    fFractionResponseClosure;       ///< Fraction of jets used for response in closure test
  Bool_t                      fFillHSparse;                   ///< Fill THnSparses
  Double_t                    fScaleShift;                    ///< Shift of the jet energy scale (fixed)
  Bool_t                      fRequireSameAcceptance;         ///< Require same acceptance type for det. level and part. level jet in response matrix
  Bool_t                      fUseStandardOutlierRejection;   ///< Use standard outlier rejection
  Bool_t                      fDebugMaxJetOutliers;           ///< Debug max jet determination for outlier rejection
  Bool_t                      fDoBkgSub;                      ///< Do background subtraction
  EJetTypeOutliers_t          fJetTypeOutliers;               ///< Jet type used for outlier detection
  TRandom                     *fSampleSplitter;               //!<! Sample splitter

  AliAnalysisTaskEmcalJetEnergyScale(const AliAnalysisTaskEmcalJetEnergyScale &);
  AliAnalysisTaskEmcalJetEnergyScale &operator=(const AliAnalysisTaskEmcalJetEnergyScale &);

  ClassDef(AliAnalysisTaskEmcalJetEnergyScale, 1);
};

}

}

std::ostream &operator<<(std::ostream &stream, const PWGJE::EMCALJetTasks::AngularityHandler &);
std::ostream &operator<<(std::ostream &stream, const PWGJE::EMCALJetTasks::AngularityHandler::AngularityBin &);
#endif // ALIANALYSISTASKEMCALJETENERGYSCALE_H
