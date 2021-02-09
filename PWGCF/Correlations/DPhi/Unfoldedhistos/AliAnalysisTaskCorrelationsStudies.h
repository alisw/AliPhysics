/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliAnalysisTaskCorrelationsStudies.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKCORRELATIONSSTUDIES_H
#define ALIANALYSISTASKCORRELATIONSSTUDIES_H

class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH3I;
class THn;
class TProfile3D;
class TObjArray;
class TList;
class TRandom3;

#include <THnSparse.h>

#include "AliCSEventCuts.h"
#include "AliCSTrackSelection.h"
#include "AliAnalysisTaskSE.h"
#include "AliDptDptCorrelations.h"
#include "AliCSPairAnalysis.h"

/// class AliAnalysisTaskCorrelationsStudies
/// \brief Analysis task for two particle correlations
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Oct 17, 2016



class AliAnalysisTaskCorrelationsStudies : public AliAnalysisTaskSE {
 public:
                                AliAnalysisTaskCorrelationsStudies();
                                AliAnalysisTaskCorrelationsStudies(const char *name);
    virtual                    ~AliAnalysisTaskCorrelationsStudies();

    virtual void                UserCreateOutputObjects();
    virtual void                NotifyRun();
    virtual void                UserExec(Option_t *option);
    virtual void                FinishTaskOutput();
    virtual void                Terminate(Option_t *);


    Bool_t                      Configure(const char *confstring);
    Bool_t                      ConfigureCorrelations(const char *confstring, const char *weighthistpattern, TString &szContainerPrefix);
    Bool_t                      ConfigureCorrelationsBinning(const char *confstring);

 private:
    const char                 *ProduceConfigurationString();
    virtual void                ProcessTracks(Bool_t simactive);
    virtual void                ProcessTrueTracks();
    Bool_t                      BuildEfficiencyProfiles();
    void                        BuildTrueRecRelation();
    void                        BuildTrueRecAccRelation();
    Int_t                       GetNoOfTruePrimaries();
    Int_t                       GetNoOfTrueParticles();
    Int_t                       GetNoOfMCRecTracks();

 private:
    /// \enum mcRecOptions
    /// \brief The options for additional rec analysis
    ///
    /// There is an additional analysis able to be done with
    /// reconstructed tracks. These are the options for it
    enum mcRecOptions {
      kNone = 0,                                              ///< no additional analysis with reconstructed tracks
      kRecWithTrue,                                           ///< additional analysis with rec tracks using true values
      kRecTruePrimaries,                                      ///< additional analysis with rec tracks but using only the ones corresponding to true primaries
      kRecWithNotAccepted,                                    ///< additional analysis with rec tracks incorporating the reconstructed which were not accepted
      kRecTruePrimariesWithTrue                               ///< additional analysis with rec tracks using only the ones corresponding to true primaries and with their true values
    };

    /// \enum trackQualityFlags
    /// \brief Flags for keeping track of track quality
    enum trackQualityFlags {
      kAccepted =           0x0001,                           ///< the track has been accepted
      kNoTrack =            0x0002,                           ///< the track has been accepted
      kFromInjectedSignal = 0x0004,                           ///< the track has been accepted
      kSamePositiveLabel =  0x0010,                           ///< the track has been accepted
      kNegativeLabel =      0x0020,                           ///< the track has been accepted
      kSameNegativeLabel =  0x0040                            ///< the track has been accepted
    };

    /// \enum onTrueEfficiencyOptions
    /// \brief Options for applying efficiency criteria on true tracks
    enum onTrueEfficiencyOptions {
      kNoEfficiencyOnTrue = 0,                                ///< no efficiency on true tracks, consider all of them within the acceptance
      kEfficiencyProfile  = 1,                                ///< apply an efficiency profile
      kOnlyReconstructed  = 2                                 ///< only consider true tracks which have been reconstructed
    };

    TList                      *fOutput;                      ///< Output histograms list
    TString                     fTaskConfigurationString;     ///< the task configuration string
    TString                     fTaskActivitiesString;        ///< the activities to perform by the task
    TString                     fEventCutsString;             ///< the event cuts string of characters
    TString                     fEfficiencyProfileToEnforce;  ///< the efficiency profile to enforce
    TString                     fEfficiencyProfileOnTrue;     ///< the efficiency profile to enforce on true data
    TString                     fAdditionalMCRecOption;       ///< the option for the additional MC rec results
    TString                     fTrackSelectionString;        ///< the track selection cuts string
    AliCSEventCuts             *fEventCuts;                   ///< the event cuts
    AliCSTrackSelection        *fTrackSelectionCuts;          ///< track cuts to select tracks

    AliCSAnalysisCutsBase::ProdPeriods      fDataPeriod;      ///< the current data period under analysis

    Bool_t                      fEnforceEfficiencyProfile;    ///< enforce the efficiency profile, only valid for MC analysis
    onTrueEfficiencyOptions     fOnTrueEfficiencyProfile;     ///< enforce efficiency on true data, only valid for MC analysis
    Bool_t                      fCorrectOnTrueEfficiency;     ///< use the efficiency profile on true data for efficiency correction on true
    Bool_t                      fDoProcessCorrelations;       ///< perform the correlation analysis data collection
    Bool_t                      fDoProcessPairAnalysis;       ///< perform the pairs analysis data collection
    mcRecOptions                fMCRecOption;                 ///< options for the additional MC rec results
    AliDptDptCorrelations       fProcessCorrelations;         ///< the processing correlations instance
    AliDptDptCorrelations       fProcessMCRecCorrelationsWithOptions; ///< the processing correlations instance for MC rec with additional options
    AliDptDptCorrelations       fProcessTrueCorrelations;     ///< the processing correlations instance for MC truth
    AliCSPairAnalysis           fProcessPairAnalysis;         ///< the processing pairs analysis instance
    AliCSPairAnalysis           fProcessTruePairAnalysis;     ///< the processing pairs analysis instance for MC truth
    TF1                        *ffEfficiencyProfile;          ///< the formula for the efficiency profile to enforce
    TRandom3                   *fRandomGenerator;             ///< the random generator object instance for enforcing the efficiency profile
    TH3F                       *fhWeightsTrack_1;             ///< the weights histogram for track one
    TH3F                       *fhWeightsTrack_2;             ///< the weights histogram for track two
    TH1F                       *fhEffCorrTrack_1;             ///< the efficiency correction for track one
    TH1F                       *fhEffCorrTrack_2;             ///< the efficiency correction for track two
    THn                        *fhPairEfficiency_PP;          ///< the plus plus pair efficiency
    THn                        *fhPairEfficiency_PM;          ///< the plus minus pair efficiency
    THn                        *fhPairEfficiency_MM;          ///< the minus minus pair efficiency
    THn                        *fhPairEfficiency_MP;          ///< the minus plus pair efficiency
    TObjArray                  *fPositiveTrackPdf;            ///< the distribution of positive tracks for simulation
    TObjArray                  *fNegativeTrackPdf;            ///< the distribution of negative tracks for simulation
    TObjArray                  *fTrueToRec;                   ///< the true to reconstructed association, via label
    TObjArray                  *fTrueToRecWrong;              ///< the true to reconstructed association, via negative label
    UInt_t                     *fMCRecFlags;                  //!<! the flags which qualify MC reconstructed tracks
    UInt_t                     *fMCTruePrimaryFlags;          //!<! the flags which qualify MC primary true tracks
    Int_t                       fMCFlagsStorageSize;          ///< the size of the MC flags storage

    TH2F                       *fhOnTrueEfficiencyProfile_1;  ///< the histogram for the efficiency profile on true data for track one
    TH2F                       *fhOnTrueEfficiencyProfile_2;  ///< the histogram for the efficiency profile on true data for track two
    TH1F                       *fhPt;                         ///< \f$ p_{T} \f$ spectrum
    TH1F                       *fhPtPos;                      ///< \f$ p_{T} \f$ spectrum for positive tracks
    TH1F                       *fhPtNeg;                      ///< \f$ p_{T} \f$ spectrum for negative tracks
    TH1F                       *fhP;                          ///< \f$ p \f$ spectrum
    TH1F                       *fhPPos;                       ///< \f$ p \f$ spectrum for positive tracks
    TH1F                       *fhPNeg;                       ///< \f$ p \f$ spectrum for negative tracks
    TH1F                       *fhTruePt;                     ///< True \f$ p_{T} \f$ spectrum
    TH1F                       *fhTruePtPos;                  ///< True \f$ p_{T} \f$ spectrum for positive tracks
    TH1F                       *fhTruePtNeg;                  ///< True \f$ p_{T} \f$ spectrum for negative tracks
    TH1F                       *fhTrueP;                      ///< True \f$ p \f$ spectrum
    TH1F                       *fhTruePPos;                   ///< True \f$ p \f$ spectrum for positive tracks
    TH1F                       *fhTruePNeg;                   ///< True \f$ p \f$ spectrum for negative tracks
    TH1F                       *fhPurePt;                     ///< \f$ p_{T} \f$ spectrum with good PID
    TH1F                       *fhPurePtPos;                  ///< \f$ p_{T} \f$ spectrum for positive tracks with good PID
    TH1F                       *fhPurePtNeg;                  ///< \f$ p_{T} \f$ spectrum for negative tracks with good PID
    TH1F                       *fhPureP;                      ///< \f$ p \f$ spectrum with good PID
    TH1F                       *fhPurePPos;                   ///< \f$ p \f$ spectrum for positive tracks with good PID
    TH1F                       *fhPurePNeg;                   ///< \f$ p \f$ spectrum for negative tracks with good PID
    TH1F                       *fhUnConstrainedPt;            ///< not constrained \f$ p_{T} \f$ spectrum
    TH1F                       *fhPtDifference;               ///< Constrained - unconstrained \f$ p_{T} \f$ spectrum
    TH1F                       *fhEtaB;                       ///< pseudorapidity spectrum before any cut
    TH1F                       *fhEtaA;                       ///< pseudorapidity spectrum after cuts
    TH1F                       *fhTrueEta;                    ///< True pseudorapidity spectrum
    TH1F                       *fhPhiB;                       ///< \f$ \phi \f$ before any cut
    TH1F                       *fhPhiA;                       ///< \f$ \phi \f$ after all cuts
    TH1F                       *fhTruePhi;                    ///< True \f$ \phi \f$
    TH2F                       *fhEtaVsPhiB;                  ///< \f$ \eta \f$ vs \f$ \phi \f$ before any cut
    TH2F                       *fhEtaVsPhiA;                  ///< \f$ \eta \f$ vs \f$ \phi \f$ after all cuts
    TH2F                       *fhTrueEtaVsPhi;               ///< True \f$ \eta \f$ vs \f$ \phi \f$
    TH2F                       *fhPtVsEtaB;                   ///< \f$ p_{T} \f$ vs \f$ \eta \f$ before any cut
    TH2F                       *fhPtVsEtaA;                   ///< \f$ p_{T} \f$ vs \f$ \eta \f$ after all cuts
    TH2F                       *fhTruePtVsEta;                ///< True\f$ p_{T} \f$ vs \f$ \eta \f$
    TProfile3D                 *fhPt3DB;                      ///< \f$ \Sigma p_{T} \f$ vs  \f$ \eta \f$, \f$ \phi \f$ and \f$ vtx_{z} \f$ before any cut
    TProfile3D                 *fhPt3DA;                      ///< \f$ \Sigma p_{T} \f$ vs  \f$ \eta \f$, \f$ \phi \f$ and \f$ vtx_{z} \f$ after all cuts
    TProfile3D                 *fhTruePt3D;                   ///< True \f$ \Sigma p_{T} \f$ vs  \f$ \eta \f$, \f$ \phi \f$ and \f$ vtx_{z} \f$
    TH3I                       *fh3Dn1B;                      ///< track density vs \f$ \eta \f$, \f$ \phi \f$ and \f$ p_{T} \f$ before any cut
    TH3I                       *fh3Dn1A;                      ///< track density vs \f$ \eta \f$, \f$ \phi \f$ and \f$ p_{T} \f$ after all cuts
    TH2I                       *fhLambdaVsMultiplicity;       ///< track lambda angle vs event multiplicity
    THnSparseI                 *fhAcceptedVsMultiplicity;     ///< number of accepted tracks vs event multiplicity (rec tracks)
    TH2I                       *fhTrueLambdaVsPrimaries;      ///< True track lambda angle vs primary tracks
    THnSparseI                 *fhTrueAcceptedVsPrimaries;    ///< True number of accepted tracks vs primary tracks
    static const Int_t          kgTHnDimension;               ///< The dimension of the THn histograms
    THnSparseI                 *fhTrueRecNotAccepted;         ///< True histograms of reconstructed tracks not accepted
    THnSparseI                 *fhTrueNotReconstructed;       ///< True histograms of not reconstructed tracks
    TH1I                       *fhTrueMultiRec;               ///< Number of true multi reconstructed tracks histogram
    TH1I                       *fhPtTrueMultiRec;             ///< pT of true multi reconstructed tracks histogram

    // NEW HISTO to be declared here

    /// Copy constructor
    /// Not allowed. Forced private.
    AliAnalysisTaskCorrelationsStudies(const AliAnalysisTaskCorrelationsStudies&); // not implemented
    /// Assignment operator
    /// Not allowed. Forced private.
    /// \return l-value reference object
    AliAnalysisTaskCorrelationsStudies& operator=(const AliAnalysisTaskCorrelationsStudies&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskCorrelationsStudies, 4);
    /// \endcond
};

#endif /* ALIANALYSISTASKCORRELATIONSSTUDIES_h */

