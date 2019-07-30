#ifndef ALICSANALYSISCUTSBASE_H
#define ALICSANALYSISCUTSBASE_H

/// \file AliCSAnalysisCutsBase.h
/// \brief Base analysis cuts support for the correlations studies tasks
///
/// Based on ideas taken from Gamma Conversion analysis under PWGGA and from AliESDTrackCuts

#include <TBits.h>
#include <TNamed.h>

/// \class AliCSAnalysisCutsBase
/// \brief Base class for correlation studies analysis cuts
///
/// Provides support for configuring a variable set of cuts
/// via a string which can be then incorporated in all histograms
/// for cut information persistence.
///
/// Based on Gamma Conversion analysis implementation within PWGGA
///
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Oct 17, 2016

class TClonesArray;
class AliInputEventHandler;
class AliMCEventHandler;

class AliCSAnalysisCutsBase : public TNamed {
public:

  /// \enum QALevel
  /// \brief The level of the QA histograms output
  enum QALevel {
    kQALevelNone,     ///< no QA histograms output
    kQALevelLight,    ///< light QA histograms output
    kQALevelHeavy     ///< full QA histograms output
  };

  /// \enum ProdPeriods
  /// \brief The supported LHC production periods
  enum ProdPeriods {
    kNoPeriod=0,      ///< no period
    /// \name 2010
    /// \brief 2010 periods
    ///@{
    kLHC10bg,         ///< pp 7 TeV (LHC10c incl 900 GeV)
    kLHC10h,          ///< PbPb 2.76TeV
    ///@}
    /// \name 2010MC
    /// \brief MC's corresponding to 2010 data
    ///@{
    kLHC10d1,         ///< anchored LHC10b pass 2
    kLHC10d2,         ///< anchored LHC10b pass 2
    kLHC10d4a,        ///< anchored LHC10c pass 2
    kLHC10d4,         ///< anchored LHC10c pass 2
    kLHC10e12,        ///< anchored LHC10c pass 2
    kLHC10e13,        ///< anchored LHC10c pass 2
    kLHC10f6a,        ///< anchored LHC10d pass 2
    kLHC10f6,         ///< anchored LHC10d pass 2
    kLHC10e20,        ///< anchored LHC10e pass 2
    kLHC10e21,        ///< anchored LHC10e pass 2
    kLHC14j4,         ///< anchored LHC10[b-g] pass 4
    kLHC11a10a,       ///< anchored LHC10h pass 2
    kLHC11a10b,       ///< anchored LHC10h pass 2
    kLHC13d2,         ///< anchored LHC10h pass 2
    kLHC13d2b,        ///< anchored LHC10h pass 2
    kLHC12a11a,       ///< anchored LHC10h pass 2
    kLHC12a11b,       ///< anchored LHC10h pass 2
    kLHC12a11c,       ///< anchored LHC10h pass 2
    kLHC12a11d,       ///< anchored LHC10h pass 2
    kLHC12a11e,       ///< anchored LHC10h pass 2
    kLHC12a11f,       ///< anchored LHC10h pass 2
    kLHC13d19,        ///< anchored LHC10h pass 2, upgrade 5.5TeV PbPb
    ///@}
    /// \name 2011
    /// \brief 2011 periods
    ///@{
    kLHC11a,          ///< pp 2.76TeV (part 7TeV)
    kLHC11b,          ///< pp 7TeV
    kLHC11cg,         ///< pp 7TeV
    kLHC11h,          ///< PbPb 2.76TeV
    ///@}
    /// \name 2011MC
    /// \brief MC's corresponding to 2011 data
    ///@{
    kLHC12a15c,       ///< anchored LHC11a pass 2 - JJ
    kLHC12f1a,        ///< anchored LHC11a pass 4
    kLHC12f1b,        ///< anchored LHC11a pass 4
    kLHC12i3,         ///< anchored LHC11a pass 4
    kLHC15g1a,        ///< anchored LHC11a pass 4 - JJ
    kLHC15g1b,        ///< anchored LHC11a pass 4 - JJ
    kLHC13e4,         ///< anchored LHC11c pass 1 - GJ
    kLHC13e5,         ///< anchored LHC11c pass 1 - JJ
    kLHC14k1a,        ///< anchored LHC11[c-d] pass 1 - JJ
    kLHC14k1b,        ///< anchored LHC11[c-d] pass 1 - JJ
    kLHC12a15f,       ///< anchored LHC11d pass 1 - JJ
    kLHC12a15g,       ///< anchored LHC11d pass 1 - GJ
    kLHC12f2a,        ///< anchored LHC11d pass 1 - JJ
    kLHC14a1a,        ///< anchored LHC11h pass 2
    kLHC14a1b,        ///< anchored LHC11h pass 2
    kLHC14a1c,        ///< anchored LHC11h pass 2
    ///@}
    /// \name 2012
    /// \brief 2012 periods
    ///@{
    kLHC12,           ///< pp 8TeV
    ///@}
    /// \name 2012MC
    /// \brief MC's corresponding to 2012 data
    ///@{
    kLHC14e2a,        ///< anchored LHC12[a-h] pass 1
    kLHC14e2b,        ///< anchored LHC12[a-h] pass 1
    kLHC14e2c,        ///< anchored LHC12[a-h] pass 1
    kLHC15h1,         ///< anchored LHC12[a-h] pass 2
    kLHC15h2,         ///< anchored LHC12[a-h] pass 2
    kLHC16c2,         ///< anchored LHC12[a-h] pass 2 - JJ
    ///@}
    /// \name 2013
    /// \brief 2013 periods
    ///@{
    kLHC13bc,         ///< pPb 5.023TeV
    kLHC13de,         ///< pPb 5.023TeV
    kLHC13f,          ///< Pbp 5.023TeV
    kLHC13g,          ///< pp 2.76TeV
    ///@}
    /// \name 2013MC
    /// \brief MC's corresponding to 2013 data
    ///@{
    kLHC13b2,         ///< anchored LHC13[b-c] pass 2
    kLHC13b2_efix,    ///< anchored LHC13[b-c] pass 2
    kLHC13e7,         ///< anchored LHC13[b-c] pass 2
    kLHC14b2,         ///< anchored LHC13[b-c] pass 2
    kLHC13b4_fix,     ///< anchored LHC13[b-c] pass 2 - JJ
    kLHC13b4_plus,    ///< anchored LHC13[b-c] pass 2 - JJ
    kLHC16c3a,        ///< anchored LHC13[d-e] pass 2 - JJ
    kLHC16c3b,        ///< anchored LHC13[d-e] pass 2 - JJ
    kLHC16c3c,        ///< anchored LHC13[d-e] pass 2 - GJ
    kLHC15g2,         ///< anchored LHC13g pass 1
    kLHC15a3a,        ///< anchored LHC13g pass 1 - JJ
    kLHC15a3a_plus,   ///< anchored LHC13g pass 1 - JJ
    kLHC15a3b,        ///< anchored LHC13g pass 1 - JJ
    kLHC15d3a,        ///< anchored LHC13g pass 1
    kLHC15d3b,        ///< anchored LHC13g pass 1
    ///@}
    /// \name 2015
    /// \brief 2015 periods
    ///@{
    kLHC15fm,         ///< pp 13 TeV
    kLHC15n,          ///< pp 5 TeV
    kLHC15oLIR,       ///< PbPb 5 TeV low interaction rate
    kLHC15oHIR,       ///< PbPb 5 TeV high interaction rate
    ///@}
    /// \name 2015MC
    /// \brief MC's corresponding to 2015 data
    ///@{
    kLHC15g3a3,       ///< anchored LHC15f pass 1
    kLHC15g3a,        ///< anchored LHC15f pass 1
    kLHC15g3c2,       ///< anchored LHC15f pass 1
    kLHC15g3c3,       ///< anchored LHC15f pass 1
    kLHC15g3,         ///< anchored LHC15f pass 1
    kLHC16a2a,        ///< anchored LHC15h pass 1
    kLHC16a2b,        ///< anchored LHC15h pass 1
    kLHC16a2c,        ///< anchored LHC15h pass 1
    kLHC15l1a2,       ///< anchored LHC15n pass 1
    kLHC15l1b2,       ///< anchored LHC15n pass 1
    kLHC15k1a1,       ///< LHC15o low IR
    kLHC15k1a2,       ///< LHC15o low IR
    kLHC15k1a3,       ///< LHC15o low IR
    kLHC16g1,         ///< anchored LHC15o pass1 - general purpose
    kLHC16g1a,        ///< anchored LHC15o pass1 - general purpose 0-10%
    kLHC16g1b,        ///< anchored LHC15o pass1 - general purpose 10-50%
    kLHC16g1c,        ///< anchored LHC15o pass1 - general purpose 50-90%
    kLHC16h4a,        ///< anchored LHC15o pass1 - injected signals 0-10%
    kLHC16h4b,        ///< anchored LHC15o pass1 - injected signals 10-50%
    kLHC16h4b2,       ///< anchored LHC15o pass1 - injected signals 10-50% with TPC gas corrections
    kLHC16h4c,        ///< anchored LHC15o pass1 - injected signals 50-90%
    kLHC16h2a,        ///< anchored LHC15o pass1 - jet-jet 0-10%
    kLHC16h2b,        ///< anchored LHC15o pass1 - jet-jet 10-50%
    kLHC16h2c,        ///< anchored LHC15o pass1 - jet-jet 50-90%
    kLHC16i3,         ///< anchored LHC15o HIJING + HF injected with electron decays
    kLHC17i2,         ///< anchored LHC15o AMPT via AliGenerators
    ///@}
    /// \name 2016
    /// \brief 2016 periods
    ///@{
    kLHC16k,          ///< pp 13 TeV
    kLHC16l,          ///< pp 13 TeV
    ///@}
    /// \name 2017
    /// \brief 2017 periods
    ///@{
    kLHC17n,          ///< XeXe 5.44 TeV
    ///@}
    /// \name 2017MC
    /// \brief MC's corresponding to 2017 data
    ///@{
    kLHC17j6,          ///< anchored LHC17n XeXe 5.44 TeV with ITS reco points
    kLHC17j7,          ///< anchored LHC17n XeXe 5.44 TeV general purpose
    ///@}
    /// \name 2018
    /// \brief 2018 periods
    ///@{
    kLHC18q,          ///< PbPb 5 TeV
    kLHC18r,          ///< PbPb 5 TeV
    ///@}
    /// \name 2018MC
    /// \brief MC's corresponding to 2018 data
    ///@{
    kLHC18l8,          ///< anchored to LHC18q/r Pb-Pb 5.02 TeV general purpose
    ///@}
    /// \name fast MC productions
    ///@{
    kLHC13f3,         ///<  PbPb, AMPT, fast generation, 2.76TeV (min. bias)
    ///@}
  };

  /// \enum EnergyValue
  /// \brief The collision energy for the production period
  enum EnergyValue {
    kUnset        = 0,  ///< not defined
    k900GeV       = 1,  ///< pp 900 GeV
    k2760GeV      = 2,  ///< pp 2.76TeV
    k5TeV         = 3,  ///< pp 5 TeV
    k7TeV         = 4,  ///< pp 7 TeV
    k8TeV         = 5,  ///< pp 8 TeV
    k13TeV        = 6,  ///< pp 13 TeV
    kpPb5TeV      = 7,  ///< pPb 5 TeV
    kpPb8TeV      = 8,  ///< pPb 8 TeV
    kPbPb2760GeV  = 9,  ///< PbPb 2.76TeV
    kPbPb5TeV     = 10, ///< PbPb 5 TeV
    kXeXe5440GeV  = 11, ///< XeXe 5.44 TeV
  };

                                AliCSAnalysisCutsBase();
                                AliCSAnalysisCutsBase(Int_t nCuts, Int_t nParams, const char *name="CS AnalysisCuts", const char * title="CS AnalysisCuts");
  virtual                      ~AliCSAnalysisCutsBase();

                                /// Gets the current activated cuts
                                /// \return the mask with the activated cuts
  const TBits                  *GetCutsActivatedMask() { return &fCutsActivatedMask; }

                                /// Initializes the cuts
                                /// \param name the name to assign to the histograms list
  virtual void                  InitCuts(const char *name) = 0;
                                /// Sets the desired level for the QA histograms output
                                /// \param level the desired QA histograms output level
  void                          SetQALevelOutput(QALevel level) { fQALevel = level; }
                                /// Get the histograms list
                                /// \return the histograms list
  virtual TList                *GetHistogramsList() { return fHistogramsList; }

                                /// Processes a potential change in the run number
                                /// Pure virtual function
  virtual void                  NotifyRun() = 0;
                                /// A new event is coming
                                /// Pure virtual function
  virtual void                  NotifyEvent() = 0;
  static void                   NotifyRunGlobal();
                                /// Prints the activated cuts mask
                                /// \param opt the print options
  virtual void                  Print(Option_t *opt = "") const { fCutsActivatedMask.Print(opt); }

  Bool_t                        InitializeCutsFromCutString(const TString cutSelectionString);
                                /// Get the event cuts associated string
                                /// \return the event cuts associated string
  const char                   *GetCutsString() const { return fCutsString; }
                                /// Get the period name corresponding to the current analysis
                                /// \return the period name of the current analysis
  static const char            *GetPeriodName() { return fgPeriodName.Data(); }
                                /// Get the period corresponding to the current analysis
                                /// \return the period code of the current analysis
  static ProdPeriods            GetGlobalPeriod() { return fgDataPeriod; }
                                /// Get the anchor period corresponding to the current analysis
                                /// \return the anchor period code of the current analysis
  static ProdPeriods            GetGlobalAnchorPeriod() { return fgAnchorPeriod; }
                                /// Is a Monte Carlo data set
                                /// \return kTRUE if it is a Monte Carlo data set, kFALSE otherwise
  static Bool_t                 IsMC() { return fgIsMC; }
                                /// Is a fast Monte Carlo data set (only MC truth)
                                /// \return kTRUE if it is a fast Monte Carlo data set, kFALSE otherwise
  static Bool_t                 IsMConlyTruth() { return fgIsMConlyTruth; }
                                /// Get the Monte Carlo event handler
                                /// \return the Monte Carlo event handler
  static AliMCEventHandler     *GetMCEventHandler() { return fgMCHandler; }
                                /// Get the Input event handler
                                /// \return the Input event handler
  static AliInputEventHandler  *GetInputEventHandler() { return fgInputHandler; }
                                /// Get the MC array of true particles
                                /// \return the MC array of true particles
  static TClonesArray          *GetMCTrueArray() { return fgMCArray; }

protected:
  Bool_t                        UpdateCutsString();
                                /// Set the value for a concrete cut
                                /// Pure virtual function
                                /// \param paramID the ID of the cut of interest
                                /// \param value the value to assign to the cut
                                /// \return kTRUE if the cut value was accepted
  virtual Bool_t                SetCutAndParams (Int_t paramID, Int_t value) = 0;
  void                          PrintCutsWithValues() const ;
                                /// Print the cut with its value
                                /// Pure virtual function
                                /// \param paramID the ID of the cut of interest
  virtual void                  PrintCutWithParams(Int_t paramID) const = 0;

private:
  static TString                GetPeriodNameFromDataFilePath();
  static Int_t                  GetCurrentRunNumber();

private:
  static TString                fgPeriodName;           ///< the period name of the ongoing analysis
  static ProdPeriods            fgDataPeriod;           ///< the global period of ongoing analysis
  static ProdPeriods            fgAnchorPeriod;         ///< the anchor period, different from #fgDataPeriod in case of MC data
  Int_t                         fNParams;               ///< the number of cuts parameters
  char                         *fCutsString;            ///< the cuts parameters associated string
protected:
  Int_t                         fNCuts;                 ///< the number of cuts
  QALevel                       fQALevel;               ///< the level of the QA histograms output
                                /// the external values for each of the cuts parameters
  Int_t                        *fParameters;            //[fNParams]
  TBits                         fCutsEnabledMask;       ///< the mask of enabled cuts
  TBits                         fCutsActivatedMask;     ///< the mask of cut activated for the ongoing event
  ProdPeriods                   fDataPeriod;            ///< the current period under analysis. Occasionally could be different from the global period
  static EnergyValue            fgEnergy;               ///< the collision energy for the analysis period
  static Bool_t                 fgIsMC;                 ///< MC flag from production information
  static Bool_t                 fgIsMConlyTruth;        ///< fast MC flag (only kinematics) from production information
  static AliInputEventHandler  *fgInputHandler;         ///< the input handler for the current ongoing analysis
  static AliMCEventHandler     *fgMCHandler;            ///< MC handler for the current ongoing analysis
  static Bool_t                 fgIsESD;                ///< kTRUE if data format is ESD, kFALSE if AOD
  static TClonesArray          *fgMCArray;              ///< the array containing true particles when MC AOD analysis

  TList                        *fHistogramsList;        ///< the list of histograms used

  /// Copy constructor
  /// Not allowed. Forced private.
  AliCSAnalysisCutsBase(const AliCSAnalysisCutsBase&);
  /// Assignment operator
  /// Not allowed. Forced private.
  /// \return l-value reference object
  AliCSAnalysisCutsBase& operator=(const AliCSAnalysisCutsBase&);

  /// \cond CLASSIMP
  ClassDef(AliCSAnalysisCutsBase,2);
  /// \endcond
};

#endif /* ALICSANALYSISCUTSBASE_H */
