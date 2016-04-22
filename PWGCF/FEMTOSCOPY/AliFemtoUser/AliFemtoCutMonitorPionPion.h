///
/// \file AliFemtoCutMonitorPionPion.h
///

#ifndef ALIFEMTOCUTMONITOR_PIONPION_H
#define ALIFEMTOCUTMONITOR_PIONPION_H

#pragma once

class TH1F;
class TH1I;
class TH2F;
class TH2I;
class AliFemtoEvent;

#include "AliFemtoCutMonitor.h"
#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoAnalysisPionLambda.h"


/// \namespace AliFemtoCutMonitorPionPion
/// \brief A set of standard cut monitors used by the pion-pion femtoscopy
///        analysis
///
/// This namespace contains cut monitors used by the AliFemtoAnalysisPionPion
/// class.
///
/// Because of the specialization for pion+pion particles, these classes
/// are not expected to be used in any other analysis, but may be used as a
/// model for other specializations.
///
/// Currently this namespace contains three cut monitors:
///   * Event
///   * Particle
///   * Pair
///
/// \author Andrew Kubera <andrew.michael.kubera@cern.ch>
///
namespace AliFemtoCutMonitorPionPion {

  /// \class AliFemtoCutMonitorPionPion::Event
  /// \brief Event cut monitor for PionPion analysis
  class Event : public AliFemtoCutMonitor {
  public:

    /**
     * Construct event cut monitor with knowledge if a passing or failing cut.
     *
     * \param passing: This will set the correct pass/fail word in histogram
     *                 titles.
     *
     * \param suffix_output: If true, this will put a _P or _F at the end of
     *   all the histogram names - required if not grouping output lists.
     */
    Event(const bool passing,
          const bool is_identical_analysis=kFALSE,
          const bool is_mc_analysis=kFALSE,
          const bool suffix_output=kFALSE);

    /**
     * Return list of Histograms to be placed in output file.
     */
    virtual TList* GetOutputList();

    /**
     * Called at beginning of event processing. Used to reset the collection
     * size counters.
     */
    virtual void EventBegin(const AliFemtoEvent*);
    virtual void EventEnd(const AliFemtoEvent*);

    /// Fill plots with event information
    virtual void Fill(const AliFemtoEvent* aEvent);

    /// Save information about the particle collection
    virtual void Fill(const AliFemtoParticleCollection *,
                      const AliFemtoParticleCollection *);

  protected:

    TH2F *fCentMult;
    TH1F *fVertexZ;
    TH2F *fVertexXY;

    TH2I *_collection_size_pass;
    TH2I *_collection_size_fail;

    TH1I *_identical_collection_size_pass;
    TH1I *_identical_collection_size_fail;

    const AliFemtoEvent *_prev_ev;
    UInt_t _prev_pion_coll_1_size;
    UInt_t _prev_pion_coll_2_size;
  };


  /// \class AliFemtoCutMonitorPionPion::Pion
  /// \brief Cut monitor containing plots of interesting pion data
  ///
  class Pion : public AliFemtoCutMonitor {
  public:
    Pion(const bool passing,
         const TString& typestr,
         const bool is_mc_analysis=kFALSE,
         const bool suffix_output=kFALSE);
    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoTrack* aEvent);

  protected:

    TH2F *fYPt;
    TH2F *fPtPhi;
    TH2F *fEtaPhi;
    TH2F *fChiTpcIts;
    TH2F *fdEdX;
    TH2F *fImpact;
    TH1F *fMinv;
  };

  /// \class AliFemtoCutMonitorPionPion::Pair
  /// \brief Cut monitor containing plots of interesting pair data
  ///
  class Pair : public AliFemtoCutMonitor {
  public:
    Pair(const bool passing,
         const TString& typestr,
         const bool is_mc_analysis=kFALSE,
         const bool suffix_output=kFALSE);

    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoPair* aEvent);

  protected:

    TH1F *fMinv;
    TH1F *fKt;
    TH2F *fDetaDphi;

    TH2F *fMCTrue_minv;
    TH2F *fMCTrue_kstar;
  };
};

#endif
