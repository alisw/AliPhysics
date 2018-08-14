///
/// \file AliFemtoCutMonitorPionLambda.h
///

#pragma once

#ifndef ALIFEMTOCUTMONITOR_PION_LAMBDA_H
#define ALIFEMTOCUTMONITOR_PION_LAMBDA_H

class TH1F;
class TH2F;
class TH2I;
class AliFemtoEvent;

#include "AliFemtoCutMonitor.h"
#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoAnalysisPionLambda.h"


/// \namespace AliFemtoCutMonitorPionLambda
/// \brief The cut monitor used by the pion-lambda femtoscopy analysis
///
/// This cut monitor is designed to be used for ALL cuts in the
/// AliFemtoAnalysisPionLambda class. This is for conveience.
///
/// Because of the specialization for pion+lambda particles, this class
/// is not expected to be used in any other analysis, but may be used as a
/// model for other specializations.
///
namespace AliFemtoPionLambdaCutMonitor {

  class Event : public AliFemtoCutMonitor {
  public:

    /**
     * Construct event cut monitor with knowledge if a passing or failing cut.
     *
     * \param passing: This will set the correct pass/fail word in histogram
     *                 titles.
     *
     * \param suffix_output: If true, this will put a _P or _F at the end of all
     *   the histogram names - required if not grouping output lists.
     */
    Event(const bool passing,
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

  private:
    Event(const Event &);
    Event& operator=(const Event &);


  protected:
    TH1F *_centrality;
    TH1F *_multiplicity;
    TH1F *_vertex_z;
    TH2F *_vertex_xy;
    TH2I *_collection_size_pass;
    TH2I *_collection_size_fail;

    const AliFemtoEvent *_prev_ev;
    UInt_t _prev_pion_coll_size;
    UInt_t _prev_lam_coll_size;
  };


  /// \class AliFemtoPionLambdaCutMonitor::Pion
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

  private:
    Pion(const Pion &);
    Pion& operator=(const Pion &);

  protected:

    TH2F *fYPt;
    TH2F *fPtPhi;
    TH2F *fEtaPhi;
    TH2F *fdEdX;
    TH1F *fMinv;
  };

  /// \class AliFemtoPionLambdaCutMonitor::Lambda
  /// \brief Cut monitor containing plots of interesting lambda data
  ///
  class Lambda : public AliFemtoCutMonitor {
  public:

    /**
     * Construct the cut monitor
     * \param passing: Determines wheter to put (PASS) or (FAIL) to
     *    the output titles
     * \param typestr: String added to the python
     * \param ltype: Lambda type which determines which pos/neg track is pion/proton daughter.
     * \param suffix_output: If True, will append _P or _F to each histogram name
     */
    Lambda(const bool passing,
           const TString& typestr,
           const AliFemtoAnalysisPionLambda::LambdaType,
           const bool is_mc_analysis=kFALSE,
           const bool suffix_output=kFALSE);

    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoV0* aEvent);

  private:
    Lambda(const Lambda &);
    Lambda& operator=(const Lambda &);


  protected:
    AliFemtoAnalysisPionLambda::LambdaType fLambdaType;

    TH1F *_minv;
    TH2F *_ypt;
    TH2F *_dedx_p_pro;
    TH2F *_dedx_p_pi;
    TH1F *fCosPointingAngle;

    TH1F *fMCTrue_minv;
    TH2F *fMCTrue_ypt;

  };

  class Pair : public AliFemtoCutMonitor {
  public:
    Pair(const bool passing,
         const TString& typestr,
         const bool is_mc_analysis=kFALSE,
         const bool suffix_output=kFALSE);

    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoPair* aEvent);

  private:
    Pair(const Pair &);
    Pair& operator=(const Pair &);

  protected:

    TH1F *_minv;
    TH1F *fKt;
    TH1F *fAvgSep_pion;
    TH1F *fAvgSep_proton;

    TH2F *fMCTrue_minv;
    TH2F *fMCTrue_kstar;

  };

/*
  ///
  AliFemtoPionLambdaCutMonitor();

  virtual TList* GetOutputList();

  /// Fill the specific
  virtual void Fill(const AliFemtoEvent* aEvent, bool pass);
  virtual void Fill(const AliFemtoTrack* aTrack, bool pass);
  virtual void Fill(const AliFemtoV0* aV0, bool pass);
  virtual void Fill(const AliFemtoPair* aPair, bool pass);
*/

/*

protected:

  //------- Event Histograms -------
  TH1F *_event_centrality_pass,
       *_event_centrality_fail;

  TH1F *_event_multiplicity_pass,
       *_event_multiplicity_fail;

  TH1F *_event_vertex_z_pass,
       *_event_vertex_z_fail;

  TH2F *_event_vertex_xy_pass,
       *_event_vertex_xy_fail;

  //------- Pion Histograms -------
  TH1F *_pion_minv_pass,
       *_pion_minv_fail;

  TH2F *_pion_dedx_pt_pass,
       *_pion_dedx_pt_fail;

  TH2F *_pion_ypt_pass,
       *_pion_ypt_fail;

  //------- Lambda Histograms -------
  TH1F *_lambda_minv_pass,
       *_lambda_minv_fail;

  TH2F *_lambda_ypt_pass,
       *_lambda_ypt_fail;

  TH2F *_lambda_dedx_pt_pro_pass,
       *_lambda_dedx_pt_pro_fail;

  TH2F *_lambda_dedx_pt_pi_pass,
       *_lambda_dedx_pt_pi_fail;

  //------- Pair Histograms -------
  TH1F  *_pair_minv_pass,
        *_pair_minv_fail;

private:

  /// To only be called from constructors as a common histogram initializer
  void _init();
*/

// #ifdef __ROOT__
//   /// \cond CLASSIMP
//   ClassDef(AliFemtoCutMonitorPionLambda, 0);
//   /// \endcond
// #endif

};

#endif
