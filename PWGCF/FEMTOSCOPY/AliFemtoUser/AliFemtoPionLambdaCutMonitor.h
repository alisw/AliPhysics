///
/// \file AliFemtoCutMonitorPionLambda.h
///

#pragma once

#ifndef ALIFEMTOCUTMONITOR_PION_LAMBDA_H
#define ALIFEMTOCUTMONITOR_PION_LAMBDA_H

class TH1F;
class TH2F;

#include "AliFemtoCutMonitor.h"
#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoAnalysisPionLambda.h"


/// \class AliFemtoCutMonitorPionLambda
/// \brief The cut monitor used by the pion-lambda femtoscopy analysis
///
/// This cut monitor is designed to be used for ALL cuts in the
/// AliFemtoAnalysisPionLambda. Instead of adding
///
/// Because of the specialization - this
///
class AliFemtoPionLambdaCutMonitor {
public:

  class Event : public AliFemtoCutMonitor {
  public:
    Event(const bool passing);
    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoEvent* aEvent);

  protected:
    TH1F *_centrality;
    TH1F *_multiplicity;
    TH1F *_vertex_z;
    TH2F *_vertex_xy;
  };

  class Pion : public AliFemtoCutMonitor {
  public:
    Pion(const bool passing, const TString& typestr);
    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoTrack* aEvent);

  protected:
    TH1F *_minv;
    TH2F *_ypt;
    TH2F *fPtPhi;
    TH2F *fEtaPhi;
  };

  class Lambda : public AliFemtoCutMonitor {
  public:
    Lambda(const bool passing,
           const TString& typestr,
           const AliFemtoAnalysisPionLambda::LambdaType);
    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoV0* aEvent);

  protected:
    AliFemtoAnalysisPionLambda::LambdaType fLambdaType;

    TH1F *_minv;
    TH2F *_ypt;
    TH2F *_dedx_pt_pro;
    TH2F *_dedx_pt_pi;
  };

  class Pair : public AliFemtoCutMonitor {
  public:
    Pair(const bool passing, const TString& typestr);
    virtual TList* GetOutputList();
    virtual void Fill(const AliFemtoPair* aEvent);

  protected:
    TH1F *_minv;
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

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCutMonitorPionLambda, 0);
  /// \endcond
#endif

};

#endif
