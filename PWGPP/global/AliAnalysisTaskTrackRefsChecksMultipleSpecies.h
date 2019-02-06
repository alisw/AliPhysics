#ifndef ALIANALYSISTASKTRACKREFSCHECKSMULTIPLESPECIES
#define ALIANALYSISTASKTRACKREFSCHECKSMULTIPLESPECIES

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskTrackRefsChecksMultipleSpecies
//*************************************************************************

class TList;
class TNtuple;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TString;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTrackRefsChecksMultipleSpecies : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskTrackRefsChecksMultipleSpecies();
  virtual ~AliAnalysisTaskTrackRefsChecksMultipleSpecies();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);
  
 private:

  AliAnalysisTaskTrackRefsChecksMultipleSpecies(const AliAnalysisTaskTrackRefsChecksMultipleSpecies &source);
  AliAnalysisTaskTrackRefsChecksMultipleSpecies& operator=(const AliAnalysisTaskTrackRefsChecksMultipleSpecies &source);
  
  TList*  fOutput;                   //! list of output histos

  // positive tracks
  TH3F*** fhdtgl_pos;                     //! histo with difference of tgLambda between consecutive layers vs 1/pt vs phi
  TH3F*** fhdtgl_prod_pos;                //! histo with difference of tgLambda with respect to value at production vs 1/pt vs phi

  TH3F*** fhdpt_pos;                     //! histo with difference of pt between consecutive layers vs 1/pt vs phi
  TH3F*** fhdpt_prod_pos;                //! histo with difference of pt with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiopt_pos;                     //! histo with ratio of pt between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiopt_prod_pos;                //! histo with ratio of pt with respect to value at production vs 1/pt vs phi

  TH3F*** fhdp_pos;                     //! histo with difference of p between consecutive layers vs 1/pt vs phi
  TH3F*** fhdp_prod_pos;                //! histo with difference of p of with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiop_pos;                     //! histo with ratio of p between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiop_prod_pos;                //! histo with ratio of p of with respect to value at production vs 1/pt vs phi

  TH3F*** fhdpz_pos;                     //! histo with difference of pz between consecutive layers vs 1/pt vs phi
  TH3F*** fhdpz_prod_pos;                //! histo with difference of pz with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiopz_pos;                     //! histo with ratio of pz between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiopz_prod_pos;                //! histo with ratio of pzwith respect to value at production vs 1/pt vs phi

  // negative tracks
  TH3F*** fhdtgl_neg;                     //! histo with difference of tgLambda between consecutive layers vs 1/pt vs phi
  TH3F*** fhdtgl_prod_neg;                //! histo with difference of tgLambda with respect to value at production vs 1/pt vs phi

  TH3F*** fhdpt_neg;                     //! histo with difference of pt between consecutive layers vs 1/pt vs phi
  TH3F*** fhdpt_prod_neg;                //! histo with difference of pt with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiopt_neg;                     //! histo with ratio of pt between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiopt_prod_neg;                //! histo with ratio of pt with respect to value at production vs 1/pt vs phi

  TH3F*** fhdp_neg;                     //! histo with difference of p between consecutive layers vs 1/pt vs phi
  TH3F*** fhdp_prod_neg;                //! histo with difference of p of with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiop_neg;                     //! histo with ratio of p between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiop_prod_neg;                //! histo with ratio of p of with respect to value at production vs 1/pt vs phi

  TH3F*** fhdpz_neg;                     //! histo with difference of pz between consecutive layers vs 1/pt vs phi
  TH3F*** fhdpz_prod_neg;                //! histo with difference of pz with respect to value at production vs 1/pt vs phi
  TH3F*** fhratiopz_neg;                     //! histo with ratio of pz between consecutive layers vs 1/pt vs phi
  TH3F*** fhratiopz_prod_neg;                //! histo with ratio of pzwith respect to value at production vs 1/pt vs phi

  TH1I*  fhNTrackRefsInITS;                    //! histo with number of trackRefs in ITS
  
  ClassDef(AliAnalysisTaskTrackRefsChecksMultipleSpecies,1);
};


#endif
