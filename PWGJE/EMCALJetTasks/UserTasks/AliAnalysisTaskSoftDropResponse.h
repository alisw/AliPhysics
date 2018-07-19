#ifndef ALIANALYSISTASKSOFTDROPRESPONSE_H
#define ALIANALYSISTASKSOFTDROPRESPONSE_H
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

//-----------------------------------------------------------------------
// Author : Kirill Lapidus, Yale University, kirill.lapidus@cern.ch
//-----------------------------------------------------------------------

#include "AliJetResponseMaker.h"
#include "AliEmcalEmbeddingQA.h"

#include "FJ_includes.h"

class AliAnalysisTaskSoftDropResponse : public AliJetResponseMaker {
 public:
  AliAnalysisTaskSoftDropResponse();
  AliAnalysisTaskSoftDropResponse(const char *name);
  virtual ~AliAnalysisTaskSoftDropResponse();
  
  void                        UserCreateOutputObjects();

  static AliAnalysisTaskSoftDropResponse * AddTaskSoftDropResponse(
      const char *ntracks1           = "Tracks",
      const char *nclusters1         = "CaloClusters",
      const char *njets1             = "Jets",
      const char *nrho1              = "Rho",
      const Double_t jetradius1      = 0.4,
      const char *ntracks2           = "MCParticles",
      const char *nclusters2         = "",
      const char *njets2             = "MCJets",
      const char *nrho2              = "",
      const Double_t    jetradius2         = 0.4,
      const Double_t    jetptcut           = 1,
      const Double_t    jetareacut         = 0.557,
      const Double_t    jetBias            = 5,
      const Int_t       biasType           = 0,   //  0 = charged, 1 = neutral, 2 = both
      const AliAnalysisTaskSoftDropResponse::MatchingType matching = AliAnalysisTaskSoftDropResponse::kGeometrical,
      const Double_t    maxDistance1       = 0.25,
      const Double_t    maxDistance2       = 0.25,
      const char *cutType            = "TPC",
      const Int_t       ptHardBin          = -999,
      const Double_t    minCent            = -999,
      const Double_t    maxCent            = -999,
      const char *taskname           = "AliAnalysisTaskSoftDropResponse",
      const Bool_t      biggerMatrix       = kFALSE,
      AliAnalysisTaskSoftDropResponse* address = 0,
      const Double_t    maxTrackPt         = 100);

  void                        SetZ2gAxis(Int_t b) { fZ2gAxis = b; }

 protected:

  Bool_t                      FillHistograms();
  void                        FillMatchingHistos(AliEmcalJet* jet1, AliEmcalJet* jet2, Double_t d, Double_t CE1, Double_t CE2);
  void                        AllocateTHnSparse();

  void                        CalculateZg(AliEmcalJet* jet, const Float_t zcut, const Float_t beta, Float_t zg);

  Int_t                       fZ2gAxis;                                 ///< add Z2g axis in matching THnSparse (default=0)

 private:
 
  void               Decluster(const fastjet::PseudoJet& jet);
  fastjet::ClusterSequence* Recluster(const AliEmcalJet* jet);
  
  void               FillZgRgVectors(const AliEmcalJet* jet);

  Int_t fNsdsteps; //!<!
  std::vector<Float_t> fZg_values; //!<! n_SD ordering
  std::vector<Float_t> fRg_values; //!<! n_SD ordering
  std::vector<Int_t>   fSDsteps_values; //!<!
 
  AliAnalysisTaskSoftDropResponse(const AliAnalysisTaskSoftDropResponse&); // not implemented
  AliAnalysisTaskSoftDropResponse &operator=(const AliAnalysisTaskSoftDropResponse&);  // not implemented

  ClassDef(AliAnalysisTaskSoftDropResponse, 1)
};
#endif