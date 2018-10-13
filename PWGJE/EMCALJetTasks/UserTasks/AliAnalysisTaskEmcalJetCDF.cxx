/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <cmath>

#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TArrayD.h>
#include <TString.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisDataContainer.h"

#include "AliEmcalList.h"

#include "AliAnalysisTaskEmcalJetCDF.h"

/// \cond CLASSIMP
ClassImp ( AliAnalysisTaskEmcalJetCDF );
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF() :
    AliAnalysisTaskEmcalJet (),
    fHistManager()
{}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF ( const char *name ) :
    AliAnalysisTaskEmcalJet ( name, kTRUE ),
    fHistManager(name)
  {
  // Standard constructor.
  SetMakeGeneralHistograms ( kTRUE );
  }

/// Destructor
AliAnalysisTaskEmcalJetCDF::~AliAnalysisTaskEmcalJetCDF() {}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetCDF::Run()
  {
  return kTRUE;
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::FillHistograms()
  {
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);

  namespace CDF = NS_AliAnalysisTaskEmcalJetCDF;
  TString histname = "", groupname = "", fullgroupname = "";

  AliJetContainer* jetCont = NULL;
  TIter next(&fJetCollArray);
  while ( (jetCont = static_cast<AliJetContainer*>(next())) ) {
    //##### EARLY VALIDITY CHECKS - BAIL OUT FAST

    // get particles connected to jets
    AliParticleContainer* fTracksCont = jetCont->GetParticleContainer();
    if (!fTracksCont) { std::cout << "*********   JET CONTAINER WITHOUT TRACKS CONTAINER   *********" << std::endl; continue; }
    TClonesArray* fTracksContArray = fTracksCont->GetArray();

    // Number of Jets found in event - accepted cuts applied by JetContainer
    Int_t fNJets_accepted = jetCont->GetNJets();

    // Multiplicity in event - accepted tracks in tracks container
    Int_t fNaccPart = fTracksCont->GetNAcceptedParticles();

    // protection
    if ( ( fNJets_accepted < 1 ) || ( fNaccPart < 1 ) ) {
      if ( fDebug > 1 ) { std::cout << "accepted (fNJets || fNPart) == 0" << std::endl; }
      continue;
      }

    if ( fDebug > 1 ) { std::cout << "fNJets = " << fNJets_accepted << " ; fNPart = " << fNaccPart << std::endl; }

    // get jet1 - if there is no leading jet there is no point to continue
    AliEmcalJet* jet1 = jetCont->GetLeadingJet(); // internaly checked for AcceptedJet
    if (!jet1) {
      if ( fDebug > 1 ) { Printf ( "Jet1 not found (did not survive cuts?)\n" ); }
      continue;
      }

//######################################################################################################
    groupname = jetCont->GetName();

    Double_t jet_pt_min = jetCont->GetMinPt();
    Double_t jet_pt_max = jetCont->GetMaxPt();

    TString jetstrmin = TString::Itoa((Int_t)jet_pt_min,10);
    TString jetstrmax = TString::Itoa((Int_t)jet_pt_max,10);

    // add to groupname the min,max pt cuts of jets in the container
    groupname = groupname + "_" + "ptbin" + "_" + jetstrmin + "_" + jetstrmax;

//######################################################################################################
//   Get histo pointers from Hist Manager
    histname = TString::Format("%s/histo1_%d", groupname.Data(), fCentBin);
    TH1D* fH1 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo2_%d", groupname.Data(), fCentBin);
    TH1D* fH2 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo3_%d", groupname.Data(), fCentBin);
    TH1D* fH3 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo4_%d", groupname.Data(), fCentBin);
    TH1D* fH4 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo4c_%d", groupname.Data(), fCentBin);
    TH1D* fH4c = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo5_%d", groupname.Data(), fCentBin);
    TH1D* fH5 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo6_%d", groupname.Data(), fCentBin);
    TH1D* fH6 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo6c_%d", groupname.Data(), fCentBin);
    TH1D* fH6c = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo7_%d", groupname.Data(), fCentBin);
    TH2D* fH7 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo7all_%d", groupname.Data(), fCentBin);
    TH2D* fH7all = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################
    histname = TString::Format("%s/histo8_%d", groupname.Data(), fCentBin);
    TH1D* fH8 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8_all_%d", groupname.Data(), fCentBin);
    TH1D* fH8_all = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8_p_%d", groupname.Data(), fCentBin);
    TH1D* fH8_p = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8_all_p_%d", groupname.Data(), fCentBin);
    TH1D* fH8_all_p = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8_pt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8_all_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8_all_pt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_all_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_all = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_p_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_p = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_all_p_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_all_p = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_pt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_all_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_all_pt = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo15_%d", groupname.Data(), fCentBin);
    TH2D* fH15 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_n70_%d", groupname.Data(), fCentBin);
    TH2D* fH15_n70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_n75_%d", groupname.Data(), fCentBin);
    TH2D* fH15_n75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_n80_%d", groupname.Data(), fCentBin);
    TH2D* fH15_n80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_n85_%d", groupname.Data(), fCentBin);
    TH2D* fH15_n85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_n90_%d", groupname.Data(), fCentBin);
    TH2D* fH15_n90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_pt70_%d", groupname.Data(), fCentBin);
    TH2D* fH15_pt70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_pt75_%d", groupname.Data(), fCentBin);
    TH2D* fH15_pt75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_pt80_%d", groupname.Data(), fCentBin);
    TH2D* fH15_pt80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_pt85_%d", groupname.Data(), fCentBin);
    TH2D* fH15_pt85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15_pt90_%d", groupname.Data(), fCentBin);
    TH2D* fH15_pt90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_%d", groupname.Data(), fCentBin);
    TH2D* fH15all = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n70_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n75_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n80_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n85_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n90_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt70_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt75_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt80_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt85_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt90_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt90 = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo20_%d", groupname.Data(), fCentBin);
    TH1D* fH20 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_n70_%d", groupname.Data(), fCentBin);
    TH1D* fH20_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_n75_%d", groupname.Data(), fCentBin);
    TH1D* fH20_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_n80_%d", groupname.Data(), fCentBin);
    TH1D* fH20_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_n85_%d", groupname.Data(), fCentBin);
    TH1D* fH20_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_n90_%d", groupname.Data(), fCentBin);
    TH1D* fH20_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_pt70_%d", groupname.Data(), fCentBin);
    TH1D* fH20_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_pt75_%d", groupname.Data(), fCentBin);
    TH1D* fH20_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fH20_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_pt85_%d", groupname.Data(), fCentBin);
    TH1D* fH20_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fH20_pt90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_%d", groupname.Data(), fCentBin);
    TH1D* fH20all = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n70_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n75_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n80_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n85_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n90_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt70_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt75_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt85_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_g_%d", groupname.Data(), fCentBin);
    TH1D* fHg = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n70_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n75_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n80_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n85_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n90_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt70_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt75_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt85_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_ptd_%d", groupname.Data(), fCentBin);
    TH1D* fHptd = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n70_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n75_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n80_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n85_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n90_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt70_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt75_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt85_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_Rjt_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n70_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n75_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n80_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n85_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n90_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt70_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt70 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt75_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt75 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt80_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt85_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt85 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt90_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt90 = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_jt_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n70_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n75_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n80_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n85_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n90_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt70_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt75_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt85_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt90 = (TH1D*)GetHistogram(histname.Data());


//######################################################################################################
    // get clusters connected to jets
    AliClusterContainer* fCaloClustersCont = jetCont->GetClusterContainer();
    // accepted clusters in cluster container
    Int_t fNaccClus = -1;
    if (fCaloClustersCont) { fNaccClus = fCaloClustersCont->GetNAcceptedClusters(); }

    fH5->Fill ( fNJets_accepted ); // Distribution of jets in events;

    UShort_t counter_part = 0; Double_t counter_pt = 0.; // counter for npart and pt recording
    UShort_t jet_n90 = -99 ; Double_t jet_pt90 = -99.99 ;
    UShort_t jet_n85 = -99 ; Double_t jet_pt85 = -99.99 ;
    UShort_t jet_n80 = -99 ; Double_t jet_pt80 = -99.99 ;
    UShort_t jet_n75 = -99 ; Double_t jet_pt75 = -99.99 ;
    UShort_t jet_n70 = -99 ; Double_t jet_pt70 = -99.99 ;

    // variables used to compute g and ptD
    Double_t g_tot = 0.;  Double_t sum_part_pt_tot  = 0.; Double_t sum_part_pt2_tot  = 0.;
    Double_t g_n90 = 0.;  Double_t sum_part_pt_n90  = 0.; Double_t sum_part_pt2_n90  = 0.;
    Double_t g_n85 = 0.;  Double_t sum_part_pt_n85  = 0.; Double_t sum_part_pt2_n85  = 0.;
    Double_t g_n80 = 0.;  Double_t sum_part_pt_n80  = 0.; Double_t sum_part_pt2_n80  = 0.;
    Double_t g_n75 = 0.;  Double_t sum_part_pt_n75  = 0.; Double_t sum_part_pt2_n75  = 0.;
    Double_t g_n70 = 0.;  Double_t sum_part_pt_n70  = 0.; Double_t sum_part_pt2_n70  = 0.;
    Double_t g_pt90 = 0.; Double_t sum_part_pt_pt90 = 0.; Double_t sum_part_pt2_pt90 = 0.;
    Double_t g_pt85 = 0.; Double_t sum_part_pt_pt85 = 0.; Double_t sum_part_pt2_pt85 = 0.;
    Double_t g_pt80 = 0.; Double_t sum_part_pt_pt80 = 0.; Double_t sum_part_pt2_pt80 = 0.;
    Double_t g_pt75 = 0.; Double_t sum_part_pt_pt75 = 0.; Double_t sum_part_pt2_pt75 = 0.;
    Double_t g_pt70 = 0.; Double_t sum_part_pt_pt70 = 0.; Double_t sum_part_pt2_pt70 = 0.;

    // **************************************************************
    //                          LEADING JETS
    // **************************************************************
    if ( jet1 ) {
      if ( fDebug > 1 ) { std::cout << "+++++++++++++++++>>>>>>>>> Leading jet found" << std::endl; jet1->Print(); }

      // vector of sorted indexes of particles in leading jet - jet1 : sorting by p_T jet constituents
      std::vector< int > jet1_sorted_idxvec = jet1->GetPtSortedTrackConstituentIndexes ( fTracksContArray );

      Double_t jet1_pt = jet1->Pt();
      UShort_t jet1_npart  = jet1->GetNumberOfTracks();
      UShort_t jet1_nconst = jet1->GetNumberOfConstituents();

      UShort_t jet1_n90  = ( UShort_t ) ( 0.90 * jet1_npart );
      Double_t jet1_pt90 = 0.90 * jet1_pt;

      UShort_t jet1_n85  = ( UShort_t ) ( 0.85 * jet1_npart );
      Double_t jet1_pt85 = 0.85 * jet1_pt;

      UShort_t jet1_n80  = ( UShort_t ) ( 0.80 * jet1_npart );
      Double_t jet1_pt80 = 0.80 * jet1_pt;

      UShort_t jet1_n75  = ( UShort_t ) ( 0.75 * jet1_npart );
      Double_t jet1_pt75 = 0.75 * jet1_pt;

      UShort_t jet1_n70  = ( UShort_t ) ( 0.70 * jet1_npart );
      Double_t jet1_pt70 = 0.70 * jet1_pt;

      fH6->Fill  ( jet1_npart );           // Multiplicity of jet1 - charged tracks
      fH6c->Fill ( jet1_nconst );          // Multiplicity of jet1 - all constituents
      fH7->Fill  ( jet1_pt, jet1_nconst ); // N(jet) vs P_{T}(jet1)

      Int_t track_idx = -999;            // index variable for tracks
      counter_part = 0; counter_pt = 0.; // reset counters
      //___________________________________________________________________________
      // parsing tracks of jet1 (leading jet) in decreasing order of Pt
      for (std::size_t i = 0; i < jet1_npart; i++ ) {
          track_idx = jet1_sorted_idxvec.at (i);
          //track = dynamic_cast<AliVParticle*>(fTracksContArray->At( track_idx ));
          AliVParticle* track = jet1->TrackAt ( track_idx, fTracksContArray );
          if (!track) { std::cout << "Parsing tracks of jet1 :: track not defined but it should!!!" << std::endl; continue; }

          Double_t dpart = jet1->DeltaR ( track );
          Double_t track_pt = track->Pt();

          fH8->Fill   ( jet1->GetZ  ( track ) );  // Momentum distribution for leading jet (FF)
          fH8xi->Fill ( jet1->GetXi ( track ) );  // Momentum distribution for leading jet (FF) xi

          Double_t z_p = CDF::Z_ptot(jet1, track);
          fH8_p->Fill   ( z_p );  // Momentum distribution for jets (FF)
          fH8xi_p->Fill ( CDF::Xi (z_p)  );  // Momentum distribution for jets (FF) xi

          Double_t z_pt = CDF::Z_pt(jet1, track);
          fH8_pt->Fill   ( z_pt );  // Momentum distribution for jets (FF)
          fH8xi_pt->Fill ( CDF::Xi (z_pt) );  // Momentum distribution for jets (FF) xi

          fH15->Fill ( dpart, track_pt, track_pt );      // <p_{T}> track vs the Distance R from Jet1
          fH20->Fill ( dpart );                          // Distribution of R in leading jet

          // fill histograms for 70% of particles with highest pt
          if ( counter_part <= jet1_n70 ) {
              fH15_n70->Fill ( dpart, track_pt, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n70->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for 75% of particles with highest pt
          if ( counter_part <= jet1_n75 ) {
              fH15_n75->Fill ( dpart, track_pt, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n75->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for 80% of particles with highest pt
          if ( counter_part <= jet1_n80 ) {
              fH15_n80->Fill ( dpart, track_pt, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n80->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for 85% of particles with highest pt
          if ( counter_part <= jet1_n85 ) {
              fH15_n85->Fill ( dpart, track_pt, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n85->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for 90% of particles with highest pt
          if ( counter_part <= jet1_n90 ) {
              fH15_n90->Fill ( dpart, track_pt, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n90->Fill ( dpart );                         // Distribution of R in leading jet
              }

          // fill histograms for particles that have first 70% of pt
          if ( counter_pt <= jet1_pt70 ) {
              fH15_pt70->Fill ( dpart, track_pt, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt70->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for particles that have first 75% of pt
          if ( counter_pt <= jet1_pt75 ) {
              fH15_pt75->Fill ( dpart, track_pt, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt75->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt80 ) {
              fH15_pt80->Fill ( dpart, track_pt, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt80->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt85 ) {
              fH15_pt85->Fill ( dpart, track_pt, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt85->Fill ( dpart );                         // Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt90 ) {
              fH15_pt90->Fill ( dpart, track_pt, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt90->Fill ( dpart );                         // Distribution of R in leading jet
              }

          ++counter_part; counter_pt += track_pt;
          } // end of loop over jet1 tracks
      } // end of jet1 (leading jet) processing

    // **************************************************************
    //                        ALL JETS
    // **************************************************************
    Double_t jet_pt = 0. ; UShort_t jet_npart = 0; UShort_t jet_nconst = 0;

    // loop over all jets
    for( auto jet : jetCont->accepted()) {
      if (!jet) { continue; }

      // vector of sorted indexes of particles in jet
      std::vector< int > jet_sorted_idxvec ;

      Int_t track_idx = -999; // index variable for tracks
      jet_pt = 0. ; jet_npart = 0; jet_nconst = 0;
      counter_part = 0; counter_pt = 0.; // reset counters

      // jet : Sorting by p_T jet constituents
      jet_sorted_idxvec.clear();
      jet_sorted_idxvec = jet->GetPtSortedTrackConstituentIndexes ( fTracksContArray );

      jet_pt = jet->Pt();
      jet_npart = jet->GetNumberOfTracks();
      jet_nconst = jet->GetNumberOfConstituents();

      // variables for g and pdt computations
      g_tot = 0.; sum_part_pt_tot = 0.; sum_part_pt2_tot = 0.;
      g_n90 = 0.; sum_part_pt_n90 = 0.; sum_part_pt2_n90 = 0.;
      g_n85 = 0.; sum_part_pt_n85 = 0.; sum_part_pt2_n85 = 0.;
      g_n80 = 0.; sum_part_pt_n80 = 0.; sum_part_pt2_n80 = 0.;
      g_n75 = 0.; sum_part_pt_n75 = 0.; sum_part_pt2_n75 = 0.;
      g_n70 = 0.; sum_part_pt_n70 = 0.; sum_part_pt2_n70 = 0.;
      g_pt90 = 0.; sum_part_pt_pt90 = 0.; sum_part_pt2_pt90 = 0.;
      g_pt85 = 0.; sum_part_pt_pt85 = 0.; sum_part_pt2_pt85 = 0.;
      g_pt80 = 0.; sum_part_pt_pt80 = 0.; sum_part_pt2_pt80 = 0.;
      g_pt75 = 0.; sum_part_pt_pt75 = 0.; sum_part_pt2_pt75 = 0.;
      g_pt70 = 0.; sum_part_pt_pt70 = 0.; sum_part_pt2_pt70 = 0.;

      // sentinels for the pt tail cuts
      jet_n90  = ( UShort_t ) ( 0.9 * jet_npart );
      jet_pt90 = 0.9 * jet_pt;

      jet_n85  = ( UShort_t ) ( 0.85 * jet_npart );
      jet_pt85 = 0.85 * jet_pt;

      jet_n80  = ( UShort_t ) ( 0.8 * jet_npart );
      jet_pt80 = 0.8 * jet_pt;

      jet_n75  = ( UShort_t ) ( 0.75 * jet_npart );
      jet_pt75 = 0.75 * jet_pt;

      jet_n70  = ( UShort_t ) ( 0.7 * jet_npart );
      jet_pt70 = 0.7 * jet_pt;

      fH1->Fill ( jet_pt );            // Pt distribution of jets
      fH2->Fill ( jet->Eta() );        // Eta distribution of jets
      fH3->Fill ( jet->Phi() );        // Phi distribution of jets
      fH4->Fill ( jet_npart );         // Multiplicity of jets
      fH4c->Fill ( jet_nconst );       // Multiplicity of jets - all constituents
      fH7all->Fill ( jet_pt, jet_nconst ); // N(jet) vs P_{T} - all jets

      // parsing all jet tracks
      for (std::size_t i = 0; i < jet_npart; i++ ) {
        track_idx = jet_sorted_idxvec.at (i);
        //track = dynamic_cast<AliVParticle*>(fTracksContArray->At( track_idx ));
        AliVParticle* track = jet->TrackAt ( track_idx, fTracksContArray );
        if (!track) { std::cout << "Parsing tracks of jets :: track not defined but it should!!!" << std::endl; continue; }

        Double_t dpart = jet->DeltaR ( track );
        Double_t track_pt = track->Pt();
        Double_t jt = CDF::Perp (*track, *jet);

        fH8_all->Fill   ( jet->GetZ  ( track ) );  // Momentum distribution for jets (FF)
        fH8xi_all->Fill ( jet->GetXi ( track ) );  // Momentum distribution for jets (FF) xi

        Double_t z_p = CDF::Z_ptot(jet,track);
        fH8_all_p->Fill   ( z_p );        // Momentum distribution for jets (FF)
        fH8xi_all_p->Fill ( CDF::Xi (z_p)  );  // Momentum distribution for jets (FF) xi

        Double_t z_pt = CDF::Z_pt(jet,track);
        fH8_all_pt->Fill   ( z_pt );       // Momentum distribution for jets (FF)
        fH8xi_all_pt->Fill ( CDF::Xi (z_pt) );  // Momentum distribution for jets (FF) xi

        fH15all->Fill ( dpart, track_pt, track_pt );    // p_{T} track vs the Distance R from jet
        fH20all->Fill ( dpart );                        // Distribution of R in leading jet

        fH_Rjt->Fill ( dpart, jt, jt );    // jt track vs dR weighted with jt
        fH_jt->Fill  ( jt );                // jt track vs dR

        // computing components for g and ptD in the jet tracks loop
        g_tot += (track_pt * dpart)/jet_pt;
        sum_part_pt_tot += TMath::Abs(track_pt);
        sum_part_pt2_tot += TMath::Power( track_pt, 2 );

  //#############################################################################################
        if ( counter_part <= jet_n90 ) {        /// N90
            fH15all_n90->Fill ( dpart, track_pt );         // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n90->Fill ( dpart );
            fH_Rjt_n90->Fill  ( dpart, jt, jt );
            fH_jt_n90->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_n90 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n90 += track_pt;
            sum_part_pt2_n90 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt90 ) {        /// PT90
            fH15all_pt90->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt90->Fill ( dpart );
            fH_Rjt_pt90->Fill  ( dpart, jt, jt );
            fH_jt_pt90->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_pt90 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt90 += track_pt;
            sum_part_pt2_pt90 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n85 ) {       /// N85
            fH15all_n85->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n85->Fill ( dpart );
            fH_Rjt_n85->Fill  ( dpart, jt, jt );
            fH_jt_n85->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_n85 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n85 += track_pt;
            sum_part_pt2_n85 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt85 ) {       /// PT85
            fH15all_pt85->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt85->Fill ( dpart );
            fH_Rjt_pt85->Fill  ( dpart, jt, jt );
            fH_jt_pt85->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_pt85 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt85 += track_pt;
            sum_part_pt2_pt85 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n80 ) {       /// N80
            fH15all_n80->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n80->Fill ( dpart );
            fH_Rjt_n80->Fill  ( dpart, jt, jt );
            fH_jt_n80->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_n80 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n80 += track_pt;
            sum_part_pt2_n80 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt80 ) {        /// PT80
            fH15all_pt80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt80->Fill ( dpart );
            fH_Rjt_pt80->Fill  ( dpart, jt, jt );
            fH_jt_pt80->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_pt80 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt80 += track_pt;
            sum_part_pt2_pt80 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n75 ) {       /// N75
            fH15all_n75->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n75->Fill ( dpart );
            fH_Rjt_n75->Fill  ( dpart, jt, jt );
            fH_jt_n75->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_n75 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n75 += track_pt;
            sum_part_pt2_n75 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt75 ) {        /// PT75
            fH15all_pt75->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt75->Fill ( dpart );
            fH_Rjt_pt75->Fill  ( dpart, jt, jt );
            fH_jt_pt75->Fill  ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_pt75 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt75 += track_pt;
            sum_part_pt2_pt75 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n70 ) {      /// N70
            fH15all_n70->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n70->Fill ( dpart );
            fH_Rjt_n70->Fill  ( dpart, jt, jt );
            fH_jt_n70->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_n70 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n70 += track_pt;
            sum_part_pt2_n70 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt70 ) {       /// PT70
            fH15all_pt70->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt70->Fill ( dpart );
            fH_Rjt_pt70->Fill  ( dpart, jt, jt );
            fH_jt_pt70->Fill   ( jt );                // jt track vs dR

            // computing components for g and ptD in the jet tracks loop
            g_pt70 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt70 += track_pt;
            sum_part_pt2_pt70 += TMath::Power( track_pt, 2 );
            }
        ++counter_part; counter_pt += track_pt;
        } // end of loop over jet tracks

      fHg->Fill     ( g_tot );
      fHg_n70->Fill ( g_n70 );    fHg_pt70->Fill ( g_pt70 );
      fHg_n75->Fill ( g_n75 );    fHg_pt75->Fill ( g_pt75 );
      fHg_n80->Fill ( g_n80 );    fHg_pt80->Fill ( g_pt80 );
      fHg_n85->Fill ( g_n85 );    fHg_pt85->Fill ( g_pt85 );
      fHg_n90->Fill ( g_n90 );    fHg_pt90->Fill ( g_pt90 );

      if ( sum_part_pt_tot > 1e-8 )
        { fHptd->Fill( TMath::Sqrt(sum_part_pt2_tot)/sum_part_pt_tot ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_tot aprox 0!!!!" << endl; }

      if ( sum_part_pt_n70 > 1e-8 )
        { fHptd_n70->Fill( TMath::Sqrt(sum_part_pt2_n70)/sum_part_pt_n70 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n70 aprox 0!!!!" << endl; }

      if ( sum_part_pt_n75 > 1e-8 )
        { fHptd_n75->Fill( TMath::Sqrt(sum_part_pt2_n75)/sum_part_pt_n75 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n75 aprox 0!!!!" << endl; }

      if ( sum_part_pt_n80 > 1e-8 )
        { fHptd_n80->Fill( TMath::Sqrt(sum_part_pt2_n80)/sum_part_pt_n80 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n80 aprox 0!!!!" << endl; }

      if ( sum_part_pt_n85 > 1e-8 )
        { fHptd_n85->Fill( TMath::Sqrt(sum_part_pt2_n85)/sum_part_pt_n85 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n85 aprox 0!!!!" << endl; }

      if ( sum_part_pt_n90 > 1e-8 )
        { fHptd_n90->Fill( TMath::Sqrt(sum_part_pt2_n90)/sum_part_pt_n90 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n90 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt70 > 1e-8 )
        { fHptd_pt70->Fill( TMath::Sqrt(sum_part_pt2_pt70)/sum_part_pt_pt70 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt70 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt75 > 1e-8 )
        { fHptd_pt75->Fill( TMath::Sqrt(sum_part_pt2_pt75)/sum_part_pt_pt75 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt75 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt80 > 1e-8 )
        { fHptd_pt80->Fill( TMath::Sqrt(sum_part_pt2_pt80)/sum_part_pt_pt80 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt80 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt85 > 1e-8 )
        { fHptd_pt85->Fill( TMath::Sqrt(sum_part_pt2_pt85)/sum_part_pt_pt85 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt85 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt90 > 1e-8 )
        { fHptd_pt90->Fill( TMath::Sqrt(sum_part_pt2_pt90)/sum_part_pt_pt90 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt90 aprox 0!!!!" << endl; }

      } // end of loopt over all jets
    } // end of loop over jet container collection

  // post data at every processing
  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  return kTRUE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::UserCreateOutputObjects()
  {
  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString histname = "", histtitle = "", groupname = "", fullgroupname = "";
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next())))
    {
    groupname = jetCont->GetName();

    Double_t jet_pt_min = jetCont->GetMinPt();
    Double_t jet_pt_max = jetCont->GetMaxPt();

    TString jetstrmin = TString::Itoa((Int_t)jet_pt_min,10);
    TString jetstrmax = TString::Itoa((Int_t)jet_pt_max,10);

    // add to groupname the min,max pt cuts of jets in the container
    groupname = groupname + "_" + "ptbin" + "_" + jetstrmin + "_" + jetstrmax;

    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++)
      {
      //=====================================================================================
      Int_t h1_nbin = 300; Double_t h1_binwidth = 1; Double_t h1_low = 0;
      Double_t h1_high = h1_low + h1_binwidth * h1_nbin; // 1GeV/bin
      histname = TString::Format("%s/histo1_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c}) (accepted);Jets", histname.Data()); // Pt distro of jets
      fHistManager.CreateTH1(histname, histtitle, h1_nbin, h1_low, h1_high);

      //=====================================================================================
      Int_t h2_nbin = 200; Double_t h2_binwidth = 0.01; Double_t h2_low = -1;
      Double_t h2_high = h2_low + h2_binwidth * h2_nbin;
      histname = TString::Format("%s/histo2_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};Jets", histname.Data());  // Eta distro of jets
      fHistManager.CreateTH1(histname, histtitle, h2_nbin, h2_low, h2_high); // 1 unit of rapidity / 100 bin

      //=====================================================================================
      Int_t h3_nbin = 126; Double_t h3_binwidth = 0.05; Double_t h3_low = 0.;
      Double_t h3_high = h3_low + h3_binwidth * h3_nbin;
      histname = TString::Format("%s/histo3_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};Jets", histname.Data()); // Phi distro of jets
      fHistManager.CreateTH1(histname, histtitle, h3_nbin, h3_low, h3_high);

      //=====================================================================================
      Int_t h4_nbin = 100; Double_t h4_binwidth = 1; Double_t h4_low = 0;
      Double_t h4_high = h4_low + h4_binwidth * h4_nbin; // 1 unit of multiplicity /bin
      histname = TString::Format("%s/histo4_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{tracks}(jet);Jets", histname.Data());  // Multiplicity of all jets; chg tracks
      fHistManager.CreateTH1(histname, histtitle, h4_nbin, h4_low, h4_high);

      histname = TString::Format("%s/histo4c_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{tracks}(jet);Jets", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h4_nbin, h4_low, h4_high);    // Multiplicity of all jets; all tracks

      histname = TString::Format("%s/histo6_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{tracks}(jet1);Jets", histname.Data()); // Multiplicity of jet1; chg tracks
      fHistManager.CreateTH1(histname, histtitle, h4_nbin, h4_low, h4_high);

      histname = TString::Format("%s/histo6c_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{tracks}(jet1);Jets", histname.Data()); // Multiplicity of jet1; all tracks
      fHistManager.CreateTH1(histname, histtitle, h4_nbin, h4_low, h4_high);
      //#####################################

      //=====================================================================================
      Int_t h5_nbin = 100; Double_t h5_binwidth = 1; Double_t h5_low = 0;
      Double_t h5_high = h5_low + h5_binwidth * h5_nbin;
      histname = TString::Format("%s/histo5_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{jets};Events", histname.Data()); // Distribution of jets in events
      fHistManager.CreateTH1(histname, histtitle, h5_nbin, h5_low, h5_high);

      //=====================================================================================
      Int_t h7_xnbin = 300; Double_t h7_xbinwidth = 1; Double_t h7_xlow = 0;
      Double_t h7_xhigh = h7_xlow + h7_xbinwidth * h7_xnbin;
      Int_t h7_ynbin = 100; Double_t h7_ybinwidth = 1; Double_t h7_ylow = 0;
      Double_t h7_yhigh = h7_ylow + h7_ybinwidth * h7_ynbin;

      histname = TString::Format("%s/histo7_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet1} (GeV/c);N_{tracks}(jet1)", histname.Data()); // N vs pt jet1
      fHistManager.CreateTH2(histname, histtitle, h7_xnbin, h7_xlow, h7_xhigh, h7_ynbin, h7_ylow, h7_yhigh);

      histname = TString::Format("%s/histo7all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/c);N_{tracks}(jets)", histname.Data()); // N vs pt all jets
      fHistManager.CreateTH2(histname, histtitle, h7_xnbin, h7_xlow, h7_xhigh, h7_ynbin, h7_ylow, h7_yhigh);
      //#####################################

      //=====================================================================================
      Int_t h8_nbin = 101; Double_t h8_binwidth = 0.01; Double_t h8_low = 0;
      Double_t h8_high = h8_low + h8_binwidth * h8_nbin;

      // Standard implementation of Z
      histname = TString::Format("%s/histo8_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;z;F(Z) = 1/N_{jet1} dN/dz", histname.Data()); // scalar z ; jet1
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      histname = TString::Format("%s/histo8_all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;z;F(Z) = 1/N_{jets} dN/dz", histname.Data()); // scalar z ; all jets
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      //########################################################
      // P_tot implementation of Z
      histname = TString::Format("%s/histo8_p_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1 P_tot;z = p_{track}/p_{jet1};F(Z) = 1/N_{jet1} dN/dz", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      histname = TString::Format("%s/histo8_all_p_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets P_tot;z = p_{track}/p_{jet};F(Z) = 1/N_{jets} dN/dz", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      //########################################################
      // Pt implementation of Z
      histname = TString::Format("%s/histo8_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1 Pt;z = p_{T,track}/p_{T,jet1};F(Z) = 1/N_{jet1} dN/dz", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      histname = TString::Format("%s/histo8_all_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets Pt;z = p_{T,track}/p_{T,jet1};F(Z) = 1/N_{jets} dN/dz", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);
      //########################################################

      //=====================================================================================
      Int_t h8xi_nbin = 140; Double_t h8xi_binwidth = 0.05; Double_t h8xi_low = 0;
      Double_t h8xi_high = h8xi_low + h8xi_binwidth * h8xi_nbin;
      histname = TString::Format("%s/histo8xi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;#xi = ln(1/z);D(#xi) = 1/N_{jet1} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      histname = TString::Format("%s/histo8xi_all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;#xi = ln(1/z);D(#xi) = 1/N_{jets} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      //########################################################
      histname = TString::Format("%s/histo8xi_p_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1 P_tot;#xi = ln(1/z);D(#xi) = 1/N_{jet1} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      histname = TString::Format("%s/histo8xi_all_p_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets P_tot;#xi = ln(1/z);D(#xi) = 1/N_{jets} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      //########################################################
      histname = TString::Format("%s/histo8xi_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1 Pt;#xi = ln(1/z);D(#xi) = 1/N_{jet1} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      histname = TString::Format("%s/histo8xi_all_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets Pt;#xi = ln(1/z);D(#xi) = 1/N_{jets} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);
      //########################################################

      //=====================================================================================
      Int_t h15_xnbin = 60; Double_t h15_xbinwidth = 0.01; Double_t h15_xlow = 0.;
      Double_t h15_xhigh = h15_xlow + h15_xbinwidth * h15_xnbin;
      Int_t h15_ynbin = 400; Double_t h15_ybinwidth = 1.; Double_t h15_ylow = 0.;
      Double_t h15_yhigh = h15_ylow + h15_ybinwidth * h15_ynbin;

      //########################################################
      histname = TString::Format("%s/histo15_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo15all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);
      //########################################################

      //=====================================================================================
      Int_t h20_nbin = 60; Double_t h20_binwidth = 0.01; Double_t h20_low = 0.;
      Double_t h20_high = h20_low + h20_binwidth * h20_nbin;

      //########################################################
      histname = TString::Format("%s/histo20_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      //########################################################
      histname = TString::Format("%s/histo20_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo20_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo20all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      //########################################################
      histname = TString::Format("%s/histo20all_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo20all_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);
      //########################################################

      //=====================================================================================
      // Distribution of girth (radial girth) g = sum_jet_parts ( r_i * ( pt_i/pt_jet ) )
      Int_t hg_nbin = 50; Double_t hg_binwidth = 0.01; Double_t hg_low = 0.;
      Double_t hg_high = hg_low + hg_binwidth * hg_nbin;

      //########################################################
      histname = TString::Format("%s/histo_g_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_g_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_g_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //=====================================================================================
      // Distribution of dispersion d pt_D = sqrt ( sum (pt_i^2) )/sum (pt_i)
      Int_t hptd_nbin = 40; Double_t hptd_binwidth = 0.05; Double_t hptd_low = 0.;
      Double_t hptd_high = hptd_low + hptd_binwidth * hptd_nbin;

      //########################################################
      histname = TString::Format("%s/histo_ptd_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_ptd_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_ptd_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);
      //########################################################

      //=====================================================================================
      Int_t h_Rjt_xnbin = 60; Double_t h_Rjt_xbinwidth = 0.01; Double_t h_Rjt_xlow = 0.;
      Double_t h_Rjt_xhigh = h_Rjt_xlow + h_Rjt_xbinwidth * h_Rjt_xnbin;
      Int_t h_Rjt_ynbin = 500; Double_t h_Rjt_ybinwidth = 0.01; Double_t h_Rjt_ylow = 0.;
      Double_t h_Rjt_yhigh = h_Rjt_ylow + h_Rjt_ybinwidth * h_Rjt_ynbin;

      //########################################################
      histname = TString::Format("%s/histo_Rjt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      //########################################################
      histname = TString::Format("%s/histo_Rjt_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      //########################################################
      histname = TString::Format("%s/histo_Rjt_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);
      //########################################################

      //=====================================================================================
      Int_t h_jt_xnbin = 500; Double_t h_jt_xbinwidth = 0.01; Double_t h_jt_xlow = 0.;
      Double_t h_jt_xhigh = h_jt_xlow + h_jt_xbinwidth * h_jt_xnbin;

      //########################################################
      histname = TString::Format("%s/histo_jt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      //########################################################
      histname = TString::Format("%s/histo_jt_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      //########################################################
      histname = TString::Format("%s/histo_jt_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);
      //########################################################

      }
      //end of loop over fNcentBins
    }
    // end of loop over jet containers

  // =========== Switch on Sumw2 for all histos ===========
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);

  // add all fHistManager content to fOutput
  TIter nexthist(fHistManager.GetListOfHistograms());
  TObject* obj = NULL;
  while ((obj = nexthist())) { fOutput->Add(obj); }

  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  }

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetCDF::ExecOnce()
  {
  AliAnalysisTaskEmcalJet::ExecOnce();
  }

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetCDF::Terminate ( Option_t * )
  {
  }

//________________________________________________________________________
TObject* AliAnalysisTaskEmcalJetCDF::GetHistogram ( const char* histName )
{
  return fHistManager.FindObject(histName);
}

//########################################################################
//   Namespace AliAnalysisTaskEmcalJetCDF
//########################################################################

//__________________________________________________________________________________________________
std::vector<Int_t> NS_AliAnalysisTaskEmcalJetCDF::SortTracksPt ( AliVEvent* event )
  {
  // Sorting by p_T (decreasing) event tracks
  Int_t entries = event->GetNumberOfTracks();

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;
  pair_list.reserve ( entries );

  for ( Int_t i_entry = 0; i_entry < entries; i_entry++ )
      {
      AliVParticle* track = event->GetTrack ( i_entry );
      if ( !track ) { std::cout << Form ("Unable to find track %d in collection %s", i_entry, event->GetName()) << std::endl ; continue; }

      pair_list.push_back ( std::make_pair ( track->Pt(), i_entry ) );
      }

  std::stable_sort ( pair_list.begin(), pair_list.end(), sort_descend() );

  // return an vector of indexes of constituents (sorted descending by pt)
  std::vector<Int_t> index_sorted_list;
  index_sorted_list.reserve ( entries );

  // populating the return object with indexes of sorted tracks
  for (auto it : pair_list) { index_sorted_list.push_back(it.second); }

  return index_sorted_list;
  }

//__________________________________________________________________________________________________
std::vector<Int_t> NS_AliAnalysisTaskEmcalJetCDF::SortTracksPt ( AliParticleContainer* trackscont )
  {
  // Sorting by p_T (decreasing) event tracks
  Int_t entries = trackscont->GetNEntries();

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;
  pair_list.reserve ( entries );

  UInt_t i_entry = 0;
  AliParticleContainer* partCont = 0;
  for(auto part : partCont->all()) {
    if (!part) {continue;}
    pair_list.push_back ( std::make_pair ( part->Pt(), i_entry++ ) );
    }

  std::stable_sort ( pair_list.begin(), pair_list.end(), sort_descend() );

  // return an vector of indexes of constituents (sorted descending by pt)
  std::vector<Int_t> index_sorted_list;
  index_sorted_list.reserve ( i_entry );

  // populating the return object with indexes of sorted tracks
  for (auto it : pair_list) { index_sorted_list.push_back(it.second); }

  return index_sorted_list;
  }


/// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
/// \param ntracks name of tracks collection
/// \param nclusters name of clusters collection
/// \param ncells name of EMCAL cell collection
/// \param tag tag name of analysis task
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF* NS_AliAnalysisTaskEmcalJetCDF::AddTaskEmcalJetCDF ( const char* ntracks, const char* nclusters, const char* ncells, const char* tag)
  {
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if ( !mgr ) { ::Error ( "AddTaskEmcalJetCDF", "No analysis manager to connect to." );  return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) { ::Error ( "AddTaskEmcalJetCDF", "This task requires an input event handler" ); return NULL; }

  enum EDataType_t { kUnknown, kESD, kAOD }; EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) { dataType = kESD; }
  else
  if (handler->InheritsFrom("AliAODInputHandler")) { dataType = kAOD; }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString suffix   ( tag );
  TString tracks   ( ntracks );
  TString clusters ( nclusters );
  TString cells    ( ncells );

  if ( tracks.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { tracks = "Tracks"; }
    else
    if ( dataType == kAOD ) { tracks = "tracks"; }
    else
      { tracks = ""; }
    }

  if ( clusters.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { clusters = "CaloClusters"; }
    else
    if ( dataType == kAOD ) { clusters = "caloClusters"; }
    else
      { clusters = ""; }
    }

  if ( cells.EqualTo("usedefault") )
    {
    if (dataType == kESD) { cells = "EMCALCells"; }
    else
    if (dataType == kAOD) { cells = "emcalCells"; }
    else
      { cells = ""; }
    }

  TString name("JetCDF");
  if (!tracks.IsNull())   { name += "_" + tracks; }
  if (!clusters.IsNull()) { name += "_" + clusters; }
  if (!cells.IsNull())    { name += "_" + cells; }
  if (!suffix.IsNull())   { name += "_" + suffix; }

  AliAnalysisTaskEmcalJetCDF* cdfTask = new AliAnalysisTaskEmcalJetCDF ( name.Data() );
  cdfTask->SetVzRange(-10,10);
  cdfTask->SetCaloCellsName(cells.Data());

  if ( tracks.EqualTo("mcparticles") ) {
      // AliMCParticleContainer* mcpartCont =
      cdfTask->AddMCParticleContainer ( tracks.Data() );
      }
  else
  if ( tracks.EqualTo("tracks") || tracks.EqualTo("Tracks") ) {
      // AliTrackContainer* trackCont =
      cdfTask->AddTrackContainer( tracks.Data() );
      }
  else
  if ( !tracks.IsNull())
    { cdfTask->AddParticleContainer(tracks.Data()); }

  cdfTask->AddClusterContainer(clusters.Data());

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask ( cdfTask );

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

  TString contname = name + "_histos";
  TString outfile (Form("%s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer ( contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data() );

  mgr->ConnectInput  ( cdfTask, 0,  cinput1 );
  mgr->ConnectOutput ( cdfTask, 1, coutput1 );

  return cdfTask;
  }

  /// Set parameters of a jet container
  /// \param jetCont AliJetContainer*
  /// \param jetptmin : min pt of jets in this container (default = 1.)
  /// \param jetptmax : max pt of jets in this container (default = 500.)
  /// \param jetareacutperc : cut jets under percentage of area given by algo radius (default = 0.)
  /// \param leadhadtype : 0 = charged, 1 = neutral, 2 = both (default = 2)
  /// \param nLeadJets : how many jets are to be considered the leading jet(s) (default = 1)
  /// \param mintrackpt : min track constituent pt to accept the jet (default = 0.15)
  /// \param maxtrackpt : max track constituent pt to accept the jet (default = 1000.)
  void NS_AliAnalysisTaskEmcalJetCDF::jetContSetParams ( AliJetContainer* jetCont, Float_t jetptmin,  Float_t jetptmax, Float_t jetareacutperc, Int_t leadhadtype, Int_t nLeadJets, Float_t mintrackpt, Float_t maxtrackpt)
    {
    if (!jetCont) { return; }
    jetCont->SetJetPtCut ( jetptmin );
    jetCont->SetJetPtCutMax ( jetptmax );
    jetCont->SetPercAreaCut ( jetareacutperc );
    jetCont->SetLeadingHadronType ( leadhadtype ); // 0 = charged, 1 = neutral, 2 = both
    jetCont->SetNLeadingJets(nLeadJets);
    jetCont->SetMinTrackPt(mintrackpt);
    jetCont->SetMaxTrackPt(maxtrackpt);
    }


// kate: indent-mode none; indent-width 2; replace-tabs on;

