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

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliEmcalList.h"

#include "AliAnalysisTaskEmcalJetCDF.h"

using std::cout;
using std::endl;

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
  TString histname = "", groupname = "";

  AliJetContainer* jetCont = NULL;
  TIter next(&fJetCollArray);

  while ((jetCont = static_cast<AliJetContainer*>(next())))
    {
    if (!jetCont) { continue; }
    groupname = jetCont->GetName();

    // Number of Jets found in event - accepted cuts applied by JetContainer
    Int_t fNJets_accepted = jetCont->GetNJets();

    // get particles connected to jets
    AliParticleContainer* fTracksCont = jetCont->GetParticleContainer();
    if (!fTracksCont) { std::cout << "*********   JET CONTAINER WITHOUT TRACKS CONTAINER   *********" << endl; continue; }
    TClonesArray* fTracksContArray = fTracksCont->GetArray();

    // Multiplicity in event - accepted tracks in tracks container
    Int_t fNaccPart = fTracksCont->GetNAcceptedParticles();

    // get clusters connected to jets
    AliClusterContainer* fCaloClustersCont = jetCont->GetClusterContainer();
    // accepted clusters in cluster container
    Int_t fNaccClus = -1;
    if (fCaloClustersCont) { fNaccClus = fCaloClustersCont->GetNAcceptedClusters(); }

    // protection
    if ( ( fNJets_accepted < 1 ) || ( fNaccPart < 1 ) )
      {
      if ( fDebug > 1 ) { std::cout << "accepted (fNJets || fNPart) == 0" << std::endl; }
      return kFALSE;
      }

    if ( fDebug > 1 )
      { std::cout << "fNJets = " << fNJets_accepted << " ; fNPart = " << fNaccPart << std::endl; }

    AliVParticle* track = NULL;
    AliEmcalJet*    jet = NULL;

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

    histname = TString::Format("%s/histo_g_%d", groupname.Data(), fCentBin);
    TH1D* fHg = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n70%d", groupname.Data(), fCentBin);
    TH1D* fHg_n70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n75%d", groupname.Data(), fCentBin);
    TH1D* fHg_n75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n80%d", groupname.Data(), fCentBin);
    TH1D* fHg_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n85%d", groupname.Data(), fCentBin);
    TH1D* fHg_n85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n90%d", groupname.Data(), fCentBin);
    TH1D* fHg_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt70%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt70 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt75%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt75 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt80%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt85%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt85 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt90%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt90 = (TH1D*)GetHistogram(histname.Data());

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

    fH5->Fill ( fNJets_accepted ); // Distribution of jets in events;

    UShort_t counter_part = 0; Double_t counter_pt = 0.; // counter for npart and pt recording

    UShort_t jet_n90 = -99 ; Double_t jet_pt90 = -99.99 ;
    UShort_t jet_n85 = -99 ; Double_t jet_pt85 = -99.99 ;
    UShort_t jet_n80 = -99 ; Double_t jet_pt80 = -99.99 ;
    UShort_t jet_n75 = -99 ; Double_t jet_pt75 = -99.99 ;
    UShort_t jet_n70 = -99 ; Double_t jet_pt70 = -99.99 ;

    // variables used to compute g and ptD
    Double_t g_tot = 0.; Double_t sum_part_pt_tot = 0.; Double_t sum_part_pt2_tot = 0.;

    Double_t g_n90 = 0.; Double_t sum_part_pt_n90 = 0.; Double_t sum_part_pt2_n90 = 0.;
    Double_t g_n85 = 0.; Double_t sum_part_pt_n85 = 0.; Double_t sum_part_pt2_n85 = 0.;
    Double_t g_n80 = 0.; Double_t sum_part_pt_n80 = 0.; Double_t sum_part_pt2_n80 = 0.;
    Double_t g_n75 = 0.; Double_t sum_part_pt_n75 = 0.; Double_t sum_part_pt2_n75 = 0.;
    Double_t g_n70 = 0.; Double_t sum_part_pt_n70 = 0.; Double_t sum_part_pt2_n70 = 0.;

    Double_t g_pt90 = 0.; Double_t sum_part_pt_pt90 = 0.; Double_t sum_part_pt2_pt90 = 0.;
    Double_t g_pt85 = 0.; Double_t sum_part_pt_pt85 = 0.; Double_t sum_part_pt2_pt85 = 0.;
    Double_t g_pt80 = 0.; Double_t sum_part_pt_pt80 = 0.; Double_t sum_part_pt2_pt80 = 0.;
    Double_t g_pt75 = 0.; Double_t sum_part_pt_pt75 = 0.; Double_t sum_part_pt2_pt75 = 0.;
    Double_t g_pt70 = 0.; Double_t sum_part_pt_pt70 = 0.; Double_t sum_part_pt2_pt70 = 0.;


    // **************************************************************
    //                          LEADING JETS
    // **************************************************************
    AliEmcalJet* jet1 = jetCont->GetLeadingJet(); // internaly checked for AcceptedJet

    if ( fDebug > 1 )
      {
      if (jet1)
        { std::cout << "+++++++++++++++++>>>>>>>>> Leading jet found" << std::endl; jet1->Print(); }
      else
        { Printf ( "Jet1 not found (did not survive cuts?)\n" ); }
      }

    if (jet1)
      {
      //vector of sorted indexes of particles in leading jet
      std::vector< int > jet1_sorted_idxvec;

      // jet1 : Sorting by p_T jet constituents
      jet1_sorted_idxvec = jet1->SortConstituentsPt ( fTracksContArray );

      Int_t track_idx = -999;            // index variable for tracks
      Double_t jet1_pt = jet1->Pt();
      UInt_t jet1_npart = jet1->GetNumberOfTracks();
      UInt_t jet1_nconst = jet1->GetNumberOfConstituents();

      UShort_t jet1_n90  = ( UShort_t ) ( 0.9 * jet1_npart );
      Double_t jet1_pt90 = 0.9 * jet1_pt;

      UShort_t jet1_n85  = ( UShort_t ) ( 0.85 * jet1_npart );
      Double_t jet1_pt85 = 0.85 * jet1_pt;

      UShort_t jet1_n80  = ( UShort_t ) ( 0.8 * jet1_npart );
      Double_t jet1_pt80 = 0.8 * jet1_pt;

      UShort_t jet1_n75  = ( UShort_t ) ( 0.75 * jet1_npart );
      Double_t jet1_pt75 = 0.75 * jet1_pt;

      UShort_t jet1_n70  = ( UShort_t ) ( 0.7 * jet1_npart );
      Double_t jet1_pt70 = 0.7 * jet1_pt;

      fH6->Fill  ( jet1_npart );        // Multiplicity of jet1 - charged tracks
      fH6c->Fill ( jet1_nconst );       // Multiplicity of jet1 - all constituents
      fH7->Fill ( jet1_pt, jet1_nconst ); // N(jet) vs P_{T}(jet1)

      counter_part = 0; counter_pt = 0.; // reset counters

      //___________________________________________________________________________
      // parsing tracks of jet1 (leading jet) in decreasing order of Pt

      for ( Size_t i = 0; i < jet1_npart; i++ )
          {
          track_idx = jet1_sorted_idxvec.at (i);
          track = jet1->TrackAt ( track_idx, fTracksContArray );
          if (!track) { std::cout << "Parsing tracks of jet1 :: track not defined but it should!!!" << std::endl; continue; }

          Double_t dpart = jet1->DeltaR ( track );
          Double_t track_pt = track->Pt();

          fH8->Fill   ( jet1->GetZ  ( track ) );  // Momentum distribution for leading jet (FF)
          fH8xi->Fill ( jet1->GetXi ( track ) );  // Momentum distribution for leading jet (FF) xi

          Double_t z_p = Z_ptot(jet1, track);
          fH8_p->Fill   ( z_p );  // Momentum distribution for jets (FF)
          fH8xi_p->Fill ( Xi (z_p)  );  // Momentum distribution for jets (FF) xi

          Double_t z_pt = Z_pt(jet1, track);
          fH8_pt->Fill   ( z_pt );  // Momentum distribution for jets (FF)
          fH8xi_pt->Fill ( Xi (z_pt) );  // Momentum distribution for jets (FF) xi

          fH15->Fill ( dpart, track_pt );          // <p_{T}> track vs the Distance R from Jet1

          fH20->Fill ( dpart );                    // Distribution of R in leading jet

          // fill histograms for 70% of particles with highest pt
          if ( counter_part <= jet1_n70 )
              {
              fH15_n70->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n70->Fill ( dpart );               // Distribution of R in leading jet
              }
          // fill histograms for 75% of particles with highest pt
          if ( counter_part <= jet1_n75 )
              {
              fH15_n75->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n75->Fill ( dpart );               //  Distribution of R in leading jet
              }
          // fill histograms for 80% of particles with highest pt
          if ( counter_part <= jet1_n80 )
              {
              fH15_n80->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n80->Fill ( dpart );               // Distribution of R in leading jet
              }
          // fill histograms for 85% of particles with highest pt
          if ( counter_part <= jet1_n85 )
              {
              fH15_n85->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n85->Fill ( dpart );               //  Distribution of R in leading jet
              }
          // fill histograms for 90% of particles with highest pt
          if ( counter_part <= jet1_n90 )
              {
              fH15_n90->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
              fH20_n90->Fill ( dpart );               //  Distribution of R in leading jet
              }

          // fill histograms for particles that have first 70% of pt
          if ( counter_pt <= jet1_pt70 )
              {
              fH15_pt70->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt70->Fill ( dpart );               //  Distribution of R in leading jet
              }
          // fill histograms for particles that have first 75% of pt
          if ( counter_pt <= jet1_pt75 )
              {
              fH15_pt75->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt75->Fill ( dpart );               //  Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt80 )
              {
              fH15_pt80->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt80->Fill ( dpart );     //  Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt85 )
              {
              fH15_pt85->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt85->Fill ( dpart );     //  Distribution of R in leading jet
              }
          // fill histograms for particles that have first 80% of pt
          if ( counter_pt <= jet1_pt90 )
              {
              fH15_pt90->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
              fH20_pt90->Fill ( dpart );     //  Distribution of R in leading jet
              }

          ++counter_part; counter_pt += track_pt;
          } // end of loop over jet1 tracks
      jet1_sorted_idxvec.clear();
      } // end of jet1 (leading jet) processing


    track = NULL; jet1 = NULL;

    // post data at every processing
    PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.

    // **************************************************************
    //                        ALL JETS
    // **************************************************************
    Double_t jet_pt = 0. ; UShort_t jet_npart = 0; UShort_t jet_nconst = 0;

    // vector of sorted indexes of particles in jet
    std::vector< int > jet_sorted_idxvec ;

    // loop over all jets
    for(auto jet : jetCont->accepted())
      {
      if (!jet) continue;
      Int_t track_idx = -999; // index variable for tracks

      // jet : Sorting by p_T jet constituents
      jet_sorted_idxvec = jet->SortConstituentsPt ( fTracksContArray );

      jet_pt = jet->Pt();
      jet_npart = jet->GetNumberOfTracks();
      jet_nconst = jet->GetNumberOfConstituents();

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

      counter_part = 0; counter_pt = 0.; // reset counters

      for ( Size_t i = 0; i < jet_npart; i++ )
        {
        track_idx = jet_sorted_idxvec.at (i);
        track = jet->TrackAt ( track_idx, fTracksContArray );
        if (!track) { std::cout << "Parsing tracks of jets :: track not defined but it should!!!" << std::endl; continue; }

        Double_t dpart = jet->DeltaR ( track );
        Double_t track_pt = track->Pt();

        fH8_all->Fill   ( jet->GetZ  ( track ) );  // Momentum distribution for jets (FF)
        fH8xi_all->Fill ( jet->GetXi ( track ) );  // Momentum distribution for jets (FF) xi

        Double_t z_p = Z_ptot(jet,track);
        fH8_all_p->Fill   ( z_p );  // Momentum distribution for jets (FF)
        fH8xi_all_p->Fill ( Xi (z_p)  );  // Momentum distribution for jets (FF) xi

        Double_t z_pt = Z_pt(jet,track);
        fH8_all_pt->Fill   ( z_pt );  // Momentum distribution for jets (FF)
        fH8xi_all_pt->Fill ( Xi (z_pt) );  // Momentum distribution for jets (FF) xi

        fH15all->Fill ( dpart, track_pt );         // p_{T} track vs the Distance R from jet
        fH20all->Fill ( dpart );                   // Distribution of R

        // computing components for g and ptD in the jet tracks loop
        g_tot += (track_pt * dpart)/jet_pt;
        sum_part_pt_tot += track_pt;
        sum_part_pt2_tot += TMath::Power( track_pt, 2 );

  //#############################################################################################
        if ( counter_part <= jet_n90 )         /// N90
            {
            fH15all_n90->Fill ( dpart, track_pt );         // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n90->Fill ( dpart );                   //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_n90 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n90 += track_pt;
            sum_part_pt2_n90 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt90 )         /// PT90
            {
            fH15all_pt90->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt90->Fill ( dpart );               //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_pt90 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt90 += track_pt;
            sum_part_pt2_pt90 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n85 )        /// N85
            {
            fH15all_n85->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n85->Fill ( dpart );      //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_n85 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n85 += track_pt;
            sum_part_pt2_n85 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt85 )        /// PT85
            {
            fH15all_pt85->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt85->Fill ( dpart );               //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_pt85 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt85 += track_pt;
            sum_part_pt2_pt85 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n80 )        /// N80
            {
            fH15all_n80->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n80->Fill ( dpart );      //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_n80 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n80 += track_pt;
            sum_part_pt2_n80 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt80 )         /// PT80
            {
            fH15all_pt80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt80->Fill ( dpart );               //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_pt80 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt80 += track_pt;
            sum_part_pt2_pt80 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n75 )        /// N75
            {
            fH15all_n75->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n75->Fill ( dpart );      //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_n75 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n75 += track_pt;
            sum_part_pt2_n75 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt75 )         /// PT75
            {
            fH15all_pt75->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt75->Fill ( dpart );               //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_pt75 += (track_pt * dpart)/jet_pt;
            sum_part_pt_pt75 += track_pt;
            sum_part_pt2_pt75 += TMath::Power( track_pt, 2 );
            }

  //#############################################################################################
        if ( counter_part <= jet_n70 )       /// N70
            {
            fH15all_n70->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
            fH20all_n70->Fill ( dpart );      //  Distribution of R in leading jet

            // computing components for g and ptD in the jet tracks loop
            g_n70 += (track_pt * dpart)/jet_pt;
            sum_part_pt_n70 += track_pt;
            sum_part_pt2_n70 += TMath::Power( track_pt, 2 );
            }

        if ( counter_pt <= jet_pt70 )        /// PT70
            {
            fH15all_pt70->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
            fH20all_pt70->Fill ( dpart );               //  Distribution of R in leading jet

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

      fHptd->Fill( TMath::Sqrt(sum_part_pt2_tot)/sum_part_pt_tot );

      fHptd->Fill( TMath::Sqrt(sum_part_pt2_n70)/sum_part_pt_n70 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_n75)/sum_part_pt_n75 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_n80)/sum_part_pt_n80 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_n85)/sum_part_pt_n85 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_n90)/sum_part_pt_n90 );

      fHptd->Fill( TMath::Sqrt(sum_part_pt2_pt70)/sum_part_pt_pt70 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_pt75)/sum_part_pt_pt75 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_pt80)/sum_part_pt_pt80 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_pt85)/sum_part_pt_pt85 );
      fHptd->Fill( TMath::Sqrt(sum_part_pt2_pt90)/sum_part_pt_pt90 );


      jet_sorted_idxvec.clear();


      }
      // end of loopt over all jets

    jet = NULL; track = NULL;

    }
    // end of loop over jet container collection

  // post data at every processing
  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  return kTRUE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::UserCreateOutputObjects()
  {
  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString histname = "", histtitle = "", groupname = "";
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next())))
    {
    groupname = jetCont->GetName();
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++)
      {
      //=====================================================================================
      Int_t h1_nbin = 200; Double_t h1_binwidth = 1; Double_t h1_low = 0;
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
      Int_t h5_nbin = 200; Double_t h5_binwidth = 1; Double_t h5_low = 0;
      Double_t h5_high = h5_low + h5_binwidth * h5_nbin;
      histname = TString::Format("%s/histo5_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{jets};Events", histname.Data()); // Distribution of jets in events
      fHistManager.CreateTH1(histname, histtitle, h5_nbin, h5_low, h5_high);

      //=====================================================================================
      Int_t h7_xnbin = 100; Double_t h7_xbinwidth = 1; Double_t h7_xlow = 0;
      Double_t h7_xhigh = h7_xlow + h7_xbinwidth * h7_xnbin;
      Int_t h7_ynbin = 200; Double_t h7_ybinwidth = 1; Double_t h7_ylow = 0;
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
      Int_t h8xi_nbin = 300; Double_t h8xi_binwidth = 0.05; Double_t h8xi_low = 0;
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
      Int_t h15_xnbin = 100; Double_t h15_xbinwidth = 0.01; Double_t h15_xlow = 0.;
      Double_t h15_xhigh = h15_xlow + h15_xbinwidth * h15_xnbin;
      Int_t h15_ynbin = 1000; Double_t h15_ybinwidth = 1.; Double_t h15_ylow = 0.;
      Double_t h15_yhigh = h15_ylow + h15_ybinwidth * h15_ynbin;

      //########################################################
      histname = TString::Format("%s/histo15_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - jet1;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo15all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // dR vs p_T track
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);
      //########################################################

      //=====================================================================================
      Int_t h20_nbin = 100; Double_t h20_binwidth = 0.01; Double_t h20_low = 0.;
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
      Int_t hg_nbin = 100; Double_t hg_binwidth = 0.01; Double_t hg_low = 0.;
      Double_t hg_high = hg_low + hg_binwidth * hg_nbin;

      //########################################################
      histname = TString::Format("%s/histo_g_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_g_n70%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n75%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n80%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n85%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n90%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_g_pt70%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt75%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt80%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt85%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt90%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //=====================================================================================
      // Distribution of dispersion d pt_D = sqrt ( sum (pt_i^2) )/sum (pt_i)
      Int_t hptd_nbin = 100; Double_t hptd_binwidth = 0.01; Double_t hptd_low = 0.;
      Double_t hptd_high = hptd_low + hptd_binwidth * hptd_nbin;

      //########################################################
      histname = TString::Format("%s/histo_ptd_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_ptd_n70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_n75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_n85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      //########################################################
      histname = TString::Format("%s/histo_ptd_pt70_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_pt75_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_pt85_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_ptd_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);
      //########################################################

      }
      //end of loop over fNcentBins
    }
    // end of loop over jet containers

  // =========== Switch on Sumw2 for all histos ===========
  for ( Int_t i = 0; i < fOutput->GetEntries(); ++i )
      {
      TH1 *h1 = dynamic_cast<TH1 *> ( fOutput->At (i) );
      if ( h1 ) { h1->Sumw2();continue; }

      TH2 *h2 = dynamic_cast<TH2 *> ( fOutput->At (i) );
      if ( h2 ) { h2->Sumw2();continue; }
      }

  TIter nexthist(fHistManager.GetListOfHistograms());
  TObject* obj = NULL;
  while ((obj = nexthist())) { fOutput->Add(obj); }

  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  }

//________________________________________________________________________
// Double_t AliAnalysisTaskEmcalJetCDF::DeltaR ( const AliVParticle *part1, const AliVParticle *part2 )
//   {
//   // Helper function to calculate the distance between two jets or a jet and
//   // particle
//   Double_t dPhi = part1->Phi() - part2->Phi();
//   Double_t dEta = part1->Eta() - part2->Eta();
//   dPhi = TVector2::Phi_mpi_pi ( dPhi );
//
//   return TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
//   }

//__________________________________________________________________________________________________
std::vector<Int_t> AliAnalysisTaskEmcalJetCDF::SortTracksPt ( AliVEvent *event ) const
  {
  //___________________________________________
  // Sorting by p_T (decreasing) event tracks
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Int_t entries = event->GetNumberOfTracks();

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;
  pair_list.reserve ( entries );

  for ( Int_t i_entry = 0; i_entry < entries; i_entry++ )
      {
      AliVParticle *track = event->GetTrack ( i_entry );

      if ( !track )
          {
          AliError ( Form ( "Unable to find track %d in collection %s", i_entry, event->GetName() ) );
          continue;
          }

      pair_list.push_back ( std::make_pair ( track->Pt(), i_entry ) );
      }

  std::stable_sort ( pair_list.begin(), pair_list.end(), sort_descend() );

  // return an vector of indexes of constituents (sorted descending by pt)
  std::vector<Int_t> index_sorted_list;
  index_sorted_list.reserve ( entries );

  for ( std::vector< std::pair <Double_t, Int_t> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it )
      {
      index_sorted_list.push_back ( ( *it ).second );
      } // populating the return object with indexes of sorted tracks

  return index_sorted_list;
  }

//__________________________________________________________________________________________________
std::vector<Int_t> AliAnalysisTaskEmcalJetCDF::SortTracksPt ( AliParticleContainer *trackscont ) const
  {
  //___________________________________________
  // Sorting by p_T (decreasing) event tracks
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;

  trackscont->ResetCurrentID();
  AliVTrack *track = NULL;
  UInt_t i_entry = 0;
  while( (track = dynamic_cast<AliVTrack*>(trackscont->GetNextAcceptParticle()) ))
    {
    i_entry++;
    pair_list.push_back ( std::make_pair ( track->Pt(), i_entry ) );
    }

  std::stable_sort ( pair_list.begin(), pair_list.end(), sort_descend() );

  // return an vector of indexes of constituents (sorted descending by pt)
  std::vector<Int_t> index_sorted_list;
  index_sorted_list.reserve ( i_entry );

  for ( std::vector< std::pair <Double_t, Int_t> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it )
      {
      index_sorted_list.push_back ( ( *it ).second );
      } // populating the return object with indexes of sorted tracks

  return index_sorted_list;
  }

//________________________________________________________________________
// Bool_t AliAnalysisTaskEmcalJetCDF::IdxInArray ( Int_t index, TArrayI &array )
//   {
//   for ( Int_t i = 0; i < array.GetSize(); i++ )
//       {
//       if ( index == array[i] ) { return kTRUE; }
//       }
//
//   return kFALSE;
//   }


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
Double_t AliAnalysisTaskEmcalJetCDF::Z_ptot( const AliEmcalJet* jet, const AliVParticle* trk)  const
{
  if (trk->P() < 1e-6) return 0.;
  return (trk != 0) ? trk->P()/ jet->P() : 0.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetCDF::Z_pt( const AliEmcalJet* jet, const AliVParticle* trk)  const
{
  if (trk->P() < 1e-6) return 0.;
  return (trk != 0) ? trk->Pt() / jet->Pt() : 0.;
}

//________________________________________________________________________
TObject* AliAnalysisTaskEmcalJetCDF::GetHistogram ( const char* histName )
{
  return fHistManager.FindObject(histName);
}


// kate: indent-mode none; indent-width 2; replace-tabs on;

