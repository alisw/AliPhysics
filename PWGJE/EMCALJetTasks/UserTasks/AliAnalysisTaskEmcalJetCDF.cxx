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
#include <iostream>

#include <TSystem.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TFile.h>
#include <TKey.h>
#include <TChain.h>
#include <TFileCollection.h>
#include <TCollection.h>
#include <THashList.h>
#include <TRegexp.h>
#include <TFileInfo.h>

#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TArrayD.h>
#include <TString.h>

#include <AliAODEvent.h>
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

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
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
    fUseAliEventCuts(kTRUE),
    fEventCuts(0),
    fEventCutList(0),
    fUseManualEventCuts(kFALSE),
    fHistManager()
{}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF ( const char *name ) :
    AliAnalysisTaskEmcalJet ( name, kTRUE ),
    fUseAliEventCuts(kTRUE),
    fEventCuts(0),
    fEventCutList(0),
    fUseManualEventCuts(kFALSE),
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

  namespace CDF = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS;
  TString histname = "", groupname = "", fullgroupname = "";

  TH2I* fSPclsvsSPDtrksBef = (TH2I*)fOutput->FindObject("fSPclsvsSPDtrksBef");
  TH2F* fMultV0onvsMultV0ofBef = (TH2F*)fOutput->FindObject("fMultV0onvsMultV0ofBef");

  // Per event QA; N.B. Event is already selected
  AliAODEvent* aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if (aod){
    // fSPclsvsSPDtrksBef
    Int_t nITSCls = aod->GetNumberOfITSClusters(0) + aod->GetNumberOfITSClusters(1);
    AliAODTracklets* aodTrkl = (AliAODTracklets*)aod->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    fSPclsvsSPDtrksBef->Fill(nITSTrkls, nITSCls);

    // fMultV0onvsMultV0ofBef
    AliAODVZERO* aodV0 = aod->GetVZEROData();
    Float_t  multV0Tot = aodV0->GetMTotV0A() + aodV0->GetMTotV0C();
    UShort_t multV0On  = aodV0->GetTriggerChargeA() + aodV0->GetTriggerChargeC();
    fMultV0onvsMultV0ofBef->Fill(multV0Tot, multV0On);
    }

  AliJetContainer* jetCont = NULL;
  TIter next(&fJetCollArray);
  while ( (jetCont = static_cast<AliJetContainer*>(next())) ) {
    //##### EARLY VALIDITY CHECKS - BAIL OUT FAST
    // get particles connected to jets
    AliParticleContainer* fTracksCont = dynamic_cast<AliParticleContainer*>(jetCont->GetParticleContainer());
    if (!fTracksCont) { std::cout << "*********   JET CONTAINER WITHOUT TRACKS CONTAINER   *********" << std::endl; continue; }
    TClonesArray* fTracksContArray = fTracksCont->GetArray();

    // Number of Jets found in event - accepted cuts applied by JetContainer
    Int_t fNJets_accepted = jetCont->GetNAcceptedJets();

    // Multiplicity in event - accepted tracks in tracks container
    Int_t fNaccPart = fTracksCont->GetNAcceptedParticles();

    // protection
    if ( ( fNJets_accepted < 1 ) || ( fNaccPart < 1 ) ) {
      if ( fDebug > 1 ) { std::cout << "accepted (fNJets || fNPart) == 0" << std::endl; }
      continue;
      }
    if ( fDebug > 1 ) { std::cout << "fNJets = " << fNJets_accepted << " ; fNPart = " << fNaccPart << std::endl; }


  // Only fill the embedding qa plots if:
  //  - We are using the embedding helper
  //  - The class has been initialized
  //  - Both jet collections are available
//   if (fEmbeddingQA.IsInitialized()) { fEmbeddingQA.RecordEmbeddedEventProperties(); }

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

    histname = TString::Format("%s/histo7all_%d", groupname.Data(), fCentBin);
    TH2D* fH7all = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################
    histname = TString::Format("%s/histo8_all_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8_all_pt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo8xi_all_pt_%d", groupname.Data(), fCentBin);
    TH1D* fH8xi_all_pt = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo15all_%d", groupname.Data(), fCentBin);
    TH2D* fH15all = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n80_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_n90_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_n90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt80_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo15all_pt90_%d", groupname.Data(), fCentBin);
    TH2D* fH15all_pt90 = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo20all_%d", groupname.Data(), fCentBin);
    TH1D* fH20all = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n80_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_n90_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo20all_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fH20all_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_g_%d", groupname.Data(), fCentBin);
    TH1D* fHg = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n80_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_n90_%d", groupname.Data(), fCentBin);
    TH1D* fHg_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_g_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fHg_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_ptd_%d", groupname.Data(), fCentBin);
    TH1D* fHptd = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n80_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_n90_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_ptd_pt90_%d", groupname.Data(), fCentBin);
    TH1D* fHptd_pt90 = (TH1D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_Rjt_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n80_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_n90_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_n90 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt80_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt80 = (TH2D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_Rjt_pt90_%d", groupname.Data(), fCentBin);
    TH2D* fH_Rjt_pt90 = (TH2D*)GetHistogram(histname.Data());
//######################################################################################################

    histname = TString::Format("%s/histo_jt_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n80_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n80 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_n90_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_n90 = (TH1D*)GetHistogram(histname.Data());

    histname = TString::Format("%s/histo_jt_pt80_%d", groupname.Data(), fCentBin);
    TH1D* fH_jt_pt80 = (TH1D*)GetHistogram(histname.Data());

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
    UShort_t jet_n80 = -99 ; Double_t jet_pt80 = -99.99 ;

    // variables used to compute g and ptD
    Double_t g_tot = 0.;  Double_t sum_part_pt_tot  = 0.; Double_t sum_part_pt2_tot  = 0.;
    Double_t g_n90 = 0.;  Double_t sum_part_pt_n90  = 0.; Double_t sum_part_pt2_n90  = 0.;
    Double_t g_n80 = 0.;  Double_t sum_part_pt_n80  = 0.; Double_t sum_part_pt2_n80  = 0.;
    Double_t g_pt90 = 0.; Double_t sum_part_pt_pt90 = 0.; Double_t sum_part_pt2_pt90 = 0.;
    Double_t g_pt80 = 0.; Double_t sum_part_pt_pt80 = 0.; Double_t sum_part_pt2_pt80 = 0.;


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
      g_n80 = 0.; sum_part_pt_n80 = 0.; sum_part_pt2_n80 = 0.;
      g_pt90 = 0.; sum_part_pt_pt90 = 0.; sum_part_pt2_pt90 = 0.;
      g_pt80 = 0.; sum_part_pt_pt80 = 0.; sum_part_pt2_pt80 = 0.;

      // sentinels for the pt tail cuts
      jet_n90  = ( UShort_t ) ( 0.9 * jet_npart );
      jet_pt90 = 0.9 * jet_pt;

      jet_n80  = ( UShort_t ) ( 0.8 * jet_npart );
      jet_pt80 = 0.8 * jet_pt;

      fH1->Fill ( jet_pt );            // Pt distribution of jets
      fH2->Fill ( jet->Eta() );        // Eta distribution of jets
      fH3->Fill ( jet->Phi() );        // Phi distribution of jets
      fH4->Fill ( jet_npart );         // Multiplicity of jets
      fH4c->Fill ( jet_nconst );       // Multiplicity of jets - all constituents
      fH7all->Fill ( jet_pt, jet_nconst ); // N(jet) vs P_{T} - all jets

      // parsing all jet tracks
      for (std::size_t i = 0; i < jet_npart; i++ ) {
        track_idx = jet_sorted_idxvec.at (i);
        AliVParticle* track = jet->TrackAt ( track_idx, fTracksContArray );
        if (!track) { std::cout << "Parsing tracks of jets :: track not defined but it should!!!" << std::endl; continue; }

        Double_t dpart = jet->DeltaR ( track );
        Double_t track_pt = track->Pt();
        Double_t jt = CDF::Perp (*track, *jet);

        // https://arxiv.org/abs/1209.6086 pag 4, (1)
        Double_t z_pt = CDF::Z_pt(jet,track);
        fH8_all_pt->Fill   ( z_pt );            // Momentum distribution for jets (FF)
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

        ++counter_part; counter_pt += track_pt;
        } // end of loop over jet tracks

      fHg->Fill     ( g_tot );
      fHg_n80->Fill ( g_n80 );    fHg_pt80->Fill ( g_pt80 );
      fHg_n90->Fill ( g_n90 );    fHg_pt90->Fill ( g_pt90 );

      if ( sum_part_pt_tot > 1e-8 )
        { fHptd->Fill( TMath::Sqrt(sum_part_pt2_tot)/sum_part_pt_tot ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_tot aprox 0!!!!" << endl; }

      if ( sum_part_pt_n80 > 1e-8 )
        { fHptd_n80->Fill( TMath::Sqrt(sum_part_pt2_n80)/sum_part_pt_n80 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n80 aprox 0!!!!" << endl; }

      if ( sum_part_pt_n90 > 1e-8 )
        { fHptd_n90->Fill( TMath::Sqrt(sum_part_pt2_n90)/sum_part_pt_n90 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_n90 aprox 0!!!!" << endl; }

      if ( sum_part_pt_pt80 > 1e-8 )
        { fHptd_pt80->Fill( TMath::Sqrt(sum_part_pt2_pt80)/sum_part_pt_pt80 ); }
      else
        { if ( fDebug > 2 ) cout << "sum_part_pt_pt80 aprox 0!!!!" << endl; }

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

  // Intialize AliEventCuts
  if (fUseAliEventCuts) {
    fEventCutList = new TList();
    fEventCutList->SetOwner();
    fEventCutList->SetName("EventCutOutput");
    
    fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
    if (fUseManualEventCuts) {
      fEventCuts.SetManualMode();
      fEventCuts.fMC = false;
      fEventCuts.SetupLHC15o();
      fEventCuts.fUseVariablesCorrelationCuts = true;
      }
    fEventCuts.AddQAplotsToList(fEventCutList);
    fOutput->Add(fEventCutList);
    }
  
  // Get the MC particle branch, in case it exists
  fGeneratorLevel = GetMCParticleContainer("mcparticles");

  // Initialize embedding QA
  const AliAnalysisTaskEmcalEmbeddingHelper* embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if ( embeddingHelper ) {
    if (fEmbeddingQA.Initialize()) { fEmbeddingQA.AddQAPlotsToList(fOutput); }
    }

  TH2I* fSPclsvsSPDtrksBef = new TH2I("fSPclsvsSPDtrksBef", "fSPclsvsSPDtrksBef;SPD N_{tracklets};SPD N_{clusters}", 1000, -0.5, 6999.5, 1000, -0.5, 24999.5);
  fOutput->Add(fSPclsvsSPDtrksBef);

  TH2F* fMultV0onvsMultV0ofBef = new TH2F("fMultV0onvsMultV0ofBef", "fMultV0onvsMultV0ofBef;V0 offline;V0 online", 1000, 0, 50000, 1000, 0, 50000);
  fOutput->Add(fMultV0onvsMultV0ofBef);

  TString histname = "", histtitle = "", groupname = "", fullgroupname = "";
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();

    Double_t jet_pt_min = jetCont->GetMinPt();
    Double_t jet_pt_max = jetCont->GetMaxPt();

    TString jetstrmin = TString::Itoa((Int_t)jet_pt_min,10);
    TString jetstrmax = TString::Itoa((Int_t)jet_pt_max,10);

    // add to groupname the min,max pt cuts of jets in the container
    groupname = groupname + "_" + "ptbin" + "_" + jetstrmin + "_" + jetstrmax;

    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
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

      //#####################################

      //=====================================================================================
      Int_t h5_nbin = 50; Double_t h5_binwidth = 1; Double_t h5_low = 0;
      Double_t h5_high = h5_low + h5_binwidth * h5_nbin;
      histname = TString::Format("%s/histo5_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;N_{jets};Events", histname.Data()); // Distribution of jets in events
      fHistManager.CreateTH1(histname, histtitle, h5_nbin, h5_low, h5_high);

      //=====================================================================================
      Int_t h7_xnbin = 200; Double_t h7_xbinwidth = 1; Double_t h7_xlow = 0;
      Double_t h7_xhigh = h7_xlow + h7_xbinwidth * h7_xnbin;
      Int_t h7_ynbin = 100; Double_t h7_ybinwidth = 1; Double_t h7_ylow = 0;
      Double_t h7_yhigh = h7_ylow + h7_ybinwidth * h7_ynbin;

      histname = TString::Format("%s/histo7all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/c);N_{tracks}(jets)", histname.Data()); // N vs pt all jets
      fHistManager.CreateTH2(histname, histtitle, h7_xnbin, h7_xlow, h7_xhigh, h7_ynbin, h7_ylow, h7_yhigh);

      //=====================================================================================
      Int_t h8_nbin = 101; Double_t h8_binwidth = 0.01; Double_t h8_low = 0;
      Double_t h8_high = h8_low + h8_binwidth * h8_nbin;

      // Pt implementation of Z
      histname = TString::Format("%s/histo8_all_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets Pt;z = p_{T,track}/p_{T,jet1};F(Z) = 1/N_{jets} dN/dz", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8_nbin, h8_low, h8_high);

      //=====================================================================================
      Int_t h8xi_nbin = 140; Double_t h8xi_binwidth = 0.05; Double_t h8xi_low = 0;
      Double_t h8xi_high = h8xi_low + h8xi_binwidth * h8xi_nbin;

      histname = TString::Format("%s/histo8xi_all_pt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets Pt;#xi = ln(1/z);D(#xi) = 1/N_{jets} dN/d#xi", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h8xi_nbin, h8xi_low, h8xi_high);

      //=====================================================================================
      Int_t h15_xnbin = 60; Double_t h15_xbinwidth = 0.01; Double_t h15_xlow = 0.;
      Double_t h15_xhigh = h15_xlow + h15_xbinwidth * h15_xnbin;
      Int_t h15_ynbin = 150; Double_t h15_ybinwidth = 1.; Double_t h15_ylow = 0.;
      Double_t h15_yhigh = h15_ylow + h15_ybinwidth * h15_ynbin;

      histname = TString::Format("%s/histo15all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //########################################################
      histname = TString::Format("%s/histo15all_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      histname = TString::Format("%s/histo15all_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;dR;#it{p}_{T,track} (GeV/c)", histname.Data()); // p_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h15_xnbin, h15_xlow, h15_xhigh, h15_ynbin, h15_ylow, h15_yhigh);

      //=====================================================================================
      Int_t h20_nbin = 60; Double_t h20_binwidth = 0.01; Double_t h20_low = 0.;
      Double_t h20_high = h20_low + h20_binwidth * h20_nbin;

      //########################################################
      histname = TString::Format("%s/histo20all_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      //########################################################
      histname = TString::Format("%s/histo20all_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      //########################################################
      histname = TString::Format("%s/histo20all_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      histname = TString::Format("%s/histo20all_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;R_{tracks};dN/dR", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, h20_nbin, h20_low, h20_high);

      //=====================================================================================
      // Distribution of girth (radial girth) g = sum_jet_parts ( r_i * ( pt_i/pt_jet ) )
      Int_t hg_nbin = 60; Double_t hg_binwidth = 0.005; Double_t hg_low = 0.;
      Double_t hg_high = hg_low + hg_binwidth * hg_nbin;

      //########################################################
      histname = TString::Format("%s/histo_g_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      //########################################################
      histname = TString::Format("%s/histo_g_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      //########################################################
      histname = TString::Format("%s/histo_g_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      histname = TString::Format("%s/histo_g_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;g;1/N_{jets} dN/dg", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hg_nbin, hg_low, hg_high);

      //=====================================================================================
      // Distribution of dispersion d pt_D = sqrt ( sum (pt_i^2) )/sum (pt_i)
      Int_t hptd_nbin = 100; Double_t hptd_binwidth = 0.01; Double_t hptd_low = 0.;
      Double_t hptd_high = hptd_low + hptd_binwidth * hptd_nbin;

      //########################################################
      histname = TString::Format("%s/histo_ptd_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      //########################################################
      histname = TString::Format("%s/histo_ptd_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);
      //########################################################

      histname = TString::Format("%s/histo_ptd_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

      histname = TString::Format("%s/histo_ptd_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s - all jets;ptd;1/N_{jets} dN/dp_{T}D", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, hptd_nbin, hptd_low, hptd_high);

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
      histname = TString::Format("%s/histo_Rjt_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      histname = TString::Format("%s/histo_Rjt_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;dR;j_{T} (GeV/c);", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH2(histname, histtitle, h_Rjt_xnbin, h_Rjt_xlow, h_Rjt_xhigh, h_Rjt_ynbin, h_Rjt_ylow, h_Rjt_yhigh);

      //########################################################
      histname = TString::Format("%s/histo_Rjt_pt80_%d", groupname.Data(), cent);
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
      histname = TString::Format("%s/histo_jt_n80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_n90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      //########################################################
      histname = TString::Format("%s/histo_jt_pt80_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);

      histname = TString::Format("%s/histo_jt_pt90_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s ;j_{T} (GeV/c);1/N_{jets} dN/dj_{T};", histname.Data()); // j_T track vs dR
      fHistManager.CreateTH1(histname, histtitle, h_jt_xnbin, h_jt_xlow, h_jt_xhigh);
      //########################################################

      } //end of loop over fNcentBins
    }// end of loop over jet containers


  // =========== Switch on Sumw2 for all histos ===========
  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);

  // add all fHistManager content to fOutput
  TIter nexthist(fHistManager.GetListOfHistograms());
  TObject* obj = NULL;
  while ((obj = nexthist())) { fOutput->Add(obj); }

  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  }

/*
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetCDF::ExecOnce() {
  AliAnalysisTaskEmcalJet::ExecOnce();
  }

/*
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetCDF::Terminate ( Option_t * ) {
  }

//________________________________________________________________________
TObject* AliAnalysisTaskEmcalJetCDF::GetHistogram ( const char* histName ) {
  return fHistManager.FindObject(histName);
}

// This function (overloading the base class) uses AliEventCuts to perform event selection.
//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::IsEventSelected() {
  if (fUseAliEventCuts) {
    if (!fEventCuts.AcceptEvent(InputEvent())) {
      PostData(1, fOutput);
      return kFALSE;
      }
    }
  else {
    Bool_t answer = AliAnalysisTaskEmcal::IsEventSelected();
    return answer;
    }
  return kTRUE;
}


//########################################################################
//   Namespace AliAnalysisTaskEmcalJetCDF
//########################################################################

//__________________________________________________________________________________________________
std::vector<Int_t> PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS::SortTracksPt ( AliVEvent* event )
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
std::vector<Int_t> PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS::SortTracksPt ( AliParticleContainer* trackscont )
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

//__________________________________________________________________________________________________
AliAnalysisTaskEmcalJetCDF* PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS::AddTaskEmcalJetCDF ( const char* ntracks, const char* nclusters, const char* ncells, const char* ntracks_mc, const char* tag)
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
  TString suffix    (tag);
  TString tracks    (ntracks);
  TString clusters  (nclusters);
  TString tracks_mc (ntracks_mc);
  TString cells     (ncells);

  if ( tracks.EqualTo("usedefault") ) {
    if ( dataType == kESD ) { tracks = "Tracks"; }
    else
    if ( dataType == kAOD ) { tracks = "tracks"; }
    else
      { tracks = ""; }
    }

  if ( clusters.EqualTo("usedefault") ) {
    if ( dataType == kESD ) { clusters = "CaloClusters"; }
    else
    if ( dataType == kAOD ) { clusters = "caloClusters"; }
    else
      { clusters = ""; }
    }

  if ( cells.EqualTo("usedefault") ) {
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
      cdfTask->AddMCParticleContainer(tracks.Data()); }
  else if ( tracks.EqualTo("tracks") || tracks.EqualTo("Tracks") ) {
      cdfTask->AddTrackContainer(tracks.Data()); }
  else if ( !tracks.IsNull()) {
      cdfTask->AddParticleContainer(tracks.Data()); }

  // Add the generator-level container, if specified
  if ( !tracks_mc.IsNull() ) {
    AliMCParticleContainer* mcpartCont = cdfTask->AddMCParticleContainer(tracks_mc.Data());
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    }

  if (!clusters.IsNull()) { cdfTask->AddClusterContainer(clusters.Data()); }

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

//__________________________________________________________________________________________________
  AliJetContainer* PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS::jetContSetParams ( AliJetContainer* jetCont, Float_t jetptmin,  Float_t jetptmax, Float_t jetareacutperc, Int_t leadhadtype, Int_t nLeadJets, Float_t mintrackpt, Float_t maxtrackpt)
    {
    if (!jetCont) { return NULL; }
    jetCont->SetJetPtCut ( jetptmin );
    jetCont->SetJetPtCutMax ( jetptmax );
    jetCont->SetPercAreaCut ( jetareacutperc );
    jetCont->SetLeadingHadronType ( leadhadtype ); // 0 = charged, 1 = neutral, 2 = both
    jetCont->SetNLeadingJets(nLeadJets);
    jetCont->SetMinTrackPt(mintrackpt);
    jetCont->SetMaxTrackPt(maxtrackpt);
    return jetCont;
    }

//__________________________________________________________________________________________________
TChain* PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetCDF_NS::CreateChain ( const char* filelist, const char* cTreeNameArg, const char* friends, UInt_t iNumFiles, UInt_t iStartWithFile ) {
TString sTreeNameArg (cTreeNameArg), treeName;

TFileCollection filecoll ("anachain","File collection for analysis"); // easy manipulation of file collections
Int_t iAddedFiles = filecoll.AddFromFile(filelist,iNumFiles,iStartWithFile);
if ( iAddedFiles < 1 ) { std::cout << "NO Files added to collection !!!" << std::endl; return NULL; }

// if cTreeNameArg is auto lets try to autodetect what type of tree we have;
// the assuption is that all files are the same and the first one is reprezentative
THashList* list =  filecoll.GetList();
if ( sTreeNameArg.EqualTo("auto") ) { // if tree name is not specified
  TRegexp tree_regex ("[aod,esd]Tree");
  TFileInfo* fileinfo = dynamic_cast<TFileInfo*>(list->At(0)); // get first fileinfo in list
  TFile file (fileinfo->GetFirstUrl()->GetFile()); // get the actual TFile
  if (file.IsZombie()) { cout << "Should not reach this message!! Error opening file" << endl; return NULL; }

  // lets parse the TFile
  TIter next(file.GetListOfKeys());
  TKey* key = NULL;
  while (( key = dynamic_cast<TKey*>(next()) )) {
    TString class_name = key->GetClassName();
    if ( ! class_name.EqualTo("TTree") ) { continue; } // searching for first TTree

    TString key_name = key->GetName();
    if ( key_name.Contains(tree_regex) ) { treeName = key_name; break;} // that is named either aodTree or esdTree
    }
  file.Close();
  }
else
  { treeName = sTreeNameArg ; } // tree name is specified

TChain* chain = new TChain (treeName.Data(),""); // lets create the chain
if ( chain->AddFileInfoList(list) == 0 ) { return NULL; } // and add file collection (THashList is a Tlist that is a TCollection)

// start treatment of friends
TChain* chainFriend = NULL;
TString sFriends (friends);
if (!sFriends.IsNull()) {
  TString friends_treename, friends_filename;
  TObjArray* arr = sFriends.Tokenize("/");
  TObjString* strobj_file = dynamic_cast<TObjString*>(arr->At(0));
  if (strobj_file) { friends_filename = strobj_file->GetString(); }
  TObjString* strobj_tree = dynamic_cast<TObjString*>(arr->At(1));
  if (strobj_tree) { friends_treename = strobj_tree->GetString(); }
  delete arr;

  if (friends_treename.IsNull()) {
    if (treeName.EqualTo("esdTree")) {
      friends_treename = "esdFriendTree"; }
    else if (treeName.EqualTo("aodTree")) {
      friends_treename = "aodTree"; }
    else {
      cout << "friends argument specified but tree name neither specified nor auto-detected (unknown tree name to associate with known friend tree name)";
      return chain; // stop processing of friends, just return the chain created so far
      }
    }

  chainFriend = new TChain(friends_treename.Data());
  TString friendinfo_for_chain = "/" + friends_filename + "/" + friends_treename;

  TIter next(list);
  TFileInfo* fileinfo = NULL;
  while (( fileinfo = dynamic_cast<TFileInfo*>(next()) )) {
    TString dirname = gSystem->DirName(fileinfo->GetFirstUrl()->GetFile());
    TString friend_for_chain = dirname + friendinfo_for_chain;
    if (chainFriend) { chainFriend->Add(friend_for_chain.Data()); }
    }

  if (chainFriend) { chain->AddFriend(chainFriend); }
  }

return chain;
}

// kate: indent-mode none; indent-width 2; replace-tabs on;

