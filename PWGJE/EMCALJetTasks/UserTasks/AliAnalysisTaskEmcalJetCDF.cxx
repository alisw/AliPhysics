#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>

#include <Rtypes.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TClonesArray.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TSeqCollection.h>
#include <TCollection.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TAttMarker.h>

#include <AliVCluster.h>
#include <AliAODCaloCluster.h>
#include <AliESDCaloCluster.h>
#include <AliVTrack.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliLog.h"

#include "AliAnalysisTaskEmcalJetCDF.h"

//using namespace std;

/// \cond CLASSIMP
ClassImp ( AliAnalysisTaskEmcalJetCDF );
/// \endcond

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF()
  : AliAnalysisTaskEmcalJet ( "AliAnalysisTaskEmcalJetCDF", kTRUE ),
    fH1( NULL ),
    fH2( NULL ),
    fH3( NULL ),
    fH4( NULL ),
    fH4c( NULL ),
    fH5( NULL ),
    fH6( NULL ),
    fH6c( NULL ),
    fH8( NULL ),
    fH8xi( NULL ),
    fH8_all( NULL ),
    fH8xi_all( NULL ),
    fH15_bin( NULL ),
    fH15_bin_n90( NULL ),
    fH15_bin_n85( NULL ),
    fH15_bin_n80( NULL ),
    fH15_bin_n75( NULL ),
    fH15_bin_n70( NULL ),
    fH15_bin_pt90( NULL ),
    fH15_bin_pt85( NULL ),
    fH15_bin_pt80( NULL ),
    fH15_bin_pt75( NULL ),
    fH15_bin_pt70( NULL ),
    fH15all_bin( NULL ),
    fH15all_bin_n90( NULL ),
    fH15all_bin_n85( NULL ),
    fH15all_bin_n80( NULL ),
    fH15all_bin_n75( NULL ),
    fH15all_bin_n70( NULL ),
    fH15all_bin_pt90( NULL ),
    fH15all_bin_pt85( NULL ),
    fH15all_bin_pt80( NULL ),
    fH15all_bin_pt75( NULL ),
    fH15all_bin_pt70( NULL ),
    fH20( NULL ),
    fH20_n90( NULL ),
    fH20_n85( NULL ),
    fH20_n80( NULL ),
    fH20_n75( NULL ),
    fH20_n70( NULL ),
    fH20_pt90( NULL ),
    fH20_pt85( NULL ),
    fH20_pt80( NULL ),
    fH20_pt75( NULL ),
    fH20_pt70( NULL ),
    fH20all( NULL ),
    fH20all_n90( NULL ),
    fH20all_n85( NULL ),
    fH20all_n80( NULL ),
    fH20all_n75( NULL ),
    fH20all_n70( NULL ),
    fH20all_pt90( NULL ),
    fH20all_pt85( NULL ),
    fH20all_pt80( NULL ),
    fH20all_pt75( NULL ),
    fH20all_pt70( NULL ),
    fHg( NULL ),
    fHg_n90( NULL ),
    fHg_n85( NULL ),
    fHg_n80( NULL ),
    fHg_n75( NULL ),
    fHg_n70( NULL ),
    fHg_pt90( NULL ),
    fHg_pt85( NULL ),
    fHg_pt80( NULL ),
    fHg_pt75( NULL ),
    fHg_pt70( NULL ),
    fHptd( NULL ),
    fHptd_n90( NULL ),
    fHptd_n85( NULL ),
    fHptd_n80( NULL ),
    fHptd_n75( NULL ),
    fHptd_n70( NULL ),
    fHptd_pt90( NULL ),
    fHptd_pt85( NULL ),
    fHptd_pt80( NULL ),
    fHptd_pt75( NULL ),
    fHptd_pt70( NULL ),
    fJetsCont ( NULL ),
    fTracksCont ( NULL ),
    fCaloClustersCont ( NULL ),
    fTracksContArray ( NULL ),
    idx_jetcont(0),
    fNJets_accepted ( 0 ),
    fNaccPart ( 0 ),
    fNaccClus ( 0 ),
    fHistManager("AliAnalysisTaskEmcalJetCDF")
  {
  // Default constructor.
  fDebug = AliLog::GetGlobalDebugLevel();
  // SetMakeGeneralHistograms ( kTRUE );
  }

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF ( const char *name )
  : AliAnalysisTaskEmcalJet ( name, kTRUE ),
    fH1( NULL ),
    fH2( NULL ),
    fH3( NULL ),
    fH4( NULL ),
    fH4c( NULL ),
    fH5( NULL ),
    fH6( NULL ),
    fH6c( NULL ),
    fH8( NULL ),
    fH8xi( NULL ),
    fH8_all( NULL ),
    fH8xi_all( NULL ),
    fH15_bin( NULL ),
    fH15_bin_n90( NULL ),
    fH15_bin_n85( NULL ),
    fH15_bin_n80( NULL ),
    fH15_bin_n75( NULL ),
    fH15_bin_n70( NULL ),
    fH15_bin_pt90( NULL ),
    fH15_bin_pt85( NULL ),
    fH15_bin_pt80( NULL ),
    fH15_bin_pt75( NULL ),
    fH15_bin_pt70( NULL ),
    fH15all_bin( NULL ),
    fH15all_bin_n90( NULL ),
    fH15all_bin_n85( NULL ),
    fH15all_bin_n80( NULL ),
    fH15all_bin_n75( NULL ),
    fH15all_bin_n70( NULL ),
    fH15all_bin_pt90( NULL ),
    fH15all_bin_pt85( NULL ),
    fH15all_bin_pt80( NULL ),
    fH15all_bin_pt75( NULL ),
    fH15all_bin_pt70( NULL ),
    fH20( NULL ),
    fH20_n90( NULL ),
    fH20_n85( NULL ),
    fH20_n80( NULL ),
    fH20_n75( NULL ),
    fH20_n70( NULL ),
    fH20_pt90( NULL ),
    fH20_pt85( NULL ),
    fH20_pt80( NULL ),
    fH20_pt75( NULL ),
    fH20_pt70( NULL ),
    fH20all( NULL ),
    fH20all_n90( NULL ),
    fH20all_n85( NULL ),
    fH20all_n80( NULL ),
    fH20all_n75( NULL ),
    fH20all_n70( NULL ),
    fH20all_pt90( NULL ),
    fH20all_pt85( NULL ),
    fH20all_pt80( NULL ),
    fH20all_pt75( NULL ),
    fH20all_pt70( NULL ),
    fHg( NULL ),
    fHg_n90( NULL ),
    fHg_n85( NULL ),
    fHg_n80( NULL ),
    fHg_n75( NULL ),
    fHg_n70( NULL ),
    fHg_pt90( NULL ),
    fHg_pt85( NULL ),
    fHg_pt80( NULL ),
    fHg_pt75( NULL ),
    fHg_pt70( NULL ),
    fHptd( NULL ),
    fHptd_n90( NULL ),
    fHptd_n85( NULL ),
    fHptd_n80( NULL ),
    fHptd_n75( NULL ),
    fHptd_n70( NULL ),
    fHptd_pt90( NULL ),
    fHptd_pt85( NULL ),
    fHptd_pt80( NULL ),
    fHptd_pt75( NULL ),
    fHptd_pt70( NULL ),
    fJetsCont ( NULL ),
    fTracksCont ( NULL ),
    fCaloClustersCont ( NULL ),
    fTracksContArray ( NULL ),
    idx_jetcont(0),
    fNJets_accepted ( 0 ),
    fNaccPart ( 0 ),
    fNaccClus ( 0 ),
    fHistManager("AliAnalysisTaskEmcalJetCDF")
  {
  // Standard constructor.
  fDebug = AliLog::GetGlobalDebugLevel();
  // SetMakeGeneralHistograms ( kTRUE );
  }

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::~AliAnalysisTaskEmcalJetCDF()
  {
  // Destructor.
  if ( fJetsCont )         { delete fJetsCont ;         fJetsCont = NULL; }
  if ( fTracksCont )       { delete fTracksCont ;       fTracksCont = NULL; }
  if ( fCaloClustersCont ) { delete fCaloClustersCont ; fCaloClustersCont = NULL; }

  if ( fTracksContArray )    { delete fTracksContArray ;    fTracksContArray = NULL; }
  if ( fCaloClustContArray ) { delete fCaloClustContArray ; fCaloClustContArray = NULL; }

  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::Run()
  {
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  return ProcessJetContainer(idx_jetcont);
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::ProcessJetContainer(Int_t idx_jet_container)
  {
  fJetsCont = GetJetContainer ( idx_jet_container );
  fJetsCont->ResetCurrentID();

  if ( !fJetsCont )
      {
      std::cout << "ERROR :: Jet Container not found!!!" << std::endl;
      return kFALSE;
      }

  // get particles connected to jets
  fTracksCont = fJetsCont->GetParticleContainer();
  if (fTracksCont)
    {
    fTracksCont->ResetCurrentID();
    fTracksCont->SetClassName ( "AliVTrack" );

    // get array of tracks
    fTracksContArray = fTracksCont->GetArray();

    // Multiplicity in event - accepted tracks in tracks container
    fNaccPart = fTracksCont->GetNAcceptedParticles();
    }

  // get clusters connected to jets
  fCaloClustersCont = fJetsCont->GetClusterContainer();
  if (fCaloClustersCont)
    {
    fCaloClustersCont->SetClassName ( "AliVCluster" );

    // get array of constituents
    fCaloClustContArray = fCaloClustersCont->GetArray();

    // accepted clusters in cluster container
    fNaccClus = fCaloClustersCont->GetNAcceptedClusters();
    }

  // Number of Jets found in event - accepted cuts applied by JetContainer
  fNJets_accepted = fJetsCont->GetNJets();

  // protection
  if ( ( fNJets_accepted < 1 ) || ( fNaccPart < 1 ) )
      {
      if ( fDebug > 1 )
          {
          std::cout << "accepted (fNJets || fNPart) == 0" << std::endl;
          }
      return kFALSE;
      }

  if ( fDebug > 1 )
      {
      std::cout << "fNJets = " << fNJets_accepted << " ; fNPart = " << fNaccPart << std::endl;
      }

  // Run analysis code here, if needed. It will be executed before FillHistograms().
  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::FillHistograms()
  {
  AliVParticle* track = NULL;
  AliEmcalJet*    jet = NULL;

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
//                        ALL JETS
// **************************************************************
  // jet propreties
  Double_t jet_pt = 0. ; UShort_t jet_npart = 0; UShort_t jet_nconst = 0;

  // vector of sorted indexes of particles in jet
  std::vector< int > jet_sorted_idxvec ;

  fJetsCont->ResetCurrentID();
  jet = fJetsCont->GetNextAcceptJet ();
  while ( jet )
    {
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

      fH20all->Fill ( dpart );                   // Distribution of R in leading jet

      fH15all->Fill ( dpart, track_pt );         // p_{T} track vs the Distance R from Jet1
      fH15all_bin->Fill ( dpart, track_pt );     // p_{T} track vs the Distance R from jet axis

      fH23all->Fill ( jet_pt, dpart );           //  Jet1 Size vs P_{T} - all jets

      // computing components for g and ptD in the jet tracks loop
      g_tot += (track_pt * dpart)/jet_pt;
      sum_part_pt_tot += track_pt;
      sum_part_pt2_tot += TMath::Power( track_pt, 2 );

//#############################################################################################
      if ( counter_part <= jet_n90 )         /// N90
          {
          fH15all_n90->Fill ( dpart, track_pt );         // p_{T} track vs the Distance R from Jet - 80% of particles
          fH15all_bin_n90->Fill ( dpart, track_pt );     // p_{T} track vs the Distance R from Jet - 80% of particles
          fH20all_n90->Fill ( dpart );                   //  Distribution of R in leading jet

          // computing components for g and ptD in the jet tracks loop
          g_n90 += (track_pt * dpart)/jet_pt;
          sum_part_pt_n90 += track_pt;
          sum_part_pt2_n90 += TMath::Power( track_pt, 2 );
          }

      if ( counter_pt <= jet_pt90 )         /// PT90
          {
          fH15all_pt90->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH15all_bin_pt90->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
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
          fH15all_bin_n85->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
          fH20all_n85->Fill ( dpart );      //  Distribution of R in leading jet

          // computing components for g and ptD in the jet tracks loop
          g_n85 += (track_pt * dpart)/jet_pt;
          sum_part_pt_n85 += track_pt;
          sum_part_pt2_n85 += TMath::Power( track_pt, 2 );
          }

      if ( counter_pt <= jet_pt85 )        /// PT85
          {
          fH15all_pt85->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH15all_bin_pt85->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
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
          fH15all_bin_n80->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
          fH20all_n80->Fill ( dpart );      //  Distribution of R in leading jet
          fH24all->Fill ( jet_pt, dpart ); //  Jet1 Size vs P_{T} - all jets - 80% of particles

          // computing components for g and ptD in the jet tracks loop
          g_n80 += (track_pt * dpart)/jet_pt;
          sum_part_pt_n80 += track_pt;
          sum_part_pt2_n80 += TMath::Power( track_pt, 2 );
          }

      if ( counter_pt <= jet_pt80 )         /// PT80
          {
          fH15all_pt80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH15all_bin_pt80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH20all_pt80->Fill ( dpart );               //  Distribution of R in leading jet
          fH25all->Fill ( jet_pt, dpart );             //  Jet1 Size vs P_{T} - all jets - 80% of Pt

          // computing components for g and ptD in the jet tracks loop
          g_pt80 += (track_pt * dpart)/jet_pt;
          sum_part_pt_pt80 += track_pt;
          sum_part_pt2_pt80 += TMath::Power( track_pt, 2 );
          }

//#############################################################################################
      if ( counter_part <= jet_n75 )        /// N75
          {
          fH15all_n75->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
          fH15all_bin_n75->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
          fH20all_n75->Fill ( dpart );      //  Distribution of R in leading jet

          // computing components for g and ptD in the jet tracks loop
          g_n75 += (track_pt * dpart)/jet_pt;
          sum_part_pt_n75 += track_pt;
          sum_part_pt2_n75 += TMath::Power( track_pt, 2 );
          }

      if ( counter_pt <= jet_pt75 )         /// PT75
          {
          fH15all_pt75->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH15all_bin_pt75->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
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
          fH15all_bin_n70->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet - 80% of particles
          fH20all_n80->Fill ( dpart );      //  Distribution of R in leading jet

          // computing components for g and ptD in the jet tracks loop
          g_n70 += (track_pt * dpart)/jet_pt;
          sum_part_pt_n70 += track_pt;
          sum_part_pt2_n70 += TMath::Power( track_pt, 2 );
          }

      if ( counter_pt <= jet_pt70 )        /// PT70
          {
          fH15all_pt70->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
          fH15all_bin_pt70->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R from Jet - 80% of pt
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

    jet = fJetsCont->GetNextAcceptJet();
    jet_sorted_idxvec.clear();
    } // end of loop over jets

  jet = NULL; track = NULL;

// **************************************************************
//                          LEADING JETS
// **************************************************************
  AliEmcalJet* jet1 = fJetsCont->GetLeadingJet(); // internaly checked for AcceptedJet

  if ( fDebug > 1 )
    {
    if (jet1)
      {
      std::cout << "+++++++++++++++++>>>>>>>>> Leading jet found" << std::endl;
      jet1->Print();
      }
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

        fH20->Fill ( dpart );                    // Distribution of R in leading jet
        fH15->Fill ( dpart, track_pt );          // <p_{T}> track vs the Distance R from Jet1
        fH15_bin->Fill ( dpart, track_pt );      // p_{T} sum track vs the Distance R from Jet1

        fH23->Fill ( jet1_pt, dpart );           //  Jet1 Size vs P_{T}(jet1)

        // fill histograms for 70% of particles with highest pt
        if ( counter_part <= jet1_n70 )
            {
            fH15_n70->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
            fH15_bin_n70->Fill ( dpart, track_pt ); // p_{T} sum track vs the Distance R from Jet1 - 80% of particles
            fH20_n70->Fill ( dpart );               // Distribution of R in leading jet
            }
        // fill histograms for 75% of particles with highest pt
        if ( counter_part <= jet1_n75 )
            {
            fH15_n75->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
            fH15_bin_n75->Fill ( dpart, track_pt ); // p_{T} sum track vs the Distance R from Jet1 - 80% of particles
            fH20_n75->Fill ( dpart );               //  Distribution of R in leading jet
            }
        // fill histograms for 80% of particles with highest pt
        if ( counter_part <= jet1_n80 )
            {
            fH15_n80->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
            fH15_bin_n80->Fill ( dpart, track_pt ); // p_{T} sum track vs the Distance R from Jet1 - 80% of particles
            fH20_n80->Fill ( dpart );               // Distribution of R in leading jet
            fH24->Fill ( jet1_pt, dpart );          // Jet1 Size vs P_{T}(jet1) - 80% of particles
            }
        // fill histograms for 85% of particles with highest pt
        if ( counter_part <= jet1_n85 )
            {
            fH15_n85->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
            fH15_bin_n85->Fill ( dpart, track_pt ); // p_{T} sum track vs the Distance R from Jet1 - 80% of particles
            fH20_n85->Fill ( dpart );               //  Distribution of R in leading jet
            }
        // fill histograms for 90% of particles with highest pt
        if ( counter_part <= jet1_n90 )
            {
            fH15_n90->Fill ( dpart, track_pt );     // <p_{T}> track vs the Distance R from Jet1 - 80% of particles
            fH15_bin_n90->Fill ( dpart, track_pt ); // p_{T} sum track vs the Distance R from Jet1 - 80% of particles
            fH20_n90->Fill ( dpart );               //  Distribution of R in leading jet
            }

        // fill histograms for particles that have first 70% of pt
        if ( counter_pt <= jet1_pt70 )
            {
            fH15_pt70->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
            fH15_bin_pt70->Fill ( dpart, track_pt ); //  p_{T} sum track vs the Distance R from Jet1 - 80% of pt
            fH20_pt70->Fill ( dpart );               //  Distribution of R in leading jet
            }
        // fill histograms for particles that have first 75% of pt
        if ( counter_pt <= jet1_pt75 )
            {
            fH15_pt75->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
            fH15_bin_pt75->Fill ( dpart, track_pt ); //  p_{T} sum track vs the Distance R from Jet1 - 80% of pt
            fH20_pt75->Fill ( dpart );               //  Distribution of R in leading jet
            }
        // fill histograms for particles that have first 80% of pt
        if ( counter_pt <= jet1_pt80 )
            {
            fH15_pt80->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
            fH15_bin_pt80->Fill ( dpart, track_pt ); //  p_{T} sum track vs the Distance R from Jet1 - 80% of pt
            fH20_pt80->Fill ( dpart );     //  Distribution of R in leading jet
            fH25->Fill ( jet1_pt, dpart ); //  Jet1 Size vs P_{T}(jet1) - 80% of Pt
            }
        // fill histograms for particles that have first 80% of pt
        if ( counter_pt <= jet1_pt85 )
            {
            fH15_pt85->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
            fH15_bin_pt85->Fill ( dpart, track_pt ); //  p_{T} sum track vs the Distance R from Jet1 - 80% of pt
            fH20_pt85->Fill ( dpart );     //  Distribution of R in leading jet
            }
        // fill histograms for particles that have first 80% of pt
        if ( counter_pt <= jet1_pt90 )
            {
            fH15_pt90->Fill ( dpart, track_pt );     //  <p_{T}> track vs the Distance R from Jet1 - 80% of pt
            fH15_bin_pt90->Fill ( dpart, track_pt ); //  p_{T} sum track vs the Distance R from Jet1 - 80% of pt
            fH20_pt90->Fill ( dpart );     //  Distribution of R in leading jet
            }

        ++counter_part; counter_pt += track_pt;
        } // end of loop over jet1 tracks
    jet1_sorted_idxvec.clear();
    } // end of jet1 (leading jet) processing


  track = NULL; jet1 = NULL;

  // post data at every processing
  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  return kTRUE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::UserCreateOutputObjects()
  {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  //  Histograms
  fOutput->SetOwner ( kTRUE );

  // Create the list of histograms. Only the list is owned.

  //____________________________________________________________________________________
  Int_t h1_nbin = 300;
  Double_t h1_binwidth = 1;
  Double_t h1_low = 0;
  Double_t h1_high = h1_low + h1_binwidth * h1_nbin; // 1GeV/bin
  fH1 = new TH1D ( "histo1", "p_{T} distribution of jets (accepted)", h1_nbin, h1_low, h1_high );
  fH1->SetStats ( kTRUE );
  fH1->GetXaxis()->SetTitle ( "p_{T,jet} in GeV/c" );
  fH1->GetYaxis()->SetTitle ( "Jets" );
  fH1->GetXaxis()->SetTitleColor ( 1 );
  fH1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH1 );

  Int_t h2_nbin = 200;
  Double_t h2_binwidth = 0.01;
  Double_t h2_low = -1;
  Double_t h2_high = h2_low + h2_binwidth * h2_nbin;
  fH2 = new TH1D ( "histo2", "#eta distribution of jets (accepted)", h2_nbin, h2_low,
                   h2_high ); // 1 unit of rapidity / 100 bin
  fH2->SetStats ( kTRUE );
  fH2->GetXaxis()->SetTitle ( "#eta_{jet}" );
  fH2->GetYaxis()->SetTitle ( "Jets" );
  fH2->GetXaxis()->SetTitleColor ( 1 );
  fH2->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH2 );

  Int_t h3_nbin = 126;
  Double_t h3_binwidth = 0.05;
  Double_t h3_low = 0.;
  Double_t h3_high = h3_low + h3_binwidth * h3_nbin;
  fH3 = new TH1D ( "histo3", "#phi distribution of jets (accepted)", h3_nbin, h3_low, h3_high );
  fH3->SetStats ( kTRUE );
  fH3->GetXaxis()->SetTitle ( "#phi_{jet}" );
  fH3->GetYaxis()->SetTitle ( "Jets" );
  fH3->GetXaxis()->SetTitleColor ( 1 );
  fH3->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH3 );

  //____________________________________________________________________________________
  Int_t h4_nbin = 100;
  Double_t h4_binwidth = 1;
  Double_t h4_low = 0;
  Double_t h4_high = h4_low + h4_binwidth * h4_nbin;
  fH4 = new TH1D ( "histo4", "Multiplicity of jets (accepted) - charged tracks", h4_nbin, h4_low, h4_high ); // 1 unit of multiplicity /bin
  fH4->SetStats ( kTRUE );
  fH4->GetXaxis()->SetTitle ( "N_{tracks}(jet)" );
  fH4->GetYaxis()->SetTitle ( "Jets" );
  fH4->GetXaxis()->SetTitleColor ( 1 );
  fH4->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH4 );

  fH4c = new TH1D ( "histo4c", "Multiplicity of jets (accepted) - all constituents", h4_nbin, h4_low, h4_high ); // 1 unit of multiplicity /bin
  fH4c->SetStats ( kTRUE );
  fH4c->GetXaxis()->SetTitle ( "N_{tracks}(jet)" );
  fH4c->GetYaxis()->SetTitle ( "Jets" );
  fH4c->GetXaxis()->SetTitleColor ( 1 );
  fH4c->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH4c );

  //____________________________________________________________________________________
  Int_t h5_nbin = 100;
  Double_t h5_binwidth = 1;
  Double_t h5_low = 0;
  Double_t h5_high = h5_low + h5_binwidth * h5_nbin;
  fH5 = new TH1D ( "histo5", "Distribution of jets in events", h5_nbin, h5_low, h5_high );
  fH5->SetStats ( kTRUE );
  fH5->GetXaxis()->SetTitle ( "N_{jets}" );
  fH5->GetYaxis()->SetTitle ( "Events" );
  fH5->GetXaxis()->SetTitleColor ( 1 );
  fH5->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH5 );

  //____________________________________________________________________________________
  Int_t h6_nbin = 100;
  Double_t h6_binwidth = 1;
  Double_t h6_low = 0;
  Double_t h6_high = h6_low + h6_binwidth * h6_nbin;
  fH6 = new TH1D ( "histo6", "Jet1 Multiplicity Distribution - charged tracks", h6_nbin, h6_low, h6_high );
  fH6->SetStats ( kTRUE );
  fH6->GetXaxis()->SetTitle ( "N_{tracks}(jet1)" );
  fH6->GetYaxis()->SetTitle ( "1/N_{jet1}" );
  fH6->GetXaxis()->SetTitleColor ( 1 );
  fH6->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH6 );

  fH6c = new TH1D ( "histo6c", "Jet1 Multiplicity Distribution - all constituents", h6_nbin, h6_low, h6_high );
  fH6c->SetStats ( kTRUE );
  fH6c->GetXaxis()->SetTitle ( "N_{tracks}(jet1)" );
  fH6c->GetYaxis()->SetTitle ( "1/N_{jet1}" );
  fH6c->GetXaxis()->SetTitleColor ( 1 );
  fH6c->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH6c );

  //____________________________________________________________________________________
  Int_t h7_nbin = 400;
  Double_t h7_binwidth = 1;
  Double_t h7_xlow = 0;
  Double_t h7_xhigh = h7_xlow + h7_binwidth * h7_nbin;
  fH7 = new TProfile ( "histo7", "N(jet1) vs P_{T} - jet1", h7_nbin, h7_xlow, h7_xhigh );
  fH7->SetStats ( kTRUE );
  fH7->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH7->GetYaxis()->SetTitle ( "<N_{tracks}(jet1)> in 1 GeV/c bin" );
  fH7->GetXaxis()->SetTitleColor ( 1 );
  fH7->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH7 );

  fH7all = new TProfile ( "histo7all", "N(jet1) vs P_{T} - all jets", h7_nbin, h7_xlow, h7_xhigh );
  fH7all->SetStats ( kTRUE );
  fH7all->GetXaxis()->SetTitle ( "p_{T}(jet) (GeV/c)" );
  fH7all->GetYaxis()->SetTitle ( "<N_{tracks}(jet1)> in 1 GeV/c bin" );
  fH7all->GetXaxis()->SetTitleColor ( 1 );
  fH7all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH7all );

  //____________________________________________________________________________________
  Int_t h8_nbin = 101;
  Double_t h8_binwidth = 0.01;
  Double_t h8_low = 0;
  Double_t h8_high = h8_low + h8_binwidth * h8_nbin;
  fH8 = new TH1D ( "histo8", "Momentum distribution for jet1 (FF)", h8_nbin, h8_low, h8_high );
  fH8->SetStats ( kTRUE );
  fH8->GetXaxis()->SetTitle ( "z = p_{T,track}/p_{T,jet1}" );
  fH8->GetYaxis()->SetTitle ( "F(Z) = 1/N_{jets1} dN/dz" );
  fH8->GetXaxis()->SetTitleColor ( 1 );
  fH8->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH8 );

  fH8_all = new TH1D ( "histo8_all", "Momentum distribution for jets (FF)", h8_nbin, h8_low, h8_high );
  fH8_all->SetStats ( kTRUE );
  fH8_all->GetXaxis()->SetTitle ( "z = p_{T,track}/p_{T,jet1}" );
  fH8_all->GetYaxis()->SetTitle ( "F(Z) = 1/N_{jets1} dN/dz" );
  fH8_all->GetXaxis()->SetTitleColor ( 1 );
  fH8_all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH8_all );

  Int_t h8xi_nbin = 300;
  Double_t h8xi_binwidth = 0.05;
  Double_t h8xi_low = 0;
  Double_t h8xi_high = h8xi_low + h8xi_binwidth * h8xi_nbin;
  fH8xi = new TH1D ( "histo8xi", "Momentum distribution for jet1 (FF)", h8xi_nbin, h8xi_low, h8xi_high );
  fH8xi->SetStats ( kTRUE );
  fH8xi->GetXaxis()->SetTitle ( "#xi = ln(1/z)" );
  fH8xi->GetYaxis()->SetTitle ( "D(#xi) = 1/N_{jets1} dN/d#xi" );
  fH8xi->GetXaxis()->SetTitleColor ( 1 );
  fH8xi->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH8xi );

  fH8xi_all = new TH1D ( "histo8xi_all", "Momentum distribution for all jets (FF)", h8xi_nbin, h8xi_low, h8xi_high );
  fH8xi_all->SetStats ( kTRUE );
  fH8xi_all->GetXaxis()->SetTitle ( "#xi = ln(1/z)" );
  fH8xi_all->GetYaxis()->SetTitle ( "D(#xi) = 1/N_{jets1} dN/d#xi" );
  fH8xi_all->GetXaxis()->SetTitleColor ( 1 );
  fH8xi_all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH8xi_all );

  //____________________________________________________________________________________
  Int_t h15_nbin = 100;
  Double_t h15_binwidth = 0.01;
  Double_t h15_xlow = 0.;
  Double_t h15_xhigh = h15_xlow + h15_binwidth * h15_nbin;

  fH15 = new TProfile ( "histo15", "<p_{T}> track in dR(jet1)", h15_nbin, h15_xlow, h15_xhigh );
  fH15->SetStats ( kTRUE );
  fH15->GetXaxis()->SetTitle ( "R" );
  fH15->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15->GetXaxis()->SetTitleColor ( 1 );
  fH15->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15 );

//=====================================================================================
  fH15_n70 = new TProfile ( "histo15_n70", "<p_{T}> track in dR(jet1) - 70% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_n70->SetStats ( kTRUE );
  fH15_n70->GetXaxis()->SetTitle ( "R" );
  fH15_n70->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n70->GetXaxis()->SetTitleColor ( 1 );
  fH15_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n70 );

  fH15_n75 = new TProfile ( "histo15_n75", "<p_{T}> track in dR(jet1) - 75% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_n75->SetStats ( kTRUE );
  fH15_n75->GetXaxis()->SetTitle ( "R" );
  fH15_n75->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n75->GetXaxis()->SetTitleColor ( 1 );
  fH15_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n75 );

  fH15_n80 = new TProfile ( "histo15_n80", "<p_{T}> track in dR(jet1) - 80% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_n80->SetStats ( kTRUE );
  fH15_n80->GetXaxis()->SetTitle ( "R" );
  fH15_n80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n80 );

  fH15_n85 = new TProfile ( "histo15_n85", "<p_{T}> track in dR(jet1) - 85% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_n85->SetStats ( kTRUE );
  fH15_n85->GetXaxis()->SetTitle ( "R" );
  fH15_n85->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n85->GetXaxis()->SetTitleColor ( 1 );
  fH15_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n85 );

  fH15_n90 = new TProfile ( "histo15_n90", "<p_{T}> track in dR(jet1) - 90% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_n90->SetStats ( kTRUE );
  fH15_n90->GetXaxis()->SetTitle ( "R" );
  fH15_n90->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n90->GetXaxis()->SetTitleColor ( 1 );
  fH15_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n90 );


//=====================================================================================
  fH15_pt70 = new TProfile ( "histo15_pt70", "<p_{T}> track in dR(jet1) - 70% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_pt70->SetStats ( kTRUE );
  fH15_pt70->GetXaxis()->SetTitle ( "R" );
  fH15_pt70->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt70 );

  fH15_pt75 = new TProfile ( "histo15_pt75", "<p_{T}> track in dR(jet1) - 75% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_pt75->SetStats ( kTRUE );
  fH15_pt75->GetXaxis()->SetTitle ( "R" );
  fH15_pt75->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt75 );

  fH15_pt80 = new TProfile ( "histo15_pt80", "<p_{T}> track in dR(jet1) - 80% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_pt80->SetStats ( kTRUE );
  fH15_pt80->GetXaxis()->SetTitle ( "R" );
  fH15_pt80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt80 );

  fH15_pt85 = new TProfile ( "histo15_pt85", "<p_{T}> track in dR(jet1) - 85% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_pt85->SetStats ( kTRUE );
  fH15_pt85->GetXaxis()->SetTitle ( "R" );
  fH15_pt85->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt85 );

  fH15_pt90 = new TProfile ( "histo15_pt90", "<p_{T}> track in dR(jet1) - 90% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_pt90->SetStats ( kTRUE );
  fH15_pt90->GetXaxis()->SetTitle ( "R" );
  fH15_pt90->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt90 );



//__________________________________________________________________________________________________
  fH15all = new TProfile ( "histo15all", "<p_{T}> track in dR(jet1)", h15_nbin, h15_xlow, h15_xhigh );
  fH15all->SetStats ( kTRUE );
  fH15all->GetXaxis()->SetTitle ( "R" );
  fH15all->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all->GetXaxis()->SetTitleColor ( 1 );
  fH15all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all );

//=====================================================================================
  fH15all_n70 = new TProfile ( "histo15all_n70", "<p_{T}> track in dR(jet1) - 70% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_n70->SetStats ( kTRUE );
  fH15all_n70->GetXaxis()->SetTitle ( "R" );
  fH15all_n70->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_n70->GetXaxis()->SetTitleColor ( 1 );
  fH15all_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_n70 );

  fH15all_n75 = new TProfile ( "histo15all_n75", "<p_{T}> track in dR(jet1) - 75% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_n75->SetStats ( kTRUE );
  fH15all_n75->GetXaxis()->SetTitle ( "R" );
  fH15all_n75->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_n75->GetXaxis()->SetTitleColor ( 1 );
  fH15all_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_n75 );

  fH15all_n80 = new TProfile ( "histo15all_n80", "<p_{T}> track in dR(jet1) - 80% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_n80->SetStats ( kTRUE );
  fH15all_n80->GetXaxis()->SetTitle ( "R" );
  fH15all_n80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15all_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_n80 );

  fH15all_n85 = new TProfile ( "histo15all_n85", "<p_{T}> track in dR(jet1) - 85% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_n85->SetStats ( kTRUE );
  fH15all_n85->GetXaxis()->SetTitle ( "R" );
  fH15all_n85->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_n85->GetXaxis()->SetTitleColor ( 1 );
  fH15all_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_n85 );

  fH15all_n90 = new TProfile ( "histo15all_n90", "<p_{T}> track in dR(jet1) - 90% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_n90->SetStats ( kTRUE );
  fH15all_n90->GetXaxis()->SetTitle ( "R" );
  fH15all_n90->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_n90->GetXaxis()->SetTitleColor ( 1 );
  fH15all_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_n90 );

//=====================================================================================
  fH15all_pt70 = new TProfile ( "histo15all_pt70", "<p_{T}> track in dR(jet1) - 70% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_pt70->SetStats ( kTRUE );
  fH15all_pt70->GetXaxis()->SetTitle ( "R" );
  fH15all_pt70->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH15all_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_pt70 );

  fH15all_pt75 = new TProfile ( "histo15all_pt75", "<p_{T}> track in dR(jet1) - 75% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_pt75->SetStats ( kTRUE );
  fH15all_pt75->GetXaxis()->SetTitle ( "R" );
  fH15all_pt75->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH15all_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_pt75 );

  fH15all_pt80 = new TProfile ( "histo15all_pt80", "<p_{T}> track in dR(jet1) - 80% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_pt80->SetStats ( kTRUE );
  fH15all_pt80->GetXaxis()->SetTitle ( "R" );
  fH15all_pt80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15all_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_pt80 );

  fH15all_pt85 = new TProfile ( "histo15all_pt85", "<p_{T}> track in dR(jet1) - 85% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_pt85->SetStats ( kTRUE );
  fH15all_pt85->GetXaxis()->SetTitle ( "R" );
  fH15all_pt85->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH15all_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_pt85 );

  fH15all_pt90 = new TProfile ( "histo15all_pt90", "<p_{T}> track in dR(jet1) - 90% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_pt90->SetStats ( kTRUE );
  fH15all_pt90->GetXaxis()->SetTitle ( "R" );
  fH15all_pt90->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15all_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH15all_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_pt90 );

//____________________________________________________________________________________
  fH15_bin = new TH1D ( "histo15_bin", "p_{T} SUM (track) in dR(jet1)", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin->SetStats ( kTRUE );
  fH15_bin->GetXaxis()->SetTitle ( " R" );
  fH15_bin->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin );


//=====================================================================================
  fH15_bin_n70 = new TH1D ( "histo15_bin_n70", "p_{T} SUM (track) in dR(jet1) - 70% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n70->SetStats ( kTRUE );
  fH15_bin_n70->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n70->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_n70->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n70 );

  fH15_bin_n75 = new TH1D ( "histo15_bin_n75", "p_{T} SUM (track) in dR(jet1) - 75% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n75->SetStats ( kTRUE );
  fH15_bin_n75->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n75->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_n75->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n75 );

  fH15_bin_n80 = new TH1D ( "histo15_bin_n80", "p_{T} SUM (track) in dR(jet1) - 80% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n80->SetStats ( kTRUE );
  fH15_bin_n80->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n80->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n80 );

  fH15_bin_n85 = new TH1D ( "histo15_bin_n85", "p_{T} SUM (track) in dR(jet1) - 85% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n85->SetStats ( kTRUE );
  fH15_bin_n85->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n85->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_n85->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n85 );

  fH15_bin_n90 = new TH1D ( "histo15_bin_n90", "p_{T} SUM (track) in dR(jet1) - 90% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n90->SetStats ( kTRUE );
  fH15_bin_n90->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n90->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_n90->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n90 );

//=====================================================================================
  fH15_bin_pt70 = new TH1D ( "histo15_bin_pt70", "p_{T} SUM (track) in dR(jet1) - 70% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt70->SetStats ( kTRUE );
  fH15_bin_pt70->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt70->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt70 );

  fH15_bin_pt75 = new TH1D ( "histo15_bin_pt75", "p_{T} SUM (track) in dR(jet1) - 75% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt75->SetStats ( kTRUE );
  fH15_bin_pt75->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt75->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt75 );

  fH15_bin_pt80 = new TH1D ( "histo15_bin_pt80", "p_{T} SUM (track) in dR(jet1) - 80% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt80->SetStats ( kTRUE );
  fH15_bin_pt80->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt80->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt80 );

  fH15_bin_pt85 = new TH1D ( "histo15_bin_pt85", "p_{T} SUM (track) in dR(jet1) - 85% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt85->SetStats ( kTRUE );
  fH15_bin_pt85->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt85->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt85 );

  fH15_bin_pt90 = new TH1D ( "histo15_bin_pt90", "p_{T} SUM (track) in dR(jet1) - 90% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt90->SetStats ( kTRUE );
  fH15_bin_pt90->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt90->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15_bin_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt90 );


//____________________________________________________________________________________
  fH15all_bin = new TH1D ( "histo15_bin_all", "p_{T} SUM (track) in dR (jet)", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin->SetStats ( kTRUE );
  fH15all_bin->GetXaxis()->SetTitle ( "R" );
  fH15all_bin->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin );

//=====================================================================================
  fH15all_bin_n70 = new TH1D ( "histo15_bin_n70_all", "p_{T} SUM (track) in dR (jet) - 70% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_n70->SetStats ( kTRUE );
  fH15all_bin_n70->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_n70->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_n70->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_n70 );

  fH15all_bin_n75 = new TH1D ( "histo15_bin_n75_all", "p_{T} SUM (track) in dR (jet) - 75% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_n75->SetStats ( kTRUE );
  fH15all_bin_n75->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_n75->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_n75->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_n75 );

  fH15all_bin_n80 = new TH1D ( "histo15_bin_n80_all", "p_{T} SUM (track) in dR (jet) - 80% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_n80->SetStats ( kTRUE );
  fH15all_bin_n80->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_n80->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_n80 );

  fH15all_bin_n85 = new TH1D ( "histo15_bin_n85_all", "p_{T} SUM (track) in dR (jet) - 85% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_n85->SetStats ( kTRUE );
  fH15all_bin_n85->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_n85->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_n85->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_n85 );

  fH15all_bin_n90 = new TH1D ( "histo15_bin_n90_all", "p_{T} SUM (track) in dR (jet) - 90% of particles", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_n90->SetStats ( kTRUE );
  fH15all_bin_n90->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_n90->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_n90->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_n90 );

//=====================================================================================
  fH15all_bin_pt70 = new TH1D ( "histo15_bin_pt70_all", "p_{T} SUM (track) in dR (jet) - 70% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_pt70->SetStats ( kTRUE );
  fH15all_bin_pt70->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_pt70->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_pt70 );

  fH15all_bin_pt75 = new TH1D ( "histo15_bin_pt75_all", "p_{T} SUM (track) in dR (jet) - 75% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_pt75->SetStats ( kTRUE );
  fH15all_bin_pt75->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_pt75->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_pt75 );

  fH15all_bin_pt80 = new TH1D ( "histo15_bin_pt80_all", "p_{T} SUM (track) in dR (jet) - 80% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_pt80->SetStats ( kTRUE );
  fH15all_bin_pt80->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_pt80->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_pt80 );

  fH15all_bin_pt85 = new TH1D ( "histo15_bin_pt85_all", "p_{T} SUM (track) in dR (jet) - 85% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_pt85->SetStats ( kTRUE );
  fH15all_bin_pt85->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_pt85->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_pt85 );

  fH15all_bin_pt90 = new TH1D ( "histo15_bin_pt90_all", "p_{T} SUM (track) in dR (jet) - 90% of Pt", h15_nbin, h15_xlow, h15_xhigh );
  fH15all_bin_pt90->SetStats ( kTRUE );
  fH15all_bin_pt90->GetXaxis()->SetTitle ( "R" );
  fH15all_bin_pt90->GetYaxis()->SetTitle ( "p_{T} SUM (track)" );
  fH15all_bin_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH15all_bin_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15all_bin_pt90 );

//____________________________________________________________________________________
  Int_t h20_nbin = 100;
  Double_t h20_binwidth = 0.01;
  Double_t h20_low = 0.;
  Double_t h20_high = h20_low + h20_binwidth * h20_nbin;
  fH20 = new TH1D ( "histo20", "dN/dR (jet1)", h20_nbin, h20_low, h20_high );
  fH20->SetStats ( kTRUE );
  fH20->GetXaxis()->SetTitle ( "R" );
  fH20->GetYaxis()->SetTitle ( "dN/dR" );
  fH20->GetXaxis()->SetTitleColor ( 1 );
  fH20->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20 );

//=====================================================================================
  fH20_n70 = new TH1D ( "histo20_n70", "dN/dR (jet1) - 70% of particles", h20_nbin, h20_low, h20_high );
  fH20_n70->SetStats ( kTRUE );
  fH20_n70->GetXaxis()->SetTitle ( "R" );
  fH20_n70->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n70->GetXaxis()->SetTitleColor ( 1 );
  fH20_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n70 );

  fH20_n75 = new TH1D ( "histo20_n75", "dN/dR (jet1) - 75% of particles", h20_nbin, h20_low, h20_high );
  fH20_n75->SetStats ( kTRUE );
  fH20_n75->GetXaxis()->SetTitle ( "R" );
  fH20_n75->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n75->GetXaxis()->SetTitleColor ( 1 );
  fH20_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n75 );

  fH20_n80 = new TH1D ( "histo20_n80", "dN/dR (jet1) - 80% of particles", h20_nbin, h20_low, h20_high );
  fH20_n80->SetStats ( kTRUE );
  fH20_n80->GetXaxis()->SetTitle ( "R" );
  fH20_n80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n80->GetXaxis()->SetTitleColor ( 1 );
  fH20_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n80 );

  fH20_n85 = new TH1D ( "histo20_n85", "dN/dR (jet1) - 85% of particles", h20_nbin, h20_low, h20_high );
  fH20_n85->SetStats ( kTRUE );
  fH20_n85->GetXaxis()->SetTitle ( "R" );
  fH20_n85->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n85->GetXaxis()->SetTitleColor ( 1 );
  fH20_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n85 );

  fH20_n90 = new TH1D ( "histo20_n90", "dN/dR (jet1) - 90% of particles", h20_nbin, h20_low, h20_high );
  fH20_n90->SetStats ( kTRUE );
  fH20_n90->GetXaxis()->SetTitle ( "R" );
  fH20_n90->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n90->GetXaxis()->SetTitleColor ( 1 );
  fH20_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n90 );

//=====================================================================================
  fH20_pt70 = new TH1D ( "histo20_pt70", "dN/dR (jet1) - 70% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt70->SetStats ( kTRUE );
  fH20_pt70->GetXaxis()->SetTitle ( "R" );
  fH20_pt70->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt70 );

  fH20_pt75 = new TH1D ( "histo20_pt75", "dN/dR (jet1) - 75% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt75->SetStats ( kTRUE );
  fH20_pt75->GetXaxis()->SetTitle ( "R" );
  fH20_pt75->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt75 );

  fH20_pt80 = new TH1D ( "histo20_pt80", "dN/dR (jet1) - 80% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt80->SetStats ( kTRUE );
  fH20_pt80->GetXaxis()->SetTitle ( "R" );
  fH20_pt80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt80 );

  fH20_pt85 = new TH1D ( "histo20_pt85", "dN/dR (jet1) - 85% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt85->SetStats ( kTRUE );
  fH20_pt85->GetXaxis()->SetTitle ( "R" );
  fH20_pt85->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt85 );

  fH20_pt90 = new TH1D ( "histo20_pt90", "dN/dR (jet1) - 90% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt90->SetStats ( kTRUE );
  fH20_pt90->GetXaxis()->SetTitle ( "R" );
  fH20_pt90->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt90 );


//____________________________________________________________________________________
  fH20all = new TH1D ( "histo20_all", "dN/dR - all jets", h20_nbin, h20_low, h20_high );
  fH20all->SetStats ( kTRUE );
  fH20all->GetXaxis()->SetTitle ( "R" );
  fH20all->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all->GetXaxis()->SetTitleColor ( 1 );
  fH20all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all );


//=====================================================================================
  fH20all_n70 = new TH1D ( "histo20_n70_all", "dN/dR - all jets - 70% of particles", h20_nbin, h20_low, h20_high );
  fH20all_n70->SetStats ( kTRUE );
  fH20all_n70->GetXaxis()->SetTitle ( "R" );
  fH20all_n70->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_n70->GetXaxis()->SetTitleColor ( 1 );
  fH20all_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_n70 );

  fH20all_n75 = new TH1D ( "histo20_n75_all", "dN/dR - all jets - 75% of particles", h20_nbin, h20_low, h20_high );
  fH20all_n75->SetStats ( kTRUE );
  fH20all_n75->GetXaxis()->SetTitle ( "R" );
  fH20all_n75->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_n75->GetXaxis()->SetTitleColor ( 1 );
  fH20all_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_n75 );

  fH20all_n80 = new TH1D ( "histo20_n80_all", "dN/dR - all jets - 80% of particles", h20_nbin, h20_low, h20_high );
  fH20all_n80->SetStats ( kTRUE );
  fH20all_n80->GetXaxis()->SetTitle ( "R" );
  fH20all_n80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_n80->GetXaxis()->SetTitleColor ( 1 );
  fH20all_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_n80 );

  fH20all_n85 = new TH1D ( "histo20_n85_all", "dN/dR - all jets - 85% of particles", h20_nbin, h20_low, h20_high );
  fH20all_n85->SetStats ( kTRUE );
  fH20all_n85->GetXaxis()->SetTitle ( "R" );
  fH20all_n85->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_n85->GetXaxis()->SetTitleColor ( 1 );
  fH20all_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_n85 );

  fH20all_n90 = new TH1D ( "histo20_n90_all", "dN/dR - all jets - 90% of particles", h20_nbin, h20_low, h20_high );
  fH20all_n90->SetStats ( kTRUE );
  fH20all_n90->GetXaxis()->SetTitle ( "R" );
  fH20all_n90->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_n90->GetXaxis()->SetTitleColor ( 1 );
  fH20all_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_n90 );

//=====================================================================================
  fH20all_pt70 = new TH1D ( "histo20_pt70_all", "dN/dR - all jets - 70% of Pt", h20_nbin, h20_low, h20_high );
  fH20all_pt70->SetStats ( kTRUE );
  fH20all_pt70->GetXaxis()->SetTitle ( "R" );
  fH20all_pt70->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_pt70->GetXaxis()->SetTitleColor ( 1 );
  fH20all_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_pt70 );

  fH20all_pt75 = new TH1D ( "histo20_pt75_all", "dN/dR - all jets - 75% of Pt", h20_nbin, h20_low, h20_high );
  fH20all_pt75->SetStats ( kTRUE );
  fH20all_pt75->GetXaxis()->SetTitle ( "R" );
  fH20all_pt75->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_pt75->GetXaxis()->SetTitleColor ( 1 );
  fH20all_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_pt75 );

  fH20all_pt80 = new TH1D ( "histo20_pt80_all", "dN/dR - all jets - 80% of Pt", h20_nbin, h20_low, h20_high );
  fH20all_pt80->SetStats ( kTRUE );
  fH20all_pt80->GetXaxis()->SetTitle ( "R" );
  fH20all_pt80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH20all_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_pt80 );

  fH20all_pt85 = new TH1D ( "histo20_pt85_all", "dN/dR - all jets - 85% of Pt", h20_nbin, h20_low, h20_high );
  fH20all_pt85->SetStats ( kTRUE );
  fH20all_pt85->GetXaxis()->SetTitle ( "R" );
  fH20all_pt85->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_pt85->GetXaxis()->SetTitleColor ( 1 );
  fH20all_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_pt85 );

  fH20all_pt90 = new TH1D ( "histo20_pt90_all", "dN/dR - all jets - 90% of Pt", h20_nbin, h20_low, h20_high );
  fH20all_pt90->SetStats ( kTRUE );
  fH20all_pt90->GetXaxis()->SetTitle ( "R" );
  fH20all_pt90->GetYaxis()->SetTitle ( "dN/dR" );
  fH20all_pt90->GetXaxis()->SetTitleColor ( 1 );
  fH20all_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20all_pt90 );

  //____________________________________________________________________________________
  Int_t h23_nbin = 400;
  Double_t h23_binwidth = 1.;
  Double_t h23_low = 0.;
  Double_t h23_high = h23_low + h23_binwidth * h23_nbin;

  fH23 = new TProfile ( "histo23", "Jet1 Size vs P_{T}(jet1)", h23_nbin, h23_low, h23_high );
  fH23->SetStats ( kTRUE );
  fH23->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH23->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH23->GetXaxis()->SetTitleColor ( 1 );
  fH23->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23 );

  fH24 = new TProfile ( "histo24", "Jet1 Size vs P_{T}(jet1) - 80% of particles", h23_nbin, h23_low, h23_high );
  fH24->SetStats ( kTRUE );
  fH24->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH24->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH24->GetXaxis()->SetTitleColor ( 1 );
  fH24->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH24 );

  fH25 = new TProfile ( "histo25", "Jet1 Size vs P_{T}(jet1) - 80% of Pt", h23_nbin, h23_low, h23_high );
  fH25->SetStats ( kTRUE );
  fH25->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH25->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH25->GetXaxis()->SetTitleColor ( 1 );
  fH25->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH25 );

  fH23all = new TProfile ( "histo23all", "Jet1 Size vs P_{T} - all jets", h23_nbin, h23_low, h23_high );
  fH23all->SetStats ( kTRUE );
  fH23all->GetXaxis()->SetTitle ( "p_{T}(jet) (GeV/c)" );
  fH23all->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH23all->GetXaxis()->SetTitleColor ( 1 );
  fH23all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23all );

  fH24all = new TProfile ( "histo24all", "Jet1 Size vs P_{T}(jet1) - 80% of particles", h23_nbin, h23_low, h23_high );
  fH24all->SetStats ( kTRUE );
  fH24all->GetXaxis()->SetTitle ( "p_{T}(jet) (GeV/c)" );
  fH24all->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH24all->GetXaxis()->SetTitleColor ( 1 );
  fH24all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH24all );

  fH25all = new TProfile ( "histo25all", "Jet1 Size vs P_{T}(jet1) - 80% of Pt", h23_nbin, h23_low, h23_high );
  fH25all->SetStats ( kTRUE );
  fH25all->GetXaxis()->SetTitle ( "p_{T}(jet) (GeV/c)" );
  fH25all->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH25all->GetXaxis()->SetTitleColor ( 1 );
  fH25all->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH25all );

//=============================================================
  Int_t hg_nbin = 100;
  Double_t hg_binwidth = 0.01;
  Double_t hg_low = 0.;
  Double_t hg_high = hg_low + hg_binwidth * hg_nbin;

  fHg = new TH1D ( "histo_g", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg->SetStats ( kTRUE );
  fHg->GetXaxis()->SetTitle ( "g" );
  fHg->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg->GetXaxis()->SetTitleColor ( 1 );
  fHg->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg );

  fHg_n90 = new TH1D ( "histo_g_n90", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_n90->SetStats ( kTRUE );
  fHg_n90->GetXaxis()->SetTitle ( "g" );
  fHg_n90->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_n90->GetXaxis()->SetTitleColor ( 1 );
  fHg_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_n90 );

  fHg_n85 = new TH1D ( "histo_g_n85", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_n85->SetStats ( kTRUE );
  fHg_n85->GetXaxis()->SetTitle ( "g" );
  fHg_n85->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_n85->GetXaxis()->SetTitleColor ( 1 );
  fHg_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_n85 );

  fHg_n80 = new TH1D ( "histo_g_n80", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_n80->SetStats ( kTRUE );
  fHg_n80->GetXaxis()->SetTitle ( "g" );
  fHg_n80->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_n80->GetXaxis()->SetTitleColor ( 1 );
  fHg_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_n80 );

  fHg_n75 = new TH1D ( "histo_g_n75", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_n75->SetStats ( kTRUE );
  fHg_n75->GetXaxis()->SetTitle ( "g" );
  fHg_n75->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_n75->GetXaxis()->SetTitleColor ( 1 );
  fHg_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_n75 );

  fHg_n70 = new TH1D ( "histo_g_n70", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_n70->SetStats ( kTRUE );
  fHg_n70->GetXaxis()->SetTitle ( "g" );
  fHg_n70->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_n70->GetXaxis()->SetTitleColor ( 1 );
  fHg_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_n70 );

  fHg_pt90 = new TH1D ( "histo_g_pt90", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_pt90->SetStats ( kTRUE );
  fHg_pt90->GetXaxis()->SetTitle ( "g" );
  fHg_pt90->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_pt90->GetXaxis()->SetTitleColor ( 1 );
  fHg_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_pt90 );

  fHg_pt85 = new TH1D ( "histo_g_pt85", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_pt85->SetStats ( kTRUE );
  fHg_pt85->GetXaxis()->SetTitle ( "g" );
  fHg_pt85->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_pt85->GetXaxis()->SetTitleColor ( 1 );
  fHg_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_pt85 );

  fHg_pt80 = new TH1D ( "histo_g_pt80", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_pt80->SetStats ( kTRUE );
  fHg_pt80->GetXaxis()->SetTitle ( "g" );
  fHg_pt80->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_pt80->GetXaxis()->SetTitleColor ( 1 );
  fHg_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_pt80 );

  fHg_pt75 = new TH1D ( "histo_g_pt75", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_pt75->SetStats ( kTRUE );
  fHg_pt75->GetXaxis()->SetTitle ( "g" );
  fHg_pt75->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_pt75->GetXaxis()->SetTitleColor ( 1 );
  fHg_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_pt75 );

  fHg_pt70 = new TH1D ( "histo_g_pt70", "dN/dg", hg_nbin, hg_low, hg_high );
  fHg_pt70->SetStats ( kTRUE );
  fHg_pt70->GetXaxis()->SetTitle ( "g" );
  fHg_pt70->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHg_pt70->GetXaxis()->SetTitleColor ( 1 );
  fHg_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHg_pt70 );


//=============================================================
  Int_t hptd_nbin = 100;
  Double_t hptd_binwidth = 0.01;
  Double_t hptd_low = 0.;
  Double_t hptd_high = hptd_low + hptd_binwidth * hptd_nbin;

  fHptd = new TH1D ( "histo_g", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd->SetStats ( kTRUE );
  fHptd->GetXaxis()->SetTitle ( "g" );
  fHptd->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd->GetXaxis()->SetTitleColor ( 1 );
  fHptd->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd );

  fHptd_n90 = new TH1D ( "histo_ptd_n90", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_n90->SetStats ( kTRUE );
  fHptd_n90->GetXaxis()->SetTitle ( "g" );
  fHptd_n90->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_n90->GetXaxis()->SetTitleColor ( 1 );
  fHptd_n90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_n90 );

  fHptd_n85 = new TH1D ( "histo_ptd_n85", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_n85->SetStats ( kTRUE );
  fHptd_n85->GetXaxis()->SetTitle ( "g" );
  fHptd_n85->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_n85->GetXaxis()->SetTitleColor ( 1 );
  fHptd_n85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_n85 );

  fHptd_n80 = new TH1D ( "histo_ptd_n80", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_n80->SetStats ( kTRUE );
  fHptd_n80->GetXaxis()->SetTitle ( "g" );
  fHptd_n80->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_n80->GetXaxis()->SetTitleColor ( 1 );
  fHptd_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_n80 );

  fHptd_n75 = new TH1D ( "histo_ptd_n75", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_n75->SetStats ( kTRUE );
  fHptd_n75->GetXaxis()->SetTitle ( "g" );
  fHptd_n75->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_n75->GetXaxis()->SetTitleColor ( 1 );
  fHptd_n75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_n75 );

  fHptd_n70 = new TH1D ( "histo_ptd_n70", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_n70->SetStats ( kTRUE );
  fHptd_n70->GetXaxis()->SetTitle ( "g" );
  fHptd_n70->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_n70->GetXaxis()->SetTitleColor ( 1 );
  fHptd_n70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_n70 );

  fHptd_pt90 = new TH1D ( "histo_ptd_pt90", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_pt90->SetStats ( kTRUE );
  fHptd_pt90->GetXaxis()->SetTitle ( "g" );
  fHptd_pt90->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_pt90->GetXaxis()->SetTitleColor ( 1 );
  fHptd_pt90->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_pt90 );

  fHptd_pt85 = new TH1D ( "histo_ptd_pt85", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_pt85->SetStats ( kTRUE );
  fHptd_pt85->GetXaxis()->SetTitle ( "g" );
  fHptd_pt85->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_pt85->GetXaxis()->SetTitleColor ( 1 );
  fHptd_pt85->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_pt85 );

  fHptd_pt80 = new TH1D ( "histo_ptd_pt80", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_pt80->SetStats ( kTRUE );
  fHptd_pt80->GetXaxis()->SetTitle ( "g" );
  fHptd_pt80->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_pt80->GetXaxis()->SetTitleColor ( 1 );
  fHptd_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_pt80 );

  fHptd_pt75 = new TH1D ( "histo_ptd_pt75", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_pt75->SetStats ( kTRUE );
  fHptd_pt75->GetXaxis()->SetTitle ( "g" );
  fHptd_pt75->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_pt75->GetXaxis()->SetTitleColor ( 1 );
  fHptd_pt75->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_pt75 );

  fHptd_pt70 = new TH1D ( "histo_ptd_pt70", "dN/dg", hptd_nbin, hptd_low, hptd_high );
  fHptd_pt70->SetStats ( kTRUE );
  fHptd_pt70->GetXaxis()->SetTitle ( "g" );
  fHptd_pt70->GetYaxis()->SetTitle ( "1/N_{jets} dN/dg" );
  fHptd_pt70->GetXaxis()->SetTitleColor ( 1 );
  fHptd_pt70->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fHptd_pt70 );


  // =========== Switch on Sumw2 for all histos ===========
  for ( Int_t i = 0; i < fOutput->GetEntries(); ++i )
      {
      TH1 *h1 = dynamic_cast<TH1 *> ( fOutput->At ( i ) );
      if ( h1 ) { h1->Sumw2();continue; }

      TProfile *hprof1 = dynamic_cast<TProfile *> ( fOutput->At ( i ) );
      if ( hprof1 ) { hprof1->Sumw2(); }
      }

  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  }

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetCDF::DeltaR ( const AliVParticle *part1, const AliVParticle *part2 )
  {
  // Helper function to calculate the distance between two jets or a jet and
  // particle
  Double_t dPhi = part1->Phi() - part2->Phi();
  Double_t dEta = part1->Eta() - part2->Eta();
  dPhi = TVector2::Phi_mpi_pi ( dPhi );

  return TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
  }

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
  Int_t entries = trackscont->GetNAcceptedParticles();

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;
  pair_list.reserve ( entries );

  trackscont->ResetCurrentID();
  for ( Int_t i_entry = 0; i_entry < entries; i_entry++ )
      {
      AliVParticle *track = trackscont->GetNextAcceptParticle ();

      if ( !track )
          {
          AliError ( Form ( "Unable to find accepted track %d in collection %s", i_entry, trackscont->GetName() ) );
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

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::IdxInArray ( Int_t index, TArrayI &array )
  {
  for ( Size_t i = 0; i < array.GetSize(); i++ )
      {
      if ( index == array[i] ) { return kTRUE; }
      }

  return kFALSE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::ExecOnce()
  {
  AliAnalysisTaskEmcalJet::ExecOnce();

  if ( fJetsCont ) { delete fJetsCont ; fJetsCont = 0; }

  if ( fTracksCont ) { delete fTracksCont ; fTracksCont = 0; }

  if ( fCaloClustersCont ) { delete fCaloClustersCont ; fCaloClustersCont = 0; }
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::Terminate ( Option_t * )
  {
  // Called once at the end of the analysis.
  // Update pointers reading them from the output slot
  fOutput = dynamic_cast<TList *> ( GetOutputData (0) );
  }

// kate: indent-mode none; indent-width 2; replace-tabs on;
