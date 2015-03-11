/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// Analysis of FF, multiplicity and pt distribution in leading jets
//
// Author: Adrian Sevcenco
//

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

using namespace std;

ClassImp ( AliAnalysisTaskEmcalJetCDF )

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF()
  : AliAnalysisTaskEmcalJet ( "AliAnalysisTaskEmcalJetCDF", kTRUE ),
    fH1 ( NULL ),
    fH2 ( NULL ),
    fH3 ( NULL ),
    fH4 ( NULL ),
    fH5 ( NULL ),
    fH6 ( NULL ),
    fH7 ( NULL ),
    fH8 ( NULL ),
    fH8xi ( NULL ),
    fH9 ( NULL ),
    fH10 ( NULL ),
    fH9_bin ( NULL ),
    fH10_bin ( NULL ),
    fH9_bin_wojet1 ( NULL ),
    fH10_bin_wojet1 ( NULL ),
    fH15 ( NULL ),
    fH15_n80 ( NULL ),
    fH15_pt80 ( NULL ),
    fH15_bin ( NULL ),
    fH15_bin_n80 ( NULL ),
    fH15_bin_pt80 ( NULL ),
    fH20 ( NULL ),
    fH20_n80 ( NULL ),
    fH20_pt80 ( NULL ),
    fH21 ( NULL ),
    fH21Toward ( NULL ),
    fH21Transverse_min ( NULL ),
    fH21Transverse_max ( NULL ),
    fH21Away ( NULL ),
    fH21_bin ( NULL ),
    fH21Toward_bin ( NULL ),
    fH21Transverse_min_bin ( NULL ),
    fH21Transverse_max_bin ( NULL ),
    fH21Away_bin ( NULL ),
    fH21_bin_wojet1 ( NULL ),
    fH21Toward_bin_wojet1 ( NULL ),
    fH21Transverse_min_bin_wojet1 ( NULL ),
    fH21Transverse_max_bin_wojet1 ( NULL ),
    fH21Away_bin_wojet1 ( NULL ),
    fH22 ( NULL ),
    fH22Toward ( NULL ),
    fH22Transverse_min ( NULL ),
    fH22Transverse_max ( NULL ),
    fH22Away ( NULL ),
    fH22_bin ( NULL ),
    fH22Toward_bin ( NULL ),
    fH22Transverse_min_bin ( NULL ),
    fH22Transverse_max_bin ( NULL ),
    fH22Away_bin ( NULL ),
    fH22_bin_wojet1 ( NULL ),
    fH22Toward_bin_wojet1 ( NULL ),
    fH22Transverse_min_bin_wojet1 ( NULL ),
    fH22Transverse_max_bin_wojet1 ( NULL ),
    fH22Away_bin_wojet1 ( NULL ),
    fH23 ( NULL ),
    fH23jet1 ( NULL ),
    fH23Toward ( NULL ),
    fH23Transverse_min ( NULL ),
    fH23Transverse_max ( NULL ),
    fH23Away ( NULL ),
    fH24 ( NULL ),
    fH25 ( NULL ),
    fH26 ( NULL ),
    fH26_n80 ( NULL ),
    fH26_pt80 ( NULL ),
    fH26jet1 ( NULL ),
    fH26jet1_n80 ( NULL ),
    fH26jet1_pt80 ( NULL ),
    fH26_bin ( NULL ),
    fH26_n80_bin ( NULL ),
    fH26_pt80_bin ( NULL ),
    fH26jet1_bin ( NULL ),
    fH26jet1_n80_bin ( NULL ),
    fH26jet1_pt80_bin ( NULL ),
    fH27 ( NULL ),
    fH27_n80 ( NULL ),
    fH27_pt80 ( NULL ),
    fH27jet1 ( NULL ),
    fH27jet1_n80 ( NULL ),
    fH27jet1_pt80 ( NULL ),
    fH27_bin ( NULL ),
    fH27_n80_bin ( NULL ),
    fH27_pt80_bin ( NULL ),
    fH27jet1_bin ( NULL ),
    fH27jet1_n80_bin ( NULL ),
    fH27jet1_pt80_bin ( NULL ),
    fH40 ( NULL ),
    fH40toward ( NULL ),
    fH40away ( NULL ),
    fH40transmin ( NULL ),
    fH40transmax ( NULL ),
    fH40_bin ( NULL ),
    fH40toward_bin ( NULL ),
    fH40away_bin ( NULL ),
    fH40transmin_bin ( NULL ),
    fH40transmax_bin ( NULL ),
    fH41 ( NULL ),
    fH41toward ( NULL ),
    fH41away ( NULL ),
    fH41transmin ( NULL ),
    fH41transmax ( NULL ),
    fH41_bin ( NULL ),
    fH41toward_bin ( NULL ),
    fH41away_bin ( NULL ),
    fH41transmin_bin ( NULL ),
    fH41transmax_bin ( NULL ),
    fJetsCont ( NULL ),
    fTracksCont ( NULL ),
    fCaloClustersCont ( NULL ),
    fTracksContArray ( NULL ),
    fJet1 ( NULL ),
    fNJets_accepted ( 0 ),
    fNaccPart ( 0 ),
    fEvPt ( -99.99 ),
    fJet1_sorted_idxvec(),
    fEvent_sorted_idxvec()
  {
  // Default constructor.
  fDebug = AliLog::GetGlobalDebugLevel();
  // SetMakeGeneralHistograms ( kTRUE );
  }

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::AliAnalysisTaskEmcalJetCDF ( const char *name )
  : AliAnalysisTaskEmcalJet ( name, kTRUE ),
    fH1 ( NULL ),
    fH2 ( NULL ),
    fH3 ( NULL ),
    fH4 ( NULL ),
    fH5 ( NULL ),
    fH6 ( NULL ),
    fH7 ( NULL ),
    fH8 ( NULL ),
    fH8xi ( NULL ),
    fH9 ( NULL ),
    fH10 ( NULL ),
    fH9_bin ( NULL ),
    fH10_bin ( NULL ),
    fH15 ( NULL ),
    fH15_n80 ( NULL ),
    fH15_pt80 ( NULL ),
    fH15_bin ( NULL ),
    fH15_bin_n80 ( NULL ),
    fH15_bin_pt80 ( NULL ),
    fH20 ( NULL ),
    fH20_n80 ( NULL ),
    fH20_pt80 ( NULL ),
    fH21 ( NULL ),
    fH21Toward ( NULL ),
    fH21Transverse_min ( NULL ),
    fH21Transverse_max ( NULL ),
    fH21Away ( NULL ),
    fH21_bin ( NULL ),
    fH21Toward_bin ( NULL ),
    fH21Transverse_min_bin ( NULL ),
    fH21Transverse_max_bin ( NULL ),
    fH21Away_bin ( NULL ),
    fH22 ( NULL ),
    fH22Toward ( NULL ),
    fH22Transverse_min ( NULL ),
    fH22Transverse_max ( NULL ),
    fH22Away ( NULL ),
    fH22_bin ( NULL ),
    fH22Toward_bin ( NULL ),
    fH22Transverse_min_bin ( NULL ),
    fH22Transverse_max_bin ( NULL ),
    fH22Away_bin ( NULL ),
    fH23 ( NULL ),
    fH23jet1 ( NULL ),
    fH23Toward ( NULL ),
    fH23Transverse_min ( NULL ),
    fH23Transverse_max ( NULL ),
    fH23Away ( NULL ),
    fH24 ( NULL ),
    fH25 ( NULL ),
    fH26 ( NULL ),
    fH26_n80 ( NULL ),
    fH26_pt80 ( NULL ),
    fH26jet1 ( NULL ),
    fH26jet1_n80 ( NULL ),
    fH26jet1_pt80 ( NULL ),
    fH26_bin ( NULL ),
    fH26_n80_bin ( NULL ),
    fH26_pt80_bin ( NULL ),
    fH26jet1_bin ( NULL ),
    fH26jet1_n80_bin ( NULL ),
    fH26jet1_pt80_bin ( NULL ),
    fH27 ( NULL ),
    fH27_n80 ( NULL ),
    fH27_pt80 ( NULL ),
    fH27jet1 ( NULL ),
    fH27jet1_n80 ( NULL ),
    fH27jet1_pt80 ( NULL ),
    fH27_bin ( NULL ),
    fH27_n80_bin ( NULL ),
    fH27_pt80_bin ( NULL ),
    fH27jet1_bin ( NULL ),
    fH27jet1_n80_bin ( NULL ),
    fH27jet1_pt80_bin ( NULL ),
    fH40 ( NULL ),
    fH40toward ( NULL ),
    fH40away ( NULL ),
    fH40transmin ( NULL ),
    fH40transmax ( NULL ),
    fH40_bin ( NULL ),
    fH40toward_bin ( NULL ),
    fH40away_bin ( NULL ),
    fH40transmin_bin ( NULL ),
    fH40transmax_bin ( NULL ),
    fH41 ( NULL ),
    fH41toward ( NULL ),
    fH41away ( NULL ),
    fH41transmin ( NULL ),
    fH41transmax ( NULL ),
    fH41_bin ( NULL ),
    fH41toward_bin ( NULL ),
    fH41away_bin ( NULL ),
    fH41transmin_bin ( NULL ),
    fH41transmax_bin ( NULL ),
    fJetsCont ( NULL ),
    fTracksCont ( NULL ),
    fCaloClustersCont ( NULL ),
    fTracksContArray ( NULL ),
    fJet1 ( NULL ),
    fNJets_accepted ( 0 ),
    fNaccPart ( 0 ),
    fEvPt ( -99.99 ),
    fJet1_sorted_idxvec(),
    fEvent_sorted_idxvec()
  {
  // Standard constructor.
  fDebug = AliLog::GetGlobalDebugLevel();
  // SetMakeGeneralHistograms ( kTRUE );
  }

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDF::~AliAnalysisTaskEmcalJetCDF()
  {
  // Destructor.
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::Run()
  {
  Int_t idx_jet_container = 0;

  fJetsCont = GetJetContainer ( idx_jet_container );

  if ( !fJetsCont )
      {
      std::cout << "ERROR :: Jet Container not found!!!" << std::endl;
      return kFALSE;
      }

  // get particles and clusters connected to jets
  fTracksCont = fJetsCont->GetParticleContainer();
  fTracksCont->SetClassName ( "AliVTrack" );

  fCaloClustersCont = fJetsCont->GetClusterContainer();
  fCaloClustersCont->SetClassName ( "AliVCluster" );

  fNJets_accepted = fJetsCont->GetNJets();          // Number of Jets found in event -
  // accepted cuts applied by
  // JetContainer
  fNaccPart = fTracksCont->GetNAcceptedParticles(); // Multiplicity in event -
  // accepted tracks in tracks
  // container

  // protection
  if ( ( fNJets_accepted < 1 ) || ( fNaccPart < 1 ) )
      {
      if ( fDebug > 2 )
          {
          std::cout << "accepted (fNJets || fNPart) < 1" << std::endl;
          }

      return kFALSE;
      }

  if ( fDebug > 1 )
      {
      std::cout << "fNJets = " << fNJets_accepted << " ; fNPart = " << fNaccPart << std::endl;
      }

  //__________________________________________________________________
  // Leading Jet
  fJet1 = fJetsCont->GetLeadingJet(); // internaly checked for AcceptedJet

  if ( !fJet1 )
      {
      if ( fDebug > 2 )
          {
          std::cout << "LEADING JET NOT FOUND " << std::endl;
          }

      return kFALSE;
      }

  if ( fDebug > 1 )
      {
      std::cout << "+++++++++++++++++>>>>>>>>> Leading jet found" << std::endl;
      fJet1->Print();
      }

  //___________________________________________
  // jet1 : Sorting by p_T jet constituents
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  fTracksContArray = fTracksCont->GetArray();
  fJet1_sorted_idxvec = fJet1->SortConstituentsPt ( fTracksContArray );

  //__________________________________________________________________
  // sorting the EVENT _ACCEPTED_ tracks by pt
  fEvent_sorted_idxvec = SortTracksPt ( fTracksCont );

  //__________________________________________________________________
  // computing of pt sum of event - acccepted tracks
  fEvPt = 0.;

  AliVTrack *track = static_cast<AliVTrack *> ( fTracksCont->GetNextAcceptParticle ( 0 ) );

  while ( track )
      {
      fEvPt += track->Pt(); // pt sum of event
      track = static_cast<AliVTrack *> ( fTracksCont->GetNextAcceptParticle() );
      }

  track = NULL;

  if ( fDebug > 2 )
      {
      std::cout << "Sum of Pt in event pt_sum_event = " << fEvPt << std::endl;
      }

  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().
  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDF::FillHistograms()
  {
  // consts used in analysis
  Double_t const kPI_3 = TMath::Pi() / 3.;
  Double_t const k2PI_3 = 2 * kPI_3;

  AliEmcalJet *jet = fJetsCont->GetNextAcceptJet ( 0 );

  while ( jet )
      {
      fH1->Fill ( jet->Pt() );             // Pt distribution of jets
      fH2->Fill ( jet->Eta() );            // Eta distribution of jets
      fH3->Fill ( jet->Phi() );            // Phi distribution of jets
      fH4->Fill ( jet->GetNumberOfTracks() ); // Multiplicity of jets

      jet = fJetsCont->GetNextAcceptJet();
      }

  jet = NULL;

  fH5->Fill ( fNJets_accepted ); // Distribution of jets in events;

  if ( fDebug > 2 )
      {
      Printf ( "CDFhistos:: end of global jet histos \n" );
      }

  Double_t jet1_pt = fJet1->Pt();
  UInt_t jet1_npart = fJet1->GetNumberOfTracks();

  Int_t jet1_n80 = ( Int_t ) ( 0.8 * jet1_npart );
  Double_t jet1_pt80 = 0.8 * jet1_pt;

  Int_t event_n80 = ( Int_t ) ( 0.8 * fNaccPart );
  Double_t event_pt80 = 0.8 * fEvPt;

  AliVParticle *jet1_trklead = fJet1->GetLeadingTrack ( fTracks );
  Double_t jet1_ptmax = jet1_trklead->Pt();

  UInt_t fNaccPart_woJet1 = fNaccPart - jet1_npart;
  Double_t fEvPt_woJet1 = fEvPt - jet1_pt;

  fH6->Fill ( jet1_npart );       // Jet1 Multiplicity Distribution ( ->Scale (1/events) )
  fH7->Fill ( jet1_pt, jet1_npart ); // N(jet1) vs P_{T}(jet1)

  Int_t counter_part = 0;
  Double_t counter_pt = 0.;
  Int_t track_idx = -999;            // index variable for tracks
  TArrayI jet1_idx_list ( jet1_npart ); // TArrayI of indexes of tracks from leading jet

  //___________________________________________________________________________
  // parsing tracks of jet1 (leading jet) in decreasing order of Pt
  AliVParticle *track = NULL;

  for ( Size_t i = 0; i < jet1_npart; i++ )
      {
      track_idx = fJet1_sorted_idxvec.at ( i );
      track = fJet1->TrackAt ( track_idx, fTracksContArray );

      jet1_idx_list.AddAt ( track_idx, i ); // fill the jet1 track index list with
      // the indexes of tracks

      Double_t dpart = fJet1->DeltaR ( track );
      Double_t track_pt = track->Pt();

      ++counter_part;
      counter_pt += track_pt;

      fH8->Fill ( fJet1->GetZ ( track ) ); // Momentum distribution for leading jet (FF)
      fH8xi->Fill ( fJet1->GetXi ( track ) ); // Momentum distribution for leading jet (FF) xi
      fH23jet1->Fill ( track_pt );      // jet1 pt distribution

      // Recomputing of radius of particles in leading jet
      fH20->Fill ( dpart );            // Distribution of R in leading jet
      fH15->Fill ( dpart, track_pt );  // <p_{T}> track vs the Distance R from Jet1
      fH15_bin->Fill ( dpart, track_pt ); // p_{T} track vs the Distance R from Jet1

      fH26jet1->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
      fH27jet1->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1

      fH26jet1_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
      fH27jet1_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1

      if ( counter_part <= jet1_n80 )        // fill histograms for 80% of particles
          {
          fH15_n80->Fill ( dpart, track_pt );  //  <p_{T}> track vs the Distance R from
          // Jet1 - 80% of particles
          fH15_bin_n80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R
          // from Jet1 - 80% of particles

          fH24->Fill ( jet1_pt, dpart ); //  Jet1 Size vs P_{T}(jet1) - 80% of particles
          fH20_n80->Fill ( dpart );   //  Distribution of R in leading jet

          fH26jet1_n80->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
          fH27jet1_n80->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1

          fH26jet1_n80_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
          fH27jet1_n80_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1
          }

      if ( counter_pt <= jet1_pt80 )          // fill histograms for 80% of pt
          {
          fH15_pt80->Fill ( dpart, track_pt );  //  <p_{T}> track vs the Distance R from
          // Jet1 - 80% of pt
          fH15_bin_pt80->Fill ( dpart, track_pt ); //  p_{T} track vs the Distance R
          // from Jet1 - 80% of pt

          fH25->Fill ( jet1_pt, dpart ); //  Jet1 Size vs P_{T}(jet1) - 80% of Pt
          fH20_pt80->Fill ( dpart );  //  Distribution of R in leading jet

          fH26jet1_pt80->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
          fH27jet1_pt80->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1

          fH26jet1_pt80_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
          fH27jet1_pt80_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1
          }
      }

  track = NULL;

  // reset counter for new usage
  counter_part = 0;
  counter_pt = 0.;
  track_idx = -999;

  // parsing accepted tracks in EVENT in decreasing order of Pt //
  for ( Size_t i = 0; i < fNaccPart; i++ ) // replace the index order by the sorted array
      {
      track_idx = fEvent_sorted_idxvec.at ( i );
      track = fTracksCont->GetNextAcceptParticle ( track_idx );

      if ( !track )
          {
          std::cout << "track not retrieved from fTracksCont" << std::endl;
          continue;
          }

      // pt of the current track
      Double_t track_pt = track->Pt();

      // dphi between leading track (max pt track from leading jet) and current
      // track
      Double_t dphi_part =
        TVector2::Phi_mpi_pi ( track->Phi() - jet1_trklead->Phi() ); // restrict the delta phi to (-pi,pi) interval
      Double_t dphi_part_absdeg = TMath::RadToDeg() * TMath::Abs ( dphi_part );

      // dphi between leading jet and current track
      Double_t dphi_part_jet1 = TVector2::Phi_mpi_pi ( track->Phi() - fJet1->Phi() ); // restrict the delta
      // phi to (-pi,pi)
      // interval
      Double_t dphi_part_jet1_absdeg = TMath::RadToDeg() * TMath::Abs ( dphi_part_jet1 );

      // dR of track from leading jet (jet1)
      Double_t dpart = fJet1->DeltaR ( track );

      ++counter_part;         // next particle
      counter_pt += track_pt; // next particle pt

      fH23->Fill ( track_pt ); //  Pt Distribution of particles in event

      fH9->Fill ( dphi_part_jet1_absdeg, counter_part ); //  N vs the Azimuthal Angle from Jet1
      fH10->Fill ( dphi_part_jet1_absdeg, counter_pt ); //  P_{T} sum vs the Azimuthal Angle from Jet1

      fH9_bin->Fill ( dphi_part_jet1_absdeg );         //  N vs the Azimuthal Angle from Jet1
      fH10_bin->Fill ( dphi_part_jet1_absdeg, track_pt ); //  P_{T} sum vs the Azimuthal Angle from Jet1

      fH26->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
      fH27->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1

      fH26_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
      fH27_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1

      if ( counter_part <= jet1_n80 )        // fill histograms for 80% of particles
          {
          fH26_n80->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
          fH27_n80->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1
          fH26_n80_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
          fH27_n80_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1
          }

      if ( counter_pt <= jet1_pt80 )          // fill histograms for 80% of pt
          {
          fH26_pt80->Fill ( dpart, counter_part ); //  N vs the Distance R from Jet1
          fH27_pt80->Fill ( dpart, counter_pt ); //  PT_{sum} vs the Distance R from Jet1
          fH26_pt80_bin->Fill ( dpart );        //  N vs the Distance R from Jet1
          fH27_pt80_bin->Fill ( dpart, track_pt ); //  PT_{sum} vs the Distance R from Jet1
          }

      // dphi track to jet1 distribution (total and per toward,away,transverse
      // regions)
      fH21->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs
      // P_{T}(jet1)
      fH22->Fill ( jet1_pt, counter_pt ); // PT_{sum}(in the event - including jet1)
      // vs P_{T}(jet1)

      fH21_bin->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
      fH22_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including
      // jet1) vs P_{T}(jet1)

      if ( ( dphi_part_jet1 > ( -1 ) * kPI_3 ) && ( dphi_part_jet1 < kPI_3 ) )
          {
          fH21Toward->Fill ( jet1_pt, counter_part ); // N (in the event - including
          // jet1) vs P_{T}(jet1)
          fH22Toward->Fill ( jet1_pt, counter_pt ); // PT_{sum}(in the event -
          // including jet1) vs P_{T}(jet1)
          fH23Toward->Fill ( track_pt );           // Pt Distribution of particles

          fH21Toward_bin->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
          fH22Toward_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event -
          // including jet1) vs P_{T}(jet1)
          }
      else if ( ( dphi_part_jet1 <= ( -1 ) * kPI_3 ) && ( dphi_part_jet1 > ( -1 ) * k2PI_3 ) )
          {
          fH21Transverse_min->Fill ( jet1_pt, counter_part ); // N (in the event -
          // including jet1) vs
          // P_{T}(jet1)
          fH22Transverse_min->Fill ( jet1_pt, counter_pt ); // PT_{sum}(in the event -
          // including jet1) vs
          // P_{T}(jet1)
          fH23Transverse_min->Fill ( track_pt );           // Pt Distribution of particles

          fH21Transverse_min_bin->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
          fH22Transverse_min_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event
          // - including jet1) vs
          // P_{T}(jet1)
          }
      else if ( ( dphi_part_jet1 >= kPI_3 ) && ( dphi_part_jet1 < k2PI_3 ) )
          {
          fH21Transverse_max->Fill ( jet1_pt, counter_part ); // N (in the event -
          // including jet1) vs
          // P_{T}(jet1)
          fH22Transverse_max->Fill ( jet1_pt, counter_pt ); // PT_{sum}(in the event -
          // including jet1) vs
          // P_{T}(jet1)
          fH23Transverse_max->Fill ( track_pt );           // Pt Distribution of particles

          fH21Transverse_max_bin->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
          fH22Transverse_max_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event
          // - including jet1) vs
          // P_{T}(jet1)
          }
      else if ( ( dphi_part_jet1 >= k2PI_3 ) || ( dphi_part_jet1 <= ( -1 ) * k2PI_3 ) )
          {
          fH21Away->Fill ( jet1_pt, counter_part ); // N (in the event - including
          // jet1) vs P_{T}(jet1)
          fH22Away->Fill ( jet1_pt, counter_pt ); // PT_{sum}(in the event - including
          // jet1) vs P_{T}(jet1)
          fH23Away->Fill ( track_pt );           // Pt Distribution of particles

          fH21Away_bin->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
          fH22Away_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event -
          // including jet1) vs P_{T}(jet1)
          }

      // Condition if track is in leading jet
      if ( !IdxInArray ( track_idx, jet1_idx_list ) )           // fill histos without tracks
          {
          // from leading jet
          fH9_bin_wojet1->Fill ( dphi_part_jet1_absdeg );         //  N vs the Azimuthal Angle from Jet1
          fH10_bin_wojet1->Fill ( dphi_part_jet1_absdeg, track_pt ); //  P_{T} sum vs
          // the Azimuthal
          // Angle from Jet1

          fH21_bin_wojet1->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
          fH22_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event -
          // including jet1) vs
          // P_{T}(jet1)

          if ( ( dphi_part_jet1 > ( -1 ) * kPI_3 ) && ( dphi_part_jet1 < kPI_3 ) )
              {
              fH21Toward_bin_wojet1->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
              fH22Toward_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event
              // - including jet1) vs
              // P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 <= ( -1 ) * kPI_3 ) && ( dphi_part_jet1 > ( -1 ) * k2PI_3 ) )
              {
              fH21Transverse_min_bin_wojet1->Fill ( jet1_pt ); // N (in the event - including jet1) vs P_{T}(jet1)
              fH22Transverse_min_bin_wojet1->Fill ( jet1_pt,
                                                    track_pt ); // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 >= kPI_3 ) && ( dphi_part_jet1 < k2PI_3 ) )
              {
              fH21Transverse_max_bin_wojet1->Fill ( jet1_pt ); // N (in the event - including jet1) vs P_{T}(jet1)
              fH22Transverse_max_bin_wojet1->Fill ( jet1_pt,
                                                    track_pt ); // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 >= k2PI_3 ) || ( dphi_part_jet1 <= ( -1 ) * k2PI_3 ) )
              {
              fH21Away_bin_wojet1->Fill ( jet1_pt );        // N (in the event - including jet1) vs P_{T}(jet1)
              fH22Away_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event -
              // including jet1) vs
              // P_{T}(jet1)
              }
          }

      // NEW UE histos
      // dphi track to leading track distribution (total and per
      // toward,away,transverse regions)
      fH40->Fill ( jet1_ptmax, counter_part ); // total particles fNPart w.r.t PTmax
      // (pt of leading particle from
      // leading jet)
      fH41->Fill ( jet1_ptmax, counter_pt ); // PTsum w.r.t PTmax

      fH40_bin->Fill ( jet1_ptmax );        // total particles fNPart w.r.t PTmax (pt of
      // leading particle from leading jet)
      fH41_bin->Fill ( jet1_ptmax, track_pt ); // PTsum w.r.t PTmax

      if ( ( dphi_part > ( -1 ) * kPI_3 ) && ( dphi_part < kPI_3 ) )
          {
          fH40toward->Fill ( jet1_ptmax, counter_part );
          fH41toward->Fill ( jet1_ptmax, counter_pt );

          fH40toward_bin->Fill ( jet1_ptmax );
          fH41toward_bin->Fill ( jet1_ptmax, track_pt );
          }
      else if ( ( dphi_part <= ( -1 ) * kPI_3 ) && ( dphi_part > ( -1 ) * k2PI_3 ) )
          {
          fH40transmin->Fill ( jet1_ptmax, counter_part );
          fH41transmin->Fill ( jet1_ptmax, counter_pt );

          fH40transmin_bin->Fill ( jet1_ptmax );
          fH41transmin_bin->Fill ( jet1_ptmax, track_pt );
          }
      else if ( ( dphi_part >= kPI_3 ) && ( dphi_part < k2PI_3 ) )
          {
          fH40transmax->Fill ( jet1_ptmax, counter_part );
          fH41transmax->Fill ( jet1_ptmax, counter_pt );

          fH40transmax_bin->Fill ( jet1_ptmax );
          fH41transmax_bin->Fill ( jet1_ptmax, track_pt );
          }
      else if ( ( dphi_part >= k2PI_3 ) || ( dphi_part <= ( -1 ) * k2PI_3 ) )
          {
          fH40away->Fill ( jet1_ptmax, counter_part );
          fH41away->Fill ( jet1_ptmax, counter_pt );

          fH40away_bin->Fill ( jet1_ptmax );
          fH41away_bin->Fill ( jet1_ptmax, track_pt );
          }
      }

  track = NULL;

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
  Int_t h1_nbin = 400;
  Double_t h1_binwidth = 1;
  Double_t h1_low = 0;
  Double_t h1_high = h1_low + h1_binwidth * h1_nbin; // 1GeV/bin
  fH1 = new TH1D ( "histo1", "Pt distribution of jets (accepted)", h1_nbin, h1_low, h1_high );
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

  Int_t h4_nbin = 100;
  Double_t h4_binwidth = 1;
  Double_t h4_low = 0;
  Double_t h4_high = h4_low + h4_binwidth * h4_nbin;
  fH4 = new TH1D ( "histo4", "Multiplicity of jets (accepted)", h4_nbin, h4_low, h4_high ); // 1 unit of multiplicity /bin
  fH4->SetStats ( kTRUE );
  fH4->GetXaxis()->SetTitle ( "N_{tracks}(jet)" );
  fH4->GetYaxis()->SetTitle ( "Jets" );
  fH4->GetXaxis()->SetTitleColor ( 1 );
  fH4->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH4 );

  Int_t h5_nbin = 100;
  Double_t h5_binwidth = 1;
  Double_t h5_low = 0;
  Double_t h5_high = h5_low + h5_binwidth * h5_nbin;
  fH5 = new TH1D ( "histo5", "Distribution of jets in events", h5_nbin, h5_low, h5_high );
  fH5->SetStats ( kTRUE );
  fH5->GetXaxis()->SetTitle ( "N_{jets}" );
  fH5->GetYaxis()->SetTitle ( "N_{events}" );
  fH5->GetXaxis()->SetTitleColor ( 1 );
  fH5->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH5 );

  Int_t h6_nbin = 100;
  Double_t h6_binwidth = 1;
  Double_t h6_low = 0;
  Double_t h6_high = h6_low + h6_binwidth * h6_nbin;
  fH6 = new TH1D ( "histo6", "Jet1 Multiplicity Distribution", h6_nbin, h6_low, h6_high );
  fH6->SetStats ( kTRUE );
  fH6->GetXaxis()->SetTitle ( "N_{tracks}(jet1)" );
  fH6->GetYaxis()->SetTitle ( "1/N(jet1)" );
  fH6->GetXaxis()->SetTitleColor ( 1 );
  fH6->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH6 );

  Int_t h7_nbin = 400;
  Double_t h7_binwidth = 1;
  Double_t h7_xlow = 0;
  Double_t h7_xhigh = h7_xlow + h7_binwidth * h7_nbin;
  fH7 = new TProfile ( "histo7", "N(jet1) vs P_{T}(jet1)", h7_nbin, h7_xlow, h7_xhigh );
  fH7->SetStats ( kTRUE );
  fH7->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH7->GetYaxis()->SetTitle ( "<N_{tracks}(jet1)> in 1 GeV/c bin" );
  fH7->GetXaxis()->SetTitleColor ( 1 );
  fH7->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH7 );

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

  //____________________________________________________________________________________
  Int_t h9_nbin = 50;
  Double_t h9_binwidth = 3.6;
  Double_t h9_xlow = 0.;
  Double_t h9_xhigh = h9_xlow + h9_binwidth * h9_nbin;
  fH9 = new TProfile ( "histo9", "N vs the Azimuthal Angle from Jet1", h9_nbin, h9_xlow, h9_xhigh );
  fH9->SetStats ( kTRUE );
  fH9->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH9->GetYaxis()->SetTitle ( "<N> in 3.6 degree bin" );
  fH9->GetXaxis()->SetTitleColor ( 1 );
  fH9->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH9 );

  fH9_bin = new TH1D ( "histo9_bin", "N vs the Azimuthal Angle from Jet1", h9_nbin, h9_xlow, h9_xhigh );
  fH9_bin->SetStats ( kTRUE );
  fH9_bin->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH9_bin->GetYaxis()->SetTitle ( "<N> in 3.6 degree bin" );
  fH9_bin->GetXaxis()->SetTitleColor ( 1 );
  fH9_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH9_bin );

  fH9_bin_wojet1 = new TH1D ( "histo9_bin_wojet1", "N vs the Azimuthal Angle from Jet1 - wo jet1 components", h9_nbin,
                              h9_xlow, h9_xhigh );
  fH9_bin_wojet1->SetStats ( kTRUE );
  fH9_bin_wojet1->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH9_bin_wojet1->GetYaxis()->SetTitle ( "<N> in 3.6 degree bin" );
  fH9_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH9_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH9_bin_wojet1 );

  //____________________________________________________________________________________
  Int_t h10_nbin = 50;
  Double_t h10_binwidth = 3.6;
  Double_t h10_xlow = 0.;
  Double_t h10_xhigh = h10_xlow + h10_binwidth * h10_nbin;
  fH10 = new TProfile ( "histo10", "P_{T} sum vs the Azimuthal Angle from Jet1", h10_nbin, h10_xlow, h10_xhigh );
  fH10->SetStats ( kTRUE );
  fH10->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH10->GetYaxis()->SetTitle ( "<P_{T} sum> in 3.6 degree bin" );
  fH10->GetXaxis()->SetTitleColor ( 1 );
  fH10->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH10 );

  fH10_bin = new TH1D ( "histo10_bin", "P_{T} sum vs the Azimuthal Angle from Jet1", h10_nbin, h10_xlow, h10_xhigh );
  fH10_bin->SetStats ( kTRUE );
  fH10_bin->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH10_bin->GetYaxis()->SetTitle ( "<P_{T} sum> in 3.6 degree bin" );
  fH10_bin->GetXaxis()->SetTitleColor ( 1 );
  fH10_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH10_bin );

  fH10_bin_wojet1 = new TH1D ( "histo10_bin_wojet1", "P_{T} sum vs the Azimuthal Angle from Jet1 - wo jet1 components",
                               h10_nbin, h10_xlow, h10_xhigh );
  fH10_bin_wojet1->SetStats ( kTRUE );
  fH10_bin_wojet1->GetXaxis()->SetTitle ( "|#Delta#phi| (degrees)" );
  fH10_bin_wojet1->GetYaxis()->SetTitle ( "<P_{T} sum> in 3.6 degree bin" );
  fH10_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH10_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH10_bin_wojet1 );

  //____________________________________________________________________________________
  Int_t h15_nbin = 100;
  Double_t h15_binwidth = 0.01;
  Double_t h15_xlow = 0.;
  Double_t h15_xhigh = h15_xlow + h15_binwidth * h15_nbin;
  fH15 = new TProfile ( "histo15", "<p_{T}> track vs the Distance R from Jet1", h15_nbin, h15_xlow, h15_xhigh );
  fH15->SetStats ( kTRUE );
  fH15->GetXaxis()->SetTitle ( "R" );
  fH15->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15->GetXaxis()->SetTitleColor ( 1 );
  fH15->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15 );

  fH15_n80 = new TProfile ( "histo15_n80", "<p_{T}> track vs the Distance R from Jet1 - 80% of particles", h15_nbin,
                            h15_xlow, h15_xhigh );
  fH15_n80->SetStats ( kTRUE );
  fH15_n80->GetXaxis()->SetTitle ( "R" );
  fH15_n80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_n80 );

  fH15_pt80 = new TProfile ( "histo15_pt80", "<p_{T}> track vs the Distance R from Jet1 - 80% of Pt", h15_nbin, h15_xlow,
                             h15_xhigh );
  fH15_pt80->SetStats ( kTRUE );
  fH15_pt80->GetXaxis()->SetTitle ( "R" );
  fH15_pt80->GetYaxis()->SetTitle ( "<p_{T}> track" );
  fH15_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_pt80 );

  //____________________________________________________________________________________
  fH15_bin = new TProfile ( "histo15_bin", "p_{T}*dR (track) vs the Distance R from Jet1", h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin->SetStats ( kTRUE );
  fH15_bin->GetXaxis()->SetTitle ( "R" );
  fH15_bin->GetYaxis()->SetTitle ( "p_{T}*dR (track)" );
  fH15_bin->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin );

  fH15_bin_n80 = new TProfile ( "histo15_bin_n80", "p_{T}*dR (track) vs the Distance R from Jet1 - 80% of particles",
                                h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_n80->SetStats ( kTRUE );
  fH15_bin_n80->GetXaxis()->SetTitle ( "R" );
  fH15_bin_n80->GetYaxis()->SetTitle ( "p_{T}*dR (track)" );
  fH15_bin_n80->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_n80 );

  fH15_bin_pt80 = new TProfile ( "histo15_bin_pt80", "pp_{T}*dR (track) vs the Distance R from Jet1 - 80% of Pt",
                                 h15_nbin, h15_xlow, h15_xhigh );
  fH15_bin_pt80->SetStats ( kTRUE );
  fH15_bin_pt80->GetXaxis()->SetTitle ( "R" );
  fH15_bin_pt80->GetYaxis()->SetTitle ( "p_{T}*dR (track)" );
  fH15_bin_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH15_bin_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH15_bin_pt80 );

  //____________________________________________________________________________________
  Int_t h20_nbin = 100;
  Double_t h20_binwidth = 0.01;
  Double_t h20_low = 0.;
  Double_t h20_high = h20_low + h20_binwidth * h20_nbin;
  fH20 = new TH1D ( "histo20", "Distribution of R in leading jet", h20_nbin, h20_low, h20_high );
  fH20->SetStats ( kTRUE );
  fH20->GetXaxis()->SetTitle ( "R" );
  fH20->GetYaxis()->SetTitle ( "dN/dR" );
  fH20->GetXaxis()->SetTitleColor ( 1 );
  fH20->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20 );

  fH20_n80 =
    new TH1D ( "histo20_n80", "Distribution of R in leading jet - 80% of particles", h20_nbin, h20_low, h20_high );
  fH20_n80->SetStats ( kTRUE );
  fH20_n80->GetXaxis()->SetTitle ( "R" );
  fH20_n80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_n80->GetXaxis()->SetTitleColor ( 1 );
  fH20_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_n80 );

  fH20_pt80 = new TH1D ( "histo20_pt80", "Distribution of R in leading jet - 80% of Pt", h20_nbin, h20_low, h20_high );
  fH20_pt80->SetStats ( kTRUE );
  fH20_pt80->GetXaxis()->SetTitle ( "R" );
  fH20_pt80->GetYaxis()->SetTitle ( "dN/dR" );
  fH20_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH20_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH20_pt80 );

  //____________________________________________________________________________________
  fH21 = new TProfile ( "histo21", "N(in the event - including jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH21->SetStats ( kTRUE );
  fH21->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21->GetXaxis()->SetTitleColor ( 1 );
  fH21->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21 );

  fH21Toward =
    new TProfile ( "histo21_toward", "N(in the event - including jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH21Toward->SetStats ( kTRUE );
  fH21Toward->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Toward->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Toward->GetXaxis()->SetTitleColor ( 1 );
  fH21Toward->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Toward );

  fH21Transverse_min = new TProfile ( "histo21_transverse_min",
                                      "N(in the event - including jet1) vs P_{T}(jet1) - transverse MIN", 200, 0., 200. );
  fH21Transverse_min->SetStats ( kTRUE );
  fH21Transverse_min->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_min->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Transverse_min->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_min->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_min );

  fH21Transverse_max = new TProfile ( "histo21_transverse_max",
                                      "N(in the event - including jet1) vs P_{T}(jet1) - transverse MAX", 200, 0., 200. );
  fH21Transverse_max->SetStats ( kTRUE );
  fH21Transverse_max->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_max->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Transverse_max->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_max->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_max );

  fH21Away = new TProfile ( "histo21_away", "N(in the event - including jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH21Away->SetStats ( kTRUE );
  fH21Away->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Away->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Away->GetXaxis()->SetTitleColor ( 1 );
  fH21Away->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Away );

  //____________________________________________________________________________________
  fH21_bin = new TH1D ( "histo21_bin", "N(in the event - including jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH21_bin->SetStats ( kTRUE );
  fH21_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21_bin->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21_bin->GetXaxis()->SetTitleColor ( 1 );
  fH21_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21_bin );

  fH21Toward_bin =
    new TH1D ( "histo21_toward_bin", "N(in the event - including jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH21Toward_bin->SetStats ( kTRUE );
  fH21Toward_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Toward_bin->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Toward_bin->GetXaxis()->SetTitleColor ( 1 );
  fH21Toward_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Toward_bin );

  fH21Transverse_min_bin = new TH1D ( "histo21_transverse_min_bin",
                                      "N(in the event - including jet1) vs P_{T}(jet1) - transverse MIN", 200, 0., 200. );
  fH21Transverse_min_bin->SetStats ( kTRUE );
  fH21Transverse_min_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_min_bin->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Transverse_min_bin->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_min_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_min_bin );

  fH21Transverse_max_bin = new TH1D ( "histo21_transverse_max_bin",
                                      "N(in the event - including jet1) vs P_{T}(jet1) - transverse MAX", 200, 0., 200. );
  fH21Transverse_max_bin->SetStats ( kTRUE );
  fH21Transverse_max_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_max_bin->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Transverse_max_bin->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_max_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_max_bin );

  fH21Away_bin = new TH1D ( "histo21_away_bin", "N(in the event - including jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH21Away_bin->SetStats ( kTRUE );
  fH21Away_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Away_bin->GetYaxis()->SetTitle ( "<N(in the event - including jet1)> in 1 GeV/c bin" );
  fH21Away_bin->GetXaxis()->SetTitleColor ( 1 );
  fH21Away_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Away_bin );

  //____________________________________________________________________________________
  fH21_bin_wojet1 = new TH1D ( "histo21_bin_wojet1", "N(in the event - without jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH21_bin_wojet1->SetStats ( kTRUE );
  fH21_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21_bin_wojet1->GetYaxis()->SetTitle ( "<N(in the event - without jet1)> in 1 GeV/c bin" );
  fH21_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH21_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21_bin_wojet1 );

  fH21Toward_bin_wojet1 =
    new TH1D ( "histo21_toward_bin_wojet1", "N(in the event - without jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH21Toward_bin_wojet1->SetStats ( kTRUE );
  fH21Toward_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Toward_bin_wojet1->GetYaxis()->SetTitle ( "<N(in the event - without jet1)> in 1 GeV/c bin" );
  fH21Toward_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH21Toward_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Toward_bin_wojet1 );

  fH21Transverse_min_bin_wojet1 =
    new TH1D ( "histo21_transverse_min_bin_wojet1", "N(in the event - without jet1) vs P_{T}(jet1) - transverse MIN", 200,
               0., 200. );
  fH21Transverse_min_bin_wojet1->SetStats ( kTRUE );
  fH21Transverse_min_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_min_bin_wojet1->GetYaxis()->SetTitle ( "<N(in the event - without jet1)> in 1 GeV/c bin" );
  fH21Transverse_min_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_min_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_min_bin_wojet1 );

  fH21Transverse_max_bin_wojet1 =
    new TH1D ( "histo21_transverse_max_bin_wojet1", "N(in the event - without jet1) vs P_{T}(jet1) - transverse MAX", 200,
               0., 200. );
  fH21Transverse_max_bin_wojet1->SetStats ( kTRUE );
  fH21Transverse_max_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Transverse_max_bin_wojet1->GetYaxis()->SetTitle ( "<N(in the event - without jet1)> in 1 GeV/c bin" );
  fH21Transverse_max_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH21Transverse_max_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Transverse_max_bin_wojet1 );

  fH21Away_bin_wojet1 =
    new TH1D ( "histo21_away_bin_wojet1", "N(in the event - without jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH21Away_bin_wojet1->SetStats ( kTRUE );
  fH21Away_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH21Away_bin_wojet1->GetYaxis()->SetTitle ( "<N(in the event - without jet1)> in 1 GeV/c bin" );
  fH21Away_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH21Away_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH21Away_bin_wojet1 );

  //__________________________________________________________________
  fH22 = new TProfile ( "histo22", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH22->SetStats ( kTRUE );
  fH22->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22->GetXaxis()->SetTitleColor ( 1 );
  fH22->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22 );

  fH22Toward =
    new TProfile ( "histo22_toward", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH22Toward->SetStats ( kTRUE );
  fH22Toward->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Toward->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Toward->GetXaxis()->SetTitleColor ( 1 );
  fH22Toward->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Toward );

  fH22Transverse_min = new TProfile (
    "histo22_transverse_min", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - transverse MIN", 200, 0., 200. );
  fH22Transverse_min->SetStats ( kTRUE );
  fH22Transverse_min->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_min->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Transverse_min->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_min->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_min );

  fH22Transverse_max = new TProfile (
    "histo22_transverse_max", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - transverse MAX", 200, 0., 200. );
  fH22Transverse_max->SetStats ( kTRUE );
  fH22Transverse_max->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_max->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Transverse_max->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_max->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_max );

  fH22Away =
    new TProfile ( "histo22_away", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH22Away->SetStats ( kTRUE );
  fH22Away->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Away->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Away->GetXaxis()->SetTitleColor ( 1 );
  fH22Away->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Away );

  //__________________________________________________________________
  fH22_bin = new TH1D ( "histo22_bin", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH22_bin->SetStats ( kTRUE );
  fH22_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22_bin->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22_bin->GetXaxis()->SetTitleColor ( 1 );
  fH22_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22_bin );

  fH22Toward_bin =
    new TH1D ( "histo22_toward_bin", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH22Toward_bin->SetStats ( kTRUE );
  fH22Toward_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Toward_bin->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Toward_bin->GetXaxis()->SetTitleColor ( 1 );
  fH22Toward_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Toward_bin );

  fH22Transverse_min_bin =
    new TH1D ( "histo22_transverse_min_bin", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - transverse MIN",
               200, 0., 200. );
  fH22Transverse_min_bin->SetStats ( kTRUE );
  fH22Transverse_min_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_min_bin->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Transverse_min_bin->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_min_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_min_bin );

  fH22Transverse_max_bin =
    new TH1D ( "histo22_transverse_max_bin", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - transverse MAX",
               200, 0., 200. );
  fH22Transverse_max_bin->SetStats ( kTRUE );
  fH22Transverse_max_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_max_bin->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Transverse_max_bin->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_max_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_max_bin );

  fH22Away_bin =
    new TH1D ( "histo22_away_bin", "PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH22Away_bin->SetStats ( kTRUE );
  fH22Away_bin->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Away_bin->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin" );
  fH22Away_bin->GetXaxis()->SetTitleColor ( 1 );
  fH22Away_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Away_bin );

  //__________________________________________________________________
  fH22_bin_wojet1 =
    new TH1D ( "histo22_bin_wojet1", "PT_{sum}(in the event - without jet1) vs P_{T}(jet1)", 200, 0., 200. );
  fH22_bin_wojet1->SetStats ( kTRUE );
  fH22_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22_bin_wojet1->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - without jet1)> in 1 GeV/c bin" );
  fH22_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH22_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22_bin_wojet1 );

  fH22Toward_bin_wojet1 = new TH1D ( "histo22_toward_bin_wojet1",
                                     "PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - toward", 200, 0., 200. );
  fH22Toward_bin_wojet1->SetStats ( kTRUE );
  fH22Toward_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Toward_bin_wojet1->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - without jet1)> in 1 GeV/c bin" );
  fH22Toward_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH22Toward_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Toward_bin_wojet1 );

  fH22Transverse_min_bin_wojet1 =
    new TH1D ( "histo22_transverse_min_bin_wojet1",
               "PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - transverse MIN", 200, 0., 200. );
  fH22Transverse_min_bin_wojet1->SetStats ( kTRUE );
  fH22Transverse_min_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_min_bin_wojet1->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - without jet1)> in 1 GeV/c bin" );
  fH22Transverse_min_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_min_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_min_bin_wojet1 );

  fH22Transverse_max_bin_wojet1 =
    new TH1D ( "histo22_transverse_max_bin_wojet1",
               "PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - transverse MAX", 200, 0., 200. );
  fH22Transverse_max_bin_wojet1->SetStats ( kTRUE );
  fH22Transverse_max_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Transverse_max_bin_wojet1->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - without jet1)> in 1 GeV/c bin" );
  fH22Transverse_max_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH22Transverse_max_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Transverse_max_bin_wojet1 );

  fH22Away_bin_wojet1 =
    new TH1D ( "histo22_away_bin_wojet1", "PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - away", 200, 0., 200. );
  fH22Away_bin_wojet1->SetStats ( kTRUE );
  fH22Away_bin_wojet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH22Away_bin_wojet1->GetYaxis()->SetTitle ( "<PT_{sum}(in the event - without jet1)> in 1 GeV/c bin" );
  fH22Away_bin_wojet1->GetXaxis()->SetTitleColor ( 1 );
  fH22Away_bin_wojet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH22Away_bin_wojet1 );

  //____________________________________________________________________________________
  fH23 = new TH1D ( "histo23", "Pt Distribution of particles", 1500, 0., 1500. );
  fH23->SetStats ( kTRUE );
  fH23->GetXaxis()->SetTitle ( "p_{T} (GeV/c)" );
  fH23->GetYaxis()->SetTitle ( "dN/dP_{T} (1/GeV/c)" );
  fH23->GetXaxis()->SetTitleColor ( 1 );
  fH23->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23 );

  fH23jet1 = new TH1D ( "histo23jet1", "Pt Distribution of particles in jet1", 1000, 0., 1000. );
  fH23jet1->SetStats ( kTRUE );
  fH23jet1->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH23jet1->GetYaxis()->SetTitle ( "dN/dP_{T}(jet1) (1/GeV/c)" );
  fH23jet1->GetXaxis()->SetTitleColor ( 1 );
  fH23jet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23jet1 );

  fH23Toward = new TH1D ( "histo23_toward", "'Toward' Pt Distribution of particles", 1500, 0., 1500. );
  fH23Toward->SetStats ( kTRUE );
  fH23Toward->GetXaxis()->SetTitle ( "p_{T} (GeV/c)" );
  fH23Toward->GetYaxis()->SetTitle ( "dN/dP_{T} (1/GeV/c)" );
  fH23Toward->GetXaxis()->SetTitleColor ( 1 );
  fH23Toward->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23Toward );

  fH23Transverse_min =
    new TH1D ( "histo23_transverse_min", "Pt Distribution of particles - Transverse MIN", 1500, 0., 1500. );
  fH23Transverse_min->SetStats ( kTRUE );
  fH23Transverse_min->GetXaxis()->SetTitle ( "p_{T} (GeV/c)" );
  fH23Transverse_min->GetYaxis()->SetTitle ( "dN/dP_{T} (1/GeV/c)" );
  fH23Transverse_min->GetXaxis()->SetTitleColor ( 1 );
  fH23Transverse_min->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23Transverse_min );

  fH23Transverse_max =
    new TH1D ( "histo23_transverse_max", "Pt Distribution of particles  - Transverse MAX", 1500, 0., 1500. );
  fH23Transverse_max->SetStats ( kTRUE );
  fH23Transverse_max->GetXaxis()->SetTitle ( "p_{T} (GeV/c)" );
  fH23Transverse_max->GetYaxis()->SetTitle ( "dN/dP_{T} (1/GeV/c)" );
  fH23Transverse_max->GetXaxis()->SetTitleColor ( 1 );
  fH23Transverse_max->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23Transverse_max );

  fH23Away = new TH1D ( "histo23_away", "'Away' Pt Distribution of particles", 1500, 0., 1500. );
  fH23Away->SetStats ( kTRUE );
  fH23Away->GetXaxis()->SetTitle ( "p_{T} (GeV/c)" );
  fH23Away->GetYaxis()->SetTitle ( "dN/dP_{T} (1/GeV/c)" );
  fH23Away->GetXaxis()->SetTitleColor ( 1 );
  fH23Away->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH23Away );

  //____________________________________________________________________________________
  Int_t h24_nbin = 400;
  Double_t h24_binwidth = 1.;
  Double_t h24_low = 0.;
  Double_t h24_high = h24_low + h24_binwidth * h24_nbin;
  fH24 = new TProfile ( "histo24", "Jet1 Size vs P_{T}(jet1) - 80% of particles", h24_nbin, h24_low, h24_high );
  fH24->SetStats ( kTRUE );
  fH24->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH24->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH24->GetXaxis()->SetTitleColor ( 1 );
  fH24->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH24 );

  Int_t h25_nbin = 400;
  Double_t h25_binwidth = 1.;
  Double_t h25_low = 0.;
  Double_t h25_high = h25_low + h25_binwidth * h25_nbin;
  fH25 = new TProfile ( "histo25", "Jet1 Size vs P_{T}(jet1) - 80% of Pt", h25_nbin, h25_low, h25_high );
  fH25->SetStats ( kTRUE );
  fH25->GetXaxis()->SetTitle ( "p_{T}(jet1) (GeV/c)" );
  fH25->GetYaxis()->SetTitle ( "<R(jet1)> in 1 GeV/c bin" );
  fH25->GetXaxis()->SetTitleColor ( 1 );
  fH25->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH25 );

  //____________________________________________________________________________________
  Int_t h26_nbin = 60;
  Double_t h26_binwidth = 0.02;
  Double_t h26_low = 0.;
  Double_t h26_high = h26_low + h26_binwidth * h26_nbin;
  fH26 = new TProfile ( "histo26", "<N> vs the Distance R from Jet1", h26_nbin, h26_low, h26_high );
  fH26->SetStats ( kTRUE );
  fH26->GetXaxis()->SetTitle ( "Distance R" );
  fH26->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26->GetXaxis()->SetTitleColor ( 1 );
  fH26->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26 );

  fH26_n80 =
    new TProfile ( "histo26_n80", "<N> vs the Distance R from Jet1 - 80% of particles", h26_nbin, h26_low, h26_high );
  fH26_n80->SetStats ( kTRUE );
  fH26_n80->GetXaxis()->SetTitle ( "Distance R" );
  fH26_n80->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26_n80->GetXaxis()->SetTitleColor ( 1 );
  fH26_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26_n80 );

  fH26_pt80 = new TProfile ( "histo26_pt80", "<N> vs the Distance R from Jet1 - 80% of Pt", h26_nbin, h26_low, h26_high );
  fH26_pt80->SetStats ( kTRUE );
  fH26_pt80->GetXaxis()->SetTitle ( "Distance R" );
  fH26_pt80->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH26_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26_pt80 );

  fH26jet1 = new TProfile ( "histo26jet1", "<N>(jet1) vs the Distance R from Jet1", h26_nbin, h26_low, h26_high );
  fH26jet1->SetStats ( kTRUE );
  fH26jet1->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1 );

  fH26jet1_n80 = new TProfile ( "histo26jet1_n80", "<N>(jet1) vs the Distance R from Jet1 - 80% of particles", h26_nbin,
                                h26_low, h26_high );
  fH26jet1_n80->SetStats ( kTRUE );
  fH26jet1_n80->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1_n80->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1_n80->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1_n80 );

  fH26jet1_pt80 =
    new TProfile ( "histo26jet1_pt80", "<N>(jet1) vs the Distance R from Jet1 - 80% of Pt", h26_nbin, h26_low, h26_high );
  fH26jet1_pt80->SetStats ( kTRUE );
  fH26jet1_pt80->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1_pt80->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1_pt80 );

  //____________________________________________________________________________________
  fH26_bin = new TH1D ( "histo26_bin", "<N> vs the Distance R from Jet1", h26_nbin, h26_low, h26_high );
  fH26_bin->SetStats ( kTRUE );
  fH26_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26_bin );

  fH26_n80_bin =
    new TH1D ( "histo26_n80_bin", "<N> vs the Distance R from Jet1 - 80% of particles", h26_nbin, h26_low, h26_high );
  fH26_n80_bin->SetStats ( kTRUE );
  fH26_n80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26_n80_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26_n80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26_n80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26_n80_bin );

  fH26_pt80_bin =
    new TH1D ( "histo26_pt80_bin", "<N> vs the Distance R from Jet1 - 80% of Pt", h26_nbin, h26_low, h26_high );
  fH26_pt80_bin->SetStats ( kTRUE );
  fH26_pt80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26_pt80_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26_pt80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26_pt80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26_pt80_bin );

  fH26jet1_bin = new TH1D ( "histo26jet1_bin", "<N>(jet1) vs the Distance R from Jet1", h26_nbin, h26_low, h26_high );
  fH26jet1_bin->SetStats ( kTRUE );
  fH26jet1_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1_bin );

  fH26jet1_n80_bin = new TH1D ( "histo26jet1_n80_bin", "<N>(jet1) vs the Distance R from Jet1 - 80% of particles",
                                h26_nbin, h26_low, h26_high );
  fH26jet1_n80_bin->SetStats ( kTRUE );
  fH26jet1_n80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1_n80_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1_n80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1_n80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1_n80_bin );

  fH26jet1_pt80_bin =
    new TH1D ( "histo26jet1_pt80_bin", "<N>(jet1) vs the Distance R from Jet1 - 80% of Pt", h26_nbin, h26_low, h26_high );
  fH26jet1_pt80_bin->SetStats ( kTRUE );
  fH26jet1_pt80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH26jet1_pt80_bin->GetYaxis()->SetTitle ( "<N> in 0.02 bin" );
  fH26jet1_pt80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH26jet1_pt80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH26jet1_pt80_bin );

  //____________________________________________________________________________________
  Int_t h27_nbin = 60;
  Double_t h27_binwidth = 0.02;
  Double_t h27_low = 0.;
  Double_t h27_high = h27_low + h27_binwidth * h27_nbin;
  fH27 = new TProfile ( "histo27", "PT_{sum} vs the Distance R from Jet1", h27_nbin, h27_low, h27_high );
  fH27->SetStats ( kTRUE );
  fH27->GetXaxis()->SetTitle ( "Distance R" );
  fH27->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27->GetXaxis()->SetTitleColor ( 1 );
  fH27->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27 );

  fH27_n80 =
    new TProfile ( "histo27_n80", "PT_{sum} vs the Distance R from Jet1 - 80% of particles", h27_nbin, h27_low, h27_high );
  fH27_n80->SetStats ( kTRUE );
  fH27_n80->GetXaxis()->SetTitle ( "Distance R" );
  fH27_n80->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27_n80->GetXaxis()->SetTitleColor ( 1 );
  fH27_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27_n80 );

  fH27_pt80 =
    new TProfile ( "histo27_pt80", "PT_{sum} vs the Distance R from Jet1 - 80% of Pt", h27_nbin, h27_low, h27_high );
  fH27_pt80->SetStats ( kTRUE );
  fH27_pt80->GetXaxis()->SetTitle ( "Distance R" );
  fH27_pt80->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH27_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27_pt80 );

  fH27jet1 = new TProfile ( "histo27jet1", "Jet1 p_{T} vs the Distance R from Jet1", h27_nbin, h27_low, h27_high );
  fH27jet1->SetStats ( kTRUE );
  fH27jet1->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1 );

  fH27jet1_n80 = new TProfile ( "histo27jet1_n80", "Jet1 p_{T} vs the Distance R from Jet1 - 80% of particles", h27_nbin,
                                h27_low, h27_high );
  fH27jet1_n80->SetStats ( kTRUE );
  fH27jet1_n80->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1_n80->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1_n80->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1_n80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1_n80 );

  fH27jet1_pt80 =
    new TProfile ( "histo27jet1_pt80", "Jet1 p_{T} vs the Distance R from Jet1 - 80% of Pt", h27_nbin, h27_low, h27_high );
  fH27jet1_pt80->SetStats ( kTRUE );
  fH27jet1_pt80->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1_pt80->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1_pt80->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1_pt80->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1_pt80 );

  //____________________________________________________________________________________
  fH27_bin = new TH1D ( "histo27_bin", "PT_{sum} vs the Distance R from Jet1", h27_nbin, h27_low, h27_high );
  fH27_bin->SetStats ( kTRUE );
  fH27_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27_bin );

  fH27_n80_bin =
    new TH1D ( "histo27_n80_bin", "PT_{sum} vs the Distance R from Jet1 - 80% of particles", h27_nbin, h27_low, h27_high );
  fH27_n80_bin->SetStats ( kTRUE );
  fH27_n80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27_n80_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27_n80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27_n80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27_n80_bin );

  fH27_pt80_bin =
    new TH1D ( "histo27_pt80_bin", "PT_{sum} vs the Distance R from Jet1 - 80% of Pt", h27_nbin, h27_low, h27_high );
  fH27_pt80_bin->SetStats ( kTRUE );
  fH27_pt80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27_pt80_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27_pt80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27_pt80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27_pt80_bin );

  fH27jet1_bin = new TH1D ( "histo27jet1_bin", "Jet1 p_{T} vs the Distance R from Jet1", h27_nbin, h27_low, h27_high );
  fH27jet1_bin->SetStats ( kTRUE );
  fH27jet1_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1_bin );

  fH27jet1_n80_bin = new TH1D ( "histo27jet1_n80_bin", "Jet1 p_{T} vs the Distance R from Jet1 - 80% of particles",
                                h27_nbin, h27_low, h27_high );
  fH27jet1_n80_bin->SetStats ( kTRUE );
  fH27jet1_n80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1_n80_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1_n80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1_n80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1_n80_bin );

  fH27jet1_pt80_bin =
    new TH1D ( "histo27jet1_pt80_bin", "Jet1 p_{T} vs the Distance R from Jet1 - 80% of Pt", h27_nbin, h27_low, h27_high );
  fH27jet1_pt80_bin->SetStats ( kTRUE );
  fH27jet1_pt80_bin->GetXaxis()->SetTitle ( "Distance R" );
  fH27jet1_pt80_bin->GetYaxis()->SetTitle ( "<PT_{sum}> (GeV/c) in 0.02 bin" );
  fH27jet1_pt80_bin->GetXaxis()->SetTitleColor ( 1 );
  fH27jet1_pt80_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH27jet1_pt80_bin );

  //____________________________________________________________________________________
  Int_t h40_nbin = 200;
  Double_t h40_binwidth = 1.;
  Double_t h40_low = 0.;
  Double_t h40_high = h40_low + h40_binwidth * h40_nbin;
  fH40 = new TProfile ( "histo40", "total particles fNPart w.r.t PTmax (pt of leading particle from jet1)", h40_nbin,
                        h40_low, h40_high );
  fH40->SetStats ( kTRUE );
  fH40->GetXaxis()->SetTitle ( "PTmax" );
  fH40->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40->GetXaxis()->SetTitleColor ( 1 );
  fH40->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40 );

  fH40toward = new TProfile ( "histo40toward", "total particles fNPart w.r.t PTmax (pt of "
                              "leading particle from jet1) - TOWARD",
                              h40_nbin, h40_low, h40_high );
  fH40toward->SetStats ( kTRUE );
  fH40toward->GetXaxis()->SetTitle ( "PTmax" );
  fH40toward->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40toward->GetXaxis()->SetTitleColor ( 1 );
  fH40toward->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40toward );

  fH40away = new TProfile ( "histo40away", "total particles fNPart w.r.t PTmax (pt of "
                            "leading particle from jet1) - AWAY",
                            h40_nbin, h40_low, h40_high );
  fH40away->SetStats ( kTRUE );
  fH40away->GetXaxis()->SetTitle ( "PTmax" );
  fH40away->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40away->GetXaxis()->SetTitleColor ( 1 );
  fH40away->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40away );

  fH40transmin = new TProfile ( "histo40transmin", "total particles fNPart w.r.t PTmax (pt of "
                                "leading particle from jet1) - TRANSMIN",
                                h40_nbin, h40_low, h40_high );
  fH40transmin->SetStats ( kTRUE );
  fH40transmin->GetXaxis()->SetTitle ( "PTmax" );
  fH40transmin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40transmin->GetXaxis()->SetTitleColor ( 1 );
  fH40transmin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40transmin );

  fH40transmax = new TProfile ( "histo40transmax", "total particles fNPart w.r.t PTmax (pt of "
                                "leading particle from jet1) - TRANSMAX",
                                h40_nbin, h40_low, h40_high );
  fH40transmax->SetStats ( kTRUE );
  fH40transmax->GetXaxis()->SetTitle ( "PTmax" );
  fH40transmax->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40transmax->GetXaxis()->SetTitleColor ( 1 );
  fH40transmax->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40transmax );

  //____________________________________________________________________________________
  fH40_bin = new TH1D ( "histo40_bin", "total particles fNPart w.r.t PTmax (pt of leading particle from jet1)", h40_nbin,
                        h40_low, h40_high );
  fH40_bin->SetStats ( kTRUE );
  fH40_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH40_bin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40_bin->GetXaxis()->SetTitleColor ( 1 );
  fH40_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40_bin );

  fH40toward_bin = new TH1D ( "histo40toward_bin", "total particles fNPart w.r.t PTmax (pt of "
                              "leading particle from jet1) - TOWARD",
                              h40_nbin, h40_low, h40_high );
  fH40toward_bin->SetStats ( kTRUE );
  fH40toward_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH40toward_bin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40toward_bin->GetXaxis()->SetTitleColor ( 1 );
  fH40toward_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40toward_bin );

  fH40away_bin = new TH1D ( "histo40away_bin", "total particles fNPart w.r.t PTmax (pt of "
                            "leading particle from jet1) - AWAY",
                            h40_nbin, h40_low, h40_high );
  fH40away_bin->SetStats ( kTRUE );
  fH40away_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH40away_bin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40away_bin->GetXaxis()->SetTitleColor ( 1 );
  fH40away_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40away_bin );

  fH40transmin_bin = new TH1D ( "histo40transmin_bin", "total particles fNPart w.r.t PTmax (pt of "
                                "leading particle from jet1) - TRANSMIN",
                                h40_nbin, h40_low, h40_high );
  fH40transmin_bin->SetStats ( kTRUE );
  fH40transmin_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH40transmin_bin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40transmin_bin->GetXaxis()->SetTitleColor ( 1 );
  fH40transmin_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40transmin_bin );

  fH40transmax_bin = new TH1D ( "histo40transmax_bin", "total particles fNPart w.r.t PTmax (pt of "
                                "leading particle from jet1) - TRANSMAX",
                                h40_nbin, h40_low, h40_high );
  fH40transmax_bin->SetStats ( kTRUE );
  fH40transmax_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH40transmax_bin->GetYaxis()->SetTitle ( "<fNPart> (accepted)" );
  fH40transmax_bin->GetXaxis()->SetTitleColor ( 1 );
  fH40transmax_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH40transmax_bin );

  //____________________________________________________________________________________
  Int_t h41_nbin = 200;
  Double_t h41_binwidth = 1.;
  Double_t h41_low = 0.;
  Double_t h41_high = h41_low + h41_binwidth * h41_nbin;
  fH41 = new TProfile ( "histo41", "PTsum w.r.t PTmax (pt of leading particle from jet1)", h41_nbin, h41_low, h41_high );
  fH41->SetStats ( kTRUE );
  fH41->GetXaxis()->SetTitle ( "PTmax" );
  fH41->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41->GetXaxis()->SetTitleColor ( 1 );
  fH41->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41 );

  fH41toward = new TProfile ( "histo41toward", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TOWARD", h41_nbin,
                              h41_low, h41_high );
  fH41toward->SetStats ( kTRUE );
  fH41toward->GetXaxis()->SetTitle ( "PTmax" );
  fH41toward->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41toward->GetXaxis()->SetTitleColor ( 1 );
  fH41toward->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41toward );

  fH41away = new TProfile ( "histo41away", "PTsum w.r.t PTmax (pt of leading particle from jet1) - AWAY", h41_nbin,
                            h41_low, h41_high );
  fH41away->SetStats ( kTRUE );
  fH41away->GetXaxis()->SetTitle ( "PTmax" );
  fH41away->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41away->GetXaxis()->SetTitleColor ( 1 );
  fH41away->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41away );

  fH41transmin = new TProfile ( "histo41transmin", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TRANSMIN",
                                h41_nbin, h41_low, h41_high );
  fH41transmin->SetStats ( kTRUE );
  fH41transmin->GetXaxis()->SetTitle ( "PTmax" );
  fH41transmin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41transmin->GetXaxis()->SetTitleColor ( 1 );
  fH41transmin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41transmin );

  fH41transmax = new TProfile ( "histo41transmax", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TRANSMAX",
                                h41_nbin, h41_low, h41_high );
  fH41transmax->SetStats ( kTRUE );
  fH41transmax->GetXaxis()->SetTitle ( "PTmax" );
  fH41transmax->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41transmax->GetXaxis()->SetTitleColor ( 1 );
  fH41transmax->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41transmax );

  //____________________________________________________________________________________
  fH41_bin =
    new TH1D ( "histo41_bin", "PTsum w.r.t PTmax (pt of leading particle from jet1)", h41_nbin, h41_low, h41_high );
  fH41_bin->SetStats ( kTRUE );
  fH41_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH41_bin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41_bin->GetXaxis()->SetTitleColor ( 1 );
  fH41_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41_bin );

  fH41toward_bin = new TH1D ( "histo41toward_bin", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TOWARD",
                              h41_nbin, h41_low, h41_high );
  fH41toward_bin->SetStats ( kTRUE );
  fH41toward_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH41toward_bin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41toward_bin->GetXaxis()->SetTitleColor ( 1 );
  fH41toward_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41toward_bin );

  fH41away_bin = new TH1D ( "histo41away_bin", "PTsum w.r.t PTmax (pt of leading particle from jet1) - AWAY", h41_nbin,
                            h41_low, h41_high );
  fH41away_bin->SetStats ( kTRUE );
  fH41away_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH41away_bin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41away_bin->GetXaxis()->SetTitleColor ( 1 );
  fH41away_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41away_bin );

  fH41transmin_bin = new TH1D ( "histo41transmin_bin", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TRANSMIN",
                                h41_nbin, h41_low, h41_high );
  fH41transmin_bin->SetStats ( kTRUE );
  fH41transmin_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH41transmin_bin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41transmin_bin->GetXaxis()->SetTitleColor ( 1 );
  fH41transmin_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41transmin_bin );

  fH41transmax_bin = new TH1D ( "histo41transmax_bin", "PTsum w.r.t PTmax (pt of leading particle from jet1) - TRANSMAX",
                                h41_nbin, h41_low, h41_high );
  fH41transmax_bin->SetStats ( kTRUE );
  fH41transmax_bin->GetXaxis()->SetTitle ( "PTmax" );
  fH41transmax_bin->GetYaxis()->SetTitle ( "PTsum (GeV/c)" );
  fH41transmax_bin->GetXaxis()->SetTitleColor ( 1 );
  fH41transmax_bin->SetMarkerStyle ( kFullCircle );
  fOutput->Add ( fH41transmax_bin );

  // =========== Switch on Sumw2 for all histos ===========
  for ( Int_t i = 0; i < fOutput->GetEntries(); ++i )
      {
      TH1 *h1 = dynamic_cast<TH1 *> ( fOutput->At ( i ) );

      if ( h1 )
          {
          h1->Sumw2();
          continue;
          }

      TProfile *hprof1 = dynamic_cast<TProfile *> ( fOutput->At ( i ) );

      if ( hprof1 )
          {
          hprof1->Sumw2();
          }
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
  typedef std::pair<Double_t, Int_t> ptidx_pair;

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
  typedef std::pair<Double_t, Int_t> ptidx_pair;

  Int_t entries = trackscont->GetNAcceptedParticles();

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;
  pair_list.reserve ( entries );

  for ( Int_t i_entry = 0; i_entry < entries; i_entry++ )
      {
      AliVParticle *track = trackscont->GetNextAcceptParticle ( i_entry );

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
      if ( index == array[i] )
          {
          return kTRUE;
          }
      }

  return kFALSE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::ExecOnce()
  {
  AliAnalysisTaskEmcalJet::ExecOnce();

  if ( fJetsCont && fJetsCont->GetArray() == 0 )
      {
      fJetsCont = 0;
      }

  if ( fTracksCont && fTracksCont->GetArray() == 0 )
      {
      fTracksCont = 0;
      }

  if ( fCaloClustersCont && fCaloClustersCont->GetArray() == 0 )
      {
      fCaloClustersCont = 0;
      }
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDF::Terminate ( Option_t * )
  {
  // Called once at the end of the analysis.
  // Update pointers reading them from the output slot
  fOutput = dynamic_cast<TList *> ( GetOutputData ( 0 ) );
  }

// kate: indent-mode none; indent-width 2; replace-tabs on;
