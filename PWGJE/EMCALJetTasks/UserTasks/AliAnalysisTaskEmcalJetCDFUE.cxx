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

#include "AliAnalysisTaskEmcalJetCDFUE.h"

using namespace std;

/// \cond CLASSIMP
ClassImp ( AliAnalysisTaskEmcalJetCDFUE );
/// \endcond

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDFUE::AliAnalysisTaskEmcalJetCDFUE()
  : AliAnalysisTaskEmcalJet ( "AliAnalysisTaskEmcalJetCDFUE", kTRUE ),
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
    fH23Toward ( NULL ),
    fH23Transverse_min ( NULL ),
    fH23Transverse_max ( NULL ),
    fH23Away ( NULL ),
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
    fCaloClustContArray ( NULL ),
    fJet1 ( NULL ),
    fNJets_accepted ( 0 ),
    fNaccPart ( 0 ),
    fNaccClus ( 0 ),
    fEvPt ( -99.99 ),
    fJet1_sorted_idxvec(),
    fEvent_sorted_idxvec()
  {
  // Default constructor.
  fDebug = AliLog::GetGlobalDebugLevel();
  // SetMakeGeneralHistograms ( kTRUE );
  }

//________________________________________________________________________
AliAnalysisTaskEmcalJetCDFUE::AliAnalysisTaskEmcalJetCDFUE ( const char *name )
  : AliAnalysisTaskEmcalJet ( name, kTRUE ),
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
    fH23Toward ( NULL ),
    fH23Transverse_min ( NULL ),
    fH23Transverse_max ( NULL ),
    fH23Away ( NULL ),
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
AliAnalysisTaskEmcalJetCDFUE::~AliAnalysisTaskEmcalJetCDFUE()
  {
  // Destructor.
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDFUE::Run()
  {
  Int_t idx_jet_container = 0;

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


  //___________________________________________
  // jet1 : Sorting by p_T jet constituents
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  fTracksContArray = fTracksCont->GetArray();
  fJet1_sorted_idxvec = fJet1->SortConstituentsPt ( fTracksContArray );

  //__________________________________________________________________
  // sorting the EVENT _ACCEPTED_ tracks by pt
  fEvent_sorted_idxvec = SortTracksPt ( fTracksCont );

  // Run analysis code here, if needed. It will be executed before FillHistograms().
  return kTRUE; // If return kFALSE FillHistogram() will NOT be executed.
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCDFUE::FillHistograms()
  {
  // consts used in analysis
  Double_t const kPI_3 = TMath::Pi() / 3.;
  Double_t const k2PI_3 = 2 * kPI_3;

  Double_t jet1_pt = fJet1->Pt();
  UInt_t jet1_npart = fJet1->GetNumberOfTracks();

  AliVParticle *jet1_trklead = fJet1->GetLeadingTrack ( fTracks );
  Double_t jet1_ptmax = jet1_trklead->Pt();

  Int_t counter_part = 0;
  Double_t counter_pt = 0.;
  Int_t track_idx = -999;            // index variable for tracks
  TArrayI jet1_idx_list ( jet1_npart ); // TArrayI of indexes of tracks from leading jet

  //___________________________________________________________________________
  AliVParticle *track = NULL;

  // reset counter for new usage
  counter_part = 0;
  counter_pt = 0.;
  track_idx = -999;

  // parsing accepted tracks in EVENT in decreasing order of Pt //
  for ( Size_t i = 0; i < fNaccPart; i++ ) // replace the index order by the sorted array
      {
      track_idx = fEvent_sorted_idxvec.at ( i );
      track = fTracksCont->GetAcceptParticle ( track_idx );

      if ( !track )
          {
          std::cout << "track not retrieved from fTracksCont" << std::endl;
          continue;
          }

      // pt of the current track
      Double_t track_pt = track->Pt();

      // dphi between leading track (max pt track from leading jet) and current track - dphi to (-pi,pi)
      Double_t dphi_part = TVector2::Phi_mpi_pi ( track->Phi() - jet1_trklead->Phi() );
//       Double_t dphi_part_absdeg = TMath::RadToDeg() * TMath::Abs ( dphi_part );

      // dphi between leading jet and current track - dphi to (-pi,pi)
      Double_t dphi_part_jet1 = TVector2::Phi_mpi_pi ( track->Phi() - fJet1->Phi() );
//       Double_t dphi_part_jet1_absdeg = TMath::RadToDeg() * TMath::Abs ( dphi_part_jet1 );

      ++counter_part;         // next particle
      counter_pt += track_pt; // next particle pt

      fH23->Fill ( track_pt ); //  Pt Distribution of particles in event

      // dphi track to jet1 distribution (total and per toward,away,transverse regions)
      fH21->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs P_{T}(jet1)
      fH21_bin->Fill ( jet1_pt );           // N (in the event - including jet1) vs P_{T}(jet1)

      fH22->Fill ( jet1_pt, counter_pt );   // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)
      fH22_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)

      if ( ( dphi_part_jet1 > ( -1 ) * kPI_3 ) && ( dphi_part_jet1 < kPI_3 ) )
          {
          fH21Toward->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs p_{T}(jet1)
          fH21Toward_bin->Fill ( jet1_pt );           // N (in the event - including jet1) vs p_{T}(jet1)

          fH22Toward->Fill ( jet1_pt, counter_pt );   // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)
          fH22Toward_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)

          fH23Toward->Fill ( track_pt );              // Pt Distribution of particles
          }
      else if ( ( dphi_part_jet1 <= ( -1 ) * kPI_3 ) && ( dphi_part_jet1 > ( -1 ) * k2PI_3 ) )
          {
          fH21Transverse_min->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs p_{T}(jet1)
          fH21Transverse_min_bin->Fill ( jet1_pt );           // N (in the event - including jet1) vs p_{T}(jet1)

          fH22Transverse_min->Fill ( jet1_pt, counter_pt );   // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)
          fH22Transverse_min_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)

          fH23Transverse_min->Fill ( track_pt );              // Pt Distribution of particles
          }
      else if ( ( dphi_part_jet1 >= kPI_3 ) && ( dphi_part_jet1 < k2PI_3 ) )
          {
          fH21Transverse_max->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs p_{T}(jet1)
          fH21Transverse_max_bin->Fill ( jet1_pt );           // N (in the event - including jet1) vs p_{T}(jet1)

          fH22Transverse_max->Fill ( jet1_pt, counter_pt );   // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)
          fH22Transverse_max_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including jet1) vs p_{T}(jet1)

          fH23Transverse_max->Fill ( track_pt );              // Pt Distribution of particles
          }
      else if ( ( dphi_part_jet1 >= k2PI_3 ) || ( dphi_part_jet1 <= ( -1 ) * k2PI_3 ) )
          {
          fH21Away->Fill ( jet1_pt, counter_part ); // N (in the event - including jet1) vs P_{T}(jet1)
          fH21Away_bin->Fill ( jet1_pt );           // N (in the event - including jet1) vs P_{T}(jet1)

          fH22Away->Fill ( jet1_pt, counter_pt );   // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)
          fH22Away_bin->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - including jet1) vs P_{T}(jet1)

          fH23Away->Fill ( track_pt );              // Pt Distribution of particles
          }


      // NEW UE histos
      // dphi track to leading track distribution (total and per toward,away,transverse regions)
      fH40->Fill ( jet1_ptmax, counter_part ); // total particles fNPart w.r.t PTmax (pt of leading particle from leading jet)
      fH40_bin->Fill ( jet1_ptmax );           // total particles fNPart w.r.t PTmax (pt of leading particle from leading jet)

      fH41->Fill ( jet1_ptmax, counter_pt );   // PTsum w.r.t PTmax
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


      // Condition if track is NOT in leading jet
      if ( !IdxInArray ( track_idx, jet1_idx_list ) )
          {
          fH21_bin_wojet1->Fill ( jet1_pt );           // N (in the event - excluding jet1) vs P_{T}(jet1)
          fH22_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - excluding jet1) vs P_{T}(jet1)

          if ( ( dphi_part_jet1 > ( -1 ) * kPI_3 ) && ( dphi_part_jet1 < kPI_3 ) )
              {
              fH21Toward_bin_wojet1->Fill ( jet1_pt );           // N (in the event - excluding jet1) vs P_{T}(jet1)
              fH22Toward_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - excluding jet1) vs P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 <= ( -1 ) * kPI_3 ) && ( dphi_part_jet1 > ( -1 ) * k2PI_3 ) )
              {
              fH21Transverse_min_bin_wojet1->Fill ( jet1_pt );           // N (in the event - excluding jet1) vs P_{T}(jet1)
              fH22Transverse_min_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - excluding jet1) vs P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 >= kPI_3 ) && ( dphi_part_jet1 < k2PI_3 ) )
              {
              fH21Transverse_max_bin_wojet1->Fill ( jet1_pt );           // N (in the event - excluding jet1) vs P_{T}(jet1)
              fH22Transverse_max_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - excluding jet1) vs P_{T}(jet1)
              }
          else if ( ( dphi_part_jet1 >= k2PI_3 ) || ( dphi_part_jet1 <= ( -1 ) * k2PI_3 ) )
              {
              fH21Away_bin_wojet1->Fill ( jet1_pt );           // N (in the event - excluding jet1) vs P_{T}(jet1)
              fH22Away_bin_wojet1->Fill ( jet1_pt, track_pt ); // PT_{sum}(in the event - excluding jet1) vs P_{T}(jet1)
              }
          }

      }

  track = NULL;

  // post data at every processing
  PostData ( 1, fOutput ); // Post data for ALL output slots > 0 here.
  return kTRUE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDFUE::UserCreateOutputObjects()
  {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  //  Histograms
  fOutput->SetOwner ( kTRUE );

  // Create the list of histograms. Only the list is owned.


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
Double_t AliAnalysisTaskEmcalJetCDFUE::DeltaR ( const AliVParticle *part1, const AliVParticle *part2 )
  {
  // Helper function to calculate the distance between two jets or a jet and
  // particle
  Double_t dPhi = part1->Phi() - part2->Phi();
  Double_t dEta = part1->Eta() - part2->Eta();
  dPhi = TVector2::Phi_mpi_pi ( dPhi );

  return TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
  }

//__________________________________________________________________________________________________
std::vector<Int_t> AliAnalysisTaskEmcalJetCDFUE::SortTracksPt ( AliVEvent *event ) const
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
std::vector<Int_t> AliAnalysisTaskEmcalJetCDFUE::SortTracksPt ( AliParticleContainer *trackscont ) const
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
      AliVParticle *track = trackscont->GetNextAcceptParticle ( );

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
Bool_t AliAnalysisTaskEmcalJetCDFUE::IdxInArray ( Int_t index, TArrayI &array )
  {
  for ( Size_t i = 0; i < array.GetSize(); i++ )
      {
      if ( index == array[i] ) { return kTRUE; }
      }

  return kFALSE;
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDFUE::ExecOnce()
  {
  AliAnalysisTaskEmcalJet::ExecOnce();

  if ( fJetsCont && fJetsCont->GetArray() == 0 ) { fJetsCont = 0; }

  if ( fTracksCont && fTracksCont->GetArray() == 0 ) { fTracksCont = 0; }

  if ( fCaloClustersCont && fCaloClustersCont->GetArray() == 0 ) { fCaloClustersCont = 0; }
  }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCDFUE::Terminate ( Option_t * )
  {
  // Called once at the end of the analysis.
  // Update pointers reading them from the output slot
  fOutput = dynamic_cast<AliEmcalList*> ( GetOutputData (0) );
  }

// kate: indent-mode none; indent-width 2; replace-tabs on;
