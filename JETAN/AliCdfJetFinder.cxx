/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to usec, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//---------------------------------------------------------------------
// Jet Finder based on CDF algortihm
// Charged jet evolution and the underlying event in proton-antiproton collisions at 1.8 TeV
// Physical Review D, vol. 65, Issue 9, id. 092002
// http://www.phys.ufl.edu/~rfield/cdf/chgjet/chgjet_intro.html
// Authors : Adrian.Sevcenco@cern.ch (adriansev@spacescience.ro )
//           Daniel.Felea@cern.ch    (dfelea@spacescience.ro)
//           Ciprian.Mihai.Mitu@cern.ch  (mcm@spacescience.ro)
//           magali.estienne@subatech.in2p3.fr &
//           alexandre.shabetai@cern.ch (Modification of the input object (reader/finder splitting))
// ** 2011
// Modified accordingly to reader/finder splitting and new handling of neutral information
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>
#include <TBits.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TVector2.h>

#include "AliAODJet.h"
#include "AliJetFinder.h"
#include "AliJetCalTrk.h"
#include "AliCdfJetFinder.h"
#include "AliCdfJetHeader.h"

ClassImp(AliCdfJetFinder)

//______________________________________________________________________________
AliCdfJetFinder::AliCdfJetFinder():
  AliJetFinder(),
  fHistos(0),
  fAODwrite(0),
  fAODtracksWrite(0),
  fAnalyseJets(0),
  fNJets(0),
  fNPart(0),
  fNInC(0),
  fNInN(0),
  fRadius(0.7),
  fMinJetParticles(1),
  fJetPtCut(0.),
  fVectParticle(NULL),
  fVectJet(NULL),
  fPtArray(NULL),
  fIdxArray(NULL)
{  
  // Default constructor
}

//______________________________________________________________________________
AliCdfJetFinder::~AliCdfJetFinder()
{
  // destructor
  Clean();

}

//______________________________________________________________________________
void AliCdfJetFinder::CreateOutputObjects(TList * const histos)
{
  // Create the list of histograms. Only the list is owned.
  fHistos = histos;

  //  gStyle->SetOptStat(11111111);

  TH1F *h1 = new TH1F ("histo1", "Pt distribution of jets", 200, 0,200);  // 1GeV/bin
  h1->SetStats(kTRUE);
  h1->GetXaxis()->SetTitle("P_{T} of jets");
  h1->GetYaxis()->SetTitle("Number of jets");
  h1->GetXaxis()->SetTitleColor(1);
  h1->SetMarkerStyle(kFullCircle);
  fHistos->Add(h1);

  TH1F *h2 = new TH1F ("histo2", "Eta distribution of jets", 240, -1.2,1.2); // 1 unit of rapidity / 100 bin
  h2->SetStats(kTRUE);
  h2->GetXaxis()->SetTitle("Eta of jets");
  h2->GetYaxis()->SetTitle("Number of jets");
  h2->GetXaxis()->SetTitleColor(1);
  h2->SetMarkerStyle(kFullCircle);
  fHistos->Add(h2);

  TH1F *h3 = new TH1F ("histo3", "Phi distribution of jets", 400, -4,4);
  h3->SetStats(kTRUE);
  h3->GetXaxis()->SetTitle("Phi of jets");
  h3->GetYaxis()->SetTitle("Number of jets");
  h3->GetXaxis()->SetTitleColor(1);
  h3->SetMarkerStyle(kFullCircle);
  fHistos->Add(h3);

  TH1F *h4 = new TH1F ("histo4", "Multiplicity of jets", 40, 0,40);  // 1 unit of multiplicity /bin
  h4->SetStats(kTRUE);
  h4->GetXaxis()->SetTitle("Particles in jets");
  h4->GetYaxis()->SetTitle("Number of jets");
  h4->GetXaxis()->SetTitleColor(1);
  h4->SetMarkerStyle(kFullCircle);
  fHistos->Add(h4);

  TH1F *h5 = new TH1F ("histo5", "Distribution of jets in events", 100, 0,100);
  h5->SetStats(kTRUE);
  h5->GetXaxis()->SetTitle("Number of jets");
  h5->GetYaxis()->SetTitle("Number of events");
  h5->GetXaxis()->SetTitleColor(1);
  h5->SetMarkerStyle(kFullCircle);
  fHistos->Add(h5);

  TH1F *h6 = new TH1F ("histo6", "Jet1 Charged Multiplicity Distribution", 30, 0,30);
  h6->SetStats(kTRUE);
  h6->GetXaxis()->SetTitle("N_{chg}");
  h6->GetYaxis()->SetTitle("Number of jets");
  h6->GetXaxis()->SetTitleColor(1);
  h6->SetMarkerStyle(kFullCircle);
  fHistos->Add(h6);

  TProfile * h7 = new TProfile ("histo7","N_{chg}(jet1) vs P_{T}(charged jet1)", 200, 0. ,200. , 0.,200. ) ;
  h7->SetStats(kTRUE);
  h7->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h7->GetYaxis()->SetTitle("<N_{chg}(jet1)> in 1 GeV/c bin");
  h7->GetXaxis()->SetTitleColor(1);
  h7->SetMarkerStyle(kFullCircle);
  fHistos->Add(h7);

  TH1F *h8 = new TH1F ("histo8", "Charge momentum distribution for leading jet", 120, 0 , 1.2);
  h8->SetStats(kTRUE);
  h8->GetXaxis()->SetTitle("Jets");
  h8->GetYaxis()->SetTitle("Particle distribution");
  h8->GetXaxis()->SetTitleColor(1);
  h8->SetMarkerStyle(kFullCircle);
  fHistos->Add(h8);

  TProfile *h9 = new TProfile ("histo9", "N_{chg} vs the Azimuthal Angle from Charged Jet1", 50 , 0. , 180. , 0 , 20 );
  h9->SetStats(kTRUE);
  h9->GetXaxis()->SetTitle("#Delta#phi (degrees)");
  h9->GetYaxis()->SetTitle("<N_{chg}> in 3.6 degree bin");
  h9->GetXaxis()->SetTitleColor(1);
  h9->SetMarkerStyle(kFullCircle);
  fHistos->Add(h9);

  TProfile *h10 = new TProfile ("histo10", "P_{T} sum vs the Azimuthal Angle from Charged Jet1", 50 , 0. , 180. , 0 , 100 );
  h10->SetStats(kTRUE);
  h10->GetXaxis()->SetTitle("#Delta#phi (degrees)");
  h10->GetYaxis()->SetTitle("<P_{T} sum> in 3.6 degree bin");
  h10->GetXaxis()->SetTitleColor(1);
  h10->SetMarkerStyle(kFullCircle);
  fHistos->Add(h10);

  TH1F *h11 = new TH1F ("histo11", " \"Transverse\" Pt Distribution ", 70, 0 , 14);
  h11->SetStats(kTRUE);
  h11->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h11->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h11->GetXaxis()->SetTitleColor(1);
  h11->SetMarkerStyle(kFullCircle);
  fHistos->Add(h11);

  TH1F *h20 = new TH1F ("histo20", "Distribution of R in leading jet", 400, 0.,4.);
  h20->SetStats(kTRUE);
  h20->GetXaxis()->SetTitle("R [formula]");
  h20->GetYaxis()->SetTitle("dN/dR");
  h20->GetXaxis()->SetTitleColor(1);
  h20->SetMarkerStyle(kFullCircle);
  fHistos->Add(h20);

  TProfile * h21 = new TProfile ("histo21","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 30. ) ;
  h21->SetStats(kTRUE);
  h21->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21->GetXaxis()->SetTitleColor(1);
  h21->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21);

  TProfile * h22 = new TProfile ("histo22","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22->SetStats(kTRUE);
  h22->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22->GetXaxis()->SetTitleColor(1);
  h22->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22);

  TProfile * h21Toward = new TProfile ("histo21_toward","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21Toward->SetStats(kTRUE);
  h21Toward->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21Toward->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21Toward->GetXaxis()->SetTitleColor(1);
  h21Toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21Toward);

  TProfile * h21Transverse = new TProfile ("histo21_transverse","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21Transverse->SetStats(kTRUE);
  h21Transverse->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21Transverse->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21Transverse->GetXaxis()->SetTitleColor(1);
  h21Transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21Transverse);

  TProfile * h21Away = new TProfile ("histo21_away","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21Away->SetStats(kTRUE);
  h21Away->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21Away->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21Away->GetXaxis()->SetTitleColor(1);
  h21Away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21Away);

  TProfile * h22Toward = new TProfile ("histo22_toward","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22Toward->SetStats(kTRUE);
  h22Toward->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22Toward->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22Toward->GetXaxis()->SetTitleColor(1);
  h22Toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22Toward);

  TProfile * h22Transverse = new TProfile ("histo22_transverse","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22Transverse->SetStats(kTRUE);
  h22Transverse->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22Transverse->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22Transverse->GetXaxis()->SetTitleColor(1);
  h22Transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22Transverse);

  TProfile * h22Away = new TProfile ("histo22_away","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22Away->SetStats(kTRUE);
  h22Away->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22Away->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22Away->GetXaxis()->SetTitleColor(1);
  h22Away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22Away);

  TH1F *h23Toward = new TH1F ("histo23_toward","'Toward' Pt Distribution of charged particles", 200, 0., 14.);
  h23Toward->SetStats(kTRUE);
  h23Toward->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23Toward->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23Toward->GetXaxis()->SetTitleColor(1);
  h23Toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23Toward);

  TH1F *h23Transverse = new TH1F ("histo23_transverse","'Transverse' Pt Distribution of charged particles", 200, 0., 14.);
  h23Transverse->SetStats(kTRUE);
  h23Transverse->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23Transverse->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23Transverse->GetXaxis()->SetTitleColor(1);
  h23Transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23Transverse);

  TH1F *h23Away = new TH1F ("histo23_away","'Away' Pt Distribution of charged particles", 200, 0., 14.);
  h23Away->SetStats(kTRUE);
  h23Away->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23Away->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23Away->GetXaxis()->SetTitleColor(1);
  h23Away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23Away);

  TProfile * h24 = new TProfile ("histo24","Jet1 Size vs P_{T}(charged jet1)", 200, 0., 50. , 0., 0.5) ;
  h24->SetStats(kTRUE);
  h24->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h24->GetYaxis()->SetTitle("<R(chgjet1)> in 1 GeV/c bin");
  h24->GetXaxis()->SetTitleColor(1);
  h24->SetMarkerStyle(kFullCircle);
  fHistos->Add(h24);

  TProfile * h25 = new TProfile ("histo25","Jet1 Size vs P_{T}(charged jet1)", 200, 0., 50. , 0., 0.5) ;
  h25->SetStats(kTRUE);
  h25->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h25->GetYaxis()->SetTitle("<R(chgjet1)> in 1 GeV/c bin");
  h25->GetXaxis()->SetTitleColor(1);
  h25->SetMarkerStyle(kFullCircle);
  fHistos->Add(h25);

  TProfile *h26 = new TProfile ("histo26", "N_{chg} vs the Distance R from Charged Jet1", 30, 0., 0.6, 0., 0.8);
  h26->SetStats(kTRUE);
  h26->GetXaxis()->SetTitle("Distance R");
  h26->GetYaxis()->SetTitle("<N_{chg}> in 0.02 bin");
  h26->GetXaxis()->SetTitleColor(1);
  h26->SetMarkerStyle(kFullCircle);
  fHistos->Add(h26);

  TProfile *h27 = new TProfile ("histo27", "N_{chg} vs the Distance R from Charged Jet1", 30, 0., 0.6, 0., 0.8);
  h27->SetStats(kTRUE);
  h27->GetXaxis()->SetTitle("Distance R");
  h27->GetYaxis()->SetTitle("<N_{chg}> in 0.02 bin");
  h27->GetXaxis()->SetTitleColor(1);
  h27->SetMarkerStyle(kFullCircle);
  fHistos->Add(h27);

  TProfile *h28 = new TProfile ("histo28", "PT_{sum} vs the Distance R from Charged Jet1", 30, 0., 0.6, 0.01, 10.);
  h28->SetStats(kTRUE);
  h28->GetXaxis()->SetTitle("Distance R");
  h28->GetYaxis()->SetTitle("<PT_{sum} (GeV/c)> in 0.02 bin");
  h28->GetXaxis()->SetTitleColor(1);
  h28->SetMarkerStyle(kFullCircle);
  fHistos->Add(h28);

  TProfile *h29 = new TProfile ("histo29", "PT_{sum} vs the Distance R from Charged Jet1", 30, 0., 0.6, 0.01, 10.);
  h29->SetStats(kTRUE);
  h29->GetXaxis()->SetTitle("Distance R");
  h29->GetYaxis()->SetTitle("<PT_{sum} (GeV/c)> in 0.02 bin");
  h29->GetXaxis()->SetTitleColor(1);
  h29->SetMarkerStyle(kFullCircle);
  fHistos->Add(h29);

}

//______________________________________________________________________________
void AliCdfJetFinder::FindJets()
{
  // Jet Algorithm:
  //  * Order all charged particles according to their PT.
  //  * Start with the highest PT particle and include in the "jet" all particles within the "radius" R = 0.7
  //     (considering each particle in the order of decreasing PT and recalculating the centroid of the jet after
  //     each new particle is added to the jet).
  //  * Go to the next highest PT particle (not already included in a jet) and include in the "jet" all particles
  //     (not already included in a jet) within the radius R =0.7.
  //  * Continue until all particles are in a "jet".
  if (fDebug) { printf("AliCDJetfinder::FindJets() %d \n", __LINE__ ); }
  AliCdfJetHeader *header = (AliCdfJetHeader*)fHeader;

  if (header)
    {
      fDebug            = header->GetDebug();
      fAODwrite         = header->IsAODwrite() ;       // write jets to AOD
      fAODtracksWrite   = header->IsAODtracksWrite() ; // write jet tracks to AOD
      fRadius           = header->GetRadius();      // get Radius from jet finder header
      fMinJetParticles  = header->GetMinPartJet (); // get minimum multiplicity of an jet
      fJetPtCut         = header->GetJetPtCut ();   // get minimum of jet pt
      fAnalyseJets      = header->GetAnalyseJets(); // get analyse jet
    }
  else
    { cout << "Header not found" << endl; return; }

  InitData();

  if (!fNPart) { 
    if (fDebug) {
      cout << "entries = 0 ; Event empty !!!" << endl ;
    }
    // no need to call clean, InitData does not 
    // create pointers if npart == 0
    return; 
  } // if event empty then exit

  FindCones();
 
  ComputeConesWeight();
 
  if (fAODwrite) { 
    if(fDebug)cout << "Writing AOD" << endl ; 
    WriteJets();
  }
 
  if (fAnalyseJets) AnalizeJets();
 
  Clean();

}

//______________________________________________________________________________
void AliCdfJetFinder::InitData()
{
  // initialisation of variables and data members

  if( fHeader->GetDebug() && fCalTrkEvent->GetNCalTrkTracks()  == 0) { cout << "No charged tracks found" << endl; }

  fNPart = fCalTrkEvent->GetNCalTrkTracks() ;

  if ( fCalTrkEvent->GetNCalTrkTracks() ) { return; } // if event empty then exit

  fVectParticle = new varContainer* [fNPart]; // container for Particles

  fPtArray  = new Double_t [fCalTrkEvent->GetNCalTrkTracks()] ; 
  fIdxArray = new Int_t    [fCalTrkEvent->GetNCalTrkTracks()] ; // index array of sorted pts

  // initialisation of momentum and index arrays
  for (  Int_t i = 0 ; i <fCalTrkEvent->GetNCalTrkTracks() ; i++ )
    {// SORTING STEP :: fPtArray with data from CalTrkTracks

      // INITIALISATION of local arrays for temporary storage
      varContainer *aParticle = new varContainer;
      aParticle->pt   = fCalTrkEvent->GetCalTrkTrack(i)->GetPt();
      aParticle->eta  = fCalTrkEvent->GetCalTrkTrack(i)->GetEta();
      aParticle->phi  = TVector2::Phi_mpi_pi ( fCalTrkEvent->GetCalTrkTrack(i)->GetPhi() ); // normalize to -pi,pi
      aParticle->njet = -999;

      fVectParticle[i] = aParticle;  // vector of Particles

      // initializing arrays
      fIdxArray [i] = -999 ;
      fPtArray [i] = aParticle->pt ;
    }

  TMath::Sort ( fNPart, fPtArray, fIdxArray ) ; // get a sorted array of indexes

}

//______________________________________________________________________________
void AliCdfJetFinder::FindCones()
{
  // parsing of particles in event and estlabish jets (label them with jet index)

  Double_t  ptSeed = 0. , etaSeed = 0. , phiSeed = 0. ; // leading particle params
  Double_t pttmp = 0. , etatmp = 0. , phitmp = 0. ; // temporary variables to be used in various calculations
  Double_t deta = 0. , dphi = 0. , dcomputed = 0. ;
  Bool_t injet = 0 ;

  fNJets = -1 ; // n jets in this event
  Int_t idxPtSort = -1 ;  // index of array of sorted pt indexes

  if (fDebug) { cout << "\n\n\n\n\n\n------------------\nBegin Event Analysis\n------------------\n\n" << endl ;}

  if(fDebug)cout << "fNPart = " << fNPart << endl;

  TBits lkupTable ( fNPart ) ;  // bit container ; 1-to-1 corespondence with fIdxArray

  while ( lkupTable.CountBits() != (UInt_t)fNPart )
    { // loop over particles in event until all flags are set
      UInt_t firstnonflagged = lkupTable.FirstNullBit() ; // set the index to the first NON flagged bit ; less conditions

      if(fDebug)cout << "\n\nfirst_non_flagged : " << firstnonflagged << endl;

      ++fNJets; // incrementing the jet counter
      if (fDebug) { printf("JET %d \n", fNJets); }

      ptSeed = 0. ; etaSeed = 0. ; phiSeed = 0. ;  // reseting leading particle params

      for (  UInt_t ipart = firstnonflagged ; ipart < (UInt_t)fNPart ; ipart++ )
	{// iteration over particles in event
	  // the loop is done over sorted array of pt
	  idxPtSort = fIdxArray[ipart] ;  // index of particle ! fIdxArray is an index list pt sorted

	  if ( lkupTable.TestBitNumber(ipart) ) { continue; } // if 4vector is already flagged skip it

	  //init computed and used vars
	  pttmp = 0. ; etatmp = 0. ; phitmp = 0. ;
	  deta = 0. ; dphi = 0. ; dcomputed = 0. ; injet = 0 ;

	  //taking info from fVectParticle ;
	  pttmp = fVectParticle[idxPtSort]->pt ;
	  etatmp = fVectParticle[idxPtSort]->eta ;
	  phitmp = fVectParticle[idxPtSort]->phi ;

	  if ( ipart == firstnonflagged )
	    {// this is first particle in event; leading particle
	      // begin the search around this particle in a fRadius

	      // CENTRE OF THE JET
	      ptSeed = pttmp ; etaSeed = etatmp ; phiSeed = phitmp ; // seeding the jet with first particle idxPtSort

	      lkupTable.SetBitNumber ( ipart ) ; // flag the index of particle in lkup_table
	      fVectParticle[idxPtSort]->njet = fNJets ; // associate particle with current jet number

	      if (fDebug) { printf("\nLeading particle :: particle index = %d ;  at sorted index %d ; in jet %d \n", idxPtSort, ipart, fNJets); }
	      if (fDebug) { printf("pt= %g ; eta= %g ; phi = %g \n", pttmp, etatmp, phitmp) ; }
	      if (fDebug) { lkupTable.Print() ;}

	      continue ; // skip to next particle
	    }

	  // condition to be in jet
	  deta = etatmp - etaSeed ;
	  dphi = TVector2::Phi_mpi_pi (phitmp - phiSeed) ; // computing dphi and normalizing to (0,2pi) interval in one step

	  dcomputed = TMath::Hypot(deta, dphi) ; // Distance(fRadius) to (eta,phi) seed

	  injet = ( ( fRadius - dcomputed ) >= 0.000000001 ) ? 1 : 0 ; // if r_computed is within jet_r in_jet == 1 else 0

	  if ( injet )
	    { // calculus of jet variables
	      lkupTable.SetBitNumber ( ipart ) ;  // flag the index of particle in lkup_table
	      fVectParticle[idxPtSort]->njet = fNJets ; // setting in particle list the associated jet

	      if (fDebug) { printf("\njet particle :: particle index = %d ; at sorted index %d ; in jet %d ; found at radius %g ;  \n", idxPtSort, ipart, fNJets, dcomputed); }
	      if (fDebug) { printf("pt= %g ; eta= %g ; phi = %g \n", pttmp, etatmp, phitmp) ; }
	      if (fDebug) { lkupTable.Print() ;}

	      continue ; // skip to next particle
	    }

	}
      // end of iteration over event; one jet definition of content ; jet parameters to be computed later
    }

}

//______________________________________________________________________________
void AliCdfJetFinder::ComputeConesWeight()
{
  // computing of jets Pt, Eta and Phi (centre of weight in (eta,phi) plane)
  // rescan the vector of particles by identify them by asociate jet number for computing of weight centre

  // JET CONTAINER
  fVectJet      = new varContainer* [fNJets]; // container for Jets

  Double_t ptJet, ptJet2 , etaJet , phiJet ; Int_t npartJet ;
  Double_t pttmp = 0. , etatmp = 0. , phitmp = 0. ; // temporary variables to be used in various calculations
  Int_t idxPtSort = -999 ;  // index of array of sorted pt indexes

  for(  Int_t jet = 0 ; jet < fNJets ; jet++ )
    {
      if (fDebug) { printf("\n\n--- Computing weight of Jet %d \n", jet ); }
      npartJet = 0 ; ptJet = 0. ; etaJet = 0. ; phiJet = 0. ; // reset variables for a new computation

      for (  Int_t ipart = 0 ; ipart < fNPart ; ipart++ )
	{// iteration over particles in event
	  // the loop is done over sorted array of pt
	  idxPtSort = fIdxArray[ipart] ;  // index of particle ! fIdxArray is an index list pt sorted

	  if ( fVectParticle[idxPtSort]->njet == jet )
	    {
	      ++npartJet; // incrementing the counter of jet particles

	      //taking info from fVectParticle ;
	      pttmp = fVectParticle[idxPtSort]->pt ;
	      etatmp = fVectParticle[idxPtSort]->eta ;
	      phitmp = TVector2::Phi_mpi_pi (fVectParticle[idxPtSort]->phi) ;

	      //      jet_new_angular_coordinate = jet_old_angular_coordinate * jet_old_pt / jet_new_pt +
	      //                                    part[i]_angular_coordinate * part[i]_pt/jet_new_pt

	      ptJet2 = ptJet + pttmp ;

	      etaJet = etaJet * ptJet / ptJet2 +  etatmp * pttmp / ptJet2 ;
	      phiJet = phiJet * ptJet / ptJet2 +  phitmp * pttmp / ptJet2 ;

	      ptJet = ptJet2 ;

	    }
	  // add a particle and recalculation of centroid
	}
      // end of 1 jet computation

      varContainer *aJet = new varContainer;  // Jet container
      aJet->pt = ptJet; aJet->eta = etaJet; aJet->phi = phiJet; aJet->njet = npartJet; // setting jet vars in container
      fVectJet[jet] = aJet;   // store the number of the jet(fNJets) and increment afterwards

      if (fDebug) { printf ("=== current jet %d : npartjet= %d ; pt_jet= %g ; eta_jet = %g ; phi_jet = %g \n\n\n",
			    jet,     npartJet,      ptJet,      etaJet,       phiJet ) ; }

    }
  //end loop over jets

}

//______________________________________________________________________________
void AliCdfJetFinder::WriteJets()  
{ 
  // Writing AOD jets and AOD tracks

  for(  Int_t jetnr = 0 ; jetnr < fNJets ; jetnr++ )
    {
      Double_t pt = 0., eta = 0., phi = 0., // jet variables
	px = 0., py = 0., pz = 0., en = 0.; // convert to 4-vector
      pt  = fVectJet[ jetnr ]->pt   ; // pt  of jet
      eta = fVectJet[ jetnr ]->eta  ; // eta of jet
      phi = fVectJet[ jetnr ]->phi  ; // phi of jet

      px = pt * TMath::Cos ( phi ) ;
      py = pt * TMath::Sin ( phi ) ;
      pz = pt / TMath::Tan ( 2.0 * TMath::ATan ( TMath::Exp ( -eta ) ) ) ;
      en = TMath::Sqrt ( px * px + py * py + pz * pz );

      AliAODJet jet (px, py, pz, en);


      if (fDebug) jet.Print("");

      if (fAODtracksWrite)
	{
	  for (  Int_t jetTrack = 0; jetTrack < fCalTrkEvent->GetNCalTrkTracks(); jetTrack++ )
	    {
	      // The first if condition below has to be checked
	      if ( fVectParticle[jetTrack]->njet == jetnr ) { jet.AddTrack(fCalTrkEvent->GetCalTrkTrack(jetTrack)->GetTrackObject()) ; }
	    }
	}
      // tracks REFs written in AOD
      AddJet(jet);

    }
  //jets vector parsed and written to AOD
}

//______________________________________________________________________________
void AliCdfJetFinder::AnalizeJets()
{
  // analyzing of jets and filling of histograms

  const Double_t kPI = TMath::Pi();
    
  //persistent pointer to histo20
  TH1F *hR = (TH1F*)fHistos->FindObject("histo20");

  Int_t   *jetsptidx = 0;     // sorted array of jets pt
  Double_t    *jetspt = 0;     // array of jets pts
  Int_t leadingjetindex = -1 ;   // index of leading jet from fVectJet
  Int_t partleadjet = 0 ; // number of particles in leading jet
  Double_t   ptleadjet = 0. ; // pt  of leading jet
  Double_t  etaleadjet = 0. ; // eta of leading jet
  Double_t  phileadjet = 0. ; // phi of leading jet

  jetsptidx = new Int_t    [fNJets] ;
  jetspt    = new Double_t [fNJets] ;

  //________________________________________________________________________________________
  //  Jet sorting and finding the leading jet that coresponds to cuts in pt and multiplicity
  //________________________________________________________________________________________

  // filing the idx_ptjets array
  if (fDebug) printf("List of unsorted jets:\n");
  for(  Int_t i = 0 ; i < fNJets ; i++ )
    {
      jetsptidx [i] = 0 ;
      jetspt    [i] = fVectJet[i]->pt ;
      if (fDebug) { cout << "   jet found: " << i << " npartjet=" << fVectJet[i]->njet << " ; jets_pt = " << jetspt[i] << endl; }
    }

  TMath::Sort ( fNJets, jetspt , jetsptidx ) ; // sorting pt of jets

  // selection of leading jet
  // looping over jets searching for __first__ one that coresponds to cuts
  for(  Int_t i = 0 ; i < fNJets ; i++ )
    {
      if ( ( fVectJet[ jetsptidx[i] ]->njet >= fMinJetParticles ) && ( fVectJet[ jetsptidx[i] ]->pt >= fJetPtCut ) )
	{
	  leadingjetindex = jetsptidx[i] ;
	  partleadjet = fVectJet[ leadingjetindex ]->njet ; // number of particles in leading jet
	  ptleadjet = fVectJet[ leadingjetindex ]->pt   ; // pt  of leading jet
	  etaleadjet = fVectJet[ leadingjetindex ]->eta  ; // eta of leading jet
	  phileadjet = fVectJet[ leadingjetindex ]->phi  ; // phi of leading jet

	  if (fDebug)
	    { printf("Leading jet %d : npart= %d ; pt= %g ; eta = %g ; phi = %g \n", leadingjetindex, partleadjet, ptleadjet, etaleadjet, phileadjet ); }

	  break ;
	}
    }
  // end of selection of leading jet



  //////////////////////////////////////////////////
  ////  Computing of values used in histograms
  //////////////////////////////////////////////////

  //___________________________________________________________________________
  // pt_sum of all particles in event
  //___________________________________________________________________________
  if (fDebug) cout << "Computing sum of pt in event" << endl ;
  Double_t ptsumevent = 0.;
  for (  Int_t i = 0 ; i< fNPart ; i++ ) { ptsumevent += fVectParticle[i]->pt ; }
  if (fDebug) printf ("Sum of all Pt in event : pt_sum_event = %g", ptsumevent) ;

  //___________________________________________________________________________
  // Filling an array with indexes of leading jet particles
  //___________________________________________________________________________
  Int_t * idxpartLJ = new Int_t [partleadjet] ;
  Int_t counterpartleadjet = 0;

  if (fDebug) cout << "Filling an array with indexes of leading jet particles" << endl;

  for( Int_t i = 0 ; i < fNPart ; i++ )
    {
      if ( fVectParticle[i]->njet == leadingjetindex )
	{  idxpartLJ[counterpartleadjet++] = i ; }
    }

  if ( (counterpartleadjet-1) > partleadjet ) { cout << " Counter_part_lead_jet > part_leadjet !!!!" << endl;}


  //___________________________________________________________________________
  // Calculus of part distribution in leading jet
  //___________________________________________________________________________
  Double_t z = 0. ;
  Double_t *zpartljet = new Double_t [ partleadjet ] ; // array of z of particles in leading jet

  if (fDebug) cout << "Entering loop of calculus of part distribution in leading jet" << endl ;

  for( Int_t j = 0 ; j < partleadjet ; j++ )
    {
      Double_t zj = fVectParticle[idxpartLJ[j]]->pt ;
      z =  zj / ptleadjet ;
      zpartljet [j] = z ;
      if (fDebug) cout << "idx_leadjet_part[j] = " << idxpartLJ[j]
		       << " p of particle = " << zj
		       << " pt lead jet = " << ptleadjet
		       << " Z = " << z << endl;
    }


  //___________________________________________________________________________
  // array of delta phi's between phi of particles and leading jet phi
  //___________________________________________________________________________
  if (fDebug) cout << "array of delta phi's between phi of particles and leading jet phi" << endl;
  Double_t dphipartLJ = 0. ;
  Double_t *dphipartljet = new Double_t [fNPart];
  for(  Int_t part = 0 ; part < fNPart ; part++ )
    {
      dphipartLJ = fVectParticle[part]->phi - phileadjet ;
      dphipartLJ = TVector2::Phi_mpi_pi (dphipartLJ) ; // restrict the delta phi to (-pi,pi) interval
      dphipartljet [part] = dphipartLJ ;
      if (fDebug) printf("part= %d ; dphi_partLJ = %g  \n", part, dphipartLJ );
    }


  //______________________________________________________________________________
  //  Pt distribution for all particles
  //______________________________________________________________________________
  TH1F * hpt = (TH1F*)fHistos->FindObject("histo11");
  if ( hpt ) { for (  Int_t i = 0 ; i < fNPart ; i++ ) { hpt->Fill( fVectParticle[i]->pt ); } }

  //___________________________________________________________________________
  // Recomputing of radius of particles in leading jet
  //___________________________________________________________________________
  if (fDebug) { printf("   Searching particles with jet index %d\n", leadingjetindex); }

  Double_t ddeta = 0. , ddphi = 0. , rpart = 0. ;

  for( Int_t j = 0 ; j < partleadjet ; j++ )
    {
      ddeta = etaleadjet - fVectParticle[idxpartLJ[j]]->eta;

      Double_t phitmp = fVectParticle[idxpartLJ[j]]->phi ;

      ddphi = TVector2::Phi_mpi_pi ( phileadjet - phitmp ) ; // restrict the delta phi to (-pi,pi) interval

      rpart = TMath::Hypot (ddeta, ddphi) ;

      if (fDebug) printf ("Particle %d with Re-Computed radius = %f ", idxpartLJ[j], rpart) ;
      if ( (rpart - fRadius) >= 0.00000001 )
	{ if (fDebug) printf ("    bigger than selected radius of %f\n", fRadius ); }
      else
	{ if (fDebug) printf ("\n") ; }

      if (hR) hR->Fill(rpart);

    }



  //_______________________________________________________________________
  // Computing of radius that contain 80% of Leading Jet ( PT and multiplicity )
  //_______________________________________________________________________
  Double_t corepartleadjet = 0.8 * partleadjet ;
  Double_t coreptleadjet = 0.8 * ptleadjet ;
  Int_t countercorepart = 0 ;
  Double_t countercorept = 0. ;
  Int_t sortedindex = -1 ;

  TProfile * hprof24 = (TProfile*)fHistos->FindObject("histo24");
  TProfile * hprof25 = (TProfile*)fHistos->FindObject("histo25");
  TProfile * hprof26 = (TProfile*)fHistos->FindObject("histo26");
  TProfile * hprof27 = (TProfile*)fHistos->FindObject("histo27");
  TProfile * hprof28 = (TProfile*)fHistos->FindObject("histo28");
  TProfile * hprof29 = (TProfile*)fHistos->FindObject("histo29");


  if ((hprof24) && (hprof25) && (hprof26) && (hprof27) && (hprof28) && (hprof29) )
    {
      for(  Int_t part = 0 ; part < fNPart ; part++ )
	{
	  Double_t pttmp = 0. ; Double_t etatmp = 0. ; Double_t phitmp = 0. ; // temporary variables
	  Double_t dpart = 0. ;
	  sortedindex = fIdxArray[part] ;

	  if ( fVectParticle [ sortedindex ]->njet == leadingjetindex )
	    {
	      pttmp = fVectParticle[sortedindex]->pt ;
	      etatmp = fVectParticle[sortedindex]->eta ;
	      phitmp = fVectParticle[sortedindex]->phi ;

	      ++countercorepart ;
	      countercorept += pttmp ;

	      dpart = TMath::Hypot ( etaleadjet - etatmp, TVector2::Phi_mpi_pi (phileadjet - phitmp) ) ;

	      if ( countercorepart <=  corepartleadjet ) { hprof24->Fill(ptleadjet, dpart); }
	      if ( countercorept <= coreptleadjet ) { hprof25->Fill(ptleadjet, dpart); }

	      if (ptleadjet >  5.) { hprof26->Fill(dpart, countercorepart); hprof28->Fill(dpart, countercorept); }
	      if (ptleadjet > 30.) { hprof27->Fill(dpart, countercorepart); hprof29->Fill(dpart, countercorept); }

	    }
	}
    }

  TH1F *hjetpt = (TH1F*)fHistos->FindObject("histo1");
  TH1F *hjeteta = (TH1F*)fHistos->FindObject("histo2");
  TH1F *hjetphi = (TH1F*)fHistos->FindObject("histo3");
  TH1F *hjetnjet = (TH1F*)fHistos->FindObject("histo4");

  for(  Int_t jet = 0 ; jet < fNJets ; jet++ )
    {
      if (hjetpt)   hjetpt   ->Fill ( fVectJet[jet]->pt   ) ;
      if (hjeteta)  hjeteta  ->Fill ( fVectJet[jet]->eta  ) ;
      if (hjetphi)  hjetphi  ->Fill ( fVectJet[jet]->phi  ) ;
      if (hjetnjet) hjetnjet ->Fill ( fVectJet[jet]->njet ) ;
    }

  TH1F *hjets = (TH1F*)fHistos->FindObject("histo5");
  if (hjets) hjets->Fill(fNJets);

  TH1F *hleadpart = (TH1F*)fHistos->FindObject("histo6");
  if (hleadpart) hleadpart->Fill(partleadjet);

  TProfile * hprof = (TProfile*)fHistos->FindObject("histo7");
  if (hprof) hprof->Fill(ptleadjet,partleadjet);

  TH1F *hMD = (TH1F*)fHistos->FindObject("histo8");
  for(  Int_t k = 0  ; k < partleadjet ; k++)
    { hMD->Fill( zpartljet[k] ); }

  TProfile * hphi = (TProfile*)fHistos->FindObject("histo9");
  for(  Int_t k = 0  ; k < partleadjet ; k++)
    { hphi->Fill( TMath::RadToDeg() * dphipartljet [k] , fNPart ) ; }

  TProfile * htpd = (TProfile*)fHistos->FindObject("histo10");
  for(  Int_t k = 0  ; k < fNPart ; k++)
    { htpd->Fill( TMath::RadToDeg() * dphipartljet [k] , ptsumevent ) ; }


  TProfile * hprof1 = (TProfile*)fHistos->FindObject("histo21");
  if (hprof1) hprof1->Fill(ptleadjet, fNPart);

  TProfile * hprof2 = (TProfile*)fHistos->FindObject("histo22");
  if (hprof2) hprof2->Fill(ptleadjet, ptsumevent);

  TProfile * hprof1toward = (TProfile*)fHistos->FindObject("histo21_toward");
  TProfile * hprof1transverse = (TProfile*)fHistos->FindObject("histo21_transverse");
  TProfile * hprof1away = (TProfile*)fHistos->FindObject("histo21_away");
  TProfile * hprof2toward = (TProfile*)fHistos->FindObject("histo22_toward");
  TProfile * hprof2transverse = (TProfile*)fHistos->FindObject("histo22_transverse");
  TProfile * hprof2away = (TProfile*)fHistos->FindObject("histo22_away");
  TH1F * hpttoward = (TH1F*)fHistos->FindObject("histo23_toward");
  TH1F * hpttransverse = (TH1F*)fHistos->FindObject("histo23_transverse");
  TH1F * hptaway = (TH1F*)fHistos->FindObject("histo23_away");

  if ( (hprof1toward) && (hprof1transverse) && (hprof1away) && (hprof2toward) && (hprof2transverse) && (hprof2away) )
    {
      for( Int_t part = 0  ; part < fNPart ; part++)
	{
	  Double_t ptpart = fVectParticle[part]->pt ; // pt of particle
	  if ( ( dphipartljet[part] >=0.) && ( dphipartljet[part] < kPI/3. ) )
	    {
	      hprof1toward->Fill( ptleadjet, fNPart );
	      hprof2toward->Fill( ptleadjet, ptsumevent);
	      hpttoward->Fill( ptpart );
	    }
	  else
	    if ( ( dphipartljet[part] >= (kPI/3.)) && ( dphipartljet[part] < (2.*kPI/3.)) )
	      {
		hprof1transverse->Fill( ptleadjet, fNPart );
		hprof2transverse->Fill( ptleadjet, ptsumevent);
		hpttransverse->Fill( ptpart );
	      }
	    else
	      if ( ( dphipartljet[part] >= ( 2.*kPI/3.)) && ( dphipartljet[part] < kPI ) )
		{
		  hprof1away->Fill( ptleadjet, fNPart );
		  hprof2away->Fill( ptleadjet, ptsumevent);
		  hptaway->Fill( ptpart );
		}
	}
    }

  delete [] dphipartljet;
  delete [] zpartljet;
  delete [] idxpartLJ;

}

//______________________________________________________________________________
void AliCdfJetFinder::Clean()
{
  // CLEANING SECTION
  for (  Int_t i = 0 ; i < fNPart ; i++ ){
    delete fVectParticle[i];
    fVectParticle[i] = 0;
  }
  delete [] fVectParticle;fVectParticle = 0;

  for (  Int_t i = 0 ; i < fNJets ; i++ ){
    delete fVectJet[i];
    fVectJet[i] = 0;
  } 
  delete [] fVectJet;fVectJet = 0;

  delete [] fPtArray;fPtArray = 0;
  delete [] fIdxArray;fIdxArray = 0;

  Reset();

}

