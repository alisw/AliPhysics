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

//---------------------------------------------------------------------
// Jet Finder based on CDF algortihm
// Charged jet evolution and the underlying event in proton-antiproton collisions at 1.8 TeV
// Physical Review D, vol. 65, Issue 9, id. 092002
// http://www.phys.ufl.edu/~rfield/cdf/chgjet/chgjet_intro.html
// Authors : Adrian.Sevcenco@cern.ch (adriansev@spacescience.ro )
//           Daniel.Felea@cern.ch    (dfelea@spacescience.ro)
//           Ciprian.Mihai.Mitu@cern.ch  (mcm@spacescience.ro)
//---------------------------------------------------------------------

/*
Changelog



*/

#include <vector>
#include <iostream>
#include <Riostream.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TBits.h>
#include <TArrayF.h>
#include "AliCdfJetFinder.h"
#include "AliCdfJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliJet.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "TProfile.h"


ClassImp ( AliCdfJetFinder )

//______________________________________________________________________________
AliCdfJetFinder::AliCdfJetFinder():
    AliJetFinder(),
    fHistos(0),
    fDebug(0),
    fFromAod(0),
    fAODwrite(0),
    fAODtracksWrite(0),
    fRefArr (NULL),
    fNJets(-9999),
    fNPart(-9999),
    fRadius(0.7),
    fMinJetParticles(1),
    fJetPtCut(0.),
    fVectParticle(NULL),
    fVectJet(NULL),
    fPtArray(NULL),
    fIdxArray(NULL)
  {  /* Constructor */  }

//______________________________________________________________________________
AliCdfJetFinder::~AliCdfJetFinder()

  {
  // destructor
  }

//______________________________________________________________________________
void AliCdfJetFinder::CreateOutputObjects(TList *histos)
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

  TProfile * h21_toward = new TProfile ("histo21_toward","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21_toward->SetStats(kTRUE);
  h21_toward->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21_toward->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21_toward->GetXaxis()->SetTitleColor(1);
  h21_toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21_toward);

  TProfile * h21_transverse = new TProfile ("histo21_transverse","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21_transverse->SetStats(kTRUE);
  h21_transverse->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21_transverse->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21_transverse->GetXaxis()->SetTitleColor(1);
  h21_transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21_transverse);

  TProfile * h21_away = new TProfile ("histo21_away","N_{chg}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0., 50. , 0., 12. ) ;
  h21_away->SetStats(kTRUE);
  h21_away->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h21_away->GetYaxis()->SetTitle("<N_{chg}(in the event - including jet1)> in 1 GeV/c bin");
  h21_away->GetXaxis()->SetTitleColor(1);
  h21_away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h21_away);

  TProfile * h22_toward = new TProfile ("histo22_toward","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22_toward->SetStats(kTRUE);
  h22_toward->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22_toward->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22_toward->GetXaxis()->SetTitleColor(1);
  h22_toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22_toward);

  TProfile * h22_transverse = new TProfile ("histo22_transverse","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22_transverse->SetStats(kTRUE);
  h22_transverse->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22_transverse->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22_transverse->GetXaxis()->SetTitleColor(1);
  h22_transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22_transverse);

  TProfile * h22_away = new TProfile ("histo22_away","PT_{sum}(in the event - including jet1) vs P_{T}(charged jet1)", 200, 0. , 50. , 0., 50. ) ;
  h22_away->SetStats(kTRUE);
  h22_away->GetXaxis()->SetTitle("P_{T} (charged jet1) (GeV/c)");
  h22_away->GetYaxis()->SetTitle("<PT_{sum}(in the event - including jet1)> in 1 GeV/c bin");
  h22_away->GetXaxis()->SetTitleColor(1);
  h22_away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h22_away);

  TH1F *h23_toward = new TH1F ("histo23_toward","'Toward' Pt Distribution of charged particles", 200, 0., 14.);
  h23_toward->SetStats(kTRUE);
  h23_toward->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23_toward->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23_toward->GetXaxis()->SetTitleColor(1);
  h23_toward->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23_toward);

  TH1F *h23_transverse = new TH1F ("histo23_transverse","'Transverse' Pt Distribution of charged particles", 200, 0., 14.);
  h23_transverse->SetStats(kTRUE);
  h23_transverse->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23_transverse->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23_transverse->GetXaxis()->SetTitleColor(1);
  h23_transverse->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23_transverse);

  TH1F *h23_away = new TH1F ("histo23_away","'Away' Pt Distribution of charged particles", 200, 0., 14.);
  h23_away->SetStats(kTRUE);
  h23_away->GetXaxis()->SetTitle("P_{T} (charged) (GeV/c)");
  h23_away->GetYaxis()->SetTitle("dN_{chg}/dP_{T} (1/GeV/c)");
  h23_away->GetXaxis()->SetTitleColor(1);
  h23_away->SetMarkerStyle(kFullCircle);
  fHistos->Add(h23_away);

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
/*
Jet Algorithm:
 * Order all charged particles according to their PT.
 * Start with the highest PT particle and include in the "jet" all particles within the "radius" R = 0.7
    (considering each particle in the order of decreasing PT and recalculating the centroid of the jet after
    each new particle is added to the jet).
 * Go to the next highest PT particle (not already included in a jet) and include in the "jet" all particles
    (not already included in a jet) within the radius R =0.7.
 * Continue until all particles are in a "jet".

*/
//______________________________________________________________________________

//______________________________________________________________________________
void AliCdfJetFinder::FindJets()
{
AliCdfJetHeader *header = (AliCdfJetHeader*)fHeader;

  if (header)
    {
    fDebug            = header->IsDebugCDF();
    fAODwrite         = header->IsAODwrite() ;       // write jets to AOD
    fAODtracksWrite   = header->IsAODtracksWrite() ; // write jet tracks to AOD
    fRadius           = header->GetRadius();      // get Radius from jet finder header
    fMinJetParticles  = header->GetMinPartJet (); // get minimum multiplicity of an jet
    fJetPtCut         = header->GetJetPtCut ();   // get minimum of jet pt
    }
  else
    { cout << "Header not found" << endl; return; }

// temporary until the other problems are resolved
fAODwrite = 0 ;

if (fAODwrite)
  {
  fFromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fFromAod) { fRefArr = fReader->GetReferences(); }
  }

fFromAod = 0 ; // disable for the moment ; only ESD reading for now

InitData();

if (!fNPart) { if (fDebug) {cout << "entries = 0 ; Event empty !!!" << endl ;} return; } // if event empty then exit

FindCones();

ComputeConesWeight();

if (fAODwrite) { cout << "Writing AOD" << endl ; WriteJets(); }

AnalizeJets();

Clean();

}

//______________________________________________________________________________
void AliCdfJetFinder::InitData()
{

  TClonesArray * vectArray = fReader->GetMomentumArray() ;
    if ( vectArray == 0 ) { cout << "Could not get the momentum array" << endl; return; }

  fNPart = vectArray->GetEntries()  ; // n particles in this event;

  if ( fNPart == 0 ) { return ; } // if event empty then exit

  fJets->SetNinput ( fNPart ) ; // number of input objects

  fVectParticle = new varContainer* [fNPart]; // container for Particles
  fVectJet      = new varContainer* [fNPart]; // container for Jets

  fPtArray  = new Double_t [fNPart] ; // momentum array
  fIdxArray = new Int_t    [fNPart] ; // index array of sorted pts

  // initialisation of momentum and index arrays
  for (  Int_t i = 0 ; i < fNPart ; i++ )
    {// SORTING STEP :: fPtArray with data from TClonesArray of TLorentzVector
    TLorentzVector * lv = (TLorentzVector*) vectArray->At(i);

   // INITIALISATION of local arrays for temporary storage
    varContainer *aParticle = new varContainer;
    aParticle->pt   = lv->Pt();
    aParticle->eta  = lv->Eta();
    aParticle->phi  = DeltaPhiNorm ( lv->Phi() ); // normalize to -pi,pi
    aParticle->njet = -999;

    fVectParticle[i] = aParticle;  // vector of Particles

    // initializing arrays
    fIdxArray [i] = -999 ;
     fPtArray [i] = aParticle->pt ;
    }

  TMath::Sort ( fNPart, fPtArray, fIdxArray ) ; // get a sorted array of indexes with TClonesArray.Size()

}


//______________________________________________________________________________
void AliCdfJetFinder::FindCones()
{
  Double_t  pt_seed = 0. , eta_seed = 0. , phi_seed = 0. ; // leading particle params
  Double_t pt_tmp = 0. , eta_tmp = 0. , phi_tmp = 0. ; // temporary variables to be used in various calculations
  Double_t deta = 0. , dphi = 0. , d_computed = 0. ;
  Bool_t in_jet = 0 ;

  fNJets = -1 ; // n jets in this event
  Int_t idxPtSort = -1 ;  // index of array of sorted pt indexes

  if (fDebug) { cout << "\n\n\n\n\n\n------------------\nBegin Event Analysis\n------------------\n\n" << endl ;}

  cout << "fNPart = " << fNPart << endl;

  TBits lkup_table ( fNPart ) ;  // bit container ; 1-to-1 corespondence with fIdxArray

  while ( lkup_table.CountBits() != (UInt_t)fNPart )
    { // loop over particles in event until all flags are set
    UInt_t first_non_flagged = lkup_table.FirstNullBit() ; // set the index to the first NON flagged bit ; less conditions

    cout << "\n\nfirst_non_flagged : " << first_non_flagged << endl;

    ++fNJets; // incrementing the jet counter
    if (fDebug) { printf("JET %d \n", fNJets); }

    pt_seed = 0. ; eta_seed = 0. ; phi_seed = 0. ;  // reseting leading particle params

    for (  UInt_t i_part = first_non_flagged ; i_part < (UInt_t)fNPart ; i_part++ )
      {// iteration over particles in event
      // the loop is done over sorted array of pt
      idxPtSort = fIdxArray[i_part] ;  // index of particle ! fIdxArray is an index list pt sorted

      if ( lkup_table.TestBitNumber(i_part) ) { continue; } // if 4vector is already flagged skip it

      //init computed and used vars
      pt_tmp = 0. ; eta_tmp = 0. ; phi_tmp = 0. ;
      deta = 0. ; dphi = 0. ; d_computed = 0. ; in_jet = 0 ;

      //taking info from fVectParticle ;
       pt_tmp = fVectParticle[idxPtSort]->pt ;
      eta_tmp = fVectParticle[idxPtSort]->eta ;
      phi_tmp = fVectParticle[idxPtSort]->phi ;

      if ( i_part == first_non_flagged )
        {// this is first particle in event; leading particle
        // begin the search around this particle in a fRadius

        // CENTRE OF THE JET
        pt_seed = pt_tmp ; eta_seed = eta_tmp ; phi_seed = phi_tmp ; // seeding the jet with first particle idxPtSort

        lkup_table.SetBitNumber ( i_part ) ; // flag the index of particle in lkup_table
        fVectParticle[idxPtSort]->njet = fNJets ; // associate particle with current jet number

        if (fDebug) { printf("\nLeading particle :: particle index = %d ;  at sorted index %d ; in jet %d \n", idxPtSort, i_part, fNJets); }
        if (fDebug) { printf("pt= %g ; eta= %g ; phi = %g \n", pt_tmp, eta_tmp, phi_tmp) ; }
        if (fDebug) { lkup_table.Print() ;}

        continue ; // skip to next particle
        }

      // condition to be in jet
      deta = eta_tmp - eta_seed ;
      dphi = DeltaPhiNorm ( phi_tmp - phi_seed ) ; // computing dphi and normalizing to (0,2pi) interval in one step

      d_computed = Distance (deta, dphi) ; // Distance(fRadius) to (eta,phi) seed

      in_jet = ( ( fRadius - d_computed ) >= 0.000000001 ) ? 1 : 0 ; // if r_computed is within jet_r in_jet == 1 else 0

      if ( in_jet )
        { // calculus of jet variables
        lkup_table.SetBitNumber ( i_part ) ;  // flag the index of particle in lkup_table
        fVectParticle[idxPtSort]->njet = fNJets ; // setting in particle list the associated jet

        if (fDebug) { printf("\njet particle :: particle index = %d ; at sorted index %d ; in jet %d ; found at radius %g ;  \n", idxPtSort, i_part, fNJets, d_computed); }
        if (fDebug) { printf("pt= %g ; eta= %g ; phi = %g \n", pt_tmp, eta_tmp, phi_tmp) ; }
        if (fDebug) { lkup_table.Print() ;}

        continue ; // skip to next particle
        }

      }
      // end of iteration over event; one jet definition of content ; jet parameters to be computed later
    }
}


//______________________________________________________________________________
void AliCdfJetFinder::ComputeConesWeight()
{
/*
      CALCULUS OF JETS Pt, Eta and Phi (centre of weight in (eta,phi) plane)
*/
  // rescan the vector of particles by identify them by asociate jet number for computing of weight centre
  // we know : fNJets = the number of jets

Double_t pt_jet , eta_jet , phi_jet ; Int_t npartJet ;
Double_t pt_tmp = 0. , eta_tmp = 0. , phi_tmp = 0. ; // temporary variables to be used in various calculations
Int_t idxPtSort = -999 ;  // index of array of sorted pt indexes

for(  Int_t jet = 0 ; jet < fNJets ; jet++ )
  {
  if (fDebug) { printf("\n\n--- Computing weight of Jet %d \n", jet ); }
  npartJet = 0 ; pt_jet = 0. ; eta_jet = 0. ; phi_jet = 0. ; // reset variables for a new computation

  for (  Int_t i_part = 0 ; i_part < fNPart ; i_part++ )
    {// iteration over particles in event
    // the loop is done over sorted array of pt
    idxPtSort = fIdxArray[i_part] ;  // index of particle ! fIdxArray is an index list pt sorted

    if ( fVectParticle[idxPtSort]->njet == jet )
      {
      ++npartJet; // incrementing the counter of jet particles

      //taking info from fVectParticle ;
       pt_tmp = fVectParticle[idxPtSort]->pt ;
      eta_tmp = fVectParticle[idxPtSort]->eta ;
      phi_tmp = DeltaPhiNorm ( fVectParticle[idxPtSort]->phi ) ;

      pt_jet += pt_tmp ;
      eta_jet = ( (pt_jet*eta_jet) + (pt_tmp*eta_tmp) )/(pt_jet + pt_tmp) ;
      phi_jet = ( (pt_jet*phi_jet) + (pt_tmp*phi_tmp) )/(pt_jet + pt_tmp) ;

      }
      // add a particle and recalculation of centroid
    }
    // end of 1 jet computation

    varContainer *aJet = new varContainer;  // Jet container
    aJet->pt = pt_jet; aJet->eta = eta_jet; aJet->phi = phi_jet; aJet->njet = npartJet; // setting jet vars in container
    fVectJet[jet] = aJet;   // store the number of the jet(fNJets) and increment afterwards

    if (fDebug) { printf ("=== current jet %d : npartjet= %d ; pt_jet= %g ; eta_jet = %g ; phi_jet = %g \n\n\n",
                                       jet,     npartJet,      pt_jet,      eta_jet,       phi_jet ) ; }

  }
  //end loop over jets

}


//______________________________________________________________________________
void AliCdfJetFinder::WriteJets()
{ // Writing AOD jets and AOD tracks

/*  for(  Int_t jet_nr = 0 ; jet_nr < fNJets ; jet_nr++ )
    {
    Double_t pt = 0., eta = 0., phi = 0., // jet variables
             px = 0., py = 0., pz = 0., en = 0.; // convert to 4-vector
    pt  = fVectJet[ jet_nr ]->pt   ; // pt  of jet
    eta = fVectJet[ jet_nr ]->eta  ; // eta of jet
    phi = fVectJet[ jet_nr ]->phi  ; // phi of jet

    px = pt * TMath::Cos ( phi ) ;
    py = pt * TMath::Sin ( phi ) ;
    pz = pt / TMath::Tan ( 2.0 * TMath::ATan ( TMath::Exp ( -eta ) ) ) ;
    en = TMath::Sqrt ( px * px + py * py + pz * pz );

    AliAODJet jet (px, py, pz, en);
    AddJet(jet);

    if (fDebug) jet.Print("");

    if (fromAod && AODtracksWrite)
      {
      for (  Int_t jet_track = 0; jet_track < fNPart; jet_track++ )
        { if ( fVectParticle[jet_track]->njet == jet_nr ) { jet.AddTrack(refs->At(jet_track)) ; } }
      }
      // tracks REFs written in AOD

    }*/

//jets vector parsed and written to AOD

}


//______________________________________________________________________________
void AliCdfJetFinder::AnalizeJets()
  {

  //persistent pointer to histo20
  TH1F *h_r = (TH1F*)fHistos->FindObject("histo20");

  Int_t   *jets_pt_idx = 0;     // sorted array of jets pt
  Double_t    *jets_pt = 0;     // array of jets pts
  Int_t leading_jet_index = -1 ;   // index of leading jet from fVectJet
  Int_t part_leadjet = 0 ; // number of particles in leading jet
  Double_t   pt_leadjet = 0. ; // pt  of leading jet
  Double_t  eta_leadjet = 0. ; // eta of leading jet
  Double_t  phi_leadjet = 0. ; // phi of leading jet

  jets_pt_idx = new Int_t    [fNJets] ;
  jets_pt     = new Double_t [fNJets] ;

//________________________________________________________________________________________
//  Jet sorting and finding the leading jet that coresponds to cuts in pt and multiplicity
//________________________________________________________________________________________

  // filing the idx_ptjets array
  if (fDebug) printf("List of unsorted jets:\n");
  for(  Int_t i = 0 ; i < fNJets ; i++ )
    {
    jets_pt_idx [i] = 0 ;
    jets_pt     [i] = fVectJet[i]->pt ;
    if (fDebug) { cout << "   jet found: " << i << " npartjet=" << fVectJet[i]->njet << " ; jets_pt = " << jets_pt[i] << endl; }
    }

  TMath::Sort ( fNJets, jets_pt , jets_pt_idx ) ; // sorting pt of jets

  // selection of leading jet
  // looping over jets searching for __first__ one that coresponds to cuts
  for(  Int_t i = 0 ; i < fNJets ; i++ )
    {
    if ( ( fVectJet[ jets_pt_idx[i] ]->njet >= fMinJetParticles ) && ( fVectJet[ jets_pt_idx[i] ]->pt >= fJetPtCut ) )
      {
      leading_jet_index = jets_pt_idx[i] ;
      part_leadjet = fVectJet[ leading_jet_index ]->njet ; // number of particles in leading jet
        pt_leadjet = fVectJet[ leading_jet_index ]->pt   ; // pt  of leading jet
       eta_leadjet = fVectJet[ leading_jet_index ]->eta  ; // eta of leading jet
       phi_leadjet = fVectJet[ leading_jet_index ]->phi  ; // phi of leading jet

      if (fDebug)
      { printf("Leading jet %d : npart= %d ; pt= %g ; eta = %g ; phi = %g \n", leading_jet_index, part_leadjet, pt_leadjet, eta_leadjet, phi_leadjet ); }

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
cout << "Computing sum of pt in event" << endl ;
Double_t pt_sum_event = 0.;
for (  Int_t i = 0 ; i< fNPart ; i++ ) { pt_sum_event += fVectParticle[i]->pt ; }
printf ("Sum of all Pt in event : pt_sum_event = %g", pt_sum_event) ;

//___________________________________________________________________________
// Filling an array with indexes of leading jet particles
//___________________________________________________________________________
Int_t * idx_leadjet_part = new Int_t [part_leadjet] ;
Int_t counter_part_lead_jet = 0;

cout << "Filling an array with indexes of leading jet particles" << endl;

for( Int_t i = 0 ; i < fNPart ; i++ )
  {
  if ( fVectParticle[i]->njet == leading_jet_index )
    {  idx_leadjet_part[counter_part_lead_jet++] = i ; }
  }

if ( (counter_part_lead_jet-1) > part_leadjet ) { cout << " Counter_part_lead_jet > part_leadjet !!!!" << endl;}


//___________________________________________________________________________
// Calculus of part distribution in leading jet
//___________________________________________________________________________
Double_t z = 0. ;
Double_t *z_part_ljet = new Double_t [ part_leadjet ] ; // array of z of particles in leading jet

cout << "Entering loop of calculus of part distribution in leading jet" << endl ;

for( Int_t j = 0 ; j < part_leadjet ; j++ )
  {
  Double_t z_j = fVectParticle[idx_leadjet_part[j]]->pt ;
  z =  z_j / pt_leadjet ;
  z_part_ljet [j] = z ;
  cout << "idx_leadjet_part[j] = " << idx_leadjet_part[j]
      << " p of particle = " << z_j
      << " pt lead jet = " << pt_leadjet
      << " Z = " << z << endl;
  }


//___________________________________________________________________________
// array of delta phi's between phi of particles and leading jet phi
//___________________________________________________________________________
cout << "array of delta phi's between phi of particles and leading jet phi" << endl;
Double_t dphi_partLJ = 0. ;
Double_t *dphi_part_ljet = new Double_t [fNPart];
for(  Int_t part = 0 ; part < fNPart ; part++ )
  {
  dphi_partLJ = fVectParticle[part]->phi - phi_leadjet ;
  dphi_partLJ = DeltaPhiNorm (dphi_partLJ) ; // restrict the delta phi to (0,pi) interval
  dphi_part_ljet [part] = dphi_partLJ ;
  printf("part= %d ; dphi_partLJ = %g  \n", part, dphi_partLJ );
  }


//______________________________________________________________________________
//  Pt distribution for all particles
//______________________________________________________________________________
TH1F * h_pt = (TH1F*)fHistos->FindObject("histo11");
if ( h_pt ) { for (  Int_t i = 0 ; i < fNPart ; i++ ) { h_pt->Fill( fVectParticle[i]->pt ); } }

//___________________________________________________________________________
// Recomputing of radius of particles in leading jet
//___________________________________________________________________________
if (fDebug) { printf("   Searching particles with jet index %d\n", leading_jet_index); }

Double_t ddeta = 0. , ddphi = 0. , r_part = 0. ;

for( Int_t j = 0 ; j < part_leadjet ; j++ )
  {
  ddeta = eta_leadjet - fVectParticle[idx_leadjet_part[j]]->eta;

  Double_t phi_tmp = fVectParticle[idx_leadjet_part[j]]->phi ;
  phi_tmp = DeltaPhiNorm (phi_tmp);

  ddphi = DeltaPhiNorm ( phi_leadjet - phi_tmp ) ; // restrict the delta phi to (-pi,pi) interval

  r_part = Distance (ddeta, ddphi) ;

  printf ("Particle %d with Re-Computed radius = %f ", idx_leadjet_part[j], r_part) ;
  if ( (r_part - fRadius) >= 0.00000001 )
    { printf ("    bigger than selected radius of %f\n", fRadius ); }
  else
    { printf ("\n") ; }

  if (h_r) h_r->Fill(r_part);

  }



//_______________________________________________________________________
// Computing of radius that contain 80% of Leading Jet ( PT and multiplicity )
//_______________________________________________________________________
Double_t core_part_leadjet = 0.8 * part_leadjet ;
Double_t core_pt_leadjet = 0.8 * pt_leadjet ;
Int_t counter_core_part = 0 ;
Double_t counter_core_pt = 0. ;
Int_t sorted_index = -1 ;

TProfile * h_prof_24 = (TProfile*)fHistos->FindObject("histo24");
TProfile * h_prof_25 = (TProfile*)fHistos->FindObject("histo25");

TProfile * h_prof_26 = (TProfile*)fHistos->FindObject("histo26");
TProfile * h_prof_27 = (TProfile*)fHistos->FindObject("histo27");
TProfile * h_prof_28 = (TProfile*)fHistos->FindObject("histo28");
TProfile * h_prof_29 = (TProfile*)fHistos->FindObject("histo29");


if ((h_prof_24) && (h_prof_25) && (h_prof_26) && (h_prof_27) && (h_prof_28) && (h_prof_29) )
{
for(  Int_t part = 0 ; part < fNPart ; part++ )
  {
  Double_t pt_tmp = 0. ; Double_t eta_tmp = 0. ; Double_t phi_tmp = 0. ; // temporary variables
  Double_t d_part = 0. ;
  sorted_index = fIdxArray[part] ;

  if ( fVectParticle [ sorted_index ]->njet == leading_jet_index )
    {
     pt_tmp = fVectParticle[sorted_index]->pt ;
    eta_tmp = fVectParticle[sorted_index]->eta ;
    phi_tmp = fVectParticle[sorted_index]->phi ;

    ++counter_core_part ;
    counter_core_pt += pt_tmp ;

    d_part = Distance ( eta_leadjet - eta_tmp, phi_leadjet - phi_tmp ) ;

    if ( counter_core_part <=  core_part_leadjet ) { h_prof_24->Fill(pt_leadjet, d_part); }
    if ( counter_core_pt <= core_pt_leadjet ) { h_prof_25->Fill(pt_leadjet, d_part); }

    if (pt_leadjet >  5.) { h_prof_26->Fill(d_part, counter_core_part); h_prof_28->Fill(d_part, counter_core_pt); }
    if (pt_leadjet > 30.) { h_prof_27->Fill(d_part, counter_core_part); h_prof_29->Fill(d_part, counter_core_pt); }

    }
  }
}










  TH1F *h_jet_pt = (TH1F*)fHistos->FindObject("histo1");
  TH1F *h_jet_eta = (TH1F*)fHistos->FindObject("histo2");
  TH1F *h_jet_phi = (TH1F*)fHistos->FindObject("histo3");
  TH1F *h_jet_njet = (TH1F*)fHistos->FindObject("histo4");

  for(  Int_t jet = 0 ; jet < fNJets ; jet++ )
    {
    if (h_jet_pt)   h_jet_pt   ->Fill ( fVectJet[ jet ]->pt   ) ;
    if (h_jet_eta)  h_jet_eta  ->Fill ( fVectJet[ jet ]->eta  ) ;
    if (h_jet_phi)  h_jet_phi  ->Fill ( fVectJet[ jet ]->phi  ) ;
    if (h_jet_njet) h_jet_njet ->Fill ( fVectJet[ jet ]->njet ) ;
    }

  TH1F *h_jets = (TH1F*)fHistos->FindObject("histo5");
  if (h_jets) h_jets->Fill(fNJets);

  TH1F *h_leadpart = (TH1F*)fHistos->FindObject("histo6");
  if (h_leadpart) h_leadpart->Fill(part_leadjet);

  TProfile * h_prof = (TProfile*)fHistos->FindObject("histo7");
  if (h_prof) h_prof->Fill(pt_leadjet,part_leadjet);

  TH1F *h_MD = (TH1F*)fHistos->FindObject("histo8");
   for(  Int_t k = 0  ; k < part_leadjet ; k++)
     { h_MD->Fill( z_part_ljet[k] ); }

  TProfile * h_phi = (TProfile*)fHistos->FindObject("histo9");
    for(  Int_t k = 0  ; k < part_leadjet ; k++)
        { h_phi->Fill( TMath::RadToDeg() * dphi_part_ljet [k] , fNPart ) ; }

  TProfile * h_tpd = (TProfile*)fHistos->FindObject("histo10");
    for(  Int_t k = 0  ; k < fNPart ; k++)
        { h_tpd->Fill( TMath::RadToDeg() * dphi_part_ljet [k] , pt_sum_event ) ; }


  TProfile * h_prof1 = (TProfile*)fHistos->FindObject("histo21");
  if (h_prof1) h_prof1->Fill(pt_leadjet, fNPart);

  TProfile * h_prof2 = (TProfile*)fHistos->FindObject("histo22");
  if (h_prof2) h_prof2->Fill(pt_leadjet, pt_sum_event);

  TProfile * h_prof1_toward = (TProfile*)fHistos->FindObject("histo21_toward");
  TProfile * h_prof1_transverse = (TProfile*)fHistos->FindObject("histo21_transverse");
  TProfile * h_prof1_away = (TProfile*)fHistos->FindObject("histo21_away");
  TProfile * h_prof2_toward = (TProfile*)fHistos->FindObject("histo22_toward");
  TProfile * h_prof2_transverse = (TProfile*)fHistos->FindObject("histo22_transverse");
  TProfile * h_prof2_away = (TProfile*)fHistos->FindObject("histo22_away");
  TH1F * h_pt_toward = (TH1F*)fHistos->FindObject("histo23_toward");
  TH1F * h_pt_transverse = (TH1F*)fHistos->FindObject("histo23_transverse");
  TH1F * h_pt_away = (TH1F*)fHistos->FindObject("histo23_away");



  if ( (h_prof1_toward) && (h_prof1_transverse) && (h_prof1_away) && (h_prof2_toward) && (h_prof2_transverse) && (h_prof2_away) )
    {
    for( Int_t part = 0  ; part < fNPart ; part++)
      {
      Double_t pt_part = fVectParticle[part]->pt ; // pt of particle
      if ( ( dphi_part_ljet[part] >=0.) && ( dphi_part_ljet[part] < pi/3. ) )
        {
        h_prof1_toward->Fill( pt_leadjet, fNPart );
        h_prof2_toward->Fill( pt_leadjet, pt_sum_event);
        h_pt_toward->Fill( pt_part );
        }
      else
      if ( ( dphi_part_ljet[part] >= (pi/3.)) && ( dphi_part_ljet[part] < (2.*pi/3.)) )
        {
        h_prof1_transverse->Fill( pt_leadjet, fNPart );
        h_prof2_transverse->Fill( pt_leadjet, pt_sum_event);
        h_pt_transverse->Fill( pt_part );
        }
      else
      if ( ( dphi_part_ljet[part] >= ( 2.*pi/3.)) && ( dphi_part_ljet[part] < pi ) )
        {
        h_prof1_away->Fill( pt_leadjet, fNPart );
        h_prof2_away->Fill( pt_leadjet, pt_sum_event);
        h_pt_away->Fill( pt_part );
        }
      }
    }







}


//______________________________________________________________________________
void AliCdfJetFinder::Clean()
{


  // CLEANING SECTION
  delete [] fVectParticle;
  delete [] fVectJet;
  delete [] fPtArray;
  delete [] fIdxArray;

//  if (z_part_ljet) delete [] z_part_ljet;
//  if (idx_leadjet_part) delete [] idx_leadjet_part;
//  if (jets_pt_idx) delete [] jets_pt_idx ;
//  if (jets_pt) delete [] jets_pt ;

//  vectArray->Delete() ; // TClonesArray of lorentz vectors

  Reset();



}




//______________________________________________________________________________
void AliCdfJetFinder::FinishRun()
{}


