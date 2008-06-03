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
//
//
//---------------------------------------------------------------------
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

using namespace std ;

struct varContainer // container for Particle variables
  { // variables of container struct
  Double_t  pt; Double_t eta; Double_t phi;
  Int_t njet; // if jets are stored in varContainer nr_jet is multiplicity of jet
  } ;

ClassImp ( AliCdfJetFinder )

//______________________________________________________________________________
AliCdfJetFinder::AliCdfJetFinder():
   AliJetFinder(),
   fHistos(0)
  {  /* Constructor */  }

//______________________________________________________________________________
AliCdfJetFinder::~AliCdfJetFinder()

  {
  // destructor
  cout << "Calling the destructor ... " << endl ;
  Reset();
  cout << "Destructor called!" << endl;
  }

//______________________________________________________________________________
void AliCdfJetFinder::CreateOutputObjects(TList *histos)
{
  // Create the list of histograms. Only the list is owned.
  fHistos = histos;

  TH1F *h1 = new TH1F ("histo1", "no jets/event", 100, 0,20);
  h1->SetStats(kTRUE);
  h1->GetXaxis()->SetTitle("N_{jets}");
  h1->GetYaxis()->SetTitle("#frac{dN}{dN_{jets}}");
  h1->GetXaxis()->SetTitleColor(1);
  h1->SetMarkerStyle(kFullCircle);
  fHistos->Add(h1);

  TH1F *h2 = new TH1F ("histo2", "no part/leading jet", 40, 0,40);
  h2->SetStats(kTRUE);
  h2->GetXaxis()->SetTitle("N_{leading}");
  h2->GetYaxis()->SetTitle("#frac{dN}{dN_{leading}}");
  h2->GetXaxis()->SetTitleColor(1);
  h2->SetMarkerStyle(kFullCircle);
  fHistos->Add(h2);

  TProfile * h3 = new TProfile ("histo3","ProfileX of (pt,npart) of leading jets", 25, 0. ,25. , 0.,25. ) ;
  h3->SetStats(kTRUE);
  h3->GetXaxis()->SetTitle("P_{t}");
  h3->GetYaxis()->SetTitle("N_{leading}");
  h3->GetXaxis()->SetTitleColor(1);
  h3->SetMarkerStyle(kFullCircle);
  fHistos->Add(h3);

  TH1F *h4 = new TH1F ("histo4", "Charge momentum distribution for leading jet", 500, 0 , 5);
  fHistos->Add(h4);

  TProfile *h5 = new TProfile ("histo5", "ProfileX of N_{charge} vs dphi leading jet", 100 , 0. , 200. , 0 , 100 );
  fHistos->Add(h5);

  TH1F *h6 = new TH1F ("histo6", " \"Transverse\" Pt Distribution ", 70, 0 , 14);
  fHistos->Add(h6);


//cout << "CreateOutputObjetcs done!" << endl ;
}

//______________________________________________________________________________
void AliCdfJetFinder::FindJets()
  {
  //cout << "Begining of FindJets ..." << endl ;
  //1) Fill 1 array with momentum and initialise another array for indexes
  //2) Run algorithm
  //   2.1) Sort momentum array
  //   2.2) loop over arrays of TLorentzVectors
  //   2.3) flag as a possible jet
  //3) fill AliAODJet list
  Bool_t debug = kFALSE;

  TRefArray *refs = 0;
  Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fromAod) { refs = fReader->GetReferences(); }

  AliCdfJetHeader *header = (AliCdfJetHeader*)fHeader;
  Double_t radius = header->GetRadius(); // get Radius from header
//  cout << "Radius is : " << radius << endl ;

  TClonesArray * vectArray = fReader->GetMomentumArray() ;
  if ( vectArray == 0 ) { cout << "Could not get the momentum array" << endl; return; }

  Int_t nPart = vectArray->GetEntries()  ; // n particles in this event;
  if ( nPart == 0 ) { if (debug) cout << "entries = 0 ; Event empty !!!" << endl ; return; }

  Int_t nJets = 0; // n jets in this event

  fJets->SetNinput ( nPart ) ; // number of input objects

  varContainer **vectParticle = new varContainer* [nPart]; // container for Particles
  varContainer **vectJet      = new varContainer* [nPart]; // container for Jets

  Double_t *ptArray = new Double_t [nPart] ; // momentum array
  Int_t   *idxArray = new Int_t    [nPart] ; // index array of sorted pts

  // initialisation of momentum and index arrays
  //  cout << "Filling idxArray && momArray for sorting " << endl;
  for ( Int_t i = 0 ; i < nPart ; i++ )
    {// SORTING STEP :: ptArray with data from TClonesArray of TLorentzVector
    TLorentzVector * lv = (TLorentzVector*) vectArray->At(i);
    // INITIALISATION of local arrays for temporary storage
    varContainer *aParticle = new varContainer;
    aParticle->pt   = lv->Pt();
    aParticle->eta  = lv->Eta();
    aParticle->phi  = lv->Phi();
    aParticle->njet = -1;
    vectParticle[i] = aParticle;

    // initializing arrays
    idxArray [i] = 0 ;
     ptArray [i] = aParticle->pt ;
    }

  TMath::Sort ( nPart, ptArray, idxArray ) ; // get a sorted array of indexes with size = size of TClonesArray

  Double_t  pt_jet = 0. ,
           eta_jet = 0. ,
           phi_jet = 0. ,
           sum_eta = 0. ,
           sum_phi = 0. ;  // jet variables

  Int_t idxPtSort = 0 ,  // index of array of sorted pt indexes
         npartJet = 0 ;  // number of particles in curent jet

  TBits lkup_table ( nPart ) ;  // bit container of size npart

  while ( lkup_table.CountBits(0) != (UInt_t)nPart )
    { // loop over particles in event until all flags are set
    UInt_t first_non_flagged = lkup_table.FirstNullBit(0) ; // set the index to the first NON flagged bit ; less conditions

    npartJet = 0 ; // reseting number of particles in jet
    pt_jet   = 0.;
    eta_jet  = 0.;
    phi_jet  = 0.;
    sum_eta  = 0.;
    sum_phi  = 0.;

    for ( UInt_t i_part = first_non_flagged ; i_part < (UInt_t)nPart ; i_part++ )
      {// iteration over particles in event
      // reseting variables

      // the loop is done over sorted array of pt
      idxPtSort = idxArray[i_part] ;  // index entry of TLorentzVector from vect_arr

      if ( lkup_table.TestBitNumber(i_part) == 1 ) { continue; } // if 4vector is already flagged skip it

      //taking info from vectParticle ;
      Double_t  pt_tmp = vectParticle[idxPtSort]->pt ;
      Double_t eta_tmp = vectParticle[idxPtSort]->eta ;
      Double_t phi_tmp = vectParticle[idxPtSort]->phi ;
      if (debug) printf("   particle %d: pt=%g\n", i_part, pt_tmp);

      // all angles are expressed in rad

      if ( TMath::Abs(eta_tmp) > 0.9 )
         {
         lkup_table.SetBitNumber ( i_part, 1 ) ; // mark particle as used to be skipped
         continue ;
         }


      if ( i_part == first_non_flagged )
        {// this is first particle in event; leading particle
        // initialise the jet variables with leading particle numbers

        npartJet++; // incrementing counter of particles in jet (nparJet == 1) (leading particle)
        pt_jet = pt_tmp ;
        eta_jet = eta_tmp ;
        phi_jet = phi_tmp ;
        if (debug) printf("  first part in jet: npartjet=%d  pt_jet=%g\n", npartJet, pt_jet);

        lkup_table.SetBitNumber ( i_part, 1 ) ; // flag the 4vector

        vectParticle[idxPtSort]->njet = nJets ; // associate particle with current jet

        continue ; // skip to next particle
        }

      // condition to be in jet
      Double_t deta = eta_jet - eta_tmp ;
      Double_t dphi = phi_jet - phi_tmp ;

      if ( dphi < -TMath::Pi() ) { dphi = -dphi - 2.0 * TMath::Pi() ; }
      if ( dphi >  TMath::Pi() ) { dphi = -dphi + 2.0 * TMath::Pi() ; }

      Double_t r_computed = TMath::Sqrt ( deta * deta + dphi * dphi ) ;

      Bool_t in_jet = ( r_computed <= radius ) ? 1 : 0 ; // if r_computed is within jet_r in_jet == 1 else 0

      if ( in_jet )
        { // calculus of jet variables
        npartJet++;  // incrementing counter of particles in jet

        pt_jet += pt_tmp ;
        sum_eta += (pt_tmp * eta_tmp) ;
        sum_phi += (pt_tmp * phi_tmp) ;
        eta_jet = sum_eta / pt_jet ;
        phi_jet = sum_phi / pt_jet ;

        lkup_table.SetBitNumber ( i_part, 1 ) ;  // flag the 4vector
        vectParticle[idxPtSort]->njet = nJets ; // setting in particle list the associated jet
        if (debug) printf("   particle added to jet: npartjet=%d  pt_jet=%g\n", npartJet, pt_jet);
        continue ; // skip to next particle
        }

      }
      // end of iteration over event; one jet definition

    varContainer *aJet = new varContainer;
    aJet->pt = pt_jet;
    aJet->eta = eta_jet;
    aJet->phi = phi_jet;
    aJet->njet = npartJet;
    vectJet[nJets++] = aJet;   // store the jet (from jet 0 to nJets-1 )
    if (debug) printf("=== current jet: npartjet=%d  ptjet=%g\n", npartJet,pt_jet);

       // writing aod information
    Double_t px = 0., py = 0., pz = 0., en = 0.; // convert to 4-vector
    px = pt_jet * TMath::Cos ( phi_jet ) ;
    py = pt_jet * TMath::Sin ( phi_jet ) ;
    pz = pt_jet / TMath::Tan ( 2.0 * TMath::ATan ( TMath::Exp ( -eta_jet ) ) ) ;
    en = TMath::Sqrt ( px * px + py * py + pz * pz );

    if (npartJet<2) continue;         // do not add jets with less than 2 particles
    if (debug) cout << "Jet 4vect : " << "px = " << px << " ; py = " << py << " ; pz = " << pz << " ; E = " << en << endl;

    AliAODJet jet (px, py, pz, en);
//    cout << "Printing jet " << endl;
    if (debug) jet.Print("");

//    cout << "Adding jet ... " ;
    AddJet(jet);
//    cout << "added \n" << endl;

//     if (fromAod)
//       {
//       for (Int_t parts = 0; parts < nPart; parts++ )
//          { if (idx_jetT[parts] == nr_jet) {jet.AddTrack(refs->At(parts));} }
//       }


    }
    // end of while loop over particles ; ends when all particles were flagged as used and all jets defined

  /////////////////////////////////
  ////   END OF EVENT PARSING   ///
  /////////////////////////////////

  Int_t   *jets_pt_idx = 0;     // sorted array of jets pt
  Double_t    *jets_pt = 0;     // array of jets pts
  Int_t leading_jet_index = -1 ;   // index of leading jet from vectJet
  Int_t part_leadjet = 0 ;         // number of particles in leading jet
  Double_t   pt_leadjet = 0. ; // pt  of leading jet
  Double_t  eta_leadjet = 0. ; // eta of leading jet
  Double_t  phi_leadjet = 0. ; // phi of leading jet
  Int_t * idx_leadjet_part = 0;

  if (nJets > 0)
    { // if there is at least one jet in event
    jets_pt_idx = new Int_t    [nJets] ;
    jets_pt     = new Double_t [nJets] ;

    // filing the idx_ptjets array
    if (debug) printf("List of unsorted jets:\n");
    for( Int_t i = 0 ; i < nJets ; i++ )
      {
      jets_pt_idx [i] = 0 ;
      jets_pt [i] = vectJet[i]->pt ;
      if (debug) cout << "   jet found: " << i << " npartjet=" << vectJet[i]->njet << "  jets_pt[i]= " << jets_pt [i] << endl;
      }
    TMath::Sort ( nJets, jets_pt , jets_pt_idx ) ; // sorting pt of jets

    // selection of leading jet with nr of particles > 1
    for( Int_t i = 0 ; i < nJets ; i++ )
      {
      if ( vectJet[ jets_pt_idx[i] ]->njet > 1 )
        {
        leading_jet_index = jets_pt_idx[i] ;
        part_leadjet = vectJet[ leading_jet_index ]->njet ; // number of particles in leading jet
          pt_leadjet = vectJet[ leading_jet_index ]->pt   ; // pt  of leading jet
         eta_leadjet = vectJet[ leading_jet_index ]->eta  ; // eta of leading jet
         phi_leadjet = vectJet[ leading_jet_index ]->phi  ; // phi of leading jet

        cout << "Leading jet: npart = " << part_leadjet
             << " ; pt  = " << pt_leadjet
             << " ; phi = " << phi_leadjet
             << " ; eta = " << eta_leadjet
             << endl ;
        break ;
        }
      }
      // end of selection of leading jet

   // Filling an array with indexes of leading jet particles
   idx_leadjet_part = new Int_t [part_leadjet] ;
   Int_t counter = 0;
   if (debug) printf("   Searching particles with jet index %d\n", leading_jet_index);
   for( Int_t i = 0 ; i <nPart ; i++ )
     {
     if ( vectParticle[i]->njet == leading_jet_index )
       {
       idx_leadjet_part[counter++] = i ;
       if (debug) cout << "   " << counter-1 << ": index=" << i << "  pt=" <<  vectParticle[i]->pt << endl;
       }
     }
   if ( (counter-1) > part_leadjet ) { cout << " Counter > part_leadjet !!!!" << endl;}

   }
    // end of loop over DETECTED jets



// Calculus of part distribution in leading jet
Double_t z = 0. ;
Double_t *z_part_ljet = new Double_t [ part_leadjet ] ; // array of z of particles in leading jet
for( Int_t j = 0 ; j < part_leadjet ; j++ )
  {
  Double_t z_j = vectParticle[idx_leadjet_part[j]]->pt ;
  z =  z_j / pt_leadjet ;
  z_part_ljet [j] = z ;
  cout << "idx_leadjet_part[j] = " << idx_leadjet_part[j]
      << " p of particle = " << z_j
      << " pt lead jet = " << pt_leadjet
      << " Z = " << z << endl;
  }


/*
Double_t dphi_partLJ = 0. ;
Double_t dphiArray [nPart];
for( Int_t part = 0 ; part < nPart ; part++ )
  {
  dphi_partLJ = phi_leadjet - phi_partT [part] ;

  if ( dphi_partLJ < -TMath::Pi() ) { dphi_partLJ = -dphi_partLJ - 2.0 * TMath::Pi() ; }
  if ( dphi_partLJ >  TMath::Pi() ) { dphi_partLJ = -dphi_partLJ + 2.0 * TMath::Pi() ; }

  dphi_partLJ = TMath::Abs(dphi_partLJ) ;

  dphiArray [part] = dphi_partLJ ;
  }
*/
  // Fill nr_jet histogram
  if (fHistos)
    {
//    printf("FILLING histograms\n");
    TH1F *h_jets = (TH1F*)fHistos->FindObject("histo1");
    if (h_jets) h_jets->Fill(nJets);

    // There may be no leading jet with more that 2 particles !
    TH1F *h_leadpart = (TH1F*)fHistos->FindObject("histo2");
    if (part_leadjet && h_leadpart) h_leadpart->Fill(part_leadjet);

    TProfile * h_prof = (TProfile*)fHistos->FindObject("histo3");
    if (part_leadjet && h_prof) h_prof->Fill(pt_leadjet,part_leadjet);

    TH1F *h_MD = (TH1F*)fHistos->FindObject("histo4");
      for( Int_t k = 0  ; k < part_leadjet ; k++)
        { h_MD->Fill( z_part_ljet[k] ); }


//     TProfile *h_phi = (TProfile*)fHistos->At(4);
//     for( Int_t k = 0  ; k < nPart ; k++)
//       { h_phi->Fill( TMath::RadToDeg() * dphiArray [k] , nPart ) ; }


//    TH1F *h_tpd = (TH1F*)fHistos->At(5);
//     for( Int_t k = 0  ; k < nPart ; k++)
//       {
//       if ( (phi_partT[k] > (TMath::Pi()/3.)) && (phi_partT[k] < (2. * TMath::Pi()/3.)) )
//         { h_tpd->Fill( pt_partT [k] ) ; }
//       }


    }


  // CLEANING SECTION
  for (Int_t i=0; i<nPart; i++) delete vectParticle[i];
  delete [] vectParticle;
  for (Int_t i=0; i<nJets; i++) delete vectJet[i];
  delete [] vectJet;
  delete [] ptArray;
  delete [] idxArray;

  if (z_part_ljet) delete [] z_part_ljet;
  if (idx_leadjet_part) delete [] idx_leadjet_part;
  if (jets_pt_idx) delete [] jets_pt_idx ;
  if (jets_pt) delete [] jets_pt ;

  vectArray->Delete() ; // TClonesArray of lorentz vectors
  }
// end of FindJets


//______________________________________________________________________________
void AliCdfJetFinder::FinishRun()
{}

