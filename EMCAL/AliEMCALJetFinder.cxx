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

/* $Id$ */
//*-- Author: Andreas Morsch (CERN)
#include <TClonesArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TAxis.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include "AliEMCALJetFinder.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "Ecommon.h"
#include "AliRun.h"
#include "AliEMCAL.h"
#include "AliHeader.h"

ClassImp(AliEMCALJetFinder)

//____________________________________________________________________________
AliEMCALJetFinder::AliEMCALJetFinder()
{
// Default constructor
    fJets  = 0;
    fNjets = 0;
    fLego  = 0;
}

AliEMCALJetFinder::AliEMCALJetFinder(const char* name, const char *title)
    : TTask(name, title)
{
// Constructor 
    fJets  = new TClonesArray("AliEMCALJet",10000);
    fNjets = 0;
    for (Int_t i=0; i<30000; i++)
    {
	fEtCell[i]  = 0.;
	fEtaCell[i] = 0.;
        fPhiCell[i] = 0.;
    }
    fLego = 0;
}




//____________________________________________________________________________
AliEMCALJetFinder::~AliEMCALJetFinder()
{
// Destructor
//
    if (fJets){
	fJets->Delete();
	delete fJets;
    }
}



#ifndef WIN32
# define jet_finder_ua1 jet_finder_ua1_
# define hf1 hf1_
# define type_of_call

#else
# define jet_finder_ua1 J
# define hf1 HF1
# define type_of_call _stdcall
#endif

extern "C" void type_of_call 
jet_finder_ua1(Int_t& ncell, Int_t& ncell_tot,
	       Float_t  etc[30000],  Float_t etac[30000],
	       Float_t  phic[30000], 
	       Float_t& min_move, Float_t& max_move, Int_t& mode, 
	       Float_t& prec_bg,  Int_t& ierror);

extern "C" void type_of_call hf1(Int_t& id, Float_t& x, Float_t& wgt);





void AliEMCALJetFinder::Find(Int_t ncell, Int_t ncell_tot, Float_t etc[30000], 
			     Float_t etac[30000], Float_t phic[30000],
			     Float_t min_move, Float_t max_move, Int_t mode, 
			     Float_t prec_bg,  Int_t   ierror)
{
// Wrapper for fortran coded jet finder
// Full argument list
    jet_finder_ua1(ncell, ncell_tot, etc, etac, phic, 
		   min_move, max_move, mode, prec_bg, ierror);
    // Write result to output
    WriteJets();
}

void AliEMCALJetFinder::Find()
{
// Wrapper for fortran coded jet finder using member data for 
// argument list

    Float_t min_move = 0;
    Float_t max_move = 0;
    Int_t   mode     = 0;
    Float_t prec_bg  = 0.;
    Int_t   ierror   = 0;

    jet_finder_ua1(fNcell, fNtot, fEtCell, fEtaCell, fPhiCell, 
		   min_move, max_move, mode, prec_bg, ierror);
    // Write result to output
    WriteJets();
}


Int_t AliEMCALJetFinder::Njets()
{
// Get number of reconstructed jets
    return EMCALJETS.njet;
}

Float_t AliEMCALJetFinder::JetEnergy(Int_t i)
{
// Get reconstructed jet energy
    return EMCALJETS.etj[i];
}

Float_t AliEMCALJetFinder::JetPhiL(Int_t i)
{
// Get reconstructed jet phi from leading particle
    return EMCALJETS.phij[0][i];
}

Float_t AliEMCALJetFinder::JetPhiW(Int_t i)
{
// Get reconstructed jet phi from weighting
    return EMCALJETS.phij[1][i];
}

Float_t  AliEMCALJetFinder::JetEtaL(Int_t i)
{
// Get reconstructed jet eta from leading particles
    return EMCALJETS.etaj[0][i];
}


Float_t  AliEMCALJetFinder::JetEtaW(Int_t i)  
{
// Get reconstructed jet eta from weighting
    return EMCALJETS.etaj[1][i];
}

void AliEMCALJetFinder::SetCellSize(Float_t eta, Float_t phi)
{
// Set grid cell size
    EMCALCELLGEO.etaCellSize = eta;
    EMCALCELLGEO.phiCellSize = phi;    
}

void AliEMCALJetFinder::SetConeRadius(Float_t par)
{
// Set jet cone radius
    EMCALJETPARAM.coneRad = par;
}

void AliEMCALJetFinder::SetEtSeed(Float_t par)
{
// Set et cut for seeds
    EMCALJETPARAM.etSeed = par;
}

void AliEMCALJetFinder::SetMinJetEt(Float_t par)
{
// Set minimum jet et
    EMCALJETPARAM.ejMin = par;
}

void AliEMCALJetFinder::SetMinCellEt(Float_t par)
{
// Set et cut per cell
    EMCALJETPARAM.etMin = par;
}


void AliEMCALJetFinder::Test()
{
//
// Test the finder call
//
    const Int_t nmax = 30000;
    Int_t ncell      = 10;
    Int_t ncell_tot  = 100;

    Float_t etc[nmax];
    Float_t etac[nmax];
    Float_t phic[nmax];
    Float_t min_move = 0;
    Float_t max_move = 0;
    Int_t   mode = 0;
    Float_t prec_bg = 0;
    Int_t   ierror = 0;

    
    Find(ncell, ncell_tot, etc, etac, phic, 
	 min_move, max_move, mode, prec_bg, ierror);

}

//
//  I/O
//	

void AliEMCALJetFinder::AddJet(const AliEMCALJet& jet)
{
    //
    // Add a jet 
    //
    TClonesArray &lrawcl = *fJets;
    new(lrawcl[fNjets++]) AliEMCALJet(jet);
}

void AliEMCALJetFinder::ResetJets()
{
    //
    // Reset Jet List 
    //
    fJets->Clear();
    fNjets = 0;
}

void AliEMCALJetFinder::WriteJets()
{
//
// Add all jets to the list
//
    const Int_t kBufferSize = 4000;
    TTree *pK = gAlice->TreeK();
    const char* file = (pK->GetCurrentFile())->GetName();
// I/O
    AliEMCAL* pEMCAL = (AliEMCAL* )gAlice->GetModule("EMCAL");
    printf("Make Branch - TreeR address %p %p\n",gAlice->TreeR(), pEMCAL);
    if (fJets && gAlice->TreeR()) {
	pEMCAL->MakeBranchInTree(gAlice->TreeR(), 
				 "Jets", 
				 &fJets, 
				 kBufferSize, 
				 file);
    }
    Int_t njet = Njets();
    for (Int_t nj=0; nj<njet; nj++)
    {
	AliEMCALJet* jet = new AliEMCALJet(JetEnergy(nj),
					   JetPhiW(nj),
					   JetEtaW(nj));
	AddJet(*jet);
	delete jet;
    }
    Int_t nev = gAlice->GetHeader()->GetEvent();
    gAlice->TreeR()->Fill();
    char hname[30];
    sprintf(hname,"TreeR%d", nev);
    gAlice->TreeR()->Write(hname);
    gAlice->TreeR()->Reset();
    ResetJets();        
}

void AliEMCALJetFinder::BookLego()
{
//
//  Book histo for discretisation
//
//
//  Get geometry parameters from 
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    AliEMCALGeometry* geom = 
	AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
    fNbinEta = geom->GetNZ();
    fNbinPhi = geom->GetNPhi();
    const Float_t  phiMin  = geom->GetArm1PhiMin()*TMath::Pi()/180.;
    const Float_t  phiMax  = geom->GetArm1PhiMax()*TMath::Pi()/180.;
    fDphi   = (phiMax-phiMin)/fNbinEta;
    fDeta   = 1.4/fNbinEta;
    fNtot   = fNbinPhi*fNbinEta;
    
//    
    fLego = new TH2F("legoH","eta-phi",
			   fNbinEta, -0.7,  0.7, 
			   fNbinPhi, phiMin, phiMax);
}

void AliEMCALJetFinder::DumpLego()
{
//
// Dump lego histo into array
//
    fNcell = 0;
    for (Int_t i = 1; i < fNbinEta; i++) {
	for (Int_t j = 1; j < fNbinPhi; j++) {
	    Float_t e    = fLego->GetBinContent(i,j);
	    TAxis* Xaxis = fLego->GetXaxis();
	    TAxis* Yaxis = fLego->GetYaxis();
	    Float_t eta  = Xaxis->GetBinCenter(i);
	    Float_t phi  = Yaxis->GetBinCenter(j);	    
	    fEtCell[fNcell]  = e;
	    fEtaCell[fNcell] = eta;
	    fPhiCell[fNcell] = phi;
	    fNcell++;
	} // phi
    } // eta
    fNcell--;
}

void AliEMCALJetFinder::ResetMap()
{
//
// Reset eta-phi array

    for (Int_t i=0; i<30000; i++)
    {
	fEtCell[i]  = 0.;
	fEtaCell[i] = 0.;
	fPhiCell[i] = 0.;
    }
}


void AliEMCALJetFinder::FillFromTracks(Int_t flag, Int_t ich)
{
//
// Fill Cells with hit information
//
//
    ResetMap();
    
    if (!fLego) BookLego();
// Reset
    if (flag == 0) fLego->Reset();
//
// Access particle information    
    Int_t npart = (gAlice->GetHeader())->GetNprimary();
    for (Int_t part=2; part<npart; part++) {
	TParticle *MPart = gAlice->Particle(part);
	Int_t mpart   = MPart->GetPdgCode();
	Int_t child1  = MPart->GetFirstDaughter();
	Float_t pT    = MPart->Pt();
	Float_t phi   = MPart->Phi();
	Float_t theta = MPart->Theta();
	Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
//	if (part == 6 || part == 7)
//	{
//	    printf("\n Simulated Jet (pt, eta, phi): %d %f %f %f", 
//		   part-5, pT, eta, phi);
//	}
	    
	if (pT == 0.) continue;
// charged or neutral 
	if (ich == 0) {
	    TParticlePDG* pdgP = MPart->GetPDG();
	    if (pdgP->Charge() == 0) continue;
	} 
// skip partons
	if (TMath::Abs(mpart) <= 6         ||
	    mpart == 21                    ||
	    mpart == 92) continue;
// acceptance cut
	if (TMath::Abs(eta) > 0.7)         continue;
	if (phi*180./TMath::Pi() > 120.)   continue;
// final state only
	if (child1 >= 0 && child1 < npart) continue;
//	printf("\n sel:%5d %5d %5d %8.2f %8.2f %8.2f", 
//	part, mpart, child1, eta, phi, pT);
	fLego->Fill(eta, phi, pT);
    } // primary loop
    DumpLego();
}

void AliEMCALJetFinder::FillFromHits(Int_t flag)
{
//
// Fill Cells with track information
//
//
    ResetMap();
    
    if (!fLego) BookLego();
    if (flag == 0) fLego->Reset();

//
// Access hit information    
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
//    AliEMCALGeometry* geom = 
//	AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
    
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
//
//   Loop over tracks
//
    Int_t nbytes = 0;

    
    for (Int_t track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	nbytes += treeH->GetEvent(track);
//
//  Loop over hits
//
	for(AliEMCALHit* mHit=(AliEMCALHit*) pEMCAL->FirstHit(-1); 
	    mHit;
	    mHit=(AliEMCALHit*) pEMCAL->NextHit()) 
	{
	    Float_t x      =    mHit->X();         // x-pos of hit
	    Float_t y      =    mHit->Y();         // y-pos
	    Float_t z      =    mHit->Z();         // z-pos
	    Float_t eloss  =    mHit->GetEnergy(); // deposited energy
//	    Int_t   index  =    mHit->GetId();     // cell index
//	    Float_t eta, phi;
//	    geom->EtaPhiFromIndex(index,  eta, phi);
	    Float_t r     = TMath::Sqrt(x*x+y*y);
	    Float_t theta = TMath::ATan2(r,z);
	    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
	    Float_t phi   = TMath::ATan2(y,x);
	    fLego->Fill(eta, phi, eloss);
//	    if (eloss > 1.) printf("\nx,y,z %f %f %f %f %f", 
//	    r, z, eta, phi, eloss);
//	    printf("\n Max %f", fLego->GetMaximum());
	} // Hit Loop
    } // Track Loop
    DumpLego();
}


void hf1(Int_t& id, Float_t& x, Float_t& wgt)
{
// dummy for hbook calls
    ;
}


