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

/*
$Log$
Revision 1.28  2002/11/20 14:13:16  morsch
- FindChargedJets() added.
- Destructor corrected.
- Geometry cuts taken from AliEMCALGeometry.

Revision 1.27  2002/11/15 17:39:10  morsch
GetPythiaParticleName removed.

Revision 1.26  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.20.4.3  2002/10/10 15:07:49  hristov
Updating VirtualMC to v3-09-02

Revision 1.25  2002/09/13 10:24:57  morsch
idem

Revision 1.24  2002/09/13 10:21:13  morsch
No cast to AliMagFCM.

Revision 1.23  2002/06/27 09:24:26  morsch
Uncomment the TH1::AddDirectory statement.

Revision 1.22  2002/05/22 13:48:43  morsch
Pdg code added to track list.

Revision 1.21  2002/04/27 07:43:08  morsch
Calculation of fDphi corrected (Renan Cabrera)

Revision 1.20  2002/03/12 01:06:23  pavlinov
Testin output from generator

Revision 1.19  2002/02/27 00:46:33  pavlinov
Added method FillFromParticles()

Revision 1.18  2002/02/21 08:48:59  morsch
Correction in FillFromHitFlaggedTrack. (Jennifer Klay)

Revision 1.17  2002/02/14 08:52:53  morsch
Major updates by Aleksei Pavlinov:
FillFromPartons, FillFromTracks, jetfinder configuration.

Revision 1.16  2002/02/08 11:43:05  morsch
SetOutputFileName(..) allows to specify an output file to which the
reconstructed jets are written. If not set output goes to input file.
Attention call Init() before processing.

Revision 1.15  2002/02/02 08:37:18  morsch
Formula for rB corrected.

Revision 1.14  2002/02/01 08:55:30  morsch
Fill map with Et and pT.

Revision 1.13  2002/01/31 09:37:36  morsch
Geometry parameters in constructor and call of SetCellSize()

Revision 1.12  2002/01/23 13:40:23  morsch
Fastidious debug print statement removed.

Revision 1.11  2002/01/22 17:25:47  morsch
Some corrections for event mixing and bg event handling.

Revision 1.10  2002/01/22 10:31:44  morsch
Some correction for bg mixing.

Revision 1.9  2002/01/21 16:28:42  morsch
Correct sign of dphi.

Revision 1.8  2002/01/21 16:05:31  morsch
- different phi-bin for hadron correction
- provisions for background mixing.

Revision 1.7  2002/01/21 15:47:18  morsch
Bending radius correctly in cm.

Revision 1.6  2002/01/21 12:53:50  morsch
authors

Revision 1.5  2002/01/21 12:47:47  morsch
Possibility to include K0long and neutrons.

Revision 1.4  2002/01/21 11:03:21  morsch
Phi propagation introduced in FillFromTracks.

Revision 1.3  2002/01/18 05:07:56  morsch
- hadronic correction
- filling of digits
- track selection upon EMCAL information

*/

//*-- Authors: Andreas Morsch   (CERN)
//*            J.L. Klay        (LBL)
//*            Aleksei Pavlinov (WSU) 

#include <stdio.h>
// From root ...
#include <TROOT.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranchElement.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TList.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TPythia6Calls.h>

// From AliRoot ...
#include "AliEMCALJetFinder.h"
#include "AliEMCALFast.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALHadronCorrection.h"
#include "AliEMCALJetMicroDst.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliMagFCM.h"
#include "AliEMCAL.h"
#include "AliHeader.h"
#include "AliPDG.h"

// Interface to FORTRAN
#include "Ecommon.h"


ClassImp(AliEMCALJetFinder)

//____________________________________________________________________________
AliEMCALJetFinder::AliEMCALJetFinder()
{
// Default constructor
    fJets             = 0;
    fNjets            = 0;
    fLego             = 0;
    fLegoB            = 0;

    fTrackList        = 0;
    fPtT              = 0;
    fEtaT             = 0;
    fPhiT             = 0;
    fPdgT             = 0;
    
    fTrackListB       = 0;
    fPtB              = 0;
    fEtaB             = 0;
    fPhiB             = 0;
    fPdgB             = 0;

    fHCorrection      = 0;
    fHadronCorrector  = 0;

    fWrite            = 0;
    fOutFileName      = 0;
    fOutFile          = 0;
    fInFile           = 0;
    fEvent            = 0;

    SetParametersForBgSubtraction();
}

AliEMCALJetFinder::AliEMCALJetFinder(const char* name, const char *title)
    : TTask(name, title)
{
// Constructor 
// Title is used in method GetFileNameForParameters();
//
    fJets  = new TClonesArray("AliEMCALJet",10000);
    fNjets = 0;
    for (Int_t i = 0; i < 30000; i++)
    {
	fEtCell[i]  = 0.;
	fEtaCell[i] = 0.;
        fPhiCell[i] = 0.;
    }
    fLego       = 0;
    fLegoB      = 0;

    fTrackList  = 0;
    fPtT        = 0;
    fEtaT       = 0;
    fPhiT       = 0;
    fPdgT       = 0;

    fTrackListB       = 0;
    fPtB        = 0;
    fEtaB       = 0;
    fPhiB       = 0;
    fPdgB       = 0;

    fHCorrection      = 0;
    fHadronCorrector  = 0;
    fBackground       = 0;
    fWrite            = 0;
    fOutFileName      = 0;
    fOutFile          = 0;
    fInFile           = 0;
    fEvent            = 0;
//
    SetPtCut();
    SetMomentumSmearing();
    SetEfficiencySim();
    SetDebug();
    SetHadronCorrection();
    SetSamplingFraction();
    SetIncludeK0andN();

    SetParametersForBgSubtraction();
}

void AliEMCALJetFinder::SetParametersForBgSubtraction
(Int_t mode, Float_t minMove, Float_t maxMove, Float_t precBg)
{
// see file /home/pavlinov/cosmic/pythia/jetFinderParamData.inc
// at WSU Linux cluster - 11-feb-2002
// These parameters must be tuned more carefull !!!
  SetMode(mode);
  SetMinMove(minMove);
  SetMaxMove(maxMove);
  SetPrecBg(precBg);
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
    delete fLego;            
    delete fLegoB;
    delete fhLegoTracks;  
    delete fhLegoEMCAL;   
    delete fhLegoHadrCorr;
    delete fhEff;         
    delete fhCellEt;      
    delete fhCellEMCALEt; 
    delete fhTrackPt;     
    delete fhTrackPtBcut; 
    delete fhChPartMultInTpc;

    delete[] fTrackList;
    delete[] fPtT;      
    delete[] fEtaT;     
    delete[] fPhiT;     
    delete[] fPdgT;     
 
    delete[] fTrackListB;
    delete[] fPtB;       
    delete[] fEtaB;      
    delete[] fPhiB;      
    delete[] fPdgB;      
}

#ifndef WIN32
# define jet_finder_ua1 jet_finder_ua1_
# define hf1 hf1_
# define type_of_call

#else
# define jet_finder_ua1 JET_FINDER_UA1
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


void AliEMCALJetFinder::Init()
{
//
// Geometry and I/O initialization
//
//
//
//  Get geometry parameters from EMCAL
//
//
//  Geometry 
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    AliEMCALGeometry* geom = 
	AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
    fNbinEta = geom->GetNZ();
    fNbinPhi = geom->GetNPhi();
    fPhiMin  = geom->GetArm1PhiMin()*TMath::Pi()/180.;
    fPhiMax  = geom->GetArm1PhiMax()*TMath::Pi()/180.;
    fEtaMin  = geom->GetArm1EtaMin();
    fEtaMax  = geom->GetArm1EtaMax();
    fDphi    = (fPhiMax-fPhiMin)/fNbinPhi;
    fDeta    = (fEtaMax-fEtaMin)/fNbinEta;
    fNtot    = fNbinPhi*fNbinEta;
//
    SetCellSize(fDeta, fDphi);
//
//  I/O
    if (fOutFileName) fOutFile = new TFile(fOutFileName, "recreate");
}

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
    if(fWrite) WriteJets();
    fEvent++;
}

void AliEMCALJetFinder::Find()
{
// Wrapper for fortran coded jet finder using member data for 
// argument list

    Float_t min_move = fMinMove;
    Float_t max_move = fMaxMove;
    Int_t   mode     = fMode;
    Float_t prec_bg  = fPrecBg;
    Int_t   ierror;

    ResetJets(); // 4-feb-2002 by PAI

    jet_finder_ua1(fNcell, fNtot, fEtCell, fEtaCell, fPhiCell, 
		   min_move, max_move, mode, prec_bg, ierror);
    fError = ierror;
    // Write result to output
    Int_t njet = Njets();
    
    for (Int_t nj=0; nj<njet; nj++)
    {
	
	fJetT[nj] = new AliEMCALJet(JetEnergy(nj),
				    JetPhiW(nj),
				    JetEtaW(nj));
    }

    FindTracksInJetCone();
    if(fWrite) WriteJets();
    fEvent++;
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
    fConeRadius = par;
}

void AliEMCALJetFinder::SetEtSeed(Float_t par)
{
// Set et cut for seeds
    EMCALJETPARAM.etSeed = par;
    fEtSeed = par;
}

void AliEMCALJetFinder::SetMinJetEt(Float_t par)
{
// Set minimum jet et
    EMCALJETPARAM.ejMin = par;
    fMinJetEt = par;
}

void AliEMCALJetFinder::SetMinCellEt(Float_t par)
{
// Set et cut per cell
    EMCALJETPARAM.etMin = par;
    fMinCellEt = par;
}

void AliEMCALJetFinder::SetPtCut(Float_t par)
{
// Set pT cut on charged tracks
    fPtCut = par;
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
    const char* file = 0;

    Int_t njet = Njets();

    for (Int_t nj = 0; nj < njet; nj++)
    {
	AddJet(*fJetT[nj]);
	delete fJetT[nj];
    }

// I/O
    if (!fOutFileName) {
//
// output written to input file
//
	AliEMCAL* pEMCAL = (AliEMCAL* )gAlice->GetModule("EMCAL");
	TTree* pK = gAlice->TreeK();
	file = (pK->GetCurrentFile())->GetName();
	if (fDebug > 1)
	    printf("Make Branch - TreeR address %p %p\n",gAlice->TreeR(), pEMCAL);
	if (fJets && gAlice->TreeR()) {
	    pEMCAL->MakeBranchInTree(gAlice->TreeR(), 
				     "EMCALJets", 
				     &fJets, 
				     kBufferSize, 
				     file);
	}
	Int_t nev = gAlice->GetHeader()->GetEvent();
	gAlice->TreeR()->Fill();
	char hname[30];
	sprintf(hname,"TreeR%d", nev);
	gAlice->TreeR()->Write(hname);
	gAlice->TreeR()->Reset();
    } else {
//
// Output written to user specified output file
//
	TTree* pK = gAlice->TreeK();
	fInFile  = pK->GetCurrentFile();

	fOutFile->cd();
	char hname[30];
	sprintf(hname,"TreeR%d", fEvent);
	TTree* treeJ = new TTree(hname, "EMCALJets");
	treeJ->Branch("EMCALJets", &fJets, kBufferSize);
	treeJ->Fill();
	treeJ->Write(hname);
	fInFile->cd();
    }
    ResetJets();        
}

void AliEMCALJetFinder::BookLego()
{
//
//  Book histo for discretisation
//

//
//  Don't add histos to the current directory
    if(fDebug) printf("\n AliEMCALJetFinder::BookLego() \n");
 
    TH2::AddDirectory(0);
    TH1::AddDirectory(0);
    gROOT->cd();
//    
//  Signal map
    fLego = new TH2F("legoH","eta-phi",
			   fNbinEta, fEtaMin, fEtaMax, 
			   fNbinPhi, fPhiMin, fPhiMax);
//
//  Background map
    fLegoB = new TH2F("legoB","eta-phi for BG event",
			   fNbinEta, fEtaMin, fEtaMax, 
			   fNbinPhi, fPhiMin, fPhiMax);

//  Tracks map
    fhLegoTracks = new TH2F("hLegoTracks","eta-phi for Tracks",
    fNbinEta, fEtaMin, fEtaMax, fNbinPhi, fPhiMin, fPhiMax);
//  EMCAL map
    fhLegoEMCAL = new TH2F("hLegoEMCAL","eta-phi for EMCAL",
    fNbinEta, fEtaMin, fEtaMax, fNbinPhi, fPhiMin, fPhiMax);
//  Hadron correction map
    fhLegoHadrCorr = new TH2F("hLegoHadrCorr","eta-phi for Hadron. Correction",
    fNbinEta, fEtaMin, fEtaMax, fNbinPhi, fPhiMin, fPhiMax);
//  Hists. for tuning jet finder 
    fhEff = new TH2F("hEff","#epsilon vs momentum ", 100,0.,20., 50,0.5,1.);

    TArrayF eTmp(1101);
    eTmp[0] = 0.0;
    for(Int_t i=1; i<=1000; i++)      eTmp[i] = 0.1*i; // step 100 mev
    for(Int_t i=1001; i<=1100; i++)   eTmp[i] = eTmp[1000] + 1.0*(i-1000); // step 1GeV

    fhCellEt  = new TH1F("hCellEt","Cell E_{T} from fLego", 
    eTmp.GetSize()-1, eTmp.GetArray()); 
    fhCellEMCALEt  = new TH1F("hCellEMCALEt","Cell E_{T} for EMCAL itself", 
    eTmp.GetSize()-1, eTmp.GetArray()); 
    fhTrackPt = new TH1F("hTrackPt","Ch.particles P_{T} ",
    eTmp.GetSize()-1, eTmp.GetArray()); 
    fhTrackPtBcut = new TH1F("hTrackPtBcut","Ch.particles P_{T} + magnetic field cut",
    eTmp.GetSize()-1, eTmp.GetArray()); 

    fhChPartMultInTpc = new TH1F("hChPartMultInTpc",
    "Charge partilcle multiplicity in |%eta|<0.9", 2000, 0, 20000);

            //! first canvas for drawing
    fHistsList=AliEMCALJetMicroDst::MoveHistsToList("Hists from AliEMCALJetFinder", kTRUE);
}

void AliEMCALJetFinder::DumpLego()
{
//
// Dump lego histo into array
//
    fNcell = 0;
    TAxis* Xaxis = fLego->GetXaxis();
    TAxis* Yaxis = fLego->GetYaxis();
    //    fhCellEt->Clear();
    Float_t e, eH;
    for (Int_t i = 1; i <= fNbinEta; i++) {
	for (Int_t j = 1; j <= fNbinPhi; j++) {
	    e = fLego->GetBinContent(i,j);
	    if (e > 0.0) {
	      Float_t eta  = Xaxis->GetBinCenter(i);
	      Float_t phi  = Yaxis->GetBinCenter(j);	    
	      fEtCell[fNcell]  = e;
	      fEtaCell[fNcell] = eta;
	      fPhiCell[fNcell] = phi;
	      fNcell++;
              fhCellEt->Fill(e);
            }
            if(fhLegoEMCAL) {
              eH = fhLegoEMCAL->GetBinContent(i,j);
              if(eH > 0.0) fhCellEMCALEt->Fill(eH);
            }
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
// Fill Cells with track information
//
    if (fDebug >= 2)
    printf("\n AliEMCALJetFinder::FillFromTracks(%i,%i) ",flag,ich);
    fNChTpc = 0;

    ResetMap();
    
    if (!fLego) BookLego();
// Reset
    if (flag == 0) fLego->Reset();
//
// Access particle information    
    Int_t npart = (gAlice->GetHeader())->GetNprimary();
    Int_t ntr   = (gAlice->GetHeader())->GetNtrack();
    printf(" : #primary particles %i # tracks %i \n", npart, ntr);
 
// Create track list
//
// 0: not selected
// 1: selected for jet finding
// 2: inside jet
// ....
    if (fTrackList) delete[] fTrackList;
    if (fPtT)       delete[] fPtT;
    if (fEtaT)      delete[] fEtaT;
    if (fPhiT)      delete[] fPhiT;
    if (fPdgT)      delete[] fPdgT;
    
    fTrackList = new Int_t  [npart];
    fPtT       = new Float_t[npart];
    fEtaT      = new Float_t[npart];
    fPhiT      = new Float_t[npart];
    fPdgT      = new Int_t[npart];

    fNt   = npart;
    fNtS  = 0;
    Float_t chTmp=0.0; // charge of current particle 
    //    Int_t idGeant;

    // this is for Pythia ??
    for (Int_t part = 0; part < npart; part++) {
	TParticle *MPart = gAlice->Particle(part);
	Int_t mpart   = MPart->GetPdgCode();
	Int_t child1  = MPart->GetFirstDaughter();
	Float_t pT    = MPart->Pt();
	Float_t p     = MPart->P();
	Float_t phi   = MPart->Phi();
	Float_t eta   = -100.;
	if(pT > 0.001) eta   = MPart->Eta();
	Float_t theta = MPart->Theta();	
        if  (fDebug>=2) { 
	   printf("ind %7i part %7i -> child1 %5i child2 %5i Status %5i\n", 
           part, mpart, child1, MPart->GetLastDaughter(), MPart->GetStatusCode());
        }
	
	if (fDebug >= 2) {
	    if (part == 6 || part == 7)
	    {
		printf("\n Simulated Jet (pt, eta, phi): %d %f %f %f\n", 
		       part-5, pT, eta, phi);
	    }

//	    if (mpart == 21)
		
//		printf("\n Simulated Jet (pt, eta, phi): %d %d %f %f %f",
//		       part, mpart, pT, eta, phi); 
	}
	
	fTrackList[part] = 0; 
	fPtT[part]       = pT; // must be change after correction for resolution !!!
	fEtaT[part]      = eta;
	fPhiT[part]      = phi;
	fPdgT[part]      = mpart;
	

	if (part < 2) continue;

	// move to fLego->Fill because hadron correction must apply 
	// if particle hit to EMCAL - 11-feb-2002
	//	if (pT == 0 || pT < fPtCut) continue;
	TParticlePDG* pdgP = 0;
// charged or neutral 
	pdgP  = MPart->GetPDG();
        chTmp = pdgP->Charge() / 3.; // 13-feb-2001!!  

	if (ich == 0) {
	    if (chTmp == 0) {
	        if (!fK0N) { 
		    continue;
		} else {
		    if (mpart != kNeutron    &&
			mpart != kNeutronBar &&
			mpart != kK0Long) continue;
		}
	    }
	}

// skip partons
	if (TMath::Abs(mpart) <= 6         ||
	    mpart == 21                    ||
	    mpart == 92) continue;

	if (TMath::Abs(eta)<=0.9) fNChTpc++;
// final state only
	if (child1 >= 0 && child1 < npart) continue;
// acceptance cut
	if (eta > fEtaMax || eta < fEtaMin)    continue;
	if (phi > fPhiMax || phi < fPhiMin )   continue;
//
	if (fDebug >= 3) 
	printf("\n=>nsel:%5d mpart %5d child1 %5d eta %6.2f phi %6.2f pT %6.2f ch %3.0f ",
	part, mpart, child1, eta, phi, pT, chTmp);
//
//
// Momentum smearing goes here ...
//
        fhTrackPt->Fill(pT);
        Float_t pw;
	if (fSmear && TMath::Abs(chTmp)) {
	    pw = AliEMCALFast::SmearMomentum(1,p);
        // p changed - take into account when calculate pT,
	// pz and so on ..  ?
            pT = (pw/p) * pT;
            if(fDebug >= 4) printf("\n Smearing : p %8.4f change to %8.4f ", p, pw);
            p  = pw;
	}
//
// Tracking Efficiency and TPC acceptance goes here ...
	Float_t eff;
	if (fEffic && TMath::Abs(chTmp)) {
	  //	    eff =  AliEMCALFast::Efficiency(1,p);
            eff = 0.95; // for testing 25-feb-2002
            if(fhEff) fhEff->Fill(p, eff);
	    if (AliEMCALFast::RandomReject(eff)) {
              if(fDebug >= 5) printf(" reject due to unefficiency ");
              continue;
            }
	}
//
// Correction of Hadronic Energy goes here ...
//
//
// phi propagation for hadronic correction

	Bool_t curls = kFALSE; // hit two the EMCAL (no curl)
	Float_t phiHC=0.0, dpH=0.0, dphi=0.0, eTdpH=0;
        if(TMath::Abs(chTmp)) {
        // hadr. correction only for charge particle 
	  dphi  = PropagatePhi(pT, chTmp, curls);
	  phiHC = phi + dphi;
	  if (fDebug >= 6) {
             printf("\n Delta phi %f pT %f ", dphi, pT);
	     if (curls) printf("\n !! Track is curling");
          }
          if(!curls) fhTrackPtBcut->Fill(pT);
	
	  if (fHCorrection && !curls) {
	      if (!fHadronCorrector)
		 Fatal("AliEMCALJetFinder",
		    "Hadronic energy correction required but not defined !");

	      dpH    = fHadronCorrector->GetEnergy(p, eta, 7);
              eTdpH  = dpH*TMath::Sin(theta);

	      if (fDebug >= 7) printf(" phi %f phiHC %f eTcorr %f\n", 
	      phi, phiHC, -eTdpH); // correction is negative
	      fLego->Fill(eta, phiHC, -eTdpH);
              fhLegoHadrCorr->Fill(eta, phiHC, eTdpH);
	  }
        }
//
//  More to do ? Just think about it !
//
	if (phi > fPhiMax || phi < fPhiMin )   continue;

        if(TMath::Abs(chTmp) ) { // charge particle
	  if (pT > fPtCut && !curls) {
	     if (fDebug >= 8) printf("Charge :  fLego->Fill(%5.2f, %5.2f, %6.2f, %d)\n",
				      eta , phi, pT, fNtS); 
             fLego->Fill(eta, phi, pT);
             fhLegoTracks->Fill(eta, phi, pT); // 20-feb for checking
	     fTrackList[part] = 1;
	     fNtS++;
          }
	} else if(ich==0 && fK0N) {
	  // case of n, nbar and K0L
	     if (fDebug >= 9) printf("Neutral :  fLego->Fill(%5.2f, %5.2f, %6.2f, %d)\n",
				      eta , phi, pT, fNtS); 
            fLego->Fill(eta, phi, pT);
	    fTrackList[part] = 1;
	    fNtS++;
        }

    } // primary loop    
    DumpLego();
}

void AliEMCALJetFinder::FillFromHits(Int_t flag)
{
//
// Fill Cells with hit information
//
//
    if (fDebug >= 2)
    printf("\n AliEMCALJetFinder::FillFromHits(%i)\n",flag);

    ResetMap();
    
    if (!fLego) BookLego();
//  Reset eta-phi maps if needed
    if (flag == 0) { // default behavior
      fLego->Reset();
      fhLegoTracks->Reset();
      fhLegoEMCAL->Reset();
      fhLegoHadrCorr->Reset();
    }
//  Initialize from background event if available
//
// Access hit information    
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
//
//   Loop over tracks
//
    Int_t nbytes = 0;
    Double_t etH = 0.0;

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

	    Float_t r      =    TMath::Sqrt(x*x+y*y);
	    Float_t theta  =    TMath::ATan2(r,z);
	    Float_t eta    =   -TMath::Log(TMath::Tan(theta/2.));
	    Float_t phi    =    TMath::ATan2(y,x);

	    if (fDebug >= 11) printf("\n Hit %f %f %f %f", x, y, z, eloss);
//	    printf("\n Hit %f %f %f %f", x, y, z, eloss);
	    
            etH = fSamplingF*eloss*TMath::Sin(theta);
	    fLego->Fill(eta, phi, etH);
	    //	    fhLegoEMCAL->Fill(eta, phi, etH);
	} // Hit Loop
    } // Track Loop
    // copy content of fLego to fhLegoEMCAL (fLego and fhLegoEMCAL are identical)
    for(Int_t i=0; i<fLego->GetSize(); i++) (*fhLegoEMCAL)[i] = (*fLego)[i]; 
    //    DumpLego(); ??
}

void AliEMCALJetFinder::FillFromDigits(Int_t flag)
{
//
// Fill Cells with digit information
//
//
    ResetMap();
    
    if (!fLego) BookLego();
    if (flag == 0) fLego->Reset();
    Int_t nbytes;
    

//
//  Connect digits    
//
    TClonesArray* digs      = new TClonesArray("AliEMCALDigit", 10000);
    TTree *treeD = gAlice->TreeD();
    TBranchElement* branchDg = (TBranchElement*)
	treeD->GetBranch("EMCAL");

    if (!branchDg) Fatal("AliEMCALJetFinder", 
			 "Reading digits requested but no digits in file !");
    
    branchDg->SetAddress(&digs);
    Int_t nent = (Int_t) branchDg->GetEntries();
//
//  Connect digitizer
//
    AliEMCALDigitizer* digr = new AliEMCALDigitizer();
    TBranchElement* branchDr = (TBranchElement*) 
	treeD->GetBranch("AliEMCALDigitizer");
    branchDr->SetAddress(&digr);
//
//
    nbytes += branchDg->GetEntry(0);
    nbytes += branchDr->GetEntry(0);
//
//  Get digitizer parameters
    Float_t towerADCped  = digr->GetTowerpedestal();
    Float_t towerADCcha  = digr->GetTowerchannel();
    Float_t preshoADCped = digr->GetPreShopedestal();
    Float_t preshoADCcha = digr->GetPreShochannel();

    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    AliEMCALGeometry* geom = 
	AliEMCALGeometry::GetInstance(pEMCAL->GetTitle(), "");
    
    if (fDebug) {
	Int_t   ndig = digs->GetEntries();
	printf("\n Number of Digits: %d %d\n", ndig, nent);
	printf("\n Parameters: %f %f %f %f\n", 
	       towerADCped, towerADCcha, preshoADCped, preshoADCcha );
	printf("\n Geometry: %d %d\n", geom->GetNEta(), geom->GetNPhi());
    }
    
//
//  Loop over digits
    AliEMCALDigit* sdg;
    TIter next(digs);
    while ((sdg = (AliEMCALDigit*)(next())))
    {
	Double_t pedestal;
	Double_t channel;
	if (sdg->GetId() > (geom->GetNZ() * geom->GetNPhi())) 
	{
	    pedestal = preshoADCped;
	    channel  = preshoADCcha; 
	} else {
	    pedestal = towerADCped;
	    channel  = towerADCcha; 
	}
	
	Float_t eta = sdg->GetEta();
	Float_t phi = sdg->GetPhi() * TMath::Pi()/180.;
	Float_t amp = (Float_t) (channel*(sdg->GetAmp())-pedestal);
	
	if (fDebug > 1) printf("\n Digit: eta %8.3f phi %8.3f amp %8.3f %8d",
			   eta, phi, amp, sdg->GetAmp());
	
	fLego->Fill(eta, phi, fSamplingF*amp);
    } // digit loop
// 
//  Dump lego hist
    DumpLego();
}


void AliEMCALJetFinder::FillFromHitFlaggedTracks(Int_t flag)
{
//
// Fill Cells with hit information
//
//
    ResetMap();
    
    if (!fLego) BookLego();
// Reset
    if (flag == 0) fLego->Reset();

// Flag charged tracks passing through TPC which 
// are associated to EMCAL Hits
    BuildTrackFlagTable();

//
// Access particle information    
    TTree *treeK = gAlice->TreeK();
    Int_t ntracks = (Int_t) treeK->GetEntries();

    if (fPtT)       delete[] fPtT;   
    if (fEtaT)      delete[] fEtaT;    
    if (fPhiT)      delete[] fPhiT;   
    if (fPdgT)      delete[] fPdgT;   
   
    fPtT       = new Float_t[ntracks];
    fEtaT      = new Float_t[ntracks];
    fPhiT      = new Float_t[ntracks];
    fPdgT      = new Int_t[ntracks];

    fNt   = ntracks;
    fNtS  = 0;
    
    for (Int_t track = 0; track < ntracks; track++) {
	  TParticle *MPart = gAlice->Particle(track);
	  Float_t pT    = MPart->Pt();
	  Float_t phi   = MPart->Phi();
	  Float_t eta   = MPart->Eta();

	if(fTrackList[track]) {
	  fPtT[track]       = pT;
	  fEtaT[track]      = eta;
	  fPhiT[track]      = phi;
	  fPdgT[track]      = MPart->GetPdgCode();
	  
	  if (track < 2) continue;	//Colliding particles?
	  if (pT == 0 || pT < fPtCut) continue;
	  fNtS++;
	  fLego->Fill(eta, phi, pT);
	}
      } // track loop
    DumpLego();
}

void AliEMCALJetFinder::FillFromParticles()
{
// 26-feb-2002 PAI - for checking all chain
// Work on particles level; accept all particle (not neutrino )
    
    Double_t PX=0, PY=0, PZ=0, E=0; // checking conservation law 
    fNChTpc = 0;

    ResetMap();
    if (!fLego) BookLego();
    fLego->Reset();
//
// Access particles information    
    Int_t npart = (gAlice->GetHeader())->GetNprimary();
    if (fDebug >= 2 || npart<=0) {
       printf(" AliEMCALJetFinder::FillFromParticles : npart %i\n", npart);
       if(npart<=0) return;
    }
    fNt   = npart;
    fNtS  = 0;
    RearrangeParticlesMemory(npart);
 
//  Go through the particles
    Int_t mpart, child1, child2, geantPdg;
    Float_t pT, phi, eta, e=0, px=0, py=0, pz=0;
    TParticle *MPart=0;
    for (Int_t part = 0; part < npart; part++) {

	fTrackList[part] = 0;

	MPart   = gAlice->Particle(part);
	mpart   = MPart->GetPdgCode();
	child1  = MPart->GetFirstDaughter();
	child2  = MPart->GetLastDaughter();
	pT      = MPart->Pt();
	phi     = MPart->Phi();
	eta     = MPart->Eta();

        px      = MPart->Px();
        py      = MPart->Py();
        pz      = MPart->Pz();
        e       = MPart->Energy();

// see pyedit in Pythia's text
        geantPdg = mpart;
	if (IsThisPartonsOrDiQuark(mpart)) continue;
        printf("%5i: %5i(%2i) px %5.1f py %5.1f pz %6.1f e %6.1f childs %5i,%5i \n", 
        part, mpart, geantPdg, px, py, pz, e, child1, child2);
	
//  exclude partons (21 - gluon, 92 - string) 
	

// exclude neutrinous also ??
	if (fDebug >= 11 && pT>0.01) 
	printf("\n part:%5d mpart %5d eta %9.2f phi %9.2f pT %9.2f ",
	part, mpart, eta, phi, pT);

	fPtT[part]       = pT;
	fEtaT[part]      = eta;
	fPhiT[part]      = phi;
	fPdgT[part]      = mpart;
	fNtS++;
	
// final state only
 	if (child1 >= 0 && child1 < npart) continue;

	//        printf("%4i -> %5i(%3i) px %6.1f py %6.1f pz %7.1f e %8.2f child1 %5i %s\n", 
	//        part, mpart, geantPdg, px, py, pz, e, child1, name.Data());
        PX += px; 
        PY += py; 
        PZ += pz;
        E  += e; 

        
	if (TMath::Abs(eta) <= 0.9) fNChTpc++;
// acceptance cut
	if (eta > fEtaMax || eta < fEtaMin)    continue;
	if (phi > fPhiMax || phi < fPhiMin )   continue;
//
        if(fK0N==0 ) { // exclude neutral hadrons
          if (mpart == kNeutron || mpart == kNeutronBar || mpart == kK0Long) continue; 
        }
	fTrackList[part] = 1;
        fLego->Fill(eta, phi, pT);

    } // primary loop
    printf("\n                PX %8.2f  PY %8.2f  PZ %8.2f  E %8.2f \n", 
    PX, PY, PZ, E);
    DumpLego();
    if(fhChPartMultInTpc) fhChPartMultInTpc->Fill(fNChTpc);
}

void AliEMCALJetFinder::FillFromPartons()
{
// function under construction - 13-feb-2002 PAI
    
    if (fDebug >= 2)
    printf("\n AliEMCALJetFinder::FillFromPartons()\n");
    // 

    ResetMap();
    if (!fLego) BookLego();
    fLego->Reset();
//
// Access particle information    
    Int_t npart = (gAlice->GetHeader())->GetNprimary();
    if (fDebug >= 2 || npart<=0)
    printf("\n AliEMCALJetFinder::FillFromPartons : npart %i\n", npart);
    fNt   = 0; // for FindTracksInJetCone
    fNtS  = 0;
 
//  Go through the partons
    Int_t statusCode=0;
    for (Int_t part = 8; part < npart; part++) {
	TParticle *MPart = gAlice->Particle(part);
	Int_t mpart   = MPart->GetPdgCode();
	//	Int_t child1  = MPart->GetFirstDaughter();
	Float_t pT    = MPart->Pt();
	//	Float_t p     = MPart->P();
	Float_t phi   = MPart->Phi();
	Float_t eta   = MPart->Eta();
	//	Float_t theta = MPart->Theta();
        statusCode    = MPart->GetStatusCode();
	
// accept partons (21 - gluon, 92 - string) 
	if (!(TMath::Abs(mpart) <= 6 || mpart == 21 ||mpart == 92)) continue;
	if (fDebug > 1 && pT>0.01) 
	printf("\n part:%5d mpart %5d status  %5d eta %8.2f phi %8.2f pT %8.2f ",
	part, mpart, statusCode, eta, phi, pT);
	//	if (fDebug >= 3) MPart->Print(); 
// accept partons before fragmentation - p.57 in Pythia manual
//        if(statusCode != 1) continue;
// acceptance cut
	if (eta > fEtaMax || eta < fEtaMin)    continue;
	if (phi > fPhiMax || phi < fPhiMin )   continue;
// final state only
//	if (child1 >= 0 && child1 < npart) continue;
//
//
        fLego->Fill(eta, phi, pT);

    } // primary loop
    DumpLego();
}

void AliEMCALJetFinder::BuildTrackFlagTable() {

// Method to generate a lookup table for TreeK particles
// which are linked to hits in the EMCAL
//
// --Author: J.L. Klay
//
// Access hit information    
    AliEMCAL* pEMCAL = (AliEMCAL*) gAlice->GetModule("EMCAL");
    
    TTree *TK = gAlice->TreeK();		// Get the number of entries in the kine tree
    Int_t nKTrks = (Int_t) TK->GetEntries();	// (Number of particles created somewhere)
    
    if(fTrackList) delete[] fTrackList;		//Make sure we get rid of the old one
    fTrackList = new Int_t[nKTrks];		//before generating a new one
    
    for (Int_t i = 0; i < nKTrks; i++) {	//Initialize members to 0
	fTrackList[i] = 0;
    }
    
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks = (Int_t) treeH->GetEntries();
//
//   Loop over tracks
//
    Int_t nbytes = 0;
    
    for (Int_t track=0; track<ntracks;track++) {
	gAlice->ResetHits();
	nbytes += treeH->GetEvent(track);
        if (pEMCAL)  {
	    
//
//  Loop over hits
//
	    for(AliEMCALHit* mHit=(AliEMCALHit*) pEMCAL->FirstHit(-1); 
		mHit;
		mHit=(AliEMCALHit*) pEMCAL->NextHit()) 
	    {
		Int_t   iTrk   =    mHit->Track();       // track number
		Int_t   idprim =    mHit->GetPrimary();  // primary particle
		
		//Determine the origin point of this particle - it made a hit in the EMCAL
		TParticle 	  *trkPart = gAlice->Particle(iTrk);
		TParticlePDG  *trkPdg  = trkPart->GetPDG();
		Int_t	   trkCode = trkPart->GetPdgCode();
		Double_t 	   trkChg;
		if (trkCode < 10000) {		//Big Ions cause problems for 
		    trkChg  = trkPdg->Charge();	//this function.  Since they aren't
		} else {				//likely to make it very far, set
		    trkChg  = 0.0;			//their charge to 0 for the Flag test
		}
		Float_t 	   vxTrk   = trkPart->Vx();
		Float_t 	   vyTrk   = trkPart->Vy();
		Float_t 	   vrTrk   = TMath::Sqrt(vxTrk*vxTrk+vyTrk*vyTrk); 
		fTrackList[iTrk] = SetTrackFlag(vrTrk,trkCode,trkChg);
		
		//Loop through the ancestry of the EMCAL entrance particles
		Int_t ancestor = trkPart->GetFirstMother();  //Get track's Mother
		while (ancestor != -1) {
		    TParticle     *ancPart = gAlice->Particle(ancestor);  //get particle info on ancestor
		    TParticlePDG  *ancPdg  = ancPart->GetPDG();
		    Int_t	      ancCode = ancPart->GetPdgCode();
		    Double_t       ancChg;
		    if (ancCode < 10000) {
			ancChg  = ancPdg->Charge();
		    } else {
			ancChg  = 0.0;
		    }
		    Float_t 	      vxAnc   = ancPart->Vx();
		    Float_t 	      vyAnc   = ancPart->Vy();
		    Float_t 	      vrAnc   = TMath::Sqrt(vxAnc*vxAnc+vyAnc*vyAnc);
		    fTrackList[ancestor] = SetTrackFlag(vrAnc,ancCode,ancChg);
		    ancestor = ancPart->GetFirstMother();           //Get the next ancestor
		}
		
		//Determine the origin point of the primary particle associated with the hit
		TParticle     *primPart = gAlice->Particle(idprim);
		TParticlePDG  *primPdg  = primPart->GetPDG();
		Int_t	   primCode = primPart->GetPdgCode();
		Double_t 	   primChg;
		if (primCode < 10000) {
		    primChg  = primPdg->Charge();
		} else {
		    primChg  = 0.0;
		}
		Float_t	   vxPrim = primPart->Vx();
		Float_t	   vyPrim = primPart->Vy();
		Float_t	   vrPrim = TMath::Sqrt(vxPrim*vxPrim+vyPrim*vyPrim);
		fTrackList[idprim] = SetTrackFlag(vrPrim,primCode,primChg);
	    } // Hit Loop
	} //if (pEMCAL)
    } // Track Loop
}

Int_t AliEMCALJetFinder
::SetTrackFlag(Float_t radius, Int_t code, Double_t charge) {

    Int_t flag    = 0; 
    Int_t parton  = 0; 
    Int_t neutral = 0;
    
    if (charge == 0) neutral = 1;
    
    if (TMath::Abs(code) <= 6   ||	
	code == 21              ||   
	code == 92) parton = 1;
    
    //It's not a parton, it's charged and it went through the TPC
    if (!parton && !neutral && radius <= 84.0) flag = 1;
    
    return flag;
}



void AliEMCALJetFinder::SaveBackgroundEvent()
{
// Saves the eta-phi lego and the tracklist
//
    if (fLegoB) {
       fLegoB->Reset();
       (*fLegoB) = (*fLegoB) + (*fLego); 
       if(fDebug) 
       printf("\n AliEMCALJetFinder::SaveBackgroundEvent() (fLegoB) %f = %f(fLego) \n", 
       fLegoB->Integral(), fLego->Integral()); 
    }
    
    if (fPtB)        delete[] fPtB;   
    if (fEtaB)       delete[] fEtaB;    
    if (fPhiB)       delete[] fPhiB;   
    if (fPdgB)       delete[] fPdgB;   
    if (fTrackListB) delete[] fTrackListB;   
   
    fPtB          = new Float_t[fNtS];
    fEtaB         = new Float_t[fNtS];
    fPhiB         = new Float_t[fNtS];
    fPdgB         = new Int_t  [fNtS];
    fTrackListB   = new Int_t  [fNtS];
    
    fNtB = 0;
    
    for (Int_t i = 0; i < fNt; i++) {
	if (!fTrackList[i]) continue;
	fPtB [fNtB]       = fPtT [i];
	fEtaB[fNtB]       = fEtaT[i];
	fPhiB[fNtB]       = fPhiT[i];
	fPdgB[fNtB]       = fPdgT[i];
	fTrackListB[fNtB] = 1;
	fNtB++;
    }
    fBackground = 1;
    printf(" fNtB %i => fNtS %i #particles %i \n", fNtB, fNtS, fNt); 
}

void AliEMCALJetFinder::InitFromBackground()
{
//
//    
    if (fDebug) printf("\n AliEMCALJetFinder::InitFromBackground() ");
    
    if (fLego) {
	fLego->Reset(); 
	(*fLego) = (*fLego) + (*fLegoB);
	if(fDebug) 
	    printf("\n AliEMCALJetFinder::SaveBackgroundEvent() (fLego) %f = %f(fLegoB) \n", 
		   fLego->Integral(), fLegoB->Integral()); 
    } else {
	printf(" => fLego undefined \n");
    }
}

    
void AliEMCALJetFinder::FindTracksInJetCone()
{
//
//  Build list of tracks inside jet cone
//
//  Loop over jets
    Int_t njet = Njets();
    for (Int_t nj = 0; nj < njet; nj++)
    {
	Float_t etaj = JetEtaW(nj);
	Float_t phij = JetPhiW(nj);	
	Int_t   nT   = 0;           // number of associated tracks
	
// Loop over particles in current event 
	for (Int_t part = 0; part < fNt; part++) {
	    if (!fTrackList[part]) continue;
	    Float_t phi      = fPhiT[part];
	    Float_t eta      = fEtaT[part];
	    Float_t dr       = TMath::Sqrt((etaj-eta)*(etaj-eta) +
					   (phij-phi)*(phij-phi));
	    if (dr < fConeRadius) {
		fTrackList[part] = nj+2;
		nT++;
	    } // < ConeRadius ?
	} // particle loop
	
// Same for background event if available
	Int_t nTB = 0;
	if (fBackground) {
	    for (Int_t part = 0; part < fNtB; part++) {
		Float_t phi      = fPhiB[part];
		Float_t eta      = fEtaB[part];
		Float_t dr       = TMath::Sqrt((etaj-eta)*(etaj-eta) +
					       (phij-phi)*(phij-phi));
		fTrackListB[part] = 1;

		if (dr < fConeRadius) {
		    fTrackListB[part] = nj+2;
		    nTB++;
		} // < ConeRadius ?
	    } // particle loop
	} // Background available ?
	
	Int_t nT0 = nT + nTB;
	printf("Total number of tracks %d\n", nT0);
	
	if (nT0 > 50) nT0 = 50;
	
	Float_t* ptT  = new Float_t[nT0];
	Float_t* etaT = new Float_t[nT0];
	Float_t* phiT = new Float_t[nT0];
	Int_t*   pdgT = new Int_t[nT0];

	Int_t iT = 0;
	Int_t j;
	
	for (Int_t part = 0; part < fNt; part++) {
	    if (fTrackList[part] == nj+2) {
		Int_t index = 0;
		for (j=iT-1; j>=0; j--) {
		    if (fPtT[part] > ptT[j]) {
			index = j+1;
			break;
		    }
		}
		for (j=iT-1; j>=index; j--) {
		    ptT [j+1]  = ptT [j];
		    etaT[j+1]  = etaT[j];
		    phiT[j+1]  = phiT[j];
		    pdgT[j+1]  = pdgT[j];
		}
		ptT [index] = fPtT [part];
		etaT[index] = fEtaT[part];
		phiT[index] = fPhiT[part];
		pdgT[index] = fPdgT[part];
		iT++;
	    } // particle associated
	    if (iT > nT0) break;
	} // particle loop
	
	if (fBackground) {
	    for (Int_t part = 0; part < fNtB; part++) {
		if (fTrackListB[part] == nj+2) {
		    Int_t index = 0;
		    for (j=iT-1; j>=0; j--) {
			if (fPtB[part] > ptT[j]) {
			    index = j+1;

			    break;
			}
		    }
		    for (j=iT-1; j>=index; j--) {
			ptT [j+1]  = ptT [j];
			etaT[j+1]  = etaT[j];
			phiT[j+1]  = phiT[j];
			pdgT[j+1]  = pdgT[j];
		    }
		    ptT [index] = fPtB [part];
		    etaT[index] = fEtaB[part];
		    phiT[index] = fPhiB[part];
		    pdgT[index] = fPdgB[part];
		    iT++;
		} // particle associated
		if (iT > nT0) break;
	    } // particle loop
	} // Background available ?

	fJetT[nj]->SetTrackList(nT0, ptT, etaT, phiT, pdgT);
	delete[] ptT;
	delete[] etaT;
	delete[] phiT;
	delete[] pdgT;
	
    } // jet loop loop
}

Float_t AliEMCALJetFinder::PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls)
{
// Propagates phi angle to EMCAL radius
//
  static Float_t b = 0.0, rEMCAL = -1.0;
  if(rEMCAL<0) {
// Get field in kGS
    b =  gAlice->Field()->SolenoidField();
// Get EMCAL radius in cm 
    rEMCAL = AliEMCALGeometry::GetInstance()->GetIPDistance();
    printf("\nMagnetic field %f rEMCAL %f ", b, rEMCAL);
  }
    Float_t dPhi = 0.;
//
//
// bending radies
// pt [Gev]
// B  [kG]
//
    Float_t rB = 3335.6 * pt / b;  // [cm]  (case of |charge|=1)
//
// check if particle is curling below EMCAL
    if (2.*rB < rEMCAL) {
	curls = kTRUE;
	return dPhi;
    }
//
// if not calculate delta phi
    Float_t phi = TMath::ACos(1.-rEMCAL*rEMCAL/(2.*rB*rB));
    dPhi = TMath::ATan2(1.-TMath::Cos(phi), TMath::Sin(phi));
    dPhi = -TMath::Sign(dPhi, charge);
//    
    return dPhi;
}

void hf1(Int_t& id, Float_t& x, Float_t& wgt)
{
// dummy for hbook calls
    ;
}

void AliEMCALJetFinder::DrawLego(Char_t *opt) 
{fLego->Draw(opt);}

void AliEMCALJetFinder::DrawLegoEMCAL(Char_t *opt) 
{fhLegoEMCAL->Draw(opt);}

void AliEMCALJetFinder::DrawHistsForTuning(Int_t mode)
{ 
  static TPaveText *varLabel=0;
  if(!fC1) {
    fC1 = new TCanvas("C1","Hists. for tunning", 0,25,600,800);
  }
  fC1->Clear();
  fC1->Divide(2,2);
  fC1->cd(1);
  gPad->SetLogy(1);
  fhCellEt->Draw();

  fC1->cd(2);
  gPad->SetLogy(1);
  fhCellEMCALEt->Draw();

  fC1->cd(3);
  gPad->SetLogy(1);
  fhTrackPt->Draw();
  fhTrackPtBcut->SetLineColor(2);
  fhTrackPtBcut->Draw("same");
 
  fC1->cd(4);
  if(!varLabel) {
    PrintParameters(1);
    varLabel = new TPaveText(0.05,0.5,0.95,0.95,"NDC");
    varLabel->SetTextAlign(12);
    varLabel->SetFillColor(19); // see TAttFill
    TString tmp(GetTitle());
    varLabel->ReadFile(GetFileNameForParameters());
  }
  varLabel->Draw();
  fC1->Update();
  if(mode) { // for saving picture to the file
    TString stmp(GetFileNameForParameters());
    stmp.ReplaceAll("_Par.txt",".ps");
    fC1->Print(stmp.Data());
  }
}

void AliEMCALJetFinder::PrintParameters(Int_t mode)
{
  FILE *file=0;
  if(mode==0) file = stdout; // output to terminal
  else {
    file = fopen(GetFileNameForParameters(),"w");
    if(file==0) file = stdout; 
  }
  fprintf(file,"====   Filling lego   ==== \n");
  fprintf(file,"Smearing          %6i  ", fSmear);
  fprintf(file,"Efficiency        %6i\n", fEffic);
  fprintf(file,"Hadr.Correct.     %6i  ", fHCorrection);
  fprintf(file,"P_{T} Cut of ch.par. %6.2f\n", fPtCut);
  fprintf(file,"====  Jet finding     ==== \n");
  fprintf(file,"Cone radius       %6.2f  ", fConeRadius);
  fprintf(file,"Seed E_{T}           %6.1f\n", fEtSeed);
  fprintf(file,"Min E_{T} of cell    %6.1f  ", fMinCellEt);
  fprintf(file,"Min E_{T} of jet     %6.1f\n", fMinJetEt);
  if(fMode) {
    fprintf(file,"====  Bg subtraction     ==== \n");
    fprintf(file,"BG subtraction  %6i  ", fMode);
    fprintf(file,"Min cone move   %6.2f\n", fMinMove);
    fprintf(file,"Max cone move   %6.2f  ", fMaxMove);
    fprintf(file,"%% change for BG %6.4f\n", fPrecBg);
  } else
  fprintf(file,"==== No Bg subtraction     ==== \n");
  if(file != stdout) fclose(file); 
}

void AliEMCALJetFinder::DrawLegos()
{ 
  if(!fC1) {
    fC1 = new TCanvas("C1","Hists. for tunning", 0,25,600,800);
  }
  fC1->Clear();
  fC1->Divide(2,2);
  gStyle->SetOptStat(111111);

  Int_t nent1, nent2, nent3, nent4;
  Double_t int1, int2, int3, int4;
  nent1 = (Int_t)fLego->GetEntries();
  int1  = fLego->Integral();
  fC1->cd(1);
  if(int1) fLego->Draw("lego");

  nent2 = (Int_t)fhLegoTracks->GetEntries();
  int2  = fhLegoTracks->Integral();
  fC1->cd(2);
  if(int2) fhLegoTracks->Draw("lego");

  nent3 = (Int_t)fhLegoEMCAL->GetEntries();
  int3  = fhLegoEMCAL->Integral();
  fC1->cd(3);
  if(int3) fhLegoEMCAL->Draw("lego");

  nent4 = (Int_t)fhLegoHadrCorr->GetEntries();
  int4  = fhLegoHadrCorr->Integral();
  fC1->cd(4);
  if(int4) fhLegoHadrCorr->Draw("lego");

  // just for checking 
  printf(" Integrals \n");
  printf("lego   %10.3f \ntrack  %10.3f \nhits   %10.3f \nHCorr   %10.3f\n--      %10.3f(must be 0)\n", 
  int1, int2, int3, int4, int1 - (int2 + int3 - int4));
}

const Char_t* AliEMCALJetFinder::GetFileNameForParameters(Char_t* dir)
{
  static TString tmp;
  if(strlen(dir)) tmp = dir;
  tmp += GetTitle();
  tmp += "_Par.txt";
  return tmp.Data();
}

void AliEMCALJetFinder::RearrangeParticlesMemory(Int_t npart)
{ // See FillFromTracks() - npart must be positive
    if (fTrackList) delete[] fTrackList;
    if (fPtT)       delete[] fPtT;
    if (fEtaT)      delete[] fEtaT;
    if (fPhiT)      delete[] fPhiT;
    if (fPdgT)      delete[] fPdgT;
    
    if(npart>0) { 
       fTrackList = new Int_t  [npart];
       fPtT       = new Float_t[npart];
       fEtaT      = new Float_t[npart];
       fPhiT      = new Float_t[npart];
       fPdgT      = new Int_t[npart];
    } else {
       printf("AliEMCALJetFinder::RearrangeParticlesMemory : npart = %d\n", npart);
    }
}

Bool_t AliEMCALJetFinder::IsThisPartonsOrDiQuark(Int_t pdg)
{
  Int_t absPdg = TMath::Abs(pdg);
  if(absPdg<=6) return kTRUE; // quarks
  if(pdg == 21) return kTRUE; // gluon 
  if(pdg == 92) return kTRUE; // string 

  // see p.51 of Pythia Manual
  // Not include diquarks with c and b quark - 4-mar-2002
  //                 ud_0,sd_0,su_0; dd_1,ud_1,uu_1;  sd_1,su_1,ss_1
  static Int_t diquark[9]={2101,3101,3201, 1103,2103,2203,  3103,3203,3303};
  for(Int_t i=0; i<9; i++) if(absPdg == diquark[i])  return kTRUE; // diquarks

  return kFALSE;
}

void AliEMCALJetFinder::FindChargedJet()
{
//
// Find jet from charged particle information only
//
    
//
//  Look for seeds
//
    Int_t njets = 0;
    Int_t part  = 0;
    Int_t nseed = 0;
  
//
//
    ResetJets();
    
//  
    for (part = 0; part < fNt; part++) {
	if (!fTrackList[part]) continue;
	if (fPtT[part] > fEtSeed) nseed++;
    }
    printf("\nThere are %d seeds (%d)\n", nseed, fNtS);
    Int_t* iSeeds = new Int_t[nseed];
    nseed = 0;
    
    for (part = 0; part < fNt; part++) {
	if (!fTrackList[part]) continue;
	if (fPtT[part] > fEtSeed) iSeeds[nseed++] =  part;
    }

//
// Loop over seeds
//
    Int_t seed = 0;
    Float_t pt;
    
    while(1){
//
// Find seed with highest pt
// 
	Float_t ptmax = -1.;
	Int_t   index = -1;
	Int_t   jndex = -1;
	for (seed = 0; seed < nseed; seed++) {
	    if ((pt = fPtT[iSeeds[seed]]) > ptmax && iSeeds[seed] != -1) {
		ptmax = pt;
		index = seed;
	    } // ptmax ?
	} // seeds 
	if (ptmax < 0.) break;
	jndex = iSeeds[index];
//
// Remove from the list   
	iSeeds[index] = -1;
	printf("\n Next Seed %d %f", jndex, ptmax);
//
// Find tracks in cone around seed
//
	Float_t phiSeed = fPhiT[jndex];
	Float_t etaSeed = fEtaT[jndex];
	Float_t eT = 0.;
	Float_t pxJ = 0.;
	Float_t pyJ = 0.;
	Float_t pzJ = 0.;
	
	for (part = 0; part < fNt; part++) {
	    if (!fTrackList[part]) continue;
	    Float_t deta = fEtaT[part] - etaSeed;
	    Float_t dphi = fPhiT[part] - phiSeed;
	    Float_t dR   = TMath::Sqrt(deta * deta + dphi * dphi);
	    if (dR < fConeRadius) {
		eT += fPtT[part];
		Float_t theta = 2. * TMath::ATan(TMath::Exp(-fEtaT[part]));
		Float_t px = fPtT[part] * TMath::Cos(fPhiT[part]);
		Float_t py = fPtT[part] * TMath::Sin(fPhiT[part]);
		Float_t pz = fPtT[part] / TMath::Tan(theta);
		pxJ += px;
		pyJ += py;
		pzJ += pz;
		//
		// if seed, remove it
		//
		for (seed = 0; seed < nseed; seed++) {
		    if (part == iSeeds[seed]) iSeeds[seed] = -1;
		} // seed ?
	    } // < cone radius
	} // particle loop

//
//      Estimate of jet direction
	Float_t phiJ   = TMath::ATan2(pyJ, pxJ);
	Float_t thetaJ = TMath::ATan2(TMath::Sqrt(pxJ * pxJ + pyJ * pyJ), pzJ);
	Float_t etaJ   = TMath::Log(TMath::Tan(thetaJ / 2.));
	Float_t ptJ    = TMath::Sqrt(pxJ * pxJ + pyJ * pyJ);
	
//
//      Sum up all energy
//
	Int_t iPhi0 = Int_t((phiJ - fPhiMin) / fDphi);
	Int_t iEta0 = Int_t((etaJ - fEtaMin) / fDeta);
	Int_t dIphi = Int_t(fConeRadius / fDphi);
	Int_t dIeta = Int_t(fConeRadius / fDeta);
	Int_t iPhi, iEta;
	Float_t sumE = 0;
	for (iPhi = iPhi0 -dIphi; iPhi < iPhi0 + dIphi; iPhi++) {
	    for (iEta = iEta0 - dIeta; iEta < iEta0 + dIeta; iEta++) {
		if (iPhi < 0 || iEta < 0) continue;
		Float_t dPhi = fPhiMin + iPhi * fDphi;
		Float_t dEta = fEtaMin + iEta * fDeta;
		if (TMath::Sqrt(dPhi * dPhi + dEta * dEta) < fConeRadius) continue;
		sumE += fLego->GetBinContent(iEta, iPhi);
	    } // eta
	} // phi
//
//
//
    	fJetT[njets++] = new AliEMCALJet(sumE, phiJ, etaJ);
	FindTracksInJetCone();
	printf("\n Jet Energy %f %f %f %f %d\n", eT, sumE, fPtT[6], fPtT[7], njets);
	printf("\n Jet Phi %f %f %f \n", phiJ, fPhiT[6], fPhiT[7]);
	printf("\n Jet Eta %f %f %f \n", etaJ, fEtaT[6], fEtaT[7]);
    } // while(1)
    EMCALJETS.njet = njets;
    if (fWrite) WriteJets();
    fEvent++;
}
