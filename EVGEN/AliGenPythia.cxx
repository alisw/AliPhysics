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
Revision 1.67  2002/12/09 15:24:09  morsch
Same trigger routine can use Pycell or Pyclus.

Revision 1.66  2002/12/09 08:22:56  morsch
UA1 jet finder (Pycell) for software triggering added.

Revision 1.65  2002/11/19 08:57:10  morsch
Configuration of pt-kick added.

Revision 1.64  2002/11/15 00:43:06  morsch
Changes for kPyJets
- initial and final state g-radiation + pt-kick default
- trigger based on parton clusters (using pyclus)
- trigger jets are stored in header.

Revision 1.63  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.52.4.4  2002/10/10 16:40:08  hristov
Updating VirtualMC to v3-09-02

Revision 1.62  2002/09/24 10:00:01  morsch
CheckTrigger() corrected.

Revision 1.61  2002/07/30 09:52:38  morsch
Call SetGammaPhiRange() and SetGammaEtaRange() in the constructor.

Revision 1.60  2002/07/19 14:49:03  morsch
Typo corrected.

Revision 1.59  2002/07/19 14:35:36  morsch
Count total number of trials. Print mean Q, x1, x2.

Revision 1.58  2002/07/17 10:04:09  morsch
SetYHard method added.

Revision 1.57  2002/05/22 13:22:53  morsch
Process kPyMbNonDiffr added.

Revision 1.56  2002/04/26 10:30:01  morsch
Option kPyBeautyPbMNR added. (N. Carrer).

Revision 1.55  2002/04/17 10:23:56  morsch
Coding Rule violations corrected.

Revision 1.54  2002/03/28 11:49:10  morsch
Pass status code in SetTrack.

Revision 1.53  2002/03/25 14:51:13  morsch
New stack-fill and count options introduced (N. Carrer).

Revision 1.51  2002/03/06 08:46:57  morsch
- Loop until np-1
- delete dyn. alloc. arrays (N. Carrer)

Revision 1.50  2002/03/03 13:48:50  morsch
Option  kPyCharmPbMNR added. Produce charm pairs in agreement with MNR
NLO calculations (Nicola Carrer).

Revision 1.49  2002/02/08 16:50:50  morsch
Add name and title in constructor.

Revision 1.48  2001/12/20 11:44:28  morsch
Add kinematic bias for direct gamma production.

Revision 1.47  2001/12/19 14:45:00  morsch
Store number of trials in header.

Revision 1.46  2001/12/19 10:36:19  morsch
Add possibility if jet kinematic biasing.

Revision 1.45  2001/11/28 08:06:52  morsch
Use fMaxLifeTime parameter.

Revision 1.44  2001/11/27 13:13:07  morsch
Maximum lifetime for long-lived particles to be put on the stack is parameter.
It can be set via SetMaximumLifetime(..).

Revision 1.43  2001/10/21 18:35:56  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.42  2001/10/15 08:21:55  morsch
Vertex truncation settings moved to AliGenMC.

Revision 1.41  2001/10/08 08:45:42  morsch
Possibility of vertex cut added.

Revision 1.40  2001/09/25 11:30:23  morsch
Pass event vertex to header.

Revision 1.39  2001/07/27 17:09:36  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.38  2001/07/13 10:58:54  morsch
- Some coded moved to AliGenMC
- Improved handling of secondary vertices.

Revision 1.37  2001/06/28 11:17:28  morsch
SetEventListRange setter added. Events in specified range are listed for
debugging. (Yuri Kharlov)

Revision 1.36  2001/03/30 07:05:49  morsch
Final print-out in finish run.
Write parton system for jet-production (preliminary solution).

Revision 1.35  2001/03/09 13:03:40  morsch
Process_t and Struc_Func_t moved to AliPythia.h

Revision 1.34  2001/02/14 15:50:40  hristov
The last particle in event marked using SetHighWaterMark

Revision 1.33  2001/01/30 09:23:12  hristov
Streamers removed (R.Brun)

Revision 1.32  2001/01/26 19:55:51  hristov
Major upgrade of AliRoot code

Revision 1.31  2001/01/17 10:54:31  hristov
Better protection against FPE

Revision 1.30  2000/12/18 08:55:35  morsch
Make AliPythia dependent generartors work with new scheme of random number generation

Revision 1.29  2000/12/04 11:22:03  morsch
Init of sRandom as in 1.15

Revision 1.28  2000/12/02 11:41:39  morsch
Use SetRandom() to initialize random number generator in constructor.

Revision 1.27  2000/11/30 20:29:02  morsch
Initialise static variable sRandom in constructor: sRandom = fRandom;

Revision 1.26  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.25  2000/10/18 19:11:27  hristov
Division by zero fixed

Revision 1.24  2000/09/18 10:41:35  morsch
Add possibility to use nuclear structure functions from PDF library V8.

Revision 1.23  2000/09/14 14:05:40  morsch
dito

Revision 1.22  2000/09/14 14:02:22  morsch
- Correct conversion from mm to cm when passing particle vertex to MC.
- Correct handling of fForceDecay == all.

Revision 1.21  2000/09/12 14:14:55  morsch
Call fDecayer->ForceDecay() at the beginning of Generate().

Revision 1.20  2000/09/06 14:29:33  morsch
Use AliPythia for event generation an AliDecayPythia for decays.
Correct handling of "nodecay" option

Revision 1.19  2000/07/11 18:24:56  fca
Coding convention corrections + few minor bug fixes

Revision 1.18  2000/06/30 12:40:34  morsch
Pythia takes care of vertex smearing. Correct conversion from Pythia units (mm) to
Geant units (cm).

Revision 1.17  2000/06/09 20:34:07  morsch
All coding rule violations except RS3 corrected

Revision 1.16  2000/05/15 15:04:20  morsch
The full event is written for fNtrack = -1
Coding rule violations corrected.

Revision 1.15  2000/04/26 10:14:24  morsch
Particles array has one entry more than pythia particle list. Upper bound of
particle loop changed to np-1 (R. Guernane, AM)

Revision 1.14  2000/04/05 08:36:13  morsch
Check status code of particles in Pythia event
to avoid double counting as partonic state and final state particle.

Revision 1.13  1999/11/09 07:38:48  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.12  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.11  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log
*/

//
// Generator using the TPythia interface (via AliPythia)
// to generate pp collisions.
// Using SetNuclei() also nuclear modifications to the structure functions
// can be taken into account. This makes, of course, only sense for the
// generation of the products of hard processes (heavy flavor, jets ...)
//
// andreas.morsch@cern.ch
//


#include "AliGenPythia.h"
#include "AliGenPythiaEventHeader.h"
#include "AliDecayerPythia.h"
#include "AliRun.h"
#include "AliPythia.h"
#include "AliPDG.h"
#include "AliConst.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TDatabasePDG.h>

 ClassImp(AliGenPythia)

AliGenPythia::AliGenPythia()
                 :AliGenMC()
{
// Default Constructor
  fParticles = 0;
  fPythia    = 0;
  fDecayer   = new AliDecayerPythia();
  SetEventListRange();
  SetJetPhiRange();
  SetJetEtaRange();
  SetJetEtRange();
  SetGammaPhiRange();
  SetGammaEtaRange();
  SetPtKick();
}

AliGenPythia::AliGenPythia(Int_t npart)
                 :AliGenMC(npart)
{
// default charm production at 5. 5 TeV
// semimuonic decay
// structure function GRVHO
//
    fName = "Pythia";
    fTitle= "Particle Generator using PYTHIA";
    fXsection  = 0.;
    fNucA1=0;
    fNucA2=0;
    SetProcess();
    SetStrucFunc();
    SetForceDecay();
    SetPtHard();
    SetYHard();
    SetEnergyCMS();
    fDecayer = new AliDecayerPythia();
    // Set random number generator 
    sRandom=fRandom;
    fFlavorSelect   = 0;
    // Produced particles  
    fParticles = new TClonesArray("TParticle",1000);
    fEventVertex.Set(3);
    SetEventListRange();
    SetJetPhiRange();
    SetJetEtaRange();
    SetJetEtRange();
    SetGammaPhiRange();
    SetGammaEtaRange();
    SetJetReconstructionMode();
    SetPtKick();
    // Options determining what to keep in the stack (Heavy flavour generation)
    fStackFillOpt = kFlavorSelection; // Keep particle with selected flavor
    fFeedDownOpt = kTRUE;             // allow feed down from higher family
    // Fragmentation on/off
    fFragmentation = kTRUE;
    // Default counting mode
    fCountMode = kCountAll;
}

AliGenPythia::AliGenPythia(const AliGenPythia & Pythia)
{
// copy constructor
    Pythia.Copy(*this);
}

AliGenPythia::~AliGenPythia()
{
// Destructor
}

void AliGenPythia::SetEventListRange(Int_t eventFirst, Int_t eventLast)
{
  // Set a range of event numbers, for which a table
  // of generated particle will be printed
  fDebugEventFirst = eventFirst;
  fDebugEventLast  = eventLast;
  if (fDebugEventLast==-1) fDebugEventLast=fDebugEventFirst;
}

void AliGenPythia::Init()
{
// Initialisation
    
    SetMC(AliPythia::Instance());
    fPythia=(AliPythia*) fgMCEvGen;
//
    fParentWeight=1./Float_t(fNpart);
//
//  Forward Paramters to the AliPythia object
    fDecayer->SetForceDecay(fForceDecay);    
    fDecayer->Init();


    fPythia->SetCKIN(3,fPtHardMin);
    fPythia->SetCKIN(4,fPtHardMax);
    fPythia->SetCKIN(7,fYHardMin);
    fPythia->SetCKIN(8,fYHardMax);
    
    if (fNucA1 > 0 && fNucA2 > 0) fPythia->SetNuclei(fNucA1, fNucA2);  
    // Fragmentation?
    if (fFragmentation) {
      fPythia->SetMSTP(111,1);
    } else {
      fPythia->SetMSTP(111,0);
    }


//  initial state radiation   
    fPythia->SetMSTP(61,fGinit);
//  final state radiation
    fPythia->SetMSTP(71,fGfinal);
//  pt - kick
    if (fPtKick > 0.) {
	fPythia->SetMSTP(91,1);
	fPythia->SetPARP(91,fPtKick);
    } else {
	fPythia->SetMSTP(91,0);
    }

 //   fPythia->SetMSTJ(1,2);
 //
    fPythia->ProcInit(fProcess,fEnergyCMS,fStrucFunc);

//  Parent and Children Selection
    switch (fProcess) 
    {
    case kPyCharm:
    case kPyCharmUnforced:
    case kPyCharmPbMNR:
	fParentSelect[0] =   411;
	fParentSelect[1] =   421;
	fParentSelect[2] =   431;
	fParentSelect[3] =  4122;
	fFlavorSelect    =  4;	
	break;
    case kPyD0PbMNR:
	fParentSelect[0] =   421;
	fFlavorSelect    =   4;	
	break;
    case kPyBeauty:
    case kPyBeautyPbMNR:
	fParentSelect[0]=  511;
	fParentSelect[1]=  521;
	fParentSelect[2]=  531;
	fParentSelect[3]= 5122;
	fParentSelect[4]= 5132;
	fParentSelect[5]= 5232;
	fParentSelect[6]= 5332;
	fFlavorSelect   = 5;	
	break;
    case kPyBeautyUnforced:
	fParentSelect[0] =  511;
	fParentSelect[1] =  521;
	fParentSelect[2] =  531;
	fParentSelect[3] = 5122;
	fParentSelect[4] = 5132;
	fParentSelect[5] = 5232;
	fParentSelect[6] = 5332;
	fFlavorSelect    = 5;	
	break;
    case kPyJpsiChi:
    case kPyJpsi:
	fParentSelect[0] = 443;
	break;
    case kPyMb:
    case kPyMbNonDiffr:
    case kPyJets:
    case kPyDirectGamma:
	break;
    }
//
//  This counts the total number of calls to Pyevnt() per run.
    fTrialsRun = 0;
    fQ         = 0.;
    fX1        = 0.;
    fX2        = 0.;    
    fNev       = 0 ;
//    
    AliGenMC::Init();
}

void AliGenPythia::Generate()
{
// Generate one event
    fDecayer->ForceDecay();

    Float_t polar[3]   =   {0,0,0};
    Float_t origin[3]  =   {0,0,0};
    Float_t p[3];
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
//
    Int_t nt=0;
    Int_t jev=0;
    Int_t j, kf;
    fTrials=0;

    //  Set collision vertex position 
    if(fVertexSmear==kPerEvent) {
	fPythia->SetMSTP(151,1);
	for (j=0;j<3;j++) {
	    fPythia->SetPARP(151+j, fOsigma[j]*10.);
	}
    } else if (fVertexSmear==kPerTrack) {
	fPythia->SetMSTP(151,0);
    }
//  event loop    
    while(1)
    {
	fPythia->Pyevnt();
	if (gAlice->GetEvNumber()>=fDebugEventFirst &&
	    gAlice->GetEvNumber()<=fDebugEventLast) fPythia->Pylist(1);
	fTrials++;
	
	fPythia->ImportParticles(fParticles,"All");

//
//
//
	Int_t i;
	
	Int_t np = fParticles->GetEntriesFast();
	if (np == 0 ) continue;
// Get event vertex and discard the event if the Z coord. is too big	
	TParticle *iparticle = (TParticle *) fParticles->At(0);
	Float_t distz = iparticle->Vz()/10.;
	if(TMath::Abs(distz)>fCutVertexZ*fOsigma[2]) continue;
//
	fEventVertex[0] = iparticle->Vx()/10.+fOrigin.At(0);
	fEventVertex[1] = iparticle->Vy()/10.+fOrigin.At(1);
	fEventVertex[2] = iparticle->Vz()/10.+fOrigin.At(2);
//
	Int_t* pParent   = new Int_t[np];
	Int_t* pSelected = new Int_t[np];
	Int_t* trackIt   = new Int_t[np];
	for (i=0; i< np; i++) {
	    pParent[i]   = -1;
	    pSelected[i] =  0;
	    trackIt[i]   =  0;
	}

	Int_t nc = 0;        // Total n. of selected particles
	Int_t nParents = 0;  // Selected parents
	Int_t nTkbles = 0;   // Trackable particles
	if (fProcess != kPyMb && fProcess != kPyJets && 
	    fProcess != kPyDirectGamma &&
	    fProcess != kPyMbNonDiffr) {
	    
	    for (i = 0; i<np; i++) {
		iparticle = (TParticle *) fParticles->At(i);
		Int_t ks = iparticle->GetStatusCode();
		kf = CheckPDGCode(iparticle->GetPdgCode());
// No initial state partons
		if (ks==21) continue;
//
// Heavy Flavor Selection
//
		// quark ?
		kf = TMath::Abs(kf);
		Int_t kfl = kf;
		// meson ?
		if  (kfl > 10) kfl/=100;
		// baryon
		if (kfl > 10) kfl/=10;
		if (kfl > 10) kfl/=10;

		Int_t ipa = iparticle->GetFirstMother()-1;
		Int_t kfMo = 0;
		
		if (ipa > -1) {
		    TParticle *  mother = (TParticle *) fParticles->At(ipa);
		    kfMo = TMath::Abs(mother->GetPdgCode());
		}
		// What to keep in Stack?
		Bool_t flavorOK = kFALSE;
		Bool_t selectOK = kFALSE;
		if (fFeedDownOpt) {
		  if (kfl >= fFlavorSelect) flavorOK = kTRUE;
		} else {
		  if (kfl > fFlavorSelect) {
		    nc = -1;
		    break;
		  }
		  if (kfl == fFlavorSelect) flavorOK = kTRUE;
		}
		switch (fStackFillOpt) {
		case kFlavorSelection:
		  selectOK = kTRUE;
		  break;
		case kParentSelection:
		  if (ParentSelected(kf) || kf <= 10) selectOK = kTRUE;
		  break;
		}
		if (flavorOK && selectOK) { 
//
// Heavy flavor hadron or quark
//
// Kinematic seletion on final state heavy flavor mesons
		    if (ParentSelected(kf) && !KinematicSelection(iparticle, 0)) 
		    {
		      continue;
		    }
		    pSelected[i] = 1;
		    if (ParentSelected(kf)) ++nParents; // Update parent count
//		    printf("\n particle (HF)  %d %d %d", i, pSelected[i], kf);
		} else {
// Kinematic seletion on decay products
		    if (fCutOnChild && ParentSelected(kfMo) && ChildSelected(kf) 
			&& !KinematicSelection(iparticle, 1))
		    {
		      continue;
		    }
//
// Decay products 
// Select if mother was selected and is not tracked

		    if (pSelected[ipa] && 
			!trackIt[ipa]  &&     // mother will be  tracked ?
			kfMo !=  5 &&         // mother is b-quark, don't store fragments          
			kfMo !=  4 &&         // mother is c-quark, don't store fragments 
			kf   != 92)           // don't store string
		    {
//
// Semi-stable or de-selected: diselect decay products:
// 
//
			if (pSelected[i] == -1 ||  fDecayer->GetLifetime(kf) > fMaxLifeTime)
			{
			    Int_t ipF = iparticle->GetFirstDaughter();
			    Int_t ipL = iparticle->GetLastDaughter();	
			    if (ipF > 0) for (j = ipF-1; j < ipL; j++) pSelected[j] = -1;
			}
//			printf("\n particle (decay)  %d %d %d", i, pSelected[i], kf);
			pSelected[i] = (pSelected[i] == -1) ? 0 : 1;
		    }
		}
		if (pSelected[i] == -1) pSelected[i] = 0;
		if (!pSelected[i]) continue;
		// Count quarks only if you did not include fragmentation
		if (fFragmentation && kf <= 10) continue;
		nc++;
// Decision on tracking
		trackIt[i] = 0;
//
// Track final state particle
		if (ks == 1) trackIt[i] = 1;
// Track semi-stable particles
		if ((ks ==1) || (fDecayer->GetLifetime(kf) > fMaxLifeTime))  trackIt[i] = 1;
// Track particles selected by process if undecayed. 
		if (fForceDecay == kNoDecay) {
		    if (ParentSelected(kf)) trackIt[i] = 1;
		} else {
		    if (ParentSelected(kf)) trackIt[i] = 0;
		}
		if (trackIt[i] == 1) ++nTkbles; // Update trackable counter
//
//

  	    } // particle selection loop
	    if (nc > 0) {
		for (i = 0; i<np; i++) {
		    if (!pSelected[i]) continue;
		    TParticle *  iparticle = (TParticle *) fParticles->At(i);
		    kf = CheckPDGCode(iparticle->GetPdgCode());
		    Int_t ks = iparticle->GetStatusCode();  
		    p[0] = iparticle->Px();
		    p[1] = iparticle->Py();
		    p[2] = iparticle->Pz();
		    origin[0] = fOrigin[0]+iparticle->Vx()/10.;
		    origin[1] = fOrigin[1]+iparticle->Vy()/10.;
		    origin[2] = fOrigin[2]+iparticle->Vz()/10.;
		    Float_t tof   = kconv*iparticle->T();
		    Int_t ipa     = iparticle->GetFirstMother()-1;
		    Int_t iparent = (ipa > -1) ? pParent[ipa] : -1;
		    SetTrack(fTrackIt*trackIt[i] ,
				     iparent, kf, p, origin, polar, tof, kPPrimary, nt, 1., ks);
		    pParent[i] = nt;
		    KeepTrack(nt); 
		} //  SetTrack loop
	    }
  	} else {
	    nc = GenerateMB();
	} // mb ?

	if (pParent)   delete[] pParent;
	if (pSelected) delete[] pSelected;
	if (trackIt)   delete[] trackIt;

	if (nc > 0) {
	  switch (fCountMode) {
	  case kCountAll:
	    // printf(" Count all \n");
	    jev += nc;
	    break;
	  case kCountParents:
	    // printf(" Count parents \n");
	    jev += nParents;
	    break;
	  case kCountTrackables:
	    // printf(" Count trackable \n");
	    jev += nTkbles;
	    break;
	  }
	    if (jev >= fNpart || fNpart == -1) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		printf("\n Trials: %i %i %i\n",fTrials, fNpart, jev);

		fQ  += fPythia->GetVINT(51);
		fX1 += fPythia->GetVINT(41);
		fX2 += fPythia->GetVINT(42);
		fTrialsRun += fTrials;
		fNev++;
		MakeHeader();
		break;
	    }
	}
    } // event loop
    SetHighWaterMark(nt);
//  adjust weight due to kinematic selection
    AdjustWeights();
//  get cross-section
    fXsection=fPythia->GetPARI(1);
}

Int_t  AliGenPythia::GenerateMB()
{
//
// Min Bias selection and other global selections
//
    Int_t i, kf, nt, iparent;
    Int_t nc = 0;
    Float_t p[3];
    Float_t polar[3]   =   {0,0,0};
    Float_t origin[3]  =   {0,0,0};
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
    
    Int_t np = fParticles->GetEntriesFast();
    Int_t* pParent = new Int_t[np];
    for (i=0; i< np; i++) pParent[i] = -1;
    if (fProcess == kPyJets || fProcess == kPyDirectGamma) {
	TParticle* jet1 = (TParticle *) fParticles->At(6);
	TParticle* jet2 = (TParticle *) fParticles->At(7);
	if (!CheckTrigger(jet1, jet2)) return 0;
    }
    
    for (i = 0; i<np; i++) {
	Int_t trackIt = 0;
	TParticle *  iparticle = (TParticle *) fParticles->At(i);
	kf = CheckPDGCode(iparticle->GetPdgCode());
	Int_t ks = iparticle->GetStatusCode();
	Int_t km = iparticle->GetFirstMother();
	if ((ks == 1  && kf!=0 && KinematicSelection(iparticle, 0)) ||
	    (ks != 1) ||
	    (fProcess == kPyJets && ks == 21 && km == 0 && i>1)) {
	    nc++;
	    if (ks == 1) trackIt = 1;
	    Int_t ipa = iparticle->GetFirstMother()-1;
	    
	    iparent = (ipa > -1) ? pParent[ipa] : -1;
	    
//
// store track information
	    p[0] = iparticle->Px();
	    p[1] = iparticle->Py();
	    p[2] = iparticle->Pz();
	    origin[0] = fOrigin[0]+iparticle->Vx()/10.;
	    origin[1] = fOrigin[1]+iparticle->Vy()/10.;
	    origin[2] = fOrigin[2]+iparticle->Vz()/10.;
	    Float_t tof=kconv*iparticle->T();
	    SetTrack(fTrackIt*trackIt, iparent, kf, p, origin, polar,
		     tof, kPPrimary, nt, 1., ks);
	    KeepTrack(nt);
	    pParent[i] = nt;
	} // select particle
    } // particle loop 

    if (pParent) delete[] pParent;
    
    printf("\n I've put %i particles on the stack \n",nc);
    return nc;
}


void AliGenPythia::FinishRun()
{
// Print x-section summary
    fPythia->Pystat(1);
    fQ  /= fNev;
    fX1 /= fNev;
    fX2 /= fNev;    
    printf("\nTotal number of Pyevnt() calls %d\n", fTrialsRun);
    printf("\nMean Q, x1, x2: %f %f %f\n", fQ, fX1, fX2);
    

}

void AliGenPythia::AdjustWeights()
{
// Adjust the weights after generation of all events
//
    TParticle *part;
    Int_t ntrack=gAlice->GetNtrack();
    for (Int_t i=0; i<ntrack; i++) {
	part= gAlice->Particle(i);
	part->SetWeight(part->GetWeight()*fKineBias);
    }
}
    
void AliGenPythia::SetNuclei(Int_t a1, Int_t a2)
{
// Treat protons as inside nuclei with mass numbers a1 and a2  
    fNucA1 = a1;
    fNucA2 = a2;
}


void AliGenPythia::MakeHeader()
{
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenPythiaEventHeader("Pythia");
//
// Event type  
    ((AliGenPythiaEventHeader*) header)->SetProcessType(fPythia->GetMSTI(1));
//
// Number of trials
    ((AliGenPythiaEventHeader*) header)->SetTrials(fTrials);
//
// Event Vertex 
    header->SetPrimaryVertex(fEventVertex);
//
// Jets that have triggered
    if (fProcess == kPyJets)
    {
	Int_t ntrig, njet;
	Float_t jets[4][10];
	GetJets(njet, ntrig, jets);
	
	for (Int_t i = 0; i < ntrig; i++) {
	    ((AliGenPythiaEventHeader*) header)->AddJet(jets[0][i], jets[1][i], jets[2][i], 
							jets[3][i]);
	}
    }
    gAlice->SetGenEventHeader(header);
}
	

Bool_t AliGenPythia::CheckTrigger(TParticle* jet1, TParticle* jet2)
{
// Check the kinematic trigger condition
//
    Double_t eta[2];
    eta[0] = jet1->Eta();
    eta[1] = jet2->Eta();
    Double_t phi[2];
    phi[0] = jet1->Phi();
    phi[1] = jet2->Phi();
    Int_t    pdg[2]; 
    pdg[0] = jet1->GetPdgCode();
    pdg[1] = jet2->GetPdgCode();    
    Bool_t   triggered = kFALSE;

    if (fProcess == kPyJets) {
	Int_t njets = 0;
	Int_t ntrig = 0;
	Float_t jets[4][10];
//
// Use Pythia clustering on parton level to determine jet axis
//
	GetJets(njets, ntrig, jets);
	
	if (ntrig) triggered = kTRUE;
//
    } else {
	Int_t ij = 0;
	Int_t ig = 1;
	if (pdg[0] == kGamma) {
	    ij = 1;
	    ig = 0;
	}
	//Check eta range first...
	if ((eta[ij] < fEtaMaxJet   && eta[ij] > fEtaMinJet) &&
	    (eta[ig] < fEtaMaxGamma && eta[ig] > fEtaMinGamma))
	{
	    //Eta is okay, now check phi range
	    if ((phi[ij] < fPhiMaxJet   && phi[ij] > fPhiMinJet) &&
		(phi[ig] < fPhiMaxGamma && phi[ig] > fPhiMinGamma))
	    {
		triggered = kTRUE;
	    }
	}
    }
    return triggered;
}
	  
AliGenPythia& AliGenPythia::operator=(const  AliGenPythia& rhs)
{
// Assignment operator
    return *this;
}

void  AliGenPythia::LoadEvent()
{
//
// Load event into Pythia Common Block
//
 

    Int_t npart = (Int_t) (gAlice->TreeK())->GetEntries(); 
   (fPythia->GetPyjets())->N = npart;

    for (Int_t part = 0; part < npart; part++) {
	TParticle *MPart = gAlice->Particle(part);
	Int_t kf     = MPart->GetPdgCode();
	Int_t ks     = MPart->GetStatusCode();
	Float_t px = MPart->Px();
	Float_t py = MPart->Py();
	Float_t pz = MPart->Pz();
	Float_t e  = MPart->Energy();
	Float_t p  = TMath::Sqrt(px * px + py * py + pz * pz);
	Float_t m  = TMath::Sqrt(e * e - p * p);
	
	
	(fPythia->GetPyjets())->P[0][part] = px;
	(fPythia->GetPyjets())->P[1][part] = py;
	(fPythia->GetPyjets())->P[2][part] = pz;
	(fPythia->GetPyjets())->P[3][part] = e;
	(fPythia->GetPyjets())->P[4][part] = m;
	
	(fPythia->GetPyjets())->K[1][part] = kf;
	(fPythia->GetPyjets())->K[0][part] = ks;
    }
}

void AliGenPythia::RecJetsUA1(Float_t eCellMin, Float_t eCellSeed, Float_t eMin, Float_t rMax, 
			      Int_t& njets, Float_t jets [4][50])
{
//
//  Calls the Pythia jet finding algorithm to find jets in the current event
//
//
//  Configure detector (EMCAL like)
//
    fPythia->SetPARU(51,2.);
    fPythia->SetMSTU(51,Int_t(96 * 2./0.7));
    fPythia->SetMSTU(52,3 * 144);
//
//  Configure Jet Finder
//  
    fPythia->SetPARU(58, eCellMin);
    fPythia->SetPARU(52, eCellSeed);
    fPythia->SetPARU(53, eMin);
    fPythia->SetPARU(54, rMax);
    fPythia->SetMSTU(54, 2);
//
//  Save jets
    Int_t n     = fPythia->GetN();

//
//  Run Jet Finder
    fPythia->Pycell(njets);
    Int_t i;
    for (i = 0; i < njets; i++) {
	Float_t px    = (fPythia->GetPyjets())->P[0][n+i];
	Float_t py    = (fPythia->GetPyjets())->P[1][n+i];
	Float_t pz    = (fPythia->GetPyjets())->P[2][n+i];
	Float_t e     = (fPythia->GetPyjets())->P[3][n+i];

	jets[0][i] = px;
	jets[1][i] = py;
	jets[2][i] = pz;
	jets[3][i] = e;
    }
}



void  AliGenPythia::GetJets(Int_t& nJets, Int_t& nJetsTrig, Float_t jets[4][10])
{
//
//  Calls the Pythia clustering algorithm to find jets in the current event
//
    Int_t n     = fPythia->GetN();
    nJets       = 0;
    nJetsTrig   = 0;
    if (fJetReconstruction == kCluster) {
//
//  Configure cluster algorithm
//    
	fPythia->SetPARU(43, 2.);
	fPythia->SetMSTU(41, 1);
//
//  Call cluster algorithm
//    
	fPythia->Pyclus(nJets);
//
//  Loading jets from common block
//
    } else {
//
//  Configure detector (EMCAL like)
//
	fPythia->SetPARU(51,2.);
	fPythia->SetMSTU(51,Int_t(96 * 2./0.7));
	fPythia->SetMSTU(52,3 * 144);
//
//  Configure Jet Finder
//  
	fPythia->SetPARU(58,  0.0);
	fPythia->SetPARU(52,  4.0);
	fPythia->SetPARU(53, 10.0);
	fPythia->SetPARU(54,  1.0);
	fPythia->SetMSTU(54,  2);
//
//  Run Jet Finder
	fPythia->Pycell(nJets);
    }

    Int_t i;
    for (i = 0; i < nJets; i++) {
	Float_t px    = (fPythia->GetPyjets())->P[0][n+i];
	Float_t py    = (fPythia->GetPyjets())->P[1][n+i];
	Float_t pz    = (fPythia->GetPyjets())->P[2][n+i];
	Float_t e     = (fPythia->GetPyjets())->P[3][n+i];
	Float_t pt    = TMath::Sqrt(px * px + py * py);
	Float_t phi   = TMath::ATan2(py,px);
	Float_t theta = TMath::ATan2(pt,pz);
	Float_t et    = e * TMath::Sin(theta);
	Float_t eta   = -TMath::Log(TMath::Tan(theta / 2.));

	if (
	    eta > fEtaMinJet && eta < fEtaMaxJet && 
	    phi > fPhiMinJet && eta < fPhiMaxJet &&
	    et  > fEtMinJet  && et  < fEtMaxJet     
	    ) 
	{
	    jets[0][nJetsTrig] = px;
	    jets[1][nJetsTrig] = py;
	    jets[2][nJetsTrig] = pz;
	    jets[3][nJetsTrig] = e;
	    nJetsTrig++;
	    
	} else {
//	    printf("\n........-Jet #%d: %10.3f %10.3f %10.3f %10.3f \n", i, pt, et, eta, phi * kRaddeg);
	}
    }
}


#ifdef never
void AliGenPythia::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGenPythia.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliGenerator::Streamer(R__b);
      R__b >> (Int_t&)fProcess;
      R__b >> (Int_t&)fStrucFunc;
      R__b >> (Int_t&)fForceDecay;
      R__b >> fEnergyCMS;
      R__b >> fKineBias;
      R__b >> fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b >> fXsection;
//      (AliPythia::Instance())->Streamer(R__b);
      R__b >> fPtHardMin;
      R__b >> fPtHardMax;
//      if (fDecayer) fDecayer->Streamer(R__b);
   } else {
      R__b.WriteVersion(AliGenPythia::IsA());
      AliGenerator::Streamer(R__b);
      R__b << (Int_t)fProcess;
      R__b << (Int_t)fStrucFunc;
      R__b << (Int_t)fForceDecay;
      R__b << fEnergyCMS;
      R__b << fKineBias;
      R__b << fTrials;
      fParentSelect.Streamer(R__b);
      fChildSelect.Streamer(R__b);
      R__b << fXsection;
//      R__b << fPythia;
      R__b << fPtHardMin;
      R__b << fPtHardMax;
      //     fDecayer->Streamer(R__b);
   }
}
#endif

