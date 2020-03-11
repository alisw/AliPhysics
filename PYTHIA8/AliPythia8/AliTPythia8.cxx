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

// Author: Andreas Morsch   27/10/2007
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TPythia8                                                                   //
//                                                                            //
// TPythia is an interface class to C++ version of Pythia 8.1                 //
// event generators, written by T.Sjostrand.                                  //
//                                                                            //
// The user is assumed to be familiar with the Pythia package.                //
// This class includes only a basic interface to Pythia8. Because Pythia8 is  //
// also written in C++, its functions/classes can be called directly from a   //
// compiled C++ script.                                                       //
// To call Pythia functions not available in this interface a dictionary must //
// be generated.                                                              //
// see $ROOTSYS/tutorials/pythia/pythia8.C for an example of use from CINT.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
/*
*------------------------------------------------------------------------------------*
 |                                                                                    |
 |  *------------------------------------------------------------------------------*  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  |   PPP   Y   Y  TTTTT  H   H  III    A      Welcome to the Lund Monte Carlo!  |  |
 |  |   P  P   Y Y     T    H   H   I    A A     This is PYTHIA version 8.100      |  |
 |  |   PPP     Y      T    HHHHH   I   AAAAA    Last date of change: 20 Oct 2007  |  |
 |  |   P       Y      T    H   H   I   A   A                                      |  |
 |  |   P       Y      T    H   H  III  A   A    Now is 27 Oct 2007 at 18:26:53    |  |
 |  |                                                                              |  |
 |  |   Main author: Torbjorn Sjostrand; CERN/PH, CH-1211 Geneva, Switzerland,     |  |
 |  |     and Department of Theoretical Physics, Lund University, Lund, Sweden;    |  |
 |  |     phone: + 41 - 22 - 767 82 27; e-mail: torbjorn@thep.lu.se                |  |
 |  |   Author: Stephen Mrenna; Computing Division, Simulations Group,             |  |
 |  |     Fermi National Accelerator Laboratory, MS 234, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 1 - 630 - 840 - 2556; e-mail: mrenna@fnal.gov                   |  |
 |  |   Author: Peter Skands; CERN/PH, CH-1211 Geneva, Switzerland,                |  |
 |  |     and Theoretical Physics Department,                                      |  |
 |  |     Fermi National Accelerator Laboratory, MS 106, Batavia, IL 60510, USA;   |  |
 |  |     phone: + 41 - 22 - 767 24 59; e-mail: skands@fnal.gov                    |  |
 |  |                                                                              |  |
 |  |   The main program reference is the 'Brief Introduction to PYTHIA 8.1',      |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, arXiv:0710.3820                     |  |
 |  |                                                                              |  |
 |  |   The main physics reference is the 'PYTHIA 6.4 Physics and Manual',         |  |
 |  |   T. Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026 [hep-ph/0603175]. |  |
 |  |                                                                              |  |
 |  |   An archive of program versions and documentation is found on the web:      |  |
 |  |   http://www.thep.lu.se/~torbjorn/Pythia.html                                |  |
 |  |                                                                              |  |
 |  |   This program is released under the GNU General Public Licence version 2.   |  |
 |  |   Please respect the MCnet Guidelines for Event Generator Authors and Users. |  |
 |  |                                                                              |  |
 |  |   Disclaimer: this program comes without any guarantees.                     |  |
 |  |   Beware of errors and use common sense when interpreting results.           |  |
 |  |                                                                              |  |
 |  |   Copyright (C) 2007 Torbjorn Sjostrand                                      |  |
 |  |                                                                              |  |
 |  |                                                                              |  |
 |  *------------------------------------------------------------------------------*  |
 |                                                                                    |
 *------------------------------------------------------------------------------------*
*/

#include "AliTPythia8.h"

#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TH1F.h"

ClassImp(AliTPythia8)

AliTPythia8*  AliTPythia8::fgInstance = 0;
char*         AliTPythia8::fgXmldocPath = 0;

//___________________________________________________________________________
AliTPythia8::AliTPythia8():
    TGenerator("AliTPythia8", "AliTPythia8"),
    fPythia(0),
    fPythiaPartonLevel(0),
    fNumberOfParticles(0),
    fSuperPositionNcoll(0),
    fLastNMPI(0),
    fLastNSuperposition(1)
{
   // Constructor
   if (fgInstance) 
      Fatal("AliTPythia8", "There's already an instance of AliTPythia8");
  
   delete fParticles; // was allocated as TObjArray in TGenerator
    
   fParticles = new TClonesArray("TParticle",50);
   if (fgXmldocPath != 0) {
     fPythia    = new Pythia8::Pythia(fgXmldocPath);
   } else {
     fPythia    = new Pythia8::Pythia();
   }
   fgInstance = this; 
}

//___________________________________________________________________________
AliTPythia8::AliTPythia8(const char *xmlDir):
    TGenerator("AliTPythia8", "AliTPythia8"),
    fPythia(0),
    fPythiaPartonLevel(0),
    fNumberOfParticles(0),
    fSuperPositionNcoll(0),
    fLastNMPI(0),
    fLastNSuperposition(1)
{
   // Constructor with an xmlDir (eg "../xmldoc"
   if (fgInstance) 
      Fatal("AliTPythia8", "There's already an instance of AliTPythia8");
  
   delete fParticles; // was allocated as TObjArray in TGenerator
    
   fParticles = new TClonesArray("TParticle",50);
   fPythia    = new Pythia8::Pythia(xmlDir);
}

//___________________________________________________________________________
AliTPythia8::~AliTPythia8()
{
   // Destructor
   if (fParticles) {
      fParticles->Delete();
      delete fParticles;
      fParticles = 0;
   }
   delete fPythia;
   
   if (fPythiaPartonLevel)
      delete fPythiaPartonLevel;
   
   if (fSuperPositionNcoll)
      delete fSuperPositionNcoll;
}

//___________________________________________________________________________
AliTPythia8* AliTPythia8::Instance()
{
   // Return an instance of AliTPythia8
   return fgInstance ? fgInstance : (fgInstance = new AliTPythia8()) ;
}

void AliTPythia8::SetPythiaSeed(UInt_t seed)
{
  //
  // set seed in PYTHIA8
  // NB. 900000000 is the maximum seed (0 is not allowed)
  //
  ReadString("Random:setSeed = on");
  ReadString(Form("Random:seed = %d", (seed % 900000000) + 1));
  
  if (fPythiaPartonLevel) {
    fPythiaPartonLevel->readString("Random:setSeed = on");
    fPythiaPartonLevel->readString(Form("Random:seed = %d", (seed % 900000000) + 1));
  }
}

void AliTPythia8::ActivateEventSuperposition(TH1 *histNcoll) 
{
	// Activates superposition mode. This means for each generated event a number is sampled from histNcoll. 
	// This number of Pythia events is generated and superimposed before hadronization and color reconnection.
	
	fSuperPositionNcoll = histNcoll;
	
	// create parton-level object
	GetPythiaPartonLevel();
}

Pythia8::Pythia* AliTPythia8::GetPythiaPartonLevel()
{
  // Return the parton level object, used when events are superpositioned
  // Basic initialization is done (can be changed by the user by retrieving the instance
  // NOTE that the standard call AliTPythia8::Instnce() gives a different object, so check where the initialization should go
  
  if (!fPythiaPartonLevel) {
    if (fgXmldocPath != 0)
      fPythiaPartonLevel    = new Pythia8::Pythia(fgXmldocPath);
    else
      fPythiaPartonLevel    = new Pythia8::Pythia();

    // special settings to avoid color reconnection and hadronization
    fPythiaPartonLevel->readString("HadronLevel:all = off");
    fPythiaPartonLevel->readString("ColourReconnection:reconnect = off");
  }
 
  return fPythiaPartonLevel;
}

//___________________________________________________________________________
Bool_t AliTPythia8::Initialize(Int_t idAin, Int_t idBin, Double_t ecms)
{
   // Initialization
   AddParticlesToPdgDataBase();
   UpdateParticleProperties();
   ReadString(Form("Beams:eCM = %13.4f", ecms));
   ReadString(Form("Beams:idA = %10d", idAin));
   ReadString(Form("Beams:idB = %10d", idBin));

   if (fPythiaPartonLevel)
     if (fPythiaPartonLevel->init() == kFALSE)
       return false;
   
   return fPythia->init();
}

//___________________________________________________________________________
void AliTPythia8::GenerateEvent()
{
  // Generate the next event
  
  if (fSuperPositionNcoll == NULL) { // default behavior
      fPythia->next();
      fLastNMPI = fPythia->info.nMPI();
  }
  else { // event superpositioning
    fLastNSuperposition = TMath::Nint(fSuperPositionNcoll->GetRandom());
    if (fLastNSuperposition <= 0)
      Fatal("AliTPythia8::GenerateEvent", "Invalid number of collisions requested (%d). Check histogram.", fLastNSuperposition);
      
    Int_t trial = 0;
    const Int_t kMaxTrial = 10;
    while (++trial <= kMaxTrial) {
      fLastNMPI = 0;
      for (Int_t i = 0; i < fLastNSuperposition; i++) {
        fPythiaPartonLevel->next();
        if (i == 0) {
          fPythia->event = fPythiaPartonLevel->event;
          fPythia->info = fPythiaPartonLevel->info;
        }
        else 
          AddEvent(fPythia, fPythiaPartonLevel);
        fLastNMPI += fPythiaPartonLevel->info.nMPI();
      }
        
      // hadronization code
      if (!fPythia->forceHadronLevel(false)) {
        Error("AliTPythia8::GenerateEvent", "Hadronization failed for %d collisions", fLastNSuperposition);
        continue;
      }
      if (!CheckEvent(fPythia->event)) {
        Error("AliTPythia8::GenerateEvent", "Checking event failed for %d collisions", fLastNSuperposition);
        continue;
      }
        
      Info("AliTPythia8::GenerateEvent", "Generated %d superimposed collisions with total %d particles (after %d trials)", fLastNSuperposition, fPythia->event.size() - 1, trial);
      break;
    }
    if (trial > kMaxTrial)
      Fatal("AliTPythia8::GenerateEvent", "Could not superimpose event after %d trials.", trial);
  }
  
  fNumberOfParticles = fPythia->event.size() - 1;
  
  ImportParticles();
}

void AliTPythia8::AddEvent(Pythia8::Pythia* target, Pythia8::Pythia* source) 
{
  //
  // add the particle lists of two pythia events
  // finalEvent += event
  // based on code by Jesper Roy Christiansen, Lund
  //

  // event content
  Pythia8::Event& finalEvent = target->event;
  Pythia8::Event& event      = source->event;

  Int_t offsetCol = finalEvent.lastColTag();
  Int_t offsetIdx = finalEvent.size() - 3;
  
//     Printf("Starting size: %d", finalEvent.size());
//     Printf("Size to add: %d", event.size());
  
  for(Int_t i = 0; i < event.size(); i++){
      if(i == 0){
          finalEvent[0].e(finalEvent[0].e() + event[0].e());
          finalEvent[0].m(finalEvent[0].m() + event[0].m());
          continue;
      } else if(i == 1 || i == 2) {
          finalEvent[i].p(finalEvent[i].p() + event[i].p());
      } else {
          Pythia8::Particle temp = event[i];
          
          // Add offset to nonzero mother, daughter and colour indices.
          // Only if they do not point to the two beam particles.
          if (temp.mother1() > 2) temp.mother1( temp.mother1() + offsetIdx );
          if (temp.mother2() > 2) temp.mother2( temp.mother2() + offsetIdx );
          if (temp.daughter1() > 2) temp.daughter1( temp.daughter1() + offsetIdx );
          if (temp.daughter2() > 2) temp.daughter2( temp.daughter2() + offsetIdx );
          if (temp.col() > 0) temp.col( temp.col() + offsetCol );
          if (temp.acol() > 0) temp.acol( temp.acol() + offsetCol );
          
          // Add the event.
          finalEvent.append(temp);
      }
  }
  
  // Read out junctions one by one.
  Pythia8::Junction tempJ;
  int begCol, endCol;
  for (int i = 0; i < event.sizeJunction(); ++i) {
      tempJ = event.getJunction(i);
      
      // Add colour offsets to all three legs.
      for (int  j = 0; j < 3; ++j) {
          begCol = tempJ.col(j);
          endCol = tempJ.endCol(j);
          if (begCol > 0) begCol += offsetCol;
          if (endCol > 0) endCol += offsetCol;
          tempJ.cols( j, begCol, endCol);
      }
      
      // Append junction to summed event.
      finalEvent.appendJunction( tempJ );
  }
//     Printf("Final size: %d", finalEvent.size());
}

Bool_t AliTPythia8::CheckEvent(Pythia8::Event event) 
{
	//
	// check if event is valid
	// needed for event superpositioning
	// based on code by Jesper Roy Christiansen, Lund
	//
    
    // Check for non-a-number.
    for (int i = 0;i < event.size(); ++i) {
        if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0.
            && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0.
            && abs(event[i].m()) >= 0.);
        else {
            cout << "Error not-a-number energ/momentum/mass" << endl;
            return false;
        }
        
        // Look for particles with not-a-number vertex/lifetime.
        if (abs(event[i].xProd()) >= 0. && abs(event[i].yProd()) >= 0.
            && abs(event[i].zProd()) >= 0.  && abs(event[i].tProd()) >= 0.
            && abs(event[i].tau()) >= 0.) ;
        else {
            cout << "Error not-a-number vertex/lifetime" << endl;
            return false;
        }
        
        // Check for negative energy.
        if (event[i].e() < 0) {
            cout << "Error negative energy" << endl;
            event.list();
            cout << "particle: " << i << endl;
            exit(1);
            return false;
        }
    }
    
    return true;
}

//___________________________________________________________________________
Int_t AliTPythia8::ImportParticles(TClonesArray *particles, Option_t *option)
{
   // Import particles from Pythia stack
   if (particles == 0) return 0;
   TClonesArray &clonesParticles = *particles;
   clonesParticles.Clear();
   Int_t nparts=0;
   Int_t i;
   Int_t ioff = 0;
   fNumberOfParticles  = fPythia->event.size();
   if (fPythia->event[0].id() == 90) {
     ioff = -1;
   }

   if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (i = 0; i < fNumberOfParticles; i++) {
	if (fPythia->event[i].id() == 90) continue;
         if (fPythia->event[i].isFinal()) {
            new(clonesParticles[nparts]) TParticle(
                fPythia->event[i].id(),
                fPythia->event[i].isFinal(),
                fPythia->event[i].mother1() + ioff,
                fPythia->event[i].mother2() + ioff,
                fPythia->event[i].daughter1() + ioff, 
                fPythia->event[i].daughter2() + ioff,
                fPythia->event[i].px(),     // [GeV/c]
                fPythia->event[i].py(),     // [GeV/c]
                fPythia->event[i].pz(),     // [GeV/c]
                fPythia->event[i].e(),      // [GeV]
                fPythia->event[i].xProd(),  // [mm]
                fPythia->event[i].yProd(),  // [mm]
                fPythia->event[i].zProd(),  // [mm]
                fPythia->event[i].tProd()); // [mm/c] 
		nparts++;
	    } // final state partice
	} // particle loop
    } else if (!strcmp(option,"All")) {
	for (i = 0; i < fNumberOfParticles; i++) {
	  if (fPythia->event[i].id() == 90) continue;
	    new(clonesParticles[nparts]) TParticle(
		fPythia->event[i].id(),
		fPythia->event[i].isFinal(),
		fPythia->event[i].mother1() + ioff,
		fPythia->event[i].mother2() + ioff,
		fPythia->event[i].daughter1() + ioff,
		fPythia->event[i].daughter2() + ioff,
		fPythia->event[i].px(),       // [GeV/c]
		fPythia->event[i].py(),       // [GeV/c]
		fPythia->event[i].pz(),       // [GeV/c]
		fPythia->event[i].e(),        // [GeV]
		fPythia->event[i].xProd(),    // [mm]
		fPythia->event[i].yProd(),    // [mm]
		fPythia->event[i].zProd(),    // [mm]
		fPythia->event[i].tProd());   // [mm/c]
     	        nparts++;
	} // particle loop	
    }
    if(ioff==-1)     fNumberOfParticles--; 
    return nparts;
}

//___________________________________________________________________________
TObjArray* AliTPythia8::ImportParticles(Option_t* /* option */)
{
   // Import particles from Pythia stack
   fParticles->Clear();
   Int_t nparts = 0;
   Int_t ioff = 0;
   fNumberOfParticles  = fPythia->event.size();
   if (fPythia->event[0].id() == 90) {
     ioff = -1;
   }


   TClonesArray &a = *((TClonesArray*)fParticles);
   for (Int_t i = 0; i < fNumberOfParticles; i++) {
     if (fPythia->event[i].id() == 90) continue;
      new(a[nparts]) TParticle(
         fPythia->event[i].id(),
         fPythia->event[i].isFinal(),
         fPythia->event[i].mother1()   + ioff,
         fPythia->event[i].mother2()   + ioff,
         fPythia->event[i].daughter1() + ioff,
         fPythia->event[i].daughter2() + ioff,
         fPythia->event[i].px(),       // [GeV/c]
         fPythia->event[i].py(),       // [GeV/c]
         fPythia->event[i].pz(),       // [GeV/c]
         fPythia->event[i].e(),        // [GeV]
         fPythia->event[i].xProd(),    // [mm]
         fPythia->event[i].yProd(),    // [mm]
         fPythia->event[i].zProd(),    // [mm]
         fPythia->event[i].tProd());   // [mm/c]
      nparts++;
   }
   if (ioff == -1) fNumberOfParticles--;
   return fParticles;
}

//___________________________________________________________________________
Int_t AliTPythia8::GetN() const
{
   // Initialization
   return (fPythia->event.size() - 1);
}

//___________________________________________________________________________
void AliTPythia8::ReadString(const char* string) const
{
   // Configuration
   fPythia->readString(string);
   
   if (fPythiaPartonLevel)
   	  fPythiaPartonLevel->readString(string);
}

//___________________________________________________________________________
void  AliTPythia8::ReadConfigFile(const char* string) const
{
  // Configuration
  fPythia->readFile(string);
}

//___________________________________________________________________________
void AliTPythia8::PrintStatistics() const
{
   // Print end of run statistics
   fPythia->stat();
}

//___________________________________________________________________________
void AliTPythia8::EventListing() const
{
   // Event listing
   fPythia->event.list();
}

//___________________________________________________________________________
void AliTPythia8::AddParticlesToPdgDataBase() const
{
   // Add some pythia specific particle code to the data base    

   TDatabasePDG *pdgDB = TDatabasePDG::Instance();
   pdgDB->AddParticle("string","string", 0, kTRUE,
                      0, 0, "QCD string", 90);
   pdgDB->AddParticle("rho_diff0", "rho_diff0", 0, kTRUE,
                      0, 0, "QCD diffr. state", 9900110);
   pdgDB->AddParticle("pi_diffr+", "pi_diffr+", 0, kTRUE,
                      0, 1, "QCD diffr. state", 9900210);
   pdgDB->AddParticle("omega_di", "omega_di", 0, kTRUE,
                      0, 0, "QCD diffr. state", 9900220);
   pdgDB->AddParticle("phi_diff","phi_diff", 0, kTRUE,
                      0, 0, "QCD diffr. state", 9900330);
   pdgDB->AddParticle("J/psi_di", "J/psi_di", 0, kTRUE,
                      0, 0, "QCD diffr. state", 9900440);
   pdgDB->AddParticle("n_diffr0","n_diffr0",0,kTRUE,
                      0, 0, "QCD diffr. state", 9902110);
   pdgDB->AddParticle("p_diffr+","p_diffr+", 0, kTRUE,
                      0, 1, "QCD diffr. state", 9902210);
}

//___________________________________________________________________________
void AliTPythia8::UpdateParticleProperties() const
{
  // Set up2date lifetimes for hadrons
  // lambda_b from PDG 2019: tau0 = 1.471 ps = 441 m/c = 0.441 mm/c
  ReadString("5122:tau0 = 4.41000e-01");
}
