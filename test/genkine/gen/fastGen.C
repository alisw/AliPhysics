// Example: generation of kinematics tree with selected properties.
// Below we select events containing the decays D* -> D0 pi, D0 -> K- pi+
// inside the barrel part of the ALICE detector (45 < theta < 135)

// To be able to compile, you can add -I<ALIROOT_INSTALL_PATH>/include to ~/.rootrc
// or use gSystem->SetIncludePath("-I<ALIROOT_INSTALL_PATH>/include")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TArrayI.h>
#include <TTree.h>

#include "AliGenerator.h"
#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythia.h"
#include "AliPythia.h"
#endif

// Forward declarations of utility functions
Float_t EtaToTheta(Float_t arg);
void GetFinalDecayProducts(Int_t ind, AliStack & stack , TArrayI & ar);

void fastGen(Int_t nev = 1, const char* filename = "galice.root")
{
  // If we want to run with Root only, we need the line below
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  AliPDG::AddParticlesToPdgDataBase();
  TDatabasePDG::Instance();
 
  //=======================================================================
  // Set Random Number seed
  gRandom->SetSeed(12345); // Set 0 to use the current time
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 

  //=======================================================================
  // Prepare the simulation environment
  
  // Create run loader
  AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
  
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);
  
  // Create stack
  rl->MakeStack();
  AliStack* stack = rl->Stack();
  
  // Create header
  AliHeader* header = rl->GetHeader();
  
  // Create and Initialize Generator
 
  // Example of charm generation taken from Config_PythiaHeavyFlavours.C
  AliGenPythia *gener = new AliGenPythia(-1);
  gener->SetEnergyCMS(14000.);
  gener->SetMomentumRange(0,999999);
  gener->SetPhiRange(0., 360.);
  gener->SetThetaRange(0.,180.);
  //  gener->SetProcess(kPyCharmppMNR); // Correct Pt distribution, wrong mult
  gener->SetProcess(kPyMb); // Correct multiplicity, wrong Pt
  gener->SetStrucFunc(kCTEQ4L);
  gener->SetPtHard(2.1,-1.0);
  gener->SetFeedDownHigherFamily(kFALSE);
  gener->SetStack(stack);
  gener->Init();

  // Go to galice.root
  rl->CdGAFile();

  // Forbid some decays. Do it after gener->Init(), because
  // the initialization of the generator includes reading of the decay table.

  AliPythia * py= AliPythia::Instance();
  py->SetMDME(737,1,0); //forbid D*+->D+ + pi0
  py->SetMDME(738,1,0);//forbid D*+->D+ + gamma

  // Forbid all D0 decays except D0->K- pi+
  for(Int_t d=747; d<=762; d++){ 
    py->SetMDME(d,1,0);
  }
  // decay 763 is D0->K- pi+
  for(Int_t d=764; d<=807; d++){ 
    py->SetMDME(d,1,0);
  }

  
  //=======================================================================
  // Main event loop
  TStopwatch timer;
  timer.Start();
  for (Int_t iev = 0; iev < nev; iev++) { // Event Loop
    
    cout <<"Event number "<< iev << endl;
    
    // Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    
    //---------------------------------------------------------------------
    // Generate event

    // Counters
    Int_t nprim = 0;  // number of primary particles
    Int_t ntrial = 0; // number of trials
    Int_t ndstar = 0; // number of D* in the event

    while(!ndstar) {// Selection of events with D*
      
      stack->Reset();
      stack->ConnectTree(rl->TreeK());
      gener->Generate();
      ntrial++;
      nprim = stack->GetNprimary();
      
      for(Int_t ipart =0; ipart < nprim; ipart++){
        TParticle * part = stack->Particle(ipart);
        if(part)    {
          
          if (TMath::Abs(part->GetPdgCode())== 413) {

	    TArrayI daughtersId;

	    GetFinalDecayProducts(ipart,*stack,daughtersId);

	    Bool_t kineOK = kTRUE;

	    Double_t thetaMin = TMath::Pi()/4;
	    Double_t thetaMax = 3*TMath::Pi()/4;

	    for (Int_t id=1; id<=daughtersId[0]; id++) {
	      TParticle * daughter = stack->Particle(daughtersId[id]);
	      if (!daughter) {
		kineOK = kFALSE;
		break;
	      }

	      Double_t theta = daughter->Theta();
	      if (theta<thetaMin || theta>thetaMax) {
		kineOK = kFALSE;
		break;
	      }
	    }

	    if (!kineOK) continue;

            part->Print();
            ndstar++;     
	    
          }
        }
      }   
    } // End of selection of events with D* 
      
    cout << "Number of particles " << nprim << endl;
    cout << "Number of trials " << ntrial << endl;
    
    // Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    
    // I/O
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");
    
  } // End of event loop
  timer.Stop();
  timer.Print();
  
  // Termination

  // Generator
  gener->FinishRun();
  // Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();
}



Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}


void GetFinalDecayProducts(Int_t ind, AliStack & stack , TArrayI & ar){

  // Recursive algorithm to get the final decay products of a particle
  //
  // ind is the index of the particle in the AliStack
  // stack is the particle stack from the generator
  // ar contains the indexes of the final decay products
  // ar[0] is the number of final decay products

  if (ind<0 || ind>stack.GetNtrack()) {
    cerr << "Invalid index of the particle " << ind << endl;
    return;
  } 
  if (ar.GetSize()==0) {
    ar.Set(10);
    ar[0] = 0;
  }

  TParticle * part = stack.Particle(ind);

  Int_t iFirstDaughter = part->GetFirstDaughter();
  if( iFirstDaughter<0) {
    // This particle is a final decay product, add its index to the array
    ar[0]++;
    if (ar.GetSize() <= ar[0]) ar.Set(ar.GetSize()+10); // resize if needed
    ar[ar[0]] = ind;
    return;
  } 

  Int_t iLastDaughter = part->GetLastDaughter();

  for (Int_t id=iFirstDaughter; id<=iLastDaughter;id++) {
    // Now search for final decay products of the daughters
    GetFinalDecayProducts(id,stack,ar);
  }
}
