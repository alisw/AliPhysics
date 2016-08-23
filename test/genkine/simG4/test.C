// Usage in compiled mode
// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");  
// gROOT->LoadMacro("test.C+");
// test()

#if !defined(__CINT__) || defined(__MAKECINT__)

// Root include files
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TParticle.h>

// AliRoot include files
#include "AliESDEvent.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"

#endif

void test() {
  
  TStopwatch timer;
  timer.Start();

  TString name;

  // Signal file, tree, and branch
  name = "AliESDs.root";
  TFile * fSig = TFile::Open(name.Data());
  TTree * tSig = (TTree*)fSig->Get("esdTree");

  AliESDEvent * esdSig = new AliESDEvent; // The signal ESD object is put here
  esdSig->ReadFromTree(tSig);

  // Run loader
  name = "galice.root";
  AliRunLoader* rlSig = AliRunLoader::Open(name.Data());

  // gAlice
  rlSig->LoadgAlice();
  gAlice = rlSig->GetAliRun();

  // Now load kinematics and event header
  rlSig->LoadKinematics();
  rlSig->LoadHeader();

  // Loop on events
  Long64_t nevSig = rlSig->GetNumberOfEvents();

  cout << nevSig << " events" << endl;

  for (Int_t iev=0; iev<nevSig; iev++) {
    cout << "Signal event " << iev << endl;

    // Get ESD
    tSig->GetEntry(iev);

    // Get MC
    rlSig->GetEvent(iev);

    // Particle stack
    AliStack * stackSig = rlSig->Stack();

    Int_t nrec = esdSig->GetNumberOfTracks();
 
 //-------------------------------------------------------------------------------------------------
    Int_t nstack = stackSig->GetNtrack();
    for(Int_t istack=0; istack < nstack; istack++){
     
      TParticle * part = stackSig->Particle(istack);

      // Loop on particles: check if the D* decay products are reconstructed
       
      if(!part) continue;
    
      if(TMath::Abs(part->GetPdgCode())== 413){
        cout<<"particle "<< istack << " is D*"<<endl;
      
       
        Int_t iDaughter1 = part->GetFirstDaughter();  //id of the Daughter = D^0
        if( iDaughter1<0) continue; 
        TParticle* daughter1 = stackSig->Particle(iDaughter1);
        if(!daughter1) continue; 
        cout<<"first daughter:  "<<daughter1->GetPdgCode()<<endl; 

        Int_t iDaughter2 = part->GetLastDaughter();  //id of the Daughter = pi+
        if( iDaughter2<0) continue; 
        TParticle* daughter2 = stackSig->Particle(iDaughter2);
        if(!daughter2) continue; 
        cout<<"last daughter: "<<daughter2->GetPdgCode()<<endl; 

        Int_t iD0 = -1;
        Int_t iPi = -1;
            
        if(TMath::Abs(daughter1->GetPdgCode())== 421){
          iD0=iDaughter1;
          iPi=iDaughter2;
        }
        else if(TMath::Abs(daughter2->GetPdgCode())== 421){
          iD0=iDaughter2; 
          iPi=iDaughter1;
        }
                   
        if (iD0<0)  continue;
 
        TParticle* secondmother = stackSig->Particle(iD0); 
        
        Int_t iD0Daughter1 = secondmother->GetFirstDaughter();
        TParticle*  D0Daughter1 = stackSig->Particle(iD0Daughter1);
        Int_t iD0Daughter2 = secondmother->GetLastDaughter();
        TParticle*  D0Daughter2 = stackSig->Particle(iD0Daughter2);
        
        for(Int_t irec=0; irec<nrec; irec++) {//loop on the ESDTree;
          AliESDtrack * track = esdSig->GetTrack(irec);
          UInt_t label = TMath::Abs(track->GetLabel());
          if(label<10000000) {  
            if(label == iPi)         cout<<label<< " We found the Pi from the D* decay"<<endl;
            if(label == iD0Daughter1) cout<<label<<" We found the K from the D0 decay"<<endl;
            if(label == iD0Daughter2) cout<<label<<" We found the Pi from the D0 decay"<<endl;
          }
        } 
      }
    }
  }

  // end loop on kine tree

  fSig->Close();

  timer.Stop();
  timer.Print();
}
