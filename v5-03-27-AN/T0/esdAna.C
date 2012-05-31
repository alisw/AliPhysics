#define esdAna_cxx
// The class definition in esdAna.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("esdAna.C")
// Root > T->Process("esdAna.C","some options")
// Root > T->Process("esdAna.C+")
//

#include "esdAna.h"
#include <TH1.h>

void esdAna::Begin(TTree *)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  
   h1 = new TH1F("hRealVertex","Primary vertex",100,-20,20);
   h3 = new TH1F("hT0vertex","T0vertex",100,-20,20);
   h2 = new TH1F("hT0start","T0 start time",100,12400,12600);
 TString option = GetOption();   
}

void esdAna::SlaveBegin(TTree *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  
     Init(tree);

   TString option = GetOption();

}

Bool_t esdAna::Process(Long64_t entry)
{

   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.

   // WARNING when a selector is used with a TChain, you must use
   //  the pointer to the current TTree to call GetEntry(entry).
   //  The entry is always the local entry number in the current tree.
   //  Assuming that fChain is the pointer to the TChain being processed,
   //  use fChain->GetTree()->GetEntry(entry).
  //  fChain->GetTree()->GetEntry(entry);
  b_ESD_fEventNumber->GetEvent(entry);
  b_ESD_fPrimaryVertex_fPosition->GetEntry(entry);
  b_ESD_fT0zVertex->GetEntry(entry);
   b_ESD_fT0timeStart->GetEntry(entry);
   printf("Processing Entry %lld %d  %f %f %f \n",entry,fEventNumber,fT0zVertex,fPrimaryVertex_fPosition[2], fT0timeStart );

   h2->Fill(fT0timeStart);
   h1->Fill(fPrimaryVertex_fPosition[2]);
   h3->Fill(fT0zVertex/2.);
   

   return kTRUE;
}

void esdAna::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   fOutput->Add(h1) ;
   fOutput->Add(h2) ;
   fOutput->Add(h3) ;
}

void esdAna::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  Float_t mean = h2->GetMean();
  printf ("mean time T0 ps %f ",mean);
  if (mean >  12600 || mean <12400 ) 
    printf (" !!!!!!!!!!-----events sample is WRONG - T0 unreal -------");  
   hfile = TFile::Open("esdAna.root","RECREATE");
   TFile * file = TFile::Open("esdAna.root", "RECREATE");
   h1->Write();
   h2->Write();
   h3->Write();
   file->Close();
}
