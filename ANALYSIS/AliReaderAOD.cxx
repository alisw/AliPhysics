#include "AliReaderAOD.h"

ClassImp(AliReaderAOD)

#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include "AliAOD.h"

Int_t AliReaderAOD::WriteAOD(AliReader* reader, const char* outfilename, Bool_t /*multcheck*/)
{
//reads tracks from runs and writes them to file
  ::Info("AliReaderAOD::Write","________________________________________________________");
  ::Info("AliReaderAOD::Write","________________________________________________________");
  ::Info("AliReaderAOD::Write","________________________________________________________");
  
  if (reader == 0x0)
   {
     ::Error("AliReaderAOD::Write","Input Reader is NULL");
     return -1;
   }
  TFile *outfile = TFile::Open(outfilename,"recreate");
  if (outfile == 0x0)
   {
     ::Error("AliReaderAOD::Write","Can not open output file %s",outfilename);
     return -1;
   }

  TTree *tree = new TTree("TAOD","Tree with tracks");
  
  TBranch *recbranch = 0x0, *simbranch = 0x0;

  AliAOD* eventsim = new AliAOD();
  AliAOD* eventrec = new AliAOD;
  
  eventsim->SetParticleClassName("AliAODParticle");
  eventrec->SetParticleClassName("AliAODParticle");
  
  if (reader->ReadsSim()) simbranch = tree->Branch("simulated","AliAOD",&eventsim,32000,99);
  if (reader->ReadsRec()) recbranch = tree->Branch("reconstructed","AliAOD",&eventrec,32000,99);

  reader->Rewind();
  while (reader->Next() == kFALSE)
   {
     
     if (reader->ReadsSim())
      {
        eventsim = reader->GetEventSim();
//        simbranch->SetAddress(&eventsim);
      }
 
     if (reader->ReadsRec())
      {
        eventrec = reader->GetEventRec();
//        recbranch->SetAddress(&eventrec);
      }
     tree->Fill();
     tree->Print();
   }
  
  ::Info("AliReaderAOD::Write","Written %d events",tree->GetEntries());
  outfile->cd();
  tree->Write();
  delete tree;
  delete outfile;
  return 0; 
}

