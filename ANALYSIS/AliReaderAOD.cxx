#include "AliReaderAOD.h"

ClassImp(AliReaderAOD)

#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include "AliAOD.h"

Int_t AliReaderAOD::WriteAOD(AliReader* reader, const char* outfilename, const char* pclassname,  Bool_t /*multcheck*/)
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
  
  
  AliAOD* eventrec = new AliAOD();//must be created before Branch is called. Otherwise clones array is not splitted
  AliAOD* eventsim = new AliAOD();//AOD together with fParticles clones array knowing exact type of particles
  
  eventrec->SetParticleClassName(pclassname);
  eventsim->SetParticleClassName(pclassname);
  
  if (reader->ReadsRec()) recbranch = tree->Branch("reconstructed","AliAOD",&eventrec,32000,99);
  if (reader->ReadsSim()) simbranch = tree->Branch("simulated","AliAOD",&eventsim,32000,99);

  delete eventsim;
  delete eventrec;
  
  reader->Rewind();
  while (reader->Next() == kFALSE)
   {
     
     if (reader->ReadsRec())
      {
        eventrec = reader->GetEventRec();
        recbranch->SetAddress(&eventrec);
      }

     if (reader->ReadsSim())
      {
        eventsim = reader->GetEventSim();
        simbranch->SetAddress(&eventsim);
      }
     eventrec->GetParticle(0)->Print();
     eventsim->GetParticle(0)->Print();
     tree->Fill();
   }
  
  ::Info("AliReaderAOD::Write","Written %d events",tree->GetEntries());
  outfile->cd();
  tree->Write();
  delete tree;
  delete outfile;
  return 0; 
}

