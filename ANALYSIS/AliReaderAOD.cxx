#include "AliReaderAOD.h"
//______________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class AliReaderAOD                                                         //
//                                                                            //
// Reader and Writer for AOD format.                                          //
// AODs are stored in a tree named by the variable fgkTreeName.               //
// There is stored 1 or 2 branches. Each of them stores AOD objects           //
// First branch is named by the variable fgkReconstructedDataBranchName       //
// ("reconstructed.") and keeps reconstructed data.                           //
// Second branch is called by the variable fgkSimulatedDataBranchName         //
// ("simulated.") and stores Monte carlo truth. If both branches are present  //
// AODs are parallel, i.e. nth particle in one branch corresponds to the nth  //
// particle in the other one.                                                 //
//                                                                            //
// Since we accept different formats of particles that are stored in AODs     //
// reader must take care of that fact: clean buffer if the next file contains //
// different particle type.                                                   //
//                                                                            //
// Piotr.Skowronski@cern.ch                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

ClassImp(AliReaderAOD)

#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include "AliAOD.h"


const TString AliReaderAOD::fgkTreeName("TAOD");
const TString AliReaderAOD::fgkReconstructedDataBranchName("reconstructed.");
const TString AliReaderAOD::fgkSimulatedDataBranchName("simulated.");

AliReaderAOD::AliReaderAOD(const Char_t* aodfilename):
 fFileName(aodfilename),
 fReadSim(kFALSE),
 fReadRec(kTRUE),
 fTree(0x0),
 fFile(0x0),
 fSimBuffer(0x0),
 fRecBuffer(0x0)
{
  //ctor
}
/********************************************************************/

AliReaderAOD::~AliReaderAOD()
{
//dtor
  if (fEventSim == fSimBuffer )
   {
    fEventSim = 0x0;
    fEventRec = 0x0;
   }
  delete fSimBuffer;
  delete fRecBuffer;
  
  delete fTree;
  delete fFile;
}
/********************************************************************/

void AliReaderAOD::Rewind()
{
//Rewinds reading
  delete fTree;
  fTree = 0x0;
  delete fFile;
  fFile = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;
}
/********************************************************************/
Int_t AliReaderAOD::ReadNext()
{
//Reads next event
  
  Info("ReadNext","Entered");
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
    {
      if (fFile == 0x0)
       {
         Int_t openfailed = OpenFile(fCurrentDir);//rl is opened here
         if (openfailed)
           {
             //Error("ReadNext","Error Occured while opening directory number %d",fCurrentDir);
             fCurrentDir++;
             continue;
           }
         fCurrentEvent = 0;
       }
      //Tree must exist because OpenFile would reuturn error in the other case
      if ( fCurrentEvent >= fTree->GetEntries() )
       {
         delete fTree;
         fTree = 0x0;
         delete fFile;
         fFile = 0x0;
         fSimBuffer = 0x0;
         fRecBuffer = 0x0;
         fCurrentDir++;
         continue;
       }
       
      Info("ReadNext","Getting event %d",fCurrentEvent);
      fTree->GetEvent(fCurrentEvent);
      Info("ReadNext","Getting event %d Done",fCurrentEvent);
      
      Int_t retval = 0;
      if (fReadRec && fReadSim)
       {
         retval = ReadRecAndSim();
       }
      else
       {
         if (fReadRec) retval = ReadRec();
         if (fReadSim) retval = ReadSim();
       } 

      fCurrentEvent++;
      if (retval == 0) fNEventsRead++;

      return retval;//success -> read one event
      
    }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array  
   
  return 1; //no more directories to read
  
  
}
/********************************************************************/

Int_t AliReaderAOD::ReadRecAndSim()
{
//Reads raconstructed and simulated data 

 Info("ReadRecAndSim","Found %d reconstructed tracks and %d simulated particles",
       fRecBuffer->GetNumberOfParticles(),fSimBuffer->GetNumberOfParticles());

 if (fCuts->GetEntriesFast() == 0x0)
  {//if there is no cuts we return pointer to the buffer
    if (fEventRec != fRecBuffer)
     {
       delete fEventRec;
       delete fEventSim;
     }
    fEventRec = fRecBuffer;//fEventRec is the pointer that the user gets when he asks about an event
    fEventSim = fSimBuffer;
  }
 else
  {//if there are cuts specified
    if ( (fEventRec == 0x0) || (fEventRec == fRecBuffer) )
     {//we need to create a new event, if it is not existing  or it is the same as branch buffer
       fEventRec = new AliAOD();
       fEventSim = new AliAOD();

       fEventRec->SetParticleClass( fRecBuffer->GetParticleClass() );
       fEventSim->SetParticleClass( fSimBuffer->GetParticleClass() );
     } 
    else
     {//or simply reset it in case it already exists
       fEventRec->Reset();
       fEventSim->Reset();
     }

    Int_t npart = fRecBuffer->GetNumberOfParticles();
    for (Int_t i = 0; i < npart; i++)
     {
       AliVAODParticle* prec = fRecBuffer->GetParticle(i);
       if (Rejected(prec)) continue;//we make cuts only on simulated data

       fEventRec->AddParticle(prec);
       fEventSim->AddParticle( fSimBuffer->GetParticle(i));
     }
  }

 Info("ReadRecAndSim","Read %d reconstructed tracks and %d simulated particles",
       fEventRec->GetNumberOfParticles(),fEventSim->GetNumberOfParticles());
 
 fTrackCounter->Fill(fEventRec->GetNumberOfParticles());
 
 return 0;
}
/********************************************************************/

Int_t AliReaderAOD::ReadRec()
{
//Reads reconstructed data only

 Info("ReadRec","Found %d reconstructed tracks",fRecBuffer->GetNumberOfParticles());

 if (fCuts->GetEntriesFast() == 0x0)
  {//if there is no cuts we return pointer to the buffer
    if (fEventRec != fRecBuffer)
     {
       delete fEventRec;
     }
    fEventRec = fRecBuffer;//fEventRec is the pointer that the user gets when he asks about an event
  }
 else
  {//if there are cuts specified
    if ( (fEventRec == 0x0) || (fEventRec == fRecBuffer) )
     {//we need to create a new event, if it is not existing  or it is the same as branch buffer
       fEventRec = new AliAOD();

       fEventRec->SetParticleClass( fRecBuffer->GetParticleClass() );
     } 
    else
     {//or simply reset it in case it already exists
       fEventRec->Reset();
     }

    Int_t npart = fRecBuffer->GetNumberOfParticles();
    for (Int_t i = 0; i < npart; i++)
     {
       AliVAODParticle* prec = fRecBuffer->GetParticle(i);
       if (Rejected(prec)) continue;//we make cuts only on simulated data

       fEventRec->AddParticle(prec);
     }
  }

 Info("ReadRec","Read %d reconstructed tracks",fEventRec->GetNumberOfParticles());
 fTrackCounter->Fill(fEventRec->GetNumberOfParticles());

 return 0;
}
/********************************************************************/

Int_t AliReaderAOD::ReadSim()
{
//Reads simulated data only

 Info("ReadSim","Found %d simulated particles",fSimBuffer->GetNumberOfParticles());

 if (fCuts->GetEntriesFast() == 0x0)
  {//if there is no cuts we return pointer to the buffer
    if (fEventSim != fSimBuffer)
     {
       delete fEventSim;
     }
    fEventSim = fSimBuffer;
  }
 else
  {//if there are cuts specified
    if ( (fEventSim == 0x0) || (fEventSim == fSimBuffer) )
     {//we need to create a new event, if it is not existing  or it is the same as branch buffer
       fEventSim = new AliAOD();

       fEventSim->SetParticleClass( fSimBuffer->GetParticleClass() );
     } 
    else
     {//or simply reset it in case it already exists
       fEventSim->Reset();
     }

    Int_t npart = fSimBuffer->GetNumberOfParticles();
    for (Int_t i = 0; i < npart; i++)
     {
       AliVAODParticle* prec = fSimBuffer->GetParticle(i);
       if (Rejected(prec)) continue;//we make cuts only on simulated data
       fEventSim->AddParticle(prec);
     }
  }

 Info("ReadSim","Read %d simulated particles",fEventSim->GetNumberOfParticles());
 fTrackCounter->Fill(fEventSim->GetNumberOfParticles());


 return 0;
}
/********************************************************************/

Int_t AliReaderAOD::OpenFile(Int_t n)
{
//opens fFile with  tree

// Info("ReadNext","Opening File %d",n);
 const TString dirname = GetDirName(n);
 if (dirname == "")
  {
    if (AliVAODParticle::GetDebug() > 2 )
      {
        Info("OpenFile","Got empty string as a directory name."); 
      }
   return 1;
  }
 
 TString filename = dirname +"/"+ fFileName;
 fFile = TFile::Open(filename.Data()); 
 if ( fFile == 0x0)
  {
    Error("OpenFile","Can't open fFile %s",filename.Data());
    return 2;
  }
 if (!fFile->IsOpen())
  {
    Error("OpenFile","Can't open fFile  %s",filename.Data());
    delete fFile;
    fFile = 0x0;
    return 3;
  }

 Info("ReadNext","File Is Opened, Getting the TREE");
 
 fTree = dynamic_cast<TTree*>(fFile->Get(fgkTreeName));
 if (fTree == 0x0)
  {
    if (AliVAODParticle::GetDebug() > 2 )
      {
        Info("ReadNext","Can not find TTree object named %s",fgkTreeName.Data());
      }
    delete fFile;
    fFile = 0x0;
    return 4;
  }

//  Info("ReadNext","Got TREE, Setting branch addresses");

  if (fReadRec)
   {
     TBranch* branch = fTree->GetBranch(fgkReconstructedDataBranchName);
     if (branch == 0x0)
      {
        Error("OpenFile","Can not find branch %s in file %s",
              fgkReconstructedDataBranchName.Data(),filename.Data());
        
        delete fTree;
        fTree = 0x0;
        delete fFile;
        fFile = 0x0;
        return 5;
      }
     fTree->SetBranchAddress(fgkReconstructedDataBranchName,&fRecBuffer);
   }
  

  if (fReadSim)
   {
     TBranch* branch = fTree->GetBranch(fgkSimulatedDataBranchName);
     if (branch == 0x0)
      {
        Error("OpenFile","Can not find branch %s in file %s",
              fgkSimulatedDataBranchName.Data(),filename.Data());
        
        delete fTree;
        fTree = 0x0;
        delete fFile;
        fFile = 0x0;
        return 6;
      }
     fTree->SetBranchAddress(fgkSimulatedDataBranchName,&fSimBuffer);
   }
//  Info("ReadNext","Got TREE, Addresses are set.");
//  Info("ReadNext","Quitting the method.");
  
  return 0;
 
}
/********************************************************************/

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

  TTree *tree = new TTree(fgkTreeName,"Tree with tracks");
  
  TBranch *recbranch = 0x0, *simbranch = 0x0;
  
  AliAOD* eventrec = new AliAOD();//must be created before Branch is called. Otherwise clones array is not splitted
  AliAOD* eventsim = new AliAOD();//AOD together with fParticles clones array knowing exact type of particles
  
  eventrec->SetParticleClassName(pclassname);
  eventsim->SetParticleClassName(pclassname);
 
  AliAOD* recbuffer = eventrec;
  AliAOD* simbuffer = eventsim;
  
  if (reader->ReadsRec()) recbranch = tree->Branch(fgkReconstructedDataBranchName,"AliAOD",&recbuffer,32000,99);
  if (reader->ReadsSim()) simbranch = tree->Branch(fgkSimulatedDataBranchName,"AliAOD",&simbuffer,32000,99);

  reader->Rewind();
  while (reader->Next() == kFALSE)
   {
     
     if (reader->ReadsRec())
      {//here we can get AOD that has different particle type
        AliAOD* event = reader->GetEventRec();
        if ( eventrec->GetParticleClass() != event->GetParticleClass() )
         {//if class type is not what what we whant we copy particles
           eventrec->CopyData(event);
           recbuffer = eventrec;
         }
        else
         {//else just pointer to event from input reader is passed
           recbuffer = event;
         } 
	recbuffer->GetParticle(0)->Print();
      }

     if (reader->ReadsSim())
      {
        AliAOD* event = reader->GetEventSim();
        if ( eventsim->GetParticleClass() != event->GetParticleClass() )
         {//if class type is not what what we whant we copy particles
           eventsim->CopyData(event);
           simbuffer = eventrec;
         }
        else
         {//else just pointer to event from input reader is passed
           simbuffer = event;
         } 
	simbuffer->GetParticle(0)->Print();
      }
     tree->Fill();
   }
  
  ::Info("AliReaderAOD::Write","Written %d events",tree->GetEntries());
  outfile->cd();
  tree->Write();

  delete eventsim;
  delete eventrec;
  
  delete tree;
  delete outfile;
  return 0; 
}

