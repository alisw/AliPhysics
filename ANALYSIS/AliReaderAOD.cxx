#include "AliReaderAOD.h"

ClassImp(AliReaderAOD)

#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include "AliAOD.h"


const TString AliReaderAOD::fgkTreeName("TAOD");
const TString AliReaderAOD::fgkRecosntructedDataBranchName("reconstructed.");
const TString AliReaderAOD::fgkSimulatedDataBranchName("simulated.");

AliReaderAOD::AliReaderAOD(const Char_t* aodfilename):
 fFileName(aodfilename),
 fReadSim(kFALSE),
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
         Int_t opened = OpenFile(fCurrentDir);//rl is opened here
         if (opened)
           {
             //Error("ReadNext","Error Occured while opening directory number %d",fCurrentDir);
             fCurrentDir++;
             continue;
           }
         fCurrentEvent = 0;
       }
      Info("ReadNext","Getting event %d",fCurrentEvent);
      fTree->GetEvent(fCurrentEvent);
      Info("ReadNext","Getting event %d Done",fCurrentEvent);
      
      //Temporary testing sollution
      fEventSim = fSimBuffer;
      fEventRec = fRecBuffer;
      
      fCurrentEvent++;
      fNEventsRead++;

      if (fTree)
       {
         if ( fCurrentEvent >= fTree->GetEntries() )
          {
            delete fTree;
            fTree = 0x0;
            delete fFile;
            fFile = 0x0;
            fSimBuffer = 0x0;
            fRecBuffer = 0x0;
            fCurrentDir++;
          } 
       }


      return 0;//success -> read one event
      
    }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array  
   
  return 1; //no more directories to read
  
  
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
    fCurrentDir++;
    delete fFile;//we have to assume there is no more ESD objects in the fFile
    fFile = 0x0;
    return 4;
  }

//  Info("ReadNext","Got TREE, Setting branch addresses");

  fTree->SetBranchAddress(fgkSimulatedDataBranchName,&fSimBuffer);
  fTree->SetBranchAddress(fgkRecosntructedDataBranchName,&fRecBuffer);
  
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
  
  if (reader->ReadsRec()) recbranch = tree->Branch(fgkRecosntructedDataBranchName,"AliAOD",&recbuffer,32000,99);
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
      }
     recbuffer->GetParticle(0)->Print();
     simbuffer->GetParticle(0)->Print();
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

