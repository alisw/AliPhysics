#include "AliHBTReaderKineTree.h"

#include <Riostream.h>
//#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliStack.h>
#include <AliHeader.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

ClassImp(AliHBTReaderKineTree)
/**********************************************************/

AliHBTReaderKineTree::
AliHBTReaderKineTree():fFileName("galice.root")
{
  fParticles = new AliHBTRun();
  fIsRead = kFALSE;
}

AliHBTReaderKineTree::
AliHBTReaderKineTree(TString& fname):fFileName(fname)
{
  fParticles = new AliHBTRun();
  fIsRead = kFALSE;
}


/**********************************************************/
AliHBTReaderKineTree::
AliHBTReaderKineTree(TObjArray* dirs,const Char_t *filename):AliHBTReader(dirs),fFileName(filename)
{
  fParticles = new AliHBTRun();
  fIsRead = kFALSE;
}
/**********************************************************/

AliHBTEvent* AliHBTReaderKineTree::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
   if (!fIsRead) 
    if(Read(fParticles,0x0))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }

   return fParticles->GetEvent(n);
 }

/********************************************************************/

Int_t AliHBTReaderKineTree::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead)
    if(Read(fParticles,0x0))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
   return fParticles->GetNumberOfEvents();
 }


/**********************************************************/
Int_t AliHBTReaderKineTree::
Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 cout<<"AliHBTReaderKineTree::Read()"<<endl;
 if (!particles) //check if an object is instatiated
   {
     Error("Read"," particles object must instatiated before passing it to the reader");
     return 1;
   }  
 particles->Reset();//clear runs == delete all old events
 
 Int_t Ndirs;//total number of directories specified in fDirs array
 Int_t Nevents; //number of events read in current directory
 Int_t totalNevents = 0; //total number of read events 
 Int_t currentdir = 0; //number of current directory name is fDirs array

 if (fDirs) //if array with directories is supplied by user
  {
    Ndirs = fDirs->GetEntries(); //get the number if directories
  }
 else
  {
    Ndirs = 0; //if the array is not supplied read only from current directory
  }

 do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
  { 
    cout<<"________________________________________________________\n";
    TFile* fileK = OpenFile(currentdir);
    TTree* treeE = (TTree*)fileK->Get("TE");
    if(treeE) 
     {
      Nevents = (Int_t)treeE->GetEntries();
      cout<<"Found "<<Nevents<<" in directory "<<GetDirName(currentdir)<<endl;
      fileK->Delete("TE");
     }
    else
     {
      Error("Read","Cannot find Header Tree (TE)");
      Nevents = 0;
     }
    
    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
//      cout<<"processing event "<<currentEvent<<" in current dir."<<endl;
      AliStack* stack = GetStack(currentEvent,fileK);
      if (!stack)
       {
         Error("Read","Can not get stack for event %d",currentEvent);
         continue;
       }
      Int_t npart = stack->GetNtrack();
      Int_t nnn = 0;
      for (Int_t i = 0;i<npart; i++)
       {
         
         TParticle * p = stack->Particle(i);
         if (p->GetFirstMother() >= 0) continue;
         
         if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                              //if not take next partilce
         
         AliHBTParticle* part = new AliHBTParticle(*p,i);
         if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track
         particles->AddParticle(totalNevents,part);//put track and particle on the run

         if ( (nnn%100) == 0) 
          {
            cout<<nnn<<"\r";
            fflush(0);
          }
         nnn++;
       }
      cout<<"Total read "<<nnn<<endl;
      delete stack;
      totalNevents++;
     }
    delete gAlice;
    gAlice = 0;
    if (fileK) 
     {
      fileK->Close();
      delete fileK;
      fileK = 0;
     }
    
    currentdir++;
  }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array
 fIsRead = kTRUE;
 return 0;
}

/**********************************************************/
AliStack* AliHBTReaderKineTree::GetStack(Int_t n, TFile* file)
{
  AliHeader *header =  new AliHeader();
  TTree* treeE = (TTree*)file->Get("TE");
  if(treeE) 
   {
    treeE->SetBranchAddress("Header", &header);
    if (!treeE->GetEntry(n)) 
     {
       Error("GetEvent","Cannot find event:%dn",n);
       delete header;
       treeE->Delete("TE");
       return 0x0;
     }
   }  
  else 
   {
     Error("GetStack","Cannot find Header Tree (TE)n");
     delete header;
     return 0x0;
   }

  AliStack* stack = header->Stack();
  if (stack) 
   {
    if (!stack->GetEvent(n)) 
     {
       delete header;
       treeE->Delete("TE");
       return 0x0;
     }
   }
  delete header;
  treeE->Delete("TE");
  return stack;
}

TFile* AliHBTReaderKineTree::OpenFile(Int_t n)
{
//opens file with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename = dirname +"/"+ fFileName;
 TFile *ret = TFile::Open(filename.Data()); 

 if ( ret == 0x0)
  {
    Error("OpenFiles","Can't open file %s",filename.Data());
    return 0x0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open file  %s",filename.Data());
    return 0x0;
  }
 
 return ret;
}
