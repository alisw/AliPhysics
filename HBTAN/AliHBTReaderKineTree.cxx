#include "AliHBTReaderKineTree.h"

#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliHeader.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

ClassImp(AliHBTReaderKineTree)
/**********************************************************/
const TString AliHBTReaderKineTree::fgkEventFolderName("HBTReaderKineTree");

AliHBTReaderKineTree::AliHBTReaderKineTree():fFileName("galice.root")
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
    AliRunLoader * rl = OpenFile(currentdir);

    if (rl == 0x0)
     {
      Error("Read","Cannot open session");
      return 2;
     }
    rl->LoadHeader();
    rl->LoadKinematics("READ");
    Nevents = rl->GetNumberOfEvents();
    
    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
//      cout<<"processing event "<<currentEvent<<" in current dir."<<endl;
      rl->GetEvent(currentEvent);
      AliStack* stack = rl->Stack();
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
//         if (p->GetFirstMother() >= 0) continue; do not apply with pythia etc
         
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
      totalNevents++;
     }
     delete rl;
     currentdir++;
  }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array
 fIsRead = kTRUE;
 return 0;
}

AliRunLoader* AliHBTReaderKineTree::OpenFile(Int_t n)
{
//opens file with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename = dirname +"/"+ fFileName;

 AliRunLoader *ret = AliRunLoader::Open(filename.Data(),fgkEventFolderName,"READ"); 

 if ( ret == 0x0)
  {
    Error("OpenFiles","Can't open session from file %s",filename.Data());
    return 0x0;
  }
 
 return ret;
}
