#include "AliHBTReaderKineTree.h"

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

AliHBTReaderKineTree::AliHBTReaderKineTree():
 fFileName("galice.root"),
 fRunLoader(0x0)
{
  //ctor
}
/**********************************************************/

AliHBTReaderKineTree::AliHBTReaderKineTree(TString& fname):
 fFileName(fname),
 fRunLoader(0x0)
{
  //ctor
}
/**********************************************************/

AliHBTReaderKineTree::AliHBTReaderKineTree(TObjArray* dirs,const Char_t *filename):
 AliHBTReader(dirs),
 fFileName(filename),
 fRunLoader(0x0)
{
  //ctor
}
/**********************************************************/

AliHBTReaderKineTree::~AliHBTReaderKineTree()
{
  //dtor
  delete fRunLoader;
}
/**********************************************************/

void AliHBTReaderKineTree::Rewind()
{
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;  
}
/**********************************************************/

Int_t AliHBTReaderKineTree::ReadNext()
{
 //Reads Kinematics Tree
  
 Info("Read","");
 if (fParticlesEvent == 0x0)  fParticlesEvent = new AliHBTEvent();
 fParticlesEvent->Reset();

 do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
  { 
    if (fRunLoader == 0x0) 
      if (OpenNextFile()) 
        { 
          fCurrentDir++;
          continue;
        }
    
    if (fCurrentEvent == fRunLoader->GetNumberOfEvents())
     {
       //read next directory
       delete fRunLoader;//close current session
       fRunLoader = 0x0;//assure pointer is null
       fCurrentDir++;//go to next dir
       continue;//directory counter is increased inside in case of error
     }
     
    Info("ReadNext","Reading Event %d",fCurrentEvent);
    
    fRunLoader->GetEvent(fCurrentEvent);
    
    AliStack* stack = fRunLoader->Stack();
    if (!stack)
     {
       Error("ReadNext","Can not get stack for event %d",fCurrentEvent);
       continue;
     }
    Int_t npart = stack->GetNtrack();
    for (Int_t i = 0;i<npart; i++)
      {
         TParticle * p = stack->Particle(i);
//         if (p->GetFirstMother() >= 0) continue; do not apply with pythia etc
         
         if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                             //if not take next partilce
         
         AliHBTParticle* part = new AliHBTParticle(*p,i);
         if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track
         fParticlesEvent->AddParticle(part);//put particle in event
      }
    Info("ReadNext","Read %d particles from event %d (event %d in dir %d).",
                     fParticlesEvent->GetNumberOfParticles(),
                     fNEventsRead,fCurrentEvent,fCurrentDir);
      
    fCurrentEvent++;
    fNEventsRead++;
    return 0;
  }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array
  
 return 1;
}
/**********************************************************/

Int_t AliHBTReaderKineTree::OpenNextFile()
{
//opens file with kine tree
 Info("OpenNextFile","________________________________________________________");
 
 const TString& dirname = GetDirName(fCurrentEvent);
 if (dirname == "")
  {
   Error("OpenNextFile","Can not get directory name");
   return 1;
  }
 TString filename = dirname +"/"+ fFileName;

 fRunLoader = AliRunLoader::Open(filename.Data(),fgkEventFolderName,"READ"); 

 if ( fRunLoader == 0x0)
  {
    Error("OpenNextFile","Can't open session from file %s",filename.Data());
    return 0x0;
  }
  
 if (fRunLoader->GetNumberOfEvents() <= 0)
  {
    Error("OpenNextFile","There is no events in this directory.");
    delete fRunLoader;
    fRunLoader = 0x0;
    return 1;
  }
  
 if (fRunLoader->LoadKinematics())
  {
    Error("OpenNextFile","Error occured while loading kinematics.");
    return 1;
  }
  
 fCurrentEvent = 0;
 return 0;
}
