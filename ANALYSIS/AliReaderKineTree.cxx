#include "AliReaderKineTree.h"
//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliReaderKineTree
//
// Reader for Kinematics
//
// Piotr.Skowronski@cern.ch
//
/////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TParticle.h>

#include <AliRunLoader.h>
#include <AliStack.h>

#include "AliAOD.h"
#include "AliAODParticle.h"

ClassImp(AliReaderKineTree)
/**********************************************************/
const TString AliReaderKineTree::fgkEventFolderName("ReaderKineTree");

AliReaderKineTree::AliReaderKineTree():
 fFileName("galice.root"),
 fRunLoader(0x0)
{
  //ctor
}
/**********************************************************/

AliReaderKineTree::AliReaderKineTree(TString& fname):
 fFileName(fname),
 fRunLoader(0x0)
{
  //ctor
}
/**********************************************************/

AliReaderKineTree::AliReaderKineTree(TObjArray* dirs,const Char_t *filename):
 AliReader(dirs),
 fFileName(filename),
 fRunLoader(0x0)
{
  //ctor
}

/**********************************************************/
AliReaderKineTree::AliReaderKineTree(const AliReaderKineTree& in):
 AliReader(in),
 fFileName(in.fFileName),
 fRunLoader(0x0)
{
  //cpy ctor
}

/**********************************************************/

AliReaderKineTree::~AliReaderKineTree()
{
  //dtor
  delete fRunLoader;
}
/**********************************************************/
AliReaderKineTree& AliReaderKineTree::operator=(const AliReaderKineTree& in)
{
//Assiment operator
  if (this == &in) return *this;
  AliReader::operator=(in);
  delete fRunLoader;
  fRunLoader = 0x0;
  return * this;
}
/**********************************************************/

void AliReaderKineTree::Rewind()
{
//Rewinds to the beginning
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;  
}
/**********************************************************/

Int_t AliReaderKineTree::ReadNext()
{
 //Reads Kinematics Tree
  
 Info("Read","");
 if (fEventSim == 0x0)  
  {
    fEventSim = new AliAOD();
  } 
 fEventSim->Reset();

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
         
         if(Rejected(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                             //if not take next partilce
         
         AliAODParticle* part = new AliAODParticle(*p,i);
         if(Rejected(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track
         fEventSim->AddParticle(part);//put particle in event
      }
    Info("ReadNext","Read %d particles from event %d (event %d in dir %d).",
                     fEventSim->GetNumberOfParticles(),
                     fNEventsRead,fCurrentEvent,fCurrentDir);
      
    fCurrentEvent++;
    fNEventsRead++;
    return 0;
  }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array
  
 return 1;
}
/**********************************************************/

Int_t AliReaderKineTree::OpenNextFile()
{
//opens file with kine tree
 Info("OpenNextFile","________________________________________________________");
 
 const TString& dirname = GetDirName(fCurrentDir);
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
    return 1;
  }
  
 if (fRunLoader->GetNumberOfEvents() <= 0)
  {
    Error("OpenNextFile","There is no events in this directory.");
    delete fRunLoader;
    fRunLoader = 0x0;
    return 2;
  }
  
 if (fRunLoader->LoadKinematics())
  {
    Error("OpenNextFile","Error occured while loading kinematics.");
    return 3;
  }
  
 fCurrentEvent = 0;
 return 0;
}
