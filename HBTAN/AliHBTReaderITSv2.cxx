#include "AliHBTReaderITSv2.h"

#include <Varargs.h>

#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliLoader.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliITStrackV2.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderITSv2)

AliHBTReaderITSv2::AliHBTReaderITSv2():
 fFileName("galice.root"),
 fRunLoader(0x0),
 fITSLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
}
/********************************************************************/

AliHBTReaderITSv2::AliHBTReaderITSv2(const Char_t* galicefilename):
 fFileName(galicefilename),
 fRunLoader(0x0),
 fITSLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
}
/********************************************************************/

AliHBTReaderITSv2::AliHBTReaderITSv2(TObjArray* dirs, const Char_t* galicefilename): 
 AliHBTReader(dirs),
 fFileName(galicefilename),
 fRunLoader(0x0),
 fITSLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
  
}
/********************************************************************/
 
void AliHBTReaderITSv2::Rewind()
{
  //rewinds reading
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;
}
/********************************************************************/

AliHBTReaderITSv2::~AliHBTReaderITSv2()
{
  //dtor
  delete fRunLoader;
}
/********************************************************************/
 
Int_t AliHBTReaderITSv2::ReadNext()
{
//reads data from next event
  
 register Int_t i = 0; //iterator
 
// AliITStrackerV2 *tracker; // ITS tracker - used for cooking labels
 TTree *tracktree = 0x0; // tree for tracks
 AliITStrackV2 *iotrack = 0x0;
 
 Double_t xk;
 Double_t par[5]; //Kalman track parameters
 Float_t phi, lam, pt;//angles and transverse momentum
 Int_t label; //label of the current track
 
 
 if (fParticlesEvent == 0x0)  fParticlesEvent = new AliHBTEvent();
 if (fTracksEvent == 0x0)  fTracksEvent = new AliHBTEvent();

 fParticlesEvent->Reset();
 fTracksEvent->Reset();
 do //do while is good even if Ndirs==0 (than read from current directory)
   {
    if (fRunLoader == 0x0) 
      if (OpenNextFile()) continue;//directory counter is increased inside in case of error

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
    
    tracktree=fITSLoader->TreeT();
    if (!tracktree) 
     {
       Error("ReadNext","Can't get a tree with ITS tracks"); 
       fCurrentEvent++;
       continue;
     }
      
    TBranch *tbranch=tracktree->GetBranch("tracks");
    if (!tbranch) 
     {
       Error("ReadNext","Can't get a branch with ITS tracks"); 
       fCurrentEvent++;
       continue;
     }

    AliStack* stack = fRunLoader->Stack();
    if (stack == 0x0)
     {
       Error("ReadNext","Can not get stack for current event",fCurrentEvent);
       fCurrentEvent++;
       continue;
     }
     
    //must be here because on the beginning conv. const. is not set yet 
    if (iotrack == 0x0) iotrack = new AliITStrackV2(); //buffer track for reading data from tree
    
    Int_t ntr = (Int_t)tracktree->GetEntries();
    
    for (i=0; i < ntr; i++) //loop over all tpc tracks
     { 
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);

       label=iotrack->GetLabel();
       if (label < 0) 
        {
          continue;
        }

       TParticle *p = stack->Particle(label);
       if(p == 0x0) continue; //if returned pointer is NULL
       if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)

       if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                           //if not take next partilce

       AliHBTParticle* part = new AliHBTParticle(*p,i);
       if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                               //if it does not delete it and take next good track

       iotrack->PropagateTo(3.,0.0028,65.19);
       iotrack->PropagateToVertex();

       iotrack->GetExternalParameters(xk,par);     //get properties of the track
       phi=TMath::ASin(par[2]) + iotrack->GetAlpha(); 
       if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
       if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
       lam=par[3]; 
       pt=1.0/TMath::Abs(par[4]);

       Double_t tpx = pt * TMath::Cos(phi); //track x coordinate of momentum
       Double_t tpy = pt * TMath::Sin(phi); //track y coordinate of momentum
       Double_t tpz = pt * lam; //track z coordinate of momentum

       Double_t mass = p->GetMass();
       Double_t tEtot = TMath::Sqrt( tpx*tpx + tpy*tpy + tpz*tpz + mass*mass);//total energy of the track

       AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(), i, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
       if(Pass(track))//check if meets all criteria of any of our cuts
                      //if it does not delete it and take next good track
        { 
         delete track;
         delete part;
         continue;
        }
        
       fParticlesEvent->AddParticle(part);
       fTracksEvent->AddParticle(track);
     }//end of loop over tracks in the event
       
    Info("ReadNext","Read %d tracks and %d particles from event %d (event %d in dir %d).",
            fParticlesEvent->GetNumberOfParticles(), fTracksEvent->GetNumberOfParticles(),
            fNEventsRead,fCurrentEvent,fCurrentDir);
     
    fCurrentEvent++;
    fNEventsRead++;
    delete iotrack;
    return 0;
   }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array

 delete iotrack;
 return 1;
}

/********************************************************************/
Int_t AliHBTReaderITSv2::OpenNextFile()
{
  //opens next file
  TString filename = GetDirName(fCurrentDir);
  if (filename.IsNull())
   {
     DoOpenError("Can not get directory name");
     return 1;
   }
  filename = filename +"/"+ fFileName;
  fRunLoader = AliRunLoader::Open(filename,AliConfig::GetDefaultEventFolderName());
  if( fRunLoader == 0x0)
   {
     DoOpenError("Can not open session.");
     return 1;
   }

  if (fRunLoader->GetNumberOfEvents() <= 0)
   {
     DoOpenError("There is no events in this directory.");
     return 1;
   }

  if (fRunLoader->LoadKinematics())
   {
     DoOpenError("Error occured while loading kinematics.");
     return 1;
   }
  fITSLoader = fRunLoader->GetLoader("ITSLoader");
  if ( fITSLoader == 0x0)
   {
     DoOpenError("Exiting due to problems with opening files.");
     return 1;
   }
   
  Info("OpenNextSession","________________________________________________________");
  Info("OpenNextSession","Found %d event(s) in directory %s",
        fRunLoader->GetNumberOfEvents(),GetDirName(fCurrentDir).Data());
  Float_t mf;
  if (fUseMagFFromRun)
   {
     if (fRunLoader->LoadgAlice())
      {
        DoOpenError("Error occured while loading AliRun.");
        return 1;
      }
     mf = fRunLoader->GetAliRun()->Field()->SolenoidField();
     Info("OpenNextSession","Setting Magnetic Field from run: B=%fT",mf/10.);
     fRunLoader->UnloadgAlice();
   }
  else
   {
     Info("OpenNextSession","Setting Own Magnetic Field: B=%fT",fMagneticField);
     if (fMagneticField == 0x0)
      {
        Fatal("OpenNextSession","Magnetic field can not be 0.");
        return 1;//pro forma
      }
     mf = fMagneticField*10.;
   }
  AliKalmanTrack::SetConvConst(1000/0.299792458/mf);

  if (fITSLoader->LoadTracks())
   {
     DoOpenError("Error occured while loading TPC tracks.");
     return 1;
   }
  
  fCurrentEvent = 0;
  return 0;
}
/********************************************************************/

void AliHBTReaderITSv2::DoOpenError( const char *va_(fmt), ...)
{
  // Does error display and clean-up in case error caught on Open Next Session

   va_list ap;
   va_start(ap,va_(fmt));
   Error("OpenNextFile", va_(fmt), ap);
   va_end(ap);
   
   delete fRunLoader;
   fRunLoader = 0x0;
   fITSLoader = 0x0;
   fCurrentDir++;
}

/********************************************************************/
