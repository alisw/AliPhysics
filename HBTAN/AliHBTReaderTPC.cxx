#include "AliHBTReaderTPC.h"

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <Varargs.h>

#include <AliRun.h>
#include <AliLoader.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliTPCtrack.h>
#include <AliTPCParam.h>
#include <AliTPCtracker.h>
#include <AliTPCLoader.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

ClassImp(AliHBTReaderTPC)
//______________________________________________
//
// class AliHBTReaderTPC
//
//reader for TPC tracking
//needs galice.root, AliTPCtracks.root, AliTPCclusters.root
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//Piotr.Skowronski@cern.ch
AliHBTReaderTPC::AliHBTReaderTPC():
 fFileName("galice.root"),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untouched
  
}

AliHBTReaderTPC::AliHBTReaderTPC(const Char_t* galicefilename):
 fFileName(galicefilename),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untouched
  
}
/********************************************************************/
AliHBTReaderTPC::AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename):
 AliHBTReader(dirs), 
 fFileName(galicefilename),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE)
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untached
  
}
/********************************************************************/

AliHBTReaderTPC::~AliHBTReaderTPC()
{
 //desctructor
   delete fRunLoader;
}
/********************************************************************/
 
void AliHBTReaderTPC::Rewind()
{
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;
}
/********************************************************************/

Int_t AliHBTReaderTPC::ReadNext()
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");

 
  TObjArray *tarray = new TObjArray(5000); //cotainer for tpc tracks
  tarray->SetOwner(); //set the ownership of the objects it contains
                      //when array is is deleted or cleared all objects 
                      //that it contains are deleted

  if (fParticlesEvent == 0x0)  fParticlesEvent = new AliHBTEvent();
  if (fTracksEvent == 0x0)  fTracksEvent = new AliHBTEvent();

  fParticlesEvent->Reset();
  fTracksEvent->Reset();

  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
      
    if (fRunLoader == 0x0) 
      if (OpenNextSession()) continue;//directory counter is increased inside in case of error

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
    TTree *tracktree = fTPCLoader->TreeT();//get the tree 
    if (!tracktree) //check if we got the tree
      {//if not return with error
         Error("ReadNext","Can't get a tree with TPC tracks !\n"); 
         continue;
      }
   
    TBranch *trackbranch=tracktree->GetBranch("tracks");//get the branch with tracks
    if (!trackbranch) ////check if we got the branch
      {//if not return with error
        Error("ReadNext","Can't get a branch with TPC tracks !\n"); 
        continue;
      }
    Int_t ntpctracks=(Int_t)tracktree->GetEntries();//get number of TPC tracks 
    Info("ReadNext","Found %d TPC tracks.",ntpctracks);
    //Copy tracks to array

    AliTPCtrack *iotrack=0;

    for (Int_t i=0; i<ntpctracks; i++) //loop over all tpc tracks
     {
       iotrack=new AliTPCtrack;   //create new tracks
       trackbranch->SetAddress(&iotrack); //tell the branch ehere to put track data from tree(file)
       tracktree->GetEvent(i); //stream track i to the iotrack
       tarray->AddLast(iotrack); //put the track in the array
     }

    Double_t xk;
    Double_t par[5];
    Float_t phi, lam, pt;//angles and transverse momentum
    Int_t label; //label of the current track

    AliStack* stack = fRunLoader->Stack();
    if (stack == 0x0)
     {
       Error("ReadNext","Can not get stack for current event",fCurrentEvent);
       fCurrentEvent++;
       continue;
     }
    stack->Particles();

    for (Int_t i=0; i<ntpctracks; i++) //loop over all good tracks
     { 
       iotrack = (AliTPCtrack*)tarray->At(i);
       label = iotrack->GetLabel();

       if (label < 0) continue;

       TParticle *p = (TParticle*)stack->Particle(label);

       if(p == 0x0) continue; //if returned pointer is NULL
       if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)

       if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                   //if not take next partilce

       AliHBTParticle* part = new AliHBTParticle(*p,i);
       if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                               //if it does not delete it and take next good track

//       iotrack->PropagateTo(3.,0.0028,65.19);
//       iotrack->PropagateToVertex(36.66,1.2e-3);//orig
       iotrack->PropagateToVertex(50.,0.0353);
       
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

       AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(),i, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
       if(Pass(track))//check if meets all criteria of any of our cuts
                    //if it does not delete it and take next good track
        { 
          delete track;
          delete part;
          continue;
        }
       fParticlesEvent->AddParticle(part);
       fTracksEvent->AddParticle(track);
      }
      
    Info("ReadNext","Read %d tracks and %d particles from event %d (event %d in dir %d).",
            fParticlesEvent->GetNumberOfParticles(), fTracksEvent->GetNumberOfParticles(),
            fNEventsRead,fCurrentEvent,fCurrentDir);
    
    fCurrentEvent++;
    fNEventsRead++;
    delete tarray;
    return 0;
   }while(fCurrentDir < GetNumberOfDirs());

  delete tarray;
  return 1;
 }
/********************************************************************/

Int_t AliHBTReaderTPC::OpenNextSession()
{
  TString filename = GetDirName(fCurrentDir);
  if (filename.IsNull())
   {
     DoOpenError("Can not get directory name");
     return 1;
   }
  filename = filename +"/"+ fFileName;
  fRunLoader = AliRunLoader::Open(filename,AliConfig::fgkDefaultEventFolderName);
  if( fRunLoader == 0x0)
   {
     DoOpenError("Can not open session.");
     return 1;
   }

  fRunLoader->LoadHeader();
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
  fTPCLoader = (AliTPCLoader*)fRunLoader->GetLoader("TPCLoader");
  if ( fTPCLoader == 0x0)
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
     fRunLoader->LoadgAlice();
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


  fRunLoader->CdGAFile();
  AliTPCParam *TPCParam= (AliTPCParam*)gDirectory->Get("75x40_100x60");
  if (!TPCParam) 
   {
    TPCParam=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
    if (!TPCParam) 
     { 
       DoOpenError("TPC parameters have not been found !\n");
       return 1;
     }
   }

  if (fTPCLoader->LoadTracks())
   {
     DoOpenError("Error occured while loading TPC tracks.");
     return 1;
   }
  
  fCurrentEvent = 0;
  return 0;
}
/********************************************************************/

void AliHBTReaderTPC::DoOpenError( const char *va_(fmt), ...)
{
  // Does error display and clean-up in case error caught on Open Next Session

   va_list ap;
   va_start(ap,va_(fmt));
   Error("OpenNextSession", va_(fmt), ap);
   va_end(ap);
   
   delete fRunLoader;
   fRunLoader = 0x0;
   fTPCLoader = 0x0;
   fCurrentDir++;
}

/********************************************************************/
/********************************************************************/
/********************************************************************/

