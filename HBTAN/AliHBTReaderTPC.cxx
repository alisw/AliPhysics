#include "AliHBTReaderTPC.h"
//______________________________________________
//
// class AliHBTReaderTPC
//
// reader for TPC tracks
// needs galice.root
// just to shut up coding conventions checker
// 
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////
#include <TTree.h>
#include <TParticle.h>
#include <TH1.h>


#include <AliRun.h>
#include <AliLoader.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliTPCtrack.h>
#include <AliTPCLoader.h>

#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTTrackPoints.h"
#include "AliHBTClusterMap.h"


ClassImp(AliHBTReaderTPC)

AliHBTReaderTPC::AliHBTReaderTPC():
 fFileName("galice.root"),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fNClustMin(0),
 fNClustMax(150),
 fNChi2PerClustMin(0.0),
 fNChi2PerClustMax(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5)
{
  //constructor
  
}
/********************************************************************/

AliHBTReaderTPC::AliHBTReaderTPC(const Char_t* galicefilename):
 fFileName(galicefilename),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE),
 fNTrackPoints(0),
 fdR(0.0),
 fNClustMin(0),
 fNClustMax(150),
 fNChi2PerClustMin(0.0),
 fNChi2PerClustMax(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5)
{
  //constructor, 
  //Defaults:
}
/********************************************************************/

AliHBTReaderTPC::AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename):
 AliHBTReader(dirs), 
 fFileName(galicefilename),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(0.0),
 fUseMagFFromRun(kTRUE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fNClustMin(0),
 fNClustMax(150),
 fNChi2PerClustMin(0.0),
 fNChi2PerClustMax(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5)
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untached
  
}
/********************************************************************/
AliHBTReaderTPC::AliHBTReaderTPC(const AliHBTReaderTPC& in):
 AliHBTReader(in),
 fFileName(in.fFileName),
 fRunLoader(0x0),
 fTPCLoader(0x0),
 fMagneticField(in.fMagneticField),
 fUseMagFFromRun(in.fUseMagFFromRun),
 fNTrackPoints(in.fNTrackPoints),
 fdR(in.fdR),
 fNClustMin(in.fNClustMin),
 fNClustMax(in.fNClustMax),
 fNChi2PerClustMin(in.fNChi2PerClustMin),
 fNChi2PerClustMax(in.fNChi2PerClustMax),
 fC00Min(in.fC00Min),
 fC00Max(in.fC00Max),
 fC11Min(in.fC11Min),
 fC11Max(in.fC11Max),
 fC22Min(in.fC22Min),
 fC22Max(in.fC22Max),
 fC33Min(in.fC33Min),
 fC33Max(in.fC33Max),
 fC44Min(in.fC44Min),
 fC44Max(in.fC44Max)
{
  //cpy constructor, 
}
/********************************************************************/

AliHBTReaderTPC::~AliHBTReaderTPC()
{
 //desctructor
   if (AliHBTParticle::GetDebug()) 
    {
      Info("~AliHBTReaderTPC","deleting Run Loader");
      AliLoader::SetDebug(AliHBTParticle::GetDebug());
    }
   
   delete fRunLoader;
   
   if (AliHBTParticle::GetDebug()) 
    {
      Info("~AliHBTReaderTPC","deleting Run Loader Done");
    }
}
/********************************************************************/

AliHBTReaderTPC& AliHBTReaderTPC::operator=(const AliHBTReaderTPC& in)
{
//Assigment operator

 delete fRunLoader;

 fFileName = in.fFileName;
 fRunLoader = 0x0;
 fTPCLoader = 0x0;
 fMagneticField = in.fMagneticField;
 fUseMagFFromRun = in.fUseMagFFromRun;
 fNTrackPoints = in.fNTrackPoints;
 fdR = in.fdR;
 fNClustMin = in.fNClustMin;
 fNClustMax = in.fNClustMax;
 fNChi2PerClustMin = in.fNChi2PerClustMin;
 fNChi2PerClustMax = in.fNChi2PerClustMax;
 fC00Min = in.fC00Min;
 fC00Max = in.fC00Max;
 fC11Min = in.fC11Min;
 fC11Max = in.fC11Max;
 fC22Min = in.fC22Min;
 fC22Max = in.fC22Max;
 fC33Min = in.fC33Min;
 fC33Max = in.fC33Max;
 fC44Min = in.fC44Min;
 fC44Max = in.fC44Max;
 return *this; 
} 
/********************************************************************/

void AliHBTReaderTPC::Rewind()
{
//Rewind reading to the beginning
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;
  if (fTrackCounter) fTrackCounter->Reset();
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
         fCurrentEvent++;//go to next dir
         continue;
      }
   
    TBranch *trackbranch=tracktree->GetBranch("tracks");//get the branch with tracks
    if (!trackbranch) ////check if we got the branch
      {//if not return with error
        Error("ReadNext","Can't get a branch with TPC tracks !\n"); 
        fCurrentEvent++;//go to next dir
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
       
       if (CheckTrack(iotrack)) continue;//Checks the cuts on track parameters cov. mtx etc
       
       TParticle *p = (TParticle*)stack->Particle(label);

       if(p == 0x0) continue; //if returned pointer is NULL
       if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)

       if(Rejected(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                   //if not take next partilce

       AliHBTParticle* part = new AliHBTParticle(*p,i);
       if(Rejected(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
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
       if(Rejected(track))//check if meets all criteria of any of our cuts
                    //if it does not delete it and take next good track
        { 
          delete track;
          delete part;
          continue;
        }
        
       if (fNTrackPoints > 0) 
        {
          AliHBTTrackPoints* tpts = new AliHBTTrackPoints(fNTrackPoints,iotrack,fdR);
          track->SetTrackPoints(tpts);
        }
       if (  fClusterMap ) 
        {
          AliHBTClusterMap* cmap = new AliHBTClusterMap(iotrack);
          track->SetClusterMap(cmap);
        }
    
       fParticlesEvent->AddParticle(part);
       fTracksEvent->AddParticle(track);
      }
      
    Info("ReadNext","Read %d tracks and %d particles from event %d (event %d in dir %d).",
            fParticlesEvent->GetNumberOfParticles(), fTracksEvent->GetNumberOfParticles(),
            fNEventsRead,fCurrentEvent,fCurrentDir);
    fTrackCounter->Fill(fTracksEvent->GetNumberOfParticles());
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
//Opens session thats from fCurrentDir 
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
Bool_t AliHBTReaderTPC::CheckTrack(AliTPCtrack* t) const
{
  //Performs check of the track
  if ( (t->GetNumberOfClusters() > fNClustMax) || (t->GetNumberOfClusters() < fNClustMin) ) return kTRUE;

  Float_t chisqpercl = t->GetChi2()/((Double_t)t->GetNumberOfClusters());
  if ( (chisqpercl < fNChi2PerClustMin) || (chisqpercl > fNChi2PerClustMax) ) return kTRUE;
  
  Double_t cc[15];
  t->GetExternalCovariance(cc);

  if ( (cc[0]  < fC00Min) || (cc[0]  > fC00Max) ) return kTRUE;
  if ( (cc[2]  < fC11Min) || (cc[2]  > fC11Max) ) return kTRUE;
  if ( (cc[5]  < fC22Min) || (cc[5]  > fC22Max) ) return kTRUE;
  if ( (cc[9]  < fC33Min) || (cc[9]  > fC33Max) ) return kTRUE;
  if ( (cc[14] < fC44Min) || (cc[14] > fC44Max) ) return kTRUE;
  
  
  return kFALSE;
  
}
/********************************************************************/

void AliHBTReaderTPC::SetNClustersRange(Int_t min,Int_t max)
{
 //sets range of Number Of Clusters that tracks have to have
 fNClustMin = min;
 fNClustMax = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetChi2PerCluserRange(Float_t min, Float_t max)
{
  //sets range of Chi2 per Cluster
  fNChi2PerClustMin = min;
  fNChi2PerClustMax = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetC00Range(Float_t min, Float_t max)
{
 //Sets range of C00 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC00Min = min;
  fC00Max = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetC11Range(Float_t min, Float_t max)
{
 //Sets range of C11 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC11Min = min;
  fC11Max = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetC22Range(Float_t min, Float_t max)
{
 //Sets range of C22 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC22Min = min;
  fC22Max = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetC33Range(Float_t min, Float_t max)
{
 //Sets range of C33 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC33Min = min;
  fC33Max = max;
}
/********************************************************************/

void AliHBTReaderTPC::SetC44Range(Float_t min, Float_t max)
{
 //Sets range of C44 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC44Min = min;
  fC44Max = max;
}

/********************************************************************/
/********************************************************************/
/********************************************************************/

/*
void (AliTPCtrack* iotrack, Double_t curv)
{
  Double_t x[5];
  iotrack->GetExternalPara
  //xk is a 
  Double_t fP4=iotrack->GetC();
  Double_t fP2=iotrack->GetEta();
  
  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1, y1=fP0, z1=fP1;
  Double_t c1=fP4*x1 - fP2, r1=sqrt(1.- c1*c1);
  Double_t c2=fP4*x2 - fP2, r2=sqrt(1.- c2*c2);
  
  fP0 += dx*(c1+c2)/(r1+r2);
  fP1 += dx*(c1+c2)/(c1*r2 + c2*r1)*fP3;

}
*/

