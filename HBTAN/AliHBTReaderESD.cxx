#include "AliHBTReaderESD.h"
//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
// class AliHBTReaderESD                                            //
//                                                                  //
// reader for ALICE Event Summary Data (ESD).                       //
//                                                                  //
// Piotr.Skowronski@cern.ch                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TPDGCode.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TKey.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliESDtrack.h>
#include <AliESD.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"
#include "AliHBTTrackPoints.h"
#include "AliHBTClusterMap.h"

ClassImp(AliHBTReaderESD)

AliHBTReaderESD::AliHBTReaderESD(const Char_t* esdfilename, const Char_t* galfilename):
 fESDFileName(esdfilename),
 fGAlFileName(galfilename),
 fFile(0x0),
 fRunLoader(0x0),
 fKeyIterator(0x0),
 fReadParticles(kFALSE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fNTPCClustMin(0),
 fNTPCClustMax(150),
 fTPCChi2PerClustMin(0.0),
 fTPCChi2PerClustMax(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5),
 fTPCC00Min(0.0),
 fTPCC00Max(10e5),
 fTPCC11Min(0.0),
 fTPCC11Max(10e5),
 fTPCC22Min(0.0),
 fTPCC22Max(10e5),
 fTPCC33Min(0.0),
 fTPCC33Max(10e5),
 fTPCC44Min(0.0),
 fTPCC44Max(10e5)

{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliESDtrack::kSPECIES))
    Fatal("AliHBTReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/
  
AliHBTReaderESD::AliHBTReaderESD(TObjArray* dirs,const Char_t* esdfilename, const Char_t* galfilename):
 AliHBTReader(dirs), 
 fESDFileName(esdfilename),
 fGAlFileName(galfilename),
 fFile(0x0),
 fRunLoader(0x0),
 fKeyIterator(0x0),
 fReadParticles(kFALSE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fNTPCClustMin(0),
 fNTPCClustMax(150),
 fTPCChi2PerClustMin(0.0),
 fTPCChi2PerClustMax(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5),
 fTPCC00Min(0.0),
 fTPCC00Max(10e5),
 fTPCC11Min(0.0),
 fTPCC11Max(10e5),
 fTPCC22Min(0.0),
 fTPCC22Max(10e5),
 fTPCC33Min(0.0),
 fTPCC33Max(10e5),
 fTPCC44Min(0.0),
 fTPCC44Max(10e5)
{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliESDtrack::kSPECIES))
    Fatal("AliHBTReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/

AliHBTReaderESD::~AliHBTReaderESD()
{
 //desctructor
  delete fRunLoader;
  delete fKeyIterator;
  delete fFile;
}
/**********************************************************/
Int_t AliHBTReaderESD::ReadNext()
{
//reads next event from fFile
  //fRunLoader is for reading Kine
  //****** Tentative particle type "concentrations"
  static const Double_t concentr[5]={0.05, 0., 0.85, 0.10, 0.05};
  
  Double_t pidtable[kNSpecies];//array used for reading pid probabilities from ESD track
  Double_t w[kNSpecies];
  Double_t mom[3];//momentum
  Double_t pos[3];//position
  
  if (AliHBTParticle::GetDebug())
    Info("ReadNext","Entered");
    
  TDatabasePDG* pdgdb = TDatabasePDG::Instance();
  if (pdgdb == 0x0)
   {
     Error("Next","Can not get PDG Database Instance.");
     return 1;
   }
  
  if (fParticlesEvent == 0x0)  fParticlesEvent = new AliHBTEvent();
  if (fTracksEvent == 0x0)  fTracksEvent = new AliHBTEvent();
  
  fParticlesEvent->Reset();
  fTracksEvent->Reset();
        
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
     if (fFile == 0x0)
      {
       fFile = OpenFile(fCurrentDir);//rl is opened here
       if (fFile == 0x0)
         {
           Error("ReadNext","Cannot get fFile for dir no. %d",fCurrentDir);
           fCurrentDir++;
           continue;
         }
       fCurrentEvent = 0;
       fKeyIterator = new TIter(fFile->GetListOfKeys());  
       fFile->Dump();
       fFile->GetListOfKeys()->Print();
      } 
     TKey* key = (TKey*)fKeyIterator->Next();
     if (key == 0x0)
      {
        if (AliHBTParticle::GetDebug() > 2 )
          {
            Info("ReadNext","No more keys.");
          }
        fCurrentDir++;
        delete fKeyIterator;
        fKeyIterator = 0x0;
        delete fFile;//we have to assume there is no more ESD objects in the fFile
        fFile = 0x0;
        delete fRunLoader;
        fRunLoader = 0x0;
        continue;
      }
     //try to read
     
     
//     TObject* esdobj = key->ReadObj();
//     if (esdobj == 0x0)
//      {
//        if (AliHBTParticle::GetDebug() > 2 )
//          {
//            Info("ReadNext","Key read NULL. Key Name is %s",key->GetName());
//            key->Dump();
//          }
//        continue;
//      }
//     esdobj->Dump();
//     AliESD* esd = dynamic_cast<AliESD*>(esdobj);
     
     TString esdname = "ESD";
     esdname+=fCurrentEvent;
     AliESD* esd = dynamic_cast<AliESD*>(fFile->Get(esdname));
     if (esd == 0x0)
      {
//        if (AliHBTParticle::GetDebug() > 2 )
//          {
//            Info("ReadNext","This key is not an AliESD object %s",key->GetName());
//          }
        if (AliHBTParticle::GetDebug() > 2 )
          {
            Info("ReadNext","Can not find AliESD object named %s",esdname.Data());
          }
        fCurrentDir++;
        delete fKeyIterator;
        fKeyIterator = 0x0;
        delete fFile;//we have to assume there is no more ESD objects in the fFile
        fFile = 0x0;
        delete fRunLoader;
        fRunLoader = 0x0;
        continue;
      }


     Float_t mf = esd->GetMagneticField(); 
     
     if ( (mf == 0.0) && (fNTrackPoints > 0) )
      {
         Error("ReadNext","Magnetic Field is 0 and Track Points Demended. Skipping to next event.");
         fCurrentEvent++;
         continue;
      }
     
     AliStack* stack = 0x0;
     if (fReadParticles && fRunLoader)
      {
        fRunLoader->GetEvent(fCurrentEvent);
        stack = fRunLoader->Stack();
      }
     Info("ReadNext","Reading Event %d",fCurrentEvent);
  
     Int_t ntr = esd->GetNumberOfTracks();
     Info("ReadNext","Found %d tracks.",ntr);
     for (Int_t i = 0;i<ntr; i++)
      {
        AliESDtrack *esdtrack = esd->GetTrack(i);
        if (esdtrack == 0x0)
         {
           Error("Next","Can not get track %d", i);
           continue;
         }
        
        //if (esdtrack->HasVertexParameters() == kFALSE) 
        if ((esdtrack->GetStatus() & AliESDtrack::kITSrefit) == kFALSE)
         {
           if (AliHBTParticle::GetDebug() > 2) 
             Info("ReadNext","Particle skipped: Data at vertex not available.");
           continue;
         }
         
        if ((esdtrack->GetStatus()&AliESDtrack::kESDpid) == kFALSE) 
         {
           if (AliHBTParticle::GetDebug() > 2) 
             Info("ReadNext","Particle skipped: PID BIT is not set.");
           continue;
         }
       
        esdtrack->GetESDpid(pidtable);
        //esdtrack->GetVertexPxPyPz(mom[0],mom[1],mom[2]); 
        esdtrack->GetPxPyPz(mom);
        //esdtrack->GetVertexXYZ(pos[0],pos[1],pos[2]);
        esdtrack->GetXYZ(pos);
        
        Double_t extx;
        Double_t extp[5];
        esdtrack->GetExternalParameters(extx,extp);
        
        Int_t charge = (extp[4] > 0)?-1:1;//if curvature=-charg/Pt is positive charge is negative
        
        //Particle from kinematics
        AliHBTParticle* particle = 0;
        Bool_t keeppart = kFALSE;
        if ( fReadParticles && stack  )
         {
           if (esdtrack->GetLabel() < 0) continue;//this is fake -  we are not able to match any track 
           TParticle *p = stack->Particle(esdtrack->GetLabel());
           if (p==0x0) 
            {
              Error("ReadNext","Can not find track with such label.");
              continue;
            }
//           if(p->GetPdgCode()<0) charge = -1;
           particle = new AliHBTParticle(*p,i);
           
         }
         
        //Here we apply Bayes' formula
        Double_t rc=0.;
        for (Int_t s=0; s<AliESDtrack::kSPECIES; s++) rc+=concentr[s]*pidtable[s];
        if (rc==0.0) 
         {
           if (AliHBTParticle::GetDebug() > 2) 
             Info("ReadNext","Particle rejected since total bayessian PID probab. is zero.");
           continue;
         }

        for (Int_t s=0; s<AliESDtrack::kSPECIES; s++) w[s]=concentr[s]*pidtable[s]/rc;

        if (AliHBTParticle::GetDebug() > 4)
         { 
           Info("ReadNext","###########################################################################");
           Info("ReadNext","Momentum: %f %f %f",mom[0],mom[1],mom[2]);
           Info("ReadNext","Position: %f %f %f",pos[0],pos[1],pos[2]);
           TString msg("Pid list got from track:");
           for (Int_t s = 0;s<kNSpecies;s++)
            {
              msg+="\n    ";
              msg+=s;
              msg+="(";
              msg+=charge*GetSpeciesPdgCode((ESpecies)s);
              msg+="): ";
              msg+=w[s];
              msg+=" (";
              msg+=pidtable[s];
              msg+=")";
            }
           Info("ReadNext","%s",msg.Data());
         }//if (AliHBTParticle::GetDebug()>4)
         
         AliHBTTrackPoints* tpts = 0x0;
         if (fNTrackPoints > 0) 
          {
            tpts = new AliHBTTrackPoints(fNTrackPoints,esdtrack,mf,fdR);
          }

         AliHBTClusterMap* cmap = 0x0; 
         if (  fClusterMap ) 
          {
            cmap = new AliHBTClusterMap(esdtrack);
          }

        for (Int_t s = 0; s<kNSpecies; s++)
         {
           Float_t pp = w[s];
           if (pp == 0.0) continue;

           Int_t pdgcode = charge*GetSpeciesPdgCode((ESpecies)s);
           if(Pass(pdgcode)) continue; //check if we are intersted with particles of this type 

           Double_t mass = pdgdb->GetParticle(pdgcode)->Mass();
           Double_t tEtot = TMath::Sqrt( mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] + mass*mass);//total energy of the track

           AliHBTParticle* track = new AliHBTParticle(pdgcode, w[s],i, 
                                                      mom[0], mom[1], mom[2], tEtot,
                                                      pos[0], pos[1], pos[2], 0.);
           //copy probabilitis of other species (if not zero)
           for (Int_t k = 0; k<kNSpecies; k++)
            {
              if (k == s) continue;
              if (w[k] == 0.0) continue;
              track->SetPIDprobability(charge*GetSpeciesPdgCode( (ESpecies)k ),w[k]);
            }
            
           if(Pass(track))//check if meets all criteria of any of our cuts
                          //if it does not delete it and take next good track
            { 
              delete track;
              continue;
            }
            
           //Single Particle cuts on cluster map and track points rather do not have sense
           if (tpts)
            {
              track->SetTrackPoints(tpts); 
            }
            
           if (cmap) 
            { 
              track->SetClusterMap(cmap);
            }
           
           fTracksEvent->AddParticle(track);
           if (particle) fParticlesEvent->AddParticle(particle);
           keeppart = kTRUE;
               
           if (AliHBTParticle::GetDebug() > 4 )
            {
              Info("ReadNext","\n\nAdding Particle with incarnation %d",pdgcode);
              track->Print();
            }
         }//for (Int_t s = 0; s<kNSpecies; s++)
         
        if (keeppart == kFALSE) 
         {
           delete particle;//particle was not stored in event
           delete tpts;
           delete cmap;
         }
        
      }//for (Int_t i = 0;i<ntr; i++)  -- loop over tracks
     
     Info("ReadNext","Read %d tracks and %d particles from event %d (event %d in dir %d).",
            fTracksEvent->GetNumberOfParticles(), fParticlesEvent->GetNumberOfParticles(),
            fNEventsRead,fCurrentEvent,fCurrentDir);
      
     fCurrentEvent++;
     fNEventsRead++;
     delete esd;
     return 0;//success -> read one event
   }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array  
   
  return 1; //no more directories to read
}
/**********************************************************/
void AliHBTReaderESD::Rewind()
{
  //rewinds reading 
  delete fKeyIterator;
  delete fFile;
  fFile = 0x0;
  fKeyIterator = 0x0;
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead = 0;
  fCurrentEvent++;
}
/**********************************************************/

TFile* AliHBTReaderESD::OpenFile(Int_t n)
{
//opens fFile with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename = dirname +"/"+ fESDFileName;
 TFile *ret = TFile::Open(filename.Data()); 

 if ( ret == 0x0)
  {
    Error("OpenFiles","Can't open fFile %s",filename.Data());
    return 0x0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open fFile  %s",filename.Data());
    return 0x0;
  }
 
 if (fReadParticles )
  {
   fRunLoader = AliRunLoader::Open(dirname +"/"+ fGAlFileName);
   if (fRunLoader == 0x0)
    {
      Error("OpenFiles","Can't get RunLoader for directory %s",dirname.Data());
      delete ret;
      return 0x0;
    }
    
   fRunLoader->LoadHeader();
   if (fRunLoader->LoadKinematics())
    {
      Error("Next","Error occured while loading kinematics.");
      delete fRunLoader;
      delete ret;
      return 0x0;
    }
  }
   
 return ret;
}
/**********************************************************/

Int_t AliHBTReaderESD::GetSpeciesPdgCode(ESpecies spec)//skowron
{
  //returns pdg code from the PID index
  //ask jura about charge
  switch (spec)
   {
     case kESDElectron:
       return kPositron;
       break;
     case kESDMuon:
       return kMuonPlus;
       break;
     case kESDPion:
       return kPiPlus;
       break;
     case kESDKaon:
       return kKPlus;
       break;
     case kESDProton:
       return kProton;
       break;
     default:
       ::Warning("GetSpeciesPdgCode","Specie with number %d is not defined.",(Int_t)spec);
       break;
   }
  return 0;
}
/********************************************************************/

void AliHBTReaderESD::SetTPCNClustersRange(Int_t min,Int_t max)
{
 //sets range of Number Of Clusters that tracks have to have
 fNTPCClustMin = min;
 fNTPCClustMax = max;
}
/********************************************************************/

void AliHBTReaderESD::SetTPCChi2PerCluserRange(Float_t min, Float_t max)
{
  //sets range of Chi2 per Cluster
  fTPCChi2PerClustMin = min;
  fTPCChi2PerClustMax = max;
}
/********************************************************************/

void AliHBTReaderESD::SetC00Range(Float_t min, Float_t max)
{
 //Sets range of C00 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC00Min = min;
  fC00Max = max;
}
/********************************************************************/

void AliHBTReaderESD::SetC11Range(Float_t min, Float_t max)
{
 //Sets range of C11 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC11Min = min;
  fC11Max = max;
}
/********************************************************************/

void AliHBTReaderESD::SetC22Range(Float_t min, Float_t max)
{
 //Sets range of C22 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC22Min = min;
  fC22Max = max;
}
/********************************************************************/

void AliHBTReaderESD::SetC33Range(Float_t min, Float_t max)
{
 //Sets range of C33 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC33Min = min;
  fC33Max = max;
}
/********************************************************************/

void AliHBTReaderESD::SetC44Range(Float_t min, Float_t max)
{
 //Sets range of C44 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC44Min = min;
  fC44Max = max;
}
