#include "AliHBTReaderESD.h"

#include <TPDGCode.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliESDtrack.h>
#include <AliESD.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

ClassImp(AliHBTReaderESD)

AliHBTReaderESD::AliHBTReaderESD(const Char_t* esdfilename):
 fParticles(new AliHBTRun()),
 fTracks(new AliHBTRun()),
 fESDFileName(esdfilename),
 fIsRead(kFALSE)
{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliESDtrack::kSPECIES))
    Fatal("AliHBTReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/
  
AliHBTReaderESD::AliHBTReaderESD(TObjArray* dirs,const Char_t* esdfilename):
 AliHBTReader(dirs), 
 fParticles(new AliHBTRun()),
 fTracks(new AliHBTRun()),
 fESDFileName(esdfilename),
 fIsRead(kFALSE)
{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliESDtrack::kSPECIES))
    Fatal("AliHBTReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/

AliHBTReaderESD::~AliHBTReaderESD()
{
 //desctructor
  delete fParticles;
  delete fTracks;
}
/********************************************************************/

AliHBTEvent* AliHBTReaderESD::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
   if (!fIsRead) 
    if(Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
   return fParticles->GetEvent(n);
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderESD::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
   if (!fIsRead) 
    if(Read(fParticles,fTracks))
     {
       Error("GetTrackEvent","Error in reading");
       return 0x0;
     }
   return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReaderESD::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) 
    if ( Read(fParticles,fTracks))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderESD::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead)
    if(Read(fParticles,fTracks))
     {
       Error("GetNumberOfTrackEvents","Error in reading");
       return 0;
     }
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/


Int_t AliHBTReaderESD::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");
  Int_t totalNevents = 0;
  Int_t currentdir = 0; //number of current directory name is fDirs array
  
  Double_t pidtable[kNSpecies];//array used for reading pid probabilities from ESD track
  Double_t w[kNSpecies];
  Double_t mom[3];//momentum
  Double_t pos[3];//position
   //****** Tentative particle type "concentrations"
  const Double_t concentr[5]={0.05, 0., 0.85, 0.10, 0.05};
  
  TDatabasePDG* pdgdb = TDatabasePDG::Instance();
  if (pdgdb == 0x0)
   {
     Error("Read","Can not get PDG Database Instance.");
     return 1;
   }
  if (!particles) //check if an object is instatiated
   {
     Error("Read"," particles object must instatiated before passing it to the reader");
     return 1;
   }
  if (!tracks)  //check if an object is instatiated
   {
     Error("Read"," tracks object must instatiated before passing it to the reader");
     return 1;
   }
  particles->Reset();//clear runs == delete all old events
  tracks->Reset();

  Int_t Ndirs;
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
     TFile* file = OpenFile(currentdir);
     if (file == 0x0)
       {
         Error("Read","Cannot get File for dir no. %d",currentdir);
         currentdir++;
         continue;
       } 
       
     for(Int_t currentEvent = 0; ;currentEvent++)//loop over all events
      {
       TString esdname;
       esdname+=currentEvent;
       AliESD* esd = dynamic_cast<AliESD*>(file->Get(esdname));
       if (esd == 0x0)
         {
           if (AliHBTParticle::fgDebug > 2 )
            {
              Info("Read","Can not find ESD object named %s.",esdname.Data());
            }
           currentdir++;
           break;//we have to assume there is no more ESD objects in the file
         } 
       
       Info("Read","Reading Event %d",currentEvent);
  
       Int_t ntr = esd->GetNumberOfTracks();
       Info("Read","Found %d tracks.",ntr);
       for (Int_t i = 0;i<ntr; i++)
        {
          AliESDtrack *esdtrack = esd->GetTrack(i);
          if (esdtrack == 0x0)
           {
             Error("Read","Can not get track %d", i);
             continue;
           }
          if ((esdtrack->GetStatus()&AliESDtrack::kESDpid) == kFALSE) 
           {
             if (AliHBTParticle::fgDebug > 2) 
               Info("Read","Particle skipped: PID BIT is not set.");
             continue;
           }

          esdtrack->GetESDpid(pidtable);
          esdtrack->GetPxPyPz(mom);
          esdtrack->GetXYZ(pos);
          Double_t rc=0.;

           //Here we apply Bayes' formula
          for (Int_t s=0; s<AliESDtrack::kSPECIES; s++) rc+=concentr[s]*pidtable[s];
          if (rc==0.0) 
           {
             if (AliHBTParticle::fgDebug > 2) 
               Info("Read","Particle rejected since total bayessian PID probab. is zero.");
             continue;
           }

          for (Int_t s=0; s<AliESDtrack::kSPECIES; s++) w[s]=concentr[s]*pidtable[s]/rc;

          if (AliHBTParticle::fgDebug > 4)
           { 
             Info("Read","###########################################################################");
             Info("Read","Momentum: %f %f %f",mom[0],mom[1],mom[2]);
             Info("Read","Position: %f %f %f",pos[0],pos[1],pos[2]);
             Info("Read","Radius: %f ,rc: %f",TMath::Hypot(pos[0],pos[1]),rc);
             TString msg("Pid list got from track:");
             for (Int_t s = 0;s<kNSpecies;s++)
              {
                msg+="\n    ";
                msg+=s;
                msg+="(";
                msg+=GetSpeciesPdgCode((ESpecies)s);
                msg+="): ";
                msg+=w[s];
                msg+=" (";
                msg+=pidtable[s];
                msg+=")";
              }
             Info("Read","%s",msg.Data());
           }//if (AliHBTParticle::fgDebug>4)
           
          for (Int_t s = 0; s<kNSpecies; s++)
           {
             if (w[s] == 0.0) continue;

             Int_t pdgcode = GetSpeciesPdgCode((ESpecies)s);
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
                track->SetPIDprobability(GetSpeciesPdgCode( (ESpecies)k ),w[k]);
              }
             if(Pass(track))//check if meets all criteria of any of our cuts
                            //if it does not delete it and take next good track
               { 
                 delete track;
                 continue;
               }
             tracks->AddParticle(totalNevents,track);
             if (AliHBTParticle::fgDebug > 4 )
              {
                Info("Read","\n\nAdding Particle with incarnation %d",pdgcode);
                track->Print();
              }
           }//for (Int_t s = 0; s<kNSpecies; s++)

        }//for (Int_t i = 0;i<ntr; i++)  -- loop over tracks
       totalNevents++;
      }//for(Int_t currentEvent = 0; ;currentEvent++) -- for over ESD in single file
      
   }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array
  
  fIsRead = kTRUE;
  return 0;
}      
/**********************************************************/
   
TFile* AliHBTReaderESD::OpenFile(Int_t n)
{
//opens file with kine tree

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
