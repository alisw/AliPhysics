#include "AliHBTReaderInternal.h"
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//                                                                       //
// class AliHBTReaderInternal                                            //
//                                                                       //
// Multi file reader for Internal Data Format                            //
//                                                                       //
// This reader reads data created by itself                              //
//   (method AliHBTReaderInternal::Write)                                //
// Data are stored in form of tree of TClonesArray of AliHBTParticle's   //
//                                                                       //
// Piotr.Skowronski@cern.ch                                              //
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TError.h>
#include <AliRun.h>
#include <AliMagF.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

ClassImp(AliHBTReaderInternal)
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal():
 fFileName(),
 fPartBranch(0x0),
 fTrackBranch(0x0),
 fTree(0x0),
 fFile(0x0),
 fPartBuffer(0x0),
 fTrackBuffer(0x0)
{
//Defalut constructor
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(const char *filename):
 fFileName(filename),
 fPartBranch(0x0),
 fTrackBranch(0x0),
 fTree(0x0),
 fFile(0x0),
 fPartBuffer(0x0),
 fTrackBuffer(0x0)
{ 
//constructor 
//filename - name of file to open
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(TObjArray* dirs, const char *filename):
 AliHBTReader(dirs),
 fFileName(filename),
 fPartBranch(0x0),
 fTrackBranch(0x0),
 fTree(0x0),
 fFile(0x0),
 fPartBuffer(0x0),
 fTrackBuffer(0x0)
{ 
//ctor
//dirs contains strings with directories to look data in
//filename - name of file to open
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(const AliHBTReaderInternal& in):
 AliHBTReader(in),
 fFileName(in.fFileName),
 fPartBranch(0x0),
 fTrackBranch(0x0),
 fTree(0x0),
 fFile(0x0),
 fPartBuffer(0x0),
 fTrackBuffer(0x0)
{
  //cpy constructor
}
/********************************************************************/

AliHBTReaderInternal& AliHBTReaderInternal::operator=(const AliHBTReaderInternal& in)
{
  //Assigment operator
  if (this == &in) return *this;
  Rewind();//close current session
  AliHBTReader::operator=((const AliHBTReader&)in);
  fFileName = in.fFileName;
  return *this;
}
/********************************************************************/
AliHBTReaderInternal::~AliHBTReaderInternal()
 {
 //desctructor
   delete fTree;
   delete fFile;
 }
/********************************************************************/
 
void AliHBTReaderInternal::Rewind()
{
  delete fTree;
  fTree = 0x0;
  delete fFile;
  fFile = 0x0;
  fCurrentDir = 0;
  fNEventsRead= 0;
}
/********************************************************************/

Int_t AliHBTReaderInternal::ReadNext()
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("ReadNext","");
  Int_t i; //iterator and some temprary values
  Int_t counter;
  AliHBTParticle* tpart = 0x0, *ttrack = 0x0;
  
  TDatabasePDG* pdgdb = TDatabasePDG::Instance();  
  if (pdgdb == 0x0)
   {
     Error("ReadNext","Can not get PDG Particles Data Base");
     return 1;
   }
   
  if (fParticlesEvent == 0x0)  fParticlesEvent = new AliHBTEvent();
  if (fTracksEvent == 0x0)  fTracksEvent = new AliHBTEvent();
  
  fParticlesEvent->Reset();
  fTracksEvent->Reset();
  
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
    if (fTree == 0x0)
     if( (i=OpenNextFile()) )
      {
        Error("ReadNext","Skipping directory due to problems with opening files. Errorcode %d",i);
        fCurrentDir++;
        continue;
      }
    if (fCurrentEvent == (Int_t)fTree->GetEntries())
     {
       delete fTree;
       fTree = 0x0;
       delete fFile;
       fFile = 0x0;
       fPartBranch = 0x0;
       fTrackBranch= 0x0;
       fCurrentDir++;
       continue;
     }
   /***************************/
   /***************************/
   /***************************/
    
        
    Info("ReadNext","Event %d",fCurrentEvent);
    fTree->GetEvent(fCurrentEvent);

    counter = 0;  
    if (fPartBranch && fTrackBranch)
     {
       Info("ReadNext","Found %d tracks in total.",fTrackBuffer->GetEntries());	
       Info("ReadNext","Found %d particles in total.",fPartBuffer->GetEntries());
       for(i = 0; i < fPartBuffer->GetEntries(); i++)
         {
           tpart = dynamic_cast<AliHBTParticle*>(fPartBuffer->At(i));
           ttrack =  dynamic_cast<AliHBTParticle*>(fTrackBuffer->At(i));

           if( tpart == 0x0 ) continue; //if returned pointer is NULL
           if( ttrack == 0x0 ) continue; //if returned pointer is NULL
           
           if (AliHBTParticle::GetDebug() > 9)
            {
              Info("ReadNext","Particle:");
              tpart->Print();
              Info("ReadNext","Track:");
              ttrack->Print();
            }
           if (ttrack->GetUID() != tpart->GetUID())
             {
               Error("ReadNext","Sth. is wrong: Track and Particle has different UID.");
               Error("ReadNext","They probobly do not correspond to each other.");
             }
           
           for (Int_t s = 0; s < ttrack->GetNumberOfPids(); s++)
            {
              //check if we are intersted with particles of this type
              //if not take next partilce
              if( Rejected(ttrack->GetNthPid(s)) ) 
               {
                 if (AliHBTParticle::GetDebug() > 9)
                  Info("ReadNext","Track Incarnation %d did not pass PID cut.",ttrack->GetNthPid(s));
                 continue; 
               }
              TParticlePDG* pdgp = pdgdb->GetParticle(ttrack->GetNthPid(s));
              if (pdgp == 0x0)//PDG part corresponding to new incarnation
               {
                 Error("ReadNext","Particle code unknown to PDG DB.");
                 continue;
               }
              
              AliHBTParticle* track = new AliHBTParticle(*ttrack);
              
              //apart of setting PDG code of an incarnation
              //it is necessary tu recalculate energy on the basis of
              //new PDG code (mass) hypothesis
              Double_t mass = pdgp->Mass();//mass of new incarnation
              Double_t tEtot = TMath::Sqrt( ttrack->Px()*ttrack->Px() + 
                                            ttrack->Py()*ttrack->Py() + 
                                            ttrack->Pz()*ttrack->Pz() + 
                                            mass*mass);//total energy of the new incarnation
              //update Energy and Calc Mass 
              track->SetMomentum(ttrack->Px(),ttrack->Py(),ttrack->Pz(),tEtot);
              track->SetCalcMass(mass);
              track->SetPdgCode(ttrack->GetNthPid(s),ttrack->GetNthPidProb(s));
              
              if( Rejected(track) )
                {
                  if (AliHBTParticle::GetDebug() > 9)
                   Info("ReadNext","Track Incarnation %d did not pass cut.",ttrack->GetNthPid(s));
                  delete track;
                  continue; 
                }
              AliHBTParticle* part = new AliHBTParticle(*tpart);//particle has only 1 incarnation (real)

              fParticlesEvent->AddParticle(part);
              fTracksEvent->AddParticle(track);

              counter++;
            }
         }
        Info("ReadNext","   Read: %d particles and tracks.",counter);
    }
   else
    {  
     if (fPartBranch)
      {
        Info("ReadNext","Found %d particles in total.",fPartBuffer->GetEntries());	
        for(i = 0; i < fPartBuffer->GetEntries(); i++)
         { 
           tpart = dynamic_cast<AliHBTParticle*>(fPartBuffer->At(i));
           if(tpart == 0x0) continue; //if returned pointer is NULL

           for (Int_t s = 0; s < tpart->GetNumberOfPids(); s++)
            {
              if( pdgdb->GetParticle(tpart->GetNthPid(s)) == 0x0 ) continue; //if particle has crazy PDG code (not known to our database)
              if( Rejected(tpart->GetNthPid(s)) ) continue; //check if we are intersted with particles of this type
              AliHBTParticle* part = new AliHBTParticle(*tpart);
              part->SetPdgCode(tpart->GetNthPid(s),tpart->GetNthPidProb(s));
              if( Rejected(part) )
                {
                  delete part;
                  continue; 
                }
              fParticlesEvent->AddParticle(part);
              counter++;
            }
         }
        Info("ReadNext","   Read: %d particles.",counter);
      }
     else Info("ReadNext","   Read: 0 particles.");

     if (fTrackBranch)
      {
        Info("ReadNext","Found %d tracks in total.",fTrackBuffer->GetEntries());	
        for(i = 0; i < fTrackBuffer->GetEntries(); i++)
         {
           ttrack =  dynamic_cast<AliHBTParticle*>(fTrackBuffer->At(i));
           if( ttrack == 0x0 ) continue; //if returned pointer is NULL

           for (Int_t s = 0; s < ttrack->GetNumberOfPids(); s++)
            {
              if( Rejected(ttrack->GetNthPid(s)) ) continue; //check if we are intersted with particles of this type
                                                         //if not take next partilce
              TParticlePDG* pdgp = pdgdb->GetParticle(ttrack->GetNthPid(s));
              if (pdgp == 0x0)//PDG part corresponding to new incarnation
               {
                 Error("ReadNext","Particle code unknown to PDG DB.");
                 continue;
               }
              AliHBTParticle* track = new AliHBTParticle(*ttrack);
              
              //apart of setting PDG code of an incarnation
              //it is necessary tu recalculate energy on the basis of
              //new PDG code (mass) hypothesis
              Double_t mass = pdgp->Mass();//mass of new incarnation
              Double_t tEtot = TMath::Sqrt( ttrack->Px()*ttrack->Px() + 
                                            ttrack->Py()*ttrack->Py() + 
                                            ttrack->Pz()*ttrack->Pz() + 
                                            mass*mass);//total energy of the new incarnation
              //update Energy and Calc Mass 
              track->SetMomentum(ttrack->Px(),ttrack->Py(),ttrack->Pz(),tEtot);
              track->SetCalcMass(mass);
              track->SetPdgCode(ttrack->GetNthPid(s),ttrack->GetNthPidProb(s));
              
              if( Rejected(track) )
                {
                  delete track;
                  continue; 
                }
              fTracksEvent->AddParticle(track);

              counter++;
            }
         }
        Info("ReadNext","   Read: %d tracks",counter);
      }
     else Info("ReadNext","   Read: 0 tracks.");
    }
    
    if (fPartBuffer) fPartBuffer->Delete();
    if (fTrackBuffer) fTrackBuffer->Delete();
    
    fCurrentEvent++;
    fNEventsRead++;
    
    return 0;

   }while(fCurrentDir < GetNumberOfDirs());

  return 1;//no more directories to read
}
/********************************************************************/

Int_t AliHBTReaderInternal::OpenNextFile()
{
  //open file in current directory
   const TString& dirname = GetDirName(fCurrentDir);
   if (dirname == "")
    {
      Error("OpenNextFile","Can not get directory name");
      return 4;
    }
   
   TString filename = dirname +"/"+ fFileName;
   fFile = TFile::Open(filename.Data());
   if ( fFile  == 0x0 ) 
     {
       Error("OpenNextFile","Can't open file named %s",filename.Data());
       return 1;
     }
   if (fFile->IsOpen() == kFALSE)
     {
       Error("OpenNextFile","Can't open filenamed %s",filename.Data());
       return 1;
     }
   
   fTree = (TTree*)fFile->Get("data");
   if (fTree == 0x0)
    {
     Error("OpenNextFile","Can not get the tree.");
     return 1;
    }

    
   fTrackBranch = fTree->GetBranch("tracks");//get the branch with tracks
   if (fTrackBranch == 0x0) ////check if we got the branch
     {//if not return with error
       Info("OpenNextFile","Can't find a branch with tracks !\n"); 
     }
   else
    {
      if (fTrackBuffer == 0x0)
        fTrackBuffer = new TClonesArray("AliHBTParticle",15000);
      fTrackBranch->SetAddress(&fTrackBuffer);
    }

   fPartBranch = fTree->GetBranch("particles");//get the branch with particles
   if (fPartBranch == 0x0) ////check if we got the branch
     {//if not return with error
       Info("OpenNextFile","Can't find a branch with particles !\n"); 
     }
   else
    {
      if (fPartBuffer == 0x0)
        fPartBuffer = new TClonesArray("AliHBTParticle",15000);
      fPartBranch->SetAddress(&fPartBuffer);
    }

   Info("OpenNextFile","________________________________________________________");
   Info("OpenNextFile","Found %d event(s) in directory %s",
         (Int_t)fTree->GetEntries(),GetDirName(fCurrentDir).Data());
   
   fCurrentEvent = 0;

   return 0; 
}
/********************************************************************/

Int_t AliHBTReaderInternal::Write(AliHBTReader* reader,const char* outfile, Bool_t multcheck)
 {
  //reads tracks from reader and writes runs to file
  //reader - provides data for writing in internal format
  //name of output file
  //multcheck - switches of checking if particle was stored with other incarnation
  // usefull e.g. when using kine data, where all particles have 100% pid prob.and saves a lot of time
  
  Int_t i,j;
  
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");

  TFile *histoOutput = TFile::Open(outfile,"recreate");
  
  if (!histoOutput->IsOpen())
   {
     ::Error("AliHBTReaderInternal::Write","File is not opened");
     return 1;
   }
    
  TTree *tracktree = new TTree("data","Tree with tracks");

  TClonesArray* pbuffer = new TClonesArray("AliHBTParticle",15000);
  TClonesArray* tbuffer = new TClonesArray("AliHBTParticle",15000);

  TClonesArray &particles = *pbuffer;
  TClonesArray &tracks = *tbuffer;
    
  TString name("Tracks");
  
  Int_t nt = reader->GetNumberOfTrackEvents();
  Int_t np = reader->GetNumberOfPartEvents();
  
  if (AliHBTParticle::GetDebug() > 0)
   ::Info("Write","Reader has %d track events and %d particles events.",nt,np);
   
  Bool_t trck = (nt > 0) ? kTRUE : kFALSE;
  Bool_t part = (np > 0) ? kTRUE : kFALSE;

  TBranch *trackbranch = 0x0, *partbranch = 0x0;
  
  if (trck) trackbranch = tracktree->Branch("tracks",&tbuffer,32000,0);
  if (part) partbranch = tracktree->Branch("particles",&pbuffer,32000,0);
  
  if ( (trck) && (part) && (np != nt))
   {
     ::Warning("AliHBTReaderInternal::Write","Number of track and particle events is different");
   }
  
  Int_t n;
  if (nt >= np ) n = nt; else n = np;
  
  if (AliHBTParticle::GetDebug() > 0)
   ::Info("Write","Will loop over %d events",n);

  for ( i =0;i< n; i++)
    {
      ::Info("AliHBTReaderInternal::Write","Event %d",i+1);
      Int_t counter = 0;
      if (trck && (i<=nt))
       { 
         AliHBTEvent* trackev = reader->GetTrackEvent(i);
         for ( j = 0; j< trackev->GetNumberOfParticles();j++)
          {
            const AliHBTParticle& t = *(trackev->GetParticle(j));
            if (multcheck)
             {
              if (FindIndex(tbuffer,t.GetUID())) 
               {
                 if (AliHBTParticle::GetDebug()>4)
                  { 
                   ::Info("Write","Track with Event UID %d already stored",t.GetUID());
                  }
                 continue; //not to write the same particles with other incarnations
               }
             }
            new (tracks[counter++]) AliHBTParticle(t);
          }
         ::Info("AliHBTReaderInternal::Write","    Tracks: %d",tracks.GetEntries());
       }else ::Info("AliHBTReaderInternal::Write","NO TRACKS");
      
      counter = 0;
      if (part && (i<=np))
       {
//        ::Warning("AliHBTReaderInternal::Write","Find index switched off!!!");

        AliHBTEvent* partev = reader->GetParticleEvent(i);
        for ( j = 0; j< partev->GetNumberOfParticles();j++)
         {
           const AliHBTParticle& part= *(partev->GetParticle(j));
            if (multcheck)
             {
              if (FindIndex(pbuffer,part.GetUID())) 
               {
                 if (AliHBTParticle::GetDebug()>4)
                  { 
                   ::Info("Write","Particle with Event UID %d already stored",part.GetUID());
                  }
                 continue; //not to write the same particles with other incarnations
               }
             } 
           new (particles[counter++]) AliHBTParticle(part);
         }
         ::Info("AliHBTReaderInternal::Write","    Particles: %d",particles.GetEntries());
       }else ::Info("AliHBTReaderInternal::Write","NO PARTICLES");

      histoOutput->cd();
      tracktree->Fill();
      tracktree->AutoSave();
      tbuffer->Delete();
      pbuffer->Delete();
     }

  histoOutput->cd();
  tracktree->Write(0,TObject::kOverwrite);
  delete tracktree;

  tbuffer->SetOwner();
  pbuffer->SetOwner();
  delete pbuffer;
  delete tbuffer;

  histoOutput->Close();
  return 0;
 }
/********************************************************************/

Bool_t AliHBTReaderInternal::FindIndex(TClonesArray* arr,Int_t idx)
{
//Checks if in the array exists already partilce with Unique ID idx
  if (arr == 0x0)
   {
     ::Error("FindIndex","Array is 0x0");
     return kTRUE;
   }
  TIter next(arr);
  AliHBTParticle* p;
  while (( p = (AliHBTParticle*)next()))
   {
     if (p->GetUID() == idx) return kTRUE;
   }
  return kFALSE;
}
