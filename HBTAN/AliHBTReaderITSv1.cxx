//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
//  class AliHBTReaderITSv1                                         //
//                                                                  //
//  Reader for ITSv1 tracks. Not maintained since v1 is not         //
//  supposed to be used                                             //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TBranch.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>

#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"
#include "AliHBTReaderITSv1.h"
#include "AliHBTRun.h"
#include "AliITSIOTrack.h"
#include "AliKalmanTrack.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliStack.h"

ClassImp(AliHBTReaderITSv1)
/********************************************************************/

AliHBTReaderITSv1::AliHBTReaderITSv1(const Char_t* tracksfilename,const Char_t* galicefilename):
 fITSTracksFileName(tracksfilename),
 fGAliceFileName(galicefilename)
{
 //ctor
}
/********************************************************************/

AliHBTReaderITSv1::AliHBTReaderITSv1(TObjArray* dirs, const Char_t* tracksfilename,const Char_t* galicefilename):
 AliHBTReader(dirs),
 fITSTracksFileName(tracksfilename),
 fGAliceFileName(galicefilename)
{
//ctor
}
/********************************************************************/

AliHBTReaderITSv1::~AliHBTReaderITSv1()
{
//dtor
}
/********************************************************************/


Int_t AliHBTReaderITSv1::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
//Reads data
 Int_t nevents = 0;
 AliITSIOTrack *iotrack=new AliITSIOTrack;
 Int_t currentdir = 0;
 Int_t ndirs;
 Int_t totalnevents = 0;
 
 if (fDirs)
  {
    ndirs = fDirs->GetEntries();
  }
 else
  {
    ndirs = 0;
  }
 
 do //do while is good even if 
  {  
   TFile* gAliceFile = OpenGAliceFile(currentdir);
   if(gAliceFile == 0x0)
    {
       Error("Read","Can not open the file with gAlice");
       delete iotrack;
       return 1;
    }
   if (gAlice->TreeE())//check if tree E exists
     {
      nevents = (Int_t)gAlice->TreeE()->GetEntries();//if yes get number of events in gAlice
      cout<<"________________________________________________________\n";
      cout<<"Found "<<nevents<<" event(s) in directory "<<GetDirName(currentdir)<<endl;
      cout<<"Setting Magnetic Field. Factor is "<<gAlice->Field()->Factor()<<endl;
      AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());
     }
    else
     {//if not return an error
       Error("Read","Can not find Header tree (TreeE) in gAlice");
       delete iotrack;
       return 4;
     }

   TFile *file = OpenTrackFile(currentdir);
   if(file == 0x0)
    {
       Error("Read","Can not open the file with ITS tracks V1");
       delete iotrack;
       return 2;
    }
    
   Int_t naccepted = 0;
   char tname[30];
   
   for (Int_t currentEvent = 0; currentEvent < nevents; currentEvent++)
    { 
      cout<<"Reading Event "<<currentEvent;
      
      sprintf(tname,"TreeT%d",currentEvent);
      file->cd(); 
      TTree *tracktree=(TTree*)file->Get(tname);
      TBranch *tbranch=tracktree->GetBranch("ITStracks");
      tbranch->SetAddress(&iotrack);
      
      gAliceFile->cd();
      gAlice->GetEvent(currentEvent);
      gAlice->Stack()->Particles();

      Int_t nentr=(Int_t)tracktree->GetEntries();
      
      cout<<".  Found "<<nentr<<" tracks.";
      fflush(0);
      
      for (Int_t i=0; i<nentr; i++) 
       {

        tracktree->GetEvent(i);
        if(!iotrack) continue;       
        Int_t label = iotrack->GetLabel();
        if (label < 0) 
         {
           continue;
         }

        TParticle *p = (TParticle*)gAlice->Stack()->Particle(label);
        if(!p)
         {
           Warning("Read","Can not get particle with label &d",label);
           continue;
         }
        if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type
                                           //if not take next partilce

        AliHBTParticle* part = new AliHBTParticle(*p,i);
        if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                //if it does not delete it and take next good track
        
        Double_t px=iotrack->GetPx();
        Double_t py=iotrack->GetPy();
        Double_t pz=iotrack->GetPz();
        Double_t mass = p->GetMass();
        Double_t tEtot = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);//total energy of the track
 
        Double_t x= iotrack->GetX();
        Double_t y= iotrack->GetY();
        Double_t z= iotrack->GetZ();
        
        AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(), i, px, py , pz, tEtot, x, y, z, 0.);
        if(Pass(track)) { delete  track;continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track

        particles->AddParticle(totalnevents,part);//put track and particle on the run
        tracks->AddParticle(totalnevents,track);
        naccepted++;
       }//end loop over tracks in the event

       totalnevents++;
       cout<<"  Accepted "<<naccepted<<" tracks"<<endl;
     }//end of loop over events in current directory
    
    gAliceFile->Close();
    delete gAliceFile;
    gAliceFile = 0;
    
    file->Close(); 
    delete file;
    file = 0;
    currentdir++;
   }while(currentdir < ndirs);//end of loop over directories specified in fDirs Obj Array


  delete iotrack;
  fIsRead = kTRUE;
  return 0;
 
 }
/********************************************************************/

TFile* AliHBTReaderITSv1::OpenTrackFile(Int_t ndir)
{
//opens files to be read for given directoru nomber in fDirs Array
   const TString& dirname = GetDirName(ndir); 
   if (dirname == "")
    {
      Error("OpenGAliceFile","Can not get directory name");
      return 0x0;
    }
   TString filename = dirname + "/" + fITSTracksFileName;

   TFile *file = TFile::Open(filename.Data());   
   if (!file)
    {
      Error("Read","Can not open file %s",filename.Data());
      return 0x0;
    }
   if (!file->IsOpen())
    {
      Error("Read","Can not open file %s",filename.Data());
      return 0x0;
    }
   
   return file;
}


/********************************************************************/
TFile* AliHBTReaderITSv1::OpenGAliceFile(Int_t ndir)
{
//Opens galice.root file
  const TString& dirname = GetDirName(ndir); 
   if (dirname == "")
    {
      Error("OpenGAliceFile","Can not get directory name");
      return 0x0;
    }
  
  TString filename = dirname + "/" + fGAliceFileName;

  TFile* gAliceFile = TFile::Open(filename.Data());
  if ( gAliceFile== 0x0)
   {
     Error("OpenFiles","Can't open file named %s",filename.Data());
     return 0x0;
   }
  if (!gAliceFile->IsOpen())
   {
     Error("OpenFiles","Can't open file named %s",filename.Data());
     return 0x0;
   }

  if (!(gAlice=(AliRun*)gAliceFile->Get("gAlice")))
   {
     Error("OpenFiles","gAlice have not been found on %s !\n",filename.Data());
     gAliceFile->Close();
     delete gAliceFile;
     return 0x0;
   }
  
  return gAliceFile;
}

/********************************************************************/
/********************************************************************/


