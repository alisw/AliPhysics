// $Id$

//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
// class AliHBTReaderESD                                            //
//                                                                  //
// Reader for ALICE Event Summary Data (ESD).                       //
//                                                                  //
// Piotr.Skowronski@cern.ch                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <AliKalmanTrack.h>
#include <AliJetEventParticles.h>
#include "AliJetParticlesReaderESD.h"

ClassImp(AliJetParticlesReaderESD)

AliJetParticlesReaderESD::AliJetParticlesReaderESD(Bool_t constrained,
						   const Char_t* esdfilename) :
  AliJetParticlesReader(),
  fConstrained(constrained),
  fESDFileName(esdfilename),
  fESD(0),
  fFile(0),
  fTree(0),
  fKeyIterator(0),
  fPassFlag(AliESDtrack::kTPCrefit)
{
  //constructor
}

/********************************************************************/
  
AliJetParticlesReaderESD::AliJetParticlesReaderESD(
                                      Bool_t constrained,
                                      TObjArray* dirs,
                                      const Char_t* esdfilename) :
  AliJetParticlesReader(dirs),
  fConstrained(constrained),
  fESDFileName(esdfilename),
  fESD(0),
  fFile(0),
  fTree(0),
  fKeyIterator(0),
  fPassFlag(AliESDtrack::kTPCrefit)
{
  //constructor
}

/********************************************************************/

AliJetParticlesReaderESD::~AliJetParticlesReaderESD()
{
  //desctructor
  if(fESD) delete fESD;
  if(fTree) delete fTree;
  if(fKeyIterator) delete fKeyIterator;
  if(fFile) delete fFile;
}

/**********************************************************/

Int_t AliJetParticlesReaderESD::ReadESD(AliESD* esd)
{
  //Reads one ESD

  if (esd == 0)
   {
     Error("ReadESD","ESD is NULL");
     return kFALSE;
   }

  /*
  TDatabasePDG* pdgdb = TDatabasePDG::Instance();
  if (pdgdb == 0)
  {
     Error("ReadESD","Can not get PDG Database Instance.");
     return kFALSE;
  }
  */
   
  Float_t mf = esd->GetMagneticField(); 
  if (mf <= 0.0)  
  {
     Error("ReadESD","Magnetic Field is 0. Skipping to next event.");
     return kFALSE;
  }
  AliKalmanTrack::SetMagneticField(mf/10.);

  Info("ReadESD","Reading Event %d",fCurrentDir*1000+fCurrentEvent);
  if((!fOwner) || (fEventParticles==0)) 
    fEventParticles = new AliJetEventParticles();

  const Int_t kntr = esd->GetNumberOfTracks();
  Info("ReadESD","Found %d tracks.",kntr);
  fEventParticles->Reset(kntr);

  TString headdesc="";
  headdesc+="Run ";
  headdesc+=esd->GetRunNumber();
  headdesc+=": Ev ";
  headdesc+=esd->GetEventNumber();
  fEventParticles->SetHeader(headdesc);

  Double_t vertexpos[3];//vertex position
  const AliESDVertex* kvertex = esd->GetVertex();
  if (kvertex == 0)
   {
     Info("ReadESD","ESD returned NULL pointer to vertex - assuming (0.0,0.0,0.0)");
     vertexpos[0] = 0.0;
     vertexpos[1] = 0.0;
     vertexpos[2] = 0.0;
   }
  else
   {
     kvertex->GetXYZ(vertexpos);
   }
  fEventParticles->SetVertex(vertexpos[0],vertexpos[1],vertexpos[2]);

  //loop over tracks
  for (Int_t i = 0;i<kntr; i++)
   {

     const AliESDtrack *kesdtrack = esd->GetTrack(i);
     if (kesdtrack == 0)
      {
        Error("ReadESD","Can't get track %d", i);
        continue;
      }
     
     ULong_t cmptest=(kesdtrack->GetStatus() & fPassFlag);
     if (cmptest!=fPassFlag)
      {
	Info("ReadNext","Particle %d skipped: %u.",i,kesdtrack->GetStatus());
	//cout << i << " "; PrintESDtrack(kesdtrack); cout << endl;
        continue;
      }

     Double_t mom[3];  //momentum
     Double_t xyz[3];  //position
     if (fConstrained) {
       if (kesdtrack->GetConstrainedChi2() > 25) continue;
       kesdtrack->GetConstrainedPxPyPz(mom);
       kesdtrack->GetConstrainedXYZ(xyz);
     } else {
       if(!kesdtrack->GetPxPyPzAt(0,mom)) continue;
       kesdtrack->GetXYZAt(0, xyz);
     }
     const Float_t kmass=kesdtrack->GetMass();
     const Float_t kp2=mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2];
     const Float_t ketot=TMath::Sqrt(kmass*kmass+kp2);
     const Float_t kpt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
     const Float_t kp=TMath::Sqrt(kp2);
     const Float_t keta=0.5*TMath::Log((kp+mom[2]+1e-30)/(kp-mom[2]+1e-30)); 
     const Float_t kphi=TMath::Pi()+TMath::ATan2(-mom[1],-mom[0]);
     //Double_t dx = xyz[0]-vertexpos[0];
     //Double_t dy = xyz[1]-vertexpos[1];
     //Float_t dca = TMath::Sqrt(dx*dx + dy*dy);
     //Float_t dz = xyz[2]-vertexpos[2];
     UInt_t index[6];
     const Int_t kncl=kesdtrack->GetITSclusters(index)
                      +kesdtrack->GetTPCclusters(NULL)
                      +kesdtrack->GetTRDclusters(NULL);
     if(IsAcceptedParticle(kpt,kphi,keta))
       fEventParticles->AddParticle(mom[0],mom[1],mom[2],ketot,i,kesdtrack->GetLabel(),kncl,kpt,kphi,keta);

   } // loop over tracks

  return kTRUE;
}

/**********************************************************/

void AliJetParticlesReaderESD::Rewind()
{
  //rewinds reading 

  if(fFile) delete fFile;
  if(fKeyIterator) delete fKeyIterator;
  fKeyIterator = 0;
  fFile = 0;
  fCurrentDir = 0;
  fNEventsRead = 0;
}

/**********************************************************/

Int_t AliJetParticlesReaderESD::ReadNext()
{
  //reads next event from fFile

  do   // is OK even if 0 dirs specified, 
    {  // in that case we try to read from "./"
      if (fFile == 0)
	{
	  fFile = OpenFile(fCurrentDir);
	  if (fFile == 0)
	    {
	      Error("ReadNext","Can't get fFile for dir no. %d",fCurrentDir);
	      fCurrentDir++;
	      continue;
	    }
     
	  fCurrentEvent = 0;
	  //fFile->GetListOfKeys()->Print();
	  
	  if(fTree) delete fTree;
	  fTree = dynamic_cast<TTree*>(fFile->Get("esdTree"));
	  if(fTree)
	    fTree->SetBranchAddress("ESD",&fESD);
	  else
	    fKeyIterator = new TIter(fFile->GetListOfKeys());  
	} 

      if(fTree)
	{
	  if(AliKalmanTrack::GetConvConst()<=0.)
	    AliKalmanTrack::SetMagneticField(0.5);
	  if(fCurrentEvent>=fTree->GetEntries())
	    {
	      fCurrentDir++;
	      delete fTree;
	      fTree = 0;
	      delete fFile;
	      fFile = 0;
	      continue;
	    }
	  fTree->GetEvent(fCurrentEvent);
	} 
      else 
	{ // "old" way via ESD objects stored in root file
	  TKey* key = (TKey*)fKeyIterator->Next();
	  if (key == 0)
	    {
	      fCurrentDir++;
	      delete fKeyIterator;
	      fKeyIterator = 0;
	      delete fFile; //we have to assume there are no more ESD objects in the fFile
	      fFile = 0;
	      continue;
	    }
	  TString esdname = "ESD";
	  esdname+=fCurrentEvent;
	  if(fESD) delete fESD;
	  fESD = dynamic_cast<AliESD*>(fFile->Get(esdname));
	  if (fESD == 0)
	    {
	      Info("ReadNext","Can't find AliESD object named %s",esdname.Data());
	      fCurrentDir++;
	      delete fKeyIterator;
	      fKeyIterator = 0;
	      delete fFile;//we have to assume there is no more ESD objects in the fFile
	      fFile = 0;
	      continue;
	    }
	}
      ReadESD(fESD);
      fCurrentEvent++;
      fNEventsRead++;
      return kTRUE;//success -> read one event
    }  while(fCurrentDir < GetNumberOfDirs());
      //end of loop over directories specified in fDirs Obj Array  
  return kFALSE; //no more directories to read
}

/**********************************************************/

TFile* AliJetParticlesReaderESD::OpenFile(Int_t n)
{
  //opens fFile with kine tree

 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can't get directory name");
   return 0;
  }
 TString filename = dirname +"/"+ fESDFileName;
 TFile *ret = TFile::Open(filename.Data()); 
 if (ret == 0)
  {
    Error("OpenFiles","Can't open fFile %s",filename.Data());
    return 0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open fFile  %s",filename.Data());
    return 0;
  }
   
 return ret;
}

/**********************************************************/

void AliJetParticlesReaderESD::PrintESDtrack(const AliESDtrack *kesdtrack) const
{
  ULong_t status=kesdtrack->GetStatus();
  cout << hex << status << dec << ": ";
  if((status & AliESDtrack::kITSin) == AliESDtrack::kITSin) cout << "ITSin ";
  if((status & AliESDtrack::kITSout) == AliESDtrack::kITSout) cout << "ITSout ";
  if((status & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) cout << "ITSrefit ";
  if((status & AliESDtrack::kITSpid) == AliESDtrack::kITSpid) cout << "ITSpid ";

  if((status & AliESDtrack::kTPCin)  == AliESDtrack::kTPCin) cout << "TPCin ";
  if((status & AliESDtrack::kTPCout) == AliESDtrack::kTPCout) cout << "TPCout ";
  if((status & AliESDtrack::kTPCrefit) == AliESDtrack::kTPCrefit) cout << "TPCrefit ";
  if((status & AliESDtrack::kTPCpid) == AliESDtrack::kTPCpid) cout << "TPCpid ";

  if((status & AliESDtrack::kTRDin) == AliESDtrack::kTRDin) cout << "TRDin ";
  if((status & AliESDtrack::kTRDout) == AliESDtrack::kTRDout) cout << "TRDout ";
  if((status & AliESDtrack::kTRDrefit) == AliESDtrack::kTRDrefit) cout << "TRDrefit ";
  if((status & AliESDtrack::kTRDpid) == AliESDtrack::kTRDpid) cout << "TRDpid ";
  //  if((status & AliESDtrack::kTRDbackup) == AliESDtrack::kTRDbackup) cout << "TRDbackup ";
  if((status & AliESDtrack::kTRDStop) == AliESDtrack::kTRDStop) cout << "TRDstop ";

  if((status & AliESDtrack::kTOFin) == AliESDtrack::kTOFin) cout << "TOFin ";
  if((status & AliESDtrack::kTOFout) == AliESDtrack::kTOFout) cout << "TOFout ";
  if((status & AliESDtrack::kTOFrefit) == AliESDtrack::kTOFrefit) cout << "TOFrefit ";
  if((status & AliESDtrack::kTOFpid) == AliESDtrack::kTOFpid) cout << "TOFpid ";

  if((status & AliESDtrack::kPHOSpid) == AliESDtrack::kPHOSpid) cout << "PHOSpid ";
  if((status & AliESDtrack::kHMPIDpid) == AliESDtrack::kHMPIDpid) cout << "HMPIDpid ";
  if((status & AliESDtrack::kEMCALpid) == AliESDtrack::kEMCALpid) cout << "EMCALpid ";
  if((status & AliESDtrack::kESDpid) == AliESDtrack::kESDpid) cout << "ESDpid ";
  if((status & AliESDtrack::kTIME) == AliESDtrack::kTIME) cout << "TIME ";
  cout << endl; 
  /*
  cout << kesdtrack->GetConstrainedChi2() 
       << " " << kesdtrack->GetITSchi2()
       << " " << kesdtrack->GetTPCchi2() 
       << " " << kesdtrack->GetTRDchi2()<< endl;
  */
}
