// $Id$

//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliJetParticlesReaderKineGoodTPC
//
// Reader for Good TPC tracks
//
// loizides@ikf.uni-frankfurt.de
//
/////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TParticle.h>
#include <TLorentzVector.h>
#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenHijingEventHeader.h>
#include "AliJetParticle.h"
#include "AliJetEventParticles.h"
#include "AliJetParticlesReaderKineGoodTPC.h"

ClassImp(AliJetParticlesReaderKineGoodTPC)

AliJetParticlesReaderKineGoodTPC::AliJetParticlesReaderKineGoodTPC() :
  AliJetParticlesReader(),
  fFileName("good_tracks_tpc"),
  fInput(0)
{
  //constructor
}

AliJetParticlesReaderKineGoodTPC::AliJetParticlesReaderKineGoodTPC(TString& fname) :
  AliJetParticlesReader(),
  fFileName(fname),
  fInput(0)
{
  //constructor
}

AliJetParticlesReaderKineGoodTPC::AliJetParticlesReaderKineGoodTPC(TObjArray* dirs, const Char_t *filename):
  AliJetParticlesReader(dirs),
  fFileName(filename),
  fInput(0)
{
  //constructor
}

AliJetParticlesReaderKineGoodTPC::~AliJetParticlesReaderKineGoodTPC()
{
  //destructor
  if(fInput) delete fInput;
}

void AliJetParticlesReaderKineGoodTPC::Rewind()
{
  //Rewinds to the beginning
  if(fInput) delete fInput;
  fInput       = 0;
  fCurrentDir  = 0;
  fNEventsRead = 0;  
}

Int_t AliJetParticlesReaderKineGoodTPC::ReadNext()
{
  //Reads good_tpc_tracks file
  if((!fOwner) || (fEventParticles == 0)) 
    fEventParticles = new AliJetEventParticles();
  else
     fEventParticles->Reset();

  while(fCurrentDir < GetNumberOfDirs())
    { 
      if(!OpenFile(fCurrentDir)) 
      { 
	delete fInput; //close current session
	fInput = 0;    //assure pointer is null
	fCurrentDir++;
	continue;
      }
    
      Info("ReadNext","Reading Event %d",fCurrentDir*1000+fCurrentEvent);

      Int_t label,code;
      Float_t px,py,pz,x,y,z;
      Int_t i=0;
      while (*fInput>>label>>code>>px>>py>>pz>>x>>y>>z)
      {
	const Float_t kp2=px*px+py*py+pz*pz;
	if(kp2<1e-3)  continue;
	const Float_t kpt=TMath::Sqrt(px*px+py*py);
	const Float_t kp=TMath::Sqrt(kp2);
	const Float_t keta=0.5*TMath::Log((kp+pz+1e-30)/(kp-pz+1e-30)); 
	const Float_t kphi=TMath::Pi()+TMath::ATan2(-py,-px);
	//cout << i << " " << label << " " << px << " " << py << " " << pz << endl;
	if(IsAcceptedParticle(kpt,kphi,keta))
	  fEventParticles->AddParticle(px,py,pz,kp,i++,label,code,kpt,kphi,keta);
      }

      fCurrentEvent++;
      fNEventsRead++;
      return kTRUE;

    }
  return kFALSE;
}

Int_t AliJetParticlesReaderKineGoodTPC::OpenFile(Int_t n)
{
  if(fInput){ //change here if you want to 
              //support more than one good file per dir
    return kFALSE;
  }

  const TString& dirname = GetDirName(n);
  if (dirname == "")
    { 
      Error("OpenNextFile","Can't get directory name with index %d",n);
      return kFALSE;
    }

  TString filename = dirname +"/"+ fFileName;
  fInput=new ifstream(filename);
#if defined(__HP_aCC) || defined(__DECCXX)
  if ( fInput->rdbuf()->is_open() == 0)
#else
  if ( fInput->is_open() == 0)
#endif
    {
      Error("OpenNextFile","Can't open session from file %s",filename.Data());
      return kFALSE;
    }

  fCurrentEvent = 0;
  return kTRUE;
}
