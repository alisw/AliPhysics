// $Id$

//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliJetParticlesReaderKine
//
// Reader for Kinematics
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
#include "AliJetParticlesReaderKine.h"

ClassImp(AliJetParticlesReaderKine)

AliJetParticlesReaderKine::AliJetParticlesReaderKine() :
  AliJetParticlesReader(),
  fFileName("galice.root"),
  fRunLoader(0),
  fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE),
  fUseTracks(kFALSE)
{
  //constructor
}

AliJetParticlesReaderKine::AliJetParticlesReaderKine(TString& fname) :
  AliJetParticlesReader(),
  fFileName(fname),
  fRunLoader(0),
  fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE),
  fUseTracks(kFALSE)
{
  //constructor
}

AliJetParticlesReaderKine::AliJetParticlesReaderKine(TObjArray* dirs, const Char_t *filename):
  AliJetParticlesReader(dirs),
  fFileName(filename),
  fRunLoader(0),
  fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE),
  fUseTracks(kFALSE)
{
  //constructor
}

AliJetParticlesReaderKine::~AliJetParticlesReaderKine()
{
  //destructor
  if(fRunLoader) delete fRunLoader;
}

void AliJetParticlesReaderKine::Rewind()
{
  //Rewinds to the beginning
  if(fRunLoader) delete fRunLoader;
  fRunLoader   = 0;
  fCurrentDir  = 0;
  fNEventsRead = 0;  
}

Int_t AliJetParticlesReaderKine::ReadNext()
{
  //Reads Kinematics Tree

  if((!fOwner) || (fEventParticles == 0)) 
    fEventParticles = new AliJetEventParticles();

  while(fCurrentDir < GetNumberOfDirs())
    { 

      if (!OpenFile(fCurrentDir)) 
      { 
	delete fRunLoader; //close current session
	fRunLoader = 0;    //assure pointer is null
	fCurrentDir++;
	continue;
      }
    
      if (fCurrentEvent == fRunLoader->GetNumberOfEvents())
	{
	  //read next directory
	  delete fRunLoader; //close current session
	  fRunLoader = 0;    //assure pointer is null
	  fCurrentDir++;     //go to next dir
	  continue; 
	}
     
      Info("ReadNext","Reading Event %d",fCurrentEvent);
      fRunLoader->GetEvent(fCurrentEvent);
      AliStack* stack = fRunLoader->Stack();
      if (!stack)
	{
	  Error("ReadNext","Can't get stack for event %d",fCurrentEvent);
	  continue;
	}

      //clear old event
      Int_t nprim = stack->GetNprimary();
      Int_t npart = nprim;
      if(fUseTracks)
	npart = stack->GetNtrack();
      fEventParticles->Reset(npart);

      TString headdesc="";
      Int_t gentype=-1;
      AliHeader *header=fRunLoader->GetHeader();
      if(!header) {
	Warning("ReadNext","Header not found in event %d",fCurrentEvent);
      } else {
	AliGenEventHeader *genheader=header->GenEventHeader();
	TString hname=genheader->GetName();
	if(hname.CompareTo("HIJINGparam")==0) {
	  gentype=1;	  ;
	  headdesc+="HIJINGparam";
	} else if(hname.CompareTo("Hijing")==0) {
	  gentype=2;
	  AliGenHijingEventHeader *hheader=(AliGenHijingEventHeader*)header->GenEventHeader();
	  headdesc+="Hijing";
	  if(!hheader) {
	    Warning("ReadNext","Hijing-Header not found in event %d",fCurrentEvent);
	  } else {
	    Int_t ntrials=hheader->Trials();
	    fEventParticles->SetTrials(ntrials);
	    fEventParticles->SetImpact(hheader->ImpactParameter());
	    fEventParticles->SetNhard(hheader->HardScatters());
	    fEventParticles->SetNpart(hheader->NwNw());

	    TLorentzVector  jet1; 
	    TLorentzVector  jet2; 
	    TLorentzVector  jetFsr1;
	    TLorentzVector  jetFsr2;
	    hheader->GetJets(jet1,jet2,jetFsr1,jetFsr2);
	    fEventParticles->AddHard(0,jet1.X(),jet1.Y(),jet1.Z(),jet1.T(),0);
	    fEventParticles->AddHard(1,jet2.X(),jet2.Y(),jet2.Z(),jet2.T(),0);
	  }
	} else if(hname.CompareTo("Pythia")==0){
	  gentype=3;
	  headdesc+="Pythia";
	  AliGenPythiaEventHeader *pheader=(AliGenPythiaEventHeader*)header->GenEventHeader();
	  if(!pheader) {
	    Warning("ReadNext","Pythia-Header not found in event %d",fCurrentEvent);
	  } else {
	    Int_t ntruq=0;
#ifndef NOUQHEADERINFO
	    ntruq=pheader->NUQTriggerJets();
	    if(ntruq){
	      Double_t x0=pheader->GetXJet();
	      Double_t y0=pheader->GetYJet();
	      if(x0==y0==-1){
		x0=y0=0.;
	      }
	      Double_t zquench[4];
	      pheader->GetZQuench(zquench);
	      fEventParticles->SetXYJet(x0,y0);
	      fEventParticles->SetZQuench(zquench);
	      for(Int_t j=0;j<ntruq;j++){
		Float_t pjet[4];
		pheader->UQJet(j,pjet);
		fEventParticles->AddUQJet(pjet);
	      }
	    }
#endif
	    //Int_t ptyp=pheader->ProcessType();
	    Int_t ntrials=pheader->Trials();
	    fEventParticles->SetTrials(ntrials);

	    Int_t ntr=pheader->NTriggerJets();
	    if(ntr){
	      for(Int_t j=0;j<ntr;j++){
		Float_t pjet[4];
		pheader->TriggerJet(j,pjet);
		fEventParticles->AddJet(pjet);
		if(!ntruq) fEventParticles->AddUQJet(pjet);
	      }
	    }
	    for(Int_t i=6;i<=7;i++){
	      TParticle *MP = stack->Particle(i);
	      if(!MP) break;
	      Int_t type=0;
	      if(MP->GetPdgCode()==21) type=2;
	      else type=1;
	      fEventParticles->AddHard(i-6,MP->Px(),MP->Py(),MP->Pz(),MP->Energy(),type);
	    }
	  }
	}
	headdesc+=" Run ";
	headdesc+=header->GetRun();
	headdesc+=": Ev ";
	headdesc+=header->GetEventNrInRun();
      }
      fEventParticles->SetHeader(headdesc);

      //get vertex
      const TParticle *kv = stack->Particle(0);
      if(kv) {
	fEventParticles->SetVertex(kv->Vx(),kv->Vy(),kv->Vz());
      }

      //loop over particles
      for (Int_t i = 0;i<npart; i++)
	{
	  TParticle *p = stack->Particle(i);
	  if(!p) continue;
	  Int_t child1 = p->GetFirstDaughter();
	  //Int_t child2 = p->GetLastDaughter();	
	  //Int_t mother = p->GetFirstMother();	   
	  //cout << child1 << " " << child2 << " " << mother << endl;
	  if((child1>=0) && (child1<nprim)) continue; 

	  //check status code depending on gentype
	  if(gentype==1){ // Hijing Param
	    if(p->GetStatusCode()!=0) continue;
	  } else if(gentype==2){ //Hijing
	    //if( p->GetStatusCode()!=0){
	    //  cout <<p->GetStatusCode() << " ";p->Print(); 
	    //}
	  } else if(gentype==3){ //Pythia
	    if(p->GetStatusCode()!=1) continue;
	    p->SetWeight(-123); //mark particle
	  }

	  //kinematic cuts
	  if(IsAcceptedParticle(p)){ //put particle in event
	    //if(p->Pt()>20){cout <<p->GetStatusCode() << " ";p->Print(); }
	    fEventParticles->AddParticle(p,i); 
	  }
	}
      fCurrentEvent++;
      fNEventsRead++;
      return kTRUE;
    }

  //end of loop over directories specified in fDirs ObjArray
  return kFALSE;
}

Int_t AliJetParticlesReaderKine::OpenFile(Int_t n)
{
  //opens file with kine tree

  if(fRunLoader){
    if(fCurrentEvent < fRunLoader->GetNumberOfEvents()) return kTRUE;
    else return kFALSE;
  }

  const TString& dirname = GetDirName(n);
  if (dirname == "")
    { 
      Error("OpenNextFile","Can't get directory name with index %d",n);
      return kFALSE;
    }

  TString filename = dirname +"/"+ fFileName;
  TFile file(filename);
  if ( file.IsOpen() == 0)
    {
      Error("OpenNextFile","Can't open session from file %s",filename.Data());
      return kFALSE;
    }
  file.Close();
  fRunLoader = AliRunLoader::Open(filename.Data()); 
  
  if (fRunLoader->GetNumberOfEvents() <= 0)
    {
      Error("OpenNextFile","There is no events in this directory.");
      delete fRunLoader;
      fRunLoader = 0;
      return kFALSE;
    }

  if (fRunLoader->LoadKinematics())
    {
      Error("OpenNextFile","Error occured while loading kinematics.");
      return kFALSE;
    }

  fCurrentEvent = 0;
  return kTRUE;
}
