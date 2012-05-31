//
// AliGenEpos.cxx
//
//  ALICE event generator based on EPOS model from Klaus Werner
//
//  Created on: Feb 28, 2009
//      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
//

#include "AliGenEpos.h"
#include "TEpos.h"
#include "TParticle.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliGenEposEventHeader.h"

ClassImp(AliGenEpos)

AliGenEpos::AliGenEpos() : AliGenMC(),
		fBmin(0),
		fBmax(10000),
		fPhiMin(0),
		fPhiMax(TMath::TwoPi()),
		fFilterModelOutput(kFALSE) {
	SetMC(new TEpos());
}

AliGenEpos::AliGenEpos(Int_t npart) : AliGenMC(npart),
		fBmin(0),
		fBmax(10000),
		fPhiMin(0),
		fPhiMax(TMath::TwoPi()),
		fFilterModelOutput(kFALSE) {
	SetMC(new TEpos());
}

void AliGenEpos::Init() {
  // Sets up TEpos
	AliGenMC::Init();
	TEpos *epos = GetTEpos();
	epos->SetLaproj(this->fZProjectile);
	epos->SetMaproj(this->fAProjectile);
	epos->SetLatarg(this->fZTarget);
	epos->SetMatarg(this->fATarget);
	epos->SetPhimin(this->fPhiMin);
	epos->SetPhimax(this->fPhiMax);
	epos->SetBminim(this->fBmin);
	epos->SetBmaxim(this->fBmax);
	epos->SetEcms(this->fEnergyCMS);
	epos->SetSplitting(kTRUE);
	GetTEpos()->Initialize();
}

void AliGenEpos::Generate() {
  // Does actual generation and output conversion
	  Float_t polar[3]    =   {0,0,0};
	  Float_t origin0[3]  =   {0,0,0};
	  Float_t origin[3]   =   {0,0,0};
	  Float_t time0 = 0.;
	  Float_t time  = 0.;
	  fNprimaries = 0;
	  Int_t nt  = 0; //output parameter for PushTrack

	  Vertex();
	  for (int j=0; j < 3; j++) origin0[j] = fVertex[j];
	  time0 = fTime;

	  // Generate one event

	  GetTEpos()->GenerateEvent();
	  AliWarning("Generated");
	  GetTEpos()->ImportParticles(&fParticles);

	  Int_t np = fParticles.GetEntriesFast();
	  AliWarning(Form("Imported %d particles", np));

	  Int_t *idsOnStack = NULL;
	  idsOnStack = new Int_t[np];
	  for (int i = 0; i < np; i++) idsOnStack[i] = 0;
	  TParticle *iparticle;

	  for (int i = 0; i < np; i++) {
		  iparticle = (TParticle *) fParticles.At(i);
		  //Bool_t isNullEntry = iparticle->GetStatusCode() == 0;
		  //Bool_t isCommentOrUnknown = iparticle->GetStatusCode() > 2;
		  Bool_t hasDecayed = iparticle->GetStatusCode() >= 2;
		  Bool_t isFinalState = iparticle->GetStatusCode() == 1;
		  Int_t imo = iparticle->GetFirstMother();
		  Bool_t  hasMother = (imo >=0);


		  if (isFinalState) {
			  origin[0] = iparticle->Vx();
			  origin[1] = iparticle->Vy();
			  origin[2] = iparticle->Vz();
			  time      = iparticle->T();
			  //doubled track with freeze out coordinates for femtoscopy
			  PushTrack(0,
					  imo>=0?idsOnStack[imo]:-1,
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
			        time,
				polar[0],polar[1],polar[2],
				hasMother ? kPDecay:kPNoProcess,nt);

		      idsOnStack[i] = nt;
		      fNprimaries++;
		      KeepTrack(nt);

		      //real track with smeared vertex
		      origin[0] += origin0[0];
		      origin[1] += origin0[1];
		      origin[2] += origin0[2];
		      time      += time0;
			  PushTrack(1,
					  nt,   //doubled track as mother
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
				time,
				polar[0],polar[1],polar[2],
				kPDecay,nt);
		      fNprimaries++;
		      KeepTrack(nt);
		  } else if(hasDecayed || !fFilterModelOutput) {
			  // don't track it and don't smear vertex
			  origin[0] = iparticle->Vx();
			  origin[1] = iparticle->Vy();
			  origin[2] = iparticle->Vz();
			  time      = iparticle->T();
			  PushTrack(0,
					  imo>=0?idsOnStack[imo]:-1,
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
				time,
				polar[0],polar[1],polar[2],
				hasMother ? kPDecay:kPNoProcess,nt);
		      idsOnStack[i] = nt;
		      fNprimaries++;
		      KeepTrack(nt);
		  } else if(fFilterModelOutput) {
			  //filtered internal model objects
			  idsOnStack[i] = -1; // to erase mother field in its dauthers
		  }
	  }
	  SetHighWaterMark(fNprimaries);
	  TArrayF eventVertex;
	  eventVertex.Set(3);
	  eventVertex[0] = origin0[0];
	  eventVertex[1] = origin0[1];
	  eventVertex[2] = origin0[2];
	// Builds the event header, to be called after each event
	  AliGenEposEventHeader* header = new AliGenEposEventHeader("EPOS");

	  header->SetNProduced(fNprimaries);
	  header->SetPrimaryVertex(eventVertex);
	  header->SetInteractionTime(time0);

	  header->SetImpactParameter(GetTEpos()->GetBimevt());
	  header->SetReactionPlaneAngle(GetTEpos()->GetPhievt());

	  header->SetHardScatters(0);
	  header->SetParticipants(GetTEpos()->GetNpjevt(), GetTEpos()->GetNtgevt());
	  header->SetCollisions(GetTEpos()->GetKolevt(), 0, 0, GetTEpos()->GetNpjevt());
	  header->SetSpectators(GetTEpos()->GetJpnevt(), GetTEpos()->GetJppevt(), GetTEpos()->GetJtnevt(), GetTEpos()->GetJtpevt());

	// Event Vertex
	  header->SetPrimaryVertex(fVertex);
	  
	  header->SetBimevt(GetTEpos()->GetBimevt());
	  header->SetPhievt(GetTEpos()->GetPhievt());
	  header->SetKolevt(GetTEpos()->GetKolevt());
	  header->SetKoievt(GetTEpos()->GetKoievt());
	  header->SetPmxevt(GetTEpos()->GetPmxevt());
	  header->SetEgyevt(GetTEpos()->GetEgyevt());
	  header->SetNpjevt(GetTEpos()->GetNpjevt());
	  header->SetNtgevt(GetTEpos()->GetNtgevt());
	  header->SetNpnevt(GetTEpos()->GetNpnevt());
	  header->SetNppevt(GetTEpos()->GetNppevt());
	  header->SetNtnevt(GetTEpos()->GetNtnevt());
	  header->SetNtpevt(GetTEpos()->GetNtpevt());
	  header->SetJpnevt(GetTEpos()->GetJpnevt());
	  header->SetJppevt(GetTEpos()->GetJppevt());
	  header->SetJtnevt(GetTEpos()->GetJtnevt());
	  header->SetJtpevt(GetTEpos()->GetJtpevt());
	  header->SetXbjevt(GetTEpos()->GetXbjevt());
	  header->SetQsqevt(GetTEpos()->GetQsqevt());
	  header->SetNglevt(GetTEpos()->GetNglevt());
	  header->SetZppevt(GetTEpos()->GetZppevt());
	  header->SetZptevt(GetTEpos()->GetZptevt());
    AddHeader(header);
    fCollisionGeometry = (AliGenEposEventHeader*)  header;
    
    delete[] idsOnStack;
}

AliGenEpos::~AliGenEpos() {
}

