/*
 * AliGenEpos.cpp
 *
 *  ALICE event generator based on EPOS model from Klaus Werner
 *
 *  Created on: Feb 28, 2009
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */

#include "AliGenEpos.h"
#include "TEpos.h"
#include "TParticle.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliGenEposEventHeader.h"

ClassImp(AliGenEpos)

AliGenEpos::AliGenEpos() : AliGenMC(),
		fBmin(0),
		fBmax(10000),
		fPhiMin(0),
		fPhiMax(2*3.1415926),
		fFilterModelOutput(kFALSE) {
	SetMC(new TEpos());
}

AliGenEpos::AliGenEpos(Int_t npart) : AliGenMC(npart),
		fBmin(0),
		fBmax(10000),
		fPhiMin(0),
		fPhiMax(2*3.1415926),
		fFilterModelOutput(kFALSE) {
	SetMC(new TEpos());
}

void AliGenEpos::Init() {
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
	GetTEpos()->Initialize();
}

void AliGenEpos::Generate() {
	  Float_t polar[3]   =   {0,0,0};
	  Float_t origin0[3]  =   {0,0,0};
	  Float_t origin[3]   =   {0,0,0};

	  Int_t nt  = 0; //output parameter for PushTrack

	  Vertex();
	  for (int j=0; j < 3; j++) origin0[j] = fVertex[j];

	  // Generate one event

	  GetTEpos()->GenerateEvent();
	  AliWarning("Generated");
	  GetTEpos()->ImportParticles(&fParticles);

	  Int_t np = fParticles.GetEntriesFast();
	  AliWarning(Form("Imported %d particles", np));

	  Int_t *idsOnStack = NULL;
	  idsOnStack = new Int_t[np];
	  TParticle *iparticle;

	  for (int i = 0; i < np; i++) {
		  iparticle = (TParticle *) fParticles.At(i);
		  //Bool_t isNullEntry = iparticle->GetStatusCode() == 0;
		  //Bool_t isCommentOrUnknown = iparticle->GetStatusCode() > 2;
		  Bool_t hasDecayed = iparticle->GetStatusCode() == 2;
		  Bool_t isFinalState = iparticle->GetStatusCode() == 1;
		  Int_t imo = iparticle->GetFirstMother();
		  Bool_t  hasMother = (imo >=0);


		  if (isFinalState) {
			  origin[0] = iparticle->Vx();
			  origin[1] = iparticle->Vy();
			  origin[2] = iparticle->Vz();
			  //doubled track with freeze out coordinates for femtoscopy
			  PushTrack(0,
					  imo>=0?idsOnStack[imo]:-1,
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
				iparticle->T(),
				polar[0],polar[1],polar[2],
				hasMother ? kPDecay:kPNoProcess,nt);

		      idsOnStack[i] = nt;
		      fNprimaries++;
		      KeepTrack(nt);

		      //real track with smeared vertex
		      origin[0] += origin0[0];
		      origin[1] += origin0[1];
		      origin[2] += origin0[2];
			  PushTrack(1,
					  nt,   //doubled track as mother
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
				iparticle->T(),
				polar[0],polar[1],polar[2],
				kPDecay,nt);
		      fNprimaries++;
		      KeepTrack(nt);
		  } else if(hasDecayed || !fFilterModelOutput) {
			  // don't track it and don't smear vertex
			  origin[0] = iparticle->Vx();
			  origin[1] = iparticle->Vy();
			  origin[2] = iparticle->Vz();
			  PushTrack(0,
					  imo>=0?idsOnStack[imo]:-1,
					  iparticle->GetPdgCode(),
				iparticle->Px(),iparticle->Py(),iparticle->Pz(),iparticle->Energy(),
				origin[0], origin[1], origin[2],
				iparticle->T(),
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

	  header->SetImpactParameter(GetTEpos()->GetBimevt());
	  header->SetReactionPlaneAngle(GetTEpos()->GetPhievt());

	  header->SetHardScatters(0);
	  header->SetParticipants(GetTEpos()->GetNpjevt(), GetTEpos()->GetNtgevt());
	  header->SetCollisions(GetTEpos()->GetKolevt(), 0, 0, GetTEpos()->GetNpjevt());
	  header->SetSpectators(GetTEpos()->GetJpnevt(), GetTEpos()->GetJppevt(), GetTEpos()->GetJtnevt(), GetTEpos()->GetJtpevt());

	// Event Vertex
	  header->SetPrimaryVertex(fVertex);
	  AddHeader(header);
	  fCollisionGeometry = (AliGenEposEventHeader*)  header;

	  delete[] idsOnStack;
}

AliGenEpos::~AliGenEpos() {
}

