/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id:  $ */
#include <Riostream.h>
#include "AliLog.h"
#include "TStarLight.h"
#include "starlight.h"
#include "upcevent.h"
#include "AliRunLoader.h"
#include "AliSLEventHeader.h"
#include "AliGenStarLight.h"

#include "TLorentzVector.h"

ClassImp(AliGenStarLight);

AliGenStarLight::AliGenStarLight()
  : AliGenMC()
  , fRapidityMotherMin( 1) // Max < Min: no cut
  , fRapidityMotherMax(-1)
  , fEtaChildMin( 1)       // Max < Min: no cut
  , fEtaChildMax(-1)
  , fSLgenerator(NULL)
  , fHeader(NULL)
{
}
//----------------------------------------------------------------------
AliGenStarLight::AliGenStarLight(Int_t npart)
  : AliGenMC(npart)
  , fRapidityMotherMin( 1) // Max < Min: no cut
  , fRapidityMotherMax(-1)
  , fEtaChildMin( 1)       // Max < Min: no cut
  , fEtaChildMax(-1)
  , fSLgenerator(new TStarLight("TStarLight",
				"StarLight UPC Event Generator",
				""))  // no config file name
  , fHeader(NULL)
{
  //
}
//----------------------------------------------------------------------
AliGenStarLight::~AliGenStarLight() {
  SafeDelete(fSLgenerator);
  SafeDelete(fHeader);
}
void AliGenStarLight::ImportConfigurationFromFile(const char* filename) {
  if (NULL == fSLgenerator) {
    AliFatal("AliGenStarLight class not constructed properly. ");
    return;
  }
  fSLgenerator->ImportConfigurationFromFile(filename);
}
void AliGenStarLight::SetParameter(const char* line) {
  if (NULL == fSLgenerator) {
    AliFatal("AliGenStarLight class not constructed properly. ");
    return;
  }
  fSLgenerator->SetParameter(line);
}
//----------------------------------------------------------------------
void AliGenStarLight::Init() {
  if (NULL == fSLgenerator) {
    AliFatal("AliGenStarLight class not constructed properly. ");
    return;
  }
  fSLgenerator->InitStarLight();
}
//----------------------------------------------------------------------
void AliGenStarLight::Generate() {
  Float_t vpos[4] = { 0, 0, 0, 0 };
  if (fVertexSmear == kPerEvent) {
    Vertex(); // get vertex
    for (Int_t i=0; i<3; ++i)
      vpos[i] = fVertex[i];
    vpos[3] = fTime;
  }

  Int_t   nt(0);     // number of tracks
  Bool_t  genOK(kFALSE);
  // generate events until all constraints are fulfilled
  for (Int_t trials=0; !genOK && trials < 100*1000; ++trials) {
    fSLgenerator->GenerateEvent();
    fSLgenerator->BoostEvent();
    fSLgenerator->ImportParticles(&fParticles, "ALL");

    TLorentzVector vSum;
    genOK = kTRUE;
    const Long64_t n(fParticles.GetEntries());
    if (n == 0) {
      AliFatal("no particles generated");
      return;
    }
    for (Long64_t i(0); i<n; ++i) {
      TParticle *partOrg(dynamic_cast<TParticle*>(fParticles.At(i)));
      if (NULL == partOrg) {
	AliFatal("NULL == partOrg");
	return;
      }
      //Reset of particle 4-vector in order to not have negative energy (to be removed at some point)
      //AliInfo(Form("Energy org=%.2f", partOrg->Energy()));
      TLorentzVector partVec;
      partVec.SetXYZM(partOrg->Px(), partOrg->Py(), partOrg->Pz(), partOrg->GetMass());
      TLorentzVector prodVec;
      prodVec.SetXYZT(vpos[0],vpos[1],vpos[2],vpos[3]);
      TParticle *part = new TParticle(partOrg->GetPdgCode(),partOrg->GetStatusCode(),-1,-1,partOrg->GetFirstDaughter(),partOrg->GetLastDaughter(),partVec,prodVec);
      //AliInfo(Form("Energy new=%.2f", part->Energy()));
      genOK =
	(part->Theta() >= fThetaMin) &&	(part->Theta() <  fThetaMax) &&
	(part->Phi()   >= fPhiMin)   &&	(part->Phi()   <  fPhiMax)   &&
	(part->Y()     >= fYMin)     &&	(part->Y()     <  fYMax)     &&
	(part->P()     >= fPMin)     &&	(part->P()     <  fPMax)     &&
	(part->Pt()    >= fPtMin)    &&	(part->Pt()    <  fPtMax);
      if (fEtaChildMin <= fEtaChildMax) // no cut if Max < Min
	genOK = genOK && (part->Eta() >= fEtaChildMin &&
			  part->Eta() <  fEtaChildMax);
      if (kFALSE == genOK)
	break;

      TLorentzVector v;
      part->Momentum(v);
      vSum += v;
    }

    if (fRapidityMotherMin <= fRapidityMotherMax) // no cut if Max < Min
      genOK = (genOK
	       ? (vSum.Rapidity() > fRapidityMotherMin &&
		  vSum.Rapidity() < fRapidityMotherMax)
	       : kFALSE);

    if (kFALSE == genOK) continue;

    fNprimaries = 0;
    for (Long64_t i(0), n(fParticles.GetEntries()); i<n; ++i) {
      const TParticle *part(dynamic_cast<TParticle*>(fParticles.At(i)));
      if (NULL == part) {
	AliFatal("NULL == part");
	return;
      }
      const Int_t   iparent(-1);
      const Float_t polar[3] = { 0, 0, 0 };
      const Float_t weight(trials+1);
      PushTrack(fTrackIt, iparent, part->GetPdgCode(),
		part->Px(), part->Py(), part->Pz(), part->Energy(),
		vpos[0],    vpos[1],    vpos[2],    vpos[3],
		polar[0],   polar[1],   polar[2],
		kPPrimary, nt, weight, part->GetStatusCode());
      //AliInfo(Form("weight=%.0f nt=%d fTrackIt=%d statusCode=%d", weight, nt, fTrackIt, part->GetStatusCode()));
      //part->Print();
      KeepTrack(nt);
      ++fNprimaries;
    }
    fParticles.Clear();
  }
  if (kFALSE == genOK)
    AliFatal("Maximum number of trials reached");

  SafeDelete(fHeader);

  fHeader = new AliSLEventHeader();
  const TArrayF vertexPosition(3, vpos);
  fHeader->SetPrimaryVertex(vertexPosition);
  fHeader->SetInteractionTime(vpos[3]);
  fHeader->SetNProduced(nt);
  fSLgenerator->ImportEventInfo(fHeader->GetEventInfo());
  AddHeader(fHeader);
  SetHighWaterMark(nt);
  AliRunLoader::Instance()->CdGAFile();
}
