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
#include "TDime.h"
#include "AliRunLoader.h"
#include "AliGenEventHeader.h"
#include "AliGenDime.h"
#include "AliGenMC.h"

ClassImp(AliGenDime);

AliGenDime::AliGenDime()
  : AliGenMC()
  , fDMgenerator(NULL) {
}

//----------------------------------------------------------------------
AliGenDime::AliGenDime(Int_t npart)
  : AliGenMC(npart)
  , fDMgenerator(new TDime())
{
//
}

//----------------------------------------------------------------------
AliGenDime::~AliGenDime() {
  if (NULL != fDMgenerator) delete fDMgenerator;
  fDMgenerator = NULL;
}

//----------------------------------------------------------------------
void AliGenDime::Init() {
  if (NULL == fDMgenerator) {
    AliFatal("AliGenDime class not constructed properly. ");
    return;
  }
  fDMgenerator->Initialize();
}

//----------------------------------------------------------------------
void AliGenDime::Generate() {

  // Generate one event -->

  const Float_t polar[3] = {0,0,0};
  Float_t vpos[4]  = {0,0,0,0};

  // Set collision vertex position
  if (fVertexSmear == kPerEvent) {
      Vertex();
      for (Int_t i = 0; i < 3; ++i) {
        vpos[i] = fVertex[i];
      }
      vpos[3] = fTime;
  }

  Int_t nt = 0; // number of tracks

  // Generate events, store into fParticles
  fDMgenerator->GenerateEvent();
  fDMgenerator->ImportParticles(&fParticles, "All");

  const Int_t   iparent = -1;
  const Float_t weight  = 1.0;

  fNprimaries = 0;

  // Loop over particles
  for (Int_t i = 0; i < fParticles.GetEntries(); ++i) {

      const TParticle* part = (TParticle*)fParticles.At(i);
      if (part == NULL) {
        AliFatal("part == NULL");
        return;
      }

      if (part->GetStatusCode() == 1) { // Push final states only

        PushTrack(fTrackIt, iparent, part->GetPdgCode(),
                part->Px(), part->Py(), part->Pz(), part->Energy(),
	        vpos[0],    vpos[1],    vpos[2],    vpos[3],
	        polar[0],   polar[1],   polar[2],
	        kPPrimary, nt, weight, part->GetStatusCode());

        AliInfo(Form("weight=%.0f nt=%d fTrackIt=%d statusCode=%d",
              weight, nt, fTrackIt, part->GetStatusCode()));
        part->Print();
        KeepTrack(nt);
        ++fNprimaries;
      }
  }
  // Clear particles for the next event
  fParticles.Clear();

  // Make header
  AliGenEventHeader* header = new AliGenEventHeader("Dime");
  const TArrayF vertexPosition(3, vpos);
  header->SetPrimaryVertex(vertexPosition);
  header->SetInteractionTime(vpos[3]);
  header->SetNProduced(nt);

  // Pass header
  AddHeader(header);
  SetHighWaterMark(nt);
  AliRunLoader::Instance()->CdGAFile();

}
