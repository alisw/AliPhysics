/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include "TMCProcess.h"
#include "TDatabasePDG.h"

#include "AliRun.h"
#include "AliVParticle.h"
#include "AliGenToyEventHeader.h"

#include "AliGenToy.h"

AliGenToy::AliGenToy(const std::string &name) :
  AliGenerator(),
  fHeader(nullptr),
  fNProduced(0),
  fEventWeight(1.),
  fEventVertex(3)
{
  SetName(name.c_str());
}

AliGenToy::~AliGenToy()
{

}

void AliGenToy::Init()
{

}

void AliGenToy::Generate()
{
  // prepare event header
  delete fHeader;
  fHeader = new AliGenToyEventHeader(GetName());

  fNProduced = 0;
  fEventWeight = 1.;
  fEventVertex.Reset();

  // run user generate
  UserGenerate();

  // finalize and add header
  fHeader->SetNProduced(fNProduced);
  fHeader->SetEventWeight(fEventWeight);
  fHeader->SetPrimaryVertex(fEventVertex);

  if (fContainer)
    fContainer->AddHeader(fHeader);
  else
    gAlice->SetGenEventHeader(fHeader);
}

void AliGenToy::SetCentrality(Double_t cent)
{
  fHeader->SetCentrality(cent);
}

Int_t AliGenToy::AddParticle(const AliVParticle &part)
{
  const Int_t done = 0;
  const Int_t parent = 0;
  const TMCProcess mech = kPPrimary;
  const Float_t weight = 1.;
  const Int_t is = 0;

  Int_t ntr;

  PushTrack(done, parent, part.PdgCode(),
            part.Px(), part.Py(), part.Pz(), part.E(),
            0., 0., 0., 0.,
            0., 0., 0.,
            mech, ntr, weight, is);
  ++fNProduced;

  return ntr;
}

Int_t AliGenToy::AddParticle(const TLorentzVector &part)
{
  const Int_t done = 0;
  const Int_t parent = 0;
  const Int_t pdg = 0;
  const TMCProcess mech = kPPrimary;
  const Float_t weight = 1.;
  const Int_t is = 0;

  Int_t ntr;

  PushTrack(done, parent, pdg,
            part.Px(), part.Py(), part.Pz(), part.E(),
            0., 0., 0., 0.,
            0., 0., 0.,
            mech, ntr, weight, is);
  ++fNProduced;

  return ntr;
}

Int_t AliGenToy::AddParticle(Double_t px, Double_t py, Double_t pz, Int_t pdg)
{
  const Int_t done = 0;
  const Int_t parent = 0;
  const TMCProcess mech = kPPrimary;
  const Float_t weight = 1.;
  const Int_t is = 0;

  const TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(pdg);
  const Float_t m = partPDG ? partPDG->Mass() : 0.;
  const Float_t e = TMath::Sqrt(px*px + py*py + pz*pz + m*m);

  Int_t ntr;

  PushTrack(done, parent, pdg,
            px, py, pz, e,
            0., 0., 0., 0.,
            0., 0., 0.,
            mech, ntr, weight, is);
  ++fNProduced;

  return ntr;
}
