// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelectionD0.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionD0.h"
#include "AliVParticle.h"
#include "AliAODRecoDecayHF2Prong.h" // libPWGHFvertexingHF
#include "TObjArray.h"
#include "THnSparse.h"
#include "TAxis.h"
#include <iostream>
#include <cerrno>
#include <memory>

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionD0)

AliDxHFEParticleSelectionD0::AliDxHFEParticleSelectionD0(const char* opt)
  : AliDxHFEParticleSelection("D0", opt)
  , fD0Properties(NULL)
{
  // constructor
  // 
  // 
  // 
  // 
}

AliDxHFEParticleSelectionD0::~AliDxHFEParticleSelectionD0()
{
  // destructor
}

int AliDxHFEParticleSelectionD0::InitControlObjects()
{
  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly
  TString name;
  const int thnSize = 4;
  const double pi=TMath::Pi();

  // TODO: very specific D0 for the moment, sort out later
  // TODO: theta?
  // 			         0    1     2   3
  // 			      mass   Pt   Phi Ptbin
  int    thnBins[thnSize] = {   200,   1000,  100, 100};
  double thnMin [thnSize] = {  1.5648,   0,    0,   0};
  double thnMax [thnSize] = {  2.1648, 100,  (2*pi), 100};

  name.Form("%s info", GetName());
  std::auto_ptr<THnSparseF> D0Properties(new THnSparseF(name, name, thnSize, thnBins, thnMin, thnMax));

  if (D0Properties.get()==NULL) {
    return -ENOMEM;
  }
  int axis=0;
  D0Properties->GetAxis(axis++)->SetTitle("mass");
  D0Properties->GetAxis(axis++)->SetTitle("Pt");
  D0Properties->GetAxis(axis++)->SetTitle("Phi"); 
  D0Properties->GetAxis(axis++)->SetTitle("Ptbin"); 

  fD0Properties=D0Properties.release();
  AddControlObject(fD0Properties);

  return AliDxHFEParticleSelection::InitControlObjects();
}

int AliDxHFEParticleSelectionD0::HistogramParticleProperties(AliVParticle* p, bool selected)
{
  /// histogram particle properties
  if (!p) return -EINVAL;

  // fill the common histograms
  AliDxHFEParticleSelection::HistogramParticleProperties(p, selected);

  // TODO: histograms for all and selected particles
  if (!selected) return 0;

  // TODO: find out which type is necessary
  AliAODRecoDecayHF2Prong* part=dynamic_cast<AliAODRecoDecayHF2Prong*>(p);
  if (part) {
    Double_t invmassD0 = part->InvMassD0();
    // TODO: use cut object to define pt bin
    Int_t ptbin=0;//cuts->PtBin(part->Pt());
    Double_t D0Stuff[] = {invmassD0,part->Pt(),part->Phi(),ptbin};
    if (fD0Properties) fD0Properties->Fill(D0Stuff);
  }

  return 0;
}

bool AliDxHFEParticleSelectionD0::IsSelected(AliVParticle* /*p*/)
{
  /// TODO: implement specific selection of D0 candidates
  return true;
}
