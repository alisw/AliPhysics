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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class creates and fills the monitor histograms for V0s              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorV0s.h"
#include "AliMonitorHisto.h"
#include "AliITSLoader.h"
#include "AliV0vertex.h"
#include "AliRunLoader.h"
#include <TFolder.h>
#include <TTree.h>
#include <TPDGCode.h>


ClassImp(AliMonitorV0s) 


//_____________________________________________________________________________
AliMonitorV0s::AliMonitorV0s()
{
// create a monitor object for V0s

}

//_____________________________________________________________________________
AliMonitorV0s::AliMonitorV0s(const AliMonitorV0s& monitor) :
  AliMonitor(monitor)
{
  Fatal("AliMonitorV0s", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorV0s& AliMonitorV0s::operator = (const AliMonitorV0s& /*monitor*/)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


//_____________________________________________________________________________
void AliMonitorV0s::CreateHistos(TFolder* folder)
{
// create the V0s monitor histograms

  fFolder = folder->AddFolder("V0s", "V0s");

  fRadius = CreateHisto1("Radius", "radius of V0 vertices", 
			 90, 0., 3., "r_{xy} [cm]", 
			 "#Delta N/N", AliMonitorHisto::kNormEvents);

  fMassK0 = CreateHisto1("MassK0", "invariant mass of K^{0} candidates", 
			 50, 0.4, 0.6, "m_{#pi^{+}#pi^{-}} [GeV/c^{2}]", 
			 "#Delta N/N", AliMonitorHisto::kNormEvents);

  fMassLambda = CreateHisto1("MassLambda", 
			     "invariant mass of #Lambda candidates", 
			     50, 1.0, 1.2, "m_{p#pi^{-}} [GeV/c^{2}]", 
			     "#Delta N/N", AliMonitorHisto::kNormEvents);

  fMassAntiLambda = CreateHisto1("MassAntiLambda", 
				 "invariant mass of #bar{#Lambda} candidates", 
				 50, 1.0, 1.2, 
				 "m_{#bar{p}#pi^{+}} [GeV/c^{2}]", 
				 "#Delta N/N", AliMonitorHisto::kNormEvents);
}


//_____________________________________________________________________________
void AliMonitorV0s::FillHistos(AliRunLoader* runLoader, 
			       AliRawReader*)
{
// fill the TPC-ITS correlation monitor histogrms

  AliITSLoader* itsLoader = (AliITSLoader*) runLoader->GetLoader("ITSLoader");
  if (!itsLoader) return;

  itsLoader->LoadV0s();
  TTree* v0s = itsLoader->TreeV0();
  if (!v0s) return;
  AliV0vertex* vertex = new AliV0vertex;
  v0s->SetBranchAddress("vertices", &vertex);

  for (Int_t i = 0; i < v0s->GetEntries(); i++) {
    v0s->GetEntry(i);
    Double_t x, y, z;
    vertex->GetXYZ(x, y, z);
    fRadius->Fill(TMath::Sqrt(x*x + y*y));
    vertex->ChangeMassHypothesis(kK0Short); 
    fMassK0->Fill(vertex->GetEffMass());
    vertex->ChangeMassHypothesis(kLambda0); 
    fMassLambda->Fill(vertex->GetEffMass());
    vertex->ChangeMassHypothesis(kLambda0Bar); 
    fMassAntiLambda->Fill(vertex->GetEffMass());
  }

  delete vertex;
  itsLoader->UnloadV0s();
}
