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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Container of the distributions for the neural network                  //
//                                                                        //
// Author:                                                                //
// Alex Wilk <wilka@uni-muenster.de>                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TROOT.h>
#include <TObject.h>
#include <TMultiLayerPerceptron.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliTRDtrack.h"

#include "AliTRDCalPIDNN.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDCalPIDNN)

//_________________________________________________________________________
AliTRDCalPIDNN::AliTRDCalPIDNN()
  :AliTRDCalPID("pid", "NN PID references for TRD")
{
  //
  //  The Default constructor
  //

  Init();

}

//_________________________________________________________________________
AliTRDCalPIDNN::AliTRDCalPIDNN(const Text_t *name, const Text_t *title) 
  :AliTRDCalPID(name,title)
{
  //
  //  The main constructor
  //
  
  Init();

}

// //_____________________________________________________________________________
// AliTRDCalPIDNN::AliTRDCalPIDNN(const AliTRDCalPIDNN &c) 
//   :TNamed(c)
//   ,fMeanChargeRatio(c.fMeanChargeRatio)
//   ,fModel(0x0)
// {
//   //
//   // Copy constructor
//   //
// 
//   if (this != &c) ((AliTRDCalPIDNN &) c).Copy(*this);
//   
// }

//_________________________________________________________________________
AliTRDCalPIDNN::~AliTRDCalPIDNN()
{
  //
  // Destructor
  //

  //if (fModel) delete fModel;
  
}

//_________________________________________________________________________
Bool_t AliTRDCalPIDNN::LoadReferences(Char_t *refFile)
{
  //
  // Read the TRD Neural Networks
  //

  // Read NN Root file  
  TFile *nnFile = TFile::Open(refFile,"READ");
  if (!nnFile || !nnFile->IsOpen()) {
    AliError(Form("Opening TRD histgram file %s failed",refFile));
    return kFALSE;
  }
  gROOT->cd();

  // Read Networks
  for (Int_t imom = 0; imom < kNMom; imom++) {
    for (Int_t iplane = 0; iplane < kNPlane; iplane++) {
      TMultiLayerPerceptron *nn = (TMultiLayerPerceptron *)
         nnFile->Get(Form("NN_Mom%d_Plane%d",imom,iplane));
      fModel->AddAt(nn,GetModelID(imom,0,iplane));
    }
  }

  nnFile->Close();
  delete nnFile;

  return kTRUE;

}

//_________________________________________________________________________
TObject *AliTRDCalPIDNN::GetModel(Int_t ip, Int_t, Int_t iplane) const
{
  //
  // Returns one selected TMultiLayerPerceptron. iType not used.
  //

  if (ip<0 || ip>= kNMom) return 0x0;
  
  AliInfo(Form("Retrive MultiLayerPerceptron for %5.2f GeV/c for plane %d" 
         ,fTrackMomentum[ip]
         ,iplane));
  
  return fModel->At(GetModelID(ip, 0, iplane));

}

//_________________________________________________________________________
Double_t AliTRDCalPIDNN::GetProbability(Int_t spec, Float_t mom, Float_t *dedx
                                      , Float_t, Int_t iplane) const
{
  //
  // Core function of AliTRDCalPID class for calculating the
  // likelihood for species "spec" (see AliTRDtrack::kNspecie) of a
  // given momentum "mom", a given dE/dx (slice "dedx") yield per
  // layer in a given layer (iplane)
  //

  if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;

  // find the interval in momentum and track segment length which applies for this data
  
  Int_t imom = 1;
  while (imom<AliTRDCalPID::kNMom-1 && mom>fTrackMomentum[imom]) imom++;
  Double_t lNN1, lNN2;
  Double_t mom1 = fTrackMomentum[imom-1], mom2 = fTrackMomentum[imom];

  TMultiLayerPerceptron *nn = 0x0;
  if(!(nn = (TMultiLayerPerceptron *) fModel->At(GetModelID(imom-1, spec, iplane/*, ilength*/)))){
    //if(!(nn = (TMultiLayerPerceptron*)fModel->At(GetModelID(imom-1, spec, iplane/*, ilength*/)))){
    AliInfo(Form("Looking for mom(%f) plane(%d)", mom-1, iplane));
    AliError(Form("NN id %d not found in DB.", GetModelID(imom-1, spec, iplane)));
    return 0.;
  }

  Double_t ddedx[AliTRDtrack::kNMLPslice];
  for (int inode=0; inode<AliTRDtrack::kNMLPslice; inode++) {
    ddedx[inode] = (Double_t) dedx[inode]
                 / (AliTRDcalibDB::Instance()->GetNumberOfTimeBins()/AliTRDtrack::kNMLPslice);
  }
  lNN1 = nn->Evaluate(spec, ddedx);
  
  if(!(nn = (TMultiLayerPerceptron*)fModel->At(GetModelID(imom, spec, iplane/*, ilength*/)))){
    //if(!(nn = (TMultiLayerPerceptron*)fModel->At(GetModelID(imom, spec, iplane/*, ilength*/)))){
    AliInfo(Form("Looking for mom(%f) plane(%d)", mom, iplane));
    AliError(Form("NN id %d not found in DB.", GetModelID(imom, spec, iplane)));
    return lNN1;
  }
  lNN2 = nn->Evaluate(spec, ddedx);
  
  // return interpolation over momentum binning
  if      (mom < fTrackMomentum[0]) {
    return lNN1;
  }
  else if (mom > fTrackMomentum[AliTRDCalPID::kNMom-1]) {
    return lNN2;
  }
  else {
    return lNN1 + (lNN2 - lNN1)*(mom - mom1)/(mom2 - mom1);
  }
  
}

//_________________________________________________________________________
void AliTRDCalPIDNN::Init()
{
  //
  // Initialization
  //

  fModel = new TObjArray(AliTRDCalPID::kNPlane * AliTRDCalPID::kNMom);
  fModel->SetOwner();
  
}

//_________________________________________________________________________
Int_t AliTRDCalPIDNN::GetModelID(Int_t mom, Int_t, Int_t plane) const
{
  
  // returns the ID of the NN distribution (66 MLPs, ID from 56 to 121)

  return /*AliPID::kSPECIES * AliTRDCalPID::kNMom + */
    plane * AliTRDCalPID::kNMom + mom;  
  
}
