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
//
// Toolkit containing various usefull things
// Usable everywhere in the hfe software package
// For more information see the cxx file
//
// Authors
//   All authors of the HFE group
//
#include <TMath.h>
#include <TParticle.h>
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TAxis.h"

#include "AliAODMCParticle.h"
#include "AliESDpid.h"
#include "AliLog.h"
#include "AliTOFPIDResponse.h"

#include "AliHFEtools.h"

ClassImp(AliHFEtools)

AliESDpid *AliHFEtools::fgDefaultPID = NULL;
Int_t AliHFEtools::fgLogLevel = 0;

//__________________________________________
AliHFEtools::AliHFEtools():
  TObject()
{
}

//__________________________________________
Double_t *AliHFEtools::MakeLinearBinning(Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for linearly binned array
  //
  Double_t *bins = new Double_t[nBins + 1];
  Double_t stepsize = (ymax - ymin) / static_cast<Double_t>(nBins);
  bins[0] = ymin;
  for(Int_t ibin = 1; ibin <= nBins; ibin++)
    bins[ibin] = bins[ibin-1] + stepsize;
  return bins;
}

//__________________________________________
Double_t *AliHFEtools::MakeLogarithmicBinning(Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for logartimically binned array
  //
  Double_t *bins = new Double_t[nBins+1];
  bins[0] = ymin;
  Double_t factor = TMath::Power(ymax/ymin, 1./nBins);
  for(Int_t ibin = 1; ibin <= nBins; ibin++){
    bins[ibin] = factor * bins[ibin-1];
  }
  return bins;
}

//_________________________________________
Bool_t AliHFEtools::BinLogAxis(TObject *o, Int_t dim){

  // 
  // converts the axis (defined by the dimension) of THx or THnSparse
  // object to Log scale. Number of bins and bin min and bin max are preserved
  //


  if(!o){
    AliError("Input histogram is null pointer");
    return kFALSE;    
  }
  
  TAxis *axis = 0x0;
  if(o->InheritsFrom("TH1")){
    axis = (dynamic_cast<TH1F*>(o))->GetXaxis();
  }
  else if(o->InheritsFrom("TH2")){
    if(0 == dim){
      axis = (dynamic_cast<TH2F*>(o))->GetXaxis();
    }
    else if(1 == dim){
      axis = (dynamic_cast<TH2F*>(o))->GetYaxis();
    }
     else{
       AliError("Only dim = 0 or 1 possible for TH2F");
     }
  }
  else if(o->InheritsFrom("THnSparse")){
    axis = (dynamic_cast<THnSparse*>(o))->GetAxis(dim);
  }
  else{
    AliError("Type of input object not recognized, please check your code or update this finction");
    return kFALSE;
  }
  if(!axis){
    AliError(Form("Axis '%d' could not be identified in the object \n", dim));
    return kFALSE;
  }
  
  Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  if(from <= 0){
    AliError(Form("Log binning not possible for this axis [min = %f]\n", from));
  }
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins+1];
  newBins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i){
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete newBins;

  return kTRUE;
}

//__________________________________________
Float_t AliHFEtools::GetRapidity(TParticle *part){
  //
  // return rapidity
  //
  Float_t rapidity;
  if(!((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)) rapidity=-999;
  else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz())));
  return rapidity;
}

//__________________________________________
Float_t AliHFEtools::GetRapidity(AliAODMCParticle *part){
  // return rapidity

  Float_t rapidity;        
  if(!((part->E() - part->Pz())*(part->E() + part->Pz())>0)) rapidity=-999; 
  else rapidity = 0.5*(TMath::Log((part->E()+part->Pz()) / (part->E()-part->Pz()))); 
  return rapidity;
}

//__________________________________________
AliESDpid* AliHFEtools::GetDefaultPID(Bool_t isMC){
  //
  // Get the default PID as singleton instance
  //
  if(!fgDefaultPID){
    fgDefaultPID = new AliESDpid;
    Double_t tres = isMC ? 80. : 130.;
    fgDefaultPID->GetTOFResponse().SetTimeResolution(tres);

    // TPC Bethe Bloch parameters
    Double_t alephParameters[5];
    if(isMC){
      // simulation
      alephParameters[0] = 2.15898e+00/50.;
      alephParameters[1] = 1.75295e+01;
      alephParameters[2] = 3.40030e-09;
      alephParameters[3] = 1.96178e+00;
      alephParameters[4] = 3.91720e+00;
    } else {
      alephParameters[0] = 0.0283086/0.97;
      //alephParameters[0] = 0.0283086;
      alephParameters[1] = 2.63394e+01;
      alephParameters[2] = 5.04114e-11;
      alephParameters[3] = 2.12543e+00;
      alephParameters[4] = 4.88663e+00;
    }
    fgDefaultPID->GetTPCResponse().SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2], alephParameters[3],alephParameters[4]);
    fgDefaultPID->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);

  }
  if(fgLogLevel){
    printf("Error - You are using the default PID: You should use the PID coming from the tender\n");
    printf("Error - Arrrrrrrrr...\n");
    printf("Error - Please rethink your program logic. Using default PID is really dangerous\n");
    printf("Error - TOF PID is adapted to Monte Carlo\n");
  }
  return fgDefaultPID;
}

//__________________________________________
void AliHFEtools::DestroyDefaultPID(){
  //
  // Destroy default PID object if existing
  //
  if(fgDefaultPID) delete fgDefaultPID;
  fgDefaultPID = NULL;
}
