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

/* $Id: AliGenTunedOnPbPb.cxx 51126 2013-08-19 13:37:49Z fnoferin $ */

// Parameterisation based on 5.5 ATeV PbPb data
// pi, K, p, neutron, K0, lambda, phi, Xi, Omega spectra, v2, v3 (no jets!)
// Author: fnoferin@cern.ch

#include <TArrayF.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TVirtualMC.h>
#include <TLorentzVector.h>

#include "AliConst.h"
#include "AliDecayer.h"
#include "AliGenEventHeaderTunedPbPb.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliGenTunedOnPbPb.h"

ClassImp(AliGenTunedOnPbPb)

  // set default parameters for 10-20% centrality
Int_t AliGenTunedOnPbPb::fgPdgInput[fgNspecies] = {211,-211,111,321,-321,2212,-2212,310,3122,-3122,333,3312,-3312,3334,-3334,2112,-2112};
Float_t AliGenTunedOnPbPb::fgMult[fgNspecies] = {450,450,450,70,70,21,21,70,20,20,8,2.4,2.4,0.4,0.4,21,21};

Float_t AliGenTunedOnPbPb::fgV3Overv2 = 6.25000000000000000e-01;
Float_t AliGenTunedOnPbPb::fgEventplane=0;

TF1 *AliGenTunedOnPbPb::fgV2 = NULL;

//_____________________________________________________________________________
AliGenTunedOnPbPb::AliGenTunedOnPbPb()
  :AliGenerator(),
  fCmin(0.),
  fCmax(100.),
  fChangeWithCentrality(kFALSE),
  fYMaxValue(2.0),
  fYlimitForFlatness(2.0),
  fYdecreaseSp(0.2),
  fYdecreaseV2(0.2)
{
  //
  // Default constructor
  //
  SetCutVertexZ();
  SetPtRange();

  for(Int_t i=0;i < fgNspecies;i++){
    fgHSpectrum[i] = NULL;
    fgHv2[i] = NULL;
  }
}

//_____________________________________________________________________________
AliGenTunedOnPbPb::~AliGenTunedOnPbPb()
{
  //
  // Standard destructor
  //
}

//_____________________________________________________________________________
void AliGenTunedOnPbPb::Init()
{
  //
  // Initialise the generator
  //

  // define histos
}


//_____________________________________________________________________________
void AliGenTunedOnPbPb::Generate()
{
  //
  // Generate one trigger
  //

  Float_t avCentr = (fCmin+fCmax)*0.5;

  Float_t centrality = avCentr;

  if(fChangeWithCentrality) centrality = fCmin + gRandom->Rndm()*(fCmax-fCmin);

  SetParameters(centrality);

  if(!fChangeWithCentrality){
    Float_t in=0;
    for(Int_t i=0;i < fgNspecies;i++){
      in=0;
      if(fgHSpectrum[i]){
        for(Int_t j=1;j<=fgHSpectrum[i]->GetNbinsX();j++){
          in += fgHSpectrum[i]->GetBinContent(j)*fgHSpectrum[i]->GetBinWidth(j);
        }
      }

      // replace n-particles with the one in input file if centralidy dependece was disable
      fgMult[i] = in;
    }
  }


  TMCProcess statusPdg[fgNspecies] = {kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary,kPPrimary};

  Float_t parV2scaling[3] = {1,0.202538,-0.00214468};

  Float_t scaleV2 = 1.0;

  TDatabasePDG *pdgD = TDatabasePDG::Instance();

  if(fChangeWithCentrality){
    parV2scaling[0] = 1. / (1 + parV2scaling[1]*avCentr + parV2scaling[2]*avCentr*avCentr); // normalize to average centrality

    scaleV2 = parV2scaling[0]*(1 + parV2scaling[1]*centrality + parV2scaling[2]*centrality*centrality); // apply the trand of v2 w.r.t. centrality
  }

  fgV2->SetParameter(2,fgV3Overv2);

  Float_t psi = gRandom->Rndm()*TMath::Pi();
  fgEventplane = psi;
  Float_t psi3 = gRandom->Rndm()*TMath::Pi()*2/3;
  Float_t psi4 = gRandom->Rndm()*TMath::Pi()*2/4;
  fgV2->SetParameter(1,psi);
  fgV2->SetParameter(3,psi3);
  fgV2->SetParameter(4,psi4);

  Int_t npart = 0;

  Float_t origin0[3];
  for (Int_t j=0;j<3;j++) origin0[j]=fOrigin[j];
  Float_t time;
  time = fTimeOrigin;
  if(fVertexSmear == kPerEvent) {
    Vertex();
    for (Int_t j=0; j < 3; j++) origin0[j] = fVertex[j];
    time = fTime;
  } // if kPerEvent

  printf("Generate event with centrality = %3.1f%c, |y|<%4.1f\n",centrality,'%',fYMaxValue);

  for(Int_t isp=0;isp < fgNspecies;isp++){
    if(! fgHSpectrum[isp]) continue;

    Int_t npartSp = Int_t(fgMult[isp]*2*fYMaxValue + gRandom->Rndm());

    printf("Total number of %i = %i\n",fgPdgInput[isp],npartSp);

    for(Int_t ipart =0; ipart < npartSp; ipart++){
      Int_t pdg = fgPdgInput[isp];

      Double_t y = gRandom->Rndm()*2*fYMaxValue - fYMaxValue;
      Double_t ytanh = TMath::TanH(y);
      Double_t pt = fgHSpectrum[isp]->GetRandom();
      Double_t mass = pdgD->GetParticle(pdg)->Mass();
      Double_t mt2 = pt*pt + mass*mass;
      Double_t pz = ytanh*TMath::Sqrt(mt2)/TMath::Sqrt(1-ytanh*ytanh);
      Double_t etot = TMath::Sqrt(mt2 + pz*pz);
      TLorentzVector tempVect(pt,0,pz,etot);
      //      Double_t eta = tempVect.PseudoRapidity();
      Double_t scaleEtaV2 = 1; // set the eta dependence
      if(TMath::Abs(y)> fYlimitForFlatness) scaleEtaV2 = 1 - fYdecreaseV2*(TMath::Abs(y) - fYlimitForFlatness);

      if(fgHv2[isp]) fgV2->SetParameter(0,fgHv2[isp]->Interpolate(pt) * scaleV2 * scaleEtaV2); 
      else fgV2->SetParameter(0,0.);
      Double_t phi = fgV2->GetRandom(-TMath::Pi(),TMath::Pi());
      Double_t px = pt*TMath::Cos(phi);
      Double_t py = pt*TMath::Sin(phi);
      Float_t p[3] = {static_cast<Float_t>(px),static_cast<Float_t>(py),static_cast<Float_t>(pz)};
      Float_t polar[3] = {0.,0.,0.};

      if(TMath::Abs(y)< fYlimitForFlatness || gRandom->Rndm() < 1 - fYdecreaseSp*(TMath::Abs(y) - fYlimitForFlatness)){// check on pseudorapidity distribution
        //	printf("%f %f\n",eta,phi - psi); // for debugging
        PushTrack(1, -1, pdg, p, origin0, polar, time, statusPdg[isp], npart, 1., 1);
        KeepTrack(npart);
        npart++;
      }
    }
  }

  TArrayF eventVertex;
  eventVertex.Set(3);
  eventVertex[0] = origin0[0];
  eventVertex[1] = origin0[1];
  eventVertex[2] = origin0[2];

  // Header
  AliGenEventHeaderTunedPbPb* header = new AliGenEventHeaderTunedPbPb("tunedOnPbPb");
  // Event Vertex
  header->SetPrimaryVertex(eventVertex);
  header->SetInteractionTime(time);
  header->SetNProduced(npart);
  header->SetCentrality(centrality);
  header->SetPsi2(psi);
  header->SetPsi3(psi3);
  header->SetPsi4(psi4);
  gAlice->SetGenEventHeader(header); 
}

void AliGenTunedOnPbPb::SetPtRange(Float_t ptmin, Float_t ptmax) {
  AliGenerator::SetPtRange(ptmin, ptmax);
}
//_____________________________________________________________________________
TH1F *AliGenTunedOnPbPb::GetMultVsCentrality(Int_t species){
  char title[100];
  snprintf(title,100,"pdg = %i;centrality;dN/dy",fgPdgInput[species]);
  TH1F *h = new TH1F("multVsCentr",title,100,0,100);

  for(Int_t i=1;i<=100;i++){
    Float_t x = i+0.5;
    SetParameters(x);
    h->SetBinContent(i,fgMult[species]);
  }

  return h;
}
//_____________________________________________________________________________
void AliGenTunedOnPbPb::SetParameters(Float_t centrality){

  if(!fgV2) fgV2 = new TF1("fv2Par","TMath::Max(0.,(1 + 2*[0]*cos(2*(x-[1])) + 2*[0]*[2]*cos(3*(x-[3])) + [0]*[2]*cos(4*(x-[4]))))",-TMath::Pi(),TMath::Pi()); // v4 is approx. 0.5*v3

  Float_t fr[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

  if(centrality < 7.5){
    fr[0] = (7.5 - centrality)/5.;
    fr[1] = (centrality-2.5)/5.;
  }
  else if(centrality < 15){
    fr[1] = (15-centrality)/7.5;
    fr[2] = (centrality-7.5)/7.5;
  }
  else if(centrality < 25){
    fr[2] = (25-centrality)/10.;
    fr[3] = (centrality-15)/10.;
  }
  else if(centrality < 35){
    fr[3] = (35-centrality)/10.;
    fr[4] = (centrality-25)/10.;
  }
  else if(centrality < 45){
    fr[4] = (45-centrality)/10.;
    fr[5] = (centrality-35)/10.;
  }
  else if(centrality < 55){
    fr[5] = (55-centrality)/10.;
    fr[6] = (centrality-45)/10.;
  }
  else if(centrality < 65){
    fr[6] = (65-centrality)/10.;
    fr[7] = (centrality-55)/10.;
  }
  else if(centrality < 75){
    fr[7] = (75-centrality)/10.;
    fr[8] = (centrality-65)/10.;
  }
  else{
    fr[8] = 1.0;
  }

  // parameters as a function of centrality
  Float_t multCent[9*fgNspecies] = {
    733.,733.,733.,109.,109.,34.0,34.0,109.,28.,28.,11.5,3.1  ,3.1  ,0.5  ,0.5  ,34.0,34.0,
    606.,606.,606.,91.0,91.0,28.0,28.0,91. ,24.,24.,9.  ,2.7  ,2.7  ,0.45 ,0.45 ,28.0,28.0,
    455.,455.,455.,68.0,68.0,21.0,21.0,68. ,20.,20.,8.  ,2.4  ,2.4  ,0.40 ,0.40 ,21.0,21.0,
    307.,307.,307.,46.0,46.0,14.5,14.5,46. ,14.,14.,5.5 ,1.5  ,1.5  ,0.2  ,0.2  ,14.5,14.5,
    201.,201.,201.,30.0,30.0,9.60,9.60,30. ,9. ,9. ,3.5 ,0.9  ,0.9  ,0.08 ,0.08 ,9.60,9.60,
    124.,124.,124.,18.3,18.3,6.10,6.10,18.3,5.1,5.1,2.2 ,0.6  ,0.6  ,0.055,0.055,6.10,6.10,
    71.0,71.0,71.0,10.2,10.2,3.60,3.60,10.2,2.6,2.6,1.4 ,0.36 ,0.36 ,0.035,0.035,3.60,3.60,
    37.0,37.0,37.0,5.10,5.10,2.00,2.00,5.10,1.5,1.5,0.50,0.020,0.020,0.015,0.015,2.00,2.00,
    17.0,17.0,17.0,2.30,2.30,0.90,0.90,2.30,0.6,0.6,0.16,0.006,0.006,0.005,0.005,0.90,0.90
  };

  Float_t v3Overv2Cent[9] = {1.2,0.82,0.625,0.5,0.45,0.4,0.37,0.3,0.3};

  fgV3Overv2 = 0;
  for(Int_t j=0;j < 9;j++)
    fgV3Overv2 += fr[j]*v3Overv2Cent[j];

  // set parameters for current centrality
  for(Int_t i=0;i < fgNspecies;i++){
    fgMult[i] = 0;

    for(Int_t j=0;j < 9;j++){
      fgMult[i] += fr[j]*multCent[i+j*fgNspecies];
    }
  }

  if(centrality > 80){
    for(Int_t i=0;i < fgNspecies;i++)
      fgMult[i] /= TMath::Log(centrality-77.);
  }
}
