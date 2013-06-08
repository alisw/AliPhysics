#include <iostream>

#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TF1.h>

#include "AliToyMCEvent.h"

#include "AliToyMCEventGeneratorSimple.h"

ClassImp(AliToyMCEventGeneratorSimple);


AliToyMCEventGeneratorSimple::AliToyMCEventGeneratorSimple()
  :AliToyMCEventGenerator()
  ,fVertexMean(0.)
  ,fVertexSigma(7.)
{
  //
  
}
//________________________________________________________________
AliToyMCEventGeneratorSimple::AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen)
  :AliToyMCEventGenerator(gen)
  ,fVertexMean(gen.fVertexMean)
  ,fVertexSigma(gen.fVertexSigma)
{
}
//________________________________________________________________
AliToyMCEventGeneratorSimple::~AliToyMCEventGeneratorSimple()
{
}
 //________________________________________________________________
AliToyMCEventGeneratorSimple& AliToyMCEventGeneratorSimple::operator = (const AliToyMCEventGeneratorSimple &gen)
{
  //assignment operator
  if (&gen == this) return *this;
  new (this) AliToyMCEventGeneratorSimple(gen);

  return *this;
}
//________________________________________________________________
void AliToyMCEventGeneratorSimple::SetParameters(Double_t vertexMean, Double_t vertexSigma) {
  fVertexMean = vertexMean;
  fVertexSigma = vertexSigma;
}
//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::Generate(Double_t time) {

  AliToyMCEvent *retEvent = new AliToyMCEvent();
  retEvent->SetT0(time);
  retEvent->SetX(0);
  retEvent->SetX(0);
  retEvent->SetZ(gRandom->Gaus(fVertexMean,fVertexSigma));

  Double_t etaCuts=.9;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  TF1 fpt("fpt",Form("x*(1+(sqrt(x*x+%f^2)-%f)/([0]*[1]))^(-[0])",mass,mass),0.4,10);
  fpt.SetParameters(7.24,0.120);
  fpt.SetNpx(10000);
  Int_t nTracks = 400; //TODO: draw from experim dist 
  
  for(Int_t iTrack=0; iTrack<nTracks; iTrack++){
    Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
    Double_t eta = gRandom->Uniform(-etaCuts, etaCuts);
    Double_t pt = fpt.GetRandom(); // momentum for f1
    Short_t sign=1;
    if(gRandom->Rndm() < 0.5){
      sign =1;
    }else{
      sign=-1;
    }
    
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
    Double_t pxyz[3];
    pxyz[0]=pt*TMath::Cos(phi);
    pxyz[1]=pt*TMath::Sin(phi);
    pxyz[2]=pt*TMath::Tan(theta);
    Double_t vertex[3]={0,0,retEvent->GetZ()};
    Double_t cv[21]={0};
    AliToyMCTrack *tempTrack = new AliToyMCTrack(vertex,pxyz,cv,sign);
    
    Bool_t trackDist = DistortTrack(*tempTrack, time);
    if(trackDist) retEvent->AddTrack(*tempTrack);
    delete tempTrack;
  }

  return retEvent;
}

//________________________________________________________________
void AliToyMCEventGeneratorSimple::RunSimulation(const Int_t nevents/*=10*/)
{
  //
  // run simple simulation with equal event spacing
  //

  if (!ConnectOutputFile()) return;

  Double_t eventTime=0.;
  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
  
  for (Int_t ievent=0; ievent<nevents; ++ievent){
    fEvent = Generate(eventTime);
    FillTree();
    delete fEvent;
    fEvent=0x0;
    eventTime+=eventSpacing;
  }

  CloseOutputFile();
}

