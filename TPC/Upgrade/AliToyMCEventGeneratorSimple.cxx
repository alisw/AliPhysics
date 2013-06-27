#include <iostream>

#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TF1.h>
#include <TStopwatch.h>
#include <AliESDtrackCuts.h>
#include <AliESDtrack.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

#include "AliToyMCEvent.h"

#include "AliToyMCEventGeneratorSimple.h"


ClassImp(AliToyMCEventGeneratorSimple);


AliToyMCEventGeneratorSimple::AliToyMCEventGeneratorSimple()
  :AliToyMCEventGenerator()
  ,fVertexMean(0.)
  ,fVertexSigma(7.)
  ,fNtracks(20)
  ,fInputFileNameESD("")
  ,fESDCuts(0x0)
{
  //
  
}
//________________________________________________________________
AliToyMCEventGeneratorSimple::AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen)
  :AliToyMCEventGenerator(gen)
  ,fVertexMean(gen.fVertexMean)
  ,fVertexSigma(gen.fVertexSigma)
  ,fNtracks(gen.fNtracks)
  ,fInputFileNameESD(gen.fInputFileNameESD)
  ,fESDCuts(gen.fESDCuts)
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
void AliToyMCEventGeneratorSimple::SetParametersSimple(Double_t vertexMean, Double_t vertexSigma) {
  fVertexMean = vertexMean;
  fVertexSigma = vertexSigma;
}
//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::Generate(Double_t time)
{
  //
  //
  //
  
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
//   Int_t nTracks = 400; //TODO: draw from experim dist
//   Int_t nTracks = 20; //TODO: draw from experim dist
  
  for(Int_t iTrack=0; iTrack<fNtracks; iTrack++){
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
void AliToyMCEventGeneratorSimple::RunSimulation(const Int_t nevents/*=10*/, const Int_t ntracks/*=400*/)
{
  //
  // run simple simulation with equal event spacing
  //

  if (!ConnectOutputFile()) return;
  //initialise the space charge. Should be done after the tree was set up
  InitSpaceCharge();
  
  fNtracks=ntracks;
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
  TStopwatch s;
  for (Int_t ievent=0; ievent<nevents; ++ievent){
    printf("Generating event %3d (%.3g)\n",ievent,eventTime);
    fEvent = Generate(eventTime);
    FillTree();
    delete fEvent;
    fEvent=0x0;
    eventTime+=eventSpacing;
  }
  s.Stop();
  s.Print();
  
  CloseOutputFile();
}

//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::GenerateESD(AliESDEvent &esdEvent, Double_t time) {

 
  AliToyMCEvent *retEvent = new AliToyMCEvent();
  retEvent->SetT0(time);
  retEvent->SetX(esdEvent.GetPrimaryVertex()->GetX());
  retEvent->SetY(esdEvent.GetPrimaryVertex()->GetY());
  retEvent->SetZ(esdEvent.GetPrimaryVertex()->GetZ());


  
  if(!fNtracks) fNtracks =  esdEvent.GetNumberOfTracks();
  Int_t nUsedTracks = 0;
  for(Int_t iTrack=0; iTrack<esdEvent.GetNumberOfTracks(); iTrack++){

    AliESDtrack *part = esdEvent.GetTrack(iTrack);
    if(!part) continue;
    if (!fESDCuts->AcceptTrack(/*(AliESDtrack*)*/part))continue;
    if(part->Pt() < 0.3) continue; //avoid tracks that fail to propagate through entire volume
    Double_t pxyz[3];
    pxyz[0]=part->Px();
    pxyz[1]=part->Py();
    pxyz[2]=part->Pz();
    Double_t vertex[3]={retEvent->GetX(),retEvent->GetY(),retEvent->GetZ()};
    
    Double_t cv[21]={0};
    Int_t sign = part->Charge();
    AliToyMCTrack *tempTrack = new AliToyMCTrack(vertex,pxyz,cv,sign);
    
    Bool_t trackDist = DistortTrack(*tempTrack, time);
    if(trackDist) {
      retEvent->AddTrack(*tempTrack);
     
      nUsedTracks++;
    }
    delete tempTrack;
    
    if(nUsedTracks >= fNtracks) break;
  }

  return retEvent;
}
//________________________________________________________________
void AliToyMCEventGeneratorSimple::RunSimulationESD(const Int_t nevents/*=10*/, const Int_t ntracks/*=400*/)
{
  //
  // run simulation using esd input with equal event spacing
  //

  if (!ConnectOutputFile()) return;
  TFile f(Form("%s%s",fInputFileNameESD.Data(),"#AliESDs.root"));
  if(!f.IsOpen()) {
    std::cout << "file "<<fInputFileNameESD.Data() <<" could not be opened" << std::endl;
    return;
  }
  TTree* esdTree = (TTree*) f.Get("esdTree");
  if(!esdTree) {
    std::cout <<"no esdTree in file " << std::endl;
    return;
  }
  InitSpaceCharge();
  AliESDEvent* esdEvent = new AliESDEvent();
  fNtracks = ntracks;
 
  esdEvent->ReadFromTree(esdTree);
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
  TStopwatch s;
  Int_t nEvents = nevents;
  fESDCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  if(nevents>esdTree->GetEntries()) nEvents = esdTree->GetEntries();
  Int_t nUsedEvents = 0;
  for (Int_t ievent=0; ievent<esdTree->GetEntries(); ++ievent){
    printf("Generating event %3d (%.3g)\n",nUsedEvents,eventTime);
    esdTree->GetEvent(ievent);
    if(esdEvent->GetNumberOfTracks()==0) {
      std::cout << " tracks == 0" << std::endl;
      continue;
    }
   
    fEvent = GenerateESD(*esdEvent, eventTime);
    if(fEvent->GetNumberOfTracks() >=10) {
      nUsedEvents++;
      FillTree();
      eventTime+=eventSpacing;
    }
    delete fEvent;
    fEvent=0x0;
    
    if(nUsedEvents>=nevents) break;
  }
  s.Stop();
  s.Print();
  
  CloseOutputFile();
}

//________________________________________________________________

void AliToyMCEventGeneratorSimple::RunSimulationBunchTrain(const Int_t nevents/*=10*/, const Int_t ntracks/*=400*/)
{
  //
  // run simple simulation with equal event spacing
  //

  //Parameters for bunc
  const Double_t abortGap = 3e-6; //
  const Double_t collFreq = 50e3;
  const Double_t bSpacing = 50e-9; //bunch spacing
  const Int_t nTrainBunches = 48;
  const Int_t nTrains = 12;
  const Double_t revFreq = 1.11e4; //revolution frequency
  const Double_t collProb = collFreq/(nTrainBunches*nTrains*revFreq);
  const Double_t trainLength = bSpacing*(nTrainBunches-1);
  const Double_t totTrainLength = nTrains*trainLength;
  const Double_t trainSpacing = (1./revFreq - abortGap - totTrainLength)/(nTrains-1); 
  Bool_t equalSpacing = kFALSE;

  TRandom3 *rand = new TRandom3();
  //rand->SetSeed();

  if (!ConnectOutputFile()) return;
  //initialise the space charge. Should be done after the tree was set up
  InitSpaceCharge();
  
  fNtracks=ntracks;
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
  TStopwatch s;
  Int_t nGeneratedEvents = 0;
  Int_t bunchCounter = 0;
  Int_t trainCounter = 0;

  //for (Int_t ievent=0; ievent<nevents; ++ievent){
  while (nGeneratedEvents<nevents){
    //  std::cout <<trainCounter << " " << bunchCounter << " "<< "eventTime " << eventTime << std::endl;
    
    if(equalSpacing)  {
      printf("Generating event %3d (%.3g)\n",nGeneratedEvents,eventTime);
      fEvent = Generate(eventTime);
      nGeneratedEvents++;
      FillTree();
      delete fEvent;
      fEvent=0x0;
      eventTime+=1./collFreq;
    }
    else{
      Int_t nCollsInCrossing = rand -> Poisson(collProb);
      for(Int_t iColl = 0; iColl<nCollsInCrossing; iColl++){
	printf("Generating event %3d (%.3g)\n",nGeneratedEvents,eventTime);
	fEvent = Generate(eventTime);
	nGeneratedEvents++;
	FillTree();
	delete fEvent;
	fEvent=0x0;

      }
      bunchCounter++;

      if(bunchCounter>=nTrainBunches){
	
	trainCounter++;
	if(trainCounter>=nTrains){
	  eventTime+=abortGap;
	  trainCounter=0;
	}
	else eventTime+=trainSpacing;

	bunchCounter=0;
      }
      else eventTime+= bSpacing;
      
    }


  }
  s.Stop();
  s.Print();
  delete rand;
  CloseOutputFile();
}
