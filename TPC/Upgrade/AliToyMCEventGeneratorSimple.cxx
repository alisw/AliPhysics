#include <iostream>

#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TF1.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

#include <AliESDtrackCuts.h>
#include <AliESDtrack.h>
#include <AliTPCLaserTrack.h>

#include "AliToyMCEvent.h"

#include "AliToyMCEventGeneratorSimple.h"

/*


*/

ClassImp(AliToyMCEventGeneratorSimple);


AliToyMCEventGeneratorSimple::AliToyMCEventGeneratorSimple()
  :AliToyMCEventGenerator()
  ,fVertexMean(0.)
  ,fVertexSigma(7.)
  ,fNtracks(20)
  ,fInputFileNameESD("")
  ,fESDCuts(0x0)
  ,fInputIndex(0)
  ,fESDEvent(0x0)
  ,fESDTree(0x0)
  ,fInputFile(0x0)
  ,fHPt(0x0)
  ,fHEta(0x0)
  ,fHMult(0x0)
  ,fHistosSet(kFALSE)
  ,fParamFile(0x0)
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
  ,fInputIndex(gen.fInputIndex)
  ,fESDEvent(gen.fESDEvent)
  ,fESDTree(gen.fESDTree)
  ,fInputFile(gen.fInputFile)
  ,fHPt(gen.fHPt)
  ,fHEta(gen.fHEta)
  ,fHMult(gen.fHMult)
  ,fHistosSet(gen.fHistosSet)
  ,fParamFile(0x0)
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
void AliToyMCEventGeneratorSimple::SetParametersToyGen(const Char_t* parfilename/*="files/params.root*/, Double_t vertexMean/*=0*/, Double_t vertexSigma/*=7.*/) {
  fVertexMean = vertexMean;
  fVertexSigma = vertexSigma;
//   fParamFile = new TFile(parfilename, "read");
  TFile f(parfilename);
  gROOT->cd();
  fHPt = (TH1F*) f.Get("hPt");
  fHEta = (TH1F*) f.Get("hEta");
  fHMult = (TH1I*) f.Get("hMult") ;
  fHistosSet = kTRUE;
  f.Close();
 
  
}
//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::Generate(Double_t time)
{
  //
  // Generate an event at 'time'
  //

  // iterate over space charge maps in case they are set
  IterateSC();
  
  AliToyMCEvent *retEvent = new AliToyMCEvent();
  retEvent->SetT0(time);
  retEvent->SetX(0);
  retEvent->SetX(0);
  Double_t vertexZ=0;
  //limit vertex to +- 10cm
  while ( TMath::Abs(vertexZ=gRandom->Gaus(fVertexMean,fVertexSigma))>10. );
  retEvent->SetZ(vertexZ);

  Double_t etaCuts=0.9;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  static TF1 fpt("fpt",Form("x*(1+(sqrt(x*x+%f^2)-%f)/([0]*[1]))^(-[0])",mass,mass),0.3,100);
  if (fpt.GetParameter(0)<1){
    printf("Set Parameters\n");
    fpt.SetParameters(7.24,0.120);
    fpt.SetNpx(200);
  }
//   Int_t nTracks = 400; //TODO: draw from experim dist
//   Int_t nTracks = 20; //TODO: draw from experim dist
   Int_t nTracksLocal = fNtracks;
  if(fHistosSet) {
    nTracksLocal = fHMult->GetRandom(); 
    if(nTracksLocal > fNtracks) nTracksLocal = fNtracks;
  }
  std::cout << "Generating " << nTracksLocal << " tracks " << std::endl;

  for(Int_t iTrack=0; iTrack<nTracksLocal; iTrack++){
    Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
    Double_t eta = 0.;
    Double_t pt = 0.;
    
    if(fHistosSet) {//draw from histograms if set
      eta = fHEta->GetRandom();
      pt = fHPt->GetRandom();
    }
    else {
      eta = gRandom->Uniform(-etaCuts, etaCuts);
      pt = fpt.GetRandom(); // momentum for f1
    }
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
    AliToyMCTrack *tempTrack = retEvent->AddTrack(vertex,pxyz,cv,sign);
    // use unique ID for track number
    // this will be used in DistortTrack track to set the cluster label
    // in one simulation the track id should be unique for performance studies
    tempTrack->SetUniqueID(fCurrentTrack++);
    DistortTrack(*tempTrack, time);
  }

  return retEvent;
}

//________________________________________________________________
void AliToyMCEventGeneratorSimple::RunSimulation(const Int_t nevents/*=10*/, const Int_t ntracks/*=400*/, const Int_t rate/*=50*/)
{
  //
  // run simple simulation with equal event spacing
  // rate in kHz
  //

  if (!ConnectOutputFile()) return;
  //initialise the space charge. Should be done after the tree was set up
  InitSpaceCharge();

  // number of tracks to simulate per interaction
  fNtracks=ntracks;
  // within one simulation the track count should be unique for effeciency studies
  // don't use the track ID 0 since the label -0 is not different from 0
  fCurrentTrack=1;
  
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./rate/1e3; //50kHz equally spaced
  TStopwatch s;
  for (Int_t ievent=0; ievent<nevents; ++ievent){
    printf("Generating event %3d (%.3g)\n",ievent,eventTime);
    fEvent = Generate(eventTime);
    SetSCScalingFactor();
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
void AliToyMCEventGeneratorSimple::RunSimulationLaser(const Int_t nevents/*=1*/)
{
  //
  // run simple simulation with equal event spacing
  //
  
  if (!ConnectOutputFile()) return;
  //initialise the space charge. Should be done after the tree was set up
  InitSpaceCharge();
  
  // within one simulation the track count should be unique for effeciency studies
  fCurrentTrack=1;
  
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./10.; //laser is running at 10Hz equally spaced
  TStopwatch s;
  for (Int_t ievent=0; ievent<nevents; ++ievent){
    printf("Generating event %3d (%.3g)\n",ievent,eventTime);
    fEvent = GenerateLaser(eventTime);
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
    AliToyMCTrack *tempTrack = retEvent->AddTrack(vertex,pxyz,cv,sign);
    // use unique ID for track number
    // this will be used in DistortTrack track to set the cluster label
    // in one simulation the track id should be unique for performance studies
    tempTrack->SetUniqueID(fCurrentTrack++);
    DistortTrack(*tempTrack, time);
    nUsedTracks++;
    
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
 
 
  esdEvent->ReadFromTree(esdTree);

  fNtracks = ntracks;
  // within one simulation the track count should be unique for effeciency studies
  fCurrentTrack=1;
  
  Double_t eventTime=0.;
  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
  TStopwatch s;
  
  fESDCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  gRandom->SetSeed();
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
  // within one simulation the track count should be unique for effeciency studies
  fCurrentTrack=1;
  
  Double_t eventTime=0.;
  //  const Double_t eventSpacing=1./50e3; //50kHz equally spaced
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
      SetSCScalingFactor();
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
        SetSCScalingFactor();
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





//________________________________________________________________
Int_t AliToyMCEventGeneratorSimple::OpenInputAndGetMaxEvents(const Int_t type, const Int_t nevents) {

  

  if(type==0) return nevents;
  
  if(type==1) {
    
    fInputFile = new TFile(Form("%s%s",fInputFileNameESD.Data(),"#AliESDs.root"),"read");
    if(!fInputFile->IsOpen()) {
      std::cout << "file "<<fInputFileNameESD.Data() <<" could not be opened" << std::endl;
      return 0;
    }
    fESDTree = (TTree*) fInputFile->Get("esdTree");
    if(!fESDTree) {
      std::cout <<"no esdTree in file " << std::endl;
      return 0;
    }
    fESDEvent = new AliESDEvent();
    fESDEvent->ReadFromTree(fESDTree);
    fESDCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);

    fInputIndex = 0;

    gRandom->SetSeed();
    return fESDTree->GetEntries();
   }

 
  std::cout << " no available input type (toymc, esd) specified" << std::endl;
  return 0;
 

}


//________________________________________________________________
void AliToyMCEventGeneratorSimple::RunSimulation2(const Bool_t equalspacing, const Int_t type, const Int_t nevents, const Int_t ntracks) {

  //type==0 simple toy
  //type==1 esd input

  Int_t nMaxEvents = OpenInputAndGetMaxEvents(type,nevents);
  if (!ConnectOutputFile()) return;
  
  InitSpaceCharge();


  fNtracks = ntracks;
  // within one simulation the track count should be unique for effeciency studies
  fCurrentTrack=1;
  
  Double_t eventTime=0.;
  TStopwatch s;
  // const Double_t eventSpacing=1./50e3; //50kHz equally spaced

  Int_t nGeneratedEvents = 0;
  
  for (Int_t ievent=0; ievent<nMaxEvents; ++ievent){
    printf("Generating event %3d (%.3g)\n",ievent,eventTime);
    
    Int_t nColls = 0;
    Double_t spacing = 0;
    GetNGeneratedEventsAndSpacing(equalspacing, nColls, spacing);

    for(Int_t iColl = 0; iColl<nColls; iColl++){
      if(type==0)    fEvent = Generate(eventTime);
      else if(type==1) fEvent = GenerateESD2(eventTime);
      if(fEvent) {
	FillTree();
	delete fEvent;
	fEvent=0x0;
	nGeneratedEvents++;
      }
    }
    eventTime+=spacing;
    if(nGeneratedEvents >= nevents) break;
  }
  s.Stop();
  s.Print();
  
  CloseOutputFile();
  CloseInputFile();


}
//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::GenerateESD2(Double_t time) {

  //test that enough tracks will pass cuts (and that there is tracks at all)
  Bool_t testEvent = kTRUE;

  // iterate over space charge maps in case they are set
  IterateSC();
  while(fESDTree->GetEvent(fInputIndex) && testEvent) {  
    
    Int_t nPassedCuts = 0;
    for(Int_t iTrack = 0; iTrack<fESDEvent->GetNumberOfTracks(); iTrack++){
      if(fESDCuts->AcceptTrack(fESDEvent->GetTrack(iTrack)) && (fESDEvent->GetTrack(iTrack)->Pt() < 0.3) ) nPassedCuts++;
    }
    if(nPassedCuts<fNtracks || (fNtracks>100 && nPassedCuts <100) ) fInputIndex++; //100 is a perhaps bit arbitrary, do we want the events were only a small number of tracks are accepted?
    else testEvent = kFALSE;
  }

  if(fInputIndex>=fESDTree->GetEntries()) return 0;
  
  

  AliToyMCEvent *retEvent = new AliToyMCEvent();
  retEvent->SetT0(time);
  retEvent->SetX(fESDEvent->GetPrimaryVertex()->GetX());
  retEvent->SetY(fESDEvent->GetPrimaryVertex()->GetY());
  retEvent->SetZ(fESDEvent->GetPrimaryVertex()->GetZ());


  
  if(!fNtracks) fNtracks =  fESDEvent->GetNumberOfTracks();
  Int_t nUsedTracks = 0;
  for(Int_t iTrack=0; iTrack<fESDEvent->GetNumberOfTracks(); iTrack++){

    AliESDtrack *part = fESDEvent->GetTrack(iTrack);
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

    AliToyMCTrack *tempTrack = retEvent->AddTrack(vertex,pxyz,cv,sign);
    // use unique ID for track number
    // this will be used in DistortTrack track to set the cluster label
    // in one simulation the track id should be unique for performance studies
    tempTrack->SetUniqueID(fCurrentTrack++);
    DistortTrack(*tempTrack, time);
    
    
    if(nUsedTracks >= fNtracks) break;
  }
  fInputIndex++;
 
  return retEvent;
}

//________________________________________________________________
AliToyMCEvent* AliToyMCEventGeneratorSimple::GenerateLaser(Double_t time)
{
  //
  // Generate an Event with laser tracks
  //

  // iterate over space charge maps in case they are set
  IterateSC();
  
  AliToyMCEvent *retEvent = new AliToyMCEvent();
  retEvent->SetEventType(AliToyMCEvent::kLaser);
  
  retEvent->SetT0(time);
  retEvent->SetX(0);
  retEvent->SetX(0);
  retEvent->SetZ(0);
  
  AliTPCLaserTrack::LoadTracks();
  TObjArray *arr=AliTPCLaserTrack::GetTracks();

  //since we have a laser track force no material budges
  Bool_t materialBudget=GetUseMaterialBudget();
  SetUseMaterialBudget(kFALSE);

  //the laser tracks have partially large inclination angles over the pads
  // -> relax the propagation contraint and switch off error messages
  SetIsLaser(kTRUE);
  Int_t debug=AliLog::GetDebugLevel("","AliToyMCEventGeneratorSimple");
  AliLog::SetClassDebugLevel("AliToyMCEventGeneratorSimple",-3);
  
  for (Int_t iTrack=0; iTrack<arr->GetEntriesFast(); ++iTrack){
    AliExternalTrackParam *track=(AliExternalTrackParam*)arr->At(iTrack);
    AliToyMCTrack *tempTrack = retEvent->AddTrack(AliToyMCTrack(*track));
    // for laser only TPC clusters exist
    tempTrack->SetUniqueID(fCurrentTrack++);
    MakeTPCClusters(*tempTrack, time);
  }

  SetIsLaser(kFALSE);
  AliLog::SetClassDebugLevel("AliToyMCEventGeneratorSimple",debug);
  SetUseMaterialBudget(materialBudget);
  return retEvent;
}

//________________________________________________________________
void AliToyMCEventGeneratorSimple::GetNGeneratedEventsAndSpacing(const Bool_t equalSpacing, Int_t &ngen, Double_t &spacing)
{

  static Int_t bunchCounter = 0;
  static Int_t trainCounter = 0;

  if(equalSpacing) {
    ngen =1;
    spacing = 1./50e3; //50kHz equally spaced
    return;
  }
  else if(!equalSpacing) {
      //parameters for bunch train
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
    Double_t time = 0;
    Int_t nCollsInCrossing = 0;
    while(nCollsInCrossing ==0){
      
      
    bunchCounter++;
    
    if(bunchCounter>=nTrainBunches){
      
      trainCounter++;
      if(trainCounter>=nTrains){
	time+=abortGap;
	trainCounter=0;
      }
      else time+=trainSpacing;

      bunchCounter=0;
    }
    else time+= bSpacing;


    nCollsInCrossing = gRandom -> Poisson(collProb);
    //std::cout << " nCollsInCrossing " << nCollsInCrossing <<  std::endl;
    }
    ngen = nCollsInCrossing;
    if(nCollsInCrossing > 1)std::cout << " nCollsInCrossing " << nCollsInCrossing <<  std::endl;
    spacing = time;
    return;

  }

}

//________________________________________________________________
Bool_t AliToyMCEventGeneratorSimple::CloseInputFile()
{
  //
  // close the input file
  //
  if (!fInputFile) return kFALSE;
  fInputFile->Close();
  delete fInputFile;
  fInputFile=0x0;

  return kTRUE;
}
