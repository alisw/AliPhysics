/*
***********************************************************
  Implementation of the AliMixingHandler class
  Contact: iarsene@cern.ch
  2015/08/07
  *********************************************************
*/

#ifndef ALIMIXINGHANDLER_H
#include "AliMixingHandler.h"
#endif

#include <iostream>
using std::cout;
using std::endl;
using std::flush;

#include <TMath.h>
#include <TTimeStamp.h>
#include <TRandom.h>

#include "AliReducedVarManager.h"
#include "AliReducedBaseTrack.h"

ClassImp(AliMixingHandler);

//_________________________________________________________________________
AliMixingHandler::AliMixingHandler() :
  TNamed(),
  fPoolDepth(100),
  fMixingThreshold(1.0),
  fDownscaleEvents(1.0),
  fDownscaleTracks(1.0),
  fPoolsLeg1("TClonesArray"),
  fPoolsLeg2("TClonesArray"),
  fNParallelCuts(0),
  fHistClassNames(""),
  fPoolSize(),
  fIsInitialized(kFALSE),
  fMixLikeSign(kTRUE),
  fCentralityLimits(),
  fEventVertexLimits(),
  fEventPlaneLimits(),
  fCentralityVariable(AliReducedVarManager::kNothing),
  fEventVertexVariable(AliReducedVarManager::kNothing),
  fEventPlaneVariable(AliReducedVarManager::kNothing),
  fHistos(0x0)
{
  // 
  // default constructor
  //
  Float_t dummyCentRange[2] = {0.0,100.};
  Float_t dummyZRange[2] = {-10.0,10.};
  Float_t dummyEPRange[2] = {(Float_t)(-0.5*TMath::Pi()),(Float_t)(0.5*TMath::Pi())};
  fCentralityLimits.Set(2,dummyCentRange);
  fEventVertexLimits.Set(2,dummyZRange);
  fEventPlaneLimits.Set(2,dummyEPRange);
  fPoolSize.Set(1);
}


//_________________________________________________________________________
AliMixingHandler::AliMixingHandler(const Char_t* name, const Char_t* title) :
  TNamed(name,title),
  fPoolDepth(100),
  fMixingThreshold(1.0),
  fDownscaleEvents(1.0),
  fDownscaleTracks(1.0),
  fPoolsLeg1("TClonesArray"),
  fPoolsLeg2("TClonesArray"),
  fNParallelCuts(0),
  fHistClassNames(""),
  fPoolSize(),
  fIsInitialized(kFALSE),
  fMixLikeSign(kTRUE),
  fCentralityLimits(),
  fEventVertexLimits(),
  fEventPlaneLimits(),
  fCentralityVariable(AliReducedVarManager::kNothing),
  fEventVertexVariable(AliReducedVarManager::kNothing),
  fEventPlaneVariable(AliReducedVarManager::kNothing),
  fHistos(0x0)
{
  //
  // Named constructor
  //
  Float_t dummyCentRange[2] = {0.0,100.};
  Float_t dummyZRange[2] = {-10.0,10.};
  Float_t dummyEPRange[2] = {(Float_t)(-0.5*TMath::Pi()),(Float_t)(0.5*TMath::Pi())};
  fCentralityLimits.Set(2,dummyCentRange);
  fEventVertexLimits.Set(2,dummyZRange);
  fEventPlaneLimits.Set(2,dummyEPRange);
  fPoolSize.Set(1);
}


//_________________________________________________________________________
AliMixingHandler::~AliMixingHandler() {
  //
  // destructor
  //
}


//_________________________________________________________________________
void AliMixingHandler::SetEventVariables(AliReducedVarManager::Variables cent, AliReducedVarManager::Variables vtx, 
                                         AliReducedVarManager::Variables ep) {
   //
   // set the event mixing variables
   //
   fCentralityVariable=cent; fEventVertexVariable=vtx; fEventPlaneVariable=ep;    
   AliReducedVarManager::SetUseVariable(fCentralityVariable);
   AliReducedVarManager::SetUseVariable(fEventVertexVariable);
   AliReducedVarManager::SetUseVariable(fEventPlaneVariable);
}


//_________________________________________________________________________
void AliMixingHandler::Init() {
  //
  // Initialization of pools
  // NOTE: The master array is a 1D array but it will be represented as a 3D array
  //       in centrality, vtxz and ep.
  //       The size of the array will be n_cent x n_z x n_ep
  //       The correct event category will be retrieved using the function FindEventCategory()
  //
  if(!fHistos) {
    cout << "AliMixingHandler::Init(): ERROR No histogram manager provided!" << endl;
    return;
  }
  if(fHistClassNames[0]=='\0') {
    cout << "AliMixingHandler::Init(): ERROR No names for the histogram classes provided!" << endl;
    return;
  }
  TObjArray* histClassArr = fHistClassNames.Tokenize(";");
  if(histClassArr->GetEntries()!=3*fNParallelCuts) {       // 3 because there is one class of histograms for each pair type: ++,+- and --
    cout << "AliMixingHandler::Init(): ERROR The number of cuts and the number of hist class names provided do not match!" << endl;
    cout << "                   hist classes: " << histClassArr->GetEntries() << ";    n-parallel cuts: " << fNParallelCuts << endl;
    return;
  }
  Int_t size = (fCentralityLimits.GetSize()-1)*(fEventVertexLimits.GetSize()-1)*(fEventPlaneLimits.GetSize()-1);
  fPoolsLeg1.Expand(size); fPoolsLeg1.SetOwner(kTRUE);
  fPoolsLeg2.Expand(size); fPoolsLeg2.SetOwner(kTRUE);
  
  fPoolSize.Set(fNParallelCuts*size);
  for(Int_t i=0;i<fNParallelCuts*size;++i) fPoolSize[i] = 0;
    
  // Initialize the random number generator for event/track downscaling
  TTimeStamp time;
  gRandom->SetSeed(time.GetNanoSec());
  
  fIsInitialized = kTRUE;
}


//_________________________________________________________________________
Bool_t AliMixingHandler::AcceptTrack() {
  //
  // Randomly decide to reject a track if a track downscale parameter has been set
  //
  if(fDownscaleTracks>1.0 && (gRandom->Rndm()>(1.0/fDownscaleTracks))) return kFALSE;
  return kTRUE;
}


//_________________________________________________________________________
void AliMixingHandler::FillEvent(TList* leg1List, TList* leg2List, Float_t* values, Int_t type) {
  //
  // Fill the leg1 and leg2 lists in the appropriate category, based on the event 
  // characteristics (centrality, vtxz, ep)
  //
  if(!fIsInitialized) Init();
  if(leg1List->GetEntries()==0 && leg2List->GetEntries()==0) return;
  
  // randomly accept/reject this event in case fDownscaleEvents is used
  if(fDownscaleEvents>1.0 && (gRandom->Rndm()>(1.0/fDownscaleEvents))) 
    return;
  
  // find the event category
  Int_t category = FindEventCategory(values[fCentralityVariable], values[fEventVertexVariable], values[fEventPlaneVariable]);
  if(category<0) return;   // event characteristics outside the defined ranges
  
  TClonesArray *leg1PoolP = static_cast<TClonesArray*>(fPoolsLeg1.At(category));
  if(!leg1PoolP) leg1PoolP = new(fPoolsLeg1[category]) TClonesArray("TList",1);
  leg1PoolP->SetOwner(kTRUE);
  TClonesArray *leg2PoolP=static_cast<TClonesArray*>(fPoolsLeg2.At(category));
  if(!leg2PoolP) leg2PoolP = new(fPoolsLeg2[category]) TClonesArray("TList",1);
  leg2PoolP->SetOwner(kTRUE);
  
  TClonesArray &leg1Pool=*leg1PoolP;
  TClonesArray &leg2Pool=*leg2PoolP;
  
  // add the leg lists to the appropriate pools
  TList *list1 = new(leg1Pool[leg1Pool.GetEntries()]) TList();
  TList *list2 = new(leg2Pool[leg2Pool.GetEntries()]) TList();
  list1->SetOwner(kTRUE); list2->SetOwner(kTRUE);
  for(Int_t it=0; it<leg1List->GetEntries(); ++it) 
    list1->Add(leg1List->At(it)->Clone());
  for(Int_t it=0; it<leg2List->GetEntries(); ++it) 
    list2->Add(leg2List->At(it)->Clone());
    
  // increment the size of the pools in this category
  ULong_t mixingMask = IncrementPoolSizes(leg1List,leg2List,category);
  
  // if full pool(s) were found then run the event mixing
  if(mixingMask) {
    RunEventMixing(leg1PoolP,leg2PoolP,mixingMask,type,values);
    ResetPoolSizes(mixingMask,category);
  }
}


//_________________________________________________________________________
Int_t AliMixingHandler::FindEventCategory(Float_t centrality, Float_t vtxz, Float_t ep) {
  //
  // Find the event category corresponding to the centrality, vtxz and ep values 
  //
  if(!fIsInitialized) Init();
  
  Int_t centBin = TMath::BinarySearch(fCentralityLimits.GetSize(), fCentralityLimits.GetArray(), centrality);
  if(centBin==-1 || centBin==fCentralityLimits.GetSize()-1) return -1;
  Int_t zBin = TMath::BinarySearch(fEventVertexLimits.GetSize(), fEventVertexLimits.GetArray(), vtxz);
  if(zBin==-1 || zBin==fEventVertexLimits.GetSize()-1) return -1;
  Int_t epBin = TMath::BinarySearch(fEventPlaneLimits.GetSize(), fEventPlaneLimits.GetArray(), ep);
  if(epBin==-1 || epBin==fEventPlaneLimits.GetSize()-1) return -1;
  return centBin*(fEventVertexLimits.GetSize()-1)*(fEventPlaneLimits.GetSize()-1) + 
         zBin*(fEventPlaneLimits.GetSize()-1) + epBin; 
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetEventPlaneBin(Int_t category) {
  //
  // Find the bin center of the centrality interval for a given category
  //
  return category%(fEventPlaneLimits.GetSize()-1);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetEventVertexBin(Int_t category) {
  //
  // Find the bin center of the centrality interval for a given category
  //
  Int_t zBin = category-GetEventPlaneBin(category);
  zBin /= (fEventPlaneLimits.GetSize()-1);
  return zBin%(fEventVertexLimits.GetSize()-1);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetCentralityBin(Int_t category) {
  //
  // Find the bin center of the centrality interval for a given category
  //
  Int_t centBin = category-GetEventPlaneBin(category);
  centBin -= GetEventVertexBin(category)*(fEventPlaneLimits.GetSize()-1);
  centBin /= (fEventPlaneLimits.GetSize()-1)*(fEventVertexLimits.GetSize()-1);
  return centBin%(fCentralityLimits.GetSize()-1);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetPoolSize(Int_t cutNumber, Float_t centrality, Float_t vtxz, Float_t ep) {
  //
  // Get the pool size for a given set of (cut,centrality,z,ep)
  //
  if(cutNumber<0 || cutNumber>fNParallelCuts) return -1;
  Int_t eventCategory = FindEventCategory(centrality,vtxz,ep);
  return GetPoolSize(cutNumber, eventCategory);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetPoolSize(Int_t cutNumber, Int_t eventCategory) {
  //
  // Get the pool size for a given set of (cut,centrality,z,ep)
  //
  if(cutNumber<0 || cutNumber>fNParallelCuts) return -1;
  if(eventCategory<0) return -1;
  
  Int_t pool = cutNumber*(fCentralityLimits.GetSize()-1)
                        *(fEventVertexLimits.GetSize()-1)
		        *(fEventPlaneLimits.GetSize()-1) 
	       + eventCategory;
  return fPoolSize[pool];
}


//_________________________________________________________________________
void AliMixingHandler::ResetPoolSizes(ULong_t mixingMask, Int_t category) {
  //
  // Reset pool sizes for a given event category and mask
  //
  Int_t nCategories = (fCentralityLimits.GetSize()-1)*(fEventVertexLimits.GetSize()-1)*(fEventPlaneLimits.GetSize()-1);
  for(Int_t icut=0; icut<fNParallelCuts; ++icut) {
    if(mixingMask&(ULong_t(1)<<icut))
	fPoolSize[icut*nCategories+category] = 0;
  }
}


//_________________________________________________________________________
ULong_t AliMixingHandler::IncrementPoolSizes(TList* list1, TList* list2, Int_t eventCategory) {
  //
  // Check which cut bits are on and increment the pool sizes accordingly
  //
  AliReducedBaseTrack* track=0x0;
  ULong_t cutsMask=0;
  
  // Check which cuts are fulfilled by the tracks in these lists
  TIter trackIter1(list1);
  for(Int_t i=0; i<list1->GetEntries();++i) {
    track=(AliReducedBaseTrack*)trackIter1();
    for(UShort_t icut=0;icut<fNParallelCuts;++icut)
      if(track->TestFlag(icut)) cutsMask |= (ULong_t(1)<<icut);
  }
  TIter trackIter2(list2);
  for(Int_t i=0; i<list2->GetEntries();++i) {
    track=(AliReducedBaseTrack*)trackIter2();
    for(UShort_t icut=0;icut<fNParallelCuts;++icut)
      if(track->TestFlag(icut)) cutsMask |= (ULong_t(1)<<icut);
  }
    
  // increment the pools for those cuts which got at least one track
  Int_t nCategories = (fCentralityLimits.GetSize()-1)*(fEventVertexLimits.GetSize()-1)*(fEventPlaneLimits.GetSize()-1);
  Bool_t fullPoolFound = kFALSE;
  for(Int_t icut=0;icut<fNParallelCuts;++icut) {
    if(cutsMask & (ULong_t(1)<<icut))
      fPoolSize[icut*nCategories+eventCategory] += 1;
    if(fPoolSize[icut*nCategories+eventCategory]==fPoolDepth) 
      fullPoolFound = kTRUE;          
  }
  
  // If a completely filled pool is found, then look for the other cuts in this event category
  // to see if any of them is above the mixing threshold (fMixingThreshold)
  ULong_t mixingMask = 0;
  for(Int_t icut=0;icut<fNParallelCuts;++icut) {
    if(fPoolSize[icut*nCategories+eventCategory]>=Int_t(fMixingThreshold*fPoolDepth)) 
      mixingMask |= (ULong_t(1)<<icut);
  }
  return mixingMask;
}


//_________________________________________________________________________
void AliMixingHandler::RunLeftoverMixing(Int_t type) {
  //
  // Run event mixing over all event categories
  //
  cout << "========================================================================" << endl;
  cout << "                            Leftover mixing " << endl;
  cout << "========================================================================" << endl;
  
  // create a mixing mask which enables all cuts
  ULong_t mixingMask = 0;
  for(Int_t i=0; i<fNParallelCuts; ++i) mixingMask |= (ULong_t(1)<<i);
  Float_t values[AliReducedVarManager::kNVars];
  
  for(Int_t icateg=0; icateg<fPoolsLeg1.GetEntries(); ++icateg) {
    TClonesArray *leg1Pool = static_cast<TClonesArray*>(fPoolsLeg1.At(icateg));
    TClonesArray *leg2Pool = static_cast<TClonesArray*>(fPoolsLeg2.At(icateg));
    Int_t centBin = GetCentralityBin(icateg);
    Int_t zBin = GetEventVertexBin(icateg);
    Int_t epBin = GetEventPlaneBin(icateg);
    //cout << "epBin/low/high :: " << epBin << "/" << fEventPlaneLimits[epBin] << "/" << fEventPlaneLimits[epBin+1] << endl;
    values[fCentralityVariable] = 0.5*(fCentralityLimits[centBin]+fCentralityLimits[centBin+1]);
    values[fEventVertexVariable] = 0.5*(fEventVertexLimits[zBin]+fEventVertexLimits[zBin+1]);
    values[fEventPlaneVariable] = 0.5*(fEventPlaneLimits[epBin]+fEventPlaneLimits[epBin+1]);
    RunEventMixing(leg1Pool,leg2Pool,mixingMask,type,values);
    ResetPoolSizes(mixingMask,icateg);
  }  // end loop over categories
}


//_________________________________________________________________________
void AliMixingHandler::RunEventMixing(TClonesArray* leg1Pool, TClonesArray* leg2Pool, ULong_t mixingMask,
				      Int_t type, Float_t* values) {
  //
  // Run event mixing
  // NOTE: The mixingMask is a bit map with bits toggled for the pools which need mixing
  //       The type is the pair candidate type. It is used in AliReducedPairInfo::CandidateType, mainly to know which mass assumption to be made for the legs
  //
  cout << "AliMixingHandler::RunEventMixing for mask " << flush;
  AliReducedVarManager::PrintBits(mixingMask,fNParallelCuts);
  cout << ";  (cent/vtx/ep): " << values[fCentralityVariable] << "/"
       << values[fEventVertexVariable] << "/" << values[fEventPlaneVariable] << endl;
  
  Int_t entries = leg1Pool->GetEntries();
  if(entries<2) return;
  
  TObjArray* histClassArr = fHistClassNames.Tokenize(";");
  
  TIter iterEv1Leg1Pool(leg1Pool);
  TIter iterEv1Leg2Pool(leg2Pool);
  ULong_t testFlags1 = 0;
  ULong_t testFlags2 = 0;
  for(Int_t iev1=0; iev1<entries; ++iev1) {                            // first event loop
    // get the lists of leg1 and leg2 tracks for the first event
    TList* ev1Leg1List = (TList*)iterEv1Leg1Pool();
    TList* ev1Leg2List = (TList*)iterEv1Leg2Pool();
    //cout << "## iev1 " << iev1 << endl;
    
    TIter iterEv2Leg1Pool(leg1Pool);
    TIter iterEv2Leg2Pool(leg2Pool);
    for(Int_t iev2=0; iev2<entries; ++iev2) {                         // second event loop 
      TList* ev2Leg1List = (TList*)iterEv2Leg1Pool();
      TList* ev2Leg2List = (TList*)iterEv2Leg2Pool();
      if(iev1==iev2) continue;
      //cout << "#### iev2 " << iev2 << endl;
      
      //loop over the ev1-leg1 list
      TIter iterLeg1(ev1Leg1List);
      AliReducedBaseTrack* ev1Leg1=0x0;
      while((ev1Leg1=(AliReducedBaseTrack*)iterLeg1())) {
	// check that this track has at least one common bit with the mixing mask
	testFlags1 = mixingMask & ev1Leg1->GetFlags();
        //cout << "###### ev1Leg1 (p/px/testFlags1): " << ev1Leg1->P() << "/" << ev1Leg1->Px() << "/";
        //AliReducedVarManager::PrintBits(testFlags1, fNParallelCuts); cout << endl;
	if(!testFlags1) continue;
	
	//loop over the ev2-leg2 list 
	TIter iterLeg2(ev2Leg2List);
	AliReducedBaseTrack* ev2Leg2=0x0;
	while((ev2Leg2=(AliReducedBaseTrack*)iterLeg2())) {
	  // check that this track has at least one common bit with the mixing mask and with ev1-leg1
	  testFlags2 = testFlags1 & ev2Leg2->GetFlags();
          //cout << "###### ev2Leg2 (p/px/testFlags1): " << ev2Leg2->P() << "/" << ev2Leg2->Px() << "/";
          //AliReducedVarManager::PrintBits(testFlags2, fNParallelCuts); cout << endl;
	  if(!testFlags2) continue;
	  
	  // fill cross-pairs (leg1 - leg2) for the enabled bits
	  AliReducedVarManager::FillPairInfoME(ev1Leg1, ev2Leg2, type, values);
          //cout << "######## cross-pair (mass): " << values[AliReducedVarManager::kMass] << endl;
	  for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
            if((testFlags2)&(ULong_t(1)<<ibit)) 
              fHistos->FillHistClass(histClassArr->At(ibit*3+1)->GetName(), values);
          }  
	}  // end loop over the ev2-leg2 list
	
	if(!fMixLikeSign) continue;
	// loop over the ev2-leg1 list
	TIter iterLeg1LS(ev2Leg1List);
	AliReducedBaseTrack* ev2Leg1=0x0;
	while((ev2Leg1=(AliReducedBaseTrack*)iterLeg1LS())) {
	  // check that this track has at least one common bit with the mixing mask and with ev1-leg1
	  testFlags2 = testFlags1 & ev2Leg1->GetFlags();
          //cout << "###### ev2Leg1 (p/px/testFlags2): " << ev2Leg1->P() << "/" << ev2Leg1->Px() << "/";
          //AliReducedVarManager::PrintBits(testFlags2, fNParallelCuts); cout << endl;
	  if(!testFlags2) continue;
	  
	  // fill like-pairs (leg1 - leg1) for the enabled bits
	  AliReducedVarManager::FillPairInfoME(ev1Leg1, ev2Leg1, type, values);
          //cout << "######## like-pair leg1-leg1 (mass): " << values[AliReducedVarManager::kMass] << endl;
	  for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
            if((testFlags2)&(ULong_t(1)<<ibit)) 
              fHistos->FillHistClass(histClassArr->At(ibit*3+0)->GetName(), values);
          }  
	}  // end loop over the ev2-leg1 list
      }  // end loop over the ev1-leg1 list
      
      if(!fMixLikeSign) continue;
      //loop over the ev1-leg2 list
      TIter iterLeg2(ev1Leg2List);
      AliReducedBaseTrack* ev1Leg2=0x0;
      while((ev1Leg2=(AliReducedBaseTrack*)iterLeg2())) {
	// check that this track has at least one common bit with the mixing mask
	testFlags1 = mixingMask & ev1Leg2->GetFlags();
        //cout << "###### ev1Leg2 (p/px/testFlags1): " << ev1Leg2->P() << "/" << ev1Leg2->Px() << "/";
        //AliReducedVarManager::PrintBits(testFlags1, fNParallelCuts); cout << endl;
	if(!testFlags1) continue;
	
	//loop over the ev2-leg2 list 
	TIter iterLeg2LS(ev2Leg2List);
	AliReducedBaseTrack* ev2Leg2=0x0;
	while((ev2Leg2=(AliReducedBaseTrack*)iterLeg2LS())) {
	  // check that this track has at least one common bit with the mixing mask and with ev1-leg1
	  testFlags2 = testFlags1 & ev2Leg2->GetFlags();
          //cout << "###### ev2Leg2 (p/px/testFlags2): " << ev2Leg2->P() << "/" << ev2Leg2->Px() << "/";
          //AliReducedVarManager::PrintBits(testFlags2, fNParallelCuts); cout << endl;
	  if(!testFlags2) continue;
	  
	  // fill like-pairs (leg2 - leg2) for the enabled bits
	  AliReducedVarManager::FillPairInfoME(ev1Leg2, ev2Leg2, type, values);
          //cout << "######## like-pair leg2-leg2 (mass): " << values[AliReducedVarManager::kMass] << endl;
	  for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
            if((testFlags2)&(ULong_t(1)<<ibit)) 
              fHistos->FillHistClass(histClassArr->At(ibit*3+2)->GetName(), values);
          }  
	}  // end loop over the ev2-leg2 list
      }  // end loop over the ev1-leg2 list
    }  // end second event loop
  }  // end first event loop
  
  // unset the mixing flags --------------------------------------
  iterEv1Leg1Pool.Reset(); 
  iterEv1Leg2Pool.Reset();
  //cout << "Unsetting bits for which mixing was performed" << endl;
  for(Int_t ie1=0; ie1<entries; ++ie1) {
    TList* leg1List = (TList*)iterEv1Leg1Pool();        
    TIter iterLeg1(leg1List);
    AliReducedBaseTrack* track=0x0;
    while((track=(AliReducedBaseTrack*)iterLeg1())) {
      //cout << "leg1 p " << track->P() << "     " << flush; 
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << "/" << flush;
      testFlags1 = mixingMask & track->GetFlags();
      if(!testFlags1) {
	//AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << endl;
	continue;
      }
      for(UShort_t ibit=0; ibit<fNParallelCuts; ++ibit) {
	if((testFlags1)&(ULong_t(1)<<ibit)) 
	  track->UnsetFlag(ibit);
      }
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << endl;
    }  // end loop over the first leg list
    
    TList* leg2List = (TList*)iterEv1Leg2Pool();
    TIter iterLeg2(leg2List);
    while((track=(AliReducedBaseTrack*)iterLeg2())) {
      //cout << "leg2 p " << track->P() << "     " << flush; 
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << "/" << flush;
      testFlags1 = mixingMask & track->GetFlags();
      if(!testFlags1) {
	//AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << endl;
	continue;
      }
      for(UShort_t ibit=0; ibit<fNParallelCuts; ++ibit) {
	if((testFlags1)&(ULong_t(1)<<ibit)) 
	  track->UnsetFlag(ibit);
      }
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << endl;
    }  // end loop over the second leg list
  }  // end loop over events
  
  // clean the tracks which don't have enabled mixing flags anymore
  AliReducedBaseTrack* track=0x0;
  // leg1 lists
  iterEv1Leg1Pool.Reset(); TList* leg1List=0x0;
  while((leg1List=(TList*)iterEv1Leg1Pool())) {
    TIter iterLeg1(leg1List);
    while((track=(AliReducedBaseTrack*)iterLeg1())) {
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << flush;
      if(!track->GetFlags()) {
         leg1List->Remove(track); delete track;
      }
    }
  }  // end while
  // leg2 lists
  iterEv1Leg2Pool.Reset(); TList* leg2List=0x0;
  while((leg2List=(TList*)iterEv1Leg2Pool())) {  
    TIter iterLeg2(leg2List);
    while((track=(AliReducedBaseTrack*)iterLeg2())) {
      //AliReducedVarManager::PrintBits(track->GetFlags(),fNParallelCuts); cout << flush;
      if(!track->GetFlags()) {
         leg2List->Remove(track); delete track;
      }
    }
  }  // end while
  
  // clean the events without any tracks left
  for(Int_t i=leg1Pool->GetEntries()-1;i>=0;--i) {
    leg1List=(TList*)leg1Pool->At(i);
    leg2List=(TList*)leg2Pool->At(i);
    if(leg1List->GetEntries()==0 && 
       leg2List->GetEntries()==0) {
       leg1Pool->RemoveAt(i); leg1Pool->Compress(); //delete leg1List;
       leg2Pool->RemoveAt(i); leg2Pool->Compress(); //delete leg2List;
    }
  }
}


//_________________________________________________________________________
void AliMixingHandler::PrintMixingLists(Int_t debugLevel) {
  //
  // Print the contents of the mixing pools
  //
  if(!fIsInitialized) Init();
  
  if(debugLevel<1) {
    cout << "Printing the event mixing lists contents. Prepare for a long dump ..." << endl;
    cout << "Centrality intervals: " << flush;
    for(Int_t icent=0; icent<fCentralityLimits.GetSize(); ++icent)
      cout << fCentralityLimits[icent] << " -- " << flush;
    cout << endl;
    cout << "Event vertex intervals: " << flush;
    for(Int_t ivtx=0; ivtx<fEventVertexLimits.GetSize(); ++ivtx)
      cout << fEventVertexLimits[ivtx] << " -- " << flush;
    cout << endl;
    cout << "Event plane intervals: " << flush;
    for(Int_t iep=0; iep<fEventPlaneLimits.GetSize(); ++iep)
      cout << fEventPlaneLimits[iep] << " -- " << flush;
    cout << endl;
    
    cout << "Mix LS pairs :: " << fMixLikeSign << endl;
    cout << "Pool depth :: " << fPoolDepth << endl;
    cout << "Mixing threshold :: " << fMixingThreshold << endl;
    cout << "Event downscale :: " << fDownscaleEvents << endl;
    cout << "Track downscale :: " << fDownscaleTracks << endl;
    cout << "No. parallel cuts :: " << fNParallelCuts << endl;
    cout << "Histogram class names :: " << fHistClassNames.Data() << endl;
  }
  
  if(debugLevel<1) return;
  
  Int_t nCategories = (fCentralityLimits.GetSize()-1)*(fEventVertexLimits.GetSize()-1)*(fEventPlaneLimits.GetSize()-1);
  AliReducedBaseTrack* track = 0x0;
  
  for(Int_t icent=0; icent<fCentralityLimits.GetSize()-1; ++icent) {
    for(Int_t iz=0; iz<fEventVertexLimits.GetSize()-1; ++iz) {
      for(Int_t iep=0; iep<fEventPlaneLimits.GetSize()-1; ++iep) {
	Int_t evCategory = FindEventCategory(fCentralityLimits[icent],fEventVertexLimits[iz],fEventPlaneLimits[iep]);
        cout << "Event category " << evCategory << " (cent/vtxz/ep) :: [" 
	     << fCentralityLimits[icent] << ";" << fCentralityLimits[icent+1] << "] -- ["
	     << fEventVertexLimits[iz] << ";" << fEventVertexLimits[iz+1] << "] -- ["
	     << fEventPlaneLimits[iep] << ";" << fEventPlaneLimits[iep+1] << "]" << endl;
	cout << "====================================================================================" << endl;     
	cout << "Event bins (cent/vtx/ep) :: " << GetCentralityBin(evCategory) << "/"
	     << GetEventVertexBin(evCategory) << "/" << GetEventPlaneBin(evCategory) << endl;
	cout << "Current pool sizes :: " << flush;
	for(Int_t icut=0;icut<fNParallelCuts;++icut) 
	  cout << fPoolSize[icut*nCategories+evCategory] << " -- " << flush;
	cout << endl;
	if(debugLevel<2) continue;
	
	TClonesArray *leg1PoolP = static_cast<TClonesArray*>(fPoolsLeg1.At(evCategory));
	if(!leg1PoolP) continue;
        TClonesArray &leg1Pool=*leg1PoolP;
	TClonesArray *leg2PoolP = static_cast<TClonesArray*>(fPoolsLeg2.At(evCategory));
	
        TIter iterLeg1Pool(leg1PoolP);
	TIter iterLeg2Pool(leg2PoolP);
	for(Int_t iev=0; iev<leg1Pool.GetEntries(); ++iev) {
	  TList* leg1List = (TList*)iterLeg1Pool();
	  TList* leg2List = (TList*)iterLeg2Pool();
	  cout << "	Event #" << iev << ";  No. of tracks (leg1/leg2) :: " 
	       << leg1List->GetEntries() << " / " << leg2List->GetEntries() << endl;
	  if(debugLevel<3) continue;
	  
	  TIter iterLeg1List(leg1List);
	  cout << "		Leg1 list" << endl;
	  for(Int_t itrack=0; itrack<leg1List->GetEntries(); ++itrack) {
	    track = (AliReducedBaseTrack*)iterLeg1List();
	    cout << "		track #" << itrack << " (p/px/py/pz/charge/flags) :: "
	         << track->P() << " / " << track->Px() << " / " 
                 << track->Py() << " / " << track->Pz() << "/" << track->Charge() << " / " << flush;
	    AliReducedVarManager::PrintBits(track->GetFlags(), fNParallelCuts);	 
	    cout << endl;
	  }  // end loop over tracks
	  
	  TIter iterLeg2List(leg2List);
	  cout << "		Leg2 list" << endl;
	  for(Int_t itrack=0; itrack<leg2List->GetEntries(); ++itrack) {
	    track = (AliReducedBaseTrack*)iterLeg2List();
	    cout << "		track #" << itrack << " (p/px/py/pz/charge/flags) :: "
	         << track->P() << " / " << track->Px() << " / " 
                 << track->Py() << " / " << track->Pz() << "/" << track->Charge() << " / " << flush;
	    AliReducedVarManager::PrintBits(track->GetFlags(), fNParallelCuts);	 
	    cout << endl;
	  }  // end loop over tracks
	  
	}  // end loop over events
      }  // end loop over event plane intervals
    }  // end loop over event vertex intervals
  }  // end loop over centrality intervals
}
