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
#include "AliReducedTrackInfo.h"

ClassImp(AliMixingHandler);

//_________________________________________________________________________
AliMixingHandler::AliMixingHandler(Int_t mixingSetup /* = kMixResonanceLegs*/) :
  TNamed(),
  fMixingSetup(mixingSetup),
  fPoolDepth(100),
  fMixingThreshold(1.0),
  fDownscaleEvents(1.0),
  fDownscaleTracks(1.0),
  fPoolsLeg1("TClonesArray"),
  fPoolsLeg2("TClonesArray"),
  fNParallelCuts(0),
  fNParallelPairCuts(0),
  fHistClassNames(""),
  fPoolSize(),
  fIsInitialized(kFALSE),
  fMixLikeSign(kTRUE),
  fVariableLimits(),
  fVariables(),
  fNMixingVariables(0),
  fHistos(0x0),
  fCrossPairsCuts(),
  fLikePairsLeg1Cuts(),
  fLikePairsLeg2Cuts()
{
  // 
  // default constructor
  //  
  Float_t dummyRange[2] = {-99999., +99999.};
  for(Int_t iVar=0; iVar<kNMaxVariables; ++iVar) {
     fVariableLimits[iVar].Set(2, dummyRange);
     fVariables[iVar] = AliReducedVarManager::kNothing;
  }
  
  fPoolSize.Set(1);
  fCrossPairsCuts.SetOwner(kTRUE);
  fLikePairsLeg1Cuts.SetOwner(kTRUE);
  fLikePairsLeg2Cuts.SetOwner(kTRUE);
}


//_________________________________________________________________________
AliMixingHandler::AliMixingHandler(const Char_t* name, const Char_t* title, Int_t mixingSetup /* = kMixResonanceLegs*/) :
  TNamed(name,title),
  fMixingSetup(mixingSetup),
  fPoolDepth(100),
  fMixingThreshold(1.0),
  fDownscaleEvents(1.0),
  fDownscaleTracks(1.0),
  fPoolsLeg1("TClonesArray"),
  fPoolsLeg2("TClonesArray"),
  fNParallelCuts(0),
  fNParallelPairCuts(0),
  fHistClassNames(""),
  fPoolSize(),
  fIsInitialized(kFALSE),
  fMixLikeSign(kTRUE),
  fVariableLimits(),
  fVariables(),
  fNMixingVariables(0),
  fHistos(0x0),
  fCrossPairsCuts(),
  fLikePairsLeg1Cuts(),
  fLikePairsLeg2Cuts()
{
  //
  // Named constructor
  //  
  Float_t dummyRange[2] = {-99999., +99999.};
  for(Int_t iVar=0; iVar<kNMaxVariables; ++iVar) {
     fVariableLimits[iVar].Set(2, dummyRange);
     fVariables[iVar] = AliReducedVarManager::kNothing;
  }
  
  fPoolSize.Set(1);
  fCrossPairsCuts.SetOwner(kTRUE);
  fLikePairsLeg1Cuts.SetOwner(kTRUE);
  fLikePairsLeg2Cuts.SetOwner(kTRUE);
}


//_________________________________________________________________________
AliMixingHandler::~AliMixingHandler() {
  //
  // destructor
  //
   fCrossPairsCuts.Clear("C");
   fLikePairsLeg1Cuts.Clear("C");
   fLikePairsLeg2Cuts.Clear("C");
}


//_________________________________________________________________________
void AliMixingHandler::AddMixingVariable(AliReducedVarManager::Variables var, Int_t nBins, const Float_t* binLims) {
   //
   // add a mixing variable
   //
   if(fNMixingVariables>=kNMaxVariables) {
      cout << "AliMixingHandler::AddMixingVariable(): ERROR Too many variables for the mixing!" << endl;
      cout << "                  Maximum number of variables: " << kNMaxVariables << endl;
      return;
   }
   fVariables[fNMixingVariables] = var;
   fVariableLimits[fNMixingVariables].Set(nBins, binLims);
   fNMixingVariables++;
   AliReducedVarManager::SetUseVariable(var);
}

//_________________________________________________________________________
void AliMixingHandler::AddMixingVariable(AliReducedVarManager::Variables var, Int_t nBins, const Double_t* binLims) {
   //
   // copy of the AddMixingVariable(AliReducedVarManager::Variables var, Int_t nBins, const Float_t* binLims) function
   //   just to support also Double array as bin limits
   //
   Float_t* bins = new Float_t[nBins];
   for(Int_t i=0;i<nBins;++i) bins[i] = binLims[i];
   AddMixingVariable(var, nBins, bins);
}

//_________________________________________________________________________
void AliMixingHandler::Init() {
  //
  // Initialization of pools
  // NOTE: The master array holding tracks is a 1D array but it will be represented as an n-dim array
  //       The size of the array will be N_1 x N_2 x ... x N_n
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
  Int_t nClassesPerCut = 0;
  if(fMixingSetup==kMixResonanceLegs) nClassesPerCut = 3;
  if(fMixingSetup==kMixCorrelation) nClassesPerCut = 1;
  if (fNParallelPairCuts>1) {
      if(histClassArr->GetEntries()!=nClassesPerCut*fNParallelCuts*fNParallelPairCuts) {
        cout << "AliMixingHandler::Init(): ERROR The number of cuts and the number of hist class names provided do not match!" << endl;
        cout << "                   hist classes: " << histClassArr->GetEntries() << ";    n-parallel cuts: " << fNParallelCuts << ";    n-parallel pair cuts: " << fNParallelPairCuts << endl;
        return;
      }
  } else {
    if(histClassArr->GetEntries()!=nClassesPerCut*fNParallelCuts) {
      cout << "AliMixingHandler::Init(): ERROR The number of cuts and the number of hist class names provided do not match!" << endl;
      cout << "                   hist classes: " << histClassArr->GetEntries() << ";    n-parallel cuts: " << fNParallelCuts << endl;
      return;
    }
  }

  Int_t size = 1;
  for(Int_t iVar = 0; iVar<fNMixingVariables; ++iVar) size *= (fVariableLimits[iVar].GetSize()-1);
  fPoolsLeg1.Expand(size); fPoolsLeg1.SetOwner(kTRUE);
  fPoolsLeg2.Expand(size); fPoolsLeg2.SetOwner(kTRUE);
  
  fPoolSize.Set(fNParallelCuts*size);
  for(Int_t i=0;i<fNParallelCuts*size;++i) fPoolSize[i] = 0;
  
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
  if (!leg1List && !leg2List) return;
  Int_t         entries1 = 0;
  if (leg1List) entries1 = leg1List->GetEntries();
  Int_t         entries2 = 0;
  if (leg2List) entries2 = leg2List->GetEntries();
  if(!entries1 && !entries2) return;
  
  // randomly accept/reject this event in case fDownscaleEvents is used
  if(fDownscaleEvents>1.0 && (gRandom->Rndm()>(1.0/fDownscaleEvents))) 
    return;
  
  // find the event category
  Int_t category = FindEventCategory(values);
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
  AliReducedTrackInfo* track = 0x0;
  if (leg1List) {
    for(Int_t it=0; it<entries1; ++it) {
      // HACK: to transmit the VZERO and TPC event plane Q vector to event mixing
      track = (AliReducedTrackInfo*)leg1List->At(it)->Clone();  
      track->SetCovMatrix(0, values[AliReducedVarManager::kVZEROQvecX+0*6+1]);
      track->SetCovMatrix(1, values[AliReducedVarManager::kVZEROQvecY+0*6+1]);
      track->SetCovMatrix(2, values[AliReducedVarManager::kVZEROQvecX+1*6+1]);
      track->SetCovMatrix(3, values[AliReducedVarManager::kVZEROQvecY+1*6+1]);
      track->SetCovMatrix(4, values[AliReducedVarManager::kTPCQvecXtree+1]);
      track->SetCovMatrix(5, values[AliReducedVarManager::kTPCQvecYtree+1]);
      list1->Add(track);
    }
  }
  if (leg2List) {
    for(Int_t it=0; it<entries2; ++it) {
      // HACK: to transmit the VZERO and TPC event plane Q vector to event mixing
      track = (AliReducedTrackInfo*)leg2List->At(it)->Clone();  
      track->SetCovMatrix(0, values[AliReducedVarManager::kVZEROQvecX+0*6+1]);
      track->SetCovMatrix(1, values[AliReducedVarManager::kVZEROQvecY+0*6+1]);
      track->SetCovMatrix(2, values[AliReducedVarManager::kVZEROQvecX+1*6+1]);
      track->SetCovMatrix(3, values[AliReducedVarManager::kVZEROQvecY+1*6+1]);
      track->SetCovMatrix(4, values[AliReducedVarManager::kTPCQvecXtree+1]);
      track->SetCovMatrix(5, values[AliReducedVarManager::kTPCQvecYtree+1]);
      list2->Add(track);  
    }
  }
    
  // increment the size of the pools in this category
  ULong_t mixingMask = IncrementPoolSizes(leg1List,leg2List,category);
  
  // if full pool(s) were found then run the event mixing
  if(mixingMask) {
    RunEventMixing(leg1PoolP,leg2PoolP,mixingMask,type,values);
    ResetPoolSizes(mixingMask,category);
  }
}


//_________________________________________________________________________
Int_t AliMixingHandler::FindEventCategory(Float_t* values) {
   //
   // Find the event category corresponding to the centrality, vtxz and ep values 
   //
   if(fNMixingVariables==0) return -1;
   if(!fIsInitialized) Init();
   
   Int_t bin[kNMaxVariables]; 
   for (Int_t i=0; i<fNMixingVariables; ++i) {
      bin[i] = TMath::BinarySearch(fVariableLimits[i].GetSize(), fVariableLimits[i].GetArray(), values[fVariables[i]]);
      if(bin[i]==-1 || bin[i]==fVariableLimits[i].GetSize()-1) return -1;      // all variables must be inside limits
   }
   
   Int_t category = 0; 
   for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) {
      Int_t tempCategory = 1;
      for(Int_t iVar2=iVar; iVar2<fNMixingVariables; ++iVar2) {
         if(iVar2==iVar) tempCategory *= bin[iVar2];
         else tempCategory *= (fVariableLimits[iVar2].GetSize()-1);
      }
      category += tempCategory;
   }
   return category;
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetBinFromCategory(Int_t iVar, Int_t category) const {
   //
   // find the bin in variable var for the n-dimensional "category"
   //
   if(fNMixingVariables==0) return -1;
   Int_t norm=1;
   for(Int_t i=fNMixingVariables-1; i>iVar; --i) norm *= (fVariableLimits[i].GetSize()-1);
   Int_t truncatedCategory = category - (category % norm);
   truncatedCategory /= norm;
   return truncatedCategory % (fVariableLimits[iVar].GetSize()-1);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetPoolSize(Int_t cutNumber, Float_t* values) {
   //
   // Get the pool size for a given set of (cut,centrality,z,ep)
   //
   if(cutNumber<0 || cutNumber>fNParallelCuts) return -1;
   Int_t eventCategory = FindEventCategory(values);
   return GetPoolSize(cutNumber, eventCategory);
}


//_________________________________________________________________________
Int_t AliMixingHandler::GetPoolSize(Int_t cutNumber, Int_t eventCategory) const {
  //
  // Get the pool size for a given set of (cut,centrality,z,ep)
  //
  if(cutNumber<0 || cutNumber>fNParallelCuts) return -1;
  if(eventCategory<0) return -1;
  
  Int_t pool = cutNumber;
  for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) pool *= (fVariableLimits[iVar].GetSize() - 1);
  pool += eventCategory;
  
  return fPoolSize[pool];
}


//_________________________________________________________________________
void AliMixingHandler::ResetPoolSizes(ULong_t mixingMask, Int_t category) {
  //
  // Reset pool sizes for a given event category and mask
  //
  Int_t nCategories = 1;
  for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) nCategories *= (fVariableLimits[iVar].GetSize() - 1);
  
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
  if (list1) {
    TIter trackIter1(list1);
    for(Int_t i=0; i<list1->GetEntries();++i) {
      track=(AliReducedBaseTrack*)trackIter1();
      for(UShort_t icut=0;icut<fNParallelCuts;++icut)
        if(track->TestFlag(icut)) cutsMask |= (ULong_t(1)<<icut);
    }
  }
  if (list2) {
    TIter trackIter2(list2);
    for(Int_t i=0; i<list2->GetEntries();++i) {
      track=(AliReducedBaseTrack*)trackIter2();
      for(UShort_t icut=0;icut<fNParallelCuts;++icut)
        if(track->TestFlag(icut)) cutsMask |= (ULong_t(1)<<icut);
    }
  }
    
  // increment the pools for those cuts which got at least one track
  Int_t nCategories = 1;
  for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) nCategories *= (fVariableLimits[iVar].GetSize() - 1);
  Bool_t fullPoolFound = kFALSE;
  for(Int_t icut=0;icut<fNParallelCuts;++icut) {
    if(cutsMask & (ULong_t(1)<<icut))
      fPoolSize[icut*nCategories+eventCategory] += 1;
    if(fPoolSize[icut*nCategories+eventCategory]==fPoolDepth) 
      fullPoolFound = kTRUE;          
  }
  if(!fullPoolFound) return 0;
  
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
  cout << "      Leftover mixing for Mixing Handler " << GetName() << endl;
  cout << "========================================================================" << endl;
  
  // create a mixing mask which enables all cuts
  ULong_t mixingMask = 0;
  for(Int_t i=0; i<fNParallelCuts; ++i) mixingMask |= (ULong_t(1)<<i);
  Float_t values[AliReducedVarManager::kNVars];
  
  for(Int_t icateg=0; icateg<fPoolsLeg1.GetEntries(); ++icateg) {
    TClonesArray *leg1Pool = static_cast<TClonesArray*>(fPoolsLeg1.At(icateg));
    TClonesArray *leg2Pool = static_cast<TClonesArray*>(fPoolsLeg2.At(icateg));
    if(!leg1Pool) continue;
    if(!leg2Pool) continue;
    
    for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) {
       Int_t bin = GetBinFromCategory(iVar, icateg);
       values[fVariables[iVar]] = 0.5*(fVariableLimits[iVar][bin] + fVariableLimits[iVar][bin+1]);
    }
    
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
    
    TIter iterEv2Leg1Pool(leg1Pool);
    TIter iterEv2Leg2Pool(leg2Pool);
    for(Int_t iev2=0; iev2<entries; ++iev2) {                         // second event loop
      TList* ev2Leg1List = (TList*)iterEv2Leg1Pool();
      TList* ev2Leg2List = (TList*)iterEv2Leg2Pool();
      if(iev1==iev2) continue;
      
      //loop over the ev1-leg1 list
      TIter iterLeg1(ev1Leg1List);
      AliReducedBaseTrack* ev1Leg1=0x0;
      while((ev1Leg1=(AliReducedBaseTrack*)iterLeg1())) {
        // check that this track has at least one common bit with the mixing mask
        testFlags1 = mixingMask & ev1Leg1->GetFlags();
        if(!testFlags1) continue;
	
        //loop over the ev2-leg2 list
        TIter iterLeg2(ev2Leg2List);
        AliReducedBaseTrack* ev2Leg2=0x0;
        while((ev2Leg2=(AliReducedBaseTrack*)iterLeg2())) {
          // check that this track has at least one common bit with the mixing mask and with ev1-leg1
          testFlags2 = testFlags1 & ev2Leg2->GetFlags();
          if(!testFlags2) continue;
	  
          // fill cross-pairs (leg1 - leg2) for the enabled bits
          if(fMixingSetup==kMixResonanceLegs) AliReducedVarManager::FillPairInfoME(ev1Leg1, ev2Leg2, type, values);
          if(fMixingSetup==kMixCorrelation)   AliReducedVarManager::FillCorrelationInfo(ev1Leg1, ev2Leg2, values);
          ULong_t pairCutMask = IsPairSelected(values, 1);
          if(!pairCutMask) continue;   // fill histograms only if pair cuts are fulfilled
          for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
            if((testFlags2)&(ULong_t(1)<<ibit)) { 
              //AliReducedVarManager::FillPairMEflow(ev1Leg1, ev2Leg2, values, ibit);
              if(fMixingSetup==kMixResonanceLegs) {
                if (fNParallelPairCuts>1) {
                  for (Int_t jbit=0; jbit<fNParallelPairCuts; jbit++) {
                    if (!((pairCutMask)&(ULong_t(1)<<jbit))) continue;
                    fHistos->FillHistClass(histClassArr->At(ibit*3+jbit*3*fNParallelCuts+1)->GetName(), values);
                  }
                } else {
                  fHistos->FillHistClass(histClassArr->At(ibit*3+1)->GetName(), values);
                }
              }
              if(fMixingSetup==kMixCorrelation) {
                Int_t pairType = (reinterpret_cast<AliReducedPairInfo*>(ev1Leg1))->PairType();
                if (fNParallelPairCuts>1) {
                  ULong_t pairCutMaskCorr = (reinterpret_cast<AliReducedPairInfo*>(ev1Leg1))->GetQualityFlags();
                  for (Int_t jbit=0; jbit<fNParallelPairCuts; jbit++) {
                    if (!((pairCutMaskCorr)&(ULong_t(1)<<jbit))) continue;
                    if (fMixLikeSign) fHistos->FillHistClass(histClassArr->At(ibit*3+jbit*fNParallelCuts+pairType)->GetName(), values);
                    else              fHistos->FillHistClass(histClassArr->At(ibit+jbit*fNParallelCuts)->GetName(), values);
                  }
                } else {
                  if (fMixLikeSign) fHistos->FillHistClass(histClassArr->At(ibit*3+pairType)->GetName(), values);
                  else              fHistos->FillHistClass(histClassArr->At(ibit)->GetName(), values);
                }
              }
            }
          }  
	}  // end loop over the ev2-leg2 list
	
	if(fMixingSetup==kMixCorrelation) continue;
	if(!fMixLikeSign) continue;
	// loop over the ev2-leg1 list
	TIter iterLeg1LS(ev2Leg1List);
	AliReducedBaseTrack* ev2Leg1=0x0;
	while((ev2Leg1=(AliReducedBaseTrack*)iterLeg1LS())) {
	  // check that this track has at least one common bit with the mixing mask and with ev1-leg1
	  testFlags2 = testFlags1 & ev2Leg1->GetFlags();
          if(!testFlags2) continue;
	  
	  // fill like-pairs (leg1 - leg1) for the enabled bits
	  AliReducedVarManager::FillPairInfoME(ev1Leg1, ev2Leg1, type, values);
      ULong_t pairCutMask = IsPairSelected(values, 0);
      if(!pairCutMask) continue;   // fill histograms only if pair cuts are fulfilled
	  for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
        if((testFlags2)&(ULong_t(1)<<ibit)) {
            //AliReducedVarManager::FillPairMEflow(ev1Leg1, ev2Leg1, values, ibit);
            if (fNParallelPairCuts>1) {
                for (Int_t jbit=0; jbit<fNParallelPairCuts; jbit++) {
                    if (!((pairCutMask)&(ULong_t(1)<<jbit))) continue;
                    fHistos->FillHistClass(histClassArr->At(ibit*3+jbit*3*fNParallelCuts+0)->GetName(), values);
                }
            } else {
                fHistos->FillHistClass(histClassArr->At(ibit*3+0)->GetName(), values);
            }
        }
      }
	}  // end loop over the ev2-leg1 list
  }  // end loop over the ev1-leg1 list
      
  if(fMixingSetup==kMixCorrelation) continue;
    if(!fMixLikeSign) continue;
    //loop over the ev1-leg2 list
    TIter iterLeg2(ev1Leg2List);
    AliReducedBaseTrack* ev1Leg2=0x0;
    while((ev1Leg2=(AliReducedBaseTrack*)iterLeg2())) {
	// check that this track has at least one common bit with the mixing mask
	    testFlags1 = mixingMask & ev1Leg2->GetFlags();
        if(!testFlags1) continue;
	
	    //loop over the ev2-leg2 list 
        TIter iterLeg2LS(ev2Leg2List);
        AliReducedBaseTrack* ev2Leg2=0x0;
        while((ev2Leg2=(AliReducedBaseTrack*)iterLeg2LS())) {
            // check that this track has at least one common bit with the mixing mask and with ev1-leg1
            testFlags2 = testFlags1 & ev2Leg2->GetFlags();
            if(!testFlags2) continue;
	  
            // fill like-pairs (leg2 - leg2) for the enabled bits
            AliReducedVarManager::FillPairInfoME(ev1Leg2, ev2Leg2, type, values);
            ULong_t pairCutMask = IsPairSelected(values, 2);
            if(!pairCutMask) continue;   // fill histograms only if pair cuts are fulfilled
            for(Int_t ibit=0; ibit<fNParallelCuts; ++ibit) {
                //AliReducedVarManager::FillPairMEflow(ev1Leg2, ev2Leg2, values, ibit);
                if((testFlags2)&(ULong_t(1)<<ibit)) {
                    if (fNParallelPairCuts>1) {
                        for (Int_t jbit=0; jbit<fNParallelPairCuts; jbit++) {
                            if (!((pairCutMask)&(ULong_t(1)<<jbit))) continue;
                            fHistos->FillHistClass(histClassArr->At(ibit*3+jbit*3*fNParallelCuts+2)->GetName(), values);
                        }
                    } else {
                        fHistos->FillHistClass(histClassArr->At(ibit*3+2)->GetName(), values);
                    }
                }
            }
        }  // end loop over the ev2-leg2 list
     }  // end loop over the ev1-leg2 list
   }  // end second event loop
 }  // end first event loop
  
  // unset the mixing flags --------------------------------------
  iterEv1Leg1Pool.Reset(); 
  iterEv1Leg2Pool.Reset();
  for(Int_t ie1=0; ie1<entries; ++ie1) {
    TList* leg1List = (TList*)iterEv1Leg1Pool();        
    TIter iterLeg1(leg1List);
    AliReducedBaseTrack* track=0x0;
    while((track=(AliReducedBaseTrack*)iterLeg1())) {
      testFlags1 = mixingMask & track->GetFlags();
      if(!testFlags1) continue;
      
      for(UShort_t ibit=0; ibit<fNParallelCuts; ++ibit) {
	if((testFlags1)&(ULong_t(1)<<ibit)) 
	  track->UnsetFlag(ibit);
      }
    }  // end loop over the first leg list
    
    TList* leg2List = (TList*)iterEv1Leg2Pool();
    TIter iterLeg2(leg2List);
    while((track=(AliReducedBaseTrack*)iterLeg2())) {
      testFlags1 = mixingMask & track->GetFlags();
      if(!testFlags1) continue;

      for(UShort_t ibit=0; ibit<fNParallelCuts; ++ibit) {
	if((testFlags1)&(ULong_t(1)<<ibit)) 
	  track->UnsetFlag(ibit);
      }
    }  // end loop over the second leg list
  }  // end loop over events
  
  // clean the tracks which don't have enabled mixing flags anymore
  AliReducedBaseTrack* track=0x0;
  // leg1 lists
  iterEv1Leg1Pool.Reset(); TList* leg1List=0x0;
  while((leg1List=(TList*)iterEv1Leg1Pool())) {
    TIter iterLeg1(leg1List);
    while((track=(AliReducedBaseTrack*)iterLeg1())) {
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
ULong_t AliMixingHandler::IsPairSelected(Float_t* values, Int_t pairType) {
   //
   // apply pair cuts
   //
   TList* cutList = 0;
   switch(pairType) {
      case 0:
         cutList = &fLikePairsLeg1Cuts;
         break;
      case 1:
         cutList = &fCrossPairsCuts;
         break;
      case 2:
         cutList = &fLikePairsLeg2Cuts;
         break;
      default:
         break;
   };
   if(!cutList) return 1;
   if(cutList->GetEntries()==0) return 1;
   
   // loop over all the cuts and make a logical AND between all of them
   ULong_t mask = 0;
   for (Int_t i=0; i<cutList->GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)cutList->At(i);
      if (cut->IsSelected(values)) mask|=(ULong_t(1)<<i);
   }
   return mask;
}


//_________________________________________________________________________
void AliMixingHandler::PrintMixingLists(Int_t debugLevel) {
   //
   // Print the contents of the mixing pools
   //
   if(!fIsInitialized) Init();
  
   cout << "Printing the event mixing lists contents. Prepare for a long dump ..." << endl;
   for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) {
      cout << "Variable #" << iVar << " (" << AliReducedVarManager::fgVariableNames[fVariables[iVar]].Data() << ") intervals: " << flush;
      for(Int_t i=0; i<fVariableLimits[iVar].GetSize(); ++i)
         cout << fVariableLimits[iVar][i] << (i<fVariableLimits[iVar].GetSize()-1 ? " -- " : "") << flush;
      cout << endl;
   }
    
   cout << "Mix LS pairs :: " << fMixLikeSign << endl;
   cout << "Pool depth :: " << fPoolDepth << endl;
   cout << "Mixing threshold :: " << fMixingThreshold << endl;
   cout << "Event downscale :: " << fDownscaleEvents << endl;
   cout << "Track downscale :: " << fDownscaleTracks << endl;
   cout << "No. parallel cuts :: " << fNParallelCuts << endl;
   cout << "Histogram class names :: " << fHistClassNames.Data() << endl;
  
   if(debugLevel<1) return;
  
   Int_t nCategories = 1;
   for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) nCategories *= (fVariableLimits[iVar].GetSize() - 1);

   AliReducedBaseTrack* track = 0x0;
   for(Int_t iCateg=0; iCateg<nCategories; ++iCateg) {
      Int_t bins[fNMixingVariables];
      for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar)  bins[iVar] = GetBinFromCategory(iVar, iCateg);
      cout << "Category #" << iCateg << ", Bins (" << flush;
      for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) 
         cout << bins[iVar] << (iVar<fNMixingVariables-1 ? "/" : "") << flush;
      cout << ")" << endl;
      for(Int_t iVar=0; iVar<fNMixingVariables; ++iVar) 
         cout << "[" << fVariableLimits[iVar][GetBinFromCategory(iVar, iCateg)] << ";" << fVariableLimits[iVar][GetBinFromCategory(iVar, iCateg)+1] << "]" << (iVar<fNMixingVariables-1 ? " -- " : "") << flush;
      cout << endl;
      cout << "====================================================================================" << endl;     
      for(Int_t icut=0;icut<fNParallelCuts;++icut) 
         cout << fPoolSize[icut*nCategories+iCateg] << (icut<fNParallelCuts-1 ? " -- " : "") << flush;
      cout << endl;
      if(debugLevel<2) continue;
      
      TClonesArray *leg1PoolP = static_cast<TClonesArray*>(fPoolsLeg1.At(iCateg));
      if(!leg1PoolP) continue;
      TClonesArray &leg1Pool=*leg1PoolP;
      TClonesArray *leg2PoolP = static_cast<TClonesArray*>(fPoolsLeg2.At(iCateg));
      
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
   }  // end loop over categories  
}
