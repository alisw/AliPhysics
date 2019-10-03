// J/psi analysis in PbPb
// author: Ionut-Cristian Arsene, i.c.arsene@gsi.de
//         2012/Apr/11

#include <iostream>
using namespace std;

#include <TTimeStamp.h>

#ifndef ALICORRELATIONREDUCEDEVENT_H
#include "AliCorrelationReducedEvent.h"
#endif

#include "DstCommonMacros.C"
using namespace DstCommonMacros;

// function prototypes
void DefineHistograms(const Char_t* histClasses, const Char_t* output);
Bool_t IsEventSelected(AliCorrelationReducedEvent* event);
Bool_t IsLegSelected(AliCorrelationReducedTrack* track, Bool_t* legQualityMap);
Bool_t IsLegSelected(AliCorrelationReducedTrack* track);
Bool_t IsLegSelectedTRD(AliCorrelationReducedTrack* track, Int_t trdCut=-1, Int_t trdMinNtr=4, Float_t eleCut=0.7);
Bool_t IsLegSelectedTOF(AliCorrelationReducedTrack* track, Int_t tofCut=-1, Float_t lowExcl=0.3, Float_t highExcl=0.93);
Bool_t IsTrackUsedForMixing(UShort_t* idxArr, UShort_t size, UShort_t idx);
Bool_t IsPairSelected(AliCorrelationReducedPair* pair);
void RunDstPbPbJpsiAnalysis(const Char_t* output, const Char_t* inputfilename,
		            Int_t howMany/* = 5000*/, Int_t offset/* = 0*/);

TFile* gSaveFile=0x0;
// event plane friend file 
const Char_t* gkEPfile="dstTree_VZERO_TPC_recentered.root";

// mixing pool depth
const Int_t gkEventMixingPoolDepth = 100;

// centrality binning for the event mixing
const Int_t gkNCentRanges = 3;
/*Double_t gkEMCentLims[gkNCentRanges+1] = { 0.0,  2.5,  5.0,  7.5, 10.0, 
                                          12.5, 15.0, 17.5, 20.0, 22.5,
					  25.0, 27.5, 30.0, 32.5, 35.0, 
					  37.5, 40.0, 50.0, 60.0, 70.0, 80.0};*/
Double_t gkEMCentLims[gkNCentRanges+1] = { 0.0,  10.0, 40.0, 80.0};
// vertex binning for the event mixing					  
const Int_t gkNVtxRanges = 1;
//Double_t gkEMVtxLims[gkNVtxRanges+1] = {-10.0, -6.0, -2.0, 2.0, 6.0, 10.0};
Double_t gkEMVtxLims[gkNVtxRanges+1] = {-10.0, 10.0};

// event plane angle binning for the event mixing
const Int_t gkNPhiEPRanges = 1;
Double_t gkEMPhiEPLims[gkNPhiEPRanges+1] = {-1.5708, 1.5708};
//Double_t gkEMPhiEPLims[gkNPhiEPRanges+1] = {-1.5708, -1.2566, -0.94248, -0.62832, -0.31415926, 0.0, 0.31415926, 0.62832, 0.94248, 1.2566, 1.5708};
/*Double_t gkEMPhiEPLims[gkNPhiEPRanges+1] = {-1.5708, -1.4137, -1.2566, -1.0996, -0.94248, -0.78540, -0.62832, -0.47124, -0.31416, -0.15708, 0.0, 
                                             0.15708, 0.31416, 0.47124, 0.62832, 0.7854,   0.94248,  1.09956,  1.2566,   1.4137,   1.5708};*/
// which event plane to be used
const Int_t gkEPused = kVZERORP + 6*1 + 1;    // VZEROC 2nd harmonic event plane
//const Int_t gkEPused = kVZERORP + 6*0 + 1;    // VZEROA 2nd harmonic event plane
//const Int_t gkEPused = kTPCRP + 1;    // TPC harmonic event plane

// define detector cuts
const Int_t gkNDetCuts = 6;
const Char_t* gkLegDetCutNames[gkNDetCuts] = {
  "TPC", 
  "TPC_TRD4_1_0.7", "TPC_TRD4_1_0.8", "TPC_TRD4_1_0.9", 
  "TPC_TOF", "TPC_TRD4_1_0.9_TOF"
};

//_________________________________________________________________________________________________
void RunDstPbPbJpsiAnalysis(const Char_t* output, const Char_t* inputfilename, 
		            Int_t howMany/* = 5000*/, Int_t offset/* = 0*/) {
  //
  // J/psi analysis in Pb-Pb
  //
  cout << "start ..." << endl;
  TTimeStamp start;
  cout << "creating chain ..." << endl;
  // create the input chains -----------------------------------------
  Long64_t entries=0;
  TChain* friendChain=0x0;
  if(gkEPfile[0]!='\0')
    friendChain = new TChain("DstFriendTree");
  TChain* chain = GetChain(inputfilename, howMany, offset, entries, friendChain, gkEPfile);
  if(!chain) return;
  if(gkEPfile[0]!='\0' && !friendChain) return;
  AliCorrelationReducedEvent* event = new AliCorrelationReducedEvent();
  chain->SetBranchAddress("Event",&event);
  AliCorrelationReducedEventFriend* eventF = 0x0;
  if(gkEPfile[0]!='\0') {
    eventF = new AliCorrelationReducedEventFriend();
    friendChain->SetBranchAddress("Event",&eventF);
  }
  
  cout << "initialize event mixing lists ..." << endl;
      
  // define histograms -----------------------------------------------
  TString histClasses = "";
  histClasses += "Event_BeforeCuts;Event_AfterCuts;VZERORP;TPCRP;";
  histClasses += "OfflineTriggers;";
  histClasses += "TrackingFlags_Before;TrackQA_JpsiLegs_BeforeCuts_ITS_TPC_TRD_TOF;";
  for(Int_t iLegCut=0; iLegCut<gkNDetCuts; ++iLegCut) {
    histClasses += Form("TrackQA_JpsiLegs_%s_AfterCuts_ITS_TPC_TRD_TOF;TrackingFlags_%s_AfterCuts;",
                        gkLegDetCutNames[iLegCut], gkLegDetCutNames[iLegCut]);
    for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx)
      histClasses += Form("PairQA_SE_%s_vtx%.1f_%.1f;", gkLegDetCutNames[iLegCut], gkEMVtxLims[iVtx], gkEMVtxLims[iVtx+1]);
  }
  
  // initialize event lists for mixing
  //const Int_t kNVarMixingRanges = (gkEventMixingType==0 ? 1 : (gkEventMixingType==1 ? gkNVtxRanges : gkNPhiEPRanges));
  gEMCategories = gkNDetCuts*gkNCentRanges*gkNPhiEPRanges*gkNVtxRanges;
  TList* list1[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges];  // event master lists
  TList* list2[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges];  // event master lists
  TList* selectedPosLegs[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges][gkEventMixingPoolDepth];
  TList* selectedNegLegs[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges][gkEventMixingPoolDepth];
  Int_t mixingPoolSize[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges] = {{{{0}}}};
  UShort_t tracksUsedForMixing[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges][100] = {{{{{0}}}}};
  UShort_t nTracksUsedForMixing[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges] = {{{{0}}}};
  Bool_t eventHasCandidates[gkNDetCuts][gkNCentRanges][gkNVtxRanges][gkNPhiEPRanges] = {{{{kFALSE}}}};
  Bool_t fullListFound = kFALSE;
  for(Int_t iDetCut=0; iDetCut<gkNDetCuts; ++iDetCut) {
    for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) 
      histClasses += Form("PairQA_ME_%s_vtx%.1f_%.1f;", gkLegDetCutNames[iDetCut], gkEMVtxLims[iVtx], gkEMVtxLims[iVtx+1]);
  }
  for(Int_t iDetCut=0; iDetCut<gkNDetCuts; ++iDetCut) {
    for(Int_t iCent=0; iCent<gkNCentRanges; ++iCent) {
      for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) {
	for(Int_t iPhi=0; iPhi<gkNPhiEPRanges; ++iPhi) {
          list1[iDetCut][iCent][iVtx][iPhi] = new TList();
          list2[iDetCut][iCent][iVtx][iPhi] = new TList();
          for(Int_t i=0; i<gkEventMixingPoolDepth; ++i) {
            selectedPosLegs[iDetCut][iCent][iVtx][iPhi][i] = new TList();
            selectedPosLegs[iDetCut][iCent][iVtx][iPhi][i]->SetOwner();
            selectedNegLegs[iDetCut][iCent][iVtx][iPhi][i] = new TList();
            selectedNegLegs[iDetCut][iCent][iVtx][iPhi][i]->SetOwner();
          }
	}  // end loop over phi bins
      }  // end loop over vtx bins
    }  // end loop over centrality bins
  }  // end loop over leg cuts
  DefineHistograms(histClasses.Data(), output);  
  
  Float_t values[kNVars];
  Int_t trackIdMap[20000] = {-1};
  TClonesArray* trackList=0x0;
  TClonesArray* pairList=0x0;
  Bool_t legsQuality[2][gkNDetCuts] = {{kFALSE}};
  Bool_t allLegCutsOr[2] = {kFALSE};
  Bool_t pairQuality[gkNDetCuts] = {kFALSE};
  
  Float_t oldVertex = -999.; Float_t oldBC=-999.; Float_t oldCentVZERO=-999.; Float_t oldCentSPD=-999.; Float_t oldCentTPC=-999.;  
  Double_t mixingTime = 0.0;
  TTimeStamp startEventLoop;
  
  for(Int_t ie=0; ie<entries; ++ie) {  
    chain->GetEntry(ie);
    friendChain->GetEntry(ie);    
    
    gCurrentEvent = event;
    gCurrentEventFriend = eventF;
    if(ie%100==0) cout << "event " << ie << endl;    
    
    FillEventInfo(event, values, eventF);
    FillHistClass("Event_BeforeCuts", values);
    // event cuts
    if(!IsEventSelected(event)) continue;
    // check wheter this event is the same as the previous
    if(TMath::Abs(event->Vertex(2)-oldVertex)<0.0001 && TMath::Abs(event->CentralityVZERO()-oldCentVZERO)<0.001) {
      cout << "This event is the a copy of the previous event" << endl;
      cout << "OLD/NEW:  vtxZ " << oldVertex << "/" << event->Vertex(2)
           << "  BC " << oldBC << "/" << event->BC()
	   << "  CentVZERO " << oldCentVZERO << "/" << event->CentralityVZERO()
	   << "  CentSPD " << oldCentSPD << "/" << event->CentralitySPD()
	   << "  CentTPC " << oldCentTPC << "/" << event->CentralityTPC() << endl;
      ++ie;
      continue;
    }
    else {
      oldVertex = event->Vertex(2); oldBC = event->BC();
      oldCentVZERO = event->CentralityVZERO(); oldCentSPD = event->CentralitySPD();
      oldCentTPC = event->CentralityTPC();
    }
  
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      FillEventOfflineTriggers(ibit, values);
      FillHistClass("OfflineTriggers", values);
    }

    // get the event vtx and centrality bins
    Float_t centVZERO = values[kCentVZERO];
    Float_t vtxZ = values[kVtxZ];
    Int_t binCent = -1;
    for(Int_t i=0; i<gkNCentRanges; ++i) {
      if(centVZERO>=gkEMCentLims[i] && centVZERO<gkEMCentLims[i+1]) {
	binCent = i; break;
      }
    }
    Int_t binVtxZ = -1;
    for(Int_t i=0; i<gkNVtxRanges; ++i) {
      if(vtxZ>=gkEMVtxLims[i] && vtxZ<gkEMVtxLims[i+1]) {
	binVtxZ = i; break;
      }
    }
    if(binCent==-1) {
      cout << "Warning: Centrality bin for this event is -1!! Something went wrong, check it out!" << endl;
      cout << "centVZERO = " << centVZERO << endl;
      continue;
    }
    if(binVtxZ==-1) {
      cout << "Warning: Vertex bin for this event is -1!! Something went wrong, check it out!" << endl;
      cout << "vtxZ = " << vtxZ << endl;
      continue;
    }
    
    Float_t phiEP = values[gkEPused];
    Int_t binPhiEP = -1;                                                                                                            
    for(Int_t i=0; i<gkNPhiEPRanges; ++i) {                                                                                         
      if(values[gkEPused]>=gkEMPhiEPLims[i] &&                                                                  
	values[gkEPused]<gkEMPhiEPLims[i+1]) {                                                                 
	  binPhiEP = i; break;                                                                                                        
      }                                                                                                                                
    }                                                                                                                                  
    if(binPhiEP==-1) {                                                                                                              
      cout << "Warning: EP Phi bin for this event is -1!! Something went wrong, check it out!" << endl;                             
      cout << "phi = " << values[gkEPused] << endl;                                                          
      continue;                                                                                                                        
    }                                                                                                                                  

    // Do the mixing if the events in the pool reached the maximum number of events
    if(fullListFound) {
      TTimeStamp startMix;
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
        for(Int_t iCent=0; iCent<gkNCentRanges; ++iCent) {
	  for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) {
	    for(Int_t iPhi=0; iPhi<gkNPhiEPRanges; ++iPhi) {
              if(mixingPoolSize[iCut][iCent][iVtx][iPhi]>=gkEventMixingPoolDepth) {
                values[kCentVZERO] = 0.5*(gkEMCentLims[iCent]+gkEMCentLims[iCent+1]);
                values[kVtxZ] = 0.5*(gkEMVtxLims[iVtx]+gkEMVtxLims[iVtx+1]);
	        values[gkEPused] = 0.5*(gkEMPhiEPLims[iPhi]+gkEMPhiEPLims[iPhi+1]);
                DoEventMixing(list1[iCut][iCent][iVtx][iPhi], list2[iCut][iCent][iVtx][iPhi], values, 
		              AliCorrelationReducedPair::kJpsiToEE, 
			      Form("PairQA_ME_%s_vtx%.1f_%.1f", gkLegDetCutNames[iCut], gkEMVtxLims[iVtx], gkEMVtxLims[iVtx+1]));
                mixingPoolSize[iCut][iCent][iVtx][iPhi] = 0;   // reset the mixing events number
              }
	    }  // end loop over phi ranges
	  }   // end loop over vertex ranges
        }   // end loop over centralities ranges
        TTimeStamp stopMix;
	mixingTime += Double_t(stopMix.GetSec())+Double_t(stopMix.GetNanoSec())/1.0e+9 - 
	             (Double_t(startMix.GetSec())+Double_t(startMix.GetNanoSec())/1.0e+9);
      }   // end loop over detector cuts 
      fullListFound = kFALSE;
    }  // end if (fullListFound)
    // reset the correct centrality, vtx and phi for this event
    values[kCentVZERO] = centVZERO;
    values[kVtxZ] = vtxZ;
    values[gkEPused] = phiEP;
    
    // make an index map for track access optimization
    AliCorrelationReducedTrack* track = 0x0;
    trackList = event->GetTracks();
    TIter nextTrack(trackList);
    for(Int_t i=0; i<20000; ++i) trackIdMap[i] = -1;
    for(Int_t it=0; it<event->NTracks(); ++it) {
      track = (AliCorrelationReducedTrack*)nextTrack();
      if(track)
        trackIdMap[track->TrackId()] = it;
    }
    
    pairList = event->GetPairs();
    TIter nextPair(pairList);
    AliCorrelationReducedPair* pair = 0x0;
    AliCorrelationReducedTrack* leg1 = 0x0;
    AliCorrelationReducedTrack* leg2 = 0x0;
    values[kNpairsSelected] = 0.0;
    
    // reset arrays
    for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
      for(Int_t iCent=0; iCent<gkNCentRanges; ++iCent) {
	for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) {
	  for(Int_t iPhi=0; iPhi<gkNPhiEPRanges; ++iPhi) {
	    nTracksUsedForMixing[iCut][iCent][iVtx][iPhi] = 0;
	    eventHasCandidates[iCut][iCent][iVtx][iPhi] = kFALSE;
	  }  // end loop over phi bins
	}  // end loop over vtx bins
      }  // end loop over centrality bins
    }  // end loop over detector cuts
        
    // loop over pairs
    for(Int_t ip=0; ip<event->NDielectrons(); ++ip) {
      pair = event->GetDielectronPair(ip);
      // pair cuts
      if(!IsPairSelected(pair)) continue;
      
      // reset flag arrays
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
	legsQuality[0][iCut] = kFALSE; legsQuality[1][iCut] = kFALSE;
	pairQuality[iCut] = kFALSE;
      }
      allLegCutsOr[0] = kFALSE; allLegCutsOr[1] = kFALSE;
      // apply cuts on the first leg
      if(trackIdMap[pair->LegId(0)]!=-1) {
        leg1 = event->GetTrack(trackIdMap[pair->LegId(0)]);
	FillTrackInfo(leg1, values);
	for(UShort_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
	  FillTrackingFlag(leg1, iflag, values);
	  FillHistClass("TrackingFlags_Before", values);
	}
	FillHistClass("TrackQA_JpsiLegs_BeforeCuts_ITS_TPC_TRD_TOF", values);
	allLegCutsOr[0] = IsLegSelected(leg1, legsQuality[0]);
      }   // end if
      
      // apply cuts on the second leg
      if(trackIdMap[pair->LegId(1)]!=-1) {
        leg2 = event->GetTrack(trackIdMap[pair->LegId(1)]);
	FillTrackInfo(leg2, values);
	for(UShort_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
	  FillTrackingFlag(leg2, iflag, values);
	  FillHistClass("TrackingFlags_Before", values);
	}
	FillHistClass("TrackQA_JpsiLegs_BeforeCuts_ITS_TPC_TRD_TOF", values);
	allLegCutsOr[1] = IsLegSelected(leg2, legsQuality[1]);
      }   // end if
      if(!allLegCutsOr[0] || !allLegCutsOr[1]) continue;
      Int_t nCutsPassed = 0;
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
	pairQuality[iCut] = legsQuality[0][iCut] && legsQuality[1][iCut];
	if(pairQuality[iCut]) ++nCutsPassed;
      }
      if(nCutsPassed==0) continue;
    
      FillPairInfo(pair, values);
                        
      // fill track leg histograms
      FillTrackInfo(leg1, values);
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
	if(pairQuality[iCut]) {
	  if(!IsTrackUsedForMixing(tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
	                           nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
			           leg1->TrackId())) {
            FillHistClass(Form("TrackQA_JpsiLegs_%s_AfterCuts_ITS_TPC_TRD_TOF", gkLegDetCutNames[iCut]), values);
            for(UShort_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
	      FillTrackingFlag(leg1, iflag, values);
	      FillHistClass(Form("TrackingFlags_%s_AfterCuts", gkLegDetCutNames[iCut]), values);
            }
	  }
	}
      }
      FillTrackInfo(leg2, values);
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
        if(pairQuality[iCut]) {
	  if(!IsTrackUsedForMixing(tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
	                           nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
			           leg2->TrackId())) {
            FillHistClass(Form("TrackQA_JpsiLegs_%s_AfterCuts_ITS_TPC_TRD_TOF", gkLegDetCutNames[iCut]), values);
            for(UShort_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
	      FillTrackingFlag(leg2, iflag, values);
	      FillHistClass(Form("TrackingFlags_%s_AfterCuts", gkLegDetCutNames[iCut]), values);
            }
	  }
	}
      }
      
      // fill single-event pair histograms
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
	if(pairQuality[iCut]) {
	  if(IsPairSelectedEM(values))
            FillHistClass(Form("PairQA_SE_%s_vtx%.1f_%.1f", gkLegDetCutNames[iCut], gkEMVtxLims[binVtxZ], gkEMVtxLims[binVtxZ+1]), values);
	}
      }
      
      values[kNpairsSelected] += 1.0;
      
      // add legs to the mixing lists
      for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
	if(pairQuality[iCut]) {
	  if(!IsTrackUsedForMixing(tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
	                           nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
			           leg1->TrackId())) {
	    if(leg1->Charge()>0) selectedPosLegs[iCut][binCent][binVtxZ][binPhiEP][mixingPoolSize[iCut][binCent][binVtxZ][binPhiEP]]->Add(leg1->Clone());
            else selectedNegLegs[iCut][binCent][binVtxZ][binPhiEP][mixingPoolSize[iCut][binCent][binVtxZ][binPhiEP]]->Add(leg1->Clone());
	    tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP][nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP]] = leg1->TrackId();
	    nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP] += 1;
	  }
	  if(!IsTrackUsedForMixing(tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
	                           nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP], 
			           leg2->TrackId())) {
	    if(leg2->Charge()>0) selectedPosLegs[iCut][binCent][binVtxZ][binPhiEP][mixingPoolSize[iCut][binCent][binVtxZ][binPhiEP]]->Add(leg2->Clone());
            else selectedNegLegs[iCut][binCent][binVtxZ][binPhiEP][mixingPoolSize[iCut][binCent][binVtxZ][binPhiEP]]->Add(leg2->Clone());
	    tracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP][nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP]] = leg2->TrackId();
	    nTracksUsedForMixing[iCut][binCent][binVtxZ][binPhiEP] += 1;
	  }
	  eventHasCandidates[iCut][binCent][binVtxZ][binPhiEP] = kTRUE;
	}  // end if(pairQuality)
      }
    }  // end loop over pairs
    
    // If this event has candidates, then store the event in the list for event mixing
    for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
      for(Int_t iCent=0; iCent<gkNCentRanges; ++iCent) {
	for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) {
	  for(Int_t iPhi=0; iPhi<gkNPhiEPRanges; ++iPhi) {
            if(!eventHasCandidates[iCut][iCent][iVtx][iPhi]) continue;
	    list1[iCut][iCent][iVtx][iPhi]->Add(selectedPosLegs[iCut][iCent][iVtx][iPhi][mixingPoolSize[iCut][iCent][iVtx][iPhi]]);
            list2[iCut][iCent][iVtx][iPhi]->Add(selectedNegLegs[iCut][iCent][iVtx][iPhi][mixingPoolSize[iCut][iCent][iVtx][iPhi]]);
            /*
	    cout << "event mixing category (cut/cent/vtx) " << iCut << "/" << iCent << "/" << iVtx << ", pool size (pos/neg): " 
                 << list1[iCut][iCent][iVtx]->GetEntries() 
                 << "/" << list2[iCut][iCent][iVtx]->GetEntries() << "; ntracks this entry (pos/neg): " 
                 << selectedPosLegs[iCut][iCent][iVtx][mixingPoolSize[iCut][iCent][iVtx]]->GetEntries() << "/" 
                 << selectedNegLegs[iCut][iCent][iVtx][mixingPoolSize[iCut][iCent][iVtx]]->GetEntries() << endl;
	    */
            mixingPoolSize[iCut][iCent][iVtx][iPhi] += 1;
	    if(mixingPoolSize[iCut][iCent][iVtx][iPhi]==gkEventMixingPoolDepth)
	      fullListFound = kTRUE;
	  }  // end loop over phi ranges
	}  // end loop over vertex ranges
      }  // end loop over centrality ranges
    }  // end loop over detector cuts
    
    FillHistClass("Event_AfterCuts", values);
    FillHistClass("TPCRP", values);
    FillHistClass("VZERORP", values);
  }  // end loop over events
  
  // Mixing the leftover events
  cout << "Leftover mixing ..." << endl;
  TTimeStamp startMix;
  for(Int_t iCut=0; iCut<gkNDetCuts; ++iCut) {
    for(Int_t iCent=0; iCent<gkNCentRanges; ++iCent) {
      for(Int_t iVtx=0; iVtx<gkNVtxRanges; ++iVtx) {
	for(Int_t iPhi=0; iPhi<gkNPhiEPRanges; ++iPhi) {
          values[kCentVZERO] = 0.5*(gkEMCentLims[iCent]+gkEMCentLims[iCent+1]);
	  values[kVtxZ] = 0.5*(gkEMVtxLims[iVtx]+gkEMVtxLims[iVtx+1]);
	  values[gkEPused] = 0.5*(gkEMPhiEPLims[iPhi]+gkEMPhiEPLims[iPhi+1]);
          DoEventMixing(list1[iCut][iCent][iVtx][iPhi], list2[iCut][iCent][iVtx][iPhi], values, 
	                AliCorrelationReducedPair::kJpsiToEE, 
			Form("PairQA_ME_%s_vtx%.1f_%.1f", gkLegDetCutNames[iCut], gkEMVtxLims[iVtx], gkEMVtxLims[iVtx+1])); // after the event mixing the lists will be cleared
          mixingPoolSize[iCut][iCent][iVtx][iPhi] = 0;   // reset the mixing events number
	}  // end loop over phi ranges
      }   // end loop over vertex ranges
    }   // end loop over centralities ranges
  }   // end loop over detector cuts 
  TTimeStamp stopMix;
  mixingTime += Double_t(stopMix.GetSec())+Double_t(stopMix.GetNanoSec())/1.0e+9 - 
	       (Double_t(startMix.GetSec())+Double_t(startMix.GetNanoSec())/1.0e+9);
  
  TTimeStamp stopEventLoop;
  
  WriteOutput(gSaveFile);
  
  cout << "Initialization time: " << startEventLoop.GetSec() - start.GetSec() << " seconds" << endl;
  cout << "Total looping time: " << stopEventLoop.GetSec() - startEventLoop.GetSec() << " seconds" << endl;
  cout << "Time spent in mixing: " << mixingTime << " seconds" << endl;
  cout << "Overall speed       : " << Double_t(stopEventLoop.GetSec() - startEventLoop.GetSec())/Double_t(entries) << " sec./event" << endl;
}


//__________________________________________________________________
Bool_t IsTrackUsedForMixing(UShort_t* idxArray, UShort_t size, UShort_t idx) {
  //
  // check whether track with idx was already added to the mixing list
  //
  for(Int_t i=0; i<size; ++i) {
    if(idxArray[i]==idx) return kTRUE;
  }
  return kFALSE;
}


//__________________________________________________________________
Bool_t IsEventSelected(AliCorrelationReducedEvent* event) {
  //
  // event selection
  //
  if(!event->IsPhysicsSelection()) return kFALSE;
  if(event->VertexNContributors()<1) return kFALSE;
  if(event->CentralityQuality()!=0) return kFALSE;
  if(event->CentralityVZERO()<-0.0001) return kFALSE;
  if(event->CentralityVZERO()>80.0) return kFALSE;
  if(TMath::Abs(event->Vertex(2))>10.0) return kFALSE;
  return kTRUE;
}


//__________________________________________________________________
Bool_t IsLegSelected(AliCorrelationReducedTrack* track, Bool_t* legQuality) {
  //
  //  track leg selection
  //  
  legQuality[0] = IsLegSelected(track);
  legQuality[1] = legQuality[0] && IsLegSelectedTRD(track, 1, 4, 0.7);
  legQuality[2] = legQuality[0] && IsLegSelectedTRD(track, 1, 4, 0.80);
  legQuality[3] = legQuality[0] && IsLegSelectedTRD(track, 1, 4, 0.90);
  legQuality[4] = legQuality[0] && IsLegSelectedTOF(track, 2);
  legQuality[5] = legQuality[3] && legQuality[4];
  Bool_t globalOr = kFALSE;
  for(Int_t i=0; i<gkNDetCuts; ++i) globalOr = globalOr || legQuality[i];
  return globalOr;
}


//__________________________________________________________________
Bool_t IsLegSelected(AliCorrelationReducedTrack* track) {
  //
  // track leg selection
  //
  if(!track->CheckTrackStatus(kTPCrefit)) return kFALSE;
  if(!track->CheckTrackStatus(kITSrefit)) return kFALSE;
  if(!(track->ITSLayerHit(0) || track->ITSLayerHit(1))) return kFALSE;
  
  if(track->Pt()<0.7) return kFALSE;
  if(TMath::Abs(track->Eta())>0.9) return kFALSE;
  
  if(track->TPCncls()<70) return kFALSE;
  
  if(track->TPCnSig(kElectron)>3.0) return kFALSE;
  if(track->TPCnSig(kElectron)<-2.0) return kFALSE;
  
  if(TMath::Abs(track->DCAz())>3.0) return kFALSE;
  if(TMath::Abs(track->DCAxy())>1.0) return kFALSE;
  
  //if(track->Pin()<2.5) {
  if(TMath::Abs(track->TPCnSig(kProton))<3.5) return kFALSE;
  if(TMath::Abs(track->TPCnSig(kPion))<3.5) return kFALSE;
  if(track->TPCnSig(kProton)<-3.5) return kFALSE;
  //if(track->TPCnSig(kProton)<-3.5) return kFALSE;
  //}
  
  return kTRUE;
}


//__________________________________________________________________
Bool_t IsLegSelectedTRD(AliCorrelationReducedTrack* track, 
		        Int_t trdCut/*=-1*/, Int_t trdMinNtr/*=4*/, Float_t eleCut/*=0.7*/) {
  //
  // TRD leg selection
  //
  if(trdCut==1 && track->TRDntracklets(0)>=trdMinNtr) { 
    if(track->TRDpid(0)<eleCut) return kFALSE;
  }
  if(trdCut==2 && track->TRDntracklets(0)>=trdMinNtr && track->Pin()>1.0) {
    if(track->TRDpid(0)<eleCut) return kFALSE;
  }
  if(trdCut==3 && track->CheckTrackStatus(kTRDpid) && track->TRDntracklets(0)>=trdMinNtr) { 
    if(track->TRDpid(0)<eleCut) return kFALSE;
  }
  if(trdCut==4 && track->CheckTrackStatus(kTRDpid) && track->TRDntracklets(1)>=trdMinNtr && track->Pin()>1.0) {
    if(track->TRDpid(0)<eleCut) return kFALSE;
  }
  return kTRUE;
}


//__________________________________________________________________
Bool_t IsLegSelectedTOF(AliCorrelationReducedTrack* track,
			Int_t tofCut/*=-1*/, Float_t lowExcl/*=0.3*/, Float_t highExcl/*=0.93*/) {
  //
  // TOF leg selection
  //
  if(tofCut==1 && track->TOFbeta()>lowExcl && track->TOFbeta()<highExcl) return kFALSE;
  if(tofCut==2 && track->CheckTrackStatus(kTOFpid) && TMath::Abs(track->TOFnSig(kElectron))>3.5) return kFALSE;
  return kTRUE;
}


//__________________________________________________________________
Bool_t IsPairSelected(AliCorrelationReducedPair* pair) {
  //
  // pair selection
  //
  if(!pair) return kFALSE;
  if(pair->CandidateId()!=AliCorrelationReducedPair::kJpsiToEE) return kFALSE;
  if(pair->PairType()!=1) return kFALSE;
  if(TMath::Abs(pair->Rapidity())>0.9) return kFALSE;
  return kTRUE;
}


//__________________________________________________________________
void DefineHistograms(const Char_t* histClasses, const Char_t* output) {
  //
  // define the histograms
  //
  cout << "Defining histograms ..." << flush;
  cout << "histogram classes: " << histClasses << endl;
  
  gSaveFile=new TFile(output,"RECREATE");
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};
  
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;
    
    // Event wise histograms
    if(classStr.Contains("Event")) {
      cout << "Event" << endl;
      AddHistClass(classStr.Data());
      AddHistogram(classStr.Data(),"RunNo","Run numbers;Run",kFALSE, kNRunBins, runHistRange[0], runHistRange[1], kRunNo);
      AddHistogram(classStr.Data(),"BC","Bunch crossing;BC",kFALSE,3000,0.,3000.,kBC);
      AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;",kFALSE,
		   2,-0.5,1.5,kIsPhysicsSelection, 0,0.0,0.0,kNothing, 0,0.0,0.0,kNothing, "off;on");
      
      AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)",kFALSE,300,-15.,15.,kVtxZ);
      AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)",kFALSE,300,-0.4,0.4,kVtxX);
      AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)",kFALSE,300,-0.4,0.4,kVtxY);
      
      AddHistogram(classStr.Data(),"VtxZ_Run_prof", "<VtxZ> vs run; Run; vtxZ (cm)", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, -20, 20., 1000., kVtxZ);
      AddHistogram(classStr.Data(),"VtxX_Run_prof", "<VtxX> vs run; Run; vtxX (cm)", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, -20, 20., 1000., kVtxX);
      AddHistogram(classStr.Data(),"VtxY_Run_prof", "<VtxY> vs run; Run; vtxY (cm)", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, -20, 20., 1000., kVtxY);
      
      AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)",kFALSE,
                   100, 0.0, 100.0, kCentVZERO);
      AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality",kFALSE,
                   100, -50.5, 49.5, kCentQuality);
      AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)",kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0.0, 100.0, kCentVZERO);      
      
      AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs",kFALSE,
                   5000,0.,5000.,kNdielectrons);
      AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
	           5000,0.,5000.,kNpairsSelected);
      AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks",kFALSE,
                   1000,0.,30000.,kNtracksTotal);
      AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks",kFALSE,
                   1000,0.,30000.,kNtracksSelected);
      AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0; tracklets", kFALSE,
                   3000, -0.5, 2999.5, kSPDntracklets);
      
      AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0., 10000., kNdielectrons);
      AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0., 10000., kNpairsSelected);
      AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0., 10000., kNtracksTotal);
      AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0., 10000., kNtracksSelected);
      AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, 0., 10000., kSPDntracklets);
      
      AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)",kFALSE,
                   300,-15.,15.,kVtxZ, 100, 0.0, 100.0, kCentVZERO);
      AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)",kFALSE,
                   300,-15.,15.,kVtxZ, 100, 0.0, 100.0, kCentSPD);
      AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)",kFALSE,
                   300,-15.,15.,kVtxZ, 100, 0.0, 100.0, kCentTPC);
      continue;
    }  // end if className contains "Event"
    
    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      cout << "OfflineTriggers" << endl;
      AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += gkOfflineTriggerNames[i]; triggerNames+=";";}
      
      AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
	           64, -0.5, 63.5, kOfflineTrigger, 2, -0.5, 1.5, kOfflineTriggerFired, 0, 0.0, 0.0, kNothing, triggerNames.Data(), "off;on");
      AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
	           64, -0.5, 63.5, kOfflineTriggerFired2, 0, 0.0, 0.0, kNothing, 0, 0.0, 0.0, kNothing, triggerNames.Data());
      AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
	           64, -0.5, 63.5, kOfflineTriggerFired2, 20, 0.0, 100.0, kCentVZERO, 0, 0.0, 0.0, kNothing, triggerNames.Data());
      AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
	           64, -0.5, 63.5, kOfflineTriggerFired2, 200, -20.0, 20.0, kVtxZ, 0, 0.0, 0.0, kNothing, triggerNames.Data());
      continue;
    }
    
    if(classStr.Contains("VZERORP")) {
      cout << "VZERORP" << endl;
      AddHistClass(classStr.Data());
      const Char_t* sname[3] = {"A","C","A&C"};
      for(Int_t ih=1; ih<2; ++ih) {
        for(Int_t iS=0; iS<3; ++iS) {
          AddHistogram(classStr.Data(), Form("QvecX_side%s_h%d_CentSPDVtxZ_prof",sname[iS],ih+1), 
                       Form("Q_{x}, side %s, harmonic %d, vs centSPD and vtxZ; Centrality (percents); VtxZ (cm); <Q_{x}>",sname[iS],ih+1), kTRUE, 
                       20, 0.0, 100.0, kCentSPD, 24, -12.0, 12.0, kVtxZ, 500, -10000.0, 10000.0, kVZEROQvecX+iS*6+ih);
	  AddHistogram(classStr.Data(), Form("QvecY_side%s_h%d_CentSPDVtxZ_prof",sname[iS],ih+1), 
                       Form("Q_{y}, side %s, harmonic %d, vs centSPD and vtxZ; Centrality (percents); VtxZ (cm); <Q_{y}>",sname[iS],ih+1), kTRUE, 
                       20, 0.0, 100.0, kCentSPD, 24, -12.0, 12.0, kVtxZ, 500, -10000.0, 10000.0, kVZEROQvecY+iS*6+ih);
	  AddHistogram(classStr.Data(), Form("QvecX_side%s_h%d_Run_prof", sname[iS], ih+1), 
		       Form("<Q_{x}>, VZERO side %s, harmonic %d, vs run; Run; <Q_{x}>", sname[iS], ih+1), kTRUE,
	               kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, -100., 100., kVZEROQvecX+iS*6+ih);
	  AddHistogram(classStr.Data(), Form("QvecY_side%s_h%d_Run_prof", sname[iS], ih+1), 
		       Form("<Q_{y}>, VZERO side %s, harmonic %d, vs run; Run; <Q_{y}>", sname[iS], ih+1), kTRUE,
	               kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, -100., 100., kVZEROQvecY+iS*6+ih);
	  AddHistogram(classStr.Data(), Form("RP_side%s_h%d_CentSPD",sname[iS],ih+1), 
                       Form("VZERO reaction plane, side %s, harmonic %d, vs centrality SPD; #Psi (rad.); centSPD (percents)",sname[iS],ih+1), kFALSE, 
                       400, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), kVZERORP+iS*6+ih, 20, 0.0, 100.0, kCentSPD);
	  AddHistogram(classStr.Data(), Form("RP_side%s_h%d_VtxZ",sname[iS],ih+1), 
                       Form("VZERO reaction plane, side %s, harmonic %d, vs vtxZ; #Psi (rad.); vtxZ (cm)",sname[iS],ih+1), kFALSE, 
                       400, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), kVZERORP+iS*6+ih, 24, -12.0, +12.0, kVtxZ);          
        }   // end loop over VZERO sides
      }   // end loop over harmonics    
      continue;
    }   // end if for the VZERO reaction plane histograms
    
    if(classStr.Contains("TPCRP")) {
      cout << "TPCRP" << endl;
      AddHistClass(classStr.Data());
      for(Int_t ih=1; ih<2; ++ih) {
	AddHistogram(classStr.Data(), Form("QvecX_TPC_h%d_Run_prof", ih+1), 
		     Form("<Q_{x}>, TPC, harmonic %d, vs run; Run; <Q_{x}>", ih+1), kTRUE,
	             kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, -100., 100., kTPCQvecX+ih);
	AddHistogram(classStr.Data(), Form("QvecY_TPC_h%d_Run_prof", ih+1), 
		     Form("<Q_{y}>, TPC, harmonic %d, vs run; Run; <Q_{y}>", ih+1), kTRUE,
	             kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 100, -100., 100., kTPCQvecY+ih);
	AddHistogram(classStr.Data(), Form("TPCRP_h%d", ih+1), 
                     Form("TPC event plane, harmonic %d; #Psi (rad.)", ih+1), kFALSE, 
                     400, -4.0/Double_t(ih+1), 4.0/Double_t(ih+1), kTPCRP+ih);
	AddHistogram(classStr.Data(), Form("TPCQvecX_h%d_CentSPDVtxZ_prof", ih+1), 
                     Form("TPC Q_{x}, harmonic %d, vs centSPD and vtxZ; Centrality (percents); VtxZ (cm); <Q_{x}>", ih+1), kTRUE, 
                     20, 0.0, 100.0, kCentSPD, 24, -12.0, 12.0, kVtxZ, 500, -600.0, 600.0, kTPCQvecX+ih);
	AddHistogram(classStr.Data(), Form("TPCQvecY_h%d_CentSPDVtxZ_prof", ih+1), 
                     Form("TPC Q_{y}, harmonic %d, vs centSPD and vtxZ; Centrality (percents); VtxZ (cm); <Q_{y}>", ih+1), kTRUE, 
                     20, 0.0, 100.0, kCentSPD, 24, -12.0, 12.0, kVtxZ, 500, -600.0, 600.0, kTPCQvecY+ih);
      }    // end loop over harmonics        
      continue;
    }    // end if for the TPC reaction plane histograms
    
    TString trkFlagNames = "";
    for(Int_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
      trkFlagNames += gkTrackingFlagNames[iflag];
      trkFlagNames += ";";
    }
    
    // Track histograms
    if(classStr.Contains("TrackingFlags")) {
      AddHistClass(classStr.Data());
      
      AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
	           kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, 0, 0.0, 0.0, kNothing, trkFlagNames.Data());
      AddHistogram(classStr.Data(), "TrackingFlags_Pt", "Tracking flags vs p_{T};p_{T} (GeV/c);", kFALSE,
	           100, 0.0, 20.0, kPt, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
      AddHistogram(classStr.Data(), "TrackingFlags_Eta", "Tracking flags vs #eta;#eta;", kFALSE,
	           30, -1.5, 1.5, kEta, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
      AddHistogram(classStr.Data(), "TrackingFlags_Phi", "Tracking flags vs #varphi;#varphi (rad.);", kFALSE,
	           60, 0.0, 6.3, kPhi, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
      AddHistogram(classStr.Data(), "TrackingFlags_CentVZERO", "Tracking flags vs centrality VZERO; centrality (%);", kFALSE,
	           20, 0.0, 100.0, kCentVZERO, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
      
      AddHistogram(classStr.Data(), "TrackingFlags_TRDntracklets", "Tracking flags vs TRD #tracklets; TRD #tracklets;", kFALSE,
	           7, -0.5, 6.5, kTRDntracklets, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());      
      AddHistogram(classStr.Data(), "TrackingFlags_TRDntrackletsPID", "Tracking flags vs TRD # pid tracklets; TRD #tracklets;", kFALSE,
	           7, -0.5, 6.5, kTRDntrackletsPID, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());      
      AddHistogram(classStr.Data(), "TrackingFlags_TRDeleProbab", "Tracking flags vs TRD electron probability; TRD probab.;", kFALSE,
	           50, 0.0, 1.0, kTRDpidProbabilities, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
      AddHistogram(classStr.Data(), "TrackingFlags_TOFbeta", "Tracking flags vs TOF #beta; TOF #beta;", kFALSE,
	           50, 0.0, 1.0, kTOFbeta, kNTrackingFlags, -0.5, kNTrackingFlags-0.5, kTrackingFlag, 0, 0.0, 0.0, kNothing, "", trkFlagNames.Data());
    }
    
    if(classStr.Contains("TrackQA")) {
      AddHistClass(classStr.Data());
    
      AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
                   1000, 0.0, 50.0, kPt);
      AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
                   1000, -1.5, 1.5, kEta);
      AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
                   1000, 0.0, 6.3, kPhi);
      AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
                   1000, -10.0, 10.0, kDcaXY);
      AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
                   1000, -10.0, 10.0, kDcaZ);

      // run dependence
      AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, 0.0, 50.0, kPt);
      AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, -1.5, 1.5, kEta);      
      AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, 0.0, 6.3, kPhi);      
      AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, -10.0, 10.0, kDcaXY);
      AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, -10.0, 10.0, kDcaZ);

      // correlations between parameters
      AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
                   300, -1.5, +1.5, kEta, 100, 0.0, 10.0, kPt);
      AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
                   300, 0.0, 6.3, kPhi, 100, 0.0, 10.0, kPt);
      AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
                   200, -1.0, +1.0, kEta, 100, 0.0, 6.3, kPhi);
      AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
                   100, 0.0, 10.0, kPt, 500, -2.0, 2.0, kDcaXY);
      AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
                   100, 0.0, 10.0, kPt, 500, -2.0, 2.0, kDcaZ);
      AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
                   100, -1.0, 1.0, kEta, 500, -2.0, 2.0, kDcaXY);
      AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
                   100, -1.0, 1.0, kEta, 500, -2.0, 2.0, kDcaZ);
      
      if(classStr.Contains("ITS")) {
        AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters;# clusters ITS", kFALSE,
                     7,-0.5,6.5,kITSncls);
	AddHistogram(classStr.Data(),"ITSncls_Run", "ITS <nclusters> vs run;run;# clusters ITS", kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 7,-0.5,6.5, kITSncls);
	AddHistogram(classStr.Data(),"ITSncls_Cent_prof", "ITS <nclusters> vs centrality; centrality;# clusters ITS", kTRUE,
                     20, 0.0, 100.0, kCentVZERO, 7,-0.5,6.5, kITSncls);
        AddHistogram(classStr.Data(),"Pt_ITSncls","ITS nclusters vs p_{T};p_{T} (GeV/c);# clusters ITS",kFALSE,
                     100, 0.0, 10.0, kPt, 7,-0.5,6.5, kITSncls);
        AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi);#eta;#phi (rad.);# clusters ITS",kTRUE,
                     192, -1.2, 1.2, kEta, 126, 0.0, 6.3, kPhi, 7, -0.5, 6.5, kITSncls);
      }  // end if ITS histograms
      
      if(classStr.Contains("TPC")) {
        AddHistogram(classStr.Data(),"TPCncls","TPC nclusters;# clusters TPC",kFALSE,
                     160,-0.5,159.5,kTPCncls);
	AddHistogram(classStr.Data(),"TPCncls_Run","TPC nclusters vs run;run;# clusters TPC",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 160,-0.5,159.5, kTPCncls);
	AddHistogram(classStr.Data(),"TPCncls_CentVZERO","TPC nclusters vs centrality;centrality;# clusters TPC",kFALSE,
                     20, 0.0, 100.0, kCentVZERO, 160,-0.5,159.5, kTPCncls);
	AddHistogram(classStr.Data(),"TPCncls_Pin","TPC nclusters vs inner param P; P (GeV/c); # clusters TPC",kFALSE,
                     200,0.0,20,kPin,160,-0.5,159.5,kTPCncls);
	AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi);#eta;#phi (rad.);# clusters TPC",kTRUE,
                     192, -1.2, 1.2, kEta, 126, 0.0, 6.3, kPhi, 160, -0.5, 159.5, kTPCncls);    
        AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P;P (GeV/c);TPC dE/dx",kFALSE,
                     1000,0.0,20.0,kPin,151,-0.5,150.5,kTPCsignal);
        AddHistogram(classStr.Data(),"TPCsignal_TPCclusters","TPC dE/dx vs. TPC nclusters;TPC #clusters;TPC dE/dx",kFALSE,
                     160,-0.5,159.5,kTPCncls,151,-0.5,150.5,kTPCsignal);
        AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P;P (GeV/c);N_{#sigma}",kFALSE,
                     1000,0.0,20.0,kPin,100,-5.0,5.0,kTPCnSig+kElectron);
	AddHistogram(classStr.Data(),"TPCnsigElectron_CentEta_prof","TPC N_{#sigma} electron vs. (#eta,centrality); #eta; centrality VZERO;N_{#sigma}",kTRUE,
                     100,-1.0,1.0,kEta, 20, 0.0, 100.0, kCentVZERO, 100,-5.0,5.0,kTPCnSig+kElectron);
	AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run; run; N_{#sigma}",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo,100,-5.0,5.0,kTPCnSig+kElectron);
      }      // end if TPC histograms
    
      if(classStr.Contains("TRD")) {
        AddHistogram(classStr.Data(),"TRDntracklets","TRD ntracklets; #tracklets TRD",kFALSE,
                     7,-0.5,6.5,kTRDntracklets);
	AddHistogram(classStr.Data(),"TRDntrackletsPID","TRD ntracklets PID; #tracklets TRD",kFALSE,
                     7,-0.5,6.5,kTRDntrackletsPID);
	AddHistogram(classStr.Data(),"TRDprobabElectron","TRD electron probability; probability",kFALSE,
                     500,0.0,1.0,kTRDpidProbabilities);
	AddHistogram(classStr.Data(),"TRDprobabPion","TRD pion probability; probability",kFALSE,
                     500,0.0,1.0,kTRDpidProbabilities+1);
        AddHistogram(classStr.Data(),"TRDntracklets_Run","TRD <ntracklets> vs run; run; #tracklets TRD",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 7,-0.5,6.5,kTRDntracklets);
        AddHistogram(classStr.Data(),"TRDntrackletsPID_Run","TRD <ntracklets PID> vs run; run; #tracklets TRD",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 7,-0.5,6.5,kTRDntrackletsPID);
        AddHistogram(classStr.Data(),"TRDprobabEle_Run","TRD <electron probability> vs run; run; probability",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 10,0.0,1.0,kTRDpidProbabilities);
	AddHistogram(classStr.Data(),"TRDprobabPion_Run","TRD <pion probability> vs run; run; probability",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 10,0.0,1.0,kTRDpidProbabilities+1);
	AddHistogram(classStr.Data(),"Eta_Phi_TRDntracklets_prof","TRD <ntracklets> vs (#eta,#phi);#eta;#phi (rad.);#tracklets TRD",kTRUE,
                     192, -1.2, 1.2, kEta, 126, 0.0, 6.3, kPhi, 7, -0.5, 6.5, kTRDntracklets);   
        AddHistogram(classStr.Data(),"TRDntracklets_cent","TRD ntracklets vs centrality; #tracklets TRD; centrality (%)",kFALSE,
                     7,-0.5,6.5,kTRDntracklets,20,0.0, 100.0, kCentSPD);
        AddHistogram(classStr.Data(),"Eta_Phi_TRDntrackletsPID_prof","TRD <ntracklets PID> vs (#eta,#phi);#eta;#phi (rad.);#tracklets TRD",kTRUE,
                     192, -1.2, 1.2, kEta, 126, 0.0, 6.3, kPhi, 7, -0.5, 6.5, kTRDntrackletsPID);
        AddHistogram(classStr.Data(),"TRDntrackletsPID_cent","TRD ntracklets PID vs centrality; #tracklets TRD; centrality (%)",kFALSE,
                     7,-0.5,6.5,kTRDntrackletsPID,20,0.0, 100.0, kCentSPD);
        AddHistogram(classStr.Data(),"TRDntracklets_TRDntrackletsPID","TRD ntracklets vs TRD ntracklets PID; #tracklets TRD; #tracklets TRD pid",kFALSE,
                     7,-0.5,6.5,kTRDntracklets,7,-0.5, 6.5, kTRDntrackletsPID);
      }   // end if TRD histograms
    
      if(classStr.Contains("TOF")) {
        AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P;P (GeV/c);#beta",kFALSE,
                     200,0.0,20.0,kP, 220,0.0,1.1,kTOFbeta);
	AddHistogram(classStr.Data(),"Eta_Phi_TOFbeta_prof","TOF <#beta> vs (#eta,#phi);#eta;#phi (rad.);TOF #beta",kTRUE,
                     192, -1.2, 1.2, kEta, 126, 0.0, 6.3, kPhi, 160, -0.5, 159.5, kTOFbeta);
        AddHistogram(classStr.Data(),"TOFnsigElectron_P","TOF N_{#sigma} electron vs. P;P (GeV/c);N_{#sigma}",kFALSE,
                     200,0.0,20.0,kP, 100,-5.0,5.0,kTOFnSig+kElectron);
	AddHistogram(classStr.Data(),"TOFnSigElectron_Run","TOF <n-#sigma_{e}> vs run number;run;TOF n-#sigma_{e}",kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 160, -0.5, 159.5, kTOFnSig+kElectron);
      }    // end if TOF histograms
      continue;
    }  // end if "TrackQA"
    
    Double_t massBinWidth = 0.04;     // *GeV/c^2
    Double_t massRange[2] = {0.0,6.0};
    const Int_t nMassBins = TMath::Nint((massRange[1]-massRange[0])/massBinWidth);
    Double_t massBinLims[nMassBins+1]; for(Int_t im=0; im<=nMassBins; ++im) massBinLims[im] = massBinWidth*im;
        
    // Histograms for pairs
    if(classStr.Contains("PairQA")) {
      AddHistClass(classStr.Data());
      
      AddHistogram(classStr.Data(), "Mass_CentVZEROCep", 
	           "Inv. mass vs (centVZERO, #Psi^{2}); centrality VZERO; #Psi^{2}; m (GeV/c^{2})", 
		   kFALSE, gkNCentRanges, gkEMCentLims, kCentVZERO, gkNPhiEPRanges, gkEMPhiEPLims, gkEPused, nMassBins, massBinLims, kMass);  
      AddHistogram(classStr.Data(), "Pt", "Pt; p_{T} (GeV/c)", kFALSE,
                     1000, 0.0, 10.0, kPairPt);
      
      if(classStr.Contains("SE")) {
        AddHistogram(classStr.Data(), "CandidateId", "Candidate id; p_{T} (GeV/c)", kFALSE,
                     AliCorrelationReducedPair::kNMaxCandidateTypes+1, -0.5, Double_t(AliCorrelationReducedPair::kNMaxCandidateTypes)+0.5, kCandidateId);
        AddHistogram(classStr.Data(), "PairType", "Pair type; pair type", kFALSE,
                     4, -0.5, 3.5, kPairType);	
        AddHistogram(classStr.Data(), "Mass", "Invariant mass; m_{inv} (GeV/c^{2})", kFALSE,
                     nMassBins, massRange[0], massRange[1], kMass);
        AddHistogram(classStr.Data(), "Rapidity", "Rapidity; y", kFALSE,
                     240, -1.2, 1.2, kPairRap);
	AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta; #eta", kFALSE,
                     240, -2.0, 2.0, kPairEta);
        AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution; #phi (rad.)", kFALSE,
                     315, 0.0, 6.3, kPairPhi);
	AddHistogram(classStr.Data(), "OpeningAngle", "Opening angle; op. angle (rad)", kFALSE,
                     1000, 0.0, 3.2, kPairOpeningAngle);
	AddHistogram(classStr.Data(), "Mass_Pt", "Pt vs invariant mass; m_{inv} (GeV/c^{2}); p_{T} (GeV/c)", kFALSE,
                     nMassBins, massRange[0], massRange[1], kMass, 100, 0.0, 10.0, kPairPt);
	AddHistogram(classStr.Data(), "Mass_Run_prof", "<Invariant mass> vs run; run;m_{inv} (GeV/c^{2})", kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, nMassBins, massRange[0], massRange[1], kMass);
        AddHistogram(classStr.Data(), "Pt_Run_prof", "<p_{T}> vs run; run;p_{T} (GeV/c)", kTRUE,
                     kNRunBins, runHistRange[0], runHistRange[1], kRunNo, 1000, 0.0, 10.0, kPairPt);
	AddHistogram(classStr.Data(), "VZEROflow_sideA_v2_Mass",
                     "Pair v2 coefficient using VZERO-A reaction plane; inv.mass (GeV/c^{2});v_{2}(#Psi_{2}^{VZERO-A})", kTRUE,
                     nMassBins, massRange[0], massRange[1], kMass, 1000, -1., 1., kPairVZEROFlowVn+0*6+1);
	AddHistogram(classStr.Data(), "VZEROflow_sideC_v2_Mass",
                     "Pair v2 coefficient using VZERO-C reaction plane; inv.mass (GeV/c^{2});v_{2}(#Psi_{2}^{VZERO-C})", kTRUE,
                     nMassBins, massRange[0], massRange[1], kMass, 1000, -1., 1., kPairVZEROFlowVn+1*6+1);
	AddHistogram(classStr.Data(), "TPCflow_v2_Mass",
                     "Pair v2 coefficient using TPC reaction plane; inv.mass (GeV/c^{2});v_{2}(#Psi_{2}^{TPC})", kTRUE,
                     nMassBins, massRange[0], massRange[1], kMass, 1000, -1., 1., kPairTPCFlowVn+1);
      }
      continue;
    }   // end if for Pair classes of histograms
    
  }   // end loop over histogram classes
}
