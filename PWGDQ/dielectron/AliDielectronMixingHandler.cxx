/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron MixingHandler                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TVectorD.h>
#include <TH1.h>
#include <TAxis.h>

#include <AliLog.h>
#include <AliVTrack.h>

#include "AliDielectron.h"
#include "AliDielectronHelper.h"
#include "AliDielectronHistos.h"
#include "AliDielectronEvent.h"

#include "AliDielectronMixingHandler.h"

ClassImp(AliDielectronMixingHandler)

AliDielectronMixingHandler::AliDielectronMixingHandler() :
  TNamed(),
  fDepth(10),
  fArrPools("TClonesArray"),
  fAxes(kMaxCuts),
  fMixType(kOSonly),
  fMixIncomplete(kTRUE),
  fMoveToSameVertex(kFALSE),
  fPID(0x0)
{
  //
  // Default Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fEventCuts[i]=0;
  }
  fAxes.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronMixingHandler::AliDielectronMixingHandler(const char* name, const char* title) :
  TNamed(name, title),
  fDepth(10),
  fArrPools("TClonesArray"),
  fAxes(kMaxCuts),
  fMixType(kOSonly),
  fMixIncomplete(kTRUE),
  fMoveToSameVertex(kFALSE),
  fPID(0x0)
{
  //
  // Named Constructor
  //
  for (Int_t i=0; i<kMaxCuts; ++i){
    fEventCuts[i]=0;
  }
  fAxes.SetOwner(kTRUE);
}

//______________________________________________
AliDielectronMixingHandler::~AliDielectronMixingHandler()
{
  //
  // Default Destructor
  //
  fAxes.Delete();
  delete fPID;
}

//________________________________________________________________
void AliDielectronMixingHandler::AddVariable(AliDielectronVarManager::ValueTypes type,
                                             Int_t nbins, Double_t min, Double_t max, Bool_t log)
{
  //
  // Add a variable to the mixing handler
  //

  // limit number of variables to kMaxCuts
  if (fAxes.GetEntriesFast()>=kMaxCuts) return;
  
  TVectorD *binLimits=0x0;
  if (!log) binLimits=AliDielectronHelper::MakeLinBinning(nbins,min,max);
  else binLimits=AliDielectronHelper::MakeLogBinning(nbins,min,max);
  if (!binLimits) return;

  Int_t size=fAxes.GetEntriesFast();
  fEventCuts[size]=(UShort_t)type;
  fAxes.Add(binLimits);
}

//________________________________________________________________
void AliDielectronMixingHandler::AddVariable(AliDielectronVarManager::ValueTypes type,
                                             const char* binLimitStr)
{
  //
  // Add a variable to the mixing handler with arbitrary binning
  //

  // limit number of variables to kMaxCuts
  if (fAxes.GetEntriesFast()>=kMaxCuts) return;
  
  TVectorD *binLimits=AliDielectronHelper::MakeArbitraryBinning(binLimitStr);
  if (!binLimits) return;
  
  Int_t size=fAxes.GetEntriesFast();
  fEventCuts[size]=(UShort_t)type;
  fAxes.Add(binLimits);
}

//______________________________________________
void AliDielectronMixingHandler::Fill(const AliVEvent *ev, AliDielectron *diele)
{
  //
  // fill event buffers and perform mixing if the pool depth is reached
  //

  //check if there are tracks available
  if (diele->GetTrackArray(0)->GetEntriesFast()==0 && diele->GetTrackArray(1)->GetEntriesFast()==0) return;

  TString dim;
  Int_t bin=FindBin(AliDielectronVarManager::GetData(),&dim);

  //add mixing bin to event data
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kMixingBin,bin);

  if (bin<0){
    AliDebug(5,Form("Bin outside range: %s",dim.Data()));
    return;
  }

  // get mixing pool, create it if it does not yet exist.
  TClonesArray *poolp=static_cast<TClonesArray*>(fArrPools.At(bin));
  if (!poolp){
    AliDebug(10,Form("New pool at %d (%s)\n",bin,dim.Data()));
    poolp=new(fArrPools[bin]) TClonesArray("AliDielectronEvent",fDepth);
  }
  TClonesArray &pool=*poolp;

  AliDebug(10,Form("new event at %d: %d",bin,pool.GetEntriesFast()));
  AliDielectronEvent *event=new(pool[pool.GetEntriesFast()]) AliDielectronEvent();
  if(ev->IsA() == AliAODEvent::Class()) event->SetAOD();
  else event->SetESD();

  event->SetProcessID(fPID);
  event->SetTracks(*diele->GetTrackArray(0), *diele->GetTrackArray(1), *diele->GetPairArray(1));
  event->SetEventData(AliDielectronVarManager::GetData());

  // check if pool depth is reached.
  if (pool.GetEntriesFast()<fDepth) return;

  // pool depth reached, do mixing
  DoMixing(pool,diele);

  // increase counter for full bins
  if (diele->fHistos) {
    diele->fHistos->Fill("Mixing","Stats",0);
    diele->fHistos->Fill("Mixing","CompletePools",bin);
  }
  
  //clear the current pool
  pool.Clear("C");
}

//______________________________________________
void AliDielectronMixingHandler::DoMixing(TClonesArray &pool, AliDielectron *diele)
{
  //
  // perform the mixing
  //

  //buffer track arrays and copy them back afterwards
  TObjArray arrTrDummy[4];
  for (Int_t i=0; i<4; ++i) arrTrDummy[i]=diele->fTracks[i];

  //buffer also global event data
  Double_t values[AliDielectronVarManager::kNMaxValues]={0};
  for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues; ++i)
    values[i]=AliDielectronVarManager::GetValue((AliDielectronVarManager::ValueTypes)i);

  for (Int_t i1=0; i1<pool.GetEntriesFast(); ++i1){
    AliDielectronEvent *ev1=static_cast<AliDielectronEvent*>(pool.At(i1));
    //use event data from the first event
    //both events are in the same mixing bin anyhow...
    AliDielectronVarManager::SetEventData(ev1->GetEventData());

    TObject *o=0x0;
    TIter ev1P(ev1->GetTrackArrayP());
    TIter ev1N(ev1->GetTrackArrayN());
    
    for (Int_t i2=i1+1; i2<pool.GetEntriesFast(); ++i2){
      //clear arryas
      diele->fTracks[0].Clear();
      diele->fTracks[1].Clear();
      diele->fTracks[2].Clear();
      diele->fTracks[3].Clear();
      
      //setup track arrays
      AliDielectronEvent *ev2=static_cast<AliDielectronEvent*>(pool.At(i2));
      ev1P.Reset();
      ev1N.Reset();
      TIter ev2P(ev2->GetTrackArrayP());
      TIter ev2N(ev2->GetTrackArrayN());

      //
      //move tracks to the same vertex (vertex of the first event), if requested
      //
      if (fMoveToSameVertex){
        const Double_t *varsFirst=ev1->GetEventData();
        const Double_t *varsMix=ev2->GetEventData();

        const Double_t vFirst[3]={varsFirst[AliDielectronVarManager::kXvPrim],
                                  varsFirst[AliDielectronVarManager::kYvPrim],
                                  varsFirst[AliDielectronVarManager::kZvPrim]};

        const Double_t vMix[3]  ={varsMix[AliDielectronVarManager::kXvPrim],
                                  varsMix[AliDielectronVarManager::kYvPrim],
                                  varsMix[AliDielectronVarManager::kZvPrim]};
                                  
        //loop over all tracks from the second event and move them to the vertex of the first
        AliVTrack *vtrack=0x0;
        while ( ( vtrack=(AliVTrack*)ev2P() ) ){
          MoveToSameVertex(vtrack, vFirst, vMix);
        }

        
        while ( ( vtrack=(AliVTrack*)ev2N() ) ){
          MoveToSameVertex(vtrack, vFirst, vMix);
        }

        
        ev2P.Reset();
        ev2N.Reset();
      }

      //mixing of ev1- ev2+ (pair type4). This is common for all mixing types
      while ( (o=ev1N()) ) diele->fTracks[1].Add(o); 
      while ( (o=ev2P()) ) diele->fTracks[2].Add(o);
      diele->FillPairArrays(1,2);
      
      if (fMixType==kAll || fMixType==kOSandLS){
        // all 4 pair arrays will be filled
        while ( (o=ev1P()) ) diele->fTracks[0].Add(o);
        while ( (o=ev2N()) ) diele->fTracks[3].Add(o);
        diele->FillPairArrays(0,2);
        diele->FillPairArrays(1,3);
        if (fMixType==kAll) diele->FillPairArrays(0,3);
      }
      
      if (fMixType==kOSonly || fMixType==kOSandLS){
        //use the pair type of ev1- ev1+ also for ev1+ ev1-
        diele->fTracks[1].Clear();
        diele->fTracks[2].Clear();
        while ( (o=ev1P()) ) diele->fTracks[1].Add(o);
        while ( (o=ev2N()) ) diele->fTracks[2].Add(o);
        diele->FillPairArrays(1,2);
      }

    }
  }

  //copy back the tracks
  for (Int_t i=0; i<4; ++i) {
    diele->fTracks[i].Clear();
    diele->fTracks[i]=arrTrDummy[i];
  }

  //set back global event values
  AliDielectronVarManager::SetEventData(values);
}

//______________________________________________
void AliDielectronMixingHandler::MixRemaining(AliDielectron *diele)
{
  //
  // mix all pools even if they are incomplete
  //

  //Check if there was any processed data and it is requested to mix incomplete bins
  if (!diele || !diele->PairArray(0) || !fMixIncomplete ) return;

  AliDielectronVarManager::SetEvent(0x0);
  for (Int_t ipool=0; ipool<fArrPools.GetSize(); ++ipool){
    TClonesArray *poolp=static_cast<TClonesArray*>(fArrPools.At(ipool));
    if (!poolp || !poolp->GetEntriesFast() || !poolp->At(0)) continue;
    //clear the arrays before the final processing"
    AliDebug(10,Form("Incomplete: Bin %d (%d)\n",ipool,poolp->GetEntriesFast()));
    diele->ClearArrays();
    DoMixing(*poolp,diele);

    // increase counter for incomplete bins
    if (diele->fHistos) {
      //buffer event data and set event data using the first event in this pool
      Double_t values[AliDielectronVarManager::kNMaxValues]={0};
      for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues; ++i)
        values[i]=AliDielectronVarManager::GetValue((AliDielectronVarManager::ValueTypes)i);
      
      AliDielectronEvent *ev1=static_cast<AliDielectronEvent*>(poolp->At(0));
      //use event data from the first event all events are in the same mixing bin anyhow...
      AliDielectronVarManager::SetEventData(ev1->GetEventData());
      
      // fill the histograms
      diele->FillHistograms(0x0, kTRUE);
      diele->fHistos->Fill("Mixing","Stats",1);
      diele->fHistos->Fill("Mixing","InCompletePools",ipool);
      diele->fHistos->Fill("Mixing","Entries_InCompletePools",poolp->GetEntriesFast());
      
      //set back global event values
      AliDielectronVarManager::SetEventData(values);
    }
    
  }
}


//______________________________________________
void AliDielectronMixingHandler::Init(const AliDielectron *diele)
{
  //
  // initialise event buffers
  //

  Int_t size=GetNumberOfBins();

  AliDebug(10,Form("Creating a pool array with size %d \n",size));

  fArrPools.Expand(size);

  //add statics histogram if we have a histogram manager
  if (diele && diele->fHistos) {
    diele->fHistos->AddClass("Mixing");
    diele->fHistos->UserHistogram("Mixing","Stats","Mixing Statistics;;#called bins",2,0,2);
    TH1* h=diele->fHistos->GetHistogram("Mixing","Stats");
    h->GetXaxis()->SetBinLabel(1,"Complete");
    h->GetXaxis()->SetBinLabel(2,"Incomplete");

    diele->fHistos->UserHistogram("Mixing","CompletePools","Mixing Statistics compete pools;bin;#fills",size,0,size);
    diele->fHistos->UserHistogram("Mixing","InCompletePools","Mixing Statistics incomplete pools;bin;#fills",size,0,size);
    diele->fHistos->UserHistogram("Mixing","Entries_InCompletePools","#entries in incomplete pools;entries;#fills",fDepth,0,fDepth);
  }

  TString values;
  for (Int_t i=0; i<fAxes.GetEntriesFast(); ++i){
    TVectorD *bins=static_cast<TVectorD*>(fAxes.At(i));
    Int_t nRows=bins->GetNrows();
    values+=Form("%s: ",AliDielectronVarManager::GetValueName(fEventCuts[i]));
    for (Int_t irow=0; irow<nRows; ++irow){
      values+=Form("%.2f, ",(*bins)[irow]);
    }
  }
  
  if (!fPID){
    fPID=TProcessID::AddProcessID();
  }

  AliDebug(10,values.Data());
}

//______________________________________________
Int_t AliDielectronMixingHandler::GetNumberOfBins() const
{
  //
  // return the number of bins this mixing handler has
  //
  Int_t size=1;
  for (Int_t i=0; i<fAxes.GetEntriesFast(); ++i)
    size*=((static_cast<TVectorD*>(fAxes.At(i)))->GetNrows()-1);
  return size;
}

//______________________________________________
Int_t AliDielectronMixingHandler::FindBin(const Double_t values[], TString *dim)
{
  //
  // bin bin in mixing stack described by 'values'
  // if the values are outside the binning range -1 is returned
  // if dim is non NULL debug info will be stored in the variable
  //

  if (fAxes.GetEntriesFast()==0) {
    if (dim) (*dim)="single bin";
    return 0;
  }
  if (dim) (*dim)="";
  Int_t sizeAdd=1;
  Int_t bin=0;
  for (Int_t i=0; i<fAxes.GetEntriesFast(); ++i){
    Double_t val=values[fEventCuts[i]];
    TVectorD *bins=static_cast<TVectorD*>(fAxes.At(i));
    Int_t nRows=bins->GetNrows();
    if ( (val<(*bins)[0]) || (val>(*bins)[nRows-1]) ) {
      return -1;
    }

    Int_t pos=TMath::BinarySearch(nRows,bins->GetMatrixArray(),val);
    bin+=sizeAdd*pos;
    if (dim) (*dim)+=Form("%s: %f (%d); ",AliDielectronVarManager::GetValueName(fEventCuts[i]),val,pos);
    sizeAdd*=(nRows-1);
  }

  return bin;
}

//______________________________________________
void AliDielectronMixingHandler::MoveToSameVertex(AliVTrack * const vtrack, const Double_t *vFirst, const Double_t* vMix)
{
  //
  // move 'track' which belongs to the vertex information of vMix to the vertex of vFirst
  //

  static Bool_t printed=kFALSE;
  
  if (vtrack->IsA()==AliESDtrack::Class()){
    AliESDtrack *track=(AliESDtrack*)vtrack;

    //get track information
    Double_t x        = track->GetX();
    Double_t alpha    = track->GetAlpha();
    Double_t param[5] = {0};
    Double_t cov[15]  = {0};

    for (Int_t i=0; i<5;  ++i) param[i]=track->GetParameter()[i];
    for (Int_t i=0; i<15; ++i) cov[i]  =track->GetCovariance()[i];

    //translation
    Double_t vt[3] = {vMix[0]-vFirst[0],vMix[1]-vFirst[1],vMix[2]-vFirst[2]};
    //rotate to the track frame
//     track->Global2LocalPosition(vt,track->GetAlpha());

    //add to track position
//     x        = x       -vt[0];
//     param[0] = param[0]-vt[1];
//     param[1] = param[1]-vt[2];
    param[1] = param[1]-vt[2];
    
    //set updated track information
    track->Set(x, alpha, param, cov);
  } else {
    //             AliAODTrack *track=(AliAODTrack*)vtrack;
    //             Double_t pos[3]={0};
    //             track->GetPosition(pos);
    //             if (pos[0]>-999.){
      //               pos[0]=pos[0]-vMix[AliDielectronVarManager::kXvPrim]+vFirst[AliDielectronVarManager::kXvPrim];
      //               pos[1]=pos[1]-vMix[AliDielectronVarManager::kYvPrim]+vFirst[AliDielectronVarManager::kYvPrim];
      //               pos[2]=pos[2]-vMix[AliDielectronVarManager::kZvPrim]+vFirst[AliDielectronVarManager::kZvPrim];
      //               track->SetPosition(pos);
//       AliError("Move To same vertex not yet implemented for AOD!");
    if (!printed) {
//       Error("AliDielectronMixingHandler::MoveToSameVertex","Move To same vertex not yet implemented for AOD!");
      printed=kTRUE;
    }
      //
      // Not that clear how to do it. In AOD track there is fPosition, fPositionAtDCA and "TRef fProdVertex"
      // where Xv(), Yv(), Zv() returns the position of fProdVertex
      //
    }
    
}
