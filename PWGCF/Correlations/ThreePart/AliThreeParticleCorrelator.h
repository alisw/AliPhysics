//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliThreeParticleCorrelator.h
/// @author Matthias Richter 
/// @author Paul Baetzing 		pbatzing@cern.ch
/// @date   2013-03-14
/// @brief  Three particle correlation steering class
///

#ifndef ALITHREEPARTICLECORRELATOR_H
#define ALITHREEPARTICLECORRELATOR_H

#include "TNamed.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliFilteredTrack.h"
#include "AliEventPoolManager.h"
#include "TObjArray.h"
#include <vector>
#include <memory>
#include <cerrno>
#include <AliAODMCParticle.h>
#include "AliCFPI0.h"
#include <TVectorT.h>
#include <AliMCEventHandler.h>
#include <AliAnalysisManager.h>
#include <TRandom3.h>
#include "THn.h"
#include "TTimeStamp.h"


using std::vector;

class TMethod;
template <class C>
class AliThreeParticleCorrelator : public TNamed {
 public:
  /// default constructor   
  AliThreeParticleCorrelator(const char* name=NULL)
    : TNamed(name?name:"AliThreeParticleCorrelator", "")
    , fCorrelations()
    , fMECorrelations()
    , fMETriggerCorrelations()
    , fMETACorrelations()
    , fMETA2Correlations()
    , factiveTriggers()
    , fAssociated()
    , fAssociatedmixed1()
    , fAssociatedmixed2()
    , fEventPoolMgr(NULL)
    , fVz(0)
    , fMultiplicity(0)
    , fRandom()
    , fMaxMixedPerEvent(-1)
    , fLeading(kTRUE)
{  
    fRandom = new TRandom3();
    TTimeStamp now;
    fRandom->SetSeed(now.GetNanoSec());
    delete gRandom;
    gRandom = fRandom;
  }
  AliThreeParticleCorrelator(const AliThreeParticleCorrelator& other)
    : TNamed(other)
    , fCorrelations(other.fCorrelations)
    , fMECorrelations(other.fMECorrelations)
    , fMETriggerCorrelations(other.fMETriggerCorrelations)
    , fMETACorrelations(other.fMETACorrelations)
    , fMETA2Correlations(other.fMETA2Correlations)
    , factiveTriggers(other.factiveTriggers)
    , fAssociated(other.fAssociated)
    , fAssociatedmixed1(other.fAssociatedmixed1)
    , fAssociatedmixed2(other.fAssociatedmixed2)
    , fEventPoolMgr(other.fEventPoolMgr)
    , fVz(other.fVz)
    , fMultiplicity(other.fMultiplicity)
    , fRandom(other.fRandom)
    , fMaxMixedPerEvent(other.fMaxMixedPerEvent)
    , fLeading(kTRUE)
  {
    if(fRandom)delete fRandom;
    fRandom = new TRandom3();
    TTimeStamp now;
    fRandom->SetSeed(now.GetNanoSec());
    delete gRandom;
    gRandom = fRandom;  
  }

  /// destructor
  virtual ~AliThreeParticleCorrelator() {if(fRandom){delete fRandom;}}

  AliThreeParticleCorrelator& operator=(const AliThreeParticleCorrelator& other){
    // assignment operator
    if (this==&other) return *this;
    this->~AliThreeParticleCorrelator();
    new (this) AliThreeParticleCorrelator(other);
    return *this;
  }

  /// add an analysis object, if ME analysis is enabled, a corresponding clone
  /// is automatically created
  C* Add(C* o, bool extramixed = true) {
    if (o) {
      fCorrelations = o;
      if (fEventPoolMgr) {
	fMECorrelations = new C(*o);
	TString name=fMECorrelations->GetName();
	name+="MEAll";
	fMECorrelations->SetName(name);
	fMECorrelations->Init();
	if(extramixed){
	  // create a clone for ME analysis
	  fMETriggerCorrelations=new C(*o);
	  name=fMETriggerCorrelations->GetName();
	  name+="METrigger";
	  fMETriggerCorrelations->SetName(name);
	  fMETriggerCorrelations->Init();
	  fMETACorrelations=new C(*o);
	  name=fMETACorrelations->GetName();
	  name+="META";
	  fMETACorrelations->SetName(name);
	  fMETACorrelations->Init();
	  fMETA2Correlations=new C(*o);
	  name=fMETA2Correlations->GetName();
	  name+="META2";
	  fMETA2Correlations->SetName(name);
	  fMETA2Correlations->Init();	  
	  fMETriggerCorrelations->SetMixedEvent(fMECorrelations);
	  fMETACorrelations->SetMixedEvent(fMECorrelations);
	  fMETA2Correlations->SetMixedEvent(fMECorrelations);
	}
	o->SetMixedEvent(fMECorrelations);
      }
    }
    return fCorrelations;
  }

  /// init the event pool manager for event mixing
  void InitEventMixing(AliEventPoolManager* mgr) {fEventPoolMgr=mgr;}
  
  /// Set the event vertex and centrality.
  void SetEventVzM(Double_t vz, Double_t m){
    fVz = vz;
    fMultiplicity = m;
  }
  
  void SetLeading(bool isleading){fLeading = isleading;}
  
  
  /// Set the maximum number of tracks kept in the mixed event pool
  void SetNMaxMixed(Int_t i){fMaxMixedPerEvent = i;}
  
  /// get ME analysis object corresponding to parameter
  C* GetCorrespondingME(C* o, int type=0) {
    if (type==0 && fMECorrelations) return fMECorrelations;
    if (type==1 && fMETriggerCorrelations) return fMETriggerCorrelations;
    if (type==2 && fMETACorrelations) return fMETACorrelations;
    if (type==3 && fMETA2Correlations) return fMETA2Correlations;
    return NULL;
  }

  //Process the Event
  virtual void Execute(TMethod */*method*/, TObjArray *arrayParticles,Int_t *error=0) {
    int err=ProcessEvent(arrayParticles);
    if (error) *error=err;
  }
  using TObject::Execute;

  int ProcessEvent(TObjArray* arrayParticles) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    //Propagate the event properties to all workers and initiate monte carlo if applicable
    if(fCorrelations->SetMultVZ(fMultiplicity,fVz)<0) return -EINVAL;
    if (!arrayParticles) return -EINVAL;
    //Create a new TObjArray that is owner to be able to cherry pick a subset of the tracks.
    //To fill with tracks and pions:
    TObjArray* Mixedparticles = NULL;
    Mixedparticles = new TObjArray();
    Mixedparticles->SetOwner();//In order for it to work in the event pool.
    
    //And for all ME workers.
    if(fEventPoolMgr)fMECorrelations->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)fMETriggerCorrelations->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)fMETACorrelations->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)fMETA2Correlations->SetMultVZ(fMultiplicity,fVz);
    //Fill the associated vector
    MakeAssociated(arrayParticles,Mixedparticles, fAssociated,true);
    //Fill the vector for all objects in fCorrelations.
    MakeTriggers(arrayParticles, Mixedparticles,fLeading);

    //Process the signal 3p and 2p correlations for all workers.
    int iResult=ProcessEvent(factiveTriggers, fAssociated,fCorrelations);
    //Event Mixing
    if (fEventPoolMgr) {
      AliEventPool* pool=fEventPoolMgr->GetEventPool(fMultiplicity, fVz);
      int EventsInPool = pool->GetCurrentNEvents();

      if (pool) {
	if (pool->IsReady()) {
	  //Correlate triggers from this event with past associated:
	  MakeAssociated(Mixedparticles, fAssociated,fMECorrelations,false);
	  MakeTriggers(Mixedparticles,fMECorrelations,fLeading,false);
	  for (int nEvent=0; nEvent<EventsInPool; nEvent++) {
	    MakeAssociated(pool->GetEvent(nEvent), fAssociatedmixed1,fMECorrelations);
	    ProcessEvent(factiveTriggers,fAssociated,fAssociatedmixed1,fMETACorrelations);
	    ProcessEvent(factiveTriggers,fAssociatedmixed1,fAssociated,fMETA2Correlations);
	    ProcessEvent(factiveTriggers, fAssociatedmixed1,fMETriggerCorrelations);
	    for(int mEvent=nEvent+1;mEvent<EventsInPool;mEvent++){
	      // all particles from different events
	      MakeAssociated(pool->GetEvent(mEvent), fAssociatedmixed2,fMECorrelations);
	      ProcessEvent(factiveTriggers,fAssociatedmixed1,fAssociatedmixed2,fMECorrelations);
	      ProcessEvent(factiveTriggers,fAssociatedmixed2,fAssociatedmixed1,fMECorrelations,false);
	    }//end inner loop
	  }//End mixed event loop
	  //Correlate associated from this event with triggers and associated from past events:
	  for (int nEvent=0; nEvent<EventsInPool; nEvent++){
	    MakeTriggers(pool->GetEvent(nEvent),fMECorrelations,fLeading);
	    MakeAssociated(pool->GetEvent(nEvent), fAssociatedmixed1,fMECorrelations);
	    //Now factiveTriggers is the same event as fAssociatedmixed1.
	    ProcessEvent(factiveTriggers,fAssociatedmixed1,fAssociated,fMETACorrelations);
	    ProcessEvent(factiveTriggers,fAssociated,fAssociatedmixed1,fMETA2Correlations);
	    ProcessEvent(factiveTriggers, fAssociated,fMETriggerCorrelations);
	    for(int mEvent=0;mEvent<EventsInPool;mEvent++){
	      if(nEvent==mEvent)continue;//dont double count
	      // all particles from different events
	      MakeAssociated(pool->GetEvent(mEvent), fAssociatedmixed2,fMECorrelations);
	      ProcessEvent(factiveTriggers,fAssociated,fAssociatedmixed2,fMECorrelations);
	      ProcessEvent(factiveTriggers,fAssociatedmixed2,fAssociated,fMECorrelations,false);
	    }//end inner loop
	  }//End mixed event loop
	}//End poolisready
	if (Mixedparticles->GetEntriesFast()>0) {//update event pool
	  EventsInPool=pool->UpdatePool(Mixedparticles);
	}//End update event Pool
      }//End ifPOOl
      
    }
    factiveTriggers.clear();
    fAssociated.clear();
    fAssociatedmixed1.clear();
    fAssociatedmixed2.clear();
    delete arrayParticles;
    return iResult;
  }
  int ProcessEvent(const std::vector<AliVParticle*>& activeTriggers, const std::vector<AliVParticle*>& associated,C* AnalysisObject) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    Double_t NAssociated = associated.size();
    Double_t weightt = 1.0;
    Double_t weighta1 = 1.0;
    Double_t weighta2 = 1.0;
    if(NAssociated==0) return 0;//No associated means we need not fill anything.
    if (activeTriggers.size()==0) return 0;//No Triggers means we need not fill anything
    for (typename std::vector<AliVParticle*>::const_iterator trigger=activeTriggers.begin(), e=activeTriggers.end(); trigger!=e; ++trigger) {
      AnalysisObject->FillTrigger(*trigger);//Fill histogram for number of triggers.
      weightt = 1.0;
      if(dynamic_cast<AliFilteredTrack*>(*trigger))weightt = dynamic_cast<AliFilteredTrack*>(*trigger)->GetEff();
      for (typename std::vector<AliVParticle*>::const_iterator assoc=associated.begin(), eassoc=associated.end(); assoc!=eassoc; ++assoc) {
	weighta1 = 1.0;
	if(dynamic_cast<AliFilteredTrack*>(*assoc))weighta1 = dynamic_cast<AliFilteredTrack*>(*assoc)->GetEff();
	for (typename std::vector<AliVParticle*>::const_iterator assoc2=assoc+1, eassoc2 = associated.end(); assoc2 !=eassoc2;++assoc2){
	  weighta2 = 1.0;
	  if(dynamic_cast<AliFilteredTrack*>(*assoc2))weighta2 = dynamic_cast<AliFilteredTrack*>(*assoc2)->GetEff();
	  AnalysisObject->Fill(*trigger,*assoc,*assoc2,weightt*weighta1*weighta2);//Fill histogram for number of triggers.
	  AnalysisObject->Fill(*trigger,*assoc2,*assoc,weightt*weighta1*weighta2);//Fill histogram for number of triggers.
	  if(trigger==activeTriggers.begin()){
	    //once per event fill the a-a 2p correlation histogram symmitrized
	    AnalysisObject->Filla(*assoc,*assoc2,weighta1*weighta2);
	    AnalysisObject->Filla(*assoc2,*assoc,weighta1*weighta2);
	  }
	}//loop over second associated
	AnalysisObject->Fill(*trigger,*assoc,weightt*weighta1);//Fill histogram for number of triggers.	        
      } // loop over first associated
    } // loop over triggers
    return 0;
  }
  int ProcessEvent(const std::vector<AliVParticle*>& activeTriggers,const std::vector<AliVParticle*>& associated, const std::vector<AliVParticle*>& associatedmixed,C* AnalysisObject,bool twop = true) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    Double_t NAssociated1 = associated.size();
    Double_t NAssociated2 = associatedmixed.size();
    Double_t weightt = 1.0;
    Double_t weighta1 = 1.0;
    Double_t weighta2 = 1.0;
    if(NAssociated1==0||NAssociated2==0) return 0;//No associated means we need not fill anything.
    if (activeTriggers.size()==0) return 0;//no triggers means nothing to be correlated
    for (typename std::vector<AliVParticle*>::const_iterator trigger=activeTriggers.begin(), e=activeTriggers.end(); trigger!=e; ++trigger) {
      AnalysisObject->FillTrigger(*trigger);//Fill histogram for number of triggers.
      weightt = 1.0;
      if(dynamic_cast<AliFilteredTrack*>(*trigger))weightt = dynamic_cast<AliFilteredTrack*>(*trigger)->GetEff();
      for (typename std::vector<AliVParticle*>::const_iterator assoc=associated.begin(), eassoc=associated.end(); assoc!=eassoc; ++assoc) {
	weighta1 = 1.0;
	if(dynamic_cast<AliFilteredTrack*>(*assoc))weighta1 = dynamic_cast<AliFilteredTrack*>(*assoc)->GetEff();
	for (typename std::vector<AliVParticle*>::const_iterator assoc2=associatedmixed.begin(), eassoc2 = associatedmixed.end(); assoc2 !=eassoc2;++assoc2){
	  weighta2 = 1.0;
	  if(dynamic_cast<AliFilteredTrack*>(*assoc2))weighta2 = dynamic_cast<AliFilteredTrack*>(*assoc2)->GetEff();
	  AnalysisObject->Fill(*trigger,*assoc,*assoc2,weightt*weighta1*weighta2);//Fill histogram for number of triggers.
	  if(trigger==activeTriggers.begin()){
	    //once per event fill the a-a 2p correlation histogram
	    AnalysisObject->Filla(*assoc,*assoc2,weighta1*weighta2);
	    AnalysisObject->Filla(*assoc2,*assoc,weighta1*weighta2);
	  }
	}//loop over second associated
	if (twop){
	  AnalysisObject->Fill(*trigger,*assoc,weightt*weighta1);//Fill histogram for number of triggers.  
	  }
	} // loop over first associated
    } // loop over triggers
    return 0;
  }
  
  void MakeTriggers(const TObjArray* arrayParticles,C* fAnalysisObject,bool fLeading,bool makehist=kFALSE) {
    /// create a particle array with reduced data objects
    /// for trigger particles containing only potential triggers.
    //if fLeading is true, it returns only the leading track in pT
    if (!arrayParticles) return ;
    factiveTriggers.clear();
    TIter nexttrigger(arrayParticles);
    TObject* otrigger=NULL;
    while ((otrigger=nexttrigger())!=NULL) {//loop over all in the array.
      AliVParticle* ptrigger=reinterpret_cast<AliVParticle*>(otrigger);
      if (!ptrigger) continue;
      if (fLeading&&(dynamic_cast<AliFilteredTrack*>(ptrigger)))if(!dynamic_cast<AliFilteredTrack*>(ptrigger)->IsLeading())continue;
      if(!fAnalysisObject->CheckTrigger(ptrigger, makehist)) continue;//check for each correlation, if it is a trigger in it. If not, reject.
      factiveTriggers.push_back(ptrigger);
    }//end while loop
    return ;
  }
  void MakeTriggers(const TObjArray* arrayParticles, TObjArray* outarrayParticles,bool fLeading) {
    /// create a particle array with reduced data objects
    /// for trigger particles containing only potential triggers.
    //if fLeading is true, it returns only the leading track in pT
    if (!arrayParticles||!outarrayParticles) return ;
    factiveTriggers.clear();
    TIter nexttrigger(arrayParticles);
    TObject* otrigger=NULL;
    while ((otrigger=nexttrigger())!=NULL) {//loop over all in the array.
      AliVParticle* ptrigger=reinterpret_cast<AliVParticle*>(otrigger);
      if (!ptrigger) continue;
      if (fLeading&&(dynamic_cast<AliFilteredTrack*>(ptrigger)))if(!dynamic_cast<AliFilteredTrack*>(ptrigger)->IsLeading())continue;
      if (!fCorrelations->CheckTrigger(ptrigger, true)) continue;//check for each correlation, if it is a trigger in it. If not, reject.
      factiveTriggers.push_back(ptrigger);
      AliFilteredTrack * outtrigger = new AliFilteredTrack(*dynamic_cast<AliFilteredTrack*>(ptrigger));
      outtrigger->SetLeading();
      outarrayParticles->Add(outtrigger);
    }//end while loop
    return ;

  }
  void MakeAssociated(const TObjArray* arrayParticles,std::vector<AliVParticle*>& associated, C* fAnalysisObject, bool makehist=kFALSE) {
    /// Create clone of particle array using a reduced track class and
    if (!arrayParticles) return ;
    associated.clear();
    TIter iassociated(arrayParticles);
    TObject* o1=NULL;
    while ((o1=iassociated())!=NULL) {
      // loop over first particle
      AliVParticle* associatedp=dynamic_cast<AliVParticle*>(o1);
      if (!associatedp) continue;
      if (!fAnalysisObject->CheckAssociated(associatedp, makehist)) continue;
      associated.push_back(associatedp);	
    } // loop over particle
    return ;
  }
  void MakeAssociated(const TObjArray* arrayParticles,TObjArray* outarrayParticles,std::vector<AliVParticle*>& associated, bool makehist=kFALSE) {
    /// Create clone of particle array using a reduced track class and
    if (!arrayParticles||!outarrayParticles) return ;
    associated.clear();
    TIter iassociated(arrayParticles);
    TObject* o1=NULL;
    while ((o1=iassociated())!=NULL) {
      // loop over first particle
      AliVParticle* associatedp=dynamic_cast<AliVParticle*>(o1);
      if (!associatedp) continue;
      if (!fCorrelations->CheckAssociated(associatedp, makehist)) continue;
      associated.push_back(associatedp);	
    } // loop over particle
    double nass = associated.size();
    Double_t fraction = 1.1;
    if(nass>fMaxMixedPerEvent&&(fMaxMixedPerEvent>0))fraction = fMaxMixedPerEvent/nass;
    for (typename std::vector<AliVParticle*>::const_iterator assoc=associated.begin(), eassoc=associated.end(); assoc!=eassoc; assoc++) {
      if(gRandom->Rndm()>fraction)continue;
      AliFilteredTrack * outass = new AliFilteredTrack(*dynamic_cast<AliFilteredTrack*>(*assoc));
      outarrayParticles->Add(outass);
    }
    return ;
  }
  void Clear(Option_t * /*option*/)
  {
    /// overloaded from TObject: cleanup
    return TObject::Clear();
  }
  void Print(Option_t */*option*/) const
  {
    /// overloaded from TObject: print info
    TObject::Print();
  }
 private:
  C* fCorrelations; //! worker object
  C* fMECorrelations; //! ME analysis clones of worker, all three particles different events, only one associated from each event
  C* fMETriggerCorrelations; //! ME analysis clones of worker, associated particles same events
  C* fMETACorrelations; //! ME analysis clones of worker, first associated particle from the same events as the trigger
  C* fMETA2Correlations; //! ME analysis clones of worker, second associated particle from the same events as the trigger
  std::vector<AliVParticle*> factiveTriggers; //!Vector to contain the triggers.
  std::vector<AliVParticle*> fAssociated;//!vector to contain the associated particles.
  std::vector<AliVParticle*> fAssociatedmixed1; //! buffer for Associated
  std::vector<AliVParticle*> fAssociatedmixed2; //! buffer for Associated
  AliEventPoolManager* fEventPoolMgr; //! event pool manager, external pointer
  Double_t fVz;//Vertex in z
  Double_t fMultiplicity;//Multiplicity in %
  TRandom3 * fRandom;//!to be able to pick mixed event tracks.
  Double_t fMaxMixedPerEvent;//limit on the size of each event in the pool
  bool     fLeading;//if true only the leading pT trigger is kept.
 ClassDef(AliThreeParticleCorrelator, 4)
};

#endif

