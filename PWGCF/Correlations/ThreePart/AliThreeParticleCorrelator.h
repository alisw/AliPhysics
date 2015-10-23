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
    , factiveTriggersME()
    , factiveTriggersMETrigger()
    , factiveTriggersMETA()
    , factiveTriggersMETA2()
    , fAssociated()
    , fAssociatedmixed1()
    , fAssociatedmixed2()
    , fEventPoolMgr(NULL)
    , fVz(0)
    , fMultiplicity(0)
    , fRandom()
{  
//     fRandom = new TRandom3();
//     fRandom->SetSeed(0);
//     delete gRandom;
//     gRandom = fRandom;
  }
  AliThreeParticleCorrelator(const AliThreeParticleCorrelator& other)
    : TNamed(other)
    , fCorrelations(other.fCorrelations)
    , fMECorrelations(other.fMECorrelations)
    , fMETriggerCorrelations(other.fMETriggerCorrelations)
    , fMETACorrelations(other.fMETACorrelations)
    , fMETA2Correlations(other.fMETA2Correlations)
    , factiveTriggers(other.factiveTriggers)
    , factiveTriggersME(other.factiveTriggersME)
    , factiveTriggersMETrigger(other.factiveTriggersMETrigger)
    , factiveTriggersMETA(other.factiveTriggersMETA)
    , factiveTriggersMETA2(other.factiveTriggersMETA2)
    , fAssociated(other.fAssociated)
    , fAssociatedmixed1(other.fAssociatedmixed1)
    , fAssociatedmixed2(other.fAssociatedmixed2)
    , fEventPoolMgr(other.fEventPoolMgr)
    , fVz(other.fVz)
    , fMultiplicity(other.fMultiplicity)
    , fRandom(other.fRandom)
  {
/*    delete gRandom;
    gRandom = fRandom;  */  
  }

  /// destructor
  virtual ~AliThreeParticleCorrelator() {delete fRandom;}

  AliThreeParticleCorrelator& operator=(const AliThreeParticleCorrelator& other){
    // assignment operator
    if (this==&other) return *this;
    this->~AliThreeParticleCorrelator();
    new (this) AliThreeParticleCorrelator(other);
    return *this;
  }

      //Class for trigger objects.
  class AliActiveTrigger {
  public:
    AliActiveTrigger(): fCorrelation(NULL), fTrigger(NULL){}
    AliActiveTrigger(C* worker, AliVParticle* trigger) : fCorrelation(worker) , fTrigger(trigger) {}
    AliActiveTrigger(const AliActiveTrigger& other): fCorrelation(other.fCorrelation) , fTrigger(other.fTrigger) {}
    ~AliActiveTrigger() {}
    AliActiveTrigger& operator=(const AliActiveTrigger& other){
      // assignment operator
      if (this==&other) return *this;
      this->~AliActiveTrigger();
      new (this) AliActiveTrigger(other);
      return *this;
    }

    //check if the associated fits the cuts for the given worker.
    bool CheckAssociated(AliVParticle* p, bool doHistogram=false) const {
      // check the conditions for an associated particle
      if (!fCorrelation) return false;
      if (!fTrigger) return false;
      return fCorrelation->CheckAssociated(p, fTrigger, doHistogram);
    }

    int Incrementtrigger() const {
      return fCorrelation->FillTrigger(fTrigger);
    }
    //3p correlation
    int Correlate(AliVParticle* p1, AliVParticle* p2, int weight) const {
      // fill
      if (!fCorrelation) return false;
      if (!fTrigger) return false;
      return fCorrelation->Fill(fTrigger, p1, p2, weight);
    }
      
    //2p correlation
    int Correlate(AliVParticle* p1) const {
      // fill 2p correlation
      if (!fCorrelation) return false;
      if (!fTrigger) return false;
      return fCorrelation->Fill(fTrigger, p1);
    }

    //to be able to compare to associated particles.
    AliVParticle* GetTrigger() const {return fTrigger;}
  private:
    C*            fCorrelation; //! correlation instance
    AliVParticle* fTrigger;     //! trigger particle

  };  
  
  /// add an analysis object, if ME analysis is enabled, a corresponding clone
  /// is automatically created
  C* Add(C* o) {
    if (o) {
      fCorrelations.push_back(o);
      if (fEventPoolMgr) {
	// create a clone for ME analysis
	fMETriggerCorrelations.push_back(new C(*o));
	TString name=fMETriggerCorrelations.back()->GetName();
	name+="METrigger";
	fMETriggerCorrelations.back()->SetName(name);
	fMETriggerCorrelations.back()->Init();
	fMETACorrelations.push_back(new C(*o));
	name=fMETACorrelations.back()->GetName();
	name+="META";
	fMETACorrelations.back()->SetName(name);
	fMETACorrelations.back()->Init();
	fMECorrelations.push_back(new C(*o));
	name=fMECorrelations.back()->GetName();
	name+="MEAll";
	fMECorrelations.back()->SetName(name);
	fMECorrelations.back()->Init();
	
	fMETA2Correlations.push_back(new C(*o));
	name=fMETA2Correlations.back()->GetName();
	name+="META2";
	fMETA2Correlations.back()->SetName(name);
	fMETA2Correlations.back()->Init();
	
	o->SetMixedEvent(fMECorrelations.back());
	fMETriggerCorrelations.back()->SetMixedEvent(fMECorrelations.back());
	fMETACorrelations.back()->SetMixedEvent(fMECorrelations.back());
	fMETA2Correlations.back()->SetMixedEvent(fMECorrelations.back());
	
      }
    }
    return fCorrelations.back();
  }

  /// init the event pool manager for event mixing
  void InitEventMixing(AliEventPoolManager* mgr) {fEventPoolMgr=mgr;}
  
  /// Set the event vertex and centrality.
  void SetEventVzM(Double_t vz, Double_t m){
    fVz = vz;
    fMultiplicity = m;
  }
  
  /// Set the efficiency weights.
  void SetWeights(TObject* weights, int type){
    for (typename vector<C*>::iterator o=fCorrelations.begin(), e=fCorrelations.end();o!=e; o++)  {
      (*o)->SetWeight(weights,type);
      //And for all ME workers.
      if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMECorrelations.begin(), e=fMECorrelations.end();o!=e; o++)(*o)->SetWeight(weights,type);
      if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETriggerCorrelations.begin(), e=fMETriggerCorrelations.end();o!=e; o++)(*o)->SetWeight(weights,type);
      if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETACorrelations.begin(), e=fMETACorrelations.end();o!=e; o++)(*o)->SetWeight(weights,type);
      if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETA2Correlations.begin(), e=fMETA2Correlations.end();o!=e; o++)(*o)->SetWeight(weights,type);
    }
  }
    
  /// get ME analysis object corresponding to parameter
  C* GetCorrespondingME(C* o, int type=0) {
    for (unsigned i=0, e=fCorrelations.size(); i<e; i++) {
      if (fCorrelations[i]==o) {
	if (type==0 && fMECorrelations.size()>i) return fMECorrelations[i];
	if (type==1 && fMETriggerCorrelations.size()>i) return fMETriggerCorrelations[i];
	if (type==2 && fMETACorrelations.size()>i) return fMETACorrelations[i];
	if (type==3 && fMETA2Correlations.size()>i) return fMETA2Correlations[i];

	break;
      }
    }
    return NULL;
  }

  //Process the Event
  virtual void Execute(TMethod */*method*/, TObjArray *arrayParticles, Int_t *error=0) {
    int err=ProcessEvent(arrayParticles);
    if (error) *error=err;
  }
  using TObject::Execute;

  int ProcessEvent(TObjArray* arrayParticles) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    //Propagate the event properties to all workers and initiate monte carlo if applicable
    for (typename vector<C*>::iterator o=fCorrelations.begin(), e=fCorrelations.end();o!=e; o++)  {
      if((*o)->SetMultVZ(fMultiplicity,fVz)<0) return -EINVAL;
    }
    
    if (!arrayParticles) return -EINVAL;

    //And for all ME workers.
    if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMECorrelations.begin(), e=fMECorrelations.end();o!=e; o++)(*o)->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETriggerCorrelations.begin(), e=fMETriggerCorrelations.end();o!=e; o++)(*o)->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETACorrelations.begin(), e=fMETACorrelations.end();o!=e; o++)(*o)->SetMultVZ(fMultiplicity,fVz);
    if(fEventPoolMgr)for (typename vector<C*>::iterator o=fMETA2Correlations.begin(), e=fMETA2Correlations.end();o!=e; o++)(*o)->SetMultVZ(fMultiplicity,fVz);
    
    //Fill the vector for all objects in fCorrelations.
    MakeTriggers(arrayParticles, factiveTriggers,fCorrelations.begin(),fCorrelations.end());

    //Fill the associated vector
//     TObjArray* associatedTracks=
    MakeAssociated(arrayParticles, fAssociated);
//     if (!associatedTracks) return -1;//The function failed somehow.
    //Process the signal 3p and 2p correlations for all workers.
    int iResult=ProcessEvent(factiveTriggers, fAssociated);

    //Event Mixing
    if (fEventPoolMgr) {
      AliEventPool* pool=fEventPoolMgr->GetEventPool(fMultiplicity, fVz);
      int EventsInPool = pool->GetCurrentNEvents();
      

      if (pool) {
	if (pool->IsReady()) {
	  //Correlate triggers from this event with past associated:
	  //Fill the vector for all objects in fCorrelationsME.
	  MakeTriggers(arrayParticles, factiveTriggersME,fMECorrelations.begin(),fMECorrelations.end());
	  //Fill the vector for all objects in fCorrelationsMETrigger.
	  MakeTriggers(arrayParticles, factiveTriggersMETrigger,fMETriggerCorrelations.begin(),fMETriggerCorrelations.end());
	  //Fill the vector for all objects in fCorrelationsMETA.
	  MakeTriggers(arrayParticles, factiveTriggersMETA,fMETACorrelations.begin(),fMETACorrelations.end());
	  //Fill the vector for all objects in fCorrelationsMETA2.
	  MakeTriggers(arrayParticles, factiveTriggersMETA2,fMETA2Correlations.begin(),fMETA2Correlations.end());	  
	  for (int nEvent=0; nEvent<EventsInPool; nEvent++) {
	    for(int mEvent=nEvent+1;mEvent<EventsInPool;mEvent++){
	      // all particles from different events
	      MakeAssociated(pool->GetEvent(nEvent), fAssociatedmixed1);
	      MakeAssociated(pool->GetEvent(mEvent), fAssociatedmixed2);
	      ProcessEvent(factiveTriggersME,fAssociatedmixed1,fAssociatedmixed2);
	      ProcessEvent(factiveTriggersME,fAssociatedmixed2,fAssociatedmixed1,false,false);
	      ProcessEvent(factiveTriggersMETA,fAssociated,fAssociatedmixed1);
	      ProcessEvent(factiveTriggersMETA2,fAssociatedmixed1,fAssociated);
	      // associated particles from same event (skips the last event)
	      ProcessEvent(factiveTriggersMETrigger, fAssociatedmixed1);
	    }//end inner loop
	    if(nEvent==EventsInPool-1){//last event.
	      // associated particles from same event (only last event)
	      MakeAssociated(pool->GetEvent(nEvent), fAssociatedmixed1);
	      ProcessEvent(factiveTriggersMETrigger, fAssociatedmixed1);
	      ProcessEvent(factiveTriggersMETA,fAssociated,fAssociatedmixed1);
	      ProcessEvent(factiveTriggersMETA2,fAssociatedmixed1,fAssociated);
	    }//end last event
	  }//End mixed event loop
	  
	  //Correlate associated from this event with triggers and associated from past events:
	  for (int nEvent=0; nEvent<EventsInPool; nEvent++){
	    MakeTriggers(pool->GetEvent(nEvent), factiveTriggersME,fMECorrelations.begin(),fMECorrelations.end());
	    //Fill the vector for all objects in fCorrelationsMETrigger.
	    MakeTriggers(pool->GetEvent(nEvent), factiveTriggersMETrigger,fMETriggerCorrelations.begin(),fMETriggerCorrelations.end());
	    //Fill the vector for all objects in fCorrelationsMETA.
	    MakeTriggers(pool->GetEvent(nEvent), factiveTriggersMETA,fMETACorrelations.begin(),fMETACorrelations.end());
	    //Fill the vector for all objects in fCorrelationsMETA.
	    MakeTriggers(pool->GetEvent(nEvent), factiveTriggersMETA2,fMETA2Correlations.begin(),fMETA2Correlations.end());	    
	    MakeAssociated(pool->GetEvent(nEvent), fAssociatedmixed1);
	    //Now factiveTriggersME* is the same event as fAssociatedmixed1.
	    ProcessEvent(factiveTriggersMETA,fAssociatedmixed1,fAssociated);
	    ProcessEvent(factiveTriggersMETA2,fAssociated,fAssociatedmixed1);
	    // associated particles from same event (skips the last event)
	    ProcessEvent(factiveTriggersMETrigger, fAssociated);
	    for(int mEvent=0;mEvent<EventsInPool;mEvent++){
	      if(nEvent==mEvent)continue;//dont double count
	      // all particles from different events
	      MakeAssociated(pool->GetEvent(mEvent), fAssociatedmixed2);
	      ProcessEvent(factiveTriggersME,fAssociated,fAssociatedmixed2);
	      ProcessEvent(factiveTriggersME,fAssociatedmixed2,fAssociated,false,false);
	    }//end inner loop
	  }//End mixed event loop
	}//End poolisready
	if (arrayParticles->GetEntriesFast()>0) {//update event pool
	  EventsInPool=pool->UpdatePool(arrayParticles);
	}//End update event Pool
      }//End ifPOOl
      
    }
    factiveTriggers.clear();
    factiveTriggersME.clear();
    factiveTriggersMETrigger.clear();
    fAssociated.clear();
    fAssociatedmixed1.clear();
    fAssociatedmixed2.clear();
    return iResult;
  }
//   template <class InputIterator>  
  int ProcessEvent(const std::vector<AliThreeParticleCorrelator::AliActiveTrigger*>& activeTriggers, const std::vector<AliVParticle*>& associated) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    Double_t NAssociated = associated.size();
    if(NAssociated==0) return 0;//No associated means we need not fill anything.
    if (activeTriggers.size()==0) return 0;//No Triggers means we need not fill anything
    for (typename std::vector<AliActiveTrigger*>::const_iterator i=activeTriggers.begin(), e=activeTriggers.end(); i!=e; i++) {
      (*i)->Incrementtrigger();//Fill histogram for number of triggers.
      for (typename std::vector<AliVParticle*>::const_iterator assoc=associated.begin(), eassoc=associated.end(); assoc!=eassoc; assoc++) {
	if (*assoc==(*i)->GetTrigger()) continue;//Trigger and associated pointer are the same particle.
	for (typename std::vector<AliVParticle*>::const_iterator assoc2=assoc+1, eassoc2 = associated.end(); assoc2 !=eassoc2;assoc2++){
	  if (!((*i)->CheckAssociated(*assoc)&& (*i)->CheckAssociated(*assoc2))){continue;}//do not fill if either is not an associated
	  if(*assoc2==(*i)->GetTrigger()) {continue;}//Do not fill if they are the same or one is the trigger
	  if(*assoc==*assoc2) continue;
// 	  (*i)->Correlate(*assoc,*assoc2, NAssociated);
// 	  (*i)->Correlate(*assoc2,*assoc, NAssociated);
	  (*i)->Correlate(*assoc,*assoc2, 1.0);
	  (*i)->Correlate(*assoc2,*assoc, 1.0);	  
	}//loop over second associated
	if ((*i)->CheckAssociated(*assoc)) (*i)->Correlate(*assoc);//2p correlation
      } // loop over first associated
    } // loop over triggers
    return 0;
  }

//   template <class InputIterator>
  int ProcessEvent(const std::vector<AliThreeParticleCorrelator::AliActiveTrigger*>& activeTriggers,const std::vector<AliVParticle*>& associated, const std::vector<AliVParticle*>& associatedmixed, bool pickone = false,bool twop = true) {
    /// Three particle correlation loop over array of AliVParticle objects
    /// Fill correlation objects of different properties 
    Double_t NAssociated1 = associated.size();
    Double_t NAssociated2 = associatedmixed.size();
    if(!(associated==associatedmixed)){//If they are the same, 1 is substracted, so add 1 if they are not and there are more combinations.
      NAssociated1+=1;
      NAssociated2+=1;}
    if(NAssociated1==0||NAssociated2==0) return 0;//No associated means we need not fill anything.
    if (activeTriggers.size()==0) return 0;//no triggers means nothing to be correlated
    if(pickone){
      Int_t NA1 = (Int_t)(NAssociated1-1)*gRandom->Rndm();//gives the index of one associated between 0 and NAssociated1-1
      Int_t NA2 = (Int_t)(NAssociated2-1)*gRandom->Rndm();//gives the index of one associatedmixed between 0 and NAssociated2-1
      for (typename std::vector<AliActiveTrigger*>::const_iterator i=activeTriggers.begin(), e=activeTriggers.end(); i!=e; i++) {
	(*i)->Incrementtrigger();//Fill histogram for number of triggers.
	if(!((*i)->CheckAssociated(associated[NA1])&&(*i)->CheckAssociated(associatedmixed[NA2])))continue;
	(*i)->Correlate(associated[NA1],associatedmixed[NA2],1);
	(*i)->Correlate(associatedmixed[NA2],associated[NA1],1);
      }
      return 0;
    }
    for (typename std::vector<AliActiveTrigger*>::const_iterator i=activeTriggers.begin(), e=activeTriggers.end(); i!=e; i++) {
      (*i)->Incrementtrigger();//Fill histogram for number of triggers.
      for (typename std::vector<AliVParticle*>::const_iterator assoc=associated.begin(), eassoc=associated.end(); assoc!=eassoc; assoc++) {
	if(*assoc==(*i)->GetTrigger()) continue;//It is possible for this to happen with associated and trigger from same event
	for (typename std::vector<AliVParticle*>::const_iterator assoc2=associatedmixed.begin(), eassoc2 = associatedmixed.end(); assoc2 !=eassoc2;assoc2++){
	  if (!((*i)->CheckAssociated(*assoc) && (*i)->CheckAssociated(*assoc2))) {continue;}//If either is not associated, do not fill.}
	  if(*assoc==*assoc2||*assoc2 == (*i)->GetTrigger()) continue;//if associated are the samme or trigger is equal to the second associated.
// 	  (*i)->Correlate(*assoc,*assoc2,(NAssociated1-1)/NAssociated2);//Fill with the number of associated in the other list.
// 	  (*i)->Correlate(*assoc2,*assoc,(NAssociated1-1)/NAssociated2);//both need be considered
	  (*i)->Correlate(*assoc,*assoc2,1.0);//Fill with the number of associated in the other list.
	}//loop over second associated
	if ((*i)->CheckAssociated(*assoc)&&twop)(*i)->Correlate(*assoc);//2p correlation
	} // loop over first associated
// 	for (typename std::vector<AliVParticle*>::const_iterator assoc2=associatedmixed.begin(), eassoc2 = associatedmixed.end(); assoc2 !=eassoc2;assoc2++){
// 	  if (dynamic_cast<AliVParticle*>(*assoc2)==i->GetTrigger()) continue;
// 	  if (i->CheckAssociated(*assoc2))i->Correlate(*assoc2);//2p correlation
// 	}//second loop so all mixed tracks get filled in 2p correlation.
    } // loop over triggers
    return 0;
  }
  
  template <class InputIterator>
  void MakeTriggers(const TObjArray* arrayParticles,std::vector<AliActiveTrigger*>& triggers,InputIterator firstAnalysisObject, InputIterator endAnalysisObject) {
    /// create a particle array with reduced data objects
    /// for trigger particles containing only potential triggers.
    if (!arrayParticles) return ;
    triggers.clear();
    TIter nexttrigger(arrayParticles);
    TObject* otrigger=NULL;
    while ((otrigger=nexttrigger())!=NULL) {//loop over all in the array.
      AliVParticle* ptrigger=reinterpret_cast<AliVParticle*>(otrigger);
      if (!ptrigger) continue;
      for (typename vector<C*>::iterator o=firstAnalysisObject,e=endAnalysisObject;o!=e; o++) {
	if (!(*o)->CheckTrigger(ptrigger, true)) continue;//check for each correlation, if it is a trigger in it. If not, reject.
	triggers.push_back(new AliActiveTrigger(*o, ptrigger));
      }
    }//end while loop
    return ;
  }
  
  void MakeAssociated(const TObjArray* arrayParticles,std::vector<AliVParticle*>& associated) {
    /// Create clone of particle array using a reduced track class and
    if (!arrayParticles) return ;
    associated.clear();
    TIter iassociated(arrayParticles);
    TObject* o1=NULL;
    while ((o1=iassociated())!=NULL) {
      // loop over first particle
      AliVParticle* associatedp=dynamic_cast<AliVParticle*>(o1);
      if (!associatedp) continue;
      // scope for local iterators
      typename vector<C*>::iterator o=fCorrelations.begin(), e=fCorrelations.end();//fine, since mixed events has the same function for ass selection as at least one signal.
      for (;o!=e; o++) if ((*o)->CheckAssociated(associatedp, NULL, false)) break;
      if (o==e) continue;// next particle if current one does not pass cuts (no break anywhere).
      associated.push_back(associatedp);	
    } // loop over particle
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
  std::vector<C*> fCorrelations; //! worker objects
  std::vector<C*> fMECorrelations; //! ME analysis clones of workers, all three particles different events, only one associated from each event
  std::vector<C*> fMETriggerCorrelations; //! ME analysis clones of workers, associated particles same events
  std::vector<C*> fMETACorrelations; //! ME analysis clones of workers, first associated particle from the same events as the trigger
  std::vector<C*> fMETA2Correlations; //! ME analysis clones of workers, second associated particle from the same events as the trigger
  std::vector<AliActiveTrigger*> factiveTriggers; //Vector to contain the triggers.
  std::vector<AliActiveTrigger*> factiveTriggersME; //Vector to contain the triggers for mixed event.
  std::vector<AliActiveTrigger*> factiveTriggersMETrigger;//Vector to contain the triggers for mixed event with both associated from the same event.
  std::vector<AliActiveTrigger*> factiveTriggersMETA;//Vector to contain the triggers for mixed event with first associated particle from the same events as the trigger
  std::vector<AliActiveTrigger*> factiveTriggersMETA2;//Vector to contain the triggers for mixed event with second associated particle from the same events as the trigger
  std::vector<AliVParticle*> fAssociated;//vector to contain the associated particles.
  std::vector<AliVParticle*> fAssociatedmixed1; //! buffer for Associated
  std::vector<AliVParticle*> fAssociatedmixed2; //! buffer for Associated
  AliEventPoolManager* fEventPoolMgr; //! event pool manager, external pointer
  Double_t fVz;//Vertex in z
  Double_t fMultiplicity;//Multiplicity in %
  TRandom3 * fRandom;//to be able to pick mixed event tracks.
  
  
  
 ClassDef(AliThreeParticleCorrelator, 2) // could not get it working for template class 
};

#endif

