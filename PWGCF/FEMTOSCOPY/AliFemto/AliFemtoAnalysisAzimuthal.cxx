////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event //
// with respect to the z position of the primary vertex and event total   //
// multiplicity and uses only events in certain reaction plane angle bin  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliFemtoAnalysisAzimuthal.h"
#include <TMath.h>
#include <string>
#include <cstdio>
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoPairCut.h"
#include "TVector2.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoAnalysisAzimuthal)
#endif

extern void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
				      AliFemtoEvent*               hbtEvent,
				      AliFemtoParticleCollection*  partCollection);


//____________________________
AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(unsigned int binsVertex, double minVertex, double maxVertex,
						       unsigned int binsMult, double minMult, double maxMult, unsigned short binsRP) 
  : 
  fFemtoParticleCut(0),
  fFlowParticleCut(0),
  fVertexZBins(binsVertex), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(binsMult) ,
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(binsRP),
  fPsi(0)
{
  //  mControlSwitch     = 0;
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexZ[0] = minVertex;
  fVertexZ[1] = maxVertex;
  fMult[0] = minMult;
  fMult[1] = maxMult;
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins,fVertexZ[0],fVertexZ[1],
										     fMultBins,fMult[0],fMult[1],
										     fRPBins,0.0,TMath::Pi());
}

//____________________________

AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) : 
  AliFemtoSimpleAnalysis(),
  fFemtoParticleCut(0),
  fFlowParticleCut(0),
  fVertexZBins(a.fVertexZBins), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(a.fMultBins) ,
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(a.fRPBins),
  fPsi(0)

{
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexZ[0] = a.fVertexZ[0]; 
  fVertexZ[1] = a.fVertexZ[1];
  fMult[0] = a.fMult[0]; 
  fMult[1] = a.fMult[1];
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins,fVertexZ[0],fVertexZ[1],
										     fMultBins,fMult[0],fMult[1],
										     fRPBins,0.0,TMath::Pi());
  // find the right event cut
  fEventCut = a.fEventCut->Clone();
  // find the right femto particle cut
  fFemtoParticleCut = a.fFemtoParticleCut->Clone();
  // find the right flow particle cut
  fFlowParticleCut = a.fFlowParticleCut->Clone();
  // find the right pair cut
  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - event cut set " << endl;
  }
  if ( fFemtoParticleCut ) {
      SetFirstParticleCut(fFemtoParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - femto particle cut set " << endl;
  }
  if ( fFlowParticleCut ) {
      SetSecondParticleCut(fFlowParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - flow particle cut set " << endl;
  }
  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - correlation function not found " << endl;
  }
  fNumEventsToMix = a.fNumEventsToMix;
  cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - analysis copied " << endl;
}

AliFemtoAnalysisAzimuthal& AliFemtoAnalysisAzimuthal::operator=(const AliFemtoAnalysisAzimuthal& a)
{
  // Assignment operator
  if (this == &a)
    return *this;

  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexZ[0] = a.fVertexZ[0]; 
  fVertexZ[1] = a.fVertexZ[1];
  fMult[0] = a.fMult[0]; 
  fMult[1] = a.fMult[1];
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexZBins,fVertexZ[0],fVertexZ[1],
										     fMultBins,fMult[0],fMult[1],
										     fRPBins,0.0,TMath::Pi());
  // find the right event cut
  fEventCut = a.fEventCut->Clone();
  // find the right femto particle cut
  fFemtoParticleCut = a.fFemtoParticleCut->Clone();
  // find the right flow particle cut
  fFlowParticleCut = a.fFlowParticleCut->Clone();
  // find the right pair cut
  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - event cut set " << endl;
  }
  if ( fFemtoParticleCut ) {
      SetFirstParticleCut(fFemtoParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - femto particle cut set " << endl;
  }
  if ( fFlowParticleCut ) {
      SetSecondParticleCut(fFlowParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - flow particle cut set " << endl;
  }
  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - correlation function not found " << endl;
  }
  fNumEventsToMix = a.fNumEventsToMix;
  cout << " AliFemtoAnalysisAzimuthal::AliFemtoAnalysisAzimuthal(const AliFemtoAnalysisAzimuthal& a) - analysis copied " << endl;

  return *this;
  
}

//____________________________
AliFemtoAnalysisAzimuthal::~AliFemtoAnalysisAzimuthal(){
  // now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete fPicoEventCollectionVectorHideAway;
}

//_________________________
void AliFemtoAnalysisAzimuthal::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Perform event processing in bins of z vertex, multiplicity and Reaction Plane angle
  // cout << " AliFemtoAnalysisAzimuthal::ProcessEvent(const AliFemtoEvent* hbtEvent) " << endl;

  //****from AliFemtoSimpleAnalysis****

  fPicoEvent=0; // we will get a new pico event, if not prevent corr. fctn to access old pico event
  fNeventsProcessed++;

  // startup for EbyE 
  fFemtoParticleCut->EventBegin(hbtEvent);
  fFlowParticleCut->EventBegin(hbtEvent);
  fPairCut->EventBegin(hbtEvent);

  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventBegin(hbtEvent);
  }

  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  if (!tmpPassEvent) 
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  if (tmpPassEvent) {
    fPicoEvent = new AliFemtoPicoEvent; // this is what we will make pairs from and put in Mixing Buffer, no memory leak. we will delete picoevents when they come out of the mixing buffer
    FillHbtParticleCollection(fFemtoParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->FirstParticleCollection());
    FillHbtParticleCollection(fFlowParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->SecondParticleCollection());
    
      // get right mixing buffer
  double vertexZ = hbtEvent->PrimVertPos().z();
  double mult = hbtEvent->UncorrectedNumberOfPrimaries();
  TVector2 tQ = GetQVector(fPicoEvent->SecondParticleCollection());
  double tPsi=tQ.Phi()/2.;
  if (tPsi > TMath::Pi()) tPsi -= TMath::Pi();
      
   fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ,mult,fPsi); 
  if (!fMixingBuffer) {
    if ( vertexZ < fVertexZ[0] ) fUnderFlowVertexZ++;
    if ( vertexZ > fVertexZ[1] ) fOverFlowVertexZ++;
    if ( mult < fMult[0] ) fUnderFlowMult++;
    if ( mult > fMult[1] ) fOverFlowMult++;
    return;
  }   
    
    //cout << "#particles in Collection: " << fPicoEvent->FirstParticleCollection()->size() << endl;
    
    //switch which allows only using events with ParticleCollections containing a minimum number of entries
    if (fPicoEvent->FirstParticleCollection()->size() >= fMinSizePartCollection ) {
      fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
    }

      //------ Make real pairs (assume identical) ------//
      MakePairs("real", fPicoEvent->FirstParticleCollection() );
      //cout << "AliFemtoAnalysisAzimuthal::ProcessEvent() - reals done ";

      //---- Make pairs for mixed events, looping over events in mixingBuffer ----//
      AliFemtoPicoEvent* storedEvent;
      AliFemtoPicoEventIterator fPicoEventIter;
      for (fPicoEventIter=MixingBuffer()->begin();fPicoEventIter!=MixingBuffer()->end();fPicoEventIter++){
        storedEvent = *fPicoEventIter;
        MakePairs("mixed",fPicoEvent->FirstParticleCollection(),
                            storedEvent->FirstParticleCollection() );
      }
      //cout << " - mixed done   " << endl;

      //--------- If mixing buffer is full, delete oldest event ---------//
      if ( MixingBufferFull() ) {
        delete MixingBuffer()->back();
        MixingBuffer()->pop_back();
      }

      //-------- Add current event (fPicoEvent) to mixing buffer --------//
      MixingBuffer()->push_front(fPicoEvent);

    }  // if ParticleCollections are big enough (mal jun2002)
    else{
      fEventCut->FillCutMonitor(hbtEvent, !tmpPassEvent);
      delete fPicoEvent;
    }

  // if currentEvent is accepted by currentAnalysis cleanup for EbyE 
  fFemtoParticleCut->EventEnd(hbtEvent);
  fFlowParticleCut->EventEnd(hbtEvent);
  fPairCut->EventEnd(hbtEvent);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventEnd(hbtEvent);
  } 
}

//_______________________________________________________________________________
void AliFemtoAnalysisAzimuthal::MakePairs(const char* typeIn, AliFemtoParticleCollection *partCollection1,
				       AliFemtoParticleCollection *partCollection2){
// Build pairs, check pair cuts, and call CFs' AddRealPair() or AddMixedPair() methods. 
// If no second particle collection is specfied, make pairs within first particle collection.

  string type = typeIn;

  //  int swpart = ((long int) partCollection1) % 2;
  int swpart = fNeventsProcessed % 2;

  AliFemtoPair* tPair = new AliFemtoPair;
  AliFemtoCorrFctnIterator tCorrFctnIter;
  AliFemtoParticleIterator tPartIter1, tPartIter2;

  AliFemtoParticleIterator tStartOuterLoop = partCollection1->begin();  // always
  AliFemtoParticleIterator tEndOuterLoop   = partCollection1->end();    // will be one less if identical
  AliFemtoParticleIterator tStartInnerLoop;
  AliFemtoParticleIterator tEndInnerLoop;
  if (partCollection2) {                         //   Two collections:
    tStartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    tEndInnerLoop   = partCollection2->end();    
  }
  else {                                         //   One collection:
    tEndOuterLoop--;                             //   Outer loop goes to next-to-last particle
    tEndInnerLoop = partCollection1->end() ;     //   Inner loop goes to last particle
  }
  for (tPartIter1=tStartOuterLoop;tPartIter1!=tEndOuterLoop;tPartIter1++) {
    if (!partCollection2){
      tStartInnerLoop = tPartIter1;
      tStartInnerLoop++;
    }
    tPair->SetTrack1(*tPartIter1);
    for (tPartIter2 = tStartInnerLoop; tPartIter2!=tEndInnerLoop;tPartIter2++) {
      tPair->SetTrack2(*tPartIter2);

//For getting the pair angle wrt EP
  if (type == "real"){
	TVector2 tQ, tQVector;
	float tPsi=0, mQx=0, mQy=0, PairAngleEP=0;
	tQ = GetQVector(partCollection1);

	mQx = tQ.X() - cos(2*(tPair->Track1()->FourMomentum().Phi()))*(tPair->Track1()->Track()->Pt()) - cos(2*(tPair->Track2()->FourMomentum().Phi()))*(tPair->Track2()->Track()->Pt());
	mQy = tQ.Y() - sin(2*(tPair->Track1()->FourMomentum().Phi()))*(tPair->Track1()->Track()->Pt()) - sin(2*(tPair->Track2()->FourMomentum().Phi()))*(tPair->Track2()->Track()->Pt());
  
	tQVector.Set(mQx,mQy);

	tPsi=tQVector.Phi()/2.;
	if (tPsi > TMath::Pi()) tPsi -= TMath::Pi();

	PairAngleEP = (tPair->EmissionAngle() - tPsi);
	tPair->SetPairAngleEP(PairAngleEP);
  }

  if (type == "mixed"){
	float tPsi1=0, tPsi2=0, mQx1=0, mQx2=0,mQy1=0, mQy2=0, px1=0, px2=0, py1=0, py2=0, PairAngleEP=0;
	TVector2 tQ1, tQ2, tQVector1, tQVector2, tP;

	tQ1 = GetQVector(partCollection1);
	tQ2 = GetQVector(partCollection2);

	mQx1 = tQ1.X() - cos(2*(tPair->Track1()->FourMomentum().Phi()))*(tPair->Track1()->Track()->Pt());
	mQx2 = tQ2.X() - cos(2*(tPair->Track2()->FourMomentum().Phi()))*(tPair->Track2()->Track()->Pt());
	mQy1 = tQ1.Y() - sin(2*(tPair->Track1()->FourMomentum().Phi()))*(tPair->Track1()->Track()->Pt());
	mQy2 = tQ2.Y() - sin(2*(tPair->Track2()->FourMomentum().Phi()))*(tPair->Track2()->Track()->Pt());
  
	tQVector1.Set(mQx1,mQy1);
	tQVector2.Set(mQx2,mQy2);

	tPsi1=tQVector1.Phi()/2.;
	if (tPsi1 > TMath::Pi()) tPsi1 -= TMath::Pi();

	tPsi2=tQVector2.Phi()/2.;
	if (tPsi2 > TMath::Pi()) tPsi2 -= TMath::Pi();

	px1 = (tPair->Track1()->Track()->Pt())*cos(tPair->Track1()->FourMomentum().Phi() - tPsi1);
	px2 = (tPair->Track2()->Track()->Pt())*cos(tPair->Track2()->FourMomentum().Phi() - tPsi2);
	py1 = (tPair->Track1()->Track()->Pt())*sin(tPair->Track1()->FourMomentum().Phi() - tPsi1);
	py2 = (tPair->Track2()->Track()->Pt())*sin(tPair->Track2()->FourMomentum().Phi() - tPsi2);

	tP.Set(px1+px2, py1+py2);
	PairAngleEP = tP.Phi();

	tPair->SetPairAngleEP(PairAngleEP);
  }

      // The following lines have to be uncommented if you want pairCutMonitors
      // they are not in for speed reasons
      // bool tmpPassPair = fPairCut->Pass(tPair);
      // fPairCut->FillCutMonitor(tPair, tmpPassPair);
      // if ( tmpPassPair )

      //---- If pair passes cut, loop over CF's and add pair to real/mixed ----//

      if (!partCollection2) {
	if (swpart) {
 	  tPair->SetTrack1(*tPartIter2);
 	  tPair->SetTrack2(*tPartIter1);
 	  swpart = 0;
 	}
 	else {
 	  tPair->SetTrack1(*tPartIter1);
 	  tPair->SetTrack2(*tPartIter2);
 	  swpart = 1;
	}
      }


      if (fPairCut->Pass(tPair)){
        for (tCorrFctnIter=fCorrFctnCollection->begin();
             tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
          AliFemtoCorrFctn* tCorrFctn = *tCorrFctnIter;
          if(type == "real")
            tCorrFctn->AddRealPair(tPair);
	  else if(type == "mixed")
            tCorrFctn->AddMixedPair(tPair);
          else
            cout << "Problem with pair type, type = " << type.c_str() << endl;
        }
      }
    }    // loop over second particle
  }      // loop over first particle

  delete tPair;
}

//_____________________________________________
TVector2 AliFemtoAnalysisAzimuthal::GetQVector(AliFemtoParticleCollection* particlecollection){
  
  TVector2 mQ;
  float mQx=0, mQy=0;

  if (!particlecollection) {
    mQ.Set(0.0, 0.0);
    return mQ;
  }

  AliFemtoParticle* flowparticle;
  AliFemtoParticleIterator pIter;
  AliFemtoParticleIterator startLoop = particlecollection->begin();
  AliFemtoParticleIterator endLoop   = particlecollection->end();
  for (pIter=startLoop;pIter!=endLoop;pIter++){
    flowparticle = *pIter;
    mQx += (cos(2*flowparticle->FourMomentum().Phi()))*(flowparticle->Track()->Pt());
    mQy += (sin(2*flowparticle->FourMomentum().Phi()))*(flowparticle->Track()->Pt());
  }
  
  mQ.Set(mQx,mQy);
  return mQ;
}

//__________________________________________________
double AliFemtoAnalysisAzimuthal::GetCurrentReactionPlane()
{
  return fPsi;
}

//_________________________
TList* AliFemtoAnalysisAzimuthal::GetOutputList()
{
  // Collect the list of output objects to be written

  TList *tOutputList = new TList();
  TList *p1Cut = fFemtoParticleCut->GetOutputList();

  TListIter nextp1(p1Cut);
  while (TObject *obj = nextp1.Next()) {
    tOutputList->Add(obj);
  }

  TList *pairCut = fPairCut->GetOutputList();

  TIter nextpair(pairCut);
  while (TObject *obj = nextpair()) {
    tOutputList->Add(obj);
  }

  TList *eventCut = fEventCut->GetOutputList();

  TIter nextevent(eventCut);
  while (TObject *obj = nextevent()) {
    tOutputList->Add(obj);
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    TList *tListCf = (*iter)->GetOutputList();
    
    TIter nextListCf(tListCf);
    while (TObject *obj = nextListCf()) {
      tOutputList->Add(obj);
    }
  }

  return tOutputList;
  
}
