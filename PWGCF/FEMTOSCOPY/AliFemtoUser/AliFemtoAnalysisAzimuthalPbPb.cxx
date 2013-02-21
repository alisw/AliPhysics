////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoAnalysisReactionPlane - Femtoscopic analysis which mixes event //
// with respect to the z position of the primary vertex and event total   //
// multiplicity and uses only events in certain reaction plane angle bin  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * Author: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *         Jorge Mercado, University of Heidelberg, jmercado@cern.ch
 *         Vera Loggins, Wayne State University, veraloggins@wayne.edu
 *
 **************************************************************************/

#include "AliFemtoAnalysisAzimuthalPbPb.h"
#include <TMath.h>
#include <string>
#include <cstdio>
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoPairCut.h"
#include "AliFemtoPairCutRadialDistanceLM.h"
#include "AliVTrack.h"
#include "TVector2.h"
#include "TTreeFormula.h"
#include "AliFemtoKinkCut.h"
#include "AliEventplane.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoPicoEventRP.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoAnalysisAzimuthalPbPb)
#endif

void FillHbtParticleCollection(       AliFemtoParticleCut*         partCut,
				      AliFemtoEvent*               hbtEvent,
				      AliFemtoPicoEventRP* 	   picoevent)
{
      AliEventplane* evpl;
      evpl = picoevent->PicoEventplane();
      *evpl = *hbtEvent->EP();
      AliFemtoTrackCut* pCut = (AliFemtoTrackCut*) partCut;
      AliFemtoTrack* pParticle;
      AliFemtoTrackIterator pIter;
      AliFemtoTrackIterator startLoop = hbtEvent->TrackCollection()->begin();
      AliFemtoTrackIterator endLoop   = hbtEvent->TrackCollection()->end();
      for (pIter=startLoop;pIter!=endLoop;pIter++){
	pParticle = *pIter;
	bool tmpPassParticle = pCut->Pass(pParticle);
	pCut->FillCutMonitor(pParticle, tmpPassParticle);
	if (tmpPassParticle){	
	  AliFemtoParticle* particle = new AliFemtoParticle(pParticle,pCut->Mass());
	  picoevent->FirstParticleCollection()->push_back(particle);
	}
      }
}

//____________________________
AliFemtoAnalysisAzimuthalPbPb::AliFemtoAnalysisAzimuthalPbPb(unsigned int binsVertex, double minVertex, double maxVertex,
						       unsigned int binsMult, double minMult, double maxMult, unsigned short binsRP) 
  : 
  fFirstParticleCut(0),
  fSecondParticleCut(0),
  fPairCutRD(0),
  fPicoEventRP(0),
  fVertexZBins(binsVertex), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(binsMult),
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(binsRP),
  fRP(0),
  fphidist(0),
  fpairphi(0),
  fRPdist(0),
  fsubRPdist(0),
  frealpsi(0),
  fmixedpsi(0)
{
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
  
  fphidist = new TH1F("fphidist","fphidist; Phi Distribution",100,-TMath::Pi(),TMath::Pi());
  fpairphi = new TH1F("fpairphi","fpairphi; Pair Phi Distribution",100,0,TMath::TwoPi());
  fRPdist = new TH1F("fRPdist","fRPdist; RP Distribution",100,0,TMath::Pi());
  fsubRPdist = new TH1F("fsubRPdist","fsubRPdist; sub RP Distribution",200,-TMath::Pi(),TMath::Pi());
  frealpsi = new TH1F("frealpsi","frealpsi; real Psi Distribution",100,0,TMath::RadToDeg()*TMath::Pi());
  fmixedpsi = new TH1F("fmixedpsi","fmixedpsi; mixed Psi Distribution",100,0,TMath::RadToDeg()*TMath::Pi());
}

//____________________________

AliFemtoAnalysisAzimuthalPbPb::AliFemtoAnalysisAzimuthalPbPb(const AliFemtoAnalysisAzimuthalPbPb& a) : 
  AliFemtoSimpleAnalysis(),
  fFirstParticleCut(0),
  fSecondParticleCut(0),
  fPairCutRD(0),
  fPicoEventRP(0),
  fVertexZBins(a.fVertexZBins), 
  fOverFlowVertexZ(0), 
  fUnderFlowVertexZ(0),
  fMultBins(a.fMultBins) ,
  fOverFlowMult(0),    
  fUnderFlowMult(0),
  fRPBins(a.fRPBins),
  fRP(0),
  fphidist(0),
  fpairphi(0),
  fRPdist(0),
  fsubRPdist(0),
  frealpsi(0),
  fmixedpsi(0) 

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
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  // find the right flow particle cut
  fSecondParticleCut = a.fSecondParticleCut->Clone();
  // find the right pair cut
  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
  }
  fNumEventsToMix = a.fNumEventsToMix;
}

AliFemtoAnalysisAzimuthalPbPb& AliFemtoAnalysisAzimuthalPbPb::operator=(const AliFemtoAnalysisAzimuthalPbPb& a)
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
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  // find the right flow particle cut
  fSecondParticleCut = a.fSecondParticleCut->Clone();
  // find the right pair cut
  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
  }
  fNumEventsToMix = a.fNumEventsToMix;

  return *this;
  
}

//____________________________
AliFemtoAnalysisAzimuthalPbPb::~AliFemtoAnalysisAzimuthalPbPb(){
  // now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete fPicoEventCollectionVectorHideAway;
}

//_________________________
void AliFemtoAnalysisAzimuthalPbPb::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Perform event processing in bins of z vertex, multiplicity and Reaction Plane angle
  //****from AliFemtoSimpleAnalysis****
  fFirstParticleCut->EventBegin(hbtEvent);
  double vertexZ = hbtEvent->PrimVertPos().z();
  double mult = hbtEvent->UncorrectedNumberOfPrimaries();
  double RP = hbtEvent->ReactionPlaneAngle();  
  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ,mult,RP); 
  if (!fMixingBuffer) {
//     cout << "no mixing buffer!!!" << endl;
    if ( vertexZ < fVertexZ[0] ) fUnderFlowVertexZ++;
    if ( vertexZ > fVertexZ[1] ) fOverFlowVertexZ++;
    if ( mult < fMult[0] ) fUnderFlowMult++;
    if ( mult > fMult[1] ) fOverFlowMult++;
    return;
  }

  // Add event to processed events
  fPicoEventRP=0; // we will get a new pico event, if not prevent corr. fctn to access old pico event
  fNeventsProcessed++;
 
  fFirstParticleCut->EventBegin(hbtEvent);
  fSecondParticleCut->EventBegin(hbtEvent);
  fPairCut->EventBegin(hbtEvent);
  fPairCutRD->EventBegin(hbtEvent);
  
  int magsign = 0;
  if(hbtEvent->MagneticField()>0) magsign = 1;
  else if(hbtEvent->MagneticField()<0) magsign = -1;
  fPairCutRD->SetMagneticFieldSign(magsign);

  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventBegin(hbtEvent);
  }
  
  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  if (!tmpPassEvent) {
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  }
  if (tmpPassEvent) {

    if (RP>0){
      fPicoEventRP = new AliFemtoPicoEventRP; // this is what we will make pairs from and put in Mixing Buffer
    // no memory leak. we will delete picoevents when they come out of the mixing buffer
    FillHbtParticleCollection(fFirstParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEventRP);
    if (fPicoEventRP->FirstParticleCollection()->size() >= fMinSizePartCollection) {
      fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
      fRPdist->Fill(RP);
      fsubRPdist->Fill(fPicoEventRP->PicoEventplane()->GetQsubRes());

        MakePairs("real", fPicoEventRP);


      //---- Make pairs for mixed events, looping over events in mixingBuffer ----//

      AliFemtoPicoEventRP* storedEvent;
      AliFemtoPicoEventIterator fPicoEventIter;
      for (fPicoEventIter=MixingBuffer()->begin();fPicoEventIter!=MixingBuffer()->end();fPicoEventIter++){
        storedEvent = (AliFemtoPicoEventRP*) *fPicoEventIter;
          MakePairs("mixed",fPicoEventRP,
                            storedEvent );
        
      }

      if ( MixingBufferFull() ) {
        delete MixingBuffer()->back();
        MixingBuffer()->pop_back();
      }

      MixingBuffer()->push_front(fPicoEventRP);

    }  // if ParticleCollections are big enough (mal jun2002)
    else{
//       cout << "here down" << endl;
      fEventCut->FillCutMonitor(hbtEvent, !tmpPassEvent);
//       cout << "and here?" << endl;
      delete fPicoEventRP;
//       cout << "and here?" << endl;
    }
  }   // if currentEvent is accepted by currentAnalysis
  fFirstParticleCut->EventEnd(hbtEvent);
  fSecondParticleCut->EventEnd(hbtEvent);
  fPairCut->EventEnd(hbtEvent);
  fPairCutRD->EventEnd(hbtEvent);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventEnd(hbtEvent);
  } 
  }
}

//_______________________________________________________________________________
void AliFemtoAnalysisAzimuthalPbPb::MakePairs(const char* typeIn, AliFemtoPicoEventRP *picoevent1,
				       AliFemtoPicoEventRP *picoevent2){
   string type = typeIn;

  int swpart = fNeventsProcessed % 2;

  AliFemtoParticleCollection* partCollection1 = picoevent1->FirstParticleCollection();
  AliFemtoParticleCollection* partCollection2 = 0;
  if (picoevent2)
    partCollection2 = picoevent2->FirstParticleCollection();
  AliFemtoPair* tPair = new AliFemtoPair;

  AliFemtoCorrFctnIterator tCorrFctnIter;

  AliFemtoParticleIterator tPartIter1, tPartIter2;

  AliFemtoParticleIterator tStartOuterLoop = partCollection1->begin();  // always
  AliFemtoParticleIterator tEndOuterLoop   = partCollection1->end();    // will be one less if identical
  AliFemtoParticleIterator tStartInnerLoop;
  AliFemtoParticleIterator tEndInnerLoop;
  if (partCollection2) {                        // Two collections:
    tStartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    tEndInnerLoop   = partCollection2->end();    //
  }
  else {                                        // One collection:
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
      
      
      //For getting the pair angle wrt EP
  if (type == "real"){
	Double_t PairAngleEP=0;
	TVector2* q=0;
	q = picoevent1->PicoEventplane()->GetQVector();

	float psi = q->Phi()/2;
	PairAngleEP = (tPair->EmissionAngle() - TMath::RadToDeg()*psi);
	while (PairAngleEP < 0) PairAngleEP += 180;
	while (PairAngleEP > 180) PairAngleEP -= 180;
	tPair->SetPairAngleEP(PairAngleEP);
	frealpsi->Fill(PairAngleEP);	
	fphidist->Fill(tPair->Track1()->FourMomentum().Phi());
	fphidist->Fill(tPair->Track2()->FourMomentum().Phi());
	fpairphi->Fill(tPair->EmissionAngle()*TMath::DegToRad());

  }

  if (type == "mixed"){

	TVector2* q1=0;
	TVector2* q2=0;
	q1 = picoevent1->PicoEventplane()->GetQVector();
	q2 = picoevent2->PicoEventplane()->GetQVector();
	Double_t PairAngleEP=0;

	
	float psi1 = q1->Phi()/2;
	float psi2 = q2->Phi()/2;
	PairAngleEP = TMath::RadToDeg()*(TMath::ATan2(((tPair->Track1()->Track()->Pt()*TMath::Sin(tPair->Track1()->FourMomentum().Phi() - psi1))+(tPair->Track2()->Track()->Pt()*TMath::Sin(tPair->Track2()->FourMomentum().Phi() - psi2))),((tPair->Track1()->Track()->Pt()*TMath::Cos(tPair->Track1()->FourMomentum().Phi() - psi1))+(tPair->Track2()->Track()->Pt()*TMath::Cos(tPair->Track2()->FourMomentum().Phi() - psi2)))));
	while (PairAngleEP < 0) PairAngleEP += 180;
	while (PairAngleEP > 180) PairAngleEP -= 180;
	fmixedpsi->Fill(PairAngleEP);
	tPair->SetPairAngleEP(PairAngleEP);
  }

      if (fPairCutRD->Pass(tPair)){
        for (tCorrFctnIter=fCorrFctnCollection->begin();
             tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
          AliFemtoCorrFctn* tCorrFctn = *tCorrFctnIter;
          if(type == "real")
            tCorrFctn->AddRealPair(tPair);
	  else if(type == "mixed") {
            tCorrFctn->AddMixedPair(tPair);
	  }
       
	}
      }
    }    // loop over second particle

  }      // loop over first particle

  delete tPair;
  
}

//_____________________________________________
TVector2 AliFemtoAnalysisAzimuthalPbPb::GetQVector(AliFemtoParticleCollection* particlecollection){
  
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
double AliFemtoAnalysisAzimuthalPbPb::GetCurrentReactionPlane()
{
  return fRP;
}

//______________________________________________________________________________
void AliFemtoAnalysisAzimuthalPbPb::SetEPhistname(char* histname)
{
  fphidist->SetName(Form("fphidist%s",histname));
  fpairphi->SetName(Form("fpairphi%s",histname));
  fRPdist->SetName(Form("fRPdist%s",histname));
  fsubRPdist->SetName(Form("fsubRPdist%s",histname));
  frealpsi->SetName(Form("frealpsi%s",histname));
  fmixedpsi->SetName(Form("fmixedpsi%s",histname));
}

//_________________________
TList* AliFemtoAnalysisAzimuthalPbPb::GetOutputList()
{
  // Collect the list of output objects to be written

  TList *tOutputList = new TList();


  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    TList *tListCf = (*iter)->GetOutputList();
    
    TIter nextListCf(tListCf);
    while (TObject *obj = nextListCf()) {
      tOutputList->Add(obj);
    }
  }

  tOutputList->Add(fphidist);
  tOutputList->Add(fpairphi);
  tOutputList->Add(fRPdist);
  tOutputList->Add(fsubRPdist);
  tOutputList->Add(frealpsi);
  tOutputList->Add(fmixedpsi);

  return tOutputList;
  
}
