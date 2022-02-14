////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSphericityEventCut - the basic cut for events.                     //
// Only cuts on event multiplicity, z-vertex position and                     //
// transverse sphericity are accepted.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoSphericityEventCut.h"
//#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoSphericityEventCut);
  /// \endcond
#endif

AliFemtoSphericityEventCut::AliFemtoSphericityEventCut() :
  AliFemtoEventCut(),
  fEventMult(),
  fVertZPos(),
  fAcceptBadVertex(false),
  fNEventsPassed(0),
  fNEventsFailed(0),
  fAcceptOnlyPhysics(0),
  fStCutMin(0.0),
  fStCutMax(1.0),
  fSelectTrigger(0)
{
  // Default constructor
  fEventMult[0] = 0;
  fEventMult[1] = 100000;
  fVertZPos[0] = -100.0;
  fVertZPos[1] = 100.0;
  fPsiEP[0] = -100.0;
  fPsiEP[1] = 100.0;
  fStCutMin = 0.0;
  fStCutMax = 1.0;
}
//------------------------------
AliFemtoSphericityEventCut::~AliFemtoSphericityEventCut()
{
  // Default destructor
}
//------------------------------
bool AliFemtoSphericityEventCut::Pass(const AliFemtoEvent* event)
{
  // Pass events if they fall within the multiplicity, z-vertex position range
  // and transverse sphericity. Fail otherwise
  //  int mult =  event->NumberOfTracks();
  int mult = (int) event->UncorrectedNumberOfPrimaries();
  double vertexZPos = event->PrimVertPos().z();

  // Double_t qxEPVZERO = 0, qyEPVZERO = 0;
  // Double_t qVZERO = -999;
  double epvzero = event->ReactionPlaneAngle();

  // cout << "AliFemtoSphericityEventCut:: epvzero:       " << fPsiEP[0] << " < " << epvzero << " < " << fPsiEP[1] << endl;
//   cout << "AliFemtoSphericityEventCut:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
//   cout << "AliFemtoSphericityEventCut:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;
//   cout << "AliFemtoSphericityEventCut:: VertexZErr: " << event->PrimVertCov()[4] << endl;

  // cout << "AliFemtoSphericityEventCut:: MagneticField: " << event->MagneticField() << endl;
  // cout << "AliFemtoSphericityEventCut:: IsCollisionCandidate: " << event->IsCollisionCandidate() << endl;
  // cout << "AliFemtoSphericityEventCut:: TriggerCluster: " << event->TriggerCluster() << endl;
  // cout << "AliFemtoSphericityEventCut:: fSelectTrigger: " << fSelectTrigger << endl;
  // cout << "AliFemtoSphericityEventCut:: " << endl;


  Int_t ParticleNumber = 0;
  Double_t SumPt = 0;
  Double_t S00=0;
  Double_t S11=0;
  Double_t S10=0;
  Double_t Lambda1 = 0;
  Double_t Lambda2 = 0;
  Double_t St = 0;

  AliFemtoTrackCollection *tracks = event->TrackCollection();


  for (AliFemtoTrackIterator iter=tracks->begin();iter!=tracks->end();iter++){

    Double_t NewPhi = (*iter)->P().Phi();
    Double_t NewPt =  (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();


    if(TMath::Abs(NewEta)>0.8 || NewPt<0.5){continue;}

    Double_t Px;
    Double_t Py;

    Px= NewPt * TMath::Cos(NewPhi);
    Py= NewPt * TMath::Sin(NewPhi);

    S00 = S00 + Px*Px/(NewPt);  // matrix elements of the transverse shpericity matrix S(i,j)
    S11 = S11 + Py*Py/(NewPt);  // i,j /in [0,1]
    S10 = S10 + Px*Py/(NewPt);
    SumPt = SumPt + NewPt;
    ParticleNumber++;

  }  	// end of track loop

  if(SumPt==0){
    return kFALSE;
  }

  S00 = S00/SumPt; // normalize
  S11 = S11/SumPt;
  S10 = S10/SumPt;

  Lambda1 = (S00 + S11 + TMath::Sqrt((S00+S11)*(S00+S11)-4.0*(S00*S11-S10*S10)))/2.0;
  Lambda2 = (S00 + S11 - TMath::Sqrt((S00+S11)*(S00+S11)-4.0*(S00*S11-S10*S10)))/2.0;

     if(Lambda1+Lambda2!=0 && ParticleNumber>2)
	{
		St = 2*Lambda2/(Lambda1+Lambda2);
	}
     else{return kFALSE;};


  //cout<<"St  = "<<St<<endl;

  if(St>fStCutMax || St<fStCutMin){
	//cout<<"Event kicked out !"<<"StCutMax= "<<fStCutMax<<"  StCutMin= "<<fStCutMin<<endl;
	//cout<<"St = "<<St<<endl;
    return kFALSE;
  }

  bool goodEvent =
    ((mult >= fEventMult[0]) &&
     (mult <= fEventMult[1]) &&
     (vertexZPos > fVertZPos[0]) &&
     (vertexZPos < fVertZPos[1]) &&
     (epvzero > fPsiEP[0]) &&
     (epvzero < fPsiEP[1]) &&
     ((!fAcceptBadVertex) || (event->ZDCParticipants() > 1.0)) &&
      ((!fSelectTrigger) || (event->TriggerCluster() == fSelectTrigger))
    );

  // cout << "AliFemtoSphericityEventCut:: goodEvent" <<goodEvent << endl;

  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
  //  cout << "AliFemtoSphericityEventCut:: return : " << goodEvent << endl;
//     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)) &&

  return goodEvent;
}
//------------------------------
AliFemtoString AliFemtoSphericityEventCut::Report()
{
  // Prepare report
  AliFemtoString report("AliFemtoSphericityEventCut Report:");
  report += Form("\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1]);
  report += Form("\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  report += Form("\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);

  return report;
}
void AliFemtoSphericityEventCut::SetAcceptBadVertex(bool b)
{
  fAcceptBadVertex = b;
}
bool AliFemtoSphericityEventCut::GetAcceptBadVertex()
{
  return fAcceptBadVertex;
}
