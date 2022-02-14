
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSphericityEventCutKK1 - the basic cut for events.                     //
// Only cuts on event multiplicity, z-vertex position and                     //
// transverse sphericity are accepted.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoSphericityEventCutKK1.h"
//#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoSphericityEventCutKK1);
  /// \endcond
#endif

AliFemtoSphericityEventCutKK1::AliFemtoSphericityEventCutKK1() :
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
  fPsiEP[0] = -1000.0;
  fPsiEP[1] = 1000.0;
  fStCutMin = 0.0;
  fStCutMax = 1.0;
} 
//------------------------------
AliFemtoSphericityEventCutKK1::~AliFemtoSphericityEventCutKK1(){
  // Default destructor
}
//------------------------------
bool AliFemtoSphericityEventCutKK1::Pass(const AliFemtoEvent* event){

  // Pass events if they fall within the multiplicity, z-vertex position range
  // and transverse sphericity. Fail otherwise
  //  int mult =  event->NumberOfTracks();
  int mult = (int) event->UncorrectedNumberOfPrimaries();
  double vertexZPos = event->PrimVertPos().z();

  // Double_t qxEPVZERO = 0, qyEPVZERO = 0;
  // Double_t qVZERO = -999;
  double epvzero = event->ReactionPlaneAngle();

//   cout << "AliFemtoSphericityEventCutKK1:: epvzero:       " << fPsiEP[0] << " < " << epvzero << " < " << fPsiEP[1] << endl;
//   cout << "AliFemtoSphericityEventCutKK1:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
 //  cout << "AliFemtoSphericityEventCutKK1:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;
//   cout << "AliFemtoSphericityEventCutKK1:: VertexZErr: " << event->PrimVertCov()[4] << endl;

  // cout << "AliFemtoSphericityEventCutKK1:: MagneticField: " << event->MagneticField() << endl;
  // cout << "AliFemtoSphericityEventCutKK1:: IsCollisionCandidate: " << event->IsCollisionCandidate() << endl;
  // cout << "AliFemtoSphericityEventCutKK1:: TriggerCluster: " << event->TriggerCluster() << endl;
  // cout << "AliFemtoSphericityEventCutKK1:: fSelectTrigger: " << fSelectTrigger << endl;
  // cout << "AliFemtoSphericityEventCutKK1:: " << endl;


  Int_t ParticleNumber = 0;
  Double_t SumPt = 0;
  Double_t S00=0; 
  Double_t S11=0;
  Double_t S10=0;
  Double_t Lambda1 = 0;
  Double_t Lambda2 = 0;
  Double_t St = 0;


 

   AliFemtoTrackCollection * tracks = event->TrackCollection(); 
   
   
  for (AliFemtoTrackIterator iter=tracks->begin();iter!=tracks->end();iter++){
  
  
    Double_t NewPhi = (*iter)->P().Phi();
    Double_t NewPt =  (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();
   
    
   if(TMath::Abs(NewEta)>0.8 || NewPt<0.5){continue;}

 //  cout<<" NewEta = " << NewEta <<" NewPt  = "<<NewPt<<endl;
    
    Double_t Px;
    Double_t Py;
    
    Px= NewPt * TMath::Cos(NewPhi);
    Py= NewPt * TMath::Sin(NewPhi);
    
    S00 = S00 + Px*Px/(NewPt);  // matrix elements of the transverse shpericity matrix S(i,j)
    S11 = S11 + Py*Py/(NewPt);  // i,j /in [0,1]
    S10 = S10 + Px*Py/(NewPt);
    SumPt = SumPt + NewPt;
    ParticleNumber++;

// cout<<"NewPt  = "<<NewPt<<endl;
    
  }  	// end of track loop

// cout<<"SumPt  = "<<SumPt<<endl;
  
    if(SumPt==0){return kFALSE;}
      
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
//	cout<<"Event kicked out !"<<"StCutMax= "<<fStCutMax<<"  StCutMin= "<<fStCutMin<<endl;
//	cout<<"St = "<<St<<endl;		  
  return kFALSE;}  


  const bool passes_ep = (fPsiEP[0] < epvzero && epvzero < fPsiEP[1]),
             passes_mult = (fEventMult[0] <= mult && mult <= fEventMult[1]),
             passes_z = (fVertZPos[0] < vertexZPos && vertexZPos < fVertZPos[1]),
             passes_zdc = ((!fAcceptBadVertex) || (event->ZDCParticipants() > 1.0)),
             passes_trigger = ((!fSelectTrigger) || (event->TriggerCluster() == fSelectTrigger));
  
  
  
     bool goodEvent = passes_ep
                        && passes_mult
                        && passes_z
                        && passes_zdc
                        && passes_trigger;


/*

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
*/
  // cout << "AliFemtoSphericityEventCutKK1:: goodEvent" <<goodEvent << endl;

  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
//    cout << "AliFemtoSphericityEventCutKK1:: return : " << goodEvent << endl;
//     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)) &&

  return (goodEvent);
}
//------------------------------
AliFemtoString AliFemtoSphericityEventCutKK1::Report(){
  // Prepare report
  string stemp;
  char ctemp[100];
  snprintf(ctemp , 100, "\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1]);
  stemp = ctemp;
  snprintf(ctemp , 100, "\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  stemp += ctemp;
  snprintf(ctemp , 100, "\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;
}
void AliFemtoSphericityEventCutKK1::SetAcceptBadVertex(bool b)
{
  fAcceptBadVertex = b;
}
bool AliFemtoSphericityEventCutKK1::GetAcceptBadVertex()
{
  return fAcceptBadVertex;
}

