////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoSpherocityEventCut - the basic cut for events.                     //
// Only cuts on event multiplicity, z-vertex position and                     //
// transverse spherocity are accepted.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoSpherocityEventCut.h"
//#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoSpherocityEventCut);
  /// \endcond
#endif

AliFemtoSpherocityEventCut::AliFemtoSpherocityEventCut():
  AliFemtoEventCut(),
  fEventMult(),
  fVertZPos(),
  fAcceptBadVertex(false),
  fNEventsPassed(0),
  fNEventsFailed(0),
  fAcceptOnlyPhysics(0),
  fSoCutMin(0.0),
  fSoCutMax(1.0),
  fSelectTrigger(0)
{
  // Default constructor
  fEventMult[0] = 0;
  fEventMult[1] = 100000;
  fVertZPos[0] = -100.0;
  fVertZPos[1] = 100.0;
  fPsiEP[0] = -1000.0;
  fPsiEP[1] = 1000.0;
}
//------------------------------
AliFemtoSpherocityEventCut::~AliFemtoSpherocityEventCut()
{
  // Default destructor
}

//------------------------------
bool AliFemtoSpherocityEventCut::Pass(const AliFemtoEvent* event)
{

  // Pass events if they fall within the multiplicity and z-vertex
  // position range. Fail otherwise
  //  int mult =  event->NumberOfTracks();
  int mult = (int) event->UncorrectedNumberOfPrimaries();
  double vertexZPos = event->PrimVertPos().z();
  double spherocity = -10;
  int MULT = 0;

  AliFemtoTrackCollection *tracks = event->TrackCollection();
  for (AliFemtoTrackIterator iter = tracks->begin(); iter != tracks->end(); iter++) {

    Double_t NewPt = (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();
    if (TMath::Abs(NewEta) > 0.8 || NewPt < 0.5) {
      continue;
    }

    MULT++;

  }
  //if(SumPt==0){return kFALSE;}
  if (MULT < 3) {
    return kFALSE;
  }

  Double_t *pxA = new Double_t[MULT]();
  Double_t *pyA = new Double_t[MULT]();

  Double_t sumapt = 0;
  Int_t counter = 0;

  AliFemtoTrackCollection *tracks2 = event->TrackCollection();
  for (AliFemtoTrackIterator iter2 = tracks2->begin(); iter2 != tracks2->end(); iter2++) {

    Double_t NewPt2 = (*iter2)->Pt();
    Double_t NewPhi2 = (*iter2)->P().Phi();
    Double_t NewEta2 = (*iter2)->P().PseudoRapidity();
    if (TMath::Abs(NewEta2) > 0.8 || NewPt2 < 0.5) {
      continue;
    }

    Double_t Px;
    Double_t Py;

    Px= NewPt2 * TMath::Cos(NewPhi2);
    Py= NewPt2 * TMath::Sin(NewPhi2);

    pxA[counter] = Px;
    pyA[counter] = Py;
    sumapt += NewPt2;
    counter++;
  }

  Double_t pFull = 0;
  Double_t Spherocity = 2;
  //Getting thrust
  for (Int_t i = 0; i < 360; ++i) {
    Double_t numerador = 0;
    Double_t phiparam  = 0;
    Double_t nx = 0;
    Double_t ny = 0;
    phiparam = ((TMath::Pi()) * i) / 180; // parametrization of the angle
    nx = TMath::Cos(phiparam);            // x component of an unitary vector n
    ny = TMath::Sin(phiparam);            // y component of an unitary vector n
    for (Int_t i1 = 0; i1 < MULT; ++i1) {
      numerador += TMath::Abs(ny * pxA[i1] - nx * pyA[i1]);//product between momentum proyection in XY plane and the unitari vector.
    }
    pFull = TMath::Power(numerador / sumapt, 2);
    if (pFull < Spherocity) {//maximization of pFull
      Spherocity = pFull;
    }
  }


  spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;

  if (pxA) {// clean up array memory used for TMath::Sort
    delete[] pxA;
    pxA = 0;
  }
  if (pyA) {// clean up array memory used for TMath::Sort
    delete[] pyA;
    pyA = 0;
  }

  if(spherocity>fSoCutMax || spherocity<fSoCutMin) {
    //cout<<" Event kicked out !"<<"SoCutMax= "<<fSoCutMax<<"  SoCutMin= "<<fSoCutMin<<endl;
    return kFALSE;
  }

  double epvzero = event->ReactionPlaneAngle();

  // cout << "AliFemtoSpherocityEventCut:: epvzero:       " << fPsiEP[0] << " < " << epvzero << " < " << fPsiEP[1] << endl;
//   cout << "AliFemtoSpherocityEventCut:: mult:       " << fEventMult[0] << " < " << mult << " < " << fEventMult[1] << endl;
//   cout << "AliFemtoSpherocityEventCut:: VertexZPos: " << fVertZPos[0] << " < " << vertexZPos << " < " << fVertZPos[1] << endl;
//   cout << "AliFemtoSpherocityEventCut:: VertexZErr: " << event->PrimVertCov()[4] << endl;

  // cout << "AliFemtoSpherocityEventCut:: MagneticField: " << event->MagneticField() << endl;
  // cout << "AliFemtoSpherocityEventCut:: IsCollisionCandidate: " << event->IsCollisionCandidate() << endl;
  // cout << "AliFemtoSpherocityEventCut:: TriggerCluster: " << event->TriggerCluster() << endl;
  // cout << "AliFemtoSpherocityEventCut:: fSelectTrigger: " << fSelectTrigger << endl;
  // cout << "AliFemtoSpherocityEventCut:: " << endl;
  bool goodEvent = ((fEventMult[0] <= mult) && (mult <= fEventMult[1])
                    && (fVertZPos[0] < vertexZPos) && (vertexZPos < fVertZPos[1])
                    && (fPsiEP[0] < epvzero) && (epvzero < fPsiEP[1])
                    && ((!fAcceptBadVertex) || (event->ZDCParticipants() > 1.0))
                    && ((!fSelectTrigger) || (event->TriggerCluster() == fSelectTrigger))
                   );

  // cout << "AliFemtoSpherocityEventCut:: goodEvent" <<goodEvent << endl;

  goodEvent ? fNEventsPassed++ : fNEventsFailed++ ;
  // cout << "AliFemtoSpherocityEventCut:: return : " << goodEvent << endl;
//     (fAcceptBadVertex || (event->PrimVertCov()[4] > -1000.0)) &&

  return goodEvent;
}
//------------------------------
AliFemtoString AliFemtoSpherocityEventCut::Report()
{
  // Prepare report
  AliFemtoString report("AliFemtoSpherocityEventCut Report:\n");
  report += Form("\nMultiplicity:\t %d-%d",fEventMult[0],fEventMult[1]);
  report += Form("\nVertex Z-position:\t %E-%E",fVertZPos[0],fVertZPos[1]);
  report += Form("\nNumber of events which passed:\t%ld  Number which failed:\t%ld",fNEventsPassed,fNEventsFailed);

  return report;
}
void AliFemtoSpherocityEventCut::SetAcceptBadVertex(bool b)
{
  fAcceptBadVertex = b;
}
bool AliFemtoSpherocityEventCut::GetAcceptBadVertex()
{
  return fAcceptBadVertex;
}
