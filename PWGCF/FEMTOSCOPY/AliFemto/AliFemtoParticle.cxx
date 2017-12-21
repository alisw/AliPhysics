///
/// \file AliFemtoParticle.cxx
///

#include "AliFemtoKink.h"
#include "AliFemtoParticle.h"
#include "AliFemtoXi.h"

double AliFemtoParticle::fgPrimPimPar0 = 9.05632e-01;
double AliFemtoParticle::fgPrimPimPar1 = -2.26737e-01;
double AliFemtoParticle::fgPrimPimPar2 = -1.03922e-01;
double AliFemtoParticle::fgPrimPipPar0 = 9.09616e-01;
double AliFemtoParticle::fgPrimPipPar1 = -9.00511e-02;
double AliFemtoParticle::fgPrimPipPar2 = -6.02940e-02;
double AliFemtoParticle::fgPrimPmPar0 = 0.;
double AliFemtoParticle::fgPrimPmPar1 = 0.;
double AliFemtoParticle::fgPrimPmPar2 = 0.;
double AliFemtoParticle::fgPrimPpPar0 = 0.;
double AliFemtoParticle::fgPrimPpPar1 = 0.;
double AliFemtoParticle::fgPrimPpPar2 = 0.;

int TpcLocalTransform(AliFmThreeVectorD &xgl,
                      int &iSector,
                      int &iPadrow,
                      float &xlocal,
                      double &ttPhi);


//_____________________
AliFemtoParticle::AliFemtoParticle() :
  fTrack(NULL),
  fV0(NULL),
  fKink(NULL),
  fXi(NULL),
  fFourMomentum(),
  fHelix(),
  fHiddenInfo(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  fHelixV0Pos(),
  fTpcV0PosEntrancePoint(),
  fTpcV0PosExitPoint(),
  fHelixV0Neg(),
  fTpcV0NegEntrancePoint(),
  fTpcV0NegExitPoint()
{
  // Default constructor
  std::fill_n(fPurity, 6, 0.0);
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoParticle &aParticle):
  fTrack(NULL),
  fV0(NULL),
  fKink(NULL),
  fXi(NULL),
  fFourMomentum(aParticle.fFourMomentum),
  fHelix(aParticle.fHelix),
  fHiddenInfo(NULL),
  fPrimaryVertex(aParticle.fPrimaryVertex),
  fSecondaryVertex(aParticle.fSecondaryVertex),
  fHelixV0Pos(aParticle.fHelixV0Pos),
  fTpcV0PosEntrancePoint(aParticle.fTpcV0PosEntrancePoint),
  fTpcV0PosExitPoint(aParticle.fTpcV0PosExitPoint),
  fHelixV0Neg(aParticle.fHelixV0Neg),
  fTpcV0NegEntrancePoint(aParticle.fTpcV0NegEntrancePoint),
  fTpcV0NegExitPoint(aParticle.fTpcV0NegExitPoint)
{
  // Copy constructor
  memcpy(fPurity, aParticle.fPurity, sizeof(fPurity));
  if (aParticle.fTrack)
    fTrack = new AliFemtoTrack(*aParticle.fTrack);
  if (aParticle.fV0)
    fV0 = new AliFemtoV0(*aParticle.fV0);
  if (aParticle.fKink)
    fKink = new AliFemtoKink(*aParticle.fKink);
  if (aParticle.fXi)
    fXi = new AliFemtoXi(*aParticle.fXi);

  if (aParticle.fHiddenInfo)
    fHiddenInfo = aParticle.HiddenInfo()->Clone();

  //CalculatePurity();
}
//_____________________
AliFemtoParticle::~AliFemtoParticle()
{
  if(fTrack)  delete fTrack;
  if(fV0)     delete fV0;
  if(fKink)   delete fKink;
  if(fXi)     delete fXi;
  if(fHiddenInfo) delete fHiddenInfo;
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoTrack *const hbtTrack, const double &mass):
  fTrack(new AliFemtoTrack(*hbtTrack)),
  fV0(NULL),
  fKink(NULL),
  fXi(NULL),
  fFourMomentum(::sqrt(hbtTrack->P().Mag2() + mass*mass), hbtTrack->P()),
  fHelix(hbtTrack->Helix()),
  fHiddenInfo(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  fHelixV0Pos(),
  fTpcV0PosEntrancePoint(),
  fTpcV0PosExitPoint(),
  fHelixV0Neg(),
  fTpcV0NegEntrancePoint(),
  fTpcV0NegExitPoint()
{
  // Constructor from normal track
  /* TO JA ODZNACZYLEM NIE WIEM DLACZEGO
  CalculateTpcExitAndEntrancePoints(&fHelix,&fPrimaryVertex,
            &fSecondaryVertex,
            &fNominalTpcEntrancePoint,
            &fNominalTpcExitPoint,
            &mNominalPosSample[0],
            &fZ[0],
            &fU[0],
            &fSect[0]);
  */
  CalculatePurity();
  if (hbtTrack->ValidHiddenInfo()) {
    fHiddenInfo = hbtTrack->GetHiddenInfo()->Clone();
  }
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoV0 *const hbtV0, const double &mass):
  fTrack(NULL),
  fV0(new AliFemtoV0(*hbtV0)),
  fKink(NULL),
  fXi(NULL),
  fFourMomentum(::sqrt(hbtV0->MomV0().Mag2() + mass*mass), hbtV0->MomV0()),
  fHelix(),
  fHiddenInfo(NULL),
  fPrimaryVertex(hbtV0->PrimaryVertex()),
  fSecondaryVertex(hbtV0->DecayVertexV0()),
  fHelixV0Pos(hbtV0->HelixPos()),
  fTpcV0PosEntrancePoint(),
  fTpcV0PosExitPoint(),
  fHelixV0Neg(hbtV0->HelixNeg()),
  fTpcV0NegEntrancePoint(),
  fTpcV0NegExitPoint()
{
  // Constructor from V0

  if (hbtV0->ValidHiddenInfo()) {
    fHiddenInfo = hbtV0->GetHiddenInfo()->Clone();
  }
  for (int ip = 0; ip < 6; ip++) fPurity[ip] = 0.0;
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoKink *const hbtKink, const double &mass) :
  fTrack(NULL),
  fV0(NULL),
  fKink(new AliFemtoKink(*hbtKink)),
  fXi(NULL),
  fFourMomentum(::sqrt(hbtKink->Parent().P().Mag2() + mass*mass), hbtKink->Parent().P()),
  fHelix(),
//   fNominalTpcExitPoint(0),
//   fNominalTpcEntrancePoint(0),
  fHiddenInfo(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  fHelixV0Pos(),
  fTpcV0PosEntrancePoint(),
  fTpcV0PosExitPoint(),
  fHelixV0Neg(),
  fTpcV0NegEntrancePoint(),
  fTpcV0NegExitPoint()
{
  // Constructor from Kink
  for (int ip = 0; ip < 6; ip++) fPurity[ip] = 0.0;
}

//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoXi *const hbtXi, const double &mass) :
  fTrack(NULL),
  fV0(NULL),
  fKink(NULL),
  fXi(new AliFemtoXi(*hbtXi)),
  fFourMomentum(::sqrt(hbtXi->MomXi().Mag2() + mass*mass), hbtXi->MomXi()),
  fHelix(),
//   fNominalTpcExitPoint(0),
//   fNominalTpcEntrancePoint(0),
  fHiddenInfo(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  fHelixV0Pos(),
  fTpcV0PosEntrancePoint(),
  fTpcV0PosExitPoint(),
  fHelixV0Neg(),
  fTpcV0NegEntrancePoint(),
  fTpcV0NegExitPoint()
{
  // Constructor from Xi
  for (int ip = 0; ip < 6; ip++) fPurity[ip] = 0.0;
}
//_____________________
AliFemtoParticle &AliFemtoParticle::operator=(const AliFemtoParticle &aParticle)
{
  // assignment operator
  if (this == &aParticle) return *this;

  delete fTrack;
  fTrack = NULL;
  if (aParticle.fTrack) {
    fTrack = new AliFemtoTrack(*aParticle.fTrack);
    CalculatePurity();
  }

  delete fV0;
  fV0 = NULL;
  if (aParticle.fV0)
    fV0 = new AliFemtoV0(*aParticle.fV0);

  delete fKink;
  fKink = NULL;
  if (aParticle.fKink)
    fKink = new AliFemtoKink(*aParticle.fKink);

  delete fXi;
  if (aParticle.fXi)
    fXi = new AliFemtoXi(*aParticle.fXi);

  fFourMomentum = aParticle.fFourMomentum;
  fHelix = aParticle.fHelix;

  fPrimaryVertex = aParticle.fPrimaryVertex;
  fSecondaryVertex = aParticle.fSecondaryVertex;

    //   for (int iter=0; iter<11; iter++)
    //     fNominalPosSample[iter] = aParticle.fNominalPosSample[iter];

    //   if (fTpcV0NegPosSample) delete fTpcV0NegPosSample;
    //   if (aParticle.fTpcV0NegPosSample) {
    //     fTpcV0NegPosSample = (AliFemtoThreeVector *) malloc(sizeof(AliFemtoThreeVector) * 11);
    //     for (int iter=0; iter<11; iter++)
    //       fTpcV0NegPosSample[iter] = aParticle.fTpcV0NegPosSample[iter];
    //   }

    //   if (fV0NegZ) delete fV0NegZ;
    //   if (aParticle.fV0NegZ) {
    //     fV0NegZ = (float *) malloc(sizeof(float) * 45);
    //     for (int iter=0; iter<11; iter++)
    //       fV0NegZ[iter] = aParticle.fV0NegZ[iter];
    //   }
    //   if (fV0NegU) delete fV0NegU;
    //   if (aParticle.fV0NegU) {
    //     fV0NegU = (float *) malloc(sizeof(float) * 45);
    //     for (int iter=0; iter<11; iter++)
    //       fV0NegU[iter] = aParticle.fV0NegU[iter];
    //   }
    //   if (fV0NegSect) delete fV0NegSect;
    //   if (aParticle.fV0NegSect) {
    //     fV0NegSect = (int *) malloc(sizeof(int) * 45);
    //     for (int iter=0; iter<11; iter++)
    //       fV0NegSect[iter] = aParticle.fV0NegSect[iter];
    //   }

//     if (fHiddenInfo) delete fHiddenInfo;
//     if(aParticle.fHiddenInfo){
//       fHiddenInfo= aParticle.fHiddenInfo->Clone();
//     }

    //   fNominalTpcEntrancePoint = aParticle.fNominalTpcEntrancePoint;
    //   fNominalTpcExitPoint     = aParticle.fNominalTpcExitPoint;

  delete fHiddenInfo;
  if (aParticle.fHiddenInfo)
    fHiddenInfo = aParticle.HiddenInfo()->Clone();

  for (int iter = 0; iter < 6; iter++)
    fPurity[iter] = aParticle.fPurity[iter];

  fHelixV0Pos = aParticle.fHelixV0Pos;
  fTpcV0PosEntrancePoint = aParticle.fTpcV0PosEntrancePoint;
  fTpcV0PosExitPoint = aParticle.fTpcV0PosExitPoint;
  fHelixV0Neg = aParticle.fHelixV0Neg;
  fTpcV0NegEntrancePoint = aParticle.fTpcV0NegEntrancePoint;
  fTpcV0NegExitPoint = aParticle.fTpcV0NegExitPoint;

  return *this;
}
// //_____________________
// const AliFemtoThreeVector& AliFemtoParticle::NominalTpcExitPoint() const{
//   // in future, may want to calculate this "on demand" only, sot this routine may get more sophisticated
//   // for now, we calculate Exit and Entrance points upon instantiation
//   return fNominalTpcExitPoint;
// }
// //_____________________
// const AliFemtoThreeVector& AliFemtoParticle::NominalTpcEntrancePoint() const{
//   // in future, may want to calculate this "on demand" only, sot this routine may get more sophisticated
//   // for now, we calculate Exit and Entrance points upon instantiation
//   return fNominalTpcEntrancePoint;
// }
//_____________________
void AliFemtoParticle::CalculatePurity()
{
  // Calculate additional parameterized purity

  double tPt = fFourMomentum.Perp();
  // pi -
  fPurity[0] = fgPrimPimPar0 * (1. - exp((tPt - fgPrimPimPar1) / fgPrimPimPar2));
  fPurity[0] *= fTrack->PidProbPion();
  // pi+
  fPurity[1] = fgPrimPipPar0 * (1. - exp((tPt - fgPrimPipPar1) / fgPrimPipPar2));
  fPurity[1] *= fTrack->PidProbPion();
  // K-
  fPurity[2] = fTrack->PidProbKaon();
  // K+
  fPurity[3] = fTrack->PidProbKaon();
  // pbar
  fPurity[4] = fTrack->PidProbProton();
  // p
  fPurity[5] = fTrack->PidProbProton();
}

double AliFemtoParticle::GetPionPurity()
{
  // Get full pion purity
  if (fTrack->Charge() > 0)
    return fPurity[1];
  else
    return fPurity[0];
}
double AliFemtoParticle::GetKaonPurity()
{
  // Get full kaon purity
  if (fTrack->Charge() > 0)
    return fPurity[3];
  else
    return fPurity[2];
}
double AliFemtoParticle::GetProtonPurity()
{
  // Get full proton purity
  if (fTrack->Charge() > 0)
    return fPurity[5];
  else
    return fPurity[4];
}

// void AliFemtoParticle::CalculateTpcExitAndEntrancePoints(AliFmPhysicalHelixD* tHelix,
//                   AliFemtoThreeVector*  PrimVert,
//                   AliFemtoThreeVector*  SecVert,
//                   AliFemtoThreeVector* tmpTpcEntrancePoint,
//                   AliFemtoThreeVector* tmpTpcExitPoint,
//                   AliFemtoThreeVector* tmpPosSample,
//                   float* tmpZ,
//                   float* tmpU,
//                   int* tmpSect){
//   // this calculates the exit point of a secondary track,
//   // either through the endcap or through the Outer Field Cage
//   // We assume the track to start at tHelix.origin-PrimaryVertex
//   // it also calculates the entrance point of the secondary track,
//   // which is the point at which it crosses the
//   // inner field cage
//   //  static AliFemtoThreeVector ZeroVec(0.,0.,0.);
//   AliFemtoThreeVector tZeroVec(0.,0.,0.);
// //   tZeroVec.SetX(tHelix->origin().x()-PrimVert->x());
// //   tZeroVec.SetY(tHelix->origin().y()-PrimVert->y());
// //   tZeroVec.SetZ(tHelix->origin().z()-PrimVert->z());
//   tZeroVec.SetX(SecVert->x()-PrimVert->x());
//   tZeroVec.SetY(SecVert->y()-PrimVert->y());
//   tZeroVec.SetZ(SecVert->z()-PrimVert->z());
//   double dip, curv, phase;
//   int h;
//   curv = tHelix->Curvature();
//   dip  = tHelix->DipAngle();
//   phase= tHelix->Phase();
//   h    = tHelix->H();

//   AliFmHelixD hel(curv,dip,phase,tZeroVec,h);

//   pairD candidates;
//   double sideLength;  // this is how much length to go to leave through sides of TPC
//   double endLength;  // this is how much length to go to leave through endcap of TPC
//   // figure out how far to go to leave through side...
//   candidates = hel.PathLength(200.0);  // bugfix MAL jul00 - 200cm NOT 2cm
//   sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

//   static AliFemtoThreeVector tWestEnd(0.,0.,200.);  // bugfix MAL jul00 - 200cm NOT 2cm
//   static AliFemtoThreeVector tEastEnd(0.,0.,-200.); // bugfix MAL jul00 - 200cm NOT 2cm
//   static AliFemtoThreeVector tEndCapNormal(0.,0.,1.0);

//   endLength = hel.PathLength(tWestEnd,tEndCapNormal);
//   if (endLength < 0.0) endLength = hel.PathLength(tEastEnd,tEndCapNormal);

//   if (endLength < 0.0) cout <<
//       "AliFemtoParticle::CalculateTpcExitAndEntrancePoints(): "
//                             << "Hey -- I cannot find an exit point out endcaps" << endl;
//   // OK, firstExitLength will be the shortest way out of the detector...
//   double firstExitLength = (endLength < sideLength) ? endLength : sideLength;
//   // now then, let's return the POSITION at which particle leaves TPC...
//   *tmpTpcExitPoint = hel.At(firstExitLength);
//   // Finally, calculate the position at which the track crosses the inner field cage
//   candidates = hel.PathLength(50.0);  // bugfix MAL jul00 - 200cm NOT 2cm

//   sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
//   //  cout << "sideLength 2 ="<<sideLength << endl;
//   *tmpTpcEntrancePoint = hel.At(sideLength);
//   // This is the secure way !
//   if (IsNaN(tmpTpcEntrancePoint->x()) ||
//       IsNaN(tmpTpcEntrancePoint->y()) ||
//       IsNaN(tmpTpcEntrancePoint->z()) ){
//     cout << "tmpTpcEntrancePoint NAN"<< endl;
//     cout << "tmpNominalTpcEntrancePoint = " <<tmpTpcEntrancePoint<< endl;
//     tmpTpcEntrancePoint->SetX(-9999.);
//     tmpTpcEntrancePoint->SetY(-9999.);
//     tmpTpcEntrancePoint->SetZ(-9999.);
//   }

//   if (IsNaN(tmpTpcExitPoint->x()) ||
//       IsNaN(tmpTpcExitPoint->y()) ||
//       IsNaN(tmpTpcExitPoint->z()) ) {
// //     cout << "tmpTpcExitPoint NAN Set at (-9999,-9999,-9999)"<< endl;
// //     cout << "tmpTpcExitPoint X= " <<tmpTpcExitPoint->x()<< endl;
// //     cout << "tmpTpcExitPoint Y= " <<tmpTpcExitPoint->y()<< endl;
// //     cout << "tmpTpcExitPoint Z= " <<tmpTpcExitPoint->z()<< endl;
//     tmpTpcExitPoint->SetX(-9999.);
//     tmpTpcExitPoint->SetY(-9999.);
//     tmpTpcExitPoint->SetZ(-9999.);
//   }


// //   if (IsNaN(tmpTpcExitPoint->x())) *tmpTpcExitPoint = AliFemtoThreeVector(-9999.,-9999.,-9999);
// //   if (IsNaN(tmpTpcEntrancetPoint->x())) *tmpTpcEntrancePoint = AliFemtoThreeVector(-9999.,-9999.,-9999);
//   //  cout << "tmpTpcEntrancePoint"<<*tmpTpcEntrancePoint << endl;

//   // 03Oct00 - mal.  OK, let's try something a little more
//   // along the lines of NA49 and E895 strategy.
//   // calculate the "nominal" position at N radii (say N=11)
//   // within the TPC, and for a pair cut
//   // use the average separation of these N
//   int irad = 0;
//   candidates = hel.PathLength(50.0);
//   sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
//   while (irad<11 && !IsNaN(sideLength)){
//     float radius = 50.0 + irad*15.0;
//     candidates = hel.PathLength(radius);
//     sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
//     tmpPosSample[irad] = hel.At(sideLength);
//     if(IsNaN(tmpPosSample[irad].x()) ||
//        IsNaN(tmpPosSample[irad].y()) ||
//        IsNaN(tmpPosSample[irad].z())
//        ){
//       cout << "tmpPosSample for radius=" << radius << " NAN"<< endl;
//       cout << "tmpPosSample=(" <<tmpPosSample[irad]<<")"<< endl;
//       tmpPosSample[irad] =  AliFemtoThreeVector(-9999.,-9999.,-9999);
//     }
//     irad++;
//     if (irad<11){
//       float radius = 50.0 + irad*15.0;
//       candidates = hel.PathLength(radius);
//       sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
//     }
//    }
//    for (int i = irad; i<11; i++)
//      {
//        tmpPosSample[i] =  AliFemtoThreeVector(-9999.,-9999.,-9999);
//      }

//   static float tRowRadius[45] = {60,64.8,69.6,74.4,79.2,84,88.8,93.6,98.8,
//         104,109.2,114.4,119.6,127.195,129.195,131.195,
//         133.195,135.195,137.195,139.195,141.195,
//         143.195,145.195,147.195,149.195,151.195,
//         153.195,155.195,157.195,159.195,161.195,
//         163.195,165.195,167.195,169.195,171.195,
//         173.195,175.195,177.195,179.195,181.195,
//         183.195,185.195,187.195,189.195};
//   int tRow,tSect,tOutOfBound;
//   double tLength,tPhi;
//   float tU;
//   AliFemtoThreeVector tPoint;
//   AliFmThreeVectorD tn(0,0,0);
//   AliFmThreeVectorD tr(0,0,0);
//   int ti =0;
//   // test to enter the loop
//   candidates =  hel.PathLength(tRowRadius[ti]);
//   tLength = (candidates.first > 0) ? candidates.first : candidates.second;
//   if (IsNaN(tLength)){
//     cout <<"tLength Init tmp NAN" << endl;
//     cout <<"padrow number= "<<ti << "not reached" << endl;
//     cout << "*** DO NOT ENTER THE LOOP***" << endl;
//     tmpSect[ti]=-1;//sector
//   }
//   // end test
//   while(ti<45 && !IsNaN(tLength)){
//     candidates =  hel.PathLength(tRowRadius[ti]);
//     tLength = (candidates.first > 0) ? candidates.first : candidates.second;
//     if (IsNaN(tLength)){
//       cout <<"tLength loop 1st NAN" << endl;
//       cout <<"padrow number=  " << ti << " not reached" << endl;
//       cout << "*** THIS IS AN ERROR SHOULDN'T  LOOP ***" << endl;
//       tmpSect[ti]=-1;//sector
//     }
//     tPoint = hel.At(tLength);
//     // Find which sector it is on
//     TpcLocalTransform(tPoint,tmpSect[ti],tRow,tU,tPhi);
//     if (IsNaN(tmpSect[ti])){
//       cout <<"***ERROR tmpSect"<< endl;
//     }
//     if (IsNaN(tRow)){
//       cout <<"***ERROR tRow"<< endl;
//     }
//     if (IsNaN(tU)){
//       cout <<"***ERROR tU"<< endl;
//     }
//     if (IsNaN(tPhi)){
//       cout <<"***ERROR tPhi"<< endl;
//     }
//     // calculate crossing plane
//     tn.SetX(cos(tPhi));
//     tn.SetY(sin(tPhi));
//     tr.SetX(tRowRadius[ti]*cos(tPhi));
//     tr.SetY(tRowRadius[ti]*sin(tPhi));
//     // find crossing point
//     tLength = hel.PathLength(tr,tn);
//     if (IsNaN(tLength)){
//       cout <<"tLength loop 2nd  NAN" << endl;
//       cout <<"padrow number=  " << ti << " not reached" << endl;
//       tmpSect[ti]=-2;//sector
//     }
//     tPoint = hel.At(tLength);
//     tmpZ[ti] = tPoint.z();
//     tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
//     if (IsNaN(tSect)){
//       cout <<"***ERROR tSect 2"<< endl;
//     }
//     if (IsNaN(tRow)){
//       cout <<"***ERROR tRow 2"<< endl;
//     }
//     if (IsNaN(tmpU[ti])){
//       cout <<"***ERROR tmpU[ti] 2"<< endl;
//     }
//     if (IsNaN(tPhi)){
//       cout <<"***ERROR tPhi 2 "<< endl;
//     }
//     if(tOutOfBound || (tmpSect[ti] == tSect && tRow!=(ti+1))){
//       tmpSect[ti]=-2;
//       //   cout << "missed once"<< endl;
//     }
//     else{
//       if(tmpSect[ti] != tSect){
//  // Try again on the other sector
//  tn.SetX(cos(tPhi));
//  tn.SetY(sin(tPhi));
//  tr.SetX(tRowRadius[ti]*cos(tPhi));
//  tr.SetY(tRowRadius[ti]*sin(tPhi));
//  // find crossing point
//  tLength = hel.PathLength(tr,tn);
//  tPoint = hel.At(tLength);
//  if (IsNaN(tLength)){
//    cout <<"tLength loop 3rd NAN" << endl;
//    cout <<"padrow number=  "<< ti << " not reached" << endl;
//    tmpSect[ti]=-1;//sector
//  }
//  tmpZ[ti] = tPoint.z();
//  tmpSect[ti] = tSect;
//  tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
//  if (IsNaN(tSect)){
//    cout <<"***ERROR tSect 3"<< endl;
//  }
//  if (IsNaN(tRow)){
//    cout <<"***ERROR tRow 3"<< endl;
//  }
//  if (IsNaN(tmpU[ti])){
//    cout <<"***ERROR tmpU[ti] 3"<< endl;
//  }
//  if (IsNaN(tPhi)){
//    cout <<"***ERROR tPhi 3 "<< endl;
//  }
//  if(tOutOfBound || tSect!= tmpSect[ti] || tRow!=(ti+1)){
//    tmpSect[ti]=-1;
//  }
//       }
//     }
//     if (IsNaN(tmpSect[ti])){
//       cout << "*******************ERROR***************************" << endl;
//       cout <<"AliFemtoParticle--Fctn tmpSect=" << tmpSect[ti] << endl;
//       cout << "*******************ERROR***************************" << endl;
//     }
//     if (IsNaN(tmpU[ti])){
//       cout << "*******************ERROR***************************" << endl;
//       cout <<"AliFemtoParticle--Fctn tmpU=" << tmpU[ti] << endl;
//       cout << "*******************ERROR***************************" << endl;
//     }
//     if (IsNaN(tmpZ[ti])){
//       cout << "*******************ERROR***************************" << endl;
//       cout <<"AliFemtoParticle--Fctn tmpZ=" << tmpZ[ti] << endl;
//       cout << "*******************ERROR***************************" << endl;
//     }
//     // If padrow ti not reached all other beyond are not reached
//     // in this case Set sector to -1
//     if (tmpSect[ti]==-1){
//       for (int tj=ti; tj<45;tj++){
//  tmpSect[tj] = -1;
//  ti=45;
//       }
//     }
//     ti++;
//     if (ti<45){
//       candidates =  hel.PathLength(tRowRadius[ti]);
//       tLength = (candidates.first > 0) ? candidates.first : candidates.second;}
//   }
// }
//_____________________
const AliFemtoThreeVector &AliFemtoParticle::TpcV0PosExitPoint() const
{
  return fTpcV0PosExitPoint;
}
//_____________________
const AliFemtoThreeVector &AliFemtoParticle::TpcV0PosEntrancePoint() const
{
  return fTpcV0PosEntrancePoint;
}
//______________________
const AliFemtoThreeVector &AliFemtoParticle::TpcV0NegExitPoint() const
{
  return fTpcV0NegExitPoint;
}
//_____________________
const AliFemtoThreeVector &AliFemtoParticle::TpcV0NegEntrancePoint() const
{
  return fTpcV0NegEntrancePoint;
}
//______________________
