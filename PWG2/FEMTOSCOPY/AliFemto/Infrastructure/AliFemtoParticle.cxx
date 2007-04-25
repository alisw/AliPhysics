/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   Particle objects are part of the PicoEvent, which is what is
 *   stored in the EventMixingBuffers
 *   A Track object gets converted to a Particle object if it
 *   passes the ParticleCut of an Analysis
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.22  2003/09/02 17:58:32  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.21  2003/05/07 15:30:43  magestro
 * Fixed bug related to finding merged hits (commit for Fabrice)
 *
 * Revision 1.20  2003/01/14 09:41:26  renault
 * changes on average separation calculation, hit shared finder and memory optimisation
 * for Z,U and Sectors variables.
 *
 * Revision 1.19  2002/12/12 17:01:49  kisiel
 * Hidden Information handling and purity calculation
 *
 * Revision 1.18  2002/11/19 23:36:00  renault
 * Enable calculation of exit/entrance separation for V0 daughters
 *
 * Revision 1.17  2001/12/14 23:11:30  fretiere
 * Add class HitMergingCut. Add class fabricesPairCut = HitMerginCut + pair purity cuts. Add TpcLocalTransform function which convert to local tpc coord (not pretty). Modify AliFemtoTrack, AliFemtoParticle, AliFemtoHiddenInfo, AliFemtoPair to handle the hit information and cope with my code
 *
 * Revision 1.16  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 * Revision 1.15  2001/04/03 21:04:36  kisiel
 *
 *
 *   Changes needed to make the Theoretical code
 *   work. The main code is the ThCorrFctn directory.
 *   The most visible change is the addition of the
 *   HiddenInfo to AliFemtoPair.
 *
 * Revision 1.14  2000/10/05 23:09:05  lisa
 * Added kT-dependent radii to mixed-event simulator AND implemented AverageSeparation Cut and CorrFctn
 *
 * Revision 1.13  2000/08/31 22:31:31  laue
 * AliFemtoAnalysis: output changed (a little bit less)
 * AliFemtoEvent: new version, members for reference mult added
 * AliFemtoIOBinary: new IO for new AliFemtoEvent version
 * AliFemtoTypes: TTree typedef to AliFemtoTTree added
 * AliFemtoVertexAnalysis: overflow and underflow added
 *
 * Revision 1.12  2000/07/23 13:52:56  laue
 * NominalExitPoint set to (-9999.,-9999.-9999.) if helix.at()
 * returns nan (not a number).
 *
 * Revision 1.11  2000/07/19 17:18:48  laue
 * Calculation of Entrance and Exit point added in AliFemtoParticle constructor
 *
 * Revision 1.10  2000/07/17 20:03:17  lisa
 * Implemented tools for addressing and assessing trackmerging
 *
 * Revision 1.9  2000/07/16 21:38:23  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers fTrack,fV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.8  2000/05/03 17:44:43  laue
 * AliFemtoEvent, AliFemtoTrack & AliFemtoV0 declared friend to AliFemtoIOBinary
 * AliFemtoParticle updated for V0 pos,neg track Id
 *
 * Revision 1.7  2000/04/03 16:21:51  laue
 * some include files changed
 * Multi track cut added
 *
 * Revision 1.6  1999/12/11 15:58:29  lisa
 * Add vertex decay position datum and accessor to AliFemtoParticle to allow pairwise cuts on seperation of V0s
 *
 * Revision 1.5  1999/09/17 22:38:02  lisa
 * first full integration of V0s into AliFemto framework
 *
 * Revision 1.4  1999/09/01 19:04:53  lisa
 * update Particle class AND add parity cf and Randys Coulomb correction
 *
 * Revision 1.3  1999/07/06 22:33:23  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.2  1999/06/29 17:50:27  fisyak
 * formal changes to account new StEvent, does not complie yet
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Infrastructure/AliFemtoParticle.h"
//#include "math_constants.h"
#ifdef __CC5__
  #include <math.h>
#else
  #include <cmath>
#endif


#include "TMath.h"
using namespace TMath;

double AliFemtoParticle::fPrimPimPar0= 9.05632e-01;
double AliFemtoParticle::fPrimPimPar1= -2.26737e-01;
double AliFemtoParticle::fPrimPimPar2= -1.03922e-01;
double AliFemtoParticle::fPrimPipPar0= 9.09616e-01;
double AliFemtoParticle::fPrimPipPar1= -9.00511e-02;
double AliFemtoParticle::fPrimPipPar2= -6.02940e-02;
double AliFemtoParticle::fPrimPmPar0= 0.;
double AliFemtoParticle::fPrimPmPar1= 0.;
double AliFemtoParticle::fPrimPmPar2= 0.;
double AliFemtoParticle::fPrimPpPar0= 0.;
double AliFemtoParticle::fPrimPpPar1= 0.;
double AliFemtoParticle::fPrimPpPar2= 0.;

int TpcLocalTransform(AliFmThreeVectorD& xgl, 
		      int& iSector, 
		      int& iPadrow, 
		      float& xlocal,
		      double& ttPhi);


//_____________________
AliFemtoParticle::AliFemtoParticle() : fTrack(0), fV0(0), fKink(0), fHiddenInfo(0) {
  /* no-op for default */
  //  cout << "Created particle " << this << endl;
}
//_____________________
AliFemtoParticle::~AliFemtoParticle(){
  //  cout << "Issuing delete for AliFemtoParticle." << endl;

  if (fTrack) delete fTrack;
  if (fV0) {
    delete[] fTpcV0NegPosSample;
    delete[] fV0NegZ;
    delete[] fV0NegU;
    delete[] fV0NegSect;
    delete fV0;
  }
  if (fKink) delete fKink;
  //  cout << "Trying to delete HiddenInfo: " << fHiddenInfo << endl;
  if (fHiddenInfo) 
    {
      //      cout << "Deleting HiddenInfo." << endl;
      delete fHiddenInfo;
    }
  //  cout << "Deleted particle " << this << endl;
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoTrack* const hbtTrack,const double& mass) : fTrack(0), fV0(0), fKink(0), fHiddenInfo(0) {
  
  
  // I know there is a better way to do this...
  fTrack = new AliFemtoTrack(*hbtTrack);
  AliFemtoThreeVector temp = hbtTrack->P();
  fFourMomentum.setVect(temp);
  double ener = ::sqrt(temp.mag2()+mass*mass);
  fFourMomentum.setE(ener);
//  fMap[0] = hbtTrack->TopologyMap(0);
 // fMap[1] = hbtTrack->TopologyMap(1);
 // fNhits = hbtTrack->NHits();
  fHelix = hbtTrack->Helix();
  //CalculateNominalTpcExitAndEntrancePoints();

 
  fPrimaryVertex.setX(0.);
  fPrimaryVertex.setY(0.);
  fPrimaryVertex.setZ(0.);
  fSecondaryVertex.setX(0.);
  fSecondaryVertex.setY(0.);
  fSecondaryVertex.setZ(0.);
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
  // ***
  fHiddenInfo= 0;
  if(hbtTrack->ValidHiddenInfo()){
    fHiddenInfo= hbtTrack->getHiddenInfo()->getParticleHiddenInfo()->clone();
  }
  // ***
  //  cout << "Created particle " << this << endl;

}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoV0* const hbtV0,const double& mass) : fTrack(0), fV0(0), fKink(0), fHiddenInfo(0) {
  fV0 = new AliFemtoV0(*hbtV0);
 //fMap[0]= 0;
  //fMap[1]= 0;
  // I know there is a better way to do this...
  AliFemtoThreeVector temp = hbtV0->momV0();
  fFourMomentum.setVect(temp);
  double ener = ::sqrt(temp.mag2()+mass*mass);
  fFourMomentum.setE(ener);
  // Calculating TpcEntrancePoint for Positive V0 daugther
  fPrimaryVertex = hbtV0->primaryVertex();
  fSecondaryVertex = hbtV0->decayVertexV0();
  fHelixV0Pos = hbtV0->HelixPos();

  fTpcV0NegPosSample = new AliFemtoThreeVector[45];//for V0Neg
  fV0NegZ = new float[45];//for V0Neg
  fV0NegU = new float[45];//for V0Neg
  fV0NegSect = new int[45];//for V0Neg
  CalculateTpcExitAndEntrancePoints(&fHelixV0Pos,&fPrimaryVertex,
				    &fSecondaryVertex,
				    &fTpcV0PosEntrancePoint,
				    &fTpcV0PosExitPoint,
				    &fNominalPosSample[0],
				    &fZ[0],
				    &fU[0],&fSect[0]);
  fHelixV0Neg = hbtV0->HelixNeg();

  CalculateTpcExitAndEntrancePoints(&fHelixV0Neg,
				    &fPrimaryVertex,
				    &fSecondaryVertex,
				    &fTpcV0NegEntrancePoint,
				    &fTpcV0NegExitPoint,
				    &fTpcV0NegPosSample[0],
				    &fV0NegZ[0],
				    &fV0NegU[0],&fV0NegSect[0]);

  // ***
  fHiddenInfo= 0;
  if(hbtV0->ValidHiddenInfo()){
    fHiddenInfo= hbtV0->getHiddenInfo()->clone();
  }
  // ***
}
//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoKink* const hbtKink,const double& mass) : fTrack(0), fV0(0), fHiddenInfo(0) {
  fKink = new AliFemtoKink(*hbtKink);
 // fMap[0]= 0;
  //fMap[1]= 0;
  // I know there is a better way to do this...
  AliFemtoThreeVector temp = hbtKink->Parent().P();
  fFourMomentum.setVect(temp);
  double ener = ::sqrt(temp.mag2()+mass*mass);
  fFourMomentum.setE(ener);
}

//_____________________
AliFemtoParticle::AliFemtoParticle(const AliFemtoXi* const hbtXi, const double& mass)  {
  fXi = new AliFemtoXi(*hbtXi);
 // fMap[0]= 0;
  //fMap[1]= 0;
  AliFemtoThreeVector temp;// = hbtXi->mMofXi;
  fFourMomentum.setVect(temp);
  double ener = ::sqrt(temp.mag2()+mass*mass);
  fFourMomentum.setE(ener);
  fHiddenInfo = 0;
}
//_____________________
const AliFemtoThreeVector& AliFemtoParticle::NominalTpcExitPoint() const{
  // in future, may want to calculate this "on demand" only, sot this routine may get more sophisticated
  // for now, we calculate Exit and Entrance points upon instantiation
  return fNominalTpcExitPoint;
}
//_____________________
const AliFemtoThreeVector& AliFemtoParticle::NominalTpcEntrancePoint() const{
  // in future, may want to calculate this "on demand" only, sot this routine may get more sophisticated
  // for now, we calculate Exit and Entrance points upon instantiation
  return fNominalTpcEntrancePoint;
}
//_____________________
void AliFemtoParticle::CalculatePurity(){
  double tPt = fFourMomentum.perp();
  // pi -
  fPurity[0] = fPrimPimPar0*(1.-exp((tPt-fPrimPimPar1)/fPrimPimPar2));
  fPurity[0] *= fTrack->PidProbPion();
  // pi+
  fPurity[1] = fPrimPipPar0*(1.-exp((tPt-fPrimPipPar1)/fPrimPipPar2));
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
  if (fTrack->Charge()>0)
    return fPurity[1];
  else
    return fPurity[0];
}
double AliFemtoParticle::GetKaonPurity()
{
  if (fTrack->Charge()>0)
    return fPurity[3];
  else
    return fPurity[2];
}
double AliFemtoParticle::GetProtonPurity()
{
  if (fTrack->Charge()>0)
    return fPurity[5];
  else
    return fPurity[4];
}

void AliFemtoParticle::CalculateTpcExitAndEntrancePoints(AliFmPhysicalHelixD* tHelix,
						       AliFemtoThreeVector*  PrimVert,
						       AliFemtoThreeVector*  SecVert,
						       AliFemtoThreeVector* tmpTpcEntrancePoint,
						       AliFemtoThreeVector* tmpTpcExitPoint,
						       AliFemtoThreeVector* tmpPosSample,
						       float* tmpZ,
						       float* tmpU,
						       int* tmpSect){
  // this calculates the exit point of a secondary track, 
  // either through the endcap or through the Outer Field Cage
  // We assume the track to start at tHelix.origin-PrimaryVertex
  // it also calculates the entrance point of the secondary track, 
  // which is the point at which it crosses the
  // inner field cage
  //  static AliFemtoThreeVector ZeroVec(0.,0.,0.);
  AliFemtoThreeVector ZeroVec(0.,0.,0.);
//   ZeroVec.setX(tHelix->origin().x()-PrimVert->x());
//   ZeroVec.setY(tHelix->origin().y()-PrimVert->y());
//   ZeroVec.setZ(tHelix->origin().z()-PrimVert->z());
  ZeroVec.setX(SecVert->x()-PrimVert->x());
  ZeroVec.setY(SecVert->y()-PrimVert->y());
  ZeroVec.setZ(SecVert->z()-PrimVert->z());
  double dip, curv, phase;
  int h;
  curv = tHelix->curvature();
  dip  = tHelix->dipAngle();
  phase= tHelix->phase();
  h    = tHelix->h();
  
  AliFmHelixD hel(curv,dip,phase,ZeroVec,h);

  pairD candidates;
  double sideLength;  // this is how much length to go to leave through sides of TPC
  double endLength;  // this is how much length to go to leave through endcap of TPC
  // figure out how far to go to leave through side...
  candidates = hel.pathLength(200.0);  // bugfix MAL jul00 - 200cm NOT 2cm
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

  static AliFemtoThreeVector WestEnd(0.,0.,200.);  // bugfix MAL jul00 - 200cm NOT 2cm
  static AliFemtoThreeVector EastEnd(0.,0.,-200.); // bugfix MAL jul00 - 200cm NOT 2cm
  static AliFemtoThreeVector EndCapNormal(0.,0.,1.0);

  endLength = hel.pathLength(WestEnd,EndCapNormal);
  if (endLength < 0.0) endLength = hel.pathLength(EastEnd,EndCapNormal);

  if (endLength < 0.0) cout << 
			 "AliFemtoParticle::CalculateTpcExitAndEntrancePoints(): "
                            << "Hey -- I cannot find an exit point out endcaps" << endl;
  // OK, firstExitLength will be the shortest way out of the detector...
  double firstExitLength = (endLength < sideLength) ? endLength : sideLength;
  // now then, let's return the POSITION at which particle leaves TPC...
  *tmpTpcExitPoint = hel.at(firstExitLength);
  // Finally, calculate the position at which the track crosses the inner field cage
  candidates = hel.pathLength(50.0);  // bugfix MAL jul00 - 200cm NOT 2cm

  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
  //  cout << "sideLength 2 ="<<sideLength << endl;
  *tmpTpcEntrancePoint = hel.at(sideLength);
  // This is the secure way !  
  if (IsNaN(tmpTpcEntrancePoint->x()) || 
      IsNaN(tmpTpcEntrancePoint->y()) || 
      IsNaN(tmpTpcEntrancePoint->z()) ){ 
    cout << "tmpTpcEntrancePoint NAN"<< endl; 
    cout << "tmpNominalTpcEntrancePoint = " <<tmpTpcEntrancePoint<< endl;
    tmpTpcEntrancePoint->setX(-9999.);
    tmpTpcEntrancePoint->setY(-9999.);
    tmpTpcEntrancePoint->setZ(-9999.);
  } 
    
  if (IsNaN(tmpTpcExitPoint->x()) || 
      IsNaN(tmpTpcExitPoint->y()) || 
      IsNaN(tmpTpcExitPoint->z()) ) {
//     cout << "tmpTpcExitPoint NAN set at (-9999,-9999,-9999)"<< endl; 
//     cout << "tmpTpcExitPoint X= " <<tmpTpcExitPoint->x()<< endl;
//     cout << "tmpTpcExitPoint Y= " <<tmpTpcExitPoint->y()<< endl;
//     cout << "tmpTpcExitPoint Z= " <<tmpTpcExitPoint->z()<< endl;
    tmpTpcExitPoint->setX(-9999.);
    tmpTpcExitPoint->setY(-9999.);
    tmpTpcExitPoint->setZ(-9999.);
  }


//   if (IsNaN(tmpTpcExitPoint->x())) *tmpTpcExitPoint = AliFemtoThreeVector(-9999.,-9999.,-9999); 
//   if (IsNaN(tmpTpcEntrancetPoint->x())) *tmpTpcEntrancePoint = AliFemtoThreeVector(-9999.,-9999.,-9999); 
  //  cout << "tmpTpcEntrancePoint"<<*tmpTpcEntrancePoint << endl;

  // 03Oct00 - mal.  OK, let's try something a little more 
  // along the lines of NA49 and E895 strategy.
  // calculate the "nominal" position at N radii (say N=11) 
  // within the TPC, and for a pair cut
  // use the average separation of these N
  int irad = 0;
  candidates = hel.pathLength(50.0);
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
  while (irad<11 && !IsNaN(sideLength)){ 
    float radius = 50.0 + irad*15.0;
    candidates = hel.pathLength(radius);
    sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    tmpPosSample[irad] = hel.at(sideLength);
    if(IsNaN(tmpPosSample[irad].x()) ||
       IsNaN(tmpPosSample[irad].y()) ||
       IsNaN(tmpPosSample[irad].z()) 
       ){
      cout << "tmpPosSample for radius=" << radius << " NAN"<< endl; 
      cout << "tmpPosSample=(" <<tmpPosSample[irad]<<")"<< endl;
      tmpPosSample[irad] =  AliFemtoThreeVector(-9999.,-9999.,-9999);
    }
    irad++;
    if (irad<11){
      float radius = 50.0 + irad*15.0;
      candidates = hel.pathLength(radius);
      sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    }
   }
   for (int i = irad; i<11; i++)
     {
       tmpPosSample[i] =  AliFemtoThreeVector(-9999.,-9999.,-9999);   
     }

  static float tRowRadius[45] = {60,64.8,69.6,74.4,79.2,84,88.8,93.6,98.8, 
				 104,109.2,114.4,119.6,127.195,129.195,131.195,
				 133.195,135.195,137.195,139.195,141.195,
				 143.195,145.195,147.195,149.195,151.195,
				 153.195,155.195,157.195,159.195,161.195,
				 163.195,165.195,167.195,169.195,171.195,
				 173.195,175.195,177.195,179.195,181.195,
				 183.195,185.195,187.195,189.195};
  int tRow,tSect,tOutOfBound;
  double tLength,tPhi;
  float tU;
  AliFemtoThreeVector tPoint;
  AliFmThreeVectorD tn(0,0,0);
  AliFmThreeVectorD tr(0,0,0);
  int ti =0;
  // test to enter the loop
  candidates =  hel.pathLength(tRowRadius[ti]);
  tLength = (candidates.first > 0) ? candidates.first : candidates.second;
  if (IsNaN(tLength)){
    cout <<"tLength Init tmp NAN" << endl;
    cout <<"padrow number= "<<ti << "not reached" << endl;
    cout << "*** DO NOT ENTER THE LOOP***" << endl;
    tmpSect[ti]=-1;//sector
  }
  // end test
  while(ti<45 && !IsNaN(tLength)){
    candidates =  hel.pathLength(tRowRadius[ti]);
    tLength = (candidates.first > 0) ? candidates.first : candidates.second;
    if (IsNaN(tLength)){
      cout <<"tLength loop 1st NAN" << endl;
      cout <<"padrow number=  " << ti << " not reached" << endl;
      cout << "*** THIS IS AN ERROR SHOULDN'T  LOOP ***" << endl;
      tmpSect[ti]=-1;//sector
    }
    tPoint = hel.at(tLength);
    // Find which sector it is on
    TpcLocalTransform(tPoint,tmpSect[ti],tRow,tU,tPhi);
    if (IsNaN(tmpSect[ti])){
      cout <<"***ERROR tmpSect"<< endl; 
    }
    if (IsNaN(tRow)){
      cout <<"***ERROR tRow"<< endl;
    }
    if (IsNaN(tU)){
      cout <<"***ERROR tU"<< endl;
    }
    if (IsNaN(tPhi)){
      cout <<"***ERROR tPhi"<< endl;
    }  
    // calculate crossing plane
    tn.setX(cos(tPhi));
    tn.setY(sin(tPhi));       
    tr.setX(tRowRadius[ti]*cos(tPhi));
    tr.setY(tRowRadius[ti]*sin(tPhi));
    // find crossing point
    tLength = hel.pathLength(tr,tn); 
    if (IsNaN(tLength)){
      cout <<"tLength loop 2nd  NAN" << endl;
      cout <<"padrow number=  " << ti << " not reached" << endl;
      tmpSect[ti]=-2;//sector
    }
    tPoint = hel.at(tLength);
    tmpZ[ti] = tPoint.z();
    tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
    if (IsNaN(tSect)){
      cout <<"***ERROR tSect 2"<< endl; 
    }
    if (IsNaN(tRow)){
      cout <<"***ERROR tRow 2"<< endl;
    }
    if (IsNaN(tmpU[ti])){
      cout <<"***ERROR tmpU[ti] 2"<< endl;
    }
    if (IsNaN(tPhi)){
      cout <<"***ERROR tPhi 2 "<< endl;
    }  
    if(tOutOfBound || (tmpSect[ti] == tSect && tRow!=(ti+1))){
      tmpSect[ti]=-2;
      //	  cout << "missed once"<< endl;
    }
    else{
      if(tmpSect[ti] != tSect){
	// Try again on the other sector
	tn.setX(cos(tPhi));
	tn.setY(sin(tPhi));       
	tr.setX(tRowRadius[ti]*cos(tPhi));
	tr.setY(tRowRadius[ti]*sin(tPhi));
	// find crossing point
	tLength = hel.pathLength(tr,tn);
	tPoint = hel.at(tLength);
	if (IsNaN(tLength)){
	  cout <<"tLength loop 3rd NAN" << endl;
	  cout <<"padrow number=  "<< ti << " not reached" << endl;
	  tmpSect[ti]=-1;//sector
	}
	tmpZ[ti] = tPoint.z();
	tmpSect[ti] = tSect;
	tOutOfBound = TpcLocalTransform(tPoint,tSect,tRow,tmpU[ti],tPhi);
	if (IsNaN(tSect)){
	  cout <<"***ERROR tSect 3"<< endl; 
	}
	if (IsNaN(tRow)){
	  cout <<"***ERROR tRow 3"<< endl;
	}
	if (IsNaN(tmpU[ti])){
	  cout <<"***ERROR tmpU[ti] 3"<< endl;
	}
	if (IsNaN(tPhi)){
	  cout <<"***ERROR tPhi 3 "<< endl;
	}  
	if(tOutOfBound || tSect!= tmpSect[ti] || tRow!=(ti+1)){
	  tmpSect[ti]=-1;
	}
      }
    }
    if (IsNaN(tmpSect[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"AliFemtoParticle--Fctn tmpSect=" << tmpSect[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    if (IsNaN(tmpU[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"AliFemtoParticle--Fctn tmpU=" << tmpU[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    if (IsNaN(tmpZ[ti])){
      cout << "*******************ERROR***************************" << endl;
      cout <<"AliFemtoParticle--Fctn tmpZ=" << tmpZ[ti] << endl;
      cout << "*******************ERROR***************************" << endl;
    }
    // If padrow ti not reached all other beyond are not reached
    // in this case set sector to -1
    if (tmpSect[ti]==-1){
      for (int tj=ti; tj<45;tj++){
	tmpSect[tj] = -1;
	ti=45;
      }
    }
    ti++;
    if (ti<45){
      candidates =  hel.pathLength(tRowRadius[ti]);
      tLength = (candidates.first > 0) ? candidates.first : candidates.second;}
  }
}
//_____________________
const AliFemtoThreeVector& AliFemtoParticle::TpcV0PosExitPoint() const{
  return fTpcV0PosExitPoint;
}
//_____________________
const AliFemtoThreeVector& AliFemtoParticle::TpcV0PosEntrancePoint() const{
  return fTpcV0PosEntrancePoint;
}
//______________________
const AliFemtoThreeVector& AliFemtoParticle::TpcV0NegExitPoint() const{
  return fTpcV0NegExitPoint;
}
//_____________________
const AliFemtoThreeVector& AliFemtoParticle::TpcV0NegEntrancePoint() const{
  return fTpcV0NegEntrancePoint;
}
//______________________
