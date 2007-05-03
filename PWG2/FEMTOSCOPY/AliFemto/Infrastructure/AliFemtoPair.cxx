/***************************************************************************
 *
 * $Id: AliFemtoPair.cc,v 1.23
 *
 * Author: Brian Laziuk, Yale University
 *         slightly modified by Mike Lisa
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    the Pair object is passed to the PairCuts for verification, and
 *    then to the AddRealPair and AddMixedPair methods of the
 *    Correlation Functions
 *
 ***************************************************************************
 * Revision 1.23  2002/09/25 19:23:25  rcwells
 * Added const to emissionAngle()
 *
 * Revision 1.22  2002/04/22 22:48:11  laue
 * corrected calculation of opening angle 
 **
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.27  2003/09/02 17:58:32  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.26  2003/01/31 19:57:15  magestro
 * Cleared up simple compiler warnings on i386_linux24
 *
 * Revision 1.25  2003/01/14 09:44:08  renault
 * corrections on average separation calculation for tracks which doesn't cross
 * all 45 padrows.
 *
 * Revision 1.24  2002/11/19 23:33:10  renault
 * Enable average separation calculation for all combinaisons of
 * V0 daughters and tracks
 *
 * Revision 1.21  2002/02/28 14:18:36  rcwells
 * Added emissionAngle function to AliFemtoPair
 *
 * Revision 1.20  2001/12/14 23:11:30  fretiere
 * Add class HitMergingCut. Add class fabricesPairCut = HitMerginCut + pair purity cuts. Add TpcLocalTransform function which convert to local tpc coord (not pretty). Modify AliFemtoTrack, AliFemtoParticle, AliFemtoHiddenInfo, AliFemtoPair to handle the hit information and cope with my code
 *
 * Revision 1.19  2001/04/25 18:05:09  perev
 * HPcorrs
 *
 * Revision 1.18  2001/04/03 21:04:36  kisiel
 *
 *
 *   Changes needed to make the Theoretical code
 *   work. The main code is the ThCorrFctn directory.
 *   The most visible change is the addition of the
 *   HiddenInfo to AliFemtoPair.
 *
 * Revision 1.17  2001/03/28 22:35:20  flierl
 * changes and bugfixes in qYKP*
 * add pairrapidity
 *
 * Revision 1.16  2001/02/15 19:23:00  rcwells
 * Fixed sign in qSideCMS
 *
 * Revision 1.15  2001/01/22 22:56:41  laue
 * Yano-Koonin-Podgoretskii Parametrisation added
 *
 * Revision 1.14  2000/12/11 21:44:30  rcwells
 * Corrected qSideCMS function
 *
 * Revision 1.13  2000/10/26 16:09:16  lisa
 * Added OpeningAngle PairCut class and method to AliFemtoPair
 *
 * Revision 1.12  2000/10/05 23:09:05  lisa
 * Added kT-dependent radii to mixed-event simulator AND implemented AverageSeparation Cut and CorrFctn
 *
 * Revision 1.11  2000/07/17 20:03:16  lisa
 * Implemented tools for addressing and assessing trackmerging
 *
 * Revision 1.10  2000/04/04 16:27:03  rcwells
 * Removed an errant cout in AliFemtoPair.cc
 *
 * Revision 1.9  2000/04/04 16:13:09  lisa
 * AliFemtoPair:quality() now returns normalized value (and so is double) and add a CorrFctn which looks at quality()
 *
 * Revision 1.8  2000/04/03 22:09:12  rcwells
 * Add member function ... quality().
 *
 * Revision 1.7  2000/02/13 21:13:33  lisa
 * changed ambiguous AliFemtoPair::fourMomentum() to fourMomentumSum() and fourMomentumDiff() and fixed related bug in QvecCorrFctn
 *
 * Revision 1.6  1999/07/29 16:16:34  lisa
 * Selemons upgrade of AliFemtoPair class
 *
 * Revision 1.5  1999/07/22 18:49:10  lisa
 * Implement idea of Fabrice to not create and delete AliFemtoPair all the time
 *
 * Revision 1.4  1999/07/12 18:57:05  lisa
 * fixed small bug in fourMomentum method of AliFemtoPair
 *
 * Revision 1.3  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.2  1999/06/29 17:50:27  fisyak
 * formal changes to account new StEvent, does not complie yet
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Infrastructure/AliFemtoPair.h"

double AliFemtoPair::fMaxDuInner = .8;
double AliFemtoPair::fMaxDzInner = 3.;
double AliFemtoPair::fMaxDuOuter = 1.4;
double AliFemtoPair::fMaxDzOuter = 3.2;


AliFemtoPair::AliFemtoPair() :
  fTrack1(0), fTrack2(0),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  kStarCalc(0),
  fNonIdParNotCalculatedGlobal(0),
  fMergingParNotCalculated(0),
  fWeightedAvSep(0),
  fFracOfMergedRow(0),
  fClosestRowAtDCA(0),
  fMergingParNotCalculatedTrkV0Pos(0),
  fFracOfMergedRowTrkV0Pos(0),
  fClosestRowAtDCATrkV0Pos(0),
  fMergingParNotCalculatedTrkV0Neg(0),
  fFracOfMergedRowTrkV0Neg(0),
  fClosestRowAtDCATrkV0Neg(0),
  fMergingParNotCalculatedV0PosV0Neg(0),
  fFracOfMergedRowV0PosV0Neg(0),
  fClosestRowAtDCAV0PosV0Neg(0),
  fMergingParNotCalculatedV0NegV0Pos(0),
  fFracOfMergedRowV0NegV0Pos(0),
  fClosestRowAtDCAV0NegV0Pos(0),
  fMergingParNotCalculatedV0PosV0Pos(0),
  fFracOfMergedRowV0PosV0Pos(0),
  fClosestRowAtDCAV0PosV0Pos(0),
  fMergingParNotCalculatedV0NegV0Neg(0),
  fFracOfMergedRowV0NegV0Neg(0),
  fClosestRowAtDCAV0NegV0Neg(0)
{
  fTrack1 = 0;
  fTrack2 = 0;
  setDefaultHalfFieldMergingPar();
}

AliFemtoPair::AliFemtoPair(AliFemtoParticle* a, AliFemtoParticle* b)
  : fTrack1(a), fTrack2(b),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  kStarCalc(0),
  fNonIdParNotCalculatedGlobal(0),
  fMergingParNotCalculated(0),
  fWeightedAvSep(0),
  fFracOfMergedRow(0),
  fClosestRowAtDCA(0),
  fMergingParNotCalculatedTrkV0Pos(0),
  fFracOfMergedRowTrkV0Pos(0),
  fClosestRowAtDCATrkV0Pos(0),
  fMergingParNotCalculatedTrkV0Neg(0),
  fFracOfMergedRowTrkV0Neg(0),
  fClosestRowAtDCATrkV0Neg(0),
  fMergingParNotCalculatedV0PosV0Neg(0),
  fFracOfMergedRowV0PosV0Neg(0),
  fClosestRowAtDCAV0PosV0Neg(0),
  fMergingParNotCalculatedV0NegV0Pos(0),
  fFracOfMergedRowV0NegV0Pos(0),
  fClosestRowAtDCAV0NegV0Pos(0),
  fMergingParNotCalculatedV0PosV0Pos(0),
  fFracOfMergedRowV0PosV0Pos(0),
  fClosestRowAtDCAV0PosV0Pos(0),
  fMergingParNotCalculatedV0NegV0Neg(0),
  fFracOfMergedRowV0NegV0Neg(0),
  fClosestRowAtDCAV0NegV0Neg(0)
{ 
  setDefaultHalfFieldMergingPar();
}

void AliFemtoPair::setDefaultHalfFieldMergingPar(){
  fMaxDuInner = 3;
  fMaxDzInner = 4.;
  fMaxDuOuter = 4.;
  fMaxDzOuter = 6.;
}
void AliFemtoPair::setDefaultFullFieldMergingPar(){
  fMaxDuInner = 0.8;
  fMaxDzInner = 3.;
  fMaxDuOuter = 1.4;
  fMaxDzOuter = 3.2;
}
void AliFemtoPair::setMergingPar(double aMaxDuInner, double aMaxDzInner,
			      double aMaxDuOuter, double aMaxDzOuter){
  fMaxDuInner = aMaxDuInner;
  fMaxDzInner = aMaxDzInner;
  fMaxDuOuter = aMaxDuOuter;
  fMaxDzOuter = aMaxDzOuter;
};

AliFemtoPair::~AliFemtoPair() {/* no-op */}

AliFemtoPair::AliFemtoPair(const AliFemtoPair &aPair):
  fTrack1(0), fTrack2(0),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  kStarCalc(0),
  fNonIdParNotCalculatedGlobal(0),
  fMergingParNotCalculated(0),
  fWeightedAvSep(0),
  fFracOfMergedRow(0),
  fClosestRowAtDCA(0),
  fMergingParNotCalculatedTrkV0Pos(0),
  fFracOfMergedRowTrkV0Pos(0),
  fClosestRowAtDCATrkV0Pos(0),
  fMergingParNotCalculatedTrkV0Neg(0),
  fFracOfMergedRowTrkV0Neg(0),
  fClosestRowAtDCATrkV0Neg(0),
  fMergingParNotCalculatedV0PosV0Neg(0),
  fFracOfMergedRowV0PosV0Neg(0),
  fClosestRowAtDCAV0PosV0Neg(0),
  fMergingParNotCalculatedV0NegV0Pos(0),
  fFracOfMergedRowV0NegV0Pos(0),
  fClosestRowAtDCAV0NegV0Pos(0),
  fMergingParNotCalculatedV0PosV0Pos(0),
  fFracOfMergedRowV0PosV0Pos(0),
  fClosestRowAtDCAV0PosV0Pos(0),
  fMergingParNotCalculatedV0NegV0Neg(0),
  fFracOfMergedRowV0NegV0Neg(0),
  fClosestRowAtDCAV0NegV0Neg(0)
{
  fTrack1 = aPair.fTrack1;
  fTrack2 = aPair.fTrack2;

  fNonIdParNotCalculated = aPair.fNonIdParNotCalculated;
  fDKSide = aPair.fDKSide;
  fDKOut = aPair.fDKOut;
  fDKLong = aPair.fDKLong;
  fCVK = aPair.fCVK;
  kStarCalc = aPair.kStarCalc;

  fNonIdParNotCalculatedGlobal = aPair.fNonIdParNotCalculatedGlobal;

  fMergingParNotCalculated = aPair.fMergingParNotCalculated;
  fWeightedAvSep = aPair.fWeightedAvSep;
  fFracOfMergedRow = aPair.fFracOfMergedRow;
  fClosestRowAtDCA = aPair.fClosestRowAtDCA;

  fMergingParNotCalculatedTrkV0Pos = aPair.fMergingParNotCalculatedTrkV0Pos;
  fFracOfMergedRowTrkV0Pos = aPair.fFracOfMergedRowTrkV0Pos;
  fClosestRowAtDCATrkV0Pos = aPair.fClosestRowAtDCATrkV0Pos;

  fMergingParNotCalculatedTrkV0Neg = aPair.fMergingParNotCalculatedTrkV0Neg;
  fFracOfMergedRowTrkV0Neg = aPair.fFracOfMergedRowTrkV0Neg;
  fClosestRowAtDCATrkV0Neg = aPair.fClosestRowAtDCATrkV0Neg;

  fMergingParNotCalculatedV0PosV0Neg = aPair.fMergingParNotCalculatedV0PosV0Neg;
  fFracOfMergedRowV0PosV0Neg = aPair.fFracOfMergedRowV0PosV0Neg;
  fClosestRowAtDCAV0PosV0Neg = aPair.fClosestRowAtDCAV0PosV0Neg;

  fMergingParNotCalculatedV0NegV0Pos = aPair.fMergingParNotCalculatedV0NegV0Pos;
  fFracOfMergedRowV0NegV0Pos = aPair.fFracOfMergedRowV0NegV0Pos;
  fClosestRowAtDCAV0NegV0Pos = aPair.fClosestRowAtDCAV0NegV0Pos;

  fMergingParNotCalculatedV0PosV0Pos = aPair.fMergingParNotCalculatedV0PosV0Pos;
  fFracOfMergedRowV0PosV0Pos = aPair.fFracOfMergedRowV0PosV0Pos;
  fClosestRowAtDCAV0PosV0Pos = aPair.fClosestRowAtDCAV0PosV0Pos;

  fMergingParNotCalculatedV0NegV0Neg = aPair.fMergingParNotCalculatedV0NegV0Neg;
  fFracOfMergedRowV0NegV0Neg = aPair.fFracOfMergedRowV0NegV0Neg;
  fClosestRowAtDCAV0NegV0Neg = aPair.fClosestRowAtDCAV0NegV0Neg;
}

AliFemtoPair& AliFemtoPair::operator=(const AliFemtoPair &aPair)
{
  if (this == &aPair)
    return *this;

  fTrack1 = aPair.fTrack1;
  fTrack2 = aPair.fTrack2;

  fNonIdParNotCalculated = aPair.fNonIdParNotCalculated;
  fDKSide = aPair.fDKSide;
  fDKOut = aPair.fDKOut;
  fDKLong = aPair.fDKLong;
  fCVK = aPair.fCVK;
  kStarCalc = aPair.kStarCalc;

  fNonIdParNotCalculatedGlobal = aPair.fNonIdParNotCalculatedGlobal;

  fMergingParNotCalculated = aPair.fMergingParNotCalculated;
  fWeightedAvSep = aPair.fWeightedAvSep;
  fFracOfMergedRow = aPair.fFracOfMergedRow;
  fClosestRowAtDCA = aPair.fClosestRowAtDCA;

  fMergingParNotCalculatedTrkV0Pos = aPair.fMergingParNotCalculatedTrkV0Pos;
  fFracOfMergedRowTrkV0Pos = aPair.fFracOfMergedRowTrkV0Pos;
  fClosestRowAtDCATrkV0Pos = aPair.fClosestRowAtDCATrkV0Pos;

  fMergingParNotCalculatedTrkV0Neg = aPair.fMergingParNotCalculatedTrkV0Neg;
  fFracOfMergedRowTrkV0Neg = aPair.fFracOfMergedRowTrkV0Neg;
  fClosestRowAtDCATrkV0Neg = aPair.fClosestRowAtDCATrkV0Neg;

  fMergingParNotCalculatedV0PosV0Neg = aPair.fMergingParNotCalculatedV0PosV0Neg;
  fFracOfMergedRowV0PosV0Neg = aPair.fFracOfMergedRowV0PosV0Neg;
  fClosestRowAtDCAV0PosV0Neg = aPair.fClosestRowAtDCAV0PosV0Neg;

  fMergingParNotCalculatedV0NegV0Pos = aPair.fMergingParNotCalculatedV0NegV0Pos;
  fFracOfMergedRowV0NegV0Pos = aPair.fFracOfMergedRowV0NegV0Pos;
  fClosestRowAtDCAV0NegV0Pos = aPair.fClosestRowAtDCAV0NegV0Pos;

  fMergingParNotCalculatedV0PosV0Pos = aPair.fMergingParNotCalculatedV0PosV0Pos;
  fFracOfMergedRowV0PosV0Pos = aPair.fFracOfMergedRowV0PosV0Pos;
  fClosestRowAtDCAV0PosV0Pos = aPair.fClosestRowAtDCAV0PosV0Pos;

  fMergingParNotCalculatedV0NegV0Neg = aPair.fMergingParNotCalculatedV0NegV0Neg;
  fFracOfMergedRowV0NegV0Neg = aPair.fFracOfMergedRowV0NegV0Neg;
  fClosestRowAtDCAV0NegV0Neg = aPair.fClosestRowAtDCAV0NegV0Neg;

  return *this;
}

//_________________
double AliFemtoPair::mInv() const
{
    double InvariantMass = abs(fTrack1->FourMomentum() + fTrack2->FourMomentum());
    return (InvariantMass);
}
//_________________
double AliFemtoPair::kT() const
{

  double  tmp = 
    (fTrack1->FourMomentum() + fTrack2->FourMomentum()).perp();
  tmp *= .5;

  return (tmp);
}
//_________________
double AliFemtoPair::rap() const
{
  // longitudinal pair rapidity : Y = 0.5 ::log( E1 + E2 + pz1 + pz2 / E1 + E2 - pz1 - pz2 )
  double  tmp = 0.5 * log (
			   (fTrack1->FourMomentum().e() + fTrack2->FourMomentum().e() + fTrack1->FourMomentum().z() + fTrack2->FourMomentum().z()) / 
			   (fTrack1->FourMomentum().e() + fTrack2->FourMomentum().e() - fTrack1->FourMomentum().z() - fTrack2->FourMomentum().z()) 
			   ) ;
  return (tmp);
}
//_________________
double AliFemtoPair::emissionAngle()const {
  double pxTotal = this->fourMomentumSum().x();
  double pyTotal = this->fourMomentumSum().y();
  double angle = atan2(pyTotal,pxTotal)*180.0/3.1415926536;
  if (angle<0.0) angle+=360.0;
  return angle;
}
//_________________
// get rid of ambiguously-named method fourMomentum() and replace it with
// fourMomentumSum() and fourMomentumDiff() - mal 13feb2000
AliFemtoLorentzVector AliFemtoPair::fourMomentumSum() const
{
  AliFemtoLorentzVector temp = fTrack1->FourMomentum()+fTrack2->FourMomentum();
  return temp;
}
AliFemtoLorentzVector AliFemtoPair::fourMomentumDiff() const
{
  AliFemtoLorentzVector temp = fTrack1->FourMomentum()-fTrack2->FourMomentum();
  return temp;
}
//__________________________________
// Yano-Koonin-Podgoretskii Parametrisation in CMS
void AliFemtoPair::qYKPCMS(double& qP, double& qT, double& q0) const
{
  ////
  // calculate momentum difference in source rest frame (= lab frame)
  ////
  AliFemtoLorentzVector l1 = fTrack1->FourMomentum() ;
  AliFemtoLorentzVector l2 = fTrack2->FourMomentum() ;
  AliFemtoLorentzVector  l ;
  // random ordering of the particles
  if ( rand()/(double)RAND_MAX > 0.50 )  
    { l = l1-l2 ; } 
  else 
    { l = l2-l1 ; } ;
  // fill momentum differences into return variables
  qP = l.z() ;
  qT = l.vect().perp() ;
  q0 = l.e() ;
}
//___________________________________
// Yano-Koonin-Podgoretskii Parametrisation in LCMS
void AliFemtoPair::qYKPLCMS(double& qP, double& qT, double& q0) const
{
  ////
  //  calculate momentum difference in LCMS : frame where pz1 + pz2 = 0
  ////
  AliFemtoLorentzVector l1 = fTrack1->FourMomentum() ;
  AliFemtoLorentzVector l2 = fTrack2->FourMomentum() ;
  // determine beta to LCMS
  double beta = (l1.z()+l2.z()) / (l1.e()+l2.e()) ;
  double beta2 =  beta*beta ;
  // unfortunately STAR Class lib knows only boost(particle) not boost(beta) :(
  // -> create particle with velocity beta and mass 1.0
  // actually this is : dummyPz = ::sqrt( (dummyMass*dummyMass*beta2) / (1-beta2) ) ; 
  double dummyPz = ::sqrt( (beta2) / (1-beta2) ) ;
  // boost in the correct direction
  if (beta>0.0) { dummyPz = -dummyPz; } ;
  // create dummy particle
  AliFemtoLorentzVector  l(0.0, 0.0, dummyPz) ; 
  double dummyMass = 1.0 ;
  l.setE(l.vect().massHypothesis(dummyMass) );
  // boost particles along the beam into a frame with velocity beta 
  AliFemtoLorentzVector l1boosted = l1.boost(l) ;
  AliFemtoLorentzVector l2boosted = l2.boost(l) ;
  // caculate the momentum difference with random ordering of the particle
  if ( rand()/(double)RAND_MAX >0.50)  
    { l = l1boosted-l2boosted ; } 
  else 
    { l = l2boosted-l1boosted ;} ;
  // fill momentum differences into return variables
  qP = l.z() ;
  qT = l.vect().perp() ;
  q0 = l.e() ;
}
//___________________________________
// Yano-Koonin-Podgoretskii Parametrisation in pair rest frame
void AliFemtoPair::qYKPPF(double& qP, double& qT, double& q0) const
{
  ////
  //  calculate momentum difference in pair rest frame : frame where (pz1 + pz2, py1 + py2, px1 + px2) = (0,0,0)
  ////
  AliFemtoLorentzVector l1 = fTrack1->FourMomentum() ;
  AliFemtoLorentzVector l2 = fTrack2->FourMomentum() ;
  // the center of gravity of the pair travels with l
  AliFemtoLorentzVector  l = l1 + l2 ; 
  l = -l ;
  l.setE(-l.e()) ;
  // boost particles  
  AliFemtoLorentzVector l1boosted = l1.boost(l) ;
  AliFemtoLorentzVector l2boosted = l2.boost(l) ;
  // caculate the momentum difference with random ordering of the particle
  if ( rand()/(double)RAND_MAX > 0.50)  
    { l = l1boosted-l2boosted ; } 
  else 
    { l = l2boosted-l1boosted ;} ;
  // fill momentum differences into return variables
  qP = l.z();
  qT = l.vect().perp();
  q0 = l.e();
}
//_________________
double AliFemtoPair::qOutCMS() const
{
    AliFemtoThreeVector tmp1 = fTrack1->FourMomentum().vect();
    AliFemtoThreeVector tmp2 = fTrack2->FourMomentum().vect();

    double dx = tmp1.x() - tmp2.x();
    double xt = tmp1.x() + tmp2.x();
    
    double dy = tmp1.y() - tmp2.y();
    double yt = tmp1.y() + tmp2.y();

    double k1 = (::sqrt(xt*xt+yt*yt));
    double k2 = (dx*xt+dy*yt);
    double tmp = k2/k1;
    return (tmp);
}
//_________________
double AliFemtoPair::qSideCMS() const
{
    AliFemtoThreeVector tmp1 = fTrack1->FourMomentum().vect();
    AliFemtoThreeVector tmp2 = fTrack2->FourMomentum().vect();

    double x1 = tmp1.x();  double y1 = tmp1.y();
    double x2 = tmp2.x();  double y2 = tmp2.y();

    double xt = x1+x2;  double yt = y1+y2;
    double k1 = ::sqrt(xt*xt+yt*yt);

    double tmp = 2.0*(x2*y1-x1*y2)/k1;
    return (tmp);
}

//_________________________
double AliFemtoPair::qLongCMS() const
{
    AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
    AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();

    double dz = tmp1.z() - tmp2.z();
    double zz = tmp1.z() + tmp2.z();

    double dt = tmp1.t() - tmp2.t();
    double tt = tmp1.t() + tmp2.t();

    double beta = zz/tt;
    double gamma = 1.0/::sqrt(1.0 - beta*beta);

    double temp = gamma*(dz - beta*dt);
    return (temp);
}

//________________________________
double AliFemtoPair::qOutPf() const
{
 AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
 AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();

    double dt = tmp1.t() - tmp2.t();
    double tt = tmp1.t() + tmp2.t();

    double xt = tmp1.x() + tmp2.x();
    double yt = tmp1.y() + tmp2.y();

    double k1 = ::sqrt(xt*xt + yt*yt);
    double bOut = k1/tt;
    double gOut = 1.0/::sqrt(1.0 - bOut*bOut);

    double temp = gOut*(this->qOutCMS() - bOut*dt);
    return (temp);
}

//___________________________________
double AliFemtoPair::qSidePf() const
{
 return(this->qSideCMS());
}

//___________________________________

double AliFemtoPair::qLongPf() const
{
 return(this->qLongCMS());
}

//___________________________________
double AliFemtoPair::qOutBf(double beta) const
{
 return(this->qOutCMS());
}

//___________________________________

double AliFemtoPair::qSideBf(double beta) const
{
 return(this->qSideCMS());
}

//___________________________________
double AliFemtoPair::qLongBf(double beta) const
{
    AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
    AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();

    double dz = tmp1.z() - tmp2.z();
    double dt = tmp1.t() + tmp2.t();

    double gamma = 1.0/::sqrt(1.0 - beta*beta);

    double temp = gamma*(dz - beta*dt);
    return (temp);
}

double AliFemtoPair::quality() const {
  unsigned long mapMask0 = 0xFFFFFF00;
  unsigned long mapMask1 = 0x1FFFFF;
  unsigned long padRow1To24Track1 = fTrack1->TopologyMap(0) & mapMask0;
  unsigned long padRow25To45Track1 = fTrack1->TopologyMap(1) & mapMask1;
  unsigned long padRow1To24Track2 = fTrack2->TopologyMap(0) & mapMask0;
  unsigned long padRow25To45Track2 = fTrack2->TopologyMap(1) & mapMask1;
  // AND logic
  unsigned long bothPads1To24 = padRow1To24Track1 & padRow1To24Track2;
  unsigned long bothPads25To45 = padRow25To45Track1 & padRow25To45Track2;
  // XOR logic
  unsigned long onePad1To24 = padRow1To24Track1 ^ padRow1To24Track2;
  unsigned long onePad25To45 = padRow25To45Track1 ^ padRow25To45Track2;
  unsigned long bitI;
  int ibits;
  int Quality = 0;
  double normQual = 0.0;
  int MaxQuality = fTrack1->NumberOfHits() + fTrack2->NumberOfHits();
  for (ibits=8;ibits<=31;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad1To24 & bitI ) {
      Quality++;
      continue;
    }
    else{
      if ( bothPads1To24 & bitI ) Quality--;
    }
  }
  for (ibits=0;ibits<=20;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad25To45 & bitI ) {
      Quality++;
      continue;
    }
    else{
      if ( bothPads25To45 & bitI ) Quality--;
    }
  }
  normQual = (double)Quality/( (double) MaxQuality );
  return ( normQual );

}

double AliFemtoPair::quality2() const {
  unsigned long mapMask0 = 0xFFFFFF00;
  unsigned long mapMask1 = 0x1FFFFF;
  unsigned long padRow1To24Track1 = fTrack1->TopologyMap(0) & mapMask0;
  unsigned long padRow25To45Track1 = fTrack1->TopologyMap(1) & mapMask1;
  unsigned long padRow1To24Track2 = fTrack2->TopologyMap(0) & mapMask0;
  unsigned long padRow25To45Track2 = fTrack2->TopologyMap(1) & mapMask1;

  // AND logic
  //unsigned long bothPads1To24 = padRow1To24Track1 & padRow1To24Track2;
  //unsigned long bothPads25To45 = padRow25To45Track1 & padRow25To45Track2;

  // XOR logic
  unsigned long onePad1To24 = padRow1To24Track1 ^ padRow1To24Track2;
  unsigned long onePad25To45 = padRow25To45Track1 ^ padRow25To45Track2;
  unsigned long bitI;
  int ibits;
  int Quality = 0;
  double normQual = 0.0;
  int MaxQuality = fTrack1->NumberOfHits() + fTrack2->NumberOfHits();
  for (ibits=8;ibits<=31;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad1To24 & bitI ) {
      Quality++;
      continue;
    }
    //else{
    //if ( bothPads1To24 & bitI ) Quality--;
    //}
  }
  for (ibits=0;ibits<=20;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad25To45 & bitI ) {
      Quality++;
      continue;
    }
    //else{
    //if ( bothPads25To45 & bitI ) Quality--;
    //}
  }
  normQual = (double)Quality/( (double) MaxQuality );
  return ( normQual );

}


double AliFemtoPair::NominalTpcExitSeparation() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcExitPoint() - fTrack2->NominalTpcExitPoint();
  return (diff.mag());
}

double AliFemtoPair::NominalTpcEntranceSeparation() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcEntrancePoint() - fTrack2->NominalTpcEntrancePoint();
  return (diff.mag());
}

double AliFemtoPair::NominalTpcAverageSeparation() const {
  AliFemtoThreeVector diff;
  double AveSep = 0.0;
  int ipt = 0;
  if (fTrack1->fNominalPosSample && fTrack2->fNominalPosSample){
  while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
	 fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
	 ipt<11
	 ){
    //  for (int ipt=0; ipt<11; ipt++){
    diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
    ipt++;
    AveSep += diff.mag();
  }
  AveSep = AveSep/(ipt+1.);
  return (AveSep);}
  else return -1;
}

double AliFemtoPair::OpeningAngle() const {
 return 57.296* fTrack1->FourMomentum().vect().angle( fTrack2->FourMomentum().vect() );
//   AliFemtoThreeVector p1 = fTrack1->FourMomentum().vect();
//   AliFemtoThreeVector p2 = fTrack2->FourMomentum().vect();
//   return 57.296*(p1.phi()-p2.phi());
//   //double dAngInv = 57.296*acos((p1.dot(p2))/(p1.mag()*p2.mag()));
//   //return (dAngInv);
}
//_________________


double AliFemtoPair::KStarFlipped() const {
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();

  AliFmThreeVectorD qwe = tP1.vect();
  qwe *= -1.; // flip it
  tP1.setVect(qwe);
  
  AliFemtoLorentzVector tSum = (tP1+fTrack2->FourMomentum());
  double tMass = abs(tSum);
  AliFmThreeVectorD tGammaBeta = (1./tMass)*tSum.vect(); 
  double tGamma = tSum.e()/tMass;
  AliFmThreeVectorD tLongMom  = ((tP1.vect()*tGammaBeta)/
			      (tGammaBeta*tGammaBeta))*tGammaBeta;
  AliFmLorentzVectorD tK(tGamma*tP1.e() - tP1.vect()*tGammaBeta,
		      tP1.vect() + (tGamma-1.)*tLongMom - tP1.e()*tGammaBeta);
//VP  tP1.vect() *= -1.; // unflip it
  return tK.vect().mag();
}

//double AliFemtoPair::CVK() const{
//const AliFemtoLorentzVector& tP1 = fTrack1->FourMomentum();
//AliFemtoLorentzVector tSum = (tP1+fTrack2->FourMomentum());
//double tMass = abs(tSum);
//AliFmThreeVectorD tGammaBeta = (1./tMass)*tSum.vect(); 
//double tGamma = tSum.e()/tMass;
//AliFmThreeVectorD tLongMom  = ((tP1.vect()*tGammaBeta)/
//		      (tGammaBeta*tGammaBeta))*tGammaBeta;
//AliFmLorentzVectorD tK(tGamma*tP1.e() - tP1.vect()*tGammaBeta,
//	      tP1.vect() + (tGamma-1.)*tLongMom - tP1.e()*tGammaBeta);
//return (tK.vect())*tGammaBeta/tK.vect().magnitude()/tGammaBeta.magnitude();
//}

double AliFemtoPair::CVKFlipped() const{
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  AliFmThreeVectorD qwe = tP1.vect();
  qwe *= -1.; // flip it
  tP1.setVect(qwe);
  
  AliFemtoLorentzVector tSum = (tP1+fTrack2->FourMomentum());
  double tMass = abs(tSum);
  AliFmThreeVectorD tGammaBeta = (1./tMass)*tSum.vect(); 
  double tGamma = tSum.e()/tMass;
  AliFmThreeVectorD tLongMom  = ((tP1.vect()*tGammaBeta)/
			      (tGammaBeta*tGammaBeta))*tGammaBeta;
  AliFmLorentzVectorD tK(tGamma*tP1.e() - tP1.vect()*tGammaBeta,
		      tP1.vect() + (tGamma-1.)*tLongMom - tP1.e()*tGammaBeta);
//VP  tP1.vect() *= -1.; // unflip it
  return (tK.vect())*tGammaBeta/tGamma;
}

double AliFemtoPair::pInv() const{
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  AliFemtoLorentzVector tP2 = fTrack2->FourMomentum();
  double tP = (tP1.px()+tP2.px())*(tP1.px()+tP2.px())+
              (tP1.py()+tP2.py())*(tP1.py()+tP2.py())+
              (tP1.pz()+tP2.pz())*(tP1.pz()+tP2.pz())-
              (tP1.e() -tP2.e() )*(tP1.e() -tP2.e() );
  return ::sqrt(fabs(tP));
}

double AliFemtoPair::qInvFlippedXY() const{
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  tP1.setX(-1.*tP1.x());
  tP1.setY(-1.*tP1.y());
  AliFemtoLorentzVector tDiff = (tP1-fTrack2->FourMomentum());
  return ( -1.* tDiff.m());
}

void AliFemtoPair::calcNonIdPar() const{ // fortran like function! faster?
  fNonIdParNotCalculated=0;
  double px1 = fTrack1->FourMomentum().vect().x();
  double py1 = fTrack1->FourMomentum().vect().y();
  double pz1 = fTrack1->FourMomentum().vect().z();
  double pE1  = fTrack1->FourMomentum().e();
  double Particle1Mass = ::sqrt(pE1*pE1 - px1*px1 - py1*py1 - pz1*pz1);
  double px2 = fTrack2->FourMomentum().vect().x();
  double py2 = fTrack2->FourMomentum().vect().y();
  double pz2 = fTrack2->FourMomentum().vect().z();
  double pE2  = fTrack2->FourMomentum().e();
  double Particle2Mass = ::sqrt(pE2*pE2 - px2*px2 - py2*py2 - pz2*pz2);

  double Px = px1+px2;
  double Py = py1+py2;
  double Pz = pz1+pz2;
  double PE = pE1+pE2;
      
  double Ptrans = Px*Px + Py*Py;
  double Mtrans = PE*PE - Pz*Pz;
  double Pinv =   ::sqrt(Mtrans - Ptrans);
  Mtrans = ::sqrt(Mtrans);
  Ptrans = ::sqrt(Ptrans);
	
  double QinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
    (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);

  double Q = (Particle1Mass*Particle1Mass - Particle2Mass*Particle2Mass)/Pinv;
  Q = sqrt ( Q*Q - QinvL);
	  
  kStarCalc = Q/2;

  // ad 1) go to LCMS
  double beta = Pz/PE;
  double gamma = PE/Mtrans;
	    
  double pz1L = gamma * (pz1 - beta * pE1);
  double pE1L = gamma * (pE1 - beta * pz1);
  
  // fill histogram for beam projection ( z - axis )
  fDKLong = pz1L;

  // ad 2) rotation px -> Pt
  double px1R = (px1*Px + py1*Py)/Ptrans;
  double py1R = (-px1*Py + py1*Px)/Ptrans;
  
  //fill histograms for side projection ( y - axis )
  fDKSide = py1R;

  // ad 3) go from LCMS to CMS
  beta = Ptrans/Mtrans;
  gamma = Mtrans/Pinv;
  
  double px1C = gamma * (px1R - beta * pE1L);
  
  // fill histogram for out projection ( x - axis )
  fDKOut  = px1C;

  fCVK = (fDKOut*Ptrans + fDKLong*Pz)/kStarCalc/::sqrt(Ptrans*Ptrans+Pz*Pz);
}


/*void AliFemtoPair::calcNonIdParGlobal() const{ // fortran like function! faster?
  fNonIdParNotCalculatedGlobal=0;
  double px1 = fTrack1->Track()->PGlobal().x();
  double py1 = fTrack1->Track()->PGlobal().y();
  double pz1 = fTrack1->Track()->PGlobal().z();
  double Particle1Mass =  fTrack1->FourMomentum().m2();
  double pE1  = ::sqrt(Particle1Mass + px1*px1 + py1*py1 + pz1*pz1);
  Particle1Mass = ::sqrt(Particle1Mass);

  double px2 = fTrack2->Track()->PGlobal().x();
  double py2 = fTrack2->Track()->PGlobal().y();
  double pz2 = fTrack2->Track()->PGlobal().z();
  double Particle2Mass =  fTrack2->FourMomentum().m2();
  double pE2  = ::sqrt(Particle2Mass + px2*px2 + py2*py2 + pz2*pz2);
  Particle2Mass = ::sqrt(Particle2Mass);

  double Px = px1+px2;
  double Py = py1+py2;
  double Pz = pz1+pz2;
  double PE = pE1+pE2;
      
  double Ptrans = Px*Px + Py*Py;
  double Mtrans = PE*PE - Pz*Pz;
  double Pinv =   ::sqrt(Mtrans - Ptrans);
  Mtrans = ::sqrt(Mtrans);
  Ptrans = ::sqrt(Ptrans);
	
  double QinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
    (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);

  double Q = (Particle1Mass*Particle1Mass - Particle2Mass*Particle2Mass)/Pinv;
  Q = sqrt ( Q*Q - QinvL);
	  
  kStarCalcGlobal = Q/2;

  // ad 1) go to LCMS
  double beta = Pz/PE;
  double gamma = PE/Mtrans;
	    
  double pz1L = gamma * (pz1 - beta * pE1);
  double pE1L = gamma * (pE1 - beta * pz1);
  
  // fill histogram for beam projection ( z - axis )
  fDKLongGlobal = pz1L;

  // ad 2) rotation px -> Pt
  double px1R = (px1*Px + py1*Py)/Ptrans;
  double py1R = (-px1*Py + py1*Px)/Ptrans;
  
  //fill histograms for side projection ( y - axis )
  fDKSideGlobal = py1R;

  // ad 3) go from LCMS to CMS
  beta = Ptrans/Mtrans;
  gamma = Mtrans/Pinv;
  
  double px1C = gamma * (px1R - beta * pE1L);
  
  // fill histogram for out projection ( x - axis )
  fDKOutGlobal  = px1C;

  fCVKGlobal = (fDKOutGlobal*Ptrans + fDKLongGlobal*Pz)/
    kStarCalcGlobal/::sqrt(Ptrans*Ptrans+Pz*Pz);
}*/



double AliFemtoPair::dcaInsideTpc() const{

  double tMinDist=NominalTpcEntranceSeparation();
  double tExit = NominalTpcExitSeparation();
  tMinDist = (tExit>tMinDist) ? tMinDist : tExit;
  double tInsideDist;
  //tMinDist = 999.;

  double rMin = 60.;
  double rMax = 190.;
  const AliFmPhysicalHelixD& tHelix1 = fTrack1->Helix();
  const AliFmPhysicalHelixD& tHelix2 = fTrack2->Helix();
  // --- One is a line and other one a helix
  //if (tHelix1.mSingularity != tHelix2.mSingularity) return -999.;
  // --- 2 lines : don't care right now
  //if (tHelix1.mSingularity)  return -999.;
  // --- 2 helix
  double dx = tHelix2.xcenter() - tHelix1.xcenter();
  double dy = tHelix2.ycenter() - tHelix1.ycenter();
  double dd = ::sqrt(dx*dx + dy*dy);
  double r1 = 1/tHelix1.curvature();
  double r2 = 1/tHelix2.curvature();
  double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);
    
  double x, y, r;
  double s;
  if (fabs(cosAlpha) < 1) {           // two solutions
    double sinAlpha = sin(acos(cosAlpha));
    x = tHelix1.xcenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
    y = tHelix1.ycenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
    r = ::sqrt(x*x+y*y);
    if( r > rMin &&  r < rMax && 
	fabs(atan2(y,x)-fTrack1->NominalTpcEntrancePoint().phi())< 0.5
	){ // first solution inside
      s = tHelix1.pathLength(x, y);
      tInsideDist=tHelix2.distance(tHelix1.at(s));
      if(tInsideDist<tMinDist) tMinDist = tInsideDist;
    }
    else{ 
      x = tHelix1.xcenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
      y = tHelix1.ycenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
      r = ::sqrt(x*x+y*y);
      if( r > rMin &&  r < rMax &&
	  fabs(atan2(y,x)-fTrack1->NominalTpcEntrancePoint().phi())< 0.5
	  ) {  // second solution inside
        s = tHelix1.pathLength(x, y);
        tInsideDist=tHelix2.distance(tHelix1.at(s));
        if(tInsideDist<tMinDist) tMinDist = tInsideDist;
      }     
    }
  }
  return tMinDist;
}

void AliFemtoPair::calcMergingPar() const{
  fMergingParNotCalculated=0;

  double tDu, tDz;
  int tN = 0;
  fFracOfMergedRow = 0.;
  fWeightedAvSep =0.;
  double tDist;
  double tDistMax = 200.;
  for(int ti=0 ; ti<45 ; ti++){
    if(fTrack1->fSect[ti]==fTrack2->fSect[ti] && fTrack1->fSect[ti]!=-1){
      tDu = fabs(fTrack1->fU[ti]-fTrack2->fU[ti]);
      tDz = fabs(fTrack1->fZ[ti]-fTrack2->fZ[ti]);
      tN++;
      if(ti<13){
	fFracOfMergedRow += (tDu<fMaxDuInner && tDz<fMaxDzInner);
	tDist = ::sqrt(tDu*tDu/fMaxDuInner/fMaxDuInner+
		     tDz*tDz/fMaxDzInner/fMaxDzInner);
	//fFracOfMergedRow += (tDu<fMaxDuInner && tDz<fMaxDzInner);
      }
      else{
	fFracOfMergedRow += (tDu<fMaxDuOuter && tDz<fMaxDzOuter);
	tDist = ::sqrt(tDu*tDu/fMaxDuOuter/fMaxDuOuter+
		     tDz*tDz/fMaxDzOuter/fMaxDzOuter);
	//fFracOfMergedRow += (tDu<fMaxDuOuter && tDz<fMaxDzOuter);
      }
      if(tDist<tDistMax){
	fClosestRowAtDCA = ti+1;
	tDistMax = tDist;
      }
      fWeightedAvSep += tDist;
    }
  }
  if(tN>0){
    fWeightedAvSep /= tN;
    fFracOfMergedRow /= tN;
  }
  else{
    fClosestRowAtDCA = -1;
    fFracOfMergedRow = -1.;
    fWeightedAvSep = -1.;
  }
}
//________________V0 daughters exit/entrance/average separation calc.
//_______1st part is a track 2nd is a V0 considering Pos daughter
double AliFemtoPair::TpcExitSeparationTrackV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcExitPoint() - fTrack2->TpcV0PosExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationTrackV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
  return (diff.mag());
}

double AliFemtoPair::TpcAverageSeparationTrackV0Pos() const {
  AliFemtoThreeVector diff;
  double AveSep = 0.0;
  int ipt = 0;
  if (fTrack1->fNominalPosSample && fTrack2->fNominalPosSample){
  while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
	 fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
	 (ipt<11)
	 ){
    diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
    ipt++;
    AveSep += diff.mag();
  }
  AveSep = AveSep/(ipt+1.);
  return (AveSep);}
  else return -1;
}
//_______1st part is a track 2nd is a V0 considering Neg daughter
double AliFemtoPair::TpcExitSeparationTrackV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcExitPoint() - fTrack2->TpcV0NegExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationTrackV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->NominalTpcEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
  return (diff.mag());
}

double AliFemtoPair::TpcAverageSeparationTrackV0Neg() const {
  AliFemtoThreeVector diff;
  double AveSep = 0.0;
  int ipt = 0;
  if (fTrack1->fNominalPosSample && fTrack2->fTpcV0NegPosSample){
  while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
	 (ipt<11)
	 ){
    diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
    ipt++;
    AveSep += diff.mag();
  }
  AveSep = AveSep/(ipt+1.);
  return (AveSep);}
  else return -1;
}

//_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Pos daughter
double AliFemtoPair::TpcExitSeparationV0PosV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0PosExitPoint() - fTrack2->TpcV0PosExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationV0PosV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0PosEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
  return (diff.mag());
}
double AliFemtoPair::TpcAverageSeparationV0PosV0Pos() const {
  AliFemtoThreeVector diff;
  double AveSep = 0.0;
  int ipt=0;
  if (fTrack1->fNominalPosSample && (fTrack2->fNominalPosSample)){
    while ((fabs(fTrack1->fNominalPosSample[ipt].x())<9999.) &&
	(fabs(fTrack1->fNominalPosSample[ipt].y())<9999.) &&
	(fabs(fTrack1->fNominalPosSample[ipt].z())<9999.) &&
	(fabs(fTrack2->fNominalPosSample[ipt].x())<9999.) &&
	(fabs(fTrack2->fNominalPosSample[ipt].y())<9999.) &&
	(fabs(fTrack2->fNominalPosSample[ipt].z())<9999.) &&
	 (ipt<11)  
	){
      diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
      ipt++;
      AveSep += diff.mag();
    }
    AveSep = AveSep/(ipt+1);
    return (AveSep);}
  else return -1;
}

//_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Neg daughter
double AliFemtoPair::TpcExitSeparationV0PosV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0PosExitPoint() - fTrack2->TpcV0NegExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationV0PosV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0PosEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
  return (diff.mag());
}
double AliFemtoPair::TpcAverageSeparationV0PosV0Neg() const {
  AliFemtoThreeVector diff;
  double AveSep = 0.0;
  int ipt = 0;
  if (fTrack1->fNominalPosSample && fTrack2->fTpcV0NegPosSample){
  while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
	 fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
	 (ipt<11)
	 ){
    diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
    ipt++;
    AveSep += diff.mag();
  }
  AveSep = AveSep/(ipt+1.);
  return (AveSep);}
  else return -1; 
}
//_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Pos daughter
// this is to check the upper case
double AliFemtoPair::TpcExitSeparationV0NegV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0NegExitPoint() - fTrack2->TpcV0PosExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationV0NegV0Pos() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0NegEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
  return (diff.mag());
}
double AliFemtoPair::TpcAverageSeparationV0NegV0Pos() const {
   AliFemtoThreeVector diff;
   double AveSep = 0.0;
   int ipt = 0;
   if ( fTrack1->fTpcV0NegPosSample &&  fTrack2->fNominalPosSample){
     while (fabs(fTrack1->fTpcV0NegPosSample[ipt].x())<9999. &&
	    fabs(fTrack1->fTpcV0NegPosSample[ipt].y())<9999. && 
	    fabs(fTrack1->fTpcV0NegPosSample[ipt].z())<9999. &&
	    fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
	    fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
	    fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
	    (ipt<11)
	    ){
       diff = fTrack1->fTpcV0NegPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
       ipt++;
       AveSep += diff.mag();
     }
     AveSep = AveSep/(ipt+1);
     return (AveSep);}
     else return -1;
}
//_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Neg daughter
double AliFemtoPair::TpcExitSeparationV0NegV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0NegExitPoint() - fTrack2->TpcV0NegExitPoint();
  return (diff.mag());
}

double AliFemtoPair::TpcEntranceSeparationV0NegV0Neg() const {
  AliFemtoThreeVector diff = fTrack1->TpcV0NegEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
  return (diff.mag());
}
double AliFemtoPair::TpcAverageSeparationV0NegV0Neg() const {
   AliFemtoThreeVector diff;
   double AveSep = 0.0;
   int ipt=0;
   if (fTrack1->fTpcV0NegPosSample && fTrack2->fTpcV0NegPosSample){
     while (fabs(fTrack1->fTpcV0NegPosSample[ipt].x())<9999. &&
	    fabs(fTrack1->fTpcV0NegPosSample[ipt].y())<9999. && 
	    fabs(fTrack1->fTpcV0NegPosSample[ipt].z())<9999. &&
	    fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
	    fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
	    fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
	    (ipt<11)
	    ){
       diff = fTrack1->fTpcV0NegPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
       ipt++;
       AveSep += diff.mag();
     }
     AveSep = AveSep/(ipt+1);
     return (AveSep);}
   else return -1;
}

//________________end V0 daughters exit/entrance/average separation calc.
void AliFemtoPair::CalcMergingParFctn(short* tmpMergingParNotCalculatedFctn,
				   float* tmpZ1,float* tmpU1,
				   float* tmpZ2,float* tmpU2,
				   int *tmpSect1,int *tmpSect2,
				   double* tmpFracOfMergedRow,
				   double* tmpClosestRowAtDCA
				   ) const{
  tmpMergingParNotCalculatedFctn=0;
  double tDu, tDz;
  int tN = 0;
  *tmpFracOfMergedRow = 0.;
  *tmpClosestRowAtDCA = 0.;
  double tDist;
  double tDistMax = 100000000.;
  for(int ti=0 ; ti<45 ; ti++){
    if(tmpSect1[ti]==tmpSect2[ti] && tmpSect1[ti]!=-1){
	tDu = fabs(tmpU1[ti]-tmpU2[ti]);
	tDz = fabs(tmpZ1[ti]-tmpZ2[ti]);
	tN++;
      if(ti<13){
	*tmpFracOfMergedRow += (tDu<fMaxDuInner && tDz<fMaxDzInner);
	tDist = ::sqrt(tDu*tDu/fMaxDuInner/fMaxDuInner+
		     tDz*tDz/fMaxDzInner/fMaxDzInner);
      }
      else{
	*tmpFracOfMergedRow += (tDu<fMaxDuOuter && tDz<fMaxDzOuter);
	tDist = ::sqrt(tDu*tDu/fMaxDuOuter/fMaxDuOuter+
		     tDz*tDz/fMaxDzOuter/fMaxDzOuter);
	}
      if(tDist<tDistMax){
	fClosestRowAtDCA = ti+1;
	tDistMax = tDist;
      }
      //fWeightedAvSep += tDist; // now, wrong but not used
    }	
  }
  if(tN>0){
    //fWeightedAvSep /= tN;
    *tmpFracOfMergedRow /= tN;
  }
  else{
    *tmpClosestRowAtDCA = -1;
    *tmpFracOfMergedRow = -1.;
    //fWeightedAvSep = -1.;
  }
}

