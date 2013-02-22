///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoPair: the Pair object is passed to the PairCuts for           //
// verification, and then to the AddRealPair and AddMixedPair methods of //
// the Correlation Functions. It holds pair-specific variables like      //
// relative momenta and has links to the particles and tracks that form  //
// the pair.                                                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include "AliFemtoPair.h"

double AliFemtoPair::fgMaxDuInner = .8;
double AliFemtoPair::fgMaxDzInner = 3.;
double AliFemtoPair::fgMaxDuOuter = 1.4;
double AliFemtoPair::fgMaxDzOuter = 3.2;


AliFemtoPair::AliFemtoPair() :
  fTrack1(0), fTrack2(0),
  fPairAngleEP(0),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  fKStarCalc(0),
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
  // Default constructor
  fTrack1 = 0;
  fTrack2 = 0;
  SetDefaultHalfFieldMergingPar();
}

AliFemtoPair::AliFemtoPair(AliFemtoParticle* a, AliFemtoParticle* b)
  : fTrack1(a), fTrack2(b),
  fPairAngleEP(0),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  fKStarCalc(0),
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
  // Construct a pair from two particles
  SetDefaultHalfFieldMergingPar();
}

void AliFemtoPair::SetDefaultHalfFieldMergingPar(){
  fgMaxDuInner = 3;
  fgMaxDzInner = 4.;
  fgMaxDuOuter = 4.;
  fgMaxDzOuter = 6.;
}
void AliFemtoPair::SetDefaultFullFieldMergingPar(){
  // Set default TPC merging parameters for STAR TPC
  fgMaxDuInner = 0.8;
  fgMaxDzInner = 3.;
  fgMaxDuOuter = 1.4;
  fgMaxDzOuter = 3.2;
}
void AliFemtoPair::SetMergingPar(double aMaxDuInner, double aMaxDzInner,
			      double aMaxDuOuter, double aMaxDzOuter)
{
  // Set TPC merging parameters for STAR TPC
  fgMaxDuInner = aMaxDuInner;
  fgMaxDzInner = aMaxDzInner;
  fgMaxDuOuter = aMaxDuOuter;
  fgMaxDzOuter = aMaxDzOuter;
}

AliFemtoPair::~AliFemtoPair() {
  // Destructor
/* no-op */
}

AliFemtoPair::AliFemtoPair(const AliFemtoPair &aPair):
  fTrack1(0), fTrack2(0),
  fPairAngleEP(0),
  fNonIdParNotCalculated(0),
  fDKSide(0),
  fDKOut(0),
  fDKLong(0),
  fCVK(0),
  fKStarCalc(0),
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
  // Copy constructor
  fTrack1 = aPair.fTrack1;
  fTrack2 = aPair.fTrack2;

  fNonIdParNotCalculated = aPair.fNonIdParNotCalculated;
  fDKSide = aPair.fDKSide;
  fDKOut = aPair.fDKOut;
  fDKLong = aPair.fDKLong;
  fCVK = aPair.fCVK;
  fKStarCalc = aPair.fKStarCalc;

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
  // Assignment operator
  if (this == &aPair)
    return *this;

  fTrack1 = aPair.fTrack1;
  fTrack2 = aPair.fTrack2;

  fNonIdParNotCalculated = aPair.fNonIdParNotCalculated;
  fDKSide = aPair.fDKSide;
  fDKOut = aPair.fDKOut;
  fDKLong = aPair.fDKLong;
  fCVK = aPair.fCVK;
  fKStarCalc = aPair.fKStarCalc;

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

//________________________
double AliFemtoPair::GetPairAngleEP() const
{
	return fPairAngleEP;
}
//_________________
double AliFemtoPair::MInv() const
{
  // invariant mass
    double tInvariantMass = abs(fTrack1->FourMomentum() + fTrack2->FourMomentum());
    return (tInvariantMass);
}
//_________________
double AliFemtoPair::KT() const
{
  // transverse momentum
  double  tmp = 
    (fTrack1->FourMomentum() + fTrack2->FourMomentum()).Perp();
  tmp *= .5;

  return (tmp);
}
//_________________
double AliFemtoPair::Rap() const
{
  // longitudinal pair rapidity : Y = 0.5 ::log( E1 + E2 + pz1 + pz2 / E1 + E2 - pz1 - pz2 )
  double  tmp = 0.5 * log (
			   (fTrack1->FourMomentum().e() + fTrack2->FourMomentum().e() + fTrack1->FourMomentum().z() + fTrack2->FourMomentum().z()) / 
			   (fTrack1->FourMomentum().e() + fTrack2->FourMomentum().e() - fTrack1->FourMomentum().z() - fTrack2->FourMomentum().z()) 
			   ) ;
  return (tmp);
}
//_________________
double AliFemtoPair::EmissionAngle() const {
  // emission angle
  double pxTotal = this->FourMomentumSum().x();
  double pyTotal = this->FourMomentumSum().y();
  double angle = atan2(pyTotal,pxTotal)*180.0/3.1415926536;
  if (angle<0.0) angle+=360.0;
  return angle;
}
//_________________
// get rid of ambiguously-named method fourMomentum() and replace it with
// fourMomentumSum() and fourMomentumDiff() - mal 13feb2000
AliFemtoLorentzVector AliFemtoPair::FourMomentumSum() const
{
  // total momentum
  AliFemtoLorentzVector temp = fTrack1->FourMomentum()+fTrack2->FourMomentum();
  return temp;
}
AliFemtoLorentzVector AliFemtoPair::FourMomentumDiff() const
{
  // momentum difference
  AliFemtoLorentzVector temp = fTrack1->FourMomentum()-fTrack2->FourMomentum();
  return temp;
}
//__________________________________
void AliFemtoPair::QYKPCMS(double& qP, double& qT, double& q0) const
{
  // Yano-Koonin-Podgoretskii Parametrisation in CMS
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
  qT = l.vect().Perp() ;
  q0 = l.e() ;
}
//___________________________________
void AliFemtoPair::QYKPLCMS(double& qP, double& qT, double& q0) const
{
  // Yano-Koonin-Podgoretskii Parametrisation in LCMS
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
  l.SetE(l.vect().MassHypothesis(dummyMass) );
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
  qT = l.vect().Perp() ;
  q0 = l.e() ;
}
//___________________________________
// Yano-Koonin-Podgoretskii Parametrisation in pair rest frame
void AliFemtoPair::QYKPPF(double& qP, double& qT, double& q0) const
{
  ////
  //  calculate momentum difference in pair rest frame : frame where (pz1 + pz2, py1 + py2, px1 + px2) = (0,0,0)
  ////
  AliFemtoLorentzVector l1 = fTrack1->FourMomentum() ;
  AliFemtoLorentzVector l2 = fTrack2->FourMomentum() ;
  // the center of gravity of the pair travels with l
  AliFemtoLorentzVector  l = l1 + l2 ; 
  l = -l ;
  l.SetE(-l.e()) ;
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
  qT = l.vect().Perp();
  q0 = l.e();
}
//_________________
double AliFemtoPair::QOutCMS() const
{
  // relative momentum out component in lab frame
    AliFemtoThreeVector tmp1 = fTrack1->FourMomentum().vect();
    AliFemtoThreeVector tmp2 = fTrack2->FourMomentum().vect();

    double dx = tmp1.x() - tmp2.x();
    double xt = tmp1.x() + tmp2.x();
    
    double dy = tmp1.y() - tmp2.y();
    double yt = tmp1.y() + tmp2.y();

    double k1 = (::sqrt(xt*xt+yt*yt));
    double k2 = (dx*xt+dy*yt);
    double tmp;

    if(k1!=0) tmp= k2/k1;
    else tmp=0;

    return (tmp);
}
//_________________
double AliFemtoPair::QSideCMS() const
{
  // relative momentum side component in lab frame
    AliFemtoThreeVector tmp1 = fTrack1->FourMomentum().vect();
    AliFemtoThreeVector tmp2 = fTrack2->FourMomentum().vect();

    double x1 = tmp1.x();  double y1 = tmp1.y();
    double x2 = tmp2.x();  double y2 = tmp2.y();

    double xt = x1+x2;  double yt = y1+y2;
    double k1 = ::sqrt(xt*xt+yt*yt);
    
    double tmp;
    if(k1!=0) tmp= 2.0*(x2*y1-x1*y2)/k1;
    else tmp=0;

    return (tmp);
}

//_________________________
double AliFemtoPair::QLongCMS() const
{
  // relative momentum component in lab frame
    AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
    AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();

    double dz = tmp1.z() - tmp2.z();
    double zz = tmp1.z() + tmp2.z();

    double dt = tmp1.t() - tmp2.t();
    double tt = tmp1.t() + tmp2.t();

    double beta = zz/tt;
    double gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));

    double temp = gamma*(dz - beta*dt);
    return (temp);
}

//________________________________
double AliFemtoPair::QOutPf() const
{
  // relative momentum out component in pair frame
  AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
  AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();
  
  double dt = tmp1.t() - tmp2.t();
  double tt = tmp1.t() + tmp2.t();
  
  double xt = tmp1.x() + tmp2.x();
  double yt = tmp1.y() + tmp2.y();
  
  double k1 = ::sqrt(xt*xt + yt*yt);
  double bOut = k1/tt;
  double gOut = 1.0/TMath::Sqrt((1.-bOut)*(1.+bOut));
  
  double temp = gOut*(QOutCMS() - bOut*dt);
  return (temp);
}

//___________________________________
double AliFemtoPair::QSidePf() const
{
  // relative momentum side component in pair frame

 return(this->QSideCMS());
}

//___________________________________

double AliFemtoPair::QLongPf() const
{
  // relative momentum long component in pair frame

  return(this->QLongCMS());
}

//___________________________________
double AliFemtoPair::QOutBf(double /* beta */) const
{
  // relative momentum out component
 return(this->QOutCMS());
}

//___________________________________

double AliFemtoPair::QSideBf(double /* beta */) const
{
  // relative momentum side component 
 return(this->QSideCMS());
}

//___________________________________
double AliFemtoPair::QLongBf(double beta) const
{
  // relative momentum long component 
    AliFemtoLorentzVector tmp1 = fTrack1->FourMomentum();
    AliFemtoLorentzVector tmp2 = fTrack2->FourMomentum();

    double dz = tmp1.z() - tmp2.z();
    double dt = tmp1.t() + tmp2.t();

    double gamma = 1.0/::sqrt((1.-beta)*(1.+beta));

    double temp = gamma*(dz - beta*dt);
    return (temp);
}

double AliFemtoPair::Quality() const {
  // Calculate split quality of the pair
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
  int tQuality = 0;
  double normQual = 0.0;
  int tMaxQuality = fTrack1->NumberOfHits() + fTrack2->NumberOfHits();
  for (ibits=8;ibits<=31;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad1To24 & bitI ) {
      tQuality++;
      continue;
    }
    else{
      if ( bothPads1To24 & bitI ) tQuality--;
    }
  }
  for (ibits=0;ibits<=20;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad25To45 & bitI ) {
      tQuality++;
      continue;
    }
    else{
      if ( bothPads25To45 & bitI ) tQuality--;
    }
  }
  normQual = (double)tQuality/( (double) tMaxQuality );
  return ( normQual );

}

double AliFemtoPair::Quality2() const {
  // second implementation of split quality
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
  int tQuality = 0;
  double normQual = 0.0;
  int tMaxQuality = fTrack1->NumberOfHits() + fTrack2->NumberOfHits();
  for (ibits=8;ibits<=31;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad1To24 & bitI ) {
      tQuality++;
      continue;
    }
    //else{
    //if ( bothPads1To24 & bitI ) tQuality--;
    //}
  }
  for (ibits=0;ibits<=20;ibits++) {
    bitI = 0;
    bitI |= 1UL<<(ibits);
    if ( onePad25To45 & bitI ) {
      tQuality++;
      continue;
    }
    //else{
    //if ( bothPads25To45 & bitI ) tQuality--;
    //}
  }
  normQual = (double)tQuality/( (double) tMaxQuality );
  return ( normQual );

}


double AliFemtoPair::NominalTpcExitSeparation() const {
  // separation at exit from STAR TPC
  AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcExitPoint() - fTrack2->Track()->NominalTpcExitPoint();
  return (diff.Mag());
}

double AliFemtoPair::NominalTpcEntranceSeparation() const {
  // separation at entrance to STAR TPC
  AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcEntrancePoint() - fTrack2->Track()->NominalTpcEntrancePoint();
  return (diff.Mag());
}

// double AliFemtoPair::NominalTpcAverageSeparation() const {
//   // average separation in STAR TPC
//   AliFemtoThreeVector diff;
//   double tAveSep = 0.0;
//   int ipt = 0;
//   if (fTrack1->fNominalPosSample && fTrack2->fNominalPosSample){
//   while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
// 	 fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
// 	 ipt<11
// 	 ){
//     //  for (int ipt=0; ipt<11; ipt++){
//     diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
//     ipt++;
//     tAveSep += diff.Mag();
//   }
//   tAveSep = tAveSep/(ipt+1.);
//   return (tAveSep);}
//   else return -1;
// }

double AliFemtoPair::OpeningAngle() const {
  // opening angle
 return 57.296* fTrack1->FourMomentum().vect().Angle( fTrack2->FourMomentum().vect() );
//   AliFemtoThreeVector p1 = fTrack1->FourMomentum().vect();
//   AliFemtoThreeVector p2 = fTrack2->FourMomentum().vect();
//   return 57.296*(p1.phi()-p2.phi());
//   //double dAngInv = 57.296*acos((p1.dot(p2))/(p1.Mag()*p2.Mag()));
//   //return (dAngInv);
}
//_________________


double AliFemtoPair::KStarFlipped() const {
  // kstar with sign flipped
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();

  AliFmThreeVectorD qwe = tP1.vect();
  qwe *= -1.; // flip it
  tP1.SetVect(qwe);
  
  AliFemtoLorentzVector tSum = (tP1+fTrack2->FourMomentum());
  double tMass = abs(tSum);
  AliFmThreeVectorD tGammaBeta = (1./tMass)*tSum.vect(); 
  double tGamma = tSum.e()/tMass;
  AliFmThreeVectorD tLongMom  = ((tP1.vect()*tGammaBeta)/
			      (tGammaBeta*tGammaBeta))*tGammaBeta;
  AliFmLorentzVectorD tK(tGamma*tP1.e() - tP1.vect()*tGammaBeta,
		      tP1.vect() + (tGamma-1.)*tLongMom - tP1.e()*tGammaBeta);
//VP  tP1.vect() *= -1.; // unflip it
  return tK.vect().Mag();
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
//return (tK.vect())*tGammaBeta/tK.vect().Magnitude()/tGammaBeta.Magnitude();
//}

double AliFemtoPair::CVKFlipped() const{
  // CVK with sign flipped
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  AliFmThreeVectorD qwe = tP1.vect();
  qwe *= -1.; // flip it
  tP1.SetVect(qwe);
  
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

double AliFemtoPair::PInv() const{
  // invariant total momentum
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  AliFemtoLorentzVector tP2 = fTrack2->FourMomentum();
  double tP = (tP1.px()+tP2.px())*(tP1.px()+tP2.px())+
              (tP1.py()+tP2.py())*(tP1.py()+tP2.py())+
              (tP1.pz()+tP2.pz())*(tP1.pz()+tP2.pz())-
              (tP1.e() -tP2.e() )*(tP1.e() -tP2.e() );
  return ::sqrt(fabs(tP));
}

double AliFemtoPair::QInvFlippedXY() const{
  // qinv with X and Y flipped
  AliFemtoLorentzVector tP1 = fTrack1->FourMomentum();
  tP1.SetX(-1.*tP1.x());
  tP1.SetY(-1.*tP1.y());
  AliFemtoLorentzVector tDiff = (tP1-fTrack2->FourMomentum());
  return ( -1.* tDiff.m());
}

void AliFemtoPair::CalcNonIdPar() const{ // fortran like function! faster?
  // Calculate generalized relative mometum 
  // Use this instead of qXYZ() function when calculating
  // anything for non-identical particles
  fNonIdParNotCalculated=0;
  double px1 = fTrack1->FourMomentum().vect().x();
  double py1 = fTrack1->FourMomentum().vect().y();
  double pz1 = fTrack1->FourMomentum().vect().z();
  double pE1  = fTrack1->FourMomentum().e();
  double tParticle1Mass = ::sqrt(pE1*pE1 - px1*px1 - py1*py1 - pz1*pz1);

  double px2 = fTrack2->FourMomentum().vect().x();
  double py2 = fTrack2->FourMomentum().vect().y();
  double pz2 = fTrack2->FourMomentum().vect().z();
  double pE2  = fTrack2->FourMomentum().e();
  double tParticle2Mass = ::sqrt(pE2*pE2 - px2*px2 - py2*py2 - pz2*pz2);

  double tPx = px1+px2;
  double tPy = py1+py2;
  double tPz = pz1+pz2;
  double tPE = pE1+pE2;
      
  double tPtrans = tPx*tPx + tPy*tPy;
  double tMtrans = tPE*tPE - tPz*tPz;
  double tPinv =   ::sqrt(tMtrans - tPtrans);
  tMtrans = ::sqrt(tMtrans);
  tPtrans = ::sqrt(tPtrans);
	
  double tQinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
    (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);

  double tQ = (tParticle1Mass*tParticle1Mass - tParticle2Mass*tParticle2Mass)/tPinv;
  tQ = sqrt ( tQ*tQ - tQinvL);
	  
  fKStarCalc = tQ/2;

  // ad 1) go to LCMS
  double beta = tPz/tPE;
  double gamma = tPE/tMtrans;
	    
  double pz1L = gamma * (pz1 - beta * pE1);
  double pE1L = gamma * (pE1 - beta * pz1);
  
  // fill histogram for beam projection ( z - axis )
  fDKLong = pz1L;

  // ad 2) rotation px -> tPt
  double px1R = (px1*tPx + py1*tPy)/tPtrans;
  double py1R = (-px1*tPy + py1*tPx)/tPtrans;
  
  //fill histograms for side projection ( y - axis )
  fDKSide = py1R;

  // ad 3) go from LCMS to CMS
  beta = tPtrans/tMtrans;
  gamma = tMtrans/tPinv;
  
  double px1C = gamma * (px1R - beta * pE1L);
  
  // fill histogram for out projection ( x - axis )
  fDKOut  = px1C;

  fCVK = (fDKOut*tPtrans + fDKLong*tPz)/fKStarCalc/::sqrt(tPtrans*tPtrans+tPz*tPz);
}


/*void AliFemtoPair::calcNonIdParGlobal() const{ // fortran like function! faster?
  fNonIdParNotCalculatedGlobal=0;
  double px1 = fTrack1->Track()->PGlobal().x();
  double py1 = fTrack1->Track()->PGlobal().y();
  double pz1 = fTrack1->Track()->PGlobal().z();
  double tParticle1Mass =  fTrack1->FourMomentum().m2();
  double pE1  = ::sqrt(tParticle1Mass + px1*px1 + py1*py1 + pz1*pz1);
  tParticle1Mass = ::sqrt(tParticle1Mass);

  double px2 = fTrack2->Track()->PGlobal().x();
  double py2 = fTrack2->Track()->PGlobal().y();
  double pz2 = fTrack2->Track()->PGlobal().z();
  double tParticle2Mass =  fTrack2->FourMomentum().m2();
  double pE2  = ::sqrt(tParticle2Mass + px2*px2 + py2*py2 + pz2*pz2);
  tParticle2Mass = ::sqrt(tParticle2Mass);

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

  double Q = (tParticle1Mass*tParticle1Mass - tParticle2Mass*tParticle2Mass)/Pinv;
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



// double AliFemtoPair::DcaInsideTpc() const{
//   // dcs inside the STAR TPC
//   double tMinDist=NominalTpcEntranceSeparation();
//   double tExit = NominalTpcExitSeparation();
//   tMinDist = (tExit>tMinDist) ? tMinDist : tExit;
//   double tInsideDist;
//   //tMinDist = 999.;

//   double rMin = 60.;
//   double rMax = 190.;
//   const AliFmPhysicalHelixD& tHelix1 = fTrack1->Helix();
//   const AliFmPhysicalHelixD& tHelix2 = fTrack2->Helix();
//   // --- One is a line and other one a helix
//   //if (tHelix1.mSingularity != tHelix2.mSingularity) return -999.;
//   // --- 2 lines : don't care right now
//   //if (tHelix1.mSingularity)  return -999.;
//   // --- 2 helix
//   double dx = tHelix2.XCenter() - tHelix1.XCenter();
//   double dy = tHelix2.YCenter() - tHelix1.YCenter();
//   double dd = ::sqrt(dx*dx + dy*dy);
//   double r1 = 1/tHelix1.Curvature();
//   double r2 = 1/tHelix2.Curvature();
//   double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);
    
//   double x, y, r;
//   double s;
//   if (fabs(cosAlpha) < 1) {           // two solutions
//     double sinAlpha = sin(acos(cosAlpha));
//     x = tHelix1.XCenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
//     y = tHelix1.YCenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
//     r = ::sqrt(x*x+y*y);
//     if( r > rMin &&  r < rMax && 
// 	fabs(atan2(y,x)-fTrack1->Track()->NominalTpcEntrancePoint().phi())< 0.5
// 	){ // first solution inside
//       s = tHelix1.PathLength(x, y);
//       tInsideDist=tHelix2.Distance(tHelix1.At(s));
//       if(tInsideDist<tMinDist) tMinDist = tInsideDist;
//     }
//     else{ 
//       x = tHelix1.XCenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
//       y = tHelix1.YCenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
//       r = ::sqrt(x*x+y*y);
//       if( r > rMin &&  r < rMax &&
// 	  fabs(atan2(y,x)-fTrack1->Track()->NominalTpcEntrancePoint().phi())< 0.5
// 	  ) {  // second solution inside
//         s = tHelix1.PathLength(x, y);
//         tInsideDist=tHelix2.Distance(tHelix1.At(s));
//         if(tInsideDist<tMinDist) tMinDist = tInsideDist;
//       }     
//     }
//   }
//   return tMinDist;
// }

// void AliFemtoPair::CalcMergingPar() const{
//   // Calculate merging factor for the pair in STAR TPC
//   fMergingParNotCalculated=0;

//   double tDu, tDz;
//   int tN = 0;
//   fFracOfMergedRow = 0.;
//   fWeightedAvSep =0.;
//   double tDist;
//   double tDistMax = 200.;
//   for(int ti=0 ; ti<45 ; ti++){
//     if(fTrack1->fSect[ti]==fTrack2->fSect[ti] && fTrack1->fSect[ti]!=-1){
//       tDu = fabs(fTrack1->fU[ti]-fTrack2->fU[ti]);
//       tDz = fabs(fTrack1->fZ[ti]-fTrack2->fZ[ti]);
//       tN++;
//       if(ti<13){
// 	fFracOfMergedRow += (tDu<fgMaxDuInner && tDz<fgMaxDzInner);
// 	tDist = ::sqrt(tDu*tDu/fgMaxDuInner/fgMaxDuInner+
// 		     tDz*tDz/fgMaxDzInner/fgMaxDzInner);
// 	//fFracOfMergedRow += (tDu<fgMaxDuInner && tDz<fgMaxDzInner);
//       }
//       else{
// 	fFracOfMergedRow += (tDu<fgMaxDuOuter && tDz<fgMaxDzOuter);
// 	tDist = ::sqrt(tDu*tDu/fgMaxDuOuter/fgMaxDuOuter+
// 		     tDz*tDz/fgMaxDzOuter/fgMaxDzOuter);
// 	//fFracOfMergedRow += (tDu<fgMaxDuOuter && tDz<fgMaxDzOuter);
//       }
//       if(tDist<tDistMax){
// 	fClosestRowAtDCA = ti+1;
// 	tDistMax = tDist;
//       }
//       fWeightedAvSep += tDist;
//     }
//   }
//   if(tN>0){
//     fWeightedAvSep /= tN;
//     fFracOfMergedRow /= tN;
//   }
//   else{
//     fClosestRowAtDCA = -1;
//     fFracOfMergedRow = -1.;
//     fWeightedAvSep = -1.;
//   }
// }
// double AliFemtoPair::TpcExitSeparationTrackV0Pos() const {
// //________________V0 daughters exit/entrance/average separation calc.
// //_______1st part is a track 2nd is a V0 considering Pos daughter
  
//   AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcExitPoint() - fTrack2->TpcV0PosExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationTrackV0Pos() const {
// //________________V0 daughters exit/entrance/average separation calc.
// //_______1st part is a track 2nd is a V0 considering Pos daughter
//   AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcAverageSeparationTrackV0Pos() const {
// //________________V0 daughters exit/entrance/average separation calc.
// //_______1st part is a track 2nd is a V0 considering Pos daughter
//   AliFemtoThreeVector diff;
//   double tAveSep = 0.0;
//   int ipt = 0;
//   if (fTrack1->fNominalPosSample && fTrack2->fNominalPosSample){
//   while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
// 	 fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
// 	 (ipt<11)
// 	 ){
//     diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
//     ipt++;
//     tAveSep += diff.Mag();
//   }
//   tAveSep = tAveSep/(ipt+1.);
//   return (tAveSep);}
//   else return -1;
// }
// double AliFemtoPair::TpcExitSeparationTrackV0Neg() const {
// //_______1st part is a track 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcExitPoint() - fTrack2->TpcV0NegExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationTrackV0Neg() const {
// //_______1st part is a track 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->Track()->NominalTpcEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcAverageSeparationTrackV0Neg() const {
// //_______1st part is a track 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff;
//   double tAveSep = 0.0;
//   int ipt = 0;
//   if (fTrack1->fNominalPosSample && fTrack2->fTpcV0NegPosSample){
//   while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
// 	 (ipt<11)
// 	 ){
//     diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
//     ipt++;
//     tAveSep += diff.Mag();
//   }
//   tAveSep = tAveSep/(ipt+1.);
//   return (tAveSep);}
//   else return -1;
// }

// double AliFemtoPair::TpcExitSeparationV0PosV0Pos() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Pos daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0PosExitPoint() - fTrack2->TpcV0PosExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationV0PosV0Pos() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Pos daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0PosEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
//   return (diff.Mag());
// }
// double AliFemtoPair::TpcAverageSeparationV0PosV0Pos() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Pos daughter
//   AliFemtoThreeVector diff;
//   double tAveSep = 0.0;
//   int ipt=0;
//   if (fTrack1->fNominalPosSample && (fTrack2->fNominalPosSample)){
//     while ((fabs(fTrack1->fNominalPosSample[ipt].x())<9999.) &&
// 	(fabs(fTrack1->fNominalPosSample[ipt].y())<9999.) &&
// 	(fabs(fTrack1->fNominalPosSample[ipt].z())<9999.) &&
// 	(fabs(fTrack2->fNominalPosSample[ipt].x())<9999.) &&
// 	(fabs(fTrack2->fNominalPosSample[ipt].y())<9999.) &&
// 	(fabs(fTrack2->fNominalPosSample[ipt].z())<9999.) &&
// 	 (ipt<11)  
// 	){
//       diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
//       ipt++;
//       tAveSep += diff.Mag();
//     }
//     tAveSep = tAveSep/(ipt+1);
//     return (tAveSep);}
//   else return -1;
// }

// double AliFemtoPair::TpcExitSeparationV0PosV0Neg() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0PosExitPoint() - fTrack2->TpcV0NegExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationV0PosV0Neg() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0PosEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
//   return (diff.Mag());
// }
// double AliFemtoPair::TpcAverageSeparationV0PosV0Neg() const {
// //_______1st part is a V0 considering Pos daughter 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff;
//   double tAveSep = 0.0;
//   int ipt = 0;
//   if (fTrack1->fNominalPosSample && fTrack2->fTpcV0NegPosSample){
//   while (fabs(fTrack1->fNominalPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack1->fNominalPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack1->fNominalPosSample[ipt].z())<9999. &&
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
// 	 fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
// 	 (ipt<11)
// 	 ){
//     diff = fTrack1->fNominalPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
//     ipt++;
//     tAveSep += diff.Mag();
//   }
//   tAveSep = tAveSep/(ipt+1.);
//   return (tAveSep);}
//   else return -1; 
// }
// double AliFemtoPair::TpcExitSeparationV0NegV0Pos() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Pos daughter
// // this is to check the upper case
//   AliFemtoThreeVector diff = fTrack1->TpcV0NegExitPoint() - fTrack2->TpcV0PosExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationV0NegV0Pos() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Pos daughter
// // this is to check the upper case
//   AliFemtoThreeVector diff = fTrack1->TpcV0NegEntrancePoint() - fTrack2->TpcV0PosEntrancePoint();
//   return (diff.Mag());
// }
// double AliFemtoPair::TpcAverageSeparationV0NegV0Pos() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Pos daughter
// // this is to check the upper case
//    AliFemtoThreeVector diff;
//    double tAveSep = 0.0;
//    int ipt = 0;
//    if ( fTrack1->fTpcV0NegPosSample &&  fTrack2->fNominalPosSample){
//      while (fabs(fTrack1->fTpcV0NegPosSample[ipt].x())<9999. &&
// 	    fabs(fTrack1->fTpcV0NegPosSample[ipt].y())<9999. && 
// 	    fabs(fTrack1->fTpcV0NegPosSample[ipt].z())<9999. &&
// 	    fabs(fTrack2->fNominalPosSample[ipt].x())<9999. &&
// 	    fabs(fTrack2->fNominalPosSample[ipt].y())<9999. && 
// 	    fabs(fTrack2->fNominalPosSample[ipt].z())<9999. &&
// 	    (ipt<11)
// 	    ){
//        diff = fTrack1->fTpcV0NegPosSample[ipt] - fTrack2->fNominalPosSample[ipt];
//        ipt++;
//        tAveSep += diff.Mag();
//      }
//      tAveSep = tAveSep/(ipt+1);
//      return (tAveSep);}
//      else return -1;
// }
// double AliFemtoPair::TpcExitSeparationV0NegV0Neg() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0NegExitPoint() - fTrack2->TpcV0NegExitPoint();
//   return (diff.Mag());
// }

// double AliFemtoPair::TpcEntranceSeparationV0NegV0Neg() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Neg daughter
//   AliFemtoThreeVector diff = fTrack1->TpcV0NegEntrancePoint() - fTrack2->TpcV0NegEntrancePoint();
//   return (diff.Mag());
// }
// double AliFemtoPair::TpcAverageSeparationV0NegV0Neg() const {
// //_______1st part is a V0 considering Neg daughter 2nd is a V0 considering Neg daughter
//    AliFemtoThreeVector diff;
//    double tAveSep = 0.0;
//    int ipt=0;
//    if (fTrack1->fTpcV0NegPosSample && fTrack2->fTpcV0NegPosSample){
//      while (fabs(fTrack1->fTpcV0NegPosSample[ipt].x())<9999. &&
// 	    fabs(fTrack1->fTpcV0NegPosSample[ipt].y())<9999. && 
// 	    fabs(fTrack1->fTpcV0NegPosSample[ipt].z())<9999. &&
// 	    fabs(fTrack2->fTpcV0NegPosSample[ipt].x())<9999. &&
// 	    fabs(fTrack2->fTpcV0NegPosSample[ipt].y())<9999. && 
// 	    fabs(fTrack2->fTpcV0NegPosSample[ipt].z())<9999. &&
// 	    (ipt<11)
// 	    ){
//        diff = fTrack1->fTpcV0NegPosSample[ipt] - fTrack2->fTpcV0NegPosSample[ipt];
//        ipt++;
//        tAveSep += diff.Mag();
//      }
//      tAveSep = tAveSep/(ipt+1);
//      return (tAveSep);}
//    else return -1;
// }

// void AliFemtoPair::CalcMergingParFctn(short* tmpMergingParNotCalculatedFctn,
// 				   float* tmpZ1,float* tmpU1,
// 				   float* tmpZ2,float* tmpU2,
// 				   int *tmpSect1,int *tmpSect2,
// 				   double* tmpFracOfMergedRow,
// 				   double* tmpClosestRowAtDCA
// 				   ) const{
// // calculate heper variables for merging 
//   tmpMergingParNotCalculatedFctn=0;
//   double tDu, tDz;
//   int tN = 0;
//   *tmpFracOfMergedRow = 0.;
//   *tmpClosestRowAtDCA = 0.;
//   double tDist;
//   double tDistMax = 100000000.;
//   for(int ti=0 ; ti<45 ; ti++){
//     if(tmpSect1[ti]==tmpSect2[ti] && tmpSect1[ti]!=-1){
// 	tDu = fabs(tmpU1[ti]-tmpU2[ti]);
// 	tDz = fabs(tmpZ1[ti]-tmpZ2[ti]);
// 	tN++;
//       if(ti<13){
// 	*tmpFracOfMergedRow += (tDu<fgMaxDuInner && tDz<fgMaxDzInner);
// 	tDist = ::sqrt(tDu*tDu/fgMaxDuInner/fgMaxDuInner+
// 		     tDz*tDz/fgMaxDzInner/fgMaxDzInner);
//       }
//       else{
// 	*tmpFracOfMergedRow += (tDu<fgMaxDuOuter && tDz<fgMaxDzOuter);
// 	tDist = ::sqrt(tDu*tDu/fgMaxDuOuter/fgMaxDuOuter+
// 		     tDz*tDz/fgMaxDzOuter/fgMaxDzOuter);
// 	}
//       if(tDist<tDistMax){
// 	fClosestRowAtDCA = ti+1;
// 	tDistMax = tDist;
//       }
//       //fWeightedAvSep += tDist; // now, wrong but not used
//     }	
//   }
//   if(tN>0){
//     //fWeightedAvSep /= tN;
//     *tmpFracOfMergedRow /= tN;
//   }
//   else{
//     *tmpClosestRowAtDCA = -1;
//     *tmpFracOfMergedRow = -1.;
//     //fWeightedAvSep = -1.;
//   }
// }

