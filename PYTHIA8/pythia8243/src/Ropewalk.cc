// Ropewalk.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// RopeRandState, RopeDipoleEnd, OverlappingRopeDipole, RopeDipole, Ropewalk,
// RopeFragPars and FlavourRope classes.

#include "Pythia8/Ropewalk.h"

namespace Pythia8 {

//==========================================================================

// OverlappingRopeDipole class.
// This class describes dipoles overlapping with a given dipole.

//--------------------------------------------------------------------------

// Constructor sets up coordinates in other dipole's rest frame.

OverlappingRopeDipole::OverlappingRopeDipole(RopeDipole* d, double m0,
  RotBstMatrix& r) : dipole(d), dir(1) {

  // Coordinates in other dipole's rest frame
  b1 = d->d1Ptr()->getParticlePtr()->vProd();
  b1.rotbst(r);
  b2 = d->d2Ptr()->getParticlePtr()->vProd();
  b2.rotbst(r);
  y1 = d->d1Ptr()->rap(m0,r);
  y2 = d->d2Ptr()->rap(m0,r);
  if (y1 < y2) dir = -1;

}

//--------------------------------------------------------------------------

// Calculate the overlap at given y and b.

bool OverlappingRopeDipole::overlap(double y, Vec4 ba, double r0) {

  if (y < min(y1, y2) || y > max(y1, y2)) return false;
  Vec4 bb  = b1 + (b2 - b1) * (y - y1) / (y2 - y1);
  Vec4 tmp = ba - bb;
  return (tmp.pT() <= 2 * r0);

}

//--------------------------------------------------------------------------

// Has the dipole been hadronized?

bool OverlappingRopeDipole::hadronized() {

  return dipole->hadronized();

}

//==========================================================================

// RopeDipole class.
// This class describes a dipoles in impact parameter space.

//--------------------------------------------------------------------------

// The RopeDipole constructor makes sure that d1 is always the colored
// end and d2 the anti-colored.

RopeDipole::RopeDipole(RopeDipoleEnd d1In, RopeDipoleEnd d2In, int iSubIn,
  Info* infoPtrIn)
  : d1(d1In), d2(d2In), iSub(iSubIn), hasRotFrom(false), hasRotTo(false),
  isHadronized(false), infoPtr(infoPtrIn) {

  // Test if d1 is colored end and d2 anti-colored.
  if (d1In.getParticlePtr()->col() == d2In.getParticlePtr()->acol()
    && d1In.getParticlePtr()->col() != 0) { }
  else { d2 = d1In, d1 = d2In; }

}

//--------------------------------------------------------------------------

// Insert an excitation on dipole, if not already there.

void RopeDipole::addExcitation(double ylab, Particle* ex) {

  pair<map<double, Particle*>::iterator, map<double, Particle*>::iterator>
    ret = excitations.equal_range(ylab);
  for (map<double, Particle*>::iterator itr = ret.first; itr != ret.second;
    ++itr)
    if (ex == itr->second) return;
  excitations.insert( make_pair(ylab,ex) );

}

//--------------------------------------------------------------------------

// Propagate the dipole itself.

void RopeDipole::propagateInit(double deltat) {

  // Dipole end momenta.
  Vec4 pcm = d1.getParticlePtr()->p();
  Vec4 pam = d2.getParticlePtr()->p();
  double mTc = sqrt(pcm.pT2() + pcm.m2Calc());
  double mTa = sqrt(pam.pT2() + pam.m2Calc());
  if (mTc == 0 || mTa == 0)
    infoPtr->errorMsg("Error in RopeDipole::propagateInit: Tried to"
      "propagate a RopeDipoleEnd with mT = 0");

  // New vertices in the lab frame.
  Vec4 newv1 = Vec4(d1.getParticlePtr()->xProd() + deltat * pcm.px() / mTc,
                d1.getParticlePtr()->yProd() + deltat * pcm.py() / mTc, 0, 0);
  Vec4 newv2 = Vec4(d2.getParticlePtr()->xProd() + deltat * pam.px() / mTa,
                d2.getParticlePtr()->yProd() + deltat * pam.py() / mTa, 0, 0);
  // Set the new vertices deep.
  d1.getParticlePtr()->vProd(newv1);
  d2.getParticlePtr()->vProd(newv2);

}

//--------------------------------------------------------------------------

// Propagate both dipole ends as well as all excitations.

void RopeDipole::propagate(double deltat, double m0) {

  // First propagate the dipole ends.
  propagateInit(deltat);
  for (map<double, Particle*>::iterator eItr = excitations.begin();
    eItr != excitations.end(); ++eItr) {

    Vec4 em = eItr->second->p();
    em.rotbst(getDipoleLabFrame());
    // Propagate excitations.

    if (em.pT() > 0.0){
      Vec4 newVert = Vec4(eItr->second->xProd() + deltat * em.px() / em.pT(),
                eItr->second->yProd() + deltat * em.py() / em.pT(), 0, 0);
      eItr->second->vProd(newVert);
    }
    else eItr->second->vProd(bInterpolateLab(eItr->first,m0));
  }

}

//--------------------------------------------------------------------------

// Put gluon excitations on the dipole.

void RopeDipole::excitationsToString(double m0, Event& event) {

  // Erase excitations below cut-off.
  map<double, Particle*>::iterator pItr = excitations.begin();
  while (pItr != excitations.end() ) {
    if (pItr->second->pAbs() < 1e-6) {
      map<double, Particle*>::iterator eraseMe = pItr;
      ++pItr;
      excitations.erase(eraseMe);
    }
    else ++pItr;
  }

  // We now colour connect the excitations to the dipole.
  // The dipole is connected like this sketch:
  // acol  (d1)  col ------ acol  (d2)  col.
  int oldcol = d1.getParticlePtr()->col();
  if (oldcol != d2.getParticlePtr()->acol()) {
    infoPtr->errorMsg("Error in Ropewalk::RopeDipole::excitationsToString: "
      "color indices do not match.");
    return;
  }
  vector<int> daughters;

  // Two cases depending on which end we should start at.
  // We always start at min rapidity and connect from there.
  if (d1.rap(m0) == minRapidity(m0)) {
    int acol = oldcol;
    for (map<double, Particle*>::iterator itr = excitations.begin();
      itr != excitations.end(); ++itr) {
      int col = event.nextColTag();
      itr->second->status(51);
      itr->second->mothers(d1.getNe(),d1.getNe());
      itr->second->cols(col,acol);
      daughters.push_back(event.append(Particle(*(itr->second))));
      acol = col;
    }
    d2.getParticlePtr()->acol(acol);
    event[d2.getNe()].acol(acol);
  }
  else {
    int acol = oldcol;
    for (map<double, Particle*>::reverse_iterator itr = excitations.rbegin();
      itr != excitations.rend(); ++itr) {
      int col = event.nextColTag();
      itr->second->status(51);
      itr->second->mothers(d1.getNe(),d1.getNe());
      itr->second->cols(col,acol);
      daughters.push_back(event.append(Particle(*(itr->second))));
      acol = col;
    }
    d2.getParticlePtr()->acol(acol);
    event[d2.getNe()].acol(acol);
  }
  bool stringEnd = false;
  if (d2.getParticlePtr()->col() == 0) stringEnd = true;

  // Update status codes and mother/daughter indexing.
  event[d1.getNe()].statusNeg();
  Particle cc1 = *d1.getParticlePtr();
  cc1.statusPos();
  cc1.mothers(d1.getNe(),d1.getNe());
  daughters.push_back(event.append(cc1));
  event[d1.getNe()].daughters( daughters[0], daughters[daughters.size() -1 ] );
  if (stringEnd) {
    event[d2.getNe()].statusNeg();
    Particle cc2 = *d2.getParticlePtr();
    cc2.statusPos();
    cc2.mothers(d2.getNe(),d2.getNe());
    int did = event.append(cc2);
    event[d2.getNe()].daughters(did,did);
  }

}

//--------------------------------------------------------------------------

// Redistribute momentum to two particles.

void RopeDipole::splitMomentum(Vec4 mom, Particle* p1, Particle* p2,
  double frac) {

  Vec4 p1new = p1->p() + frac * mom;
  Vec4 p2new = p2->p() + (1. - frac) * mom;
  p1->p(p1new);
  p2->p(p2new);
  return;

}

//--------------------------------------------------------------------------

// Recoil the dipole from adding a gluon. If the "dummy" option is set,
// the recoil will not be added, but only checked.
// Note: the gluon will not actually be added, only the recoil (if possible).

bool RopeDipole::recoil(Vec4& pg, bool dummy) {

  // Keep track of direction.
  int sign = 1;
  if (d1.rap(1.0) > d2.rap(1.0)) sign = -1;

  // Lightcone momenta after inserting the gluon.
  Particle* epaPtr = d1.getParticlePtr();
  Particle* epcPtr = d2.getParticlePtr();
  double pplus = epcPtr->pPos() + epaPtr->pPos() - pg.pPos();
  double pminus = epcPtr->pNeg() + epaPtr->pNeg() - pg.pNeg();

  // The new lightcone momenta of the dipole ends.
  double ppa = 0.0;
  double ppc = 0.0;
  double pma = 0.0;
  double pmc = 0.0;
  double mta2 = epaPtr->mT2();
  double mtc2 = epcPtr->mT2();
  double mta = sqrt(mta2);
  double mtc = sqrt(mtc2);
  if ( pplus * pminus <= pow2(mta + mtc)
    || pplus <= 0.0 || pminus <= 0.0 ) return false;

  // Calculate the new momenta.
  double sqarg = pow2(pplus * pminus - mta2 - mtc2) - 4. * mta2 * mtc2;
  if (sqarg <= 0.0 ) return false;
  if (sign > 0) {
    ppa = 0.5 * (pplus * pminus + mta2 - mtc2 + sqrt(sqarg)) / pminus;
    pma = mta2 / ppa;
    pmc = pminus - pma;
    ppc = mtc2 / pmc;
    // Check rapidity after recoil.
    if ( ppa * mtc < ppc * mta ) return false;
  }
  else {
    pma = 0.5 * (pplus * pminus + mta2 - mtc2 + sqrt(sqarg)) / pplus;
    ppa = mta2 / pma;
    ppc = pplus - ppa;
    pmc = mtc2 / ppc;
    // Check rapidity after recoil.
    if ( ppa*mtc > ppc*mta ) return false;
  }

  // Set up and store the new momenta.
  Vec4 shifta = Vec4( epaPtr->px(), epaPtr->py(),
    0.5 * (ppa - pma), 0.5 * (ppa + pma));
  Vec4 shiftc = Vec4( epcPtr->px(), epcPtr->py(),
    0.5 * (ppc - pmc), 0.5 * (ppc + pmc));
  if (!dummy) {
    epaPtr->p(shifta);
    epcPtr->p(shiftc);
  }
  return true;

}

//--------------------------------------------------------------------------

// Get the Lorentz matrix to go to the dipole rest frame.

RotBstMatrix RopeDipole::getDipoleRestFrame() {

  if (hasRotTo) return rotTo;

  RotBstMatrix r;
  r.toCMframe(d1.getParticlePtr()->p(),d2.getParticlePtr()->p());
  rotTo = r;
  hasRotTo = true;
  return rotTo;

}

//--------------------------------------------------------------------------

// Get the Lorentz matrix to go from the dipole rest frame to lab frame.

RotBstMatrix RopeDipole::getDipoleLabFrame() {

  if(hasRotFrom) return rotFrom;

  RotBstMatrix r;
  r.fromCMframe(d1.getParticlePtr()->p(),d2.getParticlePtr()->p());
  rotFrom = r;
  hasRotFrom = true;
  return rotFrom;

}
//--------------------------------------------------------------------------

// The dipole four-momentum.

Vec4 RopeDipole::dipoleMomentum() {

  Vec4 ret = d1.getParticlePtr()->p() + d2.getParticlePtr()->p();
  return ret;

}

//--------------------------------------------------------------------------

// Interpolate (linear) between dipole ends to get b position at given y.
// Here y must be given in dipole rest frame, and the resulting b-position
// will also be in the dipole rest frame.

Vec4 RopeDipole::bInterpolateDip(double y, double m0) {
  if(!hasRotTo) getDipoleRestFrame();
  Vec4 bb1 = d1.getParticlePtr()->vProd();
  bb1.rotbst(rotTo);
  Vec4 bb2 = d2.getParticlePtr()->vProd();
  bb2.rotbst(rotTo);
  double y1 = d1.rap(m0,rotTo);
  double y2 = d2.rap(m0,rotTo);
  return bb1 + y * (bb2 - bb1) / (y2 - y1);

}

//--------------------------------------------------------------------------

// Interpolate (linear) between dipole ends to get b position at given y.

Vec4 RopeDipole::bInterpolateLab(double y, double m0) {

  Vec4 bb1 = d1.getParticlePtr()->vProd();
  Vec4 bb2 = d2.getParticlePtr()->vProd();
  double y1 = d1.rap(m0);
  double y2 = d2.rap(m0);
  return bb1 + y * (bb2 - bb1) / (y2 - y1);

}

//--------------------------------------------------------------------------

// Interpolate (linear) between dipole ends to get b position in
// a given frame at given y.

Vec4 RopeDipole::bInterpolate(double y, RotBstMatrix rb, double m0) {

  Vec4 bb1 = d1.getParticlePtr()->vProd();
  Vec4 bb2 = d2.getParticlePtr()->vProd();
  bb1.rotbst(rb);
  bb2.rotbst(rb);
  double y1 = d1.rap(m0);
  double y2 = d2.rap(m0);
  return bb1 + y * (bb2 - bb1) / (y2 - y1);

}

//--------------------------------------------------------------------------

// Calculate the amount of overlapping dipoles at a given rapidity.
// Return the number of overlapping dipoles to the "left" and "right".

pair<int, int> RopeDipole::getOverlaps(double yfrac, double m0, double r0) {
  // Transform yfrac to y in dipole rest frame
  if (!hasRotTo) getDipoleRestFrame();
  double yL = d1.rap(m0,rotTo);
  double yS = d2.rap(m0,rotTo);
  double yH = yS + (yL - yS) * yfrac;
  int m = 0, n = 0;
  for (size_t i = 0; i < overlaps.size(); ++i) {
    if (overlaps[i].overlap( yfrac, bInterpolateDip(yH,m0), r0)
      && !overlaps[i].hadronized()) {
        if (overlaps[i].dir > 0) ++m;
        else                     ++n;
    }
  }
  return make_pair(m,n);

}

//==========================================================================

// Exc class.
// It is a helper class to Ropewalk, used to describe a pair of excitations
// needed for shoving. It is kept away from users, as there are a lot of
// raw pointers floating around.

//--------------------------------------------------------------------------

struct Exc {

// The constructor.
Exc(double yIn, double m0In, int iIn, int jIn, int kIn, RopeDipole* dip1In,
  RopeDipole* dip2In) : y(yIn), m0(m0In), i(iIn), j(jIn), k(kIn), pp1(NULL),
  pp2(NULL), dip1(dip1In), dip2(dip2In) { }

// Set pointers to the two excitation particles.
void setParticlePtrs(Particle* p1, Particle* p2) {
  pp1 = p1;
  pp2 = p2;
  // Set a pointer to the excitation in the dipole.
  dip1->addExcitation(y, pp1);
  dip2->addExcitation(y, pp2);
}

// Give the excitation a kick in the x and y direction,
void shove(double dpx, double dpy) {
  // The old momenta.
  Vec4 p2 = pp2->p();
  Vec4 p1 = pp1->p();
  // The new momenta, after the shove.
  double mt2new = sqrt(pow2(p2.px() - dpx)  + pow2(p2.py() - dpy));
  double e2new  = mt2new * cosh(y);
  double p2znew = mt2new * sinh(y);
  Vec4 p2new(p2.px() - dpx, p2.py() - dpy, p2znew, e2new);
  double mt1new = sqrt(pow2(p1.px() + dpx)  + pow2(p1.py() + dpy));
  double e1new  = mt1new * cosh(y);
  double p1znew = mt1new * sinh(y);
  Vec4 p1new(p1.px() + dpx, p1.py() + dpy, p1znew, e1new);
  // The differences between the two.
  Vec4 deltap1 = p1new - p1;
  Vec4 deltap2 = p2new - p2;
  // Now check if we can add these two excitations to the dipoles.
  if ( dip2->recoil(deltap2) ) {
    if ( dip1->recoil(deltap1) ) {
      pp1->p(p1new);
      pp2->p(p2new);
    } else {
      Vec4 dp2 = -deltap2;
      dip2->recoil(dp2);
    }
  }
}

// The push direction as a four vector.
Vec4 direction() {
  return dip1->bInterpolateDip(y,m0) -
    dip2->bInterpolateDip(y,m0);
}

// Member variables, slice rapidity and small cut-off mass.
double y;
double m0;

// Local particle indices.
int i, j, k;
Particle* pp1;
Particle* pp2;
RopeDipole* dip1;
RopeDipole* dip2;
};

//==========================================================================

// Ropewalk class.
// This class keeps track of all the strings making up ropes for shoving
// as well as flavour enhancement.

//--------------------------------------------------------------------------

// The Ropewalk init function sets parameters and pointers.

bool Ropewalk::init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn) {

  // Save pointers.
  infoPtr = infoPtrIn;
  rndmPtr = rndmPtrIn;

  // Parameters of the ropewalk.
  doShoving            = settings.flag("Ropewalk:doShoving");
  shoveMiniStrings     = settings.flag("Ropewalk:shoveMiniStrings");
  shoveJunctionStrings = settings.flag("Ropewalk:shoveJunctionStrings");
  shoveGluonLoops      = settings.flag("Ropewalk:shoveGluonLoops");
  limitMom             = settings.flag("Ropewalk:limitMom");
  mStringMin           = settings.parm("HadronLevel:mStringMin");
  r0                   = settings.parm("Ropewalk:r0");
  m0                   = settings.parm("Ropewalk:m0");
  pTcut                = settings.parm("Ropewalk:pTcut");
  rCutOff              = settings.parm("Ropewalk:rCutOff");
  gAmplitude           = settings.parm("Ropewalk:gAmplitude");
  gExponent            = settings.parm("Ropewalk:gExponent");
  deltay               = settings.parm("Ropewalk:deltay");
  deltat               = settings.parm("Ropewalk:deltat");
  tShove               = settings.parm("Ropewalk:tShove");
  tInit                = settings.parm("Ropewalk:tInit");
  showerCut            = settings.parm("TimeShower:pTmin");
  alwaysHighest        = settings.flag("Ropewalk:alwaysHighest");

  // Check consistency.
  if (deltat > tShove) {
    infoPtr->errorMsg("Error in Ropewalk::init: "
    "deltat cannot be larger than tShove");
    return false;
  }
  return true;

}

//--------------------------------------------------------------------------

// Calculate the average string tension of the event, in units of the default
// string tension (ie. 1 GeV/fm), using random walk in colour space.

double Ropewalk::averageKappa() {

  double kap = 0.;
  double nd = 0;
  for (DMap::iterator itr = dipoles.begin(); itr != dipoles.end(); ++itr) {

    // Getting the overlaps: m, n.
    pair<int,int> overlap = itr->second.getOverlaps( rndmPtr->flat(), m0, r0);

    // Overlaps define the number of steps taken in the random walk.
    // We define the present dipole to always point in the p-direction.
    pair<double, double> pq = select( overlap.first + 1, overlap.second,
      rndmPtr);
    double enh = 0.25 * (2. + 2. * pq.first + pq.second);
    kap += (enh > 1.0 ? enh : 1.0);
    nd  += 1.0;
  }
  return kap / nd;

}

//--------------------------------------------------------------------------

// Calculate the effective string tension a fraction yfrac in on the dipole,
// given by indices e1 and e2.

double Ropewalk::getKappaHere(int e1, int e2, double yfrac) {

  multimap< pair<int,int>, RopeDipole >::iterator
    itr = dipoles.find( make_pair(e1,e2) );
  if (itr == dipoles.end()) itr = dipoles.find( make_pair(e2,e1) );
  if (itr == dipoles.end()) return 1.0;
  RopeDipole* d = &(itr->second);
  d->hadronized(true);

  // Get quantum numbers m and n.
  pair<int, int> overlap = d->getOverlaps(yfrac, m0, r0);
  pair<double, double> pq;
  // If we are always in the highest multiplet, we need not do
  // a random walk
  if (alwaysHighest) {
    pq = make_pair(overlap.first + 1, overlap.second);
  }
  // Invoke default random walk procedure.
  else {
    pq = select(overlap.first + 1, overlap.second, rndmPtr);
  }
  // Calculate enhancement factor.
  double enh = 0.25 * (2. + 2. * pq.first + pq.second);
  return (enh > 1.0 ? enh : 1.0);

}

//--------------------------------------------------------------------------

// Calculate all overlaps of all dipoles and store as OverlappingRopeDipoles.

bool Ropewalk::calculateOverlaps() {

  // Go through all dipoles.
  for (multimap< pair<int,int>, RopeDipole>::iterator itr = dipoles.begin();
    itr != dipoles.end(); ++itr ) {
    RopeDipole* d1 = &(itr->second);
    if (d1->dipoleMomentum().m2Calc() < pow2(m0)) continue;

    // RopeDipoles rapidities in dipole rest frame.
    RotBstMatrix dipoleRestFrame = d1->getDipoleRestFrame();
    double yc1 = d1->d1Ptr()->rap(m0, dipoleRestFrame);
    double ya1 = d1->d2Ptr()->rap(m0, dipoleRestFrame);
    if (yc1 <= ya1) continue;

    // Go through all possible overlapping dipoles.
    for (multimap< pair<int,int>, RopeDipole>::iterator itr2 = dipoles.begin();
      itr2 != dipoles.end(); ++itr2) {
      RopeDipole* d2 = &(itr2->second);

      // Skip self and overlaps with miniscule dipoles.
      if (d1 == d2) continue;
      if (d2->dipoleMomentum().m2Calc() < pow2(m0)) continue;

      // Ignore if not overlapping in rapidity.
      OverlappingRopeDipole od(d2, m0, dipoleRestFrame);
      if (min(od.y1, od.y2) > yc1 || max(od.y1, od.y2) < ya1 || od.y1 == od.y2)
        continue;

      d1->addOverlappingDipole(od);

    }
  }
  return true;

}
//--------------------------------------------------------------------------

// Invoke the random walk of colour states.

pair<int, int> Ropewalk::select(int m, int n, Rndm* rndm) {

  // Initial valuesM mm and nn are step counters.
  int p = 0, q = 0;
  int mm = m, nn = n;

  // We must take all steps before terminating.
  while (mm + nn > 0) {

    // Take randomly a step in one direction.
    if (rndm->flat() < 0.5 && mm > 0) {
      --mm;

      // Calculate the step probabilities.
      double p1 = multiplicity(p + 1, q);
      double p2 = multiplicity(p, q - 1);
      double p3 = multiplicity(p - 1, q + 1);

      // Normalization.
      double sum = p1 + p2 + p3;
      p1 /= sum, p2 /= sum, p3 /= sum;

      // Select a state.
      double pick = rndm->flat();
      if      (pick < p1)      ++p;
      else if (pick < p1 + p2) --q;
      else                     --p, ++q;
    }

    // As above, but for nn.
    else if (nn > 0) {
      --nn;
      double p1 = multiplicity(p, q + 1);
      double p2 = multiplicity(p -1, q);
      double p3 = multiplicity(p + 1, q - 1);
      double sum = p1 + p2 + p3;
      p1 /= sum, p2 /= sum, p3 /= sum;
      double pick = rndm->flat();
      if      (pick < p1)      ++q;
      else if (pick < p1 + p2) --p;
      else                     --q, ++p;
    }
  }

  // Done.
  return make_pair( (p < 0 ? 0 : p), (q < 0 ? 0 : q) );

}

//--------------------------------------------------------------------------

// Shove all dipoles in the event.

void Ropewalk::shoveTheDipoles(Event& event) {

  // Possibility of some initial propagation.
  if ( tInit > 0.0) {
    for (DMap::iterator dItr = dipoles.begin(); dItr != dipoles.end(); ++dItr)
      dItr->second.propagateInit(tInit);
  }

  // The rapidity slices.
  multimap<double, RopeDipole *> rapDipoles;

  // Order the dipoles in max rapidity.
  double ymin = 0;
  double ymax = 0;
  for (DMap::iterator dItr = dipoles.begin(); dItr != dipoles.end(); ++dItr) {
    RopeDipole* dip = &(dItr->second);
    // Order the dipoles in max rapidity.
    rapDipoles.insert( make_pair(dip->maxRapidity(m0), dip) );
    // Find maximal and minimal rapidity to sample.
    if (dip->minRapidity(m0) < ymin) ymin = dip->minRapidity(m0);
    if (dip->maxRapidity(m0) > ymax) ymax = dip->maxRapidity(m0);
  }

  // Do the sampling from flat distribution.
  vector<double> rapidities;
  for (double y = ymin; y < ymax; y += deltay) rapidities.push_back(y);

  // For each value of ySample, we have a vector of excitation pairs.
  map<double, vector<Exc> > exPairs;
  for (int i = 0, N = eParticles.size(); i < N; ++i) eParticles[i].clear();
  eParticles.clear();
  for (int i = 0, N = rapidities.size(); i < N; ++i) {
  // Construct an empty vector of excitation particles.
  eParticles.push_back( vector<Particle>() );

  // Find dipoles sampled in this slice, and store them temporarily.
  double ySample = rapidities[i];
  vector<RopeDipole*> tmp;
  for (multimap<double, RopeDipole*>::iterator
    rItr = rapDipoles.lower_bound(ySample); rItr != rapDipoles.end(); ++rItr) {
    if (rItr->second->minRapidity(m0) < ySample)
      tmp.push_back(rItr->second);
  }

  // Construct excitation particles, one for each sampled dipole in this slice.
  vector<int> eraseDipoles;

  for (int j = 0, M = tmp.size(); j < M; ++j) {
    Vec4 ex;
    // Test if the dipole can bear the excitation.
    if (!tmp[j]->recoil(ex,true) ) {
      eraseDipoles.push_back(j);
    }
  }
  // Erase dipoles which could not bear an excitation.
  for (int j = 0, M = eraseDipoles.size(); j < M; ++j) {
    tmp.erase( tmp.begin() + (eraseDipoles[j]-j) );
  }
  // Add the actual excitations, but only if we can create pairs.
  if( int(tmp.size()) > 1)
    for (int j = 0, M = tmp.size(); j < M; ++j) {
      Vec4 ex;
      // We boost the excitation back from dipole rest frame.
      tmp[j]->recoil(ex,false);
      Particle pp = Particle(21, 22, 0, 0, 0, 0, 0, 0, ex);
      pp.vProd( tmp[j]->bInterpolateLab(ySample,m0) );
      eParticles[i].push_back(pp);
    }
  // Construct all pairs of possible excitations in this slice.
  exPairs[ySample] = vector<Exc>();
  for (int j = 0, M = tmp.size(); j < M; ++j)
    for (int k = 0; k < M; ++k) {
       // Don't allow a string to shove itself.
       if(j != k && tmp[j]->index() != tmp[k]->index() )
         exPairs[ySample].push_back( Exc(ySample, m0, i, j, k, tmp[j],
         tmp[k]) );
    }
  }

  // Give the excitations pointers to the excitation particles.
  for (map<double, vector<Exc> >::iterator slItr = exPairs.begin();
    slItr != exPairs.end(); ++slItr) {
    for (int i = 0, N = slItr->second.size(); i < N; ++i) {
      Exc& ep = slItr->second[i];
      ep.setParticlePtrs( &eParticles[ep.i][ep.j], &eParticles[ep.i][ep.k] );
    }
  }

  // Shoving loop.
  for (double t = tInit; t < tShove + tInit; t += deltat) {
    // For all slices.
    for (map<double, vector<Exc> >::iterator slItr = exPairs.begin();
      slItr != exPairs.end(); ++slItr)
      // For all excitation pairs.
      for (int i = 0, N = slItr->second.size(); i < N; ++i) {
        Exc& ep = slItr->second[i];
        // The direction vector is a space-time four-vector.
        Vec4 direction = ep.direction();
        // The string radius is time dependent,
        // growing with the speed of light.
        // Minimal string size is 1 / shower cut-off
        // converted to fm.
        double rt = max(t, 1. / showerCut / 5.068);
        rt = min(rt, r0 * gExponent);
        double dist = direction.pT();
        // Calculate the push, its direction and do the shoving.
        if (dist < rCutOff * rt) {
          // Gain function.
          double gain = 0.5 * deltay * deltat * gAmplitude * dist / rt / rt
                      * exp( -0.25 * dist * dist / rt / rt);
          double dpx = dist > 0.0 ? gain * direction.px() / dist: 0.0;
          double dpy = dist > 0.0 ? gain * direction.py() / dist: 0.0;
          ep.shove(dpx, dpy);
        }
      }

    // Propagate the dipoles.
    for (DMap::iterator dItr = dipoles.begin(); dItr != dipoles.end(); ++dItr)
      dItr->second.propagate(deltat, m0);
  }

  // Add the excitations to the dipoles.
  for (DMap::iterator dItr = dipoles.begin(); dItr != dipoles.end(); ++dItr) {
    RopeDipole* dip = &(dItr->second);
    if (dip->nExcitations() > 0) dip->excitationsToString( m0, event);
  }
}

//--------------------------------------------------------------------------

// Extract all dipoles from an event.

bool Ropewalk::extractDipoles(Event& event, ColConfig& colConfig) {

  // Go through all strings in the event.
  dipoles.clear();
  for (int iSub = 0; iSub < colConfig.size(); ++iSub) {

    // We can exclude loops, junctions and ministrings from the Ropewalk.
    if (colConfig[iSub].hasJunction && !shoveJunctionStrings) continue;
    if (colConfig[iSub].isClosed && !shoveGluonLoops) continue;
    if (colConfig[iSub].massExcess <= mStringMin && !shoveMiniStrings)
      continue;

    colConfig.collect(iSub,event);
    vector<int> stringPartons = colConfig[iSub].iParton;
    vector<RopeDipole> stringDipole;
    bool stringStart = true;
    RopeDipoleEnd previous;
    for (int iPar = int(stringPartons.size() - 1); iPar > -1; --iPar) {
      if (stringPartons[iPar] > 0) {
        // Ordinary particle.
        RopeDipoleEnd next( &event, stringPartons[iPar]);
        // If we are at the first parton, no dipole.
        if ( !stringStart) {
          // Get the parton placement in Event Record.
          pair<int,int> dipoleER = make_pair( stringPartons[iPar + 1],
            stringPartons[iPar] );
          RopeDipole test(previous, next, iSub, infoPtr);
          if ( limitMom && test.dipoleMomentum().pT() < pTcut)
            dipoles.insert( pair< pair<int, int>, RopeDipole>(dipoleER,
              RopeDipole( previous, next, iSub, infoPtr)) );
          else if (!limitMom)
            dipoles.insert( pair< pair<int, int>, RopeDipole>(dipoleER,
              RopeDipole( previous, next, iSub, infoPtr)));
        }
        previous = next;
        stringStart = false;
      }
      else continue;
    }
  // We have created all dipoles.
  }
  return true;

}

//==========================================================================

// RopeFragPars recalculates fragmentation parameters according to a
// changed string tension. Helper class to FlavourRope.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.

// Initial step size for a calculation.
const double RopeFragPars::DELTAA = 0.1;

// Convergence criterion for a calculation.
const double RopeFragPars::ACONV = 0.001;

// Low z cut-off in fragmentation function.
const double RopeFragPars::ZCUT = 1.0e-4;

//--------------------------------------------------------------------------

// The init function sets up initial parameters from settings.

void RopeFragPars::init(Info* infoPtrIn, Settings& settings) {

  // Info pointer.
  infoPtr = infoPtrIn;

  // The junction parameter.
  beta = settings.parm("Ropewalk:beta");

  // Initialize default values from input settings.
  const int len = 9;
  string params [len] = {"StringPT:sigma", "StringZ:aLund",
  "StringZ:aExtraDiquark","StringZ:bLund", "StringFlav:probStoUD",
  "StringFlav:probSQtoQQ", "StringFlav:probQQ1toQQ0", "StringFlav:probQQtoQ",
  "StringFlav:kappa"};
  double* variables[len] = {&sigmaIn, &aIn, &adiqIn, &bIn, &rhoIn, &xIn,
    &yIn, &xiIn, &kappaIn};
  for (int i = 0; i < len; ++i) *variables[i] = settings.parm(params[i]);

  // Insert the h = 1 case immediately.
  sigmaEff = sigmaIn, aEff = aIn, adiqEff = adiqIn, bEff = bIn,
    rhoEff = rhoIn, xEff = xIn, yEff = yIn, xiEff = xiIn, kappaEff = kappaIn;
  if (!insertEffectiveParameters(1.0)) infoPtr->errorMsg(
    "Error in RopeFragPars::init: failed to insert defaults.");

}

//--------------------------------------------------------------------------

// Return parameters at given string tension, ordered by their
// name for easy insertion in settings.

map<string,double> RopeFragPars::getEffectiveParameters(double h) {

  map<double, map<string, double> >::iterator parItr = parameters.find(h);

  // If the parameters are already calculated, return them.
  if ( parItr != parameters.end()) return parItr->second;

  // Otherwise calculate them.
  if (!calculateEffectiveParameters(h))
    infoPtr->errorMsg("Error in RopeFragPars::getEffectiveParameters:"
      " calculating effective parameters.");

  // And insert them.
  if (!insertEffectiveParameters(h))
    infoPtr->errorMsg("Error in RopeFragPars::getEffectiveParameters:"
      " inserting effective parameters.");

  // And recurse.
  return getEffectiveParameters(h);

}

//--------------------------------------------------------------------------

// Get the Fragmentation function a parameter from cache or calculate it.

double RopeFragPars::getEffectiveA(double thisb, double mT2, bool isDiquark) {

  // Check for  the trivial case.
  if (thisb == bIn) return (isDiquark ? aIn + adiqIn : aIn);

  // We order by b*mT2
  map<double, double>* aMPtr = (isDiquark ? &aDiqMap : &aMap);
  double bmT2 = mT2 * thisb;

  // Check if we have already calculated this a value before.
  map<double,double>::iterator aItr = aMPtr->find(bmT2);
  if (aItr != aMPtr->end()) return aItr->second;

  // Otherwise calculate it.
  double ae = ( isDiquark ? aEffective(aIn + adiqIn, thisb, mT2)
    : aEffective(aIn, thisb, mT2) );
  if (isDiquark) {
    double suba = getEffectiveA(thisb, mT2, false);
    aMPtr->insert( make_pair(bmT2, ae - suba) );
  }
  else aMPtr->insert(make_pair(bmT2, ae));
  return ae;

}

//--------------------------------------------------------------------------

// Calculate the effective parameters.

bool RopeFragPars::calculateEffectiveParameters(double h) {

  if (h <= 0) return false;
  double hinv = 1.0 / h;

  // Start with the easiest transformations.
  // The string tension kappa.
  kappaEff = kappaIn * h;
  // Strangeness.
  rhoEff = pow(rhoIn, hinv);
  // Strange diquarks.
  xEff = pow(xIn, hinv);
  // Spin.
  yEff = pow(yIn, hinv);
  // pT distribution width.
  sigmaEff = sigmaIn * sqrt(h);

  // Derived quantity alpha.
  double alpha = (1 + 2. * xIn * rhoIn + 9. * yIn + 6. * xIn * rhoIn * yIn
    + 3. * yIn * xIn * xIn * rhoIn * rhoIn) / (2. + rhoIn);
  double alphaEff = (1. + 2. * xEff * rhoEff + 9. * yEff
    + 6. * xEff * rhoEff * yEff + 3. * yEff * xEff * xEff * rhoEff * rhoEff)
    / (2. + rhoEff);

  // Baryons.
  xiEff = alphaEff * beta * pow( xiIn / alpha / beta, hinv);
  if (xiEff > 1.0)  xiEff = 1.0;
  if (xiEff < xiIn) xiEff = xiIn;

  // Fragmentation function b.
  bEff = (2. + rhoEff) / (2. + rhoIn) * bIn;
  if (bEff < bIn) bEff = bIn;
  if (bEff > 2.0) bEff = 2.0;

  // We calculate a for a typical particle with mT2 = 1 GeV^2.
  aEff    = getEffectiveA( bEff, 1.0, false);
  adiqEff = getEffectiveA( bEff, 1.0, true) - aEff;

  return true;

}

//--------------------------------------------------------------------------

// Insert calculated parameters in cache for later (re-)use.

bool RopeFragPars::insertEffectiveParameters(double h) {

  map<string,double> p;
  p["StringPT:sigma"]          = sigmaEff;
  p["StringZ:bLund"]           = bEff;
  p["StringFlav:probStoUD"]    = rhoEff;
  p["StringFlav:probSQtoQQ"]   = xEff;
  p["StringFlav:probQQ1toQQ0"] = yEff;
  p["StringFlav:probQQtoQ"]    = xiEff;
  p["StringZ:aLund"]           = aEff;
  p["StringZ:aExtraDiquark"]   = adiqEff;
  p["StringFlav:kappa"]        = kappaEff;

  return (parameters.insert( make_pair(h,p)).second );

}

//--------------------------------------------------------------------------

// Calculate the a parameter.

double RopeFragPars::aEffective(double aOrig, double thisb, double mT2) {

  // Calculate initial normalization constants.
  double N    = integrateFragFun(aOrig,   bIn, mT2);
  double NEff = integrateFragFun(aOrig, thisb, mT2);
  int    s    = (N < NEff) ? -1 : 1;
  double da   = DELTAA;
  double aNew = aOrig - s * da;

  // Iterate until we meet preset convergence criterion.
  do {
    // Calculate normalization with current a.
    NEff = integrateFragFun(aNew, thisb, mT2);
    if ( ((N < NEff) ? -1 : 1) != s ) {
      s = (N < NEff) ? -1 : 1;
      // If we have crossed over the solution, decrease
      // the step size and turn around.
      da /= 10.0;
    }
    aNew -= s * da;
    if (aNew < 0.0) {aNew = 0.1; break;}
    if (aNew > 2.0) {aNew = 2.0; break;}
  } while (da > ACONV);
  return aNew;

}

//--------------------------------------------------------------------------

// The Lund fragmentation function.

double RopeFragPars::fragf(double z, double a, double b, double mT2) {

  if (z < ZCUT) return 0.0;
  return pow(1 - z, a) * exp(-b * mT2 / z) / z;

}

//--------------------------------------------------------------------------

// Integral of the Lund fragmentation function with given parameter values.

double RopeFragPars::integrateFragFun(double a, double b, double mT2) {

  // Using Simpson's rule to integrate the Lund fragmentation function.
  double nextIter, nextComb;
  double thisComb = 0.0, thisIter = 0.0;
  // The target error on the integral should never be changed.
  double error = 1.0e-2;

  // 20 is the max number of iterations, 3 is min. Should not be changed.
  for (int i = 1; i < 20; ++i) {
    nextIter = trapIntegrate( a, b, mT2, thisIter, i);
    nextComb = (4.0 * nextIter - thisIter) / 3.0;
    if (i > 3 && abs(nextComb - thisComb) < error * abs(nextComb))
      return nextComb;
    thisIter = nextIter;
    thisComb = nextComb;
  }
  infoPtr->errorMsg("RopeFragPars::integrateFragFun:"
    "No convergence of frag fun integral.");
  return 0.0;

}

//--------------------------------------------------------------------------

// Helper function for integration.

double RopeFragPars::trapIntegrate( double a, double b, double mT2,
  double sOld, int n) {

  // Compute the nth correction to the integral of fragfunc between 0 and 1
  // using extended trapezoidal rule.
  if (n == 1) return 0.5 * (fragf(0.0, a, b, mT2) + fragf(1.0, a, b, mT2));
  // We want 2^(n-2) interior points (intp). Use bitwise shift to speed up.
  int intp = 1;
  intp <<= n - 2;
  double deltaz = 1.0 / double(intp);
  double z = 0.5 * deltaz;
  double sum = 0.0;
  // Do the integral.
  for (int i = 0; i < intp; ++i, z += deltaz) sum += fragf( z, a, b, mT2);
  return 0.5 * (sOld + sum / double(intp));

}

//==========================================================================

// The FlavourRope class takes care of placing a string breakup in
// the event, and assigning the string breakup effective parameters.

//--------------------------------------------------------------------------

// Change the fragmentation parameters.

bool FlavourRope::doChangeFragPar(StringFlav* flavPtr, StringZ* zPtr,
 StringPT * pTPtr, double m2Had, vector<int> iParton, int endId) {

  // The new parameters.
  map<string, double> newPar;
  if (doBuffon)
    newPar = fetchParametersBuffon(m2Had, iParton, endId);
  else
    newPar = fetchParameters(m2Had, iParton, endId);
  // Change settings to new settings.
  for (map<string, double>::iterator itr = newPar.begin(); itr!=newPar.end();
    ++itr) settingsPtr->parm( itr->first, itr->second);
  // Re-initialize flavour, z, and pT selection with new settings.
  flavPtr->init( *settingsPtr, particleDataPtr, rndmPtr, infoPtr);
  zPtr->init( *settingsPtr, *particleDataPtr, rndmPtr, infoPtr);
  pTPtr->init( *settingsPtr, particleDataPtr, rndmPtr, infoPtr);
  return true;

}

//--------------------------------------------------------------------------

// Find breakup placement and fetch effective parameters using Buffon.

map<string, double> FlavourRope::fetchParametersBuffon(double m2Had,
  vector<int> iParton, int endId) {
  // If effective string tension is set manually, use that.
  if (fixedKappa) return fp.getEffectiveParameters(h);
  if (!ePtr) {
    infoPtr->errorMsg("Error in FlavourRope::fetchParametersBuffon:"
      " Event pointer not set in FlavourRope");
    return fp.getEffectiveParameters(1.0);
  }
    if(find(hadronized.begin(),hadronized.end(),*iParton.begin()) ==
      hadronized.end()){
      hadronized.reserve(hadronized.size() + iParton.size());
      hadronized.insert(hadronized.end(),iParton.begin(),iParton.end());
    }
    // Quark string ends, default mode
   if (endId != 21){
    // Test consistency
    if(ePtr->at(*(iParton.begin())).id() != endId &&
        ePtr->at(*(iParton.end() - 1)).id() != endId) {
      infoPtr->errorMsg("Error in FlavourRope::fetchParametersBuffon:"
        " Quark end inconsistency.");
      return fp.getEffectiveParameters(1.0);
    }

      // First we must let the string vector point in the right direction
      if(ePtr->at(*(iParton.begin())).id() != endId)
       reverse(iParton.begin(),iParton.end());

      // Initialize a bit
      Vec4 hadronic4Momentum(0,0,0,0);
      double enh = 1.0;
      double dipFrac;
      vector<int>::iterator dipItr;
      // Find out when invariant mass exceeds m2Had
      for(dipItr = iParton.begin(); dipItr != iParton.end(); ++dipItr){
       double m2Big = hadronic4Momentum.m2Calc();
        if( m2Had <= m2Big){
          // Approximate the fraction we are in on the dipole, this goes
          // in three cases.
          // We are at the beginning.
          if(m2Had == 0){
            dipFrac = 0;
          }
          // We are somewhere in the first dipole
          else if(dipItr - 1  == iParton.begin()){
            dipFrac = sqrt(m2Had/m2Big);
          }
          else{
            if(ePtr->at(*(dipItr - 1)).id() != 21) {
              infoPtr->errorMsg("Error in FlavourRope::fetchParametersBuffon:"
                " Connecting partons should always be gluons.");
              return fp.getEffectiveParameters(1.0);
            }

            hadronic4Momentum -= 0.5*ePtr->at(*(dipItr -1)).p();
            double m2Small = hadronic4Momentum.m2Calc();

            dipFrac = (sqrt(m2Had) - sqrt(m2Small)) /
              (sqrt(m2Big) - sqrt(m2Small));
          }
          break;
        }
        hadronic4Momentum += ePtr->at(*dipItr).id() == 21 ?
          0.5*ePtr->at(*dipItr).p() : ePtr->at(*dipItr).p();
      }
      // If we reached the end
      // we are in a small string that should just be collapsed
      if(dipItr == iParton.end())
        return fp.getEffectiveParameters(1.0);
      // Sanity check
      if(dipFrac < 0 || dipFrac > 1) {
        infoPtr->errorMsg("Error in FlavourRope::fetchParametersBuffon:"
                " Dipole exceed with fraction less than 0 or greater than 1.");
        return fp.getEffectiveParameters(1.0);
      }
      // We now figure out at what rapidity value,
      // in the lab system, the string is breaking
      double yBreak;
      // Trivial case, just inherit
      if(dipFrac == 0)
        yBreak = ePtr->at(*dipItr).y();
      else{
        // Sanity check
        if(dipItr == iParton.begin()) {
          infoPtr->errorMsg("Error in FlavourRope::fetchParametersBuffon:"
            " We are somehow before the first dipole on a string.");
          return fp.getEffectiveParameters(1.0);
        }
        double dy = ePtr->at(*dipItr).y() - ePtr->at(*(dipItr - 1)).y();
        yBreak = ePtr->at(*(dipItr - 1)).y() + dipFrac*dy;
      }
      // Count the number of partons in the whole
      // event within deltay of breaking point
      double p = 1;
      double q = 0;

      for(int i = 0; i < ePtr->size(); ++i){
        // Don't double count partons from this
        // string (no self-overlap)
        if(find(iParton.begin(),iParton.end(), i) != iParton.end())
          continue;
        // Don't count strings that are already hadronized
        if(find(hadronized.begin(),hadronized.end(),i) != hadronized.end())
          continue;
        double pRap = ePtr->at(i).y();
        if(pRap > yBreak - rapiditySpan && pRap < yBreak + rapiditySpan ){
          // Do a "Buffon" selection to decide whether
          // two strings overlap in impact parameter space
          // given ratio of string diameter to collision diameter
          double r1 = rndmPtr->flat();
          double r2 = rndmPtr->flat();
          double theta1 = 2*M_PI*rndmPtr->flat();
          double theta2 = 2*M_PI*rndmPtr->flat();
          // Overlap?
          if(4*pow2(stringProtonRatio) > pow2(sqrt(r1)*cos(theta1) -
            sqrt(r2)*cos(theta2)) + pow2(sqrt(r1)*sin(theta1) -
            sqrt(r2)*sin(theta2))) {
              if(rndmPtr->flat() < 0.5) p += 0.5;
              else q += 0.5;
          }
        }
      }
      enh = 0.25*(2.0*p+q+2.0);

      return fp.getEffectiveParameters(enh);
    }
   // For closed gluon loops we cannot distinguish the ends.
   // Do nothing instead
   else{
      return fp.getEffectiveParameters(1.0);
   }
      return fp.getEffectiveParameters(1.0);
}

//--------------------------------------------------------------------------
// Find breakup placement and fetch effective parameters using Ropewalk.

map<string, double> FlavourRope::fetchParameters(double m2Had,
  vector<int> iParton, int endId) {
  // If effective string tension is set manually, use that.
  if (fixedKappa) return fp.getEffectiveParameters(h);
  if (!ePtr) {
    infoPtr->errorMsg("Error in FlavourRope::fetchParameters:"
      " Event pointer not set in FlavourRope");
    return fp.getEffectiveParameters(1.0);
  }
  Vec4 mom;
  int eventIndex = -1;
  // Set direction
  bool dirPos;
  if( ePtr->at(iParton[0]).id() == endId) dirPos = true;
  else if( ePtr->at(iParton[iParton.size() - 1]).id() == endId) dirPos = false;
  else {
    infoPtr->errorMsg("Error in FlavourRope::fetchParameters:"
    " Could not get string direction");
    return fp.getEffectiveParameters(1.0);
  }

  for (int i = 0, N = iParton.size(); i < N; ++i) {
    // Change to right direction
    int j = (dirPos ? i : N - 1 - i);
    // Skip the junction entry.
    if ( iParton[j] < 0) continue;
    mom += ePtr->at(iParton[j]).p();
    if ( mom.m2Calc() > m2Had) {
      eventIndex = j;
      break;
    }
  }

  // We figure out where we are on the dipole.
  // Set some values.
  double m2Here = mom.m2Calc();
  // The dipFrac signifies a fraction.
  double dipFrac = 0;
  // We are in the first dipole.
  if (eventIndex == -1 || eventIndex == 0) {
    eventIndex = 0;
    dipFrac = sqrt(m2Had / m2Here);
  }
  else {
    mom -= ePtr->at(iParton[eventIndex]).p();
    double m2Small = mom.m2Calc();
    dipFrac = (sqrt(m2Had) - sqrt(m2Small)) / (sqrt(m2Here) - sqrt(m2Small));
  }
  double enh = rwPtr->getKappaHere( iParton[eventIndex],
    iParton[eventIndex + 1], dipFrac);
  return fp.getEffectiveParameters(enh);

}


//==========================================================================

} // End namespace Pythia8
