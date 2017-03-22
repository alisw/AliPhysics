// Analysis.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// Sphericity, Thrust, ClusJet, CellJet and SlowJet classes.

#include "Analysis.h"

namespace Pythia8 {

//==========================================================================

// Sphericity class.
// This class finds sphericity-related properties of an event.

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int    Sphericity::NSTUDYMIN     = 2;

// Maximum number of times that an error warning will be printed.
const int    Sphericity::TIMESTOPRINT  = 1;

// Assign mimimum squared momentum in weight to avoid division by zero. 
const double Sphericity::P2MIN         = 1e-20;

// Second eigenvalue not too low or not possible to find eigenvectors.
const double Sphericity::EIGENVALUEMIN = 1e-10;

//--------------------------------------------------------------------------
 
// Analyze event.

bool Sphericity::analyze(const Event& event, ostream& os) {

  // Initial values, tensor and counters zero.
  eVal1 = eVal2 = eVal3 = 0.;
  eVec1 = eVec2 = eVec3 = 0.;
  double tt[4][4];
  for (int j = 1; j < 4; ++j) 
  for (int k = j; k < 4; ++k) tt[j][k] = 0.;
  int nStudy = 0;
  double denom = 0.;

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal()) {
    if (select >  2 &&  event[i].isNeutral() ) continue;
    if (select == 2 && !event[i].isVisible() ) continue;
    ++nStudy;

    // Calculate matrix to be diagonalized. Special cases for speed.
    double pNow[4];
    pNow[1] = event[i].px();
    pNow[2] = event[i].py();
    pNow[3] = event[i].pz();
    double p2Now = pNow[1]*pNow[1] + pNow[2]*pNow[2] + pNow[3]*pNow[3];
    double pWeight = 1.;
    if (powerInt == 1) pWeight = 1. / sqrt(max(P2MIN, p2Now));
    else if (powerInt == 0) pWeight = pow( max(P2MIN, p2Now), powerMod);
    for (int j = 1; j < 4; ++j)   
    for (int k = j; k < 4; ++k) tt[j][k] += pWeight * pNow[j] * pNow[k];
    denom += pWeight * p2Now;
  }

  // Very low multiplicities (0 or 1) not considered.
  if (nStudy < NSTUDYMIN) {
    if (nFew < TIMESTOPRINT) os << " PYTHIA Error in " << 
    "Sphericity::analyze: too few particles" << endl; 
    ++nFew;
    return false;
  }

  // Normalize tensor to trace = 1.
  for (int j = 1; j < 4; ++j) 
  for (int k = j; k < 4; ++k) tt[j][k] /= denom;
 
  // Find eigenvalues to matrix (third degree equation).
  double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3] 
    + tt[2][2] * tt[3][3] - pow2(tt[1][2]) - pow2(tt[1][3]) 
    - pow2(tt[2][3]) ) / 3. - 1./9.;
  double qCoefRt = sqrt( -qCoef);
  double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow2(tt[2][3]) 
    + tt[2][2] * pow2(tt[1][3]) + tt[3][3] * pow2(tt[1][2]) 
    - tt[1][1] * tt[2][2] * tt[3][3] ) 
    + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.; 
  double pTemp = max( min( rCoef / pow3(qCoefRt), 1.), -1.);
  double pCoef = cos( acos(pTemp) / 3.);
  double pCoefRt = sqrt( 3. * (1. - pow2(pCoef)) );
  eVal1 = 1./3. + qCoefRt * max( 2. * pCoef,  pCoefRt - pCoef);
  eVal3 = 1./3. + qCoefRt * min( 2. * pCoef, -pCoefRt - pCoef);
  eVal2 = 1. - eVal1 - eVal3;

  // Begin find first and last eigenvector.
  for (int iVal = 0; iVal < 2; ++iVal) {
    double eVal = (iVal == 0) ? eVal1 : eVal3;

    // If all particles are back-to-back then only first axis meaningful.
    if (iVal > 1 && eVal2 < EIGENVALUEMIN) {
      if (nBack < TIMESTOPRINT) os << " PYTHIA Error in "
      "Sphericity::analyze: particles too back-to-back" << endl; 
      ++nBack;
      return false;
    }

    // Set up matrix to diagonalize.
    double dd[4][4];
    for (int j = 1; j < 4; ++j) {
      dd[j][j] = tt[j][j] - eVal;
      for (int k = j + 1; k < 4; ++k) {
        dd[j][k] = tt[j][k]; 
        dd[k][j] = tt[j][k]; 
      }
    }

    // Find largest = pivotal element in matrix.
    int jMax = 0;
    int kMax = 0;
    double ddMax = 0.;
    for (int j = 1; j < 4; ++j) 
    for (int k = 1; k < 4; ++k) 
    if (abs(dd[j][k]) > ddMax) {
      jMax = j;
      kMax = k;
      ddMax = abs(dd[j][k]);
    }

    // Subtract one row from the other two; find new largest element. 
    int jMax2 = 0;
    ddMax = 0.;
    for (int j = 1; j < 4; ++j) 
    if ( j != jMax) {
      double pivot = dd[j][kMax] / dd[jMax][kMax];
      for (int k = 1; k < 4; ++k) {
        dd[j][k] -= pivot * dd[jMax][k];
        if (abs(dd[j][k]) > ddMax) {
          jMax2 = j;
          ddMax = abs(dd[j][k]);
	}
      } 
    }

    // Construct eigenvector. Normalize to unit length; sign irrelevant.
    int k1 = kMax + 1; if (k1 > 3) k1 -= 3;
    int k2 = kMax + 2; if (k2 > 3) k2 -= 3;
    double eVec[4];
    eVec[k1]   = -dd[jMax2][k2];    
    eVec[k2]   =  dd[jMax2][k1];    
    eVec[kMax] = (dd[jMax][k1] * dd[jMax2][k2]
      - dd[jMax][k2] * dd[jMax2][k1]) / dd[jMax][kMax];
    double length = sqrt( pow2(eVec[1]) + pow2(eVec[2])
      + pow2(eVec[3]) );
 
    // Store eigenvectors.
    if (iVal == 0) eVec1 = Vec4( eVec[1] / length,
      eVec[2] / length, eVec[3] / length, 0.);
    else eVec3 = Vec4( eVec[1] / length,
      eVec[2] / length, eVec[3] / length, 0.);
  }

  // Middle eigenvector is orthogonal to the other two; sign irrelevant.
  eVec2 = cross3( eVec1, eVec3);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Provide a listing of the info.
  
void Sphericity::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Sphericity Listing  -------- \n";
  if (powerInt !=2) os << "      Nonstandard momentum power = " 
     << fixed << setprecision(3) << setw(6) << power << "\n"; 
  os << "\n  no     lambda      e_x       e_y       e_z \n";

  // The three eigenvalues and eigenvectors.
  os << setprecision(5);
  os << "   1" << setw(11) << eVal1 << setw(11) << eVec1.px() 
     << setw(10) << eVec1.py() << setw(10) << eVec1.pz() << "\n";
  os << "   2" << setw(11) << eVal2 << setw(11) << eVec2.px() 
     << setw(10) << eVec2.py() << setw(10) << eVec2.pz() << "\n";
  os << "   3" << setw(11) << eVal3 << setw(11) << eVec3.px() 
     << setw(10) << eVec3.py() << setw(10) << eVec3.pz() << "\n";

  // Listing finished.
  os << "\n --------  End PYTHIA Sphericity Listing  ----" << endl;

}


//==========================================================================

// Thrust class.
// This class finds thrust-related properties of an event.

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int    Thrust::NSTUDYMIN    = 2;

// Maximum number of times that an error warning will be printed.
const int    Thrust::TIMESTOPRINT = 1;

// Major not too low or not possible to find major axis.
const double Thrust::MAJORMIN     = 1e-10;

//--------------------------------------------------------------------------
 
// Analyze event.

bool Thrust::analyze(const Event& event, ostream& os) {

  // Initial values and counters zero.
  eVal1 = eVal2 = eVal3 = 0.;
  eVec1 = eVec2 = eVec3 = 0.;
  int nStudy = 0;
  vector<Vec4> pOrder;
  Vec4 pSum, nRef, pPart, pFull, pMax;

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal()) {
    if (select >  2 &&  event[i].isNeutral() ) continue;
    if (select == 2 && !event[i].isVisible() ) continue;
    ++nStudy;

    // Store momenta. Use energy component for absolute momentum.
    Vec4 pNow = event[i].p();
    pNow.e(pNow.pAbs());
    pSum += pNow;
    pOrder.push_back(pNow);
  }

  // Very low multiplicities (0 or 1) not considered.
  if (nStudy < NSTUDYMIN) {
    if (nFew < TIMESTOPRINT) os << " PYTHIA Error in " << 
    "Thrust::analyze: too few particles" << endl; 
    ++nFew;
    return false;
  }

  // Try all combinations of reference vector orthogonal to two particles.
  for (int i1 = 0; i1 < nStudy - 1; ++i1) 
  for (int i2 = i1 + 1; i2 < nStudy; ++i2) {
    nRef = cross3( pOrder[i1], pOrder[i2]);
    nRef /= nRef.pAbs();
    pPart = 0.;

    // Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < nStudy; ++i) if (i != i1 && i != i2) {
      if (dot3(pOrder[i], nRef) > 0.) pPart += pOrder[i]; 
      else                            pPart -= pOrder[i];  
    }  
    for (int j = 0; j < 4; ++j) {
      if      (j == 0) pFull = pPart + pOrder[i1] + pOrder[i2];
      else if (j == 1) pFull = pPart + pOrder[i1] - pOrder[i2];
      else if (j == 2) pFull = pPart - pOrder[i1] + pOrder[i2];
      else             pFull = pPart - pOrder[i1] - pOrder[i2];
      pFull.e(pFull.pAbs());    
      if (pFull.e() > pMax.e()) pMax = pFull;
    }
  }

  // Maximum gives thrust axis and value.
  eVal1 = pMax.e() / pSum.e();
  eVec1 = pMax / pMax.e();
  eVec1.e(0.);

  // Subtract momentum along thrust axis.
  double pAbsSum = 0.;
  for (int i = 0; i < nStudy; ++i) {
    pOrder[i] -= dot3( eVec1, pOrder[i]) * eVec1;
    pOrder[i].e(pOrder[i].pAbs());
    pAbsSum += pOrder[i].e();
  }
    
  // Simpleminded major and minor axes if too little transverse left.
  if (pAbsSum < MAJORMIN * pSum.e()) {
    if ( abs(eVec1.pz()) > 0.5) eVec2 = Vec4( 1., 0., 0., 0.);
    else                        eVec2 = Vec4( 0., 0., 1., 0.); 
    eVec2 -= dot3( eVec1, eVec2) * eVec1;
    eVec2 /= eVec2.pAbs();
    eVec3  = cross3( eVec1, eVec2);
    return true;
  }

  // Try all reference vectors orthogonal to one particles.
  pMax = 0.;
  for (int i1 = 0; i1 < nStudy; ++i1) {
    nRef = cross3( pOrder[i1], eVec1);
    nRef /= nRef.pAbs();
    pPart = 0.;

    // Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < nStudy; ++i) if (i != i1) {
      if (dot3(pOrder[i], nRef) > 0.) pPart += pOrder[i]; 
      else                            pPart -= pOrder[i];  
    }  
    pFull = pPart + pOrder[i1];
    pFull.e(pFull.pAbs());    
    if (pFull.e() > pMax.e()) pMax = pFull;
    pFull = pPart - pOrder[i1];
    pFull.e(pFull.pAbs());    
    if (pFull.e() > pMax.e()) pMax = pFull;    
  }

  // Maximum gives major axis and value.
  eVal2 = pMax.e() / pSum.e();
  eVec2 = pMax / pMax.e();
  eVec2.e(0.);

  // Orthogonal direction gives minor axis, and from there value.
  eVec3 = cross3( eVec1, eVec2);
  pAbsSum = 0.;
  for (int i = 0; i < nStudy; ++i) 
    pAbsSum += abs( dot3(eVec3, pOrder[i]) );     
  eVal3 = pAbsSum / pSum.e();  

   // Done.
  return true;

}

//--------------------------------------------------------------------------

// Provide a listing of the info.
  
void Thrust::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Thrust Listing  ------------ \n"
     << "\n          value      e_x       e_y       e_z \n";

  // The thrust, major and minor values and related event axes.
  os << setprecision(5);
  os << " Thr" << setw(11) << eVal1 << setw(11) << eVec1.px() 
     << setw(10) << eVec1.py() << setw(10) << eVec1.pz() << "\n";
  os << " Maj" << setw(11) << eVal2 << setw(11) << eVec2.px() 
     << setw(10) << eVec2.py() << setw(10) << eVec2.pz() << "\n";
  os << " Min" << setw(11) << eVal3 << setw(11) << eVec3.px() 
     << setw(10) << eVec3.py() << setw(10) << eVec3.pz() << "\n";

  // Listing finished.
  os << "\n --------  End PYTHIA Thrust Listing  --------" << endl;

}

//==========================================================================

// SingleClusterJet class.
// Simple helper class to ClusterJet for a jet and its contents. 

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Assign minimal pAbs to avoid division by zero.
const double SingleClusterJet::PABSMIN  = 1e-10; 

//--------------------------------------------------------------------------
 
// Distance measures between two SingleClusterJet objects.

double dist2Fun(int measure, const SingleClusterJet& j1, 
  const SingleClusterJet& j2) {

  // JADE distance.
  if (measure == 2) return 2. * j1.pJet.e() * j2.pJet.e() 
    * (1. - dot3( j1.pJet, j2.pJet) / (j1.pAbs * j2.pAbs) );

  // Durham distance.
  if (measure == 3) return 2. * pow2( min( j1.pJet.e(), j2.pJet.e() ) ) 
    * (1. - dot3( j1.pJet, j2.pJet) / (j1.pAbs * j2.pAbs) );

  // Lund distance; "default".
  return (j1.pAbs * j2.pAbs - dot3( j1.pJet, j2.pJet)) 
    * 2. * j1.pAbs * j2.pAbs / pow2(j1.pAbs + j2.pAbs);

}  

//==========================================================================

// ClusterJet class.
// This class performs a jet clustering according to different
// distance measures: Lund, JADE or Durham.

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of times that an error warning will be printed.
const int    ClusterJet::TIMESTOPRINT   = 1;

// Assume the pi+- mass for all particles, except the photon, in one option.
const double ClusterJet::PIMASS        = 0.13957; 

// Assign minimal pAbs to avoid division by zero.
const double ClusterJet::PABSMIN        = 1e-10; 

// Initial pT/m preclustering scale as fraction of clustering one.
const double ClusterJet::PRECLUSTERFRAC = 0.1; 

// Step with which pT/m is reduced if preclustering gives too few jets.
const double ClusterJet::PRECLUSTERSTEP = 0.8;

//--------------------------------------------------------------------------
 
// Analyze event.

bool ClusterJet::analyze(const Event& event, double yScaleIn, 
  double pTscaleIn, int nJetMinIn, int nJetMaxIn, ostream& os) {

  // Input values. Initial values zero.
  yScale  = yScaleIn;
  pTscale = pTscaleIn;
  nJetMin = nJetMinIn;
  nJetMax = nJetMaxIn;
  particles.resize(0);
  jets.resize(0);
  Vec4 pSum;
  distances.clear();

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal()) {
    if (select >  2 &&  event[i].isNeutral() ) continue;
    if (select == 2 && !event[i].isVisible() ) continue;

    // Store them, possibly with modified mass => new energy.
    Vec4 pTemp = event[i].p();
    if (massSet == 0 || massSet == 1) {
      double mTemp = (massSet == 0 || event[i].id() == 22) 
        ? 0. : PIMASS; 
      double eTemp = sqrt(pTemp.pAbs2() + pow2(mTemp));
      pTemp.e(eTemp);
    }
    particles.push_back( SingleClusterJet(pTemp, i) );
    pSum += pTemp;
  }

  // Very low multiplicities not considered.
  nParticles = particles.size();
  if (nParticles < nJetMin) {
    if (nFew < TIMESTOPRINT) os << " PYTHIA Error in " << 
    "ClusterJet::analyze: too few particles" << endl; 
    ++nFew;
    return false;
  }

  // Squared maximum distance in GeV^2 for joining.
  double p2Sum = pSum.m2Calc();
  dist2Join = max( yScale * p2Sum, pow2(pTscale));
  dist2BigMin = 2. * max( dist2Join, p2Sum);

  // Do preclustering if desired and possible. 
  if (doPrecluster && nParticles > nJetMin + 2) {
    precluster();
    if (doReassign) reassign();
  }

  // If no preclustering: each particle is a starting jet.
  else for (int i = 0; i < nParticles; ++i) {
    jets.push_back( SingleClusterJet(particles[i]) );
    particles[i].daughter = i;
  }
 
  // Begin iteration towards fewer jets.
  for ( ; ; ) {
 
    // Find the two closest jets.      
    double dist2Min = dist2BigMin;
    int jMin = 0;
    int kMin = 0;
    for (int j = 0; j < int(jets.size()) - 1; ++j) 
    for (int k = j + 1; k < int(jets.size()); ++k) {
      double dist2 = dist2Fun( measure, jets[j], jets[k]); 
      if (dist2 < dist2Min) {
        dist2Min = dist2;
        jMin = j;
        kMin = k;
      }
    }

    // Stop if no pair below cut and not more jets than allowed. 
    if ( dist2Min > dist2Join  
      && (nJetMax < nJetMin || int(jets.size()) <= nJetMax) ) break;
    
    // Stop if reached minimum allowed number of jets. Else continue.
    if (int(jets.size()) <= nJetMin) break; 

    // Join two closest jets.
    jets[jMin].pJet         += jets[kMin].pJet;
    jets[jMin].pAbs          = max( PABSMIN, jets[jMin].pJet.pAbs());
    jets[jMin].multiplicity += jets[kMin].multiplicity;
    for (int i = 0; i < nParticles; ++i) 
    if (particles[i].daughter == kMin) particles[i].daughter = jMin;

    // Save the last 5 distances.
    distances.push_front(dist2Min);
    if (distances.size() > 5) distances.pop_back();

    // Move up last jet to empty slot to shrink list.
    jets[kMin]               = jets.back();
    jets.pop_back();
    int iEnd                 = jets.size();
    for (int i = 0; i < nParticles; ++i) 
    if (particles[i].daughter == iEnd) particles[i].daughter = kMin;

    // Do reassignments of particles to nearest jet if desired.
    if (doReassign) reassign();
  }

  // Order jets in decreasing energy.
  for (int j = 0; j < int(jets.size()) - 1; ++j) 
  for (int k = int(jets.size()) - 1; k > j; --k) 
  if (jets[k].pJet.e() > jets[k-1].pJet.e()) {
    swap( jets[k], jets[k-1]);
    for (int i = 0; i < nParticles; ++i) {
      if (particles[i].daughter == k) particles[i].daughter = k-1;
      else if (particles[i].daughter == k-1) particles[i].daughter = k;
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Precluster nearby particles to save computer time.
  
void ClusterJet::precluster() {

  // Begin iteration over preclustering scale.
  distPre = PRECLUSTERFRAC * sqrt(dist2Join) / PRECLUSTERSTEP;
  for ( ; ;) {
    distPre *= PRECLUSTERSTEP;
    dist2Pre = pow2(distPre);
    for (int i = 0; i < nParticles; ++i) { 
      particles[i].daughter   = -1;
      particles[i].isAssigned = false;
    }

    // Sum up low-momentum region. Jet if enough momentum.
    Vec4 pCentral;
    int multCentral = 0;
    for (int i = 0; i < nParticles; ++i) 
    if (particles[i].pAbs < 2. * distPre) {
      pCentral    += particles[i].pJet;      
      multCentral += particles[i].multiplicity;      
      particles[i].isAssigned = true;
    }
    if (pCentral.pAbs() > 2. * distPre) { 
      jets.push_back( SingleClusterJet(pCentral) );
      jets.back().multiplicity = multCentral;
      for (int i = 0; i < nParticles; ++i) 
      if (particles[i].isAssigned) particles[i].daughter = 0;
    }

    // Find fastest remaining particle until none left.
    for ( ; ;) {
      int iMax = -1;
      double pMax = 0.;
      for (int i = 0; i < nParticles; ++i) 
      if ( !particles[i].isAssigned && particles[i].pAbs > pMax) {
        iMax = i;
        pMax = particles[i].pAbs;
      }
      if (iMax == -1) break;
 
      // Sum up precluster around it according to distance function.
      Vec4 pPre;
      int multPre = 0;
      int nRemain = 0;
      for (int i = 0; i < nParticles; ++i) 
      if ( !particles[i].isAssigned) {           
        double dist2 = dist2Fun( measure, particles[iMax], 
          particles[i]); 
        if (dist2 < dist2Pre) {
          pPre += particles[i].pJet;
          ++multPre;
          particles[i].isAssigned = true;
          particles[i].daughter   = jets.size();
        } else ++nRemain;
      }
      jets.push_back( SingleClusterJet(pPre) ); 
      jets.back().multiplicity = multPre;

      // Decide whether sensible starting configuration or iterate.  
      if (int(jets.size()) + nRemain < nJetMin) break;
    }
    if (int(jets.size()) >= nJetMin) break;
  }

}

//--------------------------------------------------------------------------

// Reassign particles to nearest jet to correct misclustering.
  
void ClusterJet::reassign() {
 
  // Reset clustered momenta.
  for (int j = 0; j < int(jets.size()); ++j) {
    jets[j].pTemp        = 0.;
    jets[j].multiplicity = 0;
  }

  // Loop through particles to find closest jet.
  for (int i = 0; i < nParticles; ++i) {
    particles[i].daughter = -1;
    double dist2Min = dist2BigMin;
    int jMin = 0;
    for (int j = 0; j < int(jets.size()); ++j) {  
      double dist2 = dist2Fun( measure, particles[i], jets[j]); 
      if (dist2 < dist2Min) {
        dist2Min = dist2;
        jMin = j; 
      } 
    }  
    jets[jMin].pTemp += particles[i].pJet;
    ++jets[jMin].multiplicity;
    particles[i].daughter = jMin;
  }

  // Replace old by new jet momenta.
  for (int j = 0; j < int(jets.size()); ++j) {
    jets[j].pJet = jets[j].pTemp;
    jets[j].pAbs =  max( PABSMIN, jets[j].pJet.pAbs());
  }

  // Check that no empty clusters after reassignments.
  for ( ;  ; ) {

    // If no empty jets then done.
    int jEmpty = -1;
    for (int j = 0; j < int(jets.size()); ++j) 
      if (jets[j].multiplicity == 0) jEmpty = j;
    if (jEmpty == -1) return;

    // Find particle assigned to jet with largest distance to it.
    int iSplit = -1;
    double dist2Max = 0.;
    for (int i = 0; i < nParticles; ++i) {
      int j = particles[i].daughter;
      double dist2 = dist2Fun( measure, particles[i], jets[j]);
      if (dist2 > dist2Max) {
        iSplit = i;
        dist2Max = dist2;
      }
    } 

    // Let this particle form new jet and subtract off from existing.
    int jSplit         = particles[iSplit].daughter;    
    jets[jEmpty]       = SingleClusterJet( particles[iSplit].pJet ); 
    jets[jSplit].pJet -=  particles[iSplit].pJet;
    jets[jSplit].pAbs  =  max( PABSMIN,jets[jSplit].pJet.pAbs());
    particles[iSplit].daughter = jEmpty;
    --jets[jSplit].multiplicity;
  }      

}

//--------------------------------------------------------------------------

// Provide a listing of the info.
  
void ClusterJet::list(ostream& os) const {

  // Header.
  string method = (measure == 1) ? "Lund pT" 
        : ( (measure == 2) ? "JADE m" : "Durham kT" ) ;
  os << "\n --------  PYTHIA ClusterJet Listing, " << setw(9) <<  method 
     << " =" << fixed << setprecision(3) << setw(7) << sqrt(dist2Join) 
     << " GeV  --- \n \n  no  mult      p_x        p_y        p_z    "
     << "     e          m \n";

  // The jets.
  for (int i = 0; i < int(jets.size()); ++i) {
    os << setw(4) << i << setw(6) << jets[i].multiplicity << setw(11) 
       << jets[i].pJet.px() << setw(11) << jets[i].pJet.py() 
       << setw(11) << jets[i].pJet.pz() << setw(11) 
       << jets[i].pJet.e() << setw(11) << jets[i].pJet.mCalc() 
       << "\n";  
  }

  // Listing finished.
  os << "\n --------  End PYTHIA ClusterJet Listing  ---------------"
     << "--------" << endl;
}

//==========================================================================

// CellJet class.
// This class performs a cone jet search in (eta, phi, E_T) space.

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int CellJet::TIMESTOPRINT = 1;

//--------------------------------------------------------------------------
 
// Analyze event.

bool CellJet::analyze(const Event& event, double eTjetMinIn, 
  double coneRadiusIn, double eTseedIn, ostream& ) {

  // Input values. Initial values zero.
  eTjetMin   = eTjetMinIn;
  coneRadius = coneRadiusIn;
  eTseed     = eTseedIn;
  jets.resize(0);
  vector<SingleCell> cells;

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal()) {
    if (select >  2 &&  event[i].isNeutral() ) continue;
    if (select == 2 && !event[i].isVisible() ) continue;

    // Find particle position in (eta, phi, pT) space.
    double etaNow = event[i].eta();
    if (abs(etaNow) > etaMax) continue;
    double phiNow = event[i].phi();
    double pTnow  = event[i].pT();
    int iEtaNow   = max(1, min( nEta, 1 + int(nEta * 0.5 
      * (1. + etaNow / etaMax) ) ) );
    int iPhiNow   = max(1, min( nPhi, 1 + int(nPhi * 0.5
      * (1. + phiNow / M_PI) ) ) );
    int iCell     = nPhi * iEtaNow + iPhiNow;

    // Add pT to cell already hit or book a new cell.
    bool found = false;
    for (int j = 0; j < int(cells.size()); ++j) {
      if (iCell == cells[j].iCell) { 
        found = true;
        ++cells[j].multiplicity;
        cells[j].eTcell += pTnow; 
        continue;
      }
    }
    if (!found) {
      double etaCell = (etaMax / nEta) * (2 * iEtaNow - 1 - nEta);
      double phiCell = (M_PI / nPhi) * (2 * iPhiNow - 1 - nPhi);
      cells.push_back( SingleCell( iCell, etaCell, phiCell, pTnow, 1) );
    }
  }

  // Smear true bin content by calorimeter resolution.
  if (smear > 0 && rndmPtr != 0) 
  for (int j = 0; j < int(cells.size()); ++j) {
    double eTeConv = (smear < 2) ? 1. : cosh( cells[j].etaCell );
    double eBef = cells[j].eTcell * eTeConv; 
    double eAft = 0.;
    do eAft = eBef + resolution * sqrt(eBef) * rndmPtr->gauss();
    while (eAft < 0 || eAft > upperCut * eBef);
    cells[j].eTcell = eAft / eTeConv;
  }

  // Remove cells below threshold for seed or for use at all.
  for (int j = 0; j < int(cells.size()); ++j) { 
    if (cells[j].eTcell < eTseed)    cells[j].canBeSeed = false;
    if (cells[j].eTcell < threshold) cells[j].isUsed    = true;
  }

  // Find seed cell: the one with highest pT of not yet probed ones.
  for ( ; ; ) {
    int jMax = 0;
    double eTmax = 0.;
    for (int j = 0; j < int(cells.size()); ++j) 
    if (cells[j].canBeSeed && cells[j].eTcell > eTmax) {
      jMax = j;
      eTmax = cells[j].eTcell;
    }

    // If too small cell eT then done, else start new trial jet.  
    if (eTmax < eTseed) break;
    double etaCenterNow = cells[jMax].etaCell;
    double phiCenterNow = cells[jMax].phiCell;
    double eTjetNow     = 0.;

    //  Sum up unused cells within required distance of seed.
    for (int j = 0; j < int(cells.size()); ++j) {
      if (cells[j].isUsed) continue;
      double dEta = abs( cells[j].etaCell - etaCenterNow );
      if (dEta > coneRadius) continue;
      double dPhi = abs( cells[j].phiCell - phiCenterNow );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      if (dPhi > coneRadius) continue;
      if (pow2(dEta) + pow2(dPhi) > pow2(coneRadius)) continue;
      cells[j].isAssigned = true;
      eTjetNow += cells[j].eTcell;
    }

    // Reject cluster below minimum ET.
    if (eTjetNow < eTjetMin) {
      cells[jMax].canBeSeed = false; 
      for (int j = 0; j < int(cells.size()); ++j) 
        cells[j].isAssigned = false;

    // Else find new jet properties. 
    } else {
      double etaWeightedNow = 0.;
      double phiWeightedNow = 0.;
      int multiplicityNow   = 0;
      Vec4 pMassiveNow;
      for (int j = 0; j < int(cells.size()); ++j) 
      if (cells[j].isAssigned) {
        cells[j].canBeSeed  = false; 
        cells[j].isUsed     = true; 
        cells[j].isAssigned = false; 
        etaWeightedNow += cells[j].eTcell * cells[j].etaCell;
        double phiCell = cells[j].phiCell; 
        if (abs(phiCell - phiCenterNow) > M_PI) 
          phiCell += (phiCenterNow > 0.) ? 2. * M_PI : -2. * M_PI;
        phiWeightedNow  += cells[j].eTcell * phiCell;
        multiplicityNow += cells[j].multiplicity;
        pMassiveNow     += cells[j].eTcell * Vec4( 
           cos(cells[j].phiCell),  sin(cells[j].phiCell), 
          sinh(cells[j].etaCell), cosh(cells[j].etaCell) );
      } 
      etaWeightedNow /= eTjetNow;
      phiWeightedNow /= eTjetNow; 

      // Bookkeep new jet, in decreasing ET order.
      jets.push_back( SingleCellJet( eTjetNow, etaCenterNow, phiCenterNow,
        etaWeightedNow, phiWeightedNow, multiplicityNow, pMassiveNow) ); 
      for (int i = int(jets.size()) - 1; i > 0; --i) {
        if (jets[i-1].eTjet > jets[i].eTjet) break;
        swap( jets[i-1], jets[i]);
      }
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Provide a listing of the info.
  
void CellJet::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA CellJet Listing, eTjetMin = " 
     << fixed << setprecision(3) << setw(8) << eTjetMin 
     << ", coneRadius = " << setw(5) << coneRadius 
     << "  ------------------------------ \n \n  no    "
     << " eTjet  etaCtr  phiCtr   etaWt   phiWt mult      p_x"
     << "        p_y        p_z         e          m \n";

  // The jets.
  for (int i = 0; i < int(jets.size()); ++i) {
    os << setw(4) << i << setw(10) << jets[i].eTjet << setw(8) 
       << jets[i].etaCenter << setw(8) << jets[i].phiCenter << setw(8) 
       << jets[i].etaWeighted << setw(8) << jets[i].phiWeighted 
       << setw(5) << jets[i].multiplicity << setw(11) 
       << jets[i].pMassive.px() << setw(11) << jets[i].pMassive.py()
       << setw(11) << jets[i].pMassive.pz() << setw(11) 
       << jets[i].pMassive.e() << setw(11)
       << jets[i].pMassive.mCalc() << "\n";  
  }

  // Listing finished.
  os << "\n --------  End PYTHIA CellJet Listing  ------------------"
     << "-------------------------------------------------"
     << endl;
}

//==========================================================================

// SlowJet class.
// This class performs clustering in (y, phi, pT) space.

//--------------------------------------------------------------------------
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int    SlowJet::TIMESTOPRINT = 1;

// Assume the pi+- mass for all particles, except the photon, in one option.
const double SlowJet::PIMASS       = 0.13957; 

// Small number to avoid division by zero.
const double SlowJet::TINY         = 1e-20;

//--------------------------------------------------------------------------
 
// Set up list of particles to analyze, and initial distances.

bool SlowJet::setup(const Event& event) {

  // Initial values zero.
  clusters.resize(0);
  jets.resize(0);
  jtSize = 0;

  // Loop over final particles in the event.
  Vec4   pTemp;
  double mTemp, pT2Temp, mTTemp, yTemp, phiTemp;
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal()) {

    // Always apply selection options for visible or charged particles. 
    if      (chargedOnly &&  event[i].isNeutral() ) continue;
    else if (visibleOnly && !event[i].isVisible() ) continue;

    // Normally use built-in selection machinery.
    if (noHook) {

      // Pseudorapidity cut to describe detector range.
      if (cutInEta    && abs(event[i].eta()) > etaMax) continue;
     
      // Optionally modify mass and energy.
      pTemp = event[i].p();
      mTemp = event[i].m();
      if (modifyMass) {
        mTemp = (massSet == 0 || event[i].id() == 22) ? 0. : PIMASS; 
        pTemp.e( sqrt(pTemp.pAbs2() + mTemp*mTemp) );
      }
    
    // Alternatively pass info to SlowJetHook for decision.
    // User can also modify pTemp and mTemp.
    } else {
      pTemp = event[i].p();
      mTemp = event[i].m();
      if ( !sjHookPtr->include( i, event, pTemp, mTemp) ) continue;
    }

    // Store particle momentum, including some derived quantities.
    pT2Temp  = max( TINY*TINY, pTemp.pT2());
    mTTemp  = sqrt( mTemp*mTemp + pT2Temp);
    yTemp   = (pTemp.pz() > 0) 
            ? log( max( TINY, pTemp.e() + pTemp.pz() ) / mTTemp )
            : log( mTTemp / max( TINY, pTemp.e() - pTemp.pz() ) );
    phiTemp = pTemp.phi();
    clusters.push_back( SingleSlowJet(pTemp, pT2Temp, yTemp, phiTemp, i) );
  }

  // Resize arrays to store distances between clusters.
  origSize = clusters.size();
  clSize = origSize; 
  clLast = clSize - 1; 
  diB.resize(clSize);
  dij.resize(clSize * (clSize - 1) / 2);

  // Loop through particles and find distance to beams.
  for (int i = 0; i < clSize; ++i) {
    if (isAnti)    diB[i] = 1. / clusters[i].pT2;
    else if (isKT) diB[i] = clusters[i].pT2;
    else           diB[i] = 1.;

    // Loop through pairs and find relative distance.
    for (int j = 0; j < i; ++j) {
      dPhi = abs( clusters[i].phi - clusters[j].phi );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      dijTemp = (useStandardR) 
              ? (pow2( clusters[i].y - clusters[j].y) + dPhi*dPhi) / R2
         : 2. * (cosh( clusters[i].y - clusters[j].y) - cos(dPhi)) / R2 ;
      if (isAnti)    dijTemp /= max(clusters[i].pT2, clusters[j].pT2); 
      else if (isKT) dijTemp *= min(clusters[i].pT2, clusters[j].pT2); 
      dij[i*(i-1)/2 + j] = dijTemp;

    // End of original-particle loops.
    }
  }

  // Find first particle pair to join.
  findNext();

  // Done.
  return true;

}

//--------------------------------------------------------------------------
 
// Do one recombination step, possibly giving a jet.

bool SlowJet::doStep() {

  // Fail if no possibility to take a step.
  if (clSize == 0) return false;

  // When distance to beam is smallest the cluster is promoted to jet.
  if (jMin == -1) {

    // Store new jet if its pT is above pTMin.
    if (clusters[iMin].pT2 > pT2jetMin) {
      jets.push_back( SingleSlowJet(clusters[iMin]) );
      ++jtSize;

      // Order jets in decreasing pT.
      for (int i = jtSize - 1; i > 0; --i) {
        if (jets[i].pT2 < jets[i-1].pT2) break;
        swap( jets[i], jets[i-1]);
      } 
    }  
  } 

  // When distance between two clusters is smallest they are joined.
  else {

    // Add iMin cluster to jMin.
    clusters[jMin].p  += clusters[iMin].p;
    clusters[jMin].pT2 = max( TINY*TINY, clusters[jMin].p.pT2());
    double mTTemp  = sqrt(clusters[jMin].p.m2Calc() + clusters[jMin].pT2);
    clusters[jMin].y = (clusters[jMin].p.pz() > 0) 
      ? log( max( TINY, clusters[jMin].p.e() + clusters[jMin].p.pz() ) 
      / mTTemp ) : log( mTTemp 
      / max( TINY, clusters[jMin].p.e() - clusters[jMin].p.pz() ) );
    clusters[jMin].phi = clusters[jMin].p.phi();
    clusters[jMin].mult += clusters[iMin].mult;
    clusters[jMin].idx.insert(clusters[iMin].idx.begin(),
                              clusters[iMin].idx.end());

    // Update distances for and to new jMin.
    if (isAnti)    diB[jMin] = 1. / clusters[jMin].pT2;
    else if (isKT) diB[jMin] = clusters[jMin].pT2;
    else           diB[jMin] = 1.;
    for (int i = 0; i < clSize; ++i) if (i != jMin && i != iMin) {
      dPhi = abs( clusters[i].phi - clusters[jMin].phi );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      dijTemp = (useStandardR) 
              ? (pow2( clusters[i].y - clusters[jMin].y) + dPhi*dPhi) / R2
         : 2. * (cosh( clusters[i].y - clusters[jMin].y) - cos(dPhi)) / R2 ;
      if (isAnti)    dijTemp /= max(clusters[i].pT2, clusters[jMin].pT2); 
      else if (isKT) dijTemp *= min(clusters[i].pT2, clusters[jMin].pT2);
      if (i < jMin) dij[jMin*(jMin-1)/2 + i] = dijTemp;
      else          dij[i*(i-1)/2 + jMin]    = dijTemp;
    }
  }

  // Move up last cluster and distances to vacated position iMin.
  if (iMin < clLast) {
    clusters[iMin] = clusters[clLast];
    diB[iMin] = diB[clLast];
    for (int j = 0; j < iMin; ++j) 
      dij[iMin*(iMin-1)/2 + j] = dij[clLast*(clLast-1)/2 + j];
    for (int j = iMin + 1; j < clLast; ++j)
      dij[j*(j-1)/2 + iMin] = dij[clLast*(clLast-1)/2 + j];
  }
    
  // Shrink cluster list by one.
  clusters.pop_back();
  --clSize;
  --clLast;    

  // Find next cluster pair to join.
  findNext();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Provide a listing of the info.
  
void SlowJet::list(bool listAll, ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA SlowJet Listing, p = " << setw(2) 
     << power << ", R = " << fixed << setprecision(3) << setw(5) << R 
     << ", pTjetMin =" << setw(8) << pTjetMin << ", etaMax = " << setw(6) 
     << etaMax << "  --- \n \n  no      pTjet      y       phi   mult   "
     << "   p_x        p_y        p_z         e          m \n";

  // The jets.
  for (int i = 0; i < jtSize; ++i) {
    os << setw(4) << i << setw(11) << sqrt(jets[i].pT2) << setw(9) 
       << jets[i].y << setw(9) << jets[i].phi << setw(6) 
       << jets[i].mult << setw(11) << jets[i].p.px() << setw(11) 
       << jets[i].p.py() << setw(11) << jets[i].p.pz() << setw(11) 
       << jets[i].p.e() << setw(11) << jets[i].p.mCalc() << "\n";  
  }

  // Optionally list also clusters not yet jets.
  if (listAll && clSize > 0) {
    os << " --------  Below this line follows remaining clusters,"
       << " still pT-unordered  -------------------\n";  
    for (int i = 0; i < clSize; ++i) {
      os << setw(4) << i + jtSize << setw(11) << sqrt(clusters[i].pT2) 
         << setw(9) << clusters[i].y << setw(9) << clusters[i].phi 
         << setw(6) << clusters[i].mult << setw(11) << clusters[i].p.px() 
         << setw(11) << clusters[i].p.py() << setw(11) << clusters[i].p.pz() 
         << setw(11) << clusters[i].p.e() << setw(11) 
         << clusters[i].p.mCalc() << "\n";  
    }
  }

  // Listing finished.
  os << "\n --------  End PYTHIA SlowJet Listing  ------------------"
     << "-------------------------------------" << endl;

}

//--------------------------------------------------------------------------

// Find next cluster pair to join.
  
void SlowJet::findNext() {

  // Find smallest of diB, dij.
  if (clSize > 0) {
    iMin =  0;
    jMin = -1;
    dMin = diB[0];
    for (int i = 1; i < clSize; ++i) {
      if (diB[i] < dMin) {
        iMin = i;
        jMin = -1;
        dMin = diB[i];
      }
      for (int j = 0; j < i; ++j) {
        if (dij[i*(i-1)/2 + j] < dMin) {
          iMin = i;
          jMin = j;
          dMin = dij[i*(i-1)/2 + j];
        }
      }
    }

  // If no clusters left then instead default values.
  } else {
    iMin = -1;
    jMin = -1;
    dMin = 0.;
  } 

}

//==========================================================================

} // end namespace Pythia8
