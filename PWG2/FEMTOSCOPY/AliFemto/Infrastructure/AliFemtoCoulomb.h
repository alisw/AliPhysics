/***************************************************************************
 *
 * $Id$
 *
 * Author: Randy Wells, Ohio State, rcwells@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    This is a Coulomb correction class which
 *  1. Reads in the dat from a file
 *  2. Performs a linear interpolation in R and creates any array of interpolations
 *  3. Interpolates in eta and returns the Coulomb correction to user
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.12  2000/10/26 19:48:54  rcwells
 * Added functionality for Coulomb correction of <qInv> in 3D correltions
 *
 * Revision 1.11  2000/08/02 01:25:12  lisa
 * Add Coulomb correction capability to 3D Bertsch-Pratt CorrFctn
 *
 * Revision 1.10  2000/07/16 21:38:22  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers mTrack,mV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.9  2000/05/31 20:12:53  rcwells
 * Modified AliFemtoCoulomb for Id and Log entries
 *
 *
 **************************************************************************/

#ifndef AliFemtoCoulomb_HH
#define AliFemtoCoulomb_HH

#include <stdio.h>
#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoPair.h"
#include "Infrastructure/AliFemtoParticle.h"
#include "TH1D.h"
#include "TH3D.h"

class AliFemtoCoulomb {

public:
  AliFemtoCoulomb();
  AliFemtoCoulomb(const char *readFile, const double& radius, const double& charge);
  virtual ~AliFemtoCoulomb();

  void SetRadius(const double& radius);
  double GetRadius();
  void SetFile(const char *readFile);
  void SetChargeProduct(const double& charge);

  // These have different names so eta/Qinv don't confuse the compiler
  double CoulombCorrect(const double& eta);
  double CoulombCorrect(const double& eta, const double& radius);
  double CoulombCorrect(const AliFemtoPair* pair);
  double CoulombCorrect(const AliFemtoPair* pair, const double& radius);
  double CoulombCorrect(const double& mass, const double& charge,
		        const double& radius, const double& qInv);
  TH1D* CorrectionHistogram(const double& mass1, const double& mass2, const int& nBins, 
				    const double& low, const double& high);
#ifdef __ROOT__
  TH1D* CorrectionHistogram(const TH1D*, const double);
  TH3D* CorrectionHistogram(const TH3D*, const double);
#endif
private:
  double Eta(const AliFemtoPair* pair);                // Calculates eta
  void CreateLookupTable(const double& radius);  // Creates look-up table
  const char* fFile;                             // File to interpolate corrections from    
  double fRadius;                                // Radius from previous iteration
  double fZ1Z2;                                  // Charge product of particles
  double fEta[1000];                             // interpolated Coulomb correction table
  double fCoulomb[1000];                         // interpolated Coulomb correction table
  int fNLines;                                   // Number of Eta's in lookup-table

#ifdef __ROOT__
  ClassDef(AliFemtoCoulomb, 0)
#endif
};


#endif
