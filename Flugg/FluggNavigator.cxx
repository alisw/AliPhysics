//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

// GEANT4 tag $ Name:  $
// 
// class FluggNavigator Implementation  Paul Kent July 95/96

#include "FluggNavigator.hh"
#include "G4ios.hh"
#include <iomanip.h>

#ifdef G4GEOMETRY_DEBUG
# define G4DEBUG_NAVIGATION 1
# define G4VERBOSE 1
#endif

FluggNavigator::FluggNavigator() : 
  G4Navigator()
{
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FluggNavigator constructor" << G4endl;
#endif

  ResetStackAndState();

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FluggNavigator constructor" << G4endl;
#endif
}
    
void FluggNavigator::UpdateNavigatorHistory(const G4NavigationHistory* newNavHistory)
{
#ifdef G4GEOMETRY_DEBUG
  cout << "==> Flugg FluggNavigator::UpdateNavigatorHistory(" << newNavHistory 
       << ")" << endl;
  cout << "\t+fHistory=" << fHistory << ") ..." << G4endl;
#endif

  ResetStackAndState();
  fHistory = *newNavHistory;
  SetupHierarchy();

#ifdef G4GEOMETRY_DEBUG
  cout << "<== Flugg FluggNavigator::UpdateNavigatorHistory(" << newNavHistory 
       << ")" << endl;
#endif
}
