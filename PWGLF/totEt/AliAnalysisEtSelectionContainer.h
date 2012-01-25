//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Class for determining physics selection
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________//
#ifndef ALIANALYSISETSELECTIONCONTAINER_H
#define ALIANALYSISETSELECTIONCONTAINER_H

#include "TNamed.h"
#include <map>



class  AliPhysicsSelection; 

class AliAnalysisEtSelectionContainer : public TNamed
{
public:
  
  /** Constructor */
  AliAnalysisEtSelectionContainer();

  /** Constructor */
  AliAnalysisEtSelectionContainer(const char *name);
  
  /** Destructor */
  virtual ~AliAnalysisEtSelectionContainer();
  
  /** Return the physics selection for the current run */
  AliPhysicsSelection* GetPhysicsSelection(Int_t runNumber) { return fPhysicsSelectionMap[runNumber]; }
  
  /** Return the physics selection for the current run */
  AliPhysicsSelection* GetDefaultPhysicsSelection() { return fPhysicsSelectionMap[0]; }
  
  /** Get the map */
  std::map<int, AliPhysicsSelection*> GetPhysicsSelectionMap() const { return fPhysicsSelectionMap; }
  
  /** Add the default selection to the map */
  void AddDefaultSelection(AliPhysicsSelection *selection) { fPhysicsSelectionMap.insert(std::pair<int, AliPhysicsSelection*>(0, selection)); }
  
  /** Add a physics selection to the map */
  void AddPhysicsSelection(AliPhysicsSelection *selection, Int_t runNumber) { fPhysicsSelectionMap.insert(std::pair<int, AliPhysicsSelection*>(runNumber, selection)); }
  
  /** Copy constructor */
  AliAnalysisEtSelectionContainer(const AliAnalysisEtSelectionContainer& other);
  
  /** Assignment operator */
  AliAnalysisEtSelectionContainer& operator=(const AliAnalysisEtSelectionContainer& other);
  
private:
  
  std::map<int, AliPhysicsSelection*> fPhysicsSelectionMap; // The physics selection map

  ClassDef(AliAnalysisEtSelectionContainer, 1);
  
};

#endif // ALIANALYSISETSELECTIONCONTAINER_H
