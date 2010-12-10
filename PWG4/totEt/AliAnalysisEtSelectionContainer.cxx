//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Class for determining physics selection
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________//
#include "AliAnalysisEtSelectionContainer.h"
#include "TNamed.h"

ClassImp(AliAnalysisEtSelectionContainer)

AliAnalysisEtSelectionContainer::AliAnalysisEtSelectionContainer() : TNamed("name", "name")
,fPhysicsSelectionMap()
{

}


AliAnalysisEtSelectionContainer::AliAnalysisEtSelectionContainer(const char *name): TNamed(name, name)
,fPhysicsSelectionMap()
{

}

AliAnalysisEtSelectionContainer::AliAnalysisEtSelectionContainer(const AliAnalysisEtSelectionContainer& other): 
  TNamed(other)
  ,fPhysicsSelectionMap(other.GetPhysicsSelectionMap())
{
  // Copy constructor
}

AliAnalysisEtSelectionContainer& AliAnalysisEtSelectionContainer::operator=(const AliAnalysisEtSelectionContainer& other)
{
  // Assignment operator
  if(this != &other)
  {
    fName = other.GetName();
    fPhysicsSelectionMap = other.GetPhysicsSelectionMap();
  }
  return *this;    
}

AliAnalysisEtSelectionContainer::~AliAnalysisEtSelectionContainer()
{
  //Destructor
}
