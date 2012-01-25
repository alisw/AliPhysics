//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  for handling event selection
//  
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________//
#ifndef ALIANALYSISETSELECTIONHANDLER_H
#define ALIANALYSISETSELECTIONHANDLER_H

#include "TObject.h"
	       //#include "AliAnalysisEtSelectionContainer.h"

class AliPhysicsSelection;
class AliAnalysisEtSelectionContainer;

class AliAnalysisEtSelectionHandler : public TObject
{

public:
    AliAnalysisEtSelectionHandler();    
  
    AliAnalysisEtSelectionHandler(const char *name);
    
    virtual ~AliAnalysisEtSelectionHandler();
    
    AliPhysicsSelection* GetPhysicsSelection(Int_t runNumber);// { return fSelections->GetPhysicsSelection(runNumber); }
    
    AliPhysicsSelection* GetDefaultPhysicsSelection();// { return fSelections->GetDefaultPhysicsSelection(); }
    
    AliAnalysisEtSelectionContainer* GetSelectionContainer() const { return fSelections; }

    AliAnalysisEtSelectionHandler(const AliAnalysisEtSelectionHandler& other);
    
    AliAnalysisEtSelectionHandler& operator=(const AliAnalysisEtSelectionHandler& other);
    
    private:
  
    AliAnalysisEtSelectionContainer *fSelections; //! The selection container
  

    
    
    ClassDef(AliAnalysisEtSelectionHandler, 1);
    
    
};

#endif // ALIANALYSISETSELECTIONHANDLER_H
