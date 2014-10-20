// AliFlowEventSimpleFromTTree: Event container for flow analysis of star data
// origin:   Redmer Alexander Bertens, rbertens@cern.ch

#ifndef ALIFLOWEVENTSIMPLEFROMTTREE_H
#define ALIFLOWEVENTSIMPLEFROMTTREE_H_H

// aliroot includes
#include "AliFlowEventSimple.h"

//forward declarations
class AliFlowTTreeEvent;
class TClonesArray;
class AliFlowTrackSimpleCuts;

class AliFlowEventSimpleFromTTree : public AliFlowEventSimple {

   public:
       AliFlowEventSimpleFromTTree();
       AliFlowEventSimpleFromTTree(const AliFlowEventSimpleFromTTree& event );
       AliFlowEventSimpleFromTTree( 
               const AliFlowTTreeEvent* event,
               const TClonesArray* array,
               const AliFlowTrackSimpleCuts* rpCuts=NULL,
               const AliFlowTrackSimpleCuts* poiCuts=NULL);

       Bool_t PassesCuts(const AliFlowTTreeEvent* event);

       AliFlowEventSimpleFromTTree& operator=( const AliFlowEventSimpleFromTTree& event );
       virtual  ~AliFlowEventSimpleFromTTree() {}

       ClassDef(AliFlowEventSimpleFromTTree,1)
};

#endif
