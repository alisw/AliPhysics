#include <AliESDVertex.h>
#include <AliVertexer.h>

ClassImp(AliVertexer)

//////////////////////////////////////////////////////////////////////
// Base class for primary vertex reconstruction                     //
// AliESDVertexer is a class for full 3D primary vertex finding     //
// derived classes: AliITSVertexer                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliVertexer::AliVertexer() {
  // Default Constructor

    fCurrentVertex  = 0;
    SetDebug();
    SetFirstEvent(0);
    SetLastEvent(0);
}


//______________________________________________________________________
AliVertexer::AliVertexer(const AliVertexer &vtxr) : TObject(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliVertexer","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliVertexer& AliVertexer::operator=(const AliVertexer& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
AliVertexer::~AliVertexer() {
  // Default Destructor
  // The objects pointed by the following pointers are not owned
  // by this class and are not deleted

    fCurrentVertex  = 0;
}
