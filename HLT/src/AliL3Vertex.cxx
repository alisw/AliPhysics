//Author:        Uli Frankenfeld
//Last Modified: 07.11.2000

#include "AliL3RootTypes.h"
#include <iostream.h>
#include "AliL3Logging.h"
#include "AliL3Vertex.h"

//_____________________________________________________________
//
// AliL3Vertex
//
// stores the information of the vertex position
//


ClassImp(AliL3Vertex)
AliL3Vertex::AliL3Vertex(){
  //
  // default constructor for the AliL3Vertex class. 
  //

  //Set vertex to zero.
  SetZero();  
}

AliL3Vertex::~AliL3Vertex(){
  //
  // destructor
  //
}

void AliL3Vertex::SetZero(){
  // doit
  SetX(0);
  SetY(0);
  SetZ(0);
  SetXErr(1);
  SetYErr(1);
  SetZErr(1);
  fR=0;
  fPhi=0;

  fMWxy = 1.;
}

void AliL3Vertex::Read(AliL3VertexData *vertex){
  // doit
  SetX(vertex->fX);
  SetY(vertex->fY);
  SetZ(vertex->fZ);
  SetXErr(vertex->fXErr);
  SetYErr(vertex->fYErr);
  SetZErr(vertex->fZErr);
  fR=0;
  fPhi=0;

  fMWxy = 1.;
}

