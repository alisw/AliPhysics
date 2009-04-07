#ifndef AliGenGeVSimEventHeader_h
#define AliGenGeVSimEventHeader_h

//
// Event header for GeVSim event generator
// support event plane and elliptic flow
// in next release will suport full differential 
// directed and elliptic flow
//
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi
// 31 Oct, 2002
//

#include "AliGenEventHeader.h"


class AliGenGeVSimEventHeader: public AliGenEventHeader {

 public:
  
  //Constructors
  AliGenGeVSimEventHeader();
  AliGenGeVSimEventHeader(const char *name);
  ~AliGenGeVSimEventHeader() {}

  //Getters
  Float_t GetEventPlane() const {return fEventPlane;}
  Float_t GetEllipticFlow() const {return fEllipticFlow;}

  //Setters
  void SetEventPlane(Float_t psi);
  void SetEllipticFlow(Float_t v2);

 private:
  
  Float_t fEventPlane;       // event plane in rad.
  Float_t fEllipticFlow;     // elliptic flow (fast solution)

 public:
  ClassDef(AliGenGeVSimEventHeader, 1) // Event Header for GeVSim

};



#endif
