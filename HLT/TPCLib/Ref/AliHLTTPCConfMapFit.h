// @(#) $Id$

#ifndef ALIHLTTPC_ConfMapFit
#define ALIHLTTPC_ConfMapFit

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCConfMapTrack;
class AliHLTTPCVertex;

class AliHLTTPCConfMapFit {

 private:
  AliHLTTPCConfMapTrack *fTrack; //!
  AliHLTTPCVertex *fVertex; //!
  Double_t BFACT;
  
  static Double_t pi;

 public:
  AliHLTTPCConfMapFit (AliHLTTPCConfMapTrack *track,AliHLTTPCVertex *vertex);
  virtual ~AliHLTTPCConfMapFit() {};

  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  ClassDef(AliHLTTPCConfMapFit,1) //Conformal mapping fit class
};

#endif
