#ifndef ALIL3_ConfMapFit
#define ALIL3_ConfMapFit

#include "AliL3RootTypes.h"

class AliL3ConfMapTrack;
class AliL3Vertex;

class AliL3ConfMapFit {

 private:
  AliL3ConfMapTrack *fTrack; //!
  AliL3Vertex *fVertex; //!
  Double_t BFACT;
  Double_t bField;

  static Double_t pi;

 public:
  AliL3ConfMapFit (AliL3ConfMapTrack *track,AliL3Vertex *vertex);
  virtual ~AliL3ConfMapFit() {};

  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  ClassDef(AliL3ConfMapFit,1) //Conformal mapping fit class
};

#endif
