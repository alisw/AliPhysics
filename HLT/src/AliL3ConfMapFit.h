// @(#) $Id$

#ifndef ALIL3_ConfMapFit
#define ALIL3_ConfMapFit

class AliL3ConfMapTrack;
class AliL3Vertex;

class AliL3ConfMapFit {

 private:
  AliL3ConfMapTrack *fTrack; //!
  AliL3Vertex *fVertex; //!
  
 public:
  AliL3ConfMapFit (AliL3ConfMapTrack *track,AliL3Vertex *vertex);
  virtual ~AliL3ConfMapFit() {};

  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  ClassDef(AliL3ConfMapFit,1) //Conformal mapping fit class
};

#endif
