// @(#) $Id$

#ifndef ALIL3_ConfMapFit
#define ALIL3_ConfMapFit

class AliHLTConfMapTrack;
class AliHLTVertex;

class AliHLTConfMapFit {

 private:
  AliHLTConfMapTrack *fTrack; //!
  AliHLTVertex *fVertex; //!
  
 public:
  AliHLTConfMapFit (AliHLTConfMapTrack *track,AliHLTVertex *vertex);
  virtual ~AliHLTConfMapFit() {};

  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  ClassDef(AliHLTConfMapFit,1) //Conformal mapping fit class
};

typedef AliHLTConfMapFit AliL3ConfMapFit; // for backward compatibility

#endif
