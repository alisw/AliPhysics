// @(#) $Id$
// Original: AliHLTConfMapFit.h,v 1.5 2004/07/05 09:03:11 loizides 

#ifndef ALIHLTTPC_ConfMapFit
#define ALIHLTTPC_ConfMapFit

class AliHLTTPCConfMapTrack;
class AliHLTTPCVertex;

class AliHLTTPCConfMapFit {

 private:
  AliHLTTPCConfMapTrack *fTrack; //!
  AliHLTTPCVertex *fVertex; //!
  
 public:
  AliHLTTPCConfMapFit (AliHLTTPCConfMapTrack *track,AliHLTTPCVertex *vertex);
  virtual ~AliHLTTPCConfMapFit() {};

  // helix fit
  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  // straight line fit
  Int_t FitStraightLine();
  Int_t FitLineXY();
  Int_t FitLineSZ();
  
  ClassDef(AliHLTTPCConfMapFit,1) //Conformal mapping fit class
};

#endif
