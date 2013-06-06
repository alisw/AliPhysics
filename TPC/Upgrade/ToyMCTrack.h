#ifndef TOYMCTRACK_H
#define TOYMCTRACK_H

#include <AliExternalTrackParam.h>
#include <TClonesArray.h>
#include <AliTPCclusterMI.h>

class ToyMCTrack : public AliExternalTrackParam {
 
 public:
  ToyMCTrack();
  ToyMCTrack(Double_t x, Double_t alpha, 
	     const Double_t param[5], 
	     const Double_t covar[15]);
  ToyMCTrack(Double_t xyz[3],Double_t pxpypz[3],
	     Double_t cv[21],Short_t sign);
  ToyMCTrack(const ToyMCTrack &track);
  ToyMCTrack& operator=(const ToyMCTrack &track);
  virtual ~ToyMCTrack() {}
  
  Int_t GetNumberOfSpacePoints()     const { return fSpacePoints.GetEntriesFast(); }
  Int_t GetNumberOfDistSpacePoints() const { return fDistortedSpacePoints.GetEntriesFast(); }

  const AliTPCclusterMI* GetSpacePoint(Int_t spoint) const { return static_cast<AliTPCclusterMI*> (fSpacePoints.At(spoint)); }
  const AliTPCclusterMI* GetDistortedSpacePoint(Int_t dspoint) const { return static_cast<AliTPCclusterMI*> (fDistortedSpacePoints.At(dspoint)); }
  
  AliTPCclusterMI* AddSpacePoint(const AliTPCclusterMI &spoint);
  AliTPCclusterMI* AddDistortedSpacePoint(const AliTPCclusterMI &dspoint);

 private:

  TClonesArray fSpacePoints;
  TClonesArray fDistortedSpacePoints;
 
 
  ClassDef(ToyMCTrack,1)
};

#endif
