#ifndef TOYMCTRACK_H
#define TOYMCTRACK_H

#include <AliExternalTrackParam.h>
#include <TClonesArray.h>
#include <AliTPCclusterMI.h>

class AliToyMCTrack : public AliExternalTrackParam {
 
 public:
  AliToyMCTrack();
  AliToyMCTrack(Double_t x, Double_t alpha, 
	     const Double_t param[5], 
	     const Double_t covar[15]);
  AliToyMCTrack(Double_t xyz[3],Double_t pxpypz[3],
	     Double_t cv[21],Short_t sign);
  AliToyMCTrack(const AliToyMCTrack &track);
  AliToyMCTrack& operator=(const AliToyMCTrack &track);
  virtual ~AliToyMCTrack() {}
  
  Int_t GetNumberOfSpacePoints()     const { return fSpacePoints.GetEntriesFast(); }
  Int_t GetNumberOfDistSpacePoints() const { return fDistortedSpacePoints.GetEntriesFast(); }

  const AliTPCclusterMI* GetSpacePoint(Int_t spoint) const { return static_cast<AliTPCclusterMI*> (fSpacePoints.At(spoint)); }
  const AliTPCclusterMI* GetDistortedSpacePoint(Int_t dspoint) const { return static_cast<AliTPCclusterMI*> (fDistortedSpacePoints.At(dspoint)); }
  
  AliTPCclusterMI* AddSpacePoint(const AliTPCclusterMI &spoint);
  AliTPCclusterMI* AddDistortedSpacePoint(const AliTPCclusterMI &dspoint);

 private:

  TClonesArray fSpacePoints;
  TClonesArray fDistortedSpacePoints;
 
 
  ClassDef(AliToyMCTrack,1)
};

#endif
