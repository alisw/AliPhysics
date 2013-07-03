#ifndef TOYMCTRACK_H
#define TOYMCTRACK_H

#include <AliExternalTrackParam.h>
#include <TClonesArray.h>
#include <AliTPCclusterMI.h>
#include <AliCluster.h>

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
  Int_t GetNumberOfITSPoints() const { return fITSPoints.GetEntriesFast(); }
  Int_t GetNumberOfTRDPoints() const { return fTRDPoints.GetEntriesFast(); }

  const AliTPCclusterMI* GetSpacePoint(Int_t spoint) const { return static_cast<AliTPCclusterMI*> (fSpacePoints.At(spoint)); }
  const AliTPCclusterMI* GetDistortedSpacePoint(Int_t dspoint) const { return static_cast<AliTPCclusterMI*> (fDistortedSpacePoints.At(dspoint)); }
  const AliCluster* GetITSPoint(Int_t spoint) const { return static_cast<AliCluster*> (fITSPoints.At(spoint)); }
  const AliCluster* GetTRDPoint(Int_t spoint) const { return static_cast<AliCluster*> (fTRDPoints.At(spoint)); }

  AliTPCclusterMI* AddSpacePoint(const AliTPCclusterMI &spoint);
  AliTPCclusterMI* AddDistortedSpacePoint(const AliTPCclusterMI &dspoint);
  AliCluster* AddITSPoint(const AliCluster &spoint);
  AliCluster* AddTRDPoint(const AliCluster &spoint);
 private:

  TClonesArray fSpacePoints;
  TClonesArray fDistortedSpacePoints;
  TClonesArray fITSPoints;
  TClonesArray fTRDPoints;
 
  ClassDef(AliToyMCTrack,2)
};

#endif
