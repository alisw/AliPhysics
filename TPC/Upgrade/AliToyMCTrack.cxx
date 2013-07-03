#include "AliToyMCTrack.h"

ClassImp(AliToyMCTrack);

AliToyMCTrack::AliToyMCTrack()
  :AliExternalTrackParam()
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
  ,fITSPoints("AliCluster")
  ,fTRDPoints("AliCluster")
{
  //default constructor
}
//____________________________________________________
AliToyMCTrack::AliToyMCTrack(const AliToyMCTrack &track)
  : AliExternalTrackParam(track)
  ,fSpacePoints(track.fSpacePoints)
  ,fDistortedSpacePoints(track.fDistortedSpacePoints)
  ,fITSPoints(track.fITSPoints)
  ,fTRDPoints(track.fTRDPoints)
{
  //copy constructor
}
//_____________________________________________________
AliToyMCTrack& AliToyMCTrack::operator = (const AliToyMCTrack &track)
{
  //assignment operator
  if (&track == this) return *this;
  new (this) AliToyMCTrack(track);

  return *this;
}
//________________________________________________________________
AliToyMCTrack::AliToyMCTrack(Double_t x, Double_t alpha, 
		       const Double_t param[5], 
		       const Double_t covar[15])
  :AliExternalTrackParam(x,alpha,param,covar)
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
  ,fITSPoints("AliCluster")
  ,fTRDPoints("AliCluster")
{
  //create external track parameters from given arguments
}
//________________________________________________________________
AliToyMCTrack::AliToyMCTrack(Double_t xyz[3],Double_t pxpypz[3],
		       Double_t cv[21],Short_t sign)
  :AliExternalTrackParam(xyz,pxpypz,cv,sign)
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
  ,fITSPoints("AliCluster")
  ,fTRDPoints("AliCluster")
{
}
//________________________________________________________________
AliTPCclusterMI* AliToyMCTrack::AddSpacePoint(const AliTPCclusterMI &spoint)
{
  return new(fSpacePoints[fSpacePoints.GetEntriesFast()]) AliTPCclusterMI(spoint);
}
//________________________________________________________________
AliTPCclusterMI* AliToyMCTrack::AddDistortedSpacePoint(const AliTPCclusterMI &spoint)
{
  return new(fDistortedSpacePoints[fDistortedSpacePoints.GetEntriesFast()]) AliTPCclusterMI(spoint);
}
//________________________________________________________________
AliCluster* AliToyMCTrack::AddITSPoint(const AliCluster &spoint)
{
   return new(fITSPoints[fITSPoints.GetEntriesFast()]) AliCluster(spoint);
}
//________________________________________________________________
AliCluster* AliToyMCTrack::AddTRDPoint(const AliCluster &spoint)
{
  return new(fTRDPoints[fTRDPoints.GetEntriesFast()]) AliCluster(spoint);
}
