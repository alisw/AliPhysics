#include "ToyMCTrack.h"

ClassImp(ToyMCTrack);

ToyMCTrack::ToyMCTrack()
  :AliExternalTrackParam()
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
{
  //default constructor
}
//____________________________________________________
ToyMCTrack::ToyMCTrack(const ToyMCTrack &track)
  : AliExternalTrackParam(track)
  ,fSpacePoints(track.fSpacePoints)
  ,fDistortedSpacePoints(track.fDistortedSpacePoints)
{
  //copy constructor
}
//_____________________________________________________
ToyMCTrack& ToyMCTrack::operator = (const ToyMCTrack &track)
{
  //assignment operator
  if (&track == this) return *this;
  new (this) ToyMCTrack(track);

  return *this;
}
//________________________________________________________________
ToyMCTrack::ToyMCTrack(Double_t x, Double_t alpha, 
		       const Double_t param[5], 
		       const Double_t covar[15])
  :AliExternalTrackParam(x,alpha,param,covar)
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
{
  //create external track parameters from given arguments
}
//________________________________________________________________
ToyMCTrack::ToyMCTrack(Double_t xyz[3],Double_t pxpypz[3],
		       Double_t cv[21],Short_t sign)
  :AliExternalTrackParam(xyz,pxpypz,cv,sign)
  ,fSpacePoints("AliTPCclusterMI")
  ,fDistortedSpacePoints("AliTPCclusterMI")
{
}
//________________________________________________________________
AliTPCclusterMI* ToyMCTrack::AddSpacePoint(const AliTPCclusterMI &spoint)
{
  return new(fSpacePoints[fSpacePoints.GetEntriesFast()]) AliTPCclusterMI(spoint);
}
//________________________________________________________________
AliTPCclusterMI* ToyMCTrack::AddDistortedSpacePoint(const AliTPCclusterMI &spoint)
{
  return new(fDistortedSpacePoints[fDistortedSpacePoints.GetEntriesFast()]) AliTPCclusterMI(spoint);
}
