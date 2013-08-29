#include "AliToyMCTrack.h"
#include <AliTrackPointArray.h>

ClassImp(AliToyMCTrack);

AliToyMCTrack::AliToyMCTrack()
  :AliExternalTrackParam()
  ,fSpacePoints("AliTPCclusterMI",160)
  ,fDistortedSpacePoints("AliTPCclusterMI",160)
  ,fITSPoints("AliTrackPoint",7)
  ,fTRDPoints("AliTrackPoint",6)
{
  //default constructor
}
//________________________________________________________________
AliToyMCTrack::AliToyMCTrack(const AliToyMCTrack &track)
  : AliExternalTrackParam(track)
  ,fSpacePoints(track.fSpacePoints)
  ,fDistortedSpacePoints(track.fDistortedSpacePoints)
  ,fITSPoints(track.fITSPoints)
  ,fTRDPoints(track.fTRDPoints)
{
  //copy constructor
}

//________________________________________________________________
AliToyMCTrack::AliToyMCTrack(const AliExternalTrackParam &param)
  :AliExternalTrackParam(param)
  ,fSpacePoints("AliTPCclusterMI",160)
  ,fDistortedSpacePoints("AliTPCclusterMI",160)
  ,fITSPoints("AliTrackPoint",7)
  ,fTRDPoints("AliTrackPoint",6)
{
  //
  // initialise from external track param
  //
}
//________________________________________________________________
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
  ,fSpacePoints("AliTPCclusterMI",160)
  ,fDistortedSpacePoints("AliTPCclusterMI",160)
  ,fITSPoints("AliTrackPoint",7)
  ,fTRDPoints("AliTrackPoint",6)
{
  //create external track parameters from given arguments
}
//________________________________________________________________
AliToyMCTrack::AliToyMCTrack(Double_t xyz[3],Double_t pxpypz[3],
		       Double_t cv[21],Short_t sign)
  :AliExternalTrackParam(xyz,pxpypz,cv,sign)
  ,fSpacePoints("AliTPCclusterMI",160)
  ,fDistortedSpacePoints("AliTPCclusterMI",160)
  ,fITSPoints("AliTrackPoint",7)
  ,fTRDPoints("AliTrackPoint",6)
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
AliTrackPoint* AliToyMCTrack::AddITSPoint(const AliTrackPoint &spoint)
{
   return new(fITSPoints[fITSPoints.GetEntriesFast()]) AliTrackPoint(spoint);
}
//________________________________________________________________
AliTrackPoint* AliToyMCTrack::AddTRDPoint(const AliTrackPoint &spoint)
{
  return new(fTRDPoints[fTRDPoints.GetEntriesFast()]) AliTrackPoint(spoint);
}
