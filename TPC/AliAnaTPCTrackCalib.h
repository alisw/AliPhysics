#ifndef ALIANATPCTRACKCALIB_H
#define ALIANATPCTRACKCALIB_H

// ROOT includes
#include <TH1.h>
#include <TH2.h>

// AliRoot includes
#include "AliAnaTPCTrackBase.h"

#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliTPCseed.h>
#include <AliTPCROC.h>
#include <AliTPCcalibTracks.h>

class AliAnaTPCTrackCalib : public AliAnaTPCTrackBase {
public:
   AliAnaTPCTrackCalib();
   AliAnaTPCTrackCalib(const char *name);
   virtual ~AliAnaTPCTrackCalib() {};
   
   virtual void   CreateOutputObjects();

private:
   virtual Int_t FillTrackHistograms(Int_t nTracks, AliESDtrack* track, AliESDfriendTrack* friendTrack, AliTPCseed* seed);

   TH1I              *fNtracks;     // number of Tracks
   TH1I              *fNClusters;   // number of clusters on track
   AliTPCcalibTracks *fCalibTracks; //Analysis object
  
   ClassDef(AliAnaTPCTrackCalib, 1); // Analysis task example for TPC tracks and clusters
};

#endif
