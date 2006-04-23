#ifndef ALIESDFRIENDTRACK_H
#define ALIESDFRIENDTRACK_H

//-------------------------------------------------------------------------
//                     Class AliESDfriendTrack
//               This class contains ESD track additions
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include "AliESDtrack.h"

class AliTrackPointArray;

//_____________________________________________________________________________
class AliESDfriendTrack : public TObject {
public:
  AliESDfriendTrack();
  AliESDfriendTrack(const AliESDfriendTrack &);
  AliESDfriendTrack(const AliESDtrack &);
  virtual ~AliESDfriendTrack();

  Float_t Get1P() const {return f1P;}
  const Int_t *GetITSindices() const {return fITSindex;}
  const Int_t *GetTPCindices() const {return fTPCindex;}
  const Int_t *GetTRDindices() const {return fTRDindex;}
  const AliTrackPointArray *GetTrackPointArray() const {return fPoints;}

protected:
  Float_t f1P;                                  // 1/P (1/(GeV/c))
  Int_t fITSindex[AliESDtrack::kMaxITScluster]; // indices of the ITS clusters 
  Int_t fTPCindex[AliESDtrack::kMaxTPCcluster]; // indices of the TPC clusters
  Int_t fTRDindex[AliESDtrack::kMaxTRDcluster]; // indices of the TRD clusters

  AliTrackPointArray *fPoints; // Array of track space points in the global frame
  ClassDef(AliESDfriendTrack,1) //ESD friend track
};

#endif


