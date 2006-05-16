#ifndef ALIESDFRIENDTRACK_H
#define ALIESDFRIENDTRACK_H

//-------------------------------------------------------------------------
//                     Class AliESDfriendTrack
//               This class contains ESD track additions
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>

class AliTrackPointArray;
class AliKalmanTrack;
class TObjArrray;
//_____________________________________________________________________________
class AliESDfriendTrack : public TObject {
public:
  enum {
    kMaxITScluster=12,
    kMaxTPCcluster=160,
    kMaxTRDcluster=180
  };
  AliESDfriendTrack();
  AliESDfriendTrack(const AliESDfriendTrack &);
  virtual ~AliESDfriendTrack();

  void Set1P(Float_t p) {f1P=p;}
  void SetTrackPointArray(AliTrackPointArray *points) {
    fPoints=points;
  }
  Float_t Get1P() const  {return f1P;}
  Int_t *GetITSindices() {return fITSindex;}
  Int_t *GetTPCindices() {return fTPCindex;}
  Int_t *GetTRDindices() {return fTRDindex;}
  const AliTrackPointArray *GetTrackPointArray() const {return fPoints;}

  void SetITStrack(AliKalmanTrack *t) {fITStrack=t;}
  void SetTRDtrack(AliKalmanTrack *t) {fTRDtrack=t;}
  AliKalmanTrack *GetTRDtrack() {return fTRDtrack;}
  AliKalmanTrack *GetITStrack() {return fITStrack;}
  void AddCalibObject(TObject * calibObject); 
  TObject * GetCalibObject(Int_t index);
protected:
  Float_t f1P;                     // 1/P (1/(GeV/c))
  Int_t fITSindex[kMaxITScluster]; // indices of the ITS clusters 
  Int_t fTPCindex[kMaxTPCcluster]; // indices of the TPC clusters
  Int_t fTRDindex[kMaxTRDcluster]; // indices of the TRD clusters
  AliTrackPointArray *fPoints;//Array of track space points in the global frame
  TObjArray      *fCalibContainer; //Array of objects for calibration    
  AliKalmanTrack *fITStrack; //! pointer to the ITS track (debug purposes) 
  AliKalmanTrack *fTRDtrack; //! pointer to the TRD track (debug purposes) 

  ClassDef(AliESDfriendTrack,2) //ESD friend track
};

#endif


