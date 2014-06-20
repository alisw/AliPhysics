#ifndef ALIVVFRIENDTRACK_H
#define ALIVVFRIENDTRACK_H

//_____________________________________________________________________________
class AliVVfriendTrack {
public:
  AliVVfriendTrack();
  virtual ~AliVVfriendTrack();

  //used in calibration
  virtual AliVVTPCseed* GetTPCseed() const {return NULL;}
  virtual AliVVTRDseed* GetTRDseed() const {return NULL;}
  virtual const AliTrackPointArray *GetTrackPointArray() const {return fPoints;}
  virtual const AliExternalTrackParam * GetITSOut() const {return fITSOut;} 
  virtual const AliExternalTrackParam * GetTPCOut() const {return  fTPCOut;} 
  virtual const AliExternalTrackParam * GetTRDIn()  const {return fTRDIn;} 

private: 
  AliVVfriendTrack(const AliVVfriendTrack &);
  AliVVfriendTrack& operator=(const AliVVfriendTrack& esd);  
};

#endif


