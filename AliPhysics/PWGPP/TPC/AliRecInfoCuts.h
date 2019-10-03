#ifndef ALIRECINFOCUTS_H
#define ALIRECINFOCUTS_H

//------------------------------------------------------------------------------
// Class to keep selection cuts for reconstructed tracks. 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include "AliESDtrackCuts.h"

class AliVTrack;
class AliVEvent;
class TRootIOCtor;

class AliRecInfoCuts : public AliESDtrackCuts
{
public:
  AliRecInfoCuts(TRootIOCtor*);
  AliRecInfoCuts(const Char_t* name="", const Char_t *title ="");
  AliRecInfoCuts(const AliRecInfoCuts& obj);
  virtual ~AliRecInfoCuts() {;}
 
  // setters 
  void SetMinTPCsignalN(const Int_t min=0) 	 {fMinTPCsignalN = min;}
  void SetMaxAbsTanTheta(const Float_t max=1e99) {fMaxAbsTanTheta = max;}
  void SetMinNClustersTRD(const Int_t min=0)  {fMinNClustersTRD = min;}
  void SetMinNTrackletsTRD(const Int_t min=0)  {fMinNTrackletsTRD = min;}
  void SetTPCITSMatchingRadius(const Float_t radius=70.)  {fTPCITSMatchingRadius = radius;}
  void SetTPCTRDMatchingRadius(const Float_t radius=260.) {fTPCTRDMatchingRadius = radius;}

  // getters
  Int_t GetMinTPCsignalN()    const {return fMinTPCsignalN;}
  Float_t GetMaxAbsTanTheta() const {return fMaxAbsTanTheta;}

  // getters for selected AliESDtrackCuts data members
  Float_t GetPtMin()          const {return fPtMin;}
  Float_t GetPtMax()          const {return fPtMax;}
  Int_t GetMinNClustersTPC()  const {return fCutMinNClusterTPC;}
  Int_t GetMinNClustersITS()  const {return fCutMinNClusterITS;}
  Int_t GetMinNClustersTRD()  const {return fMinNClustersTRD;}
  Int_t GetMinNTrackletsTRD()  const {return fMinNTrackletsTRD;}
  Float_t GetTPCITSMatchingRadius()  const {return fTPCITSMatchingRadius;}
  Float_t GetTPCTRDMatchingRadius()  const {return fTPCTRDMatchingRadius;}
  Bool_t AcceptFTrack(AliVTrack *const vTrack, AliVEvent *const vEvent);
    
  // cuts init function
  void InitME();

private:
  Int_t   fMinTPCsignalN;  // min. number of TPC hits
  Float_t fMaxAbsTanTheta; // max. absolute value of tan(theta)
  Int_t   fMinNClustersTRD; // min number of TRD clusters
  Float_t fTPCITSMatchingRadius; // TPC-ITS matching radius
  Float_t fTPCTRDMatchingRadius; // TPC-TRD matching radius
  Int_t   fMinNTrackletsTRD; // min number of TRD tracklets

  ClassDef(AliRecInfoCuts, 1)
};

#endif //ALIRECINFOCUTS_H
