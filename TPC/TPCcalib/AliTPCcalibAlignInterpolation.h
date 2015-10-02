#ifndef AliTPCCALIBALIGNINTERPOLATION_H
#define AliTPCCALIBALIGNINTERPOLATION_H


#include <AliTPCcalibBase.h>

class TTreeSRedirector;
class THn; 

class AliTPCcalibAlignInterpolation : public AliTPCcalibBase {
public :
  enum EStreamFlags{ // flags to store addition data/code debugging infomation - needed in case something go wrong
    kStreamITSRefit          =0x00001,    // flag: stream ITS refit
    kStreamTRDRefit          =0x00002,    // flag: stream TRD refit
    kStreamTOFRefit          =0x00004,    // flag: stream TOF refit
    kStremInterpolation      =0x00008     // flag: stream (heavy) interpolation information at each pad-row
  };

public :
  AliTPCcalibAlignInterpolation();
  AliTPCcalibAlignInterpolation(const Text_t *name, const Text_t *title, Bool_t onTheFlyFill);
  virtual ~AliTPCcalibAlignInterpolation();
  void ProcessStandalone(const char * inputList);
  virtual void     Process(AliESDEvent *event);
  virtual void     Terminate();
  void   SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  Bool_t RefitITStrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackITS, Double_t &chi2, Double_t &npoints);
  Bool_t RefitTOFtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTOF, Double_t &chi2, Double_t &npoints);
  Bool_t RefitTRDtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTRD, Double_t &chi2, Double_t &npoints);
  void   CreateResidualHistosInterpolation(Double_t dy=5, Double_t dz=5);  
  void   CreateDistortionMapsFromFile(const char * inputFile, const char *outputFile);
  void   SetSyswatchStep(Int_t step){fSyswatchStep=step;} // step with which sys. usage is sampled
  //
  static void   FillHistogramsFromChain(const char * residualList, Double_t dy, Double_t dz, Int_t downscale);
  THn * GetHisITSDRPhi() const {return fHisITSDRPhi;}
  THn * GetHisITSTRDDRPhi() const {return fHisITSTRDDRPhi;}
  THn * GetHisITSTOFDRPhi() const {return fHisITSTOFDRPhi;}
  THn * GetHisITSDZ() const {return fHisITSDZ;}
  THn * GetHisITSTRDDZ() const {return fHisITSTRDDZ;}
  THn * GetHisITSTOFDZ() const {return fHisITSTOFDZ;}
protected:
  static Bool_t SortPointArray(AliTrackPointArray *pointarray, Int_t * sortedIndex);
  Bool_t fOnTheFlyFill;    // flag  - histogram filled on the fly?
  THn * fHisITSDRPhi;      // TPC-ITS residual histograms
  THn * fHisITSTRDDRPhi;   // TPC-ITS+TRD residual histograms
  THn * fHisITSTOFDRPhi;   // TPC-ITS_TOF residual histograms
  THn * fHisITSDZ;         // TPC-ITS residual histograms
  THn * fHisITSTRDDZ;      // TPC-ITS+TRD residual histograms
  THn * fHisITSTOFDZ;      // TPC-ITS_TOF residual histograms
public:
  TTreeSRedirector * fStreamer;  // calibration streamer 
  Int_t fStreamLevel;            // stream level - In mode 0 only basic information needed for calibration  stored (see EStreamFlags )
  Int_t fSyswatchStep;           // dump system resource information after  fSyswatchStep tracks
  Int_t fTrackCounter;           // processed track counter
  AliTPCcalibAlignInterpolation(const AliTPCcalibAlignInterpolation&); // Not implemented
  AliTPCcalibAlignInterpolation& operator=(const AliTPCcalibAlignInterpolation&); // Not implemented   
  ClassDef(AliTPCcalibAlignInterpolation,1);
};


#endif
