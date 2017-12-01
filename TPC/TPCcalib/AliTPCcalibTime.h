#ifndef ALITPCCALIBTIME_H
#define ALITPCCALIBTIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
Comments to be written here:
What do we calibrate.
  Time dependence of gain and drift velocity in order to account for changes in: temperature, pressure, gas composition.
*/


#include "AliTPCcalibBase.h"
#include "THnSparse.h"
//#include "TH1D.h"                // Temporary make code compiling for HLT in the 

class TObjArray;
class TH1F;
class TH3F;
class TH2F;
class TH1D;
class TList;
class AliVEvent;
class AliVTrack;
class AliExternalTrackParam;
class AliTPCcalibLaser;
class TGraphErrors;
class AliSplineFit;
class AliVfriendTrack;
class AliVfriendEvent;

class AliTPCcalibTime:public AliTPCcalibBase {
public:
  AliTPCcalibTime(); 
  AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeVdrift, Int_t memoryMode=2);
  virtual ~AliTPCcalibTime();
  
  virtual void           Process(AliVEvent *event);
  virtual Long64_t       Merge(TCollection *const li);
  virtual void           Analyze();
  //static Bool_t          IsLaser      (const AliESDEvent *const event) const;
  //static Bool_t          IsCosmics    (const AliESDEvent *const event) const;
  //static Bool_t          IsBeam       (const AliESDEvent *const event) const;
  void                   ProcessLaser (AliVEvent *event);
  void                   ProcessCosmic(const AliVEvent *const event);
  void                   ProcessBeam  (const AliVEvent *const event);

  Bool_t                 IsPair(const AliExternalTrackParam *tr0, const AliExternalTrackParam *tr1);
  Bool_t                 IsCross(const AliVTrack *const tr0, const AliVTrack *const tr1);
  Bool_t                 IsSame (const AliVTrack *const tr0, const AliVTrack *const tr1);
  void                   ProcessSame(const AliVTrack *const track, AliVfriendTrack *const friendTrack, const AliVEvent *const event);
  void                   ProcessAlignITS(AliVTrack *const track, const AliVfriendTrack *const friendTrack, const AliVEvent *const event, AliVfriendEvent *const vFriend);
  void                   ProcessAlignTRD(AliVTrack* const track, AliVfriendTrack *const friendTrack);
  void                   ProcessAlignTOF(AliVTrack* const track, const AliVfriendTrack *const friendTrack);

  THnSparse*    GetHistVdriftLaserA(Int_t index=1) const {return fHistVdriftLaserA[index];}
  THnSparse*    GetHistVdriftLaserC(Int_t index=1) const {return fHistVdriftLaserC[index];}
  THnSparse*    GetHistoDrift(const char* name) const;
  TObjArray*    GetHistoDrift() const;
  TGraphErrors* GetGraphDrift(const char* name);
  TObjArray*    GetGraphDrift();
  AliSplineFit* GetFitDrift(const char* name);
//  TObjArray*    GetFitDrift();
  TH1F*         GetCosmiMatchingHisto(Int_t index=0) const {return fCosmiMatchingHisto[index];}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  TObjArray* GetAlignITSTPC() const {return fAlignITSTPC;}              // alignemnt array ITS TPC match
  TObjArray* GetAlignTRDTPC() const {return fAlignTRDTPC;}              // alignemnt array TRD TPC match
  TObjArray* GetAlignTOFTPC() const {return fAlignTOFTPC;}              // alignemnt array TOF TPC match

  THnSparse * GetTPCVertexHisto(Int_t index) const { return fTPCVertex[index%12];}
  THnSparse * GetTPCVertexHistoCorrelation(Int_t index) const { return fTPCVertexCorrelation[index%5];}

  THnSparse*  GetResHistoTPCCE(Int_t index) const { return (index<5) ? fResHistoTPCCE[index]:0;}        //TPC-CE    matching map
  THnSparse*  GetResHistoTPCITS(Int_t index) const { return (index<5) ? fResHistoTPCITS[index]:0;}        //TPC-ITS    matching map
  THnSparse*  GetResHistoTPCvertex(Int_t index)      const { return (index<5) ? fResHistoTPCvertex[index]   :0;}        //TPC vertex matching map
  THnSparse*  GetResHistoTPCTRD(Int_t index)   const { return (index<5) ? fResHistoTPCTRD[index]:0;}        //TPC-TRD    matching map
  THnSparse*  GetResHistoTPCTOF(Int_t index)   const { return (index<5) ? fResHistoTPCTOF[index]:0;}        //TPC-TOF    matching map

  void        BookDistortionMaps();      // book histograms
  void        FillResHistoTPCCE(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pTPCOut );       // fill residual histo
  void        FillResHistoTPCITS(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pITSOut );       // fill residual histo
  void        FillResHistoTPC(const AliVTrack * pTrack);
  void        FillResHistoTPCTRD(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTRDIn );
  void        FillResHistoTPCTOF(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTOFIn );

  TObjArray * GetLaserArrayA() const { return fArrayLaserA;}
  TObjArray * GetLaserArrayC() const { return fArrayLaserC;}

  Int_t GetEntries() const {if (fResHistoTPCTOF[0]) return fResHistoTPCITS[0]->GetEntries(); else return 0;}

  void   SetCutTracks(Int_t maxTracks)  { fCutTracks = maxTracks; }  // set maximal number of tracks
  Int_t  GetCutTracks() const { return fCutTracks; }    // retun maximal number of tracks
  void   SetMinPt(Double_t m) { fMinPt = m; }
  Double_t GetMinPt() const { return fMinPt;}
  void   SetMinPtITSTPCalign (Double_t m) { fMinPtITSTPCalign=m; }
  Double_t GetMinPtITSTPCalign () const { return fMinPtITSTPCalign; }

  static Double_t fgResHistoMergeCut;
  static void SetResHistoMergeCut(Double_t d) {fgResHistoMergeCut=d;}

  // full reset: discard all statistics, zero histograms, start again.
  // called in online mode (HLT) after sending output for merging.
  virtual Bool_t            ResetOutputData();

protected:
  void ResetCurrent();                  // reset current values
  Int_t              fMemoryMode;       // 0 -do not fill THnSparse with residuals  1- fill only important QA THn 2 - Fill all THnsparse for calibration
  AliTPCcalibLaser * fLaser;            //! laser calibration
  //
  // current information
  //
  Float_t fDz;          //! current delta z
  
  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutMaxDz;    // maximal distance in z ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products
  Int_t   fCutTracks;   // maximal number of tracks
  Double_t   fMinPt;    //pt cut on tracks used for calibration
  Double_t   fMinPtITSTPCalign; //
 

  TH1F* fCosmiMatchingHisto[10];        // cosmic matching histogram
  //
  // distortion maps
  //
  THnSparse*  fResHistoTPCCE[5];        //TPC-TPCE matching map
  THnSparse*  fResHistoTPCITS[5];        //TPC-ITS    matching map
  THnSparse*  fResHistoTPCvertex[5];           //TPC-ITS    vertex matching map
  THnSparse*  fResHistoTPCTRD[5];        //TPC-TRD    matching map
  THnSparse*  fResHistoTPCTOF[5];        //TPC-TRD    matching map
  // laser histo
  THnSparse * fHistVdriftLaserA[3];	//Histograms for V drift from laser
  THnSparse * fHistVdriftLaserC[3];	//Histograms for V drift from laser
  TObjArray *fArrayLaserA;              //Object array of driftvelocity laserA
  TObjArray *fArrayLaserC;              //Object array of driftvelocity laserC
  //
  // TPC vertex A side C side histo
  //
  THnSparse * fTPCVertex[12];           // TPC vertex histograms A side c side - A+C -ESD
  THnSparse * fTPCVertexCorrelation[5];       // TPC vertex correlation A side C side with TPC vertex and ITS vertex     
  // DELTA Z histo
  TObjArray* fArrayDz;                  // array of DZ histograms for different triggers
  TObjArray* fAlignITSTPC;              // alignemnt array ITS TPC match
  TObjArray* fAlignTRDTPC;              // alignemnt array TRD TPC match
  TObjArray* fAlignTOFTPC;              // alignemnt array TOF TPC match
  Int_t      fTimeKalmanBin;            // width of Kalman bin - time in seconds
  Int_t    fTimeBins;			//Bins time
  Double_t fTimeStart;			//Start time
  Double_t fTimeEnd;			//End time
  Int_t    fPtBins;			//Bins pt
  Double_t fPtStart;			//Start pt
  Double_t fPtEnd;			//End pt
  Int_t    fVdriftBins;			//Bins vdrift
  Double_t fVdriftStart;		//Start vdrift
  Double_t fVdriftEnd;			//End vdrift
  Int_t    fRunBins;			//Bins run
  Double_t fRunStart;			//Start run
  Double_t fRunEnd;			//End run
  Int_t    fBinsVdrift[4];		//Bins for vdrift
  Double_t fXminVdrift[4];		//Xmax for vdrift
  Double_t fXmaxVdrift[4];		//Xmin for vdrift

private:
  AliTPCcalibTime(const AliTPCcalibTime&); 
  AliTPCcalibTime& operator=(const AliTPCcalibTime&); 

  ClassDef(AliTPCcalibTime, 10)
};

#endif


