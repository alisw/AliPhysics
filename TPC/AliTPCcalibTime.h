#ifndef ALITPCCALIBTIME_H
#define ALITPCCALIBTIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "THnSparse.h"           // Temporary
#include "TH1D.h"                // Temporary make code compiling for HLT in the 
class TObjArray;

class TH1F;
class TH3F;
class TH2F;
class THnSparse;
class TH1D;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliTPCcalibLaser;
class TGraphErrors;
class AliSplineFit;
class AliESDfriendTrack;

class AliTPCcalibTime:public AliTPCcalibBase {
public:
  AliTPCcalibTime(); 
  AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeVdrift);
  virtual ~AliTPCcalibTime();
  
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *const li);
  virtual void           Analyze();
  static Bool_t          IsLaser      (const AliESDEvent *const event);
  static Bool_t          IsCosmics    (const AliESDEvent *const event);
  static Bool_t          IsBeam       (const AliESDEvent *const event);
  void                   ProcessLaser (AliESDEvent *event);
  void                   ProcessCosmic(const AliESDEvent *const event);
  void                   ProcessBeam  (const AliESDEvent *const event);
  Bool_t                 IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1);
  Bool_t                 IsCross(AliESDtrack *const tr0, AliESDtrack *const tr1);
  Bool_t                 IsSame (AliESDtrack *const tr0, AliESDtrack *const tr1);
  void                   ProcessSame(AliESDtrack *const track, AliESDfriendTrack *const friendTrack, const AliESDEvent *const event);
  void                   ProcessAlignITS(AliESDtrack *const track, AliESDfriendTrack *const friendTrack, const AliESDEvent *const event, AliESDfriend *const ESDfriend);
  void                   ProcessAlignTRD(AliESDtrack* const track, AliESDfriendTrack *const friendTrack);
  void                   ProcessAlignTOF(AliESDtrack* const track, AliESDfriendTrack *const friendTrack);

  THnSparse*    GetHistVdriftLaserA(Int_t index=1) const {return fHistVdriftLaserA[index];};
  THnSparse*    GetHistVdriftLaserC(Int_t index=1) const {return fHistVdriftLaserC[index];};
  THnSparse*    GetHistoDrift(const char* name) const;
  TObjArray*    GetHistoDrift() const;
  TGraphErrors* GetGraphDrift(const char* name);
  TObjArray*    GetGraphDrift();
  AliSplineFit* GetFitDrift(const char* name);
//  TObjArray*    GetFitDrift();
  TH1F*         GetCosmiMatchingHisto(Int_t index=0) const {return fCosmiMatchingHisto[index];};
  
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  TObjArray* GetAlignITSTPC() const {return fAlignITSTPC;}              // alignemnt array ITS TPC match
  TObjArray* GetAlignTRDTPC() const {return fAlignTRDTPC;}              // alignemnt array TRD TPC match
  TObjArray* GetAlignTOFTPC() const {return fAlignTOFTPC;}              // alignemnt array TOF TPC match

  THnSparse*  GetResHistoTPCCE(Int_t index) const { return (index<5) ? fResHistoTPCCE[index]:0;}        //TPC-CE    matching map
  THnSparse*  GetResHistoTPCITS(Int_t index) const { return (index<5) ? fResHistoTPCITS[index]:0;}        //TPC-ITS    matching map
  THnSparse*  GetResHistoTPCvertex(Int_t index)      const { return (index<5) ? fResHistoTPCvertex[index]   :0;}        //TPC vertex matching map
  THnSparse*  GetResHistoTPCTRD(Int_t index)   const { return (index<5) ? fResHistoTPCTRD[index]:0;}        //TPC-TRD    matching map
  THnSparse*  GetResHistoTPCTOF(Int_t index)   const { return (index<5) ? fResHistoTPCTOF[index]:0;}        //TPC-TOF    matching map

  void        BookDistortionMaps();      // book histograms
  void        FillResHistoTPCCE(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pTPCOut );       // fill residual histo
  void        FillResHistoTPCITS(const AliExternalTrackParam * pTPCIn, const AliExternalTrackParam * pITSOut );       // fill residual histo
  void        FillResHistoTPC(const AliESDtrack * pTrack);
  void        FillResHistoTPCTRD(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTRDIn );
  void        FillResHistoTPCTOF(const AliExternalTrackParam * pTPCOut, const AliExternalTrackParam * pTOFIn );

private:
  void ResetCurrent();                  // reset current values

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
 
  AliTPCcalibTime(const AliTPCcalibTime&); 
  AliTPCcalibTime& operator=(const AliTPCcalibTime&); 

  TH1F* fCosmiMatchingHisto[10];
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
  // DELTA Z histo
  TObjArray* fArrayDz;                  // array of DZ histograms for different triggers
  TObjArray* fAlignITSTPC;              // alignemnt array ITS TPC match
  TObjArray* fAlignTRDTPC;              // alignemnt array TRD TPC match
  TObjArray* fAlignTOFTPC;              // alignemnt array TOF TPC match
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
  ClassDef(AliTPCcalibTime, 4); 
};

#endif


