#ifndef ALICONVERSIONTRACKCUTS_H
#define ALICONVERSIONTRACKCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: (this code is mostly copied from AliRsnTrackQuality) adapted by Svein Lindal 	*

class TH2F;
class AliESDtrackCuts;
class TList;
class AliVEvent;
class THn;
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisCuts.h"

using namespace std;

class AliConversionTrackCuts : public AliAnalysisCuts {
	
public:
  
  enum CTCuts_t {
    kPreCut = 0,
    kCutNcls,
    kCutNclsFrac,
    kCutNDF,
    kCutKinc,
    kCutDCAZ,
    kCutDCAXY,
    kCutTPCRefit,
    kNCuts
  };

  //  enum trackVal {
  //	kPt = 0, 
  //	kPhi,
  //	kEta,
  //	kNVar
  //  };
  static const char * fgkCutNames[kNCuts];

  Bool_t IsSelected(TObject * object ); 
  Bool_t IsSelected(TList * /*list*/) { return kFALSE; }
  Bool_t AcceptTrack(const AliAODTrack * track);
  Bool_t AcceptTrack(const AliESDtrack * track);
  Bool_t GetDCA(const AliAODTrack * track, Double_t dca[2]);
  Bool_t GetDCA(const AliESDtrack * track, Double_t dca[2]);


  void FillDCAHist(Float_t dcaz, Float_t dcaxy, const AliVTrack * track);
  AliConversionTrackCuts();
  AliConversionTrackCuts(TString name, TString title);
  ~AliConversionTrackCuts();
  
  void      SetEsdTrackCuts(AliESDtrackCuts * trackcuts) { fEsdTrackCuts = trackcuts; }
  void      SetDCAZmax(Double_t value)  { fDCAZmax = value*value; }
  void      SetDCAXYmax(Double_t value) { fDCAXYmax = value*value; }
  void      SetFilterBit(Int_t value)   { fFilterBit = value; }
  void      SetEvent(AliVEvent * event)  { fEvent = event; }
  void      CreateTrackEff(Bool_t create = kTRUE) { fkCreateTrackEff = create; }

  TList * CreateHistograms();
  void FillHistograms(Int_t cutIndex, const AliVTrack * track);
  virtual void   Print(const Option_t *option = "") const;

protected :

  void DefineESDCuts();


  AliESDtrackCuts * fEsdTrackCuts; //main cut
  AliESDtrackCuts * fEsdTrackCutsExtra1; //global tracks cut
  AliESDtrackCuts * fEsdTrackCutsExtra2; //others
  
  AliVEvent * fEvent;



  Int_t fFilterBit;
  
  Double_t fDCAZmax;                // maximum value for longitudinal DCA
  Double_t fDCAXYmax;                  // maximum xy value for dca
  
  Bool_t fInitialized;

  TH2F * fhPhi;
  //TH2F * fhPt;
  //TH2F * fhPhiPt;
  TH2F * fhdcaxyPt;
  TH2F * fhdcazPt;
  TH2F * fhdca;
  TH2F * fhnclpt;
  TH2F * fhnclsfpt;
  TH2F * fhEtaPhi;
  THn  * fhTrackEff;
  Bool_t fkCreateTrackEff;

  TList * fHistograms;

  AliConversionTrackCuts(const AliConversionTrackCuts&); // not implemented
  AliConversionTrackCuts& operator=(const AliConversionTrackCuts&); // not implemented

  ClassDef(AliConversionTrackCuts,5)
};

#endif
