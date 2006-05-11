//**************************************************************** 
//
//  Class for handling of ESD track quality cuts
//
//  TODO: 
//  - add functionality to save and load cuts
//  - fix the n sigma cut so it is really a n sigma cut
//

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TTree
#include "TTree.h"
#endif
#ifndef ROOT_TH2
#include "TH2.h"
#endif

#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliLog.h"

class ESDtrackQualityCuts : public TObject 
{
protected:

  //######################################################
  // esd track quality cuts
  static const Int_t fNCuts = 14;
  Char_t*            fCutNames[14];

  Int_t   fCut_MinNClusterTPC;        // min number of tpc clusters
  Int_t   fCut_MinNClusterITS;        // min number of its clusters  

  Float_t fCut_MaxChi2PerClusterTPC;  // max tpc fit chi2 per tpc cluster
  Float_t fCut_MaxChi2PerClusterITS;  // max its fit chi2 per its cluster

  Float_t fCut_MaxC11;                // max resolutions of covariance matrix diag. elements
  Float_t fCut_MaxC22;
  Float_t fCut_MaxC33;
  Float_t fCut_MaxC44;
  Float_t fCut_MaxC55;

  Float_t fCut_NsigmaToVertex;        // max number of estimated sigma from track-to-vertex
  Bool_t  fCut_SigmaToVertexRequired; // cut track if sigma from track-to-vertex could not be calculated
 
  Bool_t  fCut_AcceptKinkDaughters;   // accepting kink daughters?
  Bool_t  fCut_RequireTPCRefit;       // require TPC refit
  Bool_t  fCut_RequireITSRefit;       // require ITS refit
  
  //######################################################
  // diagnostics histograms
  Bool_t fHistogramsOn;

  TH1F** hNClustersITS;
  TH1F** hNClustersTPC;
  
  TH1F** hChi2PerClusterITS;
  TH1F** hChi2PerClusterTPC;

  TH1F** hC11;
  TH1F** hC22;
  TH1F** hC33;
  TH1F** hC44;
  TH1F** hC55;

  TH1F** hDXY;
  TH1F** hDZ;
  TH2F** hDXYvsDZ;

  TH1F** hDXYNormalized;
  TH1F** hDZNormalized;
  TH2F** hDXYvsDZNormalized;

  TH1F*  hCutStatistics;
  TH2F*  hCutCorrelation;
  

  // dummy array
  Int_t  fIdxInt[200];

public:
  ESDtrackQualityCuts();
  
  Bool_t AcceptTrack(AliESDtrack* esdTrack);
  Bool_t AcceptTrack(AliESDtrack* esdTrack, AliESDVertex* esdVtx, Double_t field);
  Bool_t AcceptTrack(AliESDtrack* esdTrack, Double_t* vtx, Double_t* vtx_res, Double_t field);
  Bool_t AcceptTrack(AliESDtrack* esdTrack, AliESDVertex* esdVtx, Float_t field)
    {return AcceptTrack(esdTrack,esdVtx, Double_t(field));}

  //######################################################
  // track quality cut setters  
  void   SetMinNClustersTPC(Int_t min=-1)          {fCut_MinNClusterTPC=min;}
  void   SetMinNClustersITS(Int_t min=-1)          {fCut_MinNClusterITS=min;}
  void   SetMaxChi2PerClusterTPC(Float_t max=1e99) {fCut_MaxChi2PerClusterTPC=max;}
  void   SetMaxChi2PerClusterITS(Float_t max=1e99) {fCut_MaxChi2PerClusterITS=max;}
  void   SetMinNsigmaToVertex(Float_t sigma=3)     {fCut_NsigmaToVertex=sigma;}
  void   SetRequireSigmaToVertex(Bool_t b=kTRUE)  { fCut_SigmaToVertexRequired = b; }
  void   SetRequireTPCRefit(Bool_t b=kFALSE)       {fCut_RequireTPCRefit=b;}
  void   SetRequireITSRefit(Bool_t b=kFALSE)       {fCut_RequireITSRefit=b;}
  void   SetAcceptKingDaughters(Bool_t b=kFALSE)   {fCut_AcceptKinkDaughters=b;}
  void   SetMaxCovDiagonalElements(Float_t c1=1e99, Float_t c2=1e99, Float_t c3=1e99, Float_t c4=1e99, Float_t c5=1e99) 
    {fCut_MaxC11=c1; fCut_MaxC22=c2; fCut_MaxC33=c3; fCut_MaxC44=c4; fCut_MaxC55=c5;}
 
  //######################################################
  void   SetHistogramsOn(Bool_t b=kFALSE) {fHistogramsOn = b;}
  void   DefineHistograms(Int_t color=1);
  void   SaveHistograms(Char_t* dir="track_selection");
  
  // void SaveQualityCuts(Char_t* file)
  // void LoadQualityCuts(Char_t* file)

  ClassDef(ESDtrackQualityCuts,0)
};

