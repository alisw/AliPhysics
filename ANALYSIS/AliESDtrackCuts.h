//
//  Class for handling of ESD track cuts.
//
//  The class manages a number of track quality cuts, a
//  track-to-vertex cut and a number of kinematic cuts. Two methods
//  can be used to figure out if an ESD track survives the cuts:
//  AcceptTrack which takes a single AliESDtrack as argument and
//  returns kTRUE/kFALSE or GetAcceptedTracks which takes an AliESDEvent
//  object and returns an TObjArray (of AliESDtracks) with the tracks
//  in the ESD that survived the cuts.
//
//
//  TODO:
//  - add functionality to save and load cuts
//  - add histograms for kinematic cut variables?
//  - upper and lower cuts for all (non-boolean) cuts
//  - update print method
//  - put comments to each variable
//

#ifndef ALIESDTRACKCUTS_H
#define ALIESDTRACKCUTS_H

#include <TF1.h>
#include <TH2.h>
#include "AliAnalysisCuts.h"

class AliESDEvent;
class AliESDtrack;
class AliLog;
class TTree;

class AliESDtrackCuts : public AliAnalysisCuts
{
public:
  enum ITSClusterRequirement { kOff = 0, kNone, kAny, kFirst, kOnlyFirst, kSecond, kOnlySecond, kBoth };
  enum Detector { kSPD = 0, kSDD, kSSD };
  
  AliESDtrackCuts(const Char_t* name = "AliESDtrackCuts", const Char_t* title = "");
  virtual ~AliESDtrackCuts();

  Bool_t IsSelected(TObject* obj)
       {return AcceptTrack((AliESDtrack*)obj);}
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  Bool_t AcceptTrack(AliESDtrack* esdTrack);
  TObjArray* GetAcceptedTracks(AliESDEvent* esd, Bool_t bTPC = kFALSE);
  Int_t CountAcceptedTracks(AliESDEvent* esd);

  static AliESDtrack* GetTPCOnlyTrack(AliESDEvent* esd, Int_t iTrack);

  virtual Long64_t Merge(TCollection* list);
  virtual void Copy(TObject &c) const;
  AliESDtrackCuts(const AliESDtrackCuts& pd);  // Copy Constructor
  AliESDtrackCuts &operator=(const AliESDtrackCuts &c);

  //######################################################
  // track quality cut setters  
  void SetMinNClustersTPC(Int_t min=-1)          {fCutMinNClusterTPC=min;}
  void SetMinNClustersITS(Int_t min=-1)          {fCutMinNClusterITS=min;}
  void SetClusterRequirementITS(Detector det, ITSClusterRequirement req = kOff) { fCutClusterRequirementITS[det] = req; }
  void SetMaxChi2PerClusterTPC(Float_t max=1e10) {fCutMaxChi2PerClusterTPC=max;}
  void SetMaxChi2PerClusterITS(Float_t max=1e10) {fCutMaxChi2PerClusterITS=max;}
  void SetRequireTPCRefit(Bool_t b=kFALSE)       {fCutRequireTPCRefit=b;}
  void SetRequireITSRefit(Bool_t b=kFALSE)       {fCutRequireITSRefit=b;}
  void SetAcceptKingDaughters(Bool_t b=kFALSE)   {fCutAcceptKinkDaughters=b;}
  void SetMaxCovDiagonalElements(Float_t c1=1e10, Float_t c2=1e10, Float_t c3=1e10, Float_t c4=1e10, Float_t c5=1e10) 
    {fCutMaxC11=c1; fCutMaxC22=c2; fCutMaxC33=c3; fCutMaxC44=c4; fCutMaxC55=c5;}

  // track to vertex cut setters
  void SetMaxNsigmaToVertex(Float_t sigma=1e10)       {fCutNsigmaToVertex = sigma; SetRequireSigmaToVertex(kTRUE);}
  void SetRequireSigmaToVertex(Bool_t b=kTRUE )       {fCutSigmaToVertexRequired = b;}
  void SetMaxDCAToVertexXY(Float_t dist=1e10)         {fCutDCAToVertexXY = dist;}
  void SetMaxDCAToVertexZ(Float_t dist=1e10)          {fCutDCAToVertexZ = dist;}
  void SetDCAToVertex2D(Bool_t b=kFALSE)              {fCutDCAToVertex2D = b;}
  
  // deprecated, will be removed in next release
  void SetMaxDCAToVertex(Float_t dist=1e10);
  void SetMinNsigmaToVertex(Float_t sigma=1e10);
  void SetDCAToVertex(Float_t dist=1e10);
  void SetDCAToVertexXY(Float_t dist=1e10);
  void SetDCAToVertexZ(Float_t dist=1e10);
  Float_t GetMinNsigmaToVertex() const;

  // getters

  Int_t   GetMinNClusterTPC()        const   { return fCutMinNClusterTPC;}
  Int_t   GetMinNClustersITS()       const   { return fCutMinNClusterITS;}
  ITSClusterRequirement GetClusterRequirementITS(Detector det) const { return fCutClusterRequirementITS[det]; }
  Float_t GetMaxChi2PerClusterTPC()  const   { return fCutMaxChi2PerClusterTPC;}
  Float_t GetMaxChi2PerClusterITS()  const   { return fCutMaxChi2PerClusterITS;}
  Bool_t  GetRequireTPCRefit()       const   { return fCutRequireTPCRefit;}
  Bool_t  GetRequireITSRefit()       const   { return fCutRequireITSRefit;}
  Bool_t  GetAcceptKingDaughters()   const   { return fCutAcceptKinkDaughters;}
  void    GetMaxCovDiagonalElements(Float_t& c1, Float_t& c2, Float_t& c3, Float_t& c4, Float_t& c5) 
      {c1 = fCutMaxC11; c2 = fCutMaxC22; c3 = fCutMaxC33; c4 = fCutMaxC44; c5 = fCutMaxC55;}
  Float_t GetMaxNsigmaToVertex()     const   { return fCutNsigmaToVertex;}
  Float_t GetMaxDCAToVertexXY()       const   { return fCutDCAToVertexXY;}
  Float_t GetMaxDCAToVertexZ()       const   { return fCutDCAToVertexZ;}
  Bool_t  GetDCAToVertex2D()         const   { return fCutDCAToVertex2D;}
  Bool_t  GetRequireSigmaToVertex( ) const   { return fCutSigmaToVertexRequired;}

  void GetPRange(Float_t& r1, Float_t& r2)   {r1=fPMin;   r2=fPMax;}
  void GetPtRange(Float_t& r1, Float_t& r2)  {r1=fPtMin;  r2=fPtMax;}
  void GetPxRange(Float_t& r1, Float_t& r2)  {r1=fPxMin;  r2=fPxMax;}
  void GetPyRange(Float_t& r1, Float_t& r2)  {r1=fPyMin;  r2=fPyMax;}
  void GetPzRange(Float_t& r1, Float_t& r2)  {r1=fPzMin;  r2=fPzMax;}
  void GetEtaRange(Float_t& r1, Float_t& r2) {r1=fEtaMin; r2=fEtaMax;}
  void GetRapRange(Float_t& r1, Float_t& r2) {r1=fRapMin; r2=fRapMax;}

  // track kinmatic cut setters
  void SetPRange(Float_t r1=0, Float_t r2=1e10)       {fPMin=r1;   fPMax=r2;}
  void SetPtRange(Float_t r1=0, Float_t r2=1e10)      {fPtMin=r1;  fPtMax=r2;}
  void SetPxRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPxMin=r1;  fPxMax=r2;}
  void SetPyRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPyMin=r1;  fPyMax=r2;}
  void SetPzRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPzMin=r1;  fPzMax=r2;}
  void SetEtaRange(Float_t r1=-1e10, Float_t r2=1e10) {fEtaMin=r1; fEtaMax=r2;}
  void SetRapRange(Float_t r1=-1e10, Float_t r2=1e10) {fRapMin=r1; fRapMax=r2;}

  //######################################################
  void SetHistogramsOn(Bool_t b=kFALSE) {fHistogramsOn = b;}
  void DefineHistograms(Int_t color=1);
  virtual Bool_t LoadHistograms(const Char_t* dir = 0);
  void SaveHistograms(const Char_t* dir = 0);
  void DrawHistograms();

  static Float_t GetSigmaToVertex(AliESDtrack* esdTrack);
  
  static void EnableNeededBranches(TTree* tree);

  // void SaveQualityCuts(Char_t* file)
  // void LoadQualityCuts(Char_t* file)

	TH1* GetDZNormalized(Int_t i) const { return fhDZNormalized[i]; }

protected:
  void Init(); // sets everything to 0
  Bool_t CheckITSClusterRequirement(ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2);
  
  enum { kNCuts = 27 };

  //######################################################
  // esd track quality cuts
  static const Char_t* fgkCutNames[kNCuts]; //! names of cuts (for internal use)

  Int_t   fCutMinNClusterTPC;         // min number of tpc clusters
  Int_t   fCutMinNClusterITS;         // min number of its clusters
  
  ITSClusterRequirement fCutClusterRequirementITS[3];  // detailed ITS cluster requirements for (SPD, SDD, SSD)

  Float_t fCutMaxChi2PerClusterTPC;   // max tpc fit chi2 per tpc cluster
  Float_t fCutMaxChi2PerClusterITS;   // max its fit chi2 per its cluster

  Float_t fCutMaxC11;                 // max cov. matrix diag. elements (res. y^2)
  Float_t fCutMaxC22;                 // max cov. matrix diag. elements (res. z^2)
  Float_t fCutMaxC33;                 // max cov. matrix diag. elements (res. sin(phi)^2)
  Float_t fCutMaxC44;                 // max cov. matrix diag. elements (res. tan(theta_dip)^2)
  Float_t fCutMaxC55;                 // max cov. matrix diag. elements (res. 1/pt^2)

  Bool_t  fCutAcceptKinkDaughters;    // accepting kink daughters?
  Bool_t  fCutRequireTPCRefit;        // require TPC refit
  Bool_t  fCutRequireITSRefit;        // require ITS refit

  // track to vertex cut
  Float_t fCutNsigmaToVertex;         // max number of estimated sigma from track-to-vertex
  Bool_t  fCutSigmaToVertexRequired;  // cut track if sigma from track-to-vertex could not be calculated
  Float_t fCutDCAToVertexXY;          // track-to-vertex cut in absolute distance in xy-plane
  Float_t fCutDCAToVertexZ;           // track-to-vertex cut in absolute distance in z-plane
  Bool_t  fCutDCAToVertex2D;          // if true a 2D DCA cut using fCutDCAToVertexXY and fCutDCAToVertexZ is made. Tracks are accepted if sqrt((DCAXY / fCutDCAToVertexXY)^2 + (DCAZ / fCutDCAToVertexZ)^2) < 1

  // esd kinematics cuts
  Float_t fPMin,   fPMax;             // definition of the range of the P
  Float_t fPtMin,  fPtMax;            // definition of the range of the Pt
  Float_t fPxMin,  fPxMax;            // definition of the range of the Px
  Float_t fPyMin,  fPyMax;            // definition of the range of the Py
  Float_t fPzMin,  fPzMax;            // definition of the range of the Pz
  Float_t fEtaMin, fEtaMax;           // definition of the range of the eta
  Float_t fRapMin, fRapMax;           // definition of the range of the y

  //######################################################
  // diagnostics histograms
  Bool_t fHistogramsOn;               // histograms on/off

  TH1F* fhNClustersITS[2];            //->
  TH1F* fhNClustersTPC[2];            //->

  TH1F* fhChi2PerClusterITS[2];       //->
  TH1F* fhChi2PerClusterTPC[2];       //->

  TH1F* fhC11[2];                     //->
  TH1F* fhC22[2];                     //->
  TH1F* fhC33[2];                     //->
  TH1F* fhC44[2];                     //->
  TH1F* fhC55[2];                     //->

  TH1F* fhDXY[2];                     //->
  TH1F* fhDZ[2];                      //->
  TH1F* fhDXYDZ[2];                   //-> absolute distance sqrt(dxy**2 + dz**2) to vertex; if 2D cut is set, normalized to given values
  TH2F* fhDXYvsDZ[2];                 //->

  TH1F* fhDXYNormalized[2];           //->
  TH1F* fhDZNormalized[2];            //->
  TH2F* fhDXYvsDZNormalized[2];       //->
  TH1F* fhNSigmaToVertex[2];          //->

  TH1F* fhPt[2];                      //-> pt of esd tracks
  TH1F* fhEta[2];                     //-> eta of esd tracks

  TF1*  ffDTheoretical;               //-> theoretical distance to vertex normalized (2d gauss)

  TH1F*  fhCutStatistics;             //-> statistics of what cuts the tracks did not survive
  TH2F*  fhCutCorrelation;            //-> 2d statistics plot

  ClassDef(AliESDtrackCuts, 4)
};


#endif
