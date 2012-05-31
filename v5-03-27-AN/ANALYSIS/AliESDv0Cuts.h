//
//  Class for handling of ESD v0 cuts.
//
//

#ifndef ALIESDV0CUTS_H
#define ALIESDV0CUTS_H

#include <TF1.h>
#include <TH2.h>
#include "AliAnalysisCuts.h"

class AliESD;
class AliESDEvent;
class AliESDVertex;
class AliESDtrack;
class AliESDv0;
class AliLog;
class TTree;

class AliESDv0Cuts : public AliAnalysisCuts
{
public:
  AliESDv0Cuts(const Char_t* name = "AliESDv0Cuts", const Char_t* title = "");
  virtual ~AliESDv0Cuts();

  Bool_t IsSelected(TObject* /*obj*/) {return kTRUE;}
  Bool_t IsSelected(TList* listObj);
  Bool_t IsSelected(TObject* const obj1, TObject* const obj2, TObject* const obj3, TObject* const obj4)
  {return AcceptV0((AliESDv0*) obj1, (AliESDtrack*) obj2, (AliESDtrack*) obj3, (const AliESDVertex*) obj4);}
  Bool_t AcceptV0(AliESDv0* const esdV0, AliESDtrack* const trackPos, AliESDtrack* const trackNeg, const AliESDVertex*  esdVertex);
  TObjArray* GetAcceptedV0s(const AliESD* esd);
  Int_t CountAcceptedV0s(const AliESD* esd);
  TObjArray* GetAcceptedV0s(const AliESDEvent* esd);
  Int_t CountAcceptedV0s(const AliESDEvent* esd);

  virtual Long64_t Merge(TCollection* list);
  virtual void Copy(TObject &c) const;
  AliESDv0Cuts(const AliESDv0Cuts& pd);  // Copy Constructor
  AliESDv0Cuts &operator=(const AliESDv0Cuts &c);

  //######################################################
  // v0 quality cut setters  
  void SetMinDcaPosToVertex(Float_t min=-1)          {fCutMinDcaPosToVertex=min;}
  void SetMinDcaNegToVertex(Float_t min=-1)          {fCutMinDcaNegToVertex=min;}
  void SetMaxChi2(Float_t max=1e10)                  {fCutMaxChi2=max;}
  void SetMaxDcaV0Daughters(Float_t max=1e10)        {fCutMaxDcaV0Daughters=max;}
  void SetMinRadius(Float_t min=-1)                  {fCutMinRadius=min;}
  void SetMaxRadius(Float_t max=1e10)                {fCutMaxRadius=max;}
  void SetMinCosinePointingAngle(Float_t min=-1)     {fCutMinCosinePointingAngle=min;}
  void SetRequireOnFlyStatus(Bool_t b=kFALSE)        {fCutRequireOnFlyStatus=b;}
  void SetMaxDcaV0ToVertex(Float_t max=1e10)         {fCutMaxDcaV0ToVertex=max;}

  // v0 kinematic cut setters
  void SetPRange(Float_t r1=0, Float_t r2=1e10)      {fPMin=r1;   fPMax=r2;}
  void SetPtRange(Float_t r1=0, Float_t r2=1e10)     {fPtMin=r1;  fPtMax=r2;}
  void SetPxRange(Float_t r1=-1e10, Float_t r2=1e10) {fPxMin=r1;  fPxMax=r2;}
  void SetPyRange(Float_t r1=-1e10, Float_t r2=1e10) {fPyMin=r1;  fPyMax=r2;}
  void SetPzRange(Float_t r1=-1e10, Float_t r2=1e10) {fPzMin=r1;  fPzMax=r2;}

  //######################################################
  void SetHistogramsOn(Bool_t b=kFALSE) {fHistogramsOn = b;}
  void DefineHistograms(Int_t color=1);
  virtual Bool_t LoadHistograms(const Char_t* dir = 0);
  void SaveHistograms(const Char_t* dir = 0);
  void DrawHistograms();

  static void EnableNeededBranches(TTree* tree);

  // void SaveQualityCuts(Char_t* file)
  // void LoadQualityCuts(Char_t* file)

protected:
  void Init(); // sets everything to 0

  enum { kNCuts = 14 };

  //######################################################
  // esd v0 quality cuts
  static const Char_t* fgkCutNames[kNCuts]; //! names of cuts (for internal use)

  Float_t fCutMinDcaPosToVertex;      // min dca of the positive daughter to the primary vertex
  Float_t fCutMinDcaNegToVertex;      // min dca of the negative daughter to the primary vertex
  Float_t fCutMaxChi2;                // max chi2
  Float_t fCutMaxDcaV0Daughters;      // max dca between the two v0 daughters
  Float_t fCutMinRadius;              // min reconstruction radius (fiducial volume)
  Float_t fCutMaxRadius;              // max reconstruction radius (fiducial volume)
  Float_t fCutMinCosinePointingAngle; // min cosine of pointing angle
  Bool_t  fCutRequireOnFlyStatus;     // require on fly status
  Float_t fCutMaxDcaV0ToVertex;       // max dca of the v0 to the primary vertex

  // v0 kinematics cuts
  Float_t fPMin,   fPMax;             // definition of the range of the P
  Float_t fPtMin,  fPtMax;            // definition of the range of the Pt
  Float_t fPxMin,  fPxMax;            // definition of the range of the Px
  Float_t fPyMin,  fPyMax;            // definition of the range of the Py
  Float_t fPzMin,  fPzMax;            // definition of the range of the Pz

  //######################################################
  // diagnostics histograms
  Bool_t fHistogramsOn;               // histograms on/off

  TH1F* fhDcaPosToVertex[2];          //->
  TH1F* fhDcaNegToVertex[2];          //->
  TH1F* fhChi2[2];                    //->
  TH1F* fhDcaV0Daughters[2];          //->
  TH1F* fhRadius[2];                  //->
  TH1F* fhCosinePointingAngle[2];     //->
  TH1F* fhOnFlyStatus[2];             //->
  TH1F* fhDcaV0ToVertex[2];           //->
  
  TH1F* fhPt[2];                      //-> pt of esd v0s

  TH1F* fhCutStatistics;              //-> statistics of what cuts the v0s did not survive
  TH2F* fhCutCorrelation;             //-> 2d statistics plot

  ClassDef(AliESDv0Cuts, 1)
};


#endif
