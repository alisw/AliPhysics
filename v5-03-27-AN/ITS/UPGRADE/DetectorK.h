#ifndef DETECTORK_H
#define DETECTORK_H

#include <TNamed.h>
#include <TList.h>
#include <TGraph.h>
#include <Riostream.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

original code of using the billoir technique was developed
for the HFT (STAR), James H. Thomas, jhthomas@lbl.gov
http://rnc.lbl.gov/~jhthomas

Changes by S. Rossegger

July 2011  - Adding the possibility of "fake calculation" and "efficiency calculation" 
             with a number of "at least correct clusters on the track"
             Done using the complete combinatorics table with 3^nLayer track outcomes.

April 2011 - Now uses the Kalman method (aliroot implementation) instead of the Billoir
             technique ... (changes by Ruben Shahoyan)

March 2011 - Changes to comply with the Alice Offline coding conventions

Feb. 2011 - Improvement in "lowest pt allowed" -> now uses helix param. for calc. of a,b

          - Adding a more sophisticaed efficiency calculation which includes
            the possibility to make chi2 cuts via Confidence Levels (method of Ruben Shahoyan)
            plus adding 'basic' detection efficiencies per layer ...

          - Adding "ITS Stand alone" tracking capabilities via 
            forward+backward tracking -> Kalman smoothing is then 
            used for the parameter estimates (important for efficiencies)

Jan. 2011 - Inclusion of ImpactParameter Plots (with vtx estimates)
            to allow comparison with ITS performances in pp data

Dec. 2010 - Translation into C++ class format 
          - Adding various Setters and Getters to build the geometry 
	    (based on cylinders) plus handling of the layer properties 


***********************************************************/

class AliExternalTrackParam; 
#include <TMatrixD.h>

class DetectorK : public TNamed {

 public:
  
  DetectorK();
  DetectorK(char *name,char *title);
  virtual ~DetectorK();

  enum {kNptBins = 100}; // less then 400 !!
 
  void AddLayer(char *name, Float_t radius, Float_t radL, Float_t phiRes=999999, Float_t zRes=999999, Float_t eff=0.99);

  void KillLayer(char *name);
  void SetRadius(char *name, Float_t radius);
  void SetRadiationLength(char *name, Float_t radL);
  void SetResolution(char *name, Float_t phiRes=999999, Float_t zRes=999999);
  void SetLayerEfficiency(char *name, Float_t eff=1.0);
  void RemoveLayer(char *name);

  Float_t GetRadius(char *name);
  Float_t GetRadiationLength(char *name);
  Float_t GetResolution(char *name, Int_t axis=0);
  Float_t GetLayerEfficiency(char *name);

  void PrintLayout(); 
  void PlotLayout(Int_t plotDead = kTRUE);
  
  void MakeAliceAllNew(Bool_t flagTPC =1,Bool_t flagMon=1);
  void MakeAliceCurrent(Int_t AlignResiduals = 0, Bool_t flagTPC =1);
  void AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip=1);
  void RemoveTPC();

  void SetBField(Float_t bfield) {fBField = bfield; }
  Float_t GetBField() const {return fBField; }
  void SetLhcUPCscale(Float_t lhcUPCscale) {fLhcUPCscale = lhcUPCscale; }
  Float_t GetLhcUPCscale() const { return fLhcUPCscale; }
  void SetParticleMass(Float_t particleMass) {fParticleMass = particleMass; }
  Float_t GetParticleMass() const { return fParticleMass; }
  void SetIntegrationTime(Float_t integrationTime) {fIntegrationTime = integrationTime; }
  Float_t GetIntegrationTime() const { return fIntegrationTime; }
  void SetMaxRadiusOfSlowDetectors(Float_t maxRadiusSlowDet) {fMaxRadiusSlowDet =  maxRadiusSlowDet; }
  Float_t GetMaxRadiusOfSlowDetectors() const { return fMaxRadiusSlowDet; }
  void SetAvgRapidity(Float_t avgRapidity) {fAvgRapidity = avgRapidity; }
  Float_t GetAvgRapidity() const { return fAvgRapidity; }
  void SetConfidenceLevel(Float_t confLevel) {fConfLevel = confLevel; }
  Float_t GetConfidenceLevel() const { return fConfLevel; }
  void SetAtLeastCorr(Int_t atLeastCorr ) {fAtLeastCorr = atLeastCorr; }
  Int_t GetAtLeastCorr() const { return fAtLeastCorr; }
  void SetAtLeastFake(Int_t atLeastFake ) {fAtLeastFake = atLeastFake; }
  Int_t GetAtLeastFake() const { return fAtLeastFake; }

  void SetdNdEtaCent(Int_t dNdEtaCent ) {fdNdEtaCent = dNdEtaCent; }
  Float_t GetdNdEtaCent() const { return fdNdEtaCent; }

  
  Float_t GetNumberOfActiveLayers() const {return fNumberOfActiveLayers; }
  Float_t GetNumberOfActiveITSLayers() const {return fNumberOfActiveITSLayers; }

  void SolveViaBilloir(Int_t flagD0=1,Int_t print=1, Bool_t allPt=1, Double_t meanPt =0.750);
  void SolveDOFminusOneAverage();

  // Helper functions
  Double_t ThetaMCS                 ( Double_t mass, Double_t RadLength, Double_t momentum ) const;
  Double_t ProbGoodHit              ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbGoodChiSqHit         ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
  Double_t ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ )   ; 
 
  // Howard W. hit distribution and convolution integral
  Double_t Dist              ( Double_t Z, Double_t radius ) ;  
  Double_t HitDensity        ( Double_t radius )   ;
  Double_t UpcHitDensity     ( Double_t radius )   ;
  Double_t IntegratedHitDensity  ( Double_t multiplicity, Double_t radius )   ;
  Double_t OneEventHitDensity    ( Double_t multiplicity, Double_t radius ) const   ;
  Double_t D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][400] ) const ;
  
  TGraph* GetGraphMomentumResolution(Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolution(Int_t axis,Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth=1);

  TGraph* GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth=1);

  TGraph* GetGraphRecoEfficiency(Int_t particle, Int_t color, Int_t linewidth=1); 
  TGraph* GetGraphRecoFakes(Int_t particle,Int_t color, Int_t linewidth);
  TGraph* GetGraphRecoPurity(Int_t particle,Int_t color, Int_t linewidth);

  TGraph* GetGraph(Int_t number, Int_t color, Int_t linewidth=1);

  void MakeStandardPlots(Bool_t add =0, Int_t color=1, Int_t linewidth=1,Bool_t onlyPionEff=0);

  // method to extend AliExternalTrackParam functionality
  Bool_t GetXatLabR(AliExternalTrackParam* tr,Double_t r,Double_t &x, Double_t bz, Int_t dir=0) const;

  Double_t* PrepareEffFakeKombinations(TMatrixD *probKomb, TMatrixD *probLay);

 protected:
 
  Int_t fNumberOfLayers;        // total number of layers in the model
  Int_t fNumberOfActiveLayers;  // number of active layers in the model
  Int_t fNumberOfActiveITSLayers;  // number of active ITS layers in the model
  TList fLayers;                // List of layer pointers
  Float_t fBField;              // Magnetic Field in Tesla
  Float_t fLhcUPCscale;         // UltraPeripheralElectrons: scale from RHIC to LHC
  Float_t fIntegrationTime;     // electronics integration time
  Float_t fConfLevel;           // Confidence Level for the tracking
  Float_t fAvgRapidity;         // rapidity of the track (= mean)
  Float_t fParticleMass;        // Particle used for tracking. Standard: mass of pion
  Double_t fMaxRadiusSlowDet;   // Maximum radius for slow detectors.  Fast detectors 
                                // and only fast detectors reside outside this radius.
  Int_t fAtLeastCorr;     // min. number of correct hits for the track to be "good"
  Int_t fAtLeastFake;     // min. number of fake hits for the track to be "fake"

  Int_t fdNdEtaCent;       // Multiplicity

  enum {kMaxNumberOfDetectors = 200};
 
  Double_t fTransMomenta[kNptBins];                          // array of transverse momenta
  Double_t fMomentumRes[kNptBins];                           // array of momentum resolution
  Double_t fResolutionRPhi[kNptBins];                        // array of rphi resolution
  Double_t fResolutionZ[kNptBins];                           // array of z resolution
  Double_t fDetPointRes[kMaxNumberOfDetectors][kNptBins];    // array of rphi resolution per layer
  Double_t fDetPointZRes[kMaxNumberOfDetectors][kNptBins];   // array of z resolution per layer
  Double_t fEfficiency[3][kNptBins];                         // efficiency for different particles
  Double_t fFake[3][kNptBins];                               // fake prob for different particles

  ClassDef(DetectorK,1);
};

#endif
