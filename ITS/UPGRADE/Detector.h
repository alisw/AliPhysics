#ifndef DETECTOR_H
#define DETECTOR_H

#include <TNamed.h>
#include <TList.h>
#include <TGraph.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

original code of using the billoir technique was developed
for the HFT (STAR), James H. Thomas, jhthomas@lbl.gov
http://rnc.lbl.gov/~jhthomas

Changes by S. Rossegger
Dec. 2010 - Translation into C++ class format 
          - Adding various Setters and Getters to build the geometry 
	    (based on cylinders) plus handling of the layer properties 

Jan. 2011 - Inclusion of ImpactParameter Plots (with vtx estimates)
            to allow comparison with ITS performances in pp data

***********************************************************/


class Detector : public TNamed {
 public:
  Detector();
  Detector(char *name,char *title);
  virtual ~Detector();

  void AddLayer(char *name, Float_t radius, Float_t radL, Float_t phiRes=99999, Float_t zRes=99999);

  void KillLayer(char *name);
  void SetRadius(char *name, Float_t radius);
  void SetRadiationLength(char *name, Float_t radL);
  void SetResolution(char *name, Float_t phiRes=99999, Float_t zRes=99999);
  void RemoveLayer(char *name);

  Float_t GetRadius(char *name);
  Float_t GetRadiationLength(char *name);
  Float_t GetResolution(char *name, Int_t axis=0);

  void PrintLayout(); 

  void MakeAliceCurrent(Int_t AlignResiduals = 0, Bool_t flagTPC =1);
  void AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip=1);

  void SetBField(Float_t bfield) {fBField = bfield; }
  Float_t GetBField() {return fBField; }
  void SetLhcUPCscale(Float_t lhcUPCscale) {fLhcUPCscale = lhcUPCscale; }
  Float_t GetLhcUPCscale() { return fLhcUPCscale; }
  void SetParticleMass(Float_t particleMass) {fParticleMass = particleMass; }
  Float_t GetParticleMass() { return fParticleMass; }
  void SetIntegrationTime(Float_t integrationTime) {fIntegrationTime = integrationTime; }
  Float_t GetIntegrationTime() { return fIntegrationTime; }
  void SetAvgRapidity(Float_t avgRapidity) {fAvgRapidity = avgRapidity; }
  Float_t GetAvgRapidity() { return fAvgRapidity; }



  void SolveViaBilloir(Int_t flagD0=1,Int_t print=1);
  // Helper functions
  Double_t ThetaMCS          ( Double_t mass, Double_t RadLength, Double_t momentum )    ;
  Double_t ProbGoodHit       ( Double_t Radius, Double_t SearchRadiusRPhi, Double_t SearchRadiusZ )   ; 
  Double_t ProbGoodChiSqHit  ( Double_t Radius, Double_t SearchRadiusRPhi, Double_t SearchRadiusZ )   ; 
  // Howard W. hit distribution and convolution integral
  Double_t Dist              ( Double_t Z, Double_t Radius ) ;  
  Double_t HitDensity        ( Double_t Radius )   ;
  Double_t UpcHitDensity     ( Double_t Radius )   ;
  Double_t IntegratedHitDensity  ( Double_t Multiplicity, Double_t Radius )   ;
  Double_t OneEventHitDensity    ( Double_t Multiplicity, Double_t Radius )   ;
  Double_t D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][400] ) ;
  
  TGraph* GetGraphMomentumResolution(Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolution(Int_t axis,Int_t color, Int_t linewidth=1);
  TGraph* GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth=1);

  TGraph* GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth=1);

  TGraph* GetGraphRecoEfficiency(Int_t particle, Int_t color, Int_t linewidth=1); 
  
  TGraph* GetGraph(Int_t number, Int_t color, Int_t linewidth=1);
 

 protected:
 
  Int_t numberOfLayers;
  TList fLayers; 
  Float_t fBField;
  Float_t fLhcUPCscale;
  Float_t fIntegrationTime;
  Float_t fAvgRapidity;
  Float_t fParticleMass; // Particle used for tracking. Standard: mass of pion

  enum {MaxNumberOfDetectors = 300};
 
  Double_t xpoint[400], ypoint[400], yprime[400], yprimeZ[400], equivalent[400] ;
  Double_t DetPointRes[MaxNumberOfDetectors][400], DetPointZRes[MaxNumberOfDetectors][400] ;
  Double_t ratio[400], PXLReference[400], efficiency[3][400];

  ClassDef(Detector,1);
};

#endif
