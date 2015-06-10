#include "DetectorK.h"
#include "AliLog.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFormula.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TText.h>
#include <TGraphErrors.h>

#include "AliExternalTrackParam.h"

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

original code of using the billoir technique was developed
for the HFT (STAR), James H. Thomas, jhthomas@lbl.gov
http://rnc.lbl.gov/~jhthomas

Changes by S. Rossegger -> see header file

***********************************************************/
Bool_t DetectorK::verboseR=0;

#define RIDICULOUS 999999 // A ridiculously large resolution (cm) to flag a dead detector

#define Luminosity    1.e27       // Luminosity of the beam (LHC HI == 1.e27, RHIC II == 8.e27 )
#define SigmaD        6.0         // Size of the interaction diamond (cm) (LHC = 6.0 cm)
#define dNdEtaMinB    1//950//660//950           // Multiplicity per unit Eta  (AuAu MinBias = 170, Central = 700)
// #define dNdEtaCent    2300//15000 //1600//2300        // Multiplicity per unit Eta  (LHC at 5.5 TeV not known)

#define CrossSectionMinB         8    // minB Cross section for event under study (PbPb MinBias ~ 8 Barns)
#define AcceptanceOfTpcAndSi     1 //1//0.60 //0.35  // Assumed geometric acceptance (efficiency) of the TPC and Si detectors
#define UPCBackgroundMultiplier  1.0   // Increase multiplicity in detector (0.0 to 1.0 * UPCRate ) (eg 1.0)
#define OtherBackground          0.0   // Increase multiplicity in detector (0.0 to 1.0 * minBias)  (eg 0.0)
#define EfficiencySearchFlag     2     // Define search method:
                                       // -> ChiSquarePlusConfLevel = 2, ChiSquare = 1, Simple = 0.  

#define PionMass                 0.139  // Mass of the Pion
#define KaonMass                 0.498  // Mass of the Kaon
#define D0Mass                   1.865  // Mass of the D0

ClassImp(TrackSol)

const double DetectorK::kPtMinFix = 0.050;
const double DetectorK::kPtMaxFix = 31.5;

//TMatrixD *probKomb; // table for efficiency kombinatorics

class ForwardLayer : public TNamed {
public:
  ForwardLayer(char *name) : TNamed(name,name) {}
  
  Float_t GetZ()         const {return zPos;}
  Float_t GetXRes()      const {return xRes;}
  Float_t GetYRes()      const {return yRes;}
  Float_t GetThickness() const {return thickness;}
  Float_t Getdensity()   const {return density;}
  Float_t GetLayerEff()  const {return eff;}

  //  void Print() {printf("  r=%3.1lf X0=%1.6lf sigPhi=%1.4lf sigZ=%1.4lf\n",radius,radL,phiRes,zRes); }
  Float_t zPos; Float_t xRes; Float_t yRes;   
  Float_t radL;
  Float_t thickness;
  Float_t density;
  Float_t eff;
  Bool_t isDead;

 ClassDef(ForwardLayer,1);
};


ClassImp(DetectorK)
DetectorK::DetectorK() 
  : TNamed("test_detector","detector"),
    fNumberOfLayers(0),
    fNumberOfActiveLayers(0),
    fNumberOfActiveITSLayers(0),
    fBField(0.5),
    fLhcUPCscale(1.0),    
    fIntegrationTime(0.02), // in ms
    fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140),    // Standard: pion mass 
    fMaxRadiusSlowDet(10.),
    fAtLeastHits(-1),     // if -1, then require hit on all ITS layers
    fAtLeastCorr(-1),     // if -1, then correct hit on all ITS layers
    fAtLeastFake(1),       // if at least x fakes, track is considered fake ...
    fMaxSeedRadius(50000),
    fptScale(10.),
  fdNdEtaCent(2300),
  kDetLayer(-1),
  fMinRadTrack(132.)
{
  //
  // default constructor
  //
  //  fLayers = new TObjArray();
  
}

DetectorK::DetectorK(char *name, char *title)
  : TNamed(name,title),
    fNumberOfLayers(0),
    fNumberOfActiveLayers(0),
    fNumberOfActiveITSLayers(0),
    fBField(0.5),
    fLhcUPCscale(1.0),
    fIntegrationTime(0.02),  // in ms
    fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140),     // Standard: pion mass
    fMaxRadiusSlowDet(10.),
    fAtLeastHits(-1),     // if -1, then require hit on all ITS layers
    fAtLeastCorr(-1),     // if -1, then correct hit on all ITS layers
    fAtLeastFake(1),       // if at least x fakes, track is considered fake ...
    fMaxSeedRadius(50000),
    fptScale(10.),
    fdNdEtaCent(2200),
  kDetLayer(-1),
  fMinRadTrack(132.)
{
  //
  // default constructor, that set the name and title
  //
  //  fLayers = new TObjArray();
}
DetectorK::~DetectorK() { // 
  // virtual destructor
  //
  //  delete fLayers;
}

void DetectorK::AddLayer(char *name, Float_t radius, Float_t radL, Float_t phiRes, Float_t zRes, Float_t eff) {
  //
  // Add additional layer to the list of layers (ordered by radius)
  // 

  CylLayerK *newLayer = (CylLayerK*) fLayers.FindObject(name);

  if (!newLayer) {
    newLayer = new CylLayerK(name);
    newLayer->radius = radius;
    newLayer->radL = radL;
    newLayer->phiRes = phiRes;
    newLayer->zRes = zRes;
    newLayer->eff = eff;

    if (newLayer->zRes==RIDICULOUS && newLayer->zRes==RIDICULOUS) 
      newLayer->isDead = kTRUE;
    else 
      newLayer->isDead = kFALSE;
  
    if (fLayers.GetEntries()==0) 
      fLayers.Add(newLayer);
    else {
      
      for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
	CylLayerK *l = (CylLayerK*)fLayers.At(i);
	if (radius<l->radius) {
	  fLayers.AddBefore(l,newLayer);
	  break;
	}
	if (radius>l->radius && (i+1)==fLayers.GetEntries() ) { 
	  // even bigger then last one
	  fLayers.Add(newLayer);
	}
      }
      
    }
    fNumberOfLayers += 1;
    if (!(newLayer->isDead)) {
      fNumberOfActiveLayers += 1;
      TString lname(newLayer->GetName());
      if ( IsITSLayer(lname) ) fNumberOfActiveITSLayers += 1;
    }


  } else {
    printf("Layer with the name %s does already exist\n",name);
  }
  

}

void DetectorK::KillLayer(char *name) {
  //
  // Marks layer as dead. Contribution only by Material Budget
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot mark as dead\n",name);
  else {
     tmp->phiRes = 999999;
     tmp->zRes = 999999;
     if (!(tmp->isDead)) {
       tmp->isDead = kTRUE;
       fNumberOfActiveLayers -= 1; 
       TString lname(tmp->GetName());
       if ( IsITSLayer(lname) ) fNumberOfActiveITSLayers -= 1;
     }     
  }
}

void DetectorK::SetRadius(char *name, Float_t radius) {
  //
  // Set layer radius [cm]
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
 

  if (!tmp) {
    printf("Layer %s not found - cannot set radius\n",name);
  } else {
      
    Float_t tmpRadL  = tmp->radL;
    Float_t tmpPhiRes = tmp->phiRes;
    Float_t tmpZRes = tmp->zRes;

    RemoveLayer(name); // so that the ordering is correct
    AddLayer(name,radius,tmpRadL,tmpPhiRes,tmpZRes);
  }
}

Float_t DetectorK::GetRadius(char *name) {
  //
  // Return layer radius [cm]
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get radius\n",name);
  else 
    return tmp->radius;

  return 0;
}

void DetectorK::SetRadiationLength(char *name, Float_t radL) {
  //
  // Set layer material [cm]
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer material\n",name);
  else {
    tmp->radL = radL;
  }
}

Float_t DetectorK::GetRadiationLength(char *name) {
  //
  // Return layer radius [cm]
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer material\n",name);
  else 
    return tmp->radL;
    
  return 0;
  
}

void DetectorK::SetResolution(char *name, Float_t phiRes, Float_t zRes) {
  //
  // Set layer resolution in [cm]
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set resolution\n",name);
  else {

    Bool_t wasDead = tmp->isDead;
    
    tmp->phiRes = phiRes;
    tmp->zRes = zRes;
    TString lname(tmp->GetName());

    if (zRes==RIDICULOUS && phiRes==RIDICULOUS) {
      tmp->isDead = kTRUE;
      if (!wasDead) {
	fNumberOfActiveLayers -= 1;
	if ( IsITSLayer(lname) ) fNumberOfActiveITSLayers -= 1;
      }
    } else {
      tmp->isDead = kFALSE;
      if (wasDead) {
	fNumberOfActiveLayers += 1;
	if ( IsITSLayer(lname) ) fNumberOfActiveITSLayers += 1;
      }
    }


  }
}

Float_t DetectorK::GetResolution(char *name, Int_t axis) {
  //
  // Return layer resolution in [cm]
  // axis = 0: resolution in rphi
  // axis = 1: resolution in z
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get resolution\n",name);
  else {
    if (axis==0) return tmp->phiRes;
    if (axis==1) return tmp->zRes;
    printf("error: axis must be either 0 or 1 (rphi or z axis)\n");
  }
  return 0;
}

void DetectorK::SetLayerEfficiency(char *name, Float_t eff) {
  //
  // Set layer efficnecy (prop that his is missed within this layer) 
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer efficiency\n",name);
  else {
    tmp->eff = eff;
  }
}

Float_t DetectorK::GetLayerEfficiency(char *name) {
  //
  // Get layer efficnecy (prop that his is missed within this layer) 
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer efficneicy\n",name);
  else 
    return tmp->eff;
    
  return 0;
  
}

void DetectorK::RemoveLayer(char *name) {
  //
  // Removes a layer from the list
  //

  CylLayerK *tmp = (CylLayerK*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot remove it\n",name);
  else {
    Bool_t wasDead = tmp->isDead;
    fLayers.Remove(tmp);
    fNumberOfLayers -= 1;
    if (!wasDead) {
      fNumberOfActiveLayers -= 1;
      TString lname(tmp->GetName());
      if ( IsITSLayer(lname) ) fNumberOfActiveITSLayers -= 1;
      
    }
  }
}


CylLayerK* DetectorK::FindLayer(char *name) const
{
  //
  // find layer by name
  //
  return (CylLayerK*) fLayers.FindObject(name);
}

CylLayerK* DetectorK::FindLayer(double r, int mode) const
{
  //
  // find layer close to radius r
  // mode = 0: closest
  // mode > 0: closest above
  // mode < 0: closest below
  //
  double drMin=-9999;
  int lrID = -1;
  int nLr = fLayers.GetEntries();
  for (Int_t i=fLayers.GetEntries(); i--;) {
    CylLayerK* tmp = (CylLayerK*)fLayers.At(i);
    double dr = tmp->radius - r;
    if (TMath::Abs(dr)<TMath::Abs(drMin)) {
      drMin = dr;
      lrID = i;
    }
  }
  if (lrID<0) return 0;  
  if (mode>0 && drMin<0) return ++lrID<nLr ? (CylLayerK*)fLayers.At(lrID) : 0;
  if (mode<0 && drMin>0) return --lrID>0   ? (CylLayerK*)fLayers.At(lrID) : 0;
  return (CylLayerK*)fLayers.At(lrID);
  //
}

Int_t DetectorK::FindLayerID(double r, int mode) const
{
  //
  // find layer ID close to radius r
  // mode = 0: closest
  // mode > 0: closest above
  // mode < 0: closest below
  //
  double drMin=-9999;
  int lrID = -1;
  int nLr = fLayers.GetEntries();
  for (Int_t i=fLayers.GetEntries(); i--;) {
    CylLayerK* tmp = (CylLayerK*)fLayers.At(i);
    double dr = tmp->radius - r;
    if (TMath::Abs(dr)<TMath::Abs(drMin)) {
      drMin = dr;
      lrID = i;
    }
  }
  if (lrID<0) return 0;  
  if (mode>0 && drMin<0) return ++lrID<nLr ? lrID : -1;
  if (mode<0 && drMin>0) return --lrID>0   ? lrID : -1;
  return lrID;
  //
}


void DetectorK::PrintLayout(Bool_t full) {
  //
  // Prints the detector layout
  //

  printf("Detector %s: \"%s\"\n",GetName(),GetTitle());
  
  if (fLayers.GetEntries()>0) 
    printf("  Name \t\t r [cm] \t  X0 \t  phi & z res [um] layerEff \n");

  CylLayerK *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (CylLayerK*)fLayers.At(i);
  
    // don't print all the tpc layers
    TString name(tmp->GetName());
    if (!full && !IsITSLayer(name) && !name.Contains("_0")) continue;

    printf("%d. %s \t %03.2f   \t%1.4f\t  ",i,
 	   tmp->GetName(), tmp->radius, tmp->radL);
    if (tmp->phiRes==RIDICULOUS) 
      printf("  -  ");
    else
      printf("%3.0f   ",tmp->phiRes*10000);
    if (tmp->zRes==RIDICULOUS) 
      printf("  -");
    else
      printf("%3.0f",tmp->zRes*10000);

    if (tmp->zRes==RIDICULOUS) 
      printf("\t  -\n");
    else 
      printf("\t%0.2f\n",tmp->eff);
    
  }
}

void DetectorK::PlotLayout(Int_t plotDead) {
  //
  // Plots the detector layout in Front view
  //

  Double_t x0=0, y0=0;

  TGraphErrors *gr = new TGraphErrors();
  gr->SetPoint(0,0,0);
  CylLayerK *lastLayer = (CylLayerK*)fLayers.At(fLayers.GetEntries()-1);  Double_t maxRad = lastLayer->radius;
  gr->SetPointError(0,maxRad,maxRad);
  gr->Draw("APE");
  

  CylLayerK *tmp = 0;
  for (Int_t i = fLayers.GetEntries()-1; i>=0; i--) {
    tmp = (CylLayerK*)fLayers.At(i);
  

    Double_t txtpos = tmp->radius;
    if ((tmp->isDead)) txtpos*=-1; //
    TText *txt = new TText(x0,txtpos,tmp->GetName());
    txt->SetTextSizePixels(5); txt->SetTextAlign(21);
    if (!tmp->isDead || plotDead) txt->Draw();

    TEllipse *layEl = new TEllipse(x0,y0,tmp->radius);
    //  layEl->SetFillColor(5);
    layEl->SetFillStyle(5001);
    layEl->SetLineStyle(tmp->isDead+1); // dashed if not active
    layEl->SetLineColor(4);
    TString name(tmp->GetName());
    if (!tmp->isDead) layEl->SetLineWidth(2);
    if (name.Contains("tpc") )  layEl->SetLineColor(29);
    if (name.Contains("trd") )  layEl->SetLineColor(30);

    if (!tmp->isDead || plotDead) layEl->Draw();
  
  }

}



void DetectorK::AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip) {
  //
  // Emulates the TPC
  // 
  // skip=1: Use every padrow, skip=2: Signal in every 2nd padrow 


  AddLayer((char*)"tpcIFC",   77.8,0.01367); // Inner Field cage
  AddLayer((char*)"tpcOFC",   254.0,0.01367); // Outer Field cage

  // % Radiation Lengths ... Average per TPC row  (i.e. total/159 )
  const int kNPassiveBound = 2;
  const Float_t radLBoubdary[kNPassiveBound] = {0.05, 0.0165};
  const Float_t rBoundary[kNPassiveBound] = {50, 70.0}; // cm

  Float_t radLPerRow = 0.000036;
  
  Float_t tpcInnerRadialPitch  =    0.75 ;    // cm
  Float_t tpcMiddleRadialPitch =    1.0  ;    // cm
  Float_t tpcOuterRadialPitch  =    1.5  ;    // cm
  //  Float_t tpcInnerPadWidth     =    0.4  ;    // cm
  //  Float_t tpcMiddlePadWidth    =    0.6   ;   // cm
  //  Float_t tpcOuterPadWidth     =    0.6   ;   // cm
  Float_t innerRows            =   63 ;
  Float_t middleRows           =   64  ;
  Float_t outerRows            =   32  ;
  Float_t tpcRows            =   (innerRows + middleRows + outerRows) ;
  Float_t rowOneRadius         =   85.2  ;    // cm
  Float_t row64Radius          =  135.1  ;    // cm
  Float_t row128Radius         =  199.2  ;    // cm                       
 
  // add boundaries between ITS and TPC
  for (int i=0;i<kNPassiveBound;i++) {
    AddLayer(Form("tpc_boundary%d",i),rBoundary[i],radLBoubdary[i]); // dummy errors
  }

  for ( Int_t k = 0 ; k < tpcRows ; k++ ) {
    
    Float_t rowRadius =0;
    if (k<innerRows) 
      rowRadius =  rowOneRadius + k*tpcInnerRadialPitch ;
    else if ( k>=innerRows && k<(innerRows+middleRows) )
      rowRadius =  row64Radius + (k-innerRows+1)*tpcMiddleRadialPitch ;
    else if (k>=(innerRows+middleRows) && k<tpcRows )
      rowRadius = row128Radius + (k-innerRows-middleRows+1)*tpcOuterRadialPitch ;

    if ( k%skip == 0 )
      AddLayer(Form("tpc_%d",k),rowRadius,radLPerRow,phiResMean,zResMean);    
    else 
      AddLayer(Form("tpc_%d",k),rowRadius,radLPerRow); // non "active" row
    
  
  }
 
}

void DetectorK::AddTRD(Float_t phiResMean, Float_t zResMean, Float_t lrEff) {
  //
  // Emulates the TRD
  // 
  const double trdX2X0=3.3e-2;
  for (int i=0;i<6;i++) AddLayer((char*)Form("trd_%d",i), 300.0+13*i ,trdX2X0, phiResMean, zResMean,
				 lrEff<1 ? lrEff : 1.0); 
 
}

void DetectorK::RemoveTPC() {

  // flag as dead, although resolution is ok ... makes live easier in the prints ... ;-)
  CylLayerK *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (CylLayerK*)fLayers.At(i);  
    TString name(tmp->GetName());
    if (name.Contains("tpc")) { RemoveLayer((char*)name.Data()); i--; }
  }
  
}


Double_t DetectorK::ThetaMCS ( Double_t mass, Double_t radLength, Double_t momentum ) const
{
  //
  // returns the Multiple Couloumb scattering angle (compare PDG boolet, 2010, equ. 27.14)
  //

  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( radLength > 0 ) theta  =  0.0136 * TMath::Sqrt(radLength) / ( beta * momentum ) * (1+0.038*TMath::Log(radLength)) ;
  return (theta) ;
}


Double_t DetectorK::ProbGoodHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Howard Wieman: http://rnc.lbl.gov/~wieman/GhostTracks.htm 
  // and http://rnc.lbl.gov/~wieman/HitFinding2D.htm
  // This is the probability of getting a good hit using 2D Gaussian distribution function and infinite search radius
  Double_t sx, sy, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusRPhi * HitDensity(radius) ;
  sy = 2 * TMath::Pi() *  searchRadiusZ    * searchRadiusZ    * HitDensity(radius) ;
  goodHit =  TMath::Sqrt(1./((1+sx)*(1+sy)))  ;
  return ( goodHit ) ;
}


Double_t DetectorK::ProbGoodChiSqHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Victor Perevoztchikov and Howard Wieman: http://rnc.lbl.gov/~wieman/HitFinding2DXsq.htm
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  Double_t sx, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusZ * HitDensity(radius) ;
  goodHit =  1./(1+sx) ;
  return ( goodHit ) ;  
}

Double_t DetectorK::ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  // Plus, in addition, taking a "confidence level" and the "layer efficiency" into account 
  // Following is correct for 2 DOF

  Double_t c = -2 *TMath::Log(fConfLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t goodHit = leff/(2*alpha) * (1 - TMath::Exp(-alpha*c));
  return ( goodHit ) ;  
}

Double_t DetectorK::ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability to not have any match to the track (see also :ProbGoodChiSqPlusConfHit:)

  Double_t c = -2 *TMath::Log(fConfLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t nullHit = (1-leff+fConfLevel*leff)*TMath::Exp(-c*(alpha-1./2));
  return ( nullHit ) ;  
}

Double_t DetectorK::HitDensity ( Double_t radius ) 
{
  // Background (0-1) is included via 'OtherBackground' which multiplies the minBias rate by a scale factor.
  // UPC electrons is a temporary kludge that is based on Kai Schweda's summary of Kai Hainken's MC results
  // See K. Hencken et al. PRC 69, 054902 (2004) and PPT slides by Kai Schweda.
  // Note that this function assumes we are working in CM and CM**2 [not meters].
  // Based on work by Yan Lu 12/20/2006, all radii and densities in centimeters or cm**2.

  //  Double_t MaxRadiusSlowDet = 0.1; //?   // Maximum radius for slow detectors.  Fast detectors 
                                        // and only fast detectors reside outside this radius.
  Double_t arealDensity = 0 ;

  if ( radius > fMaxRadiusSlowDet ) 
    {
      arealDensity  = OneEventHitDensity(fdNdEtaCent,radius)  ; // Fast detectors see central collision density (only)
      arealDensity += OtherBackground*OneEventHitDensity(dNdEtaMinB,radius)  ;  // Increase density due to background 
    }

  if (radius < fMaxRadiusSlowDet )
    { // Note that IntegratedHitDensity will always be minB one event, or more, even if integration time => zero.
      arealDensity  = OneEventHitDensity(fdNdEtaCent,radius) 
	            + IntegratedHitDensity(dNdEtaMinB,radius) 
	            + UpcHitDensity(radius) ;
      arealDensity += OtherBackground*IntegratedHitDensity(dNdEtaMinB,radius) ;  
      // Increase density due to background 
    } 

  return ( arealDensity ) ;  
}


double DetectorK::OneEventHitDensity( Double_t multiplicity, Double_t radius ) const
{
  // This is for one event at the vertex.  No smearing.

  double den   = multiplicity / (2.*TMath::Pi()*radius*radius) ; // 2 eta ?
  double tg = TMath::Tan(2*TMath::ATan(TMath::Exp(-fAvgRapidity)));
  den = den/TMath::Sqrt(1 + 1/(tg*tg));

  // double den   = multiplicity / (2.*TMath::Pi()*radius*radius) ; // 2 eta ?
  // note: surface of sphere is  '4*pi*r^2'
  //       surface of cylinder is '2*pi*r* h' 

  

  return den ;
} 


double DetectorK::IntegratedHitDensity(Double_t multiplicity, Double_t radius)
{ 
  // The integral of minBias events smeared over a gaussian vertex distribution.
  // Based on work by Yan Lu 12/20/2006, all radii in centimeters.

  Double_t zdcHz = Luminosity * 1.e-24 * CrossSectionMinB ;
  Double_t den   = zdcHz * fIntegrationTime/1000. * multiplicity * Dist(0., radius) / (2.*TMath::Pi()*radius) ;

  // Note that we do not allow the rate*time calculation to fall below one minB event at the vertex.
  if ( den < OneEventHitDensity(multiplicity,radius) )  den = OneEventHitDensity(multiplicity,radius) ;  

  return den ;
} 


double DetectorK::UpcHitDensity(Double_t radius)
{ 
  // QED electrons ...

  Double_t mUPCelectrons ;                                 ;  
  //  mUPCelectrons =  fLhcUPCscale * (1.23 - radius/6.5)      ;  // Fit to Kai Schweda summary tables at RHIC * 'scale' for LHC
  mUPCelectrons = fLhcUPCscale*5456/(radius*radius)/dNdEtaMinB;      // Fit to 'Rossegger,Sadovsky'-Alice simulation
  if ( mUPCelectrons < 0 ) mUPCelectrons =  0.0             ;  // UPC electrons fall off quickly and don't go to large R
  mUPCelectrons *= IntegratedHitDensity(dNdEtaMinB,radius) ;  // UPCs increase Mulitiplicty ~ proportional to MinBias rate
  mUPCelectrons *= UPCBackgroundMultiplier                 ;  // Allow for an external multiplier (eg 0-1) to turn off UPC

  return mUPCelectrons ;
} 


double DetectorK::Dist(double z, double r)
{
  // Convolute dEta/dZ  distribution with assumed Gaussian of vertex z distribution
  // Based on work by Howard Wieman http://rnc.lbl.gov/~wieman/HitDensityMeasuredLuminosity7.htm
  // Based on work by Yan Lu 12/20/2006, all radii and Z location in centimeters.
  Int_t    index  =  1     ;     // Start weight at 1 for Simpsons rule integration
  Int_t    nsteps =  301   ;     // NSteps must be odd for Simpson's rule to work
  double   dist   =  0.0   ;
  double   dz0    =  ( 4*SigmaD - (-4)*SigmaD ) / (nsteps-1)  ;  //cm
  double    z0    =  0.0   ;     //cm
  for(int i=0; i<nsteps; i++){
    if ( i == nsteps-1 ) index = 1 ;
    z0 = -4*SigmaD + i*dz0 ;
    dist += index * (dz0/3.) * (1/sqrt(2.*TMath::Pi())/SigmaD) * exp(-z0*z0/2./SigmaD/SigmaD) * 
      (1/sqrt((z-z0)*(z-z0) + r*r)) ;
    if ( index != 4 ) index = 4; else index = 2 ;
  }
  return dist; 
}

#define  PZero   0.861  // Momentum of back to back decay particles in the CM frame
#define  EPiZero 0.872  // Energy of the pion from a D0 decay at rest
#define  EKZero  0.993  // Energy of the Kaon from a D0 decay at rest

Double_t DetectorK::D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][400] ) const {
  // Math from Ron Longacre.  Note hardwired energy to bin conversion for PtK and PtPi.

  Double_t const1  =  pt / D0Mass ;
  Double_t const2  =  TMath::Sqrt(pt*pt+D0Mass*D0Mass) / D0Mass ;
  Double_t sum, ptPi, ptK ;
  Double_t effp, effk ;

  sum = 0.0 ;
  for ( Int_t k = 0 ; k < 360 ; k++ )   {
    
    Double_t theta = k * TMath::Pi() / 180. ;
    
    ptPi = TMath::Sqrt( 
		       PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*const2*const2 +
		       const1*const1*EPiZero*EPiZero -
		       2*PZero*TMath::Cos(theta)*const2*const1*EPiZero +
		       PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		       ) ;
    
    ptK = TMath::Sqrt( 
		      PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*const2*const2 +
		      const1*const1*EKZero*EKZero +
		      2*PZero*TMath::Cos(theta)*const2*const1*EKZero +
		      PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		      ) ;

    // JT Test Remove 100 MeV/c in pt to simulate eta!=0 decays
    Int_t pionindex = (int)((ptPi-0.1)*100.0 - 65.0*TMath::Abs(fBField)) ; 
    Int_t kaonindex = (int)((ptK -0.1)*100.0 - 65.0*TMath::Abs(fBField)) ; 
      
    if ( pionindex >= kNptBins ) pionindex = 399 ;
    if ( pionindex >= 0 )   effp = corrEfficiency[0][pionindex] ;
    if ( pionindex <  0 )   effp = (corrEfficiency[0][1]-corrEfficiency[0][0])*pionindex + corrEfficiency[0][0] ; // Extrapolate if reqd
    if ( effp < 0 )         effp = 0 ;

    if ( kaonindex >= kNptBins ) kaonindex = 399 ;
    if ( kaonindex >= 0 )   effk = corrEfficiency[1][kaonindex] ;
    if ( kaonindex <  0 )   effk = (corrEfficiency[1][1]-corrEfficiency[1][0])*kaonindex + corrEfficiency[1][0] ; // Extrapolate if reqd
    if ( effk < 0 )         effk = 0 ;

    // Note that we assume that the Kaon Decay efficiency has already been inlcuded in the kaon efficiency used here.
      
    sum += effp * effk ;
 
  }    
  
  Double_t mean =sum/360; 
  return mean ;
  
}



void DetectorK::SolveViaBilloir(Int_t flagD0,Int_t print, Bool_t allPt, Double_t meanPt, char* detLayer) {
  //
  // Solves the current geometry with the Billoir technique 
  // ( see P. Billoir, Nucl. Instr. and Meth. 225 (1984), p. 352. )
  // ABOVE IS OBSOLETE -> NOW, its uses the Aliroot Kalman technique
  //
  const float kTrackingMargin = 0.1;

  static AliExternalTrackParam probTr;   // track to propagate
  probTr.SetUseLogTermMS(kTRUE);


  Int_t nPt = kNptBins;
  // Clean up ......
  for (Int_t j=0; j<nPt; j++) {
    for (Int_t i=0; i<kMaxNumberOfDetectors; i++) {
      fDetPointRes[i][j]  = RIDICULOUS;
      fDetPointZRes[i][j] = RIDICULOUS;
    }
    fTransMomenta[j] =0;
    fMomentumRes[j] =0;
    fResolutionRPhi[j] =0;
  }
  
  if (!allPt) { // not the whole pt range -> allows a faster minimization at a defined 'meanpt'
    nPt = 3;
  }


  // Calculate track parameters using Billoirs method of matrices

  Double_t pt,tgl, lambda, deltaPoverP  ;
  Double_t charge ;
  Double_t mass[3] ;
  Int_t printOnce = 1 ;

  mass[0] = PionMass ; mass[1] = KaonMass ;  // Loop twice for the D0;  first pi then k 

  mass[2] = fParticleMass;  // third loop

  Int_t mStart =0; 
  if (!flagD0) mStart = 2; // pion and kaon is skipped -> fast mode

   

  // Prepare Probability Kombinations
  Int_t nLayer = fNumberOfActiveITSLayers;
  Int_t base = 3; // null, fake, correct

  Int_t komb = (Int_t) TMath::Power(base,nLayer);

  printf("N ITS Layers: %d\n",fNumberOfActiveITSLayers);

  TMatrixD probLay(base,fNumberOfActiveITSLayers);
  TMatrixD probKomb(komb,nLayer);
  for (Int_t num=0; num<komb; num++) {
    for (Int_t l=nLayer; l--;) {
      Int_t pow = ((Int_t)TMath::Power(base,l+1));
      probKomb(num,nLayer-1-l)=(num%pow)/((Int_t)TMath::Power(base,l));
    }
  }

  TString detLayerStr(detLayer);
  CylLayerK *theLayer = (CylLayerK*) fLayers.FindObject(detLayer);
  if (!theLayer && detLayerStr.IsNull()!=1){
    printf("Error: Layer with the name \"%s\" not found -> no detailed infos possible\n",detLayer);
    return;
  }

  for (Int_t i=0; i<fLayers.GetEntries();i++) {
    CylLayerK *l = (CylLayerK*) fLayers.At(i);
    if (detLayerStr.CompareTo(l->GetName())==0) { // is the same
      kDetLayer=i;
      break;
    }
  }

  CylLayerK *last = (CylLayerK*) fLayers.At((fLayers.GetEntries()-1));
  if (last->radius > fMinRadTrack) {
    last = 0;
    for (Int_t i=0; i<fLayers.GetEntries();i++) {
      CylLayerK *l = (CylLayerK*) fLayers.At(i);
      if (!(l->isDead) && (l->radius<fMinRadTrack)) last = l;
    }
    if (!last) {
      printf("No layer with radius < %f is found\n",fMinRadTrack);
      return;
    }
  }

  Double_t bigRad = last->radius/2 ;	// min. pt which the algorithm below could handle  
  double ptmin = ( 0.3*bigRad*TMath::Abs(fBField)*1e-2 ) + 0.005; // safety margin
  if (ptmin<kPtMinFix) ptmin = kPtMinFix;
  double ptmax = kPtMaxFix;
  double dlpt = log(ptmax/ptmin)/nPt;
  
  
  for ( Int_t massloop = mStart ; massloop < 3 ; massloop++ )  { 
    
    // PseudoRapidity OK, used as an angle
    lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity))  ; 
  

    for ( Int_t i = 0 ; i < nPt ; i++ ) { // pt loop
      //
      // Starting values based on radius of outermost layer ... log10 steps to ~20 GeV
      //   if (bigRad<61) bigRad=61; // -> min pt around 100 MeV for Bz=0.5T (don't overdo it ... ;-) )
      fTransMomenta[i] =  ptmin*TMath::Exp(dlpt*i);
      //fTransMomenta[i] = ( 0.3*bigRad*TMath::Abs(fBField)*1e-2 ) - 0.08 - (1./fptScale-0.1) + TMath::Power(10,2.3*i/nPt) / fptScale ; 
      if (!allPt) { // just 3 points around meanPt
	fTransMomenta[i] = meanPt-0.001+(Double_t)(i)*0.001;
      }
  
      // New from here ................

      // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
      // These are the EndPoint values for y, z, a, b, and d
      double bGauss = fBField*10;               // field in kgauss
      pt  =  fTransMomenta[i];                  // GeV/c
      tgl =  TMath::Tan(lambda);                // dip
      charge   = -1;                            // Assume an electron 
      enum {kY,kZ,kSnp,kTgl,kPtI};              // track parameter aliases
      enum {kY2,kYZ,kZ2,kYSnp,kZSnp,kSnp2,kYTgl,kZTgl,kSnpTgl,kTgl2,kYPtI,kZPtI,kSnpPtI,kTglPtI,kPtI2}; // cov.matrix aliases
      //
      probTr.Reset();
      double *trPars = (double*)probTr.GetParameter();
      double *trCov  = (double*)probTr.GetCovariance();
      trPars[kY] = 0;                         // start from Y = 0
      trPars[kZ] = 0;                         //            Z = 0 
      trPars[kSnp] = 0;                       //            track along X axis at the vertex
      trPars[kTgl] = tgl;      //            dip
      trPars[kPtI] = charge/pt;               //            q/pt      
      //
      // put tiny errors to propagate to the outer radius
      trCov[kY2] = trCov[kZ2] = trCov[kSnp2] = trCov[kTgl2] = trCov[kPtI2] = 1e-9;
      //
      // find max layer this track can reach
      double rmx = (TMath::Abs(fBField)>1e-5) ?  pt*100./(0.3*TMath::Abs(fBField)) : 9999;
      Int_t lastActiveLayer = -1;
      for (Int_t j=fLayers.GetEntries(); j--;) { 
	CylLayerK *l = (CylLayerK*) fLayers.At(j);
	//	printf("at lr %d r: %f vs %f, pt:%f\n",j,l->radius, 2*rmx-2.*kTrackingMargin, pt);
	if (!(l->isDead) && (l->radius < 2*(rmx-5.))) {lastActiveLayer = j; last = l; break;}
      }
      if (lastActiveLayer<0) {
	printf("No active layer with radius < %f is found, pt = %f\n",rmx, pt);
	return;
      }
      //      printf("PT=%f 2Rpt=%f Rlr=%f\n",pt,2*rmx,last->radius);
      //
      if (!PropagateToR(&probTr,last->radius + kTrackingMargin,bGauss,1)) continue;
      //if (!probTr.PropagateTo(last->radius,bGauss)) continue;
      // reset cov.matrix
      const double kLargeErr2Coord = 5*5;
      const double kLargeErr2Dir = 0.7*0.7;
      const double kLargeErr2PtI = 30.5*30.5;
      for (int ic=15;ic--;) trCov[ic] = 0.;
      trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
      trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
      trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
      probTr.CheckCovariance();
      //
      //      printf("%d - pt %lf r%lf | %lf %lf\n",massloop,fTransMomenta[i],(last->radius)/100,momentum, d);

      // Set Detector-Efficiency Storage area to unity
      fEfficiency[massloop][i] = 1.0 ;
      //
      // Back-propagate the covariance matrix along the track. 
    
      CylLayerK *layer = 0;

      /* remove      
      // find last "active layer" - start tracking at the last active layer      
      Int_t lastActiveLayer = 0;
      for (Int_t j=fLayers.GetEntries(); j--;) { 
	layer = (CylLayerK*)fLayers.At(j);
	if (!(layer->isDead)) { // is alive
	  lastActiveLayer = j;
	  break;
	}
      }
      */
      //      probTr.Print();
      for (Int_t j=lastActiveLayer+1; j--;) {  // Layer loop

	layer = (CylLayerK*)fLayers.At(j);

	if (layer->radius>fMaxSeedRadius) continue; // no seeding beyond this radius 

	TString name(layer->GetName());
	Bool_t isVertex = name.Contains("vertex");
	//
	if (!PropagateToR(&probTr,layer->radius,bGauss,-1)) exit(1);
	//	if (!probTr.PropagateTo(last->radius,bGauss)) exit(1);	//
	// rotate to frame with X axis normal to the surface
	if (!isVertex) {
	  double pos[3];
	  probTr.GetXYZ(pos);  // lab position
	  double phi = TMath::ATan2(pos[1],pos[0]);
	  if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
	  if (!probTr.Rotate(phi)) {
	    printf("Failed to rotate to the frame (phi:%+.3f)of layer at %.2f at XYZ: %+.3f %+.3f %+.3f (pt=%+.3f)\n",
		   phi,layer->radius,pos[0],pos[1],pos[2],pt);
	    
	    probTr.Print();
	    exit(1);
	  }
	}
	/*
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) 
	  {
	    printf("\nAt lr %d %s R: %f\n ",j,layer->GetName(), layer->radius);
	    probTr.Print();
	  }
	//	*/
	// save resolutions at this layer
	fDetPointRes [j][i]     =  TMath::Sqrt( probTr.GetSigmaY2() )/100  ;     // result in meters
	fDetPointZRes[j][i]     =  TMath::Sqrt( probTr.GetSigmaZ2() )/100  ;     // result in meters
	//printf(">> L%d r:%e sy: %e sz: %e\n",j,layer->radius,fDetPointRes[j][i],fDetPointZRes[j][i]);
	// End save
	//
	if (isVertex) continue;
	//
	// create fake measurement with the errors assigned to the layer
	// account for the measurement there 
	double meas[2] = {probTr.GetY(),probTr.GetZ()};
	double measErr2[3] = {layer->phiRes*layer->phiRes,0,layer->zRes*layer->zRes};
	//

	if (!probTr.Update(meas,measErr2)) {
	  printf("Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}\n",
		 meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]);
	  probTr.Print();
	  exit(1);
	}
	//printf("AfterUpdate "); probTr.Print();
	// correct for materials of this layer
	// note: if apart from MS we want also e.loss correction, the density*length should be provided as 2nd param
	if (!probTr.CorrectForMeanMaterial(layer->radL, 0, mass[massloop] , kTRUE)) {
	  printf("Failed to apply material correction, X/X0=%.4f\n",layer->radL);
	  probTr.Print();
	  exit(1);
	}
	//printf("AfterCorr "); probTr.Print();
	//
      }
    
      // Pattern recognition is done .... save values like vertex resolution etc.

      // Convert the Convariance matrix parameters into physical quantities
      // The results are propogated to the previous point but *do not* include the measurement at that point.
      //      deltaPoverP          =  TMath::Sqrt(probTr.GetSigma1Pt2())/probTr.Get1P();  // Absolute magnitude so ignore charge
      deltaPoverP          =  TMath::Sqrt(probTr.GetSigma1Pt2())/TMath::Abs(probTr.GetSigned1Pt());
      fMomentumRes[i]      =  100.* TMath::Abs( deltaPoverP );                    // results in percent
      fResolutionRPhi[i]   =  TMath::Sqrt( probTr.GetSigmaY2() ) * 1.e4;          // result in microns
      fResolutionZ[i]      =  TMath::Sqrt( probTr.GetSigmaZ2() ) * 1.e4;          // result in microns
      //      equivalent[i]  =  TMath::Sqrt(fResolutionRPhi[i]*fResolutionZ[i])           ;  // Equivalent circular radius
      //
      if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	printf("Number of active layers: %d, last Layer reached: %d\n",fNumberOfActiveLayers,lastActiveLayer);
	if (fAtLeastCorr != -1) printf("Number of combinatorics for probabilities: %d\n",komb);
	printf("Mass of tracked particle: %f (at pt=%5.0lf MeV)\n",fParticleMass,fTransMomenta[i]*1000);
	printf("Name   Radius Thickness PointResOn PointResOnZ  DetRes  DetResZ  Density Efficiency\n") ;
	//	printOnce =0;
      }
      
      // print out and efficiency calculation
      Int_t iLayActive=0;
      //      for (Int_t j=(fLayers.GetEntries()-1); j>=0; j--) {  // Layer loop
      for (Int_t j=lastActiveLayer+1; j--;) {  // Layer loop
	
	layer = (CylLayerK*)fLayers.At(j);
	
	// Convert to Meters, Tesla, and GeV
	Float_t radius = layer->radius /100;
	Float_t phiRes = layer->phiRes /100;
	Float_t zRes = layer->zRes /100;
	Float_t radLength = layer->radL;
	Float_t leff = layer->eff; // basic layer efficiency
	Bool_t isDead = layer->isDead;
	

	if ( (!isDead && radLength >0) )  { 

	    Double_t rphiError  =  TMath::Sqrt( fDetPointRes[j][i] * fDetPointRes [j][i] + 
						phiRes * phiRes ) * 100.  ; // work in cm
	    Double_t zError     =  TMath::Sqrt( fDetPointZRes[j][i] * fDetPointZRes[j][i] +
						zRes * zRes ) * 100.  ; // work in cm
	    
	    
	    Double_t layerEfficiency = 0;
	    if ( EfficiencySearchFlag == 0 )
	      layerEfficiency =  ProbGoodHit( radius*100, rphiError , zError  ) ;
	    else if ( EfficiencySearchFlag == 1 )
	      layerEfficiency =  ProbGoodChiSqHit( radius*100, rphiError , zError  ) ;
	    else if ( EfficiencySearchFlag == 2 )
	      layerEfficiency =  ProbGoodChiSqPlusConfHit( radius*100,leff, rphiError , zError  ) ;
	      
	    TString name(layer->GetName());
	    if ( IsITSLayer(name) ) {
	      probLay(2,iLayActive)= layerEfficiency ; // Pcorr
	      probLay(0,iLayActive)= ProbNullChiSqPlusConfHit( radius*100,leff, rphiError , zError  ) ; // Pnull
	      probLay(1,iLayActive)= 1 - probLay(2,iLayActive) - probLay(0,iLayActive);                 // Pfake
	      iLayActive++;    
	    }
	    if (!IsITSLayer(name) && (!name.Contains("tpc_0")) ) continue;

	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) 
	    {
	      printf("%s:\t%5.1f %9.4f %10.0f %11.0f %7.0f %8.0f %8.2f ",
		     layer->GetName(), radius*100, radLength, 
		     fDetPointRes[j][i]*1.e6, fDetPointZRes[j][i]*1.e6,
		     phiRes*1.e6, zRes*1.e6,
		     HitDensity(radius*100)) ;
	      if (!name.Contains("tpc")) 
		printf("%10.3f\n", layerEfficiency);
	      else
		printf("        -  \n");
	    }
	    
	    if (IsITSLayer(name))   fEfficiency[massloop][i] *= layerEfficiency;
	    
	    
	}
	
	if (fAtLeastCorr != -1 || fAtLeastHits) {
	  // Calculate probabilities from Kombinatorics tree ...
	  Double_t *probs = PrepareEffFakeKombinations(&probKomb, &probLay);
	  fEfficiency[massloop][i] = probs[0]; // efficiency
	  fFake[massloop][i] = probs[1];       // fake
	}

	/*
	// vertex print
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1 && radius==0) {
	printf("%s:\t -----    ----- %10.0f %11.0f \n", layer->GetName(),fDetPointRes[j][i]*1.e6, fDetPointZRes[j][i]*1.e6);
	}
	*/
      }
      if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	if (fNumberOfActiveLayers >=1500) printOnce = 0 ;
	printf("\n")  ;
      }

      if (fNumberOfActiveLayers <1500 ) {

	//      printf("Backward PtBin%d pt=%f\n",i,pt);

	// BACKWORD TRACKING +++++++++++++++++
	// number of layers is quite low ... efficiency calculation was probably nonsense 
	// Tracking outward (backword) to get reliable efficiencies from "smoothed estimates"

	// For below, see paper, NIM A262 (1987) p.444, eqs.12.
	// Equivalently, one can simply combine the forward and backward estimates. Assuming
	// pf,Cf and pb,Cb as extrapolated position estimates and errors from fwd and bwd passes one can
	// use a weighted estimate Cw = (Cf^-1 + Cb^-1)^-1,  pw = Cw (pf Cf^-1 + pb Cb^-1).
	// Surely, for the most extreme point, where one error matrices is infinite, this does not change anything.

	Bool_t doLikeAliRoot = 0; // don't do the "combined info" but do like in Aliroot

	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {	  
	  printf("- Numbers of active layer is low (%d):\n    -> \"outward\" fitting done as well to get reliable eff.estimates\n",
		 fNumberOfActiveLayers);	  
	}
	
	// RESET Covariance Matrix ( to 10 x the estimate -> as it is done in AliExternalTrackParam)
	//	mIstar.UnitMatrix(); // start with unity
	if (doLikeAliRoot) {
	  probTr.ResetCovariance(100);
	} else {
	  // cannot do complete reset, set to very large errors
	  for (int ic=15;ic--;) trCov[ic] = 0.;
	  trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
	  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
	  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
	  probTr.CheckCovariance();
	  //	  cout<<pt<<": "<<kLargeErr2Coord<<" "<<kLargeErr2Dir<<" "<<kLargeErr2PtI*trPars[kPtI]*trPars[kPtI]<<endl;
	}
	// Clean up and storing of "forward estimates"
	Double_t detPointResForw[kMaxNumberOfDetectors][kNptBins], detPointZResForw[kMaxNumberOfDetectors][kNptBins] ; 
	Double_t detPointResBwd[kMaxNumberOfDetectors][kNptBins], detPointZResBwd[kMaxNumberOfDetectors][kNptBins] ; 
	for (Int_t k=0; k<kMaxNumberOfDetectors; k++) {
	  for (Int_t l=0; l<nPt; l++) {
	    detPointResForw[k][l]  = fDetPointRes[k][l];
	    if (!doLikeAliRoot) fDetPointRes[k][l]  = RIDICULOUS;
	    detPointZResForw[k][l] = fDetPointZRes[k][l];
	    if (!doLikeAliRoot) fDetPointZRes[k][l] = RIDICULOUS;
	    detPointResBwd[k][l] = detPointZResBwd[k][l] = RIDICULOUS;
	  }
	}
 	
	// find first "active layer" - start tracking at the first active layer      
	Int_t firstActiveLayer = 0;
	for (Int_t j=0; j<=lastActiveLayer; j++) { 
	  layer = (CylLayerK*)fLayers.At(j);
	  if (!(layer->isDead)) { // is alive
	    firstActiveLayer = j;
	    break;
	  }
	}
	//probTr.Rotate(0);
   	for (Int_t j=firstActiveLayer; j<=lastActiveLayer; j++) {  // Layer loop
	  
	  layer = (CylLayerK*)fLayers.At(j);
	  //  CylLayerK *nextlayer = (CylLayerK*)fLayers.At(j+1);

	  TString name(layer->GetName());
	  Bool_t isVertex = name.Contains("vertex");
	  if (!PropagateToR(&probTr, layer->radius,bGauss,1)) exit(1);
	  //if (!probTr.PropagateTo(last->radius,bGauss))  exit(1);
	  if (!isVertex) {
	    // rotate to frame with X axis normal to the surface
	    double pos[3];
	    probTr.GetXYZ(pos);  // lab position
	    double phi = TMath::ATan2(pos[1],pos[0]);
	    if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
	    if (!probTr.Rotate(phi)) {
	      printf("Failed to rotate to the frame (phi:%+.3f)of layer at %.2f at XYZ: %+.3f %+.3f %+.3f (pt=%+.3f)\n",
		     phi,layer->radius,pos[0],pos[1],pos[2],pt);	      
	      probTr.Print();
	      exit(1);
	    }
	  }
	  /*
	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) 
	    {
	    printf("\nAt lr %d %s R: %f\n ",j,layer->GetName(), layer->radius);
	    probTr.Print();
	    }
	  */
	  //	  
	  detPointResBwd[j][i]     =  TMath::Sqrt( probTr.GetSigmaY2() )/100  ;     // result in meters
	  detPointZResBwd[j][i]    =  TMath::Sqrt( probTr.GetSigmaZ2() )/100  ;     // result in meters
	  //
	  //printf("<< L%d r:%e sy: %e sz: %e\n",j,layer->radius,fDetPointRes[j][i],fDetPointZRes[j][i]);
	  // create fake measurement with the errors assigned to the layer
	  // account for the measurement there
	  if (isVertex) continue;
	  double meas[2] = {probTr.GetY(),probTr.GetZ()};
	  double measErr2[3] = {layer->phiRes*layer->phiRes,0,layer->zRes*layer->zRes};
	  //
	  if (!probTr.Update(meas,measErr2)) {
	    printf("Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}\n",
		   meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]);
	    probTr.Print();
	    exit(1);
	  }
	  //printf("AfterUpdate "); probTr.Print();
	  // correct for materials of this layer
	  // note: if apart from MS we want also e.loss correction, the density*length should be provided as 2nd param
	  if (!probTr.CorrectForMeanMaterial(layer->radL, 0, mass[massloop] , kTRUE)) {
	    printf("Failed to apply material correction, X/X0=%.4f\n",layer->radL);
	    probTr.Print();
	    exit(1);
	  }
	  //printf("AfterCorr "); probTr.Print();
	}

	// values below NOT REALIABLE -> they do not point to the vertex but outwards !!!!!!!
	// ++++++++++++++
	// also update the values for the track position ??????
	/*
	// Pattern recognition is done .... save values like vertex resolution etc.
	
	// Invert the Matrix to recover the convariance matrix
	mIstar.Invert() ;
	// Convert the Convariance matrix parameters into physical quantities
	// The results are propogated to the previous point but *do not* include the measurement at that point.
	deltaPoverP    =  TMath::Sqrt( mIstar(4,4) ) * momentum / 0.3  ;  // Absolute magnitude so ignore charge
	fMomentumRes[i]      =  100.* TMath::Abs( deltaPoverP )             ;  // results in percent
	fResolutionRPhi[i]      =  TMath::Sqrt( mIstar(0,0) ) * 1.e6            ;  // result in microns
	fResolutionZ[i]     =  TMath::Sqrt( mIstar(1,1) ) * 1.e6            ;  // result in microns
	//      equivalent[i]  =  TMath::Sqrt(fResolutionRPhi[i]*fResolutionZ[i])           ;  // Equivalent circular radius
        */
	
	//	deltaPoverP          =  TMath::Sqrt(probTr.GetSigma1Pt2())/TMath::Abs(probTr.GetSigned1Pt());
	//	fMomentumRes[i]      =  100.* TMath::Abs( deltaPoverP );                    // results in percent
  

	

	// Weighted combination of the forward and backward estimates
	if (!doLikeAliRoot) {

	  if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) printf("\nBackward propagation estimates\n");
	  
	  for (Int_t j=lastActiveLayer+1; j--;) {  
	    //
	    fDetPointRes[j][i]  = detPointResForw[j][i]*detPointResBwd[j][i]/TMath::Sqrt((detPointResForw[j][i]*detPointResForw[j][i])     + (detPointResBwd[j][i]*detPointResBwd[j][i])); 
	    fDetPointZRes[j][i] = detPointZResForw[j][i]*detPointZResBwd[j][i]/TMath::Sqrt((detPointZResForw[j][i]*detPointZResForw[j][i]) + (detPointZResBwd[j][i]*detPointZResBwd[j][i])); 
	    //
	    layer = (CylLayerK*)fLayers.At(j);
	    
	    TString name(layer->GetName());
	    if ( layer->isDead || ( !IsITSLayer(name) && (!name.Contains("tpc_0"))) ) continue;
	    
	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) 
	    {
	      //
	      Float_t radius = layer->radius /100;
	      Float_t phiRes = layer->phiRes /100;
	      Float_t zRes = layer->zRes /100;
	      Float_t radLength = layer->radL;
	      Float_t leff = layer->eff; // basic layer efficiency
	      Double_t rphiError  =  TMath::Sqrt( detPointResBwd[j][i] * detPointResBwd[j][i] + 
						  phiRes * phiRes ) * 100.  ; // work in cm
	      Double_t zError     =  TMath::Sqrt( detPointZResBwd[j][i] * detPointZResBwd[j][i] +
						  zRes * zRes ) * 100.  ; // work in cm
	      //
	      Double_t layerEfficiency = 0;
	      if ( EfficiencySearchFlag == 0 )
		layerEfficiency =  ProbGoodHit( radius*100, rphiError , zError  ) ;
	      else if ( EfficiencySearchFlag == 1 )
		layerEfficiency =  ProbGoodChiSqHit( radius*100, rphiError , zError  ) ;
	      else if ( EfficiencySearchFlag == 2 )
		layerEfficiency =  ProbGoodChiSqPlusConfHit( radius*100,leff, rphiError , zError  ) ;
	      
	      
	      printf("%s:\t%5.1f %9.4f %10.0f %11.0f %7.0f %8.0f %8.2f ",
		     layer->GetName(), radius*100, radLength, 
		     detPointResBwd[j][i]*1.e6, detPointZResBwd[j][i]*1.e6,
		     phiRes*1.e6, zRes*1.e6,
		     HitDensity(radius*100)) ;
	      if (IsITSLayer(name)) 
		printf("%10.3f\n", layerEfficiency);
	      else
		printf("        -  \n");	    	      
	    }
	  }
	}
	// Set Detector-Efficiency Storage area to unity
	fEfficiency[massloop][i] = 1.0 ;
     
	// print out and efficiency calculation
	iLayActive=0;
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) printf("\n Combined propagation estimates\n");

	for (Int_t j=lastActiveLayer+1;j--;) {  // Layer loop
	  
	  layer = (CylLayerK*)fLayers.At(j);
	  
	  // Convert to Meters, Tesla, and GeV
	  Float_t radius = layer->radius /100;
	  Float_t phiRes = layer->phiRes /100;
	  Float_t zRes = layer->zRes /100;
	  Float_t radLength = layer->radL;
	  Float_t leff = layer->eff;
	  Bool_t isDead = layer->isDead;

	  Double_t layerEfficiency = 0;
	  if ( (!isDead && radLength >0) )  { 
	    Double_t rphiError  =  TMath::Sqrt( fDetPointRes[j][i] * fDetPointRes [j][i] + 
						phiRes * phiRes ) * 100.  ; // work in cm
	    Double_t zError     =  TMath::Sqrt( fDetPointZRes[j][i] * fDetPointZRes[j][i] +
						zRes * zRes ) * 100.  ; // work in cm
	    if ( EfficiencySearchFlag == 0 )
	      layerEfficiency =  ProbGoodHit( radius*100, rphiError , zError  ) ;
	    else if ( EfficiencySearchFlag == 1 )
	      layerEfficiency =  ProbGoodChiSqHit( radius*100, rphiError , zError  ) ;
	    else if ( EfficiencySearchFlag == 2 )
	      layerEfficiency =  ProbGoodChiSqPlusConfHit( radius*100,leff, rphiError , zError  ) ;
	      
	    TString name(layer->GetName());
	    if (IsITSLayer(name)) {
	      probLay(2,iLayActive)= layerEfficiency ; // Pcorr
	      probLay(0,iLayActive)= ProbNullChiSqPlusConfHit( radius*100,leff, rphiError , zError  ) ; // Pnull
	      probLay(1,iLayActive)= 1 - probLay(2,iLayActive) - probLay(0,iLayActive);                 // Pfake
	      iLayActive++;    
	    }
	    if (!IsITSLayer(name) && (!name.Contains("tpc_0")) ) continue;

	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) 
	    {
	      printf("%s:\t%5.1f %9.4f %10.0f %11.0f %7.0f %8.0f %8.2f ",
		     layer->GetName(), radius*100, radLength, 
		     fDetPointRes[j][i]*1.e6, fDetPointZRes[j][i]*1.e6,
		     phiRes*1.e6, zRes*1.e6,
		     HitDensity(radius*100)) ;
	      if (IsITSLayer(name)) 
		printf("%10.3f\n", layerEfficiency);
	      else
		printf("        -  \n");
	    }

	    if (massloop==2 && j==kDetLayer) 	{ // copy layer specific performances
	      fEfficProlongLay[i] = layerEfficiency;
	    }

	    if (IsITSLayer(name))   fEfficiency[massloop][i] *= layerEfficiency;



	  }
	  if (fAtLeastCorr != -1 || fAtLeastHits != -1 ) {
	    // Calculate probabilities from Kombinatorics tree ...
	    Double_t *probs = PrepareEffFakeKombinations(&probKomb, &probLay);
	    fEfficiency[massloop][i] = probs[0]; // efficiency
	    fFake[massloop][i] = probs[1];       // fake
	  }
	}
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	  printOnce = 0 ;
	  printf("\n")  ;
	}
      }      

      if (massloop==2) 	{ // copy layer specific performances
	fResolutionRPhiLay[i] = fDetPointRes[kDetLayer][i];
	fResolutionZLay[i] = fDetPointZRes[kDetLayer][i];
      }
      
    } // pt loop
    


  } // mass loop

  probTr.SetUseLogTermMS(kFALSE); // Reset of MS term usage to avoid problems since its static

  
 
}

Bool_t DetectorK::SolveTrack(TrackSol& ts) {
  //
  // Solves the current geometry for single track of given kinematics
  //
  double ptTr = ts.fPt;
  double etaTr = ts.fEta;
  double mass = ts.fMass;
  double charge = ts.fCharge;
  
  if (ptTr<0) { 
    printf("Input track is not initialized");
    return kFALSE;
  }
  
  const float kTrackingMargin = 0.1;

  static AliExternalTrackParam probTr;   // track to propagate
  probTr.SetUseLogTermMS(kTRUE);
  //
  TClonesArray &saveParInward    = ts.fTrackInw;
  TClonesArray &saveParOutwardB  = ts.fTrackOutB;
  TClonesArray &saveParOutwardA  = ts.fTrackOutA;
  TClonesArray &saveParComb      = ts.fTrackCmb;

  // Calculate track parameters using Billoirs method of matrices
  Double_t pt,lambda;
  //
  CylLayerK *last = (CylLayerK*) fLayers.At((fLayers.GetEntries()-1));
  double maxR = last->radius+kTrackingMargin*2;
  double minRad = (fMinRadTrack>0&&fMinRadTrack<maxR) ? fMinRadTrack : maxR;
  //
  if (last->radius > minRad) {
    last = 0;
    for (Int_t i=0; i<fLayers.GetEntries();i++) {
      CylLayerK *l = (CylLayerK*) fLayers.At(i);
      if (/*!(l->isDead) && */(l->radius<minRad)) last = l;
    }
    if (!last) {
      printf("No layer with radius < %f is found\n",minRad);
      return kFALSE;
    }
  }
  //  
  lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-etaTr));
  //  
  // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
  // These are the EndPoint values for y, z, a, b, and d
  double bGauss = fBField*10;               // field in kgauss
  pt  =  ptTr;
  enum {kY,kZ,kSnp,kTgl,kPtI};              // track parameter aliases
  enum {kY2,kYZ,kZ2,kYSnp,kZSnp,kSnp2,kYTgl,kZTgl,kSnpTgl,kTgl2,kYPtI,kZPtI,kSnpPtI,kTglPtI,kPtI2}; // cov.matrix aliases
  //
  probTr.Reset();
  double *trPars = (double*)probTr.GetParameter();
  double *trCov  = (double*)probTr.GetCovariance();
  trPars[kY] = 0;                         // start from Y = 0
  trPars[kZ] = 0;                         //            Z = 0 
  trPars[kSnp] = 0;                       //            track along X axis at the vertex
  trPars[kTgl] = TMath::Tan(lambda);      //            dip
  trPars[kPtI] = charge/pt;               //            q/pt      
  //
  // put tiny errors to propagate to the outer radius
  trCov[kY2] = trCov[kZ2] = trCov[kSnp2] = trCov[kTgl2] = trCov[kPtI2] = 1e-9;
  //
  // find max layer this track can reach
  double rmx = (TMath::Abs(fBField)>1e-5) ?  pt*100./(0.3*TMath::Abs(fBField)) : 9999;
  if (2*rmx-5. < minRad && minRad>0) {
    printf("Track of pt=%.3f cannot be tracked to min. r=%f\n",pt,minRad);
    return kFALSE;
  }
  Int_t lastActiveLayer = -1;
  for (Int_t j=fLayers.GetEntries(); j--;) { 
    CylLayerK *l = (CylLayerK*) fLayers.At(j);
    //	printf("at lr %d r: %f vs %f, pt:%f\n",j,l->radius, 2*rmx-2.*kTrackingMargin, pt);
    if (/*!(l->isDead) && */(l->radius <= 2*(rmx-5))) {lastActiveLayer = j; last = l; break;}
  }
  if (lastActiveLayer<0) {
    printf("No active layer with radius < %f is found, pt = %f\n",rmx, pt);
    return kFALSE;
  }
  //      printf("PT=%f 2Rpt=%f Rlr=%f\n",pt,2*rmx,last->radius);
  //
  if (!PropagateToR(&probTr,last->radius + kTrackingMargin,bGauss,1)) return kFALSE;
  //if (!probTr.PropagateTo(last->radius,bGauss)) continue;
  // reset cov.matrix
  // 
  // rotate to external layer frame
  /*
  double posL[3];
  probTr.GetXYZ(posL);  // lab position
  double phiL = TMath::ATan2(posL[1],posL[0]);
  if (!probTr.Rotate(phiL)) {
    printf("Failed to rotate to the frame (phi:%+.3f)of Extertnal layer at %.2f\n",
	   phiL,last->radius);
    probTr.Print();
    exit(1);
  }
  */
  if (!probTr.Rotate(probTr.Phi())) return kFALSE; // define large errors in track proper frame (snp=0)
  //
  const double kLargeErr2Coord = 5*5;
  const double kLargeErr2Dir = 0.7*0.7;
  const double kLargeErr2PtI = 30.5*30.5;
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  probTr.CheckCovariance();
  //
  // Back-propagate the covariance matrix along the track.   
  CylLayerK *layer = 0;
  //
  for (Int_t j=lastActiveLayer+1; j--;) {  // Layer loop
    
    layer = (CylLayerK*)fLayers.At(j);
    
    if (layer->radius>fMaxSeedRadius) continue; // no seeding beyond this radius 
    
    TString name(layer->GetName());
    Bool_t isVertex = name.Contains("vertex");
    //
    if (!PropagateToR(&probTr,layer->radius,bGauss,-1)) exit(1);
    //	if (!probTr.PropagateTo(last->radius,bGauss)) exit(1);	//
    // rotate to frame with X axis normal to the surface
    if (!isVertex) {
      double pos[3];
      probTr.GetXYZ(pos);  // lab position
      double phi = TMath::ATan2(pos[1],pos[0]);
      if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
      if (!probTr.Rotate(phi)) {
	printf("Failed to rotate to the frame (phi:%+.3f)of layer at %.2f at XYZ: %+.3f %+.3f %+.3f (pt=%+.3f)\n",
	       phi,layer->radius,pos[0],pos[1],pos[2],pt);
	
	probTr.Print();
	exit(1);
      }
    }
    // save inward parameters at this layer: before the update!
    new( saveParInward[j] ) AliExternalTrackParam(probTr);
    if (verboseR) {
      printf("SaveInw %d (%f)  ",j,layer->radius); probTr.Print();
    }    
    //
    if (!isVertex && !layer->isDead) {
      //
      // create fake measurement with the errors assigned to the layer
      // account for the measurement there 
      double meas[2] = {probTr.GetY(),probTr.GetZ()};
      double measErr2[3] = {layer->phiRes*layer->phiRes,0,layer->zRes*layer->zRes};
      //
      if (!probTr.Update(meas,measErr2)) {
	printf("Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}\n",
	       meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]);
	probTr.Print();
	exit(1);
      }
    }
    // correct for materials of this layer
    // note: if apart from MS we want also e.loss correction, the density*length should be provided as 2nd param
    if (layer->radL>0 && !probTr.CorrectForMeanMaterial(layer->radL, 0, mass , kTRUE)) {
      printf("Failed to apply material correction, X/X0=%.4f\n",layer->radL);
      probTr.Print();
      exit(1);
    }
  }
  //  
  // BACKWORD TRACKING +++++++++++++++++
  // number of layers is quite low ... efficiency calculation was probably nonsense 
  // Tracking outward (backword) to get reliable efficiencies from "smoothed estimates"
  
  // For below, see paper, NIM A262 (1987) p.444, eqs.12.
  // Equivalently, one can simply combine the forward and backward estimates. Assuming
  // pf,Cf and pb,Cb as extrapolated position estimates and errors from fwd and bwd passes one can
  // use a weighted estimate Cw = (Cf^-1 + Cb^-1)^-1,  pw = Cw (pf Cf^-1 + pb Cb^-1).
  // Surely, for the most extreme point, where one error matrices is infinite, this does not change anything.
  
  Bool_t doLikeAliRoot = 0; // don't do the "combined info" but do like in Aliroot
  
  
  // RESET Covariance Matrix ( to 10 x the estimate -> as it is done in AliExternalTrackParam)
  //	mIstar.UnitMatrix(); // start with unity
  if (doLikeAliRoot) {
    probTr.ResetCovariance(100);
  } else {
    // cannot do complete reset, set to very large errors
    for (int ic=15;ic--;) trCov[ic] = 0.;
    trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord; 
    trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
    trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
    probTr.CheckCovariance();
  }    
  // find first "active layer" - start tracking at the first active layer      
  Int_t firstActiveLayer = 0;
  for (Int_t j=0; j<=lastActiveLayer; j++) { 
    layer = (CylLayerK*)fLayers.At(j);
    if (!(layer->isDead)) { // is alive
      firstActiveLayer = j;
      break;
    }
  }
  //probTr.Rotate(0);
  for (Int_t j=0; j<=lastActiveLayer; j++) {  // Layer loop
    //
    layer = (CylLayerK*)fLayers.At(j);
    TString name(layer->GetName());
    Bool_t isVertex = name.Contains("vertex");
    if (!PropagateToR(&probTr, layer->radius,bGauss,1)) exit(1);
    //
    if (!isVertex) {
      // rotate to frame with X axis normal to the surface
      double pos[3];
      probTr.GetXYZ(pos);  // lab position
      double phi = TMath::ATan2(pos[1],pos[0]);
      if ( TMath::Abs(TMath::Abs(phi)-TMath::Pi()/2)<1e-3) phi = 0;//TMath::Sign(TMath::Pi()/2 - 1e-3,phi);
      if (!probTr.Rotate(phi)) {
	printf("Failed to rotate to the frame (phi:%+.3f)of layer at %.2f at XYZ: %+.3f %+.3f %+.3f (pt=%+.3f)\n",
	       phi,layer->radius,pos[0],pos[1],pos[2],pt);	      
	probTr.Print();
	exit(1);
      }
    }
    //
    // save outward parameters at this layer: before the update
    new( saveParOutwardB[j] ) AliExternalTrackParam(probTr);
    //
    // combined in-out prediction
    new( saveParComb[j]  ) AliExternalTrackParam(*(AliExternalTrackParam*)saveParInward[j]);
    double *covInw = (double*) ((AliExternalTrackParam*)saveParInward[j])->GetCovariance();
    double *covOut = (double*) probTr.GetCovariance();
    double *covCmb = (double*) ((AliExternalTrackParam*)saveParComb[j])->GetCovariance();
    covCmb[0] = covInw[0]*covOut[0]/(covInw[0]+covOut[0]);
    covCmb[2] = covInw[2]*covOut[2]/(covInw[2]+covOut[2]);
    covCmb[1] = 0;
    // create fake measurement with the errors assigned to the layer
    // account for the measurement there
    if (!isVertex && !layer->isDead) {
      double meas[2] = {probTr.GetY(),probTr.GetZ()};
      double measErr2[3] = {layer->phiRes*layer->phiRes,0,layer->zRes*layer->zRes};
      //
      if (!probTr.Update(meas,measErr2)) {
	printf("Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}\n",
	       meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]);
	probTr.Print();
	exit(1);
      }
    }
    // note: if apart from MS we want also e.loss correction, the density*length should be provided as 2nd param
    if (layer->radL>0 && !probTr.CorrectForMeanMaterial(layer->radL, 0, mass , kTRUE)) {
      printf("Failed to apply material correction, X/X0=%.4f\n",layer->radL);
      probTr.Print();
      exit(1);
    }
    // save outward parameters at this layer: after the update
    new( saveParOutwardA[j] ) AliExternalTrackParam(probTr);
    //
  }
  //
  probTr.SetUseLogTermMS(kFALSE); // Reset of MS term usage to avoid problems since its static
  //  
  return kTRUE;
}

Bool_t DetectorK::CalcITSEff(TrackSol& ts, Bool_t verbose)
{
  // Prepare Probability Kombinations
  Int_t nLayer = fNumberOfActiveITSLayers;
  Int_t base = 3; // null, fake, correct
  Int_t komb = (Int_t) TMath::Power(base,nLayer);
  TMatrixD probLayInw(base,fNumberOfActiveITSLayers);
  TMatrixD probLayOut(base,fNumberOfActiveITSLayers);
  TMatrixD probLayCmb(base,fNumberOfActiveITSLayers);
  TMatrixD probKomb(komb,nLayer);
  for (Int_t num=0; num<komb; num++) {
    for (Int_t l=nLayer; l--;) {
      Int_t pow = ((Int_t)TMath::Power(base,l+1));
      probKomb(num,nLayer-1-l)=(num%pow)/((Int_t)TMath::Power(base,l));
    }
  }
  int nITSAct=0, ilr=0;
  if (verbose) printf("Lr:  \t rad   x/x0   h.dens | Inw sY sZ  ->  Pr.Corr | Out sY sZ  ->  Pr.Corr | Cmb sY sZ  ->  Pr.Corr |\n");

  while (nITSAct<nLayer) {
    CylLayerK *l = (CylLayerK*) fLayers.At(ilr);
    TString name(l->GetName());
    if (l->isDead || !IsITSLayer(name) ) {ilr++; continue;}
    //
    AliExternalTrackParam* trInw = (AliExternalTrackParam*)ts.fTrackInw[ilr];
    AliExternalTrackParam* trOut = (AliExternalTrackParam*)ts.fTrackOutB[ilr];
    AliExternalTrackParam* trCmb = (AliExternalTrackParam*)ts.fTrackCmb[ilr];
    //
    double sigYInw = TMath::Sqrt(trInw->GetSigmaY2()+l->phiRes*l->phiRes);
    double sigZInw = TMath::Sqrt(trInw->GetSigmaZ2()+l->zRes*l->zRes);
    probLayInw(2,nITSAct) = ProbGoodChiSqPlusConfHit(l->radius,l->eff, sigYInw, sigZInw);// corr hit prob 						     
    probLayInw(0,nITSAct) = ProbNullChiSqPlusConfHit(l->radius,l->eff, sigYInw, sigZInw); // no hit prob 
    probLayInw(1,nITSAct) = 1.-probLayInw(2,nITSAct)-probLayInw(0,nITSAct);
    //
    double sigYOut = TMath::Sqrt(trOut->GetSigmaY2()+l->phiRes*l->phiRes);
    double sigZOut = TMath::Sqrt(trOut->GetSigmaZ2()+l->zRes*l->zRes);
    probLayOut(2,nITSAct) = ProbGoodChiSqPlusConfHit(l->radius,l->eff, sigYOut, sigZOut);// corr hit prob 
    probLayOut(0,nITSAct) = ProbNullChiSqPlusConfHit(l->radius,l->eff, sigYOut, sigZOut);// no hit prob 
    probLayOut(1,nITSAct) = 1.-probLayOut(2,nITSAct)-probLayOut(0,nITSAct);
    //
    double sigYCmb = TMath::Sqrt(trCmb->GetSigmaY2()+l->phiRes*l->phiRes);
    double sigZCmb = TMath::Sqrt(trCmb->GetSigmaZ2()+l->zRes*l->zRes);
    probLayCmb(2,nITSAct) = ProbGoodChiSqPlusConfHit(l->radius,l->eff, sigYCmb, sigZCmb); // corr hit prob 
    probLayCmb(0,nITSAct) = ProbNullChiSqPlusConfHit(l->radius,l->eff, sigYCmb, sigZCmb); // no hit prob 
    probLayCmb(1,nITSAct) = 1.-probLayCmb(2,nITSAct)-probLayCmb(0,nITSAct);
    //
    if (verbose) {
      const double kCnv=1e4;
      printf("%s:\t%5.1f %.4f %7.0f | %6.0f %6.0f -> %.3f | %6.0f %6.0f -> %.3f | %6.0f %6.0f -> %.3f\n",
	     l->GetName(),l->radius,l->radL,HitDensity(l->radius),
	     sigYInw*kCnv,sigZInw*kCnv,probLayInw(2,nITSAct),
	     sigYOut*kCnv,sigZOut*kCnv,probLayOut(2,nITSAct),
	     sigYCmb*kCnv,sigZCmb*kCnv,probLayCmb(2,nITSAct));
    }
    nITSAct++;
    ilr++;
    //
  }
  PrepareEffFakeKombinations(&probKomb,&probLayInw,(double*)ts.fProb[TrackSol::kInw]);
  PrepareEffFakeKombinations(&probKomb,&probLayOut,(double*)ts.fProb[TrackSol::kOut]);
  PrepareEffFakeKombinations(&probKomb,&probLayCmb,(double*)ts.fProb[TrackSol::kCmb]);
  if (verbose) {
    printf("Corr/Fake probs:             |    %.4f/%.4f       |     %.4f/%.4f      |     %.4f/%.4f\n",
	   ts.fProb[TrackSol::kInw][0],ts.fProb[TrackSol::kInw][1],
	   ts.fProb[TrackSol::kOut][0],ts.fProb[TrackSol::kOut][1],
	   ts.fProb[TrackSol::kCmb][0],ts.fProb[TrackSol::kCmb][1]);
  }
  //
  return kTRUE;
}

//_____________________________________________________________________
Bool_t DetectorK::ExtrapolateToR(AliExternalTrackParam* probTr, double rTgt, double mass)
{
  // propagate the track to given radius R without updates (final extrapolation in the tracking frame of cyl. at R)
  double xCurr = probTr->GetX();
  double rCurr = TMath::Sqrt(xCurr*xCurr + probTr->GetY()*probTr->GetY());
  //
  if (TMath::Abs(rCurr-rTgt)<1e-6) return kTRUE;
  //
  int dir = rCurr>rTgt ? -1 : 1; // inward or outward?
  // 
  // detemine current layer
  int lrCurr=-1,lrTgt=-1;
  //
  if (dir<0) { // inward
    for (int ilr=fNumberOfLayers;ilr--;) {
      CylLayerK *l = (CylLayerK*)fLayers.At(ilr);
      if (lrCurr<0 && l->radius<=rCurr) lrCurr = ilr;
      if (l->radius>=rTgt)  lrTgt  = ilr;
    }
  }
  else {
    for (int ilr=0;ilr<fNumberOfLayers;ilr++) {
      CylLayerK *l = (CylLayerK*)fLayers.At(ilr);
      if (lrCurr<0 && l->radius>=rCurr) lrCurr = ilr;
      if (l->radius<rTgt)  lrTgt  = ilr;
    }
  }
  //
  /*
  printf("Xcurr: %.2f Xdest: %.2f %d Icurr:%d Itgt:%d Rcurr:%.2f Rdest:%.2f\n",
	 xCurr,rTgt, dir, lrCurr, lrTgt, 
	 ((CylLayerK*)fLayers.At(lrCurr))->radius, ((CylLayerK*)fLayers.At(lrTgt))->radius );
  */
  //
  double bGauss = fBField*10;
  //
  for (int ilr=lrCurr;ilr!=lrTgt;ilr+=dir) {
    CylLayerK *l = (CylLayerK*)fLayers.At(ilr);
    if (!PropagateToR(probTr,l->radius,bGauss,dir)) return kFALSE;
    if (!probTr->CorrectForMeanMaterial(l->radL, 0, mass , kTRUE)) return kFALSE;
    if (verboseR) {
      printf("\nGot to layer %d | ",ilr); probTr->Print();
    }
  }
  // go to destination r
  if (!PropagateToR(probTr,rTgt,bGauss,dir)) return kFALSE;
  if (verboseR) {
    printf("\nGot to r= %.2f | ",rTgt); probTr->Print();
  }

  //
  return kTRUE;
}


TGraph * DetectorK::GetGraphMomentumResolution(Int_t color, Int_t linewidth) {
  //
  // returns the momentum resolution 
  //
  
  TGraph *graph = new TGraph(kNptBins, fTransMomenta, fMomentumRes);
  graph->SetTitle("Momentum Resolution .vs. Pt" ) ;
  //  graph->GetXaxis()->SetRangeUser(0.,5.0) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Momentum Resolution (%)") ;
  graph->GetYaxis()->CenterTitle();

  graph->SetMaximum(20);
  graph->SetMinimum(0.1);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph * DetectorK::GetGraphPointingResolution(Int_t axis, Int_t color, Int_t linewidth) {
 
  // Returns the pointing resolution
  // axis = 0 ... rphi pointing resolution
  // axis = 1 ... z pointing resolution
  //

  TGraph * graph =  0;

  if (axis==0) {
    graph = new TGraph ( kNptBins, fTransMomenta, fResolutionRPhi ) ;
    graph->SetTitle("R-#phi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("R-#phi Pointing Resolution (#mum)") ;
  } else {
    graph =  new TGraph ( kNptBins, fTransMomenta, fResolutionZ ) ;
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution (#mum)") ;
  }
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(1000.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineWidth(linewidth);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  
  return graph;

}

TGraph * DetectorK::GetGraphLayerInfo(Int_t plot, Int_t color, Int_t linewidth) {
 
  // Returns the pointing resolution
  // plot = 0 ... rphi pointing resolution
  // plot = 1 ... z pointing resolution
  // plot = 2 ... prolongation efficiency (outwards)
  //

   
  Double_t fDet[kNptBins]; 
  for ( Int_t i = 0 ; i < kNptBins ; i++ ) { // pt loop
    if (plot==0) 
      fDet[i] = fResolutionRPhiLay[i]*1e6; // in microns
    else if (plot==1)
      fDet[i] = fResolutionZLay[i]*1e6;    // in microns
    else 
      fDet[i] = fEfficProlongLay[i]*100;   // in percent
  }
 
  CylLayerK *l = (CylLayerK*) fLayers.At(kDetLayer);
  TGraph * graph =  0;
  graph =  new TGraph ( kNptBins, fTransMomenta, fDet ) ;
  if (plot==0) {
    graph->SetTitle(Form("R-#phi Pointing Resolution onto layer \"%s\"",(char*)l->GetName()) );
    graph->GetYaxis()->SetTitle("R-#phi Pointing Resolution (#mum)") ;
  } else if (plot==1){
    graph->SetTitle(Form("Z Pointing Resolution onto layer \"%s\"",(char*)l->GetName()) ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution (#mum)") ;
  } else {
    graph->SetTitle(Form("Prolongation efficiency onto layer \"%s\"",(char*)l->GetName()) ) ;
    graph->GetYaxis()->SetTitle("Prolongation efficiency (%)") ;
    graph->SetMinimum(0);
    graph->SetMaximum(100);
  }
  
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineWidth(linewidth);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  
  return graph;

}



TGraph * DetectorK::GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth) {
  //
  // returns the Pointing resolution (accoring to Telescope equation)
  // axis =0 ... in rphi
  // axis =1 ... in z
  //
  
  Double_t resolution[kNptBins];

  Double_t layerResolution[2];
  Double_t layerRadius[2];
  Double_t layerThickness[2];

  Int_t count =0; // search two first active layers
  printf("Telescope equation for layers:  ");
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    CylLayerK *l = (CylLayerK*)fLayers.At(i);
    if (!l->isDead && l->radius>0) {
      layerRadius[count]     = l->radius;
      layerThickness[count]  = l->radL;
      if (axis==0) {
	layerResolution[count] = l->phiRes;
      } else {
	layerResolution[count] = l->zRes;
      }
      printf("%s, ",l->GetName());
      count++;
    }
    if (count>=2) break;	
  }
  printf("\n");

  Double_t pt, momentum, thickness,aMCS ;
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 

  for ( Int_t i = 0 ; i < kNptBins ; i++ ) { 
    // Reference data as if first two layers were acting all alone 
    pt  =  fTransMomenta[i]  ;
    momentum = pt / TMath::Cos(lambda)   ;  // Total momentum
    resolution[i] =  layerResolution[0]*layerResolution[0]*layerRadius[1]*layerRadius[1] 
      +  layerResolution[1]*layerResolution[1]*layerRadius[0]*layerRadius[0] ;
    resolution[i] /= ( layerRadius[1] - layerRadius[0] ) * ( layerRadius[1] - layerRadius[0] ) ;
    thickness = layerThickness[0] / TMath::Sin(TMath::Pi()/2 - lambda) ;
    aMCS = ThetaMCS(fParticleMass, thickness, momentum) ;
    resolution[i] += layerRadius[0]*layerRadius[0]*aMCS*aMCS ;
    resolution[i] =  TMath::Sqrt(resolution[i]) * 10000.0 ;  // result in microns
  }



  TGraph* graph = new TGraph ( kNptBins, fTransMomenta, resolution ) ;
   
  if (axis==0) {
    graph->SetTitle("RPhi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("RPhi Pointing Resolution (#mum) ") ;
  } else {
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution (#mum) ") ;
  }
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineStyle(kDashed);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph * DetectorK::GetGraphRecoEfficiency(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  // particle = 3 ... D0
  //
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleEfficiency[kNptBins]; // with chosen particle mass
  Double_t kaonEfficiency[kNptBins]={0}, pionEfficiency[kNptBins]={0}, d0efficiency[kNptBins]={0}; 
  Double_t partEfficiency[2][400]={{0}};
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
	// JT Test Let the kaon decay.  If it decays inside the TPC ... then it is gone; for all decays < 130 cm.
	Double_t momentum = fTransMomenta[j] / TMath::Cos(lambda)           ;  // Total momentum at average rapidity
	if ( massloop == 1 ) { // KAON
	  doNotDecayFactor  = TMath::Exp(-130/(371*momentum/KaonMass)) ;  // Decay length for kaon is 371 cm.
	  kaonEfficiency[j] = fEfficiency[1][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;
	} else { // PION
	  doNotDecayFactor = 1.0 ;
	  pionEfficiency[j] = fEfficiency[0][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;	
	}
	partEfficiency[0][j] = pionEfficiency[j];
	partEfficiency[1][j] = kaonEfficiency[j];
      }      
    }
    
    // resulting estimate of the D0 efficiency
    for ( Int_t j = 0 ; j < kNptBins ; j++ ) {
      d0efficiency[j] = D0IntegratedEfficiency(fTransMomenta[j],partEfficiency);
    }
  } else { 
    for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
      particleEfficiency[j] = fEfficiency[2][j]* AcceptanceOfTpcAndSi;
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
    pionEfficiency[j]     *= 100;
    kaonEfficiency[j]     *= 100;
    d0efficiency[j]       *= 100;
    particleEfficiency[j] *= 100;
  }
 
  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( kNptBins, fTransMomenta, particleEfficiency ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( kNptBins, fTransMomenta, pionEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( kNptBins, fTransMomenta, kaonEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==3) {
    graph = new TGraph ( kNptBins, fTransMomenta, d0efficiency ) ;
    graph->SetLineStyle(kDashed);
  } else 
    return 0;

  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Efficiency (%)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(100)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}

TGraph * DetectorK::GetGraphRecoFakes(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  //

  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleFake[kNptBins]; // with chosen particle mass
  Double_t kaonFake[kNptBins], pionFake[kNptBins];
  Double_t partFake[2][kNptBins];
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
	// JT Test Let the kaon decay.  If it decays inside the TPC ... then it is gone; for all decays < 130 cm.
	Double_t momentum = fTransMomenta[j] / TMath::Cos(lambda)           ;  // Total momentum at average rapidity
	if ( massloop == 1 ) { // KAON
	  doNotDecayFactor  = TMath::Exp(-130/(371*momentum/KaonMass)) ;  // Decay length for kaon is 371 cm.
	  kaonFake[j] = fFake[1][j] /( doNotDecayFactor) ;
	} else { // PION
	  pionFake[j] = fFake[0][j] ;	
	}
	partFake[0][j] = pionFake[j];
	partFake[1][j] = kaonFake[j];
      }      
    }
    
  } else { 
    for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
      particleFake[j] = fFake[2][j];
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
    pionFake[j]     *= 100;
    kaonFake[j]     *= 100;
    particleFake[j] *= 100;
  }
 
  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( kNptBins, fTransMomenta, particleFake ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( kNptBins, fTransMomenta, pionFake ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( kNptBins, fTransMomenta, kaonFake ) ;
    graph->SetLineWidth(1);
  } 
  
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Fake (%)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(100)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}
TGraph * DetectorK::GetGraphRecoPurity(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  //

  //  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleFake[kNptBins]; // with chosen particle mass
  Double_t kaonFake[kNptBins], pionFake[kNptBins];
  //  Double_t partFake[2][kNptBins];
  
  if (particle != 0) {
    cout <<" not implemented"<<endl;
      
  } else { 
    for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
      particleFake[j] = fFake[2][j];
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  // Get Purity
  for ( Int_t j = 0 ; j < kNptBins ; j++ ) { 
    pionFake[j]     = (1-pionFake[j])*100;
    kaonFake[j]     = (1-kaonFake[j])*100;
    particleFake[j] = (1-particleFake[j])*100;
  }
 
  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( kNptBins, fTransMomenta, particleFake ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( kNptBins, fTransMomenta, pionFake ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( kNptBins, fTransMomenta, kaonFake ) ;
    graph->SetLineWidth(1);
  } 
  
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Purity (%)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(100)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}


TGraph* DetectorK::GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth) {
  //
  // returns the Impact Parameter d0 (convolution of pointing resolution and vtx resolution)
  // mode 0: impact parameter (convolution of pointing and vertex resolution)
  // mode 1: pointing resolution
  // mode 2: vtx resolution 
  
  
  TGraph *graph = new TGraph();

  //  TFormula vtxResRPhi("vtxRes","50-2*x"); // 50 microns at pt=0, 15 microns at pt =20 ?
  TFormula vtxResRPhi("vtxRes","35/(x+1)+10"); // 
  TFormula vtxResZ("vtxResZ","600/(x+6)+10"); // 
    
  TGraph *trackRes = GetGraphPointingResolution(axis,1);
  Double_t *pt = trackRes->GetX();
  Double_t *trRes = trackRes->GetY();
  for (Int_t ip =0; ip<trackRes->GetN(); ip++) {
    Double_t vtxRes = 0;
    if (axis==0) 
      vtxRes = vtxResRPhi.Eval(pt[ip]);
    else 
      vtxRes = vtxResZ.Eval(pt[ip]);
    
    if (mode==0)
      graph->SetPoint(ip,pt[ip],TMath::Sqrt(vtxRes*vtxRes+trRes[ip]*trRes[ip]));
    else if (mode ==1)
      graph->SetPoint(ip,pt[ip],trRes[ip]);
    else
      graph->SetPoint(ip,pt[ip],vtxRes);
  }
  
  graph->SetTitle("d_{0} r#phi resolution .vs. Pt" ) ;
  graph->GetYaxis()->SetTitle("d_{0} r#phi resolution (#mum)") ;
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph* DetectorK::GetGraph(Int_t number, Int_t color, Int_t linewidth) {
  // 
  // returns graph according to the number
  //
  switch(number) {
  case 1:
    return GetGraphPointingResolution(0,color, linewidth); // dr
  case 2:
    return GetGraphPointingResolution(1,color, linewidth); // dz
  case 3:
    return GetGraphPointingResolutionTeleEqu(0,color, linewidth); // dr - tele
  case 4:
    return GetGraphPointingResolutionTeleEqu(1,color, linewidth); // dz - tele
  case 5:
    return GetGraphMomentumResolution(color, linewidth); // pt resolution
  case 10:
    return GetGraphRecoEfficiency(0, color, linewidth);  // tracked particle
  case 11:
    return GetGraphRecoEfficiency(1, color, linewidth);  // eff. pion
  case 12:
    return GetGraphRecoEfficiency(2, color, linewidth);  // eff. kaon
  case 13: 
    return GetGraphRecoEfficiency(3, color, linewidth);  // eff. D0
  case 15:
    return GetGraphRecoFakes(0, color, linewidth);  // Fake tracked particle
  case 16:
    return GetGraphRecoFakes(1, color, linewidth);  // Fake pion
  case 17:
    return GetGraphRecoFakes(2, color, linewidth);  // Fake kaon
  default:
    printf(" Error: chosen graph number not valid\n");
  }
  return 0;

}

void DetectorK::MakeAliceAllNew(Bool_t flagTPC,Bool_t flagMon) {
  
  // All New configuration with X0 = 0.3 and resolution = 4 microns
  
  AddLayer((char*)"bpipe",2.0,0.0022); // beam pipe
  AddLayer((char*)"vertex",     0,     0); // dummy vertex for matrix calculation

  // new ideal Pixel properties?
  Double_t x0     = 0.0050;
  Double_t resRPhi = 0.0006;
  Double_t resZ   = 0.0006;
 
  if (flagMon) {
    x0     = 0.0030;
    resRPhi = 0.0004;
    resZ   = 0.0004;
  }
  
  AddLayer((char*)"ddd1",  2.2 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd2",  2.8 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd3",  3.6 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd4", 20.0 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd5", 22.0 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd6", 41.0 ,  x0, resRPhi, resZ); 
  AddLayer((char*)"ddd7", 43.0 ,  x0, resRPhi, resZ); 
 
  if (flagTPC) {
    AddTPC(0.1,0.1);                        // TPC
  }
}

void DetectorK::MakeAliceCurrent(Int_t AlignResiduals, Bool_t flagTPC) {

  // Numbers taken from 
  // 2010 JINST 5 P03003 - Alignment of the ALICE Inner Tracking System with cosmic-ray tracks
  // number for misalingment: private communication with Andrea Dainese

  AddLayer((char*)"bpipe",2.94,0.0022); // beam pipe
  AddLayer((char*)"vertex",     0,     0); // dummy vertex for matrix calculation
  AddLayer((char*)"tshld1",11.5,0.0065); // Thermal shield  // 1.3% /2
  AddLayer((char*)"tshld2",31.0,0.0065); // Thermal shield  // 1.3% /2


  if (flagTPC) {
    AddTPC(0.1,0.1);                        // TPC
  }
  // Adding the ITS - current configuration
  
  if (AlignResiduals==0) {

    AddLayer((char*)"spd1", 3.9, 0.0114, 0.0012, 0.0130);
    AddLayer((char*)"spd2", 7.6, 0.0114, 0.0012, 0.0130);
    AddLayer((char*)"sdd1",15.0, 0.0113, 0.0035, 0.0025);
    AddLayer((char*)"sdd2",23.9, 0.0126, 0.0035, 0.0025);
    AddLayer((char*)"ssd1",38.0, 0.0083, 0.0020, 0.0830);
    AddLayer((char*)"ssd2",43.0, 0.0086, 0.0020, 0.0830);

  } else if (AlignResiduals==1) {

    // tracking errors ...
    // (Additional systematic errors due to misalignments) ... 
    // itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020);  // [cm]
    // itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000);

    AddLayer((char*)"spd1", 3.9, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0010*0.0010), 
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030),
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"sdd1",15.0, 0.0113, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"sdd2",23.9, 0.0126, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"ssd1",38.0, 0.0083, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020), 
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));
    AddLayer((char*)"ssd2",43.0, 0.0086, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020),
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));   
    
  } else if (AlignResiduals==2) {
    
    // tracking errors ... PLUS ... module misalignment
    
    // itsRecoParam->SetClusterMisalErrorYBOn(0.0010,0.0030,0.0500,0.0500,0.0020,0.0020);  // [cm]
    // itsRecoParam->SetClusterMisalErrorZBOn(0.0050,0.0050,0.0050,0.0050,0.1000,0.1000);
    
    //  the ITS modules are misalignment with small gaussian smearings with
    //  sigmarphi ~ 8, 10, 10 micron in SPD, SDD, SSD
    
    AddLayer((char*)"spd1", 3.9, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0010*0.0010+0.0008*0.0008), 
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030+0.0008*0.0008),
	     TMath::Sqrt(0.0130*0.0130+0.0050*0.0050));
    AddLayer((char*)"sdd1",15.0, 0.0113, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"sdd2",23.9, 0.0126, TMath::Sqrt(0.0035*0.0035+0.0500*0.0500+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.0050*0.0050));
    AddLayer((char*)"ssd1",38.0, 0.0083, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020+0.0010*0.0010), 
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000));
    AddLayer((char*)"ssd2",43.0, 0.0086, TMath::Sqrt(0.0020*0.0020+0.0020*0.0020+0.0010*0.0010),
	     TMath::Sqrt(0.0830*0.0830+0.1000*0.1000)); 

  } else {
      
      //  the ITS modules are misalignment with small gaussian smearings with
      //  sigmarphi ~ 8, 10, 10 micron in SPD, SDD, SSD
      //  unknown in Z ????

    AddLayer((char*)"spd1", 3.9, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0008*0.0008), 
	     TMath::Sqrt(0.0130*0.0130+0.000*0.000));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0008*0.0008),
	     TMath::Sqrt(0.0130*0.0130+0.000*0.000));
    AddLayer((char*)"sdd1",15.0, 0.0113, TMath::Sqrt(0.0035*0.0035+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.000*0.000));
    AddLayer((char*)"sdd2",23.9, 0.0126, TMath::Sqrt(0.0035*0.0035+0.0010*0.0010),
	     TMath::Sqrt(0.0025*0.0025+0.000*0.000));
    AddLayer((char*)"ssd1",38.0, 0.0083, TMath::Sqrt(0.0020*0.0020+0.0010*0.0010), 
	     TMath::Sqrt(0.0830*0.0830+0.000*0.000));
    AddLayer((char*)"ssd2",43.0, 0.0086, TMath::Sqrt(0.0020*0.0020+0.0010*0.0010),
	     TMath::Sqrt(0.0830*0.0830+0.000*0.000));   
    
    
  }
  
}


void DetectorK::MakeStandardPlots(Bool_t add, Int_t color, Int_t linewidth,Bool_t onlyPionEff) {
  //
  // Produces the standard performace plots
  //
  TGraph *eff,*momRes,*pointRes;
  if (!add) {

    TCanvas *c1 = new TCanvas("c1","c1");//,100,100,500,500);  
    c1->Divide(2,2);
    
    c1->cd(1);  gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    eff = GetGraphRecoEfficiency(0,color,linewidth);
    eff->SetName(Form("grEff%d",1));
    eff->SetTitle("Efficiencies");
    eff->Draw("AL");
    eff->SetMaximum(110);
    if (!onlyPionEff) {
      eff = GetGraphRecoEfficiency(2,color,linewidth);
      eff->SetName(Form("grEff%d",2));
      eff->Draw("L");
      eff = GetGraphRecoEfficiency(3,color,linewidth);
      eff->SetName(Form("grEff%d",3));
      eff->Draw("L");
    }
    c1->cd(2); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogy();  gPad->SetLogx(); 
    momRes = GetGraphMomentumResolution(color,linewidth);
    momRes->SetName(Form("grMomRes%d",1));
    momRes->Draw("AL");
    
    c1->cd(3); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    pointRes = GetGraphPointingResolution(0,color,linewidth);
    pointRes->SetName(Form("pointRRes%d",0));
    pointRes->Draw("AL");
    //
    c1->cd(4); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    pointRes = GetGraphPointingResolution(1,color,linewidth);
    pointRes->SetName(Form("pointZRes%d",0));    
    pointRes->Draw("AL");
    
  } else {
    
    TVirtualPad *c1 = gPad->GetMother();

    c1->cd(1);
    eff = GetGraphRecoEfficiency(0,color,linewidth);
    eff->SetName(Form("grEff%dadd",1));
    eff->Draw("L");
    if (!onlyPionEff) {
      eff = GetGraphRecoEfficiency(2,color,linewidth);
      eff->SetName(Form("grEff%dadd",2));
      eff->Draw("L");
      eff = GetGraphRecoEfficiency(3,color,linewidth);
      eff->SetName(Form("grEff%dadd",3));
      eff->Draw("L");
    }
    c1->cd(2); 
    momRes = GetGraphMomentumResolution(color,linewidth);
    momRes->SetName(Form("grMomRes%dadd",1));
    momRes->Draw("L");
    
    c1->cd(3); 
    pointRes = GetGraphPointingResolution(0,color,linewidth);
    pointRes->SetName(Form("pointRRes%dadd",0));
    pointRes->Draw("L");
    
    c1->cd(4); 
    pointRes = GetGraphPointingResolution(1,color,linewidth);
    pointRes->SetName(Form("pointZRes%dadd",0));
    pointRes->Draw("L");
    
  }

}


Bool_t DetectorK::GetXatLabR(AliExternalTrackParam* tr,Double_t r,Double_t &x, Double_t bz, Int_t dir)
{
  // Get local X of the track position estimated at the radius lab radius r. 
  // The track curvature is accounted exactly
  //
  // The flag "dir" can be used to remove the ambiguity of which intersection to take (out of 2 possible)
  // 0  - take the intersection closest to the current track position
  // >0 - go along the track (increasing R)
  // <0 - go backward (decreasing R)
  //
  // special case of R=0
  if (r<kAlmost0) {x=0; return kTRUE;}
  
  const double* pars = tr->GetParameter();
  const Double_t &fy=pars[0], &sn = pars[2];
  const double kEps = 1.e-6;
  //
  double fx = tr->GetX(); 
  double crv = tr->GetC(bz);
  if (TMath::Abs(crv)>kAlmost0) {                                 // helix
    // get center of the track circle
    double tR = 1./crv;   // track radius (for the moment signed)
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double r0 = TMath::Sqrt(x0*x0+y0*y0);
    //    printf("Xc:%+e Yc:%+e tR:%e r0:%e\n",x0,y0,tR,r0);
    //
    if (r0<=kAlmost0) return kFALSE;            // the track is concentric to circle
    tR = TMath::Abs(tR);
    double tR2r0=1.,g=0,tmp=0;
    if (TMath::Abs(tR-r0)>kEps) {
      tR2r0 = tR/r0;
      g = 0.5*(r*r/(r0*tR) - tR2r0 - 1./tR2r0);
      tmp = 1.+g*tR2r0;
    }
    else {
      tR2r0 = 1.0;
      g = 0.5*r*r/(r0*tR) - 1;
      tmp = 0.5*r*r/(r0*r0);
    }
    double det = (1.-g)*(1.+g);
    if (det<0) return kFALSE;         // does not reach raduis r
    det = TMath::Sqrt(det);    
    //
    // the intersection happens in 2 points: {x0+tR*C,y0+tR*S} 
    // with C=f*c0+-|s0|*det and S=f*s0-+c0 sign(s0)*det
    // where s0 and c0 make direction for the circle center (=x0/r0 and y0/r0)
    //
    x = x0*tmp; 
    double y = y0*tmp;
    if (TMath::Abs(y0)>kAlmost0) { // when y0==0 the x,y is unique
      double dfx = tR2r0*TMath::Abs(y0)*det;
      double dfy = tR2r0*x0*TMath::Sign(det,y0);
      if (dir==0) {                    // chose the one which corresponds to smallest step 
	double delta = (x-fx)*dfx-(y-fy)*dfy; // the choice of + in C will lead to smaller step if delta<0
	if (delta<0) x += dfx;
	else         x -= dfx;
      }
      else if (dir>0) {  // along track direction: x must be > fx
	x -= dfx; // (dfx is positive)
	double dfeps = fx-x; // handle special case of very small step
	//	if (verboseR) printf("d+: x:%+e|%+e fx: %+e d %+e\n",x,x+dfx+dfx,fx,dfeps);
	if (dfeps<-kEps) return kTRUE;
	if (TMath::Abs(dfeps)<kEps &&  // are we already in right r?
	    TMath::Abs(fx*fx+fy*fy - r*r)<kEps) return fx;
	x += dfx+dfx;
	if (x-fx>0) return kTRUE;
	if (x-fx<-kEps) return kFALSE;
	x = fx; // don't move
      }
      else { // backward: x must be < fx
	x += dfx; // try the smallest step (dfx is positive)	
	double dfeps = x-fx; // handle special case of very small step
	//	if (verboseR) printf("d-: x:%+e|%+e fx: %+e d %+e\n",x,x-dfx-dfx,fx,dfeps);
	if (dfeps<-kEps) return kTRUE;
	if (TMath::Abs(dfeps)<kEps &&  // are we already in right r?
	    TMath::Abs(fx*fx+fy*fy - r*r)<kEps) return fx;
	x-=dfx+dfx;
	if (x-fx<0) return kTRUE;
	if (x-fx>kEps) return kFALSE;
	x = fx; // don't move
      }
    }
    else { // special case: track touching the circle just in 1 point
      if ( (dir>0&&x<fx) || (dir<0&&x>fx) ) return kFALSE; 
    }
  }
  else { // this is a straight track
    if (TMath::Abs(sn)>=kAlmost1) { // || to Y axis
      double det = (r-fx)*(r+fx);
      if (det<0) return kFALSE;     // does not reach raduis r
      x = fx;
      if (dir==0) return kTRUE;
      det = TMath::Sqrt(det);
      if (dir>0) {                       // along the track direction
	if (sn>0) {if (fy>det)  return kFALSE;} // track is along Y axis and above the circle
	else      {if (fy<-det) return kFALSE;} // track is against Y axis amd belo the circle
      }
      else if(dir>0) {                                    // agains track direction
	if (sn>0) {if (fy<-det) return kFALSE;} // track is along Y axis
        else if (fy>det)  return kFALSE;        // track is against Y axis
      }
    }
    else if (TMath::Abs(sn)<=kAlmost0) { // || to X axis
      double det = (r-fy)*(r+fy);
      if (det<0) return kFALSE;     // does not reach raduis r
      det = TMath::Sqrt(det);
      if (!dir) {
	x = fx>0  ? det : -det;    // choose the solution requiring the smalest step
	return kTRUE;
      }
      else if (dir>0) {                    // along the track direction
	if      (fx > det) return kFALSE;  // current point is in on the right from the circle
	else if (fx <-det) x = -det;       // on the left
	else               x =  det;       // within the circle
      }
      else {                               // against the track direction
	if      (fx <-det) return kFALSE;  
	else if (fx > det) x =  det;
	else               x = -det;
      }
    }
    else {                                 // general case of straight line
      double cs = TMath::Sqrt((1-sn)*(1+sn));
      double xsyc = fx*sn-fy*cs;
      double det = (r-xsyc)*(r+xsyc);
      if (det<0) return kFALSE;    // does not reach raduis r
      det = TMath::Sqrt(det);
      double xcys = fx*cs+fy*sn;
      double t = -xcys;
      if (dir==0) t += t>0 ? -det:det;  // chose the solution requiring the smalest step
      else if (dir>0) {                 // go in increasing fX direction. ( t+-det > 0)
	if (t>=-det) t += -det;         // take minimal step giving t>0
	else return kFALSE;             // both solutions have negative t
      }
      else {                            // go in increasing fx direction. (t+-det < 0)
	if (t<det) t -= det;            // take minimal step giving t<0
	else return kFALSE;             // both solutions have positive t
      }
      x = fx + cs*t;
    }
  }
  //
  return kTRUE;
}



Double_t* DetectorK::PrepareEffFakeKombinations(TMatrixD *probKomb, TMatrixD *probLay, double *probs) {

  if (!probLay) {  
    printf("Error: Layer tracking efficiencies not set \n");
    return 0;
  }

  TMatrixD &tProbKomb = *probKomb;
  TMatrixD &tProbLay = *probLay;


  //  Int_t base = tProbLay.GetNcols(); // 3? null, fake, correct
  Int_t nLayer = tProbKomb.GetNcols(); // nlayer? - number of ITS layers
  Int_t komb = tProbKomb.GetNrows(); // 3^nlayer? - number of kombinations

  // Fill probabilities 

  Double_t probEff =0;
  Double_t probFake =0;
  for (Int_t num=0; num<komb; num++) {
    Int_t flCorr=0, flFake=0, flNull=0; 
     for (Int_t l=0; l<nLayer; l++)  {
      if (tProbKomb(num,l)==0) 
	flNull++;
      else if (tProbKomb(num,l)==1)
	flFake++;
      else if (tProbKomb(num,l)==2)
	flCorr++;
      else 
	printf("Error: unexpected values in combinatorics table\n");
    }

    Int_t fkAtLeastHits = fAtLeastHits;
    Int_t fkAtLeastCorr = fAtLeastCorr;
    if (fAtLeastHits == -1) fkAtLeastHits = nLayer; // all hits are "correct"
    if (fAtLeastCorr == -1) fkAtLeastCorr = nLayer; // all hits are "correct"
    //
    if (flCorr+flFake < fAtLeastHits) continue;

    if (flCorr>=fkAtLeastCorr && flFake==0) { // at least correct but zero fake
      Double_t probEffLayer = 1;
      for (Int_t l=0; l<nLayer; l++) {
	probEffLayer *=  tProbLay((Int_t)tProbKomb(num,l),l);
	//	cout<<a(num,l)<<" ";
      }
      //      cout<<endl;
      probEff+=probEffLayer;
    }
 
    if (flFake>=fAtLeastFake) {
      Double_t probFakeLayer = 1;
      for (Int_t l=0; l<nLayer; l++) {
	probFakeLayer *=  tProbLay((Int_t)tProbKomb(num,l),l);
	//	cout<<a(num,l)<<" ";
      }
      //      cout<<endl;
      probFake+=probFakeLayer;
    }

  }
  if (!probs) probs = new Double_t[2];
  probs[0] = probEff; probs[1] = probFake;
  return probs;

}

//____________________________________
Bool_t DetectorK::PropagateToR(AliExternalTrackParam* trc, double r, double b, int dir) 
{
  // go to radius R
  //
  double xR = 0;
  double rr = r*r;
  int iter = 0;
  const double kTiny = 1e-6;
  //
  if (verboseR) {
    printf("Prop to %f d=%d  ",r,dir); trc->Print();
  }
  while(1) {
    //    if (!trc->GetXatLabR(r,xR,b,dir)) {
    //RRR rotate to local tracking frame
    if (!trc->Rotate(trc->Phi())) {printf("Failed to rotate to local frame %f |",trc->Phi()); trc->Print(); return kFALSE;}


    if (!GetXatLabR(trc, r ,xR, b, dir)) {
      printf("Track with pt=%f cannot reach radius %f\n",trc->Pt(),r);
      trc->Print();
      return kFALSE;
    }
    double snp = trc->GetSnpAt(xR,b);    

    if (!trc->PropagateTo(xR, b)) {printf("Failed to propagate to X=%f for R=%f snp=%f | iter=%d\n",xR,r,snp,iter); trc->Print(); return kFALSE;}
    double rcurr2 = xR*xR + trc->GetY()*trc->GetY();
    if (TMath::Abs(rcurr2-rr)<kTiny || rr<kTiny) break;
    //
    //    printf("new it%d for r=%f (xR=%f) rcurr=%f snp:%f alp:%f\n",iter, r,xR,TMath::Sqrt(rcurr2),trc->GetSnp(),trc->GetAlpha());
    if (++iter>8) {printf("Failed to propagate to R=%f after %d steps\n",r,iter); trc->Print(); return kFALSE;}
    if (verboseR) {
      printf("iter %d ",iter); trc->Print();
    }
  } 
  //
  /*
  // rotate to "sensor" frame (along intersection radius)
  if (r>kTiny) {
    double pos[3]; trc->GetXYZ(pos);
    double phi = TMath::ATan2(pos[1],pos[0]); //TMath::ASin( trc->GetSnp() );
    if (!trc->Rotate(phi)) {printf("Failed to rotate to %f to propagate to R=%f\n",phi,r); trc->Print(); return kFALSE;}
  }
  */
  if (verboseR) {
    printf("iter end "); trc->Print();
  }

  return kTRUE;
}


//_________________________________________
Bool_t DetectorK::IsITSLayer(const TString &lname)
{
  // return true for ITS layers
  return !(lname.Contains("tpc") || lname.Contains("trd"));
}
