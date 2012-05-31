#include "Detector.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFormula.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TText.h>
#include <TGraphErrors.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

original code of using the billoir technique was developed
for the HFT (STAR), James H. Thomas, jhthomas@lbl.gov
http://rnc.lbl.gov/~jhthomas

Changes by S. Rossegger -> see header file

***********************************************************/


#define RIDICULOUS 999999 // A ridiculously large resolution (cm) to flag a dead detector

#define Luminosity    1.e27       // Luminosity of the beam (LHC HI == 1.e27, RHIC II == 8.e27 )
#define SigmaD        6.0         // Size of the interaction diamond (cm) (LHC = 6.0 cm)
#define dNdEtaMinB    950         // Multiplicity per unit Eta  (AuAu MinBias = 170, Central = 700)
#define dNdEtaCent    2300        // Multiplicity per unit Eta  (LHC at 5.5 TeV not known)

#define CrossSectionMinB         8    // minB Cross section for event under study (PbPb MinBias ~ 8 Barns)
#define AcceptanceOfTpcAndSi     1//0.60 //0.35  // Assumed geometric acceptance (efficiency) of the TPC and Si detectors
#define UPCBackgroundMultiplier  1.0   // Increase multiplicity in detector (0.0 to 1.0 * UPCRate ) (eg 1.0)
#define OtherBackground          0.0   // Increase multiplicity in detector (0.0 to 1.0 * minBias)  (eg 0.0)
#define EfficiencySearchFlag     2     // Define search method:
                                       // -> ChiSquarePlusConfLevel = 2, ChiSquare = 1, Simple = 0.  

#define PionMass                 0.139  // Mass of the Pion
#define KaonMass                 0.498  // Mass of the Kaon
#define D0Mass                   1.865  // Mass of the D0



class CylLayer : public TNamed {
public:

  CylLayer(char *name) : TNamed(name,name) {}
  
  Float_t GetRadius()   const {return radius;}
  Float_t GetRadL()     const {return radL;}
  Float_t GetPhiRes()   const {return phiRes;}
  Float_t GetZRes()     const {return zRes;}
  Float_t GetLayerEff() const {return eff;}

  //  void Print() {printf("  r=%3.1lf X0=%1.6lf sigPhi=%1.4lf sigZ=%1.4lf\n",radius,radL,phiRes,zRes); }
  Float_t radius; Float_t radL; Float_t phiRes; Float_t zRes;   
  Float_t eff;
  Bool_t isDead;

 ClassDef(CylLayer,1);
};

ClassImp(Detector)
Detector::Detector() 
  : TNamed("test_detector","detector"),
    fNumberOfLayers(0),
    fNumberOfActiveLayers(0),
    fBField(0.5),
    fLhcUPCscale(2.5),    
    fIntegrationTime(0.02), // in ms
    fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140),    // Standard: pion mass 
    fMaxRadiusSlowDet(10.)
{
  //
  // default constructor
  //
  //  fLayers = new TObjArray();
  
}

Detector::Detector(char *name, char *title)
  : TNamed(name,title),
    fNumberOfLayers(0),
    fNumberOfActiveLayers(0),
    fBField(0.5),
    fLhcUPCscale(2.5),
    fIntegrationTime(0.02),  // in ms
    fConfLevel(0.0027),      // 0.27 % -> 3 sigma confidence
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140),     // Standard: pion mass
    fMaxRadiusSlowDet(10.)
{
  //
  // default constructor, that set the name and title
  //
  //  fLayers = new TObjArray();
}
Detector::~Detector() { // 
  // virtual destructor
  //
  //  delete fLayers;
}

void Detector::AddLayer(char *name, Float_t radius, Float_t radL, Float_t phiRes, Float_t zRes, Float_t eff) {
  //
  // Add additional layer to the list of layers (ordered by radius)
  // 

  CylLayer *newLayer = (CylLayer*) fLayers.FindObject(name);

  if (!newLayer) {
    newLayer = new CylLayer(name);
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
	CylLayer *l = (CylLayer*)fLayers.At(i);
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
    if (!(newLayer->isDead)) fNumberOfActiveLayers += 1;


  } else {
    printf("Layer with the name %s does already exist\n",name);
  }
  

}

void Detector::KillLayer(char *name) {
  //
  // Marks layer as dead. Contribution only by Material Budget
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot mark as dead\n",name);
  else {
     tmp->phiRes = 999999;
     tmp->zRes = 999999;
     if (!(tmp->isDead)) {
       tmp->isDead = kTRUE;
       fNumberOfActiveLayers -= 1;
     }     
  }
}

void Detector::SetRadius(char *name, Float_t radius) {
  //
  // Set layer radius [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
 

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

Float_t Detector::GetRadius(char *name) {
  //
  // Return layer radius [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get radius\n",name);
  else 
    return tmp->radius;

  return 0;
}

void Detector::SetRadiationLength(char *name, Float_t radL) {
  //
  // Set layer material [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer material\n",name);
  else {
    tmp->radL = radL;
  }
}

Float_t Detector::GetRadiationLength(char *name) {
  //
  // Return layer radius [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer material\n",name);
  else 
    return tmp->radL;
    
  return 0;
  
}

void Detector::SetResolution(char *name, Float_t phiRes, Float_t zRes) {
  //
  // Set layer resolution in [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set resolution\n",name);
  else {

    Bool_t wasDead = tmp->isDead;
    
    tmp->phiRes = phiRes;
    tmp->zRes = zRes;
    
    if (zRes==RIDICULOUS && phiRes==RIDICULOUS) {
      tmp->isDead = kTRUE;
      if (!wasDead) fNumberOfActiveLayers -= 1;
    } else {
      tmp->isDead = kFALSE;
      if (wasDead) fNumberOfActiveLayers += 1;
    }


  }
}

Float_t Detector::GetResolution(char *name, Int_t axis) {
  //
  // Return layer resolution in [cm]
  // axis = 0: resolution in rphi
  // axis = 1: resolution in z
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get resolution\n",name);
  else {
    if (axis==0) return tmp->phiRes;
    if (axis==1) return tmp->zRes;
    printf("error: axis must be either 0 or 1 (rphi or z axis)\n");
  }
  return 0;
}

void Detector::SetLayerEfficiency(char *name, Float_t eff) {
  //
  // Set layer efficnecy (prop that his is missed within this layer) 
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set layer efficiency\n",name);
  else {
    tmp->eff = eff;
  }
}

Float_t Detector::GetLayerEfficiency(char *name) {
  //
  // Get layer efficnecy (prop that his is missed within this layer) 
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot get layer efficneicy\n",name);
  else 
    return tmp->eff;
    
  return 0;
  
}

void Detector::RemoveLayer(char *name) {
  //
  // Removes a layer from the list
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot remove it\n",name);
  else {
    Bool_t wasDead = tmp->isDead;
    fLayers.Remove(tmp);
    fNumberOfLayers -= 1;
    if (!wasDead) fNumberOfActiveLayers -= 1;
  }
}


void Detector::PrintLayout() {
  //
  // Prints the detector layout
  //

  printf("Detector %s: \"%s\"\n",GetName(),GetTitle());
  
  if (fLayers.GetEntries()>0) 
    printf("  Name \t\t r [cm] \t  X0 \t  phi & z res [um]\n");

  CylLayer *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (CylLayer*)fLayers.At(i);
  
    // don't print all the tpc layers
    TString name(tmp->GetName());
    if (name.Contains("tpc") && (!name.Contains("tpc_0")) ) continue;

    printf("%d. %s \t %03.2f   \t%1.4f\t  ",i,
 	   tmp->GetName(), tmp->radius, tmp->radL);
    if (tmp->phiRes==RIDICULOUS) 
      printf("  -  ");
    else
      printf("%3.0f   ",tmp->phiRes*10000);
    if (tmp->zRes==RIDICULOUS) 
      printf("  -\n");
    else
      printf("%3.0f\n",tmp->zRes*10000);
  }
}

void Detector::PlotLayout(Int_t plotDead) {
  //
  // Plots the detector layout in Front view
  //

  Double_t x0=0, y0=0;

  TGraphErrors *gr = new TGraphErrors();
  gr->SetPoint(0,0,0);
  CylLayer *lastLayer = (CylLayer*)fLayers.At(fLayers.GetEntries()-1);  Double_t maxRad = lastLayer->radius;
  gr->SetPointError(0,maxRad,maxRad);
  gr->Draw("APE");
  

  CylLayer *tmp = 0;
  for (Int_t i = fLayers.GetEntries()-1; i>=0; i--) {
    tmp = (CylLayer*)fLayers.At(i);
  

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

    if (!tmp->isDead || plotDead) layEl->Draw();
  
  }

}



void Detector::AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip) {
  //
  // Emulates the TPC
  // 
  // skip=1: Use every padrow, skip=2: Signal in every 2nd padrow 


  AddLayer((char*)"IFC",   77.8,0.01367); // Inner Field cage
  
  // % Radiation Lengths ... Average per TPC row  (i.e. total/159 )
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

void Detector::RemoveTPC() {

  // flag as dead, although resolution is ok ... makes live easier in the prints ... ;-)
  CylLayer *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (CylLayer*)fLayers.At(i);  
    TString name(tmp->GetName());
    if (name.Contains("tpc")) { RemoveLayer((char*)name.Data()); i--; }
  }
  RemoveLayer((char*)"IFC");
  
}


Double_t Detector::ThetaMCS ( Double_t mass, Double_t radLength, Double_t momentum ) const
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


Double_t Detector::ProbGoodHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
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


Double_t Detector::ProbGoodChiSqHit ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Victor Perevoztchikov and Howard Wieman: http://rnc.lbl.gov/~wieman/HitFinding2DXsq.htm
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  Double_t sx, goodHit ;
  sx = 2 * TMath::Pi() *  searchRadiusRPhi * searchRadiusZ * HitDensity(radius) ;
  goodHit =  1./(1+sx) ;
  return ( goodHit ) ;  
}

Double_t Detector::ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ ) 
{
  // Based on work by Ruben Shahoyen 
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  // Plus, in addition, taking a "confidence level" and the "layer efficiency" into account 
  // Following is correct for 2 DOF

  Double_t gamma = -2 *TMath::Log(fConfLevel); // quantile at cut of confidence level
  Double_t alpha = (1 + 2 * TMath::Pi() * HitDensity(radius) * searchRadiusRPhi * searchRadiusZ)/2; 
  Double_t goodHit = leff/(2*alpha) * (1 - TMath::Exp(-alpha*gamma));

  return ( goodHit ) ;  
}


Double_t Detector::HitDensity ( Double_t radius ) 
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
      arealDensity  = OneEventHitDensity(dNdEtaCent,radius)  ; // Fast detectors see central collision density (only)
      arealDensity += OtherBackground*OneEventHitDensity(dNdEtaMinB,radius)  ;  // Increase density due to background 
    }

  if (radius < fMaxRadiusSlowDet )
    { // Note that IntegratedHitDensity will always be minB one event, or more, even if integration time => zero.
      arealDensity  = OneEventHitDensity(dNdEtaCent,radius) 
	            + IntegratedHitDensity(dNdEtaMinB,radius) 
	            + UpcHitDensity(radius) ;
      arealDensity += OtherBackground*IntegratedHitDensity(dNdEtaMinB,radius) ;  
      // Increase density due to background 
    } 

  return ( arealDensity ) ;  
}


double Detector::OneEventHitDensity( Double_t multiplicity, Double_t radius ) const
{
  // This is for one event at the vertex.  No smearing.
  double den   = multiplicity / (2.*TMath::Pi()*radius*radius) ; // 2 eta ?
  // note: surface of sphere is  '4*pi*r^2'
  //       surface of cylinder is '2*pi*r* h' 
  return den ;
} 


double Detector::IntegratedHitDensity(Double_t multiplicity, Double_t radius)
{ 
  // The integral of minBias events smeared over a gaussian vertex distribution.
  // Based on work by Yan Lu 12/20/2006, all radii in centimeters.

  Double_t zdcHz = Luminosity * 1.e-24 * CrossSectionMinB ;
  Double_t den   = zdcHz * fIntegrationTime/1000. * multiplicity * Dist(0., radius) / (2.*TMath::Pi()*radius) ;

  // Note that we do not allow the rate*time calculation to fall below one minB event at the vertex.
  if ( den < OneEventHitDensity(multiplicity,radius) )  den = OneEventHitDensity(multiplicity,radius) ;  

  return den ;
} 


double Detector::UpcHitDensity(Double_t radius)
{ 
  // QED electrons ...

  Double_t mUPCelectrons ;                                 ;  
  mUPCelectrons =  fLhcUPCscale * (1.23 - radius/6.5)      ;  // Fit to Kai Schweda summary tables at RHIC * 'scale' for LHC
  if ( mUPCelectrons < 0 ) mUPCelectrons =  0.0             ;  // UPC electrons fall off quickly and don't go to large R
  mUPCelectrons *= IntegratedHitDensity(dNdEtaMinB,radius) ;  // UPCs increase Mulitiplicty ~ proportional to MinBias rate
  mUPCelectrons *= UPCBackgroundMultiplier                 ;  // Allow for an external multiplier (eg 0-1) to turn off UPC

  return mUPCelectrons ;
} 


double Detector::Dist(double z, double r)
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

Double_t Detector::D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][400] ) const {
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
      
    if ( pionindex >= 400 ) pionindex = 399 ;
    if ( pionindex >= 0 )   effp = corrEfficiency[0][pionindex] ;
    if ( pionindex <  0 )   effp = (corrEfficiency[0][1]-corrEfficiency[0][0])*pionindex + corrEfficiency[0][0] ; // Extrapolate if reqd
    if ( effp < 0 )         effp = 0 ;

    if ( kaonindex >= 400 ) kaonindex = 399 ;
    if ( kaonindex >= 0 )   effk = corrEfficiency[1][kaonindex] ;
    if ( kaonindex <  0 )   effk = (corrEfficiency[1][1]-corrEfficiency[1][0])*kaonindex + corrEfficiency[1][0] ; // Extrapolate if reqd
    if ( effk < 0 )         effk = 0 ;

    // Note that we assume that the Kaon Decay efficiency has already been inlcuded in the kaon efficiency used here.
      
    sum += effp * effk ;
 
  }    
  
  Double_t mean =sum/360; 
  return mean ;
  
}


void Detector::SolveDOFminusOneAverage() {
  // 
  // Short study to address "# layers-1 efficiencies"
  // saves the means in the according arrays
  // Note: Obviously, does not work for the Telescope equation 
  //
  
  Double_t fMomentumResM[400], fResolutionRPhiM[400], fResolutionZM[400]; 
  Double_t efficiencyM[3][400];
  for (Int_t i=0; i<400; i++) {
    fMomentumResM[i] = 0;   // Momentum resolution
    fResolutionRPhiM[i] = 0;   // Resolution in R
    fResolutionZM[i] = 0; // Resolution in Z
    for (Int_t part=0; part<3; part++) 
      efficiencyM[part][i] = 0; // efficiencies
  }

  // loop over active layers in ITS (remove 1 by 1)
  Int_t nITSLayers = 0;
  CylLayer *layer =0;
  for (Int_t j=0; j<(fLayers.GetEntries()-1); j++) { 
    layer = (CylLayer*)fLayers.At(j);
    TString name(layer->GetName());
    if (name.Contains("tpc")) continue;
    if (!(layer->isDead))  {

      nITSLayers++; 
      printf("Kill Layer %s\n",name.Data());
      Double_t rRes = GetResolution((char*)name.Data(),0);
      Double_t zRes = GetResolution((char*)name.Data(),1);
      KillLayer((char*)name.Data());
      //   PrintLayout();
      SolveViaBilloir(1,0); 

      // produce sum for the mean calculation
      for (Int_t i=0; i<400; i++) {
	fMomentumResM[i] += fMomentumRes[i];   // Momentum resolution
	fResolutionRPhiM[i] += fResolutionRPhi[i];   // Resolution in R
	fResolutionZM[i] += fResolutionZ[i]; // Resolution in Z
	for (Int_t part=0; part<3; part++) 
	  efficiencyM[part][i] += fEfficiency[part][i]; // efficiencies
      }

      // "Restore" layer ...
      SetResolution((char*)name.Data(),rRes,zRes); 
      
    }
  }
  
  // save means in "std. Arrays"
  for (Int_t i=0; i<400; i++) {
    fMomentumRes[i] = fMomentumResM[i]/nITSLayers;   // Momentum resolution
    fResolutionRPhi[i] = fResolutionRPhiM[i]/nITSLayers;   // Resolution in R
    fResolutionZ[i] = fResolutionZM[i]/nITSLayers; // Resolution in Z
    for (Int_t part=0; part<3; part++) 
      fEfficiency[part][i] = efficiencyM[part][i]/nITSLayers; // efficiencies
  }


}

void Detector::SolveViaBilloir(Int_t flagD0,Int_t print, Bool_t allPt, Double_t meanPt) {
  //
  // Solves the current geometry with the Billoir technique 
  // ( see P. Billoir, Nucl. Instr. and Meth. 225 (1984), p. 352. )
  //

  Int_t nPt = 400;
  // Clean up ......
  for (Int_t i=0; i<kMaxNumberOfDetectors; i++) {
    for (Int_t j=0; j<nPt; j++) {
      fDetPointRes[i][j]  = RIDICULOUS;
      fDetPointZRes[i][j] = RIDICULOUS;
      fTransMomenta[i] =0;
      fMomentumRes[i] =0;
      fResolutionRPhi[i] =0;
    }
  }
  
  if (!allPt) { // not the whole pt range -> allows a faster minimization at a defined 'meanpt'
    nPt = 3;
  }



  // Calculate track parameters using Billoirs method of matrices

  Double_t pt, pz, lambda, momentum, rho, deltaPoverP  ;
  Double_t layerThickness, aMCS, mmm, charge ;
  Double_t vF1, vF2, vF3, vG1, vG2, vG3 ;
  Double_t dx, sigma2, sigma2Z ;
  Double_t mass[3] ;
  Int_t printOnce = 1 ;

  mass[0] = PionMass ; mass[1] = KaonMass ;  // Loop twice for the D0;  first pi then k 

  mass[2] = fParticleMass;  // third loop

  Int_t mStart =0; 
  if (!flagD0) mStart = 2; // pion and kaon is skipped -> fast mode

  for ( Int_t massloop = mStart ; massloop < 3 ; massloop++ )  { 
    
    // PseudoRapidity OK, used as an angle
    lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity))  ; 
  
    // Track moves along the x axis and deviations are in the y and z directions.
    Double_t y, z, a, b, d, e ;  
    // a = py/px = tan phi, b = pz/px = tan lambda / cos phi, d = 0.3*Ze/Abs(p)

    for ( Int_t i = 0 ; i < nPt ; i++ ) { // pt loop
 
      CylLayer *last = (CylLayer*) fLayers.At((fLayers.GetEntries()-1));

      // Convert to Meters, Tesla, and GeV
      
      // Starting values based on radius of outermost layer ... log10 steps to ~20 GeV
      Double_t bigRad = last->radius/100 /2 ;	// min. pt which the algorithm below could handle
      if (bigRad<0.61) bigRad=0.61; // -> min pt around 100 MeV for Bz=0.5T (don't overdo it ... ;-) )
      fTransMomenta[i] =  ( 0.3*bigRad*TMath::Abs(fBField) ) - 0.08 + TMath::Power(10,2.3*i/nPt) / 10.0 ; 
      if (!allPt) { // just 3 points around meanPt
	fTransMomenta[i] = meanPt-0.1+i*0.1;
      }
   
      
      // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
      // These are the EndPoint values for y, z, a, b, and d
      pt  =  fTransMomenta[i]                       ;  // GeV/c
      rho =  pt / TMath::Abs( 0.3 * fBField );  // Radius of curvature of the track (meters)
      momentum = pt / TMath::Cos(lambda)     ;  // Total momentum  
      charge   = -1                          ;  // Assume an electron 
      pz  =  pt * TMath::Tan(lambda)         ;  // GeV/c
      d   =  0.3 * charge / momentum         ;  // Its an electrons so q = -1, which makes d a negative number

      //      printf("%d - pt %lf r%lf | %lf %lf\n",massloop,fTransMomenta[i],(last->radius)/100,momentum, d);

      // Matrices for Billoir
      TMatrixD bigA(2,2); TMatrixD littleA(2,2);
      TMatrixD mIstar(5,5); TMatrixD mD(5,5); TMatrixD mDT(5,5); TMatrixD mA(5,5); TMatrixD mM(5,5); TMatrixD mWork(5,5) ; 
      TMatrixD mCovStart(5,5);
      bigA.Zero();littleA.Zero();mD.Zero();mDT.Zero();mA.Zero();mM.Zero();mWork.Zero();

      // Set Detector-Efficiency Storage area to unity
      fEfficiency[massloop][i] = 1.0 ;
      // Back-propagate the covariance matrix along the track.  I is the inverse of the covariance matrix. 
      // Start with small variations to a straight line 'information' matrix.  Tuning this matrix is 
      // equivalent to 'starting the recursion' in Billoirs paper.  It is a tricky business.  You need test
      // data and many layers for the matrices to stabilize.  Do not believe results for a small number of layers.
      // In our case, always include the TPC and this matrix will be well conditioned before it gets to the Si.


      mCovStart.UnitMatrix(); // start with unity (-> sigma estimate of 1cm)
      
      // track parametrization: (y,z,a,b,1/p) with a=py/px=tan(phi) and b = pz/px = tan(lambda)/cos(phi)
      mCovStart(0,0) = (2*2);     // 2cm variance in y
      mCovStart(1,1) = (2*2);     // 2cm variance in z
      mCovStart(2,2) = (0.1*0.1); // approx 0.1rad in angles phi and lambda
      mCovStart(3,3) = (1*1);     // approx 0.1rad in angle phi and lambda
      mCovStart(4,4) = (1*1);     // 1 GeV error in pt
      

      mIstar = mCovStart; // start with above covariance Matrix
      mIstar.Invert();
      
      // mIstar.UnitMatrix(); // start with unity (-> sigma estimate of ~ 1cm)

      // 'A' is the "angle" matrix for aMCS at each step.
      // 'D' is the "distance" matrix at each step.  It propagates the particle backwards along the track. 
      // 'M' is the "measurement" matrix. It represents the measurement resolution at each step.
    
      CylLayer *layer = 0;
      CylLayer *nextlayer = 0;
      
      // find last "active layer" - start tracking at the last active layer      
      Int_t lastActiveLayer = 0;
      for (Int_t j=(fLayers.GetEntries()-1); j>0; j--) { 
	layer = (CylLayer*)fLayers.At(j);
	if (!(layer->isDead)) { // is alive
	  lastActiveLayer = j;
	  break;
	}
      }
    
      for (Int_t j=lastActiveLayer; j>0; j--) {  // Layer loop

	layer = (CylLayer*)fLayers.At(j);
	nextlayer = (CylLayer*)fLayers.At(j-1);

	// Convert to Meters, Tesla, and GeV
	Float_t radius = layer->radius/100;
	Float_t radLength = layer->radL;

	Float_t nextRadius = nextlayer->radius /100;
	Float_t nextPhiResolution = nextlayer->phiRes /100;
	Float_t nextZResolution = nextlayer->zRes /100;
   	
	//	if ( radius >= rho*TMath::Cos(TMath::Pi()/4) )    // (meters) protect against too low pt
	//	  { // printf("pt lower bound is too low for this algorithm to take this layer into account (r=%lf)\n",radius);
	//	    continue ; }
	
	// Jims approximation of the "ideal position", breaks down when rho~radius because of the hyperbola approximation
	//	y   =  rho - TMath::Sqrt(rho*rho-radius*radius)  ; // These are 'ideal' locations and
	//	a   =  radius / ( rho - y )                      ; // not propagated values which should  
	//	b   =  rho * TMath::Tan(lambda) / ( rho - y )    ; // be done if we had data. But we don't.
	//	z   =  rho * TMath::Tan(lambda) * TMath::ATan(a) ; 

	// Helix format -> exact intersection at the layer radii
	y=0;z=0; // not used
	a = radius*radius/TMath::Sqrt(4*radius*radius*rho*rho - TMath::Power(radius,4)); // = Tan phi = y_c/x_c (in cartesian coordinates)
	b = TMath::Tan(lambda)/TMath::Cos(TMath::ATan(a));

	e   =  TMath::Sqrt( 1 + a*a + b*b )              ; 
      
	layerThickness = radLength / TMath::Sin(TMath::Pi()/2 - lambda)  ; // X0 in lambda direction
      
	aMCS =  ThetaMCS(mass[massloop], layerThickness, momentum)           ; 
	mmm =  ( 1 + a*a + b*b ) ; 
	bigA(0,0) = mmm*(1+a*a) ; bigA(0,1) = mmm*a*b     ; littleA(0,0) = mIstar(2,2)  ; littleA(0,1) = mIstar(2,3) ; 
	bigA(1,0) = mmm*a*b     ; bigA(1,1) = mmm*(1+b*b) ; littleA(1,0) = mIstar(3,2)  ; littleA(1,1) = mIstar(3,3) ; 
	littleA = bigA.Invert() + aMCS*aMCS*littleA ;
	littleA = littleA.Invert() ;
	mA(0,0) = 0.0 ;  mA(0,1) = 0.0  ;  mA(0,2) = 0.0          ;  mA(0,3) = 0.0          ;  mA(0,4) = 0.0  ; 
	mA(1,0) = 0.0 ;  mA(1,1) = 0.0  ;  mA(1,2) = 0.0          ;  mA(1,3) = 0.0          ;  mA(1,4) = 0.0  ;
	mA(2,0) = 0.0 ;  mA(2,1) = 0.0  ;  mA(2,2) = littleA(0,0) ;  mA(2,3) = littleA(0,1) ;  mA(2,4) = 0.0  ;
	mA(3,0) = 0.0 ;  mA(3,1) = 0.0  ;  mA(3,2) = littleA(1,0) ;  mA(3,3) = littleA(1,1) ;  mA(3,4) = 0.0  ;
	mA(4,0) = 0.0 ;  mA(4,1) = 0.0  ;  mA(4,2) = 0.0          ;  mA(4,3) = 0.0          ;  mA(4,4) = 0.0  ;
	mIstar = mIstar - aMCS*aMCS*mIstar*mA*mIstar ;
	dx     = ( radius - nextRadius )         ;        // Billoir works with absolute magnitude of distance
	vF3     = e * ( -1 * ( 1 + a*a ) * fBField )   ;        // Assume fBField is on Z axis, only.
	vF1     = d * ( a*vF3/(e*e) - 2*e*a*fBField )   ;
	vF2     = d * ( b*vF3/(e*e) )                  ;
	vG3     = e * ( -1 * a * b * fBField )         ;
	vG1     = d * ( a*vG3/(e*e) + e*b*fBField )     ;
	vG2     = d * ( b*vG3/(e*e) - e*a*fBField )     ;	      
	mD(0,0) = 1.0 ; mD(0,1) = 0.0 ; mD(0,2) = -1*dx+vF1*dx*dx/2 ; mD(0,3) = vF2*dx*dx/2  ;  mD(0,4) = vF3*dx*dx/2 ;
	mD(1,0) = 0.0 ; mD(1,1) = 1.0 ; mD(1,2) = vG1*dx*dx/2 ;  mD(1,3) = -1*dx+vG2*dx*dx/2 ;  mD(1,4) = vG3*dx*dx/2 ;
	mD(2,0) = 0.0 ; mD(2,1) = 0.0 ; mD(2,2) = 1-vF1*dx    ;  mD(2,3) = -1*vF2*dx         ;  mD(2,4) = -1*vF3*dx   ;
	mD(3,0) = 0.0 ; mD(3,1) = 0.0 ; mD(3,2) = -1*vG1*dx   ;  mD(3,3) = 1-vG2*dx          ;  mD(3,4) = -1*vG3*dx   ;
	mD(4,0) = 0.0 ; mD(4,1) = 0.0 ; mD(4,2) = 0.0        ;  mD(4,3) = 0.0              ;  mD(4,4) = 1.0        ;
	mDT.Transpose(mD) ;
	mIstar  = mDT.Invert() * mIstar * mD.Invert() ;
	// Prepare to save Detector efficiencies
	mWork = mIstar  ;          // Working copy of matrix 
	mWork.Invert() ;          // Invert the Matrix to recover the convariance matrix
	fDetPointRes [j-1][i]     =  TMath::Sqrt( mWork(0,0) )  ;     // result in meters
	fDetPointZRes[j-1][i]     =  TMath::Sqrt( mWork(1,1) )  ;     // result in meters
	// End save

	sigma2  = ( nextPhiResolution  * nextPhiResolution  )   ; 
	sigma2Z = ( nextZResolution    * nextZResolution )   ;
	mM(0,0) = 1/sigma2; mM(1,1) = 1/sigma2Z  ; // setting detector resolution matrix
	mIstar = mIstar + mM ;
	
      }
      
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
        
  
      if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	printf("Number of active layers: %d\n",fNumberOfActiveLayers) ;
	printf("Mass of tracked particle: %f (at pt=%5.0lf MeV)\n",fParticleMass,fTransMomenta[i]*1000);
	printf("Name   Radius Thickness PointResOn PointResOnZ  DetRes  DetResZ  Density Efficiency\n") ;
	//	printOnce =0;
      }

      // print out and efficiency calculation
      for (Int_t j=(fLayers.GetEntries()-1); j>=0; j--) {  // Layer loop
	
	layer = (CylLayer*)fLayers.At(j);
	
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
	    if (name.Contains("tpc") && (!name.Contains("tpc_0")) ) continue;

	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {


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

	    if (!name.Contains("tpc")) fEfficiency[massloop][i] *= layerEfficiency;

	}
	/*
	// vertex print
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1 && radius==0) {
	printf("%s:\t -----    ----- %10.0f %11.0f \n", layer->GetName(),fDetPointRes[j][i]*1.e6, fDetPointZRes[j][i]*1.e6);
	}
	*/
      }
      if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	if (fNumberOfActiveLayers >=15) printOnce = 0 ;
	printf("\n")  ;
      }
      



      if (fNumberOfActiveLayers <15 ) {



	// BACKWORD TRACKING +++++++++++++++++
	// number of layers is quite low ... efficiency calculation was probably nonsense 
	// Tracking outward (backword) to get reliable efficiencies from "smoothed estimates"

	// For below, see paper, NIM A262 (1987) p.444, eqs.12.
	// Equivalently, one can simply combine the forward and backward estimates. Assuming
	// pf,Cf and pb,Cb as extrapolated position estimates and errors from fwd and bwd passes one can
	// use a weighted estimate Cw = (Cf^-1 + Cb^-1)^-1,  pw = Cw (pf Cf^-1 + pb Cb^-1).
	// Surely, for the most extreme point, where one error matrices is infinite, this does not change anything.

	Bool_t doLikeAliRoot = kFALSE; // don't do the "combined info" but do like in Aliroot

	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {	  
	  printf("- Numbers of active layer is low (%d):\n    -> \"outward\" fitting done as well to get reliable eff.estimates\n",
		 fNumberOfActiveLayers);	  
	}
	
	// RESET Covariance Matrix ( to 10 x the estimate -> as it is done in AliExternalTrackParam)
	//	mIstar.UnitMatrix(); // start with unity
	if (doLikeAliRoot) {
	  Double_t factor = 10;
	  TMatrixD diagElements(5,1);
	  for (Int_t ip =0; ip<5; ip++) diagElements(ip,0) = mIstar(ip,ip)*factor;
	  mIstar.Zero();
	  for (Int_t ip =0; ip<5; ip++) mIstar(ip,ip) = diagElements(ip,0);
	} else {
	  // Complete RESET of the covariance matrix
	  mIstar = mCovStart; 
	  // starts with above covariance Matrix and uses "combinded resolutions" (all infos), see below
	  // this is the "honest" approach without biasing the results
	}
	mIstar.Invert() ; // undo the revert from above	

	// Clean up and storing of "forward estimates"
	Double_t detPointResForw[kMaxNumberOfDetectors][400], detPointZResForw[kMaxNumberOfDetectors][400] ; 
	for (Int_t k=0; k<kMaxNumberOfDetectors; k++) {
	  for (Int_t l=0; l<nPt; l++) {
	    detPointResForw[k][l]  = fDetPointRes[k][l];
	    if (!doLikeAliRoot) fDetPointRes[k][l]  = RIDICULOUS;
	    detPointZResForw[k][l] = fDetPointZRes[k][l];
	    if (!doLikeAliRoot) fDetPointZRes[k][l] = RIDICULOUS;
	  }
	}
 	
	// find first "active layer" - start tracking at the first active layer      
	Int_t firstActiveLayer = 0;
	for (Int_t j=0; j<(fLayers.GetEntries()-1); j++) { 
	  layer = (CylLayer*)fLayers.At(j);
	  if (!(layer->isDead)) { // is alive
	    firstActiveLayer = j;
	    break;
	  }
	}

   	for (Int_t j=firstActiveLayer; j<(fLayers.GetEntries()-1); j++) {  // Layer loop
	  
	  layer = (CylLayer*)fLayers.At(j);
	  nextlayer = (CylLayer*)fLayers.At(j+1);
	  
	  // Convert to Meters, Tesla, and GeV
	  Float_t radius = layer->radius/100;
	  Float_t radLength = layer->radL;
	  
	  Float_t nextRadius = nextlayer->radius /100;
	  Float_t nextPhiResolution = nextlayer->phiRes /100;
	  Float_t nextZResolution = nextlayer->zRes /100;

	  y   =  rho - TMath::Sqrt(rho*rho-radius*radius)  ; // These are 'ideal' locations and
	  a   =  radius / ( rho - y )                      ; // not propagated values which should  
	  b   =  rho * TMath::Tan(lambda) / ( rho - y )    ; // be done if we had data. But we don't.
	  z   =  rho * TMath::Tan(lambda) * TMath::ATan(a) ; 
	  e   =  TMath::Sqrt( 1 + a*a + b*b )              ; 
      
	  layerThickness = radLength / TMath::Sin(TMath::Pi()/2 - lambda)  ; 
      
	  aMCS =  ThetaMCS(mass[massloop], layerThickness, momentum)           ; 
	  mmm =  ( 1 + a*a + b*b ) ; 
	  bigA(0,0) = mmm*(1+a*a) ; bigA(0,1) = mmm*a*b     ; littleA(0,0) = mIstar(2,2)  ; littleA(0,1) = mIstar(2,3) ; 
	  bigA(1,0) = mmm*a*b     ; bigA(1,1) = mmm*(1+b*b) ; littleA(1,0) = mIstar(3,2)  ; littleA(1,1) = mIstar(3,3) ; 
	  littleA = bigA.Invert() + aMCS*aMCS*littleA ;
	  littleA = littleA.Invert() ;
	  mA(0,0) = 0.0 ;  mA(0,1) = 0.0  ;  mA(0,2) = 0.0          ;  mA(0,3) = 0.0          ;  mA(0,4) = 0.0  ; 
	  mA(1,0) = 0.0 ;  mA(1,1) = 0.0  ;  mA(1,2) = 0.0          ;  mA(1,3) = 0.0          ;  mA(1,4) = 0.0  ;
	  mA(2,0) = 0.0 ;  mA(2,1) = 0.0  ;  mA(2,2) = littleA(0,0) ;  mA(2,3) = littleA(0,1) ;  mA(2,4) = 0.0  ;
	  mA(3,0) = 0.0 ;  mA(3,1) = 0.0  ;  mA(3,2) = littleA(1,0) ;  mA(3,3) = littleA(1,1) ;  mA(3,4) = 0.0  ;
	  mA(4,0) = 0.0 ;  mA(4,1) = 0.0  ;  mA(4,2) = 0.0          ;  mA(4,3) = 0.0          ;  mA(4,4) = 0.0  ;
	  mIstar = mIstar - aMCS*aMCS*mIstar*mA*mIstar ;
	  dx     = ( radius - nextRadius )         ;        // Billoir works with absolute magnitude of distance
	  vF3     = e * ( -1 * ( 1 + a*a ) * fBField )   ;        // Assume fBField is on Z axis, only.
	  vF1     = d * ( a*vF3/(e*e) - 2*e*a*fBField )   ;
	  vF2     = d * ( b*vF3/(e*e) )                  ;
	  vG3     = e * ( -1 * a * b * fBField )         ;
	  vG1     = d * ( a*vG3/(e*e) + e*b*fBField )     ;
	  vG2     = d * ( b*vG3/(e*e) - e*a*fBField )     ;	      
	  mD(0,0) = 1.0 ; mD(0,1) = 0.0 ; mD(0,2) = -1*dx+vF1*dx*dx/2 ; mD(0,3) = vF2*dx*dx/2  ;  mD(0,4) = vF3*dx*dx/2 ;
	  mD(1,0) = 0.0 ; mD(1,1) = 1.0 ; mD(1,2) = vG1*dx*dx/2 ;  mD(1,3) = -1*dx+vG2*dx*dx/2 ;  mD(1,4) = vG3*dx*dx/2 ;
	  mD(2,0) = 0.0 ; mD(2,1) = 0.0 ; mD(2,2) = 1-vF1*dx    ;  mD(2,3) = -1*vF2*dx         ;  mD(2,4) = -1*vF3*dx   ;
	  mD(3,0) = 0.0 ; mD(3,1) = 0.0 ; mD(3,2) = -1*vG1*dx   ;  mD(3,3) = 1-vG2*dx          ;  mD(3,4) = -1*vG3*dx   ;
	  mD(4,0) = 0.0 ; mD(4,1) = 0.0 ; mD(4,2) = 0.0        ;  mD(4,3) = 0.0              ;  mD(4,4) = 1.0        ;
	  mDT.Transpose(mD) ;
	  mIstar  = mDT.Invert() * mIstar * mD.Invert() ;
	  // Prepare to save Detector efficiencies
	  mWork = mIstar  ;          // Working copy of matrix 
	  mWork.Invert() ;          // Invert the Matrix to recover the convariance matrix
	  fDetPointRes [j+1][i]     =  TMath::Sqrt( mWork(0,0) )  ;     // result in meters
	  fDetPointZRes[j+1][i]     =  TMath::Sqrt( mWork(1,1) )  ;     // result in meters
	  // End save
	  sigma2  = ( nextPhiResolution  * nextPhiResolution  )   ; 
	  sigma2Z = ( nextZResolution    * nextZResolution )   ;
	  mM(0,0) = 1/sigma2; mM(1,1) = 1/sigma2Z  ;
	  mIstar = mIstar + mM ;	
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
  	
	// Weighted combination of the forward and backward estimates
	if (!doLikeAliRoot) {
	  for (Int_t j=(fLayers.GetEntries()-1); j>=0; j--) {  
	    fDetPointRes[j][i] = 1/(1/detPointResForw[j][i] + 1/fDetPointRes[j][i]); 
	    fDetPointZRes[j][i] = 1/(1/detPointZResForw[j][i] + 1/fDetPointZRes[j][i]); 
	  }
	}
	// Set Detector-Efficiency Storage area to unity
	fEfficiency[massloop][i] = 1.0 ;
     
	// print out and efficiency calculation
	for (Int_t j=(fLayers.GetEntries()-1); j>=0; j--) {  // Layer loop
	  
	  layer = (CylLayer*)fLayers.At(j);
	
	  // Convert to Meters, Tesla, and GeV
	  Float_t radius = layer->radius /100;
	  Float_t phiRes = layer->phiRes /100;
	  Float_t zRes = layer->zRes /100;
	  Float_t radLength = layer->radL;
	  Float_t leff = layer->eff;
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
	    if (name.Contains("tpc") && (!name.Contains("tpc_0")) ) continue;

	    if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {


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

	    if (!name.Contains("tpc")) fEfficiency[massloop][i] *= layerEfficiency;
	    
	  }
	}
	if (print == 1 && fTransMomenta[i] >= meanPt && massloop == 2 && printOnce == 1) {
	  printOnce = 0 ;
	  printf("\n")  ;
	}
      }      
    } // mass loop
  } // pt loop
  
}


TGraph * Detector::GetGraphMomentumResolution(Int_t color, Int_t linewidth) {
  //
  // returns the momentum resolution 
  //
  
  TGraph *graph = new TGraph(400, fTransMomenta, fMomentumRes);
  graph->SetTitle("Momentum Resolution .vs. Pt" ) ;
  //  graph->GetXaxis()->SetRangeUser(0.,5.0) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Momentum Resolution [%]") ;
  graph->GetYaxis()->CenterTitle();

  graph->SetMaximum(20) ;
  graph->SetMinimum(0.1) ;
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph * Detector::GetGraphPointingResolution(Int_t axis, Int_t color, Int_t linewidth) {
 
  // Returns the pointing resolution
  // axis = 0 ... rphi pointing resolution
  // axis = 1 ... z pointing resolution
  //

  TGraph * graph =  0;

  if (axis==0) {
    graph = new TGraph ( 400, fTransMomenta, fResolutionRPhi ) ;
    graph->SetTitle("R-#phi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("R-#phi Pointing Resolution [#mum]") ;
  } else {
    graph =  new TGraph ( 400, fTransMomenta, fResolutionZ ) ;
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution [#mum]") ;
  }
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineWidth(linewidth);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  
  return graph;

}


TGraph * Detector::GetGraphPointingResolutionTeleEqu(Int_t axis,Int_t color, Int_t linewidth) {
  //
  // returns the Pointing resolution (accoring to Telescope equation)
  // axis =0 ... in rphi
  // axis =1 ... in z
  //
  
  Double_t resolution[400];

  Double_t layerResolution[2];
  Double_t layerRadius[2];
  Double_t layerThickness[2];

  Int_t count =0; // search two first active layers
  printf("Telescope equation for layers:  ");
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    CylLayer *l = (CylLayer*)fLayers.At(i);
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

  for ( Int_t i = 0 ; i < 400 ; i++ ) { 
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



  TGraph* graph = new TGraph ( 400, fTransMomenta, resolution ) ;
   
  if (axis==0) {
    graph->SetTitle("RPhi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("RPhi Pointing Resolution [#mum] ") ;
  } else {
    graph->SetTitle("Z Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("Z Pointing Resolution [#mum] ") ;
  }
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]") ;
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

TGraph * Detector::GetGraphRecoEfficiency(Int_t particle,Int_t color, Int_t linewidth) {
  //
  // particle = 0 ... choosen particle (setted particleMass)
  // particle = 1 ... Pion
  // particle = 2 ... Kaon
  // particle = 3 ... D0
  //
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 
  
  Double_t particleEfficiency[400]; // with chosen particle mass
  Double_t kaonEfficiency[400], pionEfficiency[400], d0efficiency[400]; 
  Double_t partEfficiency[2][400];
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < 400 ; j++ ) { 
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
    for ( Int_t j = 0 ; j < 400 ; j++ ) {
      d0efficiency[j] = D0IntegratedEfficiency(fTransMomenta[j],partEfficiency);
    }
  } else { 
    for ( Int_t j = 0 ; j < 400 ; j++ ) { 
      particleEfficiency[j] = fEfficiency[2][j]* AcceptanceOfTpcAndSi;
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( 400, fTransMomenta, particleEfficiency ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( 400, fTransMomenta, pionEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( 400, fTransMomenta, kaonEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==3) {
    graph = new TGraph ( 400, fTransMomenta, d0efficiency ) ;
    graph->SetLineStyle(kDashed);
  } else 
    return 0;

  graph->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->SetTitle("Efficiency (arbitrary units)") ;
  graph->GetYaxis()->CenterTitle();
	  
  graph->SetMinimum(0.01) ; 
  graph->SetMaximum(1.0)  ; 

  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;
}


TGraph* Detector::GetGraphImpactParam(Int_t mode, Int_t axis, Int_t color, Int_t linewidth) {
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
  graph->GetYaxis()->SetTitle("d_{0} r#phi resolution [#mum]") ;
  
  graph->SetMinimum(1) ;
  graph->SetMaximum(300.1) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetNoExponent(1) ;
  graph->GetXaxis()->SetMoreLogLabels(1) ;
  graph->GetYaxis()->CenterTitle();
  
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(linewidth);

  return graph;

}

TGraph* Detector::GetGraph(Int_t number, Int_t color, Int_t linewidth) {
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
  default:
    printf(" Error: chosen graph number not valid\n");
  }
  return 0;

}
 

void Detector::MakeAliceCurrent(Int_t AlignResiduals, Bool_t flagTPC) {

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


void Detector::MakeStandardPlots(Bool_t add, Int_t color, Int_t linewidth,Bool_t onlyPionEff) {
  //
  // Produces the standard performace plots
  //
 
  if (!add) {

    TCanvas *c1 = new TCanvas("c1","c1");  
    c1->Divide(2,2);
    
    c1->cd(1);  gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    TGraph *eff = GetGraphRecoEfficiency(1,color,linewidth);
    eff->SetTitle("Efficiencies");
    eff->Draw("AC");
    if (!onlyPionEff) {
      GetGraphRecoEfficiency(2,color,linewidth)->Draw("C");
      GetGraphRecoEfficiency(3,color,linewidth)->Draw("C");
    }
    c1->cd(2); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogy();  gPad->SetLogx(); 
    GetGraphMomentumResolution(color,linewidth)->Draw("AC");
    
    c1->cd(3); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    GetGraphPointingResolution(0,color,linewidth)->Draw("AC");
    
    c1->cd(4); gPad->SetGridx();   gPad->SetGridy(); 
    gPad->SetLogx(); 
    GetGraphPointingResolution(1,color,linewidth)->Draw("AC");

  } else {

    TVirtualPad *c1 = gPad->GetMother();

    c1->cd(1);
    GetGraphRecoEfficiency(1,color,linewidth)->Draw("C");
    if (!onlyPionEff) {
      GetGraphRecoEfficiency(2,color,linewidth)->Draw("C");
      GetGraphRecoEfficiency(3,color,linewidth)->Draw("C");
    }
    c1->cd(2); GetGraphMomentumResolution(color,linewidth)->Draw("C");
    
    c1->cd(3); GetGraphPointingResolution(0,color,linewidth)->Draw("C");
    
    c1->cd(4); GetGraphPointingResolution(1,color,linewidth)->Draw("C");
    
  }

}
