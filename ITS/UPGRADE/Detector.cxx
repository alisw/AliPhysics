#include "Detector.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFormula.h>

/***********************************************************

Fast Simulation tool for Inner Tracker Systems

original code of using the billoir technique was developed
for the HFT (STAR), James H. Thomas, jhthomas@lbl.gov
http://rnc.lbl.gov/~jhthomas

Changes by S. Rossegger
Dec. 2010 - Translation into C++ class format 
          - Adding various Setters and Getters to build the geometry 
	    (based on cylinders) plus handling of the layer properties 

***********************************************************/


#define RIDICULOUS 99999 // A ridiculously large resolution (cm) to flag a dead detector

#define Luminosity    1.e27       // Luminosity of the beam (LHC HI == 1.e27, RHIC II == 8.e27 )
#define SigmaD        6.0         // Size of the interaction diamond (cm) (LHC = 6.0 cm)
#define dNdEtaMinB    950         // Multiplicity per unit Eta  (AuAu MinBias = 170, Central = 700)
#define dNdEtaCent    2300        // Multiplicity per unit Eta  (LHC at 5.5 TeV not known)

#define CrossSectionMinB  8      // minB Cross section for event under study (PbPb MinBias ~ 8 Barns)
#define AcceptanceOfTpcAndSi     1//0.60 //0.35  // Assumed geometric acceptance (efficiency) of the TPC and Si detectors
#define UPCBackgroundMultiplier  1.0   // Increase multiplicity in detector (0.0 to 1.0 * UPCRate ) (eg 1.0)
#define OtherBackground          0.0   // Increase multiplicity in detector (0.0 to 1.0 * minBias)  (eg 0.0)
#define EfficiencySearchFlag     1     // Define search method. ChiSquare = 1, Simple = 0.  ChiSq is better

#define PionMass                 0.139  // Mass of the Pion
#define KaonMass                 0.498  // Mass of the Kaon
#define D0Mass                   1.865  // Mass of the D0



class CylLayer : public TNamed {
public:

  CylLayer(char *name) : TNamed(name,name) {}

  Float_t GetRadius() {return radius;}
  Float_t GetRadL() {return radL;}
  Float_t GetPhiRes() {return phiRes;}
  Float_t GetZRes() {return zRes;}

  //  void Print() {printf("  r=%3.1lf X0=%1.6lf sigPhi=%1.4lf sigZ=%1.4lf\n",radius,radL,phiRes,zRes); }
  Float_t radius; Float_t radL; Float_t phiRes; Float_t zRes;   
  Bool_t isDead;

 ClassDef(CylLayer,1);
};

ClassImp(Detector)
Detector::Detector() 
  : TNamed("test_detector","detector"),
    numberOfLayers(0),
    fBField(0.5),
    fLhcUPCscale(2.5),    
    fIntegrationTime(0.02), // in ms
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140)    // Standard: pion mass 
{
  //
  // default constructor
  //
  //  fLayers = new TObjArray();
}

Detector::Detector(char *name, char *title)
  : TNamed(name,title),
    numberOfLayers(0),
    fBField(0.5),
    fLhcUPCscale(2.5),
    fIntegrationTime(0.02),  // in ms
    fAvgRapidity(0.45),      // Avg rapidity, MCS calc is a function of crossing angle
    fParticleMass(0.140)     // Standard: pion mass
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

void Detector::AddLayer(char *name, Float_t radius, Float_t radL, Float_t phiRes, Float_t zRes) {
  
  CylLayer *newLayer = (CylLayer*) fLayers.FindObject(name);

  if (!newLayer) {
    newLayer = new CylLayer(name);
    newLayer->radius = radius;
    newLayer->radL = radL;
    newLayer->phiRes = phiRes;
    newLayer->zRes = zRes;

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
    numberOfLayers += 1;

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
     tmp->phiRes = 99999;
     tmp->zRes = 99999;
     tmp->isDead = kTRUE;
  }
}

void Detector::SetRadius(char *name, Float_t radius) {
  //
  // Set layer radius [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
 

  if (!tmp) {
    printf("Layer %s not found - cannot set resolution\n",name);
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
    printf("Layer %s not found - cannot set resolution\n",name);
  else 
    return tmp->radL;

  return 0;
}

void Detector::SetRadiationLength(char *name, Float_t radL) {
  //
  // Set layer radius [cm]
  //

  CylLayer *tmp = (CylLayer*) fLayers.FindObject(name);
  if (!tmp) 
    printf("Layer %s not found - cannot set resolution\n",name);
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
    printf("Layer %s not found - cannot set resolution\n",name);
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
    tmp->phiRes = phiRes;
    tmp->zRes = zRes;
    if (zRes==RIDICULOUS && phiRes==RIDICULOUS) 
      tmp->isDead = kTRUE;
    else 
      tmp->isDead = kFALSE;
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
    printf("Layer %s not found - cannot set resolution\n",name);
  else {
    if (axis==0) return tmp->phiRes;
    if (axis==1) return tmp->zRes;
    printf("error: axis must be either 0 or 1 (rphi or z axis)\n");
  }
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
    fLayers.Remove(tmp);
    numberOfLayers -= 1;
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

void Detector::AddTPC(Float_t phiResMean, Float_t zResMean, Int_t skip) {
  //
  // Emulates the TPC
  // 
  // skip=1: Use every padrow, skip=2: Signal in every 2nd padrow 
  
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

  // flag as dead, although resolution is ok ... makes live easier in the prints ... ;-)
  CylLayer *tmp = 0;
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
    tmp = (CylLayer*)fLayers.At(i);  
    TString name(tmp->GetName());
    if (name.Contains("tpc") && (!name.Contains("tpc_0")) ) tmp->isDead=kTRUE;
  }


 
}


Double_t Detector::ThetaMCS ( Double_t mass, Double_t RadLength, Double_t momentum ) 
{
  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum ) * (1+0.038*TMath::Log(RadLength)) ;
  return (theta) ;
}


Double_t Detector::ProbGoodHit ( Double_t Radius, Double_t SearchRadiusRPhi, Double_t SearchRadiusZ ) 
{
  // Based on work by Howard Wieman: http://rnc.lbl.gov/~wieman/GhostTracks.htm 
  // and http://rnc.lbl.gov/~wieman/HitFinding2D.htm
  // This is the probability of getting a good hit using 2D Gaussian distribution function and infinite search radius
  Double_t Sx, Sy, GoodHit ;
  Sx = 2 * TMath::Pi() *  SearchRadiusRPhi * SearchRadiusRPhi * HitDensity(Radius) ;
  Sy = 2 * TMath::Pi() *  SearchRadiusZ    * SearchRadiusZ    * HitDensity(Radius) ;
  GoodHit =  TMath::Sqrt(1./((1+Sx)*(1+Sy)))  ;
  return ( GoodHit ) ;
}


Double_t Detector::ProbGoodChiSqHit ( Double_t Radius, Double_t SearchRadiusRPhi, Double_t SearchRadiusZ ) 
{
  // Based on work by Victor Perevoztchikov and Howard Wieman: http://rnc.lbl.gov/~wieman/HitFinding2DXsq.htm
  // This is the probability of getting a good hit using a Chi**2 search on a 2D Gaussian distribution function
  Double_t Sx, GoodHit ;
  Sx = 2 * TMath::Pi() *  SearchRadiusRPhi * SearchRadiusZ * HitDensity(Radius) ;
  GoodHit =  1./(1+Sx) ;
  return ( GoodHit ) ;  
}


Double_t Detector::HitDensity ( Double_t Radius ) 
{
  // Background (0-1) is included via 'OtherBackground' which multiplies the minBias rate by a scale factor.
  // UPC electrons is a temporary kludge that is based on Kai Schweda's summary of Kai Hainken's MC results
  // See K. Hencken et al. PRC 69, 054902 (2004) and PPT slides by Kai Schweda.
  // Note that this function assumes we are working in CM and CM**2 [not meters].
  // Based on work by Yan Lu 12/20/2006, all radii and densities in centimeters or cm**2.

  Double_t MaxRadiusSlowDet = 10; //?   // Maximum radius for slow detectors.  Fast detectors 
                                        // and only fast detectors reside outside this radius.
  Double_t ArealDensity = 0 ;

  if ( Radius > MaxRadiusSlowDet ) 
    {
      ArealDensity  = OneEventHitDensity(dNdEtaCent,Radius)  ; // Fast detectors see central collision density (only)
      ArealDensity += OtherBackground*OneEventHitDensity(dNdEtaMinB,Radius)  ;  // Increase density due to background 
    }

  if (Radius < MaxRadiusSlowDet )
    { // Note that IntegratedHitDensity will always be minB one event, or more, even if integration time => zero.
      ArealDensity  = OneEventHitDensity(dNdEtaCent,Radius) 
	            + IntegratedHitDensity(dNdEtaMinB,Radius) 
	            + UpcHitDensity(Radius) ;
      ArealDensity += OtherBackground*IntegratedHitDensity(dNdEtaMinB,Radius) ;  
      // Increase density due to background 
    } 

  return ( ArealDensity ) ;  
}


double Detector::OneEventHitDensity( Double_t Multiplicity, Double_t Radius )
{
  // This is for one event at the vertex.  No smearing.
  double den   = Multiplicity / (2.*TMath::Pi()*Radius*Radius) ; // 2 eta ?
  // note: surface of sphere is  '4*pi*r^2'
  //       surface of cylinder is '2*pi*r* h' 
  return den ;
} 


double Detector::IntegratedHitDensity(Double_t Multiplicity, Double_t Radius)
{ 
  // The integral of minBias events smeared over a gaussian vertex distribution.
  // Based on work by Yan Lu 12/20/2006, all radii in centimeters.

  Double_t ZdcHz = Luminosity * 1.e-24 * CrossSectionMinB ;
  Double_t den   = ZdcHz * fIntegrationTime/1000. * Multiplicity * Dist(0., Radius) / (2.*TMath::Pi()*Radius) ;

  // Note that we do not allow the rate*time calculation to fall below one minB event at the vertex.
  if ( den < OneEventHitDensity(Multiplicity,Radius) )  den = OneEventHitDensity(Multiplicity,Radius) ;  

  return den ;
} 


double Detector::UpcHitDensity(Double_t Radius)
{ 
  // QED electrons ...

  Double_t UPCelectrons ;                                 ;  
  UPCelectrons =  fLhcUPCscale * (1.23 - Radius/6.5)      ;  // Fit to Kai Schweda summary tables at RHIC * 'scale' for LHC
  if ( UPCelectrons < 0 ) UPCelectrons =  0.0             ;  // UPC electrons fall off quickly and don't go to large R
  UPCelectrons *= IntegratedHitDensity(dNdEtaMinB,Radius) ;  // UPCs increase Mulitiplicty ~ proportional to MinBias rate
  UPCelectrons *= UPCBackgroundMultiplier                 ;  // Allow for an external multiplier (eg 0-1) to turn off UPC

  return UPCelectrons ;
} 


double Detector::Dist(double z, double r)
{
  // Convolute dEta/dZ  distribution with assumed Gaussian of vertex z distribution
  // Based on work by Howard Wieman http://rnc.lbl.gov/~wieman/HitDensityMeasuredLuminosity7.htm
  // Based on work by Yan Lu 12/20/2006, all radii and Z location in centimeters.
  Int_t    Index  =  1     ;     // Start weight at 1 for Simpsons rule integration
  Int_t    NSTEPS =  301   ;     // NSteps must be odd for Simpson's rule to work
  double   dist   =  0.0   ;
  double   dz0    =  ( 4*SigmaD - (-4)*SigmaD ) / (NSTEPS-1)  ;  //cm
  double    z0    =  0.0   ;     //cm
  for(int i=0; i<NSTEPS; i++){
    if ( i == NSTEPS-1 ) Index = 1 ;
    z0 = -4*SigmaD + i*dz0 ;
    dist += Index * (dz0/3.) * (1/sqrt(2.*TMath::Pi())/SigmaD) * exp(-z0*z0/2./SigmaD/SigmaD) * 
      (1/sqrt((z-z0)*(z-z0) + r*r)) ;
    if ( Index != 4 ) Index = 4; else Index = 2 ;
  }
  return dist; 
}

#define  PZero   0.861  // Momentum of back to back decay particles in the CM frame
#define  EPiZero 0.872  // Energy of the pion from a D0 decay at rest
#define  EKZero  0.993  // Energy of the Kaon from a D0 decay at rest

Double_t Detector::D0IntegratedEfficiency( Double_t pt, Double_t corrEfficiency[][400] ) {
  // Math from Ron Longacre.  Note hardwired energy to bin conversion for PtK and PtPi.

  Double_t Const1  =  pt / D0Mass ;
  Double_t Const2  =  TMath::Sqrt(pt*pt+D0Mass*D0Mass) / D0Mass ;
  Double_t sum, PtPi, PtK ;
  Double_t effp, effk ;

  sum = 0.0 ;
  for ( Int_t k = 0 ; k < 360 ; k++ )   {
    
    Double_t theta = k * TMath::Pi() / 180. ;
    
    PtPi = TMath::Sqrt( 
		       PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*Const2*Const2 +
		       Const1*Const1*EPiZero*EPiZero -
		       2*PZero*TMath::Cos(theta)*Const2*Const1*EPiZero +
		       PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		       ) ;
    
    PtK = TMath::Sqrt( 
		      PZero*PZero*TMath::Cos(theta)*TMath::Cos(theta)*Const2*Const2 +
		      Const1*Const1*EKZero*EKZero +
		      2*PZero*TMath::Cos(theta)*Const2*Const1*EKZero +
		      PZero*PZero*TMath::Sin(theta)*TMath::Sin(theta)
		      ) ;

    // JT Test Remove 100 MeV/c in pt to simulate eta!=0 decays
    Int_t pionindex = (int)((PtPi-0.1)*100.0 - 65.0*TMath::Abs(fBField)) ; 
    Int_t kaonindex = (int)((PtK -0.1)*100.0 - 65.0*TMath::Abs(fBField)) ; 
      
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




void Detector::SolveViaBilloir(Int_t flagD0,Int_t print) {

  // Calculate track parameters using Billoirs method of matrices

  
  Double_t pt, pz, lambda, momentum, rho, DeltaPoverP  ;
  Double_t layerThickness, MCS, mmm, Charge ;
  Double_t F1, F2, F3, G1, G2, G3 ;
  Double_t dx, sigma2, sigma2Z ;
  Double_t mass[3] ;
  Int_t PrintOnce = 1 ;

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

    for ( Int_t i = 0 ; i < 400 ; i++ ) { // pt loop
 
      CylLayer *last = (CylLayer*) fLayers.At((fLayers.GetEntries()-1));

      // Convert to Meters, Tesla, and GeV
      xpoint[i] =  ( 0.3*(last->radius/100)*TMath::Abs(fBField) ) 
	- 0.08 + TMath::Power(10,2.3*i/400) / 10.0 ; // Starting values based on radius in TPC ... log10 steps to ~20 GeV
   
   
      
      // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
      // These are the EndPoint values for y, z, a, b, and d
      pt  =  xpoint[i]                       ;  // GeV/c
      rho =  pt / TMath::Abs( 0.3 * fBField );  // Radius of curvature of the track (meters)
      momentum = pt / TMath::Cos(lambda)     ;  // Total momentum  
      Charge   = -1                          ;  // Assume an electron 
      pz  =  pt * TMath::Tan(lambda)         ;  // GeV/c
      d   =  0.3 * Charge / momentum         ;  // Its an electrons so q = -1, which makes d a negative number

      //      printf("%d - pt %lf r%lf | %lf %lf\n",massloop,xpoint[i],(last->radius)/100,momentum, d);

      // Matrices for Billoir
      TMatrixD BigA(2,2); TMatrixD LittleA(2,2);
      TMatrixD Istar(5,5); TMatrixD D(5,5); TMatrixD DT(5,5); TMatrixD A(5,5); TMatrixD M(5,5); TMatrixD Work(5,5) ;
    
      // Set Detector-Efficiency Storage area to unity
      efficiency[massloop][i] = 1.0 ;
      // Back-propagate the covariance matrix along the track.  I is the inverse of the covariance matrix. 
      // Start with small variations to a straight line 'information' matrix.  Tuning this matrix is 
      // equivalent to 'starting the recursion' in Billoirs paper.  It is a tricky business.  You need test
      // data and many layers for the matrices to stabilize.  Do not believe results for a small number of layers.
      // In our case, always include the TPC and this matrix will be well conditioned before it gets to the Si.

      Istar.UnitMatrix(); // start with unity

      // 'A' is the "angle" matrix for MCS at each step.
      // 'D' is the "distance" matrix at each step.  It propagates the particle backwards along the track. 
      // 'M' is the "measurement" matrix. It represents the measurement resolution at each step.
    
      CylLayer *layer = 0;
      CylLayer *nextlayer = 0;

      for (Int_t j=(fLayers.GetEntries()-1); j>0; j--) {  // Layer loop

	layer = (CylLayer*)fLayers.At(j);
	nextlayer = (CylLayer*)fLayers.At(j-1);

	//	layer->Print();
	// Convert to Meters, Tesla, and GeV
	Float_t radius = layer->radius/100;
	//	Float_t phiResolution = layer->GetPhiRes()/100;
	//	Float_t zResolution = layer->GetZRes()/100;
	Float_t radLength = layer->radL;

	Float_t nextRadius = nextlayer->radius /100;
	Float_t nextPhiResolution = nextlayer->phiRes /100;
	Float_t nextZResolution = nextlayer->zRes /100;
     
	//	printf("%d: %lf %lf %lf %lf\n",j,radius,nextRadius, radLength,nextPhiResolution);
	//	printf("%d %d %d : %lf %lf\n",j,fLayers.GetEntries()-1,0, rho, radius);
		
	if ( radius >= rho )    // (meters) protect agains too low pt
	  { printf("pt lower bound is too low for this algorithm\n");  return ; }
	
	y   =  rho - TMath::Sqrt(rho*rho-radius*radius)  ; // These are 'ideal' locations and
	a   =  radius / ( rho - y )                      ; // not propagated values which should  
	b   =  rho * TMath::Tan(lambda) / ( rho - y )    ; // be done if we had data. But we don't.
	z   =  rho * TMath::Tan(lambda) * TMath::ATan(a) ; 
	e   =  TMath::Sqrt( 1 + a*a + b*b )              ; 
      
	layerThickness = radLength / TMath::Sin(TMath::Pi()/2 - lambda)  ; 
      
	MCS =  ThetaMCS(mass[massloop], layerThickness, momentum)           ; 
	mmm =  ( 1 + a*a + b*b ) ; 
	BigA(0,0) = mmm*(1+a*a) ; BigA(0,1) = mmm*a*b     ; LittleA(0,0) = Istar(2,2)  ; LittleA(0,1) = Istar(2,3) ; 
	BigA(1,0) = mmm*a*b     ; BigA(1,1) = mmm*(1+b*b) ; LittleA(1,0) = Istar(3,2)  ; LittleA(1,1) = Istar(3,3) ; 
	LittleA = BigA.Invert() + MCS*MCS*LittleA ;
	LittleA = LittleA.Invert() ;
	A(0,0) = 0.0 ;  A(0,1) = 0.0  ;  A(0,2) = 0.0          ;  A(0,3) = 0.0          ;  A(0,4) = 0.0  ; 
	A(1,0) = 0.0 ;  A(1,1) = 0.0  ;  A(1,2) = 0.0          ;  A(1,3) = 0.0          ;  A(1,4) = 0.0  ;
	A(2,0) = 0.0 ;  A(2,1) = 0.0  ;  A(2,2) = LittleA(0,0) ;  A(2,3) = LittleA(0,1) ;  A(2,4) = 0.0  ;
	A(3,0) = 0.0 ;  A(3,1) = 0.0  ;  A(3,2) = LittleA(1,0) ;  A(3,3) = LittleA(1,1) ;  A(3,4) = 0.0  ;
	A(4,0) = 0.0 ;  A(4,1) = 0.0  ;  A(4,2) = 0.0          ;  A(4,3) = 0.0          ;  A(4,4) = 0.0  ;
	Istar = Istar - MCS*MCS*Istar*A*Istar ;
	dx     = ( radius - nextRadius )         ;        // Billoir works with absolute magnitude of distance
	F3     = e * ( -1 * ( 1 + a*a ) * fBField )   ;        // Assume fBField is on Z axis, only.
	F1     = d * ( a*F3/(e*e) - 2*e*a*fBField )   ;
	F2     = d * ( b*F3/(e*e) )                  ;
	G3     = e * ( -1 * a * b * fBField )         ;
	G1     = d * ( a*G3/(e*e) + e*b*fBField )     ;
	G2     = d * ( b*G3/(e*e) - e*a*fBField )     ;	      
	D(0,0) = 1.0 ; D(0,1) = 0.0 ; D(0,2) = -1*dx+F1*dx*dx/2 ; D(0,3) = F2*dx*dx/2  ;  D(0,4) = F3*dx*dx/2 ;
	D(1,0) = 0.0 ; D(1,1) = 1.0 ; D(1,2) = G1*dx*dx/2 ;  D(1,3) = -1*dx+G2*dx*dx/2 ;  D(1,4) = G3*dx*dx/2 ;
	D(2,0) = 0.0 ; D(2,1) = 0.0 ; D(2,2) = 1-F1*dx    ;  D(2,3) = -1*F2*dx         ;  D(2,4) = -1*F3*dx   ;
	D(3,0) = 0.0 ; D(3,1) = 0.0 ; D(3,2) = -1*G1*dx   ;  D(3,3) = 1-G2*dx          ;  D(3,4) = -1*G3*dx   ;
	D(4,0) = 0.0 ; D(4,1) = 0.0 ; D(4,2) = 0.0        ;  D(4,3) = 0.0              ;  D(4,4) = 1.0        ;
	DT.Transpose(D) ;
	Istar  = DT.Invert() * Istar * D.Invert() ;
	// Prepare to save Detector efficiencies
	Work = Istar  ;          // Working copy of matrix 
	Work.Invert() ;          // Invert the Matrix to recover the convariance matrix
	DetPointRes [j-1][i]     =  TMath::Sqrt( Work(0,0) )  ;     // result in meters
	DetPointZRes[j-1][i]     =  TMath::Sqrt( Work(1,1) )  ;     // result in meters
	// End save
	sigma2  = ( nextPhiResolution  * nextPhiResolution  )   ; 
	sigma2Z = ( nextZResolution    * nextZResolution )   ;
	M(0,0) = 1/sigma2 ;  M(0,1) = 0.0  ;        M(0,2) = 0.0  ;  M(0,3) = 0.0  ;  M(0,4) = 0.0  ; 
	M(1,0) = 0.0 ;       M(1,1) = 1/sigma2Z  ;  M(1,2) = 0.0  ;  M(1,3) = 0.0  ;  M(1,4) = 0.0  ;
	M(2,0) = 0.0 ;       M(2,1) = 0.0  ;        M(2,2) = 0.0  ;  M(2,3) = 0.0  ;  M(2,4) = 0.0  ;
	M(3,0) = 0.0 ;       M(3,1) = 0.0  ;        M(3,2) = 0.0  ;  M(3,3) = 0.0  ;  M(3,4) = 0.0  ;
	M(4,0) = 0.0 ;       M(4,1) = 0.0  ;        M(4,2) = 0.0  ;  M(4,3) = 0.0  ;  M(4,4) = 0.0  ;
	Istar = Istar + M ;
	
      }
      
      // Invert the Matrix to recover the convariance matrix
      Istar.Invert() ;
      // Convert the Convariance matrix parameters into physical quantities
      // The results are propogated to the previous point but *do not* include the measurement at that point.
      DeltaPoverP    =  TMath::Sqrt( Istar(4,4) ) * momentum / 0.3  ;  // Absolute magnitude so ignore Charge
      ypoint[i]      =  100.* TMath::Abs( DeltaPoverP )             ;  // results in percent
      yprime[i]      =  TMath::Sqrt( Istar(0,0) ) * 1.e6            ;  // result in microns
      yprimeZ[i]     =  TMath::Sqrt( Istar(1,1) ) * 1.e6            ;  // result in microns
      equivalent[i]  =  TMath::Sqrt(yprime[i]*yprimeZ[i])           ;  // Equivalent circular radius
      
      
      if (print == 1 && xpoint[i] > 0.750 && massloop == 2 && PrintOnce == 1) {
	printf("Mass of tracked particle: %f\n",fParticleMass) ;
	printf("Radius Thickness PointResOn PointResOnZ  DetRes  DetResZ  Density Efficiency\n") ;
	//	PrintOnce =0;
      }

      // print and efficiency
      for (Int_t j=(fLayers.GetEntries()-1); j>0; j--) {  // Layer loop
	
	layer = (CylLayer*)fLayers.At(j-1);
	
	// Convert to Meters, Tesla, and GeV
	Float_t radius = layer->radius /100;
	Float_t phiRes = layer->phiRes /100;
	Float_t zRes = layer->zRes /100;
	Float_t radLength = layer->radL;
	Bool_t isDead = layer->isDead;

	if ( !isDead && radLength >0 )  { 
	    Double_t rphiError  =  TMath::Sqrt( DetPointRes[j-1][i] * DetPointRes [j-1][i] + 
						phiRes * phiRes ) * 100.  ; // work in cm
	    Double_t zError     =  TMath::Sqrt( DetPointZRes[j-1][i] * DetPointZRes[j-1][i] +
						zRes * zRes ) * 100.  ; // work in cm

	    if (print == 1 && xpoint[i] > 0.750 && massloop == 2 && PrintOnce == 1) {
	      printf("%5.1f %9.4f %10.0f %11.0f %7.0f %8.0f %8.2f %10.3f\n",
		     radius*100, radLength, 
		     DetPointRes[j-1][i]*1.e6, DetPointZRes[j-1][i]*1.e6,
		     phiRes*1.e6, zRes*1.e6,
		     HitDensity(radius*100),
		     ProbGoodChiSqHit( radius*100, rphiError , zError ) ) ;
	    }

	    if ( EfficiencySearchFlag == 0 )
	      efficiency[massloop][i]    *=  ProbGoodHit( radius*100, rphiError , zError  ) ;
	    else if ( EfficiencySearchFlag == 1 )
	      efficiency[massloop][i]    *=  ProbGoodChiSqHit( radius*100, rphiError , zError  ) ;
	}
      }
      if (print == 1 && xpoint[i] > 0.750 && massloop == 2 && PrintOnce == 1) {
	PrintOnce = 0 ;
	printf("\n")  ;
      }
      
    } // mass loop
  } // pt loop
  
}


TGraph * Detector::GetGraphMomentumResolution(Int_t color, Int_t linewidth) {
  
  TGraph *graph = new TGraph(400, xpoint, ypoint);
  graph->SetMaximum(10) ;
  graph->SetMinimum(0) ;
  graph->SetTitle("Momentum Resolution .vs. Pt" ) ;
  graph->GetXaxis()->SetRangeUser(0.,5.0) ;
  graph->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]") ;
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->SetTitle("Momentum Resolution [%]") ;
  graph->GetYaxis()->CenterTitle();
  
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
    graph = new TGraph ( 400, xpoint, yprime ) ;
    graph->SetTitle("R-#phi Pointing Resolution .vs. Pt" ) ;
    graph->GetYaxis()->SetTitle("R-#phi Pointing Resolution [#mum]") ;
  } else {
    graph =  new TGraph ( 400, xpoint, yprimeZ ) ;
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
  printf("Telescope equation for layers:\n");
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
      printf("  %s \n",l->GetName());
      count++;
    }
    if (count>=2) break;	
  }
   

  Double_t pt, momentum, thickness,MCS ;
  Double_t lambda = TMath::Pi()/2.0 - 2.0*TMath::ATan(TMath::Exp(-1*fAvgRapidity)); 

  for ( Int_t i = 0 ; i < 400 ; i++ ) { 
    // Reference data as if first two layers were acting all alone 
    pt  =  xpoint[i]  ;
    momentum = pt / TMath::Cos(lambda)   ;  // Total momentum
    resolution[i] =  layerResolution[0]*layerResolution[0]*layerRadius[1]*layerRadius[1] 
      +  layerResolution[1]*layerResolution[1]*layerRadius[0]*layerRadius[0] ;
    resolution[i] /= ( layerRadius[1] - layerRadius[0] ) * ( layerRadius[1] - layerRadius[0] ) ;
    thickness = layerThickness[0] / TMath::Sin(TMath::Pi()/2 - lambda) ;
    MCS = ThetaMCS(fParticleMass, thickness, momentum) ;
    resolution[i] += layerRadius[0]*layerRadius[0]*MCS*MCS ;
    resolution[i] =  TMath::Sqrt(resolution[i]) * 10000.0 ;  // result in microns
  }



  TGraph* graph = new TGraph ( 400, xpoint, resolution ) ;
   
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
  Double_t kaonEfficiency[400], pionEfficiency[400], D0efficiency[400]; 
  Double_t partEfficiency[2][400];
  
  if (particle != 0) {
    // resulting Pion and Kaon efficiency scaled with overall efficiency
    Double_t doNotDecayFactor;
    for ( Int_t massloop = 0 ; massloop < 2 ; massloop++) { //0-pion, 1-kaon
      
      for ( Int_t j = 0 ; j < 400 ; j++ ) { 
	// JT Test Let the kaon decay.  If it decays inside the TPC ... then it is gone; for all decays < 130 cm.
	Double_t momentum = xpoint[j] / TMath::Cos(lambda)           ;  // Total momentum at average rapidity
	if ( massloop == 1 ) { // KAON
	  doNotDecayFactor  = TMath::Exp(-130/(371*momentum/KaonMass)) ;  // Decay length for kaon is 371 cm.
	  kaonEfficiency[j] = efficiency[1][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;
	} else { // PION
	  doNotDecayFactor = 1.0 ;
	  pionEfficiency[j] = efficiency[0][j] * AcceptanceOfTpcAndSi*doNotDecayFactor ;	
	}
	partEfficiency[0][j] = pionEfficiency[j];
	partEfficiency[1][j] = kaonEfficiency[j];
      }      
    }
    
    // resulting estimate of the D0 efficiency
    for ( Int_t j = 0 ; j < 400 ; j++ ) {
      D0efficiency[j] = D0IntegratedEfficiency(xpoint[j],partEfficiency);
    }
  } else { 
    for ( Int_t j = 0 ; j < 400 ; j++ ) { 
      particleEfficiency[j] = efficiency[2][j]* AcceptanceOfTpcAndSi;
      // NOTE: Decay factor (see kaon) should be included to be realiable
    }
  }

  TGraph * graph =  0;
  if (particle==0) {
    graph = new TGraph ( 400, xpoint, particleEfficiency ) ; // choosen mass
    graph->SetLineWidth(1);
  }  else if (particle==1) {
    graph = new TGraph ( 400, xpoint, pionEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==2) {
    graph = new TGraph ( 400, xpoint, kaonEfficiency ) ;
    graph->SetLineWidth(1);
  }  else if (particle ==3) {
    graph = new TGraph ( 400, xpoint, D0efficiency ) ;
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
  TFormula vtxResZ("vtxResZ","600/(x+4)+10"); // 
    
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
  AddLayer((char*)"IFC",   77.8,0.01367); // Inner Field cage

  if (flagTPC) AddTPC(0.1,0.1);             // TPC

  // Adding the ITS - current configuration
  
  if (AlignResiduals==0) {

    AddLayer((char*)"spd1", 3.9, 0.0114, 0.0012, 0.0100);
    AddLayer((char*)"spd2", 7.6, 0.0114, 0.0012, 0.0100);
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
	     TMath::Sqrt(0.0100*0.0100+0.0050*0.0050));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030),
	     TMath::Sqrt(0.0100*0.0100+0.0050*0.0050));
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
	     TMath::Sqrt(0.0100*0.0100+0.0050*0.0050));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0030*0.0030+0.0008*0.0008),
	     TMath::Sqrt(0.0100*0.0100+0.0050*0.0050));
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
	     TMath::Sqrt(0.0100*0.0100+0.000*0.000));
    AddLayer((char*)"spd2", 7.6, 0.0114, TMath::Sqrt(0.0012*0.0012+0.0008*0.0008),
	     TMath::Sqrt(0.0100*0.0100+0.000*0.000));
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
