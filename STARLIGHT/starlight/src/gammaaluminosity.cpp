///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 264                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-06-06 22:05:12 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "inputParameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "gammaaluminosity.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
photonNucleusLuminosity::photonNucleusLuminosity(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(inputParametersInstance, bbsystem)
  ,_ptBinWidthInterference(inputParametersInstance.ptBinWidthInterference())
  ,_interferenceStrength(inputParametersInstance.interferenceStrength())
  ,_protonEnergy(inputParametersInstance.protonEnergy())
  ,_beamLorentzGamma(inputParametersInstance.beamLorentzGamma())
  ,_baseFileName(inputParametersInstance.baseFileName())
  ,_maxW(inputParametersInstance.maxW())
  ,_minW(inputParametersInstance.minW())
  ,_nmbWBins(inputParametersInstance.nmbWBins())
  ,_maxRapidity(inputParametersInstance.maxRapidity())
  ,_nmbRapidityBins(inputParametersInstance.nmbRapidityBins())
  ,_productionMode(inputParametersInstance.productionMode())
  ,_beamBreakupMode(inputParametersInstance.beamBreakupMode())
  ,_interferenceEnabled(inputParametersInstance.interferenceEnabled())
  ,_maxPtInterference(inputParametersInstance.maxPtInterference())
  ,_nmbPtBinsInterference(inputParametersInstance.nmbPtBinsInterference())
{
  cout <<"Creating Luminosity Tables."<<endl;
  photonNucleusDifferentialLuminosity();
  cout <<"Luminosity Tables created."<<endl;
}


//______________________________________________________________________________
photonNucleusLuminosity::~photonNucleusLuminosity()
{ }


//______________________________________________________________________________
void photonNucleusLuminosity::photonNucleusDifferentialLuminosity()
{
  double W,dW,dY;
  double Egamma,Y;
  double testint,dndWdY;
  double csgA;
  double C;  
  int beam; 

  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);
  
  double  bwnorm,Eth;

  dW = (_wMax-_wMin)/_nWbins;
  dY  = (_yMax-(-1.0)*_yMax)/_nYbins;
    
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open(wyFileName.c_str());
  wylumfile << getbbs().beam1().Z() <<endl;
  wylumfile << getbbs().beam1().A() <<endl;
  wylumfile << getbbs().beam2().Z() <<endl;
  wylumfile << getbbs().beam2().A() <<endl;
  wylumfile << _beamLorentzGamma <<endl;
  wylumfile << _maxW <<endl;
  wylumfile << _minW <<endl;
  wylumfile << _nmbWBins <<endl;
  wylumfile << _maxRapidity <<endl;
  wylumfile << _nmbRapidityBins <<endl;
  wylumfile << _productionMode <<endl;
  wylumfile << _beamBreakupMode <<endl;
  wylumfile << _interferenceEnabled <<endl;
  wylumfile << _interferenceStrength <<endl;
  wylumfile << starlightConstants::deuteronSlopePar <<endl;
  wylumfile << _maxPtInterference <<endl;
  wylumfile << _nmbPtBinsInterference <<endl;
  
  //     Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  //Grabbing default value for C in the breit-wigner calculation
  C=getDefaultC();
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {
    W = _wMin + double(i)*dW + 0.5*dW;
    testint = testint + breitWigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  //     Write the values of Y used in the calculation to slight.txt.
  for(unsigned int i = 0; i <= _nYbins - 1; ++i) {
    Y = -1.0*_yMax + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }
    
  Eth=0.5*(((_wMin+starlightConstants::protonMass)*(_wMin +starlightConstants::protonMass)-starlightConstants::protonMass*starlightConstants::protonMass)/(_protonEnergy+sqrt(_protonEnergy*_protonEnergy-starlightConstants::protonMass*starlightConstants::protonMass))); 

  int A_1 = getbbs().beam1().A(); 
  int A_2 = getbbs().beam2().A();

  // Do this first for the case when the first beam is the photon emitter 
  // Treat pA separately with defined beams 
  // The variable beam (=1,2) defines which nucleus is the target 
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {

    W = _wMin + double(i)*dW + 0.5*dW;
    
    for(unsigned int j = 0; j <= _nYbins - 1; ++j) { 

      Y = -1.0*_yMax + double(j)*dY + 0.5*dY;

      if( A_2 == 1 && A_1 != 1 ){
        // pA, first beam is the nucleus and is in this case the target  
        Egamma = 0.5*W*exp(-Y); 
        beam = 1; 
      } else if( A_1 ==1 && A_2 != 1){
        // pA, second beam is the nucleus and is in this case the target 
        Egamma = 0.5*W*exp(Y); 
        beam = 2; 
      } else {
        Egamma = 0.5*W*exp(Y);        
        beam = 2; 
      }

      dndWdY = 0.; 

      if( Egamma > Eth && Egamma < maxPhotonEnergy() ){

	csgA=getcsgA(Egamma,W,beam);
        dndWdY = Egamma*photonFlux(Egamma,beam)*csgA*breitWigner(W,bwnorm);

      }

      wylumfile << dndWdY << endl;

    }
  }

  // Repeat the loop for the case when the second beam is the photon emitter. 
  // Don't repeat for pA
  if( !( (A_2 == 1 && A_1 != 1) || (A_1 == 1 && A_2 != 1) ) ){ 
    
    for(unsigned int i = 0; i <= _nWbins - 1; ++i) {

      W = _wMin + double(i)*dW + 0.5*dW;
    
      for(unsigned int j = 0; j <= _nYbins - 1; ++j) {

        Y = -1.0*_yMax + double(j)*dY + 0.5*dY;

        beam = 1; 
        Egamma = 0.5*W*exp(-Y);        

        dndWdY = 0.; 
 
        if( Egamma > Eth && Egamma < maxPhotonEnergy() ){

	  csgA=getcsgA(Egamma,W,beam);
          dndWdY = Egamma*photonFlux(Egamma,beam)*csgA*breitWigner(W,bwnorm);

        }

        wylumfile << dndWdY << endl;

      }
    }
  }

  wylumfile << bwnorm << endl;
  wylumfile.close();
  
  if(_interferenceEnabled==1) 
    pttablegen();
 
}


//______________________________________________________________________________
void photonNucleusLuminosity::pttablegen()
{
  //  Calculates the pt spectra for VM production with interference
  //  Follows S. Klein and J. Nystrand, Phys. Rev Lett. 84, 2330 (2000).
  //  Written by S. Klein, 8/2002
  
  //  fill in table pttable in one call
  //  Integrate over all y (using the same y values as in table yarray
  //  note that the cross section goes from ymin (<0) to ymax (>0), in numy points
  //  here,  we go from 0 to ymax in (numy/2)+1 points
  //  numy must be even.
  
  //  At each y, calculate the photon energies Egamma1 and Egamma2
  //  and the two photon-A cross sections
  
  //  loop over each p_t entry in the table.
  
  //  Then, loop over b and phi (the angle between the VM \vec{p_t} and \vec{b}
  //  and calculate the cross section at each step.
  //  Put the results in pttable

  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);

  wylumfile.open(wyFileName.c_str(),ios::app);

  double param1pt[500],param2pt[500];
  double  *ptparam1=param1pt;
  double  *ptparam2=param2pt;
  double dY=0.,Yp=0.,Egamma1=0.,Egamma2=0.,Wgp=0.,cs=0.,cvma=0.,Av=0.,tmin=0.,tmax=0.,ax=0.,bx=0.;
  double csgA=0.,t=0.,sig_ga_1=0.,sig_ga_2=0.,bmax=0.,bmin=0.,db=0.,pt=0.,sum1=0.,b=0.,A1=0.,A2=0.;
  double sumg=0.,theta=0.,amp_i_2=0.,sumint=0.;
  int NGAUSS=0,NBIN=0;
  
  double xg[16]={.0483076656877383162E0,.144471961582796493E0,
		 .239287362252137075E0, .331868602282127650E0,
		 .421351276130635345E0, .506899908932229390E0,
		 .587715757240762329E0, .663044266930215201E0,
		 .732182118740289680E0, .794483795967942407E0,
		 .849367613732569970E0, .896321155766052124E0,
		 .934906075937739689E0, .964762255587506430E0,
		 .985611511545268335E0, .997263861849481564E0};
  double ag[16]={.0965400885147278006E0, .0956387200792748594E0,
		 .0938443990808045654E0, .0911738786957638847E0,
		 .0876520930044038111E0, .0833119242269467552E0,
		 .0781938957870703065E0, .0723457941088485062E0,
		 .0658222227763618468E0, .0586840934785355471E0,
		 .0509980592623761762E0, .0428358980222266807E0,
		 .0342738629130214331E0, .0253920653092620595E0,
		 .0162743947309056706E0, .00701861000947009660E0};

  NGAUSS=16;

  //Setting input calls to variables/less calls this way.
  double Ymax=_yMax;
  int numy = _nYbins;
  double Ep = _protonEnergy;
  int ibreakup = _beamBreakupMode;
  double NPT = _nmbPtBinsInterference;
  double gamma_em = _beamLorentzGamma;
  double mass = getChannelMass();
  int beam; 
  
  //  loop over y from 0 (not -ymax) to yma
  // changed this to go from -ymax to ymax to aid asymmetric collisions
  
  dY=(2.*Ymax)/numy;
  for(int jy=1;jy<=numy;jy++){
    Yp=-Ymax+((double(jy)-0.5)*dY);
    
    // Find the photon energies.  Yp >= 0, so Egamma2 is smaller (no longer true if we integrate over all Y)
    // Use the vector meson mass for W here - neglect the width
    
    Egamma1 = 0.5*mass*exp(Yp);
    Egamma2 = 0.5*mass*exp(-Yp);
    
    //  Find the sigma(gammaA) for the two directions
    //  Photonuclear Cross Section 1
    //  Gamma-proton CM energy
    beam=2; 
    
    Wgp=sqrt(2.*Egamma1*(Ep+sqrt(Ep*Ep-starlightConstants::protonMass*
				 starlightConstants::protonMass))+starlightConstants::protonMass*starlightConstants::protonMass);
    
    // Calculate V.M.+proton cross section
    
    cs=sqrt(16.*starlightConstants::pi*vmPhotonCoupling()*slopeParameter()*
	    starlightConstants::hbarc*starlightConstants::hbarc*sigmagp(Wgp)
	    /starlightConstants::alpha);
    // Calculate V.M.+nucleus cross section
    cvma=sigma_A(cs,beam);
    
    // Calculate Av = dsigma/dt(t=0) Note Units: fm**2/Gev**2
    
    Av=(starlightConstants::alpha*cvma*cvma)/(16.*starlightConstants::pi
					      *vmPhotonCoupling()*starlightConstants::hbarc*starlightConstants::hbarc);
    
    tmin  = ((mass*mass)/(4.*Egamma1*gamma_em)*(mass*mass)/(4.*Egamma1*gamma_em));
    tmax  = tmin + 0.25;
    ax    = 0.5*(tmax-tmin);
    bx    = 0.5*(tmax+tmin);
    csgA  = 0.;
    
    for(int k=0;k<NGAUSS;k++){
      t     = sqrt(ax*xg[k]+bx);
      csgA  = csgA + ag[k]*getbbs().beam2().formFactor(t)*getbbs().beam2().formFactor(t);
      t     = sqrt(ax*(-xg[k])+bx);
      csgA  = csgA + ag[k]*getbbs().beam2().formFactor(t)*getbbs().beam2().formFactor(t);
    }
    
    csgA = 0.5*(tmax-tmin)*csgA;
    csgA = Av*csgA;
    sig_ga_1 = csgA;
	   
    // Photonuclear Cross Section 2
    beam=1; 
    
    Wgp=sqrt(2.*Egamma2*(Ep+sqrt(Ep*Ep-starlightConstants::protonMass*
				 starlightConstants::protonMass))+starlightConstants::protonMass*starlightConstants::protonMass);
    
    cs=sqrt(16.*starlightConstants::pi*vmPhotonCoupling()*slopeParameter()*
	    starlightConstants::hbarc*starlightConstants::hbarc*sigmagp(Wgp)/starlightConstants::alpha);
    
    cvma=sigma_A(cs,beam);
    
    Av=(starlightConstants::alpha*cvma*cvma)/(16.*starlightConstants::pi
					      *vmPhotonCoupling()*starlightConstants::hbarc*starlightConstants::hbarc);
    
    tmin  = (((mass*mass)/(4.*Egamma2*gamma_em))*((mass*mass)/(4.*Egamma2*gamma_em)));
    tmax  = tmin + 0.25;
    ax    = 0.5*(tmax-tmin);
    bx    = 0.5*(tmax+tmin);
    csgA  = 0.;
    
    for(int k=0;k<NGAUSS;k++){
      t     = sqrt(ax*xg[k]+bx);
      csgA  = csgA + ag[k]*getbbs().beam1().formFactor(t)*getbbs().beam1().formFactor(t);
      t     = sqrt(ax*(-xg[k])+bx);
      csgA  = csgA + ag[k]*getbbs().beam1().formFactor(t)*getbbs().beam1().formFactor(t);
    }
	   
    csgA = 0.5*(tmax-tmin)*csgA;
    csgA = Av*csgA;
    sig_ga_2 = csgA;
    
    //  Set up pttables - they find the reduction in sigma(pt)
    //  due to the nuclear form factors.
    //  Use the vector meson mass for W here - neglect width in
    //  interference calculation
    
    ptparam1=vmsigmapt(mass,Egamma1,ptparam1, 2);
    ptparam2=vmsigmapt(mass,Egamma2,ptparam2, 1);
    
    bmin = getbbs().beam1().nuclearRadius()+getbbs().beam2().nuclearRadius();
    //  if we allow for nuclear breakup, use a slightly smaller bmin
    if (ibreakup != 1) 
      bmin=0.95*bmin;

    //  set  bmax according to the smaller photon energy, following flux.f
    if (Egamma1 >=Egamma2) {
    bmax=bmin+6.*starlightConstants::hbarc*gamma_em/Egamma2;
    }
    else {
    bmax=bmin+6.*starlightConstants::hbarc*gamma_em/Egamma1;
    }    
    //  set number of bins to a reasonable number to start
    NBIN = 2000;
    db   = (bmax-bmin)/float(NBIN);
    // loop over pt
    for(int i=1;i<=NPT;i++){
      
      pt = (float(i)-0.5)*_ptBinWidthInterference;
      sum1=0.0;
      // loop over b
      for(int j=1;j<=NBIN;j++){
	
	b = bmin + (float(j)-0.5)*db;
	A1 = Egamma1*getbbs().beam1().photonDensity(Egamma1,b)*sig_ga_1*ptparam1[i];
	A2 = Egamma2*getbbs().beam2().photonDensity(Egamma2,b)*sig_ga_2*ptparam2[i];
	sumg=0.0;

	//  do this as a Gaussian integral, from 0 to pi
	for(int k=0;k<NGAUSS;k++){
	  
	  theta=xg[k]*starlightConstants::pi;
	  //  allow for a linear sum of interfering and non-interfering amplitudes
	  amp_i_2 = A1 + A2 - 2.*_interferenceStrength
	    *sqrt(A1*A2)*cos(pt*b*cos(theta)/starlightConstants::hbarc);
	  sumg  = sumg+ag[k]*amp_i_2;
	}
	//  this is dn/dpt^2
	//  The factor of 2 is because the theta integral is only from 0 to pi
	sumint=2.*sumg*b*db;
	if (ibreakup > 1)
	  sumint=sumint*getbbs().probabilityOfBreakup(b);
	sum1 = sum1 + sumint;

      }
      //  normalization is done in readDiffLum.f
      //  This is d^2sigma/dpt^2; convert to dsigma/dpt


      wylumfile << sum1*pt*_ptBinWidthInterference <<endl;
      //  end of pt loop
    }
    //  end of y loop
  }
  wylumfile.close();
}


//______________________________________________________________________________
double *photonNucleusLuminosity::vmsigmapt(double W, double Egamma, double *SIGMAPT, int beam)
{
  //
  //  This subroutine calculates the effect of the nuclear form factor
  // on the pt spectrum, for use in interference calculations
  // For an interaction with mass W and photon energy Egamma,
  // it calculates the cross section suppression SIGMAPT(PT)
  // as a function of pt.
  // The input pt values come from pttable.inc

  
  double pxmax=0.,pymax=0.,dx=0.,Epom=0.,pt=0.,px0=0.,py0=0.,sum=0.,sumy=0.;
  double py=0.,px=0.,pt1=0.,pt2=0.,f1=0.,f2=0.,q1=0.,q2=0.,norm=0.;
  int NGAUSS =0,Nxbin=0;
  double xg[16]={.0483076656877383162e0,.144471961582796493e0,
		 .239287362252137075e0, .331868602282127650e0,
		 .421351276130635345e0, .506899908932229390e0,
		 .587715757240762329e0, .663044266930215201e0,
		 .732182118740289680e0, .794483795967942407e0,
		 .849367613732569970e0, .896321155766052124e0,
		 .934906075937739689e0, .964762255587506430e0,
		 .985611511545268335e0, .997263861849481564e0};
  double ag[16]={.0965400885147278006e0, .0956387200792748594e0,
		 .0938443990808045654e0, .0911738786957638847e0,
		 .0876520930044038111e0, .0833119242269467552e0,
		 .0781938957870703065e0, .0723457941088485062e0,
		 .0658222227763618468e0, .0586840934785355471e0,
		 .0509980592623761762e0, .0428358980222266807e0,
		 .0342738629130214331e0, .0253920653092620595e0,
		 .0162743947309056706e0, .00701861000947009660e0};
  NGAUSS=16;

  //     >> Initialize
  if (beam == 1) {
  pxmax = 10.*(starlightConstants::hbarc/getbbs().beam2().nuclearRadius());
  pymax = 10.*(starlightConstants::hbarc/getbbs().beam2().nuclearRadius());
  }
  else {
  pxmax = 10.*(starlightConstants::hbarc/getbbs().beam1().nuclearRadius());
  pymax = 10.*(starlightConstants::hbarc/getbbs().beam1().nuclearRadius());
  }

  Nxbin = 500;
  
  dx = 2.*pxmax/double(Nxbin);
  Epom   = W*W/(4.*Egamma);
  
  //     >> Loop over total Pt to find distribution
      
      for(int k=1;k<=_nmbPtBinsInterference;k++){
	
              pt=_ptBinWidthInterference*(double(k)-0.5);
              
              px0 = pt;
              py0 = 0.0;
              
              //  For each total Pt, integrate over Pt1, , the photon pt
              //  The pt of the Pomeron  is the difference
              //  pt1 is
              sum=0.;
              for(int i=1;i<=Nxbin;i++){
		
		px = -pxmax + (double(i)-0.5)*dx;
		sumy=0.0;
		for(int j=0;j<NGAUSS;j++){

		  py = 0.5*pymax*xg[j]+0.5*pymax;
		  //  photon pt
		  pt1 = sqrt( px*px + py*py );
		  //  pomeron pt
		  pt2 = sqrt( (px-px0)*(px-px0) + (py-py0)*(py-py0) );
		  q1  = sqrt( ((Egamma/_beamLorentzGamma)*(Egamma/_beamLorentzGamma)) + pt1*pt1 );        
		  q2  = sqrt( ((Epom/_beamLorentzGamma)*(Epom/_beamLorentzGamma))   + pt2*pt2 );
		  
		  //  photon form factor
		  // add in phase space factor?
		  if (beam ==2) {
		  f1  = (getbbs().beam1().formFactor(q1*q1)*getbbs().beam1().formFactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  
		  //  Pomeron form factor
		  f2  = getbbs().beam2().formFactor(q2*q2)*getbbs().beam2().formFactor(q2*q2);
		  }
		  else {
		  f1  = (getbbs().beam2().formFactor(q1*q1)*getbbs().beam2().formFactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  
		  //  Pomeron form factor
		  f2  = getbbs().beam1().formFactor(q2*q2)*getbbs().beam1().formFactor(q2*q2);
		  }
		  sumy= sumy + ag[j]*f1*f2;
		  
		  //  now consider other half of py phase space - why is this split?
		  py = 0.5*pymax*(-xg[j])+0.5*pymax;
		  pt1 = sqrt( px*px + py*py );
		  pt2 = sqrt( (px-px0)*(px-px0) + (py-py0)*(py-py0) );
		  q1  = sqrt( ((Egamma/_beamLorentzGamma)*Egamma/_beamLorentzGamma) + pt1*pt1 );
		  q2  = sqrt( ((Epom/_beamLorentzGamma)*(Epom/_beamLorentzGamma))   + pt2*pt2 );
		  //  add in phase space factor?
		  if (beam ==2) {
		  f1  = (getbbs().beam1().formFactor(q1*q1)*getbbs().beam1().formFactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  
		  //  Pomeron form factor
		  f2  = getbbs().beam2().formFactor(q2*q2)*getbbs().beam2().formFactor(q2*q2);
		  }
		  else {
		  f1  = (getbbs().beam2().formFactor(q1*q1)*getbbs().beam2().formFactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  
		  //  Pomeron form factor
		  f2  = getbbs().beam1().formFactor(q2*q2)*getbbs().beam1().formFactor(q2*q2);
		  }

		  sumy= sumy + ag[j]*f1*f2;
      
		}
		//         >> This is to normalize the gaussian integration
		sumy = 0.5*pymax*sumy;
		//         >> The 2 is to account for py: 0 -- pymax
		sum  = sum + 2.*sumy*dx;
	      }
	      
	      if(k==1) norm=1./sum;
	      SIGMAPT[k]=sum*norm;
      }
      return (SIGMAPT);
}

