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
// $Rev:: 274                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2016-09-12 00:40:25 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


#include "starlightconstants.h"
#include "gammagammasingle.h"
#include "starlightconfig.h"

using namespace std;


//______________________________________________________________________________
Gammagammasingle::Gammagammasingle(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
: eventChannel(inputParametersInstance, bbsystem)
#ifdef ENABLE_PYTHIA
,_pyDecayer()
#endif
{

#ifdef ENABLE_PYTHIA
    _pyDecayer.init();
#endif
  
  //Storing inputparameters into protected members for use
  _GGsingInputnumw=inputParametersInstance.nmbWBins();
  _GGsingInputnumy=inputParametersInstance.nmbRapidityBins();
  _GGsingInputpidtest=inputParametersInstance.prodParticleType();
  _GGsingInputGamma_em=inputParametersInstance.beamLorentzGamma();
  _axionMass=inputParametersInstance.axionMass(); // AXION HACK
  cout<<"SINGLE MESON pid test: "<<_GGsingInputpidtest<<endl;
  //reading in luminosity tables
  read();
  //Now calculating crosssection
  singleCrossSection();
}


//______________________________________________________________________________
Gammagammasingle::~Gammagammasingle()
{ }


//______________________________________________________________________________
void Gammagammasingle::singleCrossSection()
{
  //This function carries out a delta function cross-section calculation. For reference, see STAR Note 243, Eq. 8
  //Multiply all _Farray[] by _f_max
  double _sigmaSum=0.,remainw=0.;//_remainwd=0.;
  int ivalw =0;//_ivalwd;
  //calculate the differential cross section and place in the sigma table
  cout << "MASS  " << getMass() << "\n"; // AXION HACK, optional
  cout << "WIDTH  " << getWidth() << "\n";// AXION HACK, optional
  _wdelta=getMass();
  for(int i=0;i<_GGsingInputnumw;i++){
    for(int j=0;j<_GGsingInputnumy;j++){
      // Eq. 1 of starnote 347
      _sigmax[i][j]=(getSpin()*2.+1.)*4*starlightConstants::pi*starlightConstants::pi*getWidth()/
	(getMass()*getMass()*getMass())*_f_max*_Farray[i][j]*starlightConstants::hbarc*starlightConstants::hbarc/100.;
    }
  }
  //find the index, i,for the value of w just less than the mass because we want to use the value from the sigma table that has w=mass

  for(int i=0;i<_GGsingInputnumw;i++){
    if(getMass()>_Warray[i]) ivalw=i;
  }

  remainw = (getMass()-_Warray[ivalw])/(_Warray[ivalw+1]-_Warray[ivalw]);
  _ivalwd = ivalw;
  _remainwd = remainw;
  //if we are interested rho pairs at threshold, the just set sigma to 100nb
  switch(_GGsingInputpidtest){
  case starlightConstants::ZOVERZ03:
    _sigmaSum =0.;
    for(int j=0;j<_GGsingInputnumy-1;j++){
                        _sigmaSum = _sigmaSum +(_Yarray[j+1]-_Yarray[j])*
			  100.0E-9*(.1/getMass())*((1.-remainw)*_f_max*
						   (_Farray[ivalw][j]+_Farray[ivalw][j])/2.+remainw*_f_max*
						   (_Farray[ivalw+1][j]+_Farray[ivalw+1][j+1])/2.);
    }
    break;
  default:
    //Sum to find the total cross-section
    _sigmaSum =0.;
    for(int j =0;j<_GGsingInputnumy-1;j++){
                        _sigmaSum = _sigmaSum+
			  (_Yarray[j+1]-_Yarray[j])*((1.-remainw)*
						   (_sigmax[ivalw][j]+_sigmax[ivalw][j+1])/2.+remainw*
						   (_sigmax[ivalw+1][j]+_sigmax[ivalw+1][j+1])/2.);
    }
  }
  // if(_sigmaSum > 0.1) cout <<"The total cross-section is: "<<_sigmaSum<<" barns."<<endl;
  // else if(_sigmaSum > 0.0001)cout <<"The total cross-section is: "<<_sigmaSum*1000<<" mb."<<endl;
  // else cout <<"The total cross-section is: "<<_sigmaSum*1000000<<" ub."<<endl;
  cout<<endl;
  if (_sigmaSum > 1.){
     cout << "Total cross section: "<<_sigmaSum<<" barn."<<endl;  
  } else if (1000.*_sigmaSum > 1.){
     cout << "Total cross section: "<<1000.*_sigmaSum<<" mb."<<endl;  
  } else if (1000000.*_sigmaSum > 1.){
    cout << "Total cross section: "<<1000000.*_sigmaSum<<" microbarn."<<endl;  
  } else if (1.E9*_sigmaSum > 1.){
    cout << "Total cross section: "<<1.E9*_sigmaSum<<" nanobarn."<<endl;  
  } else if (1.E12*_sigmaSum > 1.){
    cout << "Total cross section: "<<1.E12*_sigmaSum<<" picobarn."<<endl;  
  } else {
    cout << "Total cross section: "<<1.E15*_sigmaSum<<" femtobarn."<<endl;  
  }
  cout<<endl; 
  setTotalChannelCrossSection(_sigmaSum);
     
  return;
}


//______________________________________________________________________________
void Gammagammasingle::pickw(double &w)
{
  //This function picks a w for the 2-photon calculation. 
  double sgf=0.,signorm=0.,x=0.,remainarea=0.,remainw=0.,a=0.,b=0.,c=0.;
  int ivalw=0;
  
  double * _sigofw;
  double * sgfint;
  _sigofw = new double[starlightLimits::MAXWBINS];
  sgfint = new double[starlightLimits::MAXYBINS];
 
  if(_wdelta != 0){
    w=_wdelta;
    ivalw=_ivalwd;
    remainw=_remainwd;
  }
  else{
    //deal with the case where sigma is an array
    //_sigofw is simga integrated over y using a linear interpolation
    //sigint is the integral of sgfint, normalized
    
    //integrate sigma down to a function of just w
    for(int i=0;i<_GGsingInputnumw;i++){
      _sigofw[i]=0.;
      for(int j=0;j<_GGsingInputnumy-1;j++){
	_sigofw[i] = _sigofw[i]+(_Yarray[j+1]-_Yarray[j])*(_sigmax[i][j+1]+_sigmax[i][j])/2.;
      }
    }
    //calculate the unnormalized sgfint array
    sgfint[0]=0.;
    for(int i=0;i<_GGsingInputnumw-1;i++){
      sgf=(_sigofw[i+1]+_sigofw[i])*(_Warray[i+1]-_Warray[i])/2.;
      sgfint[i+1]=sgfint[i]+sgf;
    }
    //normalize sgfint array
    signorm=sgfint[_GGsingInputnumw-1];
    
    for(int i=0;i<_GGsingInputnumw;i++){
      sgfint[i]=sgfint[i]/signorm;
    }
    //pick a random number
    x = _randy.Rndom();
    //compare x and sgfint to find the ivalue which is just less than the random number x
    for(int i=0;i<_GGsingInputnumw;i++){
      if(x > sgfint[i]) ivalw=i;
    }
    //remainder above ivalw
    remainarea = x - sgfint[ivalw];
    
    //figure out what point corresponds to the excess area in remainarea
    c = -remainarea*signorm/(_Warray[ivalw+1]-_Warray[ivalw]);
    b = _sigofw[ivalw];
    a = (_sigofw[ivalw+1]-_sigofw[ivalw])/2.;
    if(a==0.){
      remainw = -c/b;
    }
    else{
      remainw = (-b+sqrt(b*b-4.*a*c))/(2.*a);
    }
    _ivalwd = ivalw;
    _remainwd = remainw;
    //calculate the w value
    w = _Warray[ivalw]+(_Warray[ivalw+1]-_Warray[ivalw])*remainw;
    }

  delete[] _sigofw;
  delete[] sgfint;
}


//______________________________________________________________________________
void Gammagammasingle::picky(double &y)
{
  double * sigofy;
  double * sgfint;
  sigofy = new double[starlightLimits::MAXYBINS];
  sgfint = new double[starlightLimits::MAXYBINS];
  
  double remainw =0.,remainarea=0.,remainy=0.,a=0.,b=0.,c=0.,sgf=0.,signorm=0.,x=0.;
  int ivalw=0,ivaly=0;
  
  ivalw=_ivalwd;
  remainw=_remainwd;
  //average over two colums to get y array
  for(int j=0;j<_GGsingInputnumy;j++){
    sigofy[j]=_sigmax[ivalw][j]+(_sigmax[ivalw+1][j]-_sigmax[ivalw][j])*remainw;
  }
  //calculate the unnormalized sgfint
  
  sgfint[0]=0.;
  for(int j=0;j<_GGsingInputnumy-1;j++){
    sgf = (sigofy[j+1]+sigofy[j])/2.;
    sgfint[j+1]=sgfint[j]+sgf*(_Yarray[j+1]-_Yarray[j]);
  }
  
  //normalize the sgfint array
  signorm = sgfint[_GGsingInputnumy-1];
  
  for(int j=0;j<_GGsingInputnumy;j++){
    sgfint[j]=sgfint[j]/signorm;
  }
  //pick a random number
  x = _randy.Rndom();
  //compare x and sgfint to find the ivalue which is just less then the random number x
  for(int i=0;i<_GGsingInputnumy;i++){
    if(x > sgfint[i]) 
      ivaly = i;
  }
  //remainder above ivaly
  remainarea = x - sgfint[ivaly];
  //figure what point corresponds to the leftover area in remainarea
  c = -remainarea*signorm/(_Yarray[ivaly+1]-_Yarray[ivaly]);
  b = sigofy[ivaly];
  a = (sigofy[ivaly+1]-sigofy[ivaly])/2.;
  if(a==0.){
    remainy = -c/b;
  }
  else{
    remainy = (-b + sqrt(b*b-4.*a*c))/(2.*a);
  }
  //calculate the y value
  y = _Yarray[ivaly]+(_Yarray[ivaly+1]-_Yarray[ivaly])*remainy;
  delete[] sigofy;
  delete[] sgfint;
}


//______________________________________________________________________________
void Gammagammasingle::parentMomentum(double w,double y,double &E,double &px,double &py,double &pz)
{
  //this function calculates px,py,pz,and E given w and y
  double anglepp1=0.,anglepp2=0.,ppp1=0.,ppp2=0.,E1=0.,E2=0.,signpx=0.,pt=0.;
  
  //E1 and E2 are for the 2 photons in the CM frame
  E1 = w*exp(y)/2.;
  E2 = w*exp(-y)/2.;
  //calculate px and py
  //to get x and y components-- phi is random between 0 and 2*pi
  anglepp1 = _randy.Rndom();
  anglepp2 = _randy.Rndom();
  
  ppp1 = pp1(E1);
  ppp2 = pp2(E2);
  px = ppp1*cos(2.*starlightConstants::pi*anglepp1)+ppp2*cos(2.*starlightConstants::pi*anglepp2);
  py = ppp1*sin(2.*starlightConstants::pi*anglepp1)+ppp2*sin(2.*starlightConstants::pi*anglepp2);
  //Compute vector sum Pt=Pt1+Pt2 to find pt for the produced particle
  pt = sqrt(px*px+py*py);
  //W is the mass of the produced particle (not necessarily on-mass-shell).Now compute its energy and pz
  E = sqrt(w*w+pt*pt)*cosh(y);
  pz= sqrt(w*w+pt*pt)*sinh(y);
  signpx = _randy.Rndom();
  //pick the z direction
  if(signpx > 0.5) 
    pz = -pz;	
}


//______________________________________________________________________________
double Gammagammasingle::pp1(double E)
{
  // First 'copy' of pp, for nucleus 1 form factor.  The split was needed to handle asymmetric beams.  SRK 4/2015
  //  will probably have to pass in beambeamsys? that way we can do beam1.formFactor(t) or beam2..., careful with the way sergey did it for asymmetry
  //  returns on random draw from pp(E) distribution
      
  double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
  double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
  int satisfy =0;
        
  ereds = (E/_GGsingInputGamma_em)*(E/_GGsingInputGamma_em);
  Cm = sqrt(3.)*E/_GGsingInputGamma_em;
  //the amplitude of the p_t spectrum at the maximum
  singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
  //Doing this once and then storing it as a double, for the beam 1 form factor.
  Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
  //pick a test value pp, and find the amplitude there
  x = _randy.Rndom();
  pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); //Will use nucleus #1, there should be two for symmetry//nextline
  singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
  test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

  while(satisfy==0){
    u = _randy.Rndom();
    if(u*Coef <= test){
      satisfy =1;
    }
    else{
      x =_randy.Rndom();
      pp = 5*starlightConstants::hbarc/_bbs.beam1().nuclearRadius()*x;
      singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);//Symmetry
      test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
    }
  }
  return pp;
}

//______________________________________________________________________________
double Gammagammasingle::pp2(double E)
{
  // Second 'copy' of pp, for nucleus 1 form factor.  The split was needed to handle asymmetric beams.  SRK 4/2015
  //  will probably have to pass in beambeamsys? that way we can do beam1.formFactor(t) or beam2..., careful with the way sergey did it for asymmetry
  //  returns on random draw from pp(E) distribution
      
  double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
  double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
  int satisfy =0;
        
  ereds = (E/_GGsingInputGamma_em)*(E/_GGsingInputGamma_em);
  Cm = sqrt(3.)*E/_GGsingInputGamma_em;
  //the amplitude of the p_t spectrum at the maximum
  singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
  //Doing this once and then storing it as a double, which we square later...SYMMETRY?using beam1 for now.
  Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
  //pick a test value pp, and find the amplitude there
  x = _randy.Rndom();
  pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); //Will use nucleus #1, there should be two for symmetry//nextline
  singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
  test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

  while(satisfy==0){
    u = _randy.Rndom();
    if(u*Coef <= test){
      satisfy =1;
    }
    else{
      x =_randy.Rndom();
      pp = 5*starlightConstants::hbarc/_bbs.beam2().nuclearRadius()*x;
      singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);//Symmetry
      test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
    }
  }
  return pp;
}


//______________________________________________________________________________
void Gammagammasingle::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double &px1,double &py1,double &pz1,double &px2,double &py2,double &pz2,int &iFbadevent)
{
  //     This routine decays a particle into two particles of mass mdec,
  //     taking spin into account
  
  double mdec=0.,E1=0.,E2=0.;
  double pmag,ytest=0.;
  double phi,theta,xtest,dndtheta,Ecm;
  double  betax,betay,betaz;
  
  //    set the mass of the daughter particles
  switch(_GGsingInputpidtest){ 
  case starlightConstants::ZOVERZ03:
  case starlightConstants::F2:	
    mdec = starlightConstants::pionChargedMass;
    break;
  case starlightConstants::AXION:       // AXION HACK
    mdec = 0;//axion decays to two photons, set mass of decay products to zero
    break;
  case starlightConstants::F2PRIME:
    //  decays 50% to K+/K-, 50% to K_0's
    ytest = _randy.Rndom();
    if(ytest >= 0.5){
      mdec = starlightConstants::kaonChargedMass;
    }
    else{
      mdec = 0.493677;
    }
    break;
  default :
    cout<<"No default mass selected for single photon-photon particle, expect errant results"<<endl;
  }
  
  //Calculating the momentum's magnitude
    if(W < 2*mdec){
      cout<<" ERROR: W="<<W<<endl;
      iFbadevent = 1;
      return;
    }
    pmag = sqrt(W*W/4. - mdec*mdec);
//  }
  //     pick an orientation, based on the spin
  //      phi has a flat distribution in 2*pi
  phi = _randy.Rndom()*2.*starlightConstants::pi;
  
  //     find theta, the angle between one of the outgoing particles and
  //    the beamline, in the frame of the two photons
  //this will depend on spin, F2,F2' and z/z03 all have spin 2, all other photonphoton-single mesons are handled by jetset/pythia
  //Applies to spin2 mesons.
 L300td:
  theta = starlightConstants::pi*_randy.Rndom();
  xtest = _randy.Rndom();
  dndtheta = sin(theta)*sin(theta)*sin(theta)*sin(theta)*sin(theta);
  if(xtest > dndtheta)
    goto L300td;

  //     compute unboosted momenta
  px1 = sin(theta)*cos(phi)*pmag;
  py1 = sin(theta)*sin(phi)*pmag;
  pz1 = cos(theta)*pmag;
  px2 = -px1;
  py2 = -py1;
  pz2 = -pz1;
  //        compute energies
  //Changed mass to W 11/9/2000 SRK
  Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
  E1 = sqrt(mdec*mdec+px1*px1+py1*py1+pz1*pz1);
  E2 = sqrt(mdec*mdec+px2*px2+py2*py2+pz2*pz2);

  //     Lorentz transform into the lab frame
  // betax,betay,betaz are the boost of the complete system
  betax = -(px0/Ecm);
  betay = -(py0/Ecm);
  betaz = -(pz0/Ecm);
  
  transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
  transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);
  
  
  if(iFbadevent == 1)
    return;
  
  //       change particle id from that of parent to that of daughters

  switch(_GGsingInputpidtest){
    //These decay into a pi+ pi- pair
  case starlightConstants::ZOVERZ03:
  case starlightConstants::F2:
    ipid=starlightConstants::PION;
    break;
  case starlightConstants::AXION:// AXION HACK
    ipid=starlightConstants::PHOTON;  // AXION HACK
    break;      // AXION HACK
  case starlightConstants::F2PRIME:
    if( ytest >= 0.5 )
      {
	//Decays 50/50 into k+ k- or k_s k_l
	ipid=starlightConstants::KAONCHARGE;	
      }
    else
      {
	ipid=starlightConstants::KAONNEUTRAL;
      }	
    break;
  default:
    cout<<"Rethink the daughter particles"<<endl;
  }
}


//______________________________________________________________________________
starlightConstants::event Gammagammasingle::produceEvent(int &/*ievent*/)
{
  // Not in use anymore, default event struct returned
  return starlightConstants::event();
}


//______________________________________________________________________________
upcEvent Gammagammasingle::produceEvent()
{

  //    returns the vector with the decay particles inside.
  starlightConstants::event single;
  double comenergy = 0.;
  double rapidity = 0.;
  double parentE = 0.;
  double parentmomx=0.,parentmomy=0.,parentmomz=0.;

  //this function decays particles and writes events to a file
  //zeroing out the event structure
  single._numberOfTracks=0;
  for(int i=0;i<4;i++){
    single.px[i]=0.;
    single.py[i]=0.;
    single.pz[i]=0.;
    single._fsParticle[i]=starlightConstants::UNKNOWN;
    single._charge[i]=0;
  }
  
  pickw(comenergy);
  picky(rapidity);
  parentMomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz);
  
  
  if(_GGsingInputpidtest != starlightConstants::F2 && _GGsingInputpidtest != starlightConstants::F2PRIME && _GGsingInputpidtest != starlightConstants::AXION)
  {
#ifdef ENABLE_PYTHIA
    starlightParticle particle(parentmomx,parentmomy,parentmomz, parentE, getMass(),_GGsingInputpidtest , 0);
  
    _pyDecayer.addParticle(particle);
  
    return _pyDecayer.execute();
#endif
  }


  int ievent = 0;
  int iFbadevent=0;
  starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
  double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
  double px3=0.,px4=0.,py3=0.,py4=0.,pz3=0.,pz4=0.;
  double xtest=0.,ztest=0.;
  switch(_GGsingInputpidtest){
  case starlightConstants::ZOVERZ03:
    //Decays into two pairs.
    parentmomx=parentmomx/2.;
    parentmomy=parentmomy/2.;
    parentmomz=parentmomz/2.;
    //Pair #1	
    twoBodyDecay(ipid,comenergy/2.,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
    //Pair #2
    twoBodyDecay(ipid,comenergy/2.,parentmomx,parentmomy,parentmomz,px3,py3,pz3,px4,py4,pz4,iFbadevent);
    //Now add them to vectors to be written out later.
		
    single._numberOfTracks=4;//number of tracks per event
    if (iFbadevent==0){
      xtest = _randy.Rndom();
      ztest = _randy.Rndom();
      //Assigning charges randomly.
      if (xtest<0.5){
	single._charge[0]=1;//q1=1;
	single._charge[1]=-1;//q2=-1;
      }
      else{
	single._charge[0]=-1;//q1=-1;
	single._charge[1]=1;//q2=1;
      }
      if (ztest<0.5){
	single._charge[2]=1;//q3=1;
	single._charge[3]=-1;//q4=-1;
      }
      else{
	single._charge[2]=-1;//q3=-1;
	single._charge[3]=1;//q4=1;
      }
      //Track #1
      single.px[0]=px1;
      single.py[0]=py1;
      single.pz[0]=pz1;
      single._fsParticle[0]=ipid*single._charge[0];
      //Track #2                                                                                                                      
      single.px[1]=px2;
      single.py[1]=py2;
      single.pz[1]=pz2;
      single._fsParticle[1]=ipid*single._charge[1];
      //Track #3
      single.px[2]=px3;
      single.py[2]=py3;
      single.pz[2]=pz3;
      single._fsParticle[2]=ipid*single._charge[2];
      //Track #4
      single.px[3]=px4;
      single.py[3]=py4;
      single.pz[3]=pz4;
      single._fsParticle[3]=ipid*single._charge[3];
      
      ievent=ievent+1;
    }	
    
    break;
  case starlightConstants::F2:
  case starlightConstants::F2PRIME:
    twoBodyDecay(ipid,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
    
    single._numberOfTracks=2;
    if (iFbadevent==0){
      xtest = _randy.Rndom();
      if (xtest<0.5){
	single._charge[0]=1;//q1=1;
	single._charge[1]=-1;//q2=-1;
      }
      else{
	single._charge[0]=-1;//q1=-1;
	single._charge[1]=1;//q2=1;
      }	
      //Track #1
      single.px[0]=px1;
      single.py[0]=py1;
      single.pz[0]=pz1;
      single._fsParticle[0]=ipid*single._charge[0]; 
      //Track #2
      single.px[1]=px2;
      single.py[1]=py2;
      single.pz[1]=pz2;
      single._fsParticle[1]=ipid*single._charge[1];
      ievent=ievent+1;
    }
    break;

  case starlightConstants::AXION:      // AXION HACK, start
    twoBodyDecay(ipid,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);

    single._numberOfTracks=2;
    if (iFbadevent==0){

        single._charge[0]=0;//q1=0;
        single._charge[1]=0;//q2=0;

      //Track #1
      single.px[0]=px1;
      single.py[0]=py1;
      single.pz[0]=pz1;
      single._fsParticle[0]=ipid;
      //Track #2
      single.px[1]=px2;
      single.py[1]=py2;
      single.pz[1]=pz2;
      single._fsParticle[1]=ipid;
      ievent=ievent+1;

    }
    break;  // AXION HACK, end


  default:
    break;
  }
  
  return upcEvent(single);
}


//______________________________________________________________________________
double Gammagammasingle::getMass()
{
  using namespace starlightConstants;
  double singlemass=0.;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlemass = starlightConstants::etaMass;
    break;
  case starlightConstants::ETAPRIME:
    singlemass = starlightConstants::etaPrimeMass;
    break;
  case starlightConstants::ETAC:
    singlemass = starlightConstants::etaCMass;
    break;
  case starlightConstants::F0:
    singlemass = starlightConstants::f0Mass;
    break;
  case starlightConstants::F2:
    singlemass = starlightConstants::f2Mass;
    break;
  case starlightConstants::A2:
    singlemass = starlightConstants::a2Mass;
    break;
  case starlightConstants::F2PRIME:
    singlemass = starlightConstants::f2PrimeMass;
    break;
  case starlightConstants::ZOVERZ03:
    singlemass = starlightConstants::zoverz03Mass;
    break;
  case starlightConstants::AXION: // AXION HACK
    singlemass = _axionMass;      // AXION HACK
    break; // AXION HACK
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getmass(), mass = 0."<<endl;
  }
  return singlemass;
}


//______________________________________________________________________________
double Gammagammasingle::getWidth()
{

  /* Partial widths(GAMMA(gammgamma)) taken from PDG 2014- Chinese Physics C 38, no 9, Sept. 2014.*/
  double singlewidth=0.;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlewidth = starlightConstants::etaPartialggWidth;
    break;
  case starlightConstants::ETAPRIME:
    singlewidth = starlightConstants::etaPrimePartialggWidth;
    break;
  case starlightConstants::ETAC:
    singlewidth = starlightConstants::etaCPartialggWidth;
    break;
  case starlightConstants::F0:
    singlewidth = starlightConstants::f0PartialggWidth;
    break;
  case starlightConstants::F2:
    singlewidth = starlightConstants::f2PartialggWidth;
    break;
  case starlightConstants::A2:
    singlewidth = starlightConstants::a2PartialggWidth;
    break;
  case starlightConstants::F2PRIME:
    singlewidth = starlightConstants::f2PrimePartialggWidth;
    break;
  case starlightConstants::ZOVERZ03:
    singlewidth = starlightConstants::zoverz03PartialggWidth;
    break;
  case starlightConstants::AXION: // AXION HACK
    singlewidth = 1/(64*starlightConstants::pi)*_axionMass*_axionMass*_axionMass/(1000*1000);//Fix Lambda=1000 GeV,rescaling is trivial.    // AXION HACK
    break;
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getwidth(), width = 0."<<endl;
  }
  return singlewidth; 
}


//______________________________________________________________________________
double Gammagammasingle::getSpin()
{
  double singlespin=0.5;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlespin = starlightConstants::etaSpin;
    break;
  case starlightConstants::ETAPRIME:
    singlespin = starlightConstants::etaPrimeSpin;
    break;
  case starlightConstants::ETAC:
    singlespin = starlightConstants::etaCSpin;
    break;
  case starlightConstants::F0:
    singlespin = starlightConstants::f0Spin;
    break;
  case starlightConstants::F2:
    singlespin = starlightConstants::f2Spin;
    break;
  case starlightConstants::A2:
    singlespin = starlightConstants::a2Spin;
    break;
  case starlightConstants::F2PRIME:
    singlespin = starlightConstants::f2PrimeSpin;
    break;
  case starlightConstants::ZOVERZ03:
    singlespin = starlightConstants::zoverz03Spin;
    break;
  case starlightConstants::AXION:// AXION HACK
    singlespin = starlightConstants::axionSpin;// AXION HACK
    break;// AXION HACK
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getspin(), spin = 0."<<endl;
  }
  return singlespin;
}



