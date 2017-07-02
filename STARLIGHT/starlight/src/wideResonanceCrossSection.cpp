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
// $Rev:: 286                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2017-05-03 22:32:03 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "wideResonanceCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
wideResonanceCrossSection::wideResonanceCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	: photonNucleusCrossSection(inputParametersInstance, bbsystem)//hrm
{
	_wideWmax = _wMax;
	_wideWmin = _wMin;
	_wideYmax = inputParametersInstance.maxRapidity();;
	_wideYmin = -1.0 * _wideYmax;
	_wideNumY = inputParametersInstance.nmbRapidityBins(); 
	_Ep       = inputParametersInstance.protonEnergy();
	_gamma1     = inputParametersInstance.beam1LorentzGamma();
	_gamma2     = inputParametersInstance.beam2LorentzGamma();
	_printDef = inputParametersInstance.printVM();	
}


//______________________________________________________________________________
wideResonanceCrossSection::~wideResonanceCrossSection()
{

}


//______________________________________________________________________________
void
wideResonanceCrossSection::crossSectionCalculation(const double bwnormsave)
{
	//     This subroutine calculates the cross-section assuming a wide
	//     (Breit-Wigner) resonance.

	double W,dW,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double csgA1,csgA2,csgA12,int_r,dR;
	double Eth;
	int    I,J,NW,NY,beam;
        double yprintcm,egaprint,Wgpprint,csgAprint; 
        std::vector<double> yVal(_wideNumY);
        std::vector<double> dsigdyVal(_wideNumY);
        std::vector<double> w1Val(_wideNumY);
        std::vector<double> Egamma1Val(_wideNumY);
        std::vector<double> ngamma1Val(_wideNumY);
        std::vector<double> sigGA1Val(_wideNumY);
        std::vector<double> dsigdy1Val(_wideNumY);
        std::vector<double> w2Val(_wideNumY);
        std::vector<double> Egamma2Val(_wideNumY);
        std::vector<double> ngamma2Val(_wideNumY);
        std::vector<double> sigGA2Val(_wideNumY);
        std::vector<double> dsigdy2Val(_wideNumY);

	// Get the Collider (lab) frame kinematics for printing of the differential cross section 
        double ybeam1 = acosh(_gamma1);
	double ybeam2 = -acosh(_gamma2);
        double ycm = (ybeam1+ybeam2)/2.;

	double bwnorm = bwnormsave; //used to transfer the bwnorm from the luminosity tables
        double mvm = getChannelMass();

	//gamma+nucleon threshold.
	Eth=0.5*(((_wideWmin+protonMass)*(_wideWmin+protonMass)
	          -protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
                   
	NW   = 100;
	dW   = (_wideWmax-_wideWmin)/double(NW);
  
	// NY   =  1200;
	NY   =  _wideNumY;
	dY   = (_wideYmax-_wideYmin)/double(NY);
  
	if (getBNORM()  ==  0.){
		cout<<" Using Breit-Wigner Resonance Profile."<<endl;
	}
	else{
		cout<<" Using Breit-Wigner plus direct pi+pi- profile."<<endl;
	}
  
	cout<<" Integrating over W from "<<_wideWmin<<" to "<<_wideWmax<<endl;

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();

	int_r=0.;
 
        // Do this first for the case when the first beam is the photon emitter 
        // Treat pA separately with defined beams 
        // The variable beam (=1,2) defines which nucleus is the target 
        if( ( A_1 == 1 && A_2 != 1 ) || (A_2 > 1 && A_1 > 1) || (A_1 == 1 && A_2 == 1) ){
	  for(J=0;J<=NY-1;J++){

	    y1  = _wideYmin + double(J)*dY;
	    y2  = _wideYmin + double(J+1)*dY;
	    y12 = 0.5*(y1+y2);
            yVal[J] = y12; 
	    yprintcm = y12 - ycm;
            beam = 2;
            double localdsig=0.0;
	    
	    for(I=0;I<=NW-1;I++){
    
		W = _wideWmin + double(I)*dW + 0.5*dW;
	        ega1  = 0.5*W*exp(y1);
	        ega2  = 0.5*W*exp(y2);
	        ega12 = 0.5*W*exp(y12);
	        egaprint = 0.5*W*exp(yprintcm);

	        if(ega1 < Eth || ega2 < Eth) continue;
		if(ega2 > maxPhotonEnergy()) continue;
          
		csgA1=getcsgA(ega1,W,beam);

		//         >> Middle Point                      =====>>>
		csgA12=getcsgA(ega12,W,beam);         

		//         >> Second Point                      =====>>>
		csgA2=getcsgA(ega2,W,beam);

 		// Print Point
		csgAprint=getcsgA(egaprint,W,beam);
     
		//>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
		dR  = ega1*photonFlux(ega1,beam)*csgA1;
		dR  = dR + 4.*ega12*photonFlux(ega12,beam)*csgA12;
		dR  = dR + ega2*photonFlux(ega2,beam)*csgA2;
		dR  = dR*(dY/6.)*breitWigner(W,bwnorm)*dW;
                //localdsig += dR;
                localdsig += egaprint*photonFlux(egaprint,beam)*csgAprint*dY*breitWigner(W,bwnorm)*dW;
		
		int_r = int_r+dR;  
	    }

	    //Take the printed Wgp and photon flux at the mass of the VM
	    egaprint = 0.5*mvm*exp(yprintcm);
            Wgpprint = sqrt(2.*egaprint*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
		               +starlightConstants::protonMass*starlightConstants::protonMass);
            dsigdyVal[J]  = 10.*localdsig/dY; //This is dsig/dy in millibarn 
            w1Val[J] = Wgpprint; 
            Egamma1Val[J] = egaprint;
	    ngamma1Val[J] = egaprint*photonFlux(egaprint,beam);
	    sigGA1Val[J]  = 10.*getcsgA(egaprint,mvm,beam); //This is sigma(gamma+A) in millibarn 
 	    dsigdy1Val[J] = 10.*localdsig/dY; 
            
	  }
        }  

        // Repeat the loop for the case when the second beam is the photon emitter. 
        // For pA, do it here or above
        if( ( A_2 == 1 && A_1 != 1 ) || (A_2 > 1 && A_1 > 1) || (A_1 == 1 && A_2 == 1) ){ 

	  for(J=0;J<=NY-1;J++){

	    y1  = _wideYmin + double(J)*dY;
	    y2  = _wideYmin + double(J+1)*dY;
	    y12 = 0.5*(y1+y2);
            yVal[J] = y12; 
	    yprintcm = y12 - ycm;
            beam = 1;
            double localdsig=0.0;

	    for(I=0;I<=NW-1;I++){
    
		W = _wideWmin + double(I)*dW + 0.5*dW;
    
		ega1  = 0.5*W*exp(-y1);
		ega2  = 0.5*W*exp(-y2);
		ega12 = 0.5*W*exp(-y12);
	        egaprint = 0.5*W*exp(-yprintcm);
      
		if(ega1< Eth || ega2 < Eth) continue;
		if(ega1 > maxPhotonEnergy()) continue;
          
		csgA1=getcsgA(ega1,W,beam);

		//         >> Middle Point                      =====>>>
		csgA12=getcsgA(ega12,W,beam);         

		//         >> Second Point                      =====>>>
		csgA2=getcsgA(ega2,W,beam);

 		// Print Point
		csgAprint=getcsgA(egaprint,W,beam);
      
		//>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
		dR  = ega1*photonFlux(ega1,beam)*csgA1;
		dR  = dR + 4.*ega12*photonFlux(ega12,beam)*csgA12;
		dR  = dR + ega2*photonFlux(ega2,beam)*csgA2;
		dR  = dR*(dY/6.)*breitWigner(W,bwnorm)*dW;
                //localdsig += dR;
                localdsig += egaprint*photonFlux(egaprint,beam)*csgAprint*dY*breitWigner(W,bwnorm)*dW;
      
		int_r = int_r+dR;  
	     }

	    //Take the printed Wgp and photon flux at the mass of the rho0
	    egaprint = 0.5*mvm*exp(-yprintcm);
            Wgpprint = sqrt(2.*egaprint*(_Ep+sqrt(_Ep*_Ep-starlightConstants::protonMass*starlightConstants::protonMass))
		               +starlightConstants::protonMass*starlightConstants::protonMass);
            dsigdyVal[J]  += 10.*localdsig/dY; //This is dsig/dy in millibarn 
            w2Val[J] = Wgpprint; 
            Egamma2Val[J] = egaprint;
	    ngamma2Val[J] = egaprint*photonFlux(egaprint,beam);
	    sigGA2Val[J]  = 10.*getcsgA(egaprint,mvm,beam); //This is sigma(gamma+A) in millibarn 
 	    dsigdy2Val[J] = 10.*localdsig/dY; 

	  }
        }

	cout<<endl;
	if( _impulseSelected == 1 )cout<<" Using impulse approximation. Nuclear effects removed."<<endl; 
	if (0.01*int_r > 1.){
	  cout<< " Total cross section: "<<0.01*int_r<<" barn."<<endl;
	} else if (10.*int_r > 1.){
	  cout<< " Total cross section: " <<10.*int_r<<" mb."<<endl;
        } else if (10000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000.*int_r<<" microb."<<endl;
        } else if (10000000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000000.*int_r<<" nanob."<<endl;
        } else if (1.E10*int_r > 1.){
	  cout<< " Total cross section: "<<1.E10*int_r<<" picob."<<endl;
        } else {
	  cout<< " Total cross section: " <<1.E13*int_r<<" femtob."<<endl;
        }
	cout<<endl;
	setPhotonNucleusSigma(0.01*int_r);

	//Print detailed cross section (dsig/dy), if requested  
        if( _printDef == 1 ||  _printDef == 2 ){

	  if(  _printDef == 2 ){
            printf("Printing detailed information from the vector meson cross section calculation.  \n \n");
	    printf("First column gives the rapidity in the lab frame. Second and third column give  \n");
	    printf("the gamma-nucleon center of mass energy and photon flux (k*(dn/dk))from the     \n");
	    printf("first beam (defined by BEAM_1_Z, BEAM_1_A etc.). The fourth column gives the    \n"); 
            printf("gamma+A cross section with the second beam as target. The  fifth column gives   \n");
	    printf("the contribution to the cross section (dsig/dy) for this combination. The 6th - \n");
	    printf("9th columns give the corresponding information for the opposite combination     \n");
	    printf("(the beam defined by BEAM_2_Z etc. emits the photon). The last (10th) column    \n");
	    printf("gives the total dsig/dy. \n");
	  }

	  double maxVal = 0.0; for(I=0;I<=(NY-1);I++){if(dsigdyVal[I] > maxVal)maxVal=dsigdyVal[I];}
          double scaleFactor = 0.0; 
          if( maxVal > 1.0 ){
	    scaleFactor = 1.0;//Default is millibarn 
            if( _printDef == 1){
              printf("Rapidity          dsig/dy (millibarn) \n");
	    } else if ( _printDef == 2 ){
              printf("Cross sections are in millibarn. \n \n");
	    }
	  } else if ( maxVal > 0.001 ){
	    scaleFactor = 1000.0; //This is for microbarn
            if( _printDef == 1){
              printf("Rapidity          dsig/dy (microbarn) \n");
	    } else if ( _printDef == 2 ){
              printf("Cross sections are in microbarn. \n \n");
	    }
	  } else if ( maxVal > 0.000001 ){
	    scaleFactor = 1000000.0; //This is for nanobarn
            if( _printDef == 1){
              printf("Rapidity          dsig/dy (nanobarn) \n");
            } else if ( _printDef == 2 ){
              printf("Cross sections are in nanobarn. \n \n");
	    }
	  } else if ( maxVal > 1.E-9 ){
	    scaleFactor = 1.E9;//This is for picobarn
            if( _printDef == 1){
              printf("Rapidity          dsig/dy (picobarn) \n");
            } else if ( _printDef == 2 ){
              printf("Cross sections are in picobarn. \n \n");
	    }
	  } else {
	    scaleFactor = 1.E12;//This is for femtobarn
            if( _printDef == 1){
              printf("Rapidity          dsig/dy (femtobarn) \n");
            } else if ( _printDef == 2 ){
              printf("Cross sections are in femtobarn. \n \n");
	    }
	  }

          if( _printDef == 2 ){
             printf(" y       Wgp(1) (GeV) k*(dn_1/dk)  sigma_2(gam+A)  dsig_1/dy  Wgp(2) (GeV) k*(dn_2/dk) sigma_1(gam+A)  dsig_2/dy    dsig/dy \n");
	  }
	  
	  for(J=0;J<=(NY-1);J++){
            if( _printDef == 1){ 
              printf("%+6.2f          %10.4f \n",yVal[J],scaleFactor*dsigdyVal[J]);
	    } else if ( _printDef == 2 ){
              printf("%+6.2f   %.4E   %.4E   %.4E   %9.4f     ",yVal[J],w1Val[J],ngamma1Val[J],scaleFactor*sigGA1Val[J],scaleFactor*dsigdy1Val[J]); 
              printf("%.4E   %.4E  %.4E   %9.4f    %9.4f \n",w2Val[J],ngamma2Val[J],scaleFactor*sigGA2Val[J],scaleFactor*dsigdy2Val[J],scaleFactor*dsigdyVal[J]); 
	    }
	  }
          printf("\n");
	}

}
