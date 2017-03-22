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
// $Rev:: 280                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-12-12 20:38:20 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "narrowResonanceCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
narrowResonanceCrossSection::narrowResonanceCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	:photonNucleusCrossSection(inputParametersInstance, bbsystem)
{
	_narrowYmax = inputParametersInstance.maxRapidity();
	_narrowYmin = -1.0*_narrowYmax;
	_narrowNumY = inputParametersInstance.nmbRapidityBins();
	_Ep         = inputParametersInstance.protonEnergy();	
	_printDef   = inputParametersInstance.printVM();	
}


//______________________________________________________________________________
narrowResonanceCrossSection::~narrowResonanceCrossSection()
{ }


//______________________________________________________________________________
void
narrowResonanceCrossSection::crossSectionCalculation(const double)  // _bwnormsave (unused)
{
	// This subroutine calculates the vector meson cross section assuming
	// a narrow resonance.  For reference, see STAR Note 386.
  
	double W,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double csgA1,csgA2,csgA12,int_r,dR;
	double Eth;
	int    I,J,NY,beam;
        std::vector<double> yVal(_narrowNumY);
        std::vector<double> dsigdyVal(_narrowNumY);
        std::vector<double> Egamma1Val(_narrowNumY);
        std::vector<double> ngamma1Val(_narrowNumY);
        std::vector<double> sigGA1Val(_narrowNumY);
        std::vector<double> dsigdy1Val(_narrowNumY);
        std::vector<double> Egamma2Val(_narrowNumY);
        std::vector<double> ngamma2Val(_narrowNumY);
        std::vector<double> sigGA2Val(_narrowNumY);
        std::vector<double> dsigdy2Val(_narrowNumY);

	NY   =  _narrowNumY;
	dY   = (_narrowYmax-_narrowYmin)/double(NY);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	Eth=0.5*(((W+protonMass)*(W+protonMass)-
	          protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
  
	// cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
        printf(" gamma+nucleon threshold: %e GeV \n", Eth);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();
  
 	int_r=0.;

        // Do this first for the case when the first beam is the photon emitter 
        // Treat pA separately with defined beams 
        // The variable beam (=1,2) defines which nucleus is the target 
	for(J=0;J<=(NY-1);J++){
    
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
                yVal[J] = y12; 
		
                if( A_2 == 1 && A_1 != 1 ){
                  // pA, first beam is the nucleus and is in this case the target  
                  ega1  = 0.5*W*exp(-y1);
                  ega2  = 0.5*W*exp(-y2);
                  ega12 = 0.5*W*exp(-y12);
                  beam = 1; 
                } else if( A_1 ==1 && A_2 != 1){
                  // pA, second beam is the nucleus and is in this case the target 
                  ega1  = 0.5*W*exp(y1);
                  ega2  = 0.5*W*exp(y2);
                  ega12 = 0.5*W*exp(y12);
                  beam = 2; 
                } else {
                  ega1  = 0.5*W*exp(y1);
                  ega2  = 0.5*W*exp(y2);
                  ega12 = 0.5*W*exp(y12);
                  beam = 2; 
                }
    
		if(ega1 < Eth || ega2 < Eth)   
			continue;
		if(ega2 > maxPhotonEnergy() || ega1 > maxPhotonEnergy() ) 
			continue;

		csgA1=getcsgA(ega1,W,beam);
    
		// Middle Point                      =====>>>
		csgA12=getcsgA(ega12,W,beam); 

		// Second Point                      =====>>>
		csgA2=getcsgA(ega2,W,beam);
    
		// Sum the contribution for this W,Y.
		dR  = ega1*photonFlux(ega1,beam)*csgA1;
		dR  = dR + 4.*ega12*photonFlux(ega12,beam)*csgA12;
		dR  = dR + ega2*photonFlux(ega2,beam)*csgA2;
		dR  = dR*(dY/6.);
                
                dsigdyVal[J] = 10.*dR/dY; //This is dsig/dy in millibarn
                Egamma1Val[J] = ega12;
		ngamma1Val[J] = ega12*photonFlux(ega12,beam);
		sigGA1Val[J] = 10.*csgA12; //This is sigma(gamma+A) in millibarn 
		dsigdy1Val[J] = 10.*dR/dY; 
		
		// cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12,beam)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;

		int_r = int_r+dR;
	}

        // Repeat the loop for the case when the second beam is the photon emitter. 
        // Don't repeat for pA
        if( !( (A_2 == 1 && A_1 != 1) || (A_1 == 1 && A_2 != 1) ) ){ 
	  for(J=0;J<=(NY-1);J++){
    
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
    
                beam = 1; 
		ega1  = 0.5*W*exp(-y1);
		ega2  = 0.5*W*exp(-y2);
		ega12 = 0.5*W*exp(-y12);
      
		if(ega2 < Eth || ega1 < Eth)   
			continue;
		if(ega1 > maxPhotonEnergy() || ega2 > maxPhotonEnergy() ) 
			continue;

		csgA1=getcsgA(ega1,W,beam);
    
		// Middle Point                      =====>>>
		csgA12=getcsgA(ega12,W,beam); 

		// Second Point                      =====>>>
		csgA2=getcsgA(ega2,W,beam);
    
		// Sum the contribution for this W,Y.
		dR  = ega1*photonFlux(ega1,beam)*csgA1;
		dR  = dR + 4.*ega12*photonFlux(ega12,beam)*csgA12;
		dR  = dR + ega2*photonFlux(ega2,beam)*csgA2;
		dR  = dR*(dY/6.);

		dsigdyVal[J] += 10.*dR/dY; //This is dsig/dy in millibarn
                Egamma2Val[J] = ega12;
		ngamma2Val[J] = ega12*photonFlux(ega12,beam);
		sigGA2Val[J] = 10.*csgA12; //This is sigma(gamma+A) in millibarn 
		dsigdy2Val[J] = 10.*dR/dY; 

		// cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12,beam)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;

		int_r = int_r+dR;
	  }
        }
	cout<<endl;
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
	    printf("First column gives the rapidity. Second and third column give the photon energy \n");
	    printf("and photon flux (k*(dn/dk)) from the first beam (defined by BEAM_1_Z, BEAM_1_A  \n");
	    printf("etc.). Photon energies are given in GeV in the center of mass frame. The fourth \n"); 
            printf("column gives the gamma+A cross section with the second beam as target. The      \n");
	    printf("fifth column gives the contribution to the cross section (dsig/dy) for this     \n");
	    printf("combination. The 6th - 9th columns give the corresponding information for the   \n");
	    printf("other configuration (the beam defined by BEAM_2_Z etc. emits the photon). The   \n");
	    printf("last (10th) column gives the total dsig/dy. \n");
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
  	  // cout<<" Scale Factor: "<<scaleFactor<<endl;	  

          if( _printDef == 2 ){
            printf(" y      Egamma1     k*(dn_1/dk) sigma_2(gam+A) dsig_1/dy  Egamma2     k*(dn_2/dk) sigma_1(gam+A) dsig_2/dy    dsig/dy \n");
	  }
	  
	  for(J=0;J<=(NY-1);J++){
            if( _printDef == 1){ 
              printf("%+6.2f          %10.4f \n",yVal[J],scaleFactor*dsigdyVal[J]);
	    } else if ( _printDef == 2 ){
              printf("%+5.1f   %.4E  %.4E  %.4E    %9.4f   ",yVal[J],Egamma1Val[J],ngamma1Val[J],scaleFactor*sigGA1Val[J],scaleFactor*dsigdy1Val[J]); 
              printf("%.4E  %.4E  %.4E    %9.4f    %9.4f \n",Egamma2Val[J],ngamma2Val[J],scaleFactor*sigGA2Val[J],scaleFactor*dsigdy2Val[J],scaleFactor*dsigdyVal[J]); 
	    }
	  }
          printf("\n");
	}

}
