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

// ** added code to integrate over a fixed impact parameter range, from bmin to bmax  SRK 7/2017
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 291                         $: revision of last commit
// $Author:: srklein                  $: author of last commit
// $Date:: 2017-09-11 04:55:34 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef TWOPHOTONLUMINOSITY_H
#define TWOPHOTONLUMINOSITY_H

#include "nucleus.h"
#include "beam.h"
#include "beambeamsystem.h"
#include "starlightlimits.h"


class twoPhotonLuminosity : public beamBeamSystem
{
public:
    twoPhotonLuminosity(const inputParameters& input, beam beam_1, beam beam_2);
    ~twoPhotonLuminosity();

protected:

   

private:
    struct difflumiargs
    {
	twoPhotonLuminosity *self;
        double m;
        double y;
        double res;
    };
    void twoPhotonDifferentialLuminosity(double bmina, double bmaxa);
    double D2LDMDY(double M,double Y,double &Normalize);   //simple analytic calculation, for XSEC_METHOD==0
    double D2LDMDY(double M,double Y) const;               //numerical calculation for XSEC_METHOD==1
    double D2LDMDYBRANGE(double M,double Y, double bmin, double bmax) const;        //numerical calculation for a fixed ion-ion impact parameter range
    static void * D2LDMDY_Threaded(void *a);

    double integral(double Normalize);
    double radmul(int N,double *Lower,double *Upper,int NIterMin,int NIterMax,double EPS,double *WK,int NIter,double &Result,double &ResErr,double &NFNEVL,double &Summary);
    double integrand(double N,double X[15]);
    double Nphoton(double W,double gamma,double Rho);

    double _W1; //Energy of photon #1
    double _W2; //Energy of photon #2
    double _gamma; //Gamma of the system
    
    const unsigned int _nWbins;
    const unsigned int _nYbins;
    
    const double _wMin;
    const double _yMin;
    const double _wMax;
    const double _yMax;
    const int _productionMode;
    const int _beamBreakupMode;
    const int _interferenceEnabled;
    const double _interferenceStrength;
    const double _maxPtInterference;
    const int _nmbPtBinsInterference;
    const int _xsecCalcMethod;
    const std::string _baseFileName;
    double _bmin,_bmax;  //SRK b range for peripheral collisions
};


#endif  // TWOPHOTONLUMINOSITY_H
