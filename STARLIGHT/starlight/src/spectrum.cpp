
/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

p    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "spectrum.h"
#include <cmath>
#include "beambeamsystem.h"
#include <randomgenerator.h>
#include <iostream>

spectrum::spectrum(const randomGenerator &randy, beamBeamSystem *bbs) :
	 _bMin(5.0)
        ,_bMax(128000.0)
	,_nBbins(6400)
	,_probOfBreakup(_nBbins)
        ,_beamBeamSystem(bbs)
	,_nK(10000)
        ,_fnSingle(_nK)
        ,_fnDouble(_nK)
        ,_fnSingleCumulative(_nK+1)
        ,_fnDoubleCumulative(_nK+1)
        ,_fnDoubleInt(_nK)
        ,_fnDoubleIntCumulative(_nK+1)
        ,_eGamma(_nK+1)
        ,_eGammaMin(6.0)
        ,_eGammaMax(600000.0)
        ,_zTarget(82)
        ,_aTarget(278)
        ,_hadBreakProbCalculated(false)
	,_randy(randy)
{
    _eGamma.resize(_nK+1);
    _probOfBreakup.resize(_nBbins);
}

int spectrum::generateKsingle()
{

    _fnSingle.resize(_nK);
    _fnSingleCumulative.resize(_nK+1);

    double eg_inc = exp(log(_eGammaMax/_eGammaMin)/(double)_nK);
    
    double egamma =  _eGammaMin;
    for (int i = 0; i < _nK+1; i++)
    {
        _eGamma[i] = egamma;
        egamma = egamma * eg_inc;
    }
    egamma = _eGammaMin;

    double fnorm = 0;


    if (_hadBreakProbCalculated == false)
    {
        _hadBreakProbCalculated = generateBreakupProbabilities();
    }
    double binc = exp((log(_bMax/_bMin))/(double)_nBbins);

    for (int i = 0; i < _nK; i++)
    {
        double b = _bMin;

        double bint = 0.0;

        double f1 = 0;
        double f2 = 0;

        for (int j = 0; j < _nBbins - 1; j++)
        {
            double bold = b;
            if (j == 0)
            {
		f1 = getTransformedNofe(egamma, b)*getSigma(egamma)*_probOfBreakup[j]*b;
		//std::cout << fProbOfBreakup[j] << std::endl;
            }
            else
            {
                f1 = f2;
            }
            b = b*binc;
            f2 = getTransformedNofe(egamma, b)*getSigma(egamma)*_probOfBreakup[j+1]*b;;
            bint = bint + 0.5*(f1+f2)*(b-bold);
        }
        bint = 2.0*starlightConstants::pi*bint;
        if (i == 0)
        {
            fnorm = 1.0/bint;
        }
        _fnSingle[i] = bint*(_eGamma[i+1]-_eGamma[i]);

        egamma = egamma*eg_inc;
    }

    _fnSingleCumulative[0] = 0.00;
    for (int i = 0; i < _nK; i++)
    {
        _fnSingleCumulative[i+1] = _fnSingleCumulative[i]+_fnSingle[i];
    }

    double fnormfactor = 1.00/_fnSingleCumulative[_nK];
    for (int i = 0; i < _nK; i++)
    {
        _fnSingleCumulative[i+1] = fnormfactor*_fnSingleCumulative[i+1];
    }
    
    return 0;

}

int spectrum::generateKdouble()
{
    //Quick fix for now TODO: Fix it!
    _nK = 100;

    _fnDouble.resize(_nK);
    _fnDoubleInt.resize(_nK);
    _fnDoubleIntCumulative.resize(_nK+1);
    _fnDoubleCumulative.resize(_nK+1);
    for (int i = 0; i < _nK; i++)
    {
        _fnDouble[i].resize(_nK);
        _fnDoubleCumulative[i].resize(_nK+1);
    }
    _fnDoubleCumulative[_nK].resize(_nK+1);

    double eg_inc = exp(log(_eGammaMax/_eGammaMin)/(double)_nK);
    double egamma1 =  _eGammaMin;
    double egamma2 =  _eGammaMin;

    for (int i = 0; i < _nK+1; i++)
    {
        _eGamma[i] = egamma1;
        egamma1 = egamma1 * eg_inc;
    }
    egamma1 = _eGammaMin;

    double fnorm = 0;

    if (_hadBreakProbCalculated == false)
    {
        _hadBreakProbCalculated = generateBreakupProbabilities();
    }

    double binc = exp((log(_bMax/_bMin))/(double)_nBbins);

    int nbbins = _nBbins;

    for (int i = 0; i < _nK; i++)
    {

        egamma2 = _eGammaMin;

        for (int j = 0; j < _nK; j++)
        {
            double bint = 0.0;
            double b = _bMin;
            double f1 = 0;
            double f2 = 0;

            for (int k = 0; k < nbbins - 1; k++)
            {
                double bold = b;

                if (k == 0)
                {
                  f1 = getTransformedNofe(egamma1, b) * getTransformedNofe(egamma2, b) 
                         * getSigma(egamma1) * getSigma(egamma2) *_probOfBreakup[k]*b;                }
                else
                {
                    f1 = f2;
                }
                b = b*binc;
                f2 = getTransformedNofe(egamma1, b) * getTransformedNofe(egamma2, b) 
                     * getSigma(egamma1) * getSigma(egamma2) *_probOfBreakup[k+1]*b;
                bint = bint + 0.5*(f1+f2)*(b-bold);
            }
            bint = 2.0*starlightConstants::pi*bint;
            _fnDouble[i][j] = bint * (_eGamma[i+1] - _eGamma[i]) * (_eGamma[j+1] - _eGamma[j]);
            egamma2 = egamma2 * eg_inc;
        }
        egamma1 = egamma1 * eg_inc;
    }

    for (int i = 0; i < _nK; i++)
    {
        _fnDoubleInt[i] = 0.0;
        for (int j = 0; j < _nK; j++)
        {
            _fnDoubleInt[i] = _fnDoubleInt[i] + _fnDouble[i][j];
        }
    }

    _fnDoubleIntCumulative[0] = 0.0;
    for (int i = 1; i < _nK+1; i++)
    {
        _fnDoubleIntCumulative[i] = _fnDoubleIntCumulative[i-1] + _fnDoubleInt[i-1];
    }

    fnorm = 1.0/_fnDoubleIntCumulative[_nK];
    for (int i = 0; i < _nK+1; i++)
    {
        _fnDoubleIntCumulative[i] = fnorm * _fnDoubleIntCumulative[i];
    }

    return 0;
}

double spectrum::drawKsingle()
{
    double xtest = 0;
    int itest = 0;
    double egamma = 0.0;

    xtest = _randy.Rndom();
    while (xtest > _fnSingleCumulative[itest])
    {
        itest++;
    }
    itest = itest - 1;

    if (itest >= _nK || itest < 0)
    {
        std::cerr << "ERROR: itest: " << itest << std::endl;

    }

    double delta_f = xtest - _fnSingleCumulative[itest];
    if (delta_f <= 0.0)
    {
        std::cout << "WARNING: delta_f: " << delta_f << std::endl;
        return -1;
    }
    double dE = _eGamma[itest+1] - _eGamma[itest];
    double dF = _fnSingleCumulative[itest+1] - _fnSingleCumulative[itest];

    double delta_e = delta_f*dE/dF;

    if (delta_e > (_eGamma[itest+1] - _eGamma[itest]))
    {
        std::cerr << "ERROR: delta_E: " << delta_e << std::endl;
    }
   
    egamma = _eGamma[itest] + delta_e;
    return egamma;
}

void spectrum::drawKdouble(float& egamma1, float& egamma2)
{
    double xtest1 = 0.0;
    double xtest2 = 0.0;
    int itest1 = 0;
    int itest2 = 0;

    xtest1 = _randy.Rndom();

    while (xtest1 > _fnDoubleIntCumulative[itest1])
    {
        itest1++;
    }
    itest1 = itest1 - 1;

    if (itest1 >= _nK || itest1 < 0)
    {
        std::cerr << "ERROR: itest1: " << itest1 << std::endl;
    }
    double delta_f = xtest1 - _fnDoubleIntCumulative[itest1];

    if (delta_f <= 0.0)
    {
        std::cout << "WARNING: delta_f: " << delta_f << std::endl;
    }

    double dE = _eGamma[itest1+1] - _eGamma[itest1];
    double dF = _fnDoubleIntCumulative[itest1+1] - _fnDoubleIntCumulative[itest1];

    double delta_e = delta_f*dE/dF;

    if (delta_e > (_eGamma[itest1+1] - _eGamma[itest1]))
    {
        std::cerr << "ERROR: delta_E: " << delta_e << std::endl;
    }

    egamma1 = _eGamma[itest1] + delta_e;
    
    // Second gamma

    std::vector<double> fn_second_cumulative(_nK+1);
    
    fn_second_cumulative[0] = 0.0;
    for(int i = 1; i < _nK+1; i++)
    {
       fn_second_cumulative[i] = fn_second_cumulative[i-1] + _fnDouble[itest1][i-1]; 
    }
    
    double norm_factor = 1.0/fn_second_cumulative[_nK];
    for(int i = 0; i < _nK+1; i++)
    {
      fn_second_cumulative[i] = norm_factor*fn_second_cumulative[i];
    }
    
    xtest2 = _randy.Rndom();

    while (xtest2 > fn_second_cumulative[itest2])
    {
        itest2++;
    }
    itest2 = itest2 - 1;

    if (itest2 >= _nK || itest2 < 0)
    {
        std::cerr << "ERROR: itest2: " << itest2 << std::endl;
    }
    delta_f = xtest2 - fn_second_cumulative[itest2];

    if (delta_f <= 0.0)
    {
        std::cout << "WARNING: delta_f: " << delta_f << std::endl;
    }

    dE = _eGamma[itest2+1] - _eGamma[itest2];
    dF = fn_second_cumulative[itest2+1] - fn_second_cumulative[itest2];

    delta_e = delta_f*dE/dF;

    if (delta_e > (_eGamma[itest2+1] - _eGamma[itest2]))
    {
        std::cerr << "ERROR: delta_E: " << delta_e << std::endl;
    }

    egamma2 = _eGamma[itest2] + delta_e;
    
    return;
}


bool spectrum::generateBreakupProbabilities()
{

    int nbbins = _nBbins;

    double b_min = _bMin;
    double binc = exp((log(_bMax/_bMin))/(double)_nBbins);

    double b = b_min;

    _probOfBreakup.resize(nbbins);

    for (int i = 0; i < nbbins; i++)
    {
        double bimp = b;
        double rhad = 0;
        rhad = _beamBeamSystem->probabilityOfBreakup(bimp);
        _probOfBreakup[i] = exp(-rhad);
        b = b*binc;
    }
    return true;
}

double spectrum::getFnSingle(double egamma) const
{
    double eginc = exp(log(_eGammaMax/_eGammaMin)/static_cast<double>(_nK));
    int i1 = log(egamma/_eGammaMin)/log(eginc);
    int i2 = i1 + 1;
    double fnSingle = 0.0;

    if (i1 < 0 || i2 > _nK)
    {
        std::cout << "I1, I2 out of bounds. Egamma = " << egamma << std::endl;
        std::cout << "I1, I2 = " << i1 << ", " << i2 << std::endl;
        fnSingle = 0.0;
    }
    else
    {
        double dE = _eGamma[i2] - _eGamma[i1];
        double eFrac = (egamma - _eGamma[i1])/dE;

        if (eFrac < 0.0 || eFrac > 1.0)
        {
            std::cout << "WARNING: Efrac = " << eFrac << std::endl;
        }
        fnSingle = (1.0 - eFrac)*_fnSingle[i1] + eFrac*_fnSingle[i2];
    }
    return fnSingle;
}

double spectrum::getFnDouble(double egamma1, double egamma2) const
{
    double eginc = exp(log(_eGammaMax/_eGammaMin)/static_cast<double>(_nK));
    int i1 = log(egamma1/_eGammaMin)/log(eginc);
    int i2 = i1 + 1;
    double fnDouble = 0.0;

    if (i1 < 0 || i2 > _nK)
    {
        std::cout << "I1, I2 out of bounds. Egamma1 = " << egamma1 << std::endl;
        std::cout << "I1, I2 = " << i1 << ", " << i2 << std::endl;
        fnDouble = 0.0;
        return fnDouble;
    }

    int j1 = log(egamma2/_eGammaMin)/log(eginc);
    int j2 = j1 + 1;

    if (j1 < 0 || j2 > _nK)
    {
        std::cout << "J1, J2 out of bounds. Egamma2 = " << egamma2 << std::endl;
        std::cout << "J1, J2 = " << j1 << ", " << j2 << std::endl;
        fnDouble = 0.0;
        return fnDouble;
    }

    double dE1 = _eGamma[i2] - _eGamma[i1];
    double eFrac1 = (egamma1 - _eGamma[i1])/dE1;

    if (eFrac1 < 0.0 || eFrac1 > 1.0)
    {
        std::cout << "WARNING: Efrac1 = " << eFrac1 << std::endl;
    }

    double ptemp1 = (1.0 - eFrac1)*_fnDouble[i1][j1] + eFrac1*_fnDouble[i2][j1];
    double ptemp2 = (1.0 - eFrac1)*_fnDouble[i1][j2] + eFrac1*_fnDouble[i2][j2];

    double dE2 = _eGamma[j2] - _eGamma[j1];
    double eFrac2 = (egamma2 - _eGamma[j1])/dE2;

    if (eFrac2 < 0.0 || eFrac2 > 1.0)
    {
        std::cout << "WARNING: Efrac2 = " << eFrac2 << std::endl;
    }

    fnDouble = (1.0 - eFrac2)*ptemp1 + eFrac2*ptemp2;

    return fnDouble;

}

double spectrum::getTransformedNofe(double egamma, double b)
{
   double factor = 1.0/(2.0*_beamBeamSystem->beamLorentzGamma());
   double res = factor * _beamBeamSystem->beam1().photonDensity(b, egamma*factor);
   
   return res;
}




