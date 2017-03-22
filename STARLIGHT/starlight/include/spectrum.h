
/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify
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

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include "randomgenerator.h"
class beamBeamSystem;
//class randomGenerator;

class spectrum
{
public:

    /** Spectrum must be constructed with beam-beam system, default constructor disallowed */
    spectrum(const randomGenerator &randy, beamBeamSystem *bbs);

    /**
    * Generate a table of photon energy probabilities
    * Use NK+1 logarithmic steps between Et_min and Eg_max
    */
    int generateKsingle();

    /**
    * Generate a 2-D table of photon energy probabilities
    * Use NK+1 x NK+1 logarithmic steps between Et_min and Eg_max
    */
    int generateKdouble();

    /**
    * Get the energy of a single gamma
    * @return energy of the gamma
    */
    double drawKsingle();

    /**
    * Get the energy of a single gamma
    * @param egamma1 variable passed by reference to get the energy of the frst gamma
    * @param egamma2 variable passed by reference to get the energy of the second gamma
    * @return energy of the gamma
    */
    void drawKdouble(float &egamma1, float &egamma2);

    /** Set the beam beam system */
    void setBeamBeamSystem(beamBeamSystem *bbs) {
        _beamBeamSystem = bbs;
    }

    /** Set the minimum gamma energy */
    void setMinGammaEnergy(double energy) { _eGammaMin = energy; }

    /** Set the maximum gamma energy */
    void setMaxGammaEnergy(double energy) { _eGammaMax = energy; }

    /** Set minimum impact parameter */
    void setBmin(double bmin) { _bMin = bmin; }

    /** Set maximum impact parameter */
    void setBmax(double bmax) { _bMax = bmax; }
    
protected:

    /** Generate the hadron breakup probability table */
    virtual bool generateBreakupProbabilities();

    /** Needs some explanation */
    virtual double getSigma(double /*egamma*/) const {
        return 1.05;
    }
    
    virtual double getTransformedNofe(double egamma, double b);
    
    /** Minimum impact parameter */
    double _bMin;

    /** Maximum impact parameter */
    double _bMax;
    
    /** Number of bins in impact parameter */
    int _nBbins;
    
    /** Vector containing the probability of breakup */
    std::vector<double> _probOfBreakup;
    
    /** Beam beam system */
    beamBeamSystem *_beamBeamSystem;

private:
    double getFnSingle(double egamma) const;

    double getFnDouble(double egamma1, double egamma2) const;

    /** NK */
    int _nK;

    /** Contains the 1 photon probabilities */
    std::vector<double> _fnSingle;

    /** Contains the 2 photon probabilities */
    std::vector<std::vector<double> > _fnDouble;

    /** Contains the cumulative distribution */
    std::vector<double> _fnSingleCumulative;

    /** Contains the cumulative distribution */
    std::vector<std::vector<double> > _fnDoubleCumulative;

    /**  */
    std::vector<double> _fnDoubleInt;

    /** */
    std::vector<double> _fnDoubleIntCumulative;
    
    /** Vecotr of gamma energies */
    std::vector<double> _eGamma;

    /** Min gamma energy */
    double _eGammaMin;

    /** Max gamma energy */
    double _eGammaMax;

    /** Z of target */
    int _zTarget;

    /** A of target */
    int _aTarget;

    /** Hadron breakup probability is calculated */
    bool _hadBreakProbCalculated;

    /** random number generator **/
    randomGenerator _randy;

    /** Default constructed disallowed (not implemented) */
    spectrum();

};

#endif // SPECTRUM_H
