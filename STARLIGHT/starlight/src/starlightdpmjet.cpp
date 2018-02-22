/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "starlightdpmjet.h"
#include "spectrum.h"
#include <iostream>
#include <spectrumprotonnucleus.h>
#include "starlightconfig.h"


extern "C"
{
    extern struct
        {

            double slpx[200000];
            double slpy[200000];
            double slpz[200000];
            double sle[200000];
            double slm[200000];
            int slpid[200000];
            int slcharge[200000];

        } dpmjetparticle_;

    void dt_produceevent_(float* gammaE, int* nparticles);
    void dt_getparticle_(int *ipart, int *res);
    void dt_initialise_();
}

starlightDpmJet::starlightDpmJet(const inputParameters& inputParametersInstance,randomGenerator* randy,beamBeamSystem& beamsystem ) : eventChannel(inputParametersInstance,randy,beamsystem)
        ,_spectrum(0)
        ,_doDoubleEvent(true)
	,_minGammaEnergy(6.0)
	,_maxGammaEnergy(600000.0)
	,_protonMode(false)
{

}

int starlightDpmJet::init()
{
   if(_protonMode)
   {
    _spectrum = new spectrumProtonNucleus(_randy,&_bbs);
   }
   else
   {
    _spectrum = new spectrum(_randy,&_bbs);
   }

   _spectrum->setMinGammaEnergy(_minGammaEnergy);
   _spectrum->setMaxGammaEnergy(_maxGammaEnergy);
   
   if(!_doDoubleEvent)
   {
    _spectrum->generateKsingle();
   }
   else 
   {
    _spectrum->generateKdouble();
   }

   return 0;

}


upcEvent starlightDpmJet::produceEvent()
{

    upcEvent event;

    if (!_doDoubleEvent)
    {
      int zdirection = 1;
        float gammaE = _spectrum->drawKsingle();
        event = produceSingleEvent(zdirection, gammaE);
	//        std::cout << "Gamma energy: " << gammaE << std::endl;
    }
    else
    {
        event = produceDoubleEvent();
    }

    return event;
}

upcEvent starlightDpmJet::produceSingleEvent(int zdirection, float gammaE)
{

    upcEvent event;
    event.addGamma(gammaE);

    int nParticles = 0;

    dt_produceevent_(&gammaE, &nParticles);
    

    //In which direction do we go?
    double rapidity = _bbs.beam1().rapidity()*zdirection;

    for (int i = 0; i < nParticles; i++)
    {
        starlightParticle particle(dpmjetparticle_.slpx[i], dpmjetparticle_.slpy[i], zdirection*dpmjetparticle_.slpz[i], dpmjetparticle_.sle[i], dpmjetparticle_.slm[i], dpmjetparticle_.slpid[i], dpmjetparticle_.slcharge[i]);
	vector3 boostVector(0, 0, tanh(-rapidity));
	particle.Boost(boostVector);
        event.addParticle(particle);
    }
    return event;
}

upcEvent starlightDpmJet::produceDoubleEvent()
{
    upcEvent event;

    float gammaE1 = 0.0;
    float gammaE2 = 0.0;

    _spectrum->drawKdouble(gammaE1, gammaE2);

//    std::cout << "Gamma1 energy: " << gammaE1 << std::endl;
    //std::cout << "Gamma2 energy: " << gammaE2 << std::endl;
    
    //In which direction do we go?
    int zdirection = (_randy->Rndom()) < 0.5 ? -1 : 1;

    event = produceSingleEvent(zdirection, gammaE1);

    zdirection = zdirection *-1;

    event = event + produceSingleEvent(zdirection, gammaE2);

    return event;
}


