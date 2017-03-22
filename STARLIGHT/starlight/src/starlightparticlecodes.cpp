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
// $Rev:: 90                          $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2012-06-04 16:21:17 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include "starlightparticlecodes.h"


int starlightParticleCodes::jetsetToGeant(int particleid)
{
//this function takes a jetset particle number and converts it to GEANT
    int jtogx =0;
    if(particleid==22) jtogx = 1;
    if(particleid==-11) jtogx = 2;
    if(particleid==11) jtogx = 3;
    if(particleid==12||particleid==-12||particleid==14||particleid==-14||particleid==16||particleid==-16) jtogx = 4;
    if(particleid==-13) jtogx = 5;
    if(particleid==13) jtogx = 6;
    if(particleid==111) jtogx = 7;
    if(particleid==211) jtogx = 8;
    if(particleid==-211) jtogx = 9;
    if(particleid==130) jtogx = 10;
    if(particleid==321) jtogx = 11;
    if(particleid==-321) jtogx = 12;
    if(particleid==2112) jtogx = 13;
    if(particleid==2212) jtogx = 14;
    if(particleid==-2212) jtogx = 15;
    if(particleid==310) jtogx = 16;
    if(particleid==-310) jtogx = 10;
    if(particleid==221) jtogx = 17;
    if(particleid==3122) jtogx = 18;
    if(particleid==3222) jtogx = 19;
    if(particleid==3212) jtogx = 20;
    if(particleid==3112) jtogx = 21;
    if(particleid==3322) jtogx = 22;
    if(particleid==3312) jtogx = 23;
    if(particleid==3334) jtogx = 24;
    if(particleid==-2112) jtogx = 25;
    if(particleid==-3122) jtogx = 26;
    if(particleid==-3222) jtogx = 27;
    if(particleid==-3212) jtogx = 28;
    if(particleid==-3112) jtogx = 29;
    if(particleid==-3322) jtogx = 30;
    if(particleid==-3312) jtogx = 31;
    if(particleid==-3334) jtogx = 32;
    return jtogx;
}
