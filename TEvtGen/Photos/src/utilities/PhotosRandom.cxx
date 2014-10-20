#include <iostream>
#include "PhotosRandom.h"
#include "Photos.h"
#include "Log.h"

namespace Photospp
{

bool         PhotosRandom::init    = false;
int          PhotosRandom::iseed[2]= { 1802, 9373 };
int          PhotosRandom::i97     = 96;
int          PhotosRandom::j97     = 32;
double       PhotosRandom::uran[97]= { 0.0 };
double       PhotosRandom::cran    = 362436.0  /16777216.0;
const double PhotosRandom::cdran   = 7654321.0 /16777216.0;
const double PhotosRandom::cmran   = 16777213.0/16777216.0;

void PhotosRandom::setSeed(int s1,int s2)
{
	if(s1<0 || s1>31327) Log::Fatal("PhotosRandom::setSeed(): Seed(1) out of range [0,31327]",8);
	if(s2<0 || s2>30080) Log::Fatal("PhotosRandom::setSeed(): Seed(2) out of range [0,30080]",9);
	iseed[0]=s1;
	iseed[1]=s2;
}

/*******************************************************************************
  PHORIN:   PHOton radiation  in decays RANdom number generator init

  Purpose:  Initialse PHORAN  with  the user  specified seeds in the
            array iseed.  For details  see also:  F. James  CERN DD-
            Report November 1988.

  Author(s):  B. van Eijk and F. James        Created at:  27/09/89
                                              Last Update: 22/02/90
                                         Rewritten to C++: 18/10/10
                                         by T. Przedzinski  (tprzedzi@cern.ch)
*******************************************************************************/
void PhotosRandom::initialize()
{
	long IS1,IS2,IS3,IS4,IS5;
	double S,T;

// Calculate Marsaglia and Zaman seeds (by F. James)
	IS1=(iseed[0]/177)%177+2;
	IS2= iseed[0]%177+2;
	IS3=(iseed[1]/169)%178+1;
	IS4= iseed[1]%169;
	for(int i=0;i<97;i++)
	{
		S=0.0;
		T=0.5;
		for(int j=0;j<24;j++)
		{
			IS5=( ((IS1*IS2)%179)*IS3 )%179;
			IS1=IS2;
			IS2=IS3;
			IS3=IS5;
			IS4=(53*IS4+1)%169;
			if( (IS4*IS5)%64>=32) S=S+T;
			T=0.5*T;
		}
		uran[i]=S;
	}
	init=true;
	Log::Debug(0)<<"PhotosRandom::inititalize(): seed: "<<iseed[0]<<", "<<iseed[1]<<std::endl;
}

/*******************************************************************************
  PHORAN:   PHOton radiation in decays ret number generator based
            on Marsaglia Algorithm

  Purpose:  Generate  uniformly  distributed  random numbers between
            0 and 1.  Super long period:  2**144.  See also:
            G. Marsaglia and A. Zaman,  FSU-SCR-87-50,  for seed mo-
            difications  to  this version  see:  F. James DD-Report,
            November 1988.  The generator  has  to be initialized by
            a call to PHORIN ( C++ version: initialize() ).

  Author(s):  B. van Eijk, G. Marsaglia and   Created at:  27/09/89
              A. Zaman                        Last Update: 27/09/89
                                         Rewritten to C++: 18/10/10
                                         by T. Przedzinski  (tprzedzi@cern.ch)
*******************************************************************************/
double PhotosRandom::randomReal()
{
	if(!init) Log::Fatal("PhotosRandom::randomReal(): generator not initialized",1);
	double ret=0.0;
	while(true)
	{
		ret = uran[i97]-uran[j97];
		if(ret<0.0) ret+=1.;
		uran[i97]=ret;
		i97--;
		if(i97<0) i97=96;
		j97--;
		if(j97<0) j97=96;
		cran-=cdran;
		if(cran<0.0) cran+=cmran;
		ret-=cran;
		if(ret<0.0) ret+=1.0;
		if(ret>0.0) break;
	}
	return ret;
}

} // namespace Photospp

