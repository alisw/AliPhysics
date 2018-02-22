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
// $Rev:: 298                         $: revision of last commit
// $Author:: srklein                  $: author of last commit
// $Date:: 2018-02-22 00:23:57 +0100 #$: date of last commit
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
#include "reportingUtils.h"
#include "starlightconstants.h"
#include "bessel.h"
#include "beambeamsystem.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
beamBeamSystem::beamBeamSystem(const inputParameters& inputParametersInstance,
			       const beam&            beam1,
                               const beam&            beam2)
  : _ip(&inputParametersInstance),
    _beamLorentzGamma(inputParametersInstance.beamLorentzGamma()),
    _beamLorentzGamma1(inputParametersInstance.beam1LorentzGamma()),
    _beamLorentzGamma2(inputParametersInstance.beam2LorentzGamma()),
    _beamBreakupMode (inputParametersInstance.beamBreakupMode()),
    _beam1           (beam1),
    _beam2           (beam2),
    _breakupProbabilities(0),
    _breakupImpactParameterStep(1.007),
    _breakupCutOff(10e-6)
{ 
  init();
}
  



//______________________________________________________________________________
beamBeamSystem::beamBeamSystem(const inputParameters& inputParametersInstance)
	: _beamLorentzGamma(inputParametersInstance.beamLorentzGamma()),
          _beamLorentzGamma1(inputParametersInstance.beam1LorentzGamma()),
          _beamLorentzGamma2(inputParametersInstance.beam2LorentzGamma()),
	  _beamBreakupMode (inputParametersInstance.beamBreakupMode()),
	  _beam1           (inputParametersInstance.beam1Z(),
	                    inputParametersInstance.beam1A(),
			    inputParametersInstance.productionMode(),
                            inputParametersInstance.beam1LorentzGamma()),
	  _beam2           (inputParametersInstance.beam2Z(),
	                    inputParametersInstance.beam2A(),
			    inputParametersInstance.productionMode(),
                            inputParametersInstance.beam2LorentzGamma()),
	  _breakupProbabilities(0),
	  _breakupImpactParameterStep(1.007),
	  _breakupCutOff(10e-10)
{
  init();
}



//______________________________________________________________________________
beamBeamSystem::~beamBeamSystem()
{ }

void beamBeamSystem::init()
{
   // Calculate beam gamma in CMS frame
   double rap1 = acosh(_beamLorentzGamma1);
   double rap2 = -acosh(_beamLorentzGamma2);
   
   _cmsBoost = (rap1+rap2)/2.;

   _beamLorentzGamma = cosh((rap1-rap2)/2);
   _beam1.setBeamLorentzGamma(_beamLorentzGamma);
   _beam2.setBeamLorentzGamma(_beamLorentzGamma);
   
   generateBreakupProbabilities();
}
//______________________________________________________________________________
double
beamBeamSystem::probabilityOfBreakup(const double D) const
{
	
	double bMin = (_beam1.nuclearRadius()+_beam2.nuclearRadius())/2.;
	double pOfB = 0.; // PofB = 1 means that there will be a UPC event and PofB = 0 means no UPC

	// Do pp here
	if ((_beam1.Z() == 1) && (_beam1.A() == 1) && (_beam2.Z() == 1) && (_beam2.A() == 1)) {  
		double ppslope=19.8;
		double GammaProfile = exp(-D * D / (2. * hbarc * hbarc * ppslope));
		pOfB = (1. - GammaProfile) * (1. - GammaProfile);
		return pOfB;
	}
	else if ( ( (_beam1.A() == 1) && (_beam2.A() != 1) ) || ((_beam1.A() != 1) && (_beam2.A() == 1)) ) {  
	  // This is pA
          if( _beam1.A() == 1 ){ 
            bMin = _beam2.nuclearRadius() + 0.7;
            pOfB = exp(-7.0*_beam2.rho0()*_beam2.thickness(D));
          }else if( _beam2.A() == 1 ){ 
            bMin = _beam1.nuclearRadius() + 0.7; 
            pOfB = exp(-7.0*_beam1.rho0()*_beam1.thickness(D));
          }else{
            cout<<"Some logical problem here!"<<endl;
          }
          return pOfB;
          
	}

	//use the lookup table and return...
	pOfB = 1.;
	if (D > 0.0) {             
		//Now we must determine which step number in d corresponds to this D,
		// and use appropiate Ptot(D_i)
		int i = (int)(log(D / bMin) / log(_breakupImpactParameterStep));
		if (i <= 0)
			pOfB = _breakupProbabilities[0];
		else{
			if (i >= int(_breakupProbabilities.size()-1))
				pOfB = _breakupProbabilities[_breakupProbabilities.size()-1];
			else {
				const double DLow = bMin * pow((_breakupImpactParameterStep), i);
				const double DeltaD = (_breakupImpactParameterStep-1) * DLow;
				const double DeltaP = _breakupProbabilities[i + 1] - _breakupProbabilities[i];
				pOfB   = _breakupProbabilities[i] + DeltaP * (D - DLow) / DeltaD;
			}
		}
	}

	return pOfB;
}

void
beamBeamSystem::generateBreakupProbabilities()
{

    double bMin = (_beam1.nuclearRadius()+_beam2.nuclearRadius())/2.;
    
    // Do this only for nucleus-nucleus collisions.
    // pp and pA are handled directly in probabilityOfBreakup
//    if ((_beam1.Z() != 1) && (_beam1.A() != 1) && (_beam2.Z() != 1) && _beam2.A() != 1) {
//	this change allows deuterium and tritium to be handled here
    if ((_beam1.Z() != 1) && (_beam1.A() != 1) && _beam2.A() != 1) {

        if (_beamBreakupMode == 1)
            printInfo << "Hard Sphere Break criteria. b > " << 2. * _beam1.nuclearRadius() << endl;
        if (_beamBreakupMode == 2)
            printInfo << "Requiring XnXn [Coulomb] breakup. " << endl;
        if (_beamBreakupMode == 3)
            printInfo << "Requiring 1n1n [Coulomb only] breakup. " << endl;
        if (_beamBreakupMode == 4)
            printInfo << "Requiring both nuclei to remain intact. " << endl;
        if (_beamBreakupMode == 5)
            printInfo << "Requiring no hadronic interactions. " << endl;
        if (_beamBreakupMode == 6)
            printInfo << "Requiring breakup of one or both nuclei. " << endl;
        if (_beamBreakupMode == 7)
            printInfo << "Requiring breakup of one nucleus (Xn,0n). " << endl;

	double pOfB = 0;
	double b = bMin;
        double totRad = _beam1.nuclearRadius()+_beam2.nuclearRadius();
        
	while(1)
	{
            
            if(_beamBreakupMode != 5)
            {
                if(b > (totRad*1.5))
                {
                    if(pOfB<_breakupCutOff)
                    {
//                         std::cout << "Break off b: " << b << std::endl;
//                         std::cout << "Number of PofB bins: " << _breakupProbabilities.size() << std::endl;
                        break;
                    }
                }
            }
            else
            {
                if((1-pOfB)<_breakupCutOff)
                {
//                         std::cout << "Break off b: " << b << std::endl;
//                         std::cout << "Number of PofB bins: " << _breakupProbabilities.size() << std::endl;
                        break;
                }
            }
//             std::cout << 1-pOfBreakup << std::endl;

            probabilityOfHadronBreakup(b);
            probabilityOfPhotonBreakup(b, _beamBreakupMode);

            //What was probability of photonbreakup depending upon mode selection,
            // is now done in the photonbreakupfunction
            if (_beamBreakupMode == 1) {
                if (b >_beam1.nuclearRadius()+_beam2.nuclearRadius())  // symmetry
                    _pHadronBreakup = 0;
                else
                    _pHadronBreakup = 999.;
            }
            
            b *= _breakupImpactParameterStep;
	    pOfB = exp(-1 * _pHadronBreakup) * _pPhotonBreakup;
            _breakupProbabilities.push_back(pOfB);
        } // End while(1)
    }
    
}

//______________________________________________________________________________
double
beamBeamSystem::probabilityOfHadronBreakup(const double impactparameter)
{
	//probability of hadron breakup, 
	//this is what is returned when the function is called
	double gamma = _beamLorentzGamma; 
	//input for gamma_em
	double b = impactparameter;
	int a1 = _beam1.A();  
	int a2 = _beam2.A();  

	static int IFIRSTH = 0;
	static double DELL=0., DELR=0., SIGNN=0., R1=0., A1=0., A2=0., R2=0., RHO1=0.;
	static double RHO2=0., NZ1=0., NZ2=0., NR1=0., NR2=0.,RR1=0., RR2=0., NY=0., NX=0.;
	static double AN1=0., AN2=0.;
	double RSQ=0.,Z1=0.,Z2=0.,Y=0.,X=0.,XB=0.,RPU=0.,IRUP=0.,RTU=0.;
	double IRUT=0.,T1=0.,T2=0.;
	static double DEN1[20002], DEN2[20002];
	double energy,sigmainmb;
	if (IFIRSTH != 0) goto L100;
	//Initialize
	//Integration delta x, delta z
	IFIRSTH = 1;
	DELL   = .05;
	DELR   = .01;

	// replace this with a parameterization from the particle data book  SRK 4/1025
	//use two sigma_NN's. 52mb at rhic 100gev/beam, 88mb at LHC 2.9tev/beam, gamma is in cm system
	//SIGNN = 5.2;
	//if ( gamma > 500. ) SIGNN = 8.8; 
	energy=2*gamma*0.938;   // center of mass energy, in GeV
	  // This equation is from section 50 of the particle data book, the subsection on "Total Hadronic Cross-Sections, using the parameterization for sqrt{s} > 7 GeV.
	  // only the first and second terms contribute significantly, but leave them all here for good measure
	  sigmainmb = 0.2838*pow(log(energy),2)+33.73+13.67*pow(energy,-0.412)-7.77*pow(energy,-0.5626);
	  SIGNN=sigmainmb/10.;

	//use parameter from Constants
	R1 = ( _beam1.nuclearRadius());  
        R2 = ( _beam2.nuclearRadius());
	A1 = (_beam1.woodSaxonSkinDepth()); // take values from nucleus.cpp, since this parameter may now change  // SRK Feb. 2018
        A2 = (_beam2.woodSaxonSkinDepth()); 
	//A1 = 0.535; //This is woodsaxonskindepth
        //A2 = 0.535; 
	//write(6,12)r1,a1,signn  Here is where we could probably set this up asymmetrically R2=_beam2.nuclearRadius() and RHO2=ap2=_beam2.A()
	// R2 = R1;
	RHO1 = a1;
	RHO2 = a2;
	NZ1  = ((R1+5.)/DELR);
	NR1  = NZ1;
	NZ2  = ((R2+5.)/DELR);
	NR2  = NZ2;
	RR1  = -DELR;
        RR2  = -DELR; 
	NY   = ((R1+5.)/DELL);
	NX   = 2*NY;
	// This calculates T_A(b) for beam 1 and stores it in DEN1[IR1] 
	for ( int IR1 = 1; IR1 <= NR1; IR1++) {
		DEN1[IR1] = 0.;
		RR1       = RR1+DELR;
		Z1        = -DELR/2;

		for ( int IZ1 = 1; IZ1 <= NZ1; IZ1++) {
			Z1  = Z1+DELR;
			RSQ = RR1*RR1+Z1*Z1;
			DEN1[IR1] = DEN1[IR1]+1./(1.+exp((sqrt(RSQ)-R1)/A1));
		}

		DEN1[IR1] = DEN1[IR1]*2.*DELR;
	}

	// This calculates T_A(b) for beam 2 and stores it in DEN2[IR2] 
	for ( int IR2 = 1; IR2 <= NR2; IR2++) {
		DEN2[IR2] = 0.;
		RR2       = RR2+DELR;
		Z2        = -DELR/2;

		for ( int IZ2 = 1; IZ2 <= NZ2; IZ2++) {
			Z2  = Z2+DELR;
			RSQ = RR2*RR2+Z2*Z2;
			DEN2[IR2] = DEN2[IR2]+1./(1.+exp((sqrt(RSQ)-R2)/A2));
		}

		DEN2[IR2] = DEN2[IR2]*2.*DELR;
	}

	AN1 = 0.;
	RR1 = 0.;
        RR2 = 0.;
    
	for ( int IR1 =1; IR1 <= NR1; IR1++) {
		RR1 = RR1+DELR;
		AN1 = AN1+RR1*DEN1[IR1]*DELR*2.*starlightConstants::pi;
	}
	for ( int IR2 =1; IR2 <= NR2; IR2++) {
		RR2 = RR2+DELR;
		AN2 = AN2+RR2*DEN2[IR2]*DELR*2.*starlightConstants::pi;
	}
        

	//.1 to turn mb into fm^2
	//Calculate breakup probability here
 L100:
	_pHadronBreakup = 0.;
	if ( b > 25. ) return _pHadronBreakup;
	Y = -.5*DELL;
	for ( int IY = 1; IY <= NY; IY++) {
		Y = Y+DELL;
		X = -DELL*float(NY+1);

		for ( int IX = 1; IX <=NX; IX++) {
			X = X+DELL;
			XB = b-X;
			RPU = sqrt(X*X+Y*Y);
			IRUP = (RPU/DELR)+1;
			RTU  = sqrt(XB*XB+Y*Y);
			IRUT = (RTU/DELR)+1;
			T1   = DEN2[(int)IRUT]*RHO2/AN2;
			T2   = DEN1[(int)IRUP]*RHO1/AN1;
			//Eq.6 BCW, Baltz, Chasman, White, Nucl. Inst. & Methods A 417, 1 (1998)
			_pHadronBreakup=_pHadronBreakup+2.*T1*(1.-exp(-SIGNN*T2))*DELL*DELL;
		}//for(IX)
	}//for(IY)

	return _pHadronBreakup;
}


//______________________________________________________________________________
double
beamBeamSystem::probabilityOfPhotonBreakup(const double impactparameter, const int mode)
{
	static double ee[10001], eee[162], se[10001];

	_pPhotonBreakup =0.;   //Might default the probability with a different value?
	double b = impactparameter;
	int zp = _beam1.Z();  //What about _beam2? Generic approach?
	int ap = _beam1.A();
	
	//Was initialized at the start of the function originally, been moved inward.
	double pxn=0.;
	double p1n=0.;

	//Used to be done prior to entering the function. Done properly for assymetric?
	double gammatarg = 2.*_beamLorentzGamma*_beamLorentzGamma-1.;	
	double omaxx =0.;
	//This was done prior entering the function as well
	if (_beamLorentzGamma > 500.){
		omaxx=1.E10;
	}
	else{
		omaxx=1.E7;
	}


	double e1[23]= {0.,103.,106.,112.,119.,127.,132.,145.,171.,199.,230.,235.,
	                254.,280.,300.,320.,330.,333.,373.,390.,420.,426.,440.};
	double s1[23]= {0.,12.0,11.5,12.0,12.0,12.0,15.0,17.0,28.0,33.0,
	                52.0,60.0,70.0,76.0,85.0,86.0,89.0,89.0,75.0,76.0,69.0,59.0,61.0};
	double e2[12]={0.,2000.,3270.,4100.,4810.,6210.,6600.,
	               7790.,8400.,9510.,13600.,16400.};
	double s2[12]={0.,.1266,.1080,.0805,.1017,.0942,.0844,.0841,.0755,.0827,
	               .0626,.0740};
	double e3[29]={0.,26.,28.,30.,32.,34.,36.,38.,40.,44.,46.,48.,50.,52.,55.,
	               57.,62.,64.,66.,69.,72.,74.,76.,79.,82.,86.,92.,98.,103.};
	double s3[29]={0.,30.,21.5,22.5,18.5,17.5,15.,14.5,19.,17.5,16.,14.,
	               20.,16.5,17.5,17.,15.5,18.,15.5,15.5,15.,13.5,18.,14.5,15.5,12.5,13.,
	               13.,12.};
	static double sa[161]={0.,0.,.004,.008,.013,.017,.021,.025,.029,.034,.038,.042,.046,
	                       .051,.055,.059,.063,.067,.072,.076,.08,.085,.09,.095,.1,.108,.116,
	                       .124,.132,.14,.152,.164,.176,.188,.2,.22,.24,.26,.28,.3,.32,.34,
	                       .36,.38,.4,.417,.433,.450,.467,.483,.5,.51,.516,.52,.523,.5245,
	                       .525,.5242,
	                       .5214,.518,.512,.505,.495,.482,.469,.456,.442,.428,.414,.4,.386,
	                       .370,.355,.34,.325,.310,.295,.280,.265,.25,.236,.222,.208,.194,
	                       .180,.166,
	                       .152,.138,.124,.11,.101,.095,.09,.085,.08,.076,.072,.069,.066,
	                       .063,.06,.0575,.055,.0525,.05,.04875,.0475,.04625,.045,.04375,
	                       .0425,.04125,.04,.03875,.0375,.03625,.035,.03375,.0325,.03125,.03,
	                       .02925,.0285,.02775,.027,.02625,.0255,.02475,.024,.02325,.0225,
	                       .02175,.021,.02025,.0195,.01875,.018,.01725,.0165,.01575,.015,
	                       .01425,.0135,.01275,.012,.01125,.0105,.00975,.009,.00825,.0075,
	                       .00675,.006,.00525,.0045,.00375,.003,.00225,.0015,.00075,0.};



	double sen[161]={0.,0.,.012,.025,.038,.028,.028,.038,.035,.029,.039,.035,
	                 .038,.032,.038,.041,.041,.049,.055,.061,.072,.076,.070,.067,
	                 .080,.103,.125,.138,.118,.103,.129,.155,.170,.180,.190,.200,
	                 .215,.250,.302,.310,.301,.315,.330,.355,.380,.400,.410,.420,
	                 .438,.456,.474,.492,.510,.533,.556,.578,.6,.62,.63,.638,
	                 .640,.640,.637,.631,.625,.618,.610,.600,.580,.555,.530,.505,
	                 .480,.455,.435,.410,.385,.360,.340,.320,.300,.285,.270,.255,
	                 .240,.225,.210,.180,.165,.150,.140,.132,.124,.116,.108,.100,
	                 .092,.084,.077,.071,.066,.060,.055,.051,.048,.046,.044,.042,
	                 .040,.038,.036,.034,.032,.030,.028,.027,.026,.025,.025,.025,
	                 .024,.024,.024,.024,.024,.023,.023,.023,.023,.023,.022,.022,
	                 .022,.022,.022,.021,.021,.021,.020,.020,
	                 .020,.019,.018,.017,.016,.015,.014,.013,.012,.011,.010,.009,
	                 .008,.007,.006,.005,.004,.003,.002,.001,0.};

	// gammay,p gamma,n of Armstrong begin at 265 incr 25


	double sigt[160]={0.,.4245,.4870,.5269,.4778,.4066,.3341,.2444,.2245,.2005,
	                  .1783,.1769,.1869,.1940,.2117,.2226,.2327,.2395,.2646,.2790,.2756,
	                  .2607,.2447,.2211,.2063,.2137,.2088,.2017,.2050,.2015,.2121,.2175,
	                  .2152,.1917,.1911,.1747,.1650,.1587,.1622,.1496,.1486,.1438,.1556,
	                  .1468,.1536,.1544,.1536,.1468,.1535,.1442,.1515,.1559,.1541,.1461,
	                  .1388,.1565,.1502,.1503,.1454,.1389,.1445,.1425,.1415,.1424,.1432,
	                  .1486,.1539,.1354,.1480,.1443,.1435,.1491,.1435,.1380,.1317,.1445,
	                  .1375,.1449,.1359,.1383,.1390,.1361,.1286,.1359,.1395,.1327,.1387,
	                  .1431,.1403,.1404,.1389,.1410,.1304,.1363,.1241,.1284,.1299,.1325,
	                  .1343,.1387,.1328,.1444,.1334,.1362,.1302,.1338,.1339,.1304,.1314,
	                  .1287,.1404,.1383,.1292,.1436,.1280,.1326,.1321,.1268,.1278,.1243,
	                  .1239,.1271,.1213,.1338,.1287,.1343,.1231,.1317,.1214,.1370,.1232,
	                  .1301,.1348,.1294,.1278,.1227,.1218,.1198,.1193,.1342,.1323,.1248,
	                  .1220,.1139,.1271,.1224,.1347,.1249,.1163,.1362,.1236,.1462,.1356,
	                  .1198,.1419,.1324,.1288,.1336,.1335,.1266};


	double sigtn[160]={0.,.3125,.3930,.4401,.4582,.3774,.3329,.2996,.2715,.2165,
	                   .2297,.1861,.1551,.2020,.2073,.2064,.2193,.2275,.2384,.2150,.2494,
	                   .2133,.2023,.1969,.1797,.1693,.1642,.1463,.1280,.1555,.1489,.1435,
	                   .1398,.1573,.1479,.1493,.1417,.1403,.1258,.1354,.1394,.1420,.1364,
	                   .1325,.1455,.1326,.1397,.1286,.1260,.1314,.1378,.1353,.1264,.1471,
	                   .1650,.1311,.1261,.1348,.1277,.1518,.1297,.1452,.1453,.1598,.1323,
	                   .1234,.1212,.1333,.1434,.1380,.1330,.12,.12,.12,.12,.12,.12,.12,.12,
	                   .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
	                   .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
	                   .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
	                   .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,
	                   .12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12,.12};



	static int IFIRSTP=0;


	double si1=0, g1 =0,   o1=0;
	int   ne = 0, ij =0;
	double delo=0, omax =0, gk1m=0;
	static double scon=0., zcon=0.,o0=0.;


	double x=0,y=0,eps=0,eta=0,em=0,exx=0,s=0,ictr=0,pom=0,vec=0,gk1=0;

	//  maximum energy for GDR dissocation (in target frame, in MeV)

	double omax1n=24.01;

	if (IFIRSTP != 0) goto L100;

	IFIRSTP=1;


	//This is dependenant on gold or lead....Might need to expand
	if (zp == 79)
		{


			ap=197;
			si1=540.;
			g1=4.75;

			// peak and minimum energies for GDR excitation (in MeV)
			o1=13.70;
			o0=8.1;
		}
	else
		{
			zp=82;   //assumed to be lead
			ap=208;
			si1=640.;
			g1=4.05;
			o1=13.42;
			o0=7.4;
			for(int j=1;j<=160;j++)
				{

					sa[j]=sen[j];
				}
		}
	//Part II of initialization
	delo = .05;
	//.1 to turn mb into fm^2
	scon = .1*g1*g1*si1;
	zcon = zp/(gammatarg*( pi)*( 
	                                                hbarcmev))*zp/(gammatarg*( pi)*
		             ( hbarcmev))/137.04;//alpha?

	//single neutron from GDR, Veyssiere et al. Nucl. Phys. A159, 561 (1970)
	for ( int i = 1; i <= 160; i++) {
		eee[i] = o0+.1*(i-1);
		sa[i]  = 100.*sa[i];
	}
	//See Baltz, Rhoades-Brown, and Weneser, Phys. Rev. E 54, 4233 (1996) 
	//for details of the following photo cross-sections
	eee[161]=24.1;
	ne=int((25.-o0)/delo)+1;
	//GDR any number of neutrons, Veyssiere et al., Nucl. Phys. A159, 561 (1970)
	for ( int i = 1; i <= ne; i++ ) {
		ee[i] = o0+(i-1)*delo;
		//cout<<" ee 1 "<<ee[i]<<"  "<<i<<endl;

		se[i] = scon*ee[i]*ee[i]/(((o1*o1-ee[i]*ee[i])*(o1*o1-ee[i]*ee[i]))
		                          +ee[i]*ee[i]*g1*g1);
	}
	ij = ne;   //Risky?
	//25-103 MeV, Lepretre, et al., Nucl. Phys. A367, 237 (1981)
	for ( int j = 1; j <= 27; j++ ) {
		ij = ij+1;
		ee[ij] = e3[j];
		//cout<<" ee 2 "<<ee[ij]<<"  "<<ij<<endl;

		se[ij] = .1*ap*s3[j]/208.;
	}
	//103-440 MeV, Carlos, et al., Nucl. Phys. A431, 573 (1984)
	for ( int j = 1; j <= 22; j++ ) {
		ij = ij+1;
		ee[ij] = e1[j];
		//cout<<" ee 3 "<<ee[ij]<<"  "<<ij<<endl;
		se[ij] = .1*ap*s1[j]/208.;
	}
	//440 MeV-2 GeV Armstrong et al.
	for ( int j = 9; j <= 70; j++) {
		ij = ij+1;
		ee[ij] = ee[ij-1]+25.;
		//cout<<" ee 4 "<<ee[ij]<<"  "<<ij<<endl;
		se[ij] = .1*(zp*sigt[j]+(ap-zp)*sigtn[j]);
	}
	//2-16.4 GeV Michalowski; Caldwell
	for ( int j = 1; j <= 11; j++) {
		ij = ij+1;
		ee[ij] = e2[j];
		//cout<<" ee 5 "<<ee[ij]<<"   "<<ij<<endl;
		se[ij] = .1*ap*s2[j];
	}
	//Regge paramteres
	x = .0677;
	y = .129;
	eps = .0808;
	eta = .4525;
	em = .94;
	exx = pow(10,.05);

	//Regge model for high energy
	s = .002*em*ee[ij];
	//make sure we reach LHC energies
	ictr = 100;
	if ( gammatarg > (2.*150.*150.)) ictr = 150;
	for ( int j = 1; j <= ictr; j++ ) {
		ij = ij+1;
		s = s*exx;
		ee[ij] = 1000.*.5*(s-em*em)/em;
		//cout<<" ee 6 "<<ee[ij]<<"   "<<ij<<endl;
		pom = x*pow(s,eps);
		vec = y*pow(s,(-eta));
		se[ij] = .1*.65*ap*(pom+vec);
	}
	ee[ij+1] = 99999999999.;
	//done with initaliation
	//clear counters for 1N, XN
 L100:

	p1n = 0.;
	pxn = 0.;
	//start XN calculation
	//what's the b-dependent highest energy of interest?

	omax = min(omaxx,4.*gammatarg*( hbarcmev)/b);
	if ( omax < o0 ) return _pPhotonBreakup;
	gk1m = bessel::dbesk1(ee[1]*b/(( hbarcmev)*gammatarg));
	int k = 2;
 L212:
	if (ee[k] < omax ) {
		gk1 = bessel::dbesk1(ee[k]*b/(( hbarcmev)*gammatarg));
		//Eq. 3 of BCW--NIM in Physics Research A 417 (1998) pp1-8:
		pxn=pxn+zcon*(ee[k]-ee[k-1])*.5*(se[k-1]*ee[k-1]*gk1m*gk1m+se[k]*ee[k]*gk1*gk1);
		k = k + 1;
		gk1m = gk1;
		goto L212;
	}
	//one neutron dissociation
	omax = min(omax1n,4.*gammatarg*( hbarcmev)/b);
	gk1m = bessel::dbesk1(eee[1]*b/(( hbarcmev)*gammatarg));
	k = 2;
 L102:
	if (eee[k] < omax ) {
		gk1 = bessel::dbesk1(eee[k]*b/(( hbarcmev)*gammatarg));
		//Like Eq3 but with only the one neutron out GDR photo cross section input
		p1n = p1n+zcon*(eee[k]-eee[k-1])*.5*(sa[k-1]*eee[k-1]*gk1m*gk1m+sa[k]*eee[k]*gk1*gk1);
		k = k+1;
		gk1m = gk1;
		goto L102;
	}


	if (( mode) == 1) _pPhotonBreakup = 1.;
	if (( mode) == 2) _pPhotonBreakup = (1-exp(-1*pxn))*(1-exp(-1*pxn));
	if (( mode) == 3) _pPhotonBreakup = (p1n*exp(-1*pxn))*(p1n*exp(-1*pxn));
	if (( mode) == 4) _pPhotonBreakup = exp(-2*pxn);
	if (( mode) == 5) _pPhotonBreakup = 1.;
	if (( mode) == 6) _pPhotonBreakup = (1. - exp(-2.*pxn));
	if (( mode) == 7) _pPhotonBreakup = 2.*exp(-pxn)*(1.-exp(-pxn));

	//cout<<pxn<<" "<<zcon<<" "<<ee[k]<<" "<<se[k-1]<<" "<<gk1m<<"  "<<gk1<<"  "<<k<<"  "<<ee[k+1]<< "  "<<b<< endl;

	return _pPhotonBreakup;
}
