/**
 * Main program for testing photos C++ interface.
 * Pythia events are generated first, Tauola++ used for tau decays
 * and photos used for FSR.
 *
 * @author Nadia Davidson and Tomasz Przedzinski
 * @date 10 May 2011
 */

//Pythia header files
#include "Pythia.h"
#include "HepMCInterface.h"

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

//TAUOLA header files
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCParticle.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Pythia8;
using namespace Photospp;
using namespace Tauolapp;

unsigned long NumberOfEvents = 10000;
unsigned int EventsToCheck=20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Photos (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
	//cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;
	
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
		if( (*p)->status() == 1 )
		{
			HepMC::FourVector m = (*p)->momentum();
			px+=m.px();
			py+=m.py();
			pz+=m.pz();
			e +=m.e();
			//(*p)->print();
		}
	}
  cout.precision(6);
  cout.setf(ios_base::floatfield);
	cout<<endl<<"Vector Sum: "<<px<<" "<<py<<" "<<pz<<" "<<e<<endl;
}

int main(int argc,char **argv)
{

	// Program needs at least 4 parameters
	if(argc<5)
	{
		cout<<endl<<"Usage: "<<argv[0]<<" <pythia_conf> <pythia_mode> <no_events> <tauola_mode> [ <alpha_order> <ScalarNLO_mode> ]"<<endl;
		cout<<endl<<"   eg. "<<argv[0]<<" pythia_H.conf 0 10000 4 0 0"<<endl;
		cout<<endl;
		return -1;
	}

	HepMC::I_Pythia8 ToHepMC;

	// Initialization of pythia
	Pythia pythia;
	Event& event = pythia.event;

	pythia.readString("HadronLevel:Hadronize = off");
	pythia.readString("SpaceShower:QEDshower = off");
	pythia.readString("SpaceShower:QEDshowerByL = off");
	pythia.readString("SpaceShower:QEDshowerByQ = off");
	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("PartonLevel:FSR = off");

	// Tauola is currently set to undecay taus. Otherwise, uncomment this line.
	//pythia.particleData.readString("15:mayDecay = off");

	/********************************************************
	  Read input parameters from console. List of parameters:
	  1. Pythia configuration filename
	  2. Pythia mode - e+e-@200GeV , e+e-@91.187GeV, pp@14TeV or e+e-@500GeV
	     (only e+e-@91.187GeV and e+e-@500GeV are used in this example)
	  3. Number of events
	  4. Tauola decay mode (refer to Tauola documentation)
	  5. Photos - use 1-photon mode on/off
	  6. Photos - use ScalarNLO mode on/off

	  Example where all input parameters are used:

	  ./photos_tauola_test.exe pythia_H.conf 1 100000 4 0 0
	  - use pythia_H.conf
	  - initialize using e+ e- @ 91.187GeV collisions
	  - generate 100 000 events
	  - fix TAUOLA decay to channel 4 (RHORHO_MODE)
	  - Photos is not using 1-photon mode (default option)
	  - Photos is not in ScalarNLO mode (default option)
	*********************************************************/

	// 1. Load pythia configuration file (argv[1], from console)
	if(argc>1) pythia.readFile(argv[1]);

	// 2. Initialize pythia to e+e-@91.17GeV or e+e-@500GeV collisions (argv[2], from console)
	if(atoi(argv[2])==1)      pythia.init( 11, -11, 91.187); // e+ e- collisions
	else if(atoi(argv[2])==3) pythia.init( 11, -11, 500.);   // e+ e- collisions
	else
	{
		cout<<"ERROR: Wrong Pythia mode ("<<atoi(argv[4])<<")"<<endl;
		cout<<"       Only modes '1' and '3' are used by this program."<<endl;
		return -1;
	}

	// 3. Get number of events (argv[3], from console)
	if(argc>3) NumberOfEvents=atoi(argv[3]);

	// 4. Set Tauola decay mode (argv[4], from console)
	if(argc>4)
	{
		// argv[4]=3 (tau => pi nu_tau)    for Ztautau
		// argv[4]=4 (tau => pi pi nu_tau) for Htautau
		Tauola::setSameParticleDecayMode(atoi(argv[4]));
		Tauola::setOppositeParticleDecayMode(atoi(argv[4]));
	}

	Tauola::initialize();
	Photos::initialize();

	Photos::setExponentiation(true);
	Photos::setInfraredCutOff(1.e-6);
	Photos::maxWtInterference(3.0);

	// 5. Check if we're using 1-photon mode
	if( argc>5 && atoi(argv[5]) )
	{
		Photos::setDoubleBrem(false);
		Photos::setExponentiation(false);

		// Set infrared cutoff to 10MeV for scale M_Z=91.187GeV or 500 GeV
		if(atoi(argv[2])==1) Photos::setInfraredCutOff(0.01/91.187);
		else                 Photos::setInfraredCutOff(0.01/500.);

		Photos::maxWtInterference(2.0);
	}

	// 6. Check if we're in ScalarNLO mode
	if( argc>6 )
	{
		Tauola::setEtaK0sPi(1,1,0);
    
		// Check if we are using NLO
		if(atoi(argv[6])) Photos::setMeCorrectionWtForScalar(true);
	}

	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	MC_Initialize();

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if(!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);

		if(iEvent<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation chceck BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Run TAUOLA on the event
		TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

		// Since we let Pythia decay taus, we have to undecay them first.
		t_event->undecayTaus();
		t_event->decayTaus();
		delete t_event;

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		if(iEvent<EventsToCheck)
		{
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Run MC-TESTER on the event
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		//clean up
		delete HepMCEvt;
	}
	pythia.statistics();
	MC_Finalize();
}
