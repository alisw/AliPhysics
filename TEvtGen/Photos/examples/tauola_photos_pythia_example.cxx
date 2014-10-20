/**
 * Example of use of photos C++ interface. Pythia events are
 * generated first and photos used for FSR.
 *
 * @author Nadia Davidson
 * @date 6 July 2009
 */

//pythia header files
#include "Pythia.h"
#include "HepMCInterface.h"

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

//TAUOLA header files
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

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
	HepMC::I_Pythia8 ToHepMC;

	// Initialization of pythia
	Pythia pythia;
	Event& event = pythia.event;

	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("PartonLevel:FSR = off");

	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 15");
	pythia.particleData.readString("15:mayDecay = off"); //<- uncomment for pythia+tauola

	//pythia.init( -2212, -2212, 14000.0);     //proton proton collisions
	pythia.init( 11, -11, 91.187);             //electron positron collisions

	// TAUOLA and PHOTOS initialization
	Tauola::initialize();
	Photos::initialize();

	Photos::setInfraredCutOff(0.01/91.187); // 10MeV for scale to M_Z=91.187
	//Photos::setDoubleBrem(false);
	//Photos::setExponentiation(false);

	Log::SummaryAtExit();
	cout.setf(ios::fixed);
	//Log::LogInfo(false) //To turn printing of last five events and pythia statistics off

	// Example setup - suppress processing of whole Z0 decay,
	// leaving only the Z0 -> tau+ tau- decay and whole branch starting
	// from tau- to be processed
	//Photos::suppressBremForBranch(0,23);
	//Photos::forceBremForDecay (2,23,15,-15);
	//Photos::forceBremForBranch(0,15);

	// Force mass of electron and positron to be 0.000511
	//Photos::forceMass(11,0.000511);

	// Force mass of electron and positron to be taken
	// from event record instead of being calculated from 4-vector
	//Photos::forceMassFromEventRecord(11);

	// Exclude particles with given status code from being processed
	// or taken into account during momentum conservation calculation
	//Photos::ignoreParticlesWithStatus(3);

	// Remove status code from the list of ignored status codes
	//Photos::DeIgnoreParticlesWithStatus(3);

	// Force writing history decay products for vertices
	// modified i.e. with added photons. These particles will
	// have the provided status code. Photos will ignore
	// all particles with this status code.
	//Photos::createHistoryEntries(true,3);

	MC_Initialize();

	// Begin event loop
	for (unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if(!pythia.next()) continue;

		// Convert event record to HepMC
		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);

		// Run TAUOLA on the event
		TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

		// We may want to undecay previously decayed taus.
		//t_event->undecayTaus();
		t_event->decayTaus();
		delete t_event;

		//Log::LogPhlupa(2,4);

		if(iEvent<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation chceck BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

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

		// Print out last 5 events
		if(iEvent>=NumberOfEvents-5)
		{
			Log::RedirectOutput(Log::Info());
			//pythia.event.list();
			HepMCEvt->print();
			Log::RevertOutput();
		}

		// Clean up
		delete HepMCEvt;
	}

	Log::RedirectOutput(Log::Info());
	pythia.statistics();
	Log::RevertOutput();

	MC_Finalize();
}
