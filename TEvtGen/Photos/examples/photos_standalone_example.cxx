/**
 * Example of photos usage.
 * Events are loaded from pre-generated set featuring Z0 -> tau+ tau- decays
 * and processed by photos.
 *
 * @author Tomasz Przedzinski
 * @date 17 July 2010
 */

//HepMC header files
#include "HepMC/IO_GenEvent.h"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Photospp;

int EventsToCheck=20;

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

int main()
{
	HepMC::IO_GenEvent file("photos_standalone_example.dat",std::ios::in);

	Photos::initialize();
	Photos::setInfraredCutOff(0.001/200);

	int photonAdded=0,twoAdded=0,moreAdded=0,evtCount=0;
	// Begin event loop. Generate event.
	while(true)
	{
		// Create event
		HepMC::GenEvent *HepMCEvt = new HepMC::GenEvent();
		file.fill_next_event(HepMCEvt);
		if(file.rdstate()) break;
		evtCount++;
		int buf = -HepMCEvt->particles_size();

		//cout << "BEFORE:"<<endl;
		//HepMCEvt->print();

		if(evtCount<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation chceck BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

		// Process by photos
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		if(evtCount<EventsToCheck)
		{
			checkMomentumConservationInEvent(HepMCEvt);
		}

		buf+=HepMCEvt->particles_size();
		if(buf==1)      photonAdded++;
		else if(buf==2) twoAdded++;
		else if(buf>2)  moreAdded++;

		//cout << "AFTER:"<<endl;
		//HepMCEvt->print();

		//clean up
		delete HepMCEvt;
	}

	// Print results
	cout.precision(2);
	cout.setf(ios::fixed);
	cout<<endl;
	if(evtCount==0)
	{
		cout<<"Something went wrong with the input file: photos_standalone_example.dat"<<endl;
		cout<<"No events were processed."<<endl<<endl;
		return 0;
	}
	cout<<"Summary (whole event processing):"<<endl;
	cout<<evtCount   <<"\tevents processed"<<endl;
	cout<<photonAdded<<"\ttimes one photon added to the event           \t("<<(photonAdded*100./evtCount)<<"%)"<<endl;
	cout<<twoAdded   <<"\ttimes two photons added to the event          \t("<<(twoAdded*100./evtCount)<<"%)"<<endl;
	cout<<moreAdded  <<"\ttimes more than two photons added to the event\t("<<(moreAdded*100./evtCount)<<"%)"<<endl<<endl;
	cout<<"(Contrary to results from MC-Tester, these values are technical and infrared unstable)"<<endl<<endl;
}
