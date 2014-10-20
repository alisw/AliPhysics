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

using namespace std;
using namespace Pythia8;
using namespace Photospp;

bool ShowersOn=true;
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

/* Switch Status of History Entries

   If Photos::createHistoryEntries(true,3) was called, this function changes the 
   status code of photons added by Photos and particles modified by Photos
   to 3, switching the status of history entries to 1.
   
   This results leaves all modifications performed by Photos as history entries,
   while the regular entries represent original, unmodified event.
   
   This is an example of how such operation can be performed in user analysis.
   By default, this function is not used. The example of its use is commented
   out in main event loop.
   
   NOTE: The algorithm works only on stable particles and assumes that
         there were no modifications to the order of the particles in
         which they were written to HepMC by Photos. */
void switch_history_entries_status(HepMC::GenEvent *evt)
{
  for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
    if((*p)->status()==3)
    {
      if((*p)->pdg_id()==22) continue;

      int barcode = (*p)->barcode();

      HepMC::GenVertex *v = (*p)->production_vertex();
      
      // History entries are added after photons, so we check what is the
      // position of current particle relative to photons.
      int position = 0;
      int last_photon_position = -1;
      
      for(HepMC::GenVertex::particles_out_const_iterator p2 = v->particles_out_const_begin();
          p2 != v->particles_out_const_end(); ++p2)
      {
        position++;

        if((*p2)->barcode()==barcode) break;

        if((*p2)->pdg_id()==22) { last_photon_position=position; }
      }
      
      // If particle is found prior to photons - it was already processed, so skip it
      if(last_photon_position<0) continue;
      
      position -= last_photon_position;
      HepMC::GenParticle *part = NULL;
      
      // Now, find the particle that corresponds to this history entry
      for(HepMC::GenVertex::particles_out_const_iterator p2 = v->particles_out_const_begin();
          p2 != v->particles_out_const_end(); ++p2)
      {
        --position;
        
        if     (position >  0) continue;
        else if(position == 0) part = *p2;
        else
        {
          // Give all remaining photons status 3
          if((*p2)->pdg_id()==22 ) (*p2)->set_status(3);

          // Finish if there are no more photons
          else break;
        }
      }

      // Check if this is the particle we are looking for
      if( part->pdg_id() != (*p)->pdg_id())
      {
        cout<<"switch_history_entries_status: mismatch in pdg_id of history entry"<<endl;
        cout<<"and its corresponding particle. The algorithm does not work correctly."<<endl;
        exit(-1);
      }

      // Skip this particle if its status is not 1
      if(part->status()!=1) continue;
      
      // Switch status codes of these particles
      part->set_status(3);
      (*p)->set_status(1);
    }
  }
}

int main(int argc,char **argv)
{
	// Initialization of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	//pythia.settings.listAll();

	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:FSR = off");

	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 13");
	pythia.init( 11, -11, 91.187);                           //e+ e- collisions

	MC_Initialize();

	Photos::initialize();
	//Photos::setDoubleBrem(false);
	//Photos::setExponentiation(false);

	Photos::setInfraredCutOff(0.01/91.187); // 10MeV for scale to M_Z=91.187
	Photos::maxWtInterference(3.0);
  //Photos::createHistoryEntries(true,3);

	Photos::iniInfo();
	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if (!pythia.next()) continue;

		// Convert event record to HepMC
		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

		if(iEvent<EventsToCheck)
		{
			cout<<"                                          "<<endl;
			cout<<"Momentum conservation chceck BEFORE/AFTER Photos"<<endl;
			checkMomentumConservationInEvent(HepMCEvt);
		}

		//Log::LogPhlupa(1,3);

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

    // Uncomment to turn on switching of the status code of history entries
    // with the regular entries for stable particles
    //switch_history_entries_status(HepMCEvt);

		if(iEvent<EventsToCheck)
		{
			checkMomentumConservationInEvent(HepMCEvt);
		}

		//HepMCEvt->print();

		// Run MC-TESTER on the event
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		// Print out last 5 events
		if(iEvent>=NumberOfEvents-5) HepMCEvt->print();

		// Clean up
		delete HepMCEvt;
	}
	pythia.statistics();
	MC_Finalize();
}
