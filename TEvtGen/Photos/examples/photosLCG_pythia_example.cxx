/**
 * Example of use of photos C++ interface. Pythia events are
 * generated first and photos used for FSR.
 *
 * @author LCG & Tomasz Przedzinski
 * @date 6 May 2011
 */

//pythia header files
#include "Pythia.h"
#include "HepMCInterface.h"

//PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Pythia8;
using namespace Photospp;

bool ShowersOn=true;
unsigned long NumberOfEvents = 100000;

// Calculates energy ratio between (l + bar_l) and (l + bar_l + X)
double calculate_ratio(HepMC::GenEvent *evt, double *ratio_2)
{
	double ratio = 0.0;
	for(HepMC::GenEvent::particle_const_iterator p=evt->particles_begin();p!=evt->particles_end(); p++)
	{
		// Search for Z
		if( (*p)->pdg_id() != 23 ) continue;

		// Ignore this Z if it does not have at least two daughters
		if( !(*p)->end_vertex() || (*p)->end_vertex()->particles_out_size() < 2 ) continue;

		// Sum all daughters other than photons
		double sum = 0.0;
		for(HepMC::GenVertex::particle_iterator pp = (*p)->end_vertex()->particles_begin(HepMC::children);
		    pp != (*p)->end_vertex()->particles_end(HepMC::children);
		    ++pp)
		{
		   // (*pp)->print();
		   if( (*pp)->pdg_id() != 22 ) sum += (*pp)->momentum().e();
		}

		// Calculate ratio and ratio^2
		ratio    = sum     /   (*p)->momentum().e();
		*ratio_2 = sum*sum / ( (*p)->momentum().e() * (*p)->momentum().e() );

		// Assuming only one Z decay per event
		return ratio;
	}
	return 0.0;
}

int main(int argc,char **argv)
{
	// Initialization of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	//pythia.settings.listAll();

	if(!ShowersOn)
	{
		//pythia.readString("HadronLevel:all = off");
		pythia.readString("HadronLevel:Hadronize = off");
		pythia.readString("SpaceShower:QEDshower = off");
		pythia.readString("SpaceShower:QEDshowerByL = off");
		pythia.readString("SpaceShower:QEDshowerByQ = off");
	}
	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:FSR = off");

	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 13");
	//pythia.readString("23:onIfAny = 11");
	pythia.init( 11, -11, 91.187);                           //e+ e- collisions

	Photos::initialize();

	// Turn on NLO corrections - only for PHOTOS 3.2 or higher

	//Photos::setMeCorrectionWtForZ(zNLO);
	//Photos::maxWtInterference(4.0);
	//Photos::iniInfo();

	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	double ratio_exp   = 0.0, ratio_alpha   = 0.0;
	double ratio_exp_2 = 0.0, ratio_alpha_2 = 0.0;
	double buf = 0.0;

	Photos::setDoubleBrem(true);
	Photos::setExponentiation(true);
	Photos::setInfraredCutOff(0.000001);

	Log::Info() << "---------------- First run - EXP ----------------" <<endl;

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%10000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if (!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

		//Log::LogPhlupa(1,3);

		// Call photos
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		ratio_exp   += calculate_ratio(HepMCEvt,&buf);
		ratio_exp_2 += buf;

		//HepMCEvt->print();

		// Clean up
		delete HepMCEvt;
	}

	Photos::setDoubleBrem(false);
	Photos::setExponentiation(false);
	Photos::setInfraredCutOff(1./91.187);  // that is 1 GeV for
	                                       // pythia.init( 11, -11, 91.187);

	Log::Info() << "---------------- Second run - ALPHA ORDER ----------------" <<endl;
	
	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%10000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if (!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

		//Log::LogPhlupa(1,3);

		// Call photos
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		ratio_alpha   += calculate_ratio(HepMCEvt,&buf);
		ratio_alpha_2 += buf;

		//HepMCEvt->print();

		// Clean up
		delete HepMCEvt;
	}
	
	pythia.statistics();

	ratio_exp   = ratio_exp   / (double)NumberOfEvents;
	ratio_exp_2 = ratio_exp_2 / (double)NumberOfEvents;

	ratio_alpha   = ratio_alpha   / (double)NumberOfEvents;
	ratio_alpha_2 = ratio_alpha_2 / (double)NumberOfEvents;

	double err_exp   = sqrt( (ratio_exp_2   - ratio_exp   * ratio_exp  ) / (double)NumberOfEvents );
	double err_alpha = sqrt( (ratio_alpha_2 - ratio_alpha * ratio_alpha) / (double)NumberOfEvents );
	
	cout.precision(6);
	cout.setf(ios::fixed);
	cout << "********************************************************" << endl;
	cout << "* Z -> l + bar_l + gammas                              *" << endl;
	cout << "* E(l + bar_l) / E(l + bar_l + gammas) ratio           *" << endl;
	cout << "*                                                      *" << endl;
	cout << "* PHOTOS - EXP:          " << ratio_exp   <<" +/- "<<err_exp  <<"      *" << endl;
	cout << "* PHOTOS - ALPHA ORDER:  " << ratio_alpha <<" +/- "<<err_alpha<<"      *" << endl;
	cout << "********************************************************" << endl;

}
