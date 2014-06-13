#include <stdarg.h>
#include <iostream>
#include <vector>

#include "PhotosRandom.h"
#include "PhotosEvent.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
using std::cout;
using std::endl;
using std::ios_base;

namespace Photospp
{

Photos Photos::_instance;

vector<vector<int>* >    *Photos::supBremList    = 0;
vector<vector<int>* >    *Photos::forceBremList  = 0;
vector<pair<int,double>* > *Photos::forceMassList = 0;
vector<int>              *Photos::ignoreStatusCodeList = 0;
bool Photos::isSuppressed=false;
bool Photos::massFrom4Vector=true;
double Photos::momentum_conservation_threshold   = 0.1;
bool Photos::meCorrectionWtForZ=false;
bool Photos::meCorrectionWtForW=false;
bool Photos::meCorrectionWtForScalar=false;
bool Photos::isCreateHistoryEntries=false;
int  Photos::historyEntriesStatus = 3;
double (*Photos::randomDouble)() = PhotosRandom::randomReal;

Photos::Photos()
{
	setAlphaQED           (0.00729735039);
	setInfraredCutOff     (0.01);
	setInterference       (true);
	setDoubleBrem         (true);
	setQuatroBrem         (false);
	setTopProcessRadiation(true);
	setCorrectionWtForW   (true);

	// setExponentiation(true) moved to initialize() due to communication
	// problems with Fortran under MacOS.
	phokey_.iexp = 1;
}

void Photos::initialize()
{
// Should return if already initialized?

/*******************************************************************************
  All the following parameter setters can be called after PHOINI.
  Initialization of kinematic correction against rounding errors.
  The set values will be used later if called with zero.
  Default parameter is 1 (no correction) optionally 2, 3, 4
  In case of exponentiation new version 5 is needed in most cases.
  Definition given here will be thus overwritten in such a case
  below in routine PHOCIN
*******************************************************************************/
	if(!phokey_.iexp) initializeKinematicCorrections(1);
	else              setExponentiation(true);

// Initialize status counter for warning messages
	for(int i=0;i<10;i++) phosta_.status[i]=0;
// elementary security level, should remain 1 but we may want to have a method to change.
	phosta_.ifstop=1; 

	pholun_.phlun=6; // Logical output unit for printing error messages

// Further initialization done automatically
// see places with - VARIANT A - VARIANT B - all over to switch between options

#ifndef VARIANTB
//----------- SLOWER VARIANT A, but stable ------------------
//--- it is limiting choice for small XPHCUT in fixed orer
//--- modes of operation

// Best choice is if FINT=2**N where N+1 is maximal number
// of charged daughters
// see report on overweihted events
	if(phokey_.interf) maxWtInterference(2.0);
	else               maxWtInterference(1.0);
#else

//----------- FASTER VARIANT B  ------------------
//-- it is good for tests of fixed order and small XPHCUT
//-- but is less promising for more complex cases of interference
//-- sometimes fails because of that

	if(phokey_.interf) maxWtInterference(1.8);
	else               maxWtInterference(0.0);
#endif
//------------END VARIANTS A B -----------------------

//------------------------------------------------------------------------------
// Print PHOTOS header
//------------------------------------------------------------------------------
	int                coutPrec = cout.precision(6);
	ios_base::fmtflags flags    = cout.setf(ios_base::floatfield);
	cout<<endl;
	cout<<"********************************************************************************"<<endl<<endl;

	cout<<"                            ========================="<<endl;
	cout<<"                              PHOTOS, Version:  "<<VER_MAJOR<<"."<<VER_MINOR<<endl;
	cout<<"                              Released at:  "<<DAT_DAY<<"/"<<DAT_MONTH<<"/"<<DAT_YEAR<<endl;
	cout<<"                            ========================="<<endl<<endl;

	cout<<"                     Photos QED corrections in Particle Decays"<<endl<<endl;

	cout<<"           Monte Carlo Program - by E. Barberio, B. van Eijk and Z. Was"<<endl;
	cout<<"           From version 2.09   - by P. Golonka and Z. Was"<<endl;
	cout<<"           From version 3.00   - by N. Davidson, T. Przedzinski and Z. Was"<<endl;

	cout<<"********************************************************************************"<<endl<<endl;

	cout<<"                  Internal (default) input parameters: "<<endl<<endl;
	cout<<"                    INTERF= "<<phokey_.interf<<" ISEC= " <<phokey_.isec <<" ITRE= "<<phokey_.itre
	                     <<" IEXP= "  <<phokey_.iexp  <<" IFTOP= "<<phokey_.iftop<<" IFW= " <<phokey_.ifw <<endl;
	cout<<"                    ALPHA_QED= "<<phocop_.alpha<<" XPHCUT= "<<phocop_.xphcut<<endl<<endl;

	if(phokey_.interf) cout<<"                    Option with interference is active"<<endl;
	if(phokey_.isec)   cout<<"                    Option with double photons is active"<<endl;
	if(phokey_.itre)   cout<<"                    Option with triple/quatric photons is active"<<endl;
	if(phokey_.iexp)   cout<<"                    Option with exponentiation is active EPSEXP="<<phokey_.expeps<<endl;
	if(phokey_.iftop)  cout<<"                    Emision in t tbar production is active"<<endl;
	if(phokey_.ifw)    cout<<"                    Correction wt in decay of W is active"<<endl;
	if(meCorrectionWtForZ)    cout<<"                    ME in decay of Z is active"<<endl;
	if(meCorrectionWtForW)    cout<<"                    ME in decay of W is active"<<endl;
	cout<<endl<<"          WARNING:  /HEPEVT/ is not anymore used."<<endl<<endl;
/*
	cout<<endl<<"            WARNING (1): /HEPEVT/ is not anymore the standard common block"<<endl<<endl;

	cout<<"            PHOTOS expects /HEPEVT/ to have REAL*8 variables. To change to"<<endl;
	cout<<"            REAL*4 modify its declaration in subr. PHOTOS_GET PHOTOS_SET:"<<endl;
	cout<<"                 REAL*8  d_h_phep,  d_h_vhep"<<endl;
	cout<<"            WARNING (2): check dims. of /hepevt/ /phoqed/ /ph_hepevt/."<<endl;
	cout<<"            HERE:                     d_h_nmxhep=10000 and  NMXHEP=10000"<<endl<<endl;
*/
	cout<<"********************************************************************************"<<endl;
	// Revert output stream flags and precision
	cout.precision(coutPrec);
	cout.flags    (flags);

	// Add reggeon, pomeron and its diffractive states to the list
	// of decays where bremsstrahlung is suppressed
	Photos::suppressBremForDecay(0,110);
	Photos::suppressBremForDecay(0,990);
	Photos::suppressBremForDecay(0,9902110);
	Photos::suppressBremForDecay(0,9902210);
	Photos::suppressBremForDecay(0,9900110);
	Photos::suppressBremForDecay(0,9900210);
	Photos::suppressBremForDecay(0,9900220);
	Photos::suppressBremForDecay(0,9900330);
	Photos::suppressBremForDecay(0,9900440);

  // Set suppression of all pi0 decays and K_L -> gamma e+ e- ...
  // Previously initialization in Fortran IPHEKL(i) routine and used in PHOCHK 
  // i=1 was emission allowed, i=2 was blocked 0 was when the option was used.
  // now in IPHEKL_setPi0KLnoEmmision we have only 1 to allow emissions 
  // and 2 to block.
  // Can be overriden by using 'Photos::IPHEKL_setPi0KLnoEmmision(0)'
  // method several times use Photos::forceBremForDecay() and can be 
  // over-ruled in part. 

  Photos::IPHEKL_setPi0KLnoEmission(2); // Blocks emission in  pi0  and in Kl to gamma e+ e- if parameter is !1 (enables otherwise)
  Photos::IPHQRK_setQarknoEmission (1,0);// Blocks emission from quarks if buf=1 (default); enables if buf=2
	                                 //  option 2 is not realistic and for tests only

// Initialize Marsaglia and Zaman random number generator
	PhotosRandom::initialize();
}
void Photos::iniInfo()
{
// prints infomation like Photos::initialize; may be called after reinitializations.

/*******************************************************************************
  Once  parameter setters are called after PHOINI one may want to know and write
  into output info including all reinitializations.
*******************************************************************************/


//------------------------------------------------------------------------------
// Print PHOTOS header again
//------------------------------------------------------------------------------
	int                coutPrec = cout.precision(6);
	ios_base::fmtflags flags    = cout.setf(ios_base::floatfield);
	cout<<endl;
	cout<<"********************************************************************************"<<endl<<endl;
	cout<<"                            ========================================="<<endl;
	cout<<"                            PHOTOS, information routine"<<endl;
	cout<<"                            Input parameters after reinitialization: "<<endl<<endl;
	cout<<"                            ========================================="<<endl<<endl;
	cout<<"********************************************************************************"<<endl<<endl;
	cout<<"                    INTERF= "<<phokey_.interf<<" ISEC= " <<phokey_.isec <<" ITRE= "<<phokey_.itre
	                     <<" IEXP= "  <<phokey_.iexp  <<" IFTOP= "<<phokey_.iftop<<" IFW= " <<phokey_.ifw <<endl;
	cout<<"                    ALPHA_QED= "<<phocop_.alpha<<" XPHCUT= "<<phocop_.xphcut<<endl<<endl;

	if(phokey_.interf) cout<<"                    Option with interference is active"<<endl;
	if(phokey_.isec)   cout<<"                    Option with double photons is active"<<endl;
	if(phokey_.itre)   cout<<"                    Option with triple/quatric photons is active"<<endl;
	if(phokey_.iexp)   cout<<"                    Option with exponentiation is active EPSEXP="<<phokey_.expeps<<endl;
	if(phokey_.iftop)  cout<<"                    Emision in t tbar production is active"<<endl;
	if(phokey_.ifw)    cout<<"                    Correction wt in decay of W is active"<<endl;
	if(meCorrectionWtForZ)    cout<<"                    ME in decay of Z is active"<<endl;
	if(meCorrectionWtForW)    cout<<"                    ME in decay of W is active"<<endl;
	if(meCorrectionWtForScalar)    cout<<"                    ME in decay of Scalar is active"<<endl;

	cout<<endl<<"          WARNING:  /HEPEVT/ is not anymore used."<<endl<<endl;
	// Revert output stream flags and precision
	cout.precision(coutPrec);
	cout.flags    (flags);
}

void Photos::processParticle(PhotosParticle *p)
{
	PhotosBranch b(p);
	if(!b.getSuppressionStatus()) b.process();
}

void Photos::processBranch(PhotosParticle *p)
{
	vector<PhotosParticle *> particles = p->getDecayTree();
	vector<PhotosBranch *>   branches = PhotosBranch::createBranches(particles);
	for(int i=0;i<(int)branches.size();i++) branches.at(i)->process();
}

void Photos::suppressBremForDecay(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(0);
	if(!supBremList) supBremList = new vector< vector<int>* >();
	supBremList->push_back(v);
}

void Photos::suppressBremForBranch(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(1);
	if(!supBremList) supBremList = new vector< vector<int>* >();
	supBremList->push_back(v);
}

void Photos::forceBremForDecay(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(0);
	if(!forceBremList) forceBremList = new vector< vector<int>* >();
	forceBremList->push_back(v);
}

void Photos::forceBremForBranch(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(1);
	if(!forceBremList) forceBremList = new vector< vector<int>* >();
	forceBremList->push_back(v);
}

  // Previously this functionality was encoded in FORTRAN routine
  // PHOCHK which was having some other functionality as well
void Photos::IPHEKL_setPi0KLnoEmission(int m)
{
  if(m==1)
  {
    cout << "MODOP=1 -- enables emission in pi0 to gamma e+e- : TEST " << endl ;
    cout << "MODOP=1 -- enables emission in Kl  to gamma e+e- : TEST " << endl ;
    Photos::forceBremForDecay(0,111);
    Photos::forceBremForDecay(3, 130,22,11,-11);
    Photos::forceBremForDecay(3,-130,22,11,-11);
  }
  else if(m!=1)
  {
    cout << "MODOP=2 -- blocks emission in Kl  to gamma e+e-: DEFAULT" << endl ;
    cout << "MODOP=2 -- blocks emission in pi0 to gamma e+e-: DEFAULT" << endl ;
    Photos::suppressBremForDecay(0,111);
    Photos::suppressBremForDecay(3, 130,22,11,-11);
    Photos::suppressBremForDecay(3,-130,22,11,-11);
  }
}

  // Previously this functionality was encoded in FORTRAN routine
  // PHOCHK which was having some other functionality as well
bool Photos::IPHQRK_setQarknoEmission(int MODCOR, int PDGID)
{
  static int IPHQRK_MODOP=-1;
  if(IPHQRK_MODOP==-1 && MODCOR==0){
    cout << "stop from IPHQRK_setQarknoEmission lack of initialization" << endl ;
    exit(-1);
  }
  else if (MODCOR != 0){
    IPHQRK_MODOP = MODCOR;
    if(MODCOR ==1) cout << " IPHQRK_setQarknoEmission MODOP=1 -- blocks emission from light quarks:  DEFAULT" << endl ;
    if(MODCOR !=1) cout << " IPHQRK_setQarknoEmission MODOP=2 -- emission from light quarks allowed: TEST   " << endl ;
  }
  if(IPHQRK_MODOP!=1) return true;
  
  return PDGID>9;
}

void Photos::createHistoryEntries(bool flag, int status)
{
  if(status<3)
  {
    Log::Warning()<<"Photos::createHistoryEntries: status must be >=3"<<endl;
    return;
  }

  isCreateHistoryEntries = flag;
  historyEntriesStatus   = status;
  ignoreParticlesOfStatus(status);
}

void Photos::ignoreParticlesOfStatus(int status)
{
  if(status<3)
  {
    Log::Warning()<<"Photos::ignoreParticlesOfStatus: status must be >=3"<<endl;
    return;
  }
  
  if(!ignoreStatusCodeList) ignoreStatusCodeList = new vector<int>();

  // Do not add duplicate entries to the list
  for(unsigned int i=0;i<ignoreStatusCodeList->size();i++)
    if( status==ignoreStatusCodeList->at(i) ) return;
  
  ignoreStatusCodeList->push_back(status);
}

void Photos::deIgnoreParticlesOfStatus(int status)
{
  if(!ignoreStatusCodeList) return;

  for(unsigned int i=0;i<ignoreStatusCodeList->size();i++)
  {
    if( status==ignoreStatusCodeList->at(i) )
    {
      ignoreStatusCodeList->erase(ignoreStatusCodeList->begin()+i);
      return;
    }
  }
}

bool Photos::isStatusCodeIgnored(int status)
{
  if(!ignoreStatusCodeList) return false;

  for(unsigned int i=0;i<ignoreStatusCodeList->size();i++)
    if( status==ignoreStatusCodeList->at(i) ) return true;

  return false;
}

void Photos::setRandomGenerator( double (*gen)() )
{
  if(gen==NULL) randomDouble = PhotosRandom::randomReal;
  else          randomDouble = gen;
}

void Photos::setExponentiation(bool expo)
{
	phokey_.iexp = (int) expo;
	if(expo)
	{
		setDoubleBrem(false);
		setQuatroBrem(false);
		setInfraredCutOff(0.0000001);
		initializeKinematicCorrections(5);
		phokey_.expeps=0.0001;
	}
}

void Photos::setMeCorrectionWtForW(bool corr)
{
	meCorrectionWtForW=corr;
}

void Photos::setMeCorrectionWtForZ(bool corr)
{
	meCorrectionWtForZ=corr;
}
void Photos::setMeCorrectionWtForScalar(bool corr)
{
	meCorrectionWtForScalar=corr;
}

void Photos::setStopAtCriticalError(bool stop)
{
	phosta_.ifstop=(int)stop;
	if(!stop)
	{
		Log::Info()<<"PHOTOS production mode. Elementary test of data flow from event record disabled. "<<endl
		           <<"Prior checks of the complete configuration "<<endl
		           <<"(for the particular set of input parameters) must have been done! "<<endl;
	}
}


void Photos::forceMassFromEventRecord(int pdgid)
{
  if(!forceMassList) forceMassList = new vector<pair<int,double>* >();
  forceMassList->push_back( new pair<int,double>(pdgid, -1.0) );
}

void Photos::forceMass(int pdgid, double mass)
{
  if(mass<0.0)
  {
    Log::Warning()<<"Photos::forceMass: Mass must be > 0.0"<<endl;
    return;
  }
  
  if(!forceMassList) forceMassList = new vector<pair<int,double>* >();
  forceMassList->push_back( new pair<int,double>(pdgid, mass) );
}

} // namespace Photospp
