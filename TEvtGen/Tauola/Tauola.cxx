#include <fstream>
#include <cstring>
#include <vector>
#include "Log.h"
#include "Tauola.h"
#include "TauolaEvent.h"

namespace Tauolapp
{

int    Tauola::m_pdg_id          = 15;
int    Tauola::m_firstDecayMode  = 0; 
int    Tauola::m_secondDecayMode = 0;
bool   Tauola::m_rad             = true;
double Tauola::m_rad_cut_off     = 0.001;
double Tauola::m_iniphy          = 0.1;
double Tauola::m_higgs_scalar_pseudoscalar_mix = M_PI/4;
int    Tauola::m_higgs_scalar_pseudoscalar_pdg = 35;
int    Tauola::m_helPlus  = 0;
int    Tauola::m_helMinus = 0;
double Tauola::m_wtEW     = 0.0;
double Tauola::m_wtEW0    = 0.0;
double Tauola::table11A[NS1][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table1A [NS1][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table2A [NS1][NCOS][4][4] = {{{{0.0}}}};
double Tauola::wtable11A[NS1][NCOS]      = {{0.0}};
double Tauola::wtable1A [NS1][NCOS]      = {{0.0}};
double Tauola::wtable2A [NS1][NCOS]      = {{0.0}};
double Tauola::w0table11A[NS1][NCOS]     = {{0.0}};
double Tauola::w0table1A [NS1][NCOS]     = {{0.0}};
double Tauola::w0table2A [NS1][NCOS]     = {{0.0}};

double Tauola::table11B[NS2][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table1B [NS2][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table2B [NS2][NCOS][4][4] = {{{{0.0}}}};
double Tauola::wtable11B[NS2][NCOS]      = {{0.0}};
double Tauola::wtable1B [NS2][NCOS]      = {{0.0}};
double Tauola::wtable2B [NS2][NCOS]      = {{0.0}};
double Tauola::w0table11B[NS2][NCOS]     = {{0.0}};
double Tauola::w0table1B [NS2][NCOS]     = {{0.0}};
double Tauola::w0table2B [NS2][NCOS]     = {{0.0}};

double Tauola::table11C[NS3][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table1C [NS3][NCOS][4][4] = {{{{0.0}}}};
double Tauola::table2C [NS3][NCOS][4][4] = {{{{0.0}}}};
double Tauola::wtable11C[NS3][NCOS]      = {{0.0}};
double Tauola::wtable1C [NS3][NCOS]      = {{0.0}};
double Tauola::wtable2C [NS3][NCOS]      = {{0.0}};
double Tauola::w0table11C[NS3][NCOS]     = {{0.0}};
double Tauola::w0table1C [NS3][NCOS]     = {{0.0}};
double Tauola::w0table2C [NS3][NCOS]     = {{0.0}};

double Tauola::sminA = 0;
double Tauola::smaxA = 0;

double Tauola::sminB = 0;
double Tauola::smaxB = 0;

double Tauola::sminC = 0;
double Tauola::smaxC = 0;

int    Tauola::ion[3] = {0};
double Tauola::tau_lifetime = .08711;
double Tauola::momentum_conservation_threshold   = 0.1;

Tauola::Particles Tauola::spin_correlation;

Tauola::MomentumUnits Tauola::momentumUnit = Tauola::DEFAULT_MOMENTUM;
Tauola::LengthUnits   Tauola::lengthUnit   = Tauola::DEFAULT_LENGTH;

bool   Tauola::m_is_using_decay_one        = false;
double Tauola::m_decay_one_polarization[3] = {0};
void (*Tauola::m_decay_one_boost_routine)(TauolaParticle*,TauolaParticle*) = NULL;

int     Tauola::buf_incoming_pdg_id = 0;
int     Tauola::buf_outgoing_pdg_id = 0;
double  Tauola::buf_invariant_mass_squared = -1.;
double  Tauola::buf_cosTheta               = 0.;

double  Tauola::buf_R[4][4] = {{0.0}}; //density matrix

double (*Tauola::randomDouble)()                               = Tauola::defaultRandomGenerator;
void   (*Tauola::redefineTauPlusProperties)(TauolaParticle *)  = defaultRedPlus;
void   (*Tauola::redefineTauMinusProperties)(TauolaParticle *) = defaultRedMinus;

/**************************************************************/
void Tauola::setNewCurrents(int mode)
{
  inirchl_(&mode);
}

double Tauola::defaultRandomGenerator(){
  return rand()*1./RAND_MAX;
}

void Tauola::setRandomGenerator(double (*gen)()){
  if(gen==NULL) randomDouble = defaultRandomGenerator;
  else          randomDouble = gen;
}

void Tauola::defaultRedPlus(TauolaParticle *tau)  {}
void Tauola::defaultRedMinus(TauolaParticle *tau) {}

void Tauola::setRedefineTauMinus( void (*fun)(TauolaParticle *) ){
  redefineTauMinusProperties=fun;
}

void Tauola::setRedefineTauPlus ( void (*fun)(TauolaParticle *) ){
  redefineTauPlusProperties=fun;
}

void Tauola::getBornKinematics(int *incoming_pdg_id, int *outgoing_pdg_id, double *invariant_mass_squared,double *cosTheta){
  *incoming_pdg_id = buf_incoming_pdg_id;
  *outgoing_pdg_id = buf_outgoing_pdg_id;
  *invariant_mass_squared = buf_invariant_mass_squared;
  *cosTheta               = buf_cosTheta;
  //  m_R[0][0] to be added in next step;
}

void Tauola::setUnits(MomentumUnits m, LengthUnits l){
  Tauola::momentumUnit = m;
  Tauola::lengthUnit   = l;
}

void Tauola::setTauLifetime(double t){
  tau_lifetime = t;
}

void Tauola::initialize(){
  printf("\n");
  printf(" *************************************\n");
  printf(" *     TAUOLA C++ Interface v1.1.4   *\n");
  printf(" *-----------------------------------*\n");
  printf(" *                                   *\n");
  printf(" *   (c) Nadia    Davidson,   (1,2)  *\n");
  printf(" *       Gizo     Nanava,     (3)    *\n");
  printf(" *       Tomasz   Przedzinski,(4)    *\n");
  printf(" *       Elzbieta Richter-Was,(2,4)  *\n");
  printf(" *       Zbigniew Was         (2,5)  *\n");
  printf(" *                                   *\n");
  printf(" *  1) Unimelb, Melbourne, Australia *\n");
  printf(" *  2)     INP, Krakow, Poland       *\n");
  printf(" *  3) University Bonn, Germany      *\n");
  printf(" *  4)      UJ, Krakow, Poland       *\n");
  printf(" *  5)    CERN, Geneva, Switzerland  *\n");
  printf(" *************************************\n");

  // Turn on all spin correlations
  spin_correlation.setAll(true);

  // Ininitalize tauola-fortran
  f_interface_tauolaInitialize(m_pdg_id,m_firstDecayMode,
                               m_secondDecayMode,m_rad, 
                               m_rad_cut_off, m_iniphy);

  //---------------------------------------------------------------------------
  // Initialize SANC tables
  //---------------------------------------------------------------------------
  cout<<"Reading SANC input files."<<endl;

  ifstream f("table1-1.txt");

  if(!f.is_open()){
    cout<<"File 'table1-1.txt'  missing... skipped."<<endl;
  }
  else{
    string buf;

    cout<<"Reading file 'table1-1.txt'..."<<endl;

    int dbuf1,dbuf2,dbuf3,dbufcos;
    f>>buf>>dbuf1>>dbuf2>>dbuf3>>dbufcos;

    // Check table sizes
    if(dbuf1!=NS1 || dbuf2!=NS2 || dbuf3!=NS3 || dbufcos!=NCOS)  {
      cout<<"mismatched NS1=   "<<dbuf1  <<" <--> "<<NS1<<endl; 
      cout<<"           NS2=   "<<dbuf2  <<" <--> "<<NS2<<endl; 
      cout<<"           NS3=   "<<dbuf3  <<" <--> "<<NS3<<endl; 
      cout<<"           NCOS=  "<<dbufcos<<" <--> "<<NCOS<<endl; 
      return; 
    }

    double buf1,buf2,buf3,buf4,buf5,buf6;
    f>>buf>>buf1>>buf2>>buf3>>buf4>>buf5>>buf6;

    // Set ranges
    if(sminA==0.0){
      sminA=buf1;
      smaxA=buf2;
      sminB=buf3;
      smaxB=buf4;
      sminC=buf5;
      smaxC=buf6;
    }

    // Check ranges
    if(buf1!=sminA || buf2!=smaxA || buf3!=sminB || buf4!=smaxB || buf5!=sminC || buf6!=smaxC) {
      cout<<"mismatched sminA=   "<<buf1<<" <--> "<<sminA<<endl; 
      cout<<"           smaxA=   "<<buf2<<" <--> "<<smaxA<<endl; 
      cout<<"           sminB=   "<<buf3<<" <--> "<<sminB<<endl; 
      cout<<"           smaxB=   "<<buf4<<" <--> "<<smaxB<<endl; 
      cout<<"           sminC=   "<<buf5<<" <--> "<<sminC<<endl; 
      cout<<"           smaxC=   "<<buf6<<" <--> "<<smaxC<<endl; 
      return; 
    }

    // Print out header
    while(!f.eof()){
      char head[255];
      f.getline(head,255);
      if(strcmp(head,"BeginRange1")==0) break;
      cout<<head<<endl;
    }

    // Read table
    for (int i=0;i<NS1;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table1A[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable1A[i][j];
        f>>w0table1A[i][j];
      } // for(j)
    } // for(i)

    // Find 2nd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange2")==0) break;
    }

    // Read table
    for (int i=0;i<NS2;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table1B[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable1B[i][j];
        f>>w0table1B[i][j];
      } // for(j)
    } // for(i)

    // Find 3rd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange3")==0) break;
    }

    // Read table
    for (int i=0;i<NS3;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table1C[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable1C[i][j];
        f>>w0table1C[i][j];
      } // for(j)
    } // for(i)

    // Check for proper file end
    f>>buf;
    if(buf.size() == 0 || strcmp(buf.c_str(),"End") != 0){
      cout<<"...incorrect file version or file incomplete/damaged!"<<endl;

      // In case of the error - do not use tables
      table1A[0][0][0][0] = table1B[0][0][0][0] = table1C[0][0][0][0] = 0.0;
    }
  }  // if (file is open)

  f.close();
  f.open("table2-2.txt");

  if(!f.is_open()){
    cout<<"File 'table2-2.txt'  missing... skipped."<<endl;
  }
  else{
    string buf;

    cout<<"Reading file 'table2-2.txt'..."<<endl;

    int dbuf1,dbuf2,dbuf3,dbufcos;
    f>>buf>>dbuf1>>dbuf2>>dbuf3>>dbufcos;

    // Check table sizes
    if(dbuf1!=NS1 || dbuf2!=NS2 || dbuf3!=NS3 || dbufcos!=NCOS)  {
      cout<<"mismatched NS1=   "<<dbuf1<<" <--> "<<NS1<<endl; 
      cout<<"           NS2=   "<<dbuf2<<" <--> "<<NS2<<endl; 
      cout<<"           NS3=   "<<dbuf3<<" <--> "<<NS3<<endl; 
      cout<<"           NCOS=  "<<dbufcos<<" <--> "<<NCOS<<endl; 
      return; 
    }

    double buf1,buf2,buf3,buf4,buf5,buf6;
    f>>buf>>buf1>>buf2>>buf3>>buf4>>buf5>>buf6;

    // Set ranges
    if(sminA==0.0)
    {
      sminA=buf1;
      smaxA=buf2;
      sminB=buf3;
      smaxB=buf4;
      sminC=buf5;
      smaxC=buf6;
    }

    // Check ranges
    if(buf1!=sminA || buf2!=smaxA || buf3!=sminB || buf4!=smaxB || buf5!=sminC || buf6!=smaxC) {
      cout<<"mismatched sminA=   "<<buf1<<" <--> "<<sminA<<endl; 
      cout<<"           smaxA=   "<<buf2<<" <--> "<<smaxA<<endl; 
      cout<<"           sminB=   "<<buf3<<" <--> "<<sminB<<endl; 
      cout<<"           smaxB=   "<<buf4<<" <--> "<<smaxB<<endl; 
      cout<<"           sminC=   "<<buf5<<" <--> "<<sminC<<endl; 
      cout<<"           smaxC=   "<<buf6<<" <--> "<<smaxC<<endl; 
      return; 
    }

    // Print out header
    while(!f.eof()){
      char head[255];
      f.getline(head,255);
      if(strcmp(head,"BeginRange1")==0) break;
      cout<<head<<endl;
    }

    // Read table
    for (int i=0;i<NS1;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table2A[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable2A[i][j];
        f>>w0table2A[i][j];
      } // for(j)
    } // for(i)

    // Find 2nd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange2")==0) break;
    }

    // Read table
    for (int i=0;i<NS2;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table2B[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable2B[i][j];
        f>>w0table2B[i][j];
      } // for(j)
    } // for(i)

    // Find 3rd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange3")==0) break;
    }

    // Read table
    for (int i=0;i<NS3;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table2C[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable2C[i][j];
        f>>w0table2C[i][j];
      } // for(j)
    } // for(i)

    // Check for proper file end
    f>>buf;
    if(buf.size()==0 || strcmp(buf.c_str(),"End")!=0){
      cout<<"...incorrect file version or file incomplete/damaged!"<<endl;

      // In case of the error - do not use tables
      table2A[0][0][0][0] = table2B[0][0][0][0] = table2C[0][0][0][0] = 0.0;
    }
  } // if (file is open)

  f.close();
  f.open("table11-11.txt");

  if(!f.is_open()){
    cout<<"File 'table11-11.txt' missing... skipped."<<endl;
  }
  else{
    string buf;

    cout<<"Reading file 'table11-11.txt'..."<<endl;

    int dbuf1,dbuf2,dbuf3,dbufcos;
    f>>buf>>dbuf1>>dbuf2>>dbuf3>>dbufcos;

    // Check table sizes
    if(dbuf1!=NS1 || dbuf2!=NS2 || dbuf3!=NS3 || dbufcos!=NCOS)  {
      cout<<"mismatched NS1=   "<<dbuf1<<" <--> "<<NS1<<endl; 
      cout<<"           NS2=   "<<dbuf2<<" <--> "<<NS2<<endl; 
      cout<<"           NS3=   "<<dbuf3<<" <--> "<<NS3<<endl; 
      cout<<"           NCOS=  "<<dbufcos<<" <--> "<<NCOS<<endl; 
      return; 
    }

    double buf1,buf2,buf3,buf4,buf5,buf6;
    f>>buf>>buf1>>buf2>>buf3>>buf4>>buf5>>buf6;

    // Set ranges
    if(sminA==0.0)
    {
      sminA=buf1;
      smaxA=buf2;
      sminB=buf3;
      smaxB=buf4;
      sminC=buf5;
      smaxC=buf6;
    }

    // Check ranges
    if(buf1!=sminA || buf2!=smaxA || buf3!=sminB || buf4!=smaxB || buf5!=sminC || buf6!=smaxC) {
      cout<<"mismatched sminA=   "<<buf1<<" <--> "<<sminA<<endl; 
      cout<<"           smaxA=   "<<buf2<<" <--> "<<smaxA<<endl; 
      cout<<"           sminB=   "<<buf3<<" <--> "<<sminB<<endl; 
      cout<<"           smaxB=   "<<buf4<<" <--> "<<smaxB<<endl; 
      cout<<"           sminC=   "<<buf5<<" <--> "<<sminC<<endl; 
      cout<<"           smaxC=   "<<buf6<<" <--> "<<smaxC<<endl; 
      return; 
    }

    // Print out header
    while(!f.eof()){
      char head[255];
      f.getline(head,255);
      if(strcmp(head,"BeginRange1")==0) break;
      cout<<head<<endl;
    }

    // Read table
    for (int i=0;i<NS1;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table11A[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable11A[i][j];
        f>>w0table11A[i][j];
      } // for(j)
    } // for(i)

    // Find 2nd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange2")==0) break;
    }

    // Read table
    for (int i=0;i<NS2;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table11B[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable11B[i][j];
        f>>w0table11B[i][j];
      } // for(j)
    } // for(i)

    // Find 3rd range
    while(!f.eof()){
      f>>buf;
      if(strcmp(buf.c_str(),"BeginRange3")==0) break;
    }

    // Read table
    for (int i=0;i<NS3;i++){
      for (int j=0;j<NCOS;j++){
        for (int k=0;k<4;k++){
          for (int l=0;l<4;l++){
            f>>table11C[i][j][k][l];
          } // for(l)
        } // for(k)
        f>>wtable11C[i][j];
        f>>w0table11C[i][j];
      } // for(j)
    } // for(i)

    f>>buf;
    if(buf.size()==0 || strcmp(buf.c_str(),"End")!=0){
      cout<<"...incorrect file version or file incomplete/damaged!"<<endl;

      // In case of the error - do not use tables
      table11A[0][0][0][0] = table11B[0][0][0][0] = table11C[0][0][0][0] = 0.0;
    }
  } // if (file is open)

  f.close();
  cout<<endl;
}

void Tauola::decayOne(TauolaParticle *tau, bool undecay, double polx, double poly, double polz)
{
  if(!tau) return;

  if(polx*polx+poly*poly+polz*polz>1)
  {
    Log::Warning()<<"decayOne(): ignoring wrong polarization vector: "<<polx<<" "<<poly<<" "<<polz<<endl;
    polx=poly=polz=0;
  }

  // Let the interface know that we work in the decayOne mode
  m_is_using_decay_one = true;

  m_decay_one_polarization[0] = polx;
  m_decay_one_polarization[1] = poly;
  m_decay_one_polarization[2] = polz;

  // Undecay if needed
  if(tau->hasDaughters())
  {
    if(undecay) tau->undecay();
    else
    {
      m_is_using_decay_one = false;
      return;
    }
  }

  std::vector<TauolaParticle *> list;
  list.push_back(tau);

  // Decay single tau
  TauolaParticlePair t_pair(list);
  t_pair.decayTauPair();
  t_pair.checkMomentumConservation();

  // Revert to normal mode
  m_is_using_decay_one = false;
}

void Tauola::initialise(){

  Log::Warning() <<"Deprecated routine 'Tauola::initialise'"<<endl;
  Log::Warning(0)<<"Use 'Tauola::initialize' instead."<<endl;

  initialize();
  
  // Deprecated routines:  initialise, setInitialisePhy,
  //                       f_interface_tauolaInitialise
}

bool Tauola::isUsingDecayOne()
{
  return m_is_using_decay_one;
}

bool Tauola::isUsingDecayOneBoost()
{
  return (bool) m_decay_one_boost_routine;
}

void Tauola::setBoostRoutine(void (*boost)(TauolaParticle *, TauolaParticle *))
{
  m_decay_one_boost_routine=boost;
}

void Tauola::decayOneBoost(TauolaParticle *mother, TauolaParticle *target)
{
  m_decay_one_boost_routine(mother,target);
}

const double* Tauola::getDecayOnePolarization()
{
  return m_decay_one_polarization;
}

void Tauola::setDecayingParticle(int pdg_id){
  m_pdg_id=pdg_id; 
}

int Tauola::getDecayingParticle(){
  return abs(m_pdg_id);
}

void Tauola::setSameParticleDecayMode(int firstDecayMode){
  m_firstDecayMode=firstDecayMode;
}

void Tauola::setOppositeParticleDecayMode(int secondDecayMode){
  m_secondDecayMode=secondDecayMode;
}

void Tauola::setRadiation(bool rad){
  m_rad=rad;
}

void Tauola::setRadiationCutOff(double rad_cut_off){
  m_rad_cut_off=rad_cut_off;
}


void Tauola::setInitializePhy(double iniphy){
  m_iniphy=iniphy;
}

void Tauola::setInitialisePhy(double iniphy){

  Log::Warning() <<"Deprecated routine 'Tauola::setInitialisePhy'"<<endl;
  Log::Warning(0)<<"Use 'Tauola::setInitializePhy' instead."<<endl;

  setInitializePhy(iniphy);
  
  // Deprecated routines:  initialise, setInitialisePhy,
  //                       f_interface_tauolaInitialise
}

void Tauola::setTauBr(int i, double value)
{
  if(taubra_.nchan==0)
    Log::Warning()<<"setTauBr(): run Tauola::initialize() first."<<endl;
  else if(i<1 || i>taubra_.nchan || value<0.)
    Log::Warning()<<"setTauBr(): Invalid input. Value must be >= 0 and 0 < i <= "<<taubra_.nchan<<endl;
  else taubra_.gamprt[i-1]=(float)value;
}

void Tauola::setTaukle(double bra1,double brk0, double brk0b, double brks)
{
  if(bra1<0 || bra1>1 || brk0<0 ||brk0>1 || brk0b<0 || brk0b>1 || brks<0 ||brks>1)
  {
    Log::Warning()<<"setTaukle(): variables must be in range [0,1]. Ignored."<<endl;
    return;
  }
  taukle_.bra1 =(float)bra1;
  taukle_.brk0 =(float)brk0;
  taukle_.brk0b=(float)brk0b;
  taukle_.brks =(float)brks;
}

double Tauola::getTauMass(){
  return f_getTauMass();
}

double Tauola::getHiggsScalarPseudoscalarMixingAngle(){
  return m_higgs_scalar_pseudoscalar_mix;
}

int Tauola::getHiggsScalarPseudoscalarPDG(){
  return m_higgs_scalar_pseudoscalar_pdg;
}

/** set the mixing angle. coupling: tau~(cos(phi)+isin(phi)gamma5)tau */
void Tauola::setHiggsScalarPseudoscalarMixingAngle(double angle){
  m_higgs_scalar_pseudoscalar_mix=angle;
} 

/** set the pdg code of the higgs particle which tauola should 
    treat as a scalar-pseudoscalar mix  */
void Tauola::setHiggsScalarPseudoscalarPDG(int pdg_code){

  if (particleCharge(pdg_code)!=0.0){       
    Log::Warning()<<"You want to use spin correlations of Higgs for particle of PDGID= "<<pdg_code<<endl
                  <<"This particle has charge="<<particleCharge(pdg_code)<<endl;   
    Log::Fatal("This choice is not appropriate.",0);
  }
  m_higgs_scalar_pseudoscalar_pdg=pdg_code;
} 

int Tauola::getHelPlus(){
  return m_helPlus;
}
int Tauola::getHelMinus(){
  return m_helMinus;
}
double Tauola::getEWwt(){
  return m_wtEW;
}
double Tauola::getEWwt0(){
  return m_wtEW0;
}
void Tauola::setEWwt(double wt, double wt0)
{ 
  m_wtEW  = wt;
  m_wtEW0 = wt0;
}
void Tauola::setHelicities(int minus, int plus)
{ 
  m_helMinus = minus;
  m_helPlus  = plus;
}
void Tauola::setEtaK0sPi(int eta, int k, int pi)
{  
  ion[0] = pi;
  ion[1] = k;
  ion[2] = eta;
}

void Tauola::summary()
{
  int sign_type=100;
  double pol[4] = {0};

  Log::Info()     <<"Tauola::summary(): We use old TAUOLA FORTRAN printout."<<endl;
  Log::Info(false)<<"As a consequence, there is a mismatch in printed TAUOLA version number."<<endl<<endl;

  // Print summary taken from FORTRAN TAUOLA
  dekay_(&sign_type,pol);
}

void Tauola::fill_val(int beg, int end, double* array, double value) 
{
  for (int i = beg; i < end; i++)
    array[i] = value;
}

double Tauola::particleCharge(int idhep)
{
  static double CHARGE[101] = { 0 };
  static int j=0;  
  //--
  //--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
  //--   ding to the PDG particle code...

  if(j==0) // initialization
    {   
      j=1;
      fill_val(0 ,  1, CHARGE, 0.0         );
      fill_val(1 ,  2, CHARGE,-0.3333333333);
      fill_val(2 ,  3, CHARGE, 0.6666666667);
      fill_val(3 ,  4, CHARGE,-0.3333333333);
      fill_val(4 ,  5, CHARGE, 0.6666666667);
      fill_val(5 ,  6, CHARGE,-0.3333333333);
      fill_val(6 ,  7, CHARGE, 0.6666666667);
      fill_val(7 ,  8, CHARGE,-0.3333333333);
      fill_val(8 ,  9, CHARGE, 0.6666666667);
      fill_val(9 , 11, CHARGE, 0.0         );
      fill_val(11 ,12, CHARGE,-1.0         );
      fill_val(12 ,13, CHARGE, 0.0         );
      fill_val(13 ,14, CHARGE,-1.0         );
      fill_val(14, 15, CHARGE, 0.0         );
      fill_val(15 ,16, CHARGE,-1.0         );
      fill_val(16, 17, CHARGE, 0.0         );
      fill_val(17 ,18, CHARGE,-1.0         );
      fill_val(18, 24, CHARGE, 0.0         );
      fill_val(24, 25, CHARGE, 1.0         );
      fill_val(25, 37, CHARGE, 0.0         );
      fill_val(37, 38, CHARGE, 1.0         );
      fill_val(38,101, CHARGE, 0.0         );
    }

  int idabs=abs(idhep);
  double phoch=0.0;

  //--
  //--   Charge of quark, lepton, boson etc....
  if (idabs<=100) phoch=CHARGE[idabs];
  else {
    int Q3= idabs/1000 % 10;
    int Q2= idabs/100  % 10;
    int Q1= idabs/10   % 10;
    if (Q3==0){
      //--
      //-- ...meson...
      if(Q2 % 2==0) phoch=CHARGE[Q2]-CHARGE[Q1];
      else          phoch=CHARGE[Q1]-CHARGE[Q2];
    }
    else{
      //--
      //--   ...diquarks or baryon.
      phoch=CHARGE[Q1]+CHARGE[Q2]+CHARGE[Q3];
    }
  }
  //--
  //--   Find the sign of the charge...
  if (idhep<0.0) phoch=-phoch;
  if (phoch*phoch<0.000001) phoch=0.0;
  
  return phoch;
}

} // namespace Tauolapp
