#include "TauolaHEPEVTParticle.h"

#include "Log.h"

namespace Tauolapp
{

TauolaHEPEVTParticle::~TauolaHEPEVTParticle()
{
  // Cleanup particles that do not belong to event
  for(unsigned int i=0;i<cache.size();i++)
    if(cache[i]->m_barcode<0)
      delete cache[i];
}

TauolaHEPEVTParticle::TauolaHEPEVTParticle(int pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de){
  m_px = px;
  m_py = py;
  m_pz = pz;
  m_e  = e;
  m_generated_mass = m;

  m_pdgid  = pdgid;
  m_status = status;

  m_first_mother   = ms;
  m_second_mother  = me;
  m_daughter_start = ds;
  m_daughter_end   = de;

  m_barcode = -1;
  m_event = NULL;
}

void TauolaHEPEVTParticle::undecay(){

  Log::Info()<<"TauolaHEPEVTParticle::undecay not implemented for HEPEVT"<<endl;
}

void TauolaHEPEVTParticle::setMothers(vector<TauolaParticle*> mothers){

  // If this particle has not yet been added to the event record
  // then add it to the mothers' event record
  if(m_barcode<0 && mothers.size()>0)
  {
    TauolaHEPEVTEvent *evt = ((TauolaHEPEVTParticle*)mothers[0])->m_event;
    evt->addParticle(this);
  }

  if(mothers.size()>2) Log::Fatal("TauolaHEPEVTParticle::setMothers: HEPEVT does not allow more than two mothers!");
  
  if(mothers.size()>0) m_first_mother  = mothers[0]->getBarcode();
  if(mothers.size()>1) m_second_mother = mothers[1]->getBarcode();
}

void TauolaHEPEVTParticle::setDaughters(vector<TauolaParticle*> daughters){

  // This particle must be inside some event record to be able to add daughters
  if(m_event==NULL) Log::Fatal("TauolaHEPEVTParticle::setDaughters: particle not inside event record.");

  int beg = 65535, end = -1;

  for(unsigned int i=0;i<daughters.size();i++)
  {
    int bc = daughters[i]->getBarcode();
    if(bc<0) Log::Fatal("TauolaHEPEVTParticle::setDaughters: all daughters has to be in event record first");
    if(bc<beg) beg = bc;
    if(bc>end) end = bc;
  }
  if(end == -1) beg = -1;

  m_daughter_start = beg;
  m_daughter_end   = end;
}

std::vector<TauolaParticle*> TauolaHEPEVTParticle::getMothers(){

  std::vector<TauolaParticle*> mothers;

  TauolaParticle *p1 = NULL;
  TauolaParticle *p2 = NULL;
  
  if(m_first_mother>=0)  p1 = m_event->getParticle(m_first_mother);
  if(m_second_mother>=0) p2 = m_event->getParticle(m_second_mother);

  if(p1) mothers.push_back(p1);
  if(p2) mothers.push_back(p2);

  return mothers;
}

// WARNING: this method also corrects daughter indices 
//          if such were not defined
std::vector<TauolaParticle*> TauolaHEPEVTParticle::getDaughters(){

  std::vector<TauolaParticle*> daughters;

  if(!m_event) return daughters;

  // Check if m_daughter_start and m_daughter_end are set
  // If not - try to get list of daughters from event
  if(m_daughter_end<0)
  {
    int min_d=65535, max_d=-1;
    for(int i=0;i<m_event->getParticleCount();i++)
    {
      if(m_event->getParticle(i)->isDaughterOf(this))
      {
        if(i<min_d) min_d = i;
        if(i>max_d) max_d = i;
      }
    }
    if(max_d>=0)
    {
      m_daughter_start = min_d;
      m_daughter_end   = max_d;
      m_status         = 2;
    }
  }

  // If m_daughter_end is still not set - there are no daughters
  // Otherwsie - get daughters
  if(m_daughter_end>=0)
  {
    for(int i=m_daughter_start;i<=m_daughter_end;i++)
    {
      TauolaParticle *p = m_event->getParticle(i);
      if(p==NULL)
      {
        Log::Warning()<<"TauolaHEPEVTParticle::getDaughters(): No particle with index "<<i<<endl;
        return daughters;
      }

      daughters.push_back(p);
    }
  }

  return daughters;
}

void TauolaHEPEVTParticle::checkMomentumConservation(){

  if(!m_event)           return;
  if(m_daughter_end < 0) return;

  TauolaHEPEVTParticle *buf = m_event->getParticle(m_daughter_start);

  int first_mother_idx  = buf->getFirstMotherIndex();
  int second_mother_idx = buf->getSecondMotherIndex();

  double px =0.0, py =0.0, pz =0.0, e =0.0;
  double px2=0.0, py2=0.0, pz2=0.0, e2=0.0;

  for(int i=m_daughter_start;i<=m_daughter_end;i++)
  {
    buf = m_event->getParticle(i);
    px += buf->getPx();
    py += buf->getPy();
    pz += buf->getPz();
    e  += buf->getE ();
  }

  if(first_mother_idx>=0)
  {
    buf = m_event->getParticle(first_mother_idx);
    px2 += buf->getPx();
    py2 += buf->getPy();
    pz2 += buf->getPz();
    e2  += buf->getE();
  }
  
  if(second_mother_idx>=0)
  {
    buf = m_event->getParticle(second_mother_idx);
    px2 += buf->getPx();
    py2 += buf->getPy();
    pz2 += buf->getPz();
    e2  += buf->getE();
  }
  // 3-momentum  // test HepMC style
  double dp = sqrt( (px-px2)*(px-px2) + (py-py2)*(py-py2) + (pz-pz2)*(pz-pz2) );
  // virtuality test as well.
  double m1 = sqrt( fabs( e*e   - px*px   - py*py   - pz*pz   ) );
  double m2 = sqrt( fabs( e2*e2 - px2*px2 - py2*py2 - pz2*pz2 ) );

  if( fabs(m1-m2) > 0.0001 || dp > 0.0001*(e+e2))
  {
    Log::RedirectOutput( Log::Warning()<<"Momentum not conserved in vertex: " );
    if(first_mother_idx >=0) m_event->getParticle(first_mother_idx) ->print();
    if(second_mother_idx>=0) m_event->getParticle(second_mother_idx)->print();
    for(int i=m_daughter_start;i<=m_daughter_end;i++) m_event->getParticle(i)->print();
    Log::RevertOutput();
  }
}

TauolaHEPEVTParticle * TauolaHEPEVTParticle::createNewParticle(
                        int pdg_id, int status, double mass,
                        double px, double py, double pz, double e){

  // New particles created using this method are added to cache
  // They will be deleted when this particle will be deleted

  cache.push_back(new TauolaHEPEVTParticle(pdg_id,status,px,py,pz,e,mass,-1,-1,-1,-1));
  return cache.back();
}

bool TauolaHEPEVTParticle::isDaughterOf(TauolaHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc==m_first_mother || bc==m_second_mother) return true;

  return false;
}

bool TauolaHEPEVTParticle::isMotherOf  (TauolaHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc>=m_daughter_start && bc<=m_daughter_end) return true;

  return false;
}

void TauolaHEPEVTParticle::print(){
  char buf[256];
  sprintf(buf,"P: (%2i) %6i %2i | %11.4e %11.4e %11.4e %11.4e | %11.4e | M: %2i %2i | D: %2i %2i\n",
          m_barcode, m_pdgid, m_status, m_px, m_py, m_pz, m_e, m_generated_mass,
          m_first_mother, m_second_mother,   m_daughter_start, m_daughter_end);

  cout<<buf;
}

/******** Getter and Setter methods: ***********************/

void TauolaHEPEVTParticle::setPdgID(int pdg_id){
  m_pdgid = pdg_id;
}

void TauolaHEPEVTParticle::setStatus(int status){
  m_status = status;
}

void TauolaHEPEVTParticle::setMass(double mass){
  m_generated_mass = mass;
}

int TauolaHEPEVTParticle::getPdgID(){
  return m_pdgid;
}

int TauolaHEPEVTParticle::getStatus(){
  return m_status;
}

double TauolaHEPEVTParticle::getMass(){
  return m_generated_mass;
}

inline double TauolaHEPEVTParticle::getPx(){
  return m_px;
}

inline double TauolaHEPEVTParticle::getPy(){
  return m_py;
}

double TauolaHEPEVTParticle::getPz(){
  return m_pz;
}

double TauolaHEPEVTParticle::getE(){
  return m_e;
}

void TauolaHEPEVTParticle::setPx(double px){
  m_px = px;
}

void TauolaHEPEVTParticle::setPy(double py){
  m_py = py;
}


void TauolaHEPEVTParticle::setPz(double pz){
  m_pz = pz;
}

void TauolaHEPEVTParticle::setE(double e){
  m_e = e;
}

int TauolaHEPEVTParticle::getBarcode(){
  return m_barcode;
}

void TauolaHEPEVTParticle::setBarcode(int barcode){
  m_barcode = barcode;
}

void TauolaHEPEVTParticle::setEvent(TauolaHEPEVTEvent *event){
  m_event = event;
}

int TauolaHEPEVTParticle::getFirstMotherIndex(){
  return m_first_mother;
}

int TauolaHEPEVTParticle::getSecondMotherIndex(){
  return m_second_mother;
}

int TauolaHEPEVTParticle::getDaughterRangeStart(){
  return m_daughter_start;
}

int TauolaHEPEVTParticle::getDaughterRangeEnd(){
  return m_daughter_end;
}

} // namespace Tauolapp
