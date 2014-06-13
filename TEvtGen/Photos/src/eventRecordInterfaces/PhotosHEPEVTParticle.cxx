#include "PhotosHEPEVTParticle.h"
#include "Log.h"

namespace Photospp
{

PhotosHEPEVTParticle::~PhotosHEPEVTParticle()
{
  // Cleanup particles that do not belong to event
  for(unsigned int i=0;i<cache.size();i++)
    if(cache[i]->m_barcode<0)
      delete cache[i];
}

PhotosHEPEVTParticle::PhotosHEPEVTParticle(int pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de){
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

/** Add a new daughter to this particle */
void PhotosHEPEVTParticle::addDaughter(PhotosParticle* daughter)
{
  if(!m_event) Log::Fatal("PhotosHEPEVTParticle::addDaughter - particle not in event record");

  std::vector<PhotosParticle*> mothers = daughter->getMothers();
  
  mothers.push_back( (PhotosParticle*)this );

  daughter->setMothers(mothers);
  
  int bc = daughter->getBarcode();

  if(m_daughter_end < 0)
  {
    m_daughter_start = bc;
    m_daughter_end   = bc;
  }
  // if it's in the middle of the event record
  else if(m_daughter_end != bc-1)
  {
    PhotosHEPEVTParticle *newPart = m_event->getParticle(bc);
    
    // Move all particles one spot down the list, to make place for new particle
    for(int i=bc-1;i>m_daughter_end;i--)
    {
      PhotosHEPEVTParticle *move = m_event->getParticle(i);
      move->setBarcode(i+1);
      m_event->setParticle(i+1,move);
    }

    m_daughter_end++;
    newPart->setBarcode(m_daughter_end);
    m_event->setParticle(m_daughter_end,newPart);
    
    // Now: correct all pointers before new particle
    for(int i=0;i<m_daughter_end;i++)
    {
      PhotosHEPEVTParticle *check = m_event->getParticle(i);
      int m = check->getDaughterRangeEnd();
      if(m!=-1 && m>m_daughter_end)
      {
        check->setDaughterRangeEnd(m+1);
        check->setDaughterRangeStart(check->getDaughterRangeStart()+1);
      }
    }
  }
  else m_daughter_end = bc;
}

void PhotosHEPEVTParticle::setMothers(vector<PhotosParticle*> mothers){

  // If this particle has not yet been added to the event record
  // then add it to the mothers' event record
  if(m_barcode<0 && mothers.size()>0)
  {
    PhotosHEPEVTEvent *evt = ((PhotosHEPEVTParticle*)mothers[0])->m_event;
    evt->addParticle(this);
  }

  if(mothers.size()>2) Log::Fatal("PhotosHEPEVTParticle::setMothers: HEPEVT does not allow more than two mothers!");
  
  if(mothers.size()>0) m_first_mother  = mothers[0]->getBarcode();
  if(mothers.size()>1) m_second_mother = mothers[1]->getBarcode();
}

void PhotosHEPEVTParticle::setDaughters(vector<PhotosParticle*> daughters){

  // This particle must be inside some event record to be able to add daughters
  if(m_event==NULL) Log::Fatal("PhotosHEPEVTParticle::setDaughters: particle not inside event record.");

  int beg = 65535, end = -1;

  for(unsigned int i=0;i<daughters.size();i++)
  {
    int bc = daughters[i]->getBarcode();
    if(bc<0) Log::Fatal("PhotosHEPEVTParticle::setDaughters: all daughters has to be in event record first");
    if(bc<beg) beg = bc;
    if(bc>end) end = bc;
  }
  if(end == -1) beg = -1;

  m_daughter_start = beg;
  m_daughter_end   = end;
}

std::vector<PhotosParticle*> PhotosHEPEVTParticle::getMothers(){

  std::vector<PhotosParticle*> mothers;

  PhotosParticle *p1 = NULL;
  PhotosParticle *p2 = NULL;

  // 21.XI.2013: Some generators set same mother ID in both indices if there is only one mother  
  if(m_first_mother == m_second_mother) m_second_mother = -1;

  if(m_first_mother>=0)  p1 = m_event->getParticle(m_first_mother);
  if(m_second_mother>=0) p2 = m_event->getParticle(m_second_mother);

  if(p1) mothers.push_back(p1);
  if(p2) mothers.push_back(p2);

  return mothers;
}

// WARNING: this method also corrects daughter indices 
//          if such were not defined
std::vector<PhotosParticle*> PhotosHEPEVTParticle::getDaughters(){

  std::vector<PhotosParticle*> daughters;

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
      PhotosParticle *p = m_event->getParticle(i);
      if(p==NULL)
      {
        Log::Warning()<<"PhotosHEPEVTParticle::getDaughters(): No particle with index "<<i<<endl;
        return daughters;
      }

      daughters.push_back(p);
    }
  }

  return daughters;
}

std::vector<PhotosParticle*> PhotosHEPEVTParticle::getAllDecayProducts()
{
  std::vector<PhotosParticle*> list;

  if(!hasDaughters()) // if this particle has no daughters
    return list;

  std::vector<PhotosParticle*> daughters = getDaughters();
  
  // copy daughters to list of all decay products
  list.insert(list.end(),daughters.begin(),daughters.end());
  
  // Now, get all daughters recursively, without duplicates.
  // That is, for each daughter:
  // 1)  get list of her daughters
  // 2)  for each particle on this list:
  //  a) check if it is already on the list
  //  b) if it's not, add her to the end of the list
  for(unsigned int i=0;i<list.size();i++)
  {
    std::vector<PhotosParticle*> daughters2 = list[i]->getDaughters();

    if(!list[i]->hasDaughters()) continue;
    for(unsigned int j=0;j<daughters2.size();j++)
    {
      bool add=true;
      for(unsigned int k=0;k<list.size();k++)
        if( daughters2[j]->getBarcode() == list[k]->getBarcode() )
        {
          add=false;
          break;
        }

      if(add) list.push_back(daughters2[j]);
    }
  }

  return list;
}

bool PhotosHEPEVTParticle::checkMomentumConservation(){

  if(!m_event)           return true;
  if(m_daughter_end < 0) return true;

  PhotosHEPEVTParticle *buf = m_event->getParticle(m_daughter_start);

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
    cout<<"  TO: "<<endl;
    for(int i=m_daughter_start;i<=m_daughter_end;i++) m_event->getParticle(i)->print();
    Log::RevertOutput();
    return false;
  }
  
  return true;
}

PhotosHEPEVTParticle * PhotosHEPEVTParticle::createNewParticle(
                        int pdg_id, int status, double mass,
                        double px, double py, double pz, double e){

  // New particles created using this method are added to cache
  // They will be deleted when this particle will be deleted

  cache.push_back(new PhotosHEPEVTParticle(pdg_id,status,px,py,pz,e,mass,-1,-1,-1,-1));
  return cache.back();
}

void PhotosHEPEVTParticle::createHistoryEntry()
{
  Log::Warning()<<"PhotosParticle::createHistoryEntry() not implemented for HEPEVT."<<endl;
}

void PhotosHEPEVTParticle::createSelfDecayVertex(PhotosParticle *out)
{
  Log::Warning()<<"PhotosHEPEVTParticle::createSelfDecayVertex() not implemented for HEPEVT."<<endl;
}

bool PhotosHEPEVTParticle::isDaughterOf(PhotosHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc==m_first_mother || bc==m_second_mother) return true;

  return false;
}

bool PhotosHEPEVTParticle::isMotherOf  (PhotosHEPEVTParticle *p)
{
  int bc = p->getBarcode();
  if(bc>=m_daughter_start && bc<=m_daughter_end) return true;

  return false;
}

void PhotosHEPEVTParticle::print(){
  char buf[256];
  sprintf(buf,"P: (%2i) %6i %2i | %11.4e %11.4e %11.4e %11.4e | %11.4e | M: %2i %2i | D: %2i %2i\n",
          m_barcode, m_pdgid, m_status, m_px, m_py, m_pz, m_e, m_generated_mass,
          m_first_mother, m_second_mother,   m_daughter_start, m_daughter_end);

  cout<<buf;
}

/******** Getter and Setter methods: ***********************/

void PhotosHEPEVTParticle::setPdgID(int pdg_id){
  m_pdgid = pdg_id;
}

void PhotosHEPEVTParticle::setStatus(int status){
  m_status = status;
}

void PhotosHEPEVTParticle::setMass(double mass){
  m_generated_mass = mass;
}

int PhotosHEPEVTParticle::getPdgID(){
  return m_pdgid;
}

int PhotosHEPEVTParticle::getStatus(){
  return m_status;
}

double PhotosHEPEVTParticle::getMass(){
  return m_generated_mass;
}

inline double PhotosHEPEVTParticle::getPx(){
  return m_px;
}

inline double PhotosHEPEVTParticle::getPy(){
  return m_py;
}

double PhotosHEPEVTParticle::getPz(){
  return m_pz;
}

double PhotosHEPEVTParticle::getE(){
  return m_e;
}

void PhotosHEPEVTParticle::setPx(double px){
  m_px = px;
}

void PhotosHEPEVTParticle::setPy(double py){
  m_py = py;
}


void PhotosHEPEVTParticle::setPz(double pz){
  m_pz = pz;
}

void PhotosHEPEVTParticle::setE(double e){
  m_e = e;
}

int PhotosHEPEVTParticle::getBarcode(){
  return m_barcode;
}

void PhotosHEPEVTParticle::setBarcode(int barcode){
  m_barcode = barcode;
}

void PhotosHEPEVTParticle::setEvent(PhotosHEPEVTEvent *event){
  m_event = event;
}

int PhotosHEPEVTParticle::getFirstMotherIndex(){
  return m_first_mother;
}

int PhotosHEPEVTParticle::getSecondMotherIndex(){
  return m_second_mother;
}

int PhotosHEPEVTParticle::getDaughterRangeStart(){
  return m_daughter_start;
}

int PhotosHEPEVTParticle::getDaughterRangeEnd(){
  return m_daughter_end;
}

} // namespace Photospp
