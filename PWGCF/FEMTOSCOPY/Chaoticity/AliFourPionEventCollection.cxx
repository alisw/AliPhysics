////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliFourPion 
//
//  authors: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliFourPionEventCollection.h"

AliFourPionTrackStruct::AliFourPionTrackStruct():
  fStatus(0),
  fFiltermap(0),
  fId(0),
  fPhi(0),
  fPt(0),
  fMom(0),
  fP(),
  fCharge(0),
  fEta(0),
  fMass(0),
  fDCAXY(0),
  fDCAZ(0),
  fDCA(0),
  fEaccepted(0),
  fKey(0),
  fClusterMap(0),
  fSharedMap(0),
  fX(),
  fTOFhit(0),
  fElectron(0),
  fMuon(0),
  fPion(0),
  fKaon(0),
  fProton(0),
  fLabel(0)// MC
{
  //Default constructor
}
AliFourPionTrackStruct::AliFourPionTrackStruct(const AliFourPionTrackStruct &obj)
  : fStatus(obj.fStatus),
    fFiltermap(obj.fFiltermap),
    fId(obj.fId),
    fPhi(obj.fPhi),
    fPt(obj.fPt),
    fMom(obj.fMom),
    fP(),
    fCharge(obj.fCharge),
    fEta(obj.fEta),
    fMass(obj.fMass),
    fDCAXY(obj.fDCAXY),
    fDCAZ(obj.fDCAZ),
    fDCA(obj.fDCA),
    fEaccepted(obj.fEaccepted),
    fKey(obj.fKey),
    fClusterMap(obj.fClusterMap),
    fSharedMap(obj.fSharedMap),
    fX(),
    fTOFhit(obj.fTOFhit),
    fElectron(obj.fElectron),
    fMuon(obj.fMuon),
    fPion(obj.fPion),
    fKaon(obj.fKaon),
    fProton(obj.fProton),
    fLabel(obj.fLabel)// MC
{
  // copy constructor
}
AliFourPionTrackStruct &AliFourPionTrackStruct::operator=(const AliFourPionTrackStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fStatus = obj.fStatus;
  fFiltermap = obj.fFiltermap;
  fId = obj.fId;
  fPhi = obj.fPhi;
  fPt = obj.fPt;
  fMom = obj.fMom;
  fP[0] = obj.fP[0];
  fP[1] = obj.fP[1];
  fP[2] = obj.fP[2];
  fCharge = obj.fCharge;
  fEta = obj.fEta;
  fMass = obj.fMass;
  fDCAXY = obj.fDCAXY;
  fDCAZ = obj.fDCAZ;
  fDCA = obj.fDCA;
  fEaccepted = obj.fEaccepted;
  fKey = obj.fKey;
  fClusterMap = obj.fClusterMap;
  fSharedMap = obj.fSharedMap;
  fX[0] = obj.fX[0];
  fX[1] = obj.fX[1];
  fX[2] = obj.fX[2];
  fTOFhit = obj.fTOFhit;
  fElectron = obj.fElectron;
  fMuon = obj.fMuon;
  fPion = obj.fPion;
  fKaon = obj.fKaon;
  fProton = obj.fProton;
  fLabel = obj.fLabel;// MC

  return (*this);
}
AliFourPionTrackStruct::~AliFourPionTrackStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliFourPionMCStruct::AliFourPionMCStruct():
  fPx(0),
  fPy(0),
  fPz(0),
  fPtot(0),
  fPdgCode(0),
  fMotherLabel(0)
{
  // Default constructor
}
AliFourPionMCStruct::AliFourPionMCStruct(const AliFourPionMCStruct &obj)
  : fPx(obj.fPx),
    fPy(obj.fPy),
    fPz(obj.fPz),
    fPtot(obj.fPtot),
    fPdgCode(obj.fPdgCode),
    fMotherLabel(obj.fMotherLabel)
{
  // copy constructor
}
AliFourPionMCStruct &AliFourPionMCStruct::operator=(const AliFourPionMCStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fPx = obj.fPx;
  fPy = obj.fPy;
  fPz = obj.fPz;
  fPtot = obj.fPtot;
  fPdgCode = obj.fPdgCode;
  fMotherLabel = obj.fMotherLabel;

  return (*this);
}
AliFourPionMCStruct::~AliFourPionMCStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliFourPionEventStruct::AliFourPionEventStruct():
  fFillStatus(0),
  fNtracks(0),
  fMCarraySize(0),
  fTracks(0),
  fMCtracks(0)
{
  // Default constructor
}
AliFourPionEventStruct::AliFourPionEventStruct(const AliFourPionEventStruct &obj)
  : fFillStatus(obj.fFillStatus),
    fNtracks(obj.fNtracks),
    fMCarraySize(obj.fMCarraySize),
    fTracks(obj.fTracks),
    fMCtracks(obj.fMCtracks)
{
  // copy constructor
}
AliFourPionEventStruct &AliFourPionEventStruct::operator=(const AliFourPionEventStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fFillStatus = obj.fFillStatus;
  fNtracks = obj.fNtracks;
  fMCarraySize = obj.fMCarraySize;
  fTracks = obj.fTracks;
  fMCtracks = obj.fMCtracks;
  
  return (*this);
}
AliFourPionEventStruct::~AliFourPionEventStruct()
{
  // Destructor
  if(fTracks) delete fTracks;
  if(fMCtracks) delete fMCtracks;
}

//_____________________________________________________________________________
AliFourPionEventCollection::AliFourPionEventCollection():
  fFIFO(0),
  fLimit(0),
  fMCLimit(0),
  fEvtStr(0)
{
  // Default constructor
}
AliFourPionEventCollection::AliFourPionEventCollection(Short_t a, Int_t lim, Int_t mcarraylimit, Bool_t MCcase):
  fFIFO(0),
  fLimit(0),
  fMCLimit(0),
  fEvtStr(0)
{
  
  // Main constructor
  SetBuffSize(a);
  
  fEvtStr = new AliFourPionEventStruct[fFIFO];  //allocate pointer array of type particle_event
  fLimit = lim;
  fMCLimit = mcarraylimit;

  for(Int_t ii = 0; ii < fFIFO; ii++){   //Initialize particle table pointers to NULL
    (fEvtStr + ii)->fNtracks = 0;
    (fEvtStr + ii)->fFillStatus = 0;
    (fEvtStr + ii)->fMCarraySize = 0;
    //
    (fEvtStr + ii)->fTracks = NULL;
    (fEvtStr + ii)->fTracks = new AliFourPionTrackStruct[fLimit];
    if(MCcase) (fEvtStr + ii)->fMCtracks = new AliFourPionMCStruct[fMCLimit];
    
  }
}
AliFourPionEventCollection::AliFourPionEventCollection(const AliFourPionEventCollection &obj)
  : fFIFO(obj.fFIFO),
    fLimit(obj.fLimit),
    fMCLimit(obj.fMCLimit),
    fEvtStr(obj.fEvtStr)
{
  // copy constructor
}
AliFourPionEventCollection &AliFourPionEventCollection::operator=(const AliFourPionEventCollection &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fFIFO = obj.fFIFO;
  fLimit = obj.fLimit;
  fMCLimit = obj.fMCLimit;
  fEvtStr = obj.fEvtStr;
  
  return (*this);
}
AliFourPionEventCollection::~AliFourPionEventCollection(){

    for(Int_t i = 0; i < fFIFO; i++){

	if((fEvtStr + i)->fTracks != NULL){
	  delete [] (fEvtStr + i)->fTracks;
	  delete [] (fEvtStr + i)->fMCtracks;
        }	
	
    }
    
    delete [] fEvtStr;
    //remove histos from heap

}


//_____________________________________________________________________________
void AliFourPionEventCollection::FIFOShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  
  
  for(UShort_t i=fFIFO-1 ; i > 0; i--){
    for(Int_t j=0; j<(fEvtStr + i-1)->fNtracks; j++) (fEvtStr + i)->fTracks[j] = (fEvtStr + i-1)->fTracks[j];
    for(Int_t j=0; j<(fEvtStr + i-1)->fMCarraySize; j++) (fEvtStr + i)->fMCtracks[j] = (fEvtStr + i-1)->fMCtracks[j];
    
    (fEvtStr + i)->fFillStatus = (fEvtStr + i-1)->fFillStatus;
    (fEvtStr + i)->fNtracks = (fEvtStr + i-1)->fNtracks;
    (fEvtStr + i)->fMCarraySize = (fEvtStr + i-1)->fMCarraySize;
    
  }// fifo loop


  (fEvtStr)->fNtracks=0;
  (fEvtStr)->fFillStatus=0;
  (fEvtStr)->fMCarraySize=0;
}
