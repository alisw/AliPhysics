////////////////////////////////////////////////////////////////////////////////
//
//  This class provides storage for event and track information which 
//  are used for same-event as well as mixed-event analyses in AliChaoticity 
//
//  authors: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliChaoticityEventCollection.h"

AliChaoticityTrackStruct::AliChaoticityTrackStruct():
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
AliChaoticityTrackStruct::AliChaoticityTrackStruct(const AliChaoticityTrackStruct &obj)
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
AliChaoticityTrackStruct &AliChaoticityTrackStruct::operator=(const AliChaoticityTrackStruct &obj) 
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
AliChaoticityTrackStruct::~AliChaoticityTrackStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliChaoticityPairStruct::AliChaoticityPairStruct():
  fP1(),
  fP2(),
  fE1(0),
  fE2(0),
  fCharge1(0),
  fCharge2(0),
  fIndex1(0),
  fIndex2(0),
  fQinv(0),
  fKey1(0),
  fKey2(0),
  fLabel1(0),
  fLabel2(0),
  fP1MC(),
  fP2MC()
{
  //Default constructor
}
AliChaoticityPairStruct::AliChaoticityPairStruct(const AliChaoticityPairStruct &obj)
  : fP1(),
    fP2(),
    fE1(obj.fE1),
    fE2(obj.fE2),
    fCharge1(obj.fCharge1),
    fCharge2(obj.fCharge2),
    fIndex1(obj.fIndex1),
    fIndex2(obj.fIndex2),
    fQinv(obj.fQinv),
    fKey1(obj.fKey1),
    fKey2(obj.fKey2),
    fLabel1(obj.fLabel1),
    fLabel2(obj.fLabel2),
    fP1MC(),
    fP2MC()
{
  // copy constructor
}
AliChaoticityPairStruct &AliChaoticityPairStruct::operator=(const AliChaoticityPairStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fP1[0] = obj.fP1[0];
  fP1[1] = obj.fP1[1];
  fP1[2] = obj.fP1[2];
  fP2[0] = obj.fP2[0];
  fP2[1] = obj.fP2[1];
  fP2[2] = obj.fP2[2];
  fE1 = obj.fE1;
  fE2 = obj.fE2;
  fCharge1 = obj.fCharge1;
  fCharge2 = obj.fCharge2;
  fIndex1 = obj.fIndex1;
  fIndex2 = obj.fIndex2;
  fQinv = obj.fQinv;
  fKey1 = obj.fKey1;
  fKey2 = obj.fKey2;
  fLabel1 = obj.fLabel1;
  fLabel2 = obj.fLabel2;
  fP1MC[0] = obj.fP1MC[0];
  fP1MC[1] = obj.fP1MC[1];
  fP1MC[2] = obj.fP1MC[2];
  fP2MC[0] = obj.fP2MC[0];
  fP2MC[1] = obj.fP2MC[1];
  fP2MC[2] = obj.fP2MC[2];
  
  return (*this);
}
AliChaoticityPairStruct::~AliChaoticityPairStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliChaoticityNormPairStruct::AliChaoticityNormPairStruct():
  fCharge1(0),
  fCharge2(0),
  fIndex1(0),
  fIndex2(0),
  fKey1(0),
  fKey2(0)
{
  //Default constructor
}
AliChaoticityNormPairStruct::AliChaoticityNormPairStruct(const AliChaoticityNormPairStruct &obj)
  : fCharge1(obj.fCharge1),
    fCharge2(obj.fCharge2),
    fIndex1(obj.fIndex1),
    fIndex2(obj.fIndex2),
    fKey1(obj.fKey1),
    fKey2(obj.fKey2)
{
  // copy constructor
}
AliChaoticityNormPairStruct &AliChaoticityNormPairStruct::operator=(const AliChaoticityNormPairStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fCharge1 = obj.fCharge1;
  fCharge2 = obj.fCharge2;
  fIndex1 = obj.fIndex1;
  fIndex2 = obj.fIndex2;
  fKey1 = obj.fKey1;
  fKey2 = obj.fKey2;

  return (*this);
}
AliChaoticityNormPairStruct::~AliChaoticityNormPairStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliChaoticityMCStruct::AliChaoticityMCStruct():
  fPx(0),
  fPy(0),
  fPz(0),
  fPtot(0)
{
  // Default constructor
}
AliChaoticityMCStruct::AliChaoticityMCStruct(const AliChaoticityMCStruct &obj)
  : fPx(obj.fPx),
    fPy(obj.fPy),
    fPz(obj.fPz),
    fPtot(obj.fPtot)
{
  // copy constructor
}
AliChaoticityMCStruct &AliChaoticityMCStruct::operator=(const AliChaoticityMCStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fPx = obj.fPx;
  fPy = obj.fPy;
  fPz = obj.fPz;
  fPtot = obj.fPtot;
  
  return (*this);
}
AliChaoticityMCStruct::~AliChaoticityMCStruct()
{
  // Destructor
}

//_____________________________________________________________________________
AliChaoticityEventStruct::AliChaoticityEventStruct():
  fFillStatus(0),
  fNtracks(0),
  fNpairsSE(0),
  fNpairsME(0),
  fMCarraySize(0),
  fTracks(0),
  fPairsSE(0),
  fPairsME(0),
  fMCtracks(0)
{
  // Default constructor
}
AliChaoticityEventStruct::AliChaoticityEventStruct(const AliChaoticityEventStruct &obj)
  : fFillStatus(obj.fFillStatus),
    fNtracks(obj.fNtracks),
    fNpairsSE(obj.fNpairsSE),
    fNpairsME(obj.fNpairsME),
    fMCarraySize(obj.fMCarraySize),
    fTracks(obj.fTracks),
    fPairsSE(obj.fPairsSE),
    fPairsME(obj.fPairsME),
    fMCtracks(obj.fMCtracks)
{
  // copy constructor
}
AliChaoticityEventStruct &AliChaoticityEventStruct::operator=(const AliChaoticityEventStruct &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fFillStatus = obj.fFillStatus;
  fNtracks = obj.fNtracks;
  fNpairsSE = obj.fNpairsSE;
  fNpairsME = obj.fNpairsME;
  fMCarraySize = obj.fMCarraySize;
  fTracks = obj.fTracks;
  fPairsSE = obj.fPairsSE;
  fPairsME = obj.fPairsME;
  fMCtracks = obj.fMCtracks;
  
  return (*this);
}
AliChaoticityEventStruct::~AliChaoticityEventStruct()
{
  // Destructor
  if(fTracks) delete fTracks;
  if(fPairsSE) delete fPairsSE;
  if(fPairsME) delete fPairsME;
  if(fMCtracks) delete fMCtracks;
}

//_____________________________________________________________________________
AliChaoticityEventCollection::AliChaoticityEventCollection():
  fFIFO(0),
  fLimit(0),
  fPairLimit(0),
  fMCLimit(0),
  fEvtStr(0)
{
  // Default constructor
}
AliChaoticityEventCollection::AliChaoticityEventCollection(Short_t a, Int_t lim, Int_t plimit, Int_t mcarraylimit, Bool_t MCcase):
  fFIFO(0),
  fLimit(0),
  fPairLimit(0),
  fMCLimit(0),
  fEvtStr(0)
{
  
  // Main constructor
  SetBuffSize(a);
  
  fEvtStr = new AliChaoticityEventStruct[fFIFO];  //allocate pointer array of type particle_event
  fLimit = lim;
  fPairLimit = plimit;
  fMCLimit = mcarraylimit;

  for(Int_t ii = 0; ii < fFIFO; ii++){   //Initialize particle table pointers to NULL
    (fEvtStr + ii)->fNtracks = 0;
    (fEvtStr + ii)->fNpairsSE = 0;
    (fEvtStr + ii)->fNpairsME = 0;
    (fEvtStr + ii)->fFillStatus = 0;
    (fEvtStr + ii)->fMCarraySize = 0;
    //
    (fEvtStr + ii)->fTracks = NULL;
    (fEvtStr + ii)->fTracks = new AliChaoticityTrackStruct[fLimit];
    (fEvtStr + ii)->fPairsSE = NULL;
    (fEvtStr + ii)->fPairsSE = new AliChaoticityPairStruct[fPairLimit];
    (fEvtStr + ii)->fPairsME = NULL;
    (fEvtStr + ii)->fPairsME = new AliChaoticityPairStruct[Int_t(2*fPairLimit)];
    if(MCcase) (fEvtStr + ii)->fMCtracks = new AliChaoticityMCStruct[fMCLimit];
    
  }
}
AliChaoticityEventCollection::AliChaoticityEventCollection(const AliChaoticityEventCollection &obj)
  : fFIFO(obj.fFIFO),
    fLimit(obj.fLimit),
    fPairLimit(obj.fPairLimit),
    fMCLimit(obj.fMCLimit),
    fEvtStr(obj.fEvtStr)
{
  // copy constructor
}
AliChaoticityEventCollection &AliChaoticityEventCollection::operator=(const AliChaoticityEventCollection &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;

  fFIFO = obj.fFIFO;
  fLimit = obj.fLimit;
  fPairLimit = obj.fPairLimit;
  fMCLimit = obj.fMCLimit;
  fEvtStr = obj.fEvtStr;
  
  return (*this);
}
AliChaoticityEventCollection::~AliChaoticityEventCollection(){

    for(Int_t i = 0; i < fFIFO; i++){

	if((fEvtStr + i)->fTracks != NULL){
	  delete [] (fEvtStr + i)->fTracks;
	  delete [] (fEvtStr + i)->fPairsSE;
	  delete [] (fEvtStr + i)->fPairsME;
	  delete [] (fEvtStr + i)->fMCtracks;
        }	
	
    }
    
    delete [] fEvtStr;
    //remove histos from heap

}


//_____________________________________________________________________________
void AliChaoticityEventCollection::FIFOShift(){ //Shift elements in FIFO by one and clear last element in FIFO 
  
  
  for(UShort_t i=fFIFO-1 ; i > 0; i--){
    for(Int_t j=0; j<(fEvtStr + i-1)->fNtracks; j++) (fEvtStr + i)->fTracks[j] = (fEvtStr + i-1)->fTracks[j];
    for(Int_t j=0; j<(fEvtStr + i-1)->fNpairsSE; j++) (fEvtStr + i)->fPairsSE[j] = (fEvtStr + i-1)->fPairsSE[j];
    for(Int_t j=0; j<(fEvtStr + i-1)->fMCarraySize; j++) (fEvtStr + i)->fMCtracks[j] = (fEvtStr + i-1)->fMCtracks[j];

    (fEvtStr + i)->fFillStatus = (fEvtStr + i-1)->fFillStatus;
    (fEvtStr + i)->fNtracks = (fEvtStr + i-1)->fNtracks;
    (fEvtStr + i)->fNpairsSE = (fEvtStr + i-1)->fNpairsSE;
    (fEvtStr + i)->fNpairsME = 0;
    (fEvtStr + i)->fMCarraySize = (fEvtStr + i-1)->fMCarraySize;
    
  }// fifo loop


  (fEvtStr)->fNtracks=0;
  (fEvtStr)->fNpairsSE=0;
  (fEvtStr)->fNpairsME=0;
  (fEvtStr)->fFillStatus=0;
  (fEvtStr)->fMCarraySize=0;
}
