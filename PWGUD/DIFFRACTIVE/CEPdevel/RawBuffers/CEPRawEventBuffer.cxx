/// ////////////////////////////////////////////////////////////////////////////
///
/// CEP raw event buffer
///
/// structure to hold event information
///
//______________________________________________________________________________
#include "CEPRawEventBuffer.h"
#include "CEPTrackBuffer.h"
#include "AliCEPBase.h"
#include "AliESDtrack.h"

ClassImp(CEPRawEventBuffer)

//______________________________________________________________________________
CEPRawEventBuffer::CEPRawEventBuffer()
    : TObject()
    , fEventNumber(CEPTrackBuffer::kdumval)
    , fnTracks(0)  
    , fnCaloTracks(0)  
    , fADTotalMult(0.)  
    , fADTotalTime(0.)  
    , fADTotalCharge(0.)  
    , fFMDTotalMult(0.)  
    , fV0TotalMult(0.)  
    , fV0TotalTime(0.)  
    , fV0TotalCharge(0.)  
    , fV0TotalSigWidth(0.)  
    , fEMCTotalAmplitude(0.)  
    , fEMCTotalTime(0.)  
    , fPHOSTotalAmplitude(0.)
    , fPHOSTotalTime(0.)
    , fCEPRawTracks(new TObjArray())
    , fCEPRawCaloClusterTracks(new TObjArray())
    , fADCellBuffer(0x0)
    , fV0Buffer(0x0)
    , fEMCalBuffer(0x0)
    , fPHOSBuffer(0x0)
    , fFMDBuffer(0x0)
{

}

//______________________________________________________________________________
CEPRawEventBuffer::CEPRawEventBuffer(const CEPRawEventBuffer& eb)
    : TObject(eb)
    , fEventNumber(eb.fEventNumber)
    , fnTracks(eb.fnTracks)  
    , fnCaloTracks(eb.fnCaloTracks)  
    , fADTotalMult(eb.fADTotalMult)  
    , fADTotalTime(eb.fADTotalTime)  
    , fADTotalCharge(eb.fADTotalCharge)  
    , fFMDTotalMult(eb.fFMDTotalMult)  
    , fV0TotalMult(eb.fV0TotalMult)  
    , fV0TotalTime(eb.fV0TotalTime)  
    , fV0TotalCharge(eb.fV0TotalCharge)  
    , fV0TotalSigWidth(eb.fV0TotalSigWidth)  
    , fEMCTotalAmplitude(eb.fEMCTotalAmplitude)  
    , fEMCTotalTime(eb.fEMCTotalTime)  
    , fPHOSTotalAmplitude(eb.fPHOSTotalAmplitude)  
    , fPHOSTotalTime(eb.fPHOSTotalTime)  
    , fCEPRawTracks(new TObjArray(*eb.fCEPRawTracks))
    , fCEPRawCaloClusterTracks(new TObjArray(*eb.fCEPRawCaloClusterTracks))
    , fADCellBuffer(new CEPRawADBuffer(*eb.fADCellBuffer))
    , fV0Buffer(new CEPRawV0Buffer(*eb.fV0Buffer))
    , fEMCalBuffer(new CEPRawCaloBuffer(*eb.fEMCalBuffer))
    , fPHOSBuffer(new CEPRawCaloBuffer(*eb.fPHOSBuffer))
    , fFMDBuffer(new CEPRawFMDBuffer(*eb.fFMDBuffer))
{

}

//______________________________________________________________________________
CEPRawEventBuffer::~CEPRawEventBuffer()
{
    // delete fCEPRawTracks and all the tracks it contains
    if (fCEPRawTracks) {
        fCEPRawTracks->SetOwner(kTRUE);
        fCEPRawTracks->Clear();
        delete fCEPRawTracks;
        fCEPRawTracks = 0x0;
    }
    // delete fCEPRawCaloClusterTracks and all the tracks it contains
    if (fCEPRawCaloClusterTracks) {
        fCEPRawCaloClusterTracks->SetOwner(kTRUE);
        fCEPRawCaloClusterTracks->Clear();
        delete fCEPRawCaloClusterTracks;
        fCEPRawCaloClusterTracks = 0x0;
    }    
    if (fADCellBuffer) { delete fADCellBuffer;  fADCellBuffer = 0x0; }
    if (fV0Buffer)     { delete fADCellBuffer;  fADCellBuffer = 0x0; }
    if (fEMCalBuffer)  { delete fEMCalBuffer;   fEMCalBuffer  = 0x0; }
    if (fPHOSBuffer)   { delete fPHOSBuffer;    fPHOSBuffer   = 0x0; }
    if (fFMDBuffer)    { delete fFMDBuffer;     fFMDBuffer    = 0x0; }
}

//______________________________________________________________________________
void CEPRawEventBuffer::Copy(TObject &obj) const 
{
    // interface to TOBject::Copy
    // Copies the content of this into obj!
    // bascially obj = *this

    if(this==&obj)return;
    CEPRawEventBuffer *robj = dynamic_cast<CEPRawEventBuffer*>(&obj);
    if(!robj)return; // not a CEPRawEventBuffer
    *robj = *this;
    return;
}

///
/// Assignment operator.
///
//__________________________________________________________________________
CEPRawEventBuffer & CEPRawEventBuffer::operator =(const CEPRawEventBuffer& source)  
{
    if (&source == this) return *this;

    TObject::operator=(source); // don't forget to invoke the base class' assignment operator

    Reset();
  
    // standard members
    fEventNumber = source.fEventNumber;
    fnTracks     = source.fnTracks;
    fnCaloTracks = source.fnCaloTracks;
    // AD
    fADTotalMult   = source.fADTotalMult;
    fADTotalTime   = source.fADTotalTime;
    fADTotalCharge = source.fADTotalCharge;
    // FMD
    fFMDTotalMult = source.fFMDTotalMult;
    // V0
    fV0TotalMult     = source.fV0TotalMult;
    fV0TotalTime     = source.fV0TotalTime;
    fV0TotalCharge   = source.fV0TotalCharge;
    fV0TotalSigWidth = source.fV0TotalSigWidth;
    // Calo buffers
    fEMCTotalAmplitude  = source.fEMCTotalAmplitude;
    fEMCTotalTime       = source.fEMCTotalTime;
    fPHOSTotalAmplitude = source.fPHOSTotalAmplitude;
    fPHOSTotalTime      = source.fPHOSTotalTime;

    // object members
    fADCellBuffer = source.fADCellBuffer;
    fV0Buffer     = source.fV0Buffer;
    fEMCalBuffer  = source.fEMCalBuffer;
    fPHOSBuffer   = source.fPHOSBuffer;
    fFMDBuffer    = source.fFMDBuffer;

    // Object arrays members
    if (source.fCEPRawTracks) {
        if (fCEPRawTracks) *fCEPRawTracks = *(source.fCEPRawTracks); 
        else fCEPRawTracks = new TObjArray(*(source.fCEPRawTracks));
    } else {
        fCEPRawTracks->SetOwner(kTRUE);
        fCEPRawTracks->Clear();
        delete fCEPRawTracks;
        fCEPRawTracks = 0x0;
    }

    if (source.fCEPRawCaloClusterTracks) {
        if (fCEPRawCaloClusterTracks) 
            *fCEPRawCaloClusterTracks = *(source.fCEPRawCaloClusterTracks); 
        else fCEPRawCaloClusterTracks = new TObjArray(*(source.fCEPRawCaloClusterTracks));
    } else {
        fCEPRawCaloClusterTracks->SetOwner(kTRUE);
        fCEPRawCaloClusterTracks->Clear();
        delete fCEPRawCaloClusterTracks;
        fCEPRawCaloClusterTracks = 0x0;
    }

    return *this;
}

//______________________________________________________________________________
void CEPRawEventBuffer::Reset()
{
    // reset all private variables
    fEventNumber = AliCEPBase::kdumval;
    fnTracks     = 0;
    fnCaloTracks = 0;
    // AD
    fADTotalMult   = 0.;
    fADTotalTime   = 0.;
    fADTotalCharge = 0.;
    // FMD
    fFMDTotalMult = 0.;
    // V0
    fV0TotalMult     = 0.;
    fV0TotalTime     = 0.;
    fV0TotalCharge   = 0.;
    fV0TotalSigWidth = 0.;
    // Calo buffers
    fEMCTotalAmplitude  = 0.;
    fEMCTotalTime       = 0.;
    fPHOSTotalAmplitude = 0.;
    fPHOSTotalTime      = 0.;

    // clear the track list
    fCEPRawTracks->SetOwner(kTRUE);
    fCEPRawTracks->Clear();
    // clear the calo cluster list
    fCEPRawCaloClusterTracks->SetOwner(kTRUE);
    fCEPRawCaloClusterTracks->Clear();

    fADCellBuffer->Reset();
    fV0Buffer->Reset();
    fEMCalBuffer->Reset();
    fPHOSBuffer->Reset();
    fFMDBuffer->Reset();
}

//______________________________________________________________________________
void CEPRawEventBuffer::AddTrack(CEPRawTrackBuffer* trk)
{
    // add track to next element
    fCEPRawTracks->Add(trk);
    fnTracks++;
}

//______________________________________________________________________________
void CEPRawEventBuffer::AddTrack(CEPRawCaloClusterTrack* caloTrk)
{
    // AddTrack overwritten for filling calo cluster Array
    fCEPRawCaloClusterTracks->Add(caloTrk);
    fnCaloTracks++;
}

//______________________________________________________________________________
CEPRawTrackBuffer* CEPRawEventBuffer::GetTrack(UInt_t ind)
{
    // initialize the result track
    CEPRawTrackBuffer *trk = NULL;

    if (fCEPRawTracks->GetEntries() > ind) 
    {
        trk = (CEPRawTrackBuffer*) fCEPRawTracks->At(ind);
    }

    return trk;
}

//______________________________________________________________________________
CEPRawCaloClusterTrack* CEPRawEventBuffer::GetCaloClusterTrack(UInt_t ind)
{
    // initialize the result track
    CEPRawCaloClusterTrack *caloTrk = NULL;

    if (fCEPRawCaloClusterTracks->GetEntries() > ind) 
    {
        caloTrk = (CEPRawCaloClusterTrack*) fCEPRawCaloClusterTracks->At(ind);
    }

    return caloTrk;
}

//______________________________________________________________________________
Bool_t CEPRawEventBuffer::RemoveTrack(UInt_t ind)
{
    // initialize the result track
    Bool_t done = kFALSE;
    CEPRawTrackBuffer *trk = NULL;

    if (fCEPRawTracks->GetEntries() > ind) {
        trk = (CEPRawTrackBuffer*) fCEPRawTracks->RemoveAt(ind);
        fCEPRawTracks->Compress();

        // update track counters
        fnTracks--;
        
        done = kTRUE;
    }
    return done;
}

//______________________________________________________________________________
Bool_t CEPRawEventBuffer::RemoveCaloCluster(UInt_t ind)
{
    // initialize the result track
    Bool_t done = kFALSE;
    CEPRawCaloClusterTrack *caloTrk = NULL;

    if (fCEPRawCaloClusterTracks->GetEntries() > ind) {
        caloTrk = (CEPRawCaloClusterTrack*) fCEPRawCaloClusterTracks->RemoveAt(ind);
        fCEPRawCaloClusterTracks->Compress();

        // update track counters
        fnCaloTracks--;
        
        done = kTRUE;
    }
    return done;
}

//______________________________________________________________________________
void CEPRawEventBuffer::SetEventVariables(AliESDEvent* ESDobj)
{
    // Set event number
    this->SetEventNumber( ESDobj->GetEventNumberInFile() );

    // create sub-detector buffers
    fADCellBuffer = new CEPRawADBuffer();
    fV0Buffer     = new CEPRawV0Buffer();
    fEMCalBuffer  = new CEPRawCaloBuffer();
    fPHOSBuffer   = new CEPRawCaloBuffer();
    fFMDBuffer    = new CEPRawFMDBuffer();
    // fill the buffers
    fADCellBuffer->SetADVariables(ESDobj->GetADData()); 
    fV0Buffer->SetV0Variables(ESDobj->GetVZEROData()); 
    fEMCalBuffer->SetCaloVariables(ESDobj->GetEMCALCells()); 
    fPHOSBuffer->SetCaloVariables(ESDobj->GetPHOSCells()); 
    fFMDBuffer->SetFMDVariables(ESDobj->GetFMDData()); 
    
    // AD setter
    this->SetTotalADMult(fADCellBuffer->GetADTotalMultiplicity() );
    this->SetTotalADTime(fADCellBuffer->GetADTotalTime() );
    this->SetTotalADCharge(fADCellBuffer->GetADTotalCharge() );
    // FMD setter
    this->SetTotalFMDMult(fFMDBuffer->GetFMDTotalMultiplicity() );
    // V0 setter
    this->SetTotalV0Mult(fV0Buffer->GetV0TotalMultiplicity() );
    this->SetTotalV0Time(fV0Buffer->GetV0TotalCharge() );
    this->SetTotalV0Charge(fV0Buffer->GetV0TotalCharge() );
    this->SetTotalV0SigWidth(fV0Buffer->GetV0TotalSignalWidth() );
    // Calo Setter
    this->SetTotalEMCAmplitude(fEMCalBuffer->GetCaloTotalAmplitude() );
    this->SetTotalEMCTime(fEMCalBuffer->GetCaloTotalTime() );
    this->SetTotalPHOSAmplitude(fPHOSBuffer->GetCaloTotalAmplitude() );
    this->SetTotalPHOSTime(fPHOSBuffer->GetCaloTotalTime() );

    // fill in the raw track information
    UInt_t nTracks = ESDobj->GetNumberOfTracks();
    for (UInt_t i(0); i<nTracks; i++)
    {
        // Initialize track
        AliESDtrack* aliTrk = NULL;
        CEPRawTrackBuffer* trk = new CEPRawTrackBuffer();
        // get a track from the ESD object
        aliTrk = (AliESDtrack*)ESDobj->GetTrack(i);
        // fill raw track buffer
        trk->SetTrackVariables(aliTrk);
        // add track to the CEPRawEventBuffer
        AddTrack(trk);
    }
    // fill in raw calo cluster information
    UInt_t nCaloTracks = ESDobj->GetNumberOfCaloClusters();
    for (UInt_t i(0); i<nCaloTracks; i++)
    {
        // Initialize calo cluster
        AliESDCaloCluster* aliCluster = NULL;
        CEPRawCaloClusterTrack* caloTrk = new CEPRawCaloClusterTrack();
        // get calo cluster from the ESD object
        aliCluster = (AliESDCaloCluster*)ESDobj->GetCaloCluster(i);
        // fill raw calo-cluster buffer
        caloTrk->SetCaloClusterVariables(aliCluster);
        // add track to the 
        AddTrack(caloTrk);
    }
}
  
//______________________________________________________________________________
