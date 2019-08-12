/// ////////////////////////////////////////////////////////////////////////////
///
/// CEP raw event buffer
///
/// structure to hold event information
///
//______________________________________________________________________________
#include <iostream>
#include "CEPRawEventBuffer.h"
#include "CEPTrackBuffer.h"
#include "AliCEPBase.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

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
    , fADCellBuffer(new CEPRawADBuffer())
    , fV0Buffer(new CEPRawV0Buffer())
    , fEMCalBuffer(new CEPRawCaloBuffer())
    , fPHOSBuffer(new CEPRawCaloBuffer())
    , fFMDBuffer(new CEPRawFMDBuffer())
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
    if(fADCellBuffer) { delete fADCellBuffer; fADCellBuffer = 0x0; }
    if(fV0Buffer)     { delete fV0Buffer;     fV0Buffer = 0x0;     }
    if(fEMCalBuffer)  { delete fEMCalBuffer;  fEMCalBuffer = 0x0;  }
    if(fPHOSBuffer)   { delete fPHOSBuffer;   fPHOSBuffer = 0x0;   }
    if(fFMDBuffer)    { delete fFMDBuffer;    fFMDBuffer = 0x0;    }
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
void CEPRawEventBuffer::AddCaloTrack(CEPRawCaloClusterTrack* caloTrk)
{
    // AddTrack overwritten for filling calo cluster Array
    fCEPRawCaloClusterTracks->Add(caloTrk);
    fnCaloTracks++;
}

//______________________________________________________________________________
CEPRawTrackBuffer* CEPRawEventBuffer::GetTrack(UInt_t ind)
{
    // initialize the result track
    CEPRawTrackBuffer *trk = 0x0;

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
    CEPRawCaloClusterTrack *caloTrk = 0x0;

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
    CEPRawTrackBuffer *trk = 0x0;

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
    CEPRawCaloClusterTrack *caloTrk = 0x0;

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
void CEPRawEventBuffer::SetEventVariables(AliESDEvent* ESDobj, TArrayI* TTindices)
{
    this->Reset();
    // Set event number
    this->SetEventNumber( ESDobj->GetEventNumberInFile() );

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
    AliESDVertex* vertex = (AliESDVertex*)ESDobj->GetPrimaryVertex(); 
    for (UInt_t i(0); i<nTracks; i++)
    {
        // Initialize track
        AliESDtrack* aliTrk = 0x0;
        CEPRawTrackBuffer* trk = new CEPRawTrackBuffer();
        // get a track from the ESD object
        aliTrk = (AliESDtrack*)ESDobj->GetTrack(i);
        // fill raw track buffer
        trk->SetTrackVariables(aliTrk, vertex);
        // add track to the CEPRawEventBuffer
        AddTrack(trk);
    }
    // fill in raw calo cluster information
    UInt_t nCaloTracks = ESDobj->GetNumberOfCaloClusters();
    for (UInt_t i(0); i<nCaloTracks; i++)
    {
        // Initialize calo cluster
        AliESDCaloCluster* aliCluster = 0x0;
        CEPRawCaloClusterTrack* caloTrk = new CEPRawCaloClusterTrack();
        // get calo cluster from the ESD object
        aliCluster = (AliESDCaloCluster*)ESDobj->GetCaloCluster(i);
        // fill raw calo-cluster buffer
        caloTrk->SetCaloClusterVariables(aliCluster, (AliESDCaloCells*)ESDobj->GetEMCALCells(),
                                         ESDobj, TTindices);
        // add track to the 
        AddCaloTrack(caloTrk);
    }
}
  
//______________________________________________________________________________
