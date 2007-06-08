/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the ESD class
//   This is the class to deal with during the phisical analysis of data
//   This class is generated directly by the reconstruction methods
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliESD.h"


ClassImp(AliESDRun)  
 
//______________________________________________________________________________
AliESDRun::AliESDRun() :
  fRunNumber(0),
  fPeriodNumber(0),
  fRecoVersion(0), 
  fMagneticField(0)
{
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=0.;
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=0.;
}

//______________________________________________________________________________
AliESDRun::AliESDRun(const AliESDRun &esd) :
  TObject(esd),
  fRunNumber(esd.fRunNumber),
  fPeriodNumber(esd.fPeriodNumber),
  fRecoVersion(esd.fRecoVersion),
  fMagneticField(esd.fMagneticField)
{ 
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=esd.fDiamondXY[i];
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=esd.fDiamondCovXY[i];
}

//______________________________________________________________________________
AliESDRun& AliESDRun::operator=(const AliESDRun &esd)
{ 
  if(this!=&esd) {
    TObject::operator=(esd);
    fRunNumber=esd.fRunNumber;
    fPeriodNumber=esd.fPeriodNumber;
    fRecoVersion=esd.fRecoVersion;
    fMagneticField=esd.fMagneticField;
    for (Int_t i=0; i<2; i++) fDiamondXY[i]=esd.fDiamondXY[i];
    for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=esd.fDiamondCovXY[i];
  } 
  return *this;
}

//______________________________________________________________________________
void AliESDRun::Print(const Option_t *) const
{
  printf("Mean vertex in RUN %d: X=%.4f Y=%.4f cm\n",
	 GetRunNumber(),GetDiamondX(),GetDiamondY());
  printf("Magnetic field = %f T\n",
	 GetMagneticField());
  printf("Event from reconstruction version %d \n",fRecoVersion);
}

void AliESDRun::Reset() 
{
  fRunNumber = 0;
  fPeriodNumber = 0;
  fRecoVersion = 0;
  fMagneticField = 0;
  for (Int_t i=0; i<2; i++) fDiamondXY[i]=0.;
  for (Int_t i=0; i<3; i++) fDiamondCovXY[i]=0.;
}

ClassImp(AliESDHeader)

//______________________________________________________________________________
AliESDHeader::AliESDHeader() :
  TObject(),
  fTriggerMask(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fEventType(0),
  fEventNumberInFile(0),
  fBunchCrossNumber(0),
  fTriggerCluster(0)
{
}


AliESDHeader::AliESDHeader(const AliESDHeader &header) :
  TObject(header),
  fTriggerMask(header.fTriggerMask),
  fOrbitNumber(header.fOrbitNumber),
  fTimeStamp(header.fTimeStamp),
  fEventType(header.fEventType),
  fEventNumberInFile(header.fEventNumberInFile),
  fBunchCrossNumber(header.fBunchCrossNumber),
  fTriggerCluster(header.fTriggerCluster)
{
}

AliESDHeader& AliESDHeader::operator=(const AliESDHeader &header)
{ 
  if(this!=&header) {
    TObject::operator=(header);
    fTriggerMask = header.fTriggerMask;
    fOrbitNumber = header.fOrbitNumber;
    fTimeStamp = header.fTimeStamp;
    fEventType = header.fEventType;
    fEventNumberInFile = header.fEventNumberInFile;
    fBunchCrossNumber = header.fBunchCrossNumber;
    fTriggerCluster = header.fTriggerCluster;
  } 
  return *this;
}



//______________________________________________________________________________
void AliESDHeader::Reset()
{
  fTriggerMask       = 0;
  fOrbitNumber       = 0;
  fTimeStamp         = 0;
  fEventType         = 0;
  fEventNumberInFile = 0;
  fBunchCrossNumber  = 0;
  fTriggerCluster    = 0;
}

//______________________________________________________________________________
void AliESDHeader::Print(const Option_t *) const
{
  printf("Event # %d in file Bunch crossing # %d Orbit # %d Trigger %lld \n",
	 GetEventNumberInFile(),
	 GetBunchCrossNumber(),
	 GetOrbitNumber(),
	 GetTriggerMask());
}

ClassImp(AliESDZDC)

//______________________________________________________________________________
AliESDZDC::AliESDZDC() :
  TObject(),
  fZDCN1Energy(0),
  fZDCP1Energy(0),
  fZDCN2Energy(0),
  fZDCP2Energy(0),
  fZDCEMEnergy(0),
  fZDCParticipants(0)
{
}

AliESDZDC::AliESDZDC(const AliESDZDC& zdc) :
  TObject(zdc),
  fZDCN1Energy(zdc.fZDCN1Energy),
  fZDCP1Energy(zdc.fZDCP1Energy),
  fZDCN2Energy(zdc.fZDCN2Energy),
  fZDCP2Energy(zdc.fZDCP2Energy),
  fZDCEMEnergy(zdc.fZDCEMEnergy),
  fZDCParticipants(zdc.fZDCParticipants)
{
}

AliESDZDC& AliESDZDC::operator=(const AliESDZDC&zdc)
{
  if(this!=&zdc) {
    TObject::operator=(zdc);
    fZDCN1Energy = zdc.fZDCN1Energy;
    fZDCP1Energy = zdc.fZDCP1Energy;
    fZDCN2Energy = zdc.fZDCN2Energy;
    fZDCP2Energy = zdc.fZDCP2Energy;
    fZDCEMEnergy = zdc.fZDCEMEnergy;
    fZDCParticipants = zdc.fZDCParticipants;
  } 
  return *this;
}


//______________________________________________________________________________
void AliESDZDC::Reset()
{
  fZDCN1Energy=0;
  fZDCP1Energy=0;
  fZDCN2Energy=0;
  fZDCP2Energy=0;
  fZDCEMEnergy=0;
  fZDCParticipants=0;
}

//______________________________________________________________________________
void AliESDZDC::Print(const Option_t *) const
{
}


ClassImp(AliESDTZERO)

//______________________________________________________________________________
AliESDTZERO::AliESDTZERO() :
  TObject(),
  fT0zVertex(0),
  fT0timeStart(0)   
{
  for(int i = 0;i<24;i++)fT0time[i] = fT0amplitude[i] = 0;
}

AliESDTZERO::AliESDTZERO(const AliESDTZERO &tzero ) :
  TObject(tzero),
  fT0zVertex(tzero.fT0zVertex),
  fT0timeStart(tzero.fT0timeStart)   
{
  for(int i = 0;i<24;i++){
    fT0time[i] = tzero.fT0time[i]; 
    fT0amplitude[i] = tzero.fT0amplitude[i];
  }
}

AliESDTZERO& AliESDTZERO::operator=(const AliESDTZERO& tzero){
  if(this!=&tzero) {
    TObject::operator=(tzero);
    fT0zVertex = tzero.fT0zVertex;
    fT0timeStart = tzero.fT0timeStart;   
    for(int i = 0;i<24;i++){
      fT0time[i] = tzero.fT0time[i]; 
      fT0amplitude[i] = tzero.fT0amplitude[i];
    }
  } 
  return *this;
}

//______________________________________________________________________________
void AliESDTZERO::Reset()
{
  fT0zVertex = 0;  
  fT0timeStart = 0;
  for(int i = 0;i<24;i++)fT0time[i] = fT0amplitude[i] = 0;
}

//______________________________________________________________________________
void AliESDTZERO::Print(const Option_t *) const
{
}

ClassImp(AliESDCaloTrigger)

AliESDCaloTrigger::AliESDCaloTrigger() : 
  TNamed(),
  fTriggerAmplitudes(0x0),
  fTriggerPosition(0x0)
{
}

AliESDCaloTrigger::AliESDCaloTrigger(const AliESDCaloTrigger &ctrig) : 
  TNamed(ctrig),
  fTriggerAmplitudes(ctrig.fTriggerAmplitudes),
  fTriggerPosition(ctrig.fTriggerPosition)
{
}

AliESDCaloTrigger::~AliESDCaloTrigger()
{
  delete fTriggerAmplitudes; fTriggerAmplitudes = 0;
  delete fTriggerPosition; fTriggerPosition = 0;
}

AliESDCaloTrigger& AliESDCaloTrigger::operator=(const AliESDCaloTrigger& ctrig)
{
  if(this!=&ctrig) {
    TNamed::operator=(ctrig);
    // CKB dont't want to create leak if fTriggerAmp points to 
    // somthing already, use new with placement
    if(fTriggerAmplitudes){
      fTriggerAmplitudes = new(fTriggerAmplitudes) TArrayF(*ctrig.fTriggerAmplitudes);
    }
    else{
      fTriggerAmplitudes = new TArrayF(*ctrig.fTriggerAmplitudes);
    }
    if(fTriggerPosition){
      fTriggerPosition = new(fTriggerPosition) TArrayF(*ctrig.fTriggerPosition);
    }
    else{
      fTriggerPosition = new TArrayF(*ctrig.fTriggerPosition);
    }
  } 
  return *this;
}

void AliESDCaloTrigger::Reset()
{
  
  if( fTriggerAmplitudes){  
    printf("%s %d Size %d",(char*)__FILE__,__LINE__,fTriggerAmplitudes->GetSize());
    fTriggerAmplitudes->Reset();
// delete fTriggerAmplitudes;
  }
  if( fTriggerPosition){
    fTriggerPosition->Reset();
  // delete fTriggerPosition;
  }
}



ClassImp(AliESD)

//______________________________________________________________________________
AliESD::AliESD():
  fESDObjects(new TList()),
  fESDRun(0),
  fHeader(0),
  fESDZDC(0),
  fESDFMD(0),
  fESDVZERO(0),
  fESDTZERO(0),
  fSPDVertex(0),
  fPrimaryVertex(0),
  fSPDMult(0),
  fPHOSTrigger(0),
  fEMCALTrigger(0),
  fTracks(0),
  fMuonTracks(0),
  fPmdTracks(0),
  fTrdTracks(0),
  fV0s(0),  
  fCascades(0),
  fKinks(0),
  fCaloClusters(0),
  fErrorLogs(0),
  fEMCALClusters(0), 
  fFirstEMCALCluster(-1),
  fPHOSClusters(0), 
  fFirstPHOSCluster(-1)
{
}
//______________________________________________________________________________
AliESD::AliESD(const AliESD& esd):
  TObject(esd),
  fESDObjects(new TList()),
  fESDRun(new AliESDRun(*esd.fESDRun)),
  fHeader(new AliESDHeader(*esd.fHeader)),
  fESDZDC(new AliESDZDC(*esd.fESDZDC)),
  fESDFMD(new AliESDFMD(*esd.fESDFMD)),
  fESDVZERO(new AliESDVZERO(*esd.fESDVZERO)),
  fESDTZERO(new AliESDTZERO(*esd.fESDTZERO)),
  fSPDVertex(new AliESDVertex(*esd.fSPDVertex)),
  fPrimaryVertex(new AliESDVertex(*esd.fPrimaryVertex)),
  fSPDMult(new AliMultiplicity(*esd.fSPDMult)),
  fPHOSTrigger(new AliESDCaloTrigger(*esd.fPHOSTrigger)),
  fEMCALTrigger(new AliESDCaloTrigger(*esd.fEMCALTrigger)),
  fTracks(new TClonesArray(*esd.fTracks)),
  fMuonTracks(new TClonesArray(*esd.fMuonTracks)),
  fPmdTracks(new TClonesArray(*esd.fPmdTracks)),
  fTrdTracks(new TClonesArray(*esd.fTrdTracks)),
  fV0s(new TClonesArray(*esd.fV0s)),  
  fCascades(new TClonesArray(*esd.fCascades)),
  fKinks(new TClonesArray(*esd.fKinks)),
  fCaloClusters(new TClonesArray(*esd.fCaloClusters)),
  fErrorLogs(new TClonesArray(*esd.fErrorLogs)),
  fEMCALClusters(esd.fEMCALClusters), 
  fFirstEMCALCluster(esd.fFirstEMCALCluster),
  fPHOSClusters(esd.fPHOSClusters), 
  fFirstPHOSCluster(esd.fFirstPHOSCluster)

{
  // CKB init in the constructor list and only add here ...
  AddObject(fESDRun);
  AddObject(fHeader);
  AddObject(fESDZDC);
  AddObject(fESDFMD);
  AddObject(fESDVZERO);
  AddObject(fESDTZERO);
  AddObject(fSPDVertex);
  AddObject(fPrimaryVertex);
  AddObject(fSPDMult);
  AddObject(fPHOSTrigger);
  AddObject(fEMCALTrigger);
  AddObject(fTracks);
  AddObject(fMuonTracks);
  AddObject(fPmdTracks);
  AddObject(fTrdTracks);
  AddObject(fV0s);
  AddObject(fCascades);
  AddObject(fKinks);
  AddObject(fCaloClusters);
  AddObject(fErrorLogs);

  GetStdContent();

}

//______________________________________________________________________________
AliESD & AliESD::operator=(const AliESD& source) {

  // Assignment operator

  if(&source == this) return *this;
  TObject::operator=(source);

  fESDRun = new AliESDRun(*source.fESDRun);
  fHeader = new AliESDHeader(*source.fHeader);
  fESDZDC = new AliESDZDC(*source.fESDZDC);
  fESDFMD = new AliESDFMD(*source.fESDFMD);
  fESDVZERO = new AliESDVZERO(*source.fESDVZERO);
  fESDTZERO = new AliESDTZERO(*source.fESDTZERO);
  fSPDVertex = new AliESDVertex(*source.fSPDVertex);
  fPrimaryVertex = new AliESDVertex(*source.fPrimaryVertex);
  fSPDMult = new AliMultiplicity(*source.fSPDMult);
  fPHOSTrigger = new AliESDCaloTrigger(*source.fPHOSTrigger);
  fEMCALTrigger = new AliESDCaloTrigger(*source.fEMCALTrigger);
  fTracks = new TClonesArray(*source.fTracks);
  fMuonTracks = new TClonesArray(*source.fMuonTracks);
  fPmdTracks = new TClonesArray(*source.fPmdTracks);
  fTrdTracks = new TClonesArray(*source.fTrdTracks);
  fV0s = new TClonesArray(*source.fV0s);
  fCascades = new TClonesArray(*source.fCascades);
  fKinks = new TClonesArray(*source.fKinks);
  fCaloClusters = new TClonesArray(*source.fCaloClusters);
  fErrorLogs = new TClonesArray(*source.fErrorLogs);

  // CKB this way?? or 
  // or AddObject(  fESDZDC = new AliESDZDC(*source.fESDZDC));

  fESDObjects = new TList();
  AddObject(fESDRun);
  AddObject(fHeader);
  AddObject(fESDZDC);
  AddObject(fESDFMD);
  AddObject(fESDVZERO);
  AddObject(fESDTZERO);
  AddObject(fSPDVertex);
  AddObject(fPrimaryVertex);
  AddObject(fSPDMult);
  AddObject(fPHOSTrigger);
  AddObject(fEMCALTrigger);
  AddObject(fTracks);
  AddObject(fMuonTracks);
  AddObject(fPmdTracks);
  AddObject(fTrdTracks);
  AddObject(fV0s);
  AddObject(fCascades);
  AddObject(fKinks);
  AddObject(fCaloClusters);
  AddObject(fErrorLogs);


  fEMCALClusters = source.fEMCALClusters;
  fFirstEMCALCluster = source.fFirstEMCALCluster;
  fPHOSClusters = source.fPHOSClusters;
  fFirstPHOSCluster = source.fFirstPHOSCluster;



  return *this;

}


//______________________________________________________________________________
AliESD::~AliESD()
{
  //
  // Standard destructor
  //

  delete fESDObjects;
  fESDObjects = 0;

  // everthing on the list gets deleted automatically

  /*
  fHLTConfMapTracks.Delete();
  fHLTHoughTracks.Delete();
  fMuonTracks.Delete();  
  fPmdTracks.Delete();
  fTrdTracks.Delete();
  fV0s.Delete();
  fCascades.Delete();
  fKinks.Delete();
  fCaloClusters.Delete();
  */
//   fEMCALTriggerPosition->Delete();
//   fEMCALTriggerAmplitudes->Delete();
//   fPHOSTriggerPosition->Delete();
//   fPHOSTriggerAmplitudes->Delete();
//   delete fEMCALTriggerPosition;
//   delete fEMCALTriggerAmplitudes;
//   delete fPHOSTriggerPosition;
//   delete fPHOSTriggerAmplitudes;

}

//______________________________________________________________________________
void AliESD::Reset()
{

  if(fESDRun) fESDRun->Reset();
  if(fHeader) fHeader->Reset();
  if(fESDZDC) fESDZDC->Reset();
  if(fESDFMD) fESDFMD->Clear(); // why clear.... need consistend names
  // if(fESDVZERO) fESDVZERO->; // NOT IMPLEMENTED 
  //  if(fESDVZERO) new (fESDVZERO) AliESDVZERO();
  if(fESDTZERO) fESDTZERO->Reset(); 
  // CKB no clear/reset implemented
  if(fSPDVertex){
    new (fSPDVertex) AliESDVertex();
    fSPDVertex->SetName("SPDVertex");
  }
  if(fPrimaryVertex){
    new (fPrimaryVertex) AliESDVertex();
    fPrimaryVertex->SetName("PrimaryVertex");
  }
  if(fSPDMult)new (fSPDMult) AliMultiplicity();
  if(fPHOSTrigger)fPHOSTrigger->Reset(); 
  if(fEMCALTrigger)fEMCALTrigger->Reset(); 
  if(fTracks)fTracks->Clear();
  if(fMuonTracks)fMuonTracks->Clear();
  if(fPmdTracks)fPmdTracks->Clear();
  if(fTrdTracks)fTrdTracks->Clear();
  if(fV0s)fV0s->Clear();
  if(fCascades)fCascades->Clear();
  if(fKinks)fKinks->Clear();
  if(fCaloClusters)fCaloClusters->Clear();
  if(fErrorLogs) fErrorLogs->Clear();


  fEMCALClusters=0; 
  fFirstEMCALCluster=-1; 
  fPHOSClusters=0; 
  fFirstPHOSCluster=-1; 
}

Int_t AliESD::AddV0(const AliESDv0 *v) {
  //
  // Add V0
  //
  TClonesArray &fv = *fV0s;
  Int_t idx=fV0s->GetEntriesFast();
  new(fv[idx]) AliESDv0(*v);
  return idx;
}  

//______________________________________________________________________________
void AliESD::Print(Option_t *) const 
{
  //
  // Print header information of the event
  //
  printf("ESD run information\n");
  printf("Event # in file %d Bunch crossing # %d Orbit # %d Period # %d Run # %d Trigger %lld Magnetic field %f \n",
	 GetEventNumberInFile(),
	 GetBunchCrossNumber(),
	 GetOrbitNumber(),
	 GetPeriodNumber(),
	 GetRunNumber(),
	 GetTriggerMask(),
	 GetMagneticField() );
  printf("Vertex: (%.4f +- %.4f, %.4f +- %.4f, %.4f +- %.4f) cm\n",
	   fPrimaryVertex->GetXv(), fPrimaryVertex->GetXRes(),
	   fPrimaryVertex->GetYv(), fPrimaryVertex->GetYRes(),
	   fPrimaryVertex->GetZv(), fPrimaryVertex->GetZRes());
    printf("Mean vertex in RUN: X=%.4f Y=%.4f cm\n",
	   GetDiamondX(),GetDiamondY());
    printf("SPD Multiplicity. Number of tracklets %d \n",
           fSPDMult->GetNumberOfTracklets());
  printf("Number of tracks: \n");
  printf("                 charged   %d\n", GetNumberOfTracks());
  printf("                 muon      %d\n", GetNumberOfMuonTracks());
  printf("                 pmd       %d\n", GetNumberOfPmdTracks());
  printf("                 trd       %d\n", GetNumberOfTrdTracks());
  printf("                 v0        %d\n", GetNumberOfV0s());
  printf("                 cascades  %d\n", GetNumberOfCascades());
  printf("                 kinks     %d\n", GetNumberOfKinks());
  printf("                 CaloClusters %d\n", GetNumberOfCaloClusters());
  printf("                 phos      %d\n", GetNumberOfPHOSClusters());
  printf("                 emcal     %d\n", GetNumberOfEMCALClusters());
  printf("                 FMD       %s\n", (fESDFMD ? "yes" : "no"));
  printf("                 VZERO     %s\n", (fESDVZERO ? "yes" : "no"));
}

void AliESD::SetESDfriend(const AliESDfriend *ev) {
  //
  // Attaches the complementary info to the ESD
  //
  if (!ev) return;

  Int_t ntrk=ev->GetNumberOfTracks();
 
  for (Int_t i=0; i<ntrk; i++) {
    const AliESDfriendTrack *f=ev->GetTrack(i);
    GetTrack(i)->SetFriendTrack(f);
  }
}

void AliESD::GetESDfriend(AliESDfriend *ev) const {
  //
  // Extracts the complementary info from the ESD
  //
  if (!ev) return;

  Int_t ntrk=GetNumberOfTracks();

  for (Int_t i=0; i<ntrk; i++) {
    const AliESDtrack *t=GetTrack(i);
    const AliESDfriendTrack *f=t->GetFriendTrack();
    ev->AddTrack(f);
  }
}

void AliESD::AddObject(TObject* obj) 
{
  // Add an object to the list of object.
  // Please be aware that in order to increase performance you should
  // refrain from using TObjArrays (if possible). Use TClonesArrays, instead.
  fESDObjects->AddLast(obj);
}

void AliESD::GetStdContent() 
{
  // set pointers for standard content

  fESDRun = (AliESDRun*)fESDObjects->At(kESDRun);
  fHeader = (AliESDHeader*)fESDObjects->At(kHeader);
  fESDZDC = (AliESDZDC*)fESDObjects->At(kESDZDC);
  fESDFMD = (AliESDFMD*)fESDObjects->At(kESDFMD);
  fESDVZERO = (AliESDVZERO*)fESDObjects->At(kESDVZERO);
  fESDTZERO = (AliESDTZERO*)fESDObjects->At(kESDTZERO);
  fSPDVertex = (AliESDVertex*)fESDObjects->At(kSPDVertex);
  fPrimaryVertex = (AliESDVertex*)fESDObjects->At(kPrimaryVertex);
  fSPDMult =       (AliMultiplicity*)fESDObjects->At(kSPDMult);
  fPHOSTrigger = (AliESDCaloTrigger*)fESDObjects->At(kPHOSTrigger);
  fEMCALTrigger = (AliESDCaloTrigger*)fESDObjects->At(kEMCALTrigger);
  fTracks = (TClonesArray*)fESDObjects->At(kTracks);
  fMuonTracks = (TClonesArray*)fESDObjects->At(kMuonTracks);
  fPmdTracks = (TClonesArray*)fESDObjects->At(kPmdTracks);
  fTrdTracks = (TClonesArray*)fESDObjects->At(kTrdTracks);
  fV0s = (TClonesArray*)fESDObjects->At(kV0s);
  fCascades = (TClonesArray*)fESDObjects->At(kCascades);
  fKinks = (TClonesArray*)fESDObjects->At(kKinks);
  fCaloClusters = (TClonesArray*)fESDObjects->At(kCaloClusters);
  fErrorLogs = (TClonesArray*)fESDObjects->At(kErrorLogs);

}

void AliESD::SetStdNames(){

  fSPDVertex->SetName("SPDVertex");
  fPrimaryVertex->SetName("PrimaryVertex");
  fPHOSTrigger->SetName("PHOSTrigger");
  fEMCALTrigger->SetName("EMCALTrigger");
  fTracks->SetName("Tracks");
  fMuonTracks->SetName("MuonTracks");
  fPmdTracks->SetName("PmdTracks");
  fTrdTracks->SetName("TrdTracks");
  fV0s->SetName("V0s");
  fCascades->SetName("Cascades");
  fKinks->SetName("Kinks");
  fCaloClusters->SetName("CaloClusters");

} 

void AliESD::CreateStdContent() 
{
  // create the standard AOD content and set pointers

  // create standard objects and add them to the TList of objects
  AddObject(new AliESDRun());
  AddObject(new AliESDHeader());
  AddObject(new AliESDZDC());
  AddObject(new AliESDFMD());
  AddObject(new AliESDVZERO());
  AddObject(new AliESDTZERO());
  AddObject(new AliESDVertex());
  AddObject(new AliESDVertex());
  AddObject(new AliMultiplicity());
  AddObject(new AliESDCaloTrigger());
  AddObject(new AliESDCaloTrigger());
  AddObject(new TClonesArray("AliESDtrack",0));
  AddObject(new TClonesArray("AliESDMuonTrack",0));
  AddObject(new TClonesArray("AliESDPmdTrack",0));
  AddObject(new TClonesArray("AliESDTrdTrack",0));
  AddObject(new TClonesArray("AliESDv0",0));
  AddObject(new TClonesArray("AliESDcascade",0));
  AddObject(new TClonesArray("AliESDkink",0));
  AddObject(new TClonesArray("AliESDCaloCluster",0));
  AddObject(new TClonesArray("AliRawDataErrorLog",0));

  // check the order of the indices against enum...

  // read back pointers
  GetStdContent();
  // set names
  SetStdNames();

}

void AliESD::ReadFromTree(TTree *tree){
  

  // is this really so smart that an ESDObject has a pointer to a list
  // of another ESDObject...

  fESDObjects = (TList*)((AliESD*)tree->GetTree()->GetUserInfo()->FindObject("AliESD"))->GetList(); 

  if(fESDObjects->GetEntries()<kESDListN){
    printf("%s %d AliESD::ReadFromTree() TList contains less than the standard contents %d < %d \n",(char*)__FILE__,__LINE__,fESDObjects->GetEntries(),kESDListN);
  }

  // if list is empty
  // we could still set the branch adresses based on 
  // tree->GetListOfBranches() CKB
  // or create standard list 

  // set the branch addresses
  TIter next(fESDObjects);
  TNamed *el;
  while((el=(TNamed*)next())){
    TString bname(el->GetName());

    if(bname.CompareTo("AliESDfriend")==0)
      {
	// AliESDfriend does not have a name ...
      tree->SetBranchStatus("ESDfriend.*",1);
      printf("Friend %s\n", bname.Data());
      tree->SetBranchAddress("ESDfriend.",fESDObjects->GetObjectRef(el));


    }
    else{
      printf("%s\n", bname.Data());
      tree->SetBranchAddress(bname.Data(),fESDObjects->GetObjectRef(el));
    }
  }

  GetStdContent();
}



