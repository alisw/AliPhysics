
//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//
//    UPC event containing trigger, VZERO, AD, ZDC, tracks, MC and other
//    arbitrary data in TArrayI and TArrayD containters
//_____________________________________________________________________________

#include "TH1I.h"
#include "TObjString.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TBits.h"

#include "AliUPCTrack.h"
#include "AliUPCMuonTrack.h"

#include "AliUPCEvent.h"

ClassImp(AliUPCEvent)

TClonesArray *AliUPCEvent::fgUPCTracks = 0;
TClonesArray *AliUPCEvent::fgUPCMuonTracks = 0;
TClonesArray *AliUPCEvent::fgMCParticles = 0;

//_____________________________________________________________________________
AliUPCEvent::AliUPCEvent()
 :TObject(),
  fFlags(0), fL0inputs(0), fRecoPass(0),
  fDataFilnam(0x0), fEvtNum(0), fRunNum(0),
  fVtxChi2perNDF(0), fVtxNContributors(0), fVtxTitle(0x0),
  fVtxSPDchi2perNDF(0), fVtxSPDnContributors(0), fVtxSPDtitle(0x0),
  fNTracklets(0), fNSPDfiredInner(0), fNSPDfiredOuter(0), fFOmap(0x0),
  fV0ADecision(0), fV0CDecision(0), fBBtriggerV0C(0), fBBFlagV0C(0), fBBtriggerV0A(0), fBBFlagV0A(0),
   fBBtriggerADC(0), fBBFlagADC(0), fBBtriggerADA(0), fBBFlagADA(0), fADADecision(0), fADCDecision(0),
  fZNCEnergy(0), fZPCEnergy(0), fZNAEnergy(0), fZPAEnergy(0),
  fZNCtdc(0), fZPCtdc(0), fZNAtdc(0), fZPAtdc(0),
  fZNCtdcData(0), fZPCtdcData(0), fZNAtdcData(0), fZPAtdcData(0), fZNCTime(0), fZNATime(0),fBCnumber(0),
  fUPCTracks(0x0), fNtracks(0), fUPCMuonTracks(0x0), fNmuons(0),
  fMCParticles(0x0), fNmc(0),
  fArrayInt(0x0), fArrayD(0x0)
{
  // Default constructor

  for(Int_t itrg=0; itrg<fgkNtrg; itrg++) fTrgClasses[itrg] = kFALSE;
  for(Int_t i=0; i<3; i++) {fVtxPos[i] = 0.; fVtxSPDpos[i] = 0.; fVtxMCpos[i] = 0.;}
  for(Int_t i=0; i<6; i++) {fVtxCov[i] = 0.; fVtxSPDcov[i] = 0.;}
  for(Int_t i=0; i<4; i++) {fZNATDCm[i] = -999; fZNCTDCm[i] = -999; fZPATDCm[i] = -999; fZPCTDCm[i] = -999;}

  if(!fgUPCTracks) {
    fgUPCTracks = new TClonesArray("AliUPCTrack");
    fUPCTracks = fgUPCTracks;
    fUPCTracks->SetOwner(kTRUE);
  }
  if(!fgUPCMuonTracks) {
    fgUPCMuonTracks = new TClonesArray("AliUPCMuonTrack");
    fUPCMuonTracks = fgUPCMuonTracks;
    fUPCMuonTracks->SetOwner(kTRUE);
  }

  if(!fDataFilnam) {
    fDataFilnam = new TObjString();
    fDataFilnam->SetString("");
  }

  if(!fVtxTitle) {
    fVtxTitle = new TObjString();
    fVtxTitle->SetString("");
  }

  if(!fVtxSPDtitle) {
    fVtxSPDtitle = new TObjString();
    fVtxSPDtitle->SetString("");
  }

  //if(!fUPCTracks) fUPCTracks = new TClonesArray("AliUPCTrack");
  //if(!fUPCMuonTracks) fUPCMuonTracks = new TClonesArray("AliUPCMuonTrack");
}

//_____________________________________________________________________________
AliUPCEvent::~AliUPCEvent()
{
  // destructor

  if(fDataFilnam) {delete fDataFilnam; fDataFilnam = 0x0;}
  if(fVtxTitle) {delete fVtxTitle; fVtxTitle=0x0;}
  if(fVtxSPDtitle) {delete fVtxSPDtitle; fVtxSPDtitle=0x0;}
  if(fUPCTracks) {delete fUPCTracks; fUPCTracks = 0x0;}
  if(fUPCMuonTracks) {delete fUPCMuonTracks; fUPCMuonTracks = 0x0;}
  if(fMCParticles) {delete fMCParticles; fMCParticles = 0x0;}
  if(fArrayInt) {delete fArrayInt; fArrayInt = 0x0;}
  if(fArrayD) {delete fArrayD; fArrayD = 0x0;}
}

//_____________________________________________________________________________
void AliUPCEvent::ClearEvent()
{
  // clear event variables

  Clear(); // Clear to TObject

  fBBtriggerV0C = 0;
  fBBFlagV0C = 0;
  fBBtriggerV0A = 0;
  fBBFlagV0A = 0;
  fBBtriggerADC = 0;
  fBBFlagADC = 0;
  fBBtriggerADA = 0;
  fBBFlagADA = 0;
  for(Int_t itrg=0; itrg<fgkNtrg; itrg++) fTrgClasses[itrg] = kFALSE;

  if(fUPCTracks) {fUPCTracks->Clear("C"); fNtracks = 0;}
  if(fUPCMuonTracks) {fUPCMuonTracks->Clear("C"); fNmuons = 0;}
  if(fMCParticles) {fMCParticles->Clear("C"); fNmc = 0;}
  if(fArrayInt) fArrayInt->Reset();
  if(fArrayD) fArrayD->Reset();
}

//_____________________________________________________________________________
void AliUPCEvent::SetIsESD(Bool_t isESD)
{
  //set event as created from ESD, otherwise AOD assumed

  if(!isESD) return;

  fFlags |= (1 << 0);

}

//_____________________________________________________________________________
void AliUPCEvent::SetIsMC(Bool_t isMC)
{
  // set the event as MC, initialize mc array

  if(!isMC) return;

  if(!fgMCParticles) {
    fgMCParticles = new TClonesArray("TParticle");
    fMCParticles = fgMCParticles;
    fMCParticles->SetOwner(kTRUE);
  }
  fFlags |= (1 << 1);

}

//_____________________________________________________________________________
void AliUPCEvent::SetTriggerClass(Int_t idx, Bool_t fired)
{
  // set trigger class at index idx

  if(idx < 0 || idx >= fgkNtrg) return;
  fTrgClasses[idx] = fired;
}

//_____________________________________________________________________________
void AliUPCEvent::SetInputFileName(const char *s)
{
  // set name of the input file

  fDataFilnam->Clear();
  fDataFilnam->SetString(s);

}

//_____________________________________________________________________________
void AliUPCEvent::SetPrimaryVertex(Double_t *pos, Double_t chi2, Double_t *cov, Int_t ncontrib)
{
  // set default primary vertex

  for(Int_t i=0; i<3; i++) fVtxPos[i] = pos[i];
  for(Int_t i=0; i<6; i++) fVtxCov[i] = cov[i];
  fVtxChi2perNDF = chi2;
  fVtxNContributors = ncontrib;
}

//_____________________________________________________________________________
void AliUPCEvent::SetPrimaryVertexTitle(const char *s)
{
  // set title of default primary vertex

  fVtxTitle->Clear();
  fVtxTitle->SetString(s);
}

//_____________________________________________________________________________
void AliUPCEvent::SetPrimaryVertexSPD(Double_t *pos, Double_t chi2, Double_t *cov, Int_t ncontrib)
{
  // set SPD primary vertex

  for(Int_t i=0; i<3; i++) fVtxSPDpos[i] = pos[i];
  for(Int_t i=0; i<6; i++) fVtxSPDcov[i] = cov[i];
  fVtxSPDchi2perNDF = chi2;
  fVtxSPDnContributors = ncontrib;
}

//_____________________________________________________________________________
void AliUPCEvent::SetPrimaryVertexSPDtitle(const char *s)
{
  // set title of default primary vertex

  fVtxSPDtitle->Clear();
  fVtxSPDtitle->SetString(s);
}

//_____________________________________________________________________________
void AliUPCEvent::SetPrimaryVertexMC(Float_t *pos)
{
  // set generated primary vertex

  for(Int_t i=0; i<3; i++) fVtxMCpos[i] = pos[i];
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBtriggerV0Cmask(UInt_t ibit)
{
  // fill offline fired cells in V0C setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBtriggerV0C |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBFlagV0Cmask(UInt_t ibit)
{
  // fill online fired cells in V0C setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBFlagV0C |= (1 << ibit);
}
//_____________________________________________________________________________
void AliUPCEvent::SetBBtriggerV0Amask(UInt_t ibit)
{
  // fill offline fired cells in V0C setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBtriggerV0A |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBFlagV0Amask(UInt_t ibit)
{
  // fill online fired cells in V0C setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBFlagV0A |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBtriggerADCmask(UInt_t ibit)
{
  // fill offline fired cells in ADC setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBtriggerADC |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBFlagADCmask(UInt_t ibit)
{
  // fill online fired cells in ADC setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBFlagADC |= (1 << ibit);
}
//_____________________________________________________________________________
void AliUPCEvent::SetBBtriggerADAmask(UInt_t ibit)
{
  // fill offline fired cells in ADC setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBtriggerADA |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetBBFlagADAmask(UInt_t ibit)
{
  // fill online fired cells in ADC setting the corresponding bits in the bit mask
  // AliUPCEvent::ClearEvent should be called to initialize the mask

  fBBFlagADA |= (1 << ibit);
}

//_____________________________________________________________________________
void AliUPCEvent::SetZNTDCm(Float_t *znatdcm,Float_t *znctdcm,Float_t *zpatdcm,Float_t *zpctdcm)
{
  for(Int_t i=0; i<4; i++) fZNATDCm[i] = znatdcm[i];
  for(Int_t i=0; i<4; i++) fZNCTDCm[i] = znctdcm[i];
  for(Int_t i=0; i<4; i++) fZPATDCm[i] = zpatdcm[i];
  for(Int_t i=0; i<4; i++) fZPCTDCm[i] = zpctdcm[i];

}


//_____________________________________________________________________________
AliUPCTrack *AliUPCEvent::AddTrack(void)
{
  // construct new upc central track
/*
  if(!fgUPCTracks) {
    fgUPCTracks = new TClonesArray("AliUPCTrack");
    fUPCTracks = fgUPCTracks;
    fUPCTracks->SetOwner(kTRUE);
  }
*/
  if(!fUPCTracks) return 0x0;

  AliUPCTrack *track = (AliUPCTrack*) fUPCTracks->ConstructedAt(fNtracks++);
  return track;
}

//_____________________________________________________________________________
AliUPCMuonTrack *AliUPCEvent::AddMuonTrack(void)
{
  // construct new upc muon track

  if(!fUPCMuonTracks) return 0x0;

  AliUPCMuonTrack *track = (AliUPCMuonTrack*) fUPCMuonTracks->ConstructedAt(fNmuons++);
  return track;
}

//_____________________________________________________________________________
TParticle *AliUPCEvent::AddMCParticle(void)
{
  // construct new mc TParticle

  if(!fMCParticles) return 0x0;

  TParticle *part = (TParticle*) fMCParticles->ConstructedAt(fNmc++);
  return part;
}

//_____________________________________________________________________________
Bool_t AliUPCEvent::GetIsESD(void) const
{
  // return kTRUE if event is ESD, otherwise AOD assumed

  if( fFlags & (1 << 0) ) return kTRUE;

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliUPCEvent::GetIsMC(void) const
{
  // return kTRUE if event is MC

  if( fFlags & (1 << 1) ) return kTRUE;

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliUPCEvent::GetFlagBit(UChar_t ibit) const
{
  // return kTRUE if bit is set

  if( fFlags & (1 << ibit) ) return kTRUE;

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliUPCEvent::GetTriggerClass(Int_t idx) const
{
  // returns kTRUE if trigger class at idx was fired, otherwise kFALSE;

  if(idx < 0 || idx >= fgkNtrg) return kFALSE;
  return fTrgClasses[idx];
}

//_____________________________________________________________________________
void AliUPCEvent::GetPrimaryVertex(Double_t *pos, Double_t &chi2, Double_t *cov, Int_t &ncontrib) const
{
  // get default primary vertex

  for(Int_t i=0; i<3; i++) pos[i] = fVtxPos[i];
  for(Int_t i=0; i<6; i++) cov[i] = fVtxCov[i];
  chi2 = fVtxChi2perNDF;
  ncontrib = fVtxNContributors;
}

//_____________________________________________________________________________
void AliUPCEvent::GetPrimaryVertexSPD(Double_t *pos, Double_t &chi2, Double_t *cov, Int_t &ncontrib) const
{
  // get SPD primary vertex

  for(Int_t i=0; i<3; i++) pos[i] = fVtxSPDpos[i];
  for(Int_t i=0; i<6; i++) cov[i] = fVtxSPDcov[i];
  chi2 = fVtxSPDchi2perNDF;
  ncontrib = fVtxSPDnContributors;
}

//_____________________________________________________________________________
void AliUPCEvent::GetPrimaryVertexMC(Float_t *pos) const
{
  // get generated primary vertex

  for(Int_t i=0; i<3; i++) pos[i] = fVtxMCpos[i];
}

//_____________________________________________________________________________
Bool_t AliUPCEvent::Get0SMB(void) const
{
  // get trigger input 0SMB

  if( fNSPDfiredInner + fNSPDfiredOuter >= 1 ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliUPCEvent::Get0SH1(void) const
{
  // get trigger input 0SH1

  if( fNSPDfiredOuter >= 7 ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliUPCEvent::Get0STP(void) const
{
  // get trigger input 0STP
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/0STPTrigger

  if(!fFOmap) return kFALSE;

  Int_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=0;
  Int_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=0;

  for (Int_t i=0; i<400; ++i)
    if (fFOmap->TestBitNumber(i))
      ++vPhiInner[i/20];

  for (Int_t i=400; i<1200; ++i)
    if (fFOmap->TestBitNumber(i))
      ++vPhiOuter[(i-400)/20];

  Int_t fired(0);
  for (Int_t i=0; i<10; ++i) {
    for (Int_t j=0; j<2; ++j) {
      const Int_t k(2*i+j);
      fired += (   (vPhiOuter[k]    || vPhiOuter[k+1]       || vPhiOuter[k+2])
                && (vPhiOuter[k+20] || vPhiOuter[(k+21)%40] || vPhiOuter[(k+22)%40])
                && (vPhiInner[i]    || vPhiInner[i+1])
                && (vPhiInner[i+10] || vPhiInner[(i+11)%20]));
    }
  }
  return fired;

}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNV0ChitsOffline(void) const
{
  // get number of fired cells of V0C, offline

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<32; iv++) {if( fBBtriggerV0C & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNV0ChitsOnline(void) const
{
  // get number of fired cells of V0C, online

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<32; iv++) {if( fBBFlagV0C & (1 << iv) ) nHits++;}

  return nHits;
}
//_____________________________________________________________________________
Int_t AliUPCEvent::GetNV0AhitsOffline(void) const
{
  // get number of fired cells of V0C, offline

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<32; iv++) {if( fBBtriggerV0A & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNV0AhitsOnline(void) const
{
  // get number of fired cells of V0C, online

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<32; iv++) {if( fBBFlagV0A & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNADChitsOffline(void) const
{
  // get number of fired cells of ADC, offline

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<8; iv++) {if( fBBtriggerADC & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNADChitsOnline(void) const
{
  // get number of fired cells of ADC, online

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<8; iv++) {if( fBBFlagADC & (1 << iv) ) nHits++;}

  return nHits;
}
//_____________________________________________________________________________
Int_t AliUPCEvent::GetNADAhitsOffline(void) const
{
  // get number of fired cells of ADC, offline

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<8; iv++) {if( fBBtriggerADA & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::GetNADAhitsOnline(void) const
{
  // get number of fired cells of ADC, online

  Int_t nHits = 0;
  for (UInt_t iv=0; iv<8; iv++) {if( fBBFlagADA & (1 << iv) ) nHits++;}

  return nHits;
}

//_____________________________________________________________________________
Bool_t AliUPCEvent::GetAllZDCtdc(void) const
{
  // get the final tdc of all ZDCs, return kTRUE if at least one tdc is set,
  // return kFALSE only if all tdc are empty

  if( fZNCtdc || fZPCtdc || fZNAtdc || fZPAtdc ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
AliUPCTrack *AliUPCEvent::GetTrack(Int_t iTrack) const
{
  // get upc central track

  if(!fUPCTracks) return 0x0;

  return (AliUPCTrack*) fUPCTracks->At(iTrack);
}

//_____________________________________________________________________________
AliUPCMuonTrack *AliUPCEvent::GetMuonTrack(Int_t iTrack) const
{
  // get upc muon track

  if(!fUPCMuonTracks) return 0x0;

  return (AliUPCMuonTrack*) fUPCMuonTracks->At(iTrack);
}

//_____________________________________________________________________________
TParticle *AliUPCEvent::GetMCParticle(Int_t iMC) const
{
  // get mc TParticle

  if(!fMCParticles) return 0x0;

  return (TParticle*) fMCParticles->At(iMC);
}

//_____________________________________________________________________________
Int_t AliUPCEvent::MakeArrayInt(Int_t size)
{
  //create the object of TArrayI, skip if it has been already created
  //the array allows the event to be extended for other integer variables

  if( fArrayInt ) return 999; // already initialized

  fArrayInt = new TArrayI(size);

  return 0;
}

//_____________________________________________________________________________
Int_t AliUPCEvent::MakeArrayD(Int_t size)
{
  //create the object of TArrayD, skip if it has been already created
  //the array allows the event to be extended for other double variables

  if( fArrayD ) return 999; // already initialized

  fArrayD = new TArrayD(size);

  return 0;
}

/*
//_____________________________________________________________________________
AliUPCEvent::AliUPCEvent(const AliUPCEvent &o)
 :TObject(o),
  fNTracklets(o.fNTracklets)
{

  // Copy constructor

}

//_____________________________________________________________________________
AliUPCEvent &AliUPCEvent::operator=(const AliUPCEvent &o)
{

  if(this==&o) return *this;
  AliUPCEvent::operator=(o);

  // Assignment operator
  fNTracklets = o.fNTracklets;

  return *this;

}
*/











































