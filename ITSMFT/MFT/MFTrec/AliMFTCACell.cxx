
#include "AliMFTCACell.h"

ClassImp(AliMFTCACell)

//___________________________________________________________________________
AliMFTCACell::AliMFTCACell() :
TObject(),
fGID(-1),
fTrackGID(),
fLayer(),
fHit(),
fHitp(),
fSeg(),
fStatus(0),
fUpStat(kFALSE),
fIsUsed(kFALSE),
fNNbL(0),
fNNbR(0),
fIsMerged(kFALSE),
fDetElemID()
{
  
  for (Int_t i = 0; i < 100; i++) { fNbLgid[i] = fNbRgid[i] = -1; }
  for (Int_t i = 0; i < 2; i++) {
    fTrackGID[i] = -1;
    fLayer[i] = -1;
    fDetElemID[i] = -1;
    for (Int_t j = 0; j < 3; j++) {
      fHit[i][j] = 0.;
      fHitp[i][j] = 0.;
    }
  }
  fMFTClsId[0] = fMFTClsId[1] = -1;
  
}

//___________________________________________________________________________
AliMFTCACell::AliMFTCACell(const AliMFTCACell &cell) :
TObject(cell),
fGID(cell.fGID),
fTrackGID(),
fLayer(),
fHit(),
fHitp(),
fSeg(),
fStatus(cell.fStatus),
fUpStat(cell.fUpStat),
fIsUsed(cell.fIsUsed),
fNNbL(cell.fNNbL),
fNNbR(cell.fNNbR),
fIsMerged(cell.fIsMerged),
fDetElemID()
{
  
  // copy constructor
  
  for (Int_t i = 0; i < cell.fNNbL; i++) fNbLgid[i] = cell.fNbLgid[i];
  for (Int_t i = 0; i < cell.fNNbR; i++) fNbRgid[i] = cell.fNbRgid[i];
  for (Int_t i = 0; i < 2; i++) {
    fTrackGID[i] = cell.fTrackGID[i];
    fLayer[i] = cell.fLayer[i];
    fDetElemID[i] = cell.fDetElemID[i];
    for (Int_t j = 0; j < 3; j++) {
      fHit[i][j] = cell.fHit[i][j];
      fHitp[i][j] = cell.fHitp[i][j];
    }
  }
  fSeg = new TVector3(fHit[1][0]-fHit[0][0],
                      fHit[1][1]-fHit[0][1],
                      fHit[1][2]-fHit[0][2]);
  fMFTClsId[0] = cell.fMFTClsId[0];  
  fMFTClsId[1] = cell.fMFTClsId[1];  

}

//___________________________________________________________________________
AliMFTCACell& AliMFTCACell::operator=(const AliMFTCACell& cell) 
{

  // assignment operator

  // check assignement to self
  if (this == &cell) return *this;

  TObject::operator=(cell);

  fGID = cell.fGID;
  for (Int_t i = 0; i < cell.fNNbL; i++) fNbLgid[i] = cell.fNbLgid[i];
  for (Int_t i = 0; i < cell.fNNbR; i++) fNbRgid[i] = cell.fNbRgid[i];
  for (Int_t i = 0; i < 2; i++) {
    fTrackGID[i] = cell.fTrackGID[i];
    fLayer[i] = cell.fLayer[i];
    fDetElemID[i] = cell.fDetElemID[i];
    for (Int_t j = 0; j < 3; j++) {
      fHit[i][j] = cell.fHit[i][j];
      fHitp[i][j] = cell.fHitp[i][j];
    }
  }
  fSeg = new TVector3(fHit[1][0]-fHit[0][0],
                      fHit[1][1]-fHit[0][1],
                      fHit[1][2]-fHit[0][2]);
  fStatus = cell.fStatus;
  fUpStat = cell.fUpStat;
  fIsUsed = cell.fIsUsed;
  fNNbL = cell.fNNbL;
  fNNbR = cell.fNNbR;
  fIsMerged = cell.fIsMerged;
  fMFTClsId[0] = cell.fMFTClsId[0];  
  fMFTClsId[1] = cell.fMFTClsId[1];  
 
}

//___________________________________________________________________________
void AliMFTCACell::Clear(const Option_t *) {
  
  delete fSeg;
  
  fGID = -1;
  fTrackGID[0]  = fTrackGID[1]  = -1;
  fLayer[0]     = fLayer[1]     = -1;
  fDetElemID[0] = fDetElemID[1] = -1;
  fHit[0][0] = fHit[0][1] = fHit[0][2] = 0.;
  fHit[1][0] = fHit[1][1] = fHit[1][2] = 0.;
  fHitp[0][0] = fHitp[0][1] = fHitp[0][2] = 0.;
  fHitp[1][0] = fHitp[1][1] = fHitp[1][2] = 0.;
  fStatus = 0;
  fUpStat = fIsUsed = kFALSE;
  fNNbL = fNNbR = 0;
  for (Int_t i = 0; i < 100; i++) { fNbLgid[i] = fNbRgid[i] = -1; }
  fIsMerged = kFALSE;
  
}

//___________________________________________________________________________
void AliMFTCACell::SetHits(Double_t *h1, Double_t *h2, Double_t z1, Double_t z2) {
  
  for (Int_t i = 0; i < 3; i++) {
    fHitp[0][i] = h1[i];
    fHitp[1][i] = h2[i];
    fHit[0][i] = h1[i];
    fHit[1][i] = h2[i];
  }
  /*
   // re-calculate hit projection in the median plane of a station
   fHit[0][0] = fHitp[0][0]+(fHitp[1][0]-fHitp[0][0])/(fHitp[1][2]-fHitp[0][2])*(z1-fHitp[0][2]);
   fHit[0][1] = fHitp[0][1]+(fHitp[1][1]-fHitp[0][1])/(fHitp[1][2]-fHitp[0][2])*(z1-fHitp[0][2]);
   fHit[0][2] = z1;
   fHit[1][0] = fHitp[0][0]+(fHitp[1][0]-fHitp[0][0])/(fHitp[1][2]-fHitp[0][2])*(z2-fHitp[0][2]);
   fHit[1][1] = fHitp[0][1]+(fHitp[1][1]-fHitp[0][1])/(fHitp[1][2]-fHitp[0][2])*(z2-fHitp[0][2]);
   fHit[1][2] = z2;
   */
  //fSeg = new TVector3(h2[0]-h1[0],h2[1]-h1[1],h2[2]-h1[2]);
  fSeg = new TVector3(fHit[1][0]-fHit[0][0],
                      fHit[1][1]-fHit[0][1],
                      fHit[1][2]-fHit[0][2]);
  
}

//___________________________________________________________________________
void AliMFTCACell::DrawOGL(Char_t *name) {
  
//  TEveLine *track = new TEveLine(2);
//  track->SetName(Form("%s_ID%d",name,GetGID()));
//  track->SetLineColor(1);
//  track->SetPoint(0,fHit[0][0],fHit[0][1],fHit[0][2]);
//  track->SetPoint(1,fHit[1][0],fHit[1][1],fHit[1][2]);
//  gEve->AddElement(track);
  
}
