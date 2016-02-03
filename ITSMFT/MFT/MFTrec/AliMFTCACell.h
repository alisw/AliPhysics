#ifndef AliMFTCACell_H
#define AliMFTCACell_H

#include "TObject.h"
#include "TVector3.h"

class AliMFTCACell : public TObject {
  
public:
  
  AliMFTCACell();
  ~AliMFTCACell() {};
  
  AliMFTCACell (const AliMFTCACell &cell);
  AliMFTCACell &operator=(const AliMFTCACell&);
  
  virtual void Clear(Option_t *);
  void SetGID(Int_t gid, Int_t trackid1, Int_t trackid2) {
    fGID = gid; fTrackGID[0] = trackid1; fTrackGID[1] = trackid2; }
  void SetStatus(Int_t s) { fStatus = s; }
  void SetHits(Double_t *h1, Double_t *h2, Double_t z1, Double_t z2);
  void SetLayers(Int_t iL1, Int_t iL2) { fLayer[0] = iL1; fLayer[1] = iL2; }
  void SetMFTClsId(Int_t id1, Int_t id2) { fMFTClsId[0] = id1; fMFTClsId[1] = id2; }
  void SetDetElemID(Int_t id1, Int_t id2) { fDetElemID[0] = id1; fDetElemID[1] = id2; }
  Double_t *GetHit1() { return fHit[0]; }
  Double_t *GetHit2() { return fHit[1]; }
  Double_t *GetHitp1() { return fHitp[0]; }
  Double_t *GetHitp2() { return fHitp[1]; }
  Int_t *GetLayers() { return fLayer; }
  Int_t *GetMFTClsId() { return fMFTClsId; }
  Int_t *GetDetElemID() { return fDetElemID; }
  TVector3 *GetSeg()  { return fSeg; }
  const Bool_t HasNbL() { return fNNbL > 0 ? kTRUE : kFALSE; }
  const Bool_t HasNbR() { return fNNbR > 0 ? kTRUE : kFALSE; }
  void UpdateStatus() { if (fUpStat) { fStatus++; fUpStat = kFALSE; } }
  void IncrStatus() { fUpStat = kTRUE; }
  const Int_t GetStatus() { return fStatus; }
  const Int_t GetGID() { return fGID; }
  const Int_t GetTrackGID(Int_t p) {
    if (p != 0 && p != 1) {
      printf("AliMFTCACell::GetTrackGID : wrong index.\n");
      return -1;
    } else return fTrackGID[p]; }
  void DrawOGL(Char_t *name);
  void SetUsed(Bool_t isu) { fIsUsed = isu; }
  const Bool_t IsUsed() { return fIsUsed; }
  void AddLeftNeighbour(Int_t cellgid)  { fNbLgid[fNNbL++] = cellgid; }
  void AddRightNeighbour(Int_t cellgid) { fNbRgid[fNNbR++] = cellgid; }
  const Int_t GetNNbL() { return fNNbL; }
  const Int_t GetNNbR() { return fNNbR; }
  const Int_t GetNbLgid(Int_t id) { return fNbLgid[id]; }
  const Int_t GetNbRgid(Int_t id) { return fNbRgid[id]; }
  void Reset() {
    fStatus = 1;
    fUpStat = kFALSE;
    fNNbL = fNNbR = 0;
    for (Int_t i = 0; i < 100; i++) { fNbLgid[i] = fNbRgid[i] = -1; }
  }
  const Int_t GetLength() {
    return (fLayer[1]-fLayer[0]);
  }
  void SetIsMerged() { fIsMerged = kTRUE; }
  const Bool_t IsMerged() { return fIsMerged; }
  
  void PrintCell(Option_t *opt) {
    printf("PrintCell------------------------GID = %10d \n",fGID);
    printf("Hit1: %9.4f   %9.4f   %9.4f \n",fHit[0][0],fHit[0][1],fHit[0][2]);
    printf("Hit2: %9.4f   %9.4f   %9.4f \n",fHit[1][0],fHit[1][1],fHit[1][2]);
    printf("Status: %d \n",fStatus);
    if (strcmp(opt,"FULL") == 0) {
      printf("Segment vector: \n");
      fSeg->Print();
    }
    if (strcmp(opt,"MC") == 0) {
      printf("Track:  %5d   %5d \n",fTrackGID[0],fTrackGID[1]);
    }
  }
  
private:
  
  Int_t fGID;              // Global identifier
  Int_t fTrackGID[2];      // Track GID of cell ends
  Int_t fLayer[2];         // Layer numbers of cell ends
  Double_t fHit[2][3];     // X,Y,Z values of cell ends at plane position
  Double_t fHitp[2][3];    // X,Y,Z values of cell ends
  TVector3 *fSeg;          //! Cell segment
  Int_t fStatus;           // Cell status
  Bool_t fUpStat;          // Status must be updated
  Bool_t fIsUsed;          // True if the cell has been attached to a track
  Int_t fNNbL;             // Number of neighbours at left
  Int_t fNNbR;             // Number of neighbours at right
  Int_t fNbLgid[100];      // GID of cell left neighbours
  Int_t fNbRgid[100];      // GID of cell right neighbours
  Bool_t fIsMerged;        // True if is from merged cells in the overlaps
  Int_t fDetElemID[2];     // ladders ID
  Int_t fMFTClsId[2];      // ID of MFT clusters  

  ClassDef(AliMFTCACell,1);
  
};

#endif
