#ifndef AliMFTCATrack_H
#define AliMFTCATrack_H

#include "TObject.h"
#include "TClonesArray.h"

#include "AliMFTConstants.h"

class AliMFTCACell;

//_________________________________________________________________________________
class AliMFTCATrack : public TObject {
  
public:
  
  AliMFTCATrack();
  ~AliMFTCATrack() {};
  
  AliMFTCATrack (const AliMFTCATrack &track);
  AliMFTCATrack &operator=(const AliMFTCATrack&);
  
  virtual void Clear(Option_t *);
  void AddCell(AliMFTCACell *cell);
  void SetGID(Int_t gid) { fGID = gid; }
  const Int_t GetGID() const { return fGID; }
  const Int_t GetNcells() const { return fNcells; }
  const Int_t GetCellGID(Int_t ic) const { return fCellGIDarray[ic]; }
  AliMFTCACell* GetCell(Int_t ic) const { return (AliMFTCACell*)fCells->At(ic); }
  const Int_t GetLastCellGID() const { return fCellGIDarray[fNcells-1]; }
  const Int_t GetStartLayer() const { return fStartLayer; }
  void SetStartLayer(Int_t sl) { fStartLayer = sl; }
  void SetMCflag(UChar_t mcf) { fMCflag = mcf; }
  const UChar_t GetMCflag() { return fMCflag; }
  void SetVertX(Double_t x) { fVertX = x; }
  void SetVertY(Double_t y) { fVertY = y; }
  void SetVertZ(Double_t z) { fVertZ = z; }
  void SetTheta(Double_t the) { fTheta = the; }
  void SetPhi(Double_t phi) { fPhi = phi; }
  const Double_t GetVertX() { return fVertX; }
  const Double_t GetVertY() { return fVertY; }
  const Double_t GetVertZ() { return fVertZ; }
  const Double_t GetTheta() { return fTheta; }
  const Double_t GetPhi() { return fPhi; }
  void SetChiSqX(Double_t chisq) { fChiSqX = chisq; }
  void SetChiSqY(Double_t chisq) { fChiSqY = chisq; }
  const Double_t GetChiSqX() { return fChiSqX; }
  const Double_t GetChiSqY() { return fChiSqY; }
  Double_t AddCellToChiSq(AliMFTCACell *cell);
  void SetMCindex(Int_t index) { fMCindex = index; }
  const Int_t GetMCindex() { return fMCindex; }
  void SetChargeSign(Short_t sign) { fChargeSign = sign; }
  const Short_t GetChargeSign() { return fChargeSign; }
  void SetCellGID(Int_t index, Int_t gid) { fCellGIDarray[index] = gid; };
	void EvalSignedPt();
	const Double_t GetPt() { return fPt; }

private:
  static const Int_t fNDetMax = AliMFTConstants::fNMaxPlanes;

  Int_t  fGID;                    // Track global identifier
  Int_t  fNcells;                 // Number of cells in the track
  Int_t  fStartLayer;             // Upstream start layer (RunBackward)
  Int_t  fCellGIDarray[fNDetMax]; // Array of cell global identifier
  
  TClonesArray *fCells;           //! Array of cells
  
  UChar_t fMCflag;                // MC classification of the track:
                                  // 0 = default (not set)
                                  // 1 = clean
                                  // 2 = good
                                  // 3 = fake
                                  // 4 = noise
  
  Double_t fVertX;                // x at z vertex [cm]
  Double_t fVertY;                // y at z vertex [cm]
  Double_t fVertZ;                // z vertex [cm]
  Double_t fTheta;                // theta fit [deg]
  Double_t fPhi;                  // phi fit [deg]
  Double_t fChiSqX;               // reduced ChiSq en xz
  Double_t fChiSqY;               // reduced ChiSq en yz
  Int_t    fMCindex;              // MC track index for clean tracks
  Short_t  fChargeSign;           // estimated sign of the charge
	Double_t fPt  ;									// Pt evaluated from sagitta measurement

  ClassDef(AliMFTCATrack,2);
  
};

#endif


