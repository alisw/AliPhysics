#include "TMath.h"

#include "AliMFTCATrack.h"
#include "AliMFTCACell.h"

ClassImp(AliMFTCATrack)



//___________________________________________________________________________

AliMFTCATrack::AliMFTCATrack() :
TObject(),
fGID(-1),
fNcells(0),
fStartLayer(-1),
fCellGIDarray(),
fMCflag(0),
fVertX(0.),
fVertY(0.),
fVertZ(0.),
fTheta(0.),
fPhi(0.),
fChiSqX(0.),
fChiSqY(0.),
fMCindex(-1),
fChargeSign(0.)
{
  
  fCells = new TClonesArray("AliMFTCACell", fNDetMax);
  
}

//___________________________________________________________________________

AliMFTCATrack::AliMFTCATrack(const AliMFTCATrack &track) :
TObject(track),
fGID(track.fGID),
fNcells(track.fNcells),
fStartLayer(track.fStartLayer),
fMCflag(track.fMCflag),
fVertX(track.fVertX),
fVertY(track.fVertY),
fVertZ(track.fVertZ),
fTheta(track.fTheta),
fPhi(track.fPhi),
fChiSqX(track.fChiSqX),
fChiSqY(track.fChiSqY),
fMCindex(track.fMCindex),
fChargeSign(track.fChargeSign)
{
  
  // copy constructor
  
  fCells = new TClonesArray("AliMFTCACell", fNDetMax);
  
  for (Int_t icell = 0; icell < track.fNcells; icell++)
  fCellGIDarray[icell] = track.fCellGIDarray[icell];
  
}

//___________________________________________________________________________

void AliMFTCATrack::Clear(Option_t *) {
  
}

//___________________________________________________________________________

void AliMFTCATrack::AddCell(AliMFTCACell *cell) {
  
  if (fNcells >= fNDetMax) { printf("Max number of cells in this track!\n"); return; }
  new ((*fCells)[fNcells]) AliMFTCACell(*cell);
  fCellGIDarray[fNcells] = cell->GetGID();
  fNcells++;
  
}

//___________________________________________________________________________

Double_t AliMFTCATrack::AddCellToChiSq(AliMFTCACell *cell) {
  
  const Int_t nMaxh = 100;
  Double_t xTr[nMaxh], yTr[nMaxh], zTr[nMaxh];
  Double_t a, ae, b, be, x0, xS, y0, yS, chisqx, chisqy;
  Double_t xTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t yTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t xTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t yTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t xTrErr[nMaxh], yTrErr[nMaxh];
  for (Int_t i = 0; i < nMaxh; i++) {
    xTrErr[i] = TMath::Sqrt(xTrErrDet*xTrErrDet+xTrErrMS*xTrErrMS);
    yTrErr[i] = TMath::Sqrt(yTrErrDet*yTrErrDet+yTrErrMS*yTrErrMS);
  }
  Int_t cellGID, ndof, nptr = 0;
  
  AliMFTCACell *celltr;
  
  for (Int_t iCell = 0; iCell < GetNcells(); iCell++) {
    celltr = (AliMFTCACell*)fCells->At(iCell);
    // extract hit x,y,z
    if (nptr == 0) {
      xTr[nptr] = celltr->GetHit2()[0];
      yTr[nptr] = celltr->GetHit2()[1];
      zTr[nptr] = celltr->GetHit2()[2];
      nptr++;
      xTr[nptr] = celltr->GetHit1()[0];
      yTr[nptr] = celltr->GetHit1()[1];
      zTr[nptr] = celltr->GetHit1()[2];
      nptr++;
    } else {
      xTr[nptr] = celltr->GetHit1()[0];
      yTr[nptr] = celltr->GetHit1()[1];
      zTr[nptr] = celltr->GetHit1()[2];
      nptr++;
    }
  } // END : cells loop
    // the new cell
  xTr[nptr] = cell->GetHit1()[0];
  yTr[nptr] = cell->GetHit1()[1];
  zTr[nptr] = cell->GetHit1()[2];
  nptr++;
  // linear regression
//  if (LinFit(nptr,zTr,xTr,xTrErr,a,ae,b,be)) {
//    x0 = b; xS = a;
//    if (LinFit(nptr,zTr,yTr,yTrErr,a,ae,b,be)) {
//      y0 = b; yS = a;
//      chisqx = 0.;
//      chisqy = 0.;
//      for (Int_t iptr = 0; iptr < nptr; iptr++) {
//        //printf("%d  %f  %f  %f  \n",iptr,xTr[iptr],yTr[iptr],zTr[iptr]);
//        chisqx += (xTr[iptr]-(xS*zTr[iptr]+x0))*(xTr[iptr]-(xS*zTr[iptr]+x0))/(xTrErr[iptr]*xTrErr[iptr]);
//        chisqy += (yTr[iptr]-(yS*zTr[iptr]+y0))*(yTr[iptr]-(yS*zTr[iptr]+y0))/(yTrErr[iptr]*yTrErr[iptr]);
//      }
//    }
//  }
  ndof = 2*nptr-4;
  
  return (chisqx+chisqy)/(Double_t)ndof;
  
}
