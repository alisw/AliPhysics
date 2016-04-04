#include "TMath.h"
#include <TGeoGlobalMagField.h>

#include "AliMFTCATrack.h"
#include "AliMFTCACell.h"
#include "AliMFTTrackExtrap.h"

#include "AliLog.h"

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
fChargeSign(0.),
fPt(0.)
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
fChargeSign(track.fChargeSign),
fPt(track.fPt)
{
  
  // copy constructor
  
  fCells = new TClonesArray("AliMFTCACell", fNDetMax);
  
  for (Int_t icell = 0; icell < track.fNcells; icell++)
  fCellGIDarray[icell] = track.fCellGIDarray[icell];
  
}

//___________________________________________________________________________

AliMFTCATrack& AliMFTCATrack::operator=(const AliMFTCATrack& track) 
{

  // assignment operator

  // check assignement to self
  if (this == &track) return *this;

  TObject::operator=(track);

  fGID = track.fGID;
  fNcells = track.fNcells;
  fStartLayer = track.fStartLayer;
  fCells = new TClonesArray("AliMFTCACell", fNDetMax);
  for (Int_t icell = 0; icell < track.fNcells; icell++)
  fCellGIDarray[icell] = track.fCellGIDarray[icell];
  fMCflag = track.fMCflag;
  fVertX = track.fVertX;
  fVertY = track.fVertY;
  fVertZ = track.fVertZ;
  fTheta = track.fTheta;
  fPhi = track.fPhi;
  fChiSqX = track.fChiSqX;
  fChiSqY = track.fChiSqY;
  fMCindex = track.fMCindex;
  fChargeSign = track.fChargeSign;
  fPt = track.fPt;

}

//___________________________________________________________________________

void AliMFTCATrack::Clear(Option_t *) {
  
}
//___________________________________________________________________________

/// Estimate the charge sign
void  AliMFTCATrack::EvalSignedPt(){
	const Int_t nMaxh = 100;
	Double_t xTr[nMaxh], yTr[nMaxh], zTr[nMaxh];
	Double_t r[nMaxh], u[nMaxh] , v[nMaxh];
	
	AliMFTCACell *celltr;
	Int_t nptr = 0;
	
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
	for (Int_t i = 0; i< nptr; i++) {
		AliInfo(Form("x,y,z = %f %f %f",xTr[i],yTr[i],zTr[i]));
	}
//	Double_t distL = 0., q2=0.;
//	Double_t sagitta = AliMFTTrackExtrap::Sagitta(nptr, xTr, yTr, distL,q2);
//	AliInfo(Form("sagitta = %f",sagitta));;
//	AliInfo(Form("distL = %f",distL));;
//	
//	AliInfo(Form("q2 = %f => pt = %f ",q2, 0.01/q2/2.*0.3*b[2]/10.));
//	
//	//fPt = 0.3*TMath::Abs(b[2]/10.)*distL*distL/8./sagitta;
//	fPt = 0.3*TMath::Abs(b[2]/10.)*AliMFTTrackExtrap::CircleRegression(nptr, xTr, yTr)*TMath::Sign(1.,sagitta);
// //fPt =0.01/q2/2.*0.3*b[2]/10.;
	
	
	
	
	
	for (Int_t iptr = 0; iptr < nptr; iptr++) {
		
		r[iptr] = TMath::Sqrt(yTr[iptr]*yTr[iptr]+xTr[iptr]*xTr[iptr]);
		u[iptr] = xTr[iptr]/r[iptr]/r[iptr];
		v[iptr] = yTr[iptr]/r[iptr]/r[iptr];
		
		
		AliInfo(Form("u,v,r = %f %f %f ",u[iptr],v[iptr],r[iptr]));
		
	}
	Double_t fdump, slopeX_Z, slopeY_Z;
	AliMFTTrackExtrap::LinearRegression(nptr, xTr,zTr, fdump, slopeX_Z);
	AliMFTTrackExtrap::LinearRegression(nptr, yTr,zTr, fdump, slopeY_Z);
	
	Double_t phi = TMath::ATan2(slopeY_Z,slopeX_Z);
	Double_t xS,x0;
	Double_t chi2 = AliMFTTrackExtrap::LinearRegression(nptr, u,v, x0, xS);
  AliInfo(Form("chi2 = %f x0, xS %f %f",chi2,x0,xS));;

  if (TMath::Abs(x0)<1.e-10) {
    fPt = 1.e6;
    return;
  }
  
	Double_t rx = -xS/(2.*x0);
	Double_t ry =  1./(2.*x0);
	Double_t rr = TMath::Sqrt(rx*rx+ry*ry);
	AliInfo(Form("Rx,Ry = %f %f",rx,ry));;
	AliInfo(Form("R = %f",rr));;
	AliInfo(Form("Phi = %f",phi));;
	
		Double_t zmean = 0.5 * (AliMFTConstants::DefaultPlaneZ(0) + AliMFTConstants::DefaultPlaneZ(9));
		const Double_t x[3] = {0.,0.,zmean};
		Double_t b[3] = {0.,0.,0.};
		TGeoGlobalMagField::Instance()->Field(x,b);
		AliInfo(Form("Field = %e %e %e",b[0],b[1],b[2]));
		
	fPt = 0.3*TMath::Abs(b[2]/10.)*rr/100.*TMath::Sign(1.,b[2]*rx*(-phi));
	
	
	AliInfo(Form("pt = %f",fPt));;
	
	
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
