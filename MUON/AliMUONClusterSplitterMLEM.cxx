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

// $Id$

/// \class AliMUONClusterSplitterMLEM
/// 
/// Splitter class for the MLEM algorithm...
///
/// FIXME: describe it a little bit more here...
///
/// \author Laurent Aphecetche (for the "new" C++ structure) and 
/// Alexander Zinchenko, JINR Dubna, for the hardcore of it ;-)

#include "AliMUONClusterSplitterMLEM.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMUONPad.h"
#include "AliMUONPad.h"
#include "AliMpStationType.h"
#include "AliMUONConstants.h"
#include "AliMpDEManager.h"
#include "AliMUONMathieson.h"

#include "TClonesArray.h"
#include "TH2.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TRandom.h"

ClassImp(AliMUONClusterSplitterMLEM)

const Double_t AliMUONClusterSplitterMLEM::fgkCouplMin = 1.e-3; // threshold on coupling 

//_____________________________________________________________________________
AliMUONClusterSplitterMLEM::AliMUONClusterSplitterMLEM(Int_t detElemId, 
                                                       TObjArray* fPixArray) 
: TObject(),
fPixArray(fPixArray),
fMathieson(0x0),
fDetElemId(detElemId),
fNpar(0),
fQtot(0),
fnCoupled(0)
{
  /// Constructor
  
  AliMpStationType stationType = AliMpDEManager::GetStationType(fDetElemId);
  
  Float_t kx3 = AliMUONConstants::SqrtKx3();
  Float_t ky3 = AliMUONConstants::SqrtKy3();
  Float_t pitch = AliMUONConstants::Pitch();
  
  if ( stationType == kStation1 )
  {
    kx3 = AliMUONConstants::SqrtKx3St1();
    ky3 = AliMUONConstants::SqrtKy3St1();
    pitch = AliMUONConstants::PitchSt1();
  }
  
  fMathieson = new AliMUONMathieson;
  
  fMathieson->SetPitch(pitch);
  fMathieson->SetSqrtKx3AndDeriveKx2Kx4(kx3);
  fMathieson->SetSqrtKy3AndDeriveKy2Ky4(ky3);
  
}

//_____________________________________________________________________________
AliMUONClusterSplitterMLEM::~AliMUONClusterSplitterMLEM()
{
  delete fMathieson;
}

//_____________________________________________________________________________
void 
AliMUONClusterSplitterMLEM::AddBin(TH2 *mlem, 
                                   Int_t ic, Int_t jc, Int_t mode, 
                                   Bool_t *used, TObjArray *pix)
{
  /// Add a bin to the cluster
  
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t cont1, cont = mlem->GetCellContent(jc,ic);
  AliMUONPad *pixPtr = 0;
  
  for (Int_t i=TMath::Max(ic-1,1); i<=TMath::Min(ic+1,ny); i++) {
    for (Int_t j=TMath::Max(jc-1,1); j<=TMath::Min(jc+1,nx); j++) {
      if (i != ic && j != jc) continue;
      if (used[(i-1)*nx+j-1]) continue;
      cont1 = mlem->GetCellContent(j,i);
      if (mode && cont1 > cont) continue;
      used[(i-1)*nx+j-1] = kTRUE;
      if (cont1 < 0.5) continue;
      if (pix) pix->Add(BinToPix(mlem,j,i)); 
      else {
        pixPtr = new AliMUONPad (mlem->GetXaxis()->GetBinCenter(j), 
                                   mlem->GetYaxis()->GetBinCenter(i), 0, 0, cont1);
        fPixArray->Add((TObject*)pixPtr);
      }
      AddBin(mlem, i, j, mode, used, pix); // recursive call
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONClusterSplitterMLEM::AddCluster(Int_t ic, Int_t nclust, 
                                       TMatrixD& aijcluclu, 
                                       Bool_t *used, Int_t *clustNumb, Int_t &nCoupled)
{
  /// Add a cluster to the group of coupled clusters
  
  for (Int_t i=0; i<nclust; i++) {
    if (used[i]) continue;
    if (aijcluclu(i,ic) < fgkCouplMin) continue;
    used[i] = kTRUE;
    clustNumb[nCoupled++] = i;
    AddCluster(i, nclust, aijcluclu, used, clustNumb, nCoupled);
  }
}

//_____________________________________________________________________________
TObject* 
AliMUONClusterSplitterMLEM::BinToPix(TH2 *mlem,
                                     Int_t jc, Int_t ic)
{
  /// Translate histogram bin to pixel 
  
  Double_t yc = mlem->GetYaxis()->GetBinCenter(ic);
  Double_t xc = mlem->GetXaxis()->GetBinCenter(jc);
  
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPad *pixPtr = NULL;
  
  // Compare pixel and bin positions
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPad*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 0.5) continue;
    if (TMath::Abs(pixPtr->Coord(0)-xc)<1.e-4 && TMath::Abs(pixPtr->Coord(1)-yc)<1.e-4) 
    {
      return (TObject*) pixPtr;
    }
  }
  AliError(Form(" Something wrong ??? %f %f ", xc, yc));
  return NULL;
}

//_____________________________________________________________________________
Float_t
AliMUONClusterSplitterMLEM::ChargeIntegration(Double_t x, Double_t y,
                                              const AliMUONPad& pad)
{
  /// Compute the Mathieson integral on pad area, assuming the center
  /// of the Mathieson is at (x,y)
  
  TVector2 lowerLeft(TVector2(x,y)-pad.Position()-pad.Dimensions());
  TVector2 upperRight(lowerLeft + pad.Dimensions()*2.0);
  
	return fMathieson->IntXY(lowerLeft.X(),lowerLeft.Y(),
                           upperRight.X(),upperRight.Y());
}

//_____________________________________________________________________________
void 
AliMUONClusterSplitterMLEM::Fcn1(const AliMUONCluster& cluster, 
                                    Int_t & /*fNpar*/, Double_t * /*gin*/, 
                                    Double_t &f, Double_t *par, Int_t /*iflag*/)
{
  /// Fit for one track
  
  Int_t indx, npads=0;
  Double_t charge, delta, coef=0, chi2=0, qTot = 0;
  
  for (Int_t j=0; j< cluster.Multiplicity(); ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    if ( pad->Status() != 1 ) continue;
    if ( pad->DX() > 0 ) npads++; // exclude virtual pads
    qTot += pad->Charge(); // c.fXyq[2][j];
    charge = 0;
    for (Int_t i=fNpar/3; i>=0; --i)
    { // sum over tracks
      indx = i<2 ? 2*i : 2*i+1;
      if (fNpar == 2) 
      {
        coef = 1;
      }
      else 
      {
        coef = i==fNpar/3 ? par[indx+2] : 1-coef;
      }
      coef = TMath::Max (coef, 0.);
      if ( fNpar == 8 && i < 2) 
      {
        coef = i==1 ? coef*par[indx+2] : coef - par[7];
      }
      coef = TMath::Max (coef, 0.);
      charge += ChargeIntegration(par[indx],par[indx+1],*pad);
    }
    charge *= fQtot;
    delta = charge - pad->Charge(); //c.fXyq[2][j];
    delta *= delta;
    delta /= pad->Charge(); //c.fXyq[2][j];
    chi2 += delta;
  } // for (Int_t j=0;
  f = chi2; 
  Double_t qAver = qTot/npads;
  f = chi2/qAver;
}

//_____________________________________________________________________________
Int_t 
AliMUONClusterSplitterMLEM::Fit(const AliMUONCluster& cluster,
                                Int_t iSimple, Int_t nfit, 
                                Int_t *clustFit, TObjArray **clusters, 
                                Double_t *parOk,
                                TObjArray& clusterList)
{
  /// Find selected clusters to selected pad charges
  
  //  AliDebug(2,Form("iSimple=%d nfit=%d",iSimple,nfit));
  
  TH2D *mlem = (TH2D*) gROOT->FindObject("mlem");
  Double_t xmin = mlem->GetXaxis()->GetXmin() - mlem->GetXaxis()->GetBinWidth(1);
  Double_t xmax = mlem->GetXaxis()->GetXmax() + mlem->GetXaxis()->GetBinWidth(1);
  Double_t ymin = mlem->GetYaxis()->GetXmin() - mlem->GetYaxis()->GetBinWidth(1);
  Double_t ymax = mlem->GetYaxis()->GetXmax() + mlem->GetYaxis()->GetBinWidth(1);
  Double_t step[3]={0.01,0.002,0.02}, xPad = 0, yPad = 99999;
  
  // Number of pads to use and number of virtual pads
  Int_t npads = 0, nVirtual = 0, nfit0 = nfit;
  for (Int_t i=0; i<cluster.Multiplicity(); ++i ) 
  {
    AliMUONPad* pad = cluster.Pad(i);
    if ( pad->DX() < 0 ) ++nVirtual;
    if ( pad->Status() !=1 ) continue;
    if ( pad->DX() > 0 )
    {
      ++npads;
      if (yPad > 9999) 
      { 
        xPad = pad->X();//fXyq[0][i]; 
        yPad = pad->Y();//fXyq[1][i]; 
      } 
      else 
      {
        if (pad->DY() < pad->DX() ) //fXyq[4][i] < fXyq[3][i]) 
        {
          yPad = pad->Y();//fXyq[1][i]; 
        }
        else 
        {
          xPad = pad->X();//fXyq[0][i]; 
        }
      }
    }
  }
  
  fNpar = 0;
  fQtot = 0;
  
  if (npads < 2) return 0; 
  
  // FIXME : AliWarning("Reconnect the following code for hit/track passing ?");
  
  //  Int_t tracks[3] = {-1, -1, -1};
  
  /*
   Int_t digit = 0;
   AliMUONDigit *mdig = 0;
   for (Int_t cath=0; cath<2; cath++) {  
     for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
       if (fPadIJ[0][i] != cath) continue;
       if (fPadIJ[1][i] != 1) continue;
       if (fXyq[3][i] < 0) continue; // exclude virtual pads
       digit = TMath::Nint (fXyq[5][i]);
       if (digit >= 0) mdig = fInput->Digit(cath,digit);
       else mdig = fInput->Digit(TMath::Even(cath),-digit-1);
       //if (!mdig) mdig = fInput->Digit(TMath::Even(cath),digit);
       if (!mdig) continue; // protection for cluster display
       if (mdig->Hit() >= 0) {
         if (tracks[0] < 0) {
           tracks[0] = mdig->Hit();
           tracks[1] = mdig->Track(0);
         } else if (mdig->Track(0) < tracks[1]) {
           tracks[0] = mdig->Hit();
           tracks[1] = mdig->Track(0);
         }
       }
       if (mdig->Track(1) >= 0 && mdig->Track(1) != tracks[1]) {
         if (tracks[2] < 0) tracks[2] = mdig->Track(1);
         else tracks[2] = TMath::Min (tracks[2], mdig->Track(1));
       }
     } // for (Int_t i=0;
  } // for (Int_t cath=0;
   */
  
  // Get number of pads in X and Y 
//  Int_t nInX = 0, nInY;
//  PadsInXandY(cluster,nInX, nInY);
  const Int_t kStatusToTest(1);
  
  AliMpIntPair nofPads = cluster.NofPads(kStatusToTest);
  Int_t nInX = nofPads.GetFirst();
  Int_t nInY = nofPads.GetSecond();
  //cout << " nInX and Y: " << nInX << " " << nInY << endl;
  
  Int_t nfitMax = 3; 
  nfitMax = TMath::Min (nfitMax, (npads + 1) / 3);
  if (nfitMax > 1) {
    if (nInX < 3 && nInY < 3 || nInX == 3 && nInY < 3 || nInX < 3 && nInY == 3) nfitMax = 1; // not enough pads in each direction
  }
  if (nfit > nfitMax) nfit = nfitMax;
  
  // Take cluster maxima as fitting seeds
  TObjArray *pix;
  AliMUONPad *pixPtr;
  Int_t npxclu;
  Double_t cont, cmax = 0, xseed = 0, yseed = 0, errOk[8], qq = 0;
  Double_t xyseed[3][2], qseed[3], xyCand[3][2] = {{0},{0}}, sigCand[3][2] = {{0},{0}};
  
  for (Int_t ifit=1; ifit<=nfit0; ifit++) 
  {
    cmax = 0;
    pix = clusters[clustFit[ifit-1]];
    npxclu = pix->GetEntriesFast();
    //qq = 0;
    for (Int_t clu=0; clu<npxclu; ++clu) 
    {
      pixPtr = (AliMUONPad*) pix->UncheckedAt(clu);
      cont = pixPtr->Charge();
      fQtot += cont;
      if (cont > cmax) 
      { 
        cmax = cont; 
        xseed = pixPtr->Coord(0);
        yseed = pixPtr->Coord(1);
      }
      qq += cont;
      xyCand[0][0] += pixPtr->Coord(0) * cont;
      xyCand[0][1] += pixPtr->Coord(1) * cont;
      sigCand[0][0] += pixPtr->Coord(0) * pixPtr->Coord(0) * cont;
      sigCand[0][1] += pixPtr->Coord(1) * pixPtr->Coord(1) * cont;
    }
    xyseed[ifit-1][0] = xseed;
    xyseed[ifit-1][1] = yseed;
    qseed[ifit-1] = cmax;
  } // for (Int_t ifit=1;
  
  xyCand[0][0] /= qq; // <x>
  xyCand[0][1] /= qq; // <y>
  sigCand[0][0] = sigCand[0][0]/qq - xyCand[0][0]*xyCand[0][0]; // <x^2> - <x>^2
  sigCand[0][0] = sigCand[0][0] > 0 ? TMath::Sqrt (sigCand[0][0]) : 0;
  sigCand[0][1] = sigCand[0][1]/qq - xyCand[0][1]*xyCand[0][1]; // <y^2> - <y>^2
  sigCand[0][1] = sigCand[0][1] > 0 ? TMath::Sqrt (sigCand[0][1]) : 0;
//  if (fDebug) cout << xyCand[0][0] << " " << xyCand[0][1] << " " << sigCand[0][0] << " " << sigCand[0][1] << endl;
  
  Int_t nDof, maxSeed[3];//, nMax = 0;
    Double_t fmin, chi2o = 9999, chi2n;
    
    TMath::Sort(nfit0, qseed, maxSeed, kTRUE); // in decreasing order
    
    Double_t *gin = 0, func0, func1, param[8], step0[8];
    Double_t param0[2][8]={{0},{0}}, deriv[2][8]={{0},{0}}; 
    Double_t shift[8], stepMax, derMax, parmin[8], parmax[8], func2[2], shift0;
    Double_t delta[8], scMax, dder[8], estim, shiftSave = 0;
    Int_t min, max, nCall = 0, memory[8] = {0}, nLoop, idMax = 0, iestMax = 0, nFail;
    Double_t rad, dist[3] = {0};
    
    // Try to fit with one-track hypothesis, then 2-track. If chi2/dof is 
    // lower, try 3-track (if number of pads is sufficient).
    for (Int_t iseed=0; iseed<nfit; iseed++) 
    {
      
      if (iseed) 
      { 
        for (Int_t j=0; j<fNpar; j++) 
        {
          param[j] = parOk[j]; 
        }
      } // for bounded params
      
      for (Int_t j=0; j<3; j++) 
      {
        step0[fNpar+j] = shift[fNpar+j] = step[j];
      }
      
      if (nfit == 1) 
      {
        param[fNpar] = xyCand[0][0]; // take COG
      }
      else 
      {
        param[fNpar] = xyseed[maxSeed[iseed]][0];
      }
      parmin[fNpar] = xmin; 
      parmax[fNpar++] = xmax; 
      if (nfit == 1) 
      {
        param[fNpar] = xyCand[0][1]; // take COG
      }
      else 
      {
        param[fNpar] = xyseed[maxSeed[iseed]][1];
      }
      parmin[fNpar] = ymin; 
      parmax[fNpar++] = ymax; 
      if (fNpar > 2) 
      {
        param[fNpar] = fNpar == 4 ? 0.5 : 0.3;
        parmin[fNpar] = 0; 
        parmax[fNpar++] = 1; 
      }
      if (iseed) 
      { 
        for (Int_t j=0; j<fNpar; j++) 
        {
          param0[1][j] = 0; 
        }
      }
      
      // Try new algorithm
      min = nLoop = 1; stepMax = func2[1] = derMax = 999999; nFail = 0;
      
      while (1) 
      {
        max = !min;
        Fcn1(cluster,fNpar, gin, func0, param, 1); nCall++;
        //cout << " Func: " << func0 << endl;
        
        func2[max] = func0;
        for (Int_t j=0; j<fNpar; j++) 
        {
          param0[max][j] = param[j];
          delta[j] = step0[j];
          param[j] += delta[j] / 10;
          if (j > 0) param[j-1] -= delta[j-1] / 10;
          Fcn1(cluster,fNpar, gin, func1, param, 1); nCall++;
          deriv[max][j] = (func1 - func0) / delta[j] * 10; // first derivative
                                                           //cout << j << " " << deriv[max][j] << endl;
          dder[j] = param0[0][j] != param0[1][j] ? (deriv[0][j] - deriv[1][j]) / 
            (param0[0][j] - param0[1][j]) : 0; // second derivative
        }
        param[fNpar-1] -= delta[fNpar-1] / 10;
        if (nCall > 2000) break;
        
        min = func2[0] < func2[1] ? 0 : 1;
        nFail = min == max ? 0 : nFail + 1;
        
        stepMax = derMax = estim = 0;
        for (Int_t j=0; j<fNpar; j++) 
        { 
          // Estimated distance to minimum
          shift0 = shift[j];
          if (nLoop == 1) 
          {
            shift[j] = TMath::Sign (step0[j], -deriv[max][j]); // first step
          }
          else if (TMath::Abs(deriv[0][j]) < 1.e-3 && TMath::Abs(deriv[1][j]) < 1.e-3) 
          {
            shift[j] = 0;
          }
          else if (deriv[min][j]*deriv[!min][j] > 0 && TMath::Abs(deriv[min][j]) > TMath::Abs(deriv[!min][j])
                   || TMath::Abs(deriv[0][j]-deriv[1][j]) < 1.e-3 || TMath::Abs(dder[j]) < 1.e-6) 
          {
            shift[j] = -TMath::Sign (shift[j], (func2[0]-func2[1]) * (param0[0][j]-param0[1][j]));
          }
          if (min == max) 
          { 
            if (memory[j] > 1) 
            { 
              shift[j] *= 2; 
            } 
            memory[j]++;
          }
          else 
          {
            shift[j] = dder[j] != 0 ? -deriv[min][j] / dder[j] : 0;
            memory[j] = 0;
          }
          
          if (TMath::Abs(shift[j])/step0[j] > estim) 
          { 
            estim = TMath::Abs(shift[j])/step0[j];
            iestMax = j;
          }
          
          // Too big step
          if (TMath::Abs(shift[j])/step0[j] > 10) shift[j] = TMath::Sign(10.,shift[j]) * step0[j]; // 
          
          // Failed to improve minimum
          if (min != max) 
          {
            memory[j] = 0;
            param[j] = param0[min][j];
            if (TMath::Abs(shift[j]+shift0) > 0.1*step0[j]) 
            {
              shift[j] = (shift[j] + shift0) / 2;
            }
            else 
            {
              shift[j] /= -2;
            }
          } 
          
          // Too big step
          if (TMath::Abs(shift[j]*deriv[min][j]) > func2[min]) 
          {
            shift[j] = TMath::Sign (func2[min]/deriv[min][j], shift[j]);
          }
          
          // Introduce step relaxation factor
          if (memory[j] < 3) 
          {
            scMax = 1 + 4 / TMath::Max(nLoop/2.,1.);
            if (TMath::Abs(shift0) > 0 && TMath::Abs(shift[j]/shift0) > scMax) 
            {
              shift[j] = TMath::Sign (shift0*scMax, shift[j]);
            }
          }
          param[j] += shift[j]; 
          //MLEM Check parameter limits 27-12-2004
          if (param[j] < parmin[j]) 
          { 
            shift[j] = parmin[j] - param[j]; 
            param[j] = parmin[j]; 
          } 
          else if (param[j] > parmax[j]) 
          {
            shift[j] = parmax[j] - param[j];
            param[j] = parmax[j];
          }
          //cout << " xxx " << j << " " << shift[j] << " " << param[j] << endl;
          stepMax = TMath::Max (stepMax, TMath::Abs(shift[j]/step0[j]));
          if (TMath::Abs(deriv[min][j]) > derMax) 
          {
            idMax = j;
            derMax = TMath::Abs (deriv[min][j]);
          }
        } // for (Int_t j=0; j<fNpar;
        
        if (estim < 1 && derMax < 2 || nLoop > 150) break; // minimum was found
        
        nLoop++;
        
        // Check for small step
        if (shift[idMax] == 0) 
        { 
          shift[idMax] = step0[idMax]/10; 
          param[idMax] += shift[idMax]; 
          continue; 
        }
        
        if (!memory[idMax] && derMax > 0.5 && nLoop > 10) 
        {
          if (dder[idMax] != 0 && TMath::Abs(deriv[min][idMax]/dder[idMax]/shift[idMax]) > 10) 
          {
            if (min == max) dder[idMax] = -dder[idMax];
            shift[idMax] = -deriv[min][idMax] / dder[idMax] / 10; 
            param[idMax] += shift[idMax];
            stepMax = TMath::Max (stepMax, TMath::Abs(shift[idMax])/step0[idMax]);
            if (min == max) shiftSave = shift[idMax];
          }
          if (nFail > 10) 
          {
            param[idMax] -= shift[idMax];
            shift[idMax] = 4 * shiftSave * (gRandom->Rndm(0) - 0.5);
            param[idMax] += shift[idMax];
          }
        }      
  } // while (1)
      
      fmin = func2[min];
      
      nDof = npads - fNpar + nVirtual;
      if (!nDof) nDof++;
      chi2n = fmin / nDof;
//      if (fDebug) cout << " Chi2 " << chi2n << " " << fNpar << endl;
      
      if (chi2n*1.2+1.e-6 > chi2o ) { fNpar -= 3; break; }
      
      // Save parameters and errors
      
      if (nInX == 1) {
        // One pad per direction 
        for (Int_t i=0; i<fNpar; i++) if (i == 0 || i == 2 || i == 5) param0[min][i] = xPad;
      }
      if (nInY == 1) {
        // One pad per direction 
        for (Int_t i=0; i<fNpar; i++) if (i == 1 || i == 3 || i == 6) param0[min][i] = yPad;
      }
      
      /*
       if (iseed > 0) {
         // Find distance to the nearest neighbour
         dist[0] = dist[1] = TMath::Sqrt ((param0[min][0]-param0[min][2])*
                                          (param0[min][0]-param0[min][2])
                                          +(param0[min][1]-param0[min][3])*
                                          (param0[min][1]-param0[min][3]));
         if (iseed > 1) {
           dist[2] = TMath::Sqrt ((param0[min][0]-param0[min][5])*
                                  (param0[min][0]-param0[min][5])
                                  +(param0[min][1]-param0[min][6])*
                                  (param0[min][1]-param0[min][6]));
           rad = TMath::Sqrt ((param0[min][2]-param0[min][5])*
                              (param0[min][2]-param0[min][5])
                              +(param0[min][3]-param0[min][6])*
                              (param0[min][3]-param0[min][6]));
           if (dist[2] < dist[0]) dist[0] = dist[2];
           if (rad < dist[1]) dist[1] = rad;
           if (rad < dist[2]) dist[2] = rad;
         }
         cout << dist[0] << " " << dist[1] << " " << dist[2] << endl;
         if (dist[TMath::LocMin(iseed+1,dist)] < 1.) { fNpar -= 3; break; }
       }
       */
      
      for (Int_t i=0; i<fNpar; i++) {
        parOk[i] = param0[min][i];
        //errOk[i] = fmin;
        errOk[i] = chi2n;
        // Bounded params
        parOk[i] = TMath::Max (parOk[i], parmin[i]);
        parOk[i] = TMath::Min (parOk[i], parmax[i]);
      }
      
      chi2o = chi2n;
      if (fmin < 0.1) break; // !!!???
  } // for (Int_t iseed=0; 
    
//    if (fDebug) {
//      for (Int_t i=0; i<fNpar; i++) {
//        if (i == 4 || i == 7) {
//          if (i == 7 || i == 4 && fNpar < 7) cout << parOk[i] << endl;
//          else cout << parOk[i] * (1-parOk[7]) << endl;
//          continue;
//        }
//        cout << parOk[i] << " " << errOk[i] << endl;
//      }
//    }
    nfit = (fNpar + 1) / 3;
    dist[0] = dist[1] = dist[2] = 0;
    
    if (nfit > 1) {
      // Find distance to the nearest neighbour
      dist[0] = dist[1] = TMath::Sqrt ((parOk[0]-parOk[2])*
                                       (parOk[0]-parOk[2])
                                       +(parOk[1]-parOk[3])*
                                       (parOk[1]-parOk[3]));
      if (nfit > 2) {
        dist[2] = TMath::Sqrt ((parOk[0]-parOk[5])*
                               (parOk[0]-parOk[5])
                               +(parOk[1]-parOk[6])*
                               (parOk[1]-parOk[6]));
        rad = TMath::Sqrt ((parOk[2]-parOk[5])*
                           (parOk[2]-parOk[5])
                           +(parOk[3]-parOk[6])*
                           (parOk[3]-parOk[6]));
        if (dist[2] < dist[0]) dist[0] = dist[2];
        if (rad < dist[1]) dist[1] = rad;
        if (rad < dist[2]) dist[2] = rad;
      }
    }
    
    Int_t indx;
    
    //if (!fDraw) {
    
    Double_t coef = 0;
    if (iSimple) fnCoupled = 0;
    //for (Int_t j=0; j<nfit; j++) {
    for (Int_t j=nfit-1; j>=0; j--) {
      indx = j<2 ? j*2 : j*2+1;  
      if (nfit == 1) coef = 1;
      else coef = j==nfit-1 ? parOk[indx+2] : 1-coef;
      coef = TMath::Max (coef, 0.);
      if (nfit == 3 && j < 2) coef = j==1 ? coef*parOk[indx+2] : coef - parOk[7];
      coef = TMath::Max (coef, 0.);
      
      //void AliMUONClusterFinderMLEM::AddRawCluster(Double_t x, Double_t y, 
      //                                             Double_t qTot, Double_t fmin,
      //                                             Int_t nfit, Int_t *tracks, 
      //                                             Double_t /*sigx*/, 
      //                                             Double_t /*sigy*/, 
      //                                             Double_t /*dist*/)
      
      if ( coef*fQtot >= 14 )
      {
        AliMUONCluster* cluster = new AliMUONCluster();
        
        cluster->SetCharge(coef*fQtot,coef*fQtot);
        cluster->SetPosition(TVector2(parOk[indx],parOk[indx+1]),TVector2(sigCand[0][0],sigCand[0][1]));
        cluster->SetChi2(dist[TMath::LocMin(nfit,dist)]);
        
        // FIXME: we miss some information in this cluster, as compared to 
        // the original AddRawCluster code.
        
        AliDebug(2,Form("Adding RawCluster detElemId %4d mult %2d charge %5d (xl,yl)=(%9.6g,%9.6g)",
                        fDetElemId,cluster->Multiplicity(),(Int_t)cluster->Charge(),
                        cluster->Position().X(),cluster->Position().Y()));
        
        clusterList.Add(cluster);
      }
      //      AddRawCluster (parOk[indx], // double x
      //                     parOk[indx+1], // double y
      //                     coef*qTot, // double charge
      //                     errOk[indx], // double fmin
      //                     nfit0+10*nfit+100*nMax+10000*fnCoupled, // int nfit
      //                     tracks, // int* tracks
      //                     sigCand[0][0], // double sigx
      //                     sigCand[0][1], // double sigy
      //                     dist[TMath::LocMin(nfit,dist)] // double dist
      //                     );
    }
    return nfit;
    }  


//_____________________________________________________________________________
void
AliMUONClusterSplitterMLEM::Split(const AliMUONCluster& cluster,
                                  TH2 *mlem,
                                  Double_t *coef,
                                  TObjArray& clusterList)
{
  /// The main steering function to work with clusters of pixels in anode
  /// plane (find clusters, decouple them from each other, merge them (if
  /// necessary), pick up coupled pads, call the fitting function)
  
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Int_t nPix = fPixArray->GetEntriesFast();
  
  Bool_t *used = new Bool_t[ny*nx];
  Double_t cont;
  Int_t nclust = 0, indx, indx1;
  
  for (Int_t i=0; i<ny*nx; i++) used[i] = kFALSE; 
  
  TObjArray *clusters[200]={0};
  TObjArray *pix;
  
  // Find clusters of histogram bins (easier to work in 2-D space)
  for (Int_t i=1; i<=ny; i++) 
  {
    for (Int_t j=1; j<=nx; j++) 
    {
      indx = (i-1)*nx + j - 1;
      if (used[indx]) continue;
      cont = mlem->GetCellContent(j,i);
      if (cont < 0.5) continue;
      pix = new TObjArray(20);
      used[indx] = 1;
      pix->Add(BinToPix(mlem,j,i));
      AddBin(mlem, i, j, 0, used, pix); // recursive call
      if (nclust >= 200) AliFatal(" Too many clusters !!!");
      clusters[nclust++] = pix;
    } // for (Int_t j=1; j<=nx; j++) {
    } // for (Int_t i=1; i<=ny;
//  if (fDebug) cout << nclust << endl;
  delete [] used; used = 0;
  
  // Compute couplings between clusters and clusters to pads
  Int_t npad = cluster.Multiplicity();
  
  // Exclude pads with overflows
  for (Int_t j=0; j<npad; ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    if ( pad->IsSaturated() )
    {
      pad->SetStatus(-5); 
    }
    else 
    {
      pad->SetStatus(0);
    }
  }
  
  // Compute couplings of clusters to pads
  TMatrixD aijclupad(nclust,npad);
  aijclupad = 0;
  Int_t npxclu;
  for (Int_t iclust=0; iclust<nclust; ++iclust) 
  {
    pix = clusters[iclust];
    npxclu = pix->GetEntriesFast();
    for (Int_t i=0; i<npxclu; ++i) 
    {
      indx = fPixArray->IndexOf(pix->UncheckedAt(i));
      for (Int_t j=0; j<npad; ++j) 
      {
        AliMUONPad* pad = cluster.Pad(j);
        if ( pad->Status() < 0 && pad->Status() != -5) continue;
        if (coef[j*nPix+indx] < fgkCouplMin) continue;
        aijclupad(iclust,j) += coef[j*nPix+indx];
      }
    }
  }
  
  // Compute couplings between clusters
  TMatrixD aijcluclu(nclust,nclust);
  aijcluclu = 0;
  for (Int_t iclust=0; iclust<nclust; ++iclust) 
  {
    for (Int_t j=0; j<npad; ++j) 
    {
      // Exclude overflows
      if ( cluster.Pad(j)->Status() < 0) continue;
      if (aijclupad(iclust,j) < fgkCouplMin) continue;
      for (Int_t iclust1=iclust+1; iclust1<nclust; iclust1++) 
      {
        if (aijclupad(iclust1,j) < fgkCouplMin) continue;
        aijcluclu(iclust,iclust1) += 
          TMath::Sqrt (aijclupad(iclust,j)*aijclupad(iclust1,j));
      }
    }
  }
  for (Int_t iclust=0; iclust<nclust; ++iclust) 
  {
    for (Int_t iclust1=iclust+1; iclust1<nclust; ++iclust1) 
    {
      aijcluclu(iclust1,iclust) = aijcluclu(iclust,iclust1);
    }
  }
  
  // Find groups of coupled clusters
  used = new Bool_t[nclust];
  for (Int_t i=0; i<nclust; i++) used[i] = kFALSE;
  Int_t *clustNumb = new Int_t[nclust];
  Int_t nCoupled, nForFit, minGroup[3], clustFit[3], nfit = 0;
  Double_t parOk[8];
  
  for (Int_t igroup=0; igroup<nclust; igroup++) 
  {
    if (used[igroup]) continue;
    used[igroup] = kTRUE;
    clustNumb[0] = igroup;
    nCoupled = 1;
    // Find group of coupled clusters
    AddCluster(igroup, nclust, aijcluclu, used, clustNumb, nCoupled); // recursive
    
    //    if (fDebug) {                                                                      
    //      cout << " nCoupled: " << nCoupled << endl;                                                                      
    //      for (Int_t i=0; i<nCoupled; i++) cout << clustNumb[i] << " "; cout << endl;                                                                       
    //    }
    
    fnCoupled = nCoupled;
    
    while (nCoupled > 0) 
    {
      if (nCoupled < 4) 
      {
        nForFit = nCoupled;
        for (Int_t i=0; i<nCoupled; i++) clustFit[i] = clustNumb[i];
      } 
      else 
      {
        // Too many coupled clusters to fit - try to decouple them
        // Find the lowest coupling of 1, 2, min(3,nLinks/2) pixels with 
        // all the others in the group 
        for (Int_t j=0; j<3; j++) minGroup[j] = -1;
        /*Double_t coupl =*/ MinGroupCoupl(nCoupled, clustNumb, aijcluclu, minGroup);
        
        // Flag clusters for fit
        nForFit = 0;
        while (minGroup[nForFit] >= 0 && nForFit < 3)
        {
          clustFit[nForFit] = clustNumb[minGroup[nForFit]];
          clustNumb[minGroup[nForFit]] -= 999;
          nForFit++;
        }
      } // else
      
      // Select pads for fit. 
      if (SelectPad(cluster,nCoupled, nForFit, clustNumb, clustFit, aijclupad) < 3 && nCoupled > 1) 
      {
        // Deselect pads
        for (Int_t j=0; j<npad; ++j)
        {
          AliMUONPad* pad = cluster.Pad(j);
          if ( pad->Status()==1 ) pad->SetStatus(0);
          if ( pad->Status()==-9) pad->SetStatus(-5);
        }
        // Merge the failed cluster candidates (with too few pads to fit) with 
        // the one with the strongest coupling
        Merge(cluster,nForFit, nCoupled, clustNumb, clustFit, clusters, aijcluclu, aijclupad);
      } 
      else 
      {
        // Do the fit
        nfit = Fit(cluster,0, nForFit, clustFit, clusters, parOk, clusterList);
      }
      
      // Subtract the fitted charges from pads with strong coupling and/or
      // return pads for further use
      UpdatePads(cluster,nfit, parOk);
      
      // Mark used pads
      for (Int_t j=0; j<npad; ++j) 
      {
        AliMUONPad* pad = cluster.Pad(j);
        if ( pad->Status()==1 ) pad->SetStatus(-1);
        if ( pad->Status()==-9) pad->SetStatus(-5);
      }
      
      // Sort the clusters (move to the right the used ones)
      Int_t beg = 0, end = nCoupled - 1;
      while (beg < end) 
      {
        if (clustNumb[beg] >= 0) { ++beg; continue; }
        for (Int_t j=end; j>beg; --j) 
        {
          if (clustNumb[j] < 0) continue;
          end = j - 1;
          indx = clustNumb[beg];
          clustNumb[beg] = clustNumb[j];
          clustNumb[j] = indx;
          break;
        }
        ++beg;
      }
      
      nCoupled -= nForFit;
      if (nCoupled > 3) 
      {
        // Remove couplings of used clusters
        for (Int_t iclust=nCoupled; iclust<nCoupled+nForFit;++ iclust) 
        {
          indx = clustNumb[iclust] + 999;
          for (Int_t iclust1=0; iclust1<nCoupled; ++iclust1) 
          {
            indx1 = clustNumb[iclust1];
            aijcluclu(indx,indx1) = aijcluclu(indx1,indx) = 0;
          }
        }
        
        // Update the remaining clusters couplings (exclude couplings from 
        // the used pads)
        for (Int_t j=0; j<npad; ++j) 
        {
          AliMUONPad* pad = cluster.Pad(j);
          if ( pad->Status() != -1) continue;
          for (Int_t iclust=0; iclust<nCoupled; ++iclust) 
          {
            indx = clustNumb[iclust];
            if (aijclupad(indx,j) < fgkCouplMin) continue;
            for (Int_t iclust1=iclust+1; iclust1<nCoupled; ++iclust1) 
            {
              indx1 = clustNumb[iclust1];
              if (aijclupad(indx1,j) < fgkCouplMin) continue;
              // Check this
              aijcluclu(indx,indx1) -= 
                TMath::Sqrt (aijclupad(indx,j)*aijclupad(indx1,j));
              aijcluclu(indx1,indx) = aijcluclu(indx,indx1);
            }
          }
          pad->SetStatus(-8);
        } // for (Int_t j=0; j<npad;
      } // if (nCoupled > 3)
    } // while (nCoupled > 0)
  } // for (Int_t igroup=0; igroup<nclust;
  
  for (Int_t iclust=0; iclust<nclust; iclust++)
  {
    pix = clusters[iclust]; 
    pix->Clear();
    delete pix; 
    pix = 0;
  }
  delete [] clustNumb; 
  clustNumb = 0; 
  delete [] used; 
  used = 0;

}

//_____________________________________________________________________________
void 
AliMUONClusterSplitterMLEM::Merge(const AliMUONCluster& cluster,
                                     Int_t nForFit, Int_t nCoupled, 
                                     Int_t *clustNumb, Int_t *clustFit, 
                                     TObjArray **clusters, 
                                     TMatrixD& aijcluclu, TMatrixD& aijclupad)
{
  /// Merge the group of clusters with the one having the strongest coupling with them
  
  Int_t indx, indx1, npxclu, npxclu1, imax=0;
  TObjArray *pix, *pix1;
  Double_t couplMax;
  
  for (Int_t icl=0; icl<nForFit; ++icl) 
  {
    indx = clustFit[icl];
    pix = clusters[indx];
    npxclu = pix->GetEntriesFast();
    couplMax = -1;
    for (Int_t icl1=0; icl1<nCoupled; ++icl1) 
    {
      indx1 = clustNumb[icl1];
      if (indx1 < 0) continue;
      if ( aijcluclu(indx,indx1) > couplMax) 
      {
        couplMax = aijcluclu(indx,indx1);
        imax = indx1;
      }
    } // for (Int_t icl1=0;
      // Add to it
    pix1 = clusters[imax];
    npxclu1 = pix1->GetEntriesFast();
    // Add pixels 
    for (Int_t i=0; i<npxclu; ++i) 
    { 
      pix1->Add(pix->UncheckedAt(i)); 
      pix->RemoveAt(i); 
    }
    
    //Add cluster-to-cluster couplings
    for (Int_t icl1=0; icl1<nCoupled; ++icl1) 
    {
      indx1 = clustNumb[icl1];
      if (indx1 < 0 || indx1 == imax) continue;
      aijcluclu(indx1,imax) += aijcluclu(indx,indx1);
      aijcluclu(imax,indx1) = aijcluclu(indx1,imax);
    }
    aijcluclu(indx,imax) = aijcluclu(imax,indx) = 0;
    
    //Add cluster-to-pad couplings
    for (Int_t j=0; j<cluster.Multiplicity(); ++j) 
    {
      AliMUONPad* pad = cluster.Pad(j);
      if ( pad->Status() < 0 && pad->Status() != -5 ) continue;// exclude used pads
        aijclupad(imax,j) += aijclupad(indx,j);
        aijclupad(indx,j) = 0;
    }
  } // for (Int_t icl=0; icl<nForFit;
}


//_____________________________________________________________________________
Double_t 
AliMUONClusterSplitterMLEM::MinGroupCoupl(Int_t nCoupled, Int_t *clustNumb, 
                                          TMatrixD& aijcluclu, Int_t *minGroup)
{
  /// Find group of clusters with minimum coupling to all the others
  
  Int_t i123max = TMath::Min(3,nCoupled/2); 
  Int_t indx, indx1, indx2, indx3, nTot = 0;
  Double_t *coupl1 = 0, *coupl2 = 0, *coupl3 = 0;
  
  for (Int_t i123=1; i123<=i123max; i123++) {
    
    if (i123 == 1) {
      coupl1 = new Double_t [nCoupled];
      for (Int_t i=0; i<nCoupled; i++) coupl1[i] = 0;
    }
    else if (i123 == 2) {
      nTot = nCoupled*nCoupled;
      coupl2 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl2[i] = 9999;
    } else {
      nTot = nTot*nCoupled;
      coupl3 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl3[i] = 9999;
    } // else
    
    for (Int_t i=0; i<nCoupled; i++) {
      indx1 = clustNumb[i];
      for (Int_t j=i+1; j<nCoupled; j++) {
        indx2 = clustNumb[j];
        if (i123 == 1) {
          coupl1[i] += aijcluclu(indx1,indx2);
          coupl1[j] += aijcluclu(indx1,indx2);
        } 
        else if (i123 == 2) {
          indx = i*nCoupled + j;
          coupl2[indx] = coupl1[i] + coupl1[j];
          coupl2[indx] -= 2 * (aijcluclu(indx1,indx2));
        } else {
          for (Int_t k=j+1; k<nCoupled; k++) {
            indx3 = clustNumb[k];
            indx = i*nCoupled*nCoupled + j*nCoupled + k;
            coupl3[indx] = coupl2[i*nCoupled+j] + coupl1[k];
            coupl3[indx] -= 2 * (aijcluclu(indx1,indx3)+aijcluclu(indx2,indx3));
          }
        } // else
      } // for (Int_t j=i+1;
    } // for (Int_t i=0;
  } // for (Int_t i123=1;
  
  // Find minimum coupling
  Double_t couplMin = 9999;
  Int_t locMin = 0;
  
  for (Int_t i123=1; i123<=i123max; i123++) {
    if (i123 == 1) {
      locMin = TMath::LocMin(nCoupled, coupl1);
      couplMin = coupl1[locMin];
      minGroup[0] = locMin;
      delete [] coupl1; coupl1 = 0;
    } 
    else if (i123 == 2) {
      locMin = TMath::LocMin(nCoupled*nCoupled, coupl2);
      if (coupl2[locMin] < couplMin) {
        couplMin = coupl2[locMin];
        minGroup[0] = locMin/nCoupled;
        minGroup[1] = locMin%nCoupled;
      }
      delete [] coupl2; coupl2 = 0;
    } else {
      locMin = TMath::LocMin(nTot, coupl3);
      if (coupl3[locMin] < couplMin) {
        couplMin = coupl3[locMin];
        minGroup[0] = locMin/nCoupled/nCoupled;
        minGroup[1] = locMin%(nCoupled*nCoupled)/nCoupled;
        minGroup[2] = locMin%nCoupled;
      }
      delete [] coupl3; coupl3 = 0;
    } // else
  } // for (Int_t i123=1;
  return couplMin;
}

//_____________________________________________________________________________
Int_t 
AliMUONClusterSplitterMLEM::SelectPad(const AliMUONCluster& cluster,
                                          Int_t nCoupled, Int_t nForFit, 
                                          Int_t *clustNumb, Int_t *clustFit, 
                                          TMatrixD& aijclupad)
{
  /// Select pads for fit. If too many coupled clusters, find pads giving 
  /// the strongest coupling with the rest of clusters and exclude them from the fit.
  
  Int_t npad = cluster.Multiplicity();
  Double_t *padpix = 0;
  
  if (nCoupled > 3) 
  {
    padpix = new Double_t[npad];
    for (Int_t i=0; i<npad; i++) padpix[i] = 0; 
  }
  
  Int_t nOK = 0, indx, indx1;
  for (Int_t iclust=0; iclust<nForFit; ++iclust)
  {
    indx = clustFit[iclust];
    for (Int_t j=0; j<npad; j++) 
    {
      if ( aijclupad(indx,j) < fgkCouplMin) continue;
      AliMUONPad* pad = cluster.Pad(j);
      if ( pad->Status() == -5 ) pad->SetStatus(-0); // flag overflow
      if ( pad->Status() < 0 ) continue; // exclude overflows and used pads
      if ( !pad->Status() ) 
      {
        pad->SetStatus(1);
        ++nOK; // pad to be used in fit
      }      
      if (nCoupled > 3) 
      {
        // Check other clusters
        for (Int_t iclust1=0; iclust1<nCoupled; ++iclust1) 
        {
          indx1 = clustNumb[iclust1];
          if (indx1 < 0) continue;
          if ( aijclupad(indx1,j) < fgkCouplMin ) continue;
          padpix[j] += aijclupad(indx1,j);
        }
      } // if (nCoupled > 3)
    } // for (Int_t j=0; j<npad;
  } // for (Int_t iclust=0; iclust<nForFit
  if (nCoupled < 4) return nOK;
  
  Double_t aaa = 0;
  for (Int_t j=0; j<npad; ++j) 
  {
    if (padpix[j] < fgkCouplMin) continue;
    aaa += padpix[j];
    cluster.Pad(j)->SetStatus(-1); // exclude pads with strong coupling to the other clusters
    nOK--;
  }
  delete [] padpix; 
  padpix = 0;
  return nOK;
}

//_____________________________________________________________________________
void 
AliMUONClusterSplitterMLEM::UpdatePads(const AliMUONCluster& cluster,
                                          Int_t /*nfit*/, Double_t *par)
{
  /// Subtract the fitted charges from pads with strong coupling
  
  Int_t indx;
  Double_t charge, coef=0;
  
  for (Int_t j=0; j<cluster.Multiplicity(); ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    if ( pad->Status() != 1 ) continue;
    if (fNpar != 0) 
    {
      charge = 0;
      for (Int_t i=fNpar/3; i>=0; --i) 
      { 
        // sum over tracks
        indx = i<2 ? 2*i : 2*i+1;
        if (fNpar == 2) 
        {
          coef = 1;
        }
        else 
        {
          coef = i==fNpar/3 ? par[indx+2] : 1-coef;
        }
        coef = TMath::Max (coef, 0.);
        if (fNpar == 8 && i < 2) 
        {
          coef = i==1 ? coef*par[indx+2] : coef - par[7];
        }
        coef = TMath::Max (coef, 0.);
        charge += ChargeIntegration(par[indx],par[indx+1],*pad);
      }
      charge *= fQtot;
      pad->SetCharge(pad->Charge()-charge);
    } // if (fNpar != 0)
    
    if (pad->Charge() > 6 /*fgkZeroSuppression*/) pad->SetStatus(0); 
    // return pad for further using // FIXME: remove usage of zerosuppression here
    
  } // for (Int_t j=0;
}  


