// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTransformerLUT.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Histogram.h"
#include "AliL3HistogramAdaptive.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** /class AliL3HoughTransformerLUT
//<pre>
//_____________________________________________________________
// AliL3HoughTransformerLUT
//
// Hough transformation class using Look-UP-Tables
//
//</pre>
*/

ClassImp(AliL3HoughTransformerLUT)

AliL3HoughTransformerLUT::AliL3HoughTransformerLUT() : AliL3HoughBaseTransformer()
{
  fParamSpace=0;
  fMinRow=0;
  fMaxRow=0;

  fNRows=0;
  fNEtas=0;
  fNPhi0=0;
  fSector=0;
  fSlice=0;
  fSectorRow=0;
  fZSign=0;
  fZLengthPlusOff=0;
  fTimeWidth=0;
  fPadPitch=0;
  fEtaSlice=0;

  fLUTX=0;
  fLUTY=0;
  fLUTEta=0;
  fLUTEtaReal=0;
  fLUTphi0=0;
  fLUT2sinphi0=0;
  fLUT2cosphi0=0;
  fLUTKappa=0;
  
  fLastPad=0;
  fLastIndex=0;
  fAccCharge=0;
  fX=fY=0.;
}

AliL3HoughTransformerLUT::AliL3HoughTransformerLUT(Int_t slice,Int_t patch,Int_t n_eta_segments) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  fParamSpace=0;
  fMinRow=0;
  fMaxRow=0;

  fNRows=0;
  fNEtas=0;
  fNPhi0=0;
  fSector=0;
  fSectorRow=0;
  fZSign=0;
  fZLengthPlusOff=0;
  fTimeWidth=0;
  fPadPitch=0;
  fEtaSlice=0;

  fLUTX=0;
  fLUTY=0;
  fLUTEta=0;
  fLUTEtaReal=0;
  fLUTphi0=0;
  fLUT2sinphi0=0;
  fLUT2cosphi0=0;
  fLUTKappa=0;

  fLastPad=0;
  fLastIndex=0;
  fAccCharge=0;
  fX=fY=0.;

  Init(slice,patch,n_eta_segments);
}

AliL3HoughTransformerLUT::~AliL3HoughTransformerLUT()
{
  DeleteHistograms();

  if(fNRows){
    delete[] fLUTX;
    delete[] fLUTY;
    fNRows=0;
  }

  if(fNEtas){
    delete[] fLUTEta;
    delete[] fLUTEtaReal;
    fNEtas=0;
  }
}

void AliL3HoughTransformerLUT::DeleteHistograms()
{
  if(fNPhi0){
    delete[] fLUT2sinphi0;
    delete[] fLUT2cosphi0;
    delete[] fLUTphi0;
    delete[] fLUTKappa;
    fNPhi0=0;
    fLUTphi0=0;
    fLUT2sinphi0=0;
    fLUT2cosphi0=0;
    fLUTKappa=0;
  }

  if(!fParamSpace){
    for(Int_t i=0; i<fNEtas; i++)
      {
	if(!fParamSpace[i]) continue;
	delete fParamSpace[i];
      }
    delete [] fParamSpace;
  }
}

void AliL3HoughTransformerLUT::Init(Int_t slice,Int_t patch,Int_t n_eta_segments,Int_t /*n_seqs*/) 
{
  AliL3HoughBaseTransformer::Init(slice,patch,n_eta_segments);

  //delete old LUT tables
  if(fNRows){
    fNRows=0;
    delete[] fLUTX;
    delete[] fLUTY;
  }
  if(fNEtas){
    delete[] fLUTEta;
    delete[] fLUTEtaReal;
    fNEtas=0;
  }

  //member values
  fMinRow=AliL3Transform::GetFirstRow(patch);;
  fMaxRow=AliL3Transform::GetLastRow(patch);;
  fNRows=AliL3Transform::GetNRows(patch);;
  fNEtas=n_eta_segments;
  if(fNEtas!=GetNEtaSegments()) {
    cerr << "AliL3HoughTransformerLUT::Init -> Error: Number of Etas must be equal!" << endl;
    exit(1);
  }

  AliL3Transform::Slice2Sector(slice,fMinRow,fSector,fSectorRow);
  fZSign = slice < 18 ? 1:-1; //see AliL3Transformer
  fSlice = slice;
  fZLengthPlusOff=AliL3Transform::GetZLength()+AliL3Transform::GetZOffset();
  fTimeWidth=AliL3Transform::GetZWidth();

  fPadPitch=0;
  if(fSector<AliL3Transform::GetNSectorLow())
    fPadPitch=AliL3Transform::GetPadPitchWidthLow();
  else
    fPadPitch=AliL3Transform::GetPadPitchWidthUp();  

  Float_t etamax_=GetEtaMax();
  Float_t etamin_=GetEtaMin();
  fEtaSlice=(etamax_-etamin_)/n_eta_segments;

  //lookup tables for X and Y
  fLUTX=new Float_t[fNRows];
  fLUTY=new Float_t[fNRows]; 
  for(Int_t rr=0;rr<fNRows;rr++){
    fLUTX[rr]=Float_t(AliL3Transform::Row2X(rr+fMinRow));
    fLUTY[rr]=Float_t(0.5*(AliL3Transform::GetNPads(rr+fMinRow)-1)*fPadPitch);
    //cout << rr << ": " << (Float_t)fLUTX[rr] << " " << (Float_t)fLUTY[rr] << endl;
  }

  //lookup tables for rz2s <=> etas
  //will only be used ifdef FULLLUT 
  fLUTEta=new Float_t[fNEtas];
  fLUTEtaReal=new Float_t[fNEtas];
  for(Int_t rr=0;rr<fNEtas;rr++){
    Float_t eta=etamin_+(rr+1)*fEtaSlice;
    fLUTEta[rr]=CalcRoverZ2(eta);
    fLUTEtaReal[rr]=eta-0.5*fEtaSlice;
    //cout << rr << ": " << eta << " " << fLUTEtaReal[rr] << " " << GetEta(rr,fSlice) << " - " << fLUTEta[rr] << endl;
  }
}

void AliL3HoughTransformerLUT::CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,
					     Int_t nybin,Float_t psi)
{
  //Create histograms.
  //_Only_ to be used in case of the adaptive histograms!
  //phimax is given in radians!!
  
  if(ptmin > ptmax)
    {
      cerr<<"AliL3HoughTransformerLUT::CreateHistograms: Error in ptrange "<<ptmin<<" "<<ptmax<<endl;
      return;
    }
  if(psi < 0)
    {
      cerr<<"AliL3HoughTransformerLUT::CreateHistograms: Wrong psi-angle "<<psi<<endl;
      return;
    }
  
  fParamSpace = new AliL3Histogram*[fNEtas];
  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3HistogramAdaptive(histname,ptmin,ptmax,ptres,nybin,-psi,psi);
    }

  //create lookup table for sin and cos
  fNPhi0=nybin;
  fLUTphi0=new Float_t[fNPhi0];
  fLUT2sinphi0=new Float_t[fNPhi0];
  fLUT2cosphi0=new Float_t[fNPhi0];
  fLUTKappa=new Float_t[fNPhi0];
  AliL3Histogram *hist=fParamSpace[0];
  Int_t i=0;
  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
    {
      Float_t phi0 = hist->GetBinCenterY(b);
      fLUTphi0[i]=phi0;
      fLUT2sinphi0[i]=2.*sin(phi0);
      fLUT2cosphi0[i]=2.*cos(phi0);
      fLUTKappa[i]=0.;
      i++;
      //cout << i << ": " << fLUTphi0[i] << " " << fLUT2sinphi0[i] << " " << fLUT2cosphi0[i] << endl;
    }
}

void AliL3HoughTransformerLUT::CreateHistograms(Int_t nxbin,Float_t pt_min,Int_t nybin,Float_t phimin,Float_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //pt_min = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 (degrees)
  //phi_max = maximum phi0 (degrees)
    
  Double_t x = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/pt_min;
  //Double_t torad = AliL3Transform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin/**torad*/,phimax/**torad*/);
}

void AliL3HoughTransformerLUT::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax)
{
  fParamSpace = new AliL3Histogram*[fNEtas];
  
  Char_t histname[256];
  for(Int_t i=0; i<fNEtas; i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
  
  //create lookup table for sin and cos
  fNPhi0=nybin;
  fLUTphi0=new Float_t[fNPhi0];
  fLUT2sinphi0=new Float_t[fNPhi0];
  fLUT2cosphi0=new Float_t[fNPhi0];
  fLUTKappa=new Float_t[fNPhi0];

  AliL3Histogram *hist=fParamSpace[0];
  Int_t i=0;
  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
    {
      Float_t phi0 = hist->GetBinCenterY(b);
      fLUTphi0[i]=phi0;
      fLUT2sinphi0[i]=2.*sin(phi0);
      fLUT2cosphi0[i]=2.*cos(phi0);
      fLUTKappa[i]=0.; 
      //cout << i << ": " << fLUTphi0[i] << " " << fLUT2sinphi0[i] << " " << fLUT2cosphi0[i] << endl;
      i++;
    }
}

void AliL3HoughTransformerLUT::Reset()
{
  //Reset all the histograms

  if(!fParamSpace)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformer::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<fNEtas; i++)
    fParamSpace[i]->Reset();
}

Int_t AliL3HoughTransformerLUT::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 
  
#ifdef FULLLUT 
  /* try to imitate a circuit -> should 
     go into the VHDL implementation of transformer */
  Float_t rz2=CalcRoverZ2(eta);
  return FindIndex(rz2);
#else /* optimize for speed on the computer */
  Double_t index = (eta-GetEtaMin())/fEtaSlice;
  return (Int_t)index;
#endif
}

AliL3Histogram *AliL3HoughTransformerLUT::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= fNEtas || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Double_t AliL3HoughTransformerLUT::GetEta(Int_t eta_index,Int_t slice)
{
  if(eta_index >= fNEtas || eta_index < 0){
    LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::GetEta","Index") << "Index out of range."<<ENDLOG;
    return 0.;
  }
  if(slice != fSlice){
    LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::GetEta","Index") << "Given slice does not match internal slice."<<ENDLOG;
    return 0.;
  }

  return(fLUTEtaReal[eta_index]);
}

void AliL3HoughTransformerLUT::TransformCircle()
{
  //Transform the input data with a circle HT.
  //The function loops over all the data, and transforms each pixel with the equation:
  // 
  //kappa = 2/R^2 * (y*cos(phi0) - x*sin(phi0) )
  //
  //where R^2 = x^2 + y^2
  //
  //After a day of testing it is proven that this h

  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformerLUT::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }

  Int_t lowch=GetLowerThreshold();
  Int_t highch=GetUpperThreshold();

  //Loop over the padrows:
  for(Int_t i=fMinRow, row=0; i<=fMaxRow; i++, row++)
    {
      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
	  LOG(AliL3Log::kError,"AliL3HoughTransformerLUT::TransformCircle","Data")
	    <<"AliL3HoughTransformerLUT::TransformCircle : Mismatching padrow numbering "<<i<<" != "<<(Int_t)tempPt->fRow<<ENDLOG;
	  continue;
	}


      //start a new row
      fLastPad=-1;
      //store x for this row
      fX = CalcX(row);
      //accumulate charge per etaslice
      fAccCharge=0;
      //last histogram
      fLastIndex=-1; 

      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;

	  //check threshold
	  if((Int_t)charge <=  lowch)
	    continue;
	  if ((Int_t)charge > highch)
	    charge=highch;

	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;

	  if(fLastPad!=pad){ //only update if necessary

	    //if there is accumalated charge, 
	    //put it in the old histogram
	    //using the old X and Y values
	    if(fAccCharge>0) AddCurveToHistogram(-1); //fLastIndex will be set to -1

#ifdef FULLLUT
	    fLastIndex=fNEtas-1;
#endif

	    //calculate Y for the new pad
	    fY = CalcY(pad,row);
	    //remember new pad
	    fLastPad=pad;
	  }

#ifdef FULLLUT
	  //find eta slice
	  Float_t z = CalcZ(time);

	  Float_t z2=z*z;
	  Float_t rz2 = 1 + r2/z2;
	  Int_t eta_index = FindIndex(rz2,fLastIndex);
#else
	  Float_t z = CalcZ(time);
	  Double_t r = sqrt(fX*fX+fY*fY+z*z);
	  Double_t eta = 0.5 * log((r+z)/(r-z));
	  Int_t eta_index = GetEtaIndex(eta);
#endif
	  if(eta_index < 0 || eta_index >= fNEtas){
	    LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::TransformCircle","Histograms")<<"No histograms corresponding to eta index value of "<<eta_index<<"."<<ENDLOG;
	    continue;
	  }	  

#ifndef FULLLUT
	  if(fLastIndex==-1) fLastIndex=eta_index;
#endif

	  if(fLastIndex!=eta_index){ //enter old values first
	    AddCurveToHistogram(eta_index);
	  }

	  fAccCharge+=charge;
	  //fAccCharge+=1;
	}

      //in case there is charge, store it before moving to the next row
      if(fAccCharge>0) AddCurveToHistogram(-1);
      
      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
}

void AliL3HoughTransformerLUT::Print()
{
  cout << "fSlice: " << GetSlice() << endl;
  cout << "fPatch: " << GetPatch() << endl;
  cout << "fSector: " << fSector << endl;
  cout << "fSectorRow: " << fSectorRow << endl;
  cout << "fMinRow: " << fMinRow << endl;
  cout << "fMaxRow: " << fMaxRow << endl;
  cout << "fNRows: " << fNRows << endl;
  cout << "fNEtas: " << fNEtas << endl;
  cout << "fNPhi0: " << fNPhi0 << endl;
  cout << "fZSign: " << fZSign << endl;
  cout << "fZLengthPlusOff: " << fZLengthPlusOff << endl;
  cout << "fPadPitch: " << fPadPitch << endl;
  cout << "fTimeWidth: " << fTimeWidth << endl;

  if(!fNRows) return;
  cout << "fLUTX " << fNRows << endl;
  for(Int_t i=0;i<fNRows;i++) cout << "fLUTX[" << i << "]=" << (Float_t)fLUTX[i] << endl;
  cout << "fLUTY " << fNRows << endl;
  for(Int_t i=0;i<fNRows;i++) cout << "fLUTY[" << i << "]=" << fLUTY[i] << endl;
  if(!fNEtas) return;
  cout << "fLUTEta " << fNEtas << endl;
  for(Int_t i=0;i<fNEtas;i++) cout << "fLUTEta[" << i << "]=" << fLUTEta[i] << endl;
  cout << "fLUTEtaReal " << fNEtas << endl;
  for(Int_t i=0;i<fNEtas;i++) cout << "fLUTEtaReal[" << i << "]=" << fLUTEtaReal[i] << endl;
  if(!fNPhi0) return;
  cout << "fLUTphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUTPhi0[" << i << "]=" << fLUTphi0[i] << endl;
  cout << "fLUT2sinphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2sinphi0[" << i << "]=" << fLUT2sinphi0[i] << endl;
  cout << "fLUT2cosphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2cosphi0[" << i << "]=" << fLUT2cosphi0[i] << endl;
}
