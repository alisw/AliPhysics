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

#if GCCVERSION == 3
using namespace std;
#endif

//#define FULLLUT

//_____________________________________________________________
// AliL3HoughTransformerLUT
//
// Hough transformation class using Look-UP-Tables
//

ClassImp(AliL3HoughTransformerLUT)

AliL3HoughTransformerLUT::AliL3HoughTransformerLUT() : AliL3HoughBaseTransformer()
{
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

  fDoMC = kFALSE;
}

AliL3HoughTransformerLUT::AliL3HoughTransformerLUT(Int_t slice,Int_t patch,Int_t n_eta_segments) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
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
  fDoMC=kFALSE;

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

void AliL3HoughTransformerLUT::Init(Int_t slice,Int_t patch,Int_t n_eta_segments,Int_t n_seqs) 
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
  fZSign = slice < 18 ? 1:-1;
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
  fLUTEta=new Float_t[fNEtas];
  fLUTEtaReal=new Float_t[fNEtas];
  for(Int_t rr=0;rr<fNEtas;rr++){
    Float_t eta=etamin_+(rr+1)*fEtaSlice;
    fLUTEta[rr]=CalcRoverZ2(eta);
    fLUTEtaReal[rr]=eta-0.5*fEtaSlice;
    //cout << rr << ": " << eta << " " << fLUTEtaReal[rr] << " " << GetEta(rr,fSlice) << " - " << fLUTEta[rr] << endl;
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
  Double_t torad = AliL3Transform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin*torad,phimax*torad);
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
  fNPhi0=nybin+1;

  fLUTphi0=new Float_t[fNPhi0];
  fLUT2sinphi0=new Float_t[fNPhi0];
  fLUT2cosphi0=new Float_t[fNPhi0];
  fLUTKappa=new Float_t[fNPhi0];
  Float_t diff=(ymax-ymin)/nybin;
  Float_t phi0=ymin-0.5*diff;
  for(Int_t i=0; i<fNPhi0; i++){
    phi0+=diff;
    fLUTphi0[i]=phi0;
    fLUT2sinphi0[i]=2.*sin(phi0);
    fLUT2cosphi0[i]=2.*cos(phi0);
    fLUTKappa[i]=0.;
    //cout << i << ": " << fLUTphi0[i] << " " << fLUT2sinphi0[i] << " " << fLUT2cosphi0[i] << endl;
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

Int_t AliL3HoughTransformerLUT::FindIndex(Float_t rz2, Int_t start)
{
  //could improve search through devide and conquere strategy
  
  Int_t index=start; 
  if(index==-100){
    index=0;
    while((index<fNEtas)&&(rz2<=fLUTEta[index])){
      index++;
    }
  } else {
    while((index>=0)&&(rz2>fLUTEta[index])){
      index--;
    }
    index++;
  }
  //cout << start << " - " << index << " " << ": " << rz2 << " " << fLUTEta[index] << endl;

  return index;
}

inline AliL3Histogram *AliL3HoughTransformerLUT::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= fNEtas || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Float_t AliL3HoughTransformerLUT::CalcRoverZ2(Float_t eta)
{
  Float_t e=exp(2*eta);
  Float_t ret=(e+1)/(e-1);
  ret*=ret;
  return ret;
}

Float_t AliL3HoughTransformerLUT::CalcEta(Float_t roverz2)
{
  Float_t rz=sqrt(roverz2);
  if(fZSign<0) rz=-rz;
  Float_t ret=(1+rz)/(rz-1);
  ret=0.5*log(ret);
  return ret;
}

Float_t AliL3HoughTransformerLUT::CalcX(Int_t row)
{
  return fLUTX[row];
}

Float_t AliL3HoughTransformerLUT::CalcY(Int_t pad,Int_t row)
{
  return pad*fPadPitch-fLUTY[row];
}

Float_t AliL3HoughTransformerLUT::CalcZ(Int_t time)
{
  Float_t ret=time*fTimeWidth;
  if(fZSign>0) ret=fZLengthPlusOff-ret;
  else ret=ret-fZLengthPlusOff;
  return ret;
}

Double_t AliL3HoughTransformerLUT::GetEta(Int_t eta_index,Int_t slice)
{
  if(eta_index >= fNEtas || eta_index < 0){
    //LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::GetEta","Index") << "Index out of range."<<ENDLOG;
    return 0.;
  }
  if(slice != fSlice){
    //LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::GetEta","Index") << "Given slice does not match internal slice."<<ENDLOG;
    return 0.;
  }

  //return (CalcEta(fLUTEta[eta_index])-0.5*fEtaSlice);
  return(fLUTEtaReal[eta_index]);
}

void AliL3HoughTransformerLUT::AddCurveToHistogram(Int_t new_eta_index)
{
  //Get the correct histogrampointer:
  AliL3Histogram *hist = fParamSpace[fLastIndex];
  if(!hist){
    //LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::TransformCircle","Histograms")<<"Error getting histogram in index "<<fLastIndex<<"."<<ENDLOG;
    return;
  }

  //Fill the histogram along the phirange
  for(Int_t b=0; b<fNPhi0; b++){
    hist->Fill(fLUTKappa[b],fLUTphi0[b],fAccCharge);
    //cout << kappa << " " << fLUTphi0[b] << " " << fAccCharge << endl;
  }

  fAccCharge=0;
  fLastIndex=new_eta_index;
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
  //Each pixel then transforms into a curve in the (kappa,phi0)-space. In order to find
  //which histogram in which the pixel should be transformed, the eta-value is calcluated
  //and the proper histogram index is found by GetEtaIndex(eta).

  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformerLUT::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
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

      //calculate x for this row
      Float_t x = CalcX(row);
      Float_t x2=x*x;
      Float_t r2=0;

      //start a new row
      fLastPad=-1;

      if(fAccCharge>0) cerr << "Big error " << endl;

      //accumulate charge per histogram
      fAccCharge=0;
      fLastIndex=-1;

      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;

	  //check threshold
	  if((Int_t)charge <= GetLowerThreshold() || (Int_t)charge > GetUpperThreshold())
	    continue;

	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;

	  if(fLastPad!=pad){ //only update if necessary

	    //if there is accumalated charge, put it in the right histogram
	    if(fAccCharge>0) AddCurveToHistogram(-1);

#ifdef FULLLUT
	    fLastIndex=fNEtas-1;
#endif

	    //calculate hough for the new pad
	    Float_t y = CalcY(pad,row);
	    Float_t y2 = y*y;
	    r2 = x2 + y2;
	    Float_t xr2=x/r2;
	    Float_t yr2=y/r2;

	    //AliL3Histogram *hist = fParamSpace[0];
	    //Int_t bb=hist->GetFirstYbin();
	    //if(fNPhi0-1!=hist->GetLastYbin()) cerr << "should not be" << endl;
	    for(Int_t b=0; b<fNPhi0; b++){
	      fLUTKappa[b]=(yr2*fLUT2cosphi0[b]-xr2*fLUT2sinphi0[b]);

	      //Float_t phi=atan2(y,x);
	      //Float_t phi0 = hist->GetBinCenterY(bb);
	      //Float_t R=sqrt(r2);
	      //Float_t kappa = 2*sin(phi - phi0)/R;
	      //cout << fLUTKappa[b] << " " << kappa << endl;
	      //bb++;
	      //Int_t bina=hist->FindBin(fLUTphi0[b],fLUTKappa[b]);
	      //Int_t binb=hist->FindBin(phi0,kappa);
	      //if(bina!=binb) cout << bina << " " << binb << endl;
	    }
	    
	    fLastPad=pad;
	  }

	  //find eta slice
	  Float_t z = CalcZ(time);
	  Float_t z2=z*z;

#ifdef FULLLUT
	  Float_t rz2 = 1 + r2/z2;
	  Int_t eta_index = FindIndex(rz2,fLastIndex);
#else
	  Double_t r = sqrt(r2+z2);
	  Double_t etaval = 0.5 * log((r+z)/(r-z));
	  Int_t eta_index = GetEtaIndex(etaval);
#endif
	  //cout << row << " " << (int)pad << " " << (int)time << " " << eta_index << " " << fLastIndex << endl;
	  if(eta_index < 0 || eta_index >= fNEtas){
	    //LOG(AliL3Log::kWarning,"AliL3HoughTransformerLUT::TransformCircle","Histograms")<<"No histograms corresponding to eta index value of "<<eta_index<<"."<<ENDLOG;
	    continue;
	  }	  

#ifndef FULLLUT
	  if(fLastIndex==-1) fLastIndex=eta_index;
#endif

	  if(fLastIndex!=eta_index){
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

      if(fAccCharge>0) cerr << "Big error " << endl;
}

Int_t AliL3HoughTransformerLUT::GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi)
{
  if(fDoMC)
    {
      cerr<<"AliL3HoughTransformerLUT does not provide MC information..."<<endl;
      return -1;
    }
  
  return -1;
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
