//$Id$

// Author: Constantin Loizides <mailto:loizides@fi.uib.no>
//*-- Copyright&Copy CL

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3Histogram.h"

//try to be close to VHDL version (eg. LUTs)
//VESTBO: if you switch that of, you should get your version
//of the transformer!!!
#define VHDLVERSION

//dont use dynamic objects on heap
//VESTBO: this a desperately switched on in order 
//to find the bug, which i thought hat something
//to do with a wrong heap assignment! But the
//data still gets overwritten!
#define VHDLSTATIC

#include "AliL3HoughTransformerVhdl.h"

#if GCCVERSION == 3
using namespace std;
#endif

/** \class AliL3HoughTransformerVhdl
// <pre>
//_____________________________________________________________
// AliL3HoughTransformerVhdl
//
// Hough transformation class for VHDL comparism.
//
//</pre>
*/

ClassImp(AliL3HoughTransformerVhdl)

#ifdef VHDLVERSION

AliL3HoughTransformerVhdl::AliL3HoughTransformerVhdl()
  : AliL3HoughBaseTransformer()//, /*fMinRow(0),*/ fMaxRow(0)

{
  cout << "hallo" << endl;
  fParamSpace=0;

#ifndef VHDLSTATIC
  fLUTX=0;
  fLUTY=0;
  fLUTEta=0;
  fLUTphi0=0;
  fLUT2sinphi0=0;
  fLUT2cosphi0=0;
#endif
  
  fMinRow=0;
  fMaxRow=0;

  /*
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
  */
}

AliL3HoughTransformerVhdl::AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments) 
  : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  AliL3HoughTransformerVhdl();
  cout << "wiesow hier??" << endl;
  Init(slice,patch,n_eta_segments);
}

AliL3HoughTransformerVhdl::~AliL3HoughTransformerVhdl()
{
  DeleteHistograms();

#ifndef VHDLSTATIC
  if(fNRows){
    delete[] fLUTX;
    delete[] fLUTY;
  }
  if(fNEtas) delete[] fLUTEta;
#endif
}

void AliL3HoughTransformerVhdl::Init(Int_t slice,Int_t patch,Int_t n_eta_segments)
{
  cout << "InitVhdl " << slice << " " << patch << " " << n_eta_segments << endl;
  AliL3HoughBaseTransformer::Init(slice,patch,n_eta_segments);

#ifndef VHDLSTATIC
  //delete old LUT tables
  if(fNRows){
    fNRows=0;
    delete[] fLUTX;
    delete[] fLUTY;
  }
  if(fNEtas){
    delete[] fLUTEta;
    fNEtas=0;
  }
#endif

  //set when histogram is filled
  //fLUTphi0=0;
  //fLUT2sinphi0=0;
  //fLUT2cosphi0=0;

  Int_t minrow_=AliL3Transform::GetFirstRow(patch);
  Int_t maxrow_=AliL3Transform::GetLastRow(patch);
  Int_t n_=AliL3Transform::GetNRows(patch);
  Int_t sector_=0,sectorrow_=0;
  AliL3Transform::Slice2Sector(slice,minrow_,sector_,sectorrow_);
  Float_t padpitch_=0.;
  if(sector_<AliL3Transform::GetNSectorLow())
    padpitch_=AliL3Transform::GetPadPitchWidthLow();
  else
    padpitch_=AliL3Transform::GetPadPitchWidthUp();  

  Float_t etamax_=GetEtaMax();
  Float_t etamin_=GetEtaMin();
  Float_t etaslice_=(etamax_-etamin_)/n_eta_segments;
#if 0
  //lookup tables for X and Y
#ifndef VHDLSTATIC
  fLUTX=new AliL3FFloat[n_];
  fLUTY=new AliL3FFloat[n_]; 
#endif
  for(Int_t rr=0;rr<n_;rr++){
    fLUTX[rr]=Float_t(AliL3Transform::Row2X(rr+minrow_));
    fLUTY[rr]=Float_t(0.5*(AliL3Transform::GetNPads(rr+minrow_)-1)*padpitch_);
    //VESTBO: uncomment to see values (and compare with Print function)
    //cout << rr << ": " << (Float_t)fLUTX[rr] << " " << (Float_t)fLUTY[rr] << endl;
  }

  //lookup tables for rz2s <=> etas
#ifndef VHDLSTATIC
  fLUTEta=new AliL3FFloat[n_eta_segments];
#endif
  for(Int_t rr=0;rr<n_eta_segments;rr++){
    fLUTEta[rr]=CalcRoverZ2(etamin_+(rr+1)*etaslice_);
    //VESTBO: uncomment to see values (and compare with Print function)
    //cout << rr << ": " << fLUTEta[rr] << endl;
  }
#endif /*uncomment arrays*/

  //member values
  fMinRow=minrow_;
  fMaxRow=maxrow_;
  fNRows=n_;
  fNEtas=n_eta_segments;
  fNPhi0=0;
  fSector=sector_;
  fSectorRow=sectorrow_;
  fZSign = slice < 18 ? 1:-1;
  fZLengthPlusOff=AliL3Transform::GetZLength()+AliL3Transform::GetZOffset();
  fTimeWidth=AliL3Transform::GetZWidth();
  fPadPitch=padpitch_;
  fEtaSlice=etaslice_;
}

#if 0
Float_t AliL3HoughTransformerVhdl::CalcRoverZ2(Float_t eta)
{
  Float_t e=exp(2*eta);
  Float_t ret=(e+1)/(e-1);
  ret*=ret;
  return ret;
}

Float_t AliL3HoughTransformerVhdl::CalcEta(Float_t roverz2)
{
  Float_t rz=sqrt(roverz2);
  if(fZSign>0) rz=-rz;
  Float_t ret=(1+rz)/(rz-1);
  ret=0.5*log(ret);
  return ret;
}

Float_t AliL3HoughTransformerVhdl::CalcX(Int_t row)
{
  return fLUTX[row];
}

Float_t AliL3HoughTransformerVhdl::CalcY(Int_t pad,Int_t row)
{
  return pad*fPadPitch-fLUTY[row];
}

Float_t AliL3HoughTransformerVhdl::CalcZ(Int_t time)
{
  Float_t ret=time*fTimeWidth;
  if(fZSign>0) ret=fZLengthPlusOff-ret;
  else ret=ret-fZLengthPlusOff;
  return ret;
}
#endif

void AliL3HoughTransformerVhdl::DeleteHistograms()
{
  if(!fParamSpace) return;
  for(Int_t i=0; i<GetNEtaSegments(); i++){
    if(!fParamSpace[i]) continue;
    delete fParamSpace[i];
  }
  delete[] fParamSpace;

#ifndef VHDLSTATIC
  if(fNPhi0){
    fNPhi0=0;
    delete[] fLUT2sinphi0;
    delete[] fLUT2cosphi0;
    delete[] fLUTphi0;
  }
#endif
}

void AliL3HoughTransformerVhdl::CreateHistograms(Int_t nxbin,Double_t pt_min,Int_t nybin,Double_t phimin,Double_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //pt_min = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 (degrees)
  //phi_max = maximum phi0 (degrees)
    
  AliL3FFloat x = AliL3Transform::GetBFieldValue()/pt_min;
  AliL3FFloat torad = AliL3Transform::Pi()/180;
  AliL3FFloat phimin_=phimin*torad;
  AliL3FFloat phimax_=phimax*torad;

  CreateHistograms(nxbin,-x,x,nybin,phimin_,phimax_);
}

void AliL3HoughTransformerVhdl::CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
					         Int_t nybin,Double_t ymin,Double_t ymax)
{

  AliL3FFloat xminf=xmin;
  AliL3FFloat xmaxf=xmax;
  AliL3FFloat yminf=ymin;
  AliL3FFloat ymaxf=ymax;
  /*
  cout << xminf << endl;
  cout << xmaxf << " " << ((xmaxf-xminf)/nxbin) << endl;
  cout << yminf << endl;
  cout << ymaxf << " " << ((ymaxf-yminf)/nybin) << endl;
  */


  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];

  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_vhdl_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"id",nxbin,xminf,xmaxf,nybin,yminf,ymaxf);
    }
  
  //create lookup table for sin and cos
  fNPhi0=nybin+1;
#if 0
#ifndef VHDLSTATIC
  fLUTphi0=new AliL3FFloat[fNPhi0];
  fLUT2sinphi0=new AliL3FFloat[fNPhi0];
  fLUT2cosphi0=new AliL3FFloat[fNPhi0];
#endif
  AliL3FFloat diff=(ymax-ymin)/nybin;
  AliL3FFloat phi0=ymin-0.5*diff;
  for(Int_t i=0; i<fNPhi0; i++){
    phi0+=diff;
    fLUTphi0[i]=phi0;
    fLUT2sinphi0[i]=Float_t(2*sin(phi0));
    fLUT2cosphi0[i]=Float_t(2*cos(phi0));
    //VESTBO: uncomment to see values (and compare with Print function)
    //cout << i << ": " << fLUTphi0[i] << " " << fLUT2sinphi0[i] << " " << fLUT2cosphi0[i] << endl;
  }  
#endif /*arrays*/
}

void AliL3HoughTransformerVhdl::Reset()
{
  //Reset all the histograms

  if(!fParamSpace)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fParamSpace[i]->Reset();
}

AliL3Histogram *AliL3HoughTransformerVhdl::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0) return 0;
  if(!fParamSpace[eta_index]) return 0;
  return fParamSpace[eta_index];
}

#if 0
Int_t AliL3HoughTransformerVhdl::FindIndex(Double_t rz2)
{
  Int_t index=0;
  while((index<fNEtas)&&(rz2<=fLUTEta[index])){
    index++;
    //cout << index << ": " << rz2 << " " << fLUTEta[index] << endl;
  }
  return index;
}
#endif

Int_t AliL3HoughTransformerVhdl::GetEtaIndex(Double_t eta)
{
  //AliL3FFloat rz2=CalcRoverZ2(eta);
  //return FindIndex(rz2);
  return 0;
}

Double_t AliL3HoughTransformerVhdl::GetEta(Int_t eta_index,Int_t slice){
  if(eta_index >= fNEtas || eta_index < 0){
    LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::GetEta","Index")
      <<"Index out of range."<<ENDLOG;
    return 0.;
  }
  return 0;//(CalcEta(fLUTEta[eta_index])-0.5*fEtaSlice);
}

void AliL3HoughTransformerVhdl::TransformCircle()
{
  //Transform the input data with a circle HT.
  //The function loops over all the data, and transforms each pixel with the equations:
  // 
  //kappa = 1/x * ( y/x*2*cos(phi0) - 2*sin(phi0) )

#if 0
  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformerVhdl::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }

  AliL3FFloat x=0.,y=0.,z=0.,rz2=0.;
  //Loop over the padrows:
   for(Int_t i=fMinRow, row=0; i<=fMaxRow; i++, row++){

    //Get the data on this padrow:
    AliL3DigitData *digPt = tempPt->fDigitData;
    if(i != (Int_t)tempPt->fRow){
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Padrows")
	<<"Mismatching padrow numbering "<<i<<"!="<<(Int_t)tempPt->fRow<<ENDLOG;
      continue;
    }
      
    //Loop over the data on this padrow:
    for(UInt_t j=0; j<tempPt->fNDigit; j++){
      Int_t charge = digPt[j].fCharge;
      if((Int_t)charge <= GetLowerThreshold()) continue;
      
      Int_t pad = digPt[j].fPad;
      Int_t time = digPt[j].fTime;

      x=CalcX(row);
      y=CalcY(pad,row);
      z=CalcZ(time);

      //find eta slice
      rz2=1+(x*x+y*y)/z/z;
      Int_t eta_index = FindIndex(rz2);
      if(eta_index < 0 || eta_index >= GetNEtaSegments()){
	//LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Histograms")<<"No histograms corresponding to eta index value of "<<eta_index<<"."<<ENDLOG;
	continue;
      }	  
      //Get the correct histogrampointer:
      AliL3Histogram *hist = fParamSpace[eta_index];
      if(!hist){
	//LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Histograms")<<"Error getting histogram in index "<<eta_index<<"."<<ENDLOG;
	continue;
      }

      //Fill the histogram along the phirange
      AliL3FFloat kappa,ydx=y/x;
      for(Int_t b=0; b<fNPhi0; b++){
	kappa=ydx*fLUT2cosphi0[b]-fLUT2sinphi0[b];
	kappa/=x;
	hist->Fill(kappa,fLUTphi0[b],charge);
	//cout << kappa << " " << fLUTphi0[b] << " " << charge << endl;
      }

    }
    //Move the data pointer to the next padrow:
    AliL3MemHandler::UpdateRowPointer(tempPt);
  }
#endif 
}

void AliL3HoughTransformerVhdl::Print()
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
#if 0
  if(!fNRows) return;
  cout << "fLUTX " << fNRows << endl;
  for(Int_t i=0;i<fNRows;i++) cout << "fLUTX[" << i << "]=" << (Float_t)fLUTX[i] << endl;
  cout << "fLUTY " << fNRows << endl;
  for(Int_t i=0;i<fNRows;i++) cout << "fLUTY[" << i << "]=" << fLUTY[i] << endl;
  if(!fNEtas) return;
  cout << "fLUTEta " << fNEtas << endl;
  for(Int_t i=0;i<fNEtas;i++) cout << "fLUTEta[" << i << "]=" << fLUTEta[i] << endl;
  if(!fNPhi0) return;
  cout << "fLUTphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUTPhi0[" << i << "]=" << fLUTphi0[i] << endl;
  cout << "fLUT2sinphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2sinphi0[" << i << "]=" << fLUT2sinphi0[i] << endl;
  cout << "fLUT2cosphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2cosphi0[" << i << "]=" << fLUT2cosphi0[i] << endl;
#endif 
}

//end vhdl version
#else 
//standard version

AliL3HoughTransformerVhdl::AliL3HoughTransformerVhdl()
{
  fParamSpace=0;
}

AliL3HoughTransformerVhdl::AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments) 
            : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  fParamSpace=0;
}

AliL3HoughTransformerVhdl::~AliL3HoughTransformerVhdl()
{
  DeleteHistograms();
}

void AliL3HoughTransformerVhdl::DeleteHistograms()
{
  if(!fParamSpace)
    return;
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      if(!fParamSpace[i]) continue;
      delete fParamSpace[i];
    }
  delete[] fParamSpace;
}

void AliL3HoughTransformerVhdl::CreateHistograms(Int_t nxbin,Double_t pt_min,Int_t nybin,Double_t phimin,Double_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //pt_min = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 (degrees)
  //phi_max = maximum phi0 (degrees)
    
  Double_t x = AliL3Transform::GetBFieldValue()/pt_min;
  Double_t torad = AliL3Transform::Pi()/180;

  CreateHistograms(nxbin,-1.*x,x,nybin,phimin*torad,phimax*torad);
}

void AliL3HoughTransformerVhdl::CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
					         Int_t nybin,Double_t ymin,Double_t ymax)
{

  AliL3FFloat xminf(xmin);
  AliL3FFloat xmaxf(xmax);
  AliL3FFloat yminf(ymin);
  AliL3FFloat ymaxf(ymax);

  cout << xminf << endl;
  cout << xmaxf << " " << ((xmaxf-xminf)/nxbin) << endl;
  cout << yminf << endl;
  cout << ymaxf << " " << ((ymaxf-yminf)/nybin) << endl;

  
  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];
  
  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_vhdl_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
}

void AliL3HoughTransformerVhdl::Reset()
{
  //Reset all the histograms

  if(!fParamSpace)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fParamSpace[i]->Reset();
}

AliL3Histogram *AliL3HoughTransformerVhdl::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0) return 0;
  if(!fParamSpace[eta_index]) return 0;
  return fParamSpace[eta_index];
}

Int_t AliL3HoughTransformerVhdl::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 
  AliL3FFloat etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  AliL3FFloat index = (eta-GetEtaMin())/etaslice;

  //if((Int_t)index==0) cout << index << " " << (Int_t) index << endl;
  return (Int_t)index;
}

Double_t AliL3HoughTransformerVhdl::GetEta(Int_t eta_index)
{
  if(eta_index >= GetNEtaSegments() || eta_index < 0){
    return 0;
  }
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=(Double_t)((eta_index+0.5)*eta_slice);
  if(GetSlice()>17) eta*=-1;
  return eta;
}

void AliL3HoughTransformerVhdl::TransformCircle()
{
  //Transform the input data with a circle HT.
  //The function loops over all the data, and transforms each pixel with the equations:
  // 
  //kappa = 2/R*sin(phi - phi0)
  //
  //where R = sqrt(x*x +y*y), and phi = arctan(y/x)
  //
  //Each pixel then transforms into a curve in the (kappa,phi0)-space. In order to find
  //which histogram in which the pixel should be transformed, the eta-value is calcluated
  //and the proper histogram index is found by GetEtaIndex(eta).


  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformerVhdl::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
  //Loop over the padrows:
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Padrows")
	<<"Mismatching padrow numbering."<<ENDLOG;
	  continue;
	}
      

      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  if((Int_t)charge <= GetLowerThreshold()) continue;

	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;

	  Int_t sector,row;
	  Float_t xyz[3];
	  AliL3FFloat fxyz[3];

	  //Transform data to local cartesian coordinates:
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);

	  //trunc to fixed format.
	  for(int i=0;i<3;i++) {
	    fxyz[i]=xyz[i];
	    //cout << fxyz[i];
	    xyz[i]=fxyz[i];
	  }
	  //cout << endl;

	  //Calculate the eta:
	  AliL3FFloat eta = AliL3Transform::GetEta(xyz);
	  
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = GetEtaIndex(eta);
	  if(eta_index < 0 || eta_index >= GetNEtaSegments()){
	    //LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Histograms")<<"No histograms corresponding to eta index value of "<<eta_index<<"."<<ENDLOG;
	    continue;
	  }	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      //LOG(AliL3Log::kWarning,"AliL3HoughTransformerVhdl::TransformCircle","Histograms")<<"Error getting histogram in index "<<eta_index<<"."<<ENDLOG;
	      continue;
	    }

	  //Do the transformation:
	  AliL3FFloat R = sqrt(fxyz[0]*fxyz[0] + fxyz[1]*fxyz[1]); 
	  AliL3FFloat phi = AliL3Transform::GetPhi(xyz);
	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      AliL3FFloat phi0 = hist->GetBinCenterY(b);
	      AliL3FFloat kappa = 2*sin(phi - phi0)/R;
	      //hist->Fill(kappa.GetExactVal(),phi0.GetExactVal(),charge);
	      hist->Fill(kappa,phi0,charge);
	    }
	}
      
      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
}

#endif
