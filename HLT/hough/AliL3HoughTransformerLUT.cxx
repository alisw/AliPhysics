//$Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright & copy CL

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
  fSectorRow=0;
  fZSign=0;
  fZLengthPlusOff=0;
  fTimeWidth=0;
  fPadPitch=0;
  fEtaSlice=0;

  fLUTX=0;
  fLUTY=0;
  fLUTEta=0;
  fLUTphi0=0;
  fLUT2sinphi0=0;
  fLUT2cosphi0=0;
}

AliL3HoughTransformerLUT::AliL3HoughTransformerLUT(Int_t slice,Int_t patch,Int_t n_eta_segments) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  AliL3HoughTransformerLUT();

  Init(slice,patch,n_eta_segments);
}

AliL3HoughTransformerLUT::~AliL3HoughTransformerLUT()
{
  DeleteHistograms();

#ifdef do_mc
  if(fTrackID)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fTrackID[i]) continue;
	  delete fTrackID[i];
	}
      delete [] fTrackID;
    }
#endif

  if(fNRows){
    delete[] fLUTX;
    delete[] fLUTY;
    fNRows=0;
  }
  if(fNEtas){
    delete[] fLUTEta;
    fNEtas=0;
  }
}

void AliL3HoughTransformerLUT::DeleteHistograms()
{
  if(fNPhi0){
    delete[] fLUT2sinphi0;
    delete[] fLUT2cosphi0;
    delete[] fLUTphi0;
    fNPhi0=0;
    fLUTphi0=0;
    fLUT2sinphi0=0;
    fLUT2cosphi0=0;
  }

  if(!fParamSpace){
    for(Int_t i=0; i<GetNEtaSegments(); i++)
      {
	if(!fParamSpace[i]) continue;
	delete fParamSpace[i];
      }
    delete [] fParamSpace;
  }
}

void AliL3HoughTransformerLUT::Init(Int_t slice,Int_t patch,Int_t n_eta_segments) 
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
    fNEtas=0;
  }

  //member values
  fMinRow=AliL3Transform::GetFirstRow(patch);;
  fMaxRow=AliL3Transform::GetLastRow(patch);;
  fNRows=AliL3Transform::GetNRows(patch);;
  fNEtas=n_eta_segments;

  AliL3Transform::Slice2Sector(slice,fMinRow,fSector,fSectorRow);
  fZSign = slice < 18 ? 1:-1;
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
  for(Int_t rr=0;rr<fNEtas;rr++){
    fLUTEta[rr]=CalcRoverZ2(etamin_+(rr+1)*fEtaSlice);
    //cout << rr << ": " << fLUTEta[rr] << endl;
  }
}

void AliL3HoughTransformerLUT::CreateHistograms(Int_t nxbin,Double_t pt_min,Int_t nybin,Double_t phimin,Double_t phimax)
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

void AliL3HoughTransformerLUT::CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax)
{
  fParamSpace = new AliL3Histogram*[fNEtas];
  
  Char_t histname[256];
  for(Int_t i=0; i<fNEtas; i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
  
#ifdef do_mc
  if(fDoMC)
    {
      AliL3Histogram *hist = fParamSpace[0];
      Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
      cout<<"Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fTrackID"<<endl;
      fTrackID = new TrackIndex*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fTrackID[i] = new TrackIndex[ncells];
    }
#endif

  //create lookup table for sin and cos
  fNPhi0=nybin+1;

  fLUTphi0=new Float_t[fNPhi0];
  fLUT2sinphi0=new Float_t[fNPhi0];
  fLUT2cosphi0=new Float_t[fNPhi0];
  Float_t diff=(ymax-ymin)/nybin;
  Float_t phi0=ymin-0.5*diff;
  for(Int_t i=0; i<fNPhi0; i++){
    phi0+=diff;
    fLUTphi0[i]=phi0;
    fLUT2sinphi0[i]=Float_t(2*sin(phi0));
    fLUT2cosphi0[i]=Float_t(2*cos(phi0));
    cout << i << ": " << fLUTphi0[i] << " " << fLUT2sinphi0[i] << " " << fLUT2cosphi0[i] << endl;
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
  
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fParamSpace[i]->Reset();
#ifdef do_mc
  if(fDoMC)
    {
      AliL3Histogram *hist = fParamSpace[0];
      Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	memset(fTrackID[i],0,ncells*sizeof(TrackIndex));
    }
#endif
}


Int_t AliL3HoughTransformerLUT::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 
  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliL3Histogram *AliL3HoughTransformerLUT::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
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
  if(fZSign>0) rz=-rz;
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
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=(Double_t)((eta_index+0.5)*eta_slice);
  if(slice>17) eta*=-1;
  return eta;
}

void AliL3HoughTransformerLUT::TransformCircle()
{
#if 0
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
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircle","Data")
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
	  cerr<<"AliL3HoughTransform::TransformCircle : Mismatching padrow numbering "<<i<<" "<<(Int_t)tempPt->fRow<<endl;
	  continue;
	}

      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if((Int_t)charge <= GetLowerThreshold() || (Int_t)charge > GetUpperThreshold())
	    continue;
	  Int_t sector,row;
	  Float_t xyz[3];
	  
	  //Transform data to local cartesian coordinates:
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  
	  //Calculate the eta:
	  Double_t eta = AliL3Transform::GetEta(xyz);
	  
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = GetEtaIndex(eta);
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;
	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      printf("AliL3HoughTransformer::TransformCircle : Error getting histogram in index %d\n",eta_index);
	      continue;
	    }

	  //Do the transformation:
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi = AliL3Transform::GetPhi(xyz);
	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Float_t phi0 = hist->GetBinCenterY(b);
	      Float_t kappa = 2*sin(phi - phi0)/R;
	      hist->Fill(kappa,phi0,charge);
#ifdef do_mc
	      if(fDoMC)
		{
		  Int_t bin = hist->FindBin(kappa,phi0);
		  for(Int_t t=0; t<3; t++)
		    {
		      Int_t label = digPt[j].fTrackID[t];
		      if(label < 0) break;
		      UInt_t c;
		      for(c=0; c<MaxTrack; c++)
			if(fTrackID[eta_index][bin].fLabel[c] == label || fTrackID[eta_index][bin].fNHits[c] == 0)
			  break;
		      if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
		      fTrackID[eta_index][bin].fLabel[c] = label;
		      fTrackID[eta_index][bin].fNHits[c]++;
		    }
		}
#endif
	    }
	}
      
      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
#endif //to do later
}

Int_t AliL3HoughTransformerLUT::GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi)
{
  if(!fDoMC)
    {
      cerr<<"AliL3HoughTransformer::GetTrackID : Flag switched off"<<endl;
      return -1;
    }
  
#ifdef do_mc
  if(eta_index < 0 || eta_index > GetNEtaSegments())
    {
      cerr<<"AliL3HoughTransformer::GetTrackID : Wrong etaindex "<<eta_index<<endl;
      return -1;
    }
  AliL3Histogram *hist = fParamSpace[eta_index];
  Int_t bin = hist->FindBin(kappa,psi);
  Int_t label=-1;
  Int_t max=0;
  for(UInt_t i=0; i<MaxTrack; i++)
    {
      Int_t nhits=fTrackID[eta_index][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max)
	{
	  max = nhits;
	  label = fTrackID[eta_index][bin].fLabel[i];
	}
    }
  return label;
#endif
  cout<<"AliL3HoughTransformer::GetTrackID : Compile with do_mc flag!"<<endl;
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
  if(!fNPhi0) return;
  cout << "fLUTphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUTPhi0[" << i << "]=" << fLUTphi0[i] << endl;
  cout << "fLUT2sinphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2sinphi0[" << i << "]=" << fLUT2sinphi0[i] << endl;
  cout << "fLUT2cosphi0 " << fNPhi0 << endl;
  for(Int_t i=0;i<fNPhi0;i++) cout << "fLUT2cosphi0[" << i << "]=" << fLUT2cosphi0[i] << endl;
}
