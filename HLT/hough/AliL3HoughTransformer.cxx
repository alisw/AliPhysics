//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include <math.h>
#include <TH2.h>

#include <TFile.h>
#include <TTree.h>

#include "AliTPCParam.h"
#include "AliSimDigits.h"

#include "AliL3Histogram.h"
#include "AliL3Logging.h"
#include "AliL3Defs.h"
#include "AliL3Transform.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackArray.h"
#include "AliL3DigitData.h"

ClassImp(AliL3HoughTransformer)

AliL3HoughTransformer::AliL3HoughTransformer()
{
  //Default constructor
  
}


AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Float_t *etarange)
{
  //Constructor
  
  fTransform = new AliL3Transform();

  fEtaMin = etarange[0];
  fEtaMax = etarange[1];
  fSlice = slice; 
  fPatch = patch;
  fNumOfPadRows=NRowsSlice;
}

AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Double_t *etarange,Int_t n_eta_segments)
{
  
  fTransform = new AliL3Transform();
  if(etarange)
    {
      //assumes you want to look at a given etaslice only
      fEtaMin = etarange[0];
      fEtaMax = etarange[1];
    }
  else
    {
      //transform in all etaslices
      fEtaMin = 0;
      fEtaMax = slice < 18 ? 0.9 : -0.9;
    }
  fNumEtaSegments = n_eta_segments;
  fSlice = slice;
  fPatch = patch;
  fNumOfPadRows = NRowsSlice;
  
  fNRowsInPatch = NRows[fPatch][1]-NRows[fPatch][0] + 1;
  fBinTableBounds = (fNRowsInPatch+1)*(MaxNPads+1);
  fNDigitRowData = 0;
  fDigitRowData = 0;
  fHistoPt = 0;
}


AliL3HoughTransformer::~AliL3HoughTransformer()
{
  //Destructor
  if(fBinTable)
    {
      for(Int_t i=0; i<fBinTableBounds; i++)
	delete [] fBinTable[i];
      delete [] fBinTable;
    }
  if(fEtaIndex)
    delete [] fEtaIndex;
  if(fTrackTable)
    {
      for(Int_t i=0; i<fNumEtaSegments; i++)
	delete [] fTrackTable[i];
      delete [] fTrackTable;
    }
  if(fTransform)
    delete fTransform;
  
  
}

void AliL3HoughTransformer::SetInputData(UInt_t ndigits,AliL3DigitRowData *ptr)
{
  fNDigitRowData = ndigits;
  fDigitRowData = ptr;
}

void AliL3HoughTransformer::InitTables()
{
  //Create LUT for the circle transform.
  //the actual transformation is done in TransformTables.

  AliL3Histogram *hist = fHistoPt;
  
  Int_t nbinsy = hist->GetNbinsY();
  
  fBinTable = new Int_t*[fBinTableBounds];
  Int_t index;

  Double_t etaslice = fEtaMax/fNumEtaSegments;
  
  Int_t etabounds = (MaxNPads+1)*(fNRowsInPatch+1)*(MaxNTimeBins+1);

  fEtaIndex = new Int_t[etabounds];
  for(Int_t i=0; i<etabounds; i++)
    fEtaIndex[i] = -1;
  
  fTrackTable = new UChar_t*[fNumEtaSegments];
  for(Int_t i=0; i<fNumEtaSegments; i++)
    {
      fTrackTable[i] = new UChar_t[fBinTableBounds];
      for(Int_t j=0; j<fBinTableBounds; j++)
	fTrackTable[i][j] = 0;
    }
  
  for(Int_t r=NRows[fPatch][0]; r<=NRows[fPatch][1]; r++)
    {
      Int_t prow = r - NRows[fPatch][0];
      
      for(Int_t p=0; p<fTransform->GetNPads(r); p++)
	{
	  Float_t xyz[3];
	  Int_t sector,row;
	  for(Int_t t=0; t<fTransform->GetNTimeBins(); t++)
	    {
	      fTransform->Slice2Sector(fSlice,r,sector,row);
	      fTransform->Raw2Local(xyz,sector,row,p,t);
	      Double_t eta = fTransform->GetEta(xyz);
	      if(eta < fEtaMin || eta > fEtaMax) continue;
	      Int_t ind = (prow<<17) + (p<<9) + t;
	      if(fEtaIndex[ind]>=0)
		printf("AliL3HoughTransformer::InitTables : Overlapping indexes in eta!!\n");
	      Int_t eta_index = (Int_t)(eta/etaslice);	      
	      if(eta_index < 0 || eta_index > fNumEtaSegments)
		continue;
	      fEtaIndex[ind] = eta_index;
	      
	    }
	  
	  Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	  Double_t phi_pix = fTransform->GetPhi(xyz);
	  index = (prow<<8) + p;
	  fBinTable[index] = new Int_t[nbinsy+1];
	  
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Double_t phi0 = hist->GetBinCenterY(b);
	      Double_t kappa = 2*sin(phi_pix-phi0)/r_pix;
	      Int_t bin = hist->FindBin(kappa,phi0);
	      if(fBinTable[index][b]!=0)
		printf("AliL3HoughTransformer::InitTables : Overlapping indexes %d %d b %d\n",fBinTable[index][b],index,b);
	      fBinTable[index][b] = bin;
	    }
	}

    }

}

void AliL3HoughTransformer::TransformTables(AliL3Histogram **histos,AliL3Histogram **images)
{
  //Transform all the pixels while reading, and fill the corresponding histograms.
  //Transform is done using LUT created in InitTables.
  //fTrackTable : table telling whether a specific pixel is active (nonzero):
  //fTrackTable = 0  ->  no track
  //fTrackindex = 1  ->  track present
  //fTrackindex = 2  ->  track has been removed (already found)
  //fEtaIndex : table storing the etaindex -> used to find correct histogram to fill
  //fBinTable : table storing all the bins to fill for each nonzero pixel
    
  Int_t eta_index;
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fDigitRowData;
  AliL3Histogram *hist;
  
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformTables","data")<<
	"Zero datapointer"<<ENDLOG;
      return;
    }

  Int_t out_count=0,tot_count=0;
  Int_t index,ind;
  
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      Int_t prow = i - NRows[fPatch][0];
      AliL3DigitData *bins = tempPt->fDigitData;
      for(UInt_t dig=0; dig<tempPt->fNDigit; dig++)
	{
	  index = (prow<<8) + bins[dig].fPad;
	  ind = (prow<<17) + (bins[dig].fPad<<9) + bins[dig].fTime;
	  eta_index = fEtaIndex[ind];
	  if(eta_index < 0) continue;  //pixel out of etarange
	  
	  if(fTrackTable[eta_index][index]==2) continue; //this pixel has already been removed. 
	  fTrackTable[eta_index][index] = 1; //this pixel is now marked as active (it is nonzero)
	  
	  tot_count++;
	  hist = histos[eta_index];
	  
	  if(images)
	    {
	      //Display the transformed images.
	      AliL3Histogram *image = images[eta_index];
	      Float_t xyz_local[3];
	      Int_t sector,row;
	      fTransform->Slice2Sector(fSlice,i,sector,row);
	      fTransform->Raw2Local(xyz_local,sector,row,bins[dig].fPad,bins[dig].fTime);
	      image->Fill(xyz_local[0],xyz_local[1],bins[dig].fCharge);
	    }
	  
	  if(!hist)
	    {
	      printf("Error getting histogram!\n");
	      continue;
	    }
	  for(Int_t p=hist->GetFirstYbin(); p<=hist->GetLastYbin(); p++)
	    hist->AddBinContent(fBinTable[index][p],bins[dig].fCharge);
	  
	}
      
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData) + tempPt->fNDigit*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
    }  
  
}

void AliL3HoughTransformer::WriteTables()
{
  //Write the tables to asci files.
  
  AliL3Histogram *hist = fHistoPt;
  Char_t name[100];
  sprintf(name,"histogram_table_%d.txt",fPatch);
  FILE *file = fopen(name,"w");
  
  Int_t etabounds = (MaxNPads+1)*(fNRowsInPatch+1)*(MaxNTimeBins+1);
  for(Int_t i=0; i<etabounds; i++)
    {
      if(fEtaIndex[i]<0) continue;
      fprintf(file,"%d %d\n",i,fEtaIndex[i]);
    }
  fclose(file);
  
  sprintf(name,"bin_table_%d.txt",fPatch);
  FILE *file2 = fopen(name,"w");
  for(Int_t i=0; i<fBinTableBounds; i++)
    {
      if(!fBinTable[i]) continue;
      fprintf(file2,"%d ",i);
      for(Int_t j=hist->GetFirstYbin(); j<=hist->GetLastYbin(); j++)
	{
	  fprintf(file2,"%d ",fBinTable[i][j]);
	}
      fprintf(file2,"\n");
    }
  fclose(file2);
}

/*
void AliL3HoughTransformer::InitTemplates(TH2F *hist)
{

  AliL3Digits *pixel;

  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  Int_t nbinsy = hist->GetNbinsY();

  fIndex = new Int_t*[fNDigits];
  for(Int_t i=0; i<fNDigits; i++)
    fIndex[i] = new Int_t[nbinsy+1];
    
  Int_t sector,row;
  
  for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
    {        
      
      for(pixel=(AliL3Digits*)fRowContainer[padrow].first; pixel!=0; pixel=(AliL3Digits*)pixel->nextRowPixel)
	{
	  Float_t xyz[3];
	  fTransform->Slice2Sector(fSlice,padrow,sector,row);
	  fTransform->Raw2Local(xyz,sector,row,pixel->fPad,pixel->fTime);
	  
	  Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	  
	  Double_t phi_pix = fTransform->GetPhi(xyz);
	  Int_t index = pixel->fIndex;
	  if(index >= fNDigits)
	    printf("AliL3HoughTransformer::InitTemplates : Index error! %d\n",index);
	  for(Int_t p=ymin; p<=ymax; p++)
	    {
	      
	      Double_t phi0 = hist->GetYaxis()->GetBinCenter(p);
	      Double_t kappa = 2*sin(phi_pix-phi0)/r_pix;
	      //printf("kappa %f phi0 %f\n",kappa,phi0);
	      Int_t bin = hist->FindBin(kappa,phi0);
	      if(fIndex[index][p]!=0)
		printf("AliL3HoughTransformer::InitTemplates : Overlapping indexes\n");
	      fIndex[index][p] = bin;
	    }
	}
    }
  
}


void AliL3HoughTransformer::Transform2Circle(TH2F *hist,Int_t eta_index)
{
  //Transformation is done with respect to local coordinates in slice.
  //Transform every pixel into whole phirange, using parametrisation:
  //kappa = 2*sin(phi-phi0)/R
  //Assumes you run InitTemplates first!!!!

  AliL3Digits *pix1;
    
  Int_t nbinsy = hist->GetNbinsY();

  Int_t totsignal=0,vol_index;
  if(fNumEtaSegments==1)
    eta_index=0; //only looking in one etaslice.

  for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
    {
      vol_index = eta_index*fNumOfPadRows + padrow;
      //for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
      for(pix1=(AliL3Digits*)fVolume[vol_index].first; pix1!=0; pix1=(AliL3Digits*)pix1->fNextVolumePixel)
	{
	  Float_t xyz[3];
	  fTransform->Raw2Global(xyz,2,padrow,pix1->fPad,pix1->fTime);
	  Double_t eta = fTransform->GetEta(xyz);
	  if(eta < fEtaMin || eta > fEtaMax)
	    printf("\n Eta OUT OF RANGE\n");

	  Short_t signal = pix1->fCharge;
	  Int_t index = pix1->fIndex;
	  if(index < 0) continue; //This pixel has been removed.
	  totsignal += signal;
	  for(Int_t p=0; p <= nbinsy; p++)
	    hist->AddBinContent(fIndex[index][p],signal);
	  		  
	}
      
    }
    
  printf("Total signal %d\n",totsignal);
}

void AliL3HoughTransformer::Transform2Circle(TH2F **histos,Int_t n_eta_segments,UInt_t ndigits,AliL3DigitRowData *ptr)
{
  //Transform all the pixels while reading them, and fill the corresponding histograms.
  //Everything is done in one go here. 

  Double_t etaslice = 0.9/n_eta_segments;
  Int_t eta_index;
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)ptr;
  TH2F *hist;

  Int_t out_count=0,tot_count=0;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      printf("doing row %d\n",i);
      for(UInt_t dig=0; dig<tempPt->fNDigit; dig++)
	{
	  
	  AliL3DigitData *bins = tempPt->fDigitData;
	  Float_t xyz[3];
	  Int_t sector,row;
	  fTransform->Slice2Sector(fSlice,i,sector,row);
	  fTransform->Raw2Local(xyz,sector,row,bins[dig].fPad,bins[dig].fTime);
	  Double_t eta = fTransform->GetEta(xyz);
	  eta_index = (Int_t)(eta/etaslice);
	  if(eta_index < 0 || eta_index >= n_eta_segments)
	    {
	      //printf("Eta index out of range %d\n",eta_index);
	      out_count++;
	      continue;
	    }
	  tot_count++;
	  hist = histos[eta_index];
	  if(!hist)
	    {
	      printf("Error getting histogramm!\n");
	      return;
	    }
	  
	  Int_t ymin = hist->GetYaxis()->GetFirst();
	  Int_t ymax = hist->GetYaxis()->GetLast();
	  
	  Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	  Double_t phi_pix = fTransform->GetPhi(xyz);
	  for(Int_t p=ymin; p<=ymax; p++)
	    {
	      Double_t phi0 = hist->GetYaxis()->GetBinCenter(p);
	      Double_t kappa = 2*sin(phi_pix-phi0)/r_pix;
	      //printf("kappa %f phi0 %f\n",kappa,phi0);
	      hist->Fill(kappa,phi0,bins[dig].fCharge);
	    }
	}
      
      
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData) + tempPt->fNDigit*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
      
    }
  
  printf("There were %d pixels out of range and %d inside\n", out_count,tot_count);

}


void AliL3HoughTransformer::TransformLines2Circle(TH2F *hist,AliL3TrackArray *tracks)
{

  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) {printf("AliL3HoughTransformer::TransformLines2Circle : NO TRACK IN ARRAY\n"); continue;}
     
      Double_t xmin = fTransform->Row2X(track->GetFirstRow());
      Double_t xmax = fTransform->Row2X(track->GetLastRow());
      
      Double_t a = -1./tan(track->GetPsiLine());
      Double_t b = track->GetDLine()/sin(track->GetPsiLine());
      
      Double_t ymin = a*xmin + b;
      Double_t ymax = a*xmax + b;
      
      Double_t middle_x = xmin + (xmax-xmin)/2;
      Double_t middle_y = ymin + (ymax-ymin)/2;
      
      Double_t r_middle = sqrt(middle_x*middle_x + middle_y*middle_y);
      Double_t phi = atan2(middle_y,middle_x);
      Double_t phi0 = 2*phi - track->GetPsiLine();
      Double_t kappa = 2*sin(phi-phi0)/r_middle;
      hist->Fill(kappa,phi0,track->GetWeight());
     
    }

  
}

void AliL3HoughTransformer::Transform2Line(TH2F *hist,Int_t ref_row,Int_t *rowrange,Double_t *phirange,TH2F *raw)
{
  //Transform every pixel into histogram, using parametrisation:
  //D = x*cos(psi) + y*sin(psi)

  //  printf("In Transform; rowrange %d %d ref_row %d phirange %f %f\n",rowrange[0],rowrange[1],ref_row,phirange[0],phirange[1]);

  AliL3Digits *pix1;
  
  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  
  Double_t x0 = fTransform->Row2X(ref_row);
  Double_t y0 = 0;

  Int_t sector,row;
  //for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
  
  Double_t phi_min = -10*ToRad;
  Double_t phi_max = 10*ToRad;
  Double_t delta_phi = (phi_max-phi_min)/fNPhiSegments;
  

  Int_t phi_min_index = (Int_t)((phirange[0]*ToRad-phi_min)/delta_phi);
  Int_t phi_max_index = (Int_t)((phirange[1]*ToRad-phi_min)/delta_phi);
    
  
  Int_t index;
  
  for(Int_t phi=phi_min_index; phi <= phi_max_index; phi++)
    {
      for(Int_t padrow = rowrange[0]; padrow <= rowrange[1]; padrow++)
	{
	  
	  index = phi*fNumOfPadRows + padrow;
	  //printf("Looping index %d\n",index);
	  if(index > fContainerBounds || index < 0)
	    {
	      printf("AliL3HoughTransformer::Transform2Line : index %d out of range \n",index);
	      return;
	    }
	  for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
	  //for(pix1=(AliL3Digits*)fPhiRowContainer[index].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextPhiRowPixel)
	    {
	      //printf("Transforming pixel in index %d pad %d time %d padrow %d\n",index,pix1->fPad,pix1->fTime,padrow);
	      Float_t xyz[3];
	      fTransform->Slice2Sector(fSlice,padrow,sector,row);
	      fTransform->Raw2Local(xyz,sector,row,pix1->fPad,pix1->fTime);
	      	  
    	      if(raw)
		raw->Fill(xyz[0],xyz[1],pix1->fCharge);

	      xyz[0] = xyz[0]-x0;
	      xyz[1] = xyz[1]-y0;
	      
	      //printf("calculating...");
	      for(Int_t d=xmin+1; d<=xmax; d++)
		{
		  Double_t psi = hist->GetXaxis()->GetBinCenter(d);
		  Double_t D = xyz[0]*cos(psi) + xyz[1]*sin(psi);
		  
		  Short_t signal = pix1->fCharge;
		  hist->Fill(psi,D,signal);
		  //printf("Filling histo, psi %f D %f\n",psi,D);
		}
	      //printf("done\n");
	    }
	  //printf(" done\n");
	}
      
    }
    
}


void AliL3HoughTransformer::GetPixels(Char_t *rootfile,TH2F *hist)
{

  //read data from rootfile. more or less obsolete code this.

  TFile *file = new TFile(rootfile);
  file->cd();

  AliTPCParam *param=(AliTPCParam *)file->Get("75x40_100x60");
  TTree *t=(TTree*)file->Get("TreeD_75x40_100x60");      
    
  AliSimDigits da, *digarr=&da;
  t->GetBranch("Segment")->SetAddress(&digarr);
  Stat_t num_of_entries=t->GetEntries();

  Int_t digit_counter=0;
  Float_t xyz[3];
  Double_t eta;
  
  Int_t nrows = NRows[fPatch][1] - NRows[fPatch][0] + 1;
  printf("nrows %d slice %d patch %d\n",nrows,fSlice,fPatch);
  
  if(fRowContainer)
    delete [] fRowContainer;
  fRowContainer = new AliL3HoughContainer[fNumOfPadRows];
  memset(fRowContainer,0,fNumOfPadRows*sizeof(AliL3HoughContainer));


  //fContainerBounds = (fNPhiSegments+1)*(fNumOfPadRows+1);
  //printf("Allocating %d bytes to container of size %d\n",fContainerBounds*sizeof(AliL3HoughContainer),fContainerBounds);

  
  //  if(fPhiRowContainer)
  //  delete [] fPhiRowContainer;
  //  fPhiRowContainer = new AliL3HoughContainer[fContainerBounds];
  //  memset(fPhiRowContainer,0,fContainerBounds*sizeof(AliL3HoughContainer));
  
  fContainerBounds = (fNumEtaSegments+1)*(fNumOfPadRows+1);
  if(fVolume)
    delete [] fVolume;
  fVolume = new AliL3HoughContainer[fContainerBounds];
  memset(fVolume,0,fContainerBounds*sizeof(AliL3HoughContainer));
  
  //  Double_t phi_min = -10*ToRad;
  //  Double_t phi_max = 10*ToRad;
  //  Double_t delta_phi = (phi_max-phi_min)/fNPhiSegments;
  
  Double_t eta_slice = (fEtaMax-fEtaMin)/fNumEtaSegments;
  Int_t index;
  digit_counter=0;

  //printf("\nLoading ALL pixels in slice\n\n");

  for (Int_t i=0; i<num_of_entries; i++) 
    { 
      t->GetEvent(i);
      Int_t sector; 
      Int_t row;    
      param->AdjustSectorRow(digarr->GetID(),sector,row);
      Int_t slice,padrow;
      fTransform->Sector2Slice(slice,padrow,sector,row);
      if(slice != fSlice) continue;
      if(padrow < NRows[fPatch][0]) continue;
      if(padrow > NRows[fPatch][1]) break;
      digarr->First();
      do {
	Int_t time=digarr->CurrentRow();
	Int_t pad=digarr->CurrentColumn();
	Short_t signal=digarr->CurrentDigit();
	if(time < param->GetMaxTBin()-1 && time > 0)
	  if(digarr->GetDigit(time-1,pad) <= param->GetZeroSup()
	     && digarr->GetDigit(time+1,pad) <= param->GetZeroSup()) 
	    continue;
	
	
	fTransform->Raw2Global(xyz,sector,row,pad,time);
	eta = fTransform->GetEta(xyz);
	
	if(eta < fEtaMin || eta > fEtaMax) continue;
	fTransform->Global2Local(xyz,sector);
	
	
	//phi = fTransform->GetPhi(xyz);
	if(hist)
	  hist->Fill(xyz[0],xyz[1],signal);
	
	AliL3Digits *dig = new AliL3Digits;
	dig->fIndex = digit_counter;
	digit_counter++;
	dig->fCharge = signal;
	dig->fPad = pad;
	dig->fTime = time;
	
	if(fRowContainer[padrow].first == NULL)
	  fRowContainer[padrow].first = (void*)dig;
	else
	  ((AliL3Digits*)(fRowContainer[padrow].last))->nextRowPixel=dig;
	fRowContainer[padrow].last = (void*)dig;
	
	//thisHit->etaIndex=(Int_t)((thisHit->GetEta()-fEtaMin)/etaSlice + 1);
	Int_t eta_index = (Int_t)((eta-fEtaMin)/eta_slice);
	index = eta_index*fNumOfPadRows + padrow;
	if(index > fContainerBounds || index < 0)
	  {
	    
	    //printf("AliL3HoughTransformer::GetPixels : index out of range %d %d eta_index %d padrow %d\n",index,fContainerBounds,eta_index,padrow);
	    continue;
	  }
	
	if(fVolume[index].first == NULL)
	  fVolume[index].first = (void*)dig;
	else
	  ((AliL3Digits*)(fVolume[index].last))->fNextVolumePixel = dig;
	fVolume[index].last = (void*)dig;
	
	
	//  Int_t phi_index = (Int_t)((phi-phi_min)/delta_phi);
	//  index = phi_index*fNumOfPadRows + padrow;
	//  if(phi_index > fContainerBounds || phi_index < 0)
	// {
	//  printf("AliL3HoughTransform::GetPixels : index out of range %d\n",phi_index);
	//  continue;
	//  }
	  
	//  if(fPhiRowContainer[index].first == NULL)
	//  fPhiRowContainer[index].first = (void*)dig;
	//  else
	//  ((AliL3Digits*)(fPhiRowContainer[index].last))->nextPhiRowPixel = dig;
	//  fPhiRowContainer[index].last=(void*)dig;
	
      }while (digarr->Next());
      
    }
  
  fNDigits = digit_counter;
  printf("digitcounter %d\n",digit_counter);
  printf("Allocated %d bytes to pixels\n",digit_counter*sizeof(AliL3Digits));
  file->Close();
  delete file;
  
}
*/
