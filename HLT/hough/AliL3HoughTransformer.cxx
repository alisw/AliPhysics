#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TH2.h>

#include "AliTPCParam.h"
#include "AliSimDigits.h"

#include "AliL3Defs.h"
#include "AliL3Transform.h"
#include "AliL3HoughPixel.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackArray.h"

ClassImp(AliL3HoughTransformer)

AliL3HoughTransformer::AliL3HoughTransformer()
{
  //Default constructor
  fTransform = 0;
  fEtaMin = 0;
  fEtaMax = 0;
  fSlice = 0;
  fPatch = 0;
  fRowContainer = 0;
  fPhiRowContainer = 0;
  fNumOfPadRows=0;
  fContainerBounds=0;
  fNDigits=0;
  fIndex = 0;
}


AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Float_t *etarange,Int_t phi_segments)
{
  //Constructor
  
  fTransform = new AliL3Transform();
  

  fEtaMin = etarange[0];
  fEtaMax = etarange[1];
  fSlice = slice; 
  fPatch = patch;
  fNPhiSegments = phi_segments;
  fNumOfPadRows=NRowsSlice;
}


AliL3HoughTransformer::~AliL3HoughTransformer()
{
  //Destructor
  if(fRowContainer)
    delete [] fRowContainer;
  if(fTransform)
    delete fTransform;
  if(fPhiRowContainer)
    delete [] fPhiRowContainer;
  if(fIndex)
    {
      for(Int_t i=0; i<fNDigits; i++)
	delete [] fIndex[i];
      delete [] fIndex;
    }
}

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
	  Short_t signal = pixel->fCharge;
	  Int_t index = pixel->fIndex;
	  if(index >= fNDigits)
	    printf("AliL3HoughTransformer::InitTemplates : Index error! %d\n",index);
	  for(Int_t p=ymin; p<=ymax; p++)
	    {
	      
	      Double_t phi0 = hist->GetYaxis()->GetBinCenter(p);
	      Double_t kappa = 2*sin(phi_pix-phi0)/r_pix;
	      
	      Int_t bin = hist->FindBin(kappa,phi0);
	      if(fIndex[index][p]!=0)
		printf("AliL3HoughTransformer::InitTemplates : Overlapping indexes\n");
	      fIndex[index][p] = bin;
	    }
	}
    }
  
}


void AliL3HoughTransformer::CountBins()
{
  
  Int_t middle_row = 87; //middle of the slice
  
  Double_t r_in_bundle = fTransform->Row2X(middle_row);
  //  Double_t phi_min = (fSlice*20 - 10)*ToRad;
  //Double_t phi_max = (fSlice*20 + 10)*ToRad;
  Double_t phi_min = -15*ToRad;
  Double_t phi_max = 15*ToRad;

  Double_t phi_slice = (phi_max - phi_min)/fNPhiSegments;
  Double_t min_phi0 = 10000;
  Double_t max_phi0 = 0;
  Double_t min_kappa = 100000;
  Double_t max_kappa = 0;
  
  Int_t xbin = 60;
  Int_t ybin = 60;
  Float_t xrange[2] = {-0.006 , 0.006}; //Pt 0.2->
  Float_t yrange[2] = {-0.26 , 0.26}; //slice 2 0.55->0.88
  
  TH2F *histo = new TH2F("histo","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  
  for(Int_t padrow=NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
    {
      for(Int_t pad=0; pad < fTransform->GetNPads(padrow); pad++)
	{
	  for(Int_t time = 0; time < fTransform->GetNTimeBins(); time++)
	    {
	      Float_t xyz[3];
	      Int_t sector,row;
	      fTransform->Slice2Sector(fSlice,padrow,sector,row);
	      fTransform->Raw2Global(xyz,sector,row,pad,time);
	      Double_t eta = fTransform->GetEta(xyz);
	      if(eta < fEtaMin || eta > fEtaMax) continue;
	      fTransform->Global2Local(xyz,sector);
	      Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	      Double_t phi_pix = fTransform->GetPhi(xyz);
	      
	      for(Int_t p=0; p<=fNPhiSegments; p++)
		{
		  Double_t phi_in_bundle = phi_min + p*phi_slice;
		  
		  Double_t tanPhi0 = (r_pix*sin(phi_in_bundle)-r_in_bundle*sin(phi_pix))/(r_pix*cos(phi_in_bundle)-r_in_bundle*cos(phi_pix));
		  
		  Double_t phi0 = atan(tanPhi0);
		  //  if(phi0 < 0.55 || phi0 > 0.88) continue;
		  
		  //if(phi0 < 0) phi0 = phi0 +2*Pi;
		  //Double_t kappa = sin(phi_in_bundle - phi0)*2/r_in_bundle;
		  
		  Double_t angle = phi_pix - phi0;
		  Double_t kappa = 2*sin(angle)/r_pix;
		  histo->Fill(kappa,phi0,1);
		  if(phi0 < min_phi0)
		    min_phi0 = phi0;
		  if(phi0 > max_phi0)
		    max_phi0 = phi0;
		  if(kappa < min_kappa)
		    min_kappa = kappa;
		  if(kappa > max_kappa)
		    max_kappa = kappa;
		    		  
		}
	      
	    }
	  
	}
      
    }
  Int_t count=0,bi=0;
    
  Int_t xmin = histo->GetXaxis()->GetFirst();
  Int_t xmax = histo->GetXaxis()->GetLast();
  Int_t ymin = histo->GetYaxis()->GetFirst();
  Int_t ymax = histo->GetYaxis()->GetLast();


  for(Int_t xbin=xmin+1; xbin<xmax; xbin++)
    {
      for(Int_t ybin=ymin+1; ybin<ymax; ybin++)
	{
	  bi++;
	  Int_t bin = histo->GetBin(xbin,ybin);
	  if(histo->GetBinContent(bin)>0)
	    count++;
	}
    }


  printf("Number of possible tracks in this region %d, bins %d\n",count,bi);
  printf("Phi, min %f max %f\n",min_phi0,max_phi0);
  printf("Kappa, min %f max %f\n",min_kappa,max_kappa);
  histo->Draw("box");
}


void AliL3HoughTransformer::Transform2Circle(TH2F *hist,Int_t middle_row)
{
  //Transformation is done with respect to local coordinates in slice.


  AliL3Digits *pix1;
  Int_t sector,row;
  
  //Define a common point
 
  /*
    Double_t rowdist1 = fTransform->Row2X(middle_row-1);
    Double_t rowdist2 = fTransform->Row2X(middle_row);
    Double_t r_in_bundle = rowdist1 + (rowdist1-rowdist2)/2;
  */
  
  Double_t r_in_bundle = fTransform->Row2X(middle_row);
  //Make overlap between slices
  Double_t phi_min = -15*ToRad;
  Double_t phi_max = 15*ToRad;
  
  Double_t phi_slice = (phi_max - phi_min)/fNPhiSegments;
  
  for(Int_t p=0; p <= fNPhiSegments; p++)
    {
      Double_t phi_in_bundle = phi_min + p*phi_slice;
      //printf("phi %f in slice %d patch %d middle row %f\n",phi_in_bundle/ToRad,fSlice,fPatch,r_in_bundle);
      
      for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
	//for(Int_t padrow = middle_row; padrow <= 173; padrow++)
	//for(Int_t padrow = 0; padrow <= middle_row; padrow++)
	{
	  
	  for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
	    {
	      
	      Float_t xyz[3];
	      fTransform->Slice2Sector(fSlice,padrow,sector,row);
	      fTransform->Raw2Local(xyz,sector,row,pix1->fPad,pix1->fTime);
	      //fTransform->Raw2Global(xyz,sector,row,pix1->fPad,pix1->fTime);
	      
	      Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	      
	      Double_t phi_pix = fTransform->GetPhi(xyz);
	        	      
	      //Double_t tanPhi0 = (r_pix*sin(phi_in_bundle)-r_in_bundle*sin(phi_pix))/(r_pix*cos(phi_in_bundle)-r_in_bundle*cos(phi_pix));
	      Double_t tanPhi0 = (r_in_bundle*sin(phi_pix)-r_pix*sin(phi_in_bundle))/(r_in_bundle*cos(phi_pix)-r_pix*cos(phi_in_bundle));
	      
	      Double_t phi0 = atan(tanPhi0);
	      //if(padrow > middle_row)
	      //phi0 = -phi0;
	      //if(phi0 < 0.55 || phi0 > 0.88) continue;
	      
	      Double_t angle = phi_pix - phi0;
	      Double_t kappa = 2*sin(angle)/r_pix;
	      
	      //Double_t angle = phi_in_bundle - phi0;
	      //Double_t kappa = 2*sin(angle)/r_in_bundle;

	      //if(kappa < -0.006 || kappa > 0.006) continue;
	      
	      Short_t signal = pix1->fCharge;
	      
	      hist->Fill(kappa,phi0,signal);
	      
	      
	    }
	  
	}
    }
  
}


void AliL3HoughTransformer::Transform2Circle(TH2F *hist)
{
  //Transformation is done with respect to local coordinates in slice.
  //Transform every pixel into whole phirange, using parametrisation:
  //kappa = 2*sin(phi-phi0)/R

  printf("Transforming 1 pixel only\n");

  AliL3Digits *pix1;
  Int_t sector,row;
  
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  Int_t nbinsy = hist->GetNbinsY();

  for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
    {
      
      for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
	{
	  
	  Short_t signal = pix1->fCharge;
	  Int_t index = pix1->fIndex;
	  
	  for(Int_t p=0; p <= nbinsy; p++)
	    hist->AddBinContent(fIndex[index][p],signal);
	  		  
	}
      
    }
    
  
}

/*
  void AliL3HoughTransformer::Transform2Circle(TH2F *hist)
  {
  //Transformation is done with respect to local coordinates in slice.
  //Transform every pixel into whole phirange, using parametrisation:
  //kappa = 2*sin(phi-phi0)/R
  
  printf("Transforming 1 pixel only\n");
  
  AliL3Digits *pix1;
  Int_t sector,row;
  
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  
  for(Int_t padrow = NRows[fPatch][0]; padrow <= NRows[fPatch][1]; padrow++)
  {
  
  for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
  {
  
  Float_t xyz[3];
  fTransform->Slice2Sector(fSlice,padrow,sector,row);
  fTransform->Raw2Local(xyz,sector,row,pix1->fPad,pix1->fTime);
  
  Double_t r_pix = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  
  Double_t phi_pix = fTransform->GetPhi(xyz);
  Short_t signal = pix1->fCharge;
  
  for(Int_t p=ymin+1; p<=ymax; p++)
  {
  Double_t phi0 = hist->GetYaxis()->GetBinCenter(p);
  Double_t kappa = 2*sin(phi_pix-phi0)/r_pix;
  hist->Fill(kappa,phi0,signal);
  }
  
  }
  
  }
  
  
  }
*/

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
	  //for(pix1=(AliL3Digits*)fRowContainer[padrow].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextRowPixel)
	  for(pix1=(AliL3Digits*)fPhiRowContainer[index].first; pix1!=0; pix1=(AliL3Digits*)pix1->nextPhiRowPixel)
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

  TFile *file = new TFile(rootfile);
  file->cd();

  AliTPCParam *param=(AliTPCParam *)file->Get("75x40_100x60");
  TTree *t=(TTree*)file->Get("TreeD_75x40_100x60");      
    
  AliSimDigits da, *digarr=&da;
  t->GetBranch("Segment")->SetAddress(&digarr);
  Stat_t num_of_entries=t->GetEntries();

  Int_t digit_counter=0;
  Float_t xyz[3];
  Double_t eta,phi;
  
  Int_t nrows = NRows[fPatch][1] - NRows[fPatch][0] + 1;
  printf("nrows %d slice %d patch %d\n",nrows,fSlice,fPatch);
  
  if(fRowContainer)
    delete [] fRowContainer;
  fRowContainer = new AliL3HoughContainer[fNumOfPadRows];
  memset(fRowContainer,0,fNumOfPadRows*sizeof(AliL3HoughContainer));


  fContainerBounds = (fNPhiSegments+1)*(fNumOfPadRows+1);
  printf("Allocating %d bytes to container of size %d\n",fContainerBounds*sizeof(AliL3HoughContainer),fContainerBounds);
  if(fPhiRowContainer)
    delete [] fPhiRowContainer;
  fPhiRowContainer = new AliL3HoughContainer[fContainerBounds];
  memset(fPhiRowContainer,0,fContainerBounds*sizeof(AliL3HoughContainer));

  Double_t phi_min = -10*ToRad;
  Double_t phi_max = 10*ToRad;
  Double_t delta_phi = (phi_max-phi_min)/fNPhiSegments;
  Int_t index;
  digit_counter=0;

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
	
	
	phi = fTransform->GetPhi(xyz);
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
	
	Int_t phi_index = (Int_t)((phi-phi_min)/delta_phi);
	index = phi_index*fNumOfPadRows + padrow;
	if(phi_index > fContainerBounds || phi_index < 0)
	  {
	    printf("AliL3HoughTransform::GetPixels : index out of range %d\n",phi_index);
	    continue;
	  }
		
	if(fPhiRowContainer[index].first == NULL)
	  fPhiRowContainer[index].first = (void*)dig;
	else
	  ((AliL3Digits*)(fPhiRowContainer[index].last))->nextPhiRowPixel = dig;
	fPhiRowContainer[index].last=(void*)dig;
	
      }while (digarr->Next());
      
    }
  
  fNDigits = digit_counter;
  printf("digitcounter %d\n",digit_counter);
  printf("Allocated %d bytes to pixels\n",digit_counter*sizeof(AliL3Digits));
  file->Close();
  delete file;
  
}
