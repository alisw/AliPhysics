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

/*
$Log$
Revision 1.1.4.2  2000/04/10 11:37:42  kowal2

Digits handling in a new data structure

*/

//-----------------------------------------------------------------------------
//
//
//  Author:   Marian Ivanov
//
//  Implementation of class TTPCDigitsH
//
//-----------------------------------------------------------------------------
////////////////////////////////////////////////
//  Manager class for AliTPCDigitsH             //
////////////////////////////////////////////////
 
#include <iostream.h>
#include "TMath.h"
// other include files follow here
#include "TFile.h"
//*KEEP,TFile.
#include "TROOT.h"
#include "TGraph.h"
#include "AliRun.h"
#include "AliDisplay.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TCanvasImp.h"
#include "TPaveText.h"
//*KEEP,TH1.
#include "TH1.h"
//*KEEP,TH2
#include "AliH2F.h"
//*KEEP,TF2.
#include "TF2.h"
//*KEEP,TClonesArray,T=C++.
#include "TClonesArray.h"
#include "TTree.h"
//GALICE includes
#include "GParticle.h"
#include "AliTPC.h"
#include "AliTPCParam.h"
#include "AliTPCD.h"
#include "AliTPCDigitsH.h"


R__EXTERN TSystem *  gSystem;
R__EXTERN AliRun *   gAlice;



ClassImp(AliTPCDigitsH)

AliTPCDigitsH::AliTPCDigitsH() 
{
//Begin_Html
/*
<img src="gif/TPCDigitsH.gif">
*/
//End_Html
 
  fsec = 1;
  frow = 1;
  fpad = 1;
  fTimeN = 500;
  fTimeStart = 0;
  fTimeStop = 500;  
  fOccuN    =25;
  fbDelHisto = kFALSE;
  fEventN = 0;  //event nuber connected to digit tree
  fThreshold  = 5;
  
  fParticles = 0;
  fDigits= 0; 

  fDParam  =0;
  
  fbIOState  = kFALSE;
  fbDigState = kFALSE;
}

AliTPCDigitsH::~AliTPCDigitsH() 
{

}

AliTPCDigitsH::AliTPCDigitsH(const AliTPCDigitsH &) 
{
}

AliTPCDigitsH & AliTPCDigitsH::operator = (const AliTPCDigitsH &) 
{
   return *this;
}
void  AliTPCDigitsH::SetDParam(AliTPCD * dig)
{
  if (dig!=0) 
    {
      fDParam =dig;
      fTPCParam = &(fDParam->GetParam());
    }
  else
    {
      fTPCParam = 0;
      fDParam =0;
    }
} 

AliTPCParam *&  AliTPCDigitsH::GetParam()
{
  return fTPCParam;
}


void AliTPCDigitsH::CloseFiles()
{
 if (fin) 
    {
      fin->Close();
      //try 
      // {
      delete fin; 
      //  //	}
      //  //catch(...)
      //  //	{
      //  cout<<"FIN Delete error. Contact autor of root \n";
      //	};      
      fin = 0;
    };
 if (fout) 
    {
      fout->Close();  
      //try 
      //{
      //  delete fout;
      //}
      //      catch(...)
      //{
      //  cout<<"Out Delete error. Contact autor of root \n";
      //};      
      fout =0;
    };
 fbIOState  = kFALSE;
 fbDigState = kFALSE;
}

Bool_t AliTPCDigitsH::SetIO(const char *  inFile, const char* outFile )
{
  ///   Set input and output file 
  ///   control permisions and so on
  ///   if all is OK then set flag fbIOState = kTRUE
 fbIOState = kFALSE; 
 TString s1=inFile;
 TString s2=outFile; 
 //  important ---- it close previious open file if this file exist
 if (fin) 
   {
     if (fin == fout) fout = NULL;
     fin->Close();
     //     try
     //{
     delete fin;
     //}
     //catch(...)
     //{
     //  cout<<"Fin  Delete error. Contact autor of root \n";
     //};
     fin = 0;
   };   
 //  important ---- it close previous open file if this file exist
 if (fout) 
   {
     fout->Close(); 
     //try
     //{
     delete fout;
     //}
     //catch(...)
     //{
     //  cout<<"Fout  Delete error. Contact autor of root \n";
     //};
     fout = 0;
   };

 //close the files if exist in root enviroment
 TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject((char*)inFile);
 if (file) file->Close();
 file = (TFile*)gROOT->GetListOfFiles()->FindObject((char*)outFile);
 if (file) file->Close();
 //if input file is the output file 
  if (s1 == s2 )
    {
      if (gSystem->AccessPathName((char*)inFile, kWritePermission ) == kFALSE)
         fin  = new TFile((char*)inFile,"UPDATE","trees with digits");                
      if (!(fin)) 
       {
	 cout<<"Input file couldn't be open for writing \n";
	 cout<<"Loook if the file exiest or if you have permision to write \n";        
	 return  kFALSE;
       }
      else
	{ 
	  fout = fin;
          fbIOState = kTRUE;
          return kTRUE;
	} 
    }
  //open input file 
  if (gSystem->AccessPathName((char*)inFile, kReadPermission ) == kFALSE)      
    fin = new TFile((char*)inFile,"UPDATE","trees with digits");
  else
    {
      cout<<"Input file couldn't be open\n";
      cout<<"Maybe the file is not open \n";        
      return kFALSE ;
    }          
  if (!(fin)) 
    {
      cout<<"Input file couldn't be open\n";
      cout<<"Probably not root file \n";     
      return kFALSE ;
    }

  //open output file  
  if (gSystem->AccessPathName((char*)outFile, kWritePermission) == kFALSE)    
    fout = new TFile((char *)outFile,"UPDATE");   
  else  
    fout = new TFile((char *)outFile,"NEW");     
   
  if (!(fout)) 
    {
      cout<<"Output  file couldn't be open\n";
      return kFALSE ;
    }
  //if input and output file is OK set state variable to true
  fbIOState = kTRUE;  
  return kTRUE;
}

Bool_t AliTPCDigitsH::SetEventN(Int_t EventN=0)
{
  if (!(fin))
    {
      cout<<"Warning: Input file not open !!! \n";   
      return kFALSE;
    }
  fin->cd();
  if (EventN>-1) fEventN = EventN;
  fParticles = 0;
  if (gAlice) 
    {
      delete gAlice;
      gAlice =0;
    }
  gAlice = (AliRun*) fin->Get("gAlice");   
  if (gAlice == 0 )
    cout<<" Warning : AliRun objects not found in input file \n.";
  else
    {
      gAlice->GetEvent(fEventN);
      fParticles = gAlice->Particles();
    }
  if (fParticles == 0) 
    return kFALSE;
  else
    return kTRUE;
}



Bool_t AliTPCDigitsH::SetTree(Int_t eventn  )
{
   fbDigState = kFALSE;
   ftree = 0;
   if (  fbIOState == kFALSE)
     {
       cout<<"IO files not adjusted \n FIX IT!!!";
       return kFALSE;
     }
   fin->cd();
   if (fDParam->SetTree(eventn)==0){
     fbDigState = kFALSE;
     cout<<"Input tree doesn't exist !!! \n";
     return kFALSE;
   }
   fDigits = fDParam->GetArray();
   ftree = fDParam->GetTree();
   fbDigState=kTRUE;
   return kTRUE; 
}

void AliTPCDigitsH::SetSecRowTime(Int_t sec , Int_t row , Int_t TimeN, Float_t TimeStart, Float_t TimeStop )
{
  fsec = sec;
  frow = row;
  fTimeN =TimeN;
  fTimeStart = TimeStart;
  fTimeStop = TimeStop;
}

void   AliTPCDigitsH::SetParticles(Int_t sec = -1, Int_t row = -1 ,
		       Int_t size1 = 30000,Int_t size2=300,
		       Bool_t all=kTRUE )
{
  if (sec>0) fsec = sec;
  if (row>-1) frow =row;
  char s[80];
  char sh[80];
  //  create particles histograms
  sprintf(s,"Sector %d   Row %d\n",fsec,frow);  
  sprintf(sh,"Particles%d_%d",fsec,frow);   
  fHParticles = new TH1F(sh,s,size1,4,size1); 

  sprintf(s,"Sector %d   Row %d\n",fsec,frow);  
  sprintf(sh,"All particles%d_%d",fsec,frow);   
  fHAllP = new TH1F(sh,s,200,1,25);  

  sprintf(s,"Sector %d   Row %d\n",fsec,frow);  
  sprintf(sh,"Secondary Particles%d_%d",fsec,frow);   
  fHSecondaryP = new TH1F(sh,s,200,1,25); 

  if (!(fin)) 
  {
    cout<<"Input  file not open, open file before \n";
  }
  else
{
  fin->cd();
  Int_t sectors_by_rows=(Int_t)ftree->GetEntries();
  GParticle * particle;
  //  loop over all sectors and rows
  for (Int_t n=0; n<sectors_by_rows; n++) 
    {
    if (!ftree->GetEvent(n)) continue;
    AliTPCdigit *dig=(AliTPCdigit*)fDigits->UncheckedAt(0);
    if (fsec  < dig->fSector) break;
    if (fsec != dig->fSector) continue;
    if (frow != dig->fPadRow) continue;
    
    Int_t ndigits=fDigits->GetEntriesFast();
    //loop over all digits  in sector pad
    for (Int_t ndig=0; ndig<ndigits; ndig++) 
	{
      Float_t x,y;
	dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
	fHParticles->Fill(dig->fTracks[0]);
	//     get pointer to particle information and fill All and secondary histo
	particle = (GParticle*) fParticles->UncheckedAt(dig->fTracks[0]);
      fHAllP->Fill(particle->GetKF());
      //	        Int_t id = particle->GetKF();
      x = (Float_t)particle->GetVx();
      y = (Float_t)particle->GetVy();
      if ( (x*x+y*y ) > 1 )
	  fHSecondaryP->Fill(particle->GetKF());
	if (all==kTRUE)
	  {
	    fHParticles->Fill(dig->fTracks[1]);
	    particle = (GParticle*) fParticles->UncheckedAt(dig->fTracks[1]);
	    fHAllP->Fill(particle->GetKF()); 
	    x = (Float_t)particle->GetVx();
	    y = (Float_t)particle->GetVy();
	    if ( (x*x+y*y ) > 1 )	    
	      fHSecondaryP->Fill(particle->GetKF());
	    fHParticles->Fill(dig->fTracks[2]);
	    particle = (GParticle*) fParticles->UncheckedAt(dig->fTracks[2]);
	    fHAllP->Fill(particle->GetKF());
	    x = (Float_t)particle->GetVx();
	    y = (Float_t)particle->GetVy();
	    if ( (x*x+y*y ) > 1 )	    
	      fHSecondaryP->Fill(particle->GetKF());
  	  }
  	}
      }
    //make histogram with multiplicity
      char s[80];
      char sh[80];
      if (all==kTRUE)
	sprintf(s,"Number of AliDigits over threshold  per one track in sector %d   Row %d\n (all three most important track recorded)",fsec,frow);
      else
	sprintf(s,"Number of AliDigits over threshold per sector %d   Row %d\n (only most important track)",fsec,frow);  
      sprintf(sh,"His_%d_%d",fsec,frow);   
      fHPartMultiplicity = new TH1F(sh,s,size2,1,size2);
      for (Int_t i=1;i<size1;i++)
	{
	  Int_t mul=Int_t(fHParticles->GetBinContent(i));
          if (mul>0) fHPartMultiplicity->Fill(mul);
	}
      if (fout)
	{
	  fout->cd();          
	  fHParticles->Write();     
      	  fHPartMultiplicity->Write();     
	  //fHSecondaryP->Write(); 
	  //fHAllP->Write();
	}
  }
}


void AliTPCDigitsH::Anal()
{
  if (fbIOState == kFALSE)
    {
      cout<<"Input output problem. \n Do you initialize IO files ? \n";
      return;
    }
  if (fbDigState == kFALSE)
    {
      cout<<"Input file doesn't enhalt selected tree \n";
      return;
    }  
  fout->cd();
 
//if we dont want let histogram in memory then we delete old histogram 
  if ( (fH2Digit) && (fbDelHisto == kTRUE) )  
    {
      //      try 
      //{
	  // delete fH2Digit;
      //}
      //      catch(...)
      //{
      //  cout<<"Delete error. Contact autor of root \n";
      //};
      fH2Digit = 0;
    }
  char s[80];
  char sh[80];
  sprintf(s,"Sector %d   Row %d\n",fsec,frow);  
  sprintf(sh,"h%d_%d",fsec,frow);   

  if  ( (fout) && (fbDelHisto == kFALSE) )
    {
      fH2Digit = (AliH2F *) fout->Get(sh);
      if (fH2Digit) return;     
    } 

  Int_t n_of_pads =fTPCParam->GetNPads(fsec,frow);
     
  fH2Digit = new AliH2F(sh, s, fTimeN, fTimeStart, fTimeStop, n_of_pads, 0, n_of_pads-1);
 
  if (!(fout)) 
  {
    cout<<"Input  file not open, open file before \n";
  }
  else
  {
    //fin->cd();
    Int_t sectors_by_rows=(Int_t)ftree->GetEntries();
    //loop over all sectors and rows
    for (Int_t n=0; n<sectors_by_rows; n++) {
      if (!ftree->GetEvent(n)) continue;
      AliTPCdigit *dig=(AliTPCdigit*)fDigits->UncheckedAt(0);
      if (fsec  < dig->fSector) break;
      if (fsec != dig->fSector) continue;
      if (frow != dig->fPadRow) continue;
      
      Int_t ndigits=fDigits->GetEntriesFast();
      //loop over all digits  in sector pad
      for (Int_t ndig=0; ndig<ndigits; ndig++) {
	dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
	fH2Digit->Fill(dig->fTime,dig->fPad,dig->fSignal);
      }
    }
    if (fout) fout->cd();          
    fH2Digit->Write();          
  }

}
void  AliTPCDigitsH::Draw(Option_t * opt1 ="cont1"  , Option_t * opt2 = "error",
			Option_t * opt3 = "L" )
{
  TString o1 = opt1;
  o1.ToLower();

  TString o2 = opt2;
  o2.ToLower();

  TString o3 = opt3;
  o3.ToLower();

  fcanvas  = new TCanvas("dh","Digits Histograms",700,900);
  
  if (fTitle) 
    {
      //       try 
      //{
	  //delete fTitle;
      //}
      //      catch(...)
      //{
      //  cout<<"Delete error. Contact autor of root \n";
      //};
      fTitle = 0;
    } 
  fTitle = new TPaveText(0.2,0.96,0.8,0.995);
  fTitle->AddText("Occupancy calculation for TPC");
  fTitle->Draw();
  
  fpad1 = new TPad("pad1","",0.05,0.7,0.95,0.95,21);
  fpad1->Draw();
  fpad2 = new TPad("pad2","",0.05,0.4,0.95,0.65,21);
  fpad2->Draw();
  fpad3 = new TPad("pad3","",0.05,0.05,0.95,0.35,21);
  fpad3->Draw();

  fpad1->cd();
  //  pad1->TPaveText::title.SetSize(0.1);
  if (fH2Digit) 
    { 
      fH2DigitBW->Draw(o1);
      fH2DigitBW->SetXTitle("time bin");
      fH2DigitBW->SetYTitle("pad number");
    }
  fpad2->cd();
  fH1Occu->Fit("pol0");
  if (fH1Occu) 
    {
      fH1Occu->Draw(o2);
      fH1Occu->SetXTitle("time bin");
      fH1Occu->SetYTitle("occupancy");
    }
  fpad3->cd();  
  if (fH1Digit)
    {
      fH1Digit->Draw(o3);
      fH1Digit->SetXTitle("time bin");
      fH1Digit->SetYTitle("ADC amplitude ");
    }

};

void  AliTPCDigitsH::DeleteHisto(const Text_t *namecycle)
{
  if (fout) fout->Delete(namecycle);
}
     
void  AliTPCDigitsH::SetHisto(Int_t pad = 1 )
{

  Int_t n_of_pads = fTPCParam->GetNPads(fsec,frow);
  if (pad > (n_of_pads-1)) 
    {
      cout<<"Pad number is greater then actula number of pads in thi row \n";
      cout<<"Noch einmal \n";
      return;
    }
      
  fpad = pad;
  Anal();

  //  if ( (fH1Digit) && (fbDelHisto == kTRUE)) 
  //    {
  //      try 
  //        {
  //	  // delete fH1Digit;
  //	}
  //      catch(...)
  //	{
  //	  cout<<"Delete error. Contact autor of root \n";
  //	};
  //      fH1Digit = 0;
  //    }

  char s[80];
  char sh[80];
  sprintf(s,"example sector %d   Row %d  Pad %d",fsec,frow,fpad);  
  sprintf(sh,"h%d_%d_%d",fsec,frow,fpad);   
 fH2DigitBW = new AliH2F("bw", "", fTimeN, fTimeStart, fTimeStop, n_of_pads, 0, n_of_pads-1);
  fH1Digit = new TH1F(sh,s,fTimeN,fTimeStart,fTimeStop);

  for (Int_t i = 0;i<fTimeN;i++)
    {
      Int_t index = fH2Digit->GetBin(i,pad);
      Float_t weight = fH2Digit->GetBinContent(index);
      fH1Digit->Fill(i,weight);
    };  

 
  sprintf(s,"Occupancy in sector %d   Row %d  threshold = %d",fsec,frow,fThreshold);  
  sprintf(sh,"hoccu%d_%d_%d",fsec,frow,fpad);   
  fH1Occu = new TH1F(sh,s,fOccuN,fTimeStart,fTimeStop);
  
  for (Int_t i = 0;i<fOccuN;i++)
    {
      Int_t over =0;  
      Int_t all  =0;    
      for (int itime = i*(fTimeN/fOccuN); itime<(i+1)*(fTimeN/fOccuN);itime++)
	{
	  for (Int_t ipad = 0; ipad < n_of_pads; ipad++)
	    {
	      Int_t index = fH2Digit->GetBin(itime,ipad);
	      if ( (ipad>3) && ((ipad+3)<n_of_pads)){
		all++;
		if (fH2Digit->GetBinContent(index) >fThreshold) over++ ;
	      }
              if (fH2Digit->GetBinContent(index) >fThreshold) 
		fH2DigitBW->Fill(itime,ipad,1);
	      else
		fH2DigitBW->Fill(itime,ipad,0);
	    }
	}
     Float_t occu = ((Float_t)over) /((Float_t) (all)); 
     //     Float_t time =   ((fTimeStop-fTimeStart)/fOccuN)*i+fTimeStart;
   
     fH1Occu->SetBinContent(i,occu); 
     // Int_t index = fH1Occu->GetBin(i);
     Float_t error = sqrt( ((Float_t) ((over)/25+1)) )/((Float_t)(all)/25.);
     fH1Occu->SetBinError(i,error);

    };  
 

}
  


void AliTPCDigitsH::Streamer(TBuffer & R__b)
{
  if (R__b.IsReading()) {
    //      Version_t R__v = R__b.ReadVersion();
   } else {
      R__b.WriteVersion(AliTPCDigitsH::IsA());    
   } 
}
