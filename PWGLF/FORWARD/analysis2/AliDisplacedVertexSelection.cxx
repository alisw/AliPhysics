#include "AliDisplacedVertexSelection.h"
#include <iostream>
#include <TROOT.h>
#include "AliESDEvent.h"
#include "AliESDZDC.h"
ClassImp(AliDisplacedVertexSelection)
#if 0
; // This is for Emacs
#endif

//____________________________________________________________________
AliDisplacedVertexSelection::AliDisplacedVertexSelection()
  : TObject()
{
}
//____________________________________________________________________
AliDisplacedVertexSelection::AliDisplacedVertexSelection(const AliDisplacedVertexSelection& o)
  : TObject(o)
{
}
//____________________________________________________________________
AliDisplacedVertexSelection&
AliDisplacedVertexSelection::operator=(const AliDisplacedVertexSelection& o)
{
  if (&o == this) return *this; 
  return *this;
}

//____________________________________________________________________
void
AliDisplacedVertexSelection::Output(TList* /*l*/, const char* /* name*/) const
{
}
  
//____________________________________________________________________
void
AliDisplacedVertexSelection::Print(Option_t*) const
{
#if 0
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << std::boolalpha 
	    << std::noboolalpha << std::endl;
#endif
}
//____________________________________________________________________

Double_t
AliDisplacedVertexSelection::CheckDisplacedVertex(const AliESDEvent* esd) const {
  // Code taken directly from M.Guilbaud.
  Double_t zvtx = 9999.;
  
  Float_t kZDCrefSum        =-66.5;
  Float_t kZDCrefDelta      =-2.10;
  Float_t kZDCsigmaSum      = 3.25;
  Float_t kZDCsigmaDelta    = 2.25;
  Float_t kZDCsigmaSumSat   = 2.50;
  Float_t kZDCsigmaDeltaSat = 1.35;  
  
  AliESDZDC* esdZDC = esd -> GetESDZDC();
  
  /* Double_t zdcNCEnergy  = esdZDC->GetZDCN1Energy();
     Double_t fZpcEnergy  = esdZDC->GetZDCP1Energy();
     Double_t fZnaEnergy  = esdZDC->GetZDCN2Energy();
     Double_t fZpaEnergy  = esdZDC->GetZDCP2Energy();
     Double_t fZem1Energy = esdZDC->GetZDCEMEnergy(0);
     Double_t fZem2Energy = esdZDC->GetZDCEMEnergy(1);*/

  //Double_t fzdcEn      = (esdZDC->GetZDCN1Energy()+esdZDC->GetZDCP1Energy()+esdZDC->GetZDCN2Energy()+esdZDC->GetZDCP2Energy());
  //Double_t fzemEn      = (esdZDC->GetZDCEMEnergy(0)+esdZDC->GetZDCEMEnergy(1))/8.;
  

///////////////////
//Event Selection//
///////////////////
 Double_t deltaTdc = 0;
 Double_t sumTdc   = 0;
   
for(Int_t i = 0; i < 4; ++i) 
   {
        if(esdZDC->GetZDCTDCData(10,i) != 0) 
	  {
	      Double_t  tdcCnoCorr = 0.025*(esdZDC->GetZDCTDCData(10,i)-esdZDC->GetZDCTDCData(14,i)); 
	      Double_t  tdcC       = esdZDC->GetZDCTDCCorrected(10,i); 
	       for(Int_t j = 0; j < 4; ++j) 
	          {
	               if(esdZDC->GetZDCTDCData(12,j) != 0) 
		         {
	                    Double_t   tdcAnoCorr = 0.025*(esdZDC->GetZDCTDCData(12,j)-esdZDC->GetZDCTDCData(14,j));
	                     Double_t   tdcA       = esdZDC->GetZDCTDCCorrected(12,j);
			      if(esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled)) 
			        {	    
				      deltaTdc = tdcC-tdcA;
				      sumTdc   = tdcC+tdcA;
	                        }
	                      else
			        {
				      deltaTdc = tdcCnoCorr-tdcAnoCorr;
				      sumTdc   = tdcCnoCorr+tdcAnoCorr;
				}
	                 }
	          }
          }
   }
//if(TMath::Abs(deltaTdc) > 0)
//std::cout<<deltaTdc<<"   "<<sumTdc<<std::endl;
///////////////////////
//Global Event Filter//
///////////////////////
 Bool_t zdcAccSat[21];
 // Bool_t zdcAccSatRunClass[21];
 for(Int_t k = -10; k <= 10; k++) 
   {
      zdcAccSat[k+10]         = kFALSE;
      // zdcAccSatRunClass[k+10] = kFALSE;
   }


if(deltaTdc!=0. || sumTdc!=0.)
  {

     for(Int_t k = -10; k <= 10; ++k) 
        {
	     Float_t zsat = 2.55*(Float_t)k;
		
	     if(k==0)
	       {
		     if(((deltaTdc-kZDCrefDelta)*(deltaTdc-kZDCrefDelta)/(kZDCsigmaDelta*kZDCsigmaDelta)+
		         (sumTdc  -kZDCrefSum  )*(sumTdc  -kZDCrefSum  )/(kZDCsigmaSum  *kZDCsigmaSum))<=1.0)
		       {
		           zdcAccSat[k+10] = kTRUE;
		       }	
	       }
	     else
	       {
		     if(((deltaTdc-kZDCrefDelta-zsat)*(deltaTdc-kZDCrefDelta-zsat)/(kZDCsigmaDeltaSat*kZDCsigmaDeltaSat)+
		         (sumTdc  -kZDCrefSum  -zsat)*(sumTdc  -kZDCrefSum  -zsat)/(kZDCsigmaSumSat  *kZDCsigmaSumSat  ))<=1.0)
		       {
		           zdcAccSat[k+10] = kTRUE;
		       }		
	  
	       }
	  
         }
  }
   
 for(Int_t k=-10;k<=10;++k) {
   if(zdcAccSat[k+10] && k!=0 ) {
     //std::cout<<"Displaced vertex at "<<37.5*(Float_t)k<<" cm"<<std::endl; 
     zvtx = 37.5*(Float_t)k;
   }
 }
 
 return zvtx;
 
}
//____________________________________________________________________
Double_t
AliDisplacedVertexSelection::CalculateDisplacedVertexCent(const AliESDEvent* esd) const 
{ 
  Float_t kZDCrefSum        =-66.5;
  Float_t kZDCrefDelta      =-2.10;
  Float_t kZDCsigmaSum      = 3.25;
  Float_t kZDCsigmaDelta    = 2.25;
  Float_t kZDCsigmaSumSat   = 2.50;
  Float_t kZDCsigmaDeltaSat = 1.35;  
  
  AliESDZDC* esdZDC = esd -> GetESDZDC();
  
  Int_t runnumber = esd->GetRunNumber();
  Double_t fcurrentL3 = esd->GetCurrentL3();
  Double_t fcurrentDipo = esd->GetCurrentDip();
  Double_t cent = -1;
  
  /*Double_t zdcNCEnergy  = esdZDC->GetZDCN1Energy();
  Double_t fZpcEnergy  = esdZDC->GetZDCP1Energy();
  Double_t fZnaEnergy  = esdZDC->GetZDCN2Energy();
  Double_t fZpaEnergy  = esdZDC->GetZDCP2Energy();
  Double_t fZem1Energy = esdZDC->GetZDCEMEnergy(0);
  Double_t fZem2Energy = esdZDC->GetZDCEMEnergy(1);*/
  
  Double_t fzdcEn      = (esdZDC->GetZDCN1Energy()+esdZDC->GetZDCP1Energy()+esdZDC->GetZDCN2Energy()+esdZDC->GetZDCP2Energy());
  Double_t fzemEn      = (esdZDC->GetZDCEMEnergy(0)+esdZDC->GetZDCEMEnergy(1))/8.;
  

///////////////////
//Event Selection//
///////////////////
 Double_t deltaTdc = 0;
 Double_t sumTdc   = 0;
   
for(Int_t i = 0; i < 4; ++i) 
   {
        if(esdZDC->GetZDCTDCData(10,i) != 0) 
	  {
	      Double_t  tdcCnoCorr = 0.025*(esdZDC->GetZDCTDCData(10,i)-esdZDC->GetZDCTDCData(14,i)); 
	      Double_t  tdcC       = esdZDC->GetZDCTDCCorrected(10,i); 
	       for(Int_t j = 0; j < 4; ++j) 
	          {
	               if(esdZDC->GetZDCTDCData(12,j) != 0) 
		         {
	                    Double_t   tdcAnoCorr = 0.025*(esdZDC->GetZDCTDCData(12,j)-esdZDC->GetZDCTDCData(14,j));
	                     Double_t   tdcA       = esdZDC->GetZDCTDCCorrected(12,j);
			      if(esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled)) 
			        {	    
				      deltaTdc = tdcC-tdcA;
				      sumTdc   = tdcC+tdcA;
	                        }
	                      else
			        {
				      deltaTdc = tdcCnoCorr-tdcAnoCorr;
				      sumTdc   = tdcCnoCorr+tdcAnoCorr;
				}
	                 }
	          }
          }
   }
//if(TMath::Abs(deltaTdc) > 0)
//std::cout<<deltaTdc<<"   "<<sumTdc<<std::endl;
///////////////////////
//Global Event Filter//
///////////////////////
 Bool_t zdcAccSat[21];
 // Bool_t zdcAccSatRunClass[21];
for(Int_t k = -10; k <= 10; k++) 
   {
      zdcAccSat[k+10]         = kFALSE;
      // zdcAccSatRunClass[k+10] = kFALSE;
   }


if(deltaTdc!=0. || sumTdc!=0.)
  {

     for(Int_t k = -10; k <= 10; ++k) 
        {
	     Float_t zsat = 2.55*(Float_t)k;
		
	     if(k==0)
	       {
		     if(((deltaTdc-kZDCrefDelta)*(deltaTdc-kZDCrefDelta)/(kZDCsigmaDelta*kZDCsigmaDelta)+
		         (sumTdc  -kZDCrefSum  )*(sumTdc  -kZDCrefSum  )/(kZDCsigmaSum  *kZDCsigmaSum))<=1.0)
		       {
		           zdcAccSat[k+10] = kTRUE;
		       }	
	       }
	     else
	       {
		     if(((deltaTdc-kZDCrefDelta-zsat)*(deltaTdc-kZDCrefDelta-zsat)/(kZDCsigmaDeltaSat*kZDCsigmaDeltaSat)+
		         (sumTdc  -kZDCrefSum  -zsat)*(sumTdc  -kZDCrefSum  -zsat)/(kZDCsigmaSumSat  *kZDCsigmaSumSat  ))<=1.0)
		       {
		           zdcAccSat[k+10] = kTRUE;
		       }		
	  
	       }
	  
         }
  }

 Float_t kZEMcorrPlusPlus[21]   = {0.6840,0.7879,0.8722,0.9370,0.9837,1.0137,
				   1.0292,1.0327,1.0271,1.0152,1.0000,0.9844,
				   0.9714,0.9634,0.9626,0.9708,0.9891,1.0181,
				   1.0574,1.1060,1.1617};
 Float_t kZEMcorrMoinsMoins[21] = {0.7082,0.8012,0.8809,0.9447,0.9916,1.0220,
				   1.0372,1.0395,1.0318,1.0174,1.0000,0.9830,
				   0.9699,0.9635,0.9662,0.9794,1.0034,1.0371,
				   1.0782,1.1224,1.1634};
 ///////////////////
 //ZemEn correction//
 ////////////////////

 for(Int_t k = -10; k <= 10; ++k)
   {        
     if(zdcAccSat[k+10])
       { 
	 //std::cout<<"Vertex at "<<37.5*(Float_t)k<<"cm from cent"<<std::endl;
         if(fcurrentDipo>0 && fcurrentL3>0)
	   {
	     fzemEn /= kZEMcorrPlusPlus[k+10];
	   }
	 if(fcurrentDipo<0 && fcurrentL3<0)
	   {
	     fzemEn /= kZEMcorrMoinsMoins[k+10];
	   }
       }          
   }
 
 ////////////////////////
 //Centrality selection//
 ////////////////////////
 Double_t fzdcPercentile = -1;
 Float_t slope;
 
 
 if(runnumber < 137722 && runnumber >= 137161 )
   {
     if(fzemEn > 345.)
       {
	 slope = (fzdcEn + 16992.)/(fzemEn - 345.);
	 slope += 2.23902e+02;
	 fzdcPercentile = (TMath::ATan(slope) - 1.56667)/9.49434e-05;
	 if (fzdcPercentile<0) fzdcPercentile = 0;
       }
     else fzdcPercentile = 100;                }
 else if(runnumber < 139172 && runnumber >= 137722)
   {
     if(fzemEn > 295.)
       {
	 slope = (fzdcEn + 15000.)/(fzemEn - 295.);
	 slope += 2.23660e+02;
	 fzdcPercentile = (TMath::ATan(slope) - 1.56664)/8.99571e-05;
	 if (fzdcPercentile<0) fzdcPercentile = 0;
       }                            else fzdcPercentile = 100;          }


 else if(runnumber >= 139172)
   {
     if(fzemEn > 345.)
       {
	 slope = (fzdcEn + 16992.)/(fzemEn - 345.);
	 slope += 2.04789e+02;
	 fzdcPercentile = (TMath::ATan(slope) - 1.56629)/1.02768e-04;
	 if (fzdcPercentile<0) fzdcPercentile = 0;
       }                               else fzdcPercentile = 100;           }
 else
   {
     if(fzemEn > 345.)
       {
	 slope = (fzdcEn + 16992.)/(fzemEn - 345.);
	 slope += 2.04789e+02;
	 fzdcPercentile = (TMath::ATan(slope) - 1.56629)/1.02768e-04;
	 if(fzdcPercentile<0) fzdcPercentile = 0;
       }
     else fzdcPercentile = 100;                  }
 
 if(fzdcPercentile > 0 && fzdcPercentile <100) {
   //std::cout<<fzdcPercentile<<std::endl;
   cent = fzdcPercentile;
 }
 return cent;
}
//____________________________________________________________________
//
// EOF
//
