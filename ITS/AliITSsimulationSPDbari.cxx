#include <iostream.h>
#include <TRandom.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TParticle.h>


#include "AliRun.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSMapA2.h"
#include "AliITSsimulationSPDbari.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"


ClassImp(AliITSsimulationSPDbari)
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Rocco Caliandro
// from a model developed with T. Virgili and R.A. Fini
// June 15 2000
//
// AliITSsimulationSPD is the simulation of SPDs
//
//________________________________________________________________________

AliITSsimulationSPDbari::AliITSsimulationSPDbari(){
  // constructor
  fResponse     = 0;
  fSegmentation = 0;
  fHis          = 0;
  fThresh       = 0.;
  fSigma        = 0.;
  fCouplCol     = 0.;
  fCouplRow     = 0.;
}
//_____________________________________________________________________________

AliITSsimulationSPDbari::AliITSsimulationSPDbari(AliITSsegmentation *seg, AliITSresponse *resp) {
  // constructor
      fResponse = resp;
      fSegmentation = seg;

      fResponse->Thresholds(fThresh,fSigma);
      fResponse->GetNoiseParam(fCouplCol,fCouplRow);
      
      fMapA2 = new AliITSMapA2(fSegmentation);
   
      // 
      fNPixelsZ=fSegmentation->Npz();
      fNPixelsX=fSegmentation->Npx();

}

//_____________________________________________________________________________

AliITSsimulationSPDbari::~AliITSsimulationSPDbari() { 
  // destructor

  delete fMapA2;

  if (fHis) {
     fHis->Delete(); 
     delete fHis;     
  }                
}

//__________________________________________________________________________
AliITSsimulationSPDbari::AliITSsimulationSPDbari(const AliITSsimulationSPDbari &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fMapA2 = source.fMapA2;
  this->fThresh = source.fThresh;
  this->fSigma = source.fSigma;
  this->fCouplCol = source.fCouplCol;
  this->fCouplRow = source.fCouplRow;
  this->fNPixelsX = source.fNPixelsX;
  this->fNPixelsZ = source.fNPixelsZ;
  this->fHis = source.fHis;
  return;
}

//_________________________________________________________________________
AliITSsimulationSPDbari& 
  AliITSsimulationSPDbari::operator=(const AliITSsimulationSPDbari &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fMapA2 = source.fMapA2;
  this->fThresh = source.fThresh;
  this->fSigma = source.fSigma;
  this->fCouplCol = source.fCouplCol;
  this->fCouplRow = source.fCouplRow;
  this->fNPixelsX = source.fNPixelsX;
  this->fNPixelsZ = source.fNPixelsZ;
  this->fHis = source.fHis;
  return *this;
  }
//_____________________________________________________________________________

void AliITSsimulationSPDbari::DigitiseModule(AliITSmodule *mod, Int_t module,
                                             Int_t dummy) {
  // digitize module


  TObjArray *fHits = mod->GetHits();
  Int_t nhits = fHits->GetEntriesFast();
  if (!nhits) return;


  printf("Module %d (%d hits) \n",module+1,nhits);


  Int_t  number=10000;
  Int_t    *frowpixel = new Int_t[number];
  Int_t    *fcolpixel = new Int_t[number];
  Double_t *fenepixel = new Double_t[number];

  // Array of pointers to store the track index of the digits
  // leave +1, otherwise pList crashes when col=256, row=192 
    Int_t maxNdigits = fNPixelsX*fNPixelsZ+fNPixelsX+1;
  Float_t  **pList = new Float_t* [maxNdigits];
  memset(pList,0,sizeof(Float_t*)*maxNdigits);


  // noise setting
  SetFluctuations(pList);


  // loop over hits in the module
  Int_t hitpos;
      for (hitpos=0;hitpos<nhits;hitpos++) {  
 
     HitToDigit(mod,hitpos,module,frowpixel,fcolpixel,fenepixel,pList);


  }// end loop over digits

     CreateDigit(nhits,module,pList);

  // clean memory
  delete[] frowpixel;
  delete[] fcolpixel;
  delete[] fenepixel;
  fMapA2->ClearMap();
  delete [] pList;

}
//_____________________________________________________________________________

void AliITSsimulationSPDbari::UpdateMap( Int_t row, Int_t col, Double_t ene) {
//
// updates the Map of signal, adding the energy  (ene) released by the current track
//
      Double_t signal; 
      signal = fMapA2->GetSignal(row,col);
      signal += ene;
      fMapA2->SetHit(row,col,signal);
                                         
 }
//_____________________________________________________________________________  
void AliITSsimulationSPDbari::HitToDigit(AliITSmodule *mod, Int_t hitpos, Int_t module, 
                                        Int_t *frowpixel, Int_t *fcolpixel,
					Double_t *fenepixel, Float_t **pList) {
//
//  Steering function to determine the digits associated to a given hit (hitpos)
//  The digits are created by charge sharing (ChargeSharing) and by
//  capacitive coupling (SetCoupling). At all the created digits is associated
//  the track number of the hit (ntrack)
//


   static Float_t x1l,y1l,z1l;
   Float_t x2l,y2l,z2l,etot;
   Int_t layer,r1,r2,c1,c2,row,col,npixel = 0;
   Int_t ntrack;
   Double_t ene;
   const Float_t kconv = 10000.;     // cm -> microns

   TObjArray *fHits = mod->GetHits();
   AliITShit *hit = (AliITShit*) fHits->At(hitpos);
   layer = hit->GetLayer();
   etot=hit->GetIonization();
   ntrack=hit->GetTrack();


        if (hit->GetTrackStatus()==66) {
	      hit->GetPositionL(x1l,y1l,z1l);
          // positions shifted and converted in microns 
          x1l = x1l*kconv + fSegmentation->Dx()/2.;
          z1l = z1l*kconv + fSegmentation->Dz()/2.;
        }
        else {
	      hit->GetPositionL(x2l,y2l,z2l);	      
          // positions  shifted and converted in microns
          x2l = x2l*kconv + fSegmentation->Dx()/2.;
          z2l = z2l*kconv + fSegmentation->Dz()/2.;


          // to account for 83750 effective sensitive area
          // introduced in geometry (Dz=83600)
          if (z1l>fSegmentation->Dz()) return;
          if (z2l>fSegmentation->Dz()) return;
          // to preserve for hit having x>12800
          if (x1l>fSegmentation->Dx()) return;
          if (x2l>fSegmentation->Dx()) return;

          //Get the col and row number starting from 1
          // the x direction is not inverted for the second layer!!!
	      fSegmentation->GetPadIxz(x1l, z1l, c1, r1); 
	      fSegmentation->GetPadIxz(x2l, z2l, c2, r2);

          // to account for unexpected equal entrance and 
          // exit coordinates
          if (x1l==x2l) x2l=x2l+x2l*0.000001;
          if (z1l==z2l) z2l=z2l+z2l*0.000001;

          // to account for tracks at the edge of the sensitive area
          // of SPDs
          if (x1l<0) return;
          if (x2l<0) return;
          if (z1l<0) return;
          if (z2l<0) return;
	      

	      if ((r1==r2) && (c1==c2)) 
	      {
             // no charge sharing
	         npixel = 1;		 
		     frowpixel[npixel-1] = r1;
		     fcolpixel[npixel-1] = c1;
  		     fenepixel[npixel-1] = etot;
          }
	      else {
             // charge sharing
	         ChargeSharing(x1l,z1l,x2l,z2l,c1,r1,c2,r2,etot,
		         	       npixel,frowpixel,fcolpixel,fenepixel);

          }
                  

      	  for (Int_t npix=0;npix<npixel;npix++)
	      {
		   row = frowpixel[npix];
		   col = fcolpixel[npix];
		   ene = fenepixel[npix];
		   UpdateMap(row,col,ene);                   
		   GetList(ntrack,pList,row,col); 
		   // Starting capacitive coupling effect
		   SetCoupling(row,col,ntrack,pList); 
	      }
	    x1l=x2l;
	    y1l=y2l;
	    z1l=z2l;	     	     
        }
}

//_________________________________________________________________________

void AliITSsimulationSPDbari::ChargeSharing(Float_t x1l,Float_t z1l,Float_t x2l,
                    Float_t z2l,Int_t c1,Int_t r1,Int_t c2,
				    Int_t r2,Float_t etot,
				    Int_t &npixel,Int_t *frowpixel,
				    Int_t *fcolpixel,Double_t *fenepixel){
  //
  //  Take into account the geometrical charge sharing when the track
  //  crosses more than one pixel.
  //
  //Begin_Html
  /*
  <img src="picts/ITS/barimodel_2.gif">
  </pre>
  <br clear=left>
  <font size=+2 color=red>
  <a href="mailto:Rocco.Caliandro@ba.infn.it"></a>.
  </font>
  <pre>
  */
  //End_Html


   Float_t xa,za,xb,zb,dx,dz,dtot,dm,refr,refm,refc;
   Float_t refn=0.;
   Float_t arefm, arefr, arefn, arefc, azb, az2l, axb, ax2l;
   Int_t   dirx,dirz,rb,cb;


   Int_t flag,flagrow,flagcol;
  
   Double_t epar;


   npixel = 0;
   xa = x1l;
   za = z1l;
   dx = TMath::Abs(x1l-x2l);
   dz = TMath::Abs(z1l-z2l);
   dtot = TMath::Sqrt((dx*dx)+(dz*dz));   
   dm = (x2l - x1l) / (z2l - z1l);

   dirx = (Int_t) ((x2l - x1l) / dx);
   dirz = (Int_t) ((z2l - z1l) / dz);
   
   
   // calculate the x coordinate of  the pixel in the next column    
   // and the z coordinate of  the pixel in the next row    

   Float_t xpos, zpos;

   fSegmentation->GetPadCxz(c1, r1-1, xpos, zpos); 

   Float_t xsize = fSegmentation->Dpx(0);
   Float_t zsize = fSegmentation->Dpz(r1);

   if (dirx == 1) refr = xpos+xsize/2.;
             else refr = xpos-xsize/2.;

   if (dirz == 1) refn = zpos+zsize/2.;
             else refn = zpos-zsize/2.;

   
   flag = 0;
   flagrow = 0;
   flagcol = 0;
   do
   {
       
      // calculate the x coordinate of the intersection with the pixel
      // in the next cell in row  direction

      refm = (refn - z1l)*dm + x1l;
   
      // calculate the z coordinate of the intersection with the pixel
      // in the next cell in column direction 

      refc = (refr - x1l)/dm + z1l;
      
      
      arefm = refm * dirx;
      arefr = refr * dirx;
      arefn = refn * dirz;
      arefc = refc * dirz;
            

      if ((arefm < arefr) && (arefn < arefc)){
        	 
         // the track goes in the pixel in the next cell in row direction
	     xb = refm;
	     zb = refn;
	     cb = c1;
	     rb = r1 + dirz;
	     azb = zb * dirz;
         az2l = z2l * dirz;
	     if (rb == r2) flagrow=1;
	     if (azb > az2l) {
	        zb = z2l;
	        xb = x2l;
	     }     

         // shift to the pixel in the next cell in row direction
         Float_t zsizeNext = fSegmentation->Dpz(rb);
	     refn += zsizeNext*dirz;

      }
      else {
         
         // the track goes in the pixel in the next cell in column direction
	     xb = refr;
	     zb = refc;
	     cb = c1 + dirx;
	     rb = r1;
	     axb = xb * dirx;
         ax2l = x2l * dirx;
         if (cb == c2) flagcol=1;
	     if (axb > ax2l) {
	        zb = z2l;
	        xb = x2l;
	     }

         // shift to the pixel in the next cell in column direction
         Float_t xsizeNext = fSegmentation->Dpx(cb);
	     refr += xsizeNext*dirx;
        
      }
      
      //calculate the energy lost in the crossed pixel      
      epar = TMath::Sqrt((xb-xa)*(xb-xa)+(zb-za)*(zb-za)); 
      epar = etot*(epar/dtot);

      //store row, column and energy lost in the crossed pixel


      frowpixel[npixel] = r1;
      fcolpixel[npixel] = c1;
      fenepixel[npixel] = epar;
      npixel++;
 
      // the exit point of the track is reached
      if (epar == 0) flag = 1;
      if ((r1 == r2) && (c1 == c2)) flag = 1;
      if (flag!=1) {
        r1 = rb;
        c1 = cb;
        xa = xb;
        za = zb;
      }
   
   } while (flag==0);

}
//___________________________________________________________________________
void AliITSsimulationSPDbari::SetCoupling(Int_t row, Int_t col, Int_t ntrack,
                                          Float_t **pList) {
   //
   //  Take into account the coupling between adiacent pixels.
   //  The parameters probcol and probrow are the fractions of the
   //  signal in one pixel shared in the two adjacent pixels along
   //  the column and row direction, respectively.
   //
   //Begin_Html
   /*
   <img src="picts/ITS/barimodel_3.gif">
   </pre>
   <br clear=left>
   <font size=+2 color=red>
   <a href="mailto:Rocco.Caliandro@ba.infn.it"></a>.
   </font>
   <pre>
   */
   //End_Html


   Int_t j1,j2,flag=0;
   Double_t pulse1,pulse2;
                              

   j1 = row;
   j2 = col;
  
   pulse1 = fMapA2->GetSignal(row,col);
   pulse2 = pulse1;

   for (Int_t isign=-1;isign<=1;isign+=2)
   {

// loop in row direction
      
      do
      {
         j1 += isign;
         pulse1 *= fCouplRow;                  
      
         if ((j1 < 0) || (j1 > fNPixelsZ-1) || (pulse1 < fThresh))
         { 
	       pulse1 = fMapA2->GetSignal(row,col);
	       j1 = row;
	       flag = 1;
         }
          else{                
		   UpdateMap(j1,col,pulse1);                   
		   GetList(ntrack,pList,j1,col); 
           flag = 0;
	     }
	 
      } while(flag == 0);          
      
      
// loop in column direction
      
      do
      {
         j2 += isign;
         pulse2 *= fCouplCol;                  
      
         if ((j2 < 0) || (j2 > (fNPixelsX-1)) || (pulse2 < fThresh))
         {                
	       pulse2 = fMapA2->GetSignal(row,col);
	       j2 = col;
	       flag = 1;
         }
          else{                
		   UpdateMap(row,j2,pulse2);                   
		   GetList(ntrack,pList,row,j2); 
           flag = 0;
	     }
	 
      } while(flag == 0);          
   
   }

}
//___________________________________________________________________________
void AliITSsimulationSPDbari::CreateDigit(Int_t nhits, Int_t module, Float_t
**pList) {                                   
  //
  // The pixels are fired if the energy deposited inside them is above
  // the threshold parameter ethr. Fired pixed are interpreted as digits
  // and stored in the file digitfilename.
  //

   AliITS *aliITS  = (AliITS*)gAlice->GetModule("ITS");   
 
   
   Int_t digits[3];
   Int_t tracks[3];
   Int_t hits[3];
   Float_t charges[3]; 
   Int_t gi,j1;
   

   if (nhits > 0) {
    
     for (Int_t r=1;r<=fNPixelsZ;r++) {
        for (Int_t c=1;c<=fNPixelsX;c++) {
   
           // check if the deposited energy in a pixel is above the threshold 
           Float_t signal = (Float_t) fMapA2->GetSignal(r,c);
           if ( signal > fThresh) {
	          digits[0] = r-1;  // digits starts from 0
              digits[1] = c-1;  // digits starts from 0
              digits[2] = 1;  
	          gi =r*fNPixelsX+c; // global index
	          for(j1=0;j1<3;j1++){
                 tracks[j1] = (Int_t)(*(pList[gi]+j1));
                 charges[j1] = 0;
              }
              Float_t phys = 0;        
	          aliITS->AddSimDigit(0,phys,digits,tracks,hits,charges);
	          if(pList[gi]) delete [] pList[gi];
	       }//endif of threshold condition
        }
     }// enddo on pixels
    }
    
}
//_____________________________________________________________________________

void AliITSsimulationSPDbari::GetList(Int_t label,Float_t **pList,
                                      Int_t row, Int_t col) {
  // loop over nonzero digits

  Int_t ix = col;
  Int_t iz = row;
  Int_t globalIndex;
  Float_t signal;
  Float_t highest,middle,lowest;

          
  signal=fMapA2->GetSignal(iz,ix);

  globalIndex = iz*fNPixelsX+ix; // globalIndex starts from 1


  if(!pList[globalIndex])
  {
     // 
     // Create new list (6 elements - 3 signals and 3 tracks + total sig)
     //

     pList[globalIndex] = new Float_t [6];


     // set list to -2 
     *(pList[globalIndex]) = -2.;
     *(pList[globalIndex]+1) = -2.;
     *(pList[globalIndex]+2) = -2.;
     *(pList[globalIndex]+3) =  0.;
     *(pList[globalIndex]+4) =  0.;
     *(pList[globalIndex]+5) =  0.;

     *pList[globalIndex] = (float)label;
     *(pList[globalIndex]+3) = signal;
  }
  else{


	  // check the signal magnitude
      highest = *(pList[globalIndex]+3);
      middle  = *(pList[globalIndex]+4);
      lowest  = *(pList[globalIndex]+5);


      signal -= (highest+middle+lowest);


	  //
	  //  compare the new signal with already existing list
	  //
      if(signal<lowest) return; // neglect this track

      if (signal>highest)
      {
         *(pList[globalIndex]+5) = middle;
         *(pList[globalIndex]+4) = highest;
         *(pList[globalIndex]+3) = signal;
         *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
         *(pList[globalIndex]+1) = *pList[globalIndex];
         *(pList[globalIndex]) = label;
	  }
        else if (signal>middle)
      {
         *(pList[globalIndex]+5) = middle;
         *(pList[globalIndex]+4) = signal;
         *(pList[globalIndex]+2) = *(pList[globalIndex]+1);
         *(pList[globalIndex]+1) = label;
	  }
        else
      {
         *(pList[globalIndex]+5) = signal;
         *(pList[globalIndex]+2) = label;
	  }
  }    
}
//_________________________________________________________________________ 
void AliITSsimulationSPDbari::SetFluctuations(Float_t **pList) {
  //
  //  Set the electronic noise and threshold non-uniformities to all the
  //  pixels in a detector.
  //  The parameter fSigma is the squared sum of the sigma due to noise
  //  and the sigma of the threshold distribution among pixels.
  //
  //Begin_Html
  /*
  <img src="picts/ITS/barimodel_1.gif">
  </pre>
  <br clear=left>
  <font size=+2 color=red>
  <a href="mailto:Rocco.Caliandro@ba.infn.it"></a>.
  </font>
  <pre>
  */
  //End_Html
  
  
  TRandom random; 
  Double_t signal;

  Int_t iz,ix;
  for(iz=1;iz<=fNPixelsZ;iz++){
    for(ix=1;ix<=fNPixelsX;ix++){
      signal = fSigma*random.Gaus(); 
      fMapA2->SetHit(iz,ix,signal);

      // insert in the label-signal list the pixels fired only by noise
      if ( signal > fThresh) {
        Int_t globalIndex = iz*fNPixelsX+ix; 
        pList[globalIndex] = new Float_t [6];
        *(pList[globalIndex]) = -2.;
        *(pList[globalIndex]+1) = -2.;
        *(pList[globalIndex]+2) = -2.;
        *(pList[globalIndex]+3) =  signal;
        *(pList[globalIndex]+4) =  0.;
        *(pList[globalIndex]+5) =  0.;
      }
    } // end of loop on pixels
  } // end of loop on pixels
  
 }
//____________________________________________

void AliITSsimulationSPDbari::CreateHistograms() {
  // CreateHistograms

      Int_t i;
      for(i=0;i<fNPixelsZ;i++) {
	   TString spdname("spd_");
	   Char_t candnum[4];
	   sprintf(candnum,"%d",i+1);
	   spdname.Append(candnum);
	   (*fHis)[i] = new TH1F(spdname.Data(),"SPD maps",
                              fNPixelsX,0.,(Float_t) fNPixelsX);
      }

}

//____________________________________________

void AliITSsimulationSPDbari::ResetHistograms() {
    //
    // Reset histograms for this detector
    //
    Int_t i;
    for(i=0;i<fNPixelsZ;i++ ) {
	if ((*fHis)[i])    ((TH1F*)(*fHis)[i])->Reset();
    }

}
