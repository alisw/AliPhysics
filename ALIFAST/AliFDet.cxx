////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFast Detector Class                                                 //
//                                                                        //
// to provide information of effective material (X/Xo) of the detector    //
// needed for the multiple scattering formula used in AliFTrackMaker.     // 
//                                                                        // 
// the number and dimensions of cylindrical layers of material are        //
// initialised here for the TP status and are to be updated accordingly.  //
//                                                                        //
//                                                                        //
// origin: "init_geometry" routine in "res.f" fortran by Karel Safarik    //
//         which was used to calculate the track resolution for TP.       // 
//                                                                        // 
//                                                                        //
// AliFast: E. Richter-Was and Y. Foka                                    //
//          following general structure of Makers in ATLFast              //
//          by R. Brun  and E. R. Was                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliFDet.h"
#include "TMath.h"

ClassImp(AliFDet)


//_____________________________________________________________________________
AliFDet::AliFDet(const char *name, const char *title) : TNamed(name,title)
{
   for(Int_t idDet=0; idDet<kNMaxDet; idDet++){
      fRDet[idDet]       = 0;
      fRDetSQ[idDet]     = 0;
      fThickDet[idDet]   = 0;
      fErrorRPhi[idDet]  = 0;
      fErrorZ[idDet]     = 0;
      fErrorR[idDet]     = 0;
      fIFlagDet[idDet]   = 0;
      fIFlagGas[idDet]   = 0;
   }
   fBMag       = 0;
   fConstMag   = 0;
   fNDetActive = 0;
   fNDet       = 0;
}

//_____________________________________________________________________________
void AliFDet::InitDetParam()
{
  //initialisation of the detector material to the TP status.
  //needed for multiple scattering formula used in AliFTrackMaker 
  //for track resolution calculation
  //  
  //errorRPhi, errorZ:    
  //the errors in bending and r direction are due to detector precision and alignement
  //the error  in radial direction r is due to alignement only
  //the errors are momentum dependent; for iFlagDet=2 as for TPC are calculated properly
  //
  //fErrorVertexX,fErrorVertexY, fErrorVertexZ
  //errors of vertex  
  //the vertex precision depends on particle multiplicity mult_density
  //optimistic errors for/high multiplicity
 
   Double_t   rDet[kNMaxDet];         // radius of detector material in cm
   Double_t   thickDet[kNMaxDet];     // thickness divided by X
   Double_t   errorRPhi[kNMaxDet];    // error in bending direction
   Double_t   errorZ[kNMaxDet];       // error in z direction
   Double_t   errorR[kNMaxDet];       // error in r direction,from alignement only
   Int_t      iFlagDet[kNMaxDet];     // 1: sensitive detector 
                                      // 2: errors will be calculated
   Int_t      iFlagGas[kNMaxDet];     // for gas detectors

   Int_t     nDet;

   //dummy
   nDet            = 0;
   rDet[nDet]      = 0;
   thickDet[nDet]  = 0;
   iFlagDet[nDet]  = 0;
   iFlagGas[nDet]  = 0;
   //vacum pipe
   nDet            = 1;
   rDet[nDet]      = 3.0;
   thickDet[nDet]  = 0.06/35.3;       // berylium
   iFlagDet[nDet]  = 0;               // no detection
   iFlagGas[nDet]  = 0;
   //
   nDet            = 2;
   rDet[nDet]      = 3.5;
   thickDet[nDet]  = 1.0000/30420.0;  // air
   iFlagDet[nDet]  = 0;               // no detection
   iFlagGas[nDet]  = 1;
   //
   nDet            = 3;
   rDet[nDet]      = 4.0;
   thickDet[nDet]  = 0.06/9.36;       // silicon
   errorRPhi[nDet] = 0.0015 + 0.0005;
   errorZ[nDet]    = 0.009  + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 4;
   rDet[nDet]      = 5.75;
   thickDet[nDet]  = 3.5/30420.0;       // silicon
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   //
   nDet            = 5;
   rDet[nDet]      = 7.5;
   thickDet[nDet]  = 0.06/9.36;          // silicon
   errorRPhi[nDet] = 0.0015 + 0.0005;
   errorZ[nDet]    = 0.009  + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 6;
   rDet[nDet]      = 10.75;
   thickDet[nDet]  = 6.5/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // first silicon drift
   nDet            = 7;
   rDet[nDet]      = 14.0;
   thickDet[nDet]  = 0.06/9.36;          // silicon
   errorRPhi[nDet] = 0.0025 + 0.0005;
   errorZ[nDet]    = 0.0025  + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 8;
   rDet[nDet]      = 19.0;
   thickDet[nDet]  = 10.0/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // second silicon drift
   nDet            =  9;
   rDet[nDet]      = 24.0;
   thickDet[nDet]  = 0.06/9.36;          // silicon
   errorRPhi[nDet] = 0.0025 + 0.0005;
   errorZ[nDet]    = 0.0025  + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 10;
   rDet[nDet]      = 32.0;
   thickDet[nDet]  = 16.0/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // first silicon strips
   nDet            = 11;
   rDet[nDet]      = 40.0;
   thickDet[nDet]  = 0.06/9.36;          // silicon
   errorRPhi[nDet] = 0.003 + 0.0005;
   errorZ[nDet]    = 0.100 + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 12;
   rDet[nDet]      = 42.5;
   thickDet[nDet]  = 5.0/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // second silicon strips
   nDet            = 13;
   rDet[nDet]      = 45.0;
   thickDet[nDet]  = 0.06/9.36;          // silicon
   errorRPhi[nDet] = 0.003 + 0.0005;
   errorZ[nDet]    = 0.100 + 0.0005;
   errorR[nDet]    = 0.001;
   iFlagDet[nDet]  = 1;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 14;
   rDet[nDet]      = 47.5;
   thickDet[nDet]  = 5.0/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   //
   nDet            = 15;
   rDet[nDet]      = 50.0;
   thickDet[nDet]  = 0.01;              // 1% of something ITS
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 16;
   rDet[nDet]      = 51.0;
   thickDet[nDet]  = 2.0/30420.0;        // air
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // TPC HV degrager
   nDet            = 17;
   rDet[nDet]      = 52.0;
   thickDet[nDet]  = 0.0018;           // 0.18 % of something TPC
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 18;
   rDet[nDet]      = 68.75;
   thickDet[nDet]  = 12.5/18310.0;          // CO2
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   //
   nDet            = 19;
   rDet[nDet]      = 71.25;
   thickDet[nDet]  = 12.5/18310.0;          // CO2
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // TPC inner field cage
   nDet            = 20;
   rDet[nDet]      = 78.0;
   thickDet[nDet]  = 0.0041;                // 0.41 % of something
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 0;
   //
   nDet            = 21;
   rDet[nDet]      = 83.5;
   thickDet[nDet]  = 11.0/32155.6;          // neon
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   //
   nDet            = 22;
   rDet[nDet]      = 94.5;
   thickDet[nDet]  = 11.0/32155.6;          // neon
   iFlagDet[nDet]  = 0;              
   iFlagGas[nDet]  = 1;
   // TPC
   Int_t    nPadRow = 75;
   Double_t rCurrent = 99.0;
   Double_t deltaR   =  2.0;
   for(Int_t ipad=1; ipad<nPadRow+1; ipad++){
      nDet=nDet+1;
      rCurrent = rCurrent + deltaR;
      rDet[nDet] = rCurrent;
      thickDet[nDet]  = 2.0/32155.6;        // neon
      errorRPhi[nDet] = 0.0;                //errors are momentum dependent
      errorZ[nDet]    = 0.0;                //to be calculated latter
      errorR[nDet]    = 0.0075;
      iFlagDet[nDet]  = 2;                  //means error defined latter              
      iFlagGas[nDet]  = 1;
   }

   // vertex precision
   Double_t multDensity = 3906.25;

   fErrorVertexX = 0.010/TMath::Sqrt(multDensity) + 0.00060;
   fErrorVertexY = fErrorVertexX;
   fErrorVertexZ = 0.025/TMath::Sqrt(multDensity) + 0.00075;


   // magnetic field
   fBMag = 2.0;
   fConstMag = 1.0/(fBMag*0.297792458e-3);


   // prepare more suitables variables

   Int_t nDetActive = 0;
     
   for(Int_t idDet=0; idDet<nDet+1; idDet++){
      fRDet[idDet]     = rDet[idDet]; 
      fRDetSQ[idDet]   = rDet[idDet]*rDet[idDet];
      fThickDet[idDet] = 0.0136* TMath::Sqrt(thickDet[idDet]);
      fIFlagDet[idDet]    = iFlagDet[idDet];
      fIFlagGas[idDet]    = iFlagGas[idDet];
      if(iFlagDet[idDet] > 0){
         nDetActive = nDetActive+1;
         fErrorR[idDet] = errorR[idDet]*errorR[idDet];
         if(iFlagDet[idDet] == 1){
            fErrorRPhi[idDet] = errorRPhi[idDet]*errorRPhi[idDet];
            fErrorZ[idDet] = errorZ[idDet]*errorZ[idDet];
	 }
      }
   }

   fErrorVertexX = fErrorVertexX*fErrorVertexX;
   fErrorVertexY = fErrorVertexY*fErrorVertexY;
   fErrorVertexZ = fErrorVertexZ*fErrorVertexZ;

   fNDetActive   = nDetActive;
   fNDet         = nDet;

  
}

//_____________________________________________________________________________
void AliFDet::PrintDetInfo()
{
  //to print information for the initialisation of the detector 
   printf("**************************************************************\n");
   printf("*                                                            *\n");
   printf("*                        ALICE detector                      *\n");
   printf("*                                                            *\n");
   printf("**************************************************************\n");

   for(Int_t idDet=0; idDet<fNDet+1; idDet++){
     if(fIFlagDet[idDet] == 0){
       printf("%5s %3d %8.1f %2s %10.5f %20s\n",
              "det=",idDet,fRDet[idDet],"cm",
              TMath::Power(fThickDet[idDet]/0.0136,2),
              "of X0 <---- pasive material");
     } else{
       printf("%5s %3d %8.1f %2s %10.5f %6s  %6.4f %6.4f \n",
              "det=",idDet,fRDet[idDet],"cm",
              TMath::Power(fThickDet[idDet]/0.0136,2),"of X0, errors",
              TMath::Sqrt(fErrorRPhi[idDet]),TMath::Sqrt(fErrorZ[idDet]));
     } 
   }
   printf("%20s %10.4f %10.4f %10.4f\n","vertex precision(x,y,z)",
          TMath::Sqrt(fErrorVertexX),
          TMath::Sqrt(fErrorVertexY),
          TMath::Sqrt(fErrorVertexZ));
   printf("%20s %10.4f %10s %8.1f %5s\n","magnetic field (kGauss)",fBMag,
          "(constant",fConstMag,")");

}
