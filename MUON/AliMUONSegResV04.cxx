/////////////////////////////////////////////////////
//  Segmentation and Response classes version 04   //
/////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliMUONSegResV04.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"

//___________________________________________
ClassImp(AliMUONsegmentationV04)


void AliMUONsegmentationV04::Init(AliMUONchamber* )
{
    //printf("\n Initialise segmentation v04 \n");
//
//  Fill the arrays fCx (x-contour) and fNpxS (ix-contour) for each sector
//  These arrays help in converting from real to pad co-ordinates and
//  vice versa
//
//
//  PCB module size in cm
    const Float_t dxPCB=40, dyPCB=40;
//  PCB distribution (7 rows with 1+3 segmentation regions)
    const Int_t kpcb[7][4] = {{1, 2, 2, 2}, 
			      {0, 3, 2, 2}, 
			      {0, 2, 2, 2}, 
			      {0, 0, 3, 3}, 
			      {0, 0, 2, 3}, 
			      {0, 0, 0, 4}, 
			      {0, 0, 0, 3}};
    
    
//
//                             3 3 3 | 3 3 3
//                           3 3 3 3 | 3 3 3 3
//                         3 3 3 2 2 | 2 2 3 3 3
//                       3 3 3 2 2 2 | 2 2 2 3 3 3
//                       3 3 2 2 1 1 | 1 1 2 2 3 3      
//                     3 3 2 2 1 1 1 | 1 1 1 2 2 3 3
//                     3 3 2 2 1 1 0 | 0 1 1 2 2 3 3
//                    ------------------------------
//                     3 3 2 2 1 1 0 | 0 1 1 2 2 3 3
//                     3 3 2 2 1 1 1 | 1 1 1 2 2 3 3
//                       3 3 2 2 1 1 | 1 1 2 2 3 3      
//                       3 3 3 2 2 2 | 2 2 2 3 3 3                      
//                         3 3 3 2 2 | 2 2 3 3 3
//                           3 3 3 3 | 3 3 3 3
//                             3 3 3 | 3 3 3
//
// number of pad rows per PCB
//    
    Int_t NpyPCB=Int_t(dyPCB/fDpy);
//
// maximum number of pad rows    
    fNpy=7*NpyPCB;
//
//  Calculate padsize along x
    fDpxD[fNsec-1]=fDpx;
    if (fNsec > 1) {
	for (Int_t i=fNsec-2; i>=0; i--){
	    fDpxD[i]=fDpxD[fNsec-1]/fNDiv[i];
	    //printf("\n test ---dx %d %f \n",i,fDpxD[i]);
	}
    }
//
// fill the arrays defining the pad segmentation boundaries
//
//  loop over pcb module rows
    Int_t iy=0;
    for (Int_t irow=0; irow<7; irow++) {
//  
//  loop over pads along the anode wires
	for (Int_t i=0; i<=NpyPCB; i++) {
//  iy counts the padrow
	    iy++;
//  Loop over sectors (isec=0 is the dead space surounding the beam pipe)
	    for (Int_t isec=0; isec<4; isec++) {
		if (isec==0) {
		    fNpxS[0][iy]=kpcb[irow][0]*Int_t(dxPCB/fDpxD[0]);
		    fCx[0][iy]=kpcb[irow][0]*dxPCB;
		} else {
		    fNpxS[isec][iy]=fNpxS[isec-1][iy]
			+kpcb[irow][isec]*Int_t(dxPCB/fDpxD[isec]);

		    fCx[isec][iy]=fCx[isec-1][iy]
		    +kpcb[irow][isec]*dxPCB;
		}
	    } // sectors
	} // pad raws in module
    } // PCB rows
/*
    for (Int_t iy=1; iy< fNpy; iy++) {
	    printf("\nBoundary %d %f %d %f %d %f %d %f",
		   fNpxS[0][iy], fCx[0][iy],
		   fNpxS[1][iy], fCx[1][iy],
		   fNpxS[2][iy], fCx[2][iy],
		   fNpxS[3][iy], fCx[3][iy]);
	    
    }
*/
}

void AliMUONsegmentationV04::GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)
{
    n=3;
    x[0]=(fCx[1][1]+fCx[0][1])/2/TMath::Sqrt(2.);
    y[0]=x[0];
    x[1]=(fCx[2][1]+fCx[1][1])/2/TMath::Sqrt(2.);
    y[1]=x[1];
    x[2]=(fCx[3][1]+fCx[2][1])/2/TMath::Sqrt(2.);
    y[2]=x[2];
    printf("\n 1%f %f", x[0], y[0]);
    printf("\n 2%f %f", x[1], y[1]);
    printf("\n 3%f %f\n ", x[2], y[2]);
}
















