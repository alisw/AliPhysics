/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.3  2002/11/01 20:45:53  magestro
 * Fixed bug in 3rd dimension of event collection vector, probably never encountered
 *
 * Revision 1.2  2001/11/11 18:34:13  laue
 * AliFemtoPicoEventCollectionVectorHideAway: updated for 3d grid
 * AliFemtoVertexMultAnalysis: new
 *
 * Revision 1.1  2000/07/16 21:44:11  laue
 * Collection and analysis for vertex dependent event mixing
 *
 *
 **************************************************************************/
#include "Infrastructure/AliFemtoPicoEventCollectionVectorHideAway.h"

// -----------------------------------
AliFemtoPicoEventCollectionVectorHideAway::AliFemtoPicoEventCollectionVectorHideAway(int bx, double lx, double ux,
									       int by, double ly, double uy,
									       int bz, double lz, double uz) {
    fBinsx = bx;  fMinx = lx; fMaxx = ux;
    fBinsy = by;  fMiny = ly; fMaxy = uy;
    fBinsz = bz;  fMinz = lz; fMaxz = uz;

    fBinsTot = fBinsx * fBinsy * fBinsz;
    fStepx=0;  fStepx = (fMaxx-fMinx)/fBinsx;
    fStepy=0;  fStepy = (fMaxy-fMiny)/fBinsy;
    fStepz=0;  fStepz = (fMaxz-fMinz)/fBinsz;
    

    //fCollectionVector = new AliFemtoPicoEventCollectionVector();
    fCollection = 0;
    for ( int i=0; i<fBinsTot; i++) {
	fCollection = new AliFemtoPicoEventCollection();
	fCollectionVector.push_back(fCollection);
    }
}
// -----------------------------------
AliFemtoPicoEventCollection* AliFemtoPicoEventCollectionVectorHideAway::PicoEventCollection(int ix, int iy, int iz) { 
  if ( ix<0 || ix >= fBinsx) return 0;
  if ( iy<0 || iy >= fBinsy) return 0;
  if ( iz<0 || iz >= fBinsz) return 0;
  int bin = ix + iy*fBinsx + iz*fBinsy*fBinsx; 
  cout << " AliFemtoPicoEventCollectionVectorHideAway::PicoEventCollection(...) - bin(ix,iy,iz): ";
  cout << bin << "(" << ix <<"," << iy << "," << iz <<")" << endl;
  return fCollectionVector[bin]; 
}
// -----------------------------------
AliFemtoPicoEventCollection* AliFemtoPicoEventCollectionVectorHideAway::PicoEventCollection(double x, double y, double z) {
  int ix,iy,iz;
  ix=0;iy=0;iz=0;

  ix = (int)floor( (x-fMinx)/fStepx );
  iy = (int)floor( (y-fMiny)/fStepy );
  iz = (int)floor( (z-fMinz)/fStepz );

  return PicoEventCollection( ix,iy,iz );
}

