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
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
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
										     int bz, double lz, double uz):
  fBinsTot(0),
  fBinsx(bx), fBinsy(by), fBinsz(bz),
  fMinx(lx),  fMiny(ly),  fMinz(lz),
  fMaxx(ux),  fMaxy(uy),  fMaxz(uz),
  fStepx(0),  fStepy(0),  fStepz(0),
  fCollection(0),
  fCollectionVector(0)
{
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
//___________________________________
AliFemtoPicoEventCollectionVectorHideAway::AliFemtoPicoEventCollectionVectorHideAway(const AliFemtoPicoEventCollectionVectorHideAway& aColl):
  fBinsTot(0),
  fBinsx(0), fBinsy(0), fBinsz(0),
  fMinx(0),  fMiny(0),  fMinz(0),
  fMaxx(0),  fMaxy(0),  fMaxz(0),
  fStepx(0),  fStepy(0),  fStepz(0),
  fCollection(0),
  fCollectionVector(0)
{

  fBinsTot = aColl.fBinsTot;
  fBinsx = aColl.fBinsx;
  fBinsy = aColl.fBinsy;
  fBinsz = aColl.fBinsz;
  fMinx  = aColl.fMinx;
  fMiny  = aColl.fMiny;
  fMinz  = aColl.fMinz;
  fMaxx  = aColl.fMaxx;
  fMaxy  = aColl.fMaxy;
  fMaxz  = aColl.fMaxz;
  fStepx = aColl.fStepx;
  fStepy = aColl.fStepy;
  fStepz = aColl.fStepz;
  fCollection = aColl.fCollection;

  fCollectionVector.clear();
  for (int iter=0; aColl.fCollectionVector.size();iter++){
    fCollectionVector.push_back(aColl.fCollectionVector[iter]);
  }
}
//___________________________________
AliFemtoPicoEventCollectionVectorHideAway::~AliFemtoPicoEventCollectionVectorHideAway()
{
  fCollectionVector.clear();
}
//___________________________________
AliFemtoPicoEventCollectionVectorHideAway& AliFemtoPicoEventCollectionVectorHideAway::operator=(const AliFemtoPicoEventCollectionVectorHideAway& aColl)
{
  if (this == &aColl)
    return *this;

  fBinsTot = aColl.fBinsTot;
  fBinsx = aColl.fBinsx;
  fBinsy = aColl.fBinsy;
  fBinsz = aColl.fBinsz;
  fMinx  = aColl.fMinx;
  fMiny  = aColl.fMiny;
  fMinz  = aColl.fMinz;
  fMaxx  = aColl.fMaxx;
  fMaxy  = aColl.fMaxy;
  fMaxz  = aColl.fMaxz;
  fStepx = aColl.fStepx;
  fStepy = aColl.fStepy;
  fStepz = aColl.fStepz;
  fCollection = aColl.fCollection;

  fCollectionVector.clear();

  for (int iter=0; aColl.fCollectionVector.size();iter++){
    fCollectionVector.push_back(aColl.fCollectionVector[iter]);
  }

  return *this;
}
