///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoPicoEventCollectionVectorHideAway: a helper class for         //
// managing many mixing buffers with up to three variables used for      //
// binning.                                                              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include "AliFemtoPicoEventCollectionVectorHideAway.h"

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
  // basic constructor
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
  // return mixing event collection from a given bin
  if ( ix<0 || ix >= fBinsx) return 0;
  if ( iy<0 || iy >= fBinsy) return 0;
  if ( iz<0 || iz >= fBinsz) return 0;
  int bin = ix + iy*fBinsx + iz*fBinsy*fBinsx; 
//   cout << " AliFemtoPicoEventCollectionVectorHideAway::PicoEventCollection(...) - bin(ix,iy,iz): ";
//   cout << bin << "(" << ix <<"," << iy << "," << iz <<")" << endl;
  return fCollectionVector[bin]; 
}
// -----------------------------------
AliFemtoPicoEventCollection* AliFemtoPicoEventCollectionVectorHideAway::PicoEventCollection(double x, double y, double z) {
  // return mixing event collection for given values on x, y, z axes
  int ix,iy,iz;
  ix=0;iy=0;iz=0;

  if(fStepx != 0 && fStepy != 0 && fStepz != 0)
    {
      ix = (int)floor( (x-fMinx)/fStepx );
      iy = (int)floor( (y-fMiny)/fStepy );
      iz = (int)floor( (z-fMinz)/fStepz );
    }
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
  // copy constructor
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
  // destructor
  fCollectionVector.clear();
}
//___________________________________
AliFemtoPicoEventCollectionVectorHideAway& AliFemtoPicoEventCollectionVectorHideAway::operator=(const AliFemtoPicoEventCollectionVectorHideAway& aColl)
{
  // assignment operator
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
unsigned int AliFemtoPicoEventCollectionVectorHideAway::GetBinXNumber(double x) { return (int)floor( (x-fMinx)/fStepx ); }
unsigned int AliFemtoPicoEventCollectionVectorHideAway::GetBinYNumber(double y) { return (int)floor( (y-fMiny)/fStepy ); }
unsigned int AliFemtoPicoEventCollectionVectorHideAway::GetBinZNumber(double z) { return (int)floor( (z-fMinz)/fStepz ); }
