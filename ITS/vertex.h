#ifndef VERTEX_H
#define VERTEX_H
struct Hit_tl {
   Int_t   track,n,lad,det;
   Float_t x,y,z,xl,yl,zl,xr,yr,zr;
   Float_t r,phi,eta,phir,etar;
};

void FindVertexs(Int_t evnt,Float_t frac,Float_t len);

void HitsToV(Hit_tl **spdi,Int_t &nspdi,Hit_tl **spdo,Int_t &nspdo,
	     Int_t ntracks,TTree *TH,TClonesArray *ITShits,
	     Float_t fraction,Float_t len);

Float_t vertex(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max);
Float_t vertexSlow(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max);
Float_t vertexEta(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max);

#endif
