#include <TMatrix.h>
#include <TObjArray.h>
#include <Riostream.h>
#include <Riostream.h>


#include "AliITSRad.h"


ClassImp(AliITSRad)
//______________________________________________________________________
AliITSRad::AliITSRad():
imax(0),
jmax(0),
fmrad1(0),
fmrad2(0),
fmrad3(0),
fmrad4(0),
fmrad5(0),
fmrad6(0) {
    // Default constructor

}
AliITSRad::AliITSRad(Int_t iimax, Int_t jjmax):
imax(0),
jmax(0),
fmrad1(0),
fmrad2(0),
fmrad3(0),
fmrad4(0),
fmrad5(0),
fmrad6(0){

  imax=iimax;
  jmax=jjmax;
  
  fmrad1 = new TMatrix(imax,jmax);
  fmrad2 = new TMatrix(imax,jmax);
  fmrad3 = new TMatrix(imax,jmax);
  fmrad4 = new TMatrix(imax,jmax);
  fmrad5 = new TMatrix(imax,jmax);
  fmrad6 = new TMatrix(imax,jmax);

  ifstream in("ITSlegov5.map");
  Int_t i,j;
  
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad1)(i,j);
	 }
  }
	 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad2)(i,j);
	 }
  } 
	 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad3)(i,j);
	 }
  } 
	 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad4)(i,j);
	 }
  } 
	 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad5)(i,j);
	 }
  }
  
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   in>>(*fmrad6)(i,j);
	 }
  }
	 	 
  in.close();

/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//                       Stampe provvisorie delle matrici su rad.out         
//	 
  ofstream out("rad.out");
  	 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad1)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  out<<"\n";
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad2)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  out<<"\n";
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad3)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  out<<"\n";
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad4)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  out<<"\n"; 
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad5)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  out<<"\n";
  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
	   out<<(*fmrad6)(i,j)<<" ";
	 }
	 out<<"\n";
  }
  
  out.close(); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////	 
*/
}

AliITSRad::~AliITSRad() {

  delete fmrad1;
  delete fmrad2;
  delete fmrad3;
  delete fmrad4;
  delete fmrad5;
  delete fmrad6;
  
}
