#ifndef ALIITSRAD_H
#define ALIITSRAD_H

#include <TObject.h>
#include <TMatrix.h>


class TObjArray;
class TMatrix;


//                  ITS Class to calculate the radiation lenght matrix
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
//


class AliITSRad : public TObject { 

public:

  AliITSRad(Int_t iimax, Int_t jjmax);         // class constructor
  ~AliITSRad();                                // class destructor
  
  Int_t Getimax() {return imax;}               // return the first dimension of the matrices
  Int_t Getjmax() {return jmax;}               // return the second dimension of the matrices
  
  TMatrix &GetRadMatrix1() {return *fmrad1;}   // return the radiation lengh matrix for layer 1
  TMatrix &GetRadMatrix2() {return *fmrad2;}   // return the radiation lengh matrix for layer 2
  TMatrix &GetRadMatrix3() {return *fmrad3;}   // return the radiation lengh matrix for layer 3
  TMatrix &GetRadMatrix4() {return *fmrad4;}   // return the radiation lengh matrix for layer 4
  TMatrix &GetRadMatrix5() {return *fmrad5;}   // return the radiation lengh matrix for layer 5
  TMatrix &GetRadMatrix6() {return *fmrad6;}   // return the radiation lengh matrix for layer 6
  
private:

  Int_t           imax;        // first dimension of the matrices
  Int_t           jmax;        // second dimension of the matrices
  
  TMatrix         *fmrad1;     // matrix of the radiation lenghts for layer 1
  TMatrix         *fmrad2;     // matrix of the radiation lenghts for layer 2
  TMatrix         *fmrad3;     // matrix of the radiation lenghts for layer 3
  TMatrix         *fmrad4;     // matrix of the radiation lenghts for layer 4
  TMatrix         *fmrad5;     // matrix of the radiation lenghts for layer 5
  TMatrix         *fmrad6;     // matrix of the radiation lenghts for layer 6

  ClassDef(AliITSRad, 1)
};

#endif

