#ifndef ALIITSRAD_H
#define ALIITSRAD_H

#include <TObject.h>
#include <TMatrixFfwd.h>


class TObjArray;

//                  ITS Class to calculate the radiation lenght matrix
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
//


class AliITSRad : public TObject { 

public:
    AliITSRad(); // Default constructor.
    AliITSRad(Int_t iimax, Int_t jjmax);         // class constructor
    virtual ~AliITSRad();                        // class destructor
  
  Int_t Getimax() {return imax;}               // return the first dimension of the matrices
  Int_t Getjmax() {return jmax;}               // return the second dimension of the matrices
  
  TMatrixF &GetRadMatrix1() {return *fmrad1;}   // return the radiation lengh matrix for layer 1
  TMatrixF &GetRadMatrix2() {return *fmrad2;}   // return the radiation lengh matrix for layer 2
  TMatrixF &GetRadMatrix3() {return *fmrad3;}   // return the radiation lengh matrix for layer 3
  TMatrixF &GetRadMatrix4() {return *fmrad4;}   // return the radiation lengh matrix for layer 4
  TMatrixF &GetRadMatrix5() {return *fmrad5;}   // return the radiation lengh matrix for layer 5
  TMatrixF &GetRadMatrix6() {return *fmrad6;}   // return the radiation lengh matrix for layer 6
  
private:

  AliITSRad(const AliITSRad &source); // copy constructor
  // assignment operator
  AliITSRad& operator=(const AliITSRad &source);

  Int_t           imax;        // first dimension of the matrices
  Int_t           jmax;        // second dimension of the matrices
  
  TMatrixF         *fmrad1;     // matrix of the radiation lenghts for layer 1
  TMatrixF         *fmrad2;     // matrix of the radiation lenghts for layer 2
  TMatrixF         *fmrad3;     // matrix of the radiation lenghts for layer 3
  TMatrixF         *fmrad4;     // matrix of the radiation lenghts for layer 4
  TMatrixF         *fmrad5;     // matrix of the radiation lenghts for layer 5
  TMatrixF         *fmrad6;     // matrix of the radiation lenghts for layer 6

  ClassDef(AliITSRad, 1)
};

#endif

