#ifndef AliGPCone_H
#define AliGPCone_H

#include "AliGTube.h"
const Int_t kDiv = 20;               //default number of divisions

class AliGPCone: public AliGShape {

    private:

      Double_t* fSiTab;  //Table of sin(fPhi1) .... sin(fPhil+fDphi1)
      Double_t* fCoTab;  //Table of cos(fPhi1) .... cos(fPhil+fDphi1)

    protected:

       Float_t  fPhi1;   //lower phi limit
       Float_t  fDphi1;  //range in phi
       Int_t    fNz;       //number of z segments
       Float_t* fRmin;  // pointer to array of inside radiuses
       Float_t* fRmax;  // pointer to array of outside radiuses
       Float_t* fDz;    // pointer to array of half lengths in z
       Int_t    fNdiv;     //number of divisions

    public:
        AliGPCone(); /* Default Constructor */
        AliGPCone( Text_t *name, Text_t* title,  Float_t phi1, Float_t dphi1, Int_t nz);
	AliGPCone( AliGPCone* pcone );
	AliGPCone( Text_t *name, Text_t* title,  Float_t phi1, Float_t dphi1,
	Int_t nz, Float_t Z[10], Float_t RMIN[10], Float_t RMAX[10]);
        AliGPCone( Text_t *name, Text_t* title,  Float_t *upar, Int_t np);
        virtual ~AliGPCone(); /* Destructor */

        //Float_t GetRmin2() {return fRmin2;}
        //Float_t GetRmax2() {return fRmax2;}
	        void  DefineSection(Int_t secNum, Float_t z, Float_t rmin, Float_t rmax);
        virtual void  Draw(Option_t *option);
        virtual void  DrawShape(Option_t *option); // *MENU*
        virtual Int_t GetNumberOfDivisions () const {if (fNdiv) return fNdiv; else return kDiv;}
	        void  MakeTableOfCoSin();
                void  Paint(Option_t *option);
	        void  SetNumberOfDivisions (Int_t p);
        virtual void  SetPoints(Float_t *buff);


    ClassDef(AliGPCone,1) // Polycone class
};
#endif
