#ifndef ALIGTUBE_H
#define ALIGTUBE_H

#include "AliGShape.h"

const Int_t kDivNum = 20;    //default number of divisions

class AliGTube: public AliGShape {

    protected:
        Double_t* fCoTab;       // Table of cos(fPhi1) .... cos(fPhil+fDphi1)
        Float_t   fAspectRatio; // defines  (the ellipse semi-axis in Y)/(the ellipse semi-axis in X)
        Float_t   fDz;          // half length in z
        Int_t     fNdiv;        // number of segments (precision)
        Float_t   fRmax;        // ellipse  semi-axis   in  X outside
        Float_t   fRmin;        // ellipse  semi-axis   in  X inside
        Double_t* fSiTab;       // Table of sin(fPhi1) .... sin(fPhil+fDphi1)
        
        virtual void MakeTableOfCoSin();  // Create the table of the fSiTab; fCoTab

    public:
        AliGTube(); /* Default Constructor */
        AliGTube( Text_t *name, Text_t *title, Float_t rmin, Float_t rmax, Float_t dz, Float_t aspect=1); /* Constructor*/
        AliGTube( Text_t *name, Text_t *title, Float_t rmax, Float_t dz); /* Constructor */
        AliGTube( AliGTube *tube ); 
	~AliGTube(); /* Destructor */

                void    Draw(Option_t *option);
                void    DrawShape(Option_t *option); // *MENU* 
                Float_t GetAspectRatio(){return fAspectRatio;}
                Float_t GetDz()    {return fDz;}
                Int_t   GetNdiv()  {return fNdiv;}
                Int_t   GetNumberOfDivisions () const {if (fNdiv) return fNdiv; else return kDivNum;}
	        Float_t GetRmax() {return fRmax;}
	        Float_t GetRmin() {return fRmin;}
                void    Paint(Option_t *option);
                void    PaintGLPoints(Float_t *vertex);
                void    SetDz(Float_t dz)    {fDz= dz;}
                void    SetPoints(Float_t *buff);
	        void    SetRmin(Float_t rmin) {fRmin= rmin;}
                void    SetRmax(Float_t rmax) {fRmax= rmax;}
                void    Sizeof3D() const;

    ClassDef(AliGTube,1) // Simple cone class
};

#endif
