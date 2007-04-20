#ifndef THYDJET_H
#define THYDJET_H


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THydjet                                                              //
//                                                                      //
// This class implements an interface to the Hydjet event generator.    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGenerator
#include "TGenerator.h"
#endif
class TObjArray;

class THydjet : public TGenerator {


public:

   THydjet();
   THydjet(Float_t efrm, const char *frame, Float_t aw,
           Int_t ifb, Float_t bmin, Float_t bmax, Float_t bfix, Int_t nh);
   virtual            ~THydjet();

   virtual void        Initialize();

   virtual void        GenerateEvent();

   virtual Int_t       ImportParticles(TClonesArray *particles, Option_t *option="");
   virtual TObjArray*  ImportParticles(Option_t *option="");


   //Parameters for the generation:
   
   virtual void        SetEfrm(Float_t efrm);
   virtual Float_t     GetEfrm() const;

   virtual void        SetFrame(const char *frame);
   virtual const char *GetFrame() const;

   virtual void        SetAw(Float_t aw);
   virtual Float_t     GetAw() const;

   virtual void        SetIfb(Int_t ifb);
   virtual Int_t       GetIfb() const;

   virtual void        SetBmin(Float_t bmin);
   virtual Float_t     GetBmin() const;

   virtual void        SetBmax(Float_t bmax);
   virtual Float_t     GetBmax() const;

   virtual void        SetBfix(Float_t bfix);
   virtual Float_t     GetBfix() const;

   virtual void        SetNh(Int_t nh);
   virtual Int_t       GetNh() const;



   //common HYFLOW access routines:

   virtual const void  SetYTFL(Float_t ytfl) const;
   virtual Float_t     GetYTFL() const;

   virtual const void  SetYLFL(Float_t ylfl) const;
   virtual Float_t     GetYLFL() const;

   virtual const void  SetFPART(Float_t fpart) const;
   virtual Float_t     GetFPART() const;

   //common HYJPAR access routines

   virtual const void  SetNHSEL(Int_t nhsel) const;
   virtual Int_t       GetNHSEL() const;

   virtual const void  SetPTMIN(Float_t ptmin) const;
   virtual Float_t     GetPTMIN() const;

   virtual const void  SetNJET(Int_t njet) const;
   virtual Int_t       GetNJET() const;

   // common HYFPAR access routines - read-only common:

   virtual Float_t     GetBGEN() const;

   virtual Int_t     GetNBCOL() const;

   virtual Int_t     GetNPART() const;

   virtual Int_t     GetNPYT() const;

   virtual Int_t     GetNHYD() const;


   // common LUJETS access routines - read-only common:

   virtual Int_t       GetN() const;

   virtual Int_t       GetK(Int_t key1, Int_t key2) const;

   virtual Float_t     GetP(Int_t key1, Int_t key2) const;

   virtual Float_t     GetV(Int_t key1, Int_t key2) const;


   // common HYJETS access routines - read-only common:

   virtual Int_t       GetNL() const;

   virtual Int_t       GetKL(Int_t key1, Int_t key2) const;

   virtual Float_t     GetPL(Int_t key1, Int_t key2) const;

   virtual Float_t     GetVL(Int_t key1, Int_t key2) const;

   // common PYDAT1 access routines:

   virtual void        SetMSTU(Int_t key, Int_t value);

   virtual void        SetPARU(Int_t key, Double_t value);

   virtual void        SetMSTJ(Int_t key, Int_t value);

   virtual void        SetPARJ(Int_t key, Double_t value);

   // common PYSUBS access routines:

   virtual const void  SetMSEL(Int_t msel) const;

   virtual void        SetCKIN(Int_t key, Double_t value);

   // common PYPARS access routines:

   virtual void        SetMSTP(Int_t key, Int_t value);

   // access to HYDJET routines:

   virtual void         Hydro();


   protected:

    Float_t      fEfrm;     // Energy in the centre of mass (CMS) or lab-frame (LAB)
    TString      fFrame;    // Reference frame CMS or LAB
    Float_t      fAw;       // Beam and target nucleus atomic weight
    Int_t        fIfb;       // flag of type of centrality generation
                            //    0 impact parameter fixed (bfix)
                            //    else impact parameter is generated with standard Glauber geometry
	                         //    between minimum (bmin) and maximum (bmax) values
    Float_t      fBmin;     // Minimum impact parameter in units of nucleus radius RA
                            // (i.e. minimum value in [fm] will be bmin*RA),
                            // valid only if ifb not equal to zero
    Float_t      fBmax;     // Maximum impact parameter in units of nucleus radius RA
                            // (i.e. maximum value in [fm] will be bmax*RA),
                            // valid only if ifb not equal to zero
    Float_t      fBfix;     // Fixed impact parameter in units of nucleus radius RA
                            // (i.e. fixed value in [fm] will be bfix*RA),
                            // valid only if ifb=0
    Int_t        fNh;       // Mean soft hadron multiplicity in central Pb-Pb collisions
                            // (multiplicity for other centralities and atomic numbers
                            // will be calculated automatically).

    ClassDef(THydjet,1)  //Interface to Hydjet Event Generator
};

#endif







