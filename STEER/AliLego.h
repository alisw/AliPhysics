#ifndef ALILego_H
#define ALILego_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TH2.h>

class AliLego : public TNamed  {

private:
   Float_t    fThetaMin;        //Minimum generation theta
   Float_t    fThetaMax;        //Maximum generation theta
   Float_t    fPhiMin;          //Minimum generation phi
   Float_t    fPhiMax;          //Maximum generation phi
   Float_t    fRadMin;          //Generation radius
   Float_t    fRadMax;          //Maximum tracking radius
   Float_t    fZMax;            //Maximum tracking Z
   Int_t      fNtheta;          //Number of bins in Theta
   Int_t      fNphi;            //Number of bins in Phi
   Int_t      fThetaBin;        //Current theta bin
   Int_t      fPhiBin;          //Current phi bin
   Float_t    fCurTheta;        //Current theta of track
   Float_t    fCurPhi;          //Current phi of track
   Float_t    fTotRadl;         //Total Radiation length
   Float_t    fTotAbso;         //Total absorption length
   Float_t    fTotGcm2;         //Total G/CM2 traversed
   TH2F      *fHistRadl;        //Radiation length map 
   TH2F      *fHistAbso;        //Interaction length map
   TH2F      *fHistGcm2;        //g/cm2 length map
   TH2F      *fHistReta;        //Radiation length map as a function of eta
   
public:
  AliLego();
  AliLego(const char *name, const char *title);
  virtual ~AliLego();
  virtual void  GenerateKinematics();
  virtual void  Init(Int_t ntheta,Float_t themin,Float_t themax, Int_t nphi, Float_t phimin,
                     Float_t phimax,Float_t rmin,Float_t rmax,Float_t zmax);
  Float_t       PropagateCylinder(Float_t *x, Float_t *v, Float_t r, Float_t z);
  virtual void  Run();
  virtual void  StepManager();
  
  ClassDef(AliLego,1) //Utility class to compute and draw Radiation Length Map

};

#endif

