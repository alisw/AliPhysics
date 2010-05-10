#ifndef ALITRDPIDREFMAKER_H
#define ALITRDPIDREFMAKER_H

////////////////////////////////////////////////////////////
//
// Base class for the Task to build TRD-PID reference data
// For the actual implementation please check the classes
//   - AliTRDpidRefMakerNN (Neural Networks)
//   - AliTRDpidRefMakerLQ (Multidimensional Likelihood)
//
// Authors: Alex Bercuci <A.Bercuci@gsi.de>
//          Alex Wilk    <wilka@uni-muenster.de>
//
/////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif
#ifndef ALIPID_H
#include "AliPID.h"
#endif
#ifndef ALITRDCALPID_H
#include "Cal/AliTRDCalPID.h"
#endif
#ifndef ALITRDGEOMETRY_H
#include "AliTRDgeometry.h"
#endif

class TTree;
class TObjArray;
class AliTRDReconstructor;
class AliTRDseedV1;
class AliTRDtrackInfo;
class AliTRDpidRefMaker : public AliTRDrecoTask
{
public:
  enum ETRDpidRefMakerPBins {
    k006  =  0
    ,k008 =  1
    ,k010 =  2
    ,k015 =  3
    ,k020 =  4
    ,k030 =  5
    ,k040 =  6
    ,k050 =  7
    ,k060 =  8
    ,k080 =  9
    ,k100 = 10
    ,kAll = 11
  };
  enum ETRDpidRefMakerSource {
    kV0 = 0 // use V0 as reference
   ,kMC = 1 // use MC as reference
   ,kRec= 2 // use Reconstructed PID as reference
  };

  struct AliTRDpidRefData {
    AliTRDpidRefData() : fPLbin(0xff) {
      memset(fdEdx, 0, 8*sizeof(Float_t));
    }
    virtual ~AliTRDpidRefData(){}
    UChar_t fPLbin;   // momentum / layer bin
    Float_t fdEdx[8]; // dEdx array
    ClassDef(AliTRDpidRefData, 1)  // PID data representation
  };
  struct AliTRDpidRefDataArray {
    AliTRDpidRefDataArray();
    virtual ~AliTRDpidRefDataArray();
    void    PushBack(Int_t ly, Int_t p, Float_t *dedx);
    void    Reset();

    Int_t   fNtracklets;     // number of tracklets
    AliTRDpidRefData *fData; //[fNtracklets] PID data array
  private:
    AliTRDpidRefDataArray(const AliTRDpidRefMaker::AliTRDpidRefDataArray& ref);
    AliTRDpidRefDataArray& operator=(const AliTRDpidRefMaker::AliTRDpidRefDataArray& ref);
    ClassDef(AliTRDpidRefDataArray, 1)  // track PID data representation
  };

  AliTRDpidRefMaker();
  AliTRDpidRefMaker(const char *name, const char *title);

  virtual ~AliTRDpidRefMaker();
  
  void    UserCreateOutputObjects();
  void    UserExec(Option_t *option);
  Float_t GetPthreshold() const { return fPthreshold;}

  void    SetAbundance(Float_t train);
  void    SetPthreshold(Float_t t) { fPthreshold = t;}
  void    SetRefPID(ETRDpidRefMakerSource select, AliTRDtrackInfo *t, Float_t *pid);
  void    SetSource(ETRDpidRefMakerSource pid, ETRDpidRefMakerSource momentum) {fRefPID = pid; fRefP = momentum;}


protected:
  virtual Bool_t   CheckQuality(AliTRDseedV1* trklt);
  virtual Float_t* CookdEdx(AliTRDseedV1* trklt);
  virtual void     LinkPIDdata();
  virtual void     Fill();

  AliTRDReconstructor *fReconstructor;  // reconstructor needed for recalculation the PID
  TObjArray     *fV0s;                  //! v0 array
  TTree         *fData;                 //! dEdx-P data
  TObjArray     *fInfo;                 //! list of PID info
  AliTRDpidRefDataArray *fPIDdataArray; //! pid data array
  ETRDpidRefMakerSource  fRefPID;       // reference PID source
  ETRDpidRefMakerSource  fRefP;         // reference momentum source
  UChar_t       fPIDbin;                // species bin
  Float_t       fFreq;                  // training sample relative abundance
  Float_t       fP;                     // momentum
  Float_t       fdEdx[8];               // dEdx array
  Float_t       fPID[AliPID::kSPECIES]; // pid from v0s

private:
  AliTRDpidRefMaker(const AliTRDpidRefMaker&);              // not implemented
  AliTRDpidRefMaker& operator=(const AliTRDpidRefMaker&);   // not implemented

  Float_t        fPthreshold;            // momentum threshold [GeV/c]

  ClassDef(AliTRDpidRefMaker, 3); // TRD PID reference  maker base class
};

#endif
