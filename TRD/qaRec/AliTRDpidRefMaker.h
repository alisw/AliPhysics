#ifndef ALITRDPIDREFMAKER_H
#define ALITRDPIDREFMAKER_H

////////////////////////////////////////////////////////////
//
// Base class for the Task to build TRD-PID reference data
// For the actual implementation please check the classes
//   - AliTRDpidRefMakerNN (Neural Networks)
//   - AliTRDpidRefMakerLQ (Multidimensional Likelihood)
//
// Responsible : Alex Bercuci <A.Bercuci@gsi.de>
//
/////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif
#ifndef ALIPID_H
#include "AliPID.h"
#endif
#ifndef ALITRDCALPID_H
#include "../Cal/AliTRDCalPID.h"
#endif
#ifndef ALITRDGEOMETRY_H
#include "../AliTRDgeometry.h"
#endif

class TTree;
class TObjArray;
class TEventList;
class AliTRDReconstructor;
class AliTRDseedV1;
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
  AliTRDpidRefMaker(const char *name=0, const char *title=0);

  virtual ~AliTRDpidRefMaker();
  
  void    ConnectInputData(Option_t *opt);
  void    CreateOutputObjects();
  void    Exec(Option_t *option);

  void    LoadContainer(const Char_t *InFileCont);
  void    LoadFile(const Char_t *InFile);

  void    SetAbundance(Float_t train, Float_t test);
  void    SetRefPID(void *source, Float_t *pid);
  void    SetSource(ETRDpidRefMakerSource pid, ETRDpidRefMakerSource momentum) {fRefPID = pid; fRefP = momentum;}

  void    Terminate(Option_t *);

protected:
  virtual Float_t* GetdEdx(AliTRDseedV1*) = 0;
  virtual Int_t    GetNslices() = 0;
  virtual void     Fill();

  AliTRDReconstructor *fReconstructor;  //! reconstructor needed for recalculation the PID
  TObjArray     *fV0s;                  //! v0 array
  TTree         *fData;                 //! dEdx-P data
  TEventList *fTrain[AliTRDCalPID::kNMom][AliTRDgeometry::kNlayer];          // Training list for each momentum 
  TEventList *fTest[AliTRDCalPID::kNMom][AliTRDgeometry::kNlayer];           // Test list for each momentum 
  ETRDpidRefMakerSource fRefPID;     // reference PID source
  ETRDpidRefMakerSource fRefP;       // reference momentum source
  Float_t       fTrainFreq;             //! training sample relative abundance
  Float_t       fTestFreq;              //! testing sample relative abundance
  UChar_t       fLy;                    //! TRD layer
  Float_t       fP;                     //! momentum
  Float_t       fdEdx[10];              //! dEdx array
  Float_t       fPID[AliPID::kSPECIES]; //! pid from v0s

private:
  AliTRDpidRefMaker(const AliTRDpidRefMaker&);              // not implemented
  AliTRDpidRefMaker& operator=(const AliTRDpidRefMaker&);   // not implemented

  ClassDef(AliTRDpidRefMaker, 1); // TRD PID reference  maker base class
};

#endif
