#ifndef AliBDT_H
#define AliBDT_H
#include <TNamed.h>
#include "AliMachineLearning.h"
#include "TFile.h"
#include "TH1.h"

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Boosted Decision Trees
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliBDT : public AliMachineLearning {
  
public:
  //Dummy Constructor
  AliBDT();
  
  //Standard Constructor
  AliBDT(const char * name, const char * title = "BDT");

  //Simple destructor
  ~AliBDT();
  
  //Interface to configure parameters of the machine
  void LoadModel(TString lModelName);

  double Predict(double* X, int K);
  
  //void Print(Option_t *option="");
  
private:

  //Histograms to store BDT parameters
  TH1I* fFt;  //Feature to avaliate
  TH1D* fSoS; //Split of corresponding Feature or Score on the leaf
  //-----------------------------------------------------------------------
  
  ClassDef(AliBDT, 1)
};
#endif
