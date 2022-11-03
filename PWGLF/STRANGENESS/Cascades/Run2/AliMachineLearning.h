#ifndef AliMachineLearning_H
#define AliMachineLearning_H
#include <TNamed.h>

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Machines
// virtual class: other specific classes will derive from this
// but share the interface!
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliMachineLearning : public TNamed {
  
public:
  //Simple constructor
  AliMachineLearning();
  
  //TNamed-inspired constructor
  AliMachineLearning(const char * name, const char * title = "Machine Learning");
  
  //Simple destructor
  ~AliMachineLearning();
  
  void Clear(Option_t* = "") {}; //dummy
  
  virtual void LoadModel(TString lModelName) {};

  virtual double Predict(double *X, int K) { return 0.0; };
  
  //virtual void Print(Option_t *option="") {};
  
private:
  ClassDef(AliMachineLearning, 1);
};
#endif
