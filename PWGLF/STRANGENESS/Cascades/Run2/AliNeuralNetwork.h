#ifndef AliNeuralNetwork_H
#define AliNeuralNetwork_H
#include <TNamed.h>
#include "AliMachineLearning.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Neural Network 
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliNeuralNetwork : public AliMachineLearning {
  
public:
  //Dummy Constructor
  AliNeuralNetwork();
  
  //Standard Constructor
  AliNeuralNetwork(const char * name, const char * title = "Neural Network");

  //Simple destructor
  ~AliNeuralNetwork();
  
  //Interface to configure parameters of the machine
  void LoadModel(TString lModelName);

  double Predict(double* X, int K);
  
  //void Print(Option_t *option="");
  
private:

  //Histograms to store sinaptic weights of the Neural Network
  TH2D* fW1;
  TH1D* fB1;

  TH2D* fW2;
  TH1D* fB2;

  TH2D* fW3;
  TH1D* fB3;

  TH2D* fW4;
  TH1D* fB4;

  TH2D* fW5;
  TH1D* fB5;
  //-----------------------------------------------------------------------
  
  ClassDef(AliNeuralNetwork, 1)
};
#endif
