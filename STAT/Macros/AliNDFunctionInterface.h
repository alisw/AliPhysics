#ifndef AliNDFunctionInterface_H
#define AliNDFunctionInterface_H

/// \ingroup STAT
/// \namespace AliNDFunctionInterface
/// \brief  Interface to Ndimension functional representations (THn and TMVA)
/// \authors Marian  Ivanov marian.ivanov@cern.ch
/// see example usage in

#include <map>
class TGraph;
class TH1;
class TTreeSRedirector;
class THn;

namespace AliNDFunctionInterface {
  /// generic variadic function - to get it from boost in the future
  template<typename T> vector<T> add_to_vector(vector<T> &z, T v);
  template<typename T, typename... Args> vector<T> add_to_vector(vector<T> &z, T v, Args... args);
  template<typename T, typename... Args> vector<T> make_vector(T v, Args... args);
  ///
  std::map<int, THn*> hnMapArrayInt;
  std::map<std::string, THn*> hnMapArrayName;
  // THn interpolation interface
  Double_t GetInterpolationLinear(THn * his, Double_t *xyz, Int_t  verbose);
  Double_t GetDeltaInterpolationLinear(THn * his, Double_t *xyz, Int_t dIndex, Int_t  verbose);
  // Global function interface
  Double_t GetInterpolationLinear(Int_t index, Double_t *xyz, Int_t  verbose){return GetInterpolationLinear(hnMapArrayInt[index], xyz,verbose);}
  Double_t GetDeltaInterpolationLinear(Int_t index, Double_t *xyz, Int_t dIndex, Int_t  verbose) {
      return GetDeltaInterpolationLinear(hnMapArrayInt[index], xyz, dIndex, verbose);
  }
  Double_t GetInterpolationLinear(const char *name, Double_t *xyz, Int_t  verbose){return GetInterpolationLinear(hnMapArrayName[name], xyz,verbose);}
  template<typename T, typename... Args> T EvalTHnLinear(int id, T v, Args... args);        /// variadic function evaluating THn
  /// TMVA interface
  map<int, TMVA::MethodBase *> readerMethodBase;            /// map of registered TMVA::MethodBase
  map<std::string, std::string> regressionMethodSetting;    /// map of registered TMVA regression methods
  map<std::string, TMVA::Types::EMVA> regressionMethodID;   /// map of registered TMVA regression methods
  void registerDefaultMVAMethods();                         /// example registering default methods ()
  void registerMethod(std::string method,  std::string content, TMVA::Types::EMVA id){regressionMethodSetting[method]=content; regressionMethodID[method]=id;}
  Int_t  FitMVA(TTree *tree, const char *varFit, TCut cut, const char * variableList, const char *methodList);   /// MVA regression
  Int_t  LoadMVAReader(Int_t id, const char * inputFile, const char *method, const char *dir);
  template<typename T, typename... Args> T EvalMVA(int id, T v, Args... args);        /// variadic function evaluating MVA

};


/// Helper function to create std vector in variadic function
template<typename T> vector<T> AliNDFunctionInterface::add_to_vector(vector<T> &z, T v) {z.push_back(v); return z;}
template<typename T, typename... Args> vector<T> AliNDFunctionInterface::add_to_vector(vector<T> &z, T v, Args... args) {
    z.push_back(v);return add_to_vector<T>(z, args...);
}
/// Variadic function to create an vector (boost implementation )
/// Example usage: to create vector of integers
///     auto a = AliNDFunctionInterface::make_vector<int>(1, 1, 3);
/// \tparam T
/// \tparam Args
/// \param v
/// \param args
/// \return
/// \code
///    auto a = AliNDFunctionInterface::make_vector<int>(1, 1, 3);
/// \endcode
template<typename T, typename... Args> vector<T> AliNDFunctionInterface::make_vector(T v, Args... args) {
    vector<T> z; z.push_back(v); return add_to_vector<T>(z, args...);
}

/// Variadic function to interpolate THn
/// \tparam T
/// \tparam Args
/// \param id
/// \param v
/// \param args
/// \return

template<typename T, typename... Args> T AliNDFunctionInterface::EvalTHnLinear(int id, T v, Args... args){
  auto a = make_vector<double>(v, args...);
  THn *his = hnMapArrayInt[id];
  if (his!=NULL) return GetInterpolationLinear(his, a.data(), 0);
  return 0;
};


/// Variadic function to evaluate MVA regression method registered using method ID
/// \tparam T
/// \tparam Args
/// \param id
/// \param v
/// \param args
/// \return
/// Example usage:
///    * usage in TTreeFormula. e.g calculate and visualize second derivative of regression
///\code
///   MVAInput->Draw("EvalMVA(0,fraction,Z,0)-(EvalMVA(0,fraction,Z-1,0)+EvalMVA(0,fraction,Z+1,0))*0.5:Z","Z>2&&fraction>5","");
///\edncode
///
template<typename T, typename... Args> T AliNDFunctionInterface::EvalMVA(int id, T v, Args... args){
  auto a = make_vector<float>(v, args...);
  TMVA::Event  event = TMVA::Event(a,a.size());
  TMVA::MethodBase * method= readerMethodBase[id];
  if (method==NULL) return 0;                      /// some optional verbosity needed
  return method->GetRegressionValues(&event)[0];
};


#endif

