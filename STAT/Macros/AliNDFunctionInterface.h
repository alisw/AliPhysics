#ifndef AliNDFunctionInterface_H
#define AliNDFunctionInterface_H

/// \ingroup STAT
/// \namespace AliNDFunctionInterface
/// \brief  Interface to N-dimensional functional representations (THn and TMVA)
/// \authors Marian  Ivanov marian.ivanov@cern.ch
/// see example usage in

#include <map>
class TGraph;
class TH1;
class TTreeSRedirector;
class THn;
class TObjArray;

namespace AliNDFunctionInterface {
  Int_t fVerbose=0;            /// verbosity
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
  map<int, TObjArray* > readerMethodBaseArray;              /// map of registered array of TMVA::MethodBase - used to define TMVA statistics (Mean, Median, RMS, quantile)
  map<std::string, std::string> regressionMethodSetting;    /// map of registered TMVA regression methods
  map<std::string, std::string> FactorySetting;    /// map of registered TMVA regression methods
  map<std::string, TMVA::Types::EMVA> regressionMethodID;   /// map of registered TMVA regression methods
  void registerDefaultMVAMethods();                         /// example registering default methods ()
  void registerFactory(std::string factory, std::string content){FactorySetting[factory]=content;}
  void registerMethod(std::string method,  std::string content, TMVA::Types::EMVA id){regressionMethodSetting[method]=content; regressionMethodID[method]=id;}
  Int_t  FitMVAClassification(const char * output , const char *inputTrees, const char *cuts ,const char * variableList, const char *methodList, const char * factoryString=""); /// MVA Classification
  Int_t  FitMVARegression(const char * output, TTree *tree, const char *varFit, TCut cut, const char * variables, const char *methods,const char * factoryString=""); /// new MVA Regression
  TMVA::MethodBase * LoadMVAReader(Int_t id, const char * inputFile, const char *method, const char *dir);
  Int_t  LoadMVAReaderArray(Int_t id, const char * inputFile, const char *methodMask, const char *dirMask);
  Int_t  AppendMethodToArray(Int_t index, TMVA::MethodBase * method);         /// Append method into array of methods - used e.g for bootstrap statistics
  Double_t   EvalMVAStatArray(int id, int statType, vector<float> point);     /// Evaluate statistic
  template<typename T, typename... Args> T EvalMVA(int id, T v, Args... args);        /// variadic function evaluating MVA
  template<typename T, typename... Args> T EvalMVAClasification(int id, T v, Args... args);        /// variadic function evaluating MVA
  template<typename T, typename... Args> T EvalMVAStat(int id, int statType,  T v, Args... args);        /// variadic function evaluating MVA array stat
};


/// \brief Helper function to create std vector in variadic function
template<typename T> vector<T> AliNDFunctionInterface::add_to_vector(vector<T> &z, T v) {z.push_back(v); return z;}
template<typename T, typename... Args> vector<T> AliNDFunctionInterface::add_to_vector(vector<T> &z, T v, Args... args) {
  z.push_back(v);return add_to_vector<T>(z, args...);
}
/// \brief Variadic function to create an vector (boost implementation )
/// \tparam T        - template type
/// \tparam Args     - list of arguments
/// \param v
/// \param args
/// \return
/// Example usage: to create vector of integers to create vector with 3 integers
///==============================================
/// \code
///    auto a = AliNDFunctionInterface::make_vector<int>(1, 1, 3);
/// \endcode
template<typename T, typename... Args> vector<T> AliNDFunctionInterface::make_vector(T v, Args... args) {
  vector<T> z; z.push_back(v); return add_to_vector<T>(z, args...);
}

/// Variadic function to linearly interpolate THn
/// \tparam T          - template type
/// \tparam Args       - list of arguments
/// \param id          - id of the function (e.g using  hash of the name)
/// \param v
/// \param args
/// \return interpolated value
/// Example usage:
/// ==================
/// * interpolation of 3D histogram at point xyz
/// * THn has  to be registered in map  hnMapArrayInt[id] before
///\code
/// value = AliNDFunctionInterface::EvalTHnLinear(0,x,y,z);
///\endcode
template<typename T, typename... Args> T AliNDFunctionInterface::EvalTHnLinear(int id, T v, Args... args){
  auto a = make_vector<double>(v, args...);
  THn *his = hnMapArrayInt[id];
  if (his!=NULL) return GetInterpolationLinear(his, a.data(), 0);
  return 0;
};


/// \brief Variadic function to evaluate MVA regression method registered using method ID
/// \tparam T
/// \tparam Args
/// \param id
/// \param v
/// \param args
/// \return
///
/// Example usage:
///===================
///    * usage in TTreeFormula. e.g calculate and visualize second derivative of regression
///\code
///   MVAInput->Draw("AliNDFunctionInterface::EvalMVA(0,fraction,Z,0)-(AliNDFunctionInterface::EvalMVA(0,fraction,Z-1,0)+AliNDFunctionInterface::EvalMVA(0,fraction,Z+1,0))*0.5:Z","Z>2&&fraction>5","");
///\endcode
///
template<typename T, typename... Args> T AliNDFunctionInterface::EvalMVA(int id, T v, Args... args){
  auto a = make_vector<float>(v, args...);
  TMVA::Event  event = TMVA::Event(a,a.size());
  TMVA::MethodBase * method= readerMethodBase[id];
  if (method==NULL) return 0;                      /// some optional verbosity needed
  return method->GetRegressionValues(&event)[0];
};

template<typename T, typename... Args> T AliNDFunctionInterface::EvalMVAClasification(int id, T v, Args... args){
  auto a = make_vector<float>(v, args...);
  TMVA::Event  event = TMVA::Event(a,a.size());
  TMVA::MethodBase * method= readerMethodBase[id];
  if (method==NULL) return 0;                      /// some optional verbosity needed
  return method->GetMvaValue(&event);
};


/// Template variadic function - Evaluate statistic on top of array (readerMethodBaseArray;) of MVA methods
/// To use the method - array o TMVA methods should be registered before using LoadMVAReaderArray or AppendMethodToArray
/// \tparam T            - template type - be default float
/// \tparam Args         -
/// \param id            - id of the registered array to evaluate
/// \param statType      - type of statistic (0-mean, 1-median, 2-rms)
/// \param v             -
/// \param args
/// \return
///
/// Example usage:
///===================
///   * compare mean and median statistic of array 2
///\code
///   MVAInput->Draw("AliNDFunctionInterface::EvalMVAStat(2,0,interactionRate, bz0, qmaxQASum, qmaxQASumR):AliNDFunctionInterface::EvalMVAStat(2,1,interactionRate, bz0, qmaxQASum, qmaxQASumR)","run==QA.EVS.run","");
///\endcode
template<typename T, typename... Args> T AliNDFunctionInterface::EvalMVAStat(int id, int statType,  T v, Args... args) {        /// variadic function evaluating MVA array stat
  auto a = make_vector<float>(v, args...);
  return EvalMVAStatArray(id,statType,a);
};

#endif

