#include <TGraph.h>
#include <Math/IFunction.h>

struct WrappedGraph : public ROOT::Math::IGenFunction
{
  WrappedGraph(TGraph* g) : fGraph(g) {}

  unsigned int NDim() const { return 1; }
  ROOT::Math::IGenFunction* Clone() const {return new WrappedGraph(fGraph); }
  double DoEval(double x) const { return fGraph->Eval(x); }
  TGraph* fGraph;

  
};


