#include "TCollection.h"
#include "TNamed.h"
#include <stdio.h>

// implement this function
// arguments are owned by the framework, don't delete.
// objects added to the outputs must ne NEW, so if you want to forward
// something from the input, forward a copy.
int process(TCollection* in, TCollection* out)
{
  TIter next(in);
  while (TObject* object = next())
  {
    printf("Received object: %s\n", object->GetName());
  }
  out->Add(new TNamed("result","result"));
  return 0;
}
