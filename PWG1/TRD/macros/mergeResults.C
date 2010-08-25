#include "AliTRDpwg1Helper.h"
void mergeResults(Char_t *files, Char_t *file="QAresults.root"){
  AliTRDpwg1Helper::MergeProd(file, files, 10);
}