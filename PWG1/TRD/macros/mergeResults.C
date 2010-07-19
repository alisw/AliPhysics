#include "helper.C"
void mergeResults(Char_t *files, Char_t *file="QAresults.root"){
  mergeProd(file, files, 10);
}