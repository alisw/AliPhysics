#if defined(__CLING__)
#include "AliGFWFilterTask.h"
#endif
AliGFWFilterTask *AddTaskGFWFilter(const char *name) {
  return AddFilterTask(name);
}
