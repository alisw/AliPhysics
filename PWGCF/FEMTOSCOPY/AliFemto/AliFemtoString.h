#ifndef AliFemtoString_hh
#define AliFemtoString_hh

#ifndef __CINT__

#ifndef AliFemtoString_noCint
#define AliFemtoString_noCint
#include <string>

#if !defined(ST_NO_NAMESPACES)
using std::string;
#endif

typedef string AliFemtoString; //!
#endif

#else

#ifndef AliFemtoString_yesCint
#define AliFemtoString_yesCint
class AliFemtoString; //!
#endif

#endif

#endif
