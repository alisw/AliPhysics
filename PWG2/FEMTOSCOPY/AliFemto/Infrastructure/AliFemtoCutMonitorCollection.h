#ifndef AliFemtoCutMonitorCollection_hh
#define AliFemtoCutMonitorCollection_hh


//#include <list>
#include <vector>
#if !defined(ST_NO_NAMESPACES)
using std::vector;
#endif
class AliFemtoCutMonitor;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef vector<AliFemtoCutMonitor*, allocator<AliFemtoCutMonitor*> >            AliFemtoCutMonitorCollection;
typedef vector<AliFemtoCutMonitor*, allocator<AliFemtoCutMonitor*> >::iterator  AliFemtoCutMonitorIterator;
#else
typedef vector<AliFemtoCutMonitor*>            AliFemtoCutMonitorCollection;
typedef vector<AliFemtoCutMonitor*>::iterator  AliFemtoCutMonitorIterator;
#endif

#endif
