#ifndef AliFemtoTrioFctnCollection_hh
#define AliFemtoTrioFctnCollection_hh


#include <list>
#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif
class AliFemtoTrioFctn;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<AliFemtoTrioFctn*, allocator<AliFemtoTrioFctn*> >            AliFemtoTrioFctnCollection;
typedef list<AliFemtoTrioFctn*, allocator<AliFemtoTrioFctn*> >::iterator  AliFemtoTrioFctnIterator;
#else
typedef list<AliFemtoTrioFctn*>            AliFemtoTrioFctnCollection;
typedef list<AliFemtoTrioFctn*>::iterator  AliFemtoTrioFctnIterator;
#endif

#endif
