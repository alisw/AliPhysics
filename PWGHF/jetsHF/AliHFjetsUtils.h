#define COLOR_ERROR      "\033[22;31m"   //RED
#define COLOR_WARNING    "\033[22;31;1m" //ORANGE
#define COLOR_DEBUG      "\033[22;34m"   //BLUE
#define COLOR_INFO       "\033[22;32;1m" //GREEN
#define COLOR_DEFAULT    "\033[m"        //Default

#define MSGERROR(msg)   (COLOR_ERROR   msg COLOR_DEFAULT)
#define MSGWARNING(msg) (COLOR_WARNING msg COLOR_DEFAULT)
#define MSGDEBUG(msg)   (COLOR_DEBUG   msg COLOR_DEFAULT)
#define MSGINFO(msg)    (COLOR_INFO    msg COLOR_DEFAULT)


