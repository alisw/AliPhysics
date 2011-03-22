/*
			     dateStream.c
			     ============

   Utility to simulate a DATE raw data stream using a given set of raw
   data files and a configuration file.

   Revision history:

   V01.00  4/05/2004  RD  Created
   V01.01 25/10/2005  RD  Support added for timestamp
   V01.02  4/04/2006  RD  Support for CDH
   V01.03 24/05/2006  RD  Added "Direct disk access" option
*/
#define VID "1.03"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include <cassert>
#include "event.h"

#define DESCRIPTION "DATE raw data stream simulator"
#ifdef AIX
static
#endif
char fileHandlerIdent[]= "@(#)""" __FILE__ """: """ DESCRIPTION \
                         """ """ VID """ """ \
                         """ compiled """ __DATE__ """ """ __TIME__;

#define DBG_BASE     if ( debug > 0 )
#define DBG_DETAILED if ( debug > 1 )
#define DBG_VERBOSE  if ( debug > 2 )

#ifndef TRUE
# define TRUE (0 == 0)
#endif
#ifndef FALSE
# define FALSE (0 == 1)
#endif

const char *myName;
int debug;
FILE *outF;
typedef enum { unknown, ldc, gdc } workingAsType;
typedef enum { collider, fixedTarget } workingModeType;
workingAsType workingAs;
workingModeType workingMode;
struct ldcDescriptorStruct {
  eventLdcIdType id;
  struct ldcDescriptorStruct *next;
} *ldcsHead, *ldcsTail;
void *eventsHead, *eventsTail;
struct gdcEventDescriptorStruct {
  struct ldcEventDescriptorStruct *head;
  struct ldcEventDescriptorStruct *tail;
  struct gdcEventDescriptorStruct *next;
  struct eventHeaderStruct header;
  int loaded;
} *currGdc;
struct ldcEventDescriptorStruct {
  struct equipmentEventDescriptorStruct *head;
  struct equipmentEventDescriptorStruct *tail;
  struct ldcEventDescriptorStruct *next;
  eventLdcIdType id;
  struct eventHeaderStruct header;
  int loaded;
} *currLdc;
struct equipmentEventDescriptorStruct {
  struct equipmentEventDescriptorStruct *next;
  equipmentIdType id;
  struct payloadDescriptorStruct *payload;
  struct equipmentHeaderStruct header;
} *currEvent;
struct payloadDescriptorStruct {
  struct payloadDescriptorStruct *next;
  char *fileName;
  int fileSize;
  int size;
  void *data;
} *payloadsHead, *payloadsTail;
int lineNmb;
eventGdcIdType currGdcId;
unsigned long32 currDetPattern; 
eventLdcIdType currLdcId;
equipmentIdType currEquipmentId;
int currRunNb;
int numOfLdcs;
int numOfEvents;
int createSorEor;
int handleCDH;
eventIdType oneEventDelta;
eventIdType currEventId;
int gotAliceTrigger;
int bufferData;

struct commonDataHeaderStruct *cdhRef = NULL;

void dumpPayload( const struct payloadDescriptorStruct *p ) {
  char *c;
  int i;
  int printable;
	  
  if ( p->data != NULL ) {
    for ( i = 0, c = (char *)p->data, printable = TRUE;
	  printable && i != p->size;
	  c++, i++ )
      printable = isascii( *c );
    if ( printable ) {
      printf( "       \"" );
      for ( i = 0, c = (char *)p->data; i != p->size; c++, i++ ) {
	if ( *c == '\n' )
	  printf( "\"\n       \"" );
	else
	  putchar( *c );
      }
      if ( *c != '\n' ) printf( "\"\n" );
    } else {
      long32 *v;
      for ( i = 0, v = (long32 *)p->data;
	    i+4 <= p->size;
	    v++, i += 4 ) {
	if ( i % (4*8) == 0 ) {
	  if ( i != 0 ) printf( "\n" );
	  printf( "       " );
	}
	printf( "%08x ", *v );
      }
      if ( i < p->size ) {
	int j = 0;

	printf( "\n       " );
	while ( i != p->size ) {
	  printf( "%02x ", *((char *)p->data + p->size - j - 1) & 0xff );
	  j++;
	  i++;
	}
      }
    }
    printf( "\n" );
  }
} /* End of dumpPayload */

void dumpEvents() {
  assert( workingAs == ldc || workingAs == gdc );
  if ( eventsHead != NULL ) {
    printf( "Events:\n" );
    if ( workingAs == gdc ) {
      struct gdcEventDescriptorStruct *gdc;

      for ( gdc = (struct gdcEventDescriptorStruct *)eventsHead;
	    gdc != NULL;
	    gdc = gdc->next ) {
	struct ldcEventDescriptorStruct *ldc;

	printf( " GDC (%p)\n", (void*)gdc );
	for ( ldc = gdc->head; ldc != NULL; ldc = ldc->next ) {
	  struct equipmentEventDescriptorStruct *eq;
	    
	  printf( "   LDC (%p): %d\n", (void*)ldc, ldc->id );
	  for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	    printf( "     EQUIPMENT (%p): %d PAYLOAD (%p):",
		    (void*)eq,
		    eq->id,
		    (void*)eq->payload );
	    fflush( stdout );
	    printf( "\"%s\" (%d bytes)\n",
		    eq->payload->fileName,
		    eq->payload->size );
	    dumpPayload( eq->payload );
	  }
	}
      }
    }
    if ( workingAs == ldc ) {
      struct ldcEventDescriptorStruct *ldc;

      for ( ldc = (struct ldcEventDescriptorStruct *)eventsHead;
	    ldc != NULL;
	    ldc = ldc->next ) {
	struct equipmentEventDescriptorStruct *eq;
	    
	printf( "   LDC\n" );
	for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	  printf( "     EQUIPMENT (%p): %d PAYLOAD (%p):",
		  (void*)eq,
		  eq->id,
		  (void*)eq->payload );
	  fflush( stdout );
	  printf( "\"%s\" (%d bytes)\n",
		  eq->payload->fileName,
		  eq->payload->size );
	  dumpPayload( eq->payload );
	}
      }
    }
  } else {
    printf( "Events: EMPTY\n" );
  }
} /* End of dumpEvents */

void getLine( char *line, const int maxSize ) {
  int read;
  int c;

  for ( read = 0; !feof( stdin ) && !ferror( stdin ) && read != maxSize; read++ ) {
    if ( (line[read] = getchar()) == '\n' ) break;
  }
  if ( ferror( stdin ) ) {
    fprintf( stderr,
	     "%s: failed to read configuration input errno:%d ",
	     myName, errno );
    perror( "" );
    exit( 1 );
  }
  if ( feof( stdin ) ) read--;
  if ( read == maxSize && line[read] != '\n' ) {
    fprintf( stderr,
	     "%s: Input line # %d too long (%d chars max)\n",
	     myName, lineNmb, maxSize-1 );
    exit( 1 );
  }
  line[ read ] = 0;
  DBG_VERBOSE {
    if ( !( read == 0 && feof( stdin ) ) ) {
      printf( "%d) [%3d] \"%s\"", lineNmb, read, line );
    }
  }
  for ( c = 0; c != read; c++ ) {
    if ( line[c] == '#' ) {
      line[c] = 0;
      break;
    }
  }
  DBG_VERBOSE {
    if ( read != c ) {
      printf( " => \"%s\"", line );
    }
    if ( feof( stdin ) ) printf( "<<< EOF >>>" );
    if ( ferror( stdin ) ) printf( "<<< FERROR >>>" );
    printf( "\n" );
  }
} /* End of getLine */

void handleLdc( eventLdcIdType ldcId ) {
  struct ldcDescriptorStruct *ldc;

  if ( ldcsHead != NULL ) {
    for ( ldc = ldcsHead; ldc != NULL; ldc = ldc->next ) {
      if ( ldc->id == ldcId ) {
	return;
      }
    }
  }
  if ( (ldc = (struct ldcDescriptorStruct *)malloc( sizeof( *ldc ) )) == NULL ) {
    fprintf( stderr,
	     "%s: Failed to malloc for %d bytes (struct ldcDescriptorStruct)\n",
	     myName, (int)sizeof( *ldc ) );
    exit( 1 );
  }
  ldc->id = ldcId;
  ldc->next = NULL;
  if ( ldcsHead == NULL ) {
    ldcsHead = ldcsTail = ldc;
  } else {
    ldcsTail->next = ldc;
    ldcsTail = ldc;
  }
  numOfLdcs++;
} /* End of handleLdc */

void createNewEvent() {
  assert( workingAs == ldc || workingAs == gdc );
  if ( workingAs == ldc ) {
    struct ldcEventDescriptorStruct *p;

    if ( (p = (struct ldcEventDescriptorStruct *)malloc( sizeof( *p ) ))
	   == NULL ) {
      fprintf( stderr,
	       "%s: failed to malloc for %d bytes (createNewEvent: struct ldcEventDescriptorStruct)",
	       myName, (int)sizeof( *p ) );
      perror( "" );
      exit( 1 );
    }
    p->loaded = FALSE;
    p->head = p->tail = NULL;
    p->next = NULL;
    currLdc = p;
    if ( eventsHead == NULL ) {
      eventsHead = eventsTail = p;
    } else {
      struct ldcEventDescriptorStruct *q =
	(struct ldcEventDescriptorStruct *)eventsTail;

      q->next = p;
      eventsTail = p;
    }
    p->id = currLdcId;
  } else if ( workingAs == gdc ) {
    struct gdcEventDescriptorStruct *p;

    if ( (p = (struct gdcEventDescriptorStruct *)malloc( sizeof( *p ) ))
	   == NULL ) {
      fprintf( stderr,
	       "%s: failed to malloc for %d bytes (createNewEvent: struct gdcEventDescriptorStruct)",
	       myName, (int)sizeof( *p ) );
      perror( "" );
      exit( 1 );
    }
    p->loaded = FALSE;
    p->next = NULL;
    p->head = p->tail = NULL;
    currGdc = p;
    if ( eventsHead == NULL ) {
      eventsHead = eventsTail = p;
    } else {
      struct gdcEventDescriptorStruct *q =
	(struct gdcEventDescriptorStruct *)eventsTail;

      q->next = p;
      eventsTail = p;
    }
  }
} /* End of createNewEvent */

void createNewLdcEvent() {
  struct gdcEventDescriptorStruct *gdcDesc;
  struct ldcEventDescriptorStruct *p;

  if ( (p = (struct ldcEventDescriptorStruct *)malloc( sizeof( *p ) ))
         == NULL ) {
    fprintf( stderr,
	     "%s: failed to malloc for %d bytes (createNewLdcEvent: struct ldcEventDescriptorStruct)",
	     myName, (int)sizeof( *p ) );
    perror( "" );
    exit( 1 );
  }
  p->id = currLdcId;
  p->head = p->tail = NULL;
  p->next = NULL;
  gdcDesc = (struct gdcEventDescriptorStruct *)eventsTail;
  if ( gdcDesc->head == NULL ) {
    gdcDesc->head = gdcDesc->tail = p;
  } else {
    gdcDesc->tail->next = p;
    gdcDesc->tail = p;
  }
  currLdc = p;
} /* End of createNewLdcEvent */

void loadBuffer( struct payloadDescriptorStruct * const payload ) {
  FILE *f;
  int bytesRead;

  if ( (f = fopen( payload->fileName, "r" )) == NULL ) {
    fprintf( stderr,
	     "%s: line:%d payload file \"%s\" not found or not readable, errno:%d. ",
	     myName,
	     lineNmb,
	     payload->fileName,
	     errno );
    perror( "System-dependent error " );
    exit( 1 );
  }
  if ( (payload->data = malloc( payload->size )) == NULL ) {
    fprintf( stderr,
	     "%s: line:%d Failed to malloc for payload file \"%s\" size:%d errno:%d ",
	     myName,
	     lineNmb,
	     payload->fileName,
	     payload->size,
	     errno );
    perror( "System-dependent status " );
    exit( 1 );
  }
  if ( (bytesRead = fread( payload->data, payload->fileSize, 1, f )) != 1 ) {
    fprintf( stderr,
	     "%s: line:%d Failed to read payload file \"%s\" size:%d requested:1 got:%d feof:%s ferror:%s errno:%d ",
	     myName,
	     lineNmb,
	     payload->fileName,
	     payload->size,
	     bytesRead,
	     feof(f) ? "TRUE" : "false",
	     ferror(f) ? "TRUE" : "false",
	     errno );
    perror( "System-dependent status " );
    exit( 1 );
  }
  fclose(f);
  if ( payload->size != payload->fileSize ) {
    memset( (char *)payload->data + payload->fileSize,
	    0,
	    payload->size - payload->fileSize );
  }
} /* End of loadBuffer */

void unloadBuffer( struct payloadDescriptorStruct * const payload ) {
  if ( payload->data != NULL ) {
    free( payload->data );
    payload->data = NULL;
  }
} /* End of unloadBuffer */

void unloadAllBuffers() {
  struct payloadDescriptorStruct *payload;

  for ( payload = payloadsHead; payload != NULL; payload = payload->next ) {
    unloadBuffer( payload );
  }
} /* End of unloadAllBuffers */

void loadPayload( const char *fileName ) {
  struct payloadDescriptorStruct *payload;

  for ( payload = payloadsHead; payload != NULL; payload = payload->next ) {
    if ( strcmp( fileName, payload->fileName ) == 0 )
      break;
  }
  if ( payload == NULL ) {
    FILE *f;

    if ( (payload = (struct payloadDescriptorStruct *)malloc( sizeof( *payload ) ))
	   == NULL ) {
      fprintf( stderr,
	       "%s: failed to malloc for %d bytes (loadPayload/payloadDescriptorStruct)\n",
	       myName,
	       (int)sizeof( *payload ) );
      exit( 1 );
    }
    if ( (payload->fileName = strdup( fileName )) == NULL ) {
      fprintf( stderr,
	       "%s: failed to duplicate string \"%s\" (loadPaload/fileName)\n",
	       myName,
	       fileName );
      exit( 1 );
    }
    if ( (f = fopen( fileName, "r" )) == NULL ) {
      fprintf( stderr,
	       "%s: line:%d payload file \"%s\" not found or not readable, errno:%d. ",
	       myName,
	       lineNmb,
	       fileName,
	       errno );
      perror( "System-dependent error " );
      exit( 1 );
    }
    if ( fseek( f, 0L, SEEK_END ) != 0 ) {
      fprintf( stderr,
	       "%s: line:%d Failed to seek payload file \"%s\" errno:%d ",
	       myName,
	       lineNmb,
	       fileName,
	       errno );
      perror( "System-dependent error " );
      exit( 1 );
    }
    if ( (payload->size = ftell( f )) <= 0 ) {
      fprintf( stderr,
	       "%s: line:%d Failed to get file \"%s\" size size:%d errno:%d ",
	       myName,
	       lineNmb,
	       fileName,
	       payload->size,
	       errno );
      perror( "System-dependent status " );
      exit( 1 );
    }
    payload->fileSize = payload->size;
    while ( (payload->size & 3) != 0 ) payload->size++;
    fclose( f );

    if ( bufferData ) {
      loadBuffer( payload );
    } else {
      payload->data = NULL;
    }

    payload->next = NULL;
    if ( payloadsHead == NULL ) {
      payloadsHead = payloadsTail = payload;
    } else {
      payloadsTail->next = payload;
      payloadsTail = payload;
    }
    DBG_VERBOSE {
      int b, n;

      printf( "%d)       Payload \"%s\" loaded at %p\n",
	      lineNmb,
	      fileName,
	      (void*)payload );
      if ( bufferData ) {
	if ( handleCDH &&
	     strncmp(fileName,"TRG_",4) != 0 ) {
	  struct commonDataHeaderStruct *cdh =
	    (struct commonDataHeaderStruct *)payload->data;

	  printf( " CDH: blockLenght:%d=0x%08x ",
		  cdh->cdhBlockLength, cdh->cdhBlockLength );
	  if ( cdh->cdhBlockLength < sizeof( *cdh ) ) {
	    printf( "TOO SMALL (minimum:%ld=0x%08lx)\n",
		   (unsigned long)sizeof( *cdh ),
		   (unsigned long)sizeof( *cdh ) );
	  } else {
	    printf( "version:%d=0x%x ", cdh->cdhVersion, cdh->cdhVersion );
	    if ( cdh->cdhVersion != CDH_VERSION ) {
	      printf( "EXPECTED:%d=%x (decoding may be inaccurate) ",
		      CDH_VERSION, CDH_VERSION );
	    }
	  }
	  printf( "L1TriggerMessage:0x%x", cdh->cdhL1TriggerMessage );
	  if (  cdh->cdhL1TriggerMessage != 0 ) {
	    for ( b = 0, n = 0; b != 10; b++ ) {
	      if ( (cdh->cdhL1TriggerMessage & (1<<b)) != 0 ) {
		if ( n++ != 0 )printf( "+" );
		switch (b) {
		case 0: printf( "L1SwC" ); break;
		case 1: printf( "ESR" ); break;
		case 2: printf( "RoC1" ); break;
		case 3: printf( "RoC2" ); break;
		case 4: printf( "RoC3" ); break;
		case 5: printf( "RoC4" ); break;
		case 6: printf( "ClT" ); break;
		default: printf( "spare %d", b+14 );
		}
	      }
	    }
	    printf( ">" );
	  }
	  printf( " " );
	  if ( cdh->cdhMBZ0 != 0 )
	    printf( "MBZ0:0x%x ",
		    cdh->cdhMBZ0 );
	  printf( "\n" );
	  
	  printf( "      " );
	  printf( "EventId2(orbit):%d=0x%x ",
		  cdh->cdhEventId2, cdh->cdhEventId2 );
	  printf( "EventId1(bunchCrossing):%d=0x%x ",
		  cdh->cdhEventId1, cdh->cdhEventId1 );
	  printf( "\n" );
	  
	  printf( "      " );
	  if ( cdh->cdhMBZ1 != 0 )
	    printf( "MBZ1:0x%x ",
		    cdh->cdhMBZ1 );
	  printf( "BlockAttributes:0x%x",
		  cdh->cdhBlockAttributes );
	  if ( cdh->cdhBlockAttributes != 0 ) {
	    printf( "=<" );
	    for ( b = 0, n = 0; b != 8; b++ ) {
	      if ( (cdh->cdhBlockAttributes & (1<<b)) != 0 ) {
		if ( n++ != 0 )
		  printf( "+" );
		printf( "%d", b );
	      }
	    }
	    printf( ">" );
	  }
	  printf( " " );
	  printf( "ParticipatingSubDetectors:0x%x ",
		  cdh->cdhParticipatingSubDetectors );
	  printf( "\n" );
	  printf( "      " );
	  
	  if ( cdh->cdhMBZ2 != 0 )
	    printf( "MBZ2:0x%x ",
		    cdh->cdhMBZ2 );
	  printf( "Status/Error:0x%x", cdh->cdhStatusErrorBits );
	  if ( cdh->cdhStatusErrorBits != 0 ) {
	    printf( "=<" );
	    for ( b = 0,n = 0; b != 16; b++ ) {
	      if ( (cdh->cdhStatusErrorBits & (1<<b)) != 0 ) {
		if ( n++ != 0 ) printf( "+" );
		switch (b) {
		case 0: printf( "TriggerOverLapError" ); break;
		case 1: printf( "TriggerMissingError" ); break;
		case 2: printf( "DataParityError" ); break;
		case 3: printf( "ControlParityError" ); break;
		case 4: printf( "TriggerInformationUnavailable" ); break;
		case 5: printf( "FEEError" ); break;
		case 6: printf( "HLTDecision" ); break;
		case 7: printf( "HLTPayload" ); break;
		case 8: printf( "DDGPayload" ); break;
		default: printf( "spare %d", b );
		}
	      }
	    }
	    printf( ">" );
	  }
	  printf( " " );
	  printf( "MiniEventId(bunchCrossing):%d=0x%x ",
		  cdh->cdhMiniEventId, cdh->cdhMiniEventId );
	  printf( "\n" );
	  
	  printf( "      " );
	  printf( "Trigger classes: 0x(%05x-%08x)",
		  cdh->cdhTriggerClassesHigh,
		  cdh->cdhTriggerClassesLow );
	  if ( cdh->cdhTriggerClassesHigh != 0
	       || cdh->cdhTriggerClassesLow != 0 ) {
	    printf( "=<" );
	    for ( b=0, n=0; b != 32; b++ ) {
	      if ( (cdh->cdhTriggerClassesLow & (1<<b)) != 0 ) {
		if ( n++ != 0 ) printf( "+" );
		printf( "%d", b );
	      }
	    }
	    for ( b=0; b != 18; b++ ) {
	      if ( (cdh->cdhTriggerClassesHigh & (1<<b)) != 0 ) {
		if ( n++ != 0 ) printf( "+" );
		printf( "%d", b+32 );
	      }
	    }
	    printf( ">" );
	  }
	  printf( "\n" );
	  
	  printf( "      " );
	  if ( cdh->cdhMBZ3 != 0 ) {
	    printf( "MBZ3:0x%x ",
		    cdh->cdhMBZ3 );
	  }
	  printf( "ROI:0x(%08x-%01x)", cdh->cdhRoiHigh, cdh->cdhRoiLow );
	  if ( cdh->cdhRoiHigh != 0
	       || cdh->cdhRoiLow != 0 ) {
	    printf( "=<" );
	    for ( b=0, n=0; b != 5; b++ ) {
	      if ( (cdh->cdhRoiLow & (1<<b)) != 0 ) {
		if ( n++ != 0 ) printf( "+" );
		printf( "%d", b );
	      }
	    }
	    for ( b=0; b != 32; b++ ) {
	      if ( (cdh->cdhRoiHigh & (1<<b)) != 0 ) {
		if ( n++ != 0 ) printf( "+" );
		printf( "%d", b+4 );
	      }
	    }
	    printf( ">" );
	  }
	  printf( "\n" );
	}
      }
    }
  } else {
    DBG_VERBOSE
      printf( "%d)       Payload \"%s\" already loaded at %p\n",
	      lineNmb,
	      fileName,
	      (void*)payload );
  }

  currEvent->payload = payload;
} /* End of loadPayload */

void parseEquipment( char * const line ) {
  struct equipmentEventDescriptorStruct *equipment;
  int payloadFound = FALSE;
  char *p;
  char *keyword;

  if ( (equipment =
	 (struct equipmentEventDescriptorStruct *)malloc( sizeof( *equipment ) )) == NULL ) {
    fprintf( stderr,
	     "%s: filed to malloc for %d bytes (parseEquipment/equipmentEventDescriptorStruct) errno:%d ",
	     myName,
	     (int)sizeof( *equipment ),
	     errno );
    perror( "" );
    exit( 1 );
  }
  currEvent = equipment;

  p = line;
  while ( (keyword = strtok_r( p, " \t", &p )) != NULL ) {
    DBG_VERBOSE printf( "%d)     Equipment - Keyword:\"%s\"\n",
			lineNmb,
			keyword );
    if ( strcasecmp( "id", keyword ) == 0 ) {
      char *idNum;

      if ( (idNum = strtok_r( p, " \t", &p )) == NULL ) {
	fprintf( stderr,
		 "%s: line:%d EQUIPMENT declaration, ID needed",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      if ( sscanf( idNum, "%d", &currEquipmentId ) != 1 ) {
	fprintf( stderr,
		 "%s: line:%d EQUIPMENT declaration, numeric ID needed (%s)",
		 myName,
		 lineNmb,
		 idNum );
	exit( 1 );
      }
      DBG_VERBOSE printf( "%d)     EQUIPMENT - ID:%d\n",
			  lineNmb,
			  currEquipmentId );
    } else if ( strncasecmp( "pay", keyword, 3 ) == 0 ) {
      char *fileName;

      if ( (fileName = strtok_r( p, " \t", &p )) == NULL ) {
	fprintf( stderr,
		 "%s line:%d Payload without filename found\n",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      DBG_VERBOSE printf( "%d)     Equipment - Payload:\"%s\"\n",
			  lineNmb,
			  fileName );
      if ( payloadFound ) {
	fprintf( stderr,
		 "%s line:%d Payload with multiple filenames found\n",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      loadPayload( fileName );
      payloadFound = TRUE;
    } else {
      fprintf( stderr,
	       "%s: line:%d Equipment declaration, unknown keyword \"%s\"\n",
	       myName,
	       lineNmb,
	       keyword );
      exit( 1 );
    }
  }
  if ( !payloadFound ) {
    fprintf( stderr,
	     "%s: line:%d Equipment without payload found\n",
	     myName,
	     lineNmb );
    exit( 1 );
  }

  equipment->id = currEquipmentId;
  equipment->next = NULL;
  if ( currLdc->head == NULL ) {
    currLdc->head = currLdc->tail = equipment;
  } else {
    currLdc->tail->next = equipment;
    currLdc->tail = equipment;
  }
} /* End of parseEquipment */

void parseGdc( char * const line ) {
  char *p;
  char *keyword;

  p = line;
  while ( (keyword = strtok_r( p, " \t", &p )) != NULL ) {
    if ( strcasecmp( "id", keyword ) == 0 ) {
      char *idNum;

      if ( (idNum = strtok_r( p, " \t", &p )) == NULL ) {
	fprintf( stderr,
		 "%s: line:%d GDC declaration, ID needed",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      int inCurrGdcId;
      if ( sscanf( idNum, "%d", &inCurrGdcId ) != 1 ) {
	fprintf( stderr,
		 "%s: line:%d GDC declaration, numeric ID needed (%s)",
		 myName,
		 lineNmb,
		 idNum );
	exit( 1 );
      }
      currGdcId = (eventGdcIdType)inCurrGdcId;
      DBG_VERBOSE printf( "%d)     GDC - ID:%d\n",
			  lineNmb,
			  currGdcId );
    } else if ( strcasecmp( "DetectorPattern", keyword ) == 0 ) {
      char *detPattern;

      if ( (detPattern = strtok_r( p, " \t", &p )) == NULL ) {
	fprintf( stderr,
		 "%s: line:%d GDC declaration, DetectorPattern needed",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      if ( sscanf( detPattern, "%u", &currDetPattern ) != 1 ) {
	fprintf( stderr,
		 "%s: line:%d GDC declaration, numeric DetectorPattern needed (%s)",
		 myName,
		 lineNmb,
		 detPattern );
	exit( 1 );
      }
      DBG_VERBOSE printf( "%d)     GDC - DetectorPattern:%u\n",
			  lineNmb,
			  currDetPattern );
    } else {
      fprintf( stderr,
	       "%s: line:%d GDC declaration, unknown keyword \"%s\"\n",
	       myName,
	       lineNmb,
	       keyword );
      exit( 1 );
    }  
  }
} /* End of parseGdc */

void parseLdc( char * const line ) {
  char *p;
  char *keyword;

  p = line;
  while ( (keyword = strtok_r( p, " \t", &p )) != NULL ) {
    if ( strcasecmp( "id", keyword ) == 0 ) {
      char *idNum;

      if ( (idNum = strtok_r( p, " \t", &p )) == NULL ) {
	fprintf( stderr,
		 "%s: line:%d LDC declaration, ID needed",
		 myName,
		 lineNmb );
	exit( 1 );
      }
      int inCurrLdcId;
      if ( sscanf( idNum, "%d", &inCurrLdcId ) != 1 ) {
	fprintf( stderr,
		 "%s: line:%d LDC declaration, numeric ID needed (%s)",
		 myName,
		 lineNmb,
		 idNum );
	exit( 1 );
      }
      currLdcId = (eventLdcIdType)inCurrLdcId;
      DBG_VERBOSE printf( "%d)     LDC - ID:%d\n",
			  lineNmb,
			  currLdcId );
    } else {
      fprintf( stderr,
	       "%s: line:%d LDC declaration, unknown keyword \"%s\"\n",
	       myName,
	       lineNmb,
	       keyword );
      exit( 1 );
    }  
  }
} /* End of parseLdc */

void parseRules() {
  char line[ 1025 ];

  currLdcId = HOST_ID_MIN;
  currGdcId = HOST_ID_MIN;
  currDetPattern = 0;

  for ( lineNmb = 1; !feof( stdin ); lineNmb++ ) {
    getLine( line, sizeof(line) );
    if ( strlen(line) != 0 ) {
      char *p;
      char *keyword;

      if ( (keyword = strtok_r( line, " \t", &p )) != NULL ) {
	DBG_VERBOSE printf( "%d)   Keyword:\"%s\"\n", lineNmb, keyword );
	if ( strcasecmp( "gdc", keyword ) == 0 ) {
	  if ( workingAs != gdc && workingAs != unknown ) {
	    fprintf( stderr,
		     "%s: line:%d GDC found when working in non-GDC mode (e.g. as a LDC)\n",
		     myName, lineNmb );
	    exit( 1 );
	  }
	  workingAs = gdc;
	  parseGdc( p );
	  createNewEvent();
	  currLdcId = HOST_ID_MIN;
	  currLdc = NULL;
	  currEquipmentId = 0;
	} else if ( strcasecmp( "ldc", keyword ) == 0 ) {
	  if ( workingAs != gdc && workingAs != ldc && workingAs != unknown ) {
	    fprintf( stderr,
		     "%s: line:%d LDC found when working in non-LDC/GDC mode\n",
		     myName, lineNmb );
	    exit( 1 );
	  }
	  if ( workingAs == unknown ) workingAs = ldc;
	  parseLdc( p );
	  if ( workingAs == ldc ) {
	    createNewEvent();
	    currEquipmentId = 0;
	  } else {
	    createNewLdcEvent();
	    handleLdc( currLdcId );
	    currLdcId++;
	  }
	  currEvent = NULL;
	} else if ( strncasecmp( "equ", keyword, 3 ) == 0 ) {
	  if ( workingAs == unknown
	    || (workingAs == ldc && currLdc == NULL )
	    || (workingAs == gdc && currGdc == NULL ) ) {
	    fprintf( stderr,
		     "%s: line:%d Unexpected EQUIPMENT declaration (LDC or GDC needed first)\n",
		     myName,
		     lineNmb );
	    exit( 1 );
	  }
	  parseEquipment( p );
	  currEquipmentId++;
	} else {
	  fprintf( stderr,
		   "%s: line:%d Parse error in \"%s\" unknown keyword\n",
		   myName,
		   lineNmb,
		   keyword );
	  exit( 1 );
	}
      }
    }
  } while ( !feof( stdin ) ) {}
  lineNmb -= 2;

  DBG_VERBOSE {
    printf( "End of parse: %d line%s found\n",
	    lineNmb,
	    lineNmb != 1 ? "s" : "" );
    printf( "Working as %s\n",
	    workingAs == gdc ? "GDC" :
	     workingAs == ldc ? "LDC" :
	      "UNKNOWN" );
    if ( workingAs == gdc ) {
      struct ldcDescriptorStruct *ldc;

      printf( "LDCs (%d):", numOfLdcs );
      for ( ldc = ldcsHead; ldc != NULL; ldc = ldc->next ) {
	printf( " %d", ldc->id );
      }
      printf( "\n" );
    }
    dumpEvents();
  }

  if ( workingAs == ldc ) {
    assert( ldcsHead == ldcsTail );
    assert( ldcsTail == NULL );
  }

  if ( workingAs == gdc ) {
    struct ldcDescriptorStruct *ldc;

    assert( ldcsHead != NULL );
    assert( ldcsTail != NULL );
    assert( ldcsTail->next == NULL );
    for ( ldc = ldcsHead; ldc->next != NULL; ldc = ldc->next ) {}
    assert ( ldc == ldcsTail );
  }

  if ( workingAs == unknown ) {
    DBG_VERBOSE printf( "Empty configuration: nothing to do!\n" );
    exit( 0 );
  }

  assert( (eventsHead == NULL && eventsTail == NULL)
       || (eventsHead != NULL && eventsTail != NULL) );
} /* End of parseRules */

void loadTimestamp( struct eventHeaderStruct * const ev ) {
  time_t t;

  if ( time( &t ) == (time_t)-1 ) {
    fprintf( stderr,
	     "%s: failed to get system time errno:%d (%s)\n",
	     myName, errno, strerror( errno ) );
    exit( 1 );
  }
  ev->eventTimestamp = (eventTimestampType)t;
} /* End of loadTimestamp */

void initEvent( struct eventHeaderStruct * const ev ) {
  memset( ev, 0, sizeof( *ev ) );

  ev->eventMagic = EVENT_MAGIC_NUMBER;
  ev->eventHeadSize = EVENT_HEAD_BASE_SIZE;
  ev->eventVersion = EVENT_CURRENT_VERSION;
  ev->eventRunNb = currRunNb;
  ZERO_EVENT_ID( ev->eventId );
  ZERO_TRIGGER_PATTERN( ev->eventTriggerPattern );
  ZERO_DETECTOR_PATTERN( ev->eventDetectorPattern );
  RESET_ATTRIBUTES( ev->eventTypeAttribute );
  if ( workingMode == collider )
    SET_SYSTEM_ATTRIBUTE( ev->eventTypeAttribute, ATTR_ORBIT_BC );
  ev->eventLdcId = VOID_ID;
  ev->eventGdcId = VOID_ID;
  loadTimestamp( ev );
} /* End of initEvent */

int Swap(int x)
{
   // Swap the endianess of the integer value 'x'

   return (((x & 0x000000ffU) << 24) | ((x & 0x0000ff00U) <<  8) |
           ((x & 0x00ff0000U) >>  8) | ((x & 0xff000000U) >> 24));
}

void outputEvent( const void * const ev,
		  const int size ) {
  int done;

  DBG_VERBOSE {
    const long32 * const v = (long32 *)ev; 
    printf( "Writing %d bytes @ %p (%d)\n", size, ev, *v );
  }

  // .............................Test endianess..............................
  int temp = 1;
  char* ptemp = (char*) &temp;

  if (ptemp[0]!=1) { // Mac platform: ptemp != 1..............................................................................
     int  bufSize= size; if (bufSize > (int) sizeof(eventHeaderStruct)) { bufSize = sizeof(eventHeaderStruct); }
     char* evTemp = (char*) malloc (bufSize);
     memcpy(evTemp, ev, bufSize);

     if ((bufSize % sizeof(int)) != 0) {
            fprintf( stderr, "%s: size of the input buffer ev is not multiple of 4 (size = %d)\n", myName, bufSize);
            exit( 1 );
          }
     else {
            // Invert header to evTemp.....................................................
            int* buf = (int*) evTemp; 
            for (int i=0; i < (int) (bufSize / sizeof(int)); i++, buf++) {
                 int value = Swap(*buf); 
                 memcpy(evTemp + (i * sizeof(int)), &value, sizeof(int)); 
            }

            // Write inverted header to file...............................................
            if ((done = fwrite( evTemp, bufSize, 1, outF )) != 1 ) {
                 fprintf( stderr, "%s: failed to write inverted header. event size:%d bytes, errno:%d (%s)\n", myName, size, errno, strerror( errno ) );
                 exit( 1 );
            }

            if (size > bufSize) {  // Still theraw-data payload to write (but not inverted, since it is inverted eariler).............
                if ((done = fwrite( (char*)ev + bufSize, size - bufSize, 1, outF )) != 1 ) {
                    fprintf( stderr, "%s: failed to write additional event size:%d bytes, errno:%d (%s)\n", myName, size, errno, strerror( errno ) );
                    exit( 1 );
               }
            }
     }
     free(evTemp);
  }
  else {             // Intel platform: ptemp == 1............................................................................
     if ((done = fwrite( ev, size, 1, outF )) != 1 ) {
          fprintf( stderr, "%s: failed to write event size:%d bytes, errno:%d (%s)\n", myName, size, errno, strerror( errno ) );
          exit( 1 );
     }
  }
} /* End of outputEvent */

void createSorAndEor( const int sor ) {
  unsigned char event[ 1000 ];
  struct eventHeaderStruct *ev;
  struct eventHeaderStruct sev;

  assert( workingAs == ldc || workingAs == gdc );

  if ( !createSorEor ) return;
  ev = (struct eventHeaderStruct *)event;
  initEvent( ev );
  ev->eventSize = sizeof( event );
  ev->eventType = sor ? START_OF_RUN : END_OF_RUN;
  if ( workingMode == fixedTarget )
    LOAD_RAW_EVENT_ID( ev->eventId, 0, 0, 0 );
  else
    LOAD_EVENT_ID( ev->eventId, 0, 0, 0 );
  SET_SYSTEM_ATTRIBUTE( ev->eventTypeAttribute, ATTR_P_START );

  if ( workingAs == ldc ) {
    currLdc = (struct ldcEventDescriptorStruct *)eventsHead;
  }
  if ( workingAs == gdc ) {
    initEvent( &sev );
    sev.eventGdcId = currGdcId;
    ev->eventGdcId = currGdcId;
    currGdc = (struct gdcEventDescriptorStruct *)eventsHead;
    currLdc = currGdc->head;
  }
  ev->eventLdcId = currLdc->id;

  if ( workingAs == ldc ) {
    loadTimestamp( ev );
    outputEvent( ev, ev->eventSize );
  }
  if ( workingAs == gdc ) {
    struct ldcDescriptorStruct *ldc;

    loadTimestamp( ev );

    sev.eventSize = sizeof( sev ) + numOfLdcs * ev->eventSize;
    sev.eventType = sor ? START_OF_RUN : END_OF_RUN ;
    COPY_EVENT_ID( ev->eventId, sev.eventId );
    COPY_SYSTEM_ATTRIBUTES( ev->eventTypeAttribute, sev.eventTypeAttribute );
    SET_SYSTEM_ATTRIBUTE( sev.eventTypeAttribute, ATTR_SUPER_EVENT );
    loadTimestamp( &sev );
    outputEvent( &sev, sizeof( sev ) );

    ev->eventGdcId = currGdcId;
    for ( ldc = ldcsHead; ldc != NULL; ldc = ldc->next ) {
      ev->eventLdcId = ldc->id;
      outputEvent( ev, ev->eventSize );
    }
  }

  ADD_EVENT_ID( ev->eventId, oneEventDelta );
  ev->eventSize = ev->eventSize / 2;
  ev->eventType = sor ? START_OF_RUN_FILES : END_OF_RUN_FILES;
  CLEAR_SYSTEM_ATTRIBUTE( ev->eventTypeAttribute, ATTR_P_START );
  if ( workingAs == ldc ) {
    loadTimestamp( ev );
    outputEvent( ev, ev->eventSize );
  }
  if ( workingAs == gdc ) {
    struct ldcDescriptorStruct *ldc;

    loadTimestamp( ev );

    sev.eventSize = ev->eventSize;
    sev.eventType = sor ? START_OF_RUN_FILES : END_OF_RUN_FILES;
    COPY_EVENT_ID( ev->eventId, sev.eventId );
    COPY_SYSTEM_ATTRIBUTES( ev->eventTypeAttribute, sev.eventTypeAttribute );
    CLEAR_SYSTEM_ATTRIBUTE( sev.eventTypeAttribute, ATTR_SUPER_EVENT );
    outputEvent( &sev, sizeof( sev ) );
    outputEvent( ev, ev->eventSize - sizeof( sev ) );

    sev.eventSize = sizeof( sev ) + ev->eventSize;
    sev.eventType = sor ? START_OF_RUN_FILES : END_OF_RUN_FILES;
    COPY_EVENT_ID( ev->eventId, sev.eventId );
    COPY_SYSTEM_ATTRIBUTES( ev->eventTypeAttribute, sev.eventTypeAttribute );
    SET_SYSTEM_ATTRIBUTE( sev.eventTypeAttribute, ATTR_SUPER_EVENT );

    loadTimestamp( &sev );

    ev->eventGdcId = currGdcId;
    for ( ldc = ldcsHead; ldc != NULL; ldc = ldc->next ) {
      loadTimestamp( &sev );
      outputEvent( &sev, sizeof( sev ) );
      ev->eventLdcId = ldc->id;
      outputEvent( ev, ev->eventSize );
    }
  }

  ADD_EVENT_ID( ev->eventId, oneEventDelta );
  ev->eventSize = sizeof( *ev );
  ev->eventType = sor ? START_OF_RUN : END_OF_RUN;
  SET_SYSTEM_ATTRIBUTE( ev->eventTypeAttribute, ATTR_P_END );
  if ( workingAs == ldc ) {
    loadTimestamp( ev );
    outputEvent( ev, ev->eventSize );
  }
  if ( workingAs == gdc ) {
    struct ldcDescriptorStruct *ldc;

    loadTimestamp( ev );

    sev.eventSize = sizeof( sev ) + numOfLdcs * ev->eventSize;
    sev.eventType = sor ? START_OF_RUN : END_OF_RUN;
    COPY_EVENT_ID( ev->eventId, sev.eventId );
    COPY_SYSTEM_ATTRIBUTES( ev->eventTypeAttribute, sev.eventTypeAttribute );
    SET_SYSTEM_ATTRIBUTE( sev.eventTypeAttribute, ATTR_SUPER_EVENT );
    loadTimestamp( &sev );

    outputEvent( &sev, sizeof( sev ) );

    for ( ldc = ldcsHead; ldc != NULL; ldc = ldc->next ) {
      ev->eventLdcId = ldc->id;
      outputEvent( ev, ev->eventSize );
    }
  }
} /* End of createSorEor */

void createSor() {
  createSorAndEor( TRUE );
} /* End of createSor */

void createEor() {
  createSorAndEor( FALSE );
} /* End of createEor */

void loadCdh( struct commonDataHeaderStruct * const cdh,
	             eventIdType            * const eventId,
	             equipmentIdType id ) {
  if ( !handleCDH ) return;

  // CTP raw-data does not contain CDH
  if ( id == 4352) return;

  if ( gotAliceTrigger ) {
    cdh->cdhEventId1 = EVENT_ID_GET_BUNCH_CROSSING( *eventId );
    cdh->cdhEventId2 = EVENT_ID_GET_ORBIT( *eventId );
  } else {
    cdh->cdhEventId1 = 0;
    cdh->cdhEventId2 = EVENT_ID_GET_NB_IN_RUN( *eventId );
  }
  cdh->cdhMiniEventId = cdh->cdhEventId1;
}
void decodeCDH( struct ldcEventDescriptorStruct       * const ldc,
		const struct payloadDescriptorStruct  * const payloadDesc,
	        equipmentIdType id );

void createEvent( void ) {
  assert( workingAs == ldc || workingAs == gdc );

  /* Step 1: load all buffers (if needed) and compose the GDC/LDC headers */
  if ( workingAs == gdc ) {
    struct ldcEventDescriptorStruct *ldc;

    for( ldc = currGdc->head; ldc != NULL; ldc = ldc->next ) {
      COPY_EVENT_ID( currEventId, ldc->header.eventId );
      loadTimestamp( &ldc->header );
    }
    COPY_EVENT_ID( currEventId, currGdc->header.eventId );
    loadTimestamp( &currGdc->header );

    for( ldc = currGdc->head; ldc != NULL; ldc = ldc->next ) {
      struct equipmentEventDescriptorStruct *eq;
      int n;

      for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	if ( !bufferData ) {
	  loadBuffer( eq->payload );
	  decodeCDH( ldc, eq->payload, eq->id );
	}
	loadCdh( (struct commonDataHeaderStruct*)eq->payload->data,
		 &currEventId,
		 eq->id);
      }

      if ( !currGdc->loaded ) {
	for ( n = 0; n != EVENT_TRIGGER_PATTERN_WORDS; n++ )
	  currGdc->header.eventTriggerPattern[n] |= ldc->header.eventTriggerPattern[n];
	for ( n = 0; n != EVENT_DETECTOR_PATTERN_WORDS; n++ )
	  currGdc->header.eventDetectorPattern[n] |= ldc->header.eventDetectorPattern[n];
	for ( n = 0; n != ALL_ATTRIBUTE_WORDS; n++ )
	  currGdc->header.eventTypeAttribute[n] |= ldc->header.eventTypeAttribute[n];
	currGdc->loaded = TRUE;
      }
    }
    cdhRef = NULL;
  } else if ( workingAs == ldc ) {
    struct equipmentEventDescriptorStruct *eq;

    COPY_EVENT_ID( currEventId, currLdc->header.eventId );
    loadTimestamp( &currLdc->header );

    for ( eq = currLdc->head; eq != NULL; eq = eq->next ) {
      if ( !bufferData ) {
	loadBuffer( eq->payload );
	decodeCDH( currLdc, eq->payload, eq->id );
      }
      loadCdh( (struct commonDataHeaderStruct*)eq->payload->data,
	       &currEventId,
	       eq->id);
      currLdc->loaded = TRUE;
    }
    cdhRef = NULL;
  }
  ADD_EVENT_ID( currEventId, oneEventDelta );

  /* Step 2: output the event */
  if ( workingAs == gdc ) {
    struct ldcEventDescriptorStruct *ldc;

    outputEvent( &currGdc->header, sizeof( currGdc->header ) );

    for( ldc = currGdc->head; ldc != NULL; ldc = ldc->next ) {
      struct equipmentEventDescriptorStruct *eq;

      outputEvent( &ldc->header, sizeof( ldc->header ) );

      for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	outputEvent( &eq->header, sizeof( eq->header ) );
	outputEvent( eq->payload->data, eq->payload->size );
	if ( !bufferData ) unloadBuffer( eq->payload );
      }
    }
    if ( (currGdc = currGdc->next) == NULL )
      currGdc = (struct gdcEventDescriptorStruct *)eventsHead;
  } else if ( workingAs == ldc ) {
    struct equipmentEventDescriptorStruct *eq;

    outputEvent( &currLdc->header, sizeof( currLdc->header ) );

    for ( eq = currLdc->head; eq != NULL; eq = eq->next ) {
      outputEvent( &eq->header, sizeof( eq->header ) );
      outputEvent( eq->payload->data, eq->payload->size );
      if ( !bufferData ) unloadBuffer( eq->payload );
    }
    if ( (currLdc = currLdc->next) == NULL )
      currLdc = (struct ldcEventDescriptorStruct *)eventsHead;
  }
} /* End of createEvent */

void createEvents() {
  int eventNum = 0;

  currGdc = (struct gdcEventDescriptorStruct *)eventsHead;
  currLdc = (struct ldcEventDescriptorStruct *)eventsHead;
  currEvent = NULL;

  createSor();
  for ( eventNum = 0;
	eventNum != numOfEvents && numOfEvents != 0;
	eventNum++ ) {
    createEvent();
  }
  createEor();
} /* End of createEvents */

int usage() {
  fprintf( stderr,
	   "Usage: %s [-?][-d][-i definitionFile][-o outputFile][-# numOfEvents][-s][-F|-C]\n\
   -?                  This text\n\
   -v                  Print version ID and exit\n\
   -d                  Enable debug (repeat for more verbosity)\n\
   -i definitionFile   File with the description of the events to create (default: stdin)\n\
   -o outputFile       File used to store events (default: stdout)\n\
   -# numOfEvents      Number of events to generate (default: 1 event)\n\
   -s                  Do not generate SOR/EOR files (valid only for GDCs)\n\
   -F/-C               Working in Fixed Target (F) or Collider (C) mode\n\
   -c                  Handles CDH\n\
   -D		       Direct disc access (no buffering)\n",
	   myName );
  return 1;
} /* End of usage */

void parseArgs( int argc, char **argv ) {
  int arg = 1;
  int inFileName = -1;
  int outFileName = -1;

  myName = argv[0] ;
  while ( arg < argc ) {
    if ( strcmp( "-?", argv[ arg ] ) == 0 ) {
      usage();
      exit( 0 );
    }
    if ( strcmp( "-i", argv[ arg ] ) == 0 ) {
      if ( ++arg == argc ) exit( usage() );
      inFileName = arg;
      if ( freopen( argv[arg], "r", stdin ) == NULL ){
	fprintf( stderr,
		 "%s: failed to open input definition \"%s\" errno:%d ",
		 myName, argv[arg], errno );
	perror( "" );
	exit( 1 );
      }
    } else if ( strcmp( "-v", argv[ arg ] ) == 0 ) {
      printf( "%s\n", fileHandlerIdent );
      exit( 0 );
    } else if ( strcmp( "-o", argv[ arg ] ) == 0 ) {
      if ( ++arg == argc ) exit( usage() );
      outFileName = arg;
    } else if ( strcmp( "-#", argv[ arg ] ) == 0 ) {
      int n;

      if ( ++arg == argc ) exit( usage() );
      if ( sscanf( argv[ arg ], "%d", &n ) != 1 ) exit( usage() );
      if ( n < 0 ) exit( usage() );
      numOfEvents = n;
    } else if ( strcmp( "-s", argv[ arg ] ) == 0 ) {
      createSorEor = FALSE;
    } else if ( strcmp( "-F", argv[ arg ] ) == 0 ) {
      workingMode = fixedTarget;
    } else if ( strcmp( "-C", argv[ arg ] ) == 0 ) {
      workingMode = collider;
    } else if ( strcmp( "-d", argv[ arg ] ) == 0 ) {
      debug++;
    } else if ( strcmp( "-c", argv[ arg ] ) == 0 ) {
      handleCDH = TRUE;
    } else if ( strcmp( "-D", argv[ arg ] ) == 0 ) {
      bufferData = FALSE;
    } else if ( strcmp( "-run", argv[ arg ] ) == 0 ) {
      int runnumber;
      if ( ++arg == argc ) exit( usage() );
      if ( sscanf( argv[ arg ], "%d", &runnumber ) != 1 ) exit( usage() );
      if ( runnumber < 0 ) exit( usage() );
      currRunNb = runnumber;
    } else {
      fprintf( stderr, "%s: Unknown switch \"%s\"\n", myName, argv[argc] );
      exit( usage() );
    }
    arg++;
  }

  if ( workingMode == fixedTarget )
    LOAD_RAW_EVENT_ID( oneEventDelta, 1, 0, 1 );
  else
    LOAD_EVENT_ID( oneEventDelta, 0, 0, 1 );
  ZERO_EVENT_ID( currEventId );

  DBG_VERBOSE {
    printf( "Configuration:\n" );
    printf( "  Debug level: %d\n", debug );
    printf( "  Configuration: %s\n",
	    inFileName == -1 ? "stdin" : argv[ inFileName ] );
    printf( "  Output: %s\n",
	    outFileName == -1 ? "stdout" : argv[ outFileName ] );
    printf( "  Working mode: %s\n",
	    workingMode == fixedTarget ? "fixed target" : "collider" );
    printf( "  Number of events: %d\n", numOfEvents );
    printf( "  %s SOR/EOR files\n",
	    createSorEor ? "Create" : "Do not create" );
    printf( "  CDH handling: %s\n",
	    handleCDH ? "enabled" : "disabled" );
    printf( "  data buffering: %s\n",
	    bufferData ? "enabled" : "DISABLED" );
  }

  if ( outFileName == -1 ) {
    DBG_BASE
      printf( "No more trace information from this point...\n" );
    debug = 0;
    outF = stdout;
  } else {
    if ( ( outF = fopen( argv[ outFileName ], "w" ) ) == NULL ) {
      fprintf( stderr,
	       "%s: failed to open output file \"%s\" for writing, errno:%d (%s)\n",
	       myName,
	       argv[ outFileName ],
	       errno,
	       strerror( errno ) );
      exit( 1 );
    }
    DBG_DETAILED
      printf( "Output file \"%s\" opened OK for writing\n",
	      argv[ outFileName ] );
  }
} /* End of parseArgs */
void initEquipment( struct equipmentHeaderStruct * const eq ) {
  memset( eq, 0, sizeof( *eq ) );
  RESET_ATTRIBUTES( eq->equipmentTypeAttribute );
  eq->equipmentBasicElementSize = 4;
} /* End of initEquipment */

void decodeCDH(       struct ldcEventDescriptorStruct * const ldc,
		const struct payloadDescriptorStruct  * const payloadDesc,
		      equipmentIdType id ) {
  if ( handleCDH && 
       id != 4352 ) {
    struct commonDataHeaderStruct *cdh;
    static int softwareTriggerIndicator = FALSE;
    int attr;
    int trig;

    if ( payloadDesc->size < CDH_SIZE ) {
      fprintf( stderr,
	       "%s: payload too small got:%d CDH:%d\n",
	       myName,
	       payloadDesc->size,
	       CDH_SIZE );
      exit( 1 );
    }
    if ( (cdh = (struct commonDataHeaderStruct *)payloadDesc->data) != NULL ) {
      if ( cdh->cdhVersion != CDH_VERSION ) {
	fprintf( stderr,
		 "%s: CDH version mismatch expected:%d got:%d\n",
		 myName,
		 CDH_VERSION,
		 cdh->cdhVersion );
	exit( 1 );
      }
      if ( cdhRef == NULL ) {
	cdhRef = cdh;
#define CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK (1<<CDH_TRIGGER_INFORMATION_UNAVAILABLE_BIT)
	gotAliceTrigger = (cdhRef->cdhStatusErrorBits & CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK) == 0;
	if ( gotAliceTrigger && workingMode == fixedTarget ) {
	  fprintf( stderr,
		   "%s: ALICE trigger and fixed target mode are not compatible.\n\
Either work in Collider mode or set the trigger unavailable status bit in the CDH.\n",
		   myName );
	  exit( 1 );
	}
	if ( gotAliceTrigger ) {
	  if ( (cdh->cdhL1TriggerMessage & 0x40) != 0 ) {
	    fprintf( stderr,
		     "%s: CDH is a calibration trigger (unsupported) L1TriggerMessage:0x%x\n",
		     myName, cdh->cdhL1TriggerMessage );
	    exit( 1 );
	  }
	  if ( (cdh->cdhL1TriggerMessage & 0x01) != 0 ) {
	    softwareTriggerIndicator = TRUE;
	  }
	  if ( softwareTriggerIndicator ) {
	    switch ((cdh->cdhL1TriggerMessage >> 2) & 0xF) {
	    case 0xD:
	    case 0xC:
	    case 0xB:
	    case 0xA:
	    case 0x9:
	      break;
	    case 0xF:
	      /* L1SwC bit = on, Clt bit = off, RoC[4..1] = 0xF --> END_OF_DATA */ 
	    case 0xE:
	      /* L1SwC bit = on, Clt bit = off, RoC[4..1] = 0xE0 --> START_OF_DATA */
	    case 0x8:
	      /*  L1SwC bit = on, Clt bit = off, RoC[4] = 1, but not 0xE or 0xF
		  --> SYSTEM_SOFTWARE_TRIGGER_EVENT */
	    default:
	      /*  L1SwC bit = on, Clt bit = off, RoC[4] = 0
		  --> DETECTOR_SOFTWARE_TRIGGER_EVENT */
	      fprintf( stderr,
		       "%s: CDH trigger SOD/EOD/SST/DST (unsupported) \
L1TriggerMessage:0x%x ALICETrigger:%s\n",
		       myName,
		       cdh->cdhL1TriggerMessage,
		       gotAliceTrigger ? "yes" : "no" );
	      exit( 1 );
	    }
	  }
	}
      } else {
	if ( (cdh->cdhStatusErrorBits & CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK) !=
	     (cdhRef->cdhStatusErrorBits & CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK) ) {
	  fprintf( stderr,
		   "%s: CDH coherency check failed. \
Trigger information reference:%savailable current:%savailable\n",
		   myName,
		   (cdhRef->cdhStatusErrorBits & CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK) == 0 ? "UN" : "",
		   (cdh->cdhStatusErrorBits & CDH_TRIGGER_INFORMATION_UNAVAILABLE_MASK) == 0 ? "UN" : "" );
	  exit( 1 );
	}
	if ( gotAliceTrigger ) {
	  if ( cdhRef->cdhL1TriggerMessage != cdh->cdhL1TriggerMessage ) {
	    fprintf( stderr,
		     "%s: CDH coherency check failed. \
L1 trigger message reference:0x%x current:0x%x\n",
		     myName,
		     cdhRef->cdhL1TriggerMessage,
		     cdh->cdhL1TriggerMessage );
	    exit( 1 );
	  }
	  if ( cdh->cdhParticipatingSubDetectors != cdhRef->cdhParticipatingSubDetectors ) {
	    fprintf( stderr,
		     "%s: CDH coherency check failed. \
ParticipatingSubDetectors reference:0x%x current:0x%x\n",
		     myName,
		     cdhRef->cdhParticipatingSubDetectors,
		     cdh->cdhParticipatingSubDetectors );
	    exit( 1 );
	  }
	  if ( cdh->cdhTriggerClassesLow  != cdhRef->cdhTriggerClassesLow
	       || cdh->cdhTriggerClassesHigh != cdhRef->cdhTriggerClassesHigh ) {
	    fprintf( stderr,
		     "%s: CDH coherency check failed. \
TriggerClassesHigh/Low reference:0x%x-%x current:0x%x-%x\n",
		     myName,
		     cdhRef->cdhTriggerClassesHigh, cdhRef->cdhTriggerClassesLow,
		     cdh   ->cdhTriggerClassesHigh, cdh   ->cdhTriggerClassesLow  );
	    exit( 1 );
	  }
	  if ( cdh->cdhBlockLength != 0xffffffff ) {
	    if ( (unsigned)payloadDesc->size !=  cdh->cdhBlockLength ) {
	      fprintf( stderr,
		       "%s: CDH coherency check failed. \
Payload size:%d (0x%08x) CDH block length:%d (0x%08x)\n",
		       myName,
		       payloadDesc->size, payloadDesc->size,
		       cdh->cdhBlockLength, cdh->cdhBlockLength );
	      exit( 1 );
	    }
	  }
	  if ( cdh->cdhRoiLow  != cdhRef->cdhRoiLow
	       || cdh->cdhRoiHigh != cdhRef->cdhRoiHigh ) {
	    fprintf( stderr,
		     "%s: CDH coherency check failed. \
RoiHigh/Low reference:0x%x-%x current:0x%x-%x\n",
		     myName,
		     cdhRef->cdhRoiHigh, cdhRef->cdhRoiLow,
		     cdh   ->cdhRoiHigh, cdh   ->cdhRoiLow  );
	    exit( 1 );
	  }
	}
	if ( cdh->cdhMBZ0 != 0
	     || cdh->cdhMBZ1 != 0
	     || cdh->cdhMBZ2 != 0
	     || cdh->cdhMBZ3 != 0 ) {
	  fprintf( stderr,
		   "%s: CDH check failed. MBZ0:0x%x MBZ1:0x%x MBZ2:0x%x MBZ3:0x%x\n",
		   myName,
		   cdh->cdhMBZ0, cdh->cdhMBZ1, cdh->cdhMBZ2, cdh->cdhMBZ3 );
	  exit( 1 );
	}
      }
      for ( attr = 0; attr != 8; attr++ ) {
	if ( (cdh->cdhBlockAttributes & (1<<attr)) != 0 ) {
	  SET_USER_ATTRIBUTE( ldc->header.eventTypeAttribute, attr );
	}
      }
      for ( trig = 0; trig != 32; trig++ ) {
	if ( (cdh->cdhTriggerClassesLow & (1<<trig)) != 0 ) {
	  SET_TRIGGER_IN_PATTERN( ldc->header.eventTriggerPattern,
				  trig );
	}
      }
      for ( trig = 0; trig != 18; trig++ ) {
	if ( (cdh->cdhTriggerClassesHigh & (1<<trig)) != 0 ) {
	  SET_TRIGGER_IN_PATTERN( ldc->header.eventTriggerPattern,
				  32+trig );
	}
      }
      if ( gotAliceTrigger )
	VALIDATE_TRIGGER_PATTERN( ldc->header.eventTriggerPattern );
    }
  }
} /* End of decodeCDH */

void initEvents() {
  assert( workingAs == ldc || workingAs == gdc );

  if ( workingAs == gdc ) {
    struct gdcEventDescriptorStruct *gdc;

    for ( gdc = (struct gdcEventDescriptorStruct *)eventsHead;
	  gdc != NULL;
	  gdc = gdc->next ) {
      struct ldcEventDescriptorStruct *ldc;

      initEvent( &gdc->header );
      gdc->header.eventSize = gdc->header.eventHeadSize;
      gdc->header.eventType = PHYSICS_EVENT;
      SET_SYSTEM_ATTRIBUTE( gdc->header.eventTypeAttribute, ATTR_SUPER_EVENT );
      gdc->header.eventGdcId = currGdcId;
      COPY_DETECTOR_PATTERN(&currDetPattern, gdc->header.eventDetectorPattern);
      for ( ldc = gdc->head; ldc != NULL; ldc = ldc->next ) {
	struct equipmentEventDescriptorStruct *eq;

	initEvent( &ldc->header );
	ldc->header.eventSize = ldc->header.eventHeadSize;
	ldc->header.eventType = PHYSICS_EVENT;
	ldc->header.eventGdcId = currGdcId;
	COPY_DETECTOR_PATTERN(&currDetPattern, ldc->header.eventDetectorPattern);
	ldc->header.eventLdcId = ldc->id;
	for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	  initEquipment( &eq->header );
	  eq->header.equipmentId = eq->id;
	  if ( workingMode == collider )
	    SET_SYSTEM_ATTRIBUTE( eq->header.equipmentTypeAttribute,
				  ATTR_ORBIT_BC );
	  eq->header.equipmentSize = eq->payload->size + sizeof( eq->header );
	  ldc->header.eventSize += eq->header.equipmentSize;
	  decodeCDH( ldc, eq->payload, eq->id );
	  OR_ALL_ATTRIBUTES( eq->header.equipmentTypeAttribute,
			     ldc->header.eventTypeAttribute );
	  OR_ALL_ATTRIBUTES( eq->header.equipmentTypeAttribute,
			     gdc->header.eventTypeAttribute );
	}
	gdc->header.eventSize += ldc->header.eventSize;
      }
      cdhRef = NULL;
    }

    DBG_VERBOSE {
      printf( "Headers:\n" );
      for ( gdc = (struct gdcEventDescriptorStruct *)eventsHead;
	    gdc != NULL;
	    gdc = gdc->next ) {
	struct ldcEventDescriptorStruct *ldc;
	
	printf( "   GDC:%d size:%d vers:%08x\n",
		currGdcId,
		gdc->header.eventSize,
		gdc->header.eventVersion);
	for ( ldc = gdc->head; ldc != NULL; ldc = ldc->next ) {
	  struct equipmentEventDescriptorStruct *eq;

	  printf( "      LDC:%d size:%d vers:%08x\n",
		  ldc->id,
		  ldc->header.eventSize,
		  ldc->header.eventVersion );
	  for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	    printf( "         EQ:%d size:%d %spayload:%d\n",
		    eq->id,
		    eq->header.equipmentSize,
		    eq->header.equipmentSize - sizeof( struct equipmentHeaderStruct ) == (unsigned)eq->payload->size ? "" : "-ERROR",
		    eq->payload->size );
	  }
	}
      }
    }
  } else if ( workingAs == ldc ) {
    struct ldcEventDescriptorStruct *ldc;

    for ( ldc = (struct ldcEventDescriptorStruct *)eventsHead;
	  ldc != NULL;
	  ldc = ldc->next ) {
      struct equipmentEventDescriptorStruct *eq;

      initEvent( &ldc->header );
      ldc->header.eventSize = ldc->header.eventHeadSize;
      ldc->header.eventType = PHYSICS_EVENT;
      ldc->header.eventGdcId = VOID_ID;
      ldc->header.eventLdcId = ldc->id;
      for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	initEquipment( &eq->header );
	eq->header.equipmentId = eq->id;
	if ( workingMode == collider )
	  SET_SYSTEM_ATTRIBUTE( eq->header.equipmentTypeAttribute,
				ATTR_ORBIT_BC );
	eq->header.equipmentSize = eq->payload->size + sizeof( eq->header );
	ldc->header.eventSize += eq->header.equipmentSize;
	decodeCDH( ldc, eq->payload, eq->id );
	OR_ALL_ATTRIBUTES( eq->header.equipmentTypeAttribute,
			   ldc->header.eventTypeAttribute );
      }
      cdhRef = NULL;
    }
    DBG_VERBOSE {
      printf( "Headers:\n" );
      for ( ldc = (struct ldcEventDescriptorStruct *)eventsHead;
	    ldc != NULL;
	    ldc = ldc->next ) {
	struct equipmentEventDescriptorStruct *eq;

	printf( "      LDC:%d size:%d vers:%08x\n",
		ldc->id,
		ldc->header.eventSize,
		ldc->header.eventVersion );
	for ( eq = ldc->head; eq != NULL; eq = eq->next ) {
	  printf( "         EQ:%d size:%d %spayload:%d\n",
		  eq->id,
		  eq->header.equipmentSize,
		  eq->header.equipmentSize - sizeof( struct equipmentHeaderStruct ) == (unsigned)eq->payload->size ? "" : "-ERROR",
		  eq->payload->size );
	}
      }
    }
  }
} /* End of initEvents */

void initVars() {
  debug = 0;
  workingAs = unknown;
  workingMode = fixedTarget;
  ldcsHead = ldcsTail = NULL;
  eventsHead = eventsTail = NULL;
  currGdc = NULL;
  currLdc = NULL;
  currEvent = NULL;
  payloadsHead = payloadsTail = NULL;
  currRunNb = -1;
  numOfLdcs = 0;
  numOfEvents = 1;
  createSorEor = TRUE;
  handleCDH = FALSE;
  gotAliceTrigger = TRUE;
  bufferData=TRUE;
} /* End of initVars */

int main( int argc, char **argv ) {
  initVars();
  parseArgs( argc, argv );
  parseRules();
  initEvents();
  createEvents();
  return 0;
} /* End of main */
