/* raw2date.c
**
** Take a raw data file and produce a DATE-formatted data file
**
** Revision history:
**  19/03/03 RD		Created
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef ALI_DATE

#include "event.h"

#ifndef TRUE
# define TRUE (0 == 0)
#endif
#ifndef FALSE
# define FALSE (0 == 1)
#endif

const char *whoAmI;
int numFiles = 0;
char **fileNames = NULL;
int   *dataBlockSizes = NULL;
void **dataBlocks = NULL;
int doSor = FALSE;
int doEor = FALSE;
int numLoops = 1;
struct eventHeaderStruct h;
struct eventHeaderStruct sh;
struct equipmentHeaderStruct eh;

int usage() {
  fprintf( stderr, "Usage: %s [-s][-e][-#=N] FILE [FILE]...\n",
	   whoAmI );
  return FALSE;
}

int handleArgs( const int argc, char * const * const argv ) {
  int inFile = FALSE;
  int arg;
  
  if ( argc <= 1 ) return usage();

  for ( arg = 1; arg < argc; arg++ ) {
    if ( !inFile ) {
      if ( strcmp( argv[arg], "--" ) == 0 ) {
	inFile = TRUE;
	continue;
      }
      if ( strcmp( argv[arg], "-s" ) == 0 ) {
	doSor = TRUE;
	continue;
      }
      if ( strcmp( argv[arg], "-e" ) == 0 ) {
	doEor = TRUE;
	continue;
      }
      if ( strncmp( argv[arg], "-#=", 3 ) == 0 ) {
	if ( sscanf( &argv[arg][3], "%d", &numLoops ) != 1 ) {
	  fprintf( stderr,
		   "Failed to scan \"%s\" (expected: -#=number)\n",
		   argv[arg] );
	  return usage();
	}
	continue;
      }
      if ( argv[arg][0] == '-' ) return usage();
      inFile = TRUE;
    }
    if ( (fileNames = realloc( fileNames,
			       ++numFiles * sizeof( char* ) )) == NULL ) {
      perror( "malloc failed " );
      return FALSE;
    }
    if ( (fileNames[ numFiles-1 ] = malloc( strlen(argv[arg])+1 )) == NULL ) {
      perror( "malloc failed " );
      return FALSE;
    }
    if ( (strcpy( fileNames[ numFiles-1 ], argv[arg] )) == NULL ) {
      perror( "strcpy failed " );
      return FALSE;
    }
  }

  if ( numFiles < 1 ) return usage();
  return TRUE;
}

int readData() {
  int f;

  if ( (dataBlockSizes = malloc( sizeof(int) * numFiles )) == NULL ) {
    perror( "malloc failed " );
    return FALSE;
  }
  if ( (dataBlocks = malloc( sizeof( void * ) * numFiles )) == NULL ) {
    perror( "malloc failed" );
    return FALSE;
  }
  for ( f = 0; f != numFiles; f++ ) {
    FILE *in;
    struct stat statData;

    if ( stat( fileNames[f], &statData ) != 0 ) {
      fprintf( stderr, "Cannot stat file \"%s\"", fileNames[f] );
      perror( " " );
      return FALSE;
    }

    if ( (dataBlockSizes[f] = statData.st_size) < 0 ) {
      fprintf( stderr,
	       "Stat for file \"%s\" returns size: %d\n",
	       fileNames[f], (int)statData.st_size );
      return FALSE;
    }
    if ( dataBlockSizes[f] == 0 ) {
      dataBlocks[f] = NULL;
      continue;
    }
    if ( (dataBlocks[f] = malloc( dataBlockSizes[f] )) == NULL ) {
      fprintf( stderr,
	       "Failed to malloc for file \"%s\" size:%d",
	       fileNames[f], dataBlockSizes[f] );
      perror( " " );
      return FALSE;
    }
    if ( (in = fopen( fileNames[f], "r" )) == NULL ) {
      fprintf( stderr, "Failed to open input file \"%s\"", fileNames[f] );
      perror( " " );
      return FALSE;
    }
    if ( fread( dataBlocks[f],
		dataBlockSizes[f],
		1, in ) != 1 ) {
      fprintf( stderr,
	       "Failed to read from file \"%s\"",
	       fileNames[f] );
      perror( " " );
      return FALSE;
    }
    fclose( in );
  }
  return TRUE;
}

void initStructs() {
  h.eventMagic = EVENT_MAGIC_NUMBER;
  h.eventHeadSize = EVENT_HEAD_BASE_SIZE;
  h.eventVersion = EVENT_CURRENT_VERSION;
  h.eventRunNb = 1;
  ZERO_TRIGGER_PATTERN( h.eventTriggerPattern );
  ZERO_DETECTOR_PATTERN( h.eventDetectorPattern );
  RESET_ATTRIBUTES( h.eventTypeAttribute );
  h.eventLdcId = h.eventGdcId = VOID_ID;
  memcpy( &sh, &h, sizeof( h ) );

  eh.equipmentType = 0;
  eh.equipmentId = 0;
  RESET_ATTRIBUTES( eh.equipmentTypeAttribute );
  eh.equipmentBasicElementSize = 1;
}

void dumpDummy( const int size ) {
  int i;
  char pat;

  for ( i = 0, pat = 0; i != size; i++, pat++ )
    fwrite( &pat, 1, 1, stdout );
}

void createSorEor( eventTypeType eventType ) {
  int l;
  int i;

  sh.eventSize = 10872;
  h.eventType = sh.eventType = eventType;
  LOAD_RAW_EVENT_ID( h.eventId, 1, 0, 1 );
  LOAD_RAW_EVENT_ID( sh.eventId, 1, 0, 1 );
  RESET_ATTRIBUTES( h.eventTypeAttribute );
  RESET_ATTRIBUTES( sh.eventTypeAttribute );
  SET_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_P_START );
  SET_SYSTEM_ATTRIBUTE( sh.eventTypeAttribute, ATTR_P_START );
  h.eventGdcId = sh.eventGdcId = 0;

  for ( i = 0; i != 2; i++ ) {
    h.eventSize = sh.eventSize;
    CLEAR_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_SUPER_EVENT );
    fwrite( &h, sizeof(h), 1, stdout );
    dumpDummy( h.eventSize - EVENT_HEAD_BASE_SIZE );

    SET_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_SUPER_EVENT );
    h.eventSize = sh.eventSize + EVENT_HEAD_BASE_SIZE;
  
    for ( l = 0; l != numFiles; l++ ) {
      fwrite( &h, sizeof(h), 1, stdout );
      sh.eventLdcId = l;
      fwrite( &sh, sizeof(sh), 1, stdout );
      dumpDummy( sh.eventSize - EVENT_HEAD_BASE_SIZE );
    }
    CLEAR_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_P_START );
    CLEAR_SYSTEM_ATTRIBUTE( sh.eventTypeAttribute, ATTR_P_START );
    SET_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_P_END );
    SET_SYSTEM_ATTRIBUTE( sh.eventTypeAttribute, ATTR_P_END );
  }
}

void createSor() {
  createSorEor( START_OF_RUN );
}

void createEor() {
  createSorEor( END_OF_RUN );
}

void createData( const int loopNo ) {
  int n;
  int l;
  int totSize;

  for ( totSize = EVENT_HEAD_BASE_SIZE, l = 0;
	l != numFiles;
	l++ )
    totSize += EVENT_HEAD_BASE_SIZE + sizeof( eh ) + dataBlockSizes[l];
  
  h.eventSize = totSize;
  h.eventType = sh.eventType = PHYSICS_EVENT;
  RESET_ATTRIBUTES( h.eventTypeAttribute );
  SET_SYSTEM_ATTRIBUTE( h.eventTypeAttribute, ATTR_SUPER_EVENT );
  RESET_ATTRIBUTES( sh.eventTypeAttribute );
  h.eventGdcId = sh.eventGdcId = 0;
  
  for ( n = 0; n != loopNo; n++ ) {
    LOAD_RAW_EVENT_ID( h.eventId, n, 0, n );
    LOAD_RAW_EVENT_ID( sh.eventId, n, 0, n );
    fwrite( &h, sizeof(h), 1, stdout );

    for ( l = 0; l != numFiles; l++ ) {
      sh.eventLdcId = l;
      sh.eventSize = EVENT_HEAD_BASE_SIZE + sizeof( eh ) + dataBlockSizes[l];
      eh.equipmentSize = dataBlockSizes[l];
      fwrite( &sh, sizeof(sh), 1, stdout );
      fwrite( &eh, sizeof(eh), 1, stdout );
      fwrite( dataBlocks[l], dataBlockSizes[l], 1, stdout );
    }
  }
}

void createStream() {
  if ( doSor ) createSor();
  createData( numLoops );
  if ( doEor ) createEor();
}

int main( int argc, char **argv ) {
  whoAmI = argv[0];

  if ( !handleArgs( argc, argv ) ) return 1;
  if ( !readData() ) return 1;
  initStructs();
  createStream();
  return 0;
}

#else

int main( int argc, char **argv ) {
  fprintf( stderr, "%s was compiled without DATE\n", argv[0] );
  return 1;
}

#endif
