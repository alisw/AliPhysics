      SUBROUTINE PHOTOS(ID)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:  General search routine + _GET + _SET
C.
C.    Purpose:  /HEPEVT/ is not anymore a standard at least 
C.              REAL*8 REAL*4 are in use. PHOTOS_GET and PHOTOS_SET
C.              were to be introduced.
C.              
C.
C.    Input Parameters:   ID see routine PPHOTOS_MAKE
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was                          Created at:  21/07/98
C.                                                Last Update: 21/07/98
C.
C.----------------------------------------------------------------------
      CALL PHOTOS_GET
      CALL PHOTOS_MAKE(ID)
      CALL PHOTOS_SET
      END
