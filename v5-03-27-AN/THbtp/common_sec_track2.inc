CCC   common/sec_trk_map2/
C
CCC   Second Track-to-Sector Mapping data structure:
CCC   (Used for building the reference spectrum from the
CCC    input event text file.)
CCC   (Converted to common block for ALICE application, 12/1/2000)

      Integer*4 max_trk_sec2
      Integer*4 sec_maxlen2
      Parameter (max_trk_sec2 = 150)
      Parameter (sec_maxlen2  = 28000)

      common/sec_trk_map2/
     1  stm2_sec_id(sec_maxlen2),
     2  stm2_n_trk_sec(sec_maxlen2),
     3  stm2_flag(sec_maxlen2),
     4  stm2_track_id(max_trk_sec2,sec_maxlen2)

CCC   Variable Type Declarations:

      Integer*4 stm2_sec_id       !  unique sector ID number
      Integer*4 stm2_n_trk_sec    !  Number of tracks assigned to sector
      Integer*4 stm2_flag         !  normally=0, if = 1 then more than
CCC                               !  max_trk_sec tracks could have been 
CCC                               !  assigned to this sector, however the
CCC                               !  maximum number that can be assigned is
CCC                               !  max_trk_sec.
      Integer*4 stm2_track_id     !  Foreign keys to tracks in /track/ that
CCC                               !  are assigned to this sector.

CCC   END common/sec_trk_map2/
