CCC   common/sec_trk_map/
C
CCC   Track-to-Sector Mapping data structure:
CCC   (Converted to common block for ALICE application, 12/1/2000)

      Integer*4 max_trk_sec
      Integer*4 sec_maxlen
      Parameter (max_trk_sec = 150)
      Parameter (sec_maxlen  = 28000)

      common/sec_trk_map/
     1  stm_sec_id(sec_maxlen),
     2  stm_n_trk_sec(sec_maxlen),
     3  stm_flag(sec_maxlen),
     4  stm_track_id(max_trk_sec,sec_maxlen)

CCC   Variable Type Declarations:

      Integer*4 stm_sec_id        !  unique sector ID number
      Integer*4 stm_n_trk_sec     !  Number of tracks assigned to sector
      Integer*4 stm_flag          !  normally=0, if = 1 then more than
CCC                               !  max_trk_sec tracks could have been 
CCC                               !  assigned to this sector, however the
CCC                               !  maximum number that can be assigned is
CCC                               !  max_trk_sec.
      Integer*4 stm_track_id      !  Foreign keys to tracks in /track/ that
CCC                               !  are assigned to this sector.

CCC   END common/sec_trk_map/
