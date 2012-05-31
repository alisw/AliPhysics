CCC   Second Track-to-Sector Mapping data structure:
CCC   (Used for building the reference spectrum from the
CCC    input event text file.)

      Integer*4 max_trk_sec2
      Integer*4 sec_maxlen2
      Parameter (max_trk_sec2 = 30)
      Parameter (sec_maxlen2  = 28000)

      STRUCTURE /sec_trk_map2/
         Integer*4 sec_id         !  unique sector ID number
         Integer*4 n_trk_sec      !  Number of tracks assigned to sector
         Integer*4 flag           !  normally=0, if = 1 then more than
CCC                               !  max_trk_sec tracks could have been 
CCC                               !  assigned to this sector, however the
CCC                               !  maximum number that can be assigned is
CCC                               !  max_trk_sec2.
         Integer*4 track_id(max_trk_sec2)
CCC                               !  Foreign keys to tracks in /track2/ that
CCC                               !  are assigned to this sector.
      END STRUCTURE
