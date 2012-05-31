CCC   Second Track Data Structure:
CCC   (This one is needed to build the reference spectrum from
CCC    the event text input file.)

      Integer*4 trk2_maxlen
      Parameter (trk2_maxlen = 6500)

      STRUCTURE /track2/
         Integer*4   id          ! Track ID number
         Integer*4   px_sec      ! px sector number
         Integer*4   py_sec      ! py sector number
         Integer*4   pz_sec      ! pz sector number
         Integer*4   sector      ! unique sector ID number
         Integer*4   flag        ! normally=0,if 1 indicates track assigned
CCC                              ! to sector with too many tracks, if = 1
CCC                              ! then track is not used.  See /sec_trk_map/ 
         Integer*4   out_flag    ! flag indicating track in/out of accept.
CCC                              ! non-zero for track pushed out of accept.
         Integer*4   merge_flag  ! flag indicating track is merged (not used)
         Integer*4   ge_pid      ! Geant particle ID code number
         Integer*4  start_vertex ! From input event file - track's start vrtx
         Integer*4  stop_vertex  ! From input event file - track's stop vrtx
         Integer*4   event_line  ! Line # of track in input event text file
         
         Real*4      px          ! x component of track momentum in GeV/c
         Real*4      py          ! y component of track momentum in GeV/c
         Real*4      pz          ! z component of track momentum in GeV/c
         Real*4      E           ! Total energy of track in GeV
         Real*4      pt          ! pt of track momentum in GeV/c
         Real*4      phi         ! azimuthal angle of track in degrees 
         Real*4      eta         ! pseudorapidity of track
      END STRUCTURE

