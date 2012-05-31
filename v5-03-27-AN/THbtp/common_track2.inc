CCC   common/track2/
C
CCC   Second Track Data Structure:
CCC   (This one is needed to build the reference spectrum from
CCC    the event text input file.)
CCC   (Converted to common block for ALICE application, 12/1/2000)

      Integer*4 trk2_maxlen
      Parameter (trk2_maxlen = 25000)

      common/track2/
     1  trk2_id(trk2_maxlen),
     2  trk2_px_sec(trk2_maxlen),
     3  trk2_py_sec(trk2_maxlen),
     4  trk2_pz_sec(trk2_maxlen),
     5  trk2_sector(trk2_maxlen),
     6  trk2_flag(trk2_maxlen),
     7  trk2_out_flag(trk2_maxlen),
     8  trk2_merge_flag(trk2_maxlen),
     9  trk2_ge_pid(trk2_maxlen),
     1  trk2_start_vertex(trk2_maxlen),
     2  trk2_stop_vertex(trk2_maxlen),
     3  trk2_event_line(trk2_maxlen),
     4  trk2_px(trk2_maxlen),
     5  trk2_py(trk2_maxlen),
     6  trk2_pz(trk2_maxlen),
     7  trk2_E(trk2_maxlen),
     8  trk2_pt(trk2_maxlen),
     9  trk2_phi(trk2_maxlen),
     1  trk2_eta(trk2_maxlen)

CCC   Variable Type Declarations:

      Integer*4   trk2_id          ! Track ID number
      Integer*4   trk2_px_sec      ! px sector number
      Integer*4   trk2_py_sec      ! py sector number
      Integer*4   trk2_pz_sec      ! pz sector number
      Integer*4   trk2_sector      ! unique sector ID number
      Integer*4   trk2_flag        ! normally=0,if 1 indicates track assigned
CCC                               ! to sector with too many tracks, if = 1
CCC                               ! then track is not used.  See /sec_trk2_map/ 
      Integer*4   trk2_out_flag    ! flag indicating track in/out of accept.
CCC                               ! non-zero for track pushed out of accept.
      Integer*4   trk2_merge_flag  ! flag indicating track is merged (not used)
      Integer*4   trk2_ge_pid      ! Geant particle ID code number
      Integer*4  trk2_start_vertex ! From input event file - track's start vrtx
      Integer*4  trk2_stop_vertex  ! From input event file - track's stop vrtx
      Integer*4   trk2_event_line  ! Line # of track in input event text file
      
      Real*4      trk2_px          ! x component of track momentum in GeV/c
      Real*4      trk2_py          ! y component of track momentum in GeV/c
      Real*4      trk2_pz          ! z component of track momentum in GeV/c
      Real*4      trk2_E           ! Total energy of track in GeV
      Real*4      trk2_pt          ! pt of track momentum in GeV/c
      Real*4      trk2_phi         ! azimuthal angle of track in degrees 
      Real*4      trk2_eta         ! pseudorapidity of track

CCC   END common/track2/
