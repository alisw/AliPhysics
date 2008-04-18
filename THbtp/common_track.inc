CCC   common/track/
C
CCC   Track Data Structure:
CCC   (Converted to common block for ALICE application, 12/1/2000)

      Integer*4 trk_maxlen
      Parameter (trk_maxlen = 25000)

      common/track1/
     1  trk_id(trk_maxlen),
     2  trk_px_sec(trk_maxlen),
     3  trk_py_sec(trk_maxlen),
     4  trk_pz_sec(trk_maxlen),
     5  trk_sector(trk_maxlen),
     6  trk_flag(trk_maxlen),
     7  trk_out_flag(trk_maxlen),
     8  trk_merge_flag(trk_maxlen),
     9  trk_ge_pid(trk_maxlen),
     1  trk_start_vertex(trk_maxlen),
     2  trk_stop_vertex(trk_maxlen),
     3  trk_event_line(trk_maxlen),
     4  trk_px(trk_maxlen),
     5  trk_py(trk_maxlen),
     6  trk_pz(trk_maxlen),
     7  trk_E(trk_maxlen),
     8  trk_pt(trk_maxlen),
     9  trk_phi(trk_maxlen),
     1  trk_eta(trk_maxlen)

CCC   Variable Type Declarations:

      Integer*4   trk_id          ! Track ID number
      Integer*4   trk_px_sec      ! px sector number
      Integer*4   trk_py_sec      ! py sector number
      Integer*4   trk_pz_sec      ! pz sector number
      Integer*4   trk_sector      ! unique sector ID number
      Integer*4   trk_flag        ! normally=0,if 1 indicates track assigned
CCC                               ! to sector with too many tracks, if = 1
CCC                               ! then track is not used.  See /sec_trk_map/ 
      Integer*4   trk_out_flag    ! flag indicating track in/out of accept.
CCC                               ! non-zero for track pushed out of accept.
      Integer*4   trk_merge_flag  ! flag indicating track is merged (not used)
      Integer*4   trk_ge_pid      ! Geant particle ID code number
      Integer*4  trk_start_vertex ! From input event file - track's start vrtx
      Integer*4  trk_stop_vertex  ! From input event file - track's stop vrtx
      Integer*4   trk_event_line  ! Line # of track in input event text file
      
      Real*4      trk_px          ! x component of track momentum in GeV/c
      Real*4      trk_py          ! y component of track momentum in GeV/c
      Real*4      trk_pz          ! z component of track momentum in GeV/c
      Real*4      trk_E           ! Total energy of track in GeV
      Real*4      trk_pt          ! pt of track momentum in GeV/c
      Real*4      trk_phi         ! azimuthal angle of track in degrees 
      Real*4      trk_eta         ! pseudorapidity of track

CCC   END common/track/
