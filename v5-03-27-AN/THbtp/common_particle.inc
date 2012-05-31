CCC   common/particle/
C
CCC   Geant Particle Properties:
CCC   (Converted to common block for ALICE application, 12/1/2000)

      Integer*4 part_maxlen
      Parameter (part_maxlen = 50)

      common/particle/
     1  part_id(part_maxlen),
     2  part_charge(part_maxlen),
     3  part_mass(part_maxlen),
     4  part_lifetime(part_maxlen)


CCC   Variable Type Declarations:

      Integer*4   part_id        ! Geant particle ID code number; required
CCC                              ! to be equal to the row number
      Integer*4   part_charge    ! Electric charge in units of |e|
      Real*4      part_mass      ! Rest mass in GeV/c**2
      Real*4      part_lifetime  ! Proper lifetime in sec. 

CCC   END common/particle/
