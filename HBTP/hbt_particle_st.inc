CCC   Geant Particle Properties:

      Integer*4 part_maxlen
      Parameter (part_maxlen = 50)

      STRUCTURE /particle/
         Integer*4   id          ! Geant particle ID code number; required
CCC                              ! to be equal to the row number
         Integer*4   charge      ! Electric charge in units of |e|
         Real*4      mass        ! Rest mass in GeV/c**2
         Real*4      lifetime    ! Proper lifetime in sec. 
      END STRUCTURE

