       subroutine CLEAR_FORTRAN_FILE( FILENAME )
         implicit none
         CHARACTER(LEN=30) :: FILENAME
         write(*,*) ' - (FORTRAN) FUNCTION CALLED!...'
         write(*,*) ' - (FORTRAN) CLEARING FILE ...'
         OPEN( UNIT=30, file=FILENAME, 
     +     FORM='UNFORMATTED' )
         write(*,*) ' - (FORTRAN) CLEARED! ...'
         CLOSE( UNIT=30 )
       end

       subroutine APPEND_MOINT_TO_FILE( FILENAME, DVAL,
     +                P, Q, R, S )
          implicit none
          CHARACTER(LEN=30) :: FILENAME
          DOUBLE PRECISION :: DVAL
          INTEGER :: P, Q, R, S
         OPEN( UNIT=30, file=FILENAME, ACCESS='APPEND',
     +     FORM='UNFORMATTED' )
         WRITE( 30 ) DVAL, P, Q, R, S
         CLOSE( UNIT=30 )
       end
