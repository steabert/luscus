      PROGRAM QE2LUS
      IMPLICIT NONE
      INTEGER ARGC
      INTEGER NATOM
      INTEGER IO
      DOUBLE PRECISION E
      CHARACTER*80 LINE
      CHARACTER*80 ARGV
      LOGICAL FIRST
      ARGC = IARGC()
      CALL GETARG(ARGC, ARGV)
      OPEN(UNIT=26, FILE=ARGV(1:INDEX(ARGV,'.out')-1)//'.lus',
     +     STATUS='UNKNOWN')
      NATOM = 0
      E = 0.0D0
      FIRST = .FALSE.

      OPEN(UNIT=55, IOSTAT=IO, FILE=ARGV)
      DO WHILE(IO .EQ. 0)
        READ(UNIT=55, IOSTAT=IO, FMT='(A80)') LINE
        IF (IO .EQ. 0) THEN
          IF (INDEX(LINE,'number of atoms') .NE. 0) THEN
            READ(LINE(INDEX(LINE,'=')+1:LEN_TRIM(LINE)),*) NATOM
          END IF

          IF (INDEX(LINE,'!    total energy') .NE. 0) THEN
            READ(LINE(INDEX(LINE,'=')+1:LEN_TRIM(LINE)),*) E
          END IF

          IF (INDEX(LINE,'ATOMIC_POSITIONS') .NE. 0) THEN
            IF (.NOT. FIRST) THEN
              FIRST = .TRUE.
            ELSE
              WRITE(26,1001)
            END IF
            CALL RDGEO(NATOM, IO)
            IF (E .NE. 0.0D0) THEN
              WRITE(26,1000) E
              E = 0.0D0
            END IF
          END IF
        END IF
      END DO
      CLOSE(55)
      CLOSE(26)
 1000 FORMAT(1X,'<ENERGY>',/,2X,F15.8,/,1X,'</ENERGY>')
 1001 FORMAT(1X,'<END>',/,1X,'<EDITABLE>',/,2X,'NO',/,1X,'</EDITABLE>')
      END

      SUBROUTINE RDGEO(NATOM, IO)
      IMPLICIT NONE
      INTEGER I
      INTEGER NATOM
      DOUBLE PRECISION XYZ(3)
      INTEGER IO

      CHARACTER*80 LINE

      WRITE(26,1000) NATOM
      DO 10 I = 1, NATOM
        READ(UNIT=55, IOSTAT=IO, FMT='(A80)') LINE
        IF (IO .EQ. 0) THEN
          WRITE(26,*) ' ',LINE(1:LEN_TRIM(LINE))
        END IF
 10   CONTINUE
 1000 FORMAT(2X,I5,/)
      END

