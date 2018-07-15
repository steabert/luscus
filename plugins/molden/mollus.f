      PROGRAM MOLLUS
      IMPLICIT NONE
      CHARACTER*360 IMEIN, IMEOUT
      INTEGER ARGC

      ARGC = IARGC()
      CALL GETARG(ARGC, IMEIN)

      IF (INDEX(IMEIN,'.molden') .EQ. 0) THEN
C there is no molden file!
        STOP
      END IF
      IMEOUT = IMEIN(1:INDEX(IMEIN,'.molden')-1)//'.lus'

      write(*,*) 'imein = |', imein,'|'
      write(*,*) 'imeout = |', imeout,'|'

      OPEN(UNIT=10, FILE=IMEIN, STATUS='OLD')
      OPEN(UNIT=11, FILE=IMEOUT, STATUS='UNKNOWN')

      CALL COMOLU(IMEOUT)

      END

      SUBROUTINE COMOLU(IMEOUT)
      IMPLICIT NONE
      CHARACTER*360 IMEOUT

      INTEGER IO
      CHARACTER*80 LINE
      INTEGER NGEOM
      INTEGER NATOM
      INTEGER NFREQ
      INTEGER MAXAR
      PARAMETER(MAXAR=300000000)
      INTEGER MAXNG
      PARAMETER(MAXNG=10000)
      DOUBLE PRECISION XYZ(MAXNG,3)
      DOUBLE PRECISION ENERGY(MAXNG), MAX_G(MAXNG), RMS_G(MAXNG)
      DOUBLE PRECISION FREQ(MAXNG), FR_INT(MAXNG)

      LOGICAL HAS_GE

      HAS_GE = .FALSE.

C check if file is formatted as molden file
      READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
      IF (INDEX(LINE,'MOLDEN FORMAT') .EQ. 0) RETURN

      DO WHILE(IO .EQ. 0)
        READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
        IF (IO .NE. 0) RETURN

        IF (INDEX(LINE,'N_GEO') .NE. 0) THEN
          READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
          READ(LINE,*) NGEOM
        ELSE IF (INDEX(LINE,'N_ATOMS') .NE. 0) THEN
          READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
          READ(LINE,*) NATOM
        ELSE IF (INDEX(LINE,'NATOM') .NE. 0) THEN
          READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
          READ(LINE,*) NATOM
        ELSE IF (INDEX(LINE,'N_FREQ') .NE. 0) THEN
          READ(UNIT=10, FMT='(A80)', IOSTAT=IO) LINE
          READ(LINE,*) NFREQ
        ELSE IF (INDEX(LINE,'GEOCONV') .NE. 0) THEN
C read geoconv
          CALL RDGEOC(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO)
        ELSE IF (INDEX(LINE,'GEOMETRIES') .NE. 0) THEN
          IF (INDEX(LINE,'XYZ') .NE. 0) THEN
C read geometries 'CARTESIAN'
            CALL RD_GMS(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO,
     +                  0, HAS_GE)
          ELSE IF (INDEX(LINE,'AU') .NE. 0) THEN
C read geometries 'BOHR'

            CALL RD_GMS(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO,
     +                  1, HAS_GE)
          ELSE
C read geometries 'CARTESIAN'
            CALL RD_GMS(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO,
     +                  0, HAS_GE)
          END IF
        ELSE IF (INDEX(LINE,'ATOMS') .NE. 0) THEN
          IF (INDEX(LINE,'XYZ') .NE. 0) THEN
C read geometry 'CARTESIAN'
            CALL RDGEOM(MAXNG, NATOM, XYZ, 0, IO)
          ELSE IF (INDEX(LINE,'AU') .NE. 0) THEN
C read geometrY 'BOHR'
            CALL RDGEOM(MAXNG, NATOM, XYZ, 1, IO)
          ELSE
C read geometry 'CARTESIAN'
            CALL RDGEOM(MAXNG, NATOM, XYZ, 0, IO)
          END IF
        ELSE IF (INDEX(LINE,'FR-COORD') .NE. 0) THEN
C read-geom1
          CALL RDGE1(MAXNG, NATOM, XYZ, IO)
        ELSE IF (INDEX(LINE,'CHARGE') .NE. 0) THEN
          IF (INDEX(LINE,'MULLIKEN') .NE. 0) THEN
C read charge 1
            CALL RDCHRG(MAXAR, MAXNG, NATOM, XYZ, 1, IO, IMEOUT)
          ELSE
C read charge 0
            CALL RDCHRG(MAXAR, MAXNG, NATOM, XYZ, 0, IO, IMEOUT)
          END IF
        ELSE IF (INDEX(LINE,'FREQ') .NE. 0) THEN
C read array freq
          CALL RDARAY(MAXNG, NFREQ, FREQ, IO)
        ELSE IF (INDEX(LINE,'INT') .NE. 0) THEN
C read array freq-int
          CALL RDARAY(MAXNG, NFREQ, FR_INT, IO)
        ELSE IF (INDEX(LINE,'FR-NORM-COORD') .NE. 0) THEN
C read_normod
          CALL RD_NM(MAXNG, NFREQ, NATOM, FREQ, FR_INT, IO)
        END IF





      END DO

      END

      SUBROUTINE RDGEOC(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO)
      IMPLICIT NONE
      INTEGER I
      INTEGER MAXNG, NGEOM
      DOUBLE PRECISION ENERGY(MAXNG), RMS_G(MAXNG), MAX_G(MAXNG)
      CHARACTER*80 LINE
      INTEGER IO

      DO 10 I = 1, 3
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        IF (INDEX(LINE,'energy') .NE. 0) THEN
          CALL RDARAY(MAXNG, NGEOM, ENERGY, IO)
        ELSE IF (INDEX(LINE,'max-force') .NE. 0) THEN
          CALL RDARAY(MAXNG, NGEOM, RMS_G, IO)
        ELSE IF (INDEX(LINE,'rms-force') .NE. 0) THEN
          CALL RDARAY(MAXNG, NGEOM, MAX_G, IO)
        END IF
 10   CONTINUE

      END
      SUBROUTINE RD_GMS(MAXNG, NGEOM, ENERGY, RMS_G, MAX_G, IO,
     +                  IUNIT, HAS_GE)
      IMPLICIT NONE
      INTEGER I, J
      INTEGER MAXNG, NGEOM
      DOUBLE PRECISION ENERGY(MAXNG), RMS_G(MAXNG), MAX_G(MAXNG)
      CHARACTER*80 LINE
      INTEGER IO
      INTEGER IUNIT
      LOGICAL HAS_GE

      INTEGER NATOM
      DOUBLE PRECISION X, Y, Z
      DOUBLE PRECISION C_UNIT
      CHARACTER*2 ATSYM

      IF (IUNIT .EQ. 1) THEN
        C_UNIT = 0.529177249D0
      ELSE
        C_UNIT = 1.0D0
      END IF


      DO 20 I = 1, NGEOM
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        READ(LINE,*) NATOM
        WRITE(11,1000) NATOM
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        WRITE(11,*) LINE(1:LEN_TRIM(LINE))

        DO 30 J = 1, NATOM
          READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
          READ(LINE,*) ATSYM, X, Y, Z
          WRITE(11,1111) ATSYM, X * C_UNIT, Y * C_UNIT, Z * C_UNIT
 30     CONTINUE

        IF (HAS_GE) THEN
          WRITE(11,1010) ENERGY(I)
          WRITE(11,1011) RMS_G(I)
          WRITE(11,1012) MAX_G(I)
        END IF

        WRITE(11, 1100)
C        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
 20   CONTINUE

 1000 FORMAT(2X,I5)
 1010 FORMAT(1X,'<ENERGY>',/,2X,F15.8,/,1X,'</ENERGY>')
 1011 FORMAT(1X,'<RMS_G>',/,2X,F15.8,/,1X,'</RMS_G>')
 1012 FORMAT(1X,'<MAX_G>',/,2X,F15.8,/,1X,'</MAX_G>')
 1100 FORMAT(1X,'<END>')
 1111 FORMAT(1X,A6,3(2X,F15.8))
      END

      SUBROUTINE RDGEOM(MAXAR, NATOM, XYZ, IUNIT, IO)
      IMPLICIT NONE
      INTEGER I, J

      INTEGER MAXAR, NATOM
      DOUBLE PRECISION XYZ(MAXAR,3)
      INTEGER IUNIT
      INTEGER IO

      INTEGER NNUM(NATOM), ATBR
      CHARACTER*60 ATSYM(NATOM)
      DOUBLE PRECISION C_UNIT
      CHARACTER*80 LINE
      INTEGER MAXELM
      PARAMETER(MAXELM = 103)
      CHARACTER*2 ELEM(MAXELM)
      DATA ELEM(1:2) / "H ",             "He"/
      DATA ELEM(3:10) /"Li", "Be", "B ", "C ", "N ", "O ", "F", "Ne"/
      DATA ELEM(11:18)/"Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar"/
      DATA ELEM(19:36)/"K",  "Ca",
     +              "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn",
     +              "Ga", "Ge", "As", "Se", "Br", "Kr" /
      DATA ELEM(37:54) /"Rb", "Sr",
     +              "Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
     +              "In", "Sn", "Sb", "Te", "I ","Xe" /
      DATA ELEM(55:56) / "Cs", "Ba" /
      DATA ELEM(57:70)/"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
     +                 "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"/
     
      DATA ELEM(71:86)/
     +         "Lu","Hf","Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg",
     +         "Tl", "Pb", "Bi", "Po", "At", "Rn"/
      DATA ELEM(87:88) /"Fr", "Ra" /
      DATA ELEM(89:102) / "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", 
     +                    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No" /
      DATA ELEM(103:103) /"Lr"/

      WRITE(11,1000) NATOM
 1000 FORMAT(1X,I5,/)

      IF (IUNIT .EQ. 1) THEN
        C_UNIT = 0.529177249D0
      ELSE
        C_UNIT = 1.0D0
      END IF

      DO 10 I = 1, NATOM
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        READ(LINE,*) ATSYM(I), NNUM(I), ATBR, XYZ(I,1:3)
        DO 20 J = 1, 3
          XYZ(I,J) = XYZ(I,J) * C_UNIT
 20     CONTINUE
        IF (ATBR .GT. MAXELM .OR. ATBR .LE. 0) THEN
          WRITE(11,1010) "Q ", XYZ(I,1),
     +                         XYZ(I,2),
     +                         XYZ(I,3)
        ELSE
          WRITE(11,1010) ELEM(ATBR), XYZ(I,1),
     +                               XYZ(I,2),
     +                               XYZ(I,3)
        END IF
 1010   FORMAT(1X,A2,3(2X,F15.8))

        
 10   CONTINUE
      WRITE(11,1001)
 1001 FORMAT(1X,'<ATOM>')
      DO 30 I = 1, NATOM
        WRITE(11,1020) NNUM(I)
        WRITE(11,*) ATSYM(I)(1:LEN_TRIM(ATSYM(I)))
 30   CONTINUE
      WRITE(11,1011)

 1020 FORMAT(2X,'NUMBER=',I3,1X,'NAME=',$)
 1011 FORMAT(1X,'</ATOM>')
      END

      SUBROUTINE RDCHRG(MAXAR, MAXNG, NATOM, XYZ, ISMUL, IO, IMEOUT)
      IMPLICIT NONE
      INTEGER I, J

      INTEGER MAXAR, MAXNG, NATOM
      DOUBLE PRECISION XYZ(MAXNG, 3)
      INTEGER ISMUL
      INTEGER IO
      CHARACTER*360 IMEOUT

      DOUBLE PRECISION CHARGE(NATOM)
      CHARACTER*80 LINE

      WRITE(11,1000)
 1000 FORMAT(1X,'<ATOM>')

      DO 10 I = 1, NATOM
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        READ(LINE,*) CHARGE(I)
        IF (ISMUL .EQ. 1) THEN
          WRITE(11,1010) CHARGE(I)
        ELSE
          WRITE(11,1011) CHARGE(I)
        END IF
 10   CONTINUE
 1010 FORMAT(1X,'MULLIKEN_CHARGE = ', F15.8)
 1011 FORMAT(1X,'LOPROP_CHARGE = ', F15.8)

      WRITE(11,1001)
 1001 FORMAT(1X,'</ATOM>')

      CALL MKEP(MAXAR, MAXNG, NATOM, XYZ, CHARGE, IMEOUT)

      END

      SUBROUTINE MKEP(MAXAR, MAXNG, NATOM, XYZ, CHARGE, IMEOUT)
      IMPLICIT NONE
      INTEGER I, J, K, L, M

      INTEGER MAXAR, MAXNG, NATOM
      DOUBLE PRECISION XYZ(MAXNG,3)
      DOUBLE PRECISION CHARGE(NATOM)
      CHARACTER*360 IMEOUT

      DOUBLE PRECISION MINX(3), MAXX(3)
      DOUBLE PRECISION DX(3), XV(3)
      DOUBLE PRECISION CUTOFF
      PARAMETER(CUTOFF = 1.0D-5)
      INTEGER NGRID(3)
      INTEGER N_TOT
      DOUBLE PRECISION PSI(MAXAR)
      DOUBLE PRECISION R2
      INTEGER FPOS0, FPOS1
      CHARACTER*64 TMPFIL
      DOUBLE PRECISION C_UNIT
      PARAMETER(C_UNIT=1.8897261246D0)
      DOUBLE PRECISION SIG, PI
      PARAMETER(SIG=1.20D0, PI=3.1415926535897932384626D0)

      CALL MKFN(TMPFIL)

      write(*,*) 'maxar = ', maxar
      OPEN(FILE=TMPFIL, UNIT=13, FORM='UNFORMATTED', STATUS='UNKNOWN')

      DO 10 I = 1, 3
        MINX(I) = 9.99D+99
        MAXX(I) =-9.99D+99
 10   CONTINUE

      DO 20 I = 1, NATOM
        DO 30 J = 1, 3
          IF (XYZ(I,J) .LT. MINX(J)) MINX(J) = XYZ(I,J)
          IF (XYZ(I,J) .GT. MAXX(J)) MAXX(J) = XYZ(I,J)
          XYZ(I,J) = XYZ(I,J) * C_UNIT
 30     CONTINUE
 20   CONTINUE

      N_TOT = 1.0D0
      DO 40 I = 1, 3
        DX(I) = MAXX(I) - MINX(I)
        IF (DX(I) .LT. 1.0D0) THEN
          DX(I) = 1.5D0
          MINX(I) = -0.75D0
          MAXX(I) = 0.75D0
        END IF
        MAXX(I) = MAXX(I) + 2.0D0 * DX(I)
        MINX(I) = MINX(I) - 2.0D0 * DX(I)
        write(*,*) 'maxx = ', maxx(i), 'minx = ', minx(i)
        NGRID(I) = NINT((MAXX(I) - MINX(I)) * 1.5D0)
        write(*,*) 'ngrid(i) = ', ngrid
        IF (NGRID(I) .LT. 30) THEN
          NGRID(I) = 30.0D0
        END IF
        DX(I) = (MAXX(I) - MINX(I)) / DBLE(NGRID(I))
        N_TOT = N_TOT * DBLE(NGRID(I) + 1)
 40   CONTINUE
      write(*,*) 'ngrid AAA = ', ngrid

      M = 0
      DO 50 I = 1, NGRID(1)+1
        XV(1) = MINX(1) + DX(1) * DBLE(I-1)
        DO 60 J = 1, NGRID(2)+1
          XV(2) = MINX(2) + DX(2) * DBLE(J-1)
          DO 70 K = 1, NGRID(3)+1
            XV(3) = MINX(3) + DX(3) * DBLE(K-1) 
            M = M + 1
            PSI(M) = 0.0D0
            DO 80 L = 1, NATOM
              R2 = (XV(1)-XYZ(L,1))**2 +
     +             (XV(2)-XYZ(L,2))**2 +
     +             (XV(3)-XYZ(L,3))**2
              PSI(M) = PSI(M) +
     +                 EXP(-5.0D-1*R2/SIG**2)/(SIG*DSQRT(2.0D0 * PI))
 80         CONTINUE
C            if (m .lt. 31) write(*,*) m, sqrt(r2), psi(m)
 70       CONTINUE
 60     CONTINUE
        write(*,*) 'done ', i, ' / ', ngrid(1)+1
 50   CONTINUE
      CALL FTELL(13,FPOS0)
      WRITE(13) PSI(1:M)
      CALL FTELL(13,FPOS1)
      CLOSE(13)
      CALL UNLINK(TMPFIL)

      WRITE(11,1000) N_TOT, N_TOT+1
 1000 FORMAT(1X,'<GRID>',/,
     +       1X,'N_of_MO=2 N_of_Grids=2  N_of_Points=',I8,
     +       ' Block_Size=',I8,
     +       ' N_Blocks=1 Is_cutoff=0 CutOff=1.0D-5 N_P=0')
      WRITE(11,1010) NGRID(1)+1, NGRID(2)+1, NGRID(3)+1,
     +               MINX(1), MINX(2), MINX(3),
     +               MAXX(1)-MINX(1), MAXX(2)-MINX(2), MAXX(3)-MINX(3),
     +               FPOS0, FPOS1, FPOS1-FPOS0
 1010 FORMAT(1X,'N_INDEX= 0 0 0 0 0 0 0',/,' Net=',3(1X,I5),/,
     +   1X,'Origin=',3(1X,F15.8),/,
     +   1X,'Axis_1=',1X,F15.8,2('      0.00000000'),/,
     +   1X,'Axis_2=','      0.00000000',1X,F15.8,'      0.00000000',/,
     +   1X,'Axis_3=',2('      0.00000000'),1X,F15.8,/,
     +   1X,'File_pointers = ',I8,/,
     +   1X,'File_pointers = ',I8,/,
     +   1X,'ORBOFF = ',I8,/,
     +   1X,'GridName= Dummy_surface',/,
     +   1X,'GridName= Electrostatic')
      CLOSE(11)
      OPEN(UNIT=11, FILE=IMEOUT, FORM='UNFORMATTED', ACCESS='APPEND')
C      call ftell(11,fpos0)
C      write(*,*) 'fpos0 = ', fpos0
      WRITE(11) PSI(1:M)
C      call ftell(11,fpos0)
C      write(*,*) 'fpos1 = ', fpos0

      M = 0
      DO 51 I = 1, NGRID(1)+1
        XV(1) = MINX(1) + DX(1) * DBLE(I-1)
        DO 61 J = 1, NGRID(2)+1
          XV(2) = MINX(2) + DX(2) * DBLE(J-1)
          DO 71 K = 1, NGRID(3)+1
            XV(3) = MINX(3) + DX(3) * DBLE(K-1) 
            M = M + 1
            PSI(M) = 0.0D0
            DO 81 L = 1, NATOM
              R2 = (XV(1)-XYZ(L,1))**2 +
     +             (XV(2)-XYZ(L,2))**2 +
     +             (XV(3)-XYZ(L,3))**2
             PSI(M) = PSI(M) + CHARGE(L) / DSQRT(R2)
 81         CONTINUE
C            if (m .lt. 31) write(*,*) m, psi(m)
C            write(*,*) m, psi(m)
 71       CONTINUE
 61     CONTINUE
 51   CONTINUE
      WRITE(11) PSI(1:M)
C      call ftell(11,fpos0)
C      write(*,*) 'fpos2 = ', fpos0

      CLOSE(11)
      OPEN(UNIT=11, FILE=IMEOUT, FORM='FORMATTED', ACCESS='APPEND')

      WRITE(11,1020)
 1020 FORMAT(/,1X,'</GRID>')     

      END

      SUBROUTINE MKFN(TMPFIL)
      IMPLICIT NONE
      INTEGER I
      INTEGER IC
      CHARACTER*64 TMPFIL

      CALL SRAND(INT(TIME8()))
      DO 1 I = 1, 64
 9999   CONTINUE
        IC = MOD(IRAND(),74)+48
        IF (IC .LT. 65 .AND. IC .GT. 57) THEN
          GOTO 9999
        ELSE IF (IC .LT. 97 .AND. IC .GT. 90) THEN
          GOTO 9999
        END IF
        TMPFIL(I:I) = CHAR(IC)
 1    CONTINUE
      END

      SUBROUTINE RDGE1(MAXAR, NATOM, XYZ, IO)
      IMPLICIT NONE
      INTEGER I, J
      INTEGER MAXAR, NATOM
      DOUBLE PRECISION XYZ(MAXAR, 3)
      CHARACTER*6 ATSYM
      INTEGER IO

      CHARACTER*80 LINE
      DOUBLE PRECISION B2C
      PARAMETER(B2C=5.29177249D-1)

      WRITE(11,1000) NATOM

      DO 10 I = 1, NATOM
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        READ(LINE,*) ATSYM, XYZ(I,1:3)
        WRITE(11,1010) ATSYM, XYZ(I,1)*B2C, XYZ(I,2)*B2C, XYZ(I,2)*B2C
 10   CONTINUE

 1000 FORMAT(1X,I5,/)
 1010 FORMAT(1X,A2,3(1X,F15.8))

      END
      SUBROUTINE RD_NM(MAXAR, NFREQ, NATOM, FREQ, FR_INT, IO)
      IMPLICIT NONE
      INTEGER I, J

      INTEGER MAXAR, NFREQ
      INTEGER NATOM
      DOUBLE PRECISION FREQ(MAXAR), FR_INT(MAXAR)
      INTEGER IO

      DOUBLE PRECISION XYZ(3)
      CHARACTER*80 LINE

      DO 10 I = 1, NFREQ
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        WRITE(11,1000), FREQ(I), FR_INT(I)
 1000   FORMAT(1X,'<VIBRATION>',/,1X,'FREQ=',F15.8,1X,'IR_INT=',F15.8)
        DO 20 J = 1, NATOM
          READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
          READ(LINE,*) XYZ
          WRITE(11,1010) XYZ
 1010     FORMAT(3(2X,F15.8))
 20     CONTINUE
        WRITE(11,1001)
 1001   FORMAT(1X,'</VIBRATION>')
 10   CONTINUE

      END

      SUBROUTINE RDARAY(N, M, ARRAY, IO)
      IMPLICIT NONE
      INTEGER I
      INTEGER N, M
      DOUBLE PRECISION ARRAY(N)
      CHARACTER*80 LINE
      INTEGER IO

      DO 10 I = 1, M
        READ(UNIT=10, IOSTAT=IO, FMT='(A80)') LINE
        READ(LINE,*) ARRAY(I)
 10   CONTINUE

      END

