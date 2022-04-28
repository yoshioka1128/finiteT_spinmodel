! *** MY MODIFICATIONS
! The SLATEC library has a complex error-handling system that relies on
! the XERMSG, XERPRN and other functions.  In order to make pieces of SLATEC
! easier to include in my C code (that has its own error management), I
! overwrite XERMSG so that it only only prints a simple error message.

subroutine XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)

CHARACTER*(*) LIBRAR, SUBROU, MESSG

PRINT *, "ERROR IN SLATEC '", SUBROU, "': ", MESSG

end


!!!  UPDATED D1MACH FUNCTION TAKEN FROM
!!!  http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson/f90/d1mach.f90
!!!  NB: DIFFERENT FROM THE VERSION IN BURKHARDT PORT!
!DECK D1MACH
DOUBLE PRECISION FUNCTION D1MACH (I)
IMPLICIT NONE
INTEGER :: I
DOUBLE PRECISION :: B, X
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)
!***END PROLOGUE  D1MACH
!
X = 1.0D0
B = RADIX(X)
SELECT CASE (I)
CASE (1)
D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
CASE (2)
D1MACH = HUGE(X)               ! the largest magnitude.
CASE (3)
D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
CASE (4)
D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
CASE (5)
D1MACH = LOG10(B)
CASE DEFAULT
WRITE (*, FMT = 9000)
9000     FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
STOP
END SELECT
RETURN
END


! Test function.  Used to check that arrays are passed correctly between
! C and Fortran.
subroutine set_array_to_one(pntr, ndim)

INTEGER ndim, i
DOUBLE PRECISION pntr(ndim)

DO i=1,ndim
pntr(i) = 1.
ENDDO

end subroutine



!!!  UPDATED R1MACH FUNCTION TAKEN FROM
!!!  http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson/f90/r1mach.f90
!!!  NB: DIFFERENT FROM THE VERSION IN BURKHARDT PORT!

!DECK R1MACH
REAL FUNCTION R1MACH (I)
IMPLICIT NONE
INTEGER :: I
REAL :: B, X
!***BEGIN PROLOGUE  R1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EG)
!***END PROLOGUE  R1MACH
!
X = 1.0
B = RADIX(X)
SELECT CASE (I)
CASE (1)
R1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
CASE (2)
R1MACH = HUGE(X)               ! the largest magnitude.
CASE (3)
R1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
CASE (4)
R1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
CASE (5)
R1MACH = LOG10(B)
CASE DEFAULT
WRITE (*, FMT = 9000)
9000     FORMAT ('1ERROR    1 IN R1MACH - I OUT OF BOUNDS')
STOP
END SELECT
RETURN
END
