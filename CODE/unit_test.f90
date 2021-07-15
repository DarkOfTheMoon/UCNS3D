MODULE UNIT_TEST
!> @brief
!> This module provides some basic testing subroutines
USE DECLARATION
USE LIBRARY
USE BASIS
USE MPI

IMPLICIT NONE

CONTAINS

FUNCTION TEST_BASIS_REC2D(N)
    !> @brief
    !> This subroutine tests the BASIS_REC2D function
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::ORDER,NUM_DOFS
    REAL::TEST_X,TEST_Y
    REAL,ALLOCATABLE,DIMENSION(:)::BASIS_TEST
    INTEGER::TEST_BASIS_REC2D
    
    TEST_BASIS_REC2D = 0 ! No error
    
    ORDER = 2
    NUM_DOFS = 5
    
    ALLOCATE(BASIS_TEST(NUM_DOFS)) ! 5 DoFs
    TEST_X = 5
    TEST_Y = 6
        
    BASIS_TEST = BASIS_REC2d(N,TEST_X,TEST_Y,ORDER,1,NUM_DOFS)
    
    IF (BASIS_TEST(1) /= TEST_X) THEN
        WRITE(400+N,*) 'Basis term: ', 1, 'Computed value: ', BASIS_TEST(1), 'Expected value: ', TEST_X
        TEST_BASIS_REC2D = 1 ! Error
    END IF
    IF (BASIS_TEST(2) /= TEST_Y) THEN
        WRITE(400+N,*) 'Basis term: ', 2, 'Computed value: ', BASIS_TEST(2), 'Expected value: ', TEST_X
        TEST_BASIS_REC2D = 1 ! Error
    END IF
    IF (BASIS_TEST(3) /= TEST_X*TEST_X) THEN
        WRITE(400+N,*) 'Basis term: ', 3, 'Computed value: ', BASIS_TEST(3), 'Expected value: ', TEST_X
        TEST_BASIS_REC2D = 1 ! Error
    END IF
    IF (BASIS_TEST(4) /= TEST_Y*TEST_X) THEN
        WRITE(400+N,*) 'Basis term: ', 4, 'Computed value: ', BASIS_TEST(4), 'Expected value: ', TEST_X
        TEST_BASIS_REC2D = 1 ! Error
    END IF
    IF (BASIS_TEST(5) /= TEST_Y*TEST_Y) THEN
        WRITE(400+N,*) 'Basis term: ', 5, 'Computed value: ', BASIS_TEST(5), 'Expected value: ', TEST_X
        TEST_BASIS_REC2D = 1 ! Error
    END IF
    
    DEALLOCATE(BASIS_TEST)
    
END FUNCTION TEST_BASIS_REC2D

SUBROUTINE RUN_ALL_TESTS(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    
    IF (TEST_BASIS_REC2D(N) == 0) WRITE(400+N,*) 'TEST_BASIS_REC2D passed'
    
END SUBROUTINE

END MODULE
