MODULE UNIT_TEST
!> @brief
!> This module provides some basic testing subroutines
USE BASIS
USE FLUXES
USE DG_FUNCTIONS

IMPLICIT NONE

CONTAINS

FUNCTION TEST_BASIS_REC2D(N,ORDER,NUM_DOFS)
    !> @brief
    !> Tests the BASIS_REC2D function for generic polys
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS
    REAL::TEST_X,TEST_Y
    REAL,DIMENSION(NUM_DOFS)::BASIS_TEST
    INTEGER::TEST_BASIS_REC2D
    
    TEST_BASIS_REC2D = 0 ! No error
    
    TEST_X = 5
    TEST_Y = 6
        
    BASIS_TEST = BASIS_REC2d(N,TEST_X,TEST_Y,ORDER,1,NUM_DOFS)
    
    IF (POLY == 1) THEN
        IF (BASIS_TEST(1) /= TEST_X) THEN
            WRITE(400+N,*) 'Basis term: ', 1, 'Computed value: ', BASIS_TEST(1), 'Expected value: ', TEST_X
            TEST_BASIS_REC2D = 1 ! Error
        END IF
        IF (BASIS_TEST(2) /= TEST_Y) THEN
            WRITE(400+N,*) 'Basis term: ', 2, 'Computed value: ', BASIS_TEST(2), 'Expected value: ', TEST_Y
            TEST_BASIS_REC2D = 1 ! Error
        END IF
        IF (BASIS_TEST(3) /= TEST_X*TEST_X) THEN
            WRITE(400+N,*) 'Basis term: ', 3, 'Computed value: ', BASIS_TEST(3), 'Expected value: ', TEST_X*TEST_X
            TEST_BASIS_REC2D = 1 ! Error
        END IF
        IF (BASIS_TEST(4) /= TEST_Y*TEST_X) THEN
            WRITE(400+N,*) 'Basis term: ', 4, 'Computed value: ', BASIS_TEST(4), 'Expected value: ', TEST_Y*TEST_X
            TEST_BASIS_REC2D = 1 ! Error
        END IF
        IF (BASIS_TEST(5) /= TEST_Y*TEST_Y) THEN
            WRITE(400+N,*) 'Basis term: ', 5, 'Computed value: ', BASIS_TEST(5), 'Expected value: ', TEST_Y*TEST_Y
            TEST_BASIS_REC2D = 1 ! Error
        END IF
    END IF
    
END FUNCTION TEST_BASIS_REC2D

FUNCTION TEST_BASIS_REC2D_DERIVATIVE
    IMPLICIT NONE
    INTEGER::TEST_BASIS_REC2D_DERIVATIVE
    
    TEST_BASIS_REC2D_DERIVATIVE = 0

END FUNCTION TEST_BASIS_REC2D_DERIVATIVE

FUNCTION TEST_QUADRATURE(N, ORDER)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, ORDER
    REAL::TEST_QUADRATURE_ANS, TEST_QUADRATURE_ANS_ANALYTICAL
    INTEGER::I_QP, TEST_QUADRATURE
    
    TEST_QUADRATURE = 0
    
    ! Second order function
    ! Triangle
    TEST_QUADRATURE_ANS_ANALYTICAL = 1./24 ! Integral of xy on triangle (0,0),(0,1),(1,0)
    
    VEXT(1,:) = (/ 0.,0. /); VEXT(2,:) = (/ 1.,0. /); VEXT(3,:) = (/ 0.,1. /)
    CALL QUADRATURETRIANGLE(N,ORDER)
    
    TEST_QUADRATURE_ANS = 0
    DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
        TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) * QPOINTS(2,I_QP) * WEQUA3D(I_QP)
    END DO
    TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 0.5 ! 0.5 = volume of a triangle (0,0),(0,1),(1,0)
    
    IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
        WRITE(400+N,*) "Quadrature triangle expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
        TEST_QUADRATURE = 1
    END IF
    
    ! Quad
    TEST_QUADRATURE_ANS_ANALYTICAL = 1. ! Integral of xy on quadrilateral (0,0),(0,2),(1,0),(1,2)
    
    VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,2 /); VEXT(4,:) = (/ 0,2 /)
    CALL QUADRATUREQUAD(N,ORDER)
    
    TEST_QUADRATURE_ANS = 0
    DO I_QP = 1, QP_QUAD ! QP_QUAD is global
        TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) * QPOINTS(2,I_QP) * WEQUA3D(I_QP)
    END DO
    TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 2 ! 2 = volume of a quadrilateral (0,0),(0,2),(1,0),(1,2)
    
    IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
        WRITE(400+N,*) "Quadrature quadrilateral expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
        TEST_QUADRATURE = 1
    END IF
    
    ! Fourth order function
    ! Triangle
    TEST_QUADRATURE_ANS_ANALYTICAL = 1./180 ! Integral of (x^2)(y^2) on triangle (0,0),(0,1),(1,0)
    
    VEXT(1,:) = (/ 0.,0. /); VEXT(2,:) = (/ 1.,0. /); VEXT(3,:) = (/ 0.,1. /)
    CALL QUADRATURETRIANGLE(N,ORDER)
    
    TEST_QUADRATURE_ANS = 0
    DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
        TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) ** 2 * QPOINTS(2,I_QP) ** 2 * WEQUA3D(I_QP)
    END DO
    TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 0.5 ! 0.5 = volume of a triangle (0,0),(0,1),(1,0)
    
    IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
        WRITE(400+N,*) "Quadrature triangle 4th order expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
        TEST_QUADRATURE = 1
    END IF
    
    ! Quad
    TEST_QUADRATURE_ANS_ANALYTICAL = 8./9 ! Integral of (x^2)(y^2) on quadrilateral (0,0),(0,2),(1,0),(1,2)
    
    VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,2 /); VEXT(4,:) = (/ 0,2 /)
    CALL QUADRATUREQUAD(N,ORDER)
    
    TEST_QUADRATURE_ANS = 0
    DO I_QP = 1, QP_QUAD ! QP_QUAD is global
        TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) ** 2 * QPOINTS(2,I_QP) ** 2 * WEQUA3D(I_QP)
    END DO
    TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 2 ! 2 = volume of a quadrilateral (0,0),(0,2),(1,0),(1,2)
    
    IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
        WRITE(400+N,*) "Quadrature quadrilateral 4th order expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
        TEST_QUADRATURE = 1
    END IF

END FUNCTION TEST_QUADRATURE

FUNCTION TEST_DG_SOL(N,NUM_VARIABLES,ORDER,NUM_DOFS)
!> @brief
!> Tests the DG_SOL function which returns the DG solution at a given point using the generic polynomial approximation
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS,NUM_VARIABLES
    REAL::TEST_X,TEST_Y,EXPECTED_ANS
    REAL,DIMENSION(NUM_VARIABLES)::DG_SOL_TEST
    REAL,DIMENSION(NUM_VARIABLES,NUM_DOFS+1)::TEST_U_C
    INTEGER::I_VAR,TEST_DG_SOL
        
    TEST_DG_SOL = 0

    IF (POLY == 1) THEN
        TEST_X = 2
        TEST_Y = 3
        
        EXPECTED_ANS = 139
        
        DO I_VAR = 1, NUM_VARIABLES
            TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
        END DO
        
        DG_SOL_TEST = DG_SOL(N, 1, TEST_X,TEST_Y,NUM_VARIABLES,ORDER,NUM_DOFS,TEST_U_C)
        
        IF (DG_SOL_TEST(1) /= EXPECTED_ANS) THEN
            WRITE(400+N,*) 'Expected value:', EXPECTED_ANS, 'Returned value:', DG_SOL_TEST
            TEST_DG_SOL = 1
        END IF
    END IF

END FUNCTION TEST_DG_SOL

FUNCTION TEST_DG_RHS_INTEGRAL(N, NUM_VARS, ORDER, NUM_DOFS)
!> @brief
!> Tests the integral terms of the DG RHS on a triangle (0,0)(1,0)(0,1) and a quad (0,0)(1,0)(0,1)(1,1) by comparing to the analytical integral result for generic polys. Tests for 3d element shapes, other basis polys should be added
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, NUM_VARS, ORDER, NUM_DOFS
    INTEGER::I_QP, I_DOFS, I_VAR, TEST_DG_RHS_INTEGRAL !Iterator for quadrature points, DoFs, vars
    REAL::CELL_VOLUME
    REAL,DIMENSION(NUM_VARS,NUM_DOFS+1)::TEST_DG_RHS_VOL, DG_RHS_VOL_ANALYTICAL, DG_RHS_VOL_DIFF, TEST_U_C
    
    TEST_DG_RHS_INTEGRAL = 0
    
    IF (POLY == 1) THEN
    
        ! Triangle
        CELL_VOLUME = 0.5 !Volume of a triangle (0,0),(0,1),(1,0)

        VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 0,1 /)
        CALL QUADRATURETRIANGLE(N,ORDER)
        
        DO I_VAR = 1, NUM_VARS
            TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
            DG_RHS_VOL_ANALYTICAL(I_VAR,:) = (/ 0., 41./12, 41./12, 73./30, 61./24, 53./20 /)
        END DO
        
        TEST_DG_RHS_VOL = 0
        DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
            TEST_DG_RHS_VOL = TEST_DG_RHS_VOL + DG_RHS_INTEGRAL(N,1,QPOINTS(1,I_QP),QPOINTS(2,I_QP),WEQUA3D(I_QP),NUM_VARS,ORDER,NUM_DOFS,CELL_VOLUME, DG_SOL(N,1, QPOINTS(1,I_QP),QPOINTS(2,I_QP),NUM_VARS,ORDER,NUM_DOFS,TEST_U_C), 1)
        END DO


        DG_RHS_VOL_DIFF = ABS(TEST_DG_RHS_VOL - DG_RHS_VOL_ANALYTICAL)
    !     WRITE(400+N,*) 'DIFF:', DG_RHS_VOL_DIFF, '1:', TEST_DG_RHS_VOL, '2:', DG_RHS_VOL_ANALYTICAL
            
        DO I_VAR = 1, NUM_VARS
            DO I_DOFS = 1, NUM_DOFS + 1
                IF (DG_RHS_VOL_DIFF(I_VAR,I_DOFS) > 1E-12) THEN
                    WRITE(400+N,*) I_DOFS, "DG RHS VOL TRI COMPUTED VALUE: ", TEST_DG_RHS_VOL(I_VAR,I_DOFS), "EXPECTED VALUE:", DG_RHS_VOL_ANALYTICAL(I_VAR,I_DOFS)
                    TEST_DG_RHS_INTEGRAL = 1
                END IF
            END DO
        END DO
        
        ! Quad
        CELL_VOLUME = 1 !Volume of a quadrilateral (0,0),(0,1),(1,0),(1,1)

        VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,1 /); VEXT(4,:) = (/ 0,1 /)
        CALL QUADRATUREQUAD(N,ORDER)
        
        DO I_VAR = 1, NUM_VARS
            TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
            DG_RHS_VOL_ANALYTICAL(I_VAR,:) = (/ 0., 11., 11., 77./6, 157./12, 40./3 /)
        END DO
            
        TEST_DG_RHS_VOL = 0
        DO I_QP = 1, QP_QUAD ! QP_QUAD is global
            TEST_DG_RHS_VOL = TEST_DG_RHS_VOL + DG_RHS_INTEGRAL(N,1,QPOINTS(1,I_QP),QPOINTS(2,I_QP),WEQUA3D(I_QP),NUM_VARS,ORDER,NUM_DOFS,CELL_VOLUME, DG_SOL(N,1, QPOINTS(1,I_QP),QPOINTS(2,I_QP),NUM_VARS,ORDER,NUM_DOFS,TEST_U_C), 1)
        END DO
        
        DG_RHS_VOL_DIFF = ABS(TEST_DG_RHS_VOL - DG_RHS_VOL_ANALYTICAL)
    !     WRITE(400+N,*) 'DIFF:', DG_RHS_VOL_DIFF, '1:', TEST_DG_RHS_VOL, '2:', DG_RHS_VOL_ANALYTICAL
            
        DO I_VAR = 1, NUM_VARS
            DO I_DOFS = 1, NUM_DOFS + 1
                IF (DG_RHS_VOL_DIFF(I_VAR,I_DOFS) > 1E-12) THEN
                    WRITE(400+N,*) I_DOFS, "DG RHS VOL QUAD COMPUTED VALUE: ", TEST_DG_RHS_VOL(I_VAR,I_DOFS), "EXPECTED VALUE:", DG_RHS_VOL_ANALYTICAL(I_VAR,I_DOFS)
                    TEST_DG_RHS_INTEGRAL = 1
                END IF
            END DO
        END DO
    END IF

END FUNCTION TEST_DG_RHS_INTEGRAL

FUNCTION TEST_QP_ARRAY(N)
!> @brief
!> Makes sure that QP_ARRAY, for domain integral quadrature, is stored properly in the main code. Currently only for 2 dimensions, IGQRULES == 2
    INTEGER,INTENT(IN)::N
    INTEGER::I,K,I_QP,NUM_QP,TEST_QP_ARRAY
    
    TEST_QP_ARRAY = 0
    
    IF (IGQRULES == 2) THEN
        DO I = 1, XMPIELRANK(N)
            VEXT(1:IELEM(N,I)%NONODES,1) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(1)
            VEXT(1:IELEM(N,I)%NONODES,2) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(2)
            !WRITE(400+N,*) 'VEXTX:', VEXT(1:IELEM(N,I)%NONODES,1), 'VEXTY:', VEXT(1:IELEM(N,I)%NONODES,2)
            !WRITE(400+N,*) 'X:', INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(1), 'Y:', INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(2)
            IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
                CALL QUADRATUREQUAD(N,IGQRULES)
                NUM_QP = QP_QUAD
            ELSE IF (IELEM(N,I)%ISHAPE.EQ.6)THEN
                CALL QUADRATURETRIANGLE(N,IGQRULES)
                NUM_QP = QP_LINE_N
            END IF
            
            DO I_QP = 1, NUM_QP
                IF (QP_ARRAY(I,I_QP)%X /= QPOINTS(1,I_QP)) THEN
                    WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%X RETURNED VALUE: ", QP_ARRAY(1,I_QP)%X, "EXPECTED VALUE: ", QP_ARRAY(1,I_QP)%X
                    TEST_QP_ARRAY = 1
                END IF
                IF (QP_ARRAY(I,I_QP)%Y /= QPOINTS(2,I_QP)) THEN
                    WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%Y RETURNED VALUE: ", QP_ARRAY(2,I_QP)%Y, "EXPECTED VALUE: ", QP_ARRAY(1,I_QP)%Y
                    TEST_QP_ARRAY = 1
                END IF
                IF (QP_ARRAY(I,I_QP)%QP_WEIGHT /= WEQUA3D(I_QP)) THEN
                    WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%QP_WEIGHT RETURNED VALUE: ", WEQUA3D(I_QP), "EXPECTED VALUE: ", WEQUA3D(I_QP)
                    TEST_QP_ARRAY = 1
                END IF
            END DO
        END DO
    END IF

END FUNCTION TEST_QP_ARRAY

FUNCTION TEST_TAYLOR_BASIS(N, ORDER, N_DIM, N_DOFS)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, ORDER, N_DIM, N_DOFS
    REAL,DIMENSION(N_DIM)::XYZ_IN, XCYCZC_IN, DELTA_XYZ_IN
    REAL,DIMENSION(N_DOFS)::TAYLOR_BASIS_TEST, TAYLOR_BASIS_EXPECTED
    INTEGER::TEST_TAYLOR_BASIS
    
    TEST_TAYLOR_BASIS = 0
    
    VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,1 /); VEXT(4,:) = (/ 0,1 /)
    CALL QUADRATUREQUAD(N,ORDER)
    
    XYZ_IN = (/ 1., 1. /)
    XCYCZC_IN = (/ 0.5, 0.5 /)
    DELTA_XYZ_IN = CALC_DELTA_XYZ(4, N_DIM, VEXT)
    
    TAYLOR_BASIS_TEST = TAYLOR_BASIS(XYZ_IN, XCYCZC_IN, DELTA_XYZ_IN, QPOINTS, WEQUA3D, QP_QUAD, N_DIM, N_DOFS)
    TAYLOR_BASIS_EXPECTED = (/ 1., 1., 1./3, 1./3, 1. /)
    
    IF(ANY((TAYLOR_BASIS_TEST - TAYLOR_BASIS_EXPECTED) > 1.E-12)) THEN
        WRITE(400+N,*) "Taylor basis expected:", TAYLOR_BASIS_EXPECTED, "Returned value:", TAYLOR_BASIS_TEST
        TEST_TAYLOR_BASIS = 1
    END IF

END FUNCTION TEST_TAYLOR_BASIS

FUNCTION TEST_CALC_DELTA_XYZ(N_DIM)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N_DIM
    REAL,DIMENSION(3,N_DIM)::VEXT_TEST
    REAL,DIMENSION(N_DIM)::CALC_DELTA_XYZ_TEST, MIN_TEST, MAX_TEST
    INTEGER::I_DIM, TEST_CALC_DELTA_XYZ
    
    TEST_CALC_DELTA_XYZ = 0
    
    MIN_TEST(1) = 0.0015435
    MAX_TEST(1) = 65.54154
    MIN_TEST(2) = 0.21573
    MAX_TEST(2) = 789.45613
    
    VEXT_TEST(1,:) = (/ 0.15463, MAX_TEST(2) /)
    VEXT_TEST(2,:) = (/ MIN_TEST(1), 456.1783 /)
    VEXT_TEST(3,:) = (/ MAX_TEST(1), MIN_TEST(2) /)
    CALC_DELTA_XYZ_TEST = CALC_DELTA_XYZ(3, N_DIM, VEXT_TEST)
    
    DO I_DIM = 1, N_DIM
        IF (CALC_DELTA_XYZ_TEST(I_DIM) /= 0.5 * (MAX_TEST(I_DIM) - MIN_TEST(I_DIM))) THEN
            WRITE(400+N,*) 'CALC_DELTA_XYZ SHOULD BE', 0.5 * (MAX_TEST(I_DIM) - MIN_TEST(I_DIM)), 'BUT IS', CALC_DELTA_XYZ_TEST(I_DIM)
            TEST_CALC_DELTA_XYZ = 1
        END IF
    END DO
    
END FUNCTION TEST_CALC_DELTA_XYZ

SUBROUTINE CHECK_GLOBALS(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    
    WRITE(400+N,*) 'IORDER:', IORDER, 'IDEGFREE:', IDEGFREE, 'POLY:', POLY
    
    IF (DIMENSIONA.EQ.2) THEN
        IF (NUM_DG_DOFS /= (IORDER+1) * (IORDER+2) / 2) WRITE(400+N,*) 'WARNING:NUM_DG_DOFS SHOULD BE', (IORDER+1) * (IORDER+2) / 2, 'BUT IS', NUM_DG_DOFS, 'IORDER:', IORDER
    END IF
    
    IF (IORDER == 2 .AND. QP_TRIANGLE /= 3) WRITE(400+N,*) "WARNING: QP_TRIANGLE SHOULD BE", 3, "BUT IS", QP_TRIANGLE
    IF (IORDER == 2 .AND. QP_QUAD /= 4) WRITE(400+N,*) "WARNING: QP_QUAD SHOULD BE", 4, "BUT IS", QP_QUAD

END SUBROUTINE CHECK_GLOBALS

SUBROUTINE PRINT_SOME_INFO(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM
    REAL,ALLOCATABLE,DIMENSION(:,:)::TEST_VEXT
    
    ALLOCATE(TEST_VEXT(4*XMPIELRANK(N),2))
    
    DO I_ELEM = 1, XMPIELRANK(N)
        TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,1) = INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(1)
        TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,2) = INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(2)
        WRITE(400+N,*) I_ELEM, 'VEXTX:', TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,1), 'VEXTY:', TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,2)
        !WRITE(400+N,*) 'X:', INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(1), 'Y:', INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(2)
    END DO
    
    DEALLOCATE(TEST_VEXT)
    
END SUBROUTINE PRINT_SOME_INFO

SUBROUTINE RUN_ALL_TESTS(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::NUM_VARS,ORDER,NUM_DOFS,N_DIM
    
    NUM_VARS = 1
    N_DIM = 2
    ORDER = IORDER
    
    IF (N_DIM == 2) THEN
        NUM_DOFS = ((ORDER+1) * (ORDER+2) / 2) - 1
    ELSE
        NUM_DOFS = ((ORDER+1)*(ORDER+2)*(ORDER+3)/6) - 1
    END IF
    
    CALL CHECK_GLOBALS(N)
    
    IF (TEST_BASIS_REC2D(N,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_BASIS_REC2D passed'
    IF (TEST_QP_ARRAY(N) == 0) WRITE(400+N,*) 'TEST_QP_ARRAY passed'
    IF (TEST_QUADRATURE(N,ORDER) == 0) WRITE(400+N,*) 'TEST_QUADRATURE passed'
    
    IF (POLY == 1) THEN
        IF (TEST_DG_SOL(N,NUM_VARS,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_DG_SOL passed'
        IF (TEST_DG_RHS_INTEGRAL(N,NUM_VARS,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_DG_RHS_INTEGRAL passed'
    END IF
    
    IF (TEST_CALC_DELTA_XYZ(N_DIM) == 0) WRITE(400+N,*) 'TEST_CALC_DELTA_XYZ passed'
    IF (TEST_TAYLOR_BASIS(N,ORDER,N_DIM,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_TAYLOR_BASIS passed'
    
!     CALL PRINT_SOME_INFO(N)
    
END SUBROUTINE

END MODULE
