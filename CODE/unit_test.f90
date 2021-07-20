MODULE UNIT_TEST
!> @brief
!> This module provides some basic testing subroutines
USE DECLARATION
USE LIBRARY
USE TRANSFORM
USE BASIS
USE FLUXES
USE MPI

IMPLICIT NONE

CONTAINS

FUNCTION TEST_BASIS_REC2D(N,ORDER,NUM_DOFS)
    !> @brief
    !> Tests the BASIS_REC2D function
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS
    REAL::TEST_X,TEST_Y
    REAL,DIMENSION(NUM_DOFS)::BASIS_TEST
    INTEGER::TEST_BASIS_REC2D
    
    TEST_BASIS_REC2D = 0 ! No error
    
    TEST_X = 5
    TEST_Y = 6
        
    BASIS_TEST = BASIS_REC2d(N,TEST_X,TEST_Y,ORDER,1,NUM_DOFS)
    
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
    
END FUNCTION TEST_BASIS_REC2D

FUNCTION TEST_BASIS_REC2D_DERIVATIVE
    IMPLICIT NONE
    INTEGER::TEST_BASIS_REC2D_DERIVATIVE
    
    TEST_BASIS_REC2D_DERIVATIVE = 0

END FUNCTION TEST_BASIS_REC2D_DERIVATIVE

FUNCTION TEST_DG_SOL(N,ORDER,NUM_DOFS)
!> @brief
!> Tests the DG_SOL function which returns the DG solution at a given point using the polynomial approximation
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS
    REAL::TEST_X,TEST_Y,DG_SOL_TEST
    REAL,DIMENSION(NUM_DOFS+1)::TEST_U_C
    INTEGER::TEST_DG_SOL
        
    TEST_DG_SOL = 0

    TEST_X = 2
    TEST_Y = 3
    
    TEST_U_C = (/ 2, 3, 4, 5, 6, 7 /)
    
    DG_SOL_TEST = DG_SOL(N,TEST_X,TEST_Y,ORDER,NUM_DOFS,TEST_U_C)
    
    IF (DG_SOL_TEST /= 139) THEN
        WRITE(400+N,*) 'Expected value: 139, Returned value:', DG_SOL_TEST
        TEST_DG_SOL = 1
    END IF

END FUNCTION TEST_DG_SOL

FUNCTION TEST_DG_RHS_VOL_INTEGRAL(N, ORDER, NUM_DOFS)
!> @brief
!> Tests the volume integral term of the DG RHS on a triangle (0,0)(1,0)(0,1) by comparing to the analytical integral result. Tests for other element shapes need to be added
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS
    INTEGER::I_QP !Iterator for quadrature points
    REAL::TEST_DG_RHS_VOL,CELL_VOLUME
	REAL,DIMENSION(3)::WEIGHTS_TEMP_TRI !Quadrature weights for domains
	REAL,DIMENSION(2,3)::QPOINTS_TRI !Quadrature points for domains
    REAL,DIMENSION(6)::TEST_U_C !2nd order solution value and dofs
    INTEGER::TEST_DG_RHS_VOL_INTEGRAL
    
    TEST_DG_RHS_VOL_INTEGRAL = 0
    
    CELL_VOLUME = 0.5 !Volume of a triangle (0,0),(0,1),(1,0)

    VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 0,1 /)
    CALL QUADRATURETRIANGLE(N,MIN(ORDER,6))
    WEIGHTS_TEMP_TRI = WEQUA3D(1:QP_TRIANGLE)
    QPOINTS_TRI = QPOINTS(1:2,1:QP_TRIANGLE)

    TEST_U_C = (/ 2, 3, 4, 5, 6, 7 /)
    
    DO I_QP = 1, QP_TRIANGLE   
        TEST_DG_RHS_VOL = TEST_DG_RHS_VOL+ DG_RHS_VOL_INTEGRAL(N,QPOINTS_TRI(1,I_QP),QPOINTS_TRI(2,I_QP),WEIGHTS_TEMP_TRI(I_QP),ORDER,NUM_DOFS,CELL_VOLUME,TEST_U_C)
    END DO
    
    IF (TEST_DG_RHS_VOL - 14.4583333333334 > 1E-12) THEN    
        WRITE(400+N,*) "DG RHS VOL COMPUTED VALUE: ", TEST_DG_RHS_VOL, "EXPECTED VALUE:", 14.4583333333334
        TEST_DG_RHS_VOL_INTEGRAL = 1
    END IF

END FUNCTION TEST_DG_RHS_VOL_INTEGRAL

FUNCTION TEST_QP_ARRAY(N)
!> @brief
!> Makes sure that QP_ARRAY is stored properly in the main code. Currently only for 2 dimensions, IGQRULES == 2
    INTEGER,INTENT(IN)::N
    INTEGER::I,K,I_QP,NUM_QP,KMAXE,TEST_QP_ARRAY
    
    TEST_QP_ARRAY = 0
    
    KMAXE = XMPIELRANK(N)
    
    IF (IGQRULES == 2) THEN
        DO I = 1, KMAXE
            VEXT(1:IELEM(N,I)%NONODES,1) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(1)
            VEXT(1:IELEM(N,I)%NONODES,2) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(2)
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

FUNCTION TEST_DG_RHS_SURF_INTEGRAL(N,ORDER)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER
    INTEGER::TEST_DG_RHS_SURF_INTEGRAL
    INTEGER::I_QP !Iterator for quadrature points
    REAL::TEST_DG_RHS_ANS,CELL_SURF
	REAL,DIMENSION(3)::WEIGHTS_TEMP_LINE !Quadrature weights for interfaces
	REAL,DIMENSION(2,3)::QPOINTS_LINE!Quadrature points for interfaces
    REAL,DIMENSION(6)::TEST_U_C !2nd order solution value and dofs
    
    TEST_DG_RHS_SURF_INTEGRAL = 0
    
    CELL_SURF = 0 !Volume of a triangle (0,0),(0,1),(1,0)

    VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 0,1 /)
    CALL QUADRATURELINE(N,MIN(ORDER,6))
    WEIGHTS_TEMP_LINE = WEQUA2D(1:QP_TRIANGLE)
    QPOINTS_LINE = QPOINTS2D(1:2,1:QP_LINE_N)

    TEST_U_C = (/ 2, 3, 4, 5, 6, 7 /)
    
    IF (TEST_DG_RHS_ANS - 14.4583333333334 > 1E-12) THEN    
        WRITE(400+N,*) "DG RHS SURFACE INTEGRAL COMPUTED VALUE: ", TEST_DG_RHS_ANS, "EXPECTED VALUE:", 14.4583333333334
        TEST_DG_RHS_SURF_INTEGRAL = 1
    END IF
    
END FUNCTION TEST_DG_RHS_SURF_INTEGRAL

SUBROUTINE RUN_ALL_TESTS(N,ORDER,NUM_DOFS)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS
    INTEGER::N_DIM
    
    N_DIM = 2
    
    IF (TEST_BASIS_REC2D(N,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_BASIS_REC2D passed'
    IF (TEST_DG_SOL(N,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_DG_SOL passed'
    IF (TEST_DG_RHS_VOL_INTEGRAL(N,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_DG_RHS_VOL_INTEGRAL passed'
    IF (TEST_QP_ARRAY(N) == 0) WRITE(400+N,*) 'TEST_QP_ARRAY passed'
    
END SUBROUTINE

END MODULE
