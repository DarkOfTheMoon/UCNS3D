! MODULE UNIT_TEST
! !> @brief
! !> This module provides some basic testing subroutines
! USE BASIS
! USE FLUXES
! USE DG_FUNCTIONS
! 
! IMPLICIT NONE
! 
! CONTAINS
! 
! SUBROUTINE TEST_BASIS_REC2D(N,ORDER,NUM_DOFS, ERROR_FILE_INDEX)
! !> @brief
! !> Tests the BASIS_REC2D function, currently only for generic polys
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N,ORDER,NUM_DOFS, ERROR_FILE_INDEX
!     REAL::TEST_X,TEST_Y
!     REAL,DIMENSION(NUM_DOFS)::BASIS_TEST
!     INTEGER::ERROR_YES
!     
!     ERROR_YES = 0
!     
!     TEST_X = 5
!     TEST_Y = 6
!         
!     BASIS_TEST = BASIS_REC2d(N,TEST_X,TEST_Y,ORDER,1,NUM_DOFS)
!     
!     IF (POLY == 1 .AND. ORDER > 1) THEN
!         IF (BASIS_TEST(1) /= TEST_X) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Basis term: ', 1, 'Computed value: ', BASIS_TEST(1), 'Expected value: ', TEST_X
!             ERROR_YES = 1 ! Error
!         END IF
!         IF (BASIS_TEST(2) /= TEST_Y) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Basis term: ', 2, 'Computed value: ', BASIS_TEST(2), 'Expected value: ', TEST_Y
!             ERROR_YES = 1 ! Error
!         END IF
!         IF (BASIS_TEST(3) /= TEST_X*TEST_X) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Basis term: ', 3, 'Computed value: ', BASIS_TEST(3), 'Expected value: ', TEST_X*TEST_X
!             ERROR_YES = 1 ! Error
!         END IF
!         IF (BASIS_TEST(4) /= TEST_Y*TEST_X) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Basis term: ', 4, 'Computed value: ', BASIS_TEST(4), 'Expected value: ', TEST_Y*TEST_X
!             ERROR_YES = 1 ! Error
!         END IF
!         IF (BASIS_TEST(5) /= TEST_Y*TEST_Y) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Basis term: ', 5, 'Computed value: ', BASIS_TEST(5), 'Expected value: ', TEST_Y*TEST_Y
!             ERROR_YES = 1 ! Error
!         END IF
!         IF (ERROR_YES == 0) WRITE(ERROR_FILE_INDEX,*) 'TEST_BASIS_REC2D passed'
!     ELSE
!         WRITE(ERROR_FILE_INDEX,*) 'TEST_BASIS_REC2D did not run'
!     END IF
!     
! END SUBROUTINE TEST_BASIS_REC2D
! 
! FUNCTION TEST_BASIS_REC2D_DERIVATIVE
!     IMPLICIT NONE
!     INTEGER::TEST_BASIS_REC2D_DERIVATIVE
!     
!     TEST_BASIS_REC2D_DERIVATIVE = 0
! 
! END FUNCTION TEST_BASIS_REC2D_DERIVATIVE
! 
! FUNCTION TEST_QUADRATURE(N, ORDER)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, ORDER
!     REAL::TEST_QUADRATURE_ANS, TEST_QUADRATURE_ANS_ANALYTICAL
!     INTEGER::I_QP, TEST_QUADRATURE
!     
!     TEST_QUADRATURE = 0
!     
!     ! Second order function
!     ! Triangle
!     TEST_QUADRATURE_ANS_ANALYTICAL = 65./24 !1./24 ! Integral of xy on triangle (0,0),(0,1),(1,0)
!     
!     VEXT(1,:) = (/ 2.,2. /); VEXT(2,:) = (/ 3.,2. /); VEXT(3,:) = (/ 2.,3. /)
!     CALL QUADRATURETRIANGLE(N,ORDER)
!     
!     TEST_QUADRATURE_ANS = 0
!     DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
!         TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) * QPOINTS(2,I_QP) * WEQUA3D(I_QP)
!     END DO
!     TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 0.5 ! 0.5 = volume of a triangle (0,0),(0,1),(1,0)
!     
!      IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
!         WRITE(400+N,*) "Quadrature triangle expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
!         TEST_QUADRATURE = 1
!     END IF
!     
!     ! Quad
!     TEST_QUADRATURE_ANS_ANALYTICAL = 1. ! Integral of xy on quadrilateral (0,0),(0,2),(1,0),(1,2)
!     
!     VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,2 /); VEXT(4,:) = (/ 0,2 /)
!     CALL QUADRATUREQUAD(N,ORDER)
!     
!     TEST_QUADRATURE_ANS = 0
!     DO I_QP = 1, QP_QUAD ! QP_QUAD is global
!         TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) * QPOINTS(2,I_QP) * WEQUA3D(I_QP)
!     END DO
!     TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 2 ! 2 = volume of a quadrilateral (0,0),(0,2),(1,0),(1,2)
!     
!     IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
!         WRITE(400+N,*) "Quadrature quadrilateral expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
!         TEST_QUADRATURE = 1
!     END IF
!     
!     ! Fourth order function
!     ! Triangle
!     TEST_QUADRATURE_ANS_ANALYTICAL = 1./180 ! Integral of (x^2)(y^2) on triangle (0,0),(0,1),(1,0)
!     
!     VEXT(1,:) = (/ 0.,0. /); VEXT(2,:) = (/ 1.,0. /); VEXT(3,:) = (/ 0.,1. /)
!     CALL QUADRATURETRIANGLE(N,ORDER)
!     
!     TEST_QUADRATURE_ANS = 0
!     DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
!         TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) ** 2 * QPOINTS(2,I_QP) ** 2 * WEQUA3D(I_QP)
!     END DO
!     TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 0.5 ! 0.5 = volume of a triangle (0,0),(0,1),(1,0)
!     
!     IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
!         WRITE(400+N,*) "Quadrature triangle 4th order expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
!         TEST_QUADRATURE = 1
!     END IF
!     
!     ! Quad
!     TEST_QUADRATURE_ANS_ANALYTICAL = 8./9 ! Integral of (x^2)(y^2) on quadrilateral (0,0),(0,2),(1,0),(1,2)
!     
!     VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,2 /); VEXT(4,:) = (/ 0,2 /)
!     CALL QUADRATUREQUAD(N,ORDER)
!     
!     TEST_QUADRATURE_ANS = 0
!     DO I_QP = 1, QP_QUAD ! QP_QUAD is global
!         TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS + QPOINTS(1,I_QP) ** 2 * QPOINTS(2,I_QP) ** 2 * WEQUA3D(I_QP)
!     END DO
!     TEST_QUADRATURE_ANS = TEST_QUADRATURE_ANS * 2 ! 2 = volume of a quadrilateral (0,0),(0,2),(1,0),(1,2)
!     
!     IF (ABS(TEST_QUADRATURE_ANS - TEST_QUADRATURE_ANS_ANALYTICAL) > 1.E-12) THEN
!         WRITE(400+N,*) "Quadrature quadrilateral 4th order expected value: ", TEST_QUADRATURE_ANS_ANALYTICAL, "Calculated value: ", TEST_QUADRATURE_ANS
!         TEST_QUADRATURE = 1
!     END IF
! 
! END FUNCTION TEST_QUADRATURE
! 
! SUBROUTINE TEST_DG_SOL(N, NUM_VARIABLES, ORDER, NUM_DOFS, ERROR_FILE_INDEX)
! !> @brief
! !> Tests the DG_SOL function which returns the DG solution at a given point using the generic polynomial approximation
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, ORDER, NUM_DOFS, NUM_VARIABLES, ERROR_FILE_INDEX
!     REAL::TEST_X, TEST_Y, EXPECTED_ANS
!     REAL,DIMENSION(NUM_VARIABLES,1)::DG_SOL_TEST
!     REAL,DIMENSION(NUM_VARIABLES,NUM_DOFS+1)::TEST_U_C
!     INTEGER::I_VAR
! 
!     IF (POLY == 1 .AND. ORDER > 1) THEN
!         TEST_X = 2
!         TEST_Y = 3
!         
!         EXPECTED_ANS = 139
!         
!         DO I_VAR = 1, NUM_VARIABLES
!             TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
!         END DO
!         
!         DG_SOL_TEST(:,1) = DG_SOL(N, 1, TEST_X,TEST_Y,NUM_VARIABLES,ORDER,NUM_DOFS,TEST_U_C) - DG_SOL(N, 1, TEST_X+1,TEST_Y+1,NUM_VARIABLES,ORDER,NUM_DOFS,TEST_U_C)
!         
!         WRITE(ERROR_FILE_INDEX,*) 'SHAPE(DG_SOL_TEST)', SHAPE(DG_SOL_TEST)
!         
!         IF (DG_SOL_TEST(1,1) /= EXPECTED_ANS) THEN
!             WRITE(ERROR_FILE_INDEX,*) 'Expected value:', EXPECTED_ANS, 'Returned value:', DG_SOL_TEST
!         ELSE
!             WRITE(ERROR_FILE_INDEX,*) 'TEST_DG_SOL passed'
!         END IF
!     END IF
! 
! END SUBROUTINE TEST_DG_SOL
! 
! FUNCTION TEST_DG_RHS_INTEGRAL(N, NUM_VARS, ORDER, NUM_DOFS)
! !> @brief
! !> Tests the integral terms of the DG RHS on a triangle (0,0)(1,0)(0,1) and a quad (0,0)(1,0)(0,1)(1,1) by comparing to the analytical integral result for generic polys. Tests for 3d element shapes, other basis polys should be added
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, NUM_VARS, ORDER, NUM_DOFS
!     INTEGER::I_QP, I_DOFS, I_VAR, TEST_DG_RHS_INTEGRAL !Iterator for quadrature points, DoFs, vars
!     REAL::CELL_VOLUME
!     REAL,DIMENSION(2)::CELL_CENTER
!     REAL,DIMENSION(NUM_VARS,NUM_DOFS+1)::TEST_DG_RHS_VOL, DG_RHS_VOL_ANALYTICAL, DG_RHS_VOL_DIFF, TEST_U_C
!     
!     TEST_DG_RHS_INTEGRAL = 0
!     
!     IF (POLY == 1 .AND. ORDER > 2) THEN
!     
!         ! Triangle
!         CELL_VOLUME = 0.5 !Volume of a triangle (0,0),(0,1),(1,0)
! 
!         VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 0,1 /)
!         CALL QUADRATURETRIANGLE(N,ORDER)
!         
!         CELL_CENTER = (/ 1./3, 1./3 /)
!         
!         DO I_VAR = 1, NUM_VARS
!             TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
!             DG_RHS_VOL_ANALYTICAL(I_VAR,:) = (/ 0., 41./12, 41./12, 73./30, 61./24, 53./20 /)
!         END DO
!         
!         TEST_DG_RHS_VOL = 0
!         DO I_QP = 1, QP_TRIANGLE ! QP_TRIANGLE is global
!             QPOINTS(:,I_QP) = QPOINTS(:,I_QP) - CELL_CENTER
!             
!             TEST_DG_RHS_VOL = TEST_DG_RHS_VOL + DG_RHS_INTEGRAL(N,1,QPOINTS(1,I_QP),QPOINTS(2,I_QP),WEQUA3D(I_QP),NUM_VARS,ORDER,NUM_DOFS,CELL_VOLUME, DG_SOL(N,1, QPOINTS(1,I_QP),QPOINTS(2,I_QP),NUM_VARS,ORDER,NUM_DOFS,TEST_U_C), 1)
!         END DO
! 
! 
!         DG_RHS_VOL_DIFF = ABS(TEST_DG_RHS_VOL - DG_RHS_VOL_ANALYTICAL)
!     !     WRITE(400+N,*) 'DIFF:', DG_RHS_VOL_DIFF, '1:', TEST_DG_RHS_VOL, '2:', DG_RHS_VOL_ANALYTICAL
!             
!         DO I_VAR = 1, NUM_VARS
!             DO I_DOFS = 1, NUM_DOFS + 1
!                 IF (DG_RHS_VOL_DIFF(I_VAR,I_DOFS) > 1E-12) THEN
!                     WRITE(400+N,*) I_DOFS, "DG RHS VOL TRI COMPUTED VALUE: ", TEST_DG_RHS_VOL(I_VAR,I_DOFS), "EXPECTED VALUE:", DG_RHS_VOL_ANALYTICAL(I_VAR,I_DOFS)
!                     TEST_DG_RHS_INTEGRAL = 1
!                 END IF
!             END DO
!         END DO
!         
!         ! Quad
!         CELL_VOLUME = 1 !Volume of a quadrilateral (0,0),(0,1),(1,0),(1,1)
! 
!         VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,1 /); VEXT(4,:) = (/ 0,1 /)
!         CALL QUADRATUREQUAD(N,ORDER)
!         
!         CELL_CENTER = (/ 1./2, 1./2 /)
!         DO I_VAR = 1, NUM_VARS
!             TEST_U_C(I_VAR,:) = (/ 2, 3, 4, 5, 6, 7 /)
!             DG_RHS_VOL_ANALYTICAL(I_VAR,:) = (/ 0., 11., 11., 77./6, 157./12, 40./3 /)
!         END DO
!             
!         TEST_DG_RHS_VOL = 0
!         DO I_QP = 1, QP_QUAD ! QP_QUAD is global
!             QPOINTS(:,I_QP) = QPOINTS(:,I_QP) - CELL_CENTER
!             
!             TEST_DG_RHS_VOL = TEST_DG_RHS_VOL + DG_RHS_INTEGRAL(N,1,QPOINTS(1,I_QP),QPOINTS(2,I_QP),WEQUA3D(I_QP),NUM_VARS,ORDER,NUM_DOFS,CELL_VOLUME, DG_SOL(N,1, QPOINTS(1,I_QP),QPOINTS(2,I_QP),NUM_VARS,ORDER,NUM_DOFS,TEST_U_C), 1)
!         END DO
!         
!         DG_RHS_VOL_DIFF = ABS(TEST_DG_RHS_VOL - DG_RHS_VOL_ANALYTICAL)
!     !     WRITE(400+N,*) 'DIFF:', DG_RHS_VOL_DIFF, '1:', TEST_DG_RHS_VOL, '2:', DG_RHS_VOL_ANALYTICAL
!             
!         DO I_VAR = 1, NUM_VARS
!             DO I_DOFS = 1, NUM_DOFS + 1
!                 IF (DG_RHS_VOL_DIFF(I_VAR,I_DOFS) > 1E-12) THEN
!                     WRITE(400+N,*) I_DOFS, "DG RHS VOL QUAD COMPUTED VALUE: ", TEST_DG_RHS_VOL(I_VAR,I_DOFS), "EXPECTED VALUE:", DG_RHS_VOL_ANALYTICAL(I_VAR,I_DOFS)
!                     TEST_DG_RHS_INTEGRAL = 1
!                 END IF
!             END DO
!         END DO
!     END IF
! 
! END FUNCTION TEST_DG_RHS_INTEGRAL
! 
! FUNCTION TEST_QP_ARRAY(N)
! !> @brief
! !> Makes sure that QP_ARRAY, for domain integral quadrature, is stored properly in the main code. Currently only for 2 dimensions, IGQRULES == 2
!     INTEGER,INTENT(IN)::N
!     INTEGER::I,K,I_QP,NUM_QP,TEST_QP_ARRAY
!     
!     TEST_QP_ARRAY = 0
!     
! !     IF (IGQRULES == 2) THEN
!         DO I = 1, XMPIELRANK(N)
!             VEXT(1:IELEM(N,I)%NONODES,1) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(1)
!             VEXT(1:IELEM(N,I)%NONODES,2) = INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(2)
!             !WRITE(400+N,*) 'VEXTX:', VEXT(1:IELEM(N,I)%NONODES,1), 'VEXTY:', VEXT(1:IELEM(N,I)%NONODES,2)
!             !WRITE(400+N,*) 'X:', INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(1), 'Y:', INODER(IELEM(N,I)%NODES(1:IELEM(N,I)%NONODES))%CORD(2)
!             IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
!                 CALL QUADRATUREQUAD(N,IGQRULES)
!                 NUM_QP = QP_QUAD
!             ELSE IF (IELEM(N,I)%ISHAPE.EQ.6)THEN
!                 CALL QUADRATURETRIANGLE(N,IGQRULES)
!                 NUM_QP = QP_TRIANGLE
!             END IF
!             
!             DO I_QP = 1, NUM_QP
!                 IF (QP_ARRAY(I,I_QP)%X - QPOINTS(1,I_QP) > 1E-12) THEN
!                     WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%X RETURNED VALUE: ", QP_ARRAY(I,I_QP)%X, "EXPECTED VALUE: ", QPOINTS(1,I_QP) - IELEM(N,I)%XXC
!                     TEST_QP_ARRAY = 1
!                 END IF
!                 IF (QP_ARRAY(I,I_QP)%Y /= QPOINTS(2,I_QP) > 1E-12) THEN
!                     WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%Y RETURNED VALUE: ", QP_ARRAY(I,I_QP)%Y, "EXPECTED VALUE: ", QPOINTS(2,I_QP) - IELEM(N,I)%YYC
!                     TEST_QP_ARRAY = 1
!                 END IF
!                 IF (QP_ARRAY(I,I_QP)%QP_WEIGHT /= WEQUA3D(I_QP) > 1E-12) THEN
!                     WRITE(400+N,*) "QP_ARRAY(",I,",",I_QP,")%QP_WEIGHT RETURNED VALUE: ", QP_ARRAY(I,I_QP)%QP_WEIGHT, "EXPECTED VALUE: ", WEQUA3D(I_QP)
! 
!                     TEST_QP_ARRAY = 1
!                 END IF
!             END DO
!         END DO
! 
! END FUNCTION TEST_QP_ARRAY
! 
! FUNCTION TEST_TAYLOR_BASIS(N, ORDER, N_DIM, N_DOFS)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, ORDER, N_DIM, N_DOFS
!     REAL::CELL_VOLUME
!     REAL,DIMENSION(N_DIM)::XYZ_IN, XCYCZC_IN, DELTA_XYZ_IN
!     REAL,DIMENSION(N_DOFS)::TAYLOR_BASIS_TEST, TAYLOR_BASIS_EXPECTED
!     INTEGER::I_QP,TEST_TAYLOR_BASIS
!     
!     TEST_TAYLOR_BASIS = 0
!     
!     VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 1,1 /); VEXT(4,:) = (/ 0,1 /)
!     CALL QUADRATUREQUAD(N,ORDER)
!     
!     XYZ_IN = (/ 1., 1. /)
!     XCYCZC_IN = (/ 0.5, 0.5 /)
!     CELL_VOLUME = 1.
!     DELTA_XYZ_IN = CALC_DELTA_XYZ(4, N_DIM, VEXT)
!     
!     DO I_QP = 1, QP_QUAD
!         QPOINTS(:,I_QP) = QPOINTS(:,I_QP) - XCYCZC_IN
!     END DO
!     
!     TAYLOR_BASIS_TEST = TAYLOR_BASIS(XYZ_IN - XCYCZC_IN, DELTA_XYZ_IN, QPOINTS, WEQUA3D, QP_QUAD, N_DIM, N_DOFS, CELL_VOLUME)
!     TAYLOR_BASIS_EXPECTED = (/ 1., 1., 1./3, 1./3, 1. /)
!     
!     IF(ANY((TAYLOR_BASIS_TEST - TAYLOR_BASIS_EXPECTED) > 1.E-12)) THEN
!         WRITE(400+N,*) "Taylor basis expected:", TAYLOR_BASIS_EXPECTED, "Returned value:", TAYLOR_BASIS_TEST
!         TEST_TAYLOR_BASIS = 1
!     END IF
! 
! END FUNCTION TEST_TAYLOR_BASIS
! 
! FUNCTION TEST_CALC_DELTA_XYZ(N_DIM)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N_DIM
!     REAL,DIMENSION(3,N_DIM)::VEXT_TEST
!     REAL,DIMENSION(N_DIM)::CALC_DELTA_XYZ_TEST, MIN_TEST, MAX_TEST
!     INTEGER::I_DIM, TEST_CALC_DELTA_XYZ
!     
!     TEST_CALC_DELTA_XYZ = 0
!     
!     MIN_TEST(1) = 0.0015435
!     MAX_TEST(1) = 65.54154
!     MIN_TEST(2) = 0.21573
!     MAX_TEST(2) = 789.45613
!     
!     VEXT_TEST(1,:) = (/ 0.15463, MAX_TEST(2) /)
!     VEXT_TEST(2,:) = (/ MIN_TEST(1), 456.1783 /)
!     VEXT_TEST(3,:) = (/ MAX_TEST(1), MIN_TEST(2) /)
!     CALC_DELTA_XYZ_TEST = CALC_DELTA_XYZ(3, N_DIM, VEXT_TEST)
!     
!     DO I_DIM = 1, N_DIM
!         IF (CALC_DELTA_XYZ_TEST(I_DIM) /= 0.5 * (MAX_TEST(I_DIM) - MIN_TEST(I_DIM))) THEN
!             WRITE(400+N,*) 'CALC_DELTA_XYZ SHOULD BE', 0.5 * (MAX_TEST(I_DIM) - MIN_TEST(I_DIM)), 'BUT IS', CALC_DELTA_XYZ_TEST(I_DIM)
!             TEST_CALC_DELTA_XYZ = 1
!         END IF
!     END DO
!     
! END FUNCTION TEST_CALC_DELTA_XYZ
! 
! SUBROUTINE TEST_NORMAL_LSQ(TEST_A, TEST_b, TEST_ANS, ERROR_FILE_INDEX)
!     REAL,DIMENSION(:,:),INTENT(IN)::TEST_A
!     REAL,DIMENSION(:),INTENT(IN)::TEST_b, TEST_ANS
!     INTEGER,INTENT(IN)::ERROR_FILE_INDEX
!     REAL,DIMENSION(SIZE(TEST_A,2))::NORMAL_LSQ_ANS
!     
!     NORMAL_LSQ_ANS = NORMAL_LSQ(TEST_A, TEST_b)
!     
!     IF(ALL(NORMAL_LSQ_ANS == TEST_ANS)) THEN
!         WRITE(ERROR_FILE_INDEX,*) 'TEST_NORMAL_LSQ passed'
!     ELSE
!         WRITE(ERROR_FILE_INDEX,*) 'NORMAL_LSQ returned ans:', NORMAL_LSQ_ANS, 'Expected ans:', TEST_ANS
!     END IF
!     
! END SUBROUTINE TEST_NORMAL_LSQ
! 
! SUBROUTINE CHECK_GLOBALS(N, ERROR_FILE_INDEX)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, ERROR_FILE_INDEX
!     
!     WRITE(ERROR_FILE_INDEX,*) 'IORDER:', IORDER, 'IDEGFREE:', IDEGFREE, 'POLY:', POLY, 'DG:', DG
!     
!     IF (DIMENSIONA.EQ.2) THEN
!         IF (NUM_DG_DOFS /= (IORDER+1) * (IORDER+2) / 2) WRITE(ERROR_FILE_INDEX,*) 'WARNING:NUM_DG_DOFS SHOULD BE', (IORDER+1) * (IORDER+2) / 2, 'BUT IS', NUM_DG_DOFS, 'IORDER:', IORDER
!     END IF
!     
!     IF (DG == 1) THEN
!         IF (IORDER == 1 .AND. QP_TRIANGLE /= 3) WRITE(ERROR_FILE_INDEX,*) "WARNING: QP_TRIANGLE SHOULD BE", 3, "BUT IS", QP_TRIANGLE
!         IF (IORDER == 1 .AND. QP_QUAD /= 4) WRITE(ERROR_FILE_INDEX,*) "WARNING: QP_QUAD SHOULD BE", 4, "BUT IS", QP_QUAD
!         IF (IORDER == 2 .AND. QP_TRIANGLE /= 6) WRITE(ERROR_FILE_INDEX,*) "WARNING: QP_TRIANGLE SHOULD BE", 6, "BUT IS", QP_TRIANGLE
!         IF (IORDER == 2 .AND. QP_QUAD /= 9) WRITE(ERROR_FILE_INDEX,*) "WARNING: QP_QUAD SHOULD BE", 9, "BUT IS", QP_QUAD
!     END IF
!     
!     WRITE(ERROR_FILE_INDEX,*) "NUMBEROFPOINTS, NUMBEROFPOINTS2, IGQRULES, QP_LINE_N, QP_TRIANGLE, QP_QUAD, ITESTCASE, ISCHEME, NOF_VARIABLES"
!     WRITE(ERROR_FILE_INDEX,*) NUMBEROFPOINTS, NUMBEROFPOINTS2, IGQRULES, QP_LINE_N, QP_TRIANGLE, QP_QUAD, ITESTCASE, ISCHEME, NOF_VARIABLES
! 
! END SUBROUTINE CHECK_GLOBALS
! 
! SUBROUTINE PRINT_SOME_INFO(N)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N
!     INTEGER::I_ELEM
!     REAL,DIMENSION(4*XMPIELRANK(N),2)::TEST_VEXT
!     
!     WRITE(400+N,*) 'VERTICES:'
!     DO I_ELEM = 1, XMPIELRANK(N)
!         TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,1) = INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(1)
!         TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,2) = INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(2)
!         WRITE(400+N,*) I_ELEM, IELEM(N,I_ELEM)%NODES, 'VEXTX:', TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,1), 'VEXTY:', TEST_VEXT(1:IELEM(N,I_ELEM)%NONODES,2), 'XXC:', IELEM(N,I_ELEM)%XXC, 'YYC:', IELEM(N,I_ELEM)%YYC, 'TOTVOLUME:', IELEM(N,I_ELEM)%TOTVOLUME
!         !WRITE(400+N,*) 'X:', INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(1), 'Y:', INODER(IELEM(N,I_ELEM)%NODES(1:IELEM(N,I_ELEM)%NONODES))%CORD(2)
!     END DO
!     
! END SUBROUTINE PRINT_SOME_INFO
! 
! SUBROUTINE TEST_MASS_MATRIX(ERROR_FILE_INDEX)
! !> @brief
! !> This function checks that the mass matrix is invertible for a triangle (0,0) (1,0) (0,1)
! 
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::ERROR_FILE_INDEX
!     INTEGER::I_QP, N_QP, N_DIM, N_DOFS, I_DOF, J_DOF, ERROR_YES
!     REAL::CELL_VOLUME
!     REAL,DIMENSION(DIMENSIONA)::XCYC, DELTA_XY !N_DIM
! !     REAL,DIMENSION(3)::QP_WEIGHT_TEST !N_QP
! !     REAL,DIMENSION(2,6)::QP_TEST !N_DIM, N_QP
!     REAL,DIMENSION(QP_TRIANGLE,IDEGFREE)::BASIS_VECTOR_TEST !N_QP, IDEGFREE
!     REAL,DIMENSION(IDEGFREE,IDEGFREE)::MM_TEST, IDENTITY_MATRIX !IDEGFREE, IDEGFREE
!     REAL,DIMENSION(1,IDEGFREE,IDEGFREE)::INV_MM
!     REAL,DIMENSION(IDEGFREE)::MM_TEST1 !IDEGFREE
!     
!     ERROR_YES = 0
!     N_QP = QP_TRIANGLE; N_DIM = 2; N_DOFS = IDEGFREE
!     
!     IDENTITY_MATRIX = 0.
!     DO I_DOF = 1, N_DOFS
!         DO J_DOF = 1, N_DOFS
!             IF (I_DOF == J_DOF) IDENTITY_MATRIX(I_DOF, J_DOF) = 1.
!         END DO
!     END DO
!                 
! !     CELL_VOLUME = 6.560397999963079E-004
! !     XCYC = (/ 0.228697584059357, 0.420303461038695 /)
! !     QP_TEST = RESHAPE((/ 0.225734681133017, 0.429656077863782, 0.220108049013362, 0.410513961552453, 0.240250022031694, 0.420740343699850 /), SHAPE(QP_TEST))
! !     QP_WEIGHT_TEST = (/ 0.333333333333333, 0.333333333333333, 0.333333333333333 /)
! !     DELTA_XY = (/ 2.014197301833204E-002, 1.914211631132953E-002 /)
!     
!     CELL_VOLUME = 1./2
!     XCYC = (/ 1./3, 1./3 /)
! !     QP_TEST = RESHAPE((/ 2./3, 1./6, 1./6, 2./3, 1./6, 1./6 /), SHAPE(QP_TEST))
! !     write(400+N,*) QP_TEST(1,1),QP_TEST(2,1),QP_TEST(1,2),QP_TEST(2,2),QP_TEST(1,3),QP_TEST(2,3)
! !     QP_WEIGHT_TEST = (/ 1./3, 1./3, 1./3 /)
!     VEXT(1,:) = (/ 0,0 /); VEXT(2,:) = (/ 1,0 /); VEXT(3,:) = (/ 0,1 /)
!     CALL QUADRATURETRIANGLE(N,IGQRULES)
!     
!     DELTA_XY = (/ 1./2, 1./2 /)
!     
!     DO I_QP = 1, N_QP
!         BASIS_VECTOR_TEST(I_QP,:) = BASIS_REC2D(N,QPOINTS(1,I_QP),QPOINTS(2,I_QP),IORDER,1,N_DOFS)
!     END DO
!     
!     MM_TEST = 0
!     DO I_DOF = 1, N_DOFS
!         DO J_DOF = 1, N_DOFS
!             DO I_QP = 1, N_QP
!                 MM_TEST(I_DOF,J_DOF) = MM_TEST(I_DOF,J_DOF) + BASIS_VECTOR_TEST(I_QP,I_DOF) * BASIS_VECTOR_TEST(I_QP,J_DOF) * WEQUA3D(I_QP)
!             END DO
!         END DO
!     END DO
!     
!     MM_TEST = MM_TEST * CELL_VOLUME
!     CALL COMPMASSINV(RESHAPE(MM_TEST,(/1,IDEGFREE,IDEGFREE/)),INV_MM,IDEGFREE,1)
!     
!     IF (ANY(MATMUL(MM_TEST, INV_MM(1,:,:)) - IDENTITY_MATRIX > 1E-12)) THEN
!         WRITE(ERROR_FILE_INDEX,*) 'MM_TEST', MM_TEST
!         WRITE(ERROR_FILE_INDEX,*) 'INV MM', INV_MM
!         WRITE(ERROR_FILE_INDEX,*) 'MM * INV_MM', MATMUL(MM_TEST, INV_MM(1,:,:))
!         ERROR_YES = 1
!     END IF
!     
!     MM_TEST1 = 0.
!     DO I_DOF = 1, N_DOFS
!         DO I_QP = 1, N_QP
!             MM_TEST1(I_DOF) = MM_TEST1(I_DOF) + BASIS_VECTOR_TEST(I_QP,I_DOF) * WEQUA3D(I_QP)
!         END DO
!     END DO
!     
!     IF (ANY(MM_TEST > 1E12)) THEN
!         WRITE(ERROR_FILE_INDEX,*) 'MM_TEST1', MM_TEST1 !Should be zero for every element
!         ERROR_YES = 1
!     END IF
!     
!     IF(ERROR_YES == 0) WRITE(ERROR_FILE_INDEX,*) 'TEST_MASS_MATRIX PASSED'
! 
! END SUBROUTINE TEST_MASS_MATRIX
! 
! SUBROUTINE RUN_ALL_TESTS(N)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N
!     INTEGER::NUM_VARS,ORDER,NUM_DOFS,N_DIM,ERROR_FILE_INDEX
!     REAL,DIMENSION(5,3)::TEST_MATRIX
!     REAL,DIMENSION(5)::TEST_MATRIX_B
!     REAL,DIMENSION(3)::TEST_MATRIX_LSQ_ANS
!     
!     NUM_VARS = 1
!     N_DIM = 2
!     ORDER = IORDER
!     ERROR_FILE_INDEX = 400 + N
!     
!     IF (N_DIM == 2) THEN
!         NUM_DOFS = ((ORDER+1) * (ORDER+2) / 2) - 1
!     ELSE
!         NUM_DOFS = ((ORDER+1)*(ORDER+2)*(ORDER+3)/6) - 1
!     END IF
!     
!     CALL CHECK_GLOBALS(N, ERROR_FILE_INDEX)
!     
!     CALL TEST_BASIS_REC2D(N,ORDER,NUM_DOFS, ERROR_FILE_INDEX)
!     IF (TEST_QP_ARRAY(N) == 0) WRITE(400+N,*) 'TEST_QP_ARRAY passed'
!     IF (TEST_QUADRATURE(N,ORDER) == 0) WRITE(400+N,*) 'TEST_QUADRATURE passed'
!     
!     IF (POLY == 1) THEN
!         CALL TEST_DG_SOL(N, NUM_VARS, ORDER, NUM_DOFS, ERROR_FILE_INDEX)
!         IF (TEST_DG_RHS_INTEGRAL(N,NUM_VARS,ORDER,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_DG_RHS_INTEGRAL passed'
!     END IF
!     
!     IF (TEST_CALC_DELTA_XYZ(N_DIM) == 0) WRITE(400+N,*) 'TEST_CALC_DELTA_XYZ passed'
!     IF (TEST_TAYLOR_BASIS(N,ORDER,N_DIM,NUM_DOFS) == 0) WRITE(400+N,*) 'TEST_TAYLOR_BASIS passed'
!     
!     TEST_MATRIX = RESHAPE((/ 1., 1., 1., 1., 1., 76., 92., 68., 86., 54., 48., 92., 82., 68., 70. /), SHAPE(TEST_MATRIX))
!     TEST_MATRIX_B = RESHAPE((/ 43., 90., 64., 69., 50. /), SHAPE(TEST_MATRIX_B))
!     TEST_MATRIX_LSQ_ANS = (/ -42.4, 0.639, 0.799 /)
!     
!     CALL TEST_NORMAL_LSQ(TEST_MATRIX, TEST_MATRIX_B, TEST_MATRIX_LSQ_ANS, ERROR_FILE_INDEX)
!     
! !     IF (TEST_QR_DECOMPOSITION(TEST_MATRIX,ERROR_FILE_INDEX) == 0) WRITE(ERROR_FILE_INDEX,*) 'TEST_QR_DECOMPOSITION passed'
!     
!     CALL TEST_MASS_MATRIX(ERROR_FILE_INDEX)
!     CALL PRINT_SOME_INFO(N)
!     
! END SUBROUTINE
! 
! END MODULE
