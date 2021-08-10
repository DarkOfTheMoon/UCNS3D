MODULE DG_FUNCTIONS

USE BASIS

IMPLICIT NONE

CONTAINS

FUNCTION DG_SOL(N, X_IN, Y_IN, NUM_VARIABLES, IORDER, IDEGFREE, U_C_VALDG)
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, IDEGFREE: number of basis terms

    INTEGER,INTENT(IN)::N,IORDER,IDEGFREE,NUM_VARIABLES
    REAL,INTENT(IN)::X_IN,Y_IN ! Coordinates of the point where the solution is requested
    REAL,DIMENSION(:,:),INTENT(IN)::U_C_VALDG
    REAL,DIMENSION(IDEGFREE)::BASIS_TEMP
    INTEGER::I_DOF, I_VAR
    REAL,DIMENSION(NUM_VARIABLES)::DG_SOL

    IF(ALL(SHAPE(U_C_VALDG) /= (/ NUM_VARIABLES,IDEGFREE+1 /))) THEN
        WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
        STOP
    END IF
    
    BASIS_TEMP = BASIS_REC2D(N,X_IN,Y_IN,IORDER,0,IDEGFREE)

    DO I_VAR = 1, NUM_VARIABLES
        DG_SOL = U_C_VALDG(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(:), U_C_VALDG(I_VAR,2:))
    END DO

END FUNCTION DG_SOL

FUNCTION DG_RHS_INTEGRAL(N,QP_X,QP_Y,QP_WEIGHT,NUM_VARS,ORDER,NUM_DOFS,CELL_VOL_OR_SURF,FLUX_TERM,VOL_OR_SURF)
!> @brief
!> Calculates the volume or surface integral term in the DG RHS for scalar linear advection with speed = 1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N,ORDER,NUM_VARS,NUM_DOFS,VOL_OR_SURF
    REAL,INTENT(IN)::QP_X,QP_Y,QP_WEIGHT,CELL_VOL_OR_SURF
    REAL,DIMENSION(:),INTENT(IN)::FLUX_TERM
    INTEGER::I_VAR
    REAL,DIMENSION(NUM_VARS,NUM_DOFS+1)::DG_RHS_INTEGRAL
    
    IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
        WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
        STOP
    END IF
    
    IF (VOL_OR_SURF == 1) THEN ! VOLUME INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(I_VAR,1) = 0
            DG_RHS_INTEGRAL(I_VAR,2:) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * (BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,0,NUM_DOFS,1) + BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,0,NUM_DOFS,2))
        END DO
    ELSE IF (VOL_OR_SURF == 2) THEN ! SURFACE INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(I_VAR,1) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF
            DG_RHS_INTEGRAL(I_VAR,2:) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * BASIS_REC2D(N,QP_X,QP_Y,ORDER,0,NUM_DOFS)
        END DO
    END IF

END FUNCTION DG_RHS_INTEGRAL

SUBROUTINE LUO_LSQ_RECONSTRUCT(N)
INTEGER,INTENT(IN)::N
    
END SUBROUTINE LUO_LSQ_RECONSTRUCT

END MODULE DG_FUNCTIONS