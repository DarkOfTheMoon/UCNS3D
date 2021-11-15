MODULE DG_FUNCTIONS

USE BASIS
USE DECLARATION
USE DERIVATIVES

IMPLICIT NONE

CONTAINS

FUNCTION DG_SOL(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
INTEGER,INTENT(IN)::N
INTEGER::I_DOF, I_VAR
REAL,DIMENSION(Nof_VARIABLES)::DG_SOL

    NUMBER = IELEM(N,ICONSIDERED)%IORDER
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)

DO I_VAR = 1, NOF_VARIABLES
    DG_SOL(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:NUMBER_OF_DOG), U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1))
END DO

END FUNCTION DG_SOL

FUNCTION DG_SOL_I_CENTER(N, I)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    INTEGER,INTENT(IN)::N, I
    INTEGER::I_VAR
    REAL,DIMENSION(NOF_VARIABLES)::DG_SOL_I_CENTER

    NUMBER = IELEM(N, I)%IORDER
    X1 = IELEM(N, I)%XXC
    Y1 = IELEM(N, I)%YYC
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, I, NUMBER_OF_DOG)

    DO I_VAR = 1, NOF_VARIABLES
        DG_SOL_I_CENTER(I_VAR) = U_C(I)%VALDG(1, I_VAR, 1) + DOT_PRODUCT(BASIS_TEMP, U_C(I)%VALDG(1, I_VAR, 2:))
    END DO

END FUNCTION DG_SOL_I_CENTER

FUNCTION DG_SOLFACE(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given surface point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
INTEGER::I_DOF, I_VAR
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:nof_Variables)::DG_SOLFACE

!     IF(ALL(SHAPE(U_C_VALDG) /= (/ nof_Variables, NUMBER_OF_DOG+1 /))) THEN
!         WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
!         STOP
!     END IF

! X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
! Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
X1=ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1) - 1.0d0 / 3.0d0
Y1=ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2) - 1.0d0 / 3.0d0
NUMBER=IELEM(N,ICONSIDERED)%IORDER

BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)

DO I_VAR = 1, NOF_VARIABLES
    DG_SOLFACE(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:NUMBER_OF_DOG), U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1))
END DO

END FUNCTION DG_SOLFACE

FUNCTION DG_SURF_FLUX(N)
!> @brief
!> Calculates the RHS flux term to be integrated in the DG formulation
IMPLICIT NONE
REAL,DIMENSION(NUMBER_OF_DOG+1,NOF_VARIABLES)::DG_SURF_FLUX
INTEGER::I_VAR
INTEGER,INTENT(IN)::N

! X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
! Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
X1=ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1) - 1.0d0 / 3.0d0
Y1=ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2) - 1.0d0 / 3.0d0
NUMBER=IELEM(N,ICONSIDERED)%IORDER

DO I_VAR = 1, NOF_VARIABLES
    DG_SURF_FLUX(1,I_VAR) = HLLCFLUX(I_VAR) * WEQUA2D(POINTX) * IELEM(N,ICONSIDERED)%SURF(FACEX)
    DG_SURF_FLUX(2:NUMBER_OF_DOG+1,I_VAR) = HLLCFLUX(I_VAR) * WEQUA2D(POINTX) * IELEM(N,ICONSIDERED)%SURF(FACEX) * BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
END DO

END FUNCTION DG_SURF_FLUX

FUNCTION DG_CELL_AVG(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(NOF_VARIABLES)::DG_CELL_AVG
REAL,DIMENSION(NOF_VARIABLES)::FLUX_TERM
INTEGER::I,NQP,I_QP,I_VAR

IF (IELEM(N,ICONSIDERED)%ISHAPE == 5) NQP = QP_TRIANGLE * 2
IF (IELEM(N,ICONSIDERED)%ISHAPE == 6) NQP = QP_TRIANGLE
NUMBER=IELEM(N,ICONSIDERED)%IORDER
NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE
           
DO I_VAR = 1, NOF_VARIABLES
    DG_CELL_AVG(I_VAR) = 0.0d0
      
    DO I_QP = 1, NQP 
        X1=QP_ARRAY(ICONSIDERED,I_QP)%X
        Y1=QP_ARRAY(ICONSIDERED,I_QP)%Y
            
        FLUX_TERM = DG_SOL(N)
        
        IF (ITESTCASE.EQ.3)THEN
            LEFTV=DG_SOL(N)
            CALL PRIM2CONS2D(N)
            FLUX_TERM=FLUXEVAL2D(LEFTV)
        END IF
            
        DG_CELL_AVG(I_VAR) = DG_CELL_AVG(I_VAR) + FLUX_TERM(I_VAR) * QP_ARRAY(ICONSIDERED,I_QP)%QP_WEIGHT!* IELEM(N,ICONSIDERED)%TOTVOLUME
    END DO
    DG_CELL_AVG(I_VAR) = DG_CELL_AVG(I_VAR) / IELEM(N,ICONSIDERED)%TOTVOLUME
END DO

END FUNCTION DG_CELL_AVG

FUNCTION DG_VOL_INTEGRAL(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(IDEGFREE+1,NOF_VARIABLES)::DG_VOL_INTEGRAL
REAL,DIMENSION(NOF_VARIABLES)::FLUX_TERM
REAL,DIMENSION(DIMENSIONA)::DER_JAC
INTEGER::I_DOF,NQP,I_QP,I_VAR

IF (IELEM(N,ICONSIDERED)%ISHAPE == 5) NQP = QP_TRIANGLE*2
IF (IELEM(N,ICONSIDERED)%ISHAPE == 6) NQP = QP_TRIANGLE
NUMBER=IELEM(N,ICONSIDERED)%IORDER
NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE
           
DO I_VAR = 1, NOF_VARIABLES
    DG_VOL_INTEGRAL(1,I_VAR) = 0.0d0
    DG_VOL_INTEGRAL(2:,I_VAR)= 0.0D0

      DO I_DOF=1,NUMBER_OF_DOG
      
        DO I_QP = 1, NQP 
            X1=ILOCAL_RECON3(ICONSIDERED)%VOL_QPOINTS(I_QP, 1)!QP_ARRAY(ICONSIDERED,I_QP)%X
            Y1=ILOCAL_RECON3(ICONSIDERED)%VOL_QPOINTS(I_QP, 2)!QP_ARRAY(ICONSIDERED,I_QP)%Y
            
            FLUX_TERM=DG_SOL(N)
            
            IF (ITESTCASE.EQ.3)THEN
                LEFTV=DG_SOL(N)
                CALL PRIM2CONS2D(N)
                FLUX_TERM=FLUXEVAL2D(LEFTV)
            END IF
            
            DER_JAC = (/ DF2DX(X1,Y1,I_DOF), DF2DY(X1,Y1,I_DOF) /)
            DER_JAC = MATMUL(DER_JAC, ILOCAL_RECON3(ICONSIDERED)%INVCCJAC)
            
            DG_VOL_INTEGRAL(I_DOF+1,I_VAR) = DG_VOL_INTEGRAL(I_DOF+1,I_VAR) + FLUX_TERM(I_VAR) * ILOCAL_RECON3(ICONSIDERED)%VOL_QPWEIGHTS(I_QP) * (DER_JAC(1) + DER_JAC(2)) * IELEM(N,ICONSIDERED)%TOTVOLUME
        END DO
     END DO
END DO

END FUNCTION DG_VOL_INTEGRAL


SUBROUTINE RECONSTRUCT_DG(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I_FACE, I_ELEM, I_QP
        
!$OMP DO
DO I_ELEM = 1, XMPIELRANK(N)
    DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA
            !SOMEWHERE PRESTORED THE GUASSIAN QUADRATURE POINTS FOR YOUR SIDES IN ANOTHER SUBROUTINE AND YOU BUILD ONLY THE cleft STATES AND VOLUME INTEGRAL
            
        DO I_QP = 1, QP_LINE_N
            
            FACEX=I_FACE
            POINTX=I_QP
            ICONSIDERED=I_ELEM
            NUMBER_OF_DOG=IELEM(N,I_ELEM)%IDEGFREE
            ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(:, FACEX, POINTX) = DG_SOLFACE(N)
        
!             write(700+n,*)"look here","cell",i_elem,"face",facex,"point",pointx
!             write(700+n,*)ILOCAL_RECON3(I_ELEM)%ULEFT_DG(1:NOF_VARIABLES, I_FACE, I_QP)
                
                !STORE IT HERE (ILOCAL_RECON3(I)%ULEFT_DG(1,1:FACES,1:NGP) ! you need to allocate it in memory
        END DO
    END DO
END DO
!$OMP END DO

END SUBROUTINE RECONSTRUCT_DG

FUNCTION CALC_DELTA_XYZ(NUM_NODES, NUM_DIMS, NODES_IN)
!> @brief
!> Calculates the "delta x/y/z" as in Luo 2012 eq 3.12
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NUM_NODES, NUM_DIMS
    REAL,DIMENSION(:,:),INTENT(IN)::NODES_IN ! (NODE, DIMENSION)
    INTEGER::I_NODES, I_DIM
    REAL::XYZ_MAX,XYZ_MIN
    REAL,DIMENSION(NUM_DIMS)::CALC_DELTA_XYZ
    
    DO I_DIM = 1, NUM_DIMS
        XYZ_MAX = -1E12
        XYZ_MIN = 1E12
        
        DO I_NODES = 1, NUM_NODES
            IF (NODES_IN(I_NODES,I_DIM) > XYZ_MAX) XYZ_MAX = NODES_IN(I_NODES,I_DIM)
            IF (NODES_IN(I_NODES,I_DIM) < XYZ_MIN) XYZ_MIN = NODES_IN(I_NODES,I_DIM)
        END DO
        
        CALC_DELTA_XYZ(I_DIM) = 0.5 * ABS(XYZ_MAX - XYZ_MIN)
    END DO
    
END FUNCTION CALC_DELTA_XYZ

SUBROUTINE ALLOCATE_DG_MASS_MATRICES(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::KMAXE
    
    KMAXE = XMPIELRANK(N)
    
    ALLOCATE(MASS_MATRIX_CENTERS(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); MASS_MATRIX_CENTERS(N,:,:,:) = ZERO
    ALLOCATE(INV_MASS_MATRIX(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); INV_MASS_MATRIX(N,:,:,:) = ZERO
END SUBROUTINE ALLOCATE_DG_MASS_MATRICES

SUBROUTINE PRESTORE_AND_ALLOCATE_DG(N, I)
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY, mass matrix
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, I
    INTEGER::K, I_QP, COUNT_1, DECOMP_BOOL
        
!     DO I = 1, XMPIELRANK(N)
        
        DECOMP_BOOL = 0

        !Store delta xyz (normalization factor from Luo 2012)
        IELEM(N,I)%DELTA_XYZ = CALC_DELTA_XYZ(IELEM(N,I)%NONODES, DIMENSIONA, NODES_LIST)
        
        IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
            ALLOCATE(ILOCAL_RECON3(I)%VOL_QPOINTS(QP_TRIANGLE*2, DIMENSIONA))
            ALLOCATE(ILOCAL_RECON3(I)%VOL_QPWEIGHTS(QP_TRIANGLE*2))
        ELSE
            ALLOCATE(ILOCAL_RECON3(I)%VOL_QPOINTS(QP_TRIANGLE, DIMENSIONA))
            ALLOCATE(ILOCAL_RECON3(I)%VOL_QPWEIGHTS(QP_TRIANGLE))
        END IF
        
        DO K = 1,IELEM(N,I)%NONODES
            NODES_LIST(k,1:DIMENSIONA)=INODER(IELEM(N,I)%NODES(K))%CORD(1:DIMENSIONA)
            VEXT(k,1:DIMENSIONA)=NODES_LIST(k,1:DIMENSIONA)
            VEXT(k,1:DIMENSIONA)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:DIMENSIONA)-ILOCAL_RECON3(I)%VEXT_REF(1:DIMENSIONA)) ! Transforming to reference space
        END DO
        
        ALLOCATE(IELEM(N,I)%TAYLOR_INTEGRAL(IDEGFREE - DIMENSIONA))
        IELEM(N,I)%TAYLOR_INTEGRAL = ZERO
    
        !Store volume quadrature points and Taylor integral
        SELECT CASE(ielem(n,i)%ishape)
        CASE(5) ! Quad
            IF (DECOMP_BOOL == 1) THEN
                ELTYPE=IELEM(N,I)%ISHAPE ! Needed for DECOMPOSE2
                CALL DECOMPOSE2
        
                COUNT_1=0
                DO K=1,ELEM_DEC
                    VEXT(1:3,1:DIMENSIONA)=ELEM_LISTD(k,1:3,1:DIMENSIONA)
                
                    CALL QUADRATUREtriangle(N,IGQRULES)
                    
                    VOLTEMP = TRIANGLEVOLUME(N) / IELEM(N,I)%totvolume
                    
                    DO I_QP = 1, QP_Triangle
!                         COUNT_1=COUNT_1+1
!                         QP_ARRAY(I,COUNT_1)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
!                         QP_ARRAY(I,COUNT_1)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
!                         QP_ARRAY(I,COUNT_1)%QP_WEIGHT = WEQUA3D(I_QP) * VOLTEMP
                        ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1:DIMENSIONA) = QPOINTS(1:DIMENSIONA, I_QP) - 0.50D0 ! Center of reference element
                        ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP) = WEQUA3D(I_QP) * VOLTEMP
                        
                        IELEM(N,I)%TAYLOR_INTEGRAL(1) = IELEM(N,I)%TAYLOR_INTEGRAL(1) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP + QP_TRIANGLE * (K - 1), 1) / IELEM(N,I)%DELTA_XYZ(1)) ** 2 / 2 * ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP)
                        IELEM(N,I)%TAYLOR_INTEGRAL(2) = IELEM(N,I)%TAYLOR_INTEGRAL(2) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP + QP_TRIANGLE * (K - 1), 2) / IELEM(N,I)%DELTA_XYZ(2)) ** 2 / 2 * ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP)
                        IELEM(N,I)%TAYLOR_INTEGRAL(3) = IELEM(N,I)%TAYLOR_INTEGRAL(3) + ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP + QP_TRIANGLE * (K - 1), 1) / IELEM(N,I)%DELTA_XYZ(1) *  ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP + QP_TRIANGLE * (K - 1), 2) / IELEM(N,I)%DELTA_XYZ(2) * ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP)
                    END DO
                END DO
            ELSE ! DECOMP_BOOL /= 1
                CALL QUADRATUREQUAD(N, IGQRULES)
                
                 DO I_QP = 1, QP_QUAD
!                     QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - 0.50D0 !CENTER OF REF ELEM IELEM(N,I)%XXC
!                     QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - 0.50D0 !CENTER OF REF ELEM IELEM(N,I)%YYC
!                     QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP)! * IELEM(N,I)%totvolume
                    ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1:DIMENSIONA) = QPOINTS(1:DIMENSIONA, I_QP) - 0.50D0 ! Center of reference element
                    ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP) = WEQUA3D(I_QP)
                    
                    IELEM(N,I)%TAYLOR_INTEGRAL(1) = IELEM(N,I)%TAYLOR_INTEGRAL(1) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1) / IELEM(N,I)%DELTA_XYZ(1)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                    IELEM(N,I)%TAYLOR_INTEGRAL(2) = IELEM(N,I)%TAYLOR_INTEGRAL(2) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 2) / IELEM(N,I)%DELTA_XYZ(2)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                    IELEM(N,I)%TAYLOR_INTEGRAL(3) = IELEM(N,I)%TAYLOR_INTEGRAL(3) + ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1) / IELEM(N,I)%DELTA_XYZ(1) *  ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 2) / IELEM(N,I)%DELTA_XYZ(2) * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                END DO
            END IF
        CASE(6) !Tri
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            
            DO I_QP = 1, QP_Triangle
!                 QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - 1.0D0 / 3.0D0 !CENTER OF REF ELEM IELEM(N,I)%XXC
!                 QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - 1.0D0 / 3.0D0 !CENTER OF REF ELEM IELEM(N,I)%YYC
!                 QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP)! * IELEM(N,I)%totvolume
                ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1:DIMENSIONA) = QPOINTS(1:DIMENSIONA, I_QP) - 1.0D0 / 3.0D0 ! Center of reference element
                ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP) = WEQUA3D(I_QP)
                
                IELEM(N,I)%TAYLOR_INTEGRAL(1) = IELEM(N,I)%TAYLOR_INTEGRAL(1) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1) / IELEM(N,I)%DELTA_XYZ(1)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME !No multiplication because cancelled out
                IELEM(N,I)%TAYLOR_INTEGRAL(2) = IELEM(N,I)%TAYLOR_INTEGRAL(2) + (ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 2) / IELEM(N,I)%DELTA_XYZ(2)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME !No multiplication because cancelled out
                IELEM(N,I)%TAYLOR_INTEGRAL(3) = IELEM(N,I)%TAYLOR_INTEGRAL(3) + ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1) / IELEM(N,I)%DELTA_XYZ(1) *  ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1) / IELEM(N,I)%DELTA_XYZ(2) * WEQUA3D(I_QP) ! * CELL_VOLUME !No multiplication because cancelled out
            END DO
!             IELEM(N,I)%TAYLOR_INTEGRAL(:) = IELEM(N,I)%TAYLOR_INTEGRAL(:) * 2.0d0 ! Divided by volume
        END SELECT
        
        ALLOCATE(ILOCAL_RECON3(I)%ULEFT_DG(NOF_VARIABLES, IELEM(N,I)%IFCA, QP_LINE_N))
!     END DO !1, XMPIELRANK(N)
    
    CALL ASS_MASS_MATRIX(N, I)
    
END SUBROUTINE PRESTORE_AND_ALLOCATE_DG
        
SUBROUTINE PRESTORE_AND_ALLOCATE_DG_SURF_QPOINTS  
!> @brief
!> Prestores SURF_QPOINTS
    IMPLICIT NONE
    INTEGER::I, K, I_QP, I_FACE, NND, IDUMMY
    
    DO I = 1, XMPIELRANK(N)
        ICONSIDERED = I
        ALLOCATE(ILOCAL_RECON3(I)%SURF_QPOINTS(IELEM(N,I)%IFCA, QP_LINE_N, DIMENSIONA))
        DO I_FACE = 1, IELEM(N,I)%IFCA
            IDUMMY=0
            !GAUSSIAN POINTS FIXED
            IF ((IPERIODICITY.EQ.1).AND.(IELEM(N,I)%INTERIOR.EQ.1))THEN	
                IF (IELEM(N,I)%IBOUNDS(I_FACE).GT.0)THEN	!CHECK FOR BOUNDARIES
                    IF (IBOUND(N,IELEM(N,I)%IBOUNDS(I_FACE))%ICODE.EQ.5)THEN	!PERIODIC
                        IDUMMY=1
                    END IF
                END IF
        
                NND=2
                IF (IDUMMY.EQ.0)THEN
                    DO K=1,NND
                        VEXT(K,1:2)=INODER(IELEM(N,I)%NODES_FACES(I_FACE,K))%CORD(1:DIMS)
                    END DO
                ELSE
                    FACEX=I_FACE;
                    CALL COORDINATES_FACE_PERIOD2D1(N,ICONSIDERED,FACEX)
                END IF
                CALL QUADRATURELINE(N,IGQRULES)	  
            ELSE
                NND=2
                DO K=1,NND
                    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(I_FACE,K))%CORD(1:dims)
                END DO
                CALL QUADRATURELINE(N,IGQRULES)
            END IF

            DO I_QP = 1, QP_LINE_N
                ILOCAL_RECON3(I)%SURF_QPOINTS(I_FACE,I_QP,1) = QPOINTS2D(1,I_QP) - IELEM(N,I)%XXC
                ILOCAL_RECON3(I)%SURF_QPOINTS(I_FACE,I_QP,2) = QPOINTS2D(2,I_QP) - IELEM(N,I)%YYC
            END DO
        END DO
    END DO
END SUBROUTINE PRESTORE_AND_ALLOCATE_DG_SURF_QPOINTS


SUBROUTINE ASS_MASS_MATRIX(N, I)
!> @brief
!> Assembles the mass matrix
!> REQUIRES: Globals: IELEM, QP_QUAD, QP_TRIANGLE, MASS_MATRIX
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, I
    INTEGER::I_QP, N_QP, I_DOF, J_DOF, KMAXE
    REAL,DIMENSION(IDEGFREE)::BASIS_VECTOR
    REAL::INTEG_TEST,PHX,INTEG_MM
    
!     KMAXE = XMPIELRANK(N)
    
!     DO I_ELEM = 1, KMAXE
        SELECT CASE(IELEM(N,I)%ISHAPE)
        CASE(5) ! Quadrilateral
            N_QP = QP_QUAD!QP_TRIANGLE*2
        CASE(6) ! Triangle
            N_QP = QP_TRIANGLE
        END SELECT
        
        NUMBER_OF_DOG = IELEM(N,I)%IDEGFREE
            
        DO I_DOF = 1, NUM_DG_DOFS       
            DO J_DOF = 1, NUM_DG_DOFS
                INTEG_MM = ZERO
                
                    DO I_QP = 1, N_QP
                        IXX = I; x1=ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 1); y1=ILOCAL_RECON3(I)%VOL_QPOINTS(I_QP, 2)
                        BASIS_VECTOR = BASIS_REC2D(N,X1,Y1,IORDER,IXX,NUMBER_OF_DOG)
                        
                            IF (I_DOF==1.and.J_DOF==1)THEN
                                PHX = 1
                                
                            ELSE IF (I_DOF==1.and.J_DOF/=1)THEN
                                PHX = BASIS_VECTOR(J_DOF-1)
                                
                            ELSE IF (I_DOF/=1.and.J_DOF==1)THEN
                                PHX = BASIS_VECTOR(I_DOF-1) 
                                
                            ELSE IF (I_DOF/=1.and.J_DOF/=1)THEN
                                PHX = BASIS_VECTOR(I_DOF-1)*BASIS_VECTOR(J_DOF-1)
                                
                            END IF
                        INTEG_MM = INTEG_MM + PHX * ILOCAL_RECON3(I)%VOL_QPWEIGHTS(I_QP)
                    END DO
            
                MASS_MATRIX_CENTERS(N, I, I_DOF, J_DOF) = MASS_MATRIX_CENTERS(N, I, I_DOF, J_DOF) + INTEG_MM
            
            END DO
        END DO
!     END DO !I_ELEM = 1, KMAXE

    MASS_MATRIX_CENTERS(N, I, :, :) = MASS_MATRIX_CENTERS(N, I, :, :) * IELEM(N,I)%TOTVOLUME
    INV_MASS_MATRIX(N,I,:,:) = INVERSE_MATRIX(MASS_MATRIX_CENTERS(N,I,:,:), NUM_DG_DOFS)
    
!     DO I_ELEM = 1, KMAXE
!         WRITE(300+N,*) I_ELEM
!         DO I_QP = 1, N_QP
!             WRITE(300+N,*) 'QP_ARRAY', QP_ARRAY(I_ELEM,I_QP)%X, QP_ARRAY(I_ELEM,I_QP)%Y, QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!             WRITE(300+N,*) 'BASIS,', BASIS_REC2D(N,QP_ARRAY(I_ELEM,I_QP)%X,QP_ARRAY(I_ELEM,I_QP)%Y,IORDER,I_ELEM,IDEGFREE)
!         END DO
!         WRITE(300+N,*) 'XYZ', IELEM(N,I_ELEM)%NODES
!         WRITE(300+N,*) 'DELTAXYZ', IELEM(N,I_ELEM)%DELTA_XYZ
!         WRITE(300+N,*) 'MMC', MASS_MATRIX_CENTERS(N,I_ELEM,:,:)
!         WRITE(300+N,*) 'Inverse,', INV_MASS_MATRIX(N,I_ELEM,:,:)
!         WRITE(300+N,*) 'Identity', MATMUL(MASS_MATRIX_CENTERS(N,I_ELEM,:,:),INV_MASS_MATRIX(N,I_ELEM,:,:))
!     END DO
    
END SUBROUTINE ASS_MASS_MATRIX

SUBROUTINE COMPMASSINV(totalMM,invMM,N_DOFS)
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
 
integer :: i,j,k,l,m,irow,P,kmaxe
real:: big,dum
real,DIMENSION(N_DOFS,N_DOFS)::a,b
integer,INTENT(IN)::N_DOFS
REAL,DIMENSION(:,:,:),INTENT(IN)::totalMM
REAL,DIMENSION(:,:,:),INTENT(INOUT)::invMM
kmaxe=xmpielrank(n)
DO P=1,kmaxe

a(:,:)=totalMM(P,:,:)
b(:,:)=zero

do i = 1,N_DOFS
    do j = 1,N_DOFS
        b(i,j) = 0.0
    end do
    b(i,i) = 1.0
end do 

do i = 1,N_DOFS   
   big = a(i,i)
   do j = i,N_DOFS
     if (a(j,i).gt.big) then
       big = a(j,i)
       irow = j
     end if
   end do
   ! interchange lines i with irow for both a() and b() matrices
   if (big.gt.a(i,i)) then
     do k = 1,N_DOFS
       dum = a(i,k)                      ! matrix a()
       a(i,k) = a(irow,k)
       a(irow,k) = dum
       dum = b(i,k)                 ! matrix b()
       b(i,k) = b(irow,k)
       b(irow,k) = dum
     end do
   end if
   ! divide all entries in line i from a(i,j) by the value a(i,i); 
   ! same operation for the identity matrix
   dum = a(i,i)
   do j = 1,N_DOFS
     a(i,j) = a(i,j)/dum
     b(i,j) = b(i,j)/dum
   end do
   ! make zero all entries in the column a(j,i); same operation for indent()
   do j = i+1,N_DOFS
     dum = a(j,i)
     do k = 1,N_DOFS
       a(j,k) = a(j,k) - dum*a(i,k)
       b(j,k) = b(j,k) - dum*b(i,k)               
            
     end do
   end do
end do
  
 do i = 1,N_DOFS-1
   do j = i+1,N_DOFS
     dum = a(i,j)
     do l = 1,N_DOFS
       a(i,l) = a(i,l)-dum*a(j,l)
       b(i,l) = b(i,l)-dum*b(j,l)
     end do
   end do
 end do
 
 invMM(P,:,:)=b(:,:)
  
END DO
 
END SUBROUTINE COMPMASSINV

FUNCTION INVERSE_MATRIX(MATRIX_IN, N_DOFS)
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
    INTEGER,INTENT(IN)::N_DOFS
    REAL,DIMENSION(N_DOFS,N_DOFS),INTENT(IN)::MATRIX_IN
    integer :: i,j,k,l,m,irow
    real:: big,dum
    real,DIMENSION(N_DOFS,N_DOFS)::a,b
    REAL,DIMENSION(N_DOFS,N_DOFS)::INVERSE_MATRIX

    a(:,:)=MATRIX_IN
    b(:,:)=zero

    do i = 1,N_DOFS
        do j = 1,N_DOFS
            b(i,j) = 0.0
        end do
        b(i,i) = 1.0
    end do 

    do i = 1,N_DOFS   
    big = a(i,i)
    do j = i,N_DOFS
        if (a(j,i).gt.big) then
        big = a(j,i)
        irow = j
        end if
    end do
    ! interchange lines i with irow for both a() and b() matrices
    if (big.gt.a(i,i)) then
        do k = 1,N_DOFS
        dum = a(i,k)                      ! matrix a()
        a(i,k) = a(irow,k)
        a(irow,k) = dum
        dum = b(i,k)                 ! matrix b()
        b(i,k) = b(irow,k)
        b(irow,k) = dum
        end do
    end if
    ! divide all entries in line i from a(i,j) by the value a(i,i); 
    ! same operation for the identity matrix
    dum = a(i,i)
    do j = 1,N_DOFS
        a(i,j) = a(i,j)/dum
        b(i,j) = b(i,j)/dum
    end do
    ! make zero all entries in the column a(j,i); same operation for indent()
    do j = i+1,N_DOFS
        dum = a(j,i)
        do k = 1,N_DOFS
        a(j,k) = a(j,k) - dum*a(i,k)
        b(j,k) = b(j,k) - dum*b(i,k)               
                
        end do
    end do
    end do
    
    do i = 1,N_DOFS-1
    do j = i+1,N_DOFS
        dum = a(i,j)
        do l = 1,N_DOFS
        a(i,l) = a(i,l)-dum*a(j,l)
        b(i,l) = b(i,l)-dum*b(j,l)
        end do
    end do
    end do
    
    INVERSE_MATRIX = b(:,:)
    
END FUNCTION INVERSE_MATRIX

FUNCTION NORMAL_LSQ(A, b)
    REAL,DIMENSION(:,:),INTENT(IN)::A
    REAL,DIMENSION(:),INTENT(IN)::b
    REAL,DIMENSION(SIZE(A,1))::NORMAL_LSQ
    
    NORMAL_LSQ = MATMUL(MATMUL(INVERSE_MATRIX(MATMUL(TRANSPOSE(A), A), SIZE(A,2)), TRANSPOSE(A)), b)
    
END FUNCTION NORMAL_LSQ

SUBROUTINE LUO_LSQ_RECONSTRUCT(N)
!> @brief
!> This function reconstructs an approximation of NUM_DG_RECONSTRUCT_DOFS
!> REQUIRES: IELEM, U_C as globals
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_FACE, I_QP, I_VAR
    REAL,DIMENSION(NUM_DG_RECONSTRUCT_DOFS-1)::BASIS_TEMP
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_RECONSTRUCT_DOFS)::RECONSTRUCTED_SOL
    
    !$OMP DO
    DO I_ELEM = 1, XMPIELRANK(N)
        ICONSIDERED = I_ELEM
        NUMBER = IELEM(N,ICONSIDERED)%IORDER
        
        RECONSTRUCTED_SOL(:, 1:NUM_DG_DOFS) = U_C(ICONSIDERED)%VALDG(1, :, :)
        RECONSTRUCTED_SOL(:, NUM_DG_DOFS+1:NUM_DG_RECONSTRUCT_DOFS) = LUO_LSQ_RECONSTRUCT_SOL(N)
        
        DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA ! Set ULEFT_DG
            FACEX = I_FACE
            DO I_QP = 1, QP_LINE_N
                POINTX = I_QP
                
                X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
                Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
                BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER + 1, ICONSIDERED, NUM_DG_RECONSTRUCT_DOFS-1)
                
                DO I_VAR = 1, NOF_VARIABLES
                    ILOCAL_RECON3(I_ELEM)%ULEFT_DG(I_VAR, I_FACE, I_QP) = RECONSTRUCTED_SOL(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:NUM_DG_RECONSTRUCT_DOFS-1), RECONSTRUCTED_SOL(I_VAR, 2:NUM_DG_RECONSTRUCT_DOFS))
                END DO
                
!                 write(700+n,*)"look here","cell",i_elem,"face",facex,"point",pointx
!                 write(700+n,*)ILOCAL_RECON3(I_ELEM)%ULEFT_DG(:, I_FACE, I_QP)
            END DO
        END DO
    END DO
    !$OMP END DO
    
END SUBROUTINE LUO_LSQ_RECONSTRUCT

FUNCTION LUO_LSQ_RECONSTRUCT_SOL(N)
!> @brief
!> Returns the higher dofs reconstructed solution using the Luo LSQ method
    IMPLICIT NONE
    INTEGER, INTENT(IN)::N
    INTEGER::NEIGHBOR_INDEX, I_DIM, I_FACE, I_DG_RECONSTRUCT_DOF, I_VAR, NEIGHBOR_TYPE
    REAL,DIMENSION(NUM_DG_RECONSTRUCT_DOFS - 1)::BASIS_TEMP
    REAL,DIMENSION(NOF_VARIABLES)::DG_SOL_TEMP, DG_SOL_I_J
    REAL,DIMENSION(DIMENSIONA)::NEIGHBOR_DELTA_XYZ
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_DOFS)::NEIGHBOR_U
    REAL,DIMENSION(NOF_VARIABLES, IELEM(N, ICONSIDERED)%IFCA * NUM_DG_DOFS, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS)::LHS_MATRIX ! See eq. 3.19 of Luo 2012
    REAL,DIMENSION(NOF_VARIABLES, IELEM(N, ICONSIDERED)%IFCA * NUM_DG_DOFS)::RHS_DG_RECONSTRUCT ! See eq. 3.19 of Luo 2012
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS)::LUO_LSQ_RECONSTRUCT_SOL
    
    LHS_MATRIX = ZERO
    
    DO I_FACE = 1, IELEM(N, ICONSIDERED)%IFCA ! Construct LHS, RHS matrices
        NEIGHBOR_INDEX = IELEM(N,ICONSIDERED)%INEIGH(I_FACE)
        OFFSET = NUM_DG_DOFS * (I_FACE - 1) + 1 ! Incrementing LHS and RHS matrix indices for each face
    
        IF (IELEM(N, ICONSIDERED)%INTERIOR == 0)THEN ! Element is interior
            NEIGHBOR_TYPE = 0
        ELSE ! Element is not interior
            IF (IELEM(N,ICONSIDERED)%INEIGHB(I_FACE) == N) THEN ! Same CPU
                IF (IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE).GT.0) THEN ! Boundary
                    IF (IBOUND(N,IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE))%ICODE.EQ.5) THEN ! PERIODIC IN MY CPU
                        NEIGHBOR_TYPE = 0
                    ELSE
                        NEIGHBOR_TYPE = 1
                    END IF
                ELSE ! Not boundary
                    NEIGHBOR_TYPE = 0
                END IF
            ELSE ! Other CPU
                IF (IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE).GT.0) THEN ! Boundary
                    IF (IBOUND(N,IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE))%ICODE.EQ.5) THEN ! PERIODIC IN OTHER CPU
                        NEIGHBOR_TYPE = 2
                    END IF
                ELSE ! Not boundary
                    NEIGHBOR_TYPE = 2
                END IF
            END IF
        END IF
        
        SELECT CASE(NEIGHBOR_TYPE)
        CASE(0) ! Interior or same CPU (periodic or not boundary)
            X1 = IELEM(N, NEIGHBOR_INDEX)%XXC
            Y1 = IELEM(N, NEIGHBOR_INDEX)%YYC
            BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER + 1, NEIGHBOR_INDEX, NUM_DG_RECONSTRUCT_DOFS - 1)
            NEIGHBOR_U = U_C(NEIGHBOR_INDEX)%VALDG(1,:,:)
            NEIGHBOR_DELTA_XYZ = IELEM(N, NEIGHBOR_INDEX)%DELTA_XYZ
            DG_SOL_TEMP = DG_SOL_I_CENTER(N, NEIGHBOR_INDEX) ! To call dg_sol without changing ICONSIDERED
        CASE(1) ! Same CPU boundary and not periodic
        CASE(2) ! Other CPU
            ! BASIS_TEMP (basis at center of neighboring element)
            BASIS_TEMP = IEXSOLHIR(ICONSIDERED)%BASIS_NEIGHBOR_CENTER(I_FACE, :)
            ! DELTA_XYZ
            NEIGHBOR_DELTA_XYZ = IEXSOLHIR(ICONSIDERED)%DELTA_XYZ(I_FACE, :)
            ! NEIGHBOR_U
            NEIGHBOR_U = IEXSOLHIR(ICONSIDERED)%SOL_DG(I_FACE,1:NOF_VARIABLES,1:NUM_DG_DOFS)
            DO I_VAR = 1, NOF_VARIABLES
                DG_SOL_TEMP(I_VAR) = NEIGHBOR_U(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP, NEIGHBOR_U(I_VAR, 2:))
            END DO
        END SELECT
        
        DO I_DG_RECONSTRUCT_DOF = 1, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS
            LHS_MATRIX(:, OFFSET, I_DG_RECONSTRUCT_DOF) = BASIS_TEMP(NUMBER_OF_DOG + I_DG_RECONSTRUCT_DOF) ! Higher order terms
        END DO
        
        ! 2D, first order terms
        LHS_MATRIX(:,OFFSET+1,1) = BASIS_TEMP(1) ! x, x term
        LHS_MATRIX(:,OFFSET+1,DIMENSIONA+1) = BASIS_TEMP(2) ! x, xy term
        LHS_MATRIX(:,OFFSET+2,2) = BASIS_TEMP(2) ! y, y term
        LHS_MATRIX(:,OFFSET+2,DIMENSIONA+1) = BASIS_TEMP(1) ! y, xy term
        
        DO I_VAR = 1, NOF_VARIABLES
            DG_SOL_I_J(I_VAR) = U_C(ICONSIDERED)%VALDG(1, I_VAR, 1) + DOT_PRODUCT(BASIS_TEMP, U_C(ICONSIDERED)%VALDG(1, I_VAR, 2:)) ! DG solution using values at i and basis at center of j
        END DO
        
        DG_SOL_TEMP = DG_SOL_TEMP - DG_SOL_I_J
        DO I_VAR = 1, NOF_VARIABLES
            RHS_DG_RECONSTRUCT(I_VAR,OFFSET) =  DG_SOL_TEMP(I_VAR) ! 1st term
            DO I_DIM = 1, DIMENSIONA ! 1st order terms
                RHS_DG_RECONSTRUCT(I_VAR, OFFSET+I_DIM) = IELEM(N, ICONSIDERED)%DELTA_XYZ(I_DIM) / NEIGHBOR_DELTA_XYZ(I_DIM) * NEIGHBOR_U(I_VAR,I_DIM+1) - U_C(ICONSIDERED)%VALDG(1,I_VAR,I_DIM+1)
            END DO
        END DO
    END DO
    
    DO I_VAR = 1, NOF_VARIABLES
        LUO_LSQ_RECONSTRUCT_SOL(I_VAR, :) = NORMAL_LSQ(LHS_MATRIX(I_VAR,:,:), RHS_DG_RECONSTRUCT(I_VAR,:)) ! Solve for reconstructed solution
    END DO
    
    WRITE(700+N,*) LHS_MATRIX, 'RHS_DG_RECONSTRUCT', RHS_DG_RECONSTRUCT, 'LUO_LSQ_RECONSTRUCT_SOL', LUO_LSQ_RECONSTRUCT_SOL
    
END FUNCTION LUO_LSQ_RECONSTRUCT_SOL

END MODULE DG_FUNCTIONS
