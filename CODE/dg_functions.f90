MODULE DG_FUNCTIONS

USE BASIS
USE DECLARATION
USE DERIVATIVES

IMPLICIT NONE

CONTAINS


! FUNCTION DG_RHS_INTEGRAL(N, I_ELEM, QP_X, QP_Y, QP_WEIGHT, NUM_VARS, ORDER, NUM_DOFS, CELL_VOL_OR_SURF, FLUX_TERM, VOL_OR_SURF)
! !> @brief
! !> Calculates the volume or surface integral term in the DG RHS for scalar linear advection with speed = 1
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::N, I_ELEM, ORDER, NUM_VARS, NUM_DOFS, VOL_OR_SURF
!     REAL,INTENT(IN)::QP_X, QP_Y, QP_WEIGHT, CELL_VOL_OR_SURF
!     REAL,DIMENSION(:),INTENT(IN)::FLUX_TERM
!     INTEGER::I_VAR
!     REAL,DIMENSION(NUM_DOFS+1,NUM_VARS)::DG_RHS_INTEGRAL
!     
!     IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
!         WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
!         STOP
!     END IF
!     
!     IF (VOL_OR_SURF == 1) THEN ! VOLUME INTEGRAL
!         DO I_VAR = 1, NUM_VARS
!             DG_RHS_INTEGRAL(1,I_VAR) = 0
!             DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * (BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,1) + BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,2)) ! For linear advection with speed = 1
!         END DO
!     ELSE IF (VOL_OR_SURF == 2) THEN ! SURFACE INTEGRAL
!         DO I_VAR = 1, NUM_VARS
!             DG_RHS_INTEGRAL(1,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF
!             DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * BASIS_REC2D(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS)
!         END DO
!     END IF
! 
! END FUNCTION DG_RHS_INTEGRAL


FUNCTION DG_SOL(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
INTEGER,INTENT(IN)::N
INTEGER::I_DOF, I_VAR
REAL,DIMENSION(Nof_VARIABLES)::DG_SOL

!     IF(ALL(SHAPE(U_C_VALDG) /= (/ NUM_VARIABLES, NUM_DOFS+1 /))) THEN
!         WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
!         STOP
!     END IF
NUMBER=IELEM(N,ICONSIDERED)%IORDER
BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)

DO I_VAR = 1, NOF_VARIABLES
    DG_SOL(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:NUMBER_OF_DOG), U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1))
END DO

END FUNCTION DG_SOL


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

X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
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
INTEGER::I
INTEGER,INTENT(IN)::N

X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
NUMBER=IELEM(N,ICONSIDERED)%IORDER

DO I = 1, NOF_VARIABLES
    DG_SURF_FLUX(1,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)
    DG_SURF_FLUX(2:NUMBER_OF_DOG+1,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
END DO

END FUNCTION DG_SURF_FLUX


FUNCTION DG_VOL_INTEGRAL(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(IDEGFREE+1,NOF_VARIABLES)::DG_VOL_INTEGRAL
REAL,DIMENSION(NOF_VARIABLES)::FLUX_TERM
INTEGER::I,J,K,NQP,I_QP,I_VAR
REAL::PH,INTEG
    
!     IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
!         WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
!         STOP
!     END IF
    

IF (IELEM(N,ICONSIDERED)%ISHAPE == 5) NQP = QP_TRIANGLE*2
IF (IELEM(N,ICONSIDERED)%ISHAPE == 6) NQP = QP_TRIANGLE
    NUMBER=IELEM(N,ICONSIDERED)%IORDER
    NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE
           
DO I_VAR = 1, NOF_VARIABLES
    DG_VOL_INTEGRAL(1,I_VAR) = 0.0d0
    DG_VOL_INTEGRAL(2:,I_VAR)= 0.0D0

      DO I=1,NUMBER_OF_DOG
      
        DO I_QP = 1, NQP 
            X1=QP_ARRAY(ICONSIDERED,I_QP)%X
            Y1=QP_ARRAY(ICONSIDERED,I_QP)%Y
             
            FLUX_TERM=DG_SOL(N)
            
            IF (ITESTCASE.EQ.3)THEN
                LEFTV=DG_SOL(N)
                CALL PRIM2CONS2D(N)
                FLUX_TERM=FLUXEVAL2D(LEFTV)
             
            END IF
             
              DG_VOL_INTEGRAL(I+1,I_VAR) = DG_VOL_INTEGRAL(I+1,I_VAR)+FLUX_TERM(I_VAR)* QP_ARRAY(ICONSIDERED,I_QP)%QP_WEIGHT *(DF2DX(X1,Y1,I)+DF2DY(X1,Y1,I))!* IELEM(N,ICONSIDERED)%TOTVOLUME

!               DG_VOL_INTEGRAL(2:,I_VAR) = DG_VOL_INTEGRAL(2:,I_VAR)+FLUX_TERM(I_VAR)* QP_ARRAY(ICONSIDERED,I_QP)%QP_WEIGHT * IELEM(N,ICONSIDERED)%TOTVOLUME*(BASIS_REC2D_DERIVATIVE(N,x1,y1,number,ICONSIDERED,NUMBER_OF_DOG,1) + BASIS_REC2D_DERIVATIVE(N,x1,y1,number,ICONSIDERED,NUMBER_OF_DOG,2))
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
        
            write(700+n,*)"look here","cell",i_elem,"face",facex,"point",pointx
            write(700+n,*)ILOCAL_RECON3(I_ELEM)%ULEFT_DG(1:NOF_VARIABLES, I_FACE, I_QP)
                
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

SUBROUTINE PRESTORE_AND_ALLOCATE_DG
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY, SURF_QPOINTS, mass matrix
    IMPLICIT NONE
    INTEGER::I, K, I_QP, N_QP, I_FACE,nnd,iqp,idummy,COUNT_1
    
	ALLOCATE(QP_ARRAY(XMPIELRANK(N),QP_Triangle*2)); !Allocates for 2D
    
    DO I = 1, XMPIELRANK(N)    
        iconsidered=i
        !Store volume quadrature points
        ELTYPE=IELEM(N,I)%ISHAPE
        DO K = 1,IELEM(N,I)%NONODES
            NODES_LIST(k,1:2)=INODER(IELEM(N,I)%NODES(K))%CORD(1:2)
            VEXT(k,1:2)=NODES_LIST(k,1:2)
        END DO
        
        CALL DECOMPOSE2
        
        !Store delta xyz (normalization factor from Luo 2012)
        IELEM(N,I)%DELTA_XYZ = CALC_DELTA_XYZ(IELEM(N,I)%NONODES, DIMENSIONA, NODES_LIST)
    
        SELECT CASE(ielem(n,i)%ishape)
        CASE(5)
            COUNT_1=0
             DO K=1,ELEM_DEC
                 VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            
                CALL QUADRATUREtriangle(N,IGQRULES)
!                CALL QUADRATUREQUAD(N,IGQRULES)
                
                VOLTEMP=TRIANGLEVOLUME(N)
                
                DO I_QP = 1, QP_Triangle
                    COUNT_1=COUNT_1+1
                    QP_ARRAY(I,COUNT_1)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
                    QP_ARRAY(I,COUNT_1)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
                    QP_ARRAY(I,COUNT_1)%QP_WEIGHT = WEQUA3D(I_QP) * VOLTEMP
                END DO
                
           END DO
            
!             DO I_QP = 1, QP_QUAD!N_QP
!                 QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
!                 QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
!                 QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP) * QUADVOLUME(N)
!             END DO
            
        CASE(6)
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            N_QP = QP_Triangle
            
            DO I_QP = 1, N_QP
                QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
                QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
                QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP) * IELEM(N,I)%totvolume
            END DO
        END SELECT
                

        
        ALLOCATE(ILOCAL_RECON3(I)%SURF_QPOINTS(IELEM(N,I)%IFCA, QP_LINE_N, DIMENSIONA))
        !Store surface quadrature points
        DO I_FACE = 1, IELEM(N,I)%IFCA
            IDUMMY=0
            !GAUSSIAN POINTS FIXED
                if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                IF (IELEM(N,I)%IBOUNDS(I_FACE).GT.0)THEN	!CHECK FOR BOUNDARIES
                    if (ibound(n,ielem(n,i)%ibounds(I_FACE))%icode.eq.5)then	!PERIODIC
                        IDUMMY=1
                    END IF
                END IF
        
                IQP=QP_LINE
                NND=2
                IF (IDUMMY.EQ.0)THEN
                    DO K=1,NND
                        VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(I_FACE,K))%CORD(1:dims)
                    END DO
                ELSE
                    facex=I_FACE;
                    CALL coordinates_face_PERIOD2D1(n,iconsidered,facex)
                    
                    
                END IF
                CALL QUADRATURELINE(N,IGQRULES)	  
            ELSE
                IQP=QP_LINE
                NND=2
                DO K=1,NND
                    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(I_FACE,K))%CORD(1:dims)
                END DO
                CALL QUADRATURELINE(N,IGQRULES)
            END IF
        
        
        
        
        
        
        
        
        
!             VEXT(1,1:2) = inoder(IELEM(N,I)%NODES_FACES(I_FACE,1))%CORD(1:2)  !COPY THE COORDINATE OF THE FIRST NODE OF THID EDGE
!             VEXT(2,1:2) = inoder(IELEM(N,I)%NODES_FACES(I_FACE,2))%CORD(1:2)  !COPY THE COORDINATE OF THE SECOND NODE OF THID EDGE
!             CALL QUADRATURELINE(N,IGQRULES)

            DO I_QP = 1, QP_LINE_N
                ILOCAL_RECON3(I)%SURF_QPOINTS(I_FACE,I_QP,1) = QPOINTS2D(1,I_QP) - IELEM(N,I)%XXC
                ILOCAL_RECON3(I)%SURF_QPOINTS(I_FACE,I_QP,2) = QPOINTS2D(2,I_QP) - IELEM(N,I)%YYC
            END DO
        END DO
        
        ALLOCATE(ILOCAL_RECON3(I)%ULEFT_DG(NOF_VARIABLES, IELEM(N,I)%IFCA, QP_LINE_N))
    END DO
    
    CALL ASS_MASS_MATRIX(N)
    
END SUBROUTINE PRESTORE_AND_ALLOCATE_DG


SUBROUTINE ASS_MASS_MATRIX(N)
!> @brief
!> Assembles the mass matrix
!> REQUIRES: Globals: IELEM, QP_QUAD, QP_TRIANGLE, MASS_MATRIX
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_QP, N_QP, I_DOF, J_DOF, KMAXE
    REAL,DIMENSION(IDEGFREE)::BASIS_VECTOR
    REAL::INTEG_TEST,PHX,INTEG_MM
    
    KMAXE = XMPIELRANK(N)
    
    ALLOCATE(MASS_MATRIX_CENTERS(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); MASS_MATRIX_CENTERS(N,:,:,:) = ZERO; MASS_MATRIX_CENTERS(N,:,1,1) = 0.0D0
    
    DO I_ELEM = 1, KMAXE
        SELECT CASE(IELEM(N,I_ELEM)%ISHAPE)
        CASE(5) ! Quadrilateral
            N_QP = QP_TRIANGLE*2
        CASE(6) ! Triangle
            N_QP = QP_TRIANGLE
        END SELECT
        
    NUMBER_OF_DOG = IELEM(N,I_ELEM)%IDEGFREE
        
    DO I_DOF = 1, NUM_DG_DOFS       
        DO J_DOF = 1, NUM_DG_DOFS
            INTEG_MM = ZERO
            
                DO I_QP = 1, N_QP
                    IXX = I_ELEM; x1=QP_ARRAY(I_ELEM,I_QP)%X; y1=QP_ARRAY(I_ELEM,I_QP)%Y
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
                    INTEG_MM = INTEG_MM + PHX*QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
                END DO
        
        MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF, J_DOF) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF, J_DOF) + INTEG_MM
        
        END DO
    END DO

   
END DO


    ALLOCATE(INV_MASS_MATRIX(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); INV_MASS_MATRIX(N,:,:,:) = ZERO
    
    CALL COMPMASSINV(MASS_MATRIX_CENTERS(N,:,:,:), INV_MASS_MATRIX(N,:,:,:), NUM_DG_DOFS)
    
    DO I_ELEM = 1, KMAXE
        WRITE(300+N,*) I_ELEM
!         DO I_QP = 1, N_QP
!             WRITE(300+N,*) 'QP_ARRAY', QP_ARRAY(I_ELEM,I_QP)%X, QP_ARRAY(I_ELEM,I_QP)%Y, QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!             WRITE(300+N,*) 'BASIS,', BASIS_REC2D(N,QP_ARRAY(I_ELEM,I_QP)%X,QP_ARRAY(I_ELEM,I_QP)%Y,IORDER,I_ELEM,IDEGFREE)
!         END DO
        WRITE(300+N,*) 'XYZ', IELEM(N,I_ELEM)%NODES
!         WRITE(300+N,*) 'DELTAXYZ', IELEM(N,I_ELEM)%DELTA_XYZ
        WRITE(300+N,*) 'MMC', MASS_MATRIX_CENTERS(N,I_ELEM,:,:)
        WRITE(300+N,*) 'Inverse,', INV_MASS_MATRIX(N,I_ELEM,:,:)
        WRITE(300+N,*) 'Identity', MATMUL(MASS_MATRIX_CENTERS(N,I_ELEM,:,:),INV_MASS_MATRIX(N,I_ELEM,:,:))
    END DO
    
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

FUNCTION INVERSE_MATRIX(MATRIX_IN, N_COLS)
    INTEGER,INTENT(IN)::N_COLS
    REAL,DIMENSION(N_COLS,N_COLS),INTENT(IN)::MATRIX_IN
    REAL,DIMENSION(1,N_COLS,N_COLS)::DUMMY_MATRIX_IN, DUMMY_INVERSE_MATRIX
    REAL,DIMENSION(N_COLS,N_COLS)::INVERSE_MATRIX
    
    DUMMY_MATRIX_IN(1,:,:) = MATRIX_IN
    CALL COMPMASSINV(DUMMY_MATRIX_IN, DUMMY_INVERSE_MATRIX, N_COLS)
    INVERSE_MATRIX = DUMMY_INVERSE_MATRIX(1,:,:)
    
END FUNCTION INVERSE_MATRIX

FUNCTION NORMAL_LSQ(A, b)
    REAL,DIMENSION(:,:),INTENT(IN)::A
    REAL,DIMENSION(:),INTENT(IN)::b
    REAL,DIMENSION(SIZE(A,1))::NORMAL_LSQ
    
    NORMAL_LSQ = MATMUL(MATMUL(INVERSE_MATRIX(MATMUL(TRANSPOSE(A), A), SIZE(A,2)), TRANSPOSE(A)), b)
    
END FUNCTION NORMAL_LSQ

! FUNCTION LUO_LSQ_RECONSTRUCT(N, N_DIM, ORDER, DEGFREE)
! !> @brief
! !> This function reconstructs an approximation of NUM_DG_RECONSTRUCT_DOFS
! !> REQUIRES: IELEM, U_C as globals
!     INTEGER,INTENT(IN)::N, N_DIM, ORDER, DEGFREE
!     INTEGER::I_ELEM, I_FACE, NEIGHBOR_INDEX, I_DIM, KMAXE
!     REAL,DIMENSION(DEGFREE)::BASIS_TEMP
!     REAL,DIMENSION(IELEM(N, I_ELEM)%IFCA * 3, NUM_DG_RECONSTRUCT_DOFS)::LHS_MATRIX ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(IELEM(N, I_ELEM)%IFCA * 3)::RHS_DG_RECONSTRUCT ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(XMPIELRANK(N),NUM_DG_RECONSTRUCT_DOFS)::LUO_LSQ_RECONSTRUCT ! NUM_DG_RECONSTRUCT_DOFS
!     
!     KMAXE = XMPIELRANK(N)
!     
!     DO I_ELEM = 1, KMAXE
!         LHS_MATRIX = ZERO
!         IF (IELEM(N, I_ELEM)%INTERIOR == 0)THEN ! Element is interior
!             DO I_FACE = 1, IELEM(N, I_ELEM)%IFCA
!                 NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(I_FACE)
!                 OFFSET = (N_DIM + 1) * (I_FACE - 1)
!             
!                 BASIS_TEMP = BASIS_REC2D(N,IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, IORDER,I_ELEM,IDEGFREE)
!                 
!                 DO I_DG_RECONSTRUCT_DOF = 1, NUM_DG_RECONSTRUCT_DOFS
!                     LHS_MATRIX(OFFSET+1,I_DG_RECONSTRUCT_DOF) = BASIS_TEMP(N_DIM+I_DG_RECONSTRUCT_DOF-1)
!                 END DO
!                 
!                 LHS_MATRIX(OFFSET+2,1) = BASIS_TEMP(1)
!                 LHS_MATRIX(OFFSET+2,N_DIM+1) = BASIS_TEMP(2)
!                 LHS_MATRIX(OFFSET+3,2) = BASIS_TEMP(2)
!                 LHS_MATRIX(OFFSET+3,N_DIM+1) = BASIS_TEMP(1)
!                 
!                 RHS_DG_RECONSTRUCT(OFFSET+1) = DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, ORDER, DEGFREE, U_C(NEIGHBOR_INDEX)%VALDG(1,:,:)) - DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, ORDER, DEGFREE, U_C(I_ELEM)%VALDG(1,:,:))
!                 DO I_DIM = 1, N_DIM
!                     RHS_DG_RECONSTRUCT(OFFSET+I_DIM+1) = IELEM(N,I_ELEM)%DELTA_XYZ(I_DIM) / IELEM(N,NEIGHBOR_INDEX)%DELTA_XYZ(I_DIM) * U_C(I_ELEM)%VALDG(1,:,I_DIM+1) - U_C(NEIGHBOR_INDEX)%VALDG(1,:,I_DIM+1)
!                 END DO
!                 
!             END DO
!         END IF
!         
!         LUO_LSQ_RECONSTRUCT(I_ELEM,:) = NORMAL_LSQ(LHS_MATRIX, RHS_DG_RECONSTRUCT)
!         
!     END DO
!     
! END FUNCTION LUO_LSQ_RECONSTRUCT

END MODULE DG_FUNCTIONS
