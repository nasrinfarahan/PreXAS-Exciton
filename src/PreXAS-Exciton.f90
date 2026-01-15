!===============================================================================
! PreXAS-Exciton: Pump-probe XAS from valence/core BSE eigenstates
! Computes x-ray absorption cross section of optically excited material
! Based on methodology in Phys. Rev. B 110, 235126 (2024)
!===============================================================================
program PreXAS

  implicit none
  ! Mathematical constants
  complex, parameter :: z = (0.0d0, 1.0d0)
  
  ! Physical parameters (modify in input.cfg or here)
  double precision, parameter :: h = 0.489805d0  ! Broadening Γ (eV)
  double precision, parameter :: m = 0.489805d0  
  double precision :: xI = 5.06833921667d0      ! Pump energy ω_o (eV)
  double precision :: xI_min = 4.1, xI_max = 5.1, stepI = 0.1d0
  double precision :: xF_min = 257.0d0, xF_max = 266.0d0, stepF = 0.01d0
  
  ! Polarizations (a-axis example: modify as needed)
  ! eps stands for polarization of x-ray pulse in x,y,z direction
  ! eps_opt  stands for polarization of optical pulse in x,y,z direction
  double precision :: eps_opt_x = 0.5d0, eps_opt_y = -0.866025404d0, eps_opt_z = 0.0d0
  double precision :: eps_x = 0.5d0, eps_y = -0.866025404d0, eps_z = 0.0d0
  
  ! Dimensions (from exciting calculations)
  integer, parameter :: lambdaI_max = 30, lambdaF_max = 30
  integer, parameter :: conduction_max = 22, valence_max = 16, core_max = 2, kpoints = 64. ,conduction_min = 17
  
  ! Arrays for BSE eigenvectors (Re/Im parts)
  real*8, dimension(lambdaI_max, conduction_min:conduction_max, 1:valence_max, kpoints, 2) :: B
  real*8, dimension(lambdaF_max, conduction_min:conduction_max, 1:core_max, kpoints, 2) :: X
  
  ! Momentum matrix elements
  complex(8), allocatable :: pmcov(:,:,:,:), pmvv(:,:,:,:)
  
  ! Signal and weights
  complex(8), allocatable :: ReSignal(:,:), ImSignal(:,:)
  complex(8), allocatable :: Reweight(:), Imweight(:)
  
  ! K-points
  real*8, allocatable :: kpoint(:,:)
  
  ! Working variables
  ! ic = conduction, iv= valence, iy= core, index= kpoints number
  integer :: i, j, k, ic, iv, iy, index, io, lambdaI, lambdaF
  integer :: conduction, valence, core
  ! kcx to kcz is a brief defination of kpoints in three dimension x, y,z in conduction bands 
  ! and respectively the same for kvx to kvz for valence bands and core states
  double precision :: kcx, kcy, kcz, kvx, kvy, kvz
  ! Re and Im means repectively real part and imaginary part of a complex number
  ! BSE means valence-excitations and Xray means core-excitations
  ! ene is brief of energy corresponds to BSE and Xray
  double precision :: ReBSE, ImBSE, ReXray, ImXray, eneBSE, eneXray, deltaE
  ! RePhi, ImPhi, is the real and imaginary part of momentum matrix elemetns based on basis functions (we have it in x,y,z direction)
  ! this momentum matrix element is between core to valence
  double precision :: RePhi_x, ImPhi_x, RePhi_y, ImPhi_y, RePhi_z, ImPhi_z
  ! reFe, imFe, is the real and imaginary part of momentum matrix elemetns based on basis functions (we have it in x,y,z direction)
  ! this momentum matrix element is between valence to conduction
  double precision :: reFe_x, imFe_x, reFe_y, imFe_y, reFe_z, imFe_z
  ! real and imaginary part of multipling BSE and Xray complex number, give us Recoeff, Imcoeff, coeff is a brief of coefficient
  complex*8 :: Recoeff, Imcoeff
  double precision :: sigma, xF
  INTEGER :: rows, numberOfRows, numberOfRowz, rowz, numberofExcitonicWave, cloumns

  double precision :: dummy ! we define dummy to avoid reading useless coloumns in read files
  integer :: dummy_int, rank
  integer(kind=8) :: recl
  integer, allocatable :: index2rank(:)
  integer :: Bunit, Xunit, Funit
  character(500) :: dummyString, BSEfile, Xrayfile


 
!==============================================================================
!   111111     PPPPPP   U    U   M     M   PPPPPP
!   11         P    P   U    U   MM   MM   P    P
!   11         PPPPPP   U    U   M M M M   PPPPPP
!   11         P        U    U   M  M  M   P
!   11  @      P        UUUUUU   M     M   P
!==============================================================================
    allocate(kpoint(1:3,1:kpoints)) !allocate kpoints in three diemension x, y, z and number of kpoints
    
! read the k-points file
    open(5,file='KPOINTS_QMT001_EXC.OUT', status='old', action='read')
! write it again to make sure read it correctly and avoid using unusable coloumns in previous file produceed by exciting
    open(unit = 7, file = 'kpoints.txt', status = 'old')
    DO k =1,kpoints
        read(5,*,iostat=io) dummy_int, kpoint(1,k),kpoint(2,k), kpoint(3,k), dummy, dummy, dummy, dummy, dummy
        write(7,*,iostat=io) k, kpoint(1,k), kpoint(2,k), kpoint(3,k)
    END DO
    rewind(7)
print*, 'Number of k-points and kpoint in each direction print out in kpoints.txt:', k
            
! read the BSE files for differents lambdaI
! read all files for BSE and try to skip 13 first lines which is not useful
DO i = 1, lambdaI_max
    10 format('PUMP/EXCITON_EVEC/EXEVEC_BSE-singlet-TDA-BAR_SCR-full_QMT001_LAMBDA',i6.6,'.OUT') 
    write(BSEfile,10)i
    Bunit = i + 520
    open(Bunit, file = BSEfile, action='read', status='old')

    if (i == 1) then
        Do rows = 1 , 3
            READ(Bunit,*)
        END DO
! read number of rows in BSE files from header
        READ(Bunit,'(A38,I8)') dummyString, numberOfRows
! skip lines 5 to 13 in file 1 and read line 4 to get the rows in the file
        Do rows = 5 , 13
            READ(Bunit,*)
        END DO
    else
        ! skip the first 13 lines in all other files
        Do rows = 1 , 13
            READ(Bunit,*)
        END DO
    ENDIF
END DO 

print*, 'Number of rows in BSE:', numberOfRows

! for each file in BSE we need to make k-points equal by kpoint file which we read it already (kpoints from files and kcx to kcz in conduction band 
! and kvx to kvz in valence band from BSE file must be equal to connect real and imaginary part of eigenvector to different lambda and kpoint number)
B = 0d0
DO i = 1, lambdaI_max
    DO j= 1, numberOfRows
        Bunit = i + 520
        READ(Bunit,*,iostat=io) dummy, lambdaI, ic, iv, kcx, kcy, kcz, kvx, kvy, kvz, dummy, ReBSE, ImBSE
        DO k = 1, kpoints
            if ((ABS(kpoint(1,k)-kcx).lt.1d-6).and.(ABS(kpoint(2,k)-kcy).lt.1d-6).and.(ABS(kpoint(3,k)-kcz).lt.1d-6))then
            index = k
            end if
        END DO 
        B(lambdaI,ic,iv,index,1) = ReBSE
        B(lambdaI,ic,iv,index,2) = ImBSE
    END DO
END DO 
print*, 'coincidence of conduction kpoint number with valence kpoint numbers in BSE'        

!==============================================================================
!   222222     PPPPPP   RRRRR    OOOOO   BBBBBB   EEEEEE
!       22     P    P   R    R   O   O   B    B   E
!   222222     PPPPPP   RRRRR    O   O   BBBBBB   EEEE
!   22         P        R   R    O   O   B    B   E
!   222222 @   P        R    R   OOOOO   BBBBBB   EEEEEE
!==============================================================================
! read X-ray files for differents lambdaF
do i = 1, lambdaF_max
    20 format('PROBE/EXCITON_EVEC/EXEVEC_BSE-singlet-TDA-BAR_SCR-full_QMT001_LAMBDA',i6.6,'.OUT') 
    write(Xrayfile,20)i
    Xunit = i + 200 
    open(Xunit, file = Xrayfile, action='read', status='old')

    if (i == 1) then
        ! skip first 13 lines in file 1 and read line 4 to get the number of rows in the file
        Do rowz = 1 , 3
            READ(Xunit,*)
        END DO

        READ(Xunit,'(A38,I8)') dummyString, numberOfRowz
! skip lines 5 to 13 in file 1 and read line 4 to get the rows in the file
        Do rowz = 5 , 13
            READ(Xunit,*)
        END DO
    else
! skip the first 13 lines in all other files
        Do rowz = 1 , 13
            READ(Xunit,*)
        END DO
    ENDIF
END DO 
print*, 'Number of rows in x-ray:', numberOfRowz

! for each file in BSE we need to make k-points equal by kpoint file which we read it already (kpoints from files and kcx to kcz in conduction band 
! and kvx to kvz in valence band from BSE file must be equal to connect real and imaginary part of eigenvector to different lambda and kpoint number)
X = 0d0
DO i = 1, lambdaF_max
    DO j= 1, numberOfRowz
        Xunit = i + 200
        READ(Xunit,*,iostat=io) dummy, lambdaF, ic, iy, kcx, kcy, kcz, kvx, kvy, kvz, dummy, ReXray, ImXray
        DO k = 1, kpoints
            if ((ABS(kpoint(1,k)-kcx).lt.1d-6).and.(ABS(kpoint(2,k)-kcy).lt.1d-6).and.(ABS(kpoint(3,k)-kcz).lt.1d-6))then
            index = k
            end if
        END DO 
        X(lambdaF,ic,iy,index,1) = ReXray
        X(lambdaF,ic,iy,index,2) = ImXray
    END DO
END DO 
print*, 'coincidence of conduction kpoint number with core kpoint numbers in Xray'               
!print out some real and imaginary part of X-ray to test code
print*,'A',  X(7,17,2,19,1),X(7,17,2,19,2)
print*,'A',  X(4,18,2,55,1),X(4,18,2,55,2)
print*,'A',  X(2,20,2,19,1),X(2,20,2,19,2) 
print*,'A',  X(1,18,1,19,1),X(1,18,1,19,2) 
print*,'A',  X(3,18,1,19,1),X(3,18,1,19,2) 

! SO far, I read the BSE eigenvectors after x-ray pulse and get them for the exact k-points 
! (cause we want to get those eigenvectors which excite from core to valence at the exact k-point)
!----------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------- 




!----------------------------------------------------------------------------------------
! \sum_{k,lambda_F, lambda_I, c, v, y} (A^{\lambda_F}_{ck,yk}*A^{\lambda_I}_{ck,vk})^2  
!----------------------------------------------------------------------------------------
!OPEN (4, file = 'NEWRESULT.OUT', status='replace', action='readwrite')


DO index = 1 , kpoints
    DO lambdaF = 1, lambdaF_max
        DO lambdaI = 1, lambdaI_max
        Recoeff= 0d0
        Imcoeff= 0d0
            DO ic = 17, conduction_max
                DO iv = 1, valence_max
                    DO iy = 1, core_max
                        Recoeff = Recoeff + (B(lambdaI,ic,iv,index,1)*X(lambdaF,ic,iy,index,1))+&
                        (B(lambdaI,ic,iv,index,2)*X(lambdaF,ic,iy,index,2))
                        Imcoeff = Imcoeff + (B(lambdaI,ic,iv,index,2)*X(lambdaF,ic,iy,index,1))-&
                        (B(lambdaI,ic,iv,index,1)*X(lambdaF,ic,iy,index,2))
                        !WRITE(4,*,iostat=io) index, lambdaI, lambdaF, iv, iy, ic, Recoeff + z*Imcoeff
                    END DO ! iy core
                END DO ! iv valence
            END DO ! ic conduction
        END DO ! number of lambda_I
    END DO ! number of lambda_F
END DO ! kpoints

print*, 'print multiply eigenvectors for different-\lambda'
 


!==============================================================================
!   333333     PPPPPP   M     M   AAAAA   TTTTTT
!       33     P    P   MM   MM   A   A     TT
!   333333     PPPPPP   M M M M   AAAAA     TT
!       33     P        M  M  M   A   A     TT
!   333333 @   P        M     M   A   A     TT
!------------------------------------------------------------------------------
!   \sum_{k, v, y} <\phi_vk|\phi_yk> momentum matrix element between core to valence
!==============================================================================


! try to read nonreadable momentum matrix element by exciting code that we modified the code to produce
allocate(index2rank(index))
index2rank = -1
DO index = 1, kpoints
    index2rank(index) = rank
END DO
        
allocate(pmcov(1:core_max,1:valence_max,3,1:kpoints))
DO index = 1, kpoints
    DO iy = 1, core_max
        DO iv = 1, valence_max
            if (rank == index2rank(index)) then
                inquire(iolength=recl) pmcov(iy,iv,:,index)
                open(700,File='PROBE/PMATCOV.OUT',Action='read',Form='UNFORMATTED',Access='DIRECT',Status='old',Recl=recl)
                read(700,rec=index) pmcov(iy,iv,:,index)
            end if !rank
        END DO !core
    END DO !valence
END DO !index



!------------------------------------------------------------------------------------------------------------------------------
! Signal = \sum_{K,lambda_F, lambda_I} \sum_{v, y} <\phi_vk|\phi_yk> * \sum_{c} A^{*\lambda_F}_{ck,yk}*A^{\lambda_I}_{ck,vk}) 
!------------------------------------------------------------------------------------------------------------------------------
! allocate imaginary and real part of produced signal
allocate(ReSignal(1:lambdaF_max,1:lambdaI_max)) 
allocate(ImSignal(1:lambdaF_max,1:lambdaI_max)) 

!open files related Final energy and initial energy to calculate differences (E_F - E_I)(ev)
OPEN (1, file = 'BSE.energy.OUT', status='old', action='read') ! Energy of excite an elcetron after optical pulse
OPEN (2, file = 'X-ray.energy.OUT', status='old', action='read') ! Energy of excite an elcetron after x-ray pulse
OPEN (9, file = 'Signalnew.OUT', status='replace', action='readwrite') ! open a file to write the signal 

DO lambdaF = 1, lambdaF_max
READ(2,*,end=5000,iostat=io) eneXray
    DO lambdaI = 1, lambdaI_max
    READ(1,*,end=6000,iostat=io) eneBSE
    deltaE = eneXray-eneBSE 
    ImSignal(lambdaF,lambdaI) = 0d0
    ReSignal(lambdaF,lambdaI) = 0d0      
        DO index = 1 , kpoints
            DO iy = 1, core_max
                DO iv = 1, valence_max
                    if (rank == index2rank(index)) then
                        inquire(iolength=recl) pmcov(iy,iv,:,index)
                        open(3,File='pmcov.OUT',Action='readwrite',Form='FORMATTED',Access='STREAM',Status='old')
                        write(3,*) pmcov(iy,iv,:,index)
                        !REAL(REAL(pmcov(iy,iv,1,index))), REAL(REAL(pmcov(iy,iv,2,index))), &
                        !REAL(REAL(pmcov(iy,iv,3,index))), REAL(AIMAG(pmcov(iy,iv,1,index))),&
                        !REAL(AIMAG(pmcov(iy,iv,2,index))), REAL(AIMAG(pmcov(iy,iv,3,index)))
                    end if !rank
                    RePhi_x = ((REAL(REAL(pmcov(iy,iv,1,index))))*(eps_x))
                    ImPhi_x = ((REAL(AIMAG(pmcov(iy,iv,1,index))))*(eps_x))
                    RePhi_y = ((REAL(REAL(pmcov(iy,iv,2,index))))*(eps_y))
                    ImPhi_y = ((REAL(AIMAG(pmcov(iy,iv,2,index))))*(eps_y))
                    RePhi_z = ((REAL(REAL(pmcov(iy,iv,3,index))))*(eps_z))
                    ImPhi_z = ((REAL(AIMAG(pmcov(iy,iv,3,index))))*(eps_z))
                    Recoeff= 0d0
                    Imcoeff= 0d0
                    DO ic = 17, conduction_max
                        Recoeff = Recoeff + (B(lambdaI,ic,iv,index,1)*X(lambdaF,ic,iy,index,1))+&
                        (B(lambdaI,ic,iv,index,2)*X(lambdaF,ic,iy,index,2))
                        Imcoeff = Imcoeff + (B(lambdaI,ic,iv,index,2)*X(lambdaF,ic,iy,index,1))-&
                        (B(lambdaI,ic,iv,index,1)*X(lambdaF,ic,iy,index,2))
                        ReSignal(lambdaF,lambdaI) = (ReSignal(lambdaF,lambdaI))+((REAL(Recoeff)*(RePhi_x+RePhi_y+RePhi_z))-&
                        ((REAL(Imcoeff)*(ImPhi_x+ImPhi_y+ImPhi_z))))
                        ImSignal(lambdaF,lambdaI) = (ImSignal(lambdaF,lambdaI))+((REAL(Recoeff)*(ImPhi_x+ImPhi_y+ImPhi_z))+&
                        ((REAL(Imcoeff)*(RePhi_x+RePhi_y+RePhi_z))))
                    END DO !ic conduction band 
                END DO !iy core
            END DO !iv valence
        END DO ! index --> k-points
        WRITE(9,*,iostat=io) lambdaF, lambdaI, deltaE, eneBSE, &
        ((abs(ReSignal(lambdaF,lambdaI)+ z*ImSignal(lambdaF,lambdaI)))**2d0)             
    END DO !lambda_I
    6000 continue
    rewind(1)
END DO !lambda_F
5000 continue
rewind(2)
print*, 'print out the signal'
!----------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------- 





!-----------------------------------------------------------------------------------------
! \sum_{k, c, v} <\phi_ck|\phi_vk>  momentum matrix element between valence to conduction
!-----------------------------------------------------------------------------------------   
allocate(pmvv(1:valence_max,17:conduction_max,3,1:kpoints))
DO index = 1, kpoints
    DO iv = 1, valence_max
        DO ic = 17 , conduction_max
            if (rank == index2rank(index)) then
                inquire(iolength=recl) pmvv(iv,ic,:,index)
                open(70,File='PUMP/PMATVV.OUT',Action='read',Form='UNFORMATTED',Access='DIRECT',Status='old',Recl=recl)
                read(70,rec=index) pmvv(iv,ic,:,index)
            end if !rank
        END DO !valence
    END DO ! conduction
END DO !index

! Calculate weigth for the signal by multiplying the coefficient into momentum matrix element
OPEN (90, file = 'Weight.OUT', status='replace', action='readwrite') !open file to calculate weight


!lambdaI_max = 120

allocate(Reweight(1:lambdaI_max)) 
allocate(Imweight(1:lambdaI_max)) 

! try to produce a sigma which is lorentzian broadening over the signal

!xI = 5.06833921667 !this attribute you change to get different signal in different excitation energies
xI = xI_min
Do while (xI < xI_max + stepI)
xI = xI + stepI
    DO lambdaI = 1, lambdaI_max
    READ(1,*,end=60000,iostat=io) eneBSE
    Imweight(lambdaI) = 0d0
    Reweight(lambdaI) = 0d0      
        DO index = 1 , kpoints
            DO iv = 1, valence_max
                DO ic = 17 , conduction_max
                    if (rank == index2rank(index)) then
                        inquire(iolength=recl) pmvv(iv,ic,:,index)
                        open(300,File='pmvv.OUT',Action='readwrite',Form='FORMATTED',Access='STREAM',Status='old')
                        write(300,*) pmvv(iv,ic,:,index)
                        !REAL(REAL(pmvv(iv,ic,1,index))), REAL(REAL(pmvv(iv,ic,2,index))), &
                        !REAL(REAL(pmvv(iv,ic,3,index))), REAL(AIMAG(pmvv(iv,ic,1,index))),&
                        !REAL(AIMAG(pmvv(iv,ic,2,index))), REAL(AIMAG(pmvv(iv,ic,3,index)))
                    end if !rank
                    reFe_x = ((REAL(REAL(pmvv(iv,ic,1,index))))*(eps_opt_x))
                    imFe_x = ((REAL(AIMAG(pmvv(iv,ic,1,index))))*(eps_opt_x))
                    reFe_y = ((REAL(REAL(pmvv(iv,ic,2,index))))*(eps_opt_y))
                    imFe_y = ((REAL(AIMAG(pmvv(iv,ic,2,index))))*(eps_opt_y))
                    reFe_z = ((REAL(REAL(pmvv(iv,ic,3,index))))*(eps_opt_z))
                    imFe_z = ((REAL(AIMAG(pmvv(iv,ic,3,index))))*(eps_opt_z))
                        Reweight(lambdaI) = Reweight(lambdaI)+((B(lambdaI,ic,iv,index,1)*(reFe_x+reFe_y+reFe_z))-&
                        ((B(lambdaI,ic,iv,index,2)*(imFe_x+imFe_y+imFe_z))))
                        Imweight(lambdaI) = Imweight(lambdaI)+((REAL(B(lambdaI,ic,iv,index,1))*(imFe_x+imFe_y+imFe_z))+&
                        ((REAL(B(lambdaI,ic,iv,index,2))*(reFe_x+reFe_y+reFe_z)))) 
                END DO !iv valence
            END DO !ic conduction band
        END DO ! index --> k-points
        WRITE(90,*,iostat=io) lambdaI, eneBSE, ((abs(Reweight(lambdaI)+ z*Imweight(lambdaI)))**2d0)
        !(abs(h-z*(xI-(eneBSE))))**2d0)        
    END DO !lambda_I
    60000 continue
END DO
rewind(10)
print*, 'weight'


!==============================================================================
!   4444444     SSSSSS   IIIIII   GGGGGG   M     M   AAAAA
!   44   44     S          II     G        MM   MM   A   A
!   44444444    SSSSSS     II     G  GGG   M M M M   AAAAA
!        44          S     II     G    G   M  M  M   A   A
!        44 @   SSSSSS   IIIIII   GGGGGG   M     M   A   A
!------------------------------------------------------------------------------
!   4. SIGMA
!==============================================================================
!FINAL LOOP
! So far we calculate the absolute value of weight and the signal 
! Now, we want to add the broadening for both weight and signal based on text file and multiply each other
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!Sigma = \sum_{lambda_F, lambda_I}  {|Weight(lambdaI,c,v,k)| * (1)/(\Gamma^2+ (\omega_I-eneBSE)^2)} * {|Signal(lambda_F,lambda_I,v,y,k)| *(1)/(\Gamma^2+ (\omega_F-DeltaE)^2)}   
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

OPEN (8, file = 'sigmanew.OUT', status='replace', action='readwrite') !open file for the main sigma file which calculate cross-section

! try to produce a sigma which is lorentzian broadening over the signal

xF = xF_min
Do while (xF < xF_max + stepF)
xF = xF + stepF
sigma = 0d0
    DO lambdaF = 1, lambdaF_max
    READ(2,*,end=7000,iostat=io) eneXray
        DO lambdaI = 1, lambdaI_max
        READ(1,*,end=8000,iostat=io) eneBSE
        deltaE = eneXray-eneBSE
        !sigma = sigma+((((abs(ReSignal(lambdaF,lambdaI)+z*ImSignal(lambdaF,lambdaI)))**2d0)/((h**2d0+((xF-deltaE)**2d0))))*&
        !(((abs(Reweight(lambdaI)+z*Imweight(lambdaI)))**2d0)/((h**2d0+(((xI-(eneBSE)))**2d0)))))
        sigma = sigma+((((abs(ReSignal(lambdaF,lambdaI)+z*ImSignal(lambdaF,lambdaI)))**2d0)/&
        (abs(h-z*(xF-deltaE)-z*(xI-(eneBSE))))**2d0)*&
        (((abs(Reweight(lambdaI)+z*Imweight(lambdaI)))**2d0)/(abs(m-z*(xI-(eneBSE))))**2d0))
        WRITE(8,*) xF, xI, sigma
        END DO !lambda_I
        8000 continue
        rewind(1)
    END DO !lambda_F
    7000 continue
    rewind(2)
END DO !xF x-ray energy






close(1)
close(2)
close(3) 
close(4)
close(5)
close(6)
close(7)
close(8)
close(10)
close(700)
    
    
DO i = 1, 30
    Bunit = i + 520
    close(Bunit)
END DO


do i = 1, 30
    Xunit = i + 10
    close(Xunit)
end do   
            
End Program PreXAS
