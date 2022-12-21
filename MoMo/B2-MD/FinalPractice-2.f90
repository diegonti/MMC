program main

    ! Variable declaration
    implicit none
    real, parameter :: pi = 4.*atan(1.)
    real, parameter :: Na = 6.0221408e23
    real, parameter :: kb = 1.380649e-23
    real :: L,dens,a, E, cutoff, dt, Ekin,Epot, T,s,nu,Tinst,P,pt_mod,msd
    real :: ru_time,ru_dens,ru_dist,ru_temp,ru_E,ru_press,ru_vel,ru_mom, eps,sigma,MM,mass
    real, allocatable,dimension(:,:) :: r,rold,v,F
    real, dimension(3) :: pt_vec
    character(len=10) :: cell
    character(len=25) :: fileTraj,fileData,file1,file2,fileV1,fileV2
    integer :: M, N, i,j,k, Nsteps, frames_to_RDF
    logical :: addpbc

    ! Initial parameters/inputs                             (216)         (256)             (250)
    N = 256                ! Choose the right number, N=M**3 for sc, N=4M**3 for fcc, N=2M**3 for bcc
    dens = 0.8             ! (0.05; 0.1; 0.2; 0.4; 0.6; 0.8)
    cell = "fcc"           ! Specify unit cell (sc, fcc, bcc)

    sigma = 3.4     ! Amstrongs
    eps = 0.998     ! kJ/mol
    MM = 40         ! g/mol
    mass = MM/Na    ! g

    L = (N/dens)**(1./3.)
    if (cell=="sc") then; M = (N)**(1./3.)
    else if (cell=="fcc") then; M = (N/4.)**(1./3.)
    else if (cell=="bcc") then; M = (N/2.)**(1./3.)
    end if

    a = L/M              ! Lattice parameter
    dt = 0.0001          ! timestep
    cutoff = 0.4*L       ! Cutoff for interactions
    nu = 0.1             ! Thermostat nu
    addpbc = .true.      ! To add or not PBC
    frames_to_RDF = 1000 ! Frames to save positions to compute later RDF

    ! Name of files the data is outputed
    file1 = "initial.xyz"                   ! Initial coordinates
    file2 = "final.xyz"                     ! Final coordinates
    fileV1 = "initialV.dat"                 ! Initial velocities
    fileV2 = "finalV.dat"                   ! Final velocities
    fileData = "thermodynamics_raw.dat"     ! Observables data
    fileTraj = "trajectory.xyz"             ! Saved trajectories

    ! Conversion factors between reduced and real units
    ru_time = sqrt(mass*Na*(sigma*1e-10)**2 / (eps*1e6))    ! t in seconds
    ru_dist = sigma                                         ! distance in Ang
    ru_dens = 1e24 * MM / (Na*sigma**3)                     ! density in g/mL
    ru_E = eps                                              ! energy in kJ/mol
    ru_press = 1e3*eps/((sigma*1e-10)**3 * Na)              ! pressure in Pa
    ru_temp = 1e3*eps/(kb*Na)                               ! temperature in K
    ru_vel = ru_dist/ru_time *1e-10                         ! velocity in Ang/s
    ru_mom = m*ru_vel

    ! Print parameters
    print*, "Input parameters (reduced units)"
    print*, "Unit Cell:           ", cell
    print*, "Number of particles: ", N
    print*, "Density:             ", dens
    print*, "Box length:          ", L
    print*, "Lattice parameter:   ", a
    print*

    print*, "Input parameters (real units)"
    print*, "Unit Cell:           ", cell
    print*, "Number of particles: ", N
    print*, "Density:             ", dens*ru_dens
    print*, "Box length:          ", L*ru_dist
    print*, "Lattice parameter:   ", a*ru_dist
    print*


    allocate(r(N,3))
    allocate(v(N,3))
    allocate(F(N,3))
    allocate(rold(N,3))

    ! Opening output files
    open(1,file=file1,status="replace")
    open(2,file=file2,status="replace")
    open(3,file=fileData,status="replace")
    write(3,"(8a15)") "# Time","Epot","Ekin","E","Momentum","Tinst","Pressure","MSD"
    open(4,file=fileV1,status="replace")
    open(5,file=fileV2,status="replace")
    open(8,file=fileTraj,status="replace")

    ! Initializing positions (r) and velocities (v)
    call initializePositions(N,dens,r,cell)
    v = 0

    call writePOS(r*ru_dist,1)     ! Saves initial geometry
    call writePOS(v*ru_vel,4)      ! Saves initial velocities


    ! Main MD loop to disorder the configuration !!!!!!!!!!!!!!!!!!
    Nsteps = 10000      ! Number of disordering steps
    T = 100             ! Reduced temperature (melting)
    Tinst = T
    s = sqrt(T)         ! Sigma (themrostat)
    print*, "Starting Melting at T = ",T*ru_temp
    do i=1,Nsteps
        call velocityVerlet(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
!        call euler(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
!        call verlet(r,rold,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)

        call thermostat(v,nu,s)
!        call pbcN(r,L)

!        call momentum(v,pt_vec,pt_mod)
!        Tinst = 2.*Ekin/(3.*N)
!        write(3,*) i*dt,Epot,Ekin, Epot+Ekin, pt_mod, Tinst,P
    end do


    ! Main MD production loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Nsteps = 500000     ! Number of simulation steps
    T = 2.0             ! Reduced temperature
    Tinst = T
    s = sqrt(T)         ! Sigma (thermostat)
    rold = r            ! Saving initial coordinates (MSD)

    print*, "Starting Production Cycles at T = ",T*ru_temp
    write(*,"(a)",advance='no') "Completed (%): "
    do i=1,Nsteps
        call velocityVerlet(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
!        call euler(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
!        call verlet(r,rold,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)

        call thermostat(v,nu,s)
!        call pbcN(r,L)

        call momentum(v,pt_vec,pt_mod)
        Tinst = 2.*Ekin/(3.*N)
        call getMSD(r,rold,L,msd)

        write(3,*) i*dt*ru_time,Epot*ru_E,Ekin*ru_E, (Epot+Ekin)*ru_E, pt_mod*ru_mom, Tinst*ru_temp,P*ru_press,msd*ru_dist**2

        ! Saves last snapshots of trajectoy (to compute RDF correctly)
        if ((i>=(Nsteps-frames_to_RDF)) .and. (i<=Nsteps)) then
            call writePOS(r*ru_dist,8)
        end if

        if (mod(i,Nsteps/10)==0) then
            write(*,'(1x,i0)',advance='no') (100*i)/Nsteps
        end if
    end do
    print*,

    call writePOS(r*ru_dist,2)  ! Saves final position
    call writePOS(v*ru_vel,5)  ! Saves final velocities


    deallocate(r)
    deallocate(v)
    deallocate(F)
    deallocate(rold)



    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!      Subroutines      !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine boxmuller(S, X1,X2, XOUT1, XOUT2)
        implicit none
        real, intent(in) :: x1,x2, s
        real, intent(out) :: xout1, xout2
        real pi
        pi = 4.*atan(1.)

        XOUT1=s*sqrt(-2.*(log(1.-x1)))*cos(2.*pi*x2)
        XOUT2=s*sqrt(-2.*(log(1.-x1)))*sin(2.*pi*x2)

        end subroutine boxmuller

        ! Andersen Thermostat
        subroutine thermostat(v,nu,s)
        implicit none
        real, intent(inout),dimension(:,:) :: v
        real, intent(in) :: nu,s
        real :: x1,x2,xout1,xout2
        integer :: N,i,d

        N = size(v,dim=1)
        do i=1,N
            if (rand()<nu) then
                do d=1,3
                    x1=rand();x2=rand()
                    call boxmuller(s,x1,x2,xout1,xout2)
                    v(i,d) = xout1
                end do
            end if
        end do
        end subroutine


        ! To choose and run the specified initial position
        subroutine initializePositions(N,dens,r,cell)
        implicit none
        real, intent(inout),dimension(:,:) :: r
        integer, intent(in) :: N
        real, intent(in) :: dens
        character(len=10), intent(in) :: cell

        if ((index(cell,"fcc")==1) .OR. (index(cell,"FCC")==1))  then
            call initFCC(N,dens,r)
        else if((index(cell,"sc")==1) .OR. (index(cell,"SC")==1)) then
            call initSC(N,dens,r)
        else if ((index(cell,"bcc")==1) .OR. (index(cell,"BCC")==1))  then
            call initBCC(N,dens,r)
        end if

        end subroutine initializePositions

        ! Initial position and velocities of two particles
        subroutine init2(L,r,v)
        implicit none
        real, intent(in) :: L
        real, intent(inout),dimension(:,:) :: r,v

        r = reshape( (/L/2.-0.75,0.,0., L/2+0.75,0.,0./),shape(r), order=(/2,1/))
        v = reshape( (/0.,0.,0., 0.,0.,0./),shape(v), order=(/2,1/))

        end subroutine

        ! Initial positions for SC cell
        subroutine initSC(N,dens,r)
        implicit none
        real, intent(inout),dimension(:,:) :: r
        integer, intent(in) :: N
        real, intent(in) :: dens
        real :: L,a
        integer :: p, i,j,k, M

        L = (N/dens)**(1./3.)
        M = N**(1./3.)
        a = L/M
        p = 1
        do i=0,M-1
            do j=0,M-1
                do k=0,M-1
                    r(p, :) = [i,j,k]
                    p=p+1
                end do
            end do
        end do
        r = r*a
        end subroutine initSC

        ! Initial positions for FCC cell
        subroutine initFCC(N,dens,r)
        implicit none
        real, intent(inout),dimension(:,:) :: r
        integer, intent(in) :: N
        real, intent(in) :: dens
        real, dimension(4,3) :: ucell
        real :: L,a
        integer :: i,j,k, p,at,M

        L = (N/dens)**(1./3.)
        M = (N/4.)**(1./3.)
        a = L/M

        ucell = reshape( (/0.,0.,0., 0.,0.5,0.5,  0.5,0.,0.5,  0.5,0.5,0./),shape(ucell), order=(/2,1/))
        p = 1
        do i=0,M-1
            do j=0,M-1
                do k=0,M-1
                    do at=1,size(ucell, dim=1)
                        r(p, :) = ucell(at,:) + (/i,j,k/)
                        p=p+1
                    end do
                end do
            end do
        end do
        r = r*a

        end subroutine initFCC

        ! Initial positions for BCC cell
        subroutine initBCC(N,dens,r)
        implicit none
        real, intent(inout),dimension(:,:) :: r
        integer, intent(in) :: N
        real, intent(in) :: dens
        real, dimension(2,3) :: ucell
        real :: L,a
        integer :: i,j,k, p,at,M

        L = (N/dens)**(1./3.)
        M = (N/2.)**(1./3.)
        a = L/M

        ucell = reshape( (/0.,0.,0., 0.5,0.5,0.5/),shape(ucell), order=(/2,1/))
        p = 1
        do i=0,M-1
            do j=0,M-1
                do k=0,M-1
                    do at=1,size(ucell, dim=1)
                        r(p, :) = ucell(at,:) + (/i,j,k/)
                        p=p+1
                    end do
                end do
            end do
        end do
        r = r*a
        end subroutine initBCC

        ! Initial velocities as Bimodal distribution
        subroutine initBimodal(T,v)
        implicit none
        real,intent(inout),dimension(:,:) :: v
        real,intent(in) :: T
        real :: vi
        integer :: N
        N = size(v, dim=1)
        v=0
        if (mod(N,2)/=0) then
            print*, "Number of particles (N) should be multiple of 2."
        end if
        vi = sqrt(T)
        v(1:N:2,:) = -vi
        v(2:N:2,:) = +vi

        end subroutine initBimodal

        ! Initialize velocities at zero
        subroutine initZero(v)
        implicit none
        real,intent(inout),dimension(:,:) :: v
        v = 0
        end subroutine initZero


        ! Writes current positions in the specified file (indexed)
        subroutine writePOS(r,fileN)
        implicit none
        real, intent(in),dimension(:,:) :: r
        integer, intent(in) :: fileN
        integer :: i, N

        N = size(r, dim=1)

        write(fileN,*) N
        write(fileN,*)
        do i= 1,N
            write(fileN,*) "Ar", r(i,:)
        end do

        end subroutine writePOS


        ! Computes truncated Lenard-Jones Energy and Forces with PBC
        subroutine getForces(r,cutoff,dens,Tinst,addpbc,E,F,P)
        implicit none
        real, intent(in),dimension(:,:) :: r
        real, intent(in) :: cutoff,dens,Tinst
        logical, intent(in) :: addpbc
        real, intent(out) :: E,P
        real, intent(inout), dimension(:,:) :: F
        real :: d2,d6,d8,d12,d14, cf2,cf6,cf12, L,V
        real, dimension(3) :: rij,Fij
        integer :: N,i,j

        N = size(r, dim=1)
        L = (N/dens)**(1./3.)
        cf2 = cutoff*cutoff
        F = 0
        E = 0
        P = (dens*Tinst)
        V = L**3
        do i = 1,N
            do j=i+1,N
                rij = r(i,:) - r(j,:)

                if (addpbc) then
                    call pbc(rij,L)
                end if


                d2 = sum(rij**2)
                if (d2<cf2) then
                    d6=d2*d2*d2; d12=d6*d6
                    d8=d6*d2; d14=d8*d6
                    cf6=cf2*cf2*cf2; cf12=cf6*cf6

                    Fij = (48./d14 - 24./d8)*rij

                    ! F(i) = -F(j)
                    F(i,:) = F(i,:) + Fij
                    F(j,:) = F(j,:) - Fij

                    !Energy: E = E(r) - E(cf) --> smoother truncation
                    E = E + 4.*(1./d12 - 1./d6) - 4.*(1./cf12 - 1./cf6)

                    ! Pressure
                    P = P + (1./(3.*V))*(dot_product(rij,Fij))

                end if
            end do
        end do
        end subroutine getForces


        ! Periodic Boundary Conditions for 1 particle
        subroutine pbc(rij,L)
        implicit none
        real, intent(inout),dimension(3) :: rij
        real, intent(in) :: L
        integer :: d,i

        do d=1,3
            if (rij(d)>(L/2.)) then
                rij(d) = rij(d) - L
            else if (rij(d)<-(L/2.)) then
                rij(d) = rij(d) + L
            end if
        end do

        end subroutine pbc

        ! do pbc pero <L <0 ??

        ! Periodic Boundary Conditions for a matrix
        subroutine pbcN(r,L)
        implicit none
        real, intent(inout),dimension(:,:) :: r
        real, intent(in) :: L
        integer :: d,i

        do i=1,size(r,dim=1)
            do d=1,3
                if (r(i,d)>(L)) then
                    r(i,d) = r(i,d) - L
                else if (r(i,d)<0) then
                    r(i,d) = r(i,d) + L
                end if
            end do
        end do

        end subroutine pbcN

        ! Implementation of the Verlet integrator
        subroutine verlet(r,rold,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
        implicit none
        logical, intent(in) :: addpbc
        real, intent(in) :: dt, cutoff, dens,Tinst
        real, intent(inout), dimension(:,:) :: F
        real, intent(inout),dimension(:,:) :: r, rold,v
        real, intent(out) :: Ekin,Epot,P
        real, dimension(:,:), allocatable :: roldtemp

        allocate(roldtemp(size(r, dim=1),3))

        call getForces(r,cutoff,dens,Tinst,addpbc,Epot,F,P)
        roldtemp = r
        r = 2*r - rold + F*dt*dt
        call pbcN(r,L)
        v = (r - rold)/(2*dt)
        rold = roldtemp

        call kinetc(v,Ekin)

        deallocate(roldtemp)

        end subroutine verlet

        ! Implementation of the Velocity-Verlet integrator
        subroutine velocityVerlet(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
        implicit none
        logical, intent(in) :: addpbc
        real, intent(in) :: dt, cutoff, dens,Tinst
        real, intent(inout), dimension(:,:) :: F
        real, intent(inout),dimension(:,:) :: r, v
        real,intent(out) :: Ekin,Epot,P

        call getForces(r,cutoff,dens,Tinst,addpbc,Epot,F,P)
        r = r + v*dt + 0.5*F* dt*dt
        call pbcN(r,L)
        v = v + 0.5*F*dt
        call getForces(r,cutoff,dens,Tinst,addpbc,Epot,F,P)
        v = v + 0.5*F*dt
        call kinetc(v,Ekin)

        end subroutine velocityVerlet


        ! Implementation of the Euler integrator
        subroutine euler(r,v,dt,F,Ekin,Epot,P,cutoff,dens,Tinst,addpbc)
        implicit none
        logical, intent(in) :: addpbc
        real, intent(in) :: dt, cutoff, dens,Tinst
        real, intent(inout), dimension(:,:) :: F
        real, intent(inout),dimension(:,:) :: r, v
        real,intent(out) :: Ekin,Epot,P

        call getForces(r,cutoff,dens,Tinst,addpbc,Epot,F,P)

        r = r + v*dt + 0.5*F*dt*dt
        v = v + F*dt
        call pbcN(r,L)
        call kinetc(v,Ekin)

        end subroutine euler

        ! Kinetic Energy Calculation
        subroutine kinetc(v,Ekin)
        implicit none
        real, intent(in),dimension(:,:) :: v
        real, intent(out) :: Ekin
        integer :: N,i

        Ekin = 0
        N = size(v, dim=1)
        do i=1,N
            Ekin = Ekin + sum((v(i,:)**2)/2.)
        end do

        end subroutine kinetc

        ! Compute the total momentum
        subroutine momentum(v,pt_vec,pt_mod)
        implicit none
        real, intent(in),dimension(:,:) :: v
        real, intent(out),dimension(3) :: pt_vec
        real, intent(out) :: pt_mod
        integer :: N,i
        N = size(v, dim=1)

        pt_vec = 0
        do i=1,N
            pt_vec = pt_vec + v(i,:)
        end do
        pt_mod = sqrt(sum(pt_vec**2))

        end subroutine momentum

        subroutine getMSD(r,rold,L,msd)
        implicit none
        real, intent(in) :: L
        real, intent(in),dimension(:,:) :: r,rold
        real,intent(out) :: msd
        real, dimension(3) :: rij
        integer :: N, i
        N = size(r, dim=1)
        msd = 0.
        do i=1,N
            rij = r(i,:) - rold(i,:)
            call pbc(rij,L)
            msd = msd + sum(rij**2)
        end do
        msd = msd/N
        end subroutine


        ! Test function to print a variable array (r,v,f,...)
        subroutine test(var)
        implicit none
        real, intent(in), dimension(:,:) :: var
        integer :: i

        do i=1,size(var,dim=1)
            print*, var(i,:)
        end do

        end subroutine test


end program main


