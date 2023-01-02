program practice3

    ! Declaration of variables
    implicit none
    real :: T, r1279, t_start, t_finish, tMC_start, tMC_finish
    integer :: dE,E,M
    integer(kind=8) :: accepted, not_accepted, Nsteps
    integer :: L,N,z, nMCS,nmeas, i,rnd,time
    character(len=:),allocatable :: file_name
    integer, allocatable, dimension(:) :: s
    integer, allocatable, dimension(:,:) :: nbr, shift
    real, allocatable, dimension(:) :: table

    call cpu_time(t_start)

    ! INPUTS
    L = 100                     ! 2D Lattice size
    N = L**2                    ! Number of spins
    z = 4                       ! Coordination number (2D square)
    T = 2.00                    ! Temperature
    nMCS = 1e6                  ! Number of MC Steps
    Nsteps = N*nMCS             ! Total steps (# of attempted flips)
    nmeas = 10                  ! Measure each nmeas MCS
    call setr1279(3333)         ! Randomizer seed

    ! Allocating memory of arrays
    allocate(nbr(z,N))
    allocate(shift(2,L))
    allocate(s(N))
    allocate(table(-2*z:+2*z))

    ! Parametrized table of exponentials exp(-dE/T)
    call set_table(z,table,T)


    call init_spins(s)              ! Initializing spin array at random

    call pbc(shift,L)               ! Applying PBC
    call init_nbr(nbr,shift)        ! Setting neighbors

    call energy(s,nbr,E)            ! Computing initial Energy
    call magnetization(s,M)         ! Computing initial Magnetization

    file_name = "out.dat"           ! Output file
    open(1,file=trim(file_name),status="replace")
    write(1,*) "# L = ", L
    write(1,"(3a11)") "# Time","E","M"

    ! Main MonteCarlo loop
    call cpu_time(tMC_start)
    accepted = 0_8; not_accepted = 0_8
    do time=1,nMCS
        do i=1,N
            rnd = mod(int(N*r1279()),N)+1    ! Generating random number (index)
            call Ediff(rnd,s,nbr,dE)         ! Computing energy difference

            ! Acceptance of step and updating system
            if (dE < 0) then
                E = E + dE              ! Updating Energy
                M = M - 2*s(rnd)        ! Updating Magnetization
                s(rnd) = -s(rnd)        ! Flipping Spin
                accepted = accepted + 1_8

            else
                if ((dE>0) .and. (r1279()<table(dE))) then
                    E = E + dE              ! Updating Energy
                    M = M - 2*s(rnd)        ! Updating Magnetization
                    s(rnd) = -s(rnd)        ! Flipping Spin
                    accepted = accepted + 1_8
                else
                    not_accepted = not_accepted + 1_8
                end if
            end if
        end do

        ! Writing data into output file
        if (mod(time,nmeas) == 0) then
            write(1,*) time,E,M
        end if
    end do
    call cpu_time(tMC_finish)

    ! Printing
    print*,"Number of spins:   ", N
    print*,"Last Energy:       ", E
    print*,"Last Magnetization:", M
    print*,"Accepted steps     ", accepted
    print*,"Total steps (N*nMCS)",Nsteps

    print*,"MonteCarlo time(attempted flips/s):", Nsteps/(tMC_finish-tMC_start)


    ! Deallocating memory of arrays
    deallocate(nbr)
    deallocate(shift)
    deallocate(s)
    deallocate(table)

    call cpu_time(t_finish)
    print*
    print*, "Process finished. Time(s): ", t_finish-t_start


    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!      Subroutines      !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialization of neighbors matrix for a square lattice
    subroutine init_nbr(nbr,shift)
    implicit none
    integer, intent(in), dimension(:,:) :: shift
    integer, intent(inout), dimension(:,:) :: nbr
    integer :: L, i, x,y
    L = size(shift,dim=2)

    i = 0
    do y=1,L
        do x=1,L
            i = i + 1
            nbr(1,i) = shift(2,x) + L*(y-1)
            nbr(2,i) = shift(1,x) + L*(y-1)
            nbr(3,i) = x + L*(shift(2,y)-1)
            nbr(4,i) = x + L*(shift(1,y)-1)
        end do
    end do

    end subroutine init_nbr


    ! Applies PBC (shift)
    subroutine pbc(shift,L)
    implicit none
    integer, intent(in) :: L
    integer, intent(inout), dimension(:,:) :: shift
    integer :: i

    do i=1,L
        shift(1,i) = i-1
        shift(2,i) = i+1
    end do
    shift(1,1) = L
    shift(2,L) = 1

    end subroutine pbc


    ! Creates table of energies
    subroutine set_table(z,table,T)
    implicit none
    real, intent(in) :: T
    real,intent(inout),dimension(-2*z:2*z) :: table
    integer :: n,z

    do n=-2*z,2*z
        table(n) = exp(-n/T)
    end do

    end subroutine set_table


    ! Initialize randomly the spin array
    subroutine init_spins(s)
    implicit none
    integer, intent(inout), dimension(:) :: s
    real :: r1279
    integer :: N,i
    N = size(s)

    do i=1,N
        s(i) = 2*mod(int(2*r1279()),2) - 1
    end do

    end subroutine init_spins


    ! Compute energy difference
    subroutine Ediff(rnd_i,s,nbr,dE)
    implicit none
    integer,intent(in) :: rnd_i
    integer, intent(in), dimension(:) :: s
    integer, intent(in), dimension(:,:) :: nbr
    integer, intent(out) :: dE
    integer :: h,k

    z = size(nbr, dim=1)
    h=0
    do k=1,z
        h = h + s(nbr(k,rnd_i))
    end do
    dE = 2*s(rnd_i)*h

    end subroutine


    ! Compute total energy
    subroutine energy(s,nbr,E)
    implicit none
    integer, intent(in), dimension(:) :: s
    integer, intent(in), dimension(:,:) :: nbr
    integer, intent(out) :: E
    integer :: E1,E2
    integer :: N,z, i,k

    z = size(nbr, dim=1)
    N = size(nbr, dim=2)

    E1 = 0
    do i=1,N
        E2 = 0
        do k=1,z
            E2 = E2 + s(nbr(k,i))
        end do
        E1 = E1 + s(i)*E2
    end do

    E = -E1/2

    end subroutine energy


    ! Compute magnetization
    subroutine magnetization(s,M)
    implicit none
    integer, intent(in), dimension(:) :: s
    integer, intent(out) :: M

    M = sum(s)

    end subroutine magnetization


    ! Test print for 2D arrays
    subroutine test(array)
    implicit none
    integer, intent(in),dimension(:,:) :: array
    integer :: N,i
    N = size(array, dim=1)

    do i=1,N
        print*, array(i,:)
    end do

    end subroutine test


end program practice3
