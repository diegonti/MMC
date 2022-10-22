program molecules
    ! With an .xyz file from user input, opens the file and returns
    ! Molar mass and geometric information (max and min distances, hydroxyl and carbonyl groups)
    implicit none

    !Variable declaration
    real, dimension(100) :: x,y,z               !Coordinates arrays
    character(len=1), dimension(100) :: label   !Atom Label array
    character(len=100) :: path, molec           !Files
    integer :: n,i,j,selectGroup                !Number of atoms and indexes
    real :: molarMass, r, maxDist, minDist, maxDistNoH, minDistNoH          !Distances
    integer, dimension(2) :: maxLocat, minLocat, maxLocatNoH, minLocatNoH   !Location
    character(len=2) :: maxL1, maxL2, minL1, minL2, maxL1NoH, maxL2NoH, minL1NoH, minL2NoH !Format
    real, allocatable, dimension(:,:) :: distances,distancesNoH       !Matrix with distances

    !Main Program
    print*, "Input the molecule filename from the molecules folder (Ex: paracetamol.xyz):"
    read(5,*) molec                     !Read user input
    path = "./molecules/"//molec        !Gets path of molecule
    print *, "path: ",path

    open(1, file = path, status="old")
        read(1,*) n                 !Number of atoms
        read(1,*)                   !Energy
        allocate(distances(n,n))    !Allocate matrix with n
        allocate(distancesNoH(n,n))
        do i = 1,n          !Get Coordinates
            read(1,*) label(i), x(i), y(i), z(i)
        end do
    close(1)

    !Loop to get distances
    do i=1,n
        do j=1,n
            if (i/=j) then    !For avoiding ii distances (atom with itself)
                r = sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
                distances(i,j) = r
                distancesNoH(i,j) = r
                if ((label(i)=="H") .OR. (label(j)=="H")) then !to dismiss hydrogen atoms
                    distancesNoH(i,j) = 2
                end if
            else
                distances(i,j) = 2 !Adding 2 as a mean value. If leaved empty, 0 or a random number would interfere with max and min results.
                distancesNoH(i,j) = 2
            end if
        end do
    end do

    !Get Max and Min Distances and Location (with Hydrogens)
    maxDist = maxval(distances)
    maxLocat = maxloc(distances)
    write(maxL1, "(I2)") maxLocat(1)
    write(maxL2, "(I2)") maxLocat(2)

    minDist = minval(distances)
    minLocat = minloc(distances)
    write(minL1, "(I2)") minLocat(1)
    write(minL2, "(I2)") minLocat(2)

    !Get Max and Min Distances and Location (without Hydrogens)
    maxDistNoH = maxval(distancesNoH)
    maxLocatNoH = maxloc(distancesNoH)
    write(maxL1NoH, "(I2)") maxLocatNoH(1)
    write(maxL2NoH, "(I2)") maxLocatNoH(2)

    minDistNoH = minval(distancesNoH)
    minLocatNoH = minloc(distancesNoH)
    write(minL1NoH, "(I2)") minLocatNoH(1)
    write(minL2NoH, "(I2)") minLocatNoH(2)

    !Printing Results
    print*
    print*, "Molar Mass (g/mol): ", molarMass(label)
    print*
    print*,"Considering H atoms..."
    print*,"Max distance (Ang): ",maxDist,"Between atoms: ",label(maxLocat(1)),maxL1,"-",label(maxLocat(2)),maxL2
    print*,"Min distance (Ang): ",minDist,"Between atoms: ",label(minLocat(1)),minL1,"-",label(minLocat(2)),minL2
    print*
    print*,"NOT considering H atoms..."
    print*,"Max distance (Ang): ",maxDistNoH,"Between atoms: ",label(maxLocatNoH(1)),maxL1NoH,"-",label(maxLocatNoH(2)),maxL2NoH
    print*,"Min distance (Ang): ",minDistNoH,"Between atoms: ",label(minLocatNoH(1)),minL1NoH,"-",label(minLocatNoH(2)),minL2NoH
    print*
    print*,"Number of OH (hydroxyl):", selectGroup(distances,label,"O","H",0.94,1.05,n)
    print*,"Number of CO (carbonyl):", selectGroup(distances,label,"C","O",1.20,1.30,n)
    print*

    deallocate(distances)
    deallocate(distancesNoH)
end program molecules
