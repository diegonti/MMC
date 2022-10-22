function molarMass(labels)
    !Returns molarMass of a given array containing the atom labels of the molecule
    implicit none

    !Constants: Molar masses (g/mol)
    real,parameter :: mC=12.011
    real,parameter :: mH=1.00784
    real,parameter :: mO=15.999
    real,parameter :: mN=14.0067

    integer :: nC, nH, nO, nN
    character(len=1), dimension(100) :: labels
    real :: molarMass

    nC = count(labels=="C")
    nH = count(labels=="H")
    nO = count(labels=="O")
    nN = count(labels=="N")

    molarMass = nC*mC+nH*mH+nO*mO+nN*mN

end function molarMass