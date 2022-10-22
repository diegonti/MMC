function selectGroup(distances,label,a1,a2, dmin,dmax, n)
    !Returns number of ocurrences in array of a selected group of two atoms (a1,a2) and a range (dmin,dmax)
    integer :: n
    real, dimension(n,n) :: distances
    character(len=1), dimension(100) :: label
    character(len=1) :: a1,a2
    real :: dmin,dmax
    integer :: selectGroup

    selectGroup = 0
    do i=1,n
        do j=1,n

            if ((label(i)==a1) .AND. (label(j)==a2)) then
                if ((distances(i,j) >= dmin) .AND. (distances(i,j) <= dmax)) then
                   selectGroup = selectGroup + 1
                end if
            end if
        end do
    end do

end function selectGroup