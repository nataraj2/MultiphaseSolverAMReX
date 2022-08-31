module fileio
  use precision
  implicit none
  
  ! File index
  integer :: fileindex
  integer, dimension(128) :: iunits
  
contains
  
  ! ====================== !
  ! File index management: !
  !   - open a file        !
  !   - add it to the list !
  ! ====================== !
  integer function iopen()
    implicit none

    integer, save :: icall = 1
    integer :: i

    if (icall .eq. 1) then
       fileindex = 1
       icall = 0
    end if
    iunits(fileindex) = 0
    do i=1,fileindex
       if (iunits(i) .eq. 0) exit
    end do
    if (i .eq. fileindex) then
       fileindex = fileindex + 1
       if (fileindex .ge. 128) stop "iopen: maximum units number exceeded"
    end if
    iunits(i) = 1
    iopen = i + 10
    return
  end function iopen
  
  ! ======================= !
  ! File index management:  !
  !   - close a file        !
  !   - remove it from list !
  ! ======================= !
  integer function iclose(iu)
    implicit none

    integer :: iu

    iu = iu - 10
    if (iu .gt. 0 .and. iu .lt. fileindex) then
       iunits(iu) = 0
       iclose = iu + 10
    else
       iclose = -1
    end if
    return
  end function iclose
  
end module fileio
