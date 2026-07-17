!> Global registry for restart read/write callbacks used by statistics modules.
module modrestart_registry
  use modglobal, only: output_prefix, startfile, lwarmstart
  use modmpi, only: cmyid, myid

  implicit none
  private

  integer, parameter :: max_entries = 128
  integer, parameter :: fmt_version = 1
  character(len=*), parameter :: magic = 'DALES_STAT_RESTART'

  abstract interface
    subroutine restart_writer(iunit)
      integer, intent(in) :: iunit
    end subroutine restart_writer

    subroutine restart_reader(iunit)
      integer, intent(in) :: iunit
    end subroutine restart_reader
  end interface

  type restart_entry
    character(len=64) :: name = ''
    procedure(restart_writer), pointer, nopass :: writer => null()
    procedure(restart_reader), pointer, nopass :: reader => null()
  end type restart_entry

  type(restart_entry), save :: entries(max_entries)
  integer, save :: nentries = 0

  public :: register_restart_handlers
  public :: run_restart_writers
  public :: run_restart_readers
  public :: run_cleanups

contains

  subroutine register_restart_handlers(name, reader, writer)
    character(len=*), intent(in) :: name
    procedure(restart_reader) :: reader
    procedure(restart_writer) :: writer

    integer :: idx

    idx = find_entry(trim(name))
    if (idx == 0) then
      if (nentries >= max_entries) return
      nentries = nentries + 1
      idx = nentries
      entries(idx)%name = trim(name)
    end if

    entries(idx)%reader => reader
    entries(idx)%writer => writer
  end subroutine register_restart_handlers

  subroutine run_restart_writers(restart_name)
    character(len=*), intent(in) :: restart_name

    integer :: i, iunit, nwriters
    character(len=50) :: state_name, latest_name

    state_name = restart_name
    if (len_trim(state_name) >= 5) state_name(5:5) = 'r'

        open(newunit=iunit, file=trim(output_prefix)//trim(state_name), form='unformatted', &
          status='replace', action='write')

    nwriters = 0
    do i = 1, nentries
      if (associated(entries(i)%writer)) nwriters = nwriters + 1
    end do

    write(iunit) magic, fmt_version, nwriters

    do i = 1, nentries
      if (.not. associated(entries(i)%writer)) cycle
      write(iunit) entries(i)%name
      call entries(i)%writer(iunit)
    end do

    close(iunit)

    latest_name = state_name
    if (len_trim(latest_name) >= 13) then
      latest_name(6:13) = '_latest_'
      call system('ln -s -f '//trim(state_name)//' '//trim(output_prefix)//trim(latest_name))
    end if
  end subroutine run_restart_writers

  subroutine run_restart_readers
    integer :: i, iunit, nstored, version
    character(len=50) :: state_name
    character(len=64) :: name
    character(len=len(magic)) :: file_magic

    if (.not. lwarmstart) return

    state_name = startfile
    if (len_trim(state_name) >= 5) state_name(5:5) = 'r'
    if (len_trim(state_name) >= 21) state_name(14:21) = cmyid

        open(newunit=iunit, file=trim(output_prefix)//trim(state_name), form='unformatted', &
          status='old', action='read')

    read(iunit) file_magic, version, nstored
    if (file_magic /= magic .or. version /= fmt_version) then
      close(iunit)
      return
    end if

    do i = 1, nstored
      read(iunit) name
      call dispatch_reader(trim(name), iunit)
    end do

    close(iunit)

    if (myid == 0) then
      write(*,*) 'Restarted statistics state from ', trim(state_name)
    end if
  end subroutine run_restart_readers

  subroutine run_cleanups(restart_name)
    character(len=*), intent(in) :: restart_name

    integer :: i, iunit, nwriters
    character(len=50) :: state_name, latest_name

    state_name = restart_name
    if (len_trim(state_name) >= 5) state_name(5:5) = 'r'

        open(newunit=iunit, file=trim(output_prefix)//trim(state_name), form='unformatted', &
          status='replace', action='write')

    nwriters = 0
    do i = 1, nentries
      if (associated(entries(i)%writer)) nwriters = nwriters + 1
    end do

    write(iunit) magic, fmt_version, nwriters

    do i = nentries, 1, -1
      if (.not. associated(entries(i)%writer)) cycle
      write(iunit) entries(i)%name
      call entries(i)%writer(iunit)
    end do

    close(iunit)

    latest_name = state_name
    if (len_trim(latest_name) >= 13) then
      latest_name(6:13) = '_latest_'
      call system('ln -s -f '//trim(state_name)//' '//trim(output_prefix)//trim(latest_name))
    end if
  end subroutine run_cleanups

  integer function find_entry(name)
    character(len=*), intent(in) :: name
    integer :: i

    find_entry = 0
    do i = 1, nentries
      if (trim(entries(i)%name) == trim(name)) then
        find_entry = i
        return
      end if
    end do
  end function find_entry

  subroutine dispatch_reader(name, iunit)
    character(len=*), intent(in) :: name
    integer, intent(in) :: iunit

    integer :: idx

    idx = find_entry(name)
    if (idx == 0) return
    if (.not. associated(entries(idx)%reader)) return

    call entries(idx)%reader(iunit)
  end subroutine dispatch_reader

end module modrestart_registry
