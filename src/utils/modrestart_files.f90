!> Global registry for restart read/write callbacks and tagged restart field I/O.
module modrestart_files
  use iso_fortran_env, only: real32
  use netcdf
  use modglobal, only: output_prefix, nsv, lrestart_netcdf
  use modlogging, only: finish
  use modprecision, only: field_r
  use modsurfdata, only: isurf

  use modrestart_writers, only: write_restart_field, read_restart_field, handle_netcdf_status

  implicit none
  private    
  integer, parameter :: max_entries = 128
  integer, parameter :: fmt_version = 1
  integer, parameter :: nc_name_len = 128
  integer, parameter :: field_nc_type = merge(nf90_float, nf90_double, field_r == real32)
  character(len=*), parameter :: magic = 'DALES_STAT_RESTART'

  character(len=*), parameter :: modname = 'modrestart_files'

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
    character :: char_identifier = ' '
    integer :: unit = -1
    logical :: open = .false.
    procedure(restart_writer), pointer, nopass :: writer => null()
    procedure(restart_reader), pointer, nopass :: reader => null()
  end type restart_entry

  type(restart_entry), save :: entries(max_entries)
  integer, save :: nentries = 0

  public :: register_restart_handlers
  public :: open_restart_file
  public :: close_restart_file
  public :: run_restart_writers
  public :: run_restart_writer
  public :: run_restart_reader
  public :: run_restart_readers

contains

  subroutine register_restart_handlers(name, character_identifier, reader, writer)
    character(len=*), intent(in) :: name
    character, intent(in) :: character_identifier
    procedure(restart_reader), optional :: reader
    procedure(restart_writer), optional :: writer

    integer :: idx

    idx = find_entry(name)
    if (idx == 0) then
      if (nentries >= max_entries) return
      nentries = nentries + 1
      idx = nentries
      entries(idx)%name = name
      entries(idx)%char_identifier = character_identifier
    end if

    if (present(reader)) entries(idx)%reader => reader
    if (present(writer)) entries(idx)%writer => writer
  end subroutine register_restart_handlers

  subroutine run_restart_writer(entry)
    type(restart_entry), intent(in) :: entry

    integer :: idx

    if (.not. associated(entry%writer)) then
      call finish(modname//'/run_restart_writer', 'restart writer not registered: ', entry%name)
    end if

    call entry%writer(entry%unit)
  end subroutine run_restart_writer

  subroutine run_restart_reader(entry)
    type(restart_entry), intent(in) :: entry

    integer :: idx

    if (.not. associated(entry%reader)) then
      call finish(modname//'/run_restart_reader', 'restart reader not registered: ', entry%name)
    end if

    call entry%reader(entry%unit)
  end subroutine run_restart_reader

  subroutine run_restart_writers(restart_name)
    use modprecision, only: longint
    character(len=*), intent(in) :: restart_name
    character(len=64) :: name

    integer :: i
    integer(kind=longint) :: fmt_version_l
    fmt_version_l = int(fmt_version, kind=longint)

    name = restart_name

    do i = 1, nentries
        name(5:5) = entries(i)%char_identifier

        call open_restart_file(name, entries(i)%unit, status='replace', action='write')

        ! call write_restart_field(entries(i)%unit, 'fmt_version', fmt_version_l)
        call run_restart_writer(entries(i))

        call close_restart_file(entries(i)%unit)

        call create_latest_symlink(name)
    end do
  end subroutine run_restart_writers
  subroutine run_restart_readers(restart_name)
    character(len=*), intent(in) :: restart_name
    character(len=64) :: name

    integer :: i

    do i = 1, nentries
        name = restart_name
        name(5:5) = entries(i)%char_identifier

        call open_restart_file(name, entries(i)%unit, action='read')

        call run_restart_reader(entries(i))
        call close_restart_file(entries(i)%unit)
    end do
  end subroutine run_restart_readers

  subroutine open_restart_file(name, iunit, status, action)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: iunit
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: action

    integer :: idx
    character(len=16) :: file_status, file_action
    character(len=nc_name_len) :: storage_name

    file_status = 'unknown'
    if (present(status)) file_status = status
    file_action = 'readwrite'
    if (present(action)) file_action = action

    storage_name = restart_storage_name(name)
    if (lrestart_netcdf) then
      call open_restart_file_netcdf(storage_name, iunit, file_status, file_action)
    else
      open(newunit=iunit, file=trim(output_prefix)//storage_name, form='unformatted', status=file_status, action=file_action)
    end if
  end subroutine open_restart_file

  subroutine close_restart_file(iunit)
    integer, intent(inout) :: iunit

    if (lrestart_netcdf) then
      call handle_netcdf_status(modname//'/close_restart_file', nf90_close(iunit), 'restart file')
    else
      close(iunit)
    end if

    iunit = -1
  end subroutine close_restart_file

  subroutine create_latest_symlink(file_name)
    character(len=*), intent(in) :: file_name

    character(len=64) :: latest_name
    character(len=nc_name_len) :: source_name, target_name

    latest_name = file_name
    latest_name(6:13) = '_latest_'

    source_name = restart_storage_name(file_name)
    target_name = restart_storage_name(latest_name)
    call system('ln -s -f '//trim(source_name)//' '//trim(output_prefix)//trim(target_name))
  end subroutine create_latest_symlink

  subroutine open_restart_file_netcdf(file_name, iunit, file_status, file_action)
    use netcdf, only: nf90_noerr
    character(len=*), intent(in) :: file_name
    integer, intent(inout) :: iunit
    character(len=*), intent(in) :: file_status
    character(len=*), intent(in) :: file_action

    integer :: status, mode

    if (file_status == 'replace') then
      status = nf90_create(trim(output_prefix)//trim(file_name), ior(nf90_clobber, nf90_netcdf4), iunit)
      call handle_netcdf_status(modname//'/open_restart_file', status, file_name)
      return
    end if

    mode = nf90_write
    if (file_action == 'read') mode = nf90_nowrite

    status = nf90_open(trim(output_prefix)//trim(file_name), mode, iunit)
    if (status == nf90_noerr .and. file_status /= 'old' .and. file_action /= 'read') then
      status = nf90_create(trim(output_prefix)//trim(file_name), ior(nf90_clobber, nf90_netcdf4), iunit)
    end if
    call handle_netcdf_status(modname//'/open_restart_file', status, file_name)
  end subroutine open_restart_file_netcdf

  function restart_storage_name(name) result(storage_name)
    character(len=*), intent(in) :: name
    character(len=nc_name_len) :: storage_name
    integer :: nlen

    storage_name = name
    if (lrestart_netcdf) then
      nlen = len_trim(name)
      if (nlen >= 3 .and. name(nlen-2:nlen) == '.nc') then
        storage_name = name
      else
        storage_name = trim(name)//'.nc'
      end if
    end if
  end function restart_storage_name

  integer function find_entry(name)
    character(len=*), intent(in) :: name
    integer :: i

    find_entry = 0
    do i = 1, nentries
      if (entries(i)%name == name) then
        find_entry = i
        return
      end if
    end do
  end function find_entry

end module modrestart_files
