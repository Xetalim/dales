!> Global registry for restart read/write callbacks and tagged restart field I/O.
module modrestart_registry
  use iso_fortran_env, only: real32
  use netcdf
  use modglobal, only: output_prefix, nsv, lrestart_netcdf
  use modlogging, only: finish
  use modprecision, only: field_r
  use modsurfdata, only: isurf

  implicit none
  private

  integer, parameter :: max_entries = 128
  integer, parameter :: fmt_version = 1
  integer, parameter :: tag_len = 64
  integer, parameter :: nc_name_len = 128
  integer, parameter :: field_nc_type = merge(nf90_float, nf90_double, field_r == real32)
  character(len=*), parameter :: magic = 'DALES_STAT_RESTART'

  abstract interface
    subroutine restart_writer(iunit)
      integer, intent(in) :: iunit
    end subroutine restart_writer

    subroutine restart_reader(iunit)
      integer, intent(in) :: iunit
    end subroutine restart_reader
  end interface

  interface write_restart_field
    procedure :: write_restart_field_real_s
    procedure :: write_restart_field_real_r1
    procedure :: write_restart_field_real_r2
    procedure :: write_restart_field_real_r3
    procedure :: write_restart_field_real_r4
    procedure :: write_restart_field_longint_s
  end interface

  interface read_restart_field
    procedure :: read_restart_field_real_s
    procedure :: read_restart_field_real_r1
    procedure :: read_restart_field_real_r2
    procedure :: read_restart_field_real_r3
    procedure :: read_restart_field_real_r4
    procedure :: read_restart_field_longint_s
  end interface

  type restart_entry
    character(len=64) :: name = ''
    procedure(restart_writer), pointer, nopass :: writer => null()
    procedure(restart_reader), pointer, nopass :: reader => null()
    logical :: include_in_cleanup = .true.
  end type restart_entry

  type(restart_entry), save :: entries(max_entries)
  integer, save :: nentries = 0

  type restart_file_entry
    character(len=64) :: name = ''
    integer :: unit = -1
    logical :: open = .false.
  end type restart_file_entry

  type(restart_file_entry), save :: files(max_entries)
  integer, save :: nfiles = 0

  public :: register_restart_handlers
  public :: open_restart_file
  public :: close_restart_file
  public :: write_restart_field
  public :: read_restart_field
  public :: run_cleanups
  public :: run_restart_writers
  public :: run_restart_writer
  public :: run_restart_reader

contains

  subroutine register_restart_handlers(name, reader, writer, include_in_cleanup)
    character(len=*), intent(in) :: name
    procedure(restart_reader), optional :: reader
    procedure(restart_writer), optional :: writer
    logical, intent(in), optional :: include_in_cleanup

    integer :: idx

    idx = find_entry(name)
    if (idx == 0) then
      if (nentries >= max_entries) return
      nentries = nentries + 1
      idx = nentries
      entries(idx)%name = name
    end if

    if (present(reader)) entries(idx)%reader => reader
    if (present(writer)) entries(idx)%writer => writer
    if (present(include_in_cleanup)) entries(idx)%include_in_cleanup = include_in_cleanup
  end subroutine register_restart_handlers

  subroutine run_restart_writer(name, iunit)
    character(len=*), intent(in) :: name
    integer, intent(in) :: iunit

    integer :: idx

    idx = find_entry(name)
    if (idx == 0) then
      call finish('modrestart_registry/run_restart_writer', 'unknown restart handler: ', name)
    end if
    if (.not. associated(entries(idx)%writer)) then
      call finish('modrestart_registry/run_restart_writer', 'restart writer not registered: ', name)
    end if

    call entries(idx)%writer(iunit)
  end subroutine run_restart_writer

  subroutine run_restart_reader(name, iunit)
    character(len=*), intent(in) :: name
    integer, intent(in) :: iunit

    integer :: idx

    idx = find_entry(name)
    if (idx == 0) then
      call finish('modrestart_registry/run_restart_reader', 'unknown restart handler: ', name)
    end if
    if (.not. associated(entries(idx)%reader)) then
      call finish('modrestart_registry/run_restart_reader', 'restart reader not registered: ', name)
    end if

    call entries(idx)%reader(iunit)
  end subroutine run_restart_reader

  subroutine run_restart_writers(restart_name)
    character(len=*), intent(in) :: restart_name

    integer :: i, iunit
    character(len=64) :: file_name
    logical :: should_write
    logical :: replace_file

    do i = 1, nentries
      if (entries(i)%include_in_cleanup) cycle
      if (.not. associated(entries(i)%writer)) cycle

      call resolve_restart_file(restart_name, entries(i)%name, file_name, should_write, replace_file)
      if (.not. should_write) cycle

      if (replace_file) then
        call open_restart_file(file_name, iunit, status='replace')
      else
        call open_restart_file(file_name, iunit)
      end if

      call entries(i)%writer(iunit)
      call close_restart_file(file_name)
      call run_cleanups(file_name)
    end do
  end subroutine run_restart_writers

  subroutine open_restart_file(name, iunit, status, action)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: iunit
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: action

    integer :: idx
    character(len=16) :: file_status, file_action
    character(len=nc_name_len) :: storage_name

    idx = find_file(name)
    if (idx /= 0 .and. files(idx)%open) then
      call finish('modrestart_registry/open_restart_file', 'restart file already open: ', name)
    end if

    if (idx == 0) then
      if (nfiles >= max_entries) then
        call finish('modrestart_registry/open_restart_file', 'too many open restart files')
      end if
      nfiles = nfiles + 1
      idx = nfiles
      files(idx)%name = name
    end if

    file_status = 'unknown'
    if (present(status)) file_status = status
    file_action = 'readwrite'
    if (present(action)) file_action = action

    storage_name = restart_storage_name(name)
    if (lrestart_netcdf) then
      call open_restart_file_netcdf(storage_name, iunit, file_status, file_action)
    else
      open(unit=iunit, file=trim(output_prefix)//storage_name, form='unformatted', status=file_status, action=file_action)
    end if

    files(idx)%unit = iunit
    files(idx)%open = .true.
  end subroutine open_restart_file

  subroutine close_restart_file(name)
    character(len=*), intent(in) :: name

    integer :: idx

    idx = find_file(name)
    if (idx == 0) then
      call finish('modrestart_registry/close_restart_file', 'unknown restart file handle: ', name)
    end if
    if (.not. files(idx)%open) then
      call finish('modrestart_registry/close_restart_file', 'restart file is not open: ', name)
    end if

    if (lrestart_netcdf) then
      call handle_netcdf_status('modrestart_registry/close_restart_file', nf90_close(files(idx)%unit), name)
    else
      close(files(idx)%unit)
    end if

    files(idx)%unit = -1
    files(idx)%open = .false.
  end subroutine close_restart_file

  subroutine write_restart_tag(iunit, tag)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag

    character(len=tag_len) :: stored_tag

    stored_tag = tag
    write(iunit) stored_tag
  end subroutine write_restart_tag

  subroutine read_restart_tag(iunit, tag)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag

    character(len=tag_len) :: stored_tag
    integer :: ios

    stored_tag = ''
    read(iunit, iostat=ios) stored_tag
    if (ios /= 0) then
      call finish('modrestart_registry/read_restart_field', 'restart file ended before tagged field ', tag, ' could be read')
    end if

    if (stored_tag /= tag) then
      call finish('modrestart_registry/read_restart_field', 'expected tagged restart field ', tag, ' but found ', stored_tag, &
                  '. The restart file is missing field tags or the field order is different.')
    end if
  end subroutine read_restart_tag

  subroutine write_restart_field_real_s(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(in) :: field
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, field_nc_type, varid)
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_real_s

  subroutine write_restart_field_real_r1(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(in) :: field(:)
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, field_nc_type, varid, [size(field,1)])
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_real_r1

  subroutine write_restart_field_real_r2(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(in) :: field(:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, field_nc_type, varid, [size(field,1), size(field,2)])
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_real_r2

  subroutine write_restart_field_real_r3(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(in) :: field(:,:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, field_nc_type, varid, [size(field,1), size(field,2), size(field,3)])
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_real_r3

  subroutine write_restart_field_real_r4(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(in) :: field(:,:,:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, field_nc_type, varid, [size(field,1), size(field,2), size(field,3), size(field,4)])
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_real_r4

  subroutine write_restart_field_longint_s(iunit, tag, field)
    use modprecision, only : longint

    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    integer(kind=longint), intent(in) :: field
    integer :: varid

    if (lrestart_netcdf) then
      call define_restart_var_netcdf(iunit, tag, nf90_int64, varid)
      call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_put_var(iunit, varid, field), tag)
    else
      call write_restart_tag(iunit, tag)
      write(iunit) field
    end if
  end subroutine write_restart_field_longint_s

  subroutine read_restart_field_real_s(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(inout) :: field
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_real_s

  subroutine read_restart_field_real_r1(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(inout) :: field(:)
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_real_r1

  subroutine read_restart_field_real_r2(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(inout) :: field(:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_real_r2

  subroutine read_restart_field_real_r3(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(inout) :: field(:,:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_real_r3

  subroutine read_restart_field_real_r4(iunit, tag, field)
    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    real(kind=field_r), intent(inout) :: field(:,:,:,:)
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_real_r4

  subroutine read_restart_field_longint_s(iunit, tag, field)
    use modprecision, only : longint

    integer, intent(in) :: iunit
    character(len=*), intent(in) :: tag
    integer(kind=longint), intent(inout) :: field
    integer :: varid

    if (lrestart_netcdf) then
      call lookup_restart_var_netcdf(iunit, tag, varid)
      call handle_netcdf_status('modrestart_registry/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_longint_s

  subroutine run_cleanups(file_name)
    character(len=*), intent(in) :: file_name
    use modprecision, only : longint

    integer :: i, iunit, nwriters
    character(len=64) :: cleanup_name
    integer(kind=longint) :: fmt_version_l, nwriters_l

    call create_latest_symlink(file_name)

    if (file_name(5:5) /= 'd') return

    cleanup_name = file_name
    cleanup_name(5:5) = 'r'

    if (lrestart_netcdf) then
      nwriters = 0
      do i = 1, nentries
        if (entries(i)%include_in_cleanup .and. associated(entries(i)%writer)) nwriters = nwriters + 1
      end do

      call open_restart_file(cleanup_name, iunit, status='replace', action='write')

      fmt_version_l = int(fmt_version, kind=longint)
      nwriters_l = int(nwriters, kind=longint)
      call write_restart_field(iunit, 'fmt_version', fmt_version_l)
      call write_restart_field(iunit, 'num_cleanup_writers', nwriters_l)

      do i = nentries, 1, -1
        if (.not. entries(i)%include_in_cleanup) cycle
        if (.not. associated(entries(i)%writer)) cycle
        call entries(i)%writer(iunit)
      end do

      call close_restart_file(cleanup_name)
      call create_latest_symlink(cleanup_name)
      return
    end if

    call open_restart_file(cleanup_name, iunit, status='replace', action='write')

    nwriters = 0
    do i = 1, nentries
      if (entries(i)%include_in_cleanup .and. associated(entries(i)%writer)) nwriters = nwriters + 1
    end do

    write(iunit) magic, fmt_version, nwriters

    do i = nentries, 1, -1
      if (.not. entries(i)%include_in_cleanup) cycle
      if (.not. associated(entries(i)%writer)) cycle
      write(iunit) entries(i)%name
      call entries(i)%writer(iunit)
    end do

    call close_restart_file(cleanup_name)
    call create_latest_symlink(cleanup_name)
  end subroutine run_cleanups

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
    character(len=*), intent(in) :: file_name
    integer, intent(inout) :: iunit
    character(len=*), intent(in) :: file_status
    character(len=*), intent(in) :: file_action

    integer :: status, mode

    if (file_status == 'replace') then
      status = nf90_create(trim(output_prefix)//trim(file_name), nf90_clobber, iunit)
      call handle_netcdf_status('modrestart_registry/open_restart_file', status, file_name)
      return
    end if

    mode = nf90_write
    if (file_action == 'read') mode = nf90_nowrite

    status = nf90_open(trim(output_prefix)//trim(file_name), mode, iunit)
    if (status == nf90_enoent .and. file_status /= 'old' .and. file_action /= 'read') then
      status = nf90_create(trim(output_prefix)//trim(file_name), nf90_clobber, iunit)
    end if
    call handle_netcdf_status('modrestart_registry/open_restart_file', status, file_name)
  end subroutine open_restart_file_netcdf

  subroutine define_restart_var_netcdf(ncid, tag, xtype, varid, dims)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: tag
    integer, intent(in) :: xtype
    integer, intent(out) :: varid
    integer, intent(in), optional :: dims(:)

    integer :: status, i, rank
    integer :: dimids(4)
    character(len=64) :: dim_name

    status = nf90_redef(ncid)
    if (status /= nf90_noerr .and. status /= nf90_eindefine) then
      call handle_netcdf_status('modrestart_registry/write_restart_field', status, tag)
    end if

    rank = 0
    if (present(dims)) rank = size(dims)

    do i = 1, rank
      write(dim_name, '(A,"_",I1)') trim(tag), i
      status = nf90_def_dim(ncid, trim(dim_name), dims(i), dimids(i))
      if (status == nf90_enameinuse) then
        status = nf90_inq_dimid(ncid, trim(dim_name), dimids(i))
      end if
      call handle_netcdf_status('modrestart_registry/write_restart_field', status, dim_name)
    end do

    if (rank == 0) then
      status = nf90_def_var(ncid, trim(tag), xtype, varid)
    else
      status = nf90_def_var(ncid, trim(tag), xtype, dimids(1:rank), varid)
    end if
    call handle_netcdf_status('modrestart_registry/write_restart_field', status, tag)

    call handle_netcdf_status('modrestart_registry/write_restart_field', nf90_enddef(ncid), tag)
  end subroutine define_restart_var_netcdf

  subroutine lookup_restart_var_netcdf(ncid, tag, varid)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: tag
    integer, intent(out) :: varid

    integer :: status

    status = nf90_inq_varid(ncid, trim(tag), varid)
    if (status == nf90_enotvar) then
      call finish('modrestart_registry/read_restart_field', 'expected tagged restart field ', tag, &
                  ' but no NetCDF variable with that name was found')
    end if
    call handle_netcdf_status('modrestart_registry/read_restart_field', status, tag)
  end subroutine lookup_restart_var_netcdf

  subroutine handle_netcdf_status(context, status, name)
    character(len=*), intent(in) :: context
    integer, intent(in) :: status
    character(len=*), intent(in) :: name

    if (status /= nf90_noerr) then
      call finish(context, 'netcdf error for ', name, ': ', nf90_strerror(status))
    end if
  end subroutine handle_netcdf_status

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

  subroutine resolve_restart_file(restart_name, handler_name, file_name, should_write, replace_file)
    character(len=*), intent(in) :: restart_name
    character(len=*), intent(in) :: handler_name
    character(len=64), intent(out) :: file_name
    logical, intent(out) :: should_write
    logical, intent(out) :: replace_file

    file_name = restart_name
    should_write = .true.
    replace_file = .false.

    select case (handler_name)
    case ('state')
      replace_file = .true.
    case ('scalar')
      file_name(5:5) = 's'
      should_write = nsv > 0
    case ('surface')
      file_name(5:5) = 'l'
      should_write = isurf == 1 .or. isurf == 11
    case default
      call finish('modrestart_registry/run_restart_writers', 'unknown main restart handler: ', handler_name)
    end select
  end subroutine resolve_restart_file

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

  integer function find_file(name)
    character(len=*), intent(in) :: name
    integer :: i

    find_file = 0
    do i = 1, nfiles
      if (files(i)%name == name) then
        find_file = i
        return
      end if
    end do
  end function find_file

end module modrestart_registry
