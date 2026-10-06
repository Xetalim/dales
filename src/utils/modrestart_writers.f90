!> Global registry for restart read/write callbacks and tagged restart field I/O.
module modrestart_writers
  use iso_fortran_env, only: real32
  use netcdf
  use modglobal, only: lrestart_netcdf
  use modlogging, only: finish
  use modprecision, only: field_r

  implicit none
  private

  public :: write_restart_field, read_restart_field, handle_netcdf_status

  integer, parameter :: max_entries = 128
  integer, parameter :: fmt_version = 1
  integer, parameter :: tag_len = 64
  integer, parameter :: field_nc_type = merge(nf90_float, nf90_double, field_r == real32)
  
  character(len=*), parameter :: modname = 'modrestart_writers'

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

contains
  subroutine handle_netcdf_status(context, status, name)
    character(len=*), intent(in) :: context
    integer, intent(in) :: status
    character(len=*), intent(in) :: name

    if (status /= nf90_noerr) then
      call finish(context, 'netcdf error for ', name, ': ', nf90_strerror(status))
    end if
  end subroutine handle_netcdf_status
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
      call finish(modname//'/read_restart_field', 'restart file ended before tagged field ', tag, ' could be read')
    end if

    if (stored_tag /= tag) then
      call finish(modname//'/read_restart_field', 'expected tagged restart field ', tag, ' but found ', stored_tag, &
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/write_restart_field', nf90_put_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
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
      call handle_netcdf_status(modname//'/read_restart_field', nf90_get_var(iunit, varid, field), tag)
    else
      call read_restart_tag(iunit, tag)
      read(iunit) field
    end if
  end subroutine read_restart_field_longint_s

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
      call handle_netcdf_status(modname//'/write_restart_field', status, tag)
    end if

    rank = 0
    if (present(dims)) rank = size(dims)

    do i = 1, rank
      write(dim_name, '(A,"_",I1)') trim(tag), i
      status = nf90_def_dim(ncid, trim(dim_name), dims(i), dimids(i))
      if (status == nf90_enameinuse) then
        status = nf90_inq_dimid(ncid, trim(dim_name), dimids(i))
      end if
      call handle_netcdf_status(modname//'/write_restart_field', status, dim_name)
    end do

    if (rank == 0) then
      status = nf90_def_var(ncid, trim(tag), xtype, varid)
    else
      status = nf90_def_var(ncid, trim(tag), xtype, dimids(1:rank), varid)
    end if
    call handle_netcdf_status(modname//'/write_restart_field', status, tag)

    call handle_netcdf_status(modname//'/write_restart_field', nf90_enddef(ncid), tag)
  end subroutine define_restart_var_netcdf

  subroutine lookup_restart_var_netcdf(ncid, tag, varid)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: tag
    integer, intent(out) :: varid

    integer :: status

    status = nf90_inq_varid(ncid, trim(tag), varid)
    if (status == nf90_enotvar) then
      call finish(modname//'/read_restart_field', 'expected tagged restart field ', tag, &
                  ' but no NetCDF variable with that name was found')
    end if
    call handle_netcdf_status(modname//'/read_restart_field', status, tag)
  end subroutine lookup_restart_var_netcdf

end module modrestart_writers
