module modinputchecking

use modprecision, only: field_r
use ieee_exceptions

interface check_grid_variable
    module procedure check_grid_variable_1d_int
    module procedure check_grid_variable_1d_real
    module procedure check_grid_variable_2d_real
    module procedure check_grid_variable_2d_int
    module procedure check_grid_variable_3d_real
    module procedure check_grid_variable_3d_int
end interface check_grid_variable

private :: check_grid_variable_1d_int
private :: check_grid_variable_1d_real
private :: check_grid_variable_2d_real
private :: check_grid_variable_2d_int
private :: check_grid_variable_3d_real
private :: check_grid_variable_3d_int


contains
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_1d_int( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    INTEGER, INTENT(IN) ::  valid_max  !< upper bound of valid range for var
    INTEGER, INTENT(IN) ::  valid_min  !< lower bound of valid range for var

    INTEGER, DIMENSION(:), INTENT(IN) ::  var  !< target variable


    !
    !-- Check if the input is within allowed bounds.
    DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
        IF ( ( var(i) < valid_min  .OR.  var(i) > valid_max ) )  THEN
            WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                       &
                                        ' for index (i) = ', i, ' set to',   &
                                        var(i), ' valid range is [', valid_min, valid_max, '].'
            call abort
        ENDIF
    ENDDO

END SUBROUTINE check_grid_variable_1d_int


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_1d_real( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(field_r), INTENT(IN) ::  valid_max  !< upper bound of valid range for var
    REAL(field_r), INTENT(IN) ::  valid_min  !< lower bound of valid range for var

    REAL(field_r), DIMENSION(:), INTENT(IN) ::  var  !< target variable


    !
    !-- Check if the input is within allowed bounds.
    DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
        IF ( ( var(i) < valid_min  .OR.  var(i) > valid_max ) )  THEN
            WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                       &
                                        ' for index (i) = ', i, ' set to',   &
                                        var(i), ' valid range is [', valid_min, valid_max, '].'
            call abort
        ENDIF
    ENDDO

END SUBROUTINE check_grid_variable_1d_real


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_2d_real( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(field_r), INTENT(IN) ::  valid_max   !< upper bound of valid range for var
    REAL(field_r), INTENT(IN) ::  valid_min   !< lower bound of valid range for var

    REAL(field_r), DIMENSION(:,:), INTENT(IN) ::  var   !< target variable


    !
    !-- Check if the input is within allowed bounds
    DO j = LBOUND( var, 2 ), UBOUND( var, 2 )
        DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
            IF ( ( var(i,j) < valid_min  .OR.  var(i,j) > valid_max ) )  THEN
                WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                    &
                                        ' for grid cell (i,j) = ', i,j,      &
                                        ' set to', var(i,j), ' valid range is [',                  &
                                        valid_min, valid_max, '].'
                call abort
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE check_grid_variable_2d_real

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_2d_int( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    integer, INTENT(IN) ::  valid_max   !< upper bound of valid range for var
    integer, INTENT(IN) ::  valid_min   !< lower bound of valid range for var

    integer, DIMENSION(:,:), INTENT(IN) ::  var   !< target variable


    !
    !-- Check if the input is within allowed bounds
    DO j = LBOUND( var, 2 ), UBOUND( var, 2 )
        DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
            IF ( ( var(i,j) < valid_min  .OR.  var(i,j) > valid_max ) )  THEN
                WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                    &
                                        ' for grid cell (i,j) = ', i,j,      &
                                        ' set to', var(i,j), ' valid range is [',                  &
                                        valid_min, valid_max, '].'
                    call abort
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE check_grid_variable_2d_int
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_3d_real( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    REAL(field_r), INTENT(IN) ::  valid_max   !< upper bound of valid range for var
    REAL(field_r), INTENT(IN) ::  valid_min   !< lower bound of valid range for var

    REAL(field_r), DIMENSION(:,:,:), INTENT(IN) ::  var   !< target variable


    !
    !-- Check if the input is within allowed bounds
    DO k = LBOUND( var, 3), UBOUND( var, 3)
        DO j = LBOUND( var, 2 ), UBOUND( var, 2 )
            DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
                IF ( ( var(i,j,k) < valid_min  .OR.  var(i,j,k) > valid_max ) )  THEN
                    WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                    &
                                            ' for grid cell (i,j,k) = ', i,j,k,      &
                                            ' set to', var(i,j,k), ' valid range is [',                  &
                                            valid_min, valid_max, '].'
                    call abort
                ENDIF
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE check_grid_variable_3d_real

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the given variable is in the predefined valid range.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_grid_variable_3d_int( varname, var, valid_min, valid_max )

    CHARACTER(LEN=*), INTENT(IN) ::  varname  !< variable name in file

    integer, INTENT(IN) ::  valid_max   !< upper bound of valid range for var
    integer, INTENT(IN) ::  valid_min   !< lower bound of valid range for var

    integer, DIMENSION(:,:,:), INTENT(IN) ::  var   !< target variable


    !
    !-- Check if the input is within allowed bounds
    DO k = LBOUND( var, 3), UBOUND( var, 3)
        DO j = LBOUND( var, 2 ), UBOUND( var, 2 )
            DO  i = LBOUND( var, 1 ), UBOUND( var, 1 )
                IF ( ( var(i,j,k) < valid_min  .OR.  var(i,j,k) > valid_max ) )  THEN
                    WRITE( 6, * ) 'Input variable ' // TRIM( varname ) //                    &
                                            ' for grid cell (i,j,k) = ', i,j,k,      &
                                            ' set to', var(i,j,k), ' valid range is [',                  &
                                            valid_min, valid_max, '].'
                    call abort
                ENDIF
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE check_grid_variable_3d_int
end module modinputchecking