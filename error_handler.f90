MODULE Error_Logging
    !> Module designed to handle error and warning logging
    !! This module provides various subroutines to manage, update and print 
    !! lists of warnings and errors ocurring within code

    !> @author G. A. O. Seaford
    !! @version 1.1
    !! @date 11th November 2024

    IMPLICIT NONE

    !> New type for messages
    TYPE :: Message
        CHARACTER(LEN=:), ALLOCATABLE :: text
    END TYPE Message

    TYPE(Message), ALLOCATABLE, DIMENSION(:) :: error_log
    TYPE(Message), ALLOCATABLE, DIMENSION(:) :: warning_log

    INTEGER, PARAMETER :: error_unit = 0 
    INTEGER :: num_errors = 0
    INTEGER :: num_warnings = 0

    !> Flags to control printing behavior
    LOGICAL :: with_warnings = .FALSE.
    LOGICAL :: with_errors = .TRUE.

    CONTAINS

    SUBROUTINE Initialise_Log()
        !> Initializes the error and warning logs by allocating initial space.

        IF (.NOT. ALLOCATED(error_log)) THEN
            ALLOCATE(error_log(20))
            IF (.NOT. ALLOCATED(error_log)) THEN
                WRITE (error_unit, *) "Error allocating warning log."
                STOP 1
            END IF
            num_errors = 0
        END IF

        IF (.NOT. ALLOCATED(warning_log)) THEN
            ALLOCATE(warning_log(20))
            IF (.NOT. ALLOCATED(warning_log)) THEN
                WRITE (error_unit, *) "Error allocating message log."
                STOP 1
            END IF
            num_warnings = 0
        END IF
    END SUBROUTINE Initialise_Log

    SUBROUTINE Add_Message(log, msg, count)
        !> Subroutine to add a message to a specified log
        !! Allocates memory for log if not allocated, and 
        !! contains process to allocate additional memeory where necessary
        !! @param[inout] log    Array for message storage
        !! @param[in]    msg    Message to add to storage array
        !> @param[inout] count  Number of messages in log.

        TYPE(Message), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: log
        INTEGER :: err_code, count
        CHARACTER(LEN=*), INTENT(IN) :: msg
        
        IF (.NOT. ALLOCATED(log)) THEN
            !> Preallocate with 20 messages
            ALLOCATE(log(20), STAT=err_code) 
            IF (err_code /= 0) THEN
                WRITE (error_unit, *) "Error allocating message log."
                RETURN
            END IF
            count = 0
        END IF
        
        count = count + 1

        IF (count > SIZE(log)) THEN
            ! Reallocate with more space
            ALLOCATE(log(SIZE(log) + 10), STAT=err_code)
            IF (err_code /= 0) THEN
                WRITE (error_unit, *) "Error reallocating message log."
                RETURN
            END IF
        END IF
        
        log(count)%text = msg

    END SUBROUTINE Add_Message

    SUBROUTINE Add_Error_Message(err_msg)
        !> Subroutine to add an error message to error log
        !! @param[in] err_msg      Error message to be added to the error log

        CHARACTER(LEN=*), INTENT(IN) :: err_msg
        CALL Add_Message(error_log, err_msg, num_errors)

    END SUBROUTINE Add_Error_Message

    SUBROUTINE Add_Warning_Message(war_msg)
        !> Subroutine to add a warning message to warning log
        !! @param[in] war_msg      Warning message to be added to the error log

        CHARACTER(LEN=*), INTENT(IN) :: war_msg
        CALL Add_Message(warning_log, war_msg, num_warnings)

    END SUBROUTINE Add_Warning_Message

    SUBROUTINE Print_Errors()
        !> Prints all logged error messages by default.
        !! Can be disabled via the --without_errors argument.
        !! If no messages are logged, a default message is printed.

        INTEGER :: err_indx

        IF (.NOT. with_errors) THEN 
            RETURN
        END IF

        IF (ALLOCATED(error_log) .AND. num_errors>0) THEN
            PRINT *, "The following errors were encountered:"
            PRINT *, ""

            DO err_indx = 1, num_errors
                IF (TRIM(error_log(err_indx)%text) /= "") THEN
                    WRITE (*, *) TRIM(error_log(err_indx)%text)
                END IF
            END DO

        ELSE
            PRINT *, "No errors were logged."
            PRINT *, ""
        END IF
    END SUBROUTINE Print_Errors

    SUBROUTINE Print_Warnings()
        !> Prints all logged warning messages if the --with_warnings argument is passed.
        !! Disabled by default as warnings are non-critical.
        !! If no messages are logged, a default message is printed.
        INTEGER :: war_indx

        IF (.NOT. with_warnings) THEN 
            RETURN
        END IF

        IF (ALLOCATED(warning_log) .AND. num_warnings>0) THEN
            PRINT *, "The following warnings were encountered:"
            PRINT *, ""

            DO war_indx = 1, num_warnings
                IF (TRIM(warning_log(war_indx)%text) /= "") THEN
                    WRITE (*, *) TRIM(warning_log(war_indx)%text)
                END IF
            END DO
            PRINT *, ""

        ELSE
            PRINT *, "No warnings were logged."
            PRINT *, ""
        END IF

    END SUBROUTINE Print_Warnings

END MODULE Error_Logging
