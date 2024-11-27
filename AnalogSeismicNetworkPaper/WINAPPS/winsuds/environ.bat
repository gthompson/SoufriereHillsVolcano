@ECHO OFF

ECHO Setting up Win-SUDS environment...

REM - Edit these next two paths as appropriate for your system:
SET SUDS_HOME=c:\suds
SET PC_SUDS=c:\pcsuds

REM - SCSITape host and drive identifiers...
SET SCSI_HOST=0
SET SCSI_TAPE=5

SET SUDS_ETC=%SUDS_HOME%\etc
SET SUDS_BIN=%SUDS_HOME%\bin

PATH=%PATH%;%SUDS_BIN%;%PC_SUDS%

