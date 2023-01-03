@ECHO OFF
REM - Add FAST.exe directory to Windows search path ---------------------------------
ECHO * Actual Windows search path is:
ECHO * ------------------------------
ECHO * %PATH%
ECHO *
	set PATH=%PATH%D:\University_Flensburg\Third semester\WEC Project\DLC\Week 10\1.2_V10.8_04122022
ECHO * New Windows search path including FAST directory is:
ECHO * ----------------------------------------------------
ECHO * %PATH%
ECHO *
REM - Call FAST ---------------------------------------------------------------------
	openfast_x64 IEA-15-255-RWT-UMaineSemi.fst
PAUSE
EXIT /B
