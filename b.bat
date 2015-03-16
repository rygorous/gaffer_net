@echo off
cl /W4 /O2 /nologo /EHsc main.cpp && main
if errorlevel 1 goto end
fc /b delta_data_realnew.bin output.bin
:end