@echo off
cl /O2 /nologo /EHsc main.cpp && main
if errorlevel 1 goto end
fc /b delta_data_new.bin output.bin
:end