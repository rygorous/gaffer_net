@echo off
cl /O2 /nologo /EHsc main.cpp && main
fc /b delta_data.bin output.bin