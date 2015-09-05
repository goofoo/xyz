@echo off
for /r %%I in (*.c;*.h;*.cpp;*.hpp) do (
	clang-format -i "%%I"
)
exit
