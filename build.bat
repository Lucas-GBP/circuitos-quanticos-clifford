@echo off

REM Cria a pasta build se não existir
if not exist build (
  mkdir build
)

REM Vai para a pasta build
cd build

REM Configura e compila o projeto
cmake -S . -B ./build
if %errorlevel% neq 0 exit /b %errorlevel%
cmake --build ./build --config Release

REM Executa o programa
./build/CircuitosClifford.out