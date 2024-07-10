# Cria a pasta build se não existir
if [ ! -d "build" ]; then
  mkdir build
fi

cmake -S . -B ./build
(cd build && make)

# Executa o programa
./build/CircuitosClifford.out