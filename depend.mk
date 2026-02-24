./build/main.o : ./build/env.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/mp_base.o : ./build/kinds.o 
./build/mp_global.o : ./build/io_global.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/kinds.o : 
