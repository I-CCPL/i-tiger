./build/main.o : ./build/io_input.o ./build/env.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/io_input.o : ./build/io_global.o ./build/io_w90.o 
./build/io_w90.o : ./build/io_global.o ./build/kinds.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/kinds.o : ./build/io_param.o 
