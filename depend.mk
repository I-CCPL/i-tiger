./build/w90_chk_dump.o : ./build/io_global.o ./build/io_w90.o ./build/kinds.o 
./build/main.o : ./build/io_input.o ./build/env.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/io_input.o : ./build/io_global.o ./build/w90_chk_dump.o ./build/io_w90.o 
./build/io_w90.o : ./build/mp_base.o ./build/io_global.o ./build/mat3x3.o ./build/kinds.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/kinds.o : ./build/io_param.o 
./build/mat3x3.o : ./build/kinds.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
