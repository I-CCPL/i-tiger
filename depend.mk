./build/main.o : ./build/io_w90.o ./build/io_input.o ./build/env.o 
./build/check_Hermicity.o : ./build/kinds.o 
./build/dump_vec_io.o : ./build/io_global.o ./build/kinds.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/io_input.o : ./build/io_w90.o ./build/io_global.o 
./build/io_w90.o : ./build/lin_eig_H.o ./build/dump_vec_io.o ./build/mp_base.o ./build/io_global.o ./build/mat3x3.o ./build/kinds.o 
./build/kinds.o : ./build/io_param.o 
./build/lin_eig_H.o : ./build/mp_base.o ./build/kinds.o 
./build/mat3x3.o : ./build/kinds.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
