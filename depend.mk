./build/main.o : ./build/w90.o ./build/system.o ./build/io_input.o ./build/env.o ./build/kinds.o 
./build/cell.o : ./build/lin_mat3x3.o ./build/kinds.o 
./build/io_input.o : ./build/w90.o ./build/io_global.o 
./build/w90.o : ./build/unique.o ./build/lin_eig_H.o ./build/dump_vec_io.o ./build/cell.o ./build/system.o ./build/mp_base.o ./build/lin_mat3x3.o ./build/io_global.o ./build/kinds.o 
./build/check_Hermicity.o : ./build/kinds.o 
./build/constants.o : ./build/kinds.o 
./build/dump_vec_io.o : ./build/io_global.o ./build/kinds.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/kinds.o : ./build/io_param.o 
./build/lin_eig_H.o : ./build/mp_base.o ./build/kinds.o 
./build/lin_mat3x3.o : ./build/kinds.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
./build/unique.o : ./build/kinds.o 
