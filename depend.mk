./build/OAM.o : ./build/kpoints.o ./build/system.o ./build/io_input.o ./build/constants.o ./build/kinds.o 
./build/itg_R.o : ./build/write_matrix.o ./build/fft.o ./build/constants.o ./build/R_vector.o ./build/wannier90.o ./build/system.o ./build/kinds.o 
./build/itg_k.o : ./build/io_input.o ./build/write_matrix.o ./build/lin_eig_H.o ./build/fft.o ./build/kpoints.o ./build/itg_R.o ./build/system.o ./build/kinds.o 
./build/itg_q.o : ./build/write_matrix.o ./build/wannier90.o 
./build/main.o : ./build/itg_k.o ./build/itg_R.o ./build/itg_q.o ./build/io_input.o ./build/env.o 
./build/velocity.o : ./build/derivation.o ./build/fft.o ./build/kpoints.o ./build/R_vector.o ./build/write_matrix.o ./build/system.o ./build/constants.o ./build/kinds.o 
./build/check_Hermicity.o : ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/dump_vec_io.o : ./build/io_global.o ./build/kinds.o 
./build/w90_chk_dum.o : ./build/system.o ./build/dump_vec_io.o 
./build/warning.o : ./build/io_global.o 
./build/write_matrix.o : ./build/kpoints.o ./build/system.o ./build/kinds.o 
./build/R_vector.o : ./build/wannier90.o ./build/constants.o ./build/algo_unique.o ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/derivation.o : ./build/constants.o ./build/R_vector.o ./build/wannier90.o ./build/system.o ./build/kinds.o 
./build/fft.o : ./build/kpoints.o ./build/R_vector.o ./build/wannier90.o ./build/system.o ./build/io_global.o ./build/constants.o ./build/kinds.o 
./build/grid.o : ./build/kinds.o 
./build/io_input.o : ./build/kpoints.o ./build/wannier90.o ./build/char.o ./build/mp_base.o ./build/io_global.o ./build/kinds.o 
./build/kpoints.o : ./build/constants.o ./build/system.o ./build/kinds.o 
./build/system.o : ./build/lin_mat.o ./build/kinds.o 
./build/w90_base.o : ./build/io_input.o 
./build/w90_chk.o : ./build/lin_eig_H.o ./build/wannier90.o ./build/system.o ./build/mp_base.o ./build/io_global.o ./build/kinds.o 
./build/w90_eig.o : ./build/mp_base.o ./build/io_global.o 
./build/w90_mmn.o : ./build/constants.o ./build/system.o ./build/algo_unique.o ./build/kinds.o ./build/io_global.o ./build/mp_base.o 
./build/wannier90.o : ./build/kpoints.o ./build/io_global.o ./build/kinds.o 
./build/algo_unique.o : ./build/kinds.o 
./build/constants.o : ./build/kinds.o 
./build/env.o : ./build/mp_global.o ./build/io_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/kinds.o : ./build/io_param.o 
./build/lin_eig_H.o : ./build/mp_base.o ./build/kpoints.o ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/lin_mat.o : ./build/kinds.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
