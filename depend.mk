./build/Berry.o : ./build/io_input.o ./build/system.o ./build/constants.o ./build/kinds.o 
./build/NLO.o : ./build/delta_funct.o ./build/io_output.o ./build/io_global.o ./build/mp_base.o ./build/kpoints.o ./build/system.o ./build/constants.o ./build/io_input.o ./build/kinds.o 
./build/OAM.o : ./build/system.o ./build/io_input.o ./build/constants.o ./build/kinds.o 
./build/io_output.o : ./build/io_input.o ./build/kpoints.o ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/itg_R.o : ./build/io_input.o ./build/io_global.o ./build/fft.o ./build/constants.o ./build/R_vector.o ./build/wannier90.o ./build/system.o ./build/kinds.o 
./build/itg_k.o : ./build/io_output.o ./build/mp_base.o ./build/io_global.o ./build/NLO.o ./build/lin_eig_H.o ./build/io_input.o ./build/fft.o ./build/constants.o ./build/kpoints.o ./build/itg_R.o ./build/system.o ./build/kinds.o 
./build/itg_q.o : ./build/write_matrix.o ./build/io_global.o ./build/wannier90.o ./build/io_input.o 
./build/main.o : ./build/itg_k.o ./build/itg_R.o ./build/itg_q.o ./build/kpoints.o ./build/io_input.o ./build/io_global.o ./build/env.o 
./build/velocity.o : ./build/io_input.o ./build/fft.o ./build/kpoints.o ./build/R_vector.o ./build/system.o ./build/constants.o ./build/kinds.o 
./build/check_Hermicity.o : ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/check_nan.o : ./build/kinds.o 
./build/debug_R.o : ./build/itg_R.o ./build/io_input.o ./build/io_global.o 
./build/debug_k.o : ./build/kpoints.o ./build/itg_k.o ./build/mp_global.o ./build/io_input.o ./build/io_global.o 
./build/debug_q.o : ./build/system.o ./build/itg_q.o ./build/wannier90.o ./build/io_input.o ./build/io_global.o ./build/kinds.o 
./build/dump_vec_io.o : ./build/io_global.o ./build/kinds.o 
./build/module_check.o : ./build/kinds.o 
./build/w90_chk_dum.o : ./build/system.o ./build/dump_vec_io.o 
./build/warning.o : ./build/io_global.o 
./build/write_matrix.o : ./build/kpoints.o ./build/system.o ./build/kinds.o 
./build/R_vector.o : ./build/wannier90.o ./build/constants.o ./build/algo_unique.o ./build/system.o ./build/io_global.o ./build/kinds.o 
./build/delta_funct.o : ./build/constants.o ./build/kinds.o 
./build/fft.o : ./build/kpoints.o ./build/R_vector.o ./build/wannier90.o ./build/system.o ./build/io_input.o ./build/io_global.o ./build/constants.o ./build/kinds.o 
./build/grid.o : ./build/kinds.o 
./build/io_input.o : ./build/system.o ./build/kpoints.o ./build/wannier90.o ./build/char.o ./build/mp_base.o ./build/io_global.o ./build/kinds.o 
./build/k_mp.o : ./build/kpoints.o ./build/kinds.o ./build/mp_global.o 
./build/kpoints.o : ./build/constants.o ./build/system.o ./build/kinds.o 
./build/system.o : ./build/lin_vec.o ./build/lin_mat.o ./build/kinds.o 
./build/w90_base.o : ./build/constants.o ./build/io_input.o 
./build/w90_chk.o : ./build/io_output.o ./build/lin_eig_H.o ./build/wannier90.o ./build/system.o ./build/mp_base.o ./build/io_global.o ./build/constants.o ./build/kinds.o 
./build/w90_eig.o : ./build/mp_base.o ./build/io_global.o 
./build/w90_mmn.o : ./build/constants.o ./build/system.o ./build/algo_unique.o ./build/kinds.o ./build/io_global.o ./build/mp_base.o 
./build/wannier90.o : ./build/kpoints.o ./build/io_global.o ./build/kinds.o 
./build/algo_unique.o : ./build/kinds.o 
./build/clock.o : ./build/io_global.o ./build/mp_global.o ./build/kinds.o 
./build/constants.o : ./build/kinds.o 
./build/env.o : ./build/char.o ./build/kinds.o ./build/mp_global.o ./build/io_global.o 
./build/errore.o : ./build/io_global.o ./build/mp_global.o 
./build/io_global.o : ./build/mp_base.o ./build/io_param.o 
./build/kinds.o : ./build/io_param.o 
./build/lin_eig_H.o : ./build/constants.o ./build/mp_base.o ./build/io_global.o ./build/kinds.o 
./build/lin_mat.o : ./build/kinds.o 
./build/lin_vec.o : ./build/kinds.o 
./build/mp_base.o : ./build/mp_global.o ./build/kinds.o 
./build/mp_global.o : ./build/io_param.o 
./build/mp_ops.o : ./build/kinds.o ./build/mp_global.o 
