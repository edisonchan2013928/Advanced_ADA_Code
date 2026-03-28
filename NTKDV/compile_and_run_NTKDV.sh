#Compile the program
g++ -c shortest_path.cpp -w -o shortest_path.o
g++ -c KAF.cpp -w -o KAF.o
g++ -c TIC.cpp -w -o TIC.o
g++ -c TIC_phi.cpp -w -o TIC_phi.o
g++ -c TIS.cpp -w -o TIS.o
g++ -c alg_NTKDV.cpp -w -o alg_NTKDV.o
g++ main.cpp -O3 -o main shortest_path.o KAF.o TIC.o TIC_phi.o TIS.o alg_NTKDV.o
exit #Remove this line if you want to run the code

#our_model.network_fileName = argv[1];
#our_model.out_NTKDV_fileName = argv[2];
#our_model.method = atoi(argv[3]);
#our_model.lixel_reg_length = atoi(argv[4]);
#our_model.k_type = atoi(argv[5]);
#our_model.bandwidth_s = atof(argv[6]);
#our_model.bandwidth_t = atof(argv[7]);
#our_model.T = atoi(argv[8]);

dir="./NTKDV_Datasets/"
out_dir="./Results/"
dataset="New_York_taxi"
lixel_reg_length=10
k_type=2
bandwidth_s=1000
bandwidth_t=14
T=32

method=1 #RQS
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_bs"$bandwidth_s"_bt"$bandwidth_t $method $lixel_reg_length $k_type $bandwidth_s $bandwidth_t $T

method=2 #SPS
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_bs"$bandwidth_s"_bt"$bandwidth_t $method $lixel_reg_length $k_type $bandwidth_s $bandwidth_t $T

method=3 #ADA_{TIC}^{B+}
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_bs"$bandwidth_s"_bt"$bandwidth_t $method $lixel_reg_length $k_type $bandwidth_s $bandwidth_t $T

method=4 #ADA_{TIS}^{B+}
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_bs"$bandwidth_s"_bt"$bandwidth_t $method $lixel_reg_length $k_type $bandwidth_s $bandwidth_t $T