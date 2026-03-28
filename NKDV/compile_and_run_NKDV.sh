#Compile the program
g++ -c shortest_path.cpp -w -o shortest_path.o
g++ -c KAF.cpp -w -o KAF.o
g++ -c alg_NKDV.cpp -w -o alg_NKDV.o
g++ main.cpp -O3 -o main shortest_path.o KAF.o alg_NKDV.o
exit #Remove this line if you want to run the code

#These are the parameters for calling our code.
#our_model.network_fileName = argv[1];
#our_model.out_NKDV_fileName = argv[2];
#our_model.method = atoi(argv[3]);
#our_model.lixel_reg_length = atoi(argv[4]);
#our_model.k_type = atoi(argv[5]);
#our_model.bandwidth = atof(argv[6]);

dir="./NKDV_Datasets/"
out_dir="./Results/"
k_type=2
lixel_reg_length=10
bandwidth=1000
dataset="New_York_taxi"

method=1 #RQS
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_b"$bandwidth $method $lixel_reg_length $k_type $bandwidth

method=2 #SPS
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_b"$bandwidth $method $lixel_reg_length $k_type $bandwidth

method=3 #ADA
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_b"$bandwidth $method $lixel_reg_length $k_type $bandwidth

method=6 #ADA^{B+}
./main $dir$dataset"/"$dataset $out_dir$dataset"_m"$method"_k"$k_type"_l"$lixel_reg_length"_b"$bandwidth $method $lixel_reg_length $k_type $bandwidth