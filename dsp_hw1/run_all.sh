# Create the 5 models
#g++ train.cpp -o train
#./train 1 model_init.txt seq_model_01.txt model_01.txt
#./train 1 model_init.txt seq_model_02.txt model_02.txt
#./train 1 model_init.txt seq_model_03.txt model_03.txt
#./train 1 model_init.txt seq_model_04.txt model_04.txt
#./train 1 model_init.txt seq_model_05.txt model_05.txt

# Start testing 
#g++ test.cpp -o test
./test modellist.txt testing_data1.txt result1.txt
./test modellist.txt testing_data2.txt result2.txt

# Output accuracy 
#g++ accuracy.cpp -o accuracy
./accuracy result1.txt testing_answer.txt acc.txt