# Create the 5 models
g++ train.cpp -o train
./train 10 model_init.txt seq_model_01.txt model_01.txt
./train 10 model_init.txt seq_model_02.txt model_02.txt
./train 10 model_init.txt seq_model_03.txt model_03.txt
./train 10 model_init.txt seq_model_04.txt model_04.txt
./train 10 model_init.txt seq_model_05.txt model_05.txt

# Start testing 
#g++ test.cpp -o test
#./test modellist.txt testing_data.txt result.txt

# Output accuracy 
#g++ accuracy.cpp -o accuracy
#./accuracy result.txt testing_answer.txt