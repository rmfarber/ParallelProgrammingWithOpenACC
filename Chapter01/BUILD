pgc++ -std=c++11 -acc -ta=multicore,nvidia -Minfo=accel accFill_ex1.cpp -o accFill_ex1

pgc++ -std=c++11 -acc -ta=multicore,nvidia -Minfo=accel accFill_ex2.cpp -o accFill_ex2

pgc++ -std=c++11 -acc -O3 -ta=multicore,nvidia -Minfo=accel accTask.cpp -o accTask

echo "Building accCounter.acc"
pgc++ -DUSE_ATOMIC -std=c++11 -acc -ta=multicore,nvidia:cc35 -Minfo=accel accCounter.cpp -o accCounter_w_atomic

echo "Building accCounter.seq"
pgc++ -std=c++11 -acc -ta=multicore,nvidia:cc35 -Minfo=accel accCounter.cpp -o accCounter_no_atomic

echo "Building accParaCounter
pgc++ -std=c++11 -acc -ta=multicore,nvidia:cc35 -Minfo=accel accParaCounter.cpp -o accParaCounter

echo "Building accParaRNG
pgc++ -std=c++11 -acc -ta=multicore,nvidia:cc35 -Minfo=accel accParaRNG.cpp -o accParaRNG
