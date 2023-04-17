CC    = g++
SRC1  = src/script_cpp/angle_calculation.cpp
EXE1  = bin/angle_calculation

all: $(EXE1) 

$(EXE1): $(SRC1)
	mkdir -p bin
	$(CC) -std=c++17 -lstdc++fs $(SRC1) -o $(EXE1)

docker_start:
	docker build -t my_rna_assessment .
	docker run -it my_rna_assessment

run:
	bin/angle_calculation -d decoys_data  -o output.csv -R -p -f -t
