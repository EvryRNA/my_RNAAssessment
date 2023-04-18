CC    = g++
SRC1  = src/script_cpp/angle_calculation.cpp
EXE1  = bin/angle_calculation
LDFLAGS=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

all: $(EXE1) 

$(EXE1): $(SRC1)
	mkdir -p bin
	#$(CC) -std=c++17 -lstdc++fs $(SRC1) -o $(EXE1)
	$(CC) -std=c++17 -L $(LDFLAGS) $(SRC1) -o $(EXE1)

docker_start:
	docker build -t my_rna_assessment .
	docker run -it my_rna_assessment

run:
	bin/angle_calculation -d decoys_data  -o output.csv -R -p -f -t
	#bin/angle_calculation -d decoys_data/3D2V.pdb  -o output.csv -R -p -f -t
	#bin/angle_calculation -d ../my_RNAAssessment_updated/data/3d_robot_dataset/1BYIA/1BYIA_decoy0_1.pdb -o output.csv -p -f -t