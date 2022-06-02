# my_RNAAssessment

## Dependencies

You can see all dependencies (tools or packages used) with their corresponding versions in the [DEPENDENCY](https://github.com/FranGASTRIN/my_RNAAssessment/blob/main/DEPENDENCY) file.

## Dockerfile

**Requirement :** This Dockerfile use the image of `tensorflow/tensorflow`. Make sure to have this image before building an image with this Dockerfile :
```
sudo docker pull tensorflow/tensorflow
```

To generate the Docker image from the provided Dockerfile, use the following command line when you are on this repository:
```
sudo docker build -t my_rnaa .
```

From this image, you can now run your container using:
```
sudo docker run -it --rm my_rnaa
```

:warning: The `--rm` option causes the container and its contents to be deleted when it is stopped.

Once the container is running, 2 volumes (`data` and `results`) are created at the pathway `/var/lib/docker/volumes`. It allows the communication of files and directories between the host machine and the containe.
For example, to transfer a folder from the host machine to the container, you will need to enter in a second terminal a command like:
```
sudo cp -r /pathway/of/my/folder/ /var/lib/docker/volumes/random_volume_name/_data/.
```

*Help :* `random_volume_name` *can be find using* `sudo ls /var/lib/docker/volumes`.

Your folder will be copied in one of the 2 volume directories of your container.
You can also recover the files or folders contained in this container by copying them to one of these directories, then executing the following type of command in the second terminal :
```
sudo cp -r /var/lib/docker/volumes/random_volume_name/_data/folder /pathway/of/my/host/machine/folder/ 
```

For more information about the use of Docker :
- [How to install Docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-20-04-fr)
- [Start to use and understand Docker](https://takacsmark.com/dockerfile-tutorial-by-example-dockerfile-best-practices-2018/)

## Programs

### Downloader.py

Allows the user to download a list of PDB files, if he has one, in the format:
```
1BTKA
2F46A 
4TXRA
```
or at least:
```
1BTK
2F46
4TXR
```
using the command :
```
python3 Downloader.py your_dataset_repository_pathway/ PDB_list.txt
```
:warning: `your_dataset_repository_pathway/` must be an existing repository !
### Install angle_calculation & distance_calculation

Since they are C++ files, they must be compiled. Enter the following command, when you are in `my_RNAAssessment/`, to do it:
```
make
```
Than, run them without any argument or with the option `-h` to get their help.

## How to use
### angle_calculation

Basically, if you want to have the list of all the torsion angle values from a list of protein PDB files (Psi/Phi), you can simply run the command :
```
./angle_calculation -d dataset_pathway/ -l list_of_PDB.txt -o output_file_name
```

To apply this program to a list of RNA PDB files (eta/theta), you can use the previous command line adding to it at the end the option `-R` :
```
./angle_calculation -d dataset_pathway/ -l list_RNA_PDB.txt -o output_RNA_file -R
```
You can also access to the pseudotorsion angle values **eta'** and **theta'** using options `-Ra` or `-RA` instead of option `-R`.

There are also few more options that we leave it to the user to discover them by himself with the help option `-h`, like `-O` option to add the Omega torsion angle values for proteins, just to name one of them. 
