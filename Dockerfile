# Requires the latest stable version of Tensorflow
FROM tensorflow/tensorflow

# Installs tools to clone Git repositories, write/read files, and 
# installs Python scikit-learn package to run python scripts
RUN apt update -y && \
    apt install vim -y && \
    apt install git -y && \
    pip3 install sklearn && \
    pip3 install seaborn

# To have a clean environment to work     
WORKDIR home

# Clone Git repository in it
RUN git clone https://github.com/FranGASTRIN/my_RNAAssessment.git

# 2 empty directories (connection between virtual machine and host)
VOLUME /home/data
VOLUME /home/results

# To make sure to have the bash command prompt
CMD ["/bin/bash"]
