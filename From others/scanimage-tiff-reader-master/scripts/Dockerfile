FROM ubuntu
MAINTAINER Nathan Clack <nathan@vidriotech.com>
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install software-properties-common
RUN apt-add-repository "deb http://llvm.org/apt/trusty/ llvm-toolchain-trusty-3.8 main"
RUN add-apt-repository ppa:staticfloat/juliareleases
RUN add-apt-repository ppa:staticfloat/julia-deps
RUN add-apt-repository ppa:george-edison55/cmake-3.x
RUN apt-get -y install clang-3.8 lldb-3.8 git wget bzip2 vim julia make cmake openssh-server
RUN wget http://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh
RUN bash Anaconda2-4.0.0-Linux-x86_64.sh -bf
RUN export PATH="/root/anaconda2/bin:$PATH"
RUN mkdir /root/.ssh
RUN export CC=clang-3.8
RUN export CXX=clang++-3.8

RUN mkdir -p /var/run/sshd
RUN echo 'root:screencast' | chpasswd
RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd


RUN useradd -m -d /home/jenkins -s /bin/sh jenkins
RUN echo "jenkins:jenkins" | chpasswd

# EXPOSE 22
# CMD ["/usr/sbin/sshd", "-D"]
