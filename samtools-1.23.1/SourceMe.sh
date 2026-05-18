prefix=/home/liaojing/tttt/samtools-1.23.1
export PATH=$prefix/bin:$PATH
export MANPATH=$prefix/man:$MANPATH
export LD_LIBRARY_PATH=$prefix/lib:$prefix/lib64:$LD_LIBRARY_PATH
export CFLAGS="-I$prefix/include"
export LDFLAGS="-L$prefix/lib"


export C_INCLUDE_PATH=$C_INCLUDE_PATH:$prefix/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$prefix/include

