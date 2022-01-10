#!/bin/bash
# CentOS 7.6.1810

err_msg() { echo "$@";} >&2

usage() {
    err_msg \
"Usage: $(basename "$0")
Short options:          Description
        -s              Samtools install flag.
        -o              Install or download path.
        -h              Usage.
"
}

# Get option
while getopts "o:sh" opt; do
    case $opt in
        s)
            samtools="install"
            ;;
        o)
            outdir=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            err_msg "$@ is not valid option"
            exit 0
            ;;
    esac
done

# check option
[ $# -eq 0 ] && echo -e \
"None of option is inserted
Try 'bash ./install.sh -h' for more information." && exit 0

if  [[ $samtools -eq "install"  ]] && [[ -v outdir ]]; then
    yum install -y wget bzip2 gcc-c++ make ncurses-devel zlib-devel bzip2-devel xz-devel && \
    mkdir  && \
    TOOLS=~/tools && \
    wget "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2" -P $TOOLS --no-check-certificate && \
    tar xjf $TOOLS/samtools-1.9.tar.bz2 -C $TOOLS && \
    cd $TOOLS/samtools-1.9 && \
    bash ./configure && \
    make && \
    make install && \
    rm -rf $TOOLS/samtools-1.9.tar.bz2
fi

