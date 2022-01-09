# CentOS 7.6.1810

while getopts -l samtools flag
do
    case "${flag}" in
        samtools) 
            ;;
        \?)	# 지정이외 옵션은 이 문자로 할당된다.
            echo $@ is not valid option
            exit 0
            ;;  
    esac
done

yum install -y wget bzip2 gcc-c++ make ncurses-devel zlib-devel bzip2-devel xz-devel && \
mkdir ~/tools && \
TOOLS=~/tools && \
wget "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2" -P $TOOLS --no-check-certificate && \
tar xjf $TOOLS/samtools-1.9.tar.bz2 -C $TOOLS && \
cd $TOOLS/samtools-1.9 && \
bash ./configure && \
make && \
make install && \
rm -rf $TOOLS/samtools-1.9.tar.bz2
