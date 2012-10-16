#! /bin/sh
# wrapper script for allpathslg PrepareAllPathsInputs
version=VERSION
pkg=allpathslg-$version
path=/usr/lib/$pkg
prog=PrepareAllPathsInputs.pl

# see how we were called
#prog=`basename $0`
#if [ $? != 0 -o -z "$prog" ]; then
#    echo "$pkg: unable to determine program name: $0"
#    exit 2
#fi

# see if we are on a 64-bit machine
#if [ "`uname -m`" = x86_64 ]; then
#    path="$path-64"
#fi

# make sure executable exists under path
if [ ! -x "$path/$prog" ]; then
    echo "$pkg: program does not exist under $path: $prog"
    exit 2
fi

# set up environment
PATH=$path:$PATH
export PATH

# execute program (should be in path)
exec $prog "$@"
