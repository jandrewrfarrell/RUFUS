if [ -z "$nCPUs" ]; then
    nCPUs=$(grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
fi
pref=$(basename $0 .sh)
DIR=../bin
JF="$DIR/jellyfish"
[ -n "$VALGRIND" ] && JF="valgrind $JF"
SRCDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS/cloud/jellyfish-MODIFIED-merge
BUILDDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS/cloud/jellyfish-MODIFIED-merge

check () {
    cut -d\  -f 2 $1 | xargs md5sum | sort -k2,2 | diff -w $DIFFFLAGS $1 -
}

ENABLE_RUBY_BINDING=""
RUBY=""
ENABLE_PYTHON_BINDING=""
PYTHON=""
ENABLE_PERL_BINDING=""
PERL=""

if [ -n "$DEBUG" ]; then
    set -x;
    DIFFFLAGS="-y"
fi

set -e
