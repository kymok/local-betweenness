# usage: compute_all.sh input_dir output_dir

function julia () { /Applications/Julia-1.3.app/Contents/Resources/julia/bin/julia "$@"; }

BASEDIR=$(dirname "$0")

if [ ! -d "$1" ]; then
	echo "\"$1\" is not a directory"
	exit 1
fi

if [ ! -d "$2" ]; then
	echo "\"$2\" is not a directory"
	exit 1
fi

for file in $1/*.graphml;do
	julia $BASEDIR/betweenness.jl --input $file --output $2
done
