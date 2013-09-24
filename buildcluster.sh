
case $1 in
	debug)
		echo "debug version"
		gcc -lm -g -Wall -D DEBUG -std=c99 -pedantic cluster.c -o cluster2
		;;
	*)
		echo "optimised version"
		gcc -lm -o3 -Wall -std=c99 -pedantic cluster.c -o cluster2
		;;
esac



